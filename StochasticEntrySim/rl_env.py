"""
Phase 1 -- OrionEntryEnv: a Gymnasium environment around the existing Orion
entry simulator, for training and benchmarking an RL surrogate controller
against the traditional predictor-corrector.

Design (see docs/rl_mdp_spec.md):
  * Controlled segment: entry interface -> drogue deploy (~7.6 km).
  * Action: 1-D continuous bank command, a in [-1, 1] -> sigma_cmd = a * pi/2.
  * Observation: ~13 normalized features (altitude, velocity, gamma, heading,
    range-to-go, errors, energy, current bank, heat/g margins, hdot).
  * Reward: terminal-dominated (miss at drogue) + dense heat/g/fuel shaping.
  * Transition: the real physics -- bank actuator -> roll PD -> 12-jet RCS ->
    roll dynamics -> 3DOF EOM -- via step_closed_loop_milestone1. The ONLY
    substitution vs run_sim.py is the guidance decision.

Two controller modes:
  * "policy"   (default): the env's action drives the bank (for RL).
  * "baseline": the real continuous predictor-corrector drives the bank
    (action ignored). Used to VALIDATE the env reproduces run_sim.py and to
    characterize the classical baseline (Phase 3).

The gymnasium import is optional so the file is usable locally for validation;
in Colab, `pip install gymnasium` gives the real base class + spaces.
"""

from __future__ import annotations

import math
from typing import Any, Dict, Optional, Tuple

import numpy as np

import constants
import mission_config
import telemetry
import AtmosphereModel
from control import (
    SimpleBankGuidance,
    CapsuleControlStack,
    CapsuleControlConfig,
    GuidanceScheduler,
    BasicObservationProvider,
    BankActuator,
    BankActuatorLimits,
    RollTorqueController,
    ReentryState,
)
from ReactionControl import build_orion_cm_rcs_12
from math_3d import (
    nominal_aero_forces_from_state,
    nominal_heating_envelope_from_state,
    IntervalSupervisorConfig,
    annotate_nominal_state_with_interval_supervisor,
)
from point_math_3d import make_initial_capsule_attitude, step_closed_loop_milestone1, aero_forces


# --- Optional gymnasium (fallback shim keeps the env usable without it) ------
try:
    import gymnasium as gym
    from gymnasium import spaces
    _HAVE_GYM = True
except Exception:  # pragma: no cover - exercised only when gym is absent
    _HAVE_GYM = False

    class _Box:
        def __init__(self, low, high, shape, dtype=np.float32):
            self.low = np.asarray(low, dtype=dtype) if np.ndim(low) else np.full(shape, low, dtype=dtype)
            self.high = np.asarray(high, dtype=dtype) if np.ndim(high) else np.full(shape, high, dtype=dtype)
            self.shape = tuple(shape)
            self.dtype = dtype

        def sample(self):
            return np.random.uniform(self.low, self.high).astype(self.dtype)

    class _Spaces:
        Box = _Box

    spaces = _Spaces()  # type: ignore

    class _EnvBase:
        metadata: Dict[str, Any] = {}

    class gym:  # type: ignore
        Env = _EnvBase


# Observation feature order (kept explicit for traceability / logging).
OBS_NAMES = (
    "altitude_norm", "velocity_norm", "gamma_norm", "sin_chi", "cos_chi",
    "range_to_go_norm", "heading_error_norm", "cross_track_norm",
    "specific_energy_norm", "sigma_norm", "heat_margin", "g_margin", "hdot_norm",
)

# Extra interval (worst-case) features appended when interval_obs is enabled.
# These come from the interval supervisor's guaranteed bounds, so the policy
# can see how much uncertainty / worst-case margin it is carrying.
INTERVAL_OBS_NAMES = (
    "heat_margin_hi",   # (q_lim - qdot_hi) / q_lim  -- worst-case heat margin
    "g_margin_hi",      # (n_lim - n_hi)   / n_lim   -- worst-case g margin
    "box_width_norm",   # scalar summary of the state-box uncertainty
)

# Normalization scales (denominators).
_NORM = {
    "alt": 120_000.0, "vel": 11_000.0, "gamma": math.pi / 2.0,
    "range": 2_500_000.0, "heading": math.pi, "cross": 500_000.0,
    "energy": 6.0e7, "sigma": math.pi / 2.0, "hdot": 300.0,
}


class PolicyGuidance:
    """Guidance stand-in whose bank command is set externally (the RL action).

    Implements the minimal surface CapsuleControlStack needs: a `params`
    attribute, compute_sigma_cmd(), reset(), get_last_debug().
    """

    def __init__(self, params: Optional[dict] = None,
                 target_phi_rad: float = 0.0, target_lam_rad: float = 0.0):
        self.params = params
        self.target_phi_rad = float(target_phi_rad)
        self.target_lam_rad = float(target_lam_rad)
        self._action_sigma_rad = 0.0

    def set_action(self, sigma_rad: float) -> None:
        self._action_sigma_rad = float(sigma_rad)

    def compute_sigma_cmd(self, t_s, state, obs) -> float:
        return float(self._action_sigma_rad)

    def reset(self) -> None:
        self._action_sigma_rad = 0.0

    def get_last_debug(self) -> Dict[str, Any]:
        return {"predictor_mode": "rl_policy",
                "chosen_sigma_cmd_deg": math.degrees(self._action_sigma_rad)}


class OrionEntryEnv(gym.Env):
    """Gymnasium env for Orion entry guidance (entry -> drogue deploy)."""

    metadata = {"render_modes": []}

    def __init__(
        self,
        config_path: str = "configs/default.json",
        controller_mode: str = "policy",      # "policy" | "baseline"
        max_episode_steps: int = 1400,
        gload_limit_g: float = 8.0,
        target_radius_km: float = 25.0,
        reward_weights: Optional[telemetry.RewardWeights] = None,
        fast_atmosphere: bool = True,
        dispersion: Optional[bool] = None,
        w_progress_per_km: float = 0.05,
        interval_obs: bool = False,
        interval_reward: bool = False,
        interval_shield: bool = False,
        shield_horizon: int = 1,
        seed: Optional[int] = None,
    ):
        super().__init__()
        if controller_mode not in ("policy", "baseline"):
            raise ValueError("controller_mode must be 'policy' or 'baseline'")
        # --- interval-arithmetic integration (3 independent, ablatable knobs) ---
        #   interval_obs    : append worst-case interval features to the obs
        #   interval_reward : penalize the interval UPPER BOUND on heat / g
        #   interval_shield : veto unsafe bank commands via 1-step reachability
        # Any of them turns on the per-step interval propagation (use_intervals).
        self.interval_obs = bool(interval_obs)
        self.interval_reward = bool(interval_reward)
        self.interval_shield = bool(interval_shield)
        self.shield_horizon = max(1, int(shield_horizon))
        self.use_intervals = self.interval_obs or self.interval_reward or self.interval_shield
        # Dense progress shaping: reward per km of range-to-target closed each
        # step. Turns the otherwise-sparse (terminal-only) reward into a dense
        # "warmer/colder" signal -- without it PPO can't credit-assign a single
        # end-of-episode miss across the ~540 bank decisions in an episode.
        self.w_progress_per_km = float(w_progress_per_km)

        # Fast tabulated atmosphere for RL throughput (~100x fewer PyMSIS
        # calls). Built once, module-global, shared across all env instances.
        # Disable for final high-fidelity eval to recover lat/lon variation.
        if fast_atmosphere and not AtmosphereModel.atmosphere_table_active():
            AtmosphereModel.build_atmosphere_table()
        self.controller_mode = controller_mode
        self.max_episode_steps = int(max_episode_steps)
        self.gload_limit_g = float(gload_limit_g)
        self.target_radius_km = float(target_radius_km)
        self.reward_weights = reward_weights or telemetry.RewardWeights()

        # --- config + fixed parameters ---
        self.cfg = mission_config.load_config(config_path)
        self.params = {
            "mass_kg": float(constants.CAPSULE_MASS_KG),
            "ref_area_m2": float(constants.CAPSULE_REFERENCE_AREA_M2),
            "aero_model": "orion_cm_trim",
            "CD_const": 1.15, "CL_const": 0.28,
        }
        mission_config.apply_aero_to_params(self.cfg, self.params)
        self._radius_earth = float(constants.RADIUS_EARTH)
        # Nominal IC (config units) + entry-interface state vector.
        self._ic_nominal = mission_config.nominal_initial_conditions(self.cfg)
        self._x0 = mission_config.initial_state_vector(self.cfg, self._radius_earth)
        # IC dispersion: explicit env flag overrides the config's `enabled`.
        self.dispersion = (mission_config.dispersion_enabled(self.cfg)
                           if dispersion is None else bool(dispersion))
        self._ic_drawn = dict(self._ic_nominal)   # IC actually used this episode
        tgt = mission_config.mission_targets(self.cfg)
        self.target_phi_rad = float(tgt["target_phi_rad"])
        self.target_lam_rad = float(tgt["target_lam_rad"])
        self.cos_gamma_gate = float(tgt["cos_gamma_termination_gate"])
        self.Izz = float(constants.CAPSULE_IZZ_KGM2)
        self.dt_s = float(constants.ENTRY_DT_S)
        self.guidance_period_s = float(constants.GUIDANCE_PERIOD_S)
        self.n_substeps = max(1, int(round(self.guidance_period_s / self.dt_s)))
        # Drogue-deploy altitude = end of the controlled segment.
        self.drogue_alt_m = float(mission_config.build_cpas_config(self.cfg).drogue_deploy_alt_m)
        self.heat_rate_limit = float(constants.HEAT_RATE_LIMIT_DEFAULT)

        # Interval supervisor config (only used when use_intervals). Mirrors the
        # one run_sim.py builds so the RL interval bounds match the classical
        # controller's. Candidate bank set the shield can fall back to.
        self.supervisor_cfg = self._build_supervisor_cfg() if self.use_intervals else None
        self._shield_candidates_rad = sorted(
            {math.radians(d) for d in (0.0, 25.0, 45.0, 70.0, 90.0)}
            | {-math.radians(d) for d in (25.0, 45.0, 70.0, 90.0)})

        # --- spaces (obs grows when interval features are appended) ---
        self.obs_names = OBS_NAMES + (INTERVAL_OBS_NAMES if self.interval_obs else ())
        self.action_space = spaces.Box(-1.0, 1.0, shape=(1,), dtype=np.float32)
        self.observation_space = spaces.Box(
            -np.inf, np.inf, shape=(len(self.obs_names),), dtype=np.float32)

        # --- build the control stack (guidance depends on mode) ---
        self._build_control_stack()

        self._rng = np.random.default_rng(seed)
        self._reset_episode_state()

    # ------------------------------------------------------------------
    def _build_control_stack(self) -> None:
        if self.controller_mode == "baseline":
            self.guidance = SimpleBankGuidance(
                target_phi_rad=self.target_phi_rad,
                target_lam_rad=self.target_lam_rad,
                params=dict(self.params),
                velocity_enable_mps=40_000.0,
                use_continuous_predictor=True,
                use_bank_reversal=True,
            )
        else:
            self.guidance = PolicyGuidance(
                params=dict(self.params),
                target_phi_rad=self.target_phi_rad,
                target_lam_rad=self.target_lam_rad,
            )

        self.obs_provider = BasicObservationProvider(
            target_phi_rad=self.target_phi_rad,
            target_lam_rad=self.target_lam_rad,
            params=self.params,
        )
        self.bank_actuator = BankActuator(BankActuatorLimits(
            sigma_rate_max_rps=math.radians(20.0),
            sigma_accel_max_rps2=math.radians(40.0),
        ))
        self.roll_controller = RollTorqueController(
            kp_Nm_per_rad=float(constants.ROLL_KP_NM_PER_RAD),
            kd_Nm_per_rad_s=float(constants.ROLL_KD_NM_PER_RAD_S),
            max_torque_Nm=float(constants.ROLL_CMD_MAX_TORQUE_NM),
            sigma_deadband_rad=float(constants.ROLL_SIGMA_DEADBAND_RAD),
            rate_deadband_rad_s=float(constants.ROLL_RATE_DEADBAND_RAD_S),
        )
        self.control = CapsuleControlStack(
            cfg=CapsuleControlConfig(guidance_period_s=self.guidance_period_s),
            scheduler=GuidanceScheduler(period_s=self.guidance_period_s),
            guidance=self.guidance,
            obs_provider=self.obs_provider,
            bank_actuator=self.bank_actuator,
            roll_controller=self.roll_controller,
        )
        self.rcs = build_orion_cm_rcs_12()

    def _control_step_fn(self, t_s, x_state, sigma_actual_rad, roll_rate_rad_s):
        r_m, phi, lam, V, gamma, chi = [float(v) for v in x_state]
        state = ReentryState(
            r_m=r_m, phi_rad=phi, lam_rad=lam, V_mps=V, gamma_rad=gamma, chi_rad=chi,
            sigma_actual_rad=float(sigma_actual_rad), roll_rate_rad_s=float(roll_rate_rad_s),
            sigma_cmd_rad=0.0, sigma_target_rad=0.0,
        )
        Dr, Dth, Dph, Lr, Lth, Lph, _, _ = aero_forces(
            r=r_m, V=V, gamma=gamma, chi=chi, params=self.params,
            sigma_actual_rad=float(sigma_actual_rad),
        )
        drag = math.sqrt(Dr * Dr + Dth * Dth + Dph * Dph)
        lift = math.sqrt(Lr * Lr + Lth * Lth + Lph * Lph)
        return self.control.step(
            t_s=float(t_s), dt_s=float(self.dt_s), state=state,
            lift_N=float(lift), drag_N=float(drag), mass_kg=float(self.params["mass_kg"]),
        )

    def _reset_episode_state(self) -> None:
        self.x = list(self._x0)
        self.att = make_initial_capsule_attitude(0.0)
        self.heat_shield = None
        self.t_s = 0.0
        self.step_count = 0
        self.rcs_fuel_kg = 0.0
        self.peak_g = 0.0
        self.peak_qdot = 0.0
        self._last_qdot = 0.0
        self._last_g = 0.0
        self._q_nom = 0.0           # nominal dynamic pressure (for the g bound)
        # --- interval state (worst-case tube), maintained when use_intervals ---
        self.x_interval = None              # current state box (None -> inflate)
        self.interval_heat_shield = None    # interval heat-shield accumulator
        self.interval_active = self.use_intervals
        self._qdot_hi = 0.0         # worst-case heat rate this step (W/m^2)
        self._g_hi = 0.0            # worst-case load factor this step (g)
        self._box_width = 0.0       # scalar uncertainty summary
        self.peak_qdot_hi = 0.0
        self.peak_g_hi = 0.0
        self.interval_violations = 0   # steps whose worst-case breached a limit
        self.shield_interventions = 0  # steps where the shield overrode the policy
        self.shield_steps = 0          # shield decision points evaluated
        # Range-to-target at the start of the episode, for progress shaping.
        self._prev_range_km = self._great_circle_miss_km()

    # ------------------------------------------------------------------
    def reset(self, *, seed: Optional[int] = None, options: Optional[dict] = None):
        if seed is not None:
            self._rng = np.random.default_rng(seed)
        # Resolve the initial condition for this episode (config units), in
        # priority order:
        #   1. options["ic_override"] -> exact IC (frozen test-set replay)
        #   2. dispersion on          -> sample around nominal via the RNG
        #   3. otherwise              -> nominal IC
        options = options or {}
        use_dispersion = options.get("dispersion", self.dispersion)
        if options.get("ic_override"):
            self._ic_drawn = {**self._ic_nominal, **dict(options["ic_override"])}
        elif use_dispersion:
            self._ic_drawn = mission_config.sample_dispersed_ic(self.cfg, self._rng)
        else:
            self._ic_drawn = dict(self._ic_nominal)
        self._x0 = mission_config.initial_state_vector_from_ic(
            self._ic_drawn, self._radius_earth)

        self.control.reset()
        self.rcs.reset()
        self._reset_episode_state()
        # Prime heat/g from the initial state so the first obs has margins.
        self._update_aux_metrics()
        return self._build_obs(), self._info(terminal=False)

    def _alt(self) -> float:
        return float(self.x[0] - constants.RADIUS_EARTH)

    def _update_aux_metrics(self) -> None:
        """Refresh heat shield, instantaneous heat rate, and g-load."""
        heat = nominal_heating_envelope_from_state(
            x_nominal=list(self.x), dt_s=float(self.dt_s), heat_shield=self.heat_shield)
        self.heat_shield = heat["heat_shield"]
        self._last_qdot = float(heat["qdot_total_stag"].hi)
        aero = nominal_aero_forces_from_state(
            x=self.x, sigma_rad=float(self.att.sigma_rel_rad), params=self.params)
        rho = float(aero["rho_kgm3"])
        V = float(self.x[3])
        self._q_nom = 0.5 * rho * V * V        # nominal dynamic pressure (Pa)
        d = telemetry.compute_derived_metrics(
            r_m=float(self.x[0]), V_mps=V, altitude_m=self._alt(),
            T_K=None, rho_kgm3=rho,
            drag_mag_N=float(aero["drag_mag_N"]), lift_mag_N=float(aero["lift_mag_N"]),
            mass_kg=float(self.params["mass_kg"]))
        self._last_g = float(d["load_factor_g"])
        self.peak_g = max(self.peak_g, self._last_g)
        self.peak_qdot = max(self.peak_qdot, self._last_qdot)

    # ================= interval-arithmetic integration ===================
    def _build_supervisor_cfg(self) -> IntervalSupervisorConfig:
        """Interval supervisor config, mirroring run_sim.py so the RL worst-case
        bounds are computed exactly like the classical controller's."""
        return IntervalSupervisorConfig(
            r_half_width_m=10.0,
            phi_half_width_rad=math.radians(0.01),
            lam_half_width_rad=math.radians(0.01),
            V_half_width_mps=20.0,
            gamma_half_width_rad=math.radians(0.15),
            chi_half_width_rad=math.radians(0.20),
            min_altitude_m=0.0,
            max_altitude_m=130_000.0,
            max_speed_mps=80_000.0,
            max_dynamic_pressure_pa=5.0e7,
            include_heating=True,
            heat_rate_limit=float(constants.HEAT_RATE_LIMIT_DEFAULT),
            heat_load_limit=float(constants.HEAT_LOAD_LIMIT_DEFAULT),
            interval_recenter_enabled=bool(constants.INTERVAL_RECENTER_ENABLED),
            interval_recenter_cadence_s=float(constants.INTERVAL_RECENTER_CADENCE_S),
            interval_recenter_width_thresholds=dict(constants.INTERVAL_RECENTER_WIDTH_THRESHOLDS),
            interval_box_split_enabled=bool(constants.INTERVAL_BOX_SPLIT_ENABLED),
            interval_box_split_max_depth=int(constants.INTERVAL_BOX_SPLIT_MAX_DEPTH),
            interval_box_split_width_thresholds=dict(constants.INTERVAL_BOX_SPLIT_WIDTH_THRESHOLDS),
            interval_denominator_safety_V_mps=float(constants.INTERVAL_DENOMINATOR_SAFETY_V_MPS),
            interval_denominator_safety_cos_gamma=float(constants.INTERVAL_DENOMINATOR_SAFETY_COS_GAMMA),
            interval_denominator_safety_cos_phi=float(constants.INTERVAL_DENOMINATOR_SAFETY_COS_PHI),
        )

    def _annotate(self, x_nominal, sigma_rad, x_interval, heat_shield):
        """One interval-propagation step. Returns the IntervalAnnotationResult."""
        return annotate_nominal_state_with_interval_supervisor(
            x_nominal_old=list(x_nominal), params=self.params,
            sigma_actual_after_rad=float(sigma_rad), x_interval_old=x_interval,
            supervisor_cfg=self.supervisor_cfg, dt_s=float(self.dt_s),
            heat_shield=heat_shield, t_s=float(self.t_s))

    def _metrics_from_annotation(self, ann) -> Tuple[float, float, float]:
        """Worst-case (qdot_hi [W/m^2], g_hi [g], box_width_scalar) from one
        annotation. The heat bound is rigorous (supervisor); the g bound scales
        the nominal load factor by the dynamic-pressure upper bound, assuming the
        aero coefficient is ~constant over the (small) state box."""
        qdot_hi = (float(ann.heating_qdot_max_interval.hi)
                   if ann.heating_qdot_max_interval is not None else self._last_qdot)
        q_hi = float(ann.q_interval.hi) if ann.q_interval is not None else self._q_nom
        if self._q_nom > 1e-9 and self._last_g > 0.0:
            g_hi = self._last_g * (q_hi / self._q_nom)
        else:
            g_hi = self._last_g
        w = ann.state_widths_new or {}
        # Scalar uncertainty summary: velocity + altitude box widths, normalized.
        box_width = float(w.get("V", 0.0)) / 200.0 + float(w.get("r", 0.0)) / 2000.0
        return qdot_hi, max(g_hi, self._last_g), box_width

    def _interval_propagate(self, x_nominal_old, sigma_rad) -> None:
        """Advance the interval tube one substep alongside the nominal state and
        refresh the worst-case metrics. Deactivates intervals on failure."""
        if not self.use_intervals or not self.interval_active:
            return
        ann = self._annotate(x_nominal_old, sigma_rad, self.x_interval,
                             self.interval_heat_shield)
        # Only a numerical failure kills the tube. A heat/g violation is a real
        # worst-case signal we want to keep feeding to the obs / reward / shield.
        if (not ann.interval_valid) or getattr(ann, "interval_numerical_failure", False):
            self.interval_active = False
            return
        self.x_interval = [type(iv)(iv.lo, iv.hi) for iv in ann.x_interval_new]
        self.interval_heat_shield = getattr(ann, "heat_shield", self.interval_heat_shield)
        self._qdot_hi, self._g_hi, self._box_width = self._metrics_from_annotation(ann)
        self.peak_qdot_hi = max(self.peak_qdot_hi, self._qdot_hi)
        self.peak_g_hi = max(self.peak_g_hi, self._g_hi)
        if self._qdot_hi > self.heat_rate_limit or self._g_hi > self.gload_limit_g:
            self.interval_violations += 1

    def _bank_is_safe(self, sigma_rad) -> bool:
        """1-step reachability check: would commanding this bank keep the
        worst-case heat rate and g within limits over the next substep?"""
        if not self.interval_active:
            return True   # no valid tube -> cannot veto
        ann = self._annotate(self.x, sigma_rad, self.x_interval,
                             self.interval_heat_shield)
        if not ann.interval_valid:
            return False
        qdot_hi, g_hi, _ = self._metrics_from_annotation(ann)
        return qdot_hi <= self.heat_rate_limit and g_hi <= self.gload_limit_g

    def _shield(self, sigma_cmd_rad: float) -> Tuple[float, bool]:
        """Interval safety shield. If the policy's bank passes the 1-step
        reachability check, use it. Otherwise project to the nearest candidate
        bank that is safe; if none is safe, pick the candidate with the lowest
        worst-case heat (most conservative). Returns (bank, intervened)."""
        self.shield_steps += 1
        if self._bank_is_safe(sigma_cmd_rad):
            return sigma_cmd_rad, False
        # Search candidates ordered by closeness to the policy's command.
        cands = sorted(self._shield_candidates_rad, key=lambda c: abs(c - sigma_cmd_rad))
        for c in cands:
            if self._bank_is_safe(c):
                self.shield_interventions += 1
                return c, True
        # Nothing certified safe -> minimize worst-case heat as a fallback.
        best_c, best_q = sigma_cmd_rad, float("inf")
        for c in cands:
            if not self.interval_active:
                break
            ann = self._annotate(self.x, c, self.x_interval, self.interval_heat_shield)
            q_hi, _, _ = self._metrics_from_annotation(ann)
            if q_hi < best_q:
                best_q, best_c = q_hi, c
        self.shield_interventions += 1
        return best_c, True

    def _build_obs(self) -> np.ndarray:
        r, phi, lam, V, gamma, chi = [float(v) for v in self.x]
        obs_d = self.obs_provider.observe(
            ReentryState(r_m=r, phi_rad=phi, lam_rad=lam, V_mps=V, gamma_rad=gamma,
                         chi_rad=chi, sigma_actual_rad=float(self.att.sigma_rel_rad),
                         roll_rate_rad_s=float(self.att.omega_b_rad_s[2]),
                         sigma_cmd_rad=0.0, sigma_target_rad=0.0), self.t_s)
        alt = self._alt()
        spec_e = 0.5 * V * V + constants.STANDARD_GRAVITY_MPS2 * alt
        hdot = V * math.sin(gamma)
        feats = np.array([
            alt / _NORM["alt"],
            V / _NORM["vel"],
            gamma / _NORM["gamma"],
            math.sin(chi), math.cos(chi),
            float(obs_d.get("range_to_go_m", 0.0)) / _NORM["range"],
            float(obs_d.get("heading_error_rad", 0.0)) / _NORM["heading"],
            float(obs_d.get("cross_track_error_m", 0.0)) / _NORM["cross"],
            spec_e / _NORM["energy"],
            float(self.att.sigma_rel_rad) / _NORM["sigma"],
            (self.heat_rate_limit - self._last_qdot) / self.heat_rate_limit,
            (self.gload_limit_g - self._last_g) / self.gload_limit_g,
            hdot / _NORM["hdot"],
        ], dtype=np.float32)
        if self.interval_obs:
            interval_feats = np.array([
                (self.heat_rate_limit - self._qdot_hi) / self.heat_rate_limit,
                (self.gload_limit_g - self._g_hi) / self.gload_limit_g,
                min(self._box_width, 5.0),   # clipped scalar uncertainty
            ], dtype=np.float32)
            feats = np.concatenate([feats, interval_feats])
        return feats

    # ------------------------------------------------------------------
    def step(self, action):
        self._intervened_this_step = False
        if self.controller_mode == "policy":
            a = float(np.asarray(action).reshape(-1)[0])
            a = max(-1.0, min(1.0, a))
            sigma_cmd = a * (math.pi / 2.0)
            # Interval safety shield: veto/clip an unsafe bank before it flies.
            if self.interval_shield:
                sigma_cmd, self._intervened_this_step = self._shield(sigma_cmd)
            self.guidance.set_action(sigma_cmd)
        # baseline mode: action ignored; the predictor computes sigma.

        terminated = False
        truncated = False
        outcome = "running"
        range_to_go_m = 0.0

        for _ in range(self.n_substeps):
            x_old = list(self.x)
            sr = step_closed_loop_milestone1(
                t_s=float(self.t_s), x_trans=list(self.x), att=self.att,
                params=self.params, control_step_fn=self._control_step_fn,
                rcs=self.rcs, dt_s=float(self.dt_s), Izz_kgm2=float(self.Izz))
            self.x = list(sr.x_new)
            self.att = sr.att_new
            self.t_s += self.dt_s
            self.step_count_substep = getattr(self, "step_count_substep", 0) + 1

            # RCS fuel for this substep
            fired = int(sr.roll_step.num_fired_internal_steps) * float(constants.RCS_INTERNAL_DT_S)
            tsec = float(len(sr.roll_step.fire_cmd.active_names())) * fired
            self.rcs_fuel_kg += telemetry.rcs_propellant_kg(tsec)

            self._update_aux_metrics()
            # Advance the interval tube alongside the nominal state.
            self._interval_propagate(x_old, float(sr.sigma_actual_after_rad))

            alt = self._alt()
            if alt <= self.drogue_alt_m:
                terminated = True; outcome = "drogue"; break
            if abs(math.cos(self.x[4])) < self.cos_gamma_gate:
                terminated = True; outcome = "cos_gamma_fail"; break
            if (not all(math.isfinite(v) for v in self.x)) or self.peak_g > 50.0:
                terminated = True; outcome = "diverged"; break

        self.step_count += 1
        if not terminated and self.step_count >= self.max_episode_steps:
            truncated = True; outcome = "timeout"

        reward, miss_km = self._reward(terminated, outcome)
        info = self._info(terminal=terminated or truncated, outcome=outcome, miss_km=miss_km)
        return self._build_obs(), float(reward), bool(terminated), bool(truncated), info

    # ------------------------------------------------------------------
    def _great_circle_miss_km(self) -> float:
        R = constants.RADIUS_EARTH / 1000.0
        p1, l1 = float(self.x[1]), float(self.x[2])
        p2, l2 = self.target_phi_rad, self.target_lam_rad
        dp, dl = p2 - p1, l2 - l1
        a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
        return float(2 * R * math.asin(min(1.0, math.sqrt(a))))

    def _reward(self, terminated: bool, outcome: str) -> Tuple[float, float]:
        miss_km = self._great_circle_miss_km()

        # Dense per-step shaping (heat / g penalties -- only bite near limits).
        # interval_reward penalizes the worst-case (interval upper bound) instead
        # of the nominal point, so robustness to bounded uncertainty is learned.
        qdot_pen = self._qdot_hi if (self.interval_reward and self.use_intervals) else self._last_qdot
        g_pen = self._g_hi if (self.interval_reward and self.use_intervals) else self._last_g
        ing = telemetry.composite_reward(
            range_to_go_m=miss_km * 1000.0,
            qdot_w_m2=qdot_pen, load_factor_g=g_pen,
            rcs_fuel_rate_kg_s=0.0, weights=self.reward_weights)
        r = ing["rew_heat_term"] + ing["rew_gload_term"]

        # Dense progress shaping: reward range-to-target closed since last step
        # (positive when getting closer, negative when drifting away). This is
        # the signal that actually drives targeting; without it the reward is
        # effectively terminal-only and PPO stalls at a large miss.
        r += self.w_progress_per_km * (self._prev_range_km - miss_km)
        self._prev_range_km = miss_km

        if terminated and outcome == "drogue":
            w = self.reward_weights
            r += -w.w_range * (miss_km * 1000.0)              # terminal miss
            if miss_km < self.target_radius_km:
                r += 50.0                                     # on-target bonus
            if self.peak_g < self.gload_limit_g:
                r += 20.0                                     # survivability bonus
        elif terminated and outcome in ("cos_gamma_fail", "diverged"):
            r += -500.0                                       # hard failure
        return float(r), float(miss_km)

    def _info(self, terminal: bool, outcome: str = "running", miss_km: float = float("nan")) -> Dict[str, Any]:
        return {
            "t_s": float(self.t_s), "altitude_m": self._alt(),
            "V_mps": float(self.x[3]), "miss_km": float(miss_km),
            "peak_g": float(self.peak_g), "peak_qdot_MWm2": float(self.peak_qdot) / 1e6,
            "rcs_fuel_kg": float(self.rcs_fuel_kg), "outcome": outcome,
            "success": bool(terminal and outcome == "drogue" and
                            self._great_circle_miss_km() < self.target_radius_km),
            "controller_mode": self.controller_mode,
            "dispersed": bool(self.dispersion),
            "ic": dict(self._ic_drawn),
            # --- interval diagnostics (meaningful only when use_intervals) ---
            "use_intervals": bool(self.use_intervals),
            "interval_active": bool(self.interval_active),
            "peak_qdot_hi_MWm2": float(self.peak_qdot_hi) / 1e6,
            "peak_g_hi": float(self.peak_g_hi),
            "interval_violations": int(self.interval_violations),
            "shield_interventions": int(self.shield_interventions),
            "shield_steps": int(self.shield_steps),
            "shield_intervention_rate": (float(self.shield_interventions) / self.shield_steps
                                         if self.shield_steps > 0 else 0.0),
        }


# ============================================================================
# Validation: run the env in BASELINE mode and confirm it reproduces run_sim.py
# (the predictor-corrector should land near the target with survivable entry g).
# ============================================================================
if __name__ == "__main__":
    import sys
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else "configs/default.json"
    print(f"[validate] OrionEntryEnv baseline mode on {cfg_path}")
    env = OrionEntryEnv(config_path=cfg_path, controller_mode="baseline")
    obs, info = env.reset(seed=42)
    print(f"  obs dim = {obs.shape[0]} ({', '.join(OBS_NAMES)})")
    done = False
    steps = 0
    while not done:
        obs, rew, term, trunc, info = env.step(env.action_space.sample())  # ignored in baseline
        done = term or trunc
        steps += 1
        if steps % 100 == 0:
            print(f"  t={info['t_s']:.0f}s alt={info['altitude_m']/1000:.1f}km "
                  f"V={info['V_mps']:.0f} peak_g={info['peak_g']:.1f}")
    print("-" * 60)
    print(f"  outcome      : {info['outcome']}")
    print(f"  MISS         : {info['miss_km']:.1f} km   (run_sim.py baseline ~4.5 km)")
    print(f"  entry peak g : {info['peak_g']:.1f} g")
    print(f"  peak qdot    : {info['peak_qdot_MWm2']:.2f} MW/m^2")
    print(f"  RCS fuel     : {info['rcs_fuel_kg']:.1f} kg")
    print(f"  guidance steps: {steps}")
