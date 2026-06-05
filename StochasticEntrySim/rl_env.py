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
from math_3d import nominal_aero_forces_from_state, nominal_heating_envelope_from_state
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
        seed: Optional[int] = None,
    ):
        super().__init__()
        if controller_mode not in ("policy", "baseline"):
            raise ValueError("controller_mode must be 'policy' or 'baseline'")

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
        self._x0 = mission_config.initial_state_vector(self.cfg, float(constants.RADIUS_EARTH))
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

        # --- spaces ---
        self.action_space = spaces.Box(-1.0, 1.0, shape=(1,), dtype=np.float32)
        self.observation_space = spaces.Box(-np.inf, np.inf, shape=(len(OBS_NAMES),), dtype=np.float32)

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

    # ------------------------------------------------------------------
    def reset(self, *, seed: Optional[int] = None, options: Optional[dict] = None):
        if seed is not None:
            self._rng = np.random.default_rng(seed)
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
        d = telemetry.compute_derived_metrics(
            r_m=float(self.x[0]), V_mps=float(self.x[3]), altitude_m=self._alt(),
            T_K=None, rho_kgm3=float(aero["rho_kgm3"]),
            drag_mag_N=float(aero["drag_mag_N"]), lift_mag_N=float(aero["lift_mag_N"]),
            mass_kg=float(self.params["mass_kg"]))
        self._last_g = float(d["load_factor_g"])
        self.peak_g = max(self.peak_g, self._last_g)
        self.peak_qdot = max(self.peak_qdot, self._last_qdot)

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
        return feats

    # ------------------------------------------------------------------
    def step(self, action):
        if self.controller_mode == "policy":
            a = float(np.asarray(action).reshape(-1)[0])
            a = max(-1.0, min(1.0, a))
            self.guidance.set_action(a * (math.pi / 2.0))
        # baseline mode: action ignored; the predictor computes sigma.

        terminated = False
        truncated = False
        outcome = "running"
        range_to_go_m = 0.0

        for _ in range(self.n_substeps):
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
        # Dense per-step shaping (heat / g / fuel).
        ing = telemetry.composite_reward(
            range_to_go_m=self._great_circle_miss_km() * 1000.0,
            qdot_w_m2=self._last_qdot, load_factor_g=self._last_g,
            rcs_fuel_rate_kg_s=0.0, weights=self.reward_weights)
        r = ing["rew_heat_term"] + ing["rew_gload_term"]

        miss_km = self._great_circle_miss_km()
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
