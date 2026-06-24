"""
Orion Capsule Parachute Assembly System (CPAS), a high fidelity model.

The model captures the full deployment sequence the real Orion CM uses:

    STOWED to FBC_JETTISONED to DROGUE to PILOT to MAIN to LANDED

and inside each chute phase it models the discrete pyrotechnic reefing stages
plus the mortar and line stretch transient at deployment.

Features implemented:

  1. Forward Bay Cover jettison. Fires before the drogue and sheds 150 kg.
  2. Mortar deployment and line stretch. About 0.6 s of zero drag, then a short
     opening shock transient above the reefed drag.
  3. Reefing pyrotechnics. Drogues use 2 stages (50 percent then 100 percent),
     mains use 3 stages (13 then 35 then 100 percent).
  4. Chute pendulum swing. A self sustaining oscillator that adds horizontal
     velocity, used only when the cluster is degraded.
  5. Asymmetric inflation, also called squidding. Each operational main has a
     random chance per step to briefly collapse to about 20 percent area.
  6. Wind drift under chutes. A configurable east and north wind applied to the
     position update.
  7. A 2 of 3 main failure mode. Set by num_mains_operational, or by a random
     failure trigger.
  8. RCS roll for splashdown. Exposes a target bank during the drogue and main
     phases so the existing roll controller keeps the capsule oriented for
     impact.

Integration contract:

Each simulation step the caller invokes CPAS.step(t_s, alt_m, V_mps, mach=None)
and applies the returned CPASOutput to its dynamics:

    drag_area_cdA_m2          add a force q_bar times drag_area_cdA_m2 opposite V
    lift_scale                multiply the aero CL by this scalar
    mass_shed_kg              subtract from params['mass_kg'] once
    wind_east_mps, north_mps  add to the position rate under chutes
    pendulum_lateral_v_mps    add to horizontal velocity
    target_sigma_rad          guidance or RCS roll target, when not NaN

The module has zero dependency on the rest of the simulator, so it stays unit
testable on its own. See verify_cpas.py.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional, Tuple


# =============================================================================
# Enums
# =============================================================================

class CPASPhase(str, Enum):
    """Discrete deployment phases of the CPAS state machine."""
    STOWED = "stowed"
    DROGUE = "drogue"
    PILOT = "pilot"
    MAIN = "main"
    LANDED = "landed"


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class CPASConfig:
    """Tunable CPAS parameters. Defaults reflect public Orion specs and
    enable every modeled feature."""

    # -------- Forward Bay Cover --------
    enable_fbc: bool = True
    fbc_jettison_alt_m: float = 8500.0       # ~28,000 ft, ~900 m above drogue
    fbc_mass_shed_kg: float = 150.0

    # -------- Drogue chutes (2) --------
    drogue_deploy_alt_m: float = 7600.0
    drogue_deploy_max_speed_mps: float = 250.0
    drogue_deploy_max_mach: float = 0.9
    drogue_cd: float = 0.50
    drogue_ref_area_m2: float = 76.0
    drogue_line_stretch_s: float = 0.6
    drogue_opening_shock_s: float = 0.3
    drogue_opening_shock_factor: float = 2.0   # multiplier on reefed drag
    drogue_reefed_area_fraction: float = 0.50
    drogue_reefed_duration_s: float = 8.0

    # -------- Pilot chutes (3) --------
    pilot_deploy_alt_m: float = 2900.0
    pilot_cd: float = 0.60
    pilot_ref_area_m2: float = 26.0
    pilot_open_time_s: float = 1.5

    # -------- Main chutes (3) --------
    main_deploy_alt_m: float = 2400.0
    main_cd: float = 0.85
    main_ref_area_m2: float = 2950.0
    main_line_stretch_s: float = 0.6
    main_opening_shock_s: float = 0.3
    main_opening_shock_factor: float = 2.5
    main_reef1_area_fraction: float = 0.13
    main_reef1_duration_s: float = 8.0
    main_reef2_area_fraction: float = 0.35
    main_reef2_duration_s: float = 8.0

    # -------- Pendulum dynamics --------
    # In real Orion, the pendulum is an *emergent failure mode* that only
    # occurs when a main chute is lost (Edwards et al. 2024). With all
    # three mains healthy the cluster is stable and there is no swing.
    # Edwards reports "pendulum limit cycle oscillations exceeding 25 deg"
    # in 2-of-3 main failure cases, which we model as a Van der Pol-style
    # self-sustaining oscillator: amplitude is pumped up below the target
    # and bled off above it, producing a steady ~25 deg swing for the
    # full main-chute descent rather than an exponential decay.
    enable_pendulum: bool = True
    pendulum_period_s: float = 8.0
    pendulum_target_amplitude_rad: float = math.radians(25.0)
    pendulum_limit_cycle_strength: float = 0.30   # Van der Pol mu (1/s)
    pendulum_riser_length_m: float = 50.0

    # -------- Squidding (asymmetric inflation) --------
    enable_squidding: bool = True
    squid_probability_per_second: float = 0.005   # per main, per second
    squid_duration_s: float = 1.5
    squid_area_fraction: float = 0.20

    # -------- Failure mode --------
    # Stochastic main-chute failure prior from:
    #   Fuqua, B. C. (2010). "Development of the Initial Main Parachute
    #   Failure Probability for the CxP Orion CEV CPAS." ENRE 655 Class
    #   Project, NASA Safety & Mission Assurance Space Shuttle &
    #   Exploration Analysis Section.
    # Beta(alpha=0.2, beta=82.7) prior built from Apollo (45 demands),
    # Soyuz (99 demands), military supply-drop (96,988 demands), and SRB
    # post-Challenger recovery data. Mean per-main per-mission failure
    # probability = alpha / (alpha + beta) = 0.2 / 82.9 ~= 2.41e-3.
    enable_stochastic_failure: bool = True
    main_failure_probability_per_chute: float = 0.00241

    # Used when enable_stochastic_failure is False. Lets you deterministically
    # force 2-of-3 (or fewer) scenarios for testing the pendulum behavior.
    num_mains_operational: int = 3                # 0..3 when stochastic OFF

    # -------- Wind drift under chutes --------
    enable_wind: bool = True
    wind_east_mps: float = 5.0
    wind_north_mps: float = 2.0

    # -------- RCS roll for splashdown --------
    enable_splash_orientation: bool = True
    splash_orientation_sigma_rad: float = 0.0     # target bank during chutes

    # -------- Behavioral flags carried over from v1 --------
    zero_lift_after_drogue: bool = True


# =============================================================================
# Runtime state
# =============================================================================

@dataclass
class CPASState:
    """Mutable state of the CPAS sequencer. One instance per simulation."""

    phase: CPASPhase = CPASPhase.STOWED

    # FBC
    fbc_jettisoned: bool = False
    fbc_jettison_t_s: float = math.nan

    # Phase timestamps
    drogue_deploy_t_s: float = math.nan
    pilot_deploy_t_s: float = math.nan
    main_deploy_t_s: float = math.nan
    landed_t_s: float = math.nan

    # Pendulum (1-D for simplicity; oscillates in a fixed azimuth direction)
    pendulum_angle_rad: float = 0.0
    pendulum_rate_rad_s: float = 0.0
    pendulum_initialized: bool = False

    # Per-main state (each independent)
    mains_operational_mask: Tuple[bool, bool, bool] = (True, True, True)
    mains_squidding_mask: Tuple[bool, bool, bool] = (False, False, False)
    squid_end_t_s: Tuple[float, float, float] = (math.nan, math.nan, math.nan)

    # Cumulative mass shed (for the dynamics)
    mass_shed_kg: float = 0.0

    # Last reefing stage we emitted -- used to detect disreef transitions
    last_reefing_stage: int = 0


# =============================================================================
# Per-step output
# =============================================================================

@dataclass
class CPASOutput:
    """What CPAS contributes to the translational dynamics this step."""

    phase: str
    # Drag CD*A from currently deployed chutes (with shock/reefing factored in)
    drag_area_cdA_m2: float
    # Multiplier on the aero lift coefficient
    lift_scale: float
    # Fraction of full inflation for the leading chute stage [0..1]
    open_fraction: float
    # True once mains are at full inflation
    force_vertical: bool

    # New high-fidelity outputs
    fbc_jettisoned: bool
    mass_shed_kg: float                 # cumulative shed mass this step
    mass_shed_delta_kg: float           # NEW shed this step (one-shot)
    reefing_stage: int                  # 0 stowed/extending, 1/2/3 stages
    num_mains_active: int               # operational chutes currently inflated
    num_mains_squidding: int            # currently squidding
    pendulum_angle_rad: float
    pendulum_rate_rad_s: float
    pendulum_lateral_v_mps: float       # horizontal velocity contribution
    pendulum_azimuth_rad: float         # direction of pendulum swing
    wind_east_mps: float
    wind_north_mps: float
    target_sigma_rad: float             # for RCS roll-for-splash (NaN if none)

    events: List[str] = field(default_factory=list)


# =============================================================================
# CPAS class
# =============================================================================

class CPAS:
    """Stateful Orion CPAS sequencer with reefing, pendulum, squidding,
    wind drift and failure-mode support."""

    def __init__(self, config: Optional[CPASConfig] = None,
                 rng: Optional[random.Random] = None,
                 pendulum_azimuth_rad: Optional[float] = None):
        """
        Build the CPAS sequencer and roll the random outcomes for this descent.

        Parameters:
            config:               tunable CPAS parameters. Defaults are used when
                                  None.
            rng:                  a random.Random for the failure and squidding
                                  draws. A fixed seed source is used when None.
            pendulum_azimuth_rad: fixed swing direction for the pendulum. A
                                  random direction is chosen when None.

        At construction each main is independently failed against the configured
        probability, or set from num_mains_operational when stochastic failure is
        off, and the pendulum constants are precomputed for the hot path.
        """
        self.config: CPASConfig = config or CPASConfig()
        self.state: CPASState = CPASState()
        self._rng = rng or random.Random(0)
        # The pendulum swings in a fixed azimuth direction (lat-lon frame).
        # If not given explicitly, pick a random one once.
        if pendulum_azimuth_rad is None:
            pendulum_azimuth_rad = self._rng.uniform(0.0, 2.0 * math.pi)
        self._pendulum_azimuth_rad: float = float(pendulum_azimuth_rad)

        # --- Failure-mode roll at construction time ---
        # When stochastic mode is enabled (default), each main is rolled
        # independently against the Fuqua 2010 Beta-prior mean failure
        # probability. Otherwise the deterministic num_mains_operational
        # knob is honored, so the user can force a 2-of-3 scenario for
        # testing the pendulum failure mode.
        if self.config.enable_stochastic_failure:
            p_fail = float(self.config.main_failure_probability_per_chute)
            mask = tuple(self._rng.random() >= p_fail for _ in range(3))
        else:
            n_ops = max(0, min(3, int(self.config.num_mains_operational)))
            mask = tuple(i < n_ops for i in range(3))
        self.state.mains_operational_mask = mask

        # Stash any "main_X_failed" events to emit on the first step()
        # so the events log records the determined failure pattern.
        self._pending_init_events: List[str] = [
            f"main{i + 1}_failed_at_init"
            for i in range(3) if not mask[i]
        ]

        # --- Pre-compute pendulum constants (don't depend on t/alt/V) ---
        # Hoisted out of the per-step _integrate_pendulum hot path.
        _omega = 2.0 * math.pi / max(1e-3, self.config.pendulum_period_s)
        _A_half = 0.5 * float(self.config.pendulum_target_amplitude_rad)
        self._pend_omega2 = _omega * _omega
        self._pend_A_half_sq = _A_half * _A_half
        self._pend_mu = float(self.config.pendulum_limit_cycle_strength)
        self._pend_L = float(self.config.pendulum_riser_length_m)

    # ------------------------------------------------------------------
    # Public step
    # ------------------------------------------------------------------
    def step(self,
             t_s: float,
             alt_m: float,
             V_mps: float,
             mach: Optional[float] = None,
             dt_s: float = 0.25) -> CPASOutput:
        """Advance the CPAS state machine and return all contributions."""
        cfg = self.config
        st = self.state
        events: List[str] = []
        mass_shed_delta_kg = 0.0

        # Surface any pending init-time events (e.g. main_X_failed_at_init)
        # on the very first step so the events log records them.
        if self._pending_init_events:
            events.extend(self._pending_init_events)
            self._pending_init_events = []

        # --- (a) Forward bay cover jettison (precedes the drogues) ---
        if cfg.enable_fbc and not st.fbc_jettisoned:
            if alt_m <= cfg.fbc_jettison_alt_m:
                st.fbc_jettisoned = True
                st.fbc_jettison_t_s = float(t_s)
                mass_shed_delta_kg = float(cfg.fbc_mass_shed_kg)
                st.mass_shed_kg += mass_shed_delta_kg
                events.append("fbc_jettison")

        # --- (b) Phase transitions (monotonic, one-shot) ---
        if st.phase == CPASPhase.STOWED and self._can_deploy_drogue(alt_m, V_mps, mach):
            st.phase = CPASPhase.DROGUE
            st.drogue_deploy_t_s = float(t_s)
            events.append("drogue_deploy")
        if st.phase == CPASPhase.DROGUE and alt_m <= cfg.pilot_deploy_alt_m:
            st.phase = CPASPhase.PILOT
            st.pilot_deploy_t_s = float(t_s)
            events.append("pilot_deploy")
            events.append("drogue_cut")  # drogues are pyro-cut the instant pilots fire
        if st.phase == CPASPhase.PILOT and alt_m <= cfg.main_deploy_alt_m:
            st.phase = CPASPhase.MAIN
            st.main_deploy_t_s = float(t_s)
            events.append("main_deploy")
            self._init_pendulum_if_needed()
        if st.phase == CPASPhase.MAIN and alt_m <= 0.0:
            st.phase = CPASPhase.LANDED
            st.landed_t_s = float(t_s)
            events.append("landed")

        # --- (c) Compute the multi-stage opening profile -------------
        drag_area, open_frac, reef_stage = self._drag_area_and_stage(t_s)

        # Emit a disreef event whenever the stage advances.
        if reef_stage > st.last_reefing_stage:
            if st.phase == CPASPhase.DROGUE and reef_stage == 2:
                events.append("drogue_disreef")
            elif st.phase == CPASPhase.MAIN and reef_stage == 2:
                events.append("main_disreef_1")
            elif st.phase == CPASPhase.MAIN and reef_stage == 3:
                events.append("main_disreef_2")
            st.last_reefing_stage = reef_stage
        elif reef_stage < st.last_reefing_stage:
            # Reset when entering a new phase (e.g., drogue->pilot->main)
            st.last_reefing_stage = reef_stage
        lift_scale = 1.0 if st.phase == CPASPhase.STOWED else (
            1.0 if not cfg.zero_lift_after_drogue else 0.0)

        # --- (d) Squidding (random per-main brief collapses) ---------
        n_active, n_squid, squid_events = self._update_squidding(t_s, dt_s)
        for ev in squid_events:
            events.append(ev)

        # Scale the main-chute contribution by (n_active / n_design)
        # to model the 2-of-3 failure mode AND active squidding.
        # _update_squidding already returned the per-frame counts; reuse them
        # instead of re-walking the masks.
        if st.phase == CPASPhase.MAIN:
            n_healthy = n_active - n_squid
            cluster_fraction = (n_healthy + n_squid * cfg.squid_area_fraction) / 3.0
            drag_area *= cluster_fraction

        # --- (e) Pendulum integration --------------------------------
        # Restricted to MAIN/LANDED phases — pendulum is a main-chute
        # phenomenon (Edwards 2024). Initial amplitude is set in
        # _init_pendulum_if_needed() at the main_deploy transition.
        if cfg.enable_pendulum and st.phase in (CPASPhase.MAIN, CPASPhase.LANDED):
            self._integrate_pendulum(dt_s)
        # Lateral velocity contribution: v_lat = L * theta_dot * cos(theta).
        # Use the cached riser length so the hot path stays allocation-free.
        if cfg.enable_pendulum and st.phase in (CPASPhase.MAIN, CPASPhase.LANDED):
            pendulum_v = self._pend_L * st.pendulum_rate_rad_s * math.cos(st.pendulum_angle_rad)
        else:
            pendulum_v = 0.0

        # --- (f) Wind drift (only applies when chutes are out) -------
        if cfg.enable_wind and st.phase != CPASPhase.STOWED:
            wind_e = cfg.wind_east_mps
            wind_n = cfg.wind_north_mps
        else:
            wind_e = 0.0
            wind_n = 0.0

        # --- (g) RCS roll-for-splash target --------------------------
        if cfg.enable_splash_orientation and st.phase in (
                CPASPhase.DROGUE, CPASPhase.PILOT, CPASPhase.MAIN, CPASPhase.LANDED):
            target_sigma = float(cfg.splash_orientation_sigma_rad)
        else:
            target_sigma = math.nan

        # --- (h) Force-vertical flag (kept from v1) -------------------
        force_vertical = (st.phase == CPASPhase.MAIN and reef_stage >= 3) \
                         or st.phase == CPASPhase.LANDED

        return CPASOutput(
            phase=st.phase.value,
            drag_area_cdA_m2=float(drag_area),
            lift_scale=float(lift_scale),
            open_fraction=float(open_frac),
            force_vertical=force_vertical,
            fbc_jettisoned=bool(st.fbc_jettisoned),
            mass_shed_kg=float(st.mass_shed_kg),
            mass_shed_delta_kg=float(mass_shed_delta_kg),
            reefing_stage=int(reef_stage),
            num_mains_active=int(n_active),
            num_mains_squidding=int(n_squid),
            pendulum_angle_rad=float(st.pendulum_angle_rad),
            pendulum_rate_rad_s=float(st.pendulum_rate_rad_s),
            pendulum_lateral_v_mps=float(pendulum_v),
            pendulum_azimuth_rad=float(self._pendulum_azimuth_rad),
            wind_east_mps=float(wind_e),
            wind_north_mps=float(wind_n),
            target_sigma_rad=float(target_sigma),
            events=events,
        )

    # ------------------------------------------------------------------
    # Convenience accessors
    # ------------------------------------------------------------------
    def summary(self) -> dict:
        """
        Return a snapshot of the sequencer state for logging.

        Returns:
            A dictionary with the current phase, the jettison and deploy
            timestamps, the total mass shed, the per main operational mask, and
            the pendulum swing direction in degrees.
        """
        return {
            "phase": self.state.phase.value,
            "fbc_jettisoned": bool(self.state.fbc_jettisoned),
            "fbc_jettison_t_s": self.state.fbc_jettison_t_s,
            "drogue_deploy_t_s": self.state.drogue_deploy_t_s,
            "pilot_deploy_t_s": self.state.pilot_deploy_t_s,
            "main_deploy_t_s": self.state.main_deploy_t_s,
            "landed_t_s": self.state.landed_t_s,
            "mass_shed_kg_total": self.state.mass_shed_kg,
            "mains_operational_mask": list(self.state.mains_operational_mask),
            "pendulum_azimuth_deg": math.degrees(self._pendulum_azimuth_rad),
        }

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _can_deploy_drogue(self, alt_m, V_mps, mach):
        """
        Report whether the drogue may deploy at the current condition.

        Parameters:
            alt_m: geometric altitude in meters.
            V_mps: speed in meters per second.
            mach:  Mach number, or None when not available.

        Returns:
            True only when the altitude is at or below the deploy altitude and
            the speed and Mach are within their deploy limits.
        """
        cfg = self.config
        if alt_m > cfg.drogue_deploy_alt_m:
            return False
        if V_mps > cfg.drogue_deploy_max_speed_mps:
            return False
        if mach is not None and mach > cfg.drogue_deploy_max_mach:
            return False
        return True

    def _drag_area_and_stage(self, t_s):
        """
        Compute the effective chute CD*A right now plus the integer
        reefing stage (0 = stowed/line-stretch, 1/2 = reefed stages,
        3 = fully open). Returns (drag_area, open_fraction, reefing_stage).
        """
        cfg = self.config
        st = self.state

        if st.phase == CPASPhase.DROGUE:
            t_since = t_s - st.drogue_deploy_t_s
            stage, open_fraction = self._opening_profile_drogue(t_since)
            drag_area = open_fraction * cfg.drogue_cd * cfg.drogue_ref_area_m2
            return drag_area, open_fraction, stage

        if st.phase == CPASPhase.PILOT:
            # Drogues are pyro-cut the instant pilots fire (Edwards 2024,
            # Fuqua 2010). Only the pilot chutes provide drag now -- they
            # in turn extract the mains from the forward bay.
            t_since_pilot = t_s - st.pilot_deploy_t_s
            pilot_ramp = max(0.0, min(1.0, t_since_pilot / max(1e-6, cfg.pilot_open_time_s)))
            pilot_part = pilot_ramp * cfg.pilot_cd * cfg.pilot_ref_area_m2
            return (pilot_part, pilot_ramp, 3)

        if st.phase == CPASPhase.MAIN:
            t_since = t_s - st.main_deploy_t_s
            stage, open_fraction = self._opening_profile_main(t_since)
            drag_area = open_fraction * cfg.main_cd * cfg.main_ref_area_m2
            return drag_area, open_fraction, stage

        if st.phase == CPASPhase.LANDED:
            return (cfg.main_cd * cfg.main_ref_area_m2, 1.0, 3)

        return (0.0, 0.0, 0)

    def _opening_profile_drogue(self, t_since: float):
        """Drogue 2-stage opening:
          [0, line_stretch)            zero drag (chute extending)
          [..., +shock_duration)       opening shock burst (>reefed)
          [..., reefed_duration)       reefed at drogue_reefed_area_fraction
          [reefed_duration, inf)       fully open
        """
        cfg = self.config
        ls = cfg.drogue_line_stretch_s
        sh = cfg.drogue_opening_shock_s
        rd = cfg.drogue_reefed_duration_s
        reefed = cfg.drogue_reefed_area_fraction
        if t_since < ls:
            return (0, 0.0)
        if t_since < ls + sh:
            return (1, reefed * cfg.drogue_opening_shock_factor)
        if t_since < ls + rd:
            return (1, reefed)
        return (2, 1.0)

    def _opening_profile_main(self, t_since: float):
        """Main 3-stage opening:
          [0, line_stretch)               zero drag
          [..., +shock_duration)          opening shock burst
          [..., line_stretch + reef1_dur) reef stage 1 (~13%)
          [..., +reef2_dur)               reef stage 2 (~35%)
          [..., inf)                      fully open
        """
        cfg = self.config
        ls = cfg.main_line_stretch_s
        sh = cfg.main_opening_shock_s
        d1 = cfg.main_reef1_duration_s
        d2 = cfg.main_reef2_duration_s
        f1 = cfg.main_reef1_area_fraction
        f2 = cfg.main_reef2_area_fraction
        if t_since < ls:
            return (0, 0.0)
        if t_since < ls + sh:
            return (1, f1 * cfg.main_opening_shock_factor)
        if t_since < ls + d1:
            return (1, f1)
        if t_since < ls + d1 + d2:
            return (2, f2)
        return (3, 1.0)

    def _update_squidding(self, t_s: float, dt_s: float):
        """Independently update squidding state for each operational main."""
        cfg = self.config
        st = self.state
        events: List[str] = []
        if st.phase != CPASPhase.MAIN or not cfg.enable_squidding:
            n_active = sum(1 for i in range(3)
                           if st.mains_operational_mask[i])
            return (n_active, 0, events)

        new_squid = list(st.mains_squidding_mask)
        new_end = list(st.squid_end_t_s)
        p_step = 1.0 - math.exp(-cfg.squid_probability_per_second * float(dt_s))

        for i in range(3):
            if not st.mains_operational_mask[i]:
                continue
            if new_squid[i]:
                # End squidding if duration expired
                if t_s >= new_end[i]:
                    new_squid[i] = False
                    new_end[i] = math.nan
                    events.append(f"main{i+1}_squid_end")
            else:
                # Roll for new squid event
                if self._rng.random() < p_step:
                    new_squid[i] = True
                    new_end[i] = t_s + cfg.squid_duration_s
                    events.append(f"main{i+1}_squid_start")

        st.mains_squidding_mask = tuple(new_squid)
        st.squid_end_t_s = tuple(new_end)

        n_active = sum(1 for i in range(3) if st.mains_operational_mask[i])
        n_squid = sum(1 for i in range(3) if st.mains_squidding_mask[i])
        return (n_active, n_squid, events)

    def _init_pendulum_if_needed(self):
        """
        Kick the pendulum at main_deploy ONLY when the cluster is
        degraded (one or more mains lost). With all three mains healthy
        the cluster is stable -- per Edwards et al. 2024, the planar
        pendulum is an emergent 2-of-3 failure-mode behavior and is not
        observed in nominal 3-chute descents.

        For the healthy-3 case we deliberately leave pendulum_initialized
        = False so the integrator early-returns every step, instead of
        running the Van der Pol math on a permanent (0, 0) state.
        """
        if self.state.pendulum_initialized:
            return
        n_active = sum(1 for ok in self.state.mains_operational_mask if ok)
        if n_active >= 3:
            return  # stable cluster, integrator stays a no-op
        # 2-of-3 (or worse): kick the pendulum so the limit-cycle solver
        # has something to lock onto immediately.
        self.state.pendulum_angle_rad = float(self.config.pendulum_target_amplitude_rad)
        self.state.pendulum_rate_rad_s = 0.0
        self.state.pendulum_initialized = True

    def _integrate_pendulum(self, dt_s: float):
        """
        Van der Pol-style self-sustaining limit-cycle oscillator:

            theta_ddot + mu*(theta^2 - A^2)*theta_dot + omega^2*sin(theta) = 0

        Below |theta| < A the damping term is negative (pumps energy in);
        above it is positive (bleeds energy out). Result: amplitude
        converges to A and holds there. Initial amplitude is set in
        _init_pendulum_if_needed() at the main_deploy transition (and
        only when one or more mains have failed).
        """
        st = self.state
        if not st.pendulum_initialized:
            return
        # Pendulum coefficients are constant across the descent — they are
        # precomputed in __init__ and cached on self.
        omega2 = self._pend_omega2
        A_half_sq = self._pend_A_half_sq
        mu = self._pend_mu
        theta = st.pendulum_angle_rad
        theta_dot = st.pendulum_rate_rad_s
        damping_coeff = mu * (theta * theta - A_half_sq)
        theta_ddot = -damping_coeff * theta_dot - omega2 * math.sin(theta)
        # Semi-implicit Euler
        theta_dot += theta_ddot * float(dt_s)
        theta += theta_dot * float(dt_s)
        st.pendulum_angle_rad = theta
        st.pendulum_rate_rad_s = theta_dot


# =============================================================================
# Module helpers
# =============================================================================

def apply_horizontal_drift_to_state(x_next, cpas_out, dt_s, R_earth_m):
    """
    In-place lat/lon update from CPAS wind + pendulum lateral velocity.

    Both drivers used to open-code this inline (same five lines twice). The
    capsule's east/north velocity contributions come from:
      - constant wind from the config (wind_east_mps, wind_north_mps)
      - pendulum swing at azimuth pendulum_azimuth_rad

    Operates only when chutes are deployed (cpas_out.phase != "stowed").
    Mutates x_next[1] (phi) and x_next[2] (lambda) and returns nothing.
    """
    if cpas_out.phase == "stowed":
        return
    phi = float(x_next[1])
    v_lat_pend = float(cpas_out.pendulum_lateral_v_mps)
    az = float(cpas_out.pendulum_azimuth_rad)
    v_east = float(cpas_out.wind_east_mps) + v_lat_pend * math.sin(az)
    v_north = float(cpas_out.wind_north_mps) + v_lat_pend * math.cos(az)
    x_next[1] += dt_s * v_north / R_earth_m
    x_next[2] += dt_s * v_east / (R_earth_m * max(1e-6, math.cos(phi)))


def terminal_velocity_mps(config: CPASConfig,
                          mass_kg: float,
                          rho_kgm3: float,
                          g_mps2: float = 9.81,
                          n_mains_active: Optional[int] = None) -> float:
    """Analytic terminal velocity under fully-open mains.

    If n_mains_active is given, scale by that fraction of the 3-chute cluster
    (so e.g. 2-of-3 mains gives a higher terminal velocity).
    """
    n = config.num_mains_operational if n_mains_active is None else n_mains_active
    cd_a = config.main_cd * config.main_ref_area_m2 * (float(n) / 3.0)
    if cd_a <= 0.0 or rho_kgm3 <= 0.0:
        return float("inf")
    return math.sqrt(2.0 * mass_kg * g_mps2 / (rho_kgm3 * cd_a))
