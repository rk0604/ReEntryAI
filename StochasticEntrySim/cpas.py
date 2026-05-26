"""
Orion Capsule Parachute Assembly System (CPAS) model.

Provides a self-contained state machine and force model for the three-stage
parachute sequence used by the Orion Crew Module for terminal descent:

    STOWED  ->  DROGUE  ->  PILOT  ->  MAIN  ->  LANDED

Public deployment specs are based on NASA Orion reference materials and
the Strahan et al. "Orion Entry Flight Control Stability and Performance"
paper. Approximate trigger conditions (defaults below are tunable):

    Drogues (2):  deploy at ~7.6 km altitude, Mach <= 0.9
                  total reference area ~76 m^2, CD ~0.5
                  reefed open over ~8 s
    Pilots  (3):  deploy at ~2.9 km altitude (drogues stay attached)
                  small chutes whose job is to pull out the mains
    Mains   (3):  deploy at ~2.4 km altitude
                  total reference area ~280 m^2, CD ~0.85
                  full open over ~8 s, terminal sink rate ~7.6 m/s

Integration model
-----------------
Each simulation step the caller invokes CPAS.step(...) with the current
altitude, speed and (optional) Mach number. CPAS returns a CPASOutput
describing what to add to the translational dynamics for this step:

    drag_area_cdA_m2 : add a chute drag force F = q_bar * drag_area_cdA_m2
                       opposite to the velocity vector
    lift_scale       : multiplier for the aero lift coefficient
                       (1.0 stowed, 0.0 once any chute is deployed)
    force_vertical   : True once mains are fully open. The capsule is
                       hanging from the chute lines and no longer flying
                       aerodynamically; recommend overriding gamma to
                       -pi/2 and freezing heading chi.

The module has no dependency on the rest of the simulator — it only does
arithmetic on the inputs you give it. That keeps it unit-testable in
isolation (see verify_cpas.py).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional


class CPASPhase(str, Enum):
    """Discrete deployment phases of the CPAS state machine."""
    STOWED = "stowed"
    DROGUE = "drogue"
    PILOT = "pilot"
    MAIN = "main"
    LANDED = "landed"


@dataclass
class CPASConfig:
    """Tunable CPAS parameters. Defaults reflect public Orion specs."""

    # --- Drogue chutes (2 of them) ---
    drogue_deploy_alt_m: float = 7600.0
    drogue_deploy_max_speed_mps: float = 250.0
    drogue_deploy_max_mach: float = 0.9
    drogue_cd: float = 0.50
    drogue_ref_area_m2: float = 76.0
    drogue_open_time_s: float = 8.0

    # --- Pilot chutes (3 of them, each ~3 m diameter) ---
    pilot_deploy_alt_m: float = 2900.0
    pilot_cd: float = 0.60
    pilot_ref_area_m2: float = 26.0    # 3 chutes x ~8.7 m^2 each
    pilot_open_time_s: float = 1.5

    # --- Main chutes (3 of them, each ~35 m / 116 ft diameter) ---
    # Disk area per chute: pi * (35.4/2)^2 ~= 984 m^2. 3 chutes -> ~2950 m^2.
    # With CD ~0.85 this yields V_terminal ~7.6 m/s for a 8500 kg capsule
    # at sea-level density.
    main_deploy_alt_m: float = 2400.0
    main_cd: float = 0.85
    main_ref_area_m2: float = 2950.0
    main_open_time_s: float = 8.0

    # --- Behavioral flags ---
    # Kill aerodynamic lift the moment any chute is deployed.
    zero_lift_after_drogue: bool = True
    # Force the velocity vector vertical once mains are fully open
    # (capsule hangs from chute lines, no longer flying).
    force_vertical_under_mains: bool = True


@dataclass
class CPASState:
    """Mutable state of the CPAS sequencer. One instance per simulation run."""

    phase: CPASPhase = CPASPhase.STOWED
    drogue_deploy_t_s: float = math.nan
    pilot_deploy_t_s: float = math.nan
    main_deploy_t_s: float = math.nan
    landed_t_s: float = math.nan


@dataclass
class CPASOutput:
    """What CPAS contributes to the translational dynamics this step."""

    phase: str
    # Effective drag CD*A from currently deployed chutes (m^2).
    # Drag force = 0.5 * rho * V^2 * drag_area_cdA_m2, opposite to velocity.
    drag_area_cdA_m2: float
    # Multiplier on the aero lift coefficient (1.0 = normal, 0.0 = killed).
    lift_scale: float
    # If True, the dynamics should force gamma to -pi/2 (capsule on mains).
    force_vertical: bool
    # Fraction of full inflation for the leading chute stage [0..1].
    open_fraction: float
    # Events this step (e.g. "drogue_deploy", "main_deploy").
    events: List[str] = field(default_factory=list)


class CPAS:
    """
    Orion CPAS state machine.

    Usage:
        cpas = CPAS()
        for k in range(num_steps):
            out = cpas.step(t_s=t, alt_m=alt, V_mps=V, mach=M)
            # Add chute drag to dynamics:
            drag_chute_N = 0.5 * rho * V**2 * out.drag_area_cdA_m2
            # Optionally zero out aero lift, force vertical, etc.
    """

    def __init__(self, config: Optional[CPASConfig] = None):
        self.config: CPASConfig = config or CPASConfig()
        self.state: CPASState = CPASState()

    # ------------------------------------------------------------------
    # Trigger predicates -- separated out for readability and testing.
    # ------------------------------------------------------------------
    def _should_deploy_drogue(self, alt_m: float, V_mps: float,
                              mach: Optional[float]) -> bool:
        cfg = self.config
        if alt_m > cfg.drogue_deploy_alt_m:
            return False
        if V_mps > cfg.drogue_deploy_max_speed_mps:
            return False
        if mach is not None and mach > cfg.drogue_deploy_max_mach:
            return False
        return True

    def _should_deploy_pilot(self, alt_m: float) -> bool:
        return alt_m <= self.config.pilot_deploy_alt_m

    def _should_deploy_main(self, alt_m: float) -> bool:
        return alt_m <= self.config.main_deploy_alt_m

    # ------------------------------------------------------------------
    # Main step
    # ------------------------------------------------------------------
    def step(self,
             t_s: float,
             alt_m: float,
             V_mps: float,
             mach: Optional[float] = None) -> CPASOutput:
        """Advance the CPAS state machine and return what to apply this step."""
        cfg = self.config
        st = self.state
        events: List[str] = []

        # --- One-shot phase transitions (monotonic) ---
        if st.phase == CPASPhase.STOWED:
            if self._should_deploy_drogue(alt_m, V_mps, mach):
                st.phase = CPASPhase.DROGUE
                st.drogue_deploy_t_s = float(t_s)
                events.append("drogue_deploy")

        if st.phase == CPASPhase.DROGUE:
            if self._should_deploy_pilot(alt_m):
                st.phase = CPASPhase.PILOT
                st.pilot_deploy_t_s = float(t_s)
                events.append("pilot_deploy")

        if st.phase == CPASPhase.PILOT:
            if self._should_deploy_main(alt_m):
                st.phase = CPASPhase.MAIN
                st.main_deploy_t_s = float(t_s)
                events.append("main_deploy")

        if st.phase == CPASPhase.MAIN:
            if alt_m <= 0.0:
                st.phase = CPASPhase.LANDED
                st.landed_t_s = float(t_s)
                events.append("landed")

        # --- Compute effective drag area, lift scale, vertical flag ---
        drag_area = 0.0
        lift_scale = 1.0
        force_vertical = False
        open_fraction = 0.0

        if st.phase == CPASPhase.DROGUE:
            open_fraction = _ramp(t_s, st.drogue_deploy_t_s, cfg.drogue_open_time_s)
            drag_area = open_fraction * cfg.drogue_cd * cfg.drogue_ref_area_m2
            lift_scale = 0.0 if cfg.zero_lift_after_drogue else 1.0

        elif st.phase == CPASPhase.PILOT:
            # Drogues remain attached and fully open; pilots add a little drag.
            drogue_full = cfg.drogue_cd * cfg.drogue_ref_area_m2
            open_fraction = _ramp(t_s, st.pilot_deploy_t_s, cfg.pilot_open_time_s)
            pilot_partial = open_fraction * cfg.pilot_cd * cfg.pilot_ref_area_m2
            drag_area = drogue_full + pilot_partial
            lift_scale = 0.0

        elif st.phase == CPASPhase.MAIN:
            # Mains replace drogues / pilots once they take over.
            open_fraction = _ramp(t_s, st.main_deploy_t_s, cfg.main_open_time_s)
            drag_area = open_fraction * cfg.main_cd * cfg.main_ref_area_m2
            lift_scale = 0.0
            force_vertical = (cfg.force_vertical_under_mains and open_fraction >= 1.0)

        elif st.phase == CPASPhase.LANDED:
            drag_area = cfg.main_cd * cfg.main_ref_area_m2
            lift_scale = 0.0
            force_vertical = cfg.force_vertical_under_mains
            open_fraction = 1.0

        return CPASOutput(
            phase=st.phase.value,
            drag_area_cdA_m2=drag_area,
            lift_scale=lift_scale,
            force_vertical=force_vertical,
            open_fraction=open_fraction,
            events=events,
        )

    # ------------------------------------------------------------------
    # Convenience accessors
    # ------------------------------------------------------------------
    def summary(self) -> dict:
        """Return a JSON-friendly snapshot of all deployment timestamps."""
        return {
            "phase": self.state.phase.value,
            "drogue_deploy_t_s": self.state.drogue_deploy_t_s,
            "pilot_deploy_t_s": self.state.pilot_deploy_t_s,
            "main_deploy_t_s": self.state.main_deploy_t_s,
            "landed_t_s": self.state.landed_t_s,
        }


def _ramp(t_s: float, deploy_t_s: float, open_time_s: float) -> float:
    """Linear ramp from 0 to 1 over open_time_s starting at deploy_t_s."""
    if math.isnan(deploy_t_s):
        return 0.0
    if open_time_s <= 0.0:
        return 1.0
    return max(0.0, min(1.0, (t_s - deploy_t_s) / open_time_s))


def terminal_velocity_mps(config: CPASConfig,
                          mass_kg: float,
                          rho_kgm3: float,
                          g_mps2: float = 9.81) -> float:
    """
    Analytic terminal velocity under fully-open main chutes.

    V_terminal = sqrt( 2 * m * g / (rho * CD * A) )
    """
    cd_a = config.main_cd * config.main_ref_area_m2
    if cd_a <= 0.0 or rho_kgm3 <= 0.0:
        return float("inf")
    return math.sqrt(2.0 * mass_kg * g_mps2 / (rho_kgm3 * cd_a))
