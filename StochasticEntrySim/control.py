"""
control.py

Milestone-1 capsule control stack for atmospheric entry.

What this file does
-------------------
This module owns the high-level control path:

    raw vehicle state
        -> observation features
        -> 1 Hz guidance update
        -> smoothed bank target
        -> roll torque command
        -> physical RCS / attitude layer elsewhere
        -> actual bank angle comes back from physics

This is the key architectural change from the old version:
the control stack no longer directly "flies" sigma itself.

Instead:
- guidance computes an ideal bank command sigma_cmd
- the bank actuator smooths that into sigma_target
- a roll torque controller compares sigma_target against the
  actual bank angle being flown by the capsule
- the physical RCS + attitude dynamics layer is responsible
  for turning that torque command into real motion

So in milestone 1, this file outputs a torque command.
It does not directly integrate the true capsule attitude.

Main flow
---------
1. Observe the current state and build guidance features.
2. Let guidance refresh sigma_cmd only at the guidance cadence.
3. Smooth sigma_cmd into sigma_target.
4. Compare sigma_target against sigma_actual from physics.
5. Output tau_roll_cmd for the RCS / rigid-body layer.

Notes
-----
- This file keeps your existing guidance structure and observation layer.
- The old "achieved sigma" logic has been split:
    * sigma_target is owned here
    * sigma_actual is owned by the rigid-body physics layer
- Aerodynamic acceleration is still computed because it is useful telemetry.
  It is no longer the main gate for whether control is "real".
"""

import math
import AtmosphereModel
import constants

from dataclasses import dataclass, field
from typing import Dict, Optional, Protocol


# ------------------------------------------------------------------------------
# Small math helpers
# ------------------------------------------------------------------------------

def wrap_to_pi(angle_rad: float) -> float:
    """
    Wrap an angle into [-pi, pi].

    This is important because sigma is an angle. Once you start smoothing
    commands or comparing actual vs target bank, you want the *shortest*
    angular error, not a raw difference that might jump across ±pi.
    """
    while angle_rad > math.pi:
        angle_rad -= 2.0 * math.pi
    while angle_rad < -math.pi:
        angle_rad += 2.0 * math.pi
    return angle_rad


def clamp(value: float, lower: float, upper: float) -> float:
    """
    Clamp a scalar into [lower, upper].
    """
    return max(lower, min(upper, value))


def sign_no_zero(x: float) -> float:
    """
    Return +1 for positive, -1 for negative, 0 for zero.
    """
    if x > 0.0:
        return 1.0
    if x < 0.0:
        return -1.0
    return 0.0


# ------------------------------------------------------------------------------
# Core state and config
# ------------------------------------------------------------------------------

@dataclass
class ReentryState:
    """
    Minimal control-facing reentry state.

    Primary translational states
    ----------------------------
    r_m
        Radial distance from Earth center, meters.

    phi_rad
        Geocentric latitude, radians.

    lam_rad
        Longitude, radians.

    V_mps
        Speed, meters per second.

    gamma_rad
        Flight-path angle, radians.

    chi_rad
        Heading angle, radians.

    New milestone-1 attitude/control fields
    ---------------------------------------
    sigma_actual_rad
        Actual bank angle currently being flown by the capsule.
        This comes back from the rigid-body / RCS layer, not from
        the high-level controller.

    roll_rate_rad_s
        Actual roll rate of the capsule.

    sigma_cmd_rad
        Latest ideal bank command issued by guidance.

    sigma_target_rad
        Smoothed target bank angle after actuator shaping.

    Why both sigma_cmd and sigma_target exist
    -----------------------------------------
    sigma_cmd is what guidance *wants*.
    sigma_target is what the command shaper *asks the rigid body to track*.
    sigma_actual is what the rigid body is *actually doing*.
    """
    r_m: float
    phi_rad: float
    lam_rad: float
    V_mps: float
    gamma_rad: float
    chi_rad: float

    sigma_actual_rad: float = 0.0
    roll_rate_rad_s: float = 0.0

    sigma_cmd_rad: float = 0.0
    sigma_target_rad: float = 0.0

    # ------------------------------------------------------------------
    # Backward-compatibility aliases
    # ------------------------------------------------------------------
    # These properties help older code that still expects sigma_rad and
    # sigma_dot_rps from the previous control model.

    @property
    def sigma_rad(self) -> float:
        """
        Compatibility alias for the actual flown bank angle.
        """
        return self.sigma_actual_rad

    @sigma_rad.setter
    def sigma_rad(self, value: float) -> None:
        self.sigma_actual_rad = float(value)

    @property
    def sigma_dot_rps(self) -> float:
        """
        Compatibility alias for the actual roll rate.
        """
        return self.roll_rate_rad_s

    @sigma_dot_rps.setter
    def sigma_dot_rps(self, value: float) -> None:
        self.roll_rate_rad_s = float(value)


@dataclass
class SimConfig:
    """
    Basic simulation timing configuration.

    dt_s
        Physics integration step.

    max_steps
        Max number of simulation steps.

    guidance_period_s
        How often guidance is allowed to refresh sigma_cmd.
    """
    dt_s: float = constants.ENTRY_DT_S
    max_steps: int = 2000
    guidance_period_s: float = constants.GUIDANCE_PERIOD_S


# ------------------------------------------------------------------------------
# Protocols
# ------------------------------------------------------------------------------

class GuidanceLaw(Protocol):
    """
    Structural interface for any bank-guidance module.
    """
    def reset(self) -> None:
        ...

    def compute_sigma_cmd(
        self,
        t_s: float,
        state: ReentryState,
        obs: Dict[str, float],
    ) -> float:
        ...


class ObservationProvider(Protocol):
    """
    Structural interface for any module that converts raw vehicle state
    into guidance features.
    """
    def reset(self) -> None:
        ...

    def observe(self, state: ReentryState, t_s: float) -> Dict[str, float]:
        ...


# ------------------------------------------------------------------------------
# Guidance law
# ------------------------------------------------------------------------------

@dataclass
class SimpleBankGuidance:
    """
    Simple capsule-style bank guidance for atmospheric entry.

    Behavior
    --------
    1. Steering is disabled above a velocity threshold.
    2. Crossrange error selects bank sign.
    3. Downrange error adjusts bank magnitude.
    4. Output is the ideal bank command sigma_cmd.

    This guidance law is still intentionally simple.
    It is a good first stage because it gives you a clean sigma_cmd
    without mixing guidance logic with rigid-body attitude control.
    """
    max_bank_deg: float = 70.0
    min_bank_deg: float = 40.0
    downrange_gain: float = 1e-6
    velocity_enable_mps: float = 6000.0

    def reset(self) -> None:
        """
        No persistent guidance memory yet.
        """
        pass

    def compute_sigma_cmd(
        self,
        t_s: float,
        state: ReentryState,
        obs: Dict[str, float],
    ) -> float:
        """
        Return the ideal commanded bank angle in radians.

        The velocity gate is preserved from your earlier version:
        above a threshold speed, guidance simply commands zero bank.
        """
        if state.V_mps > self.velocity_enable_mps:
            return 0.0

        cross_err = obs.get("crossrange_error_m", 0.0)
        downrange_err = obs.get("downrange_error_m", 0.0)

        # Crossrange decides which side to bank toward.
        if cross_err > 0.0:
            bank_sign = 1.0
        elif cross_err < 0.0:
            bank_sign = -1.0
        else:
            bank_sign = 0.0

        base_mag = math.radians(self.min_bank_deg)
        max_mag = math.radians(self.max_bank_deg)

        mag_adjust = abs(downrange_err) * self.downrange_gain
        mag = clamp(base_mag + mag_adjust, base_mag, max_mag)

        sigma_cmd = bank_sign * mag
        return float(wrap_to_pi(sigma_cmd))


# ------------------------------------------------------------------------------
# Guidance scheduler
# ------------------------------------------------------------------------------

class GuidanceScheduler:
    """
    Gate guidance updates so sigma_cmd refreshes more slowly than the
    physics loop.

    Example:
    dt = 0.25 s and guidance period = 1.0 s
    means guidance only issues a fresh sigma_cmd once every 4 physics steps.
    """

    def __init__(self, period_s: float):
        self.period_s = float(period_s)
        self._last_update_time: Optional[float] = None

    def reset(self) -> None:
        self._last_update_time = None

    def should_update(self, t_s: float) -> bool:
        """
        Return True only when guidance is allowed to update.

        The first call always updates so the command hold is initialized
        immediately at the beginning of a trajectory.
        """
        if self._last_update_time is None:
            self._last_update_time = t_s
            return True

        elapsed = t_s - self._last_update_time
        if elapsed >= self.period_s:
            self._last_update_time = t_s
            return True

        return False


# ------------------------------------------------------------------------------
# Aero authority utility
# ------------------------------------------------------------------------------

def compute_aero_accel_mps2(
    lift_N: float,
    drag_N: float,
    mass_kg: float,
) -> float:
    """
    Return aerodynamic acceleration magnitude.

    Formula
    -------
    sqrt(L^2 + D^2) / m

    Why this is kept
    ----------------
    In the old file this was part of the steering gate.
    In milestone 1 it is mainly useful as telemetry and debugging:
    it tells you how much the atmosphere is contributing relative
    to the RCS-driven attitude motion.
    """
    if mass_kg <= 0.0:
        raise ValueError("mass_kg must be positive")

    total_force = math.sqrt(lift_N * lift_N + drag_N * drag_N)
    return total_force / mass_kg


# ------------------------------------------------------------------------------
# Observation layer
# ------------------------------------------------------------------------------

class BasicObservationProvider:
    """
    Build guidance features from the raw reentry state.

    Current approximation
    ---------------------
    - latitude error -> crossrange
    - longitude error -> downrange

    This is not yet a full great-circle targeting model, but it is
    a good first observation layer for closing the guidance loop.
    """

    def __init__(
        self,
        target_phi_rad: float,
        target_lam_rad: float,
        params: Optional[dict] = None,
    ):
        self.target_phi = float(target_phi_rad)
        self.target_lam = float(target_lam_rad)
        self.params = params or {}

    def reset(self) -> None:
        pass

    def observe(self, state: ReentryState, t_s: float) -> Dict[str, float]:
        """
        Return a feature dictionary for guidance.

        Returned features include:
        - altitude
        - dynamic pressure
        - crossrange error
        - downrange error
        - speed
        - flight-path angle
        - actual sigma
        - actual roll rate
        """
        altitude_m = state.r_m - constants.RADIUS_EARTH
        altitude_m = max(0.0, altitude_m)

        atm = AtmosphereModel.US_Standard_ATM(altitude_m)
        rho = atm["rho_kgm3"] if atm["rho_kgm3"] is not None else 0.0

        dynamic_pressure_pa = 0.5 * rho * state.V_mps * state.V_mps

        R = constants.RADIUS_EARTH
        dphi = state.phi_rad - self.target_phi
        dlam = state.lam_rad - self.target_lam

        downrange_error_m = R * dlam * math.cos(self.target_phi)
        crossrange_error_m = R * dphi

        return {
            "altitude_m": float(altitude_m),
            "dynamic_pressure_pa": float(dynamic_pressure_pa),
            "crossrange_error_m": float(crossrange_error_m),
            "downrange_error_m": float(downrange_error_m),
            "speed_mps": float(state.V_mps),
            "flight_path_angle_rad": float(state.gamma_rad),
            "sigma_actual_rad": float(state.sigma_actual_rad),
            "roll_rate_rad_s": float(state.roll_rate_rad_s),
        }


# ------------------------------------------------------------------------------
# Bank command shaper
# ------------------------------------------------------------------------------

@dataclass
class BankActuatorLimits:
    """
    Limits for smoothing sigma_cmd into sigma_target.

    sigma_rate_max_rps
        Max target bank rate.

    sigma_accel_max_rps2
        Max target bank angular acceleration.
    """
    sigma_rate_max_rps: float
    sigma_accel_max_rps2: float


class BankActuator:
    """
    Smooth sigma_cmd into sigma_target using simple second-order limits.

    Important design change
    -----------------------
    In the old version, this class tried to produce the achieved bank angle
    directly, which made the control file act like the capsule attitude
    dynamics were already solved.

    In the revised version, this class only produces a smooth *target* bank
    angle that the rigid-body / RCS layer should try to track.
    """

    def __init__(self, limits: BankActuatorLimits):
        self.limits = limits
        self.sigma_target_rad = 0.0
        self.sigma_target_dot_rps = 0.0

    def reset(
        self,
        sigma_target0_rad: float = 0.0,
        sigma_target_dot0_rps: float = 0.0,
    ) -> None:
        self.sigma_target_rad = float(wrap_to_pi(sigma_target0_rad))
        self.sigma_target_dot_rps = float(sigma_target_dot0_rps)

    def step(
        self,
        dt_s: float,
        sigma_cmd_rad: float,
    ) -> tuple[float, float]:
        """
        Advance the command shaper by one step.

        Output
        ------
        sigma_target_rad
            Smoothed target bank angle.

        sigma_target_dot_rps
            Rate of the smoothed target.

        Why this exists
        ---------------
        Guidance should be allowed to make large discrete bank decisions
        at 1 Hz, but the rigid-body system should be asked to track a
        smoother target rather than an instantaneous step.
        """
        if dt_s <= 0.0:
            raise ValueError("dt_s must be positive")

        sigma_cmd_rad = wrap_to_pi(sigma_cmd_rad)
        prev_error = wrap_to_pi(sigma_cmd_rad - self.sigma_target_rad)

        # If target already matches command, stop moving.
        if abs(prev_error) <= 1.0e-12:
            self.sigma_target_rad = sigma_cmd_rad
            self.sigma_target_dot_rps = 0.0
            return self.sigma_target_rad, self.sigma_target_dot_rps

        # Simple bang-bang target acceleration toward the command.
        sigma_ddot_cmd = sign_no_zero(prev_error) * self.limits.sigma_accel_max_rps2

        self.sigma_target_dot_rps += sigma_ddot_cmd * dt_s
        self.sigma_target_dot_rps = clamp(
            self.sigma_target_dot_rps,
            -self.limits.sigma_rate_max_rps,
            self.limits.sigma_rate_max_rps,
        )

        self.sigma_target_rad = wrap_to_pi(
            self.sigma_target_rad + self.sigma_target_dot_rps * dt_s
        )

        # Overshoot protection:
        # if we crossed the command in this step, clamp exactly to it.
        new_error = wrap_to_pi(sigma_cmd_rad - self.sigma_target_rad)
        crossed = (
            (prev_error > 0.0 and new_error < 0.0)
            or (prev_error < 0.0 and new_error > 0.0)
            or abs(new_error) <= 1.0e-12
        )

        if crossed:
            self.sigma_target_rad = sigma_cmd_rad
            self.sigma_target_dot_rps = 0.0

        return self.sigma_target_rad, self.sigma_target_dot_rps


# ------------------------------------------------------------------------------
# Roll torque controller
# ------------------------------------------------------------------------------

@dataclass
class RollTorqueController:
    """
    PD controller that converts bank-tracking error into a roll torque command.

    Formula
    -------
    tau_roll_cmd = kp * (sigma_target - sigma_actual) - kd * roll_rate

    Interpretation
    --------------
    - If the capsule is not yet at the desired bank angle, command roll torque.
    - If the capsule is already rolling too fast, the damping term reduces torque.
    """

    kp_Nm_per_rad: float = constants.ROLL_KP_NM_PER_RAD
    kd_Nm_per_rad_s: float = constants.ROLL_KD_NM_PER_RAD_S
    max_torque_Nm: float = constants.ROLL_CMD_MAX_TORQUE_NM

    sigma_deadband_rad: float = constants.ROLL_SIGMA_DEADBAND_RAD
    rate_deadband_rad_s: float = constants.ROLL_RATE_DEADBAND_RAD_S

    def reset(self) -> None:
        """
        No persistent memory yet.
        """
        pass

    def step(
        self,
        sigma_target_rad: float,
        sigma_actual_rad: float,
        roll_rate_rad_s: float,
    ) -> float:
        """
        Return the commanded body roll torque, in N m.

        Important:
        This controller compares the *target* bank angle against the *actual*
        bank angle coming back from the rigid-body / RCS layer.
        """
        sigma_error = wrap_to_pi(sigma_target_rad - sigma_actual_rad)

        # Deadband prevents torque chatter when the capsule is already
        # very close to target and not rotating meaningfully.
        if (
            abs(sigma_error) < self.sigma_deadband_rad
            and abs(roll_rate_rad_s) < self.rate_deadband_rad_s
        ):
            return 0.0

        tau_roll_cmd = (
            self.kp_Nm_per_rad * sigma_error
            - self.kd_Nm_per_rad_s * roll_rate_rad_s
        )

        tau_roll_cmd = clamp(
            tau_roll_cmd,
            -self.max_torque_Nm,
            self.max_torque_Nm,
        )
        return float(tau_roll_cmd)


# ------------------------------------------------------------------------------
# Control stack config and debug output
# ------------------------------------------------------------------------------

@dataclass
class CapsuleControlConfig:
    """
    Configuration shared across the full capsule control stack.
    """
    guidance_period_s: float = constants.GUIDANCE_PERIOD_S


@dataclass
class CapsuleControlOutput:
    """
    Per-step control telemetry.

    This is the most useful object to log when validating the new flow.

    obs
        What guidance saw.

    sigma_cmd_rad
        Ideal bank command from guidance.

    sigma_target_rad
        Smoothed bank target sent toward the attitude loop.

    sigma_actual_rad
        Actual bank angle currently being flown.

    roll_rate_rad_s
        Actual roll rate currently being flown.

    tau_roll_cmd_Nm
        Torque request sent to the physical RCS / attitude layer.
    """
    obs: Dict[str, float] = field(default_factory=dict)

    guidance_updated: bool = False

    sigma_cmd_rad: float = 0.0
    sigma_target_rad: float = 0.0
    sigma_target_dot_rps: float = 0.0

    sigma_actual_rad: float = 0.0
    roll_rate_rad_s: float = 0.0

    tau_roll_cmd_Nm: float = 0.0
    aero_accel_mps2: float = 0.0


# ------------------------------------------------------------------------------
# Full stack
# ------------------------------------------------------------------------------

class CapsuleControlStack:
    """
    Full high-level control chain for milestone 1.

    Flow inside step
    ----------------
    1. Build observations.
    2. Decide if guidance may refresh sigma_cmd.
    3. Hold sigma_cmd between guidance updates.
    4. Smooth sigma_cmd into sigma_target.
    5. Read sigma_actual and roll_rate from the current state.
    6. Compute roll torque command to track sigma_target.
    7. Return telemetry for the rigid-body / RCS layer.

    Important design choice
    -----------------------
    This class does NOT update sigma_actual itself.

    Why?
    Because actual bank angle is now a *physical result* of:
    - RCS jet firing
    - rigid-body rotational dynamics
    - capsule inertia
    - quaternion / roll-state integration

    That physics belongs in the attitude / RCS layer, not here.
    """

    def __init__(
        self,
        cfg: CapsuleControlConfig,
        scheduler: GuidanceScheduler,
        guidance: GuidanceLaw,
        obs_provider: ObservationProvider,
        bank_actuator: BankActuator,
        roll_controller: RollTorqueController,
    ):
        self.cfg = cfg
        self.scheduler = scheduler
        self.guidance = guidance
        self.obs_provider = obs_provider
        self.bank_actuator = bank_actuator
        self.roll_controller = roll_controller

        # Held guidance command between 1 Hz refreshes.
        self._sigma_cmd_hold_rad = 0.0

        # Latest telemetry for logging / notebook inspection.
        self.debug = CapsuleControlOutput()

    def reset(self) -> None:
        """
        Reset all stateful pieces so a new trajectory starts cleanly.
        """
        self.scheduler.reset()
        self.guidance.reset()
        self.obs_provider.reset()
        self.bank_actuator.reset()
        self.roll_controller.reset()

        self._sigma_cmd_hold_rad = 0.0
        self.debug = CapsuleControlOutput()

    def step(
        self,
        t_s: float,
        dt_s: float,
        state: ReentryState,
        lift_N: float,
        drag_N: float,
        mass_kg: float,
    ) -> CapsuleControlOutput:
        """
        Advance the high-level control stack by one step.

        Inputs
        ------
        t_s, dt_s
            Current time and step size.

        state
            Current control-facing state. In particular, this must already
            contain sigma_actual_rad and roll_rate_rad_s from the physics layer.

        lift_N, drag_N, mass_kg
            Used only to compute aero acceleration telemetry.

        Output
        ------
        CapsuleControlOutput

        Most important line
        -------------------
        The roll controller is the key new bridge:

            tau_roll_cmd = roll_controller.step(
                sigma_target_rad,
                sigma_actual_rad,
                roll_rate_rad_s
            )

        That is where the old sigma-only logic becomes a physically meaningful
        torque request for the RCS / attitude layer.
        """
        if dt_s <= 0.0:
            raise ValueError("dt_s must be positive")

        # 1. Build observations from the raw vehicle state.
        obs = self.obs_provider.observe(state, t_s)

        # 2. Refresh sigma_cmd only at the guidance cadence.
        guidance_updated = self.scheduler.should_update(t_s)
        if guidance_updated:
            self._sigma_cmd_hold_rad = float(
                self.guidance.compute_sigma_cmd(t_s, state, obs)
            )

        # 3. Smooth sigma_cmd into sigma_target.
        sigma_target_rad, sigma_target_dot_rps = self.bank_actuator.step(
            dt_s=dt_s,
            sigma_cmd_rad=self._sigma_cmd_hold_rad,
        )

        # 4. Read actual bank motion from the physics-owned state.
        sigma_actual_rad = float(state.sigma_actual_rad)
        roll_rate_rad_s = float(state.roll_rate_rad_s)

        # 5. Compute the roll torque needed to track sigma_target.
        tau_roll_cmd_Nm = self.roll_controller.step(
            sigma_target_rad=sigma_target_rad,
            sigma_actual_rad=sigma_actual_rad,
            roll_rate_rad_s=roll_rate_rad_s,
        )

        # 6. Compute aero acceleration only as useful telemetry.
        aero_accel_mps2 = compute_aero_accel_mps2(
            lift_N=lift_N,
            drag_N=drag_N,
            mass_kg=mass_kg,
        )

        # 7. Mirror command/target values back into the state so logs,
        # notebooks, and surrounding code can inspect them easily.
        state.sigma_cmd_rad = float(self._sigma_cmd_hold_rad)
        state.sigma_target_rad = float(sigma_target_rad)

        self.debug = CapsuleControlOutput(
            obs=dict(obs),
            guidance_updated=bool(guidance_updated),
            sigma_cmd_rad=float(self._sigma_cmd_hold_rad),
            sigma_target_rad=float(sigma_target_rad),
            sigma_target_dot_rps=float(sigma_target_dot_rps),
            sigma_actual_rad=float(sigma_actual_rad),
            roll_rate_rad_s=float(roll_rate_rad_s),
            tau_roll_cmd_Nm=float(tau_roll_cmd_Nm),
            aero_accel_mps2=float(aero_accel_mps2),
        )

        return self.debug


# ------------------------------------------------------------------------------
# Backward-compatibility aliases
# ------------------------------------------------------------------------------
# These aliases let older code keep importing the old names while you migrate
# the rest of the simulator to the milestone-1 architecture.

SigmaControlConfig = CapsuleControlConfig
SigmaControlDebug = CapsuleControlOutput
SigmaControlStack = CapsuleControlStack