"""
ReactionControl.py

Milestone-1 physical RCS model for the entry capsule

What this file does
-------------------
This module replaces the old timing-only RCS placeholder with a physical,
body-frame thruster model that can:

- define fixed thrusters on the capsule body
- compute force from each thruster in body coordinates
- compute torque from each thruster using r x F
- allocate roll thrusters from a commanded roll torque
- respect the fixed Euler step dt = 0.25 s
- return a real body-frame wrench for the attitude integrator

Milestone-1 scope
-----------------
This file is intentionally scoped to the first useful step:
- roll control is active
- pitch and yaw geometry is already represented in the layout
- pitch and yaw allocation can be added later without replacing the file

Important design choice
-----------------------
This module works in the BODY frame only

That is deliberate:
- thrusters are fixed in body coordinates
- force and torque are easiest to compute in body coordinates
- world-frame pose conversion should stay in math_3d.py, where attitude
  and frame transforms already belong

Numerical behavior
------------------
The physics loop runs at dt = 0.25 s.
To stay compatible with that fixed-step loop, minimum pulse width is modeled
as one full simulation step.

So each step, a selected thruster is either:
- OFF for the whole step
- ON for the whole step

A small commanded torque is handled through a simple per-channel duty
accumulator, which effectively creates pulse-width modulation across steps.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import constants


# Small math helpers ------------------------------------------------------------------

def clamp(value: float, lower: float, upper: float) -> float:
    """
    Clamp a scalar into [lower, upper].
    """
    return max(lower, min(upper, value))


def normalize(vec: np.ndarray, eps: float = 1.0e-12) -> np.ndarray:
    """
    Return vec normalized to unit length.

    Raises
    ------
    ValueError if the vector norm is too small.
    """
    arr = np.asarray(vec, dtype=float).reshape(3)
    n = float(np.linalg.norm(arr))
    if n < eps:
        raise ValueError("Cannot normalize near-zero vector")
    return arr / n

# Core thruster and wrench data models ------------------------------------------------------------------

@dataclass
class RCSThruster:
    """
    One body-fixed RCS thruster.

    Parameters
    ----------
    name
        Human-readable thruster identifier.

    position_b_m
        Thruster location in the BODY frame, meters, relative to the chosen
        body origin or center of mass reference.

    direction_b_unit
        Unit thrust direction in the BODY frame.

    max_thrust_N
        Maximum thrust magnitude in Newtons.

    min_pulse_s
        Minimum pulse width this thruster can physically produce.
        For milestone 1, this is typically one whole simulation step.
    """
    name: str
    position_b_m: np.ndarray
    direction_b_unit: np.ndarray
    max_thrust_N: float
    min_pulse_s: float = constants.ORION_CM_RCS_MIN_PULSE_S

    def __post_init__(self) -> None:
        self.position_b_m = np.asarray(self.position_b_m, dtype=float).reshape(3)
        self.direction_b_unit = normalize(self.direction_b_unit)
        self.max_thrust_N = float(self.max_thrust_N)
        self.min_pulse_s = float(self.min_pulse_s)

        if self.max_thrust_N <= 0.0:
            raise ValueError(f"Thruster {self.name} must have positive max_thrust_N")
        if self.min_pulse_s <= 0.0:
            raise ValueError(f"Thruster {self.name} must have positive min_pulse_s")

    def body_force(self, duty: float) -> np.ndarray:
        """
        Return body-frame force vector for a duty in [0, 1].
        """
        duty = clamp(float(duty), 0.0, 1.0)
        return self.direction_b_unit * (self.max_thrust_N * duty)

    def body_torque(self, duty: float) -> np.ndarray:
        """
        Return body-frame torque vector for a duty in [0, 1].

        torque = r x F
        """
        return np.cross(self.position_b_m, self.body_force(duty))


@dataclass
class ThrusterFireCommand:
    """
    Per-step thruster command.

    duty_by_name maps thruster name -> duty over the current full simulation step.
    In milestone 1 the duties are usually either 0.0 or 1.0.
    """
    duty_by_name: Dict[str, float] = field(default_factory=dict)

    def is_empty(self) -> bool:
        return len(self.duty_by_name) == 0

    def active_names(self) -> List[str]:
        return [name for name, duty in self.duty_by_name.items() if duty > 0.0]


@dataclass
class RCSWrench:
    """
    Resulting body-frame wrench from a fire command.
    """
    force_b_N: np.ndarray
    torque_b_Nm: np.ndarray


@dataclass
class RollChannelStepResult:
    """
    Convenience debug packet for one roll allocation step.
    """
    tau_roll_cmd_Nm: float
    tau_roll_capacity_Nm: float
    requested_duty: float
    fired_this_step: bool
    fire_cmd: ThrusterFireCommand
    wrench: RCSWrench


# Duty accumulation for fixed-step minimum pulse logic ------------------------------------------------------------------
class PulseAccumulator:
    """
    Convert a fractional duty request into whole-step ON/OFF firing.

    Why this exists
    ---------------
    The integrator uses a fixed step of dt = 0.25 s.
    If minimum pulse width is also 0.25 s, then inside one step the jet cannot
    fire for "half a step" physically in this simple model

    So we accumulate requested duty across steps:
    - if requested duty is large enough, fire this step
    - otherwise, carry fractional demand forward
    """

    def __init__(self) -> None:
        self._accum = 0.0

    def reset(self) -> None:
        self._accum = 0.0

    def step(self, duty_request: float) -> bool:
        """
        Return True if the channel should fire for this entire simulation step
        """
        duty_request = clamp(float(duty_request), 0.0, 1.0)
        self._accum += duty_request

        if self._accum >= 1.0:
            self._accum -= 1.0
            return True

        return False

    @property
    def residual(self) -> float:
        return self._accum

# Physical RCS system ------------------------------------------------------------------

class CapsuleRCSSystem:
    """
    Physical body-frame RCS system for the entry capsule.

    Layout assumptions
    ------------------
    - Thrusters are fixed in body coordinates.
    - Body-frame force and torque are computed exactly from the thruster layout
    - The roll allocator uses symmetric tangential jets
    - Axial jets are stored now for future pitch/yaw control, even though
      milestone 1 does not allocate them yet

    Main API
    --------
    reset()
    allocate_roll_step(tau_roll_cmd_Nm, dt_s)
    body_wrench(fire_cmd)
    step_roll_channel(tau_roll_cmd_Nm, dt_s)
    """

    def __init__(
        self,
        thrusters: Dict[str, RCSThruster],
        roll_pos_names: Iterable[str],
        roll_neg_names: Iterable[str],
        axial_names: Optional[Iterable[str]] = None,
    ):
        self.thrusters: Dict[str, RCSThruster] = dict(thrusters)
        self.roll_pos_names: List[str] = list(roll_pos_names)
        self.roll_neg_names: List[str] = list(roll_neg_names)
        self.axial_names: List[str] = list(axial_names) if axial_names is not None else []

        if not self.roll_pos_names:
            raise ValueError("roll_pos_names cannot be empty")
        if not self.roll_neg_names:
            raise ValueError("roll_neg_names cannot be empty")

        missing = [
            name
            for name in (self.roll_pos_names + self.roll_neg_names + self.axial_names)
            if name not in self.thrusters
        ]
        if missing:
            raise ValueError(f"Unknown thruster names in layout: {missing}")

        # Separate accumulators for positive-roll and negative-roll channels
        self._roll_pos_pwm = PulseAccumulator()
        self._roll_neg_pwm = PulseAccumulator()

    # Lifecycle

    def reset(self) -> None:
        """
        Reset internal PWM state for a new trajectory.
        """
        self._roll_pos_pwm.reset()
        self._roll_neg_pwm.reset()

    # Basic access helpers

    def get_thruster(self, name: str) -> RCSThruster:
        return self.thrusters[name]

    def names(self) -> List[str]:
        return list(self.thrusters.keys())

    # Wrench evaluation ------------------------------------------------------------------

    def body_wrench(self, fire_cmd: ThrusterFireCommand) -> RCSWrench:
        """
        Sum the body-frame force and torque from a fire command.
        """
        force_b = np.zeros(3, dtype=float)
        torque_b = np.zeros(3, dtype=float)

        for name, duty in fire_cmd.duty_by_name.items():
            thr = self.thrusters[name]
            force_b += thr.body_force(duty)
            torque_b += thr.body_torque(duty)

        return RCSWrench(force_b_N=force_b, torque_b_Nm=torque_b)

    # Roll-channel capability ------------------------------------------------------------------

    def roll_torque_capacity_Nm(self, positive: bool = True) -> float:
        """
        Return the magnitude of the total z-axis roll torque produced when the
        whole positive-roll or negative-roll channel fires at full thrust.
        """
        names = self.roll_pos_names if positive else self.roll_neg_names

        tau_z = 0.0
        for name in names:
            thr = self.thrusters[name]
            tau_z += thr.body_torque(duty=1.0)[2]

        return abs(float(tau_z))

    # Roll allocator ------------------------------------------------------------------

    def allocate_roll_step(
        self,
        tau_roll_cmd_Nm: float,
        dt_s: float,
    ) -> ThrusterFireCommand:
        """
        Allocate body-fixed roll thrusters for one simulation step.

        Inputs
        ------
        tau_roll_cmd_Nm
            Commanded roll torque about body z.

        dt_s
            Current simulation step size.

        Behavior
        --------
        - Positive command -> positive-roll channel
        - Negative command -> negative-roll channel
        - Requested duty is abs(cmd) / channel_capacity
        - Minimum pulse width is modeled as one whole step
        - If the channel fires, all jets in that channel fire for the full step

        This is intentionally simple and stable for milestone 1.
        """
        if dt_s <= 0.0:
            raise ValueError("dt_s must be positive")

        fire = ThrusterFireCommand()

        if abs(tau_roll_cmd_Nm) <= 1.0e-12:
            return fire

        sign_positive = tau_roll_cmd_Nm > 0.0
        channel_names = self.roll_pos_names if sign_positive else self.roll_neg_names
        channel_capacity = self.roll_torque_capacity_Nm(positive=sign_positive)

        if channel_capacity <= 1.0e-12:
            return fire

        requested_duty = clamp(abs(float(tau_roll_cmd_Nm)) / channel_capacity, 0.0, 1.0)

        # Milestone-1 minimum pulse logic:
        # one step is the smallest firing unit.
        pwm = self._roll_pos_pwm if sign_positive else self._roll_neg_pwm
        fire_this_step = pwm.step(requested_duty)

        if not fire_this_step:
            return fire

        for name in channel_names:
            thr = self.thrusters[name]

            # If the thruster min pulse is larger than dt, then this simple
            # model is not consistent and should fail loudly.
            if thr.min_pulse_s - dt_s > 1.0e-12:
                raise ValueError(
                    f"Thruster {name} min_pulse_s={thr.min_pulse_s} exceeds dt_s={dt_s}"
                )

            fire.duty_by_name[name] = 1.0

        return fire

    def step_roll_channel(
        self,
        tau_roll_cmd_Nm: float,
        dt_s: float,
    ) -> RollChannelStepResult:
        """
        Convenience wrapper that:
        1. allocates roll thrusters
        2. computes the resulting wrench
        3. returns a debug-rich result packet
        """
        sign_positive = tau_roll_cmd_Nm > 0.0
        tau_capacity = 0.0
        requested_duty = 0.0

        if abs(tau_roll_cmd_Nm) > 1.0e-12:
            tau_capacity = self.roll_torque_capacity_Nm(positive=sign_positive)
            if tau_capacity > 1.0e-12:
                requested_duty = clamp(abs(float(tau_roll_cmd_Nm)) / tau_capacity, 0.0, 1.0)

        fire_cmd = self.allocate_roll_step(
            tau_roll_cmd_Nm=tau_roll_cmd_Nm,
            dt_s=dt_s,
        )
        wrench = self.body_wrench(fire_cmd)

        return RollChannelStepResult(
            tau_roll_cmd_Nm=float(tau_roll_cmd_Nm),
            tau_roll_capacity_Nm=float(tau_capacity),
            requested_duty=float(requested_duty),
            fired_this_step=not fire_cmd.is_empty(),
            fire_cmd=fire_cmd,
            wrench=wrench,
        )


# Orion-inspired milestone-1 12-thruster layout
def build_orion_cm_rcs_12(
    ring_radius_m: float = constants.CAPSULE_RCS_RING_RADIUS_M,
    ring_z_m: float = constants.CAPSULE_RCS_RING_Z_M,
    thrust_N: float = constants.ORION_CM_RCS_THRUST_N,
    min_pulse_s: float = constants.ORION_CM_RCS_MIN_PULSE_S,
) -> CapsuleRCSSystem:
    """
    Build a symmetric 12-thruster Orion-inspired crew-module layout.

    Layout
    ------
    Four pods around the shoulder ring:
    - one axial jet per pod
    - one positive-roll tangential jet per pod
    - one negative-roll tangential jet per pod

    Naming
    ------
    POD0_AX
    POD0_RPOS
    POD0_RNEG
    ...
    POD3_AX
    POD3_RPOS
    POD3_RNEG

    Body-axis convention
    --------------------
    x_b, y_b, z_b form a right-handed body frame.
    z_b is the roll axis used by milestone-1 bank control.

    Geometry
    --------
    Pod azimuths are:
    0 deg, 90 deg, 180 deg, 270 deg

    Tangential direction at pod azimuth psi is:
        t_hat = [-sin(psi), cos(psi), 0]

    This makes the RPOS channel produce positive z torque and the RNEG
    channel produce negative z torque.
    """
    thrusters: Dict[str, RCSThruster] = {}
    roll_pos_names: List[str] = []
    roll_neg_names: List[str] = []
    axial_names: List[str] = []

    azimuths = (
        0.0,
        0.5 * math.pi,
        1.0 * math.pi,
        1.5 * math.pi,
    )

    for pod_idx, psi in enumerate(azimuths):
        c = math.cos(psi)
        s = math.sin(psi)

        position_b = np.array([ring_radius_m * c, ring_radius_m * s, ring_z_m], dtype=float)

        axial_dir = np.array([0.0, 0.0, 1.0], dtype=float)
        tangential_dir = np.array([-s, c, 0.0], dtype=float)

        ax_name = f"POD{pod_idx}_AX"
        rp_name = f"POD{pod_idx}_RPOS"
        rn_name = f"POD{pod_idx}_RNEG"

        thrusters[ax_name] = RCSThruster(
            name=ax_name,
            position_b_m=position_b,
            direction_b_unit=axial_dir,
            max_thrust_N=thrust_N,
            min_pulse_s=min_pulse_s,
        )
        axial_names.append(ax_name)

        thrusters[rp_name] = RCSThruster(
            name=rp_name,
            position_b_m=position_b,
            direction_b_unit=tangential_dir,
            max_thrust_N=thrust_N,
            min_pulse_s=min_pulse_s,
        )
        roll_pos_names.append(rp_name)

        thrusters[rn_name] = RCSThruster(
            name=rn_name,
            position_b_m=position_b,
            direction_b_unit=-tangential_dir,
            max_thrust_N=thrust_N,
            min_pulse_s=min_pulse_s,
        )
        roll_neg_names.append(rn_name)

    return CapsuleRCSSystem(
        thrusters=thrusters,
        roll_pos_names=roll_pos_names,
        roll_neg_names=roll_neg_names,
        axial_names=axial_names,
    )


# compatibility helpers

@dataclass
class LegacyPulseResult:
    """
    Small compatibility packet for older logs or notebooks that still expect
    something pulse-like from an RCS object.

    This is not used by the new milestone-1 control stack directly
    """
    fired: bool
    active_names: List[str]
    requested_duty: float
    residual_pos: float
    residual_neg: float


class LegacyPulseAdapter:
    """
    Thin adapter that exposes a quantize-like interface on top of the new
    physical roll-channel allocator.

    Use this only if some surrounding code still expects an object with a
    quantize_on_time-like flavor during migration.

    New code should use:
        rcs.step_roll_channel(...)
    directly.
    """

    def __init__(self, rcs_system: CapsuleRCSSystem):
        self.rcs = rcs_system

    def reset(self) -> None:
        self.rcs.reset()

    def quantize_tau_roll(
        self,
        tau_roll_cmd_Nm: float,
        dt_s: float,
    ) -> LegacyPulseResult:
        result = self.rcs.step_roll_channel(
            tau_roll_cmd_Nm=tau_roll_cmd_Nm,
            dt_s=dt_s,
        )
        return LegacyPulseResult(
            fired=result.fired_this_step,
            active_names=result.fire_cmd.active_names(),
            requested_duty=result.requested_duty,
            residual_pos=self.rcs._roll_pos_pwm.residual,
            residual_neg=self.rcs._roll_neg_pwm.residual,
        )


# ------------------------------------------------------------------------------
# Example milestone-1 usage
# ------------------------------------------------------------------------------
#
# rcs = build_orion_cm_rcs_12()
#
# # Every 0.25 s physics step:
# roll_step = rcs.step_roll_channel(
#     tau_roll_cmd_Nm=ctrl_out.tau_roll_cmd_Nm,
#     dt_s=constants.ENTRY_DT_S,
# )
#
# F_rcs_b = roll_step.wrench.force_b_N
# tau_rcs_b = roll_step.wrench.torque_b_Nm
#
# Important:
# - control.py owns sigma_cmd and sigma_target
# - ReactionControl.py turns tau_roll_cmd into a real body wrench
# - math_3d.py integrates roll rate, attitude, and sigma_actual