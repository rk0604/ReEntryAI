"""
Reaction control system (RCS) model for the capsule roll channel.

This module turns a requested roll torque into real, timed thruster firings and
the body force and torque those firings produce. It models the parts of a small
thruster that matter for roll control: a minimum on time once a jet starts, a
minimum off time before it can fire again, a command lag, and the accumulated
demand that has not yet been delivered.

Main pieces:
    RCSThruster        one jet, with a body position, a thrust direction, a
                       thrust level, and pulse timing.
    TimedPulseChannel  converts a requested average duty into a real on or off
                       jet state, honoring the timing rules above.
    CapsuleRCSSystem   the full set of jets plus the positive and negative roll
                       channels, with the allocator that decides which jets fire.
    build_orion_cm_rcs_12  builds a symmetric four pod, twelve jet layout.

The roll axis is the body z axis. Positive roll and negative roll are handled by
two independent timed channels.
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass, field
from typing import Deque, Dict, Iterable, List, Optional

import numpy as np
import constants


# Small math helpers

def clamp(value: float, lower: float, upper: float) -> float:
    """Return value held inside the closed range [lower, upper]."""
    return max(lower, min(upper, value))


def normalize(vec: np.ndarray, eps: float = 1.0e-12) -> np.ndarray:
    """
    Return the unit length version of a three component vector.

    Parameters:
        vec: a three element array like value.
        eps: the smallest length treated as nonzero.

    Returns:
        The input divided by its length. Raises ValueError when the input is too
        close to zero to define a direction.
    """
    # Convert the input into a flat three component float array.
    arr = np.asarray(vec, dtype=float).reshape(3)

    # Compute the Euclidean length so the direction can be normalized.
    n = float(np.linalg.norm(arr))

    # A near zero vector cannot define a thrust direction.
    if n < eps:
        raise ValueError("Cannot normalize near zero vector")

    # Return the unit direction.
    return arr / n


# Core thruster and wrench data models

@dataclass
class RCSThruster:
    """
    One reaction control thruster.

    Fields:
        name:             human readable thruster name.
        position_b_m:     thruster location in body coordinates, in meters.
        direction_b_unit: unit thrust direction in body coordinates.
        max_thrust_N:     maximum thrust level in newtons.
        min_pulse_s:      minimum time the jet stays on once it starts firing.
        min_off_pulse_s:  minimum time the jet stays off before it can refire.
        lag_s:            delay before a command reaches the jet logic.
    """

    name: str
    position_b_m: np.ndarray
    direction_b_unit: np.ndarray
    max_thrust_N: float
    min_pulse_s: float = constants.ORION_CM_RCS_MIN_PULSE_S
    min_off_pulse_s: float = constants.ORION_CM_RCS_MIN_OFF_PULSE_S
    lag_s: float = constants.ORION_CM_RCS_LAG_S

    def __post_init__(self) -> None:
        """Normalize the stored vectors and validate the thrust and timing values."""
        # Store the body position as a three component float vector.
        self.position_b_m = np.asarray(self.position_b_m, dtype=float).reshape(3)

        # Normalize the thrust direction once so later force scaling is correct.
        self.direction_b_unit = normalize(self.direction_b_unit)

        # Store all timing and thrust values as floats.
        self.max_thrust_N = float(self.max_thrust_N)
        self.min_pulse_s = float(self.min_pulse_s)
        self.min_off_pulse_s = float(self.min_off_pulse_s)
        self.lag_s = float(self.lag_s)

        # Thrust magnitude must be positive.
        if self.max_thrust_N <= 0.0:
            raise ValueError(f"Thruster {self.name} must have positive max_thrust_N")

        # Minimum on time must be positive.
        if self.min_pulse_s <= 0.0:
            raise ValueError(f"Thruster {self.name} must have positive min_pulse_s")

        # Minimum off time must be zero or greater.
        if self.min_off_pulse_s < 0.0:
            raise ValueError(f"Thruster {self.name} must have nonnegative min_off_pulse_s")

        # Lag must be zero or greater.
        if self.lag_s < 0.0:
            raise ValueError(f"Thruster {self.name} must have nonnegative lag_s")

    def body_force(self, duty: float) -> np.ndarray:
        """
        Return the body frame force for a given duty.

        Parameters:
            duty: fraction of full thrust in the range zero to one.

        Returns:
            The force vector, the unit direction scaled by max thrust times duty.
        """
        # Clamp the duty so force never exceeds the physical maximum.
        duty = clamp(float(duty), 0.0, 1.0)

        # Force equals thrust magnitude times unit direction.
        return self.direction_b_unit * (self.max_thrust_N * duty)

    def body_torque(self, duty: float) -> np.ndarray:
        """
        Return the body frame torque for a given duty.

        Parameters:
            duty: fraction of full thrust in the range zero to one.

        Returns:
            The torque vector, position crossed with force.
        """
        # Torque is position crossed with force, computed by hand. np.cross has
        # large per call overhead and, called on the order of a million times
        # per episode, dominated the runtime. The explicit form is far faster.
        r = self.position_b_m
        f = self.body_force(duty)
        return np.array([
            r[1] * f[2] - r[2] * f[1],
            r[2] * f[0] - r[0] * f[2],
            r[0] * f[1] - r[1] * f[0],
        ])


@dataclass
class ThrusterFireCommand:
    """
    A firing command for one integration sub step.

    Fields:
        duty_by_name: maps each thruster name to its duty for this sub step.
    """

    duty_by_name: Dict[str, float] = field(default_factory=dict)

    def is_empty(self) -> bool:
        """Return True when no thruster is active this sub step."""
        return len(self.duty_by_name) == 0

    def active_names(self) -> List[str]:
        """Return the names of every thruster with a positive commanded duty."""
        return [name for name, duty in self.duty_by_name.items() if duty > 0.0]


@dataclass
class RCSWrench:
    """
    The net body force and torque from a set of firings.

    Fields:
        force_b_N:   net body force in newtons.
        torque_b_Nm: net body torque in newton meters.
    """

    force_b_N: np.ndarray
    torque_b_Nm: np.ndarray


@dataclass
class RollChannelStepResult:
    """
    Debug rich result of stepping the roll channel for one sub step.

    Fields:
        tau_roll_cmd_Nm:      commanded roll torque that entered the allocator.
        tau_roll_capacity_Nm: total available torque for the commanded direction.
        requested_duty:       requested average duty before pulse timing.
        fired_this_step:      True when any roll thruster fired this sub step.
        fire_cmd:             the low level thruster command.
        wrench:               the real body force and torque produced.
        roll_pos_backlog_s:   undelivered positive roll demand, in seconds.
        roll_neg_backlog_s:   undelivered negative roll demand, in seconds.
        roll_pos_is_on:       True when the positive channel is firing.
        roll_neg_is_on:       True when the negative channel is firing.
    """

    tau_roll_cmd_Nm: float
    tau_roll_capacity_Nm: float
    requested_duty: float
    fired_this_step: bool
    fire_cmd: ThrusterFireCommand
    wrench: RCSWrench
    roll_pos_backlog_s: float
    roll_neg_backlog_s: float
    roll_pos_is_on: bool
    roll_neg_is_on: bool


# Timed pulse channel model

class TimedPulseChannel:
    """
    Convert a requested average duty into a real timed jet state.

    This models command lag, minimum on time, minimum off time, and the
    accumulated demand that has not yet been delivered. One channel drives one
    roll direction.
    """

    def __init__(
        self,
        min_pulse_s: float,
        min_off_pulse_s: float,
        lag_s: float,
    ) -> None:
        """
        Build a timed pulse channel.

        Parameters:
            min_pulse_s:     minimum on time once firing begins.
            min_off_pulse_s: minimum off time after turning off.
            lag_s:           command lag before demand reaches the firing logic.
        """
        # Store timing values as floats.
        self.min_pulse_s = float(min_pulse_s)
        self.min_off_pulse_s = float(min_off_pulse_s)
        self.lag_s = float(lag_s)

        # Validate physical timing values.
        if self.min_pulse_s <= 0.0:
            raise ValueError("min_pulse_s must be positive")
        if self.min_off_pulse_s < 0.0:
            raise ValueError("min_off_pulse_s must be nonnegative")
        if self.lag_s < 0.0:
            raise ValueError("lag_s must be nonnegative")

        # Create an internal delay line that is sized on first use.
        self._delay_steps = 0
        self._delay_queue: Deque[float] = deque()

        # Initialize dynamic state.
        self.reset()

    def reset(self) -> None:
        """Reset the firing state, locks, backlog, and delay line to idle."""
        # backlog_s stores requested on time not yet delivered by the jet.
        self._backlog_s = 0.0

        # is_on marks whether the channel is firing during the current sub step.
        self._is_on = False

        # on_lock_remaining_s enforces the minimum on time once firing begins.
        self._on_lock_remaining_s = 0.0

        # off_lock_remaining_s enforces the minimum off time after turn off.
        self._off_lock_remaining_s = 0.0

        # Delay line contents are reset to zero demand.
        if self._delay_steps > 0:
            self._delay_queue = deque([0.0] * self._delay_steps, maxlen=self._delay_steps)
        else:
            self._delay_queue = deque()

    def _ensure_delay_queue(self, dt_s: float) -> None:
        """
        Validate the step size and size the lag delay line to match it.

        Parameters:
            dt_s: the sub step size in seconds. It must be positive, small
                  enough to resolve the timing values, and an exact divisor of
                  the lag.
        """
        # dt_s must be positive because the timing logic advances forward in time.
        if dt_s <= 0.0:
            raise ValueError("dt_s must be positive")

        # The internal model needs dt small enough to resolve on, off, and lag.
        smallest_positive_time = min(
            self.min_pulse_s,
            self.min_off_pulse_s if self.min_off_pulse_s > 0.0 else self.min_pulse_s,
            self.lag_s if self.lag_s > 0.0 else self.min_pulse_s,
        )

        if dt_s - smallest_positive_time > 1.0e-12:
            raise ValueError("dt_s is too large for the timed pulse channel model")

        # Convert lag into an integer number of delayed sub steps.
        new_delay_steps = int(round(self.lag_s / dt_s))

        # The chosen dt must resolve the configured lag exactly.
        if abs(new_delay_steps * dt_s - self.lag_s) > 1.0e-12:
            raise ValueError("lag_s must be an exact multiple of dt_s")

        # Only rebuild the queue if the step size changed or it is uninitialized.
        if new_delay_steps != self._delay_steps or len(self._delay_queue) != new_delay_steps:
            self._delay_steps = new_delay_steps
            self._delay_queue = deque([0.0] * self._delay_steps, maxlen=self._delay_steps)

    def _delayed_request(self, duty_request: float) -> float:
        """
        Push a new request through the lag delay line and return the active one.

        Parameters:
            duty_request: the freshly requested average duty.

        Returns:
            The request that becomes active now, which is the oldest delayed
            request, or the input itself when there is no lag.
        """
        # Without lag the current request reaches the channel immediately.
        if self._delay_steps == 0:
            return duty_request

        # With lag the oldest delayed request becomes active now.
        delayed = self._delay_queue.popleft()
        self._delay_queue.append(duty_request)
        return float(delayed)

    def step(self, duty_request: float, dt_s: float) -> bool:
        """
        Advance the channel one sub step and report whether the jet fires.

        Applies the lag, accumulates demand, and enforces the minimum on time,
        minimum off time, and the threshold to start a real pulse.

        Parameters:
            duty_request: requested average duty in the range zero to one.
            dt_s:         sub step size in seconds.

        Returns:
            True when the jet is on during this sub step, otherwise False.
        """
        # Configure and validate timing resolution for this sub step.
        self._ensure_delay_queue(float(dt_s))

        # Clamp the requested average duty into the physical range.
        duty_request = clamp(float(duty_request), 0.0, 1.0)

        # Apply the lag model before demand reaches the firing logic.
        delayed_request = self._delayed_request(duty_request)

        # Convert delayed average duty into requested on time for this sub step.
        self._backlog_s += delayed_request * float(dt_s)

        # If the channel is currently on then this sub step fires.
        if self._is_on:
            # Delivered on time reduces the remaining backlog.
            self._backlog_s = max(0.0, self._backlog_s - float(dt_s))

            # While the minimum on lock remains the channel must stay on.
            if self._on_lock_remaining_s > 0.0:
                self._on_lock_remaining_s = max(0.0, self._on_lock_remaining_s - float(dt_s))
                return True

            # After the minimum on time the channel stays on while backlog remains.
            if self._backlog_s > 1.0e-15:
                return True

            # No remaining demand, so the channel turns off after this fired sub step.
            self._is_on = False
            self._off_lock_remaining_s = self.min_off_pulse_s
            return True

        # If the channel is off then the minimum off lock must expire first.
        if self._off_lock_remaining_s > 0.0:
            self._off_lock_remaining_s = max(0.0, self._off_lock_remaining_s - float(dt_s))
            return False

        # Start firing only when enough accumulated on time justifies a real pulse.
        if self._backlog_s + 1.0e-15 >= self.min_pulse_s:
            self._is_on = True

            # The minimum on lock starts now and one sub step is consumed immediately.
            self._on_lock_remaining_s = max(0.0, self.min_pulse_s - float(dt_s))

            # This sub step delivers real on time immediately.
            self._backlog_s = max(0.0, self._backlog_s - float(dt_s))

            return True

        # Otherwise remain off and keep storing requested on time in backlog.
        return False

    @property
    def backlog_s(self) -> float:
        """Return the remaining undelivered on time, for logging."""
        return float(self._backlog_s)

    @property
    def is_on(self) -> bool:
        """Return the current firing state, for logging."""
        return bool(self._is_on)


# Physical RCS system

class CapsuleRCSSystem:
    """
    The full body frame thruster layout and the roll channel timing logic.

    Thruster forces and torques are computed exactly in body coordinates. The
    positive and negative roll channels are timed independently.
    """

    def __init__(
        self,
        thrusters: Dict[str, RCSThruster],
        roll_pos_names: Iterable[str],
        roll_neg_names: Iterable[str],
        axial_names: Optional[Iterable[str]] = None,
    ):
        """
        Build the RCS system from a thruster layout and channel membership.

        Parameters:
            thrusters:      map from thruster name to RCSThruster.
            roll_pos_names: names of the jets that make positive roll torque.
            roll_neg_names: names of the jets that make negative roll torque.
            axial_names:    optional names of axial jets, kept for future use.
        """
        # Copy the thruster layout and named channel membership.
        self.thrusters: Dict[str, RCSThruster] = dict(thrusters)
        self.roll_pos_names: List[str] = list(roll_pos_names)
        self.roll_neg_names: List[str] = list(roll_neg_names)
        self.axial_names: List[str] = list(axial_names) if axial_names is not None else []

        # The positive roll channel must exist.
        if not self.roll_pos_names:
            raise ValueError("roll_pos_names cannot be empty")

        # The negative roll channel must exist.
        if not self.roll_neg_names:
            raise ValueError("roll_neg_names cannot be empty")

        # Every named thruster must exist in the layout.
        missing = [
            name
            for name in (self.roll_pos_names + self.roll_neg_names + self.axial_names)
            if name not in self.thrusters
        ]
        if missing:
            raise ValueError(f"Unknown thruster names in layout: {missing}")

        # Build one timed pulse channel for positive roll.
        pos_timing = self._channel_timing_from_names(self.roll_pos_names)
        self._roll_pos_channel = TimedPulseChannel(
            min_pulse_s=pos_timing["min_pulse_s"],
            min_off_pulse_s=pos_timing["min_off_pulse_s"],
            lag_s=pos_timing["lag_s"],
        )

        # Build one timed pulse channel for negative roll.
        neg_timing = self._channel_timing_from_names(self.roll_neg_names)
        self._roll_neg_channel = TimedPulseChannel(
            min_pulse_s=neg_timing["min_pulse_s"],
            min_off_pulse_s=neg_timing["min_off_pulse_s"],
            lag_s=neg_timing["lag_s"],
        )

    def _channel_timing_from_names(self, names: Iterable[str]) -> Dict[str, float]:
        """
        Read the shared timing values for a roll channel.

        Parameters:
            names: thruster names that belong to one channel.

        Returns:
            A dict with min_pulse_s, min_off_pulse_s, and lag_s. Every thruster
            in the channel must share these values, otherwise ValueError is
            raised.
        """
        # Pull timing values from the first thruster in the channel.
        name_list = list(names)
        first = self.thrusters[name_list[0]]

        min_pulse_s = float(first.min_pulse_s)
        min_off_pulse_s = float(first.min_off_pulse_s)
        lag_s = float(first.lag_s)

        # All thrusters in the same roll channel must share the same timing.
        for name in name_list[1:]:
            thr = self.thrusters[name]

            if abs(float(thr.min_pulse_s) - min_pulse_s) > 1.0e-12:
                raise ValueError("All thrusters in a roll channel must share min_pulse_s")

            if abs(float(thr.min_off_pulse_s) - min_off_pulse_s) > 1.0e-12:
                raise ValueError("All thrusters in a roll channel must share min_off_pulse_s")

            if abs(float(thr.lag_s) - lag_s) > 1.0e-12:
                raise ValueError("All thrusters in a roll channel must share lag_s")

        return {
            "min_pulse_s": min_pulse_s,
            "min_off_pulse_s": min_off_pulse_s,
            "lag_s": lag_s,
        }

    def reset(self) -> None:
        """Reset the positive and negative roll channel timing states."""
        self._roll_pos_channel.reset()
        self._roll_neg_channel.reset()

    def get_thruster(self, name: str) -> RCSThruster:
        """Return one thruster by name."""
        return self.thrusters[name]

    def names(self) -> List[str]:
        """Return all thruster names in the layout."""
        return list(self.thrusters.keys())

    def body_wrench(self, fire_cmd: ThrusterFireCommand) -> RCSWrench:
        """
        Sum the body force and torque from every active thruster.

        Parameters:
            fire_cmd: the per sub step firing command.

        Returns:
            The net body wrench for the sub step.
        """
        # Initialize zero net force and zero net torque in body coordinates.
        force_b = np.zeros(3, dtype=float)
        torque_b = np.zeros(3, dtype=float)

        # Sum every active thruster contribution.
        for name, duty in fire_cmd.duty_by_name.items():
            thr = self.thrusters[name]
            force_b += thr.body_force(duty)
            torque_b += thr.body_torque(duty)

        # Return the net body wrench for the current sub step.
        return RCSWrench(force_b_N=force_b, torque_b_Nm=torque_b)

    def roll_torque_capacity_Nm(self, positive: bool = True) -> float:
        """
        Return the total roll torque available in one direction.

        Parameters:
            positive: True for the positive roll channel, False for negative.

        Returns:
            The magnitude of the summed z axis torque when every jet in the
            channel fires at full thrust.
        """
        # Choose the positive or negative roll channel.
        names = self.roll_pos_names if positive else self.roll_neg_names

        # Sum the z axis torque when every jet in the channel fires at full thrust.
        tau_z = 0.0
        for name in names:
            thr = self.thrusters[name]
            tau_z += thr.body_torque(duty=1.0)[2]

        # Capacity is the magnitude of the signed z torque sum.
        return abs(float(tau_z))

    def allocate_roll_step(
        self,
        tau_roll_cmd_Nm: float,
        dt_s: float,
    ) -> ThrusterFireCommand:
        """
        Decide which jets fire this sub step for a commanded roll torque.

        Parameters:
            tau_roll_cmd_Nm: commanded roll torque. Its sign selects the channel.
            dt_s:            sub step size in seconds.

        Returns:
            A firing command. When the active channel fires, every jet in that
            channel is set to full duty.
        """
        # The time step must be positive for timing state propagation.
        if dt_s <= 0.0:
            raise ValueError("dt_s must be positive")

        # The default output is no active thruster.
        fire = ThrusterFireCommand()

        # Compute the requested average duty for the commanded direction.
        requested_duty = 0.0
        sign_positive = tau_roll_cmd_Nm > 0.0

        if abs(tau_roll_cmd_Nm) > 1.0e-12:
            channel_capacity = self.roll_torque_capacity_Nm(positive=sign_positive)

            if channel_capacity > 1.0e-12:
                requested_duty = clamp(abs(float(tau_roll_cmd_Nm)) / channel_capacity, 0.0, 1.0)

        # Feed the proper request into the positive and negative timed channels.
        if sign_positive and requested_duty > 0.0:
            fire_pos = self._roll_pos_channel.step(requested_duty, float(dt_s))
            fire_neg = self._roll_neg_channel.step(0.0, float(dt_s))
        elif (not sign_positive) and requested_duty > 0.0:
            fire_pos = self._roll_pos_channel.step(0.0, float(dt_s))
            fire_neg = self._roll_neg_channel.step(requested_duty, float(dt_s))
        else:
            fire_pos = self._roll_pos_channel.step(0.0, float(dt_s))
            fire_neg = self._roll_neg_channel.step(0.0, float(dt_s))

        # If the positive roll channel is active then all positive roll jets fire fully.
        if fire_pos:
            for name in self.roll_pos_names:
                fire.duty_by_name[name] = 1.0

        # If the negative roll channel is active then all negative roll jets fire fully.
        if fire_neg:
            for name in self.roll_neg_names:
                fire.duty_by_name[name] = 1.0

        return fire

    def step_roll_channel(
        self,
        tau_roll_cmd_Nm: float,
        dt_s: float,
    ) -> RollChannelStepResult:
        """
        Allocate firings and compute the resulting wrench for one sub step.

        Parameters:
            tau_roll_cmd_Nm: commanded roll torque.
            dt_s:            sub step size in seconds.

        Returns:
            A RollChannelStepResult with the command, the available capacity,
            the requested duty, the firing command, the produced wrench, and
            channel timing telemetry.
        """
        # Determine which signed capacity applies to the current command.
        sign_positive = tau_roll_cmd_Nm > 0.0
        tau_capacity = 0.0
        requested_duty = 0.0

        if abs(tau_roll_cmd_Nm) > 1.0e-12:
            tau_capacity = self.roll_torque_capacity_Nm(positive=sign_positive)

            if tau_capacity > 1.0e-12:
                requested_duty = clamp(abs(float(tau_roll_cmd_Nm)) / tau_capacity, 0.0, 1.0)

        # Allocate real timed thruster firing for this sub step.
        fire_cmd = self.allocate_roll_step(
            tau_roll_cmd_Nm=tau_roll_cmd_Nm,
            dt_s=dt_s,
        )

        # Convert the active thruster set into a physical body wrench.
        wrench = self.body_wrench(fire_cmd)

        # Return debug rich timing and wrench data.
        return RollChannelStepResult(
            tau_roll_cmd_Nm=float(tau_roll_cmd_Nm),
            tau_roll_capacity_Nm=float(tau_capacity),
            requested_duty=float(requested_duty),
            fired_this_step=not fire_cmd.is_empty(),
            fire_cmd=fire_cmd,
            wrench=wrench,
            roll_pos_backlog_s=float(self._roll_pos_channel.backlog_s),
            roll_neg_backlog_s=float(self._roll_neg_channel.backlog_s),
            roll_pos_is_on=bool(self._roll_pos_channel.is_on),
            roll_neg_is_on=bool(self._roll_neg_channel.is_on),
        )


# Orion inspired twelve thruster layout

def build_orion_cm_rcs_12(
    ring_radius_m: float = constants.CAPSULE_RCS_RING_RADIUS_M,
    ring_z_m: float = constants.CAPSULE_RCS_RING_Z_M,
    thrust_N: float = constants.ORION_CM_RCS_THRUST_N,
    min_pulse_s: float = constants.ORION_CM_RCS_MIN_PULSE_S,
    min_off_pulse_s: float = constants.ORION_CM_RCS_MIN_OFF_PULSE_S,
    lag_s: float = constants.ORION_CM_RCS_LAG_S,
) -> CapsuleRCSSystem:
    """
    Build a symmetric four pod, twelve jet crew module RCS layout.

    Each pod sits on the shoulder ring and carries three jets: one axial jet
    along body z, one positive roll jet, and one negative roll jet. The roll
    jets point tangentially so they make a torque about the body z axis.

    Parameters:
        ring_radius_m:   radius of the shoulder ring from the body origin.
        ring_z_m:        height of the ring above the body origin.
        thrust_N:        thrust of each jet in newtons.
        min_pulse_s:     minimum on time for each jet.
        min_off_pulse_s: minimum off time for each jet.
        lag_s:           command lag for each jet.

    Returns:
        A CapsuleRCSSystem with the twelve jets and the two roll channels.
    """
    # Build a symmetric crew module style layout with four pods. Each pod has
    # one axial jet, one positive roll jet, and one negative roll jet.

    thrusters: Dict[str, RCSThruster] = {}
    roll_pos_names: List[str] = []
    roll_neg_names: List[str] = []
    axial_names: List[str] = []

    # Pod azimuths are equally spaced around the shoulder ring.
    azimuths = tuple(constants.CAPSULE_RCS_POD_AZIMUTHS_RAD)

    for pod_idx, psi in enumerate(azimuths):
        # Convert azimuth into planar cosine and sine values.
        c = math.cos(psi)
        s = math.sin(psi)

        # Pod position lies on the shoulder ring in body coordinates.
        position_b = np.array([ring_radius_m * c, ring_radius_m * s, ring_z_m], dtype=float)

        # The axial jet points along positive body z.
        axial_dir = np.array([0.0, 0.0, 1.0], dtype=float)

        # The tangential direction at this azimuth makes roll torque about body z.
        tangential_dir = np.array([-s, c, 0.0], dtype=float)

        # Build consistent pod thruster names.
        ax_name = f"POD{pod_idx}_AX"
        rp_name = f"POD{pod_idx}_RPOS"
        rn_name = f"POD{pod_idx}_RNEG"

        # Store one axial thruster for future pitch, yaw, or translation use.
        thrusters[ax_name] = RCSThruster(
            name=ax_name,
            position_b_m=position_b,
            direction_b_unit=axial_dir,
            max_thrust_N=thrust_N,
            min_pulse_s=min_pulse_s,
            min_off_pulse_s=min_off_pulse_s,
            lag_s=lag_s,
        )
        axial_names.append(ax_name)

        # The positive roll jet uses the positive tangential direction.
        thrusters[rp_name] = RCSThruster(
            name=rp_name,
            position_b_m=position_b,
            direction_b_unit=tangential_dir,
            max_thrust_N=thrust_N,
            min_pulse_s=min_pulse_s,
            min_off_pulse_s=min_off_pulse_s,
            lag_s=lag_s,
        )
        roll_pos_names.append(rp_name)

        # The negative roll jet uses the opposite tangential direction.
        thrusters[rn_name] = RCSThruster(
            name=rn_name,
            position_b_m=position_b,
            direction_b_unit=-tangential_dir,
            max_thrust_N=thrust_N,
            min_pulse_s=min_pulse_s,
            min_off_pulse_s=min_off_pulse_s,
            lag_s=lag_s,
        )
        roll_neg_names.append(rn_name)

    # Return the full physical RCS system.
    return CapsuleRCSSystem(
        thrusters=thrusters,
        roll_pos_names=roll_pos_names,
        roll_neg_names=roll_neg_names,
        axial_names=axial_names,
    )


# Compatibility helpers

@dataclass
class LegacyPulseResult:
    """
    Simple result packet for older logging code.

    Fields:
        fired:          True when any jet fired this step.
        active_names:   names of the jets that fired.
        requested_duty: the requested average duty.
        residual_pos:   leftover positive roll demand as a pulse fraction.
        residual_neg:   leftover negative roll demand as a pulse fraction.
    """

    fired: bool
    active_names: List[str]
    requested_duty: float
    residual_pos: float
    residual_neg: float


class LegacyPulseAdapter:
    """
    Thin adapter that exposes older pulse style telemetry on top of the new
    timed channel model.
    """

    def __init__(self, rcs_system: CapsuleRCSSystem):
        """Store the physical RCS system this adapter wraps."""
        self.rcs = rcs_system

    def reset(self) -> None:
        """Reset the physical channel timing states."""
        self.rcs.reset()

    def quantize_tau_roll(
        self,
        tau_roll_cmd_Nm: float,
        dt_s: float,
    ) -> LegacyPulseResult:
        """
        Step the real timed roll channel and return older style telemetry.

        Parameters:
            tau_roll_cmd_Nm: commanded roll torque.
            dt_s:            sub step size in seconds.

        Returns:
            A LegacyPulseResult with the firing flag, active jet names, the
            requested duty, and the leftover demand on each channel expressed as
            a fraction of the minimum pulse.
        """
        # Step the real timed roll channel model.
        result = self.rcs.step_roll_channel(
            tau_roll_cmd_Nm=tau_roll_cmd_Nm,
            dt_s=dt_s,
        )

        # Convert remaining backlog into a pulse fraction style residual.
        pos_min_pulse = max(constants.ORION_CM_RCS_MIN_PULSE_S, 1.0e-12)
        neg_min_pulse = max(constants.ORION_CM_RCS_MIN_PULSE_S, 1.0e-12)

        return LegacyPulseResult(
            fired=result.fired_this_step,
            active_names=result.fire_cmd.active_names(),
            requested_duty=result.requested_duty,
            residual_pos=float(self.rcs._roll_pos_channel.backlog_s / pos_min_pulse),
            residual_neg=float(self.rcs._roll_neg_channel.backlog_s / neg_min_pulse),
        )
