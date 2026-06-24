"""
control.py

Milestone 1 capsule control stack for atmospheric entry.

This module owns the high level control path:

raw vehicle state
to observation features
to guidance update
to smoothed bank target
to roll torque command
to physical RCS and attitude layer elsewhere
to actual bank angle returned from physics

The control stack does not directly fly sigma itself

Instead guidance computes sigma_cmd.
The bank actuator smooths sigma_cmd into sigma_target
The roll torque controller compares sigma_target against the actual
bank angle sigma_actual.
The physical RCS plus roll dynamics layer is responsible for turning
that torque command into real motion

This revised version keeps that architecture intact but upgrades the
guidance law into a short horizon heat aware predictor corrector

Important milestone 1 approximation note

This is not a literal full reproduction of Lu.
The current implementation uses a discrete candidate search over bank
magnitudes as the predictor corrector mechanism because that fits the
current codebase cleanly and preserves the current guidance interface.

The predictor corrector uses the corrected spherical Earth 3DOF model
and the current interval heat shield model as the predictive
feasibility layer.
"""

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Protocol

import AtmosphereModel
import constants
from math_3d import (
    IntervalSupervisorConfig,
    rollout_predictor_corrector_candidate,
    rollout_nominal_to_terminal,
)


# Small math helpers ------------------------------------------------------------

def wrap_to_pi(angle_rad: float) -> float:
    """
    Wrap an angle into the closed interval from negative pi to pi.
    """
    while angle_rad > math.pi:
        angle_rad -= 2.0 * math.pi
    while angle_rad < -math.pi:
        angle_rad += 2.0 * math.pi
    return angle_rad


def wrap_longitude(lam_rad: float) -> float:
    """
    Wrap longitude into the interval from negative pi to pi.
    """
    return wrap_to_pi(lam_rad)


def clamp(value: float, lower: float, upper: float) -> float:
    """
    Clamp a scalar into the closed interval from lower to upper.
    """
    return max(lower, min(upper, value))


def sign_no_zero(x: float) -> float:
    """
    Return plus one for positive, minus one for negative, and zero for zero.
    """
    if x > 0.0:
        return 1.0
    if x < 0.0:
        return -1.0
    return 0.0


def clamp_unit(x: float) -> float:
    """
    Clamp a scalar into the interval from negative one to one.
    This is used before inverse trigonometric functions.
    """
    return clamp(x, -1.0, 1.0)


def chi_to_psi_north_rad(chi_rad: float) -> float:
    """
    Convert the codebase heading convention into a north referenced azimuth.

    The current codebase keeps chi measured from east toward north.
    That convention is preserved here and is not silently changed.

    chi = 0 means eastbound
    chi = pi / 2 means northbound
    chi = pi means westbound
    chi = negative pi / 2 means southbound

    psi is the more common azimuth measured from north toward east.
    """
    return wrap_to_pi(0.5 * math.pi - chi_rad)


def psi_north_to_chi_rad(psi_rad: float) -> float:
    """
    Convert a north referenced azimuth into the codebase chi convention.
    """
    return wrap_to_pi(0.5 * math.pi - psi_rad)


def spherical_range_and_course(
    phi_rad: float,
    lam_rad: float,
    target_phi_rad: float,
    target_lam_rad: float,
    radius_m: float,
) -> Dict[str, float]:
    """
    Compute great circle geometry from the current location to the target.

    Returned values

    central_angle_rad
    Great circle arc angle between current point and target.

    range_to_go_m
    Great circle arc length on the spherical Earth model.

    target_azimuth_north_rad
    Initial great circle azimuth from the current point to the target.
    This is measured from north toward east.

    target_course_chi_rad
    The same direction expressed in the codebase chi convention.

    north_error_m and east_error_m
    Local tangent plane components of the target direction.
    """
    # Compute wrapped longitude difference first so path geometry is continuous.
    dlam = wrap_longitude(target_lam_rad - lam_rad)

    # Precompute trigonometric terms for the spherical cosine law.
    sin_phi = math.sin(phi_rad)
    cos_phi = math.cos(phi_rad)
    sin_target_phi = math.sin(target_phi_rad)
    cos_target_phi = math.cos(target_phi_rad)

    # Compute the central angle between the two surface points.
    cos_delta = (
        sin_phi * sin_target_phi
        + cos_phi * cos_target_phi * math.cos(dlam)
    )
    central_angle_rad = math.acos(clamp_unit(cos_delta))

    # Convert the arc angle into a surface range using the spherical Earth radius.
    range_to_go_m = radius_m * central_angle_rad

    # If the current point is already at the target, return a zero direction.
    if central_angle_rad <= 1.0e-12:
        target_azimuth_north_rad = 0.0
        target_course_chi_rad = 0.0
        north_error_m = 0.0
        east_error_m = 0.0
    else:
        # Compute the initial great circle azimuth from the current point.
        y = math.sin(dlam) * cos_target_phi
        x = (
            cos_phi * sin_target_phi
            - sin_phi * cos_target_phi * math.cos(dlam)
        )
        target_azimuth_north_rad = wrap_to_pi(math.atan2(y, x))

        # Convert the north referenced azimuth into the codebase chi convention.
        target_course_chi_rad = psi_north_to_chi_rad(target_azimuth_north_rad)

        # Resolve the range into local north and east components.
        north_error_m = range_to_go_m * math.cos(target_azimuth_north_rad)
        east_error_m = range_to_go_m * math.sin(target_azimuth_north_rad)

    return {
        "central_angle_rad": float(central_angle_rad),
        "range_to_go_m": float(range_to_go_m),
        "target_azimuth_north_rad": float(target_azimuth_north_rad),
        "target_course_chi_rad": float(target_course_chi_rad),
        "north_error_m": float(north_error_m),
        "east_error_m": float(east_error_m),
    }

# Core state and config ------------------------------------------------------------

@dataclass
class ReentryState:
    """
    Minimal control facing reentry state.

    Primary translational states

    r_m
    Radial distance from Earth center in meters.

    phi_rad
    Geocentric latitude in radians.

    lam_rad
    Longitude in radians.

    V_mps
    Speed in meters per second.

    gamma_rad
    Flight path angle in radians.

    chi_rad
    Heading angle in the local horizontal plane.
    This codebase keeps chi measured from east toward north.

    Attitude and control fields

    sigma_actual_rad
    Actual bank angle currently being flown by the capsule.
    This comes back from the rigid body plus RCS layer.

    roll_rate_rad_s
    Actual roll rate of the capsule.

    sigma_cmd_rad
    Latest ideal bank command issued by guidance.

    sigma_target_rad
    Smoothed target bank angle after actuator shaping.
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

    @property
    def sigma_rad(self) -> float:
        """
        Compatibility alias for the actual flown bank angle.
        """
        return self.sigma_actual_rad

    @sigma_rad.setter
    def sigma_rad(self, value: float) -> None:
        """Set the actual flown bank angle through the compatibility alias."""
        self.sigma_actual_rad = float(value)

    @property
    def sigma_dot_rps(self) -> float:
        """
        Compatibility alias for the actual roll rate.
        """
        return self.roll_rate_rad_s

    @sigma_dot_rps.setter
    def sigma_dot_rps(self, value: float) -> None:
        """Set the actual roll rate through the compatibility alias."""
        self.roll_rate_rad_s = float(value)


@dataclass
class SimConfig:
    """
    Basic simulation timing configuration.
    """
    dt_s: float = constants.ENTRY_DT_S
    max_steps: int = 2000
    guidance_period_s: float = constants.GUIDANCE_PERIOD_S


# Predictor corrector support models ------------------------------------------------------------
@dataclass
class GuidancePredictionContext:
    """
    Runtime context passed into guidance before a guidance update.

    This gives the predictor corrector access to the live interval state
    and the live heat shield state without changing the public control
    stack call signature.
    """
    x_interval_old: Optional[List[Any]] = None
    heat_shield: Optional[Any] = None
    supervisor_cfg: Optional[IntervalSupervisorConfig] = None
    trajectory_id: str = ""
    guidance_cycle_index: int = -1
    interval_active: bool = True


# Protocols still a work in progress

class GuidanceLaw(Protocol):
    """
    Structural interface for any bank guidance module.
    """

    def reset(self) -> None:
        """Reset any internal guidance state. Part of the interface."""
        ...

    def compute_sigma_cmd(
        self,
        t_s: float,
        state: ReentryState,
        obs: Dict[str, float],
    ) -> float:
        """Return the commanded bank angle in radians. Part of the interface."""
        ...


class ObservationProvider(Protocol):
    """
    Structural interface for any module that converts raw vehicle state
    into guidance features.
    """

    def reset(self) -> None:
        """Reset any internal observation state. Part of the interface."""
        ...

    def observe(self, state: ReentryState, t_s: float) -> Dict[str, float]:
        """Return the guidance feature dictionary for a state. Part of the interface."""
        ...


# Guidance law ------------------------------------------------------------
@dataclass
class SimpleBankGuidance:
    """
    Heat aware short horizon predictor corrector guidance.

    This class provides the guidance layer that outputs sigma_cmd for the
    existing control stack.

    The guidance logic follows a short horizon candidate search structure
    inspired by predictor corrector style entry guidance methods such as Lu,
    while remaining a discrete candidate implementation for milestone 1.

    At each guidance update, the algorithm computes the current spherical
    target geometry from the observation state. It then determines the bank
    sign from the heading offset logic while preserving the current chi
    convention used elsewhere in the simulation.

    A set of candidate bank magnitudes is generated within the allowed bank
    limits. Each candidate is propagated forward through the corrected
    translational model over a short horizon. During the same rollout, the
    interval heat shield state is also propagated so that thermal feasibility
    is evaluated together with path shaping.

    Each rollout is scored using a geometry tracking cost together with a
    heat feasibility penalty. Candidates that improve target geometry while
    remaining thermally feasible are preferred. If all candidates violate the
    heat constraint, the guidance selects the least severe option according
    to the combined cost.

    This implementation is informed by Lu style predictor corrector guidance
    concepts and related entry guidance literature, but it does not perform
    the continuous scalar corrector solve described in Lu Eq 19 through Eq 25.
    Instead, it uses a discrete bank candidate evaluation over a short
    prediction horizon.
    """

    # Spherical target location used by both observation and predictor logic.
    target_phi_rad: float = 0.0
    target_lam_rad: float = 0.0

    # Aerodynamic and mass parameters used by the predictive rollout.
    # If this is left as None, the class builds a practical fallback set.
    params: Optional[Dict[str, float]] = None

    # Allowed bank range in degrees.
    max_bank_deg: float = 70.0
    min_bank_deg: float = 40.0

    # This gain is retained for compatibility with the old class even though
    # the short horizon rollout now handles most of the correction logic.
    downrange_gain: float = 0.0

    # Guidance can be disabled above a threshold if desired.
    velocity_enable_mps: float = 6000.0

    # Simple heading deadband schedule retained from the older version.
    low_speed_deadband_deg: float = 3.0
    high_speed_deadband_deg: float = 12.0
    low_speed_mps: float = 1500.0
    high_speed_mps: float = 7800.0

    # Prediction horizon settings (short interval/heat screen).
    predictor_horizon_steps: int = constants.PREDICTOR_CORRECTOR_HORIZON_STEPS
    prediction_dt_s: float = constants.ENTRY_DT_S

    # --- Continuous predictor-corrector (Lu 2014) ---
    # When enabled, the bank-angle MAGNITUDE is found by a golden-section
    # search that drives the predicted terminal miss distance to a minimum,
    # using a long nominal rollout all the way to drogue-deploy altitude.
    # The bank SIGN still comes from the existing heading-offset logic.
    # Set to False to fall back to the legacy discrete-bin selector.
    use_continuous_predictor: bool = True
    # Altitude at which the entry phase ends and the long rollout stops.
    terminal_prediction_altitude_m: float = 7600.0
    # Coarse timestep for the long forward prediction (s).
    long_prediction_dt_s: float = 2.0
    # Hard cap on long-rollout Euler steps (safety against skip-out).
    long_prediction_max_steps: int = 600
    # Upper bound of the bank-magnitude search (deg).
    sigma_search_max_deg: float = 90.0
    # Golden-section iterations (≈ this many + 2 rollout evals per cycle).
    sigma_search_iters: int = 16

    # --- Bank-reversal logic (Apollo / Lu 2014) ---
    # The bank MAGNITUDE comes from the predictor; the SIGN is set by a
    # reversal logic with a velocity-dependent deadband. The sign is HELD
    # until the heading error crosses the deadband on the opposite side,
    # then reversed -- hysteresis that prevents the sign chatter a naive
    # sign(heading_error) rule produces. Without this the bank sign flips at
    # every zero-crossing, slewing the actuator through large reversals and
    # bleeding range erratically. Set False to fall back to the naive rule.
    use_bank_reversal: bool = True

    # Candidate bank grid in degrees.
    candidate_bank_deg: List[float] = field(
        default_factory=lambda: list(constants.PREDICTOR_CANDIDATE_BANK_DEG)
    )

    # Cost weights for the geometry portion of the predictor.
    weight_range: float = constants.PREDICTOR_COST_WEIGHT_RANGE
    weight_heading: float = constants.PREDICTOR_COST_WEIGHT_HEADING
    weight_cross_track: float = constants.PREDICTOR_COST_WEIGHT_CROSS_TRACK
    weight_heat: float = constants.PREDICTOR_COST_WEIGHT_HEAT

    # Heat feasibility limits.
    heat_rate_limit: Optional[float] = constants.HEAT_RATE_LIMIT_DEFAULT
    heat_load_limit: Optional[float] = constants.HEAT_LOAD_LIMIT_DEFAULT

    # Large penalties used when a candidate violates heat feasibility or
    # interval validity.
    infeasible_penalty: float = constants.PREDICTOR_HEAT_INFEASIBLE_PENALTY
    invalid_interval_penalty: float = constants.PREDICTOR_INVALID_INTERVAL_PENALTY

    # Optional zero bank candidate.
    include_zero_bank_candidate: bool = True

    # Internal state used for live context handoff and debug output.
    _prediction_context: GuidancePredictionContext = field(
        default_factory=GuidancePredictionContext,
        init=False,
        repr=False,
    )
    _last_debug: Dict[str, Any] = field(
        default_factory=dict,
        init=False,
        repr=False,
    )
    # Held bank sign for the reversal logic (0 = not yet set).
    _current_bank_sign: float = field(default=0.0, init=False, repr=False)

    def reset(self) -> None:
        """
        Reset internal context and last debug packet.
        """
        self._prediction_context = GuidancePredictionContext()
        self._last_debug = {}
        self._current_bank_sign = 0.0

    def set_prediction_context(
        self,
        x_interval_old: Optional[List[Any]] = None,
        heat_shield: Optional[Any] = None,
        supervisor_cfg: Optional[IntervalSupervisorConfig] = None,
        trajectory_id: str = "",
        guidance_cycle_index: int = -1,
        interval_active: bool = True,
    ) -> None:
        """
        Store the live interval and heat context for the next guidance update.

        The notebook or outer simulation loop should call this before the
        control stack step when guidance is supposed to evaluate heat aware
        candidate rollouts using the current interval state.
        """
        self._prediction_context = GuidancePredictionContext(
            x_interval_old=x_interval_old,
            heat_shield=heat_shield,
            supervisor_cfg=supervisor_cfg,
            trajectory_id=str(trajectory_id),
            guidance_cycle_index=int(guidance_cycle_index),
            interval_active=bool(interval_active),
        )

    def get_last_debug(self) -> Dict[str, Any]:
        """
        Return the most recent guidance debug dictionary.
        """
        return dict(self._last_debug)

    def heading_deadband_rad(self, speed_mps: float) -> float:
        """
        Interpolate a velocity dependent heading deadband.
        """
        # If the schedule is malformed, fall back to the high speed value.
        if self.high_speed_mps <= self.low_speed_mps:
            return math.radians(self.high_speed_deadband_deg)

        # Compute interpolation coordinate across the speed schedule.
        alpha = (speed_mps - self.low_speed_mps) / (self.high_speed_mps - self.low_speed_mps)

        # Clamp the interpolation coordinate to the valid range.
        alpha = clamp(alpha, 0.0, 1.0)

        # Linearly interpolate the deadband in degrees.
        deadband_deg = (
            self.low_speed_deadband_deg
            + alpha * (self.high_speed_deadband_deg - self.low_speed_deadband_deg)
        )

        # Convert the scheduled deadband into radians.
        return math.radians(deadband_deg)

    def _default_rollout_params(self) -> Dict[str, float]:
        """
        Build a conservative fallback parameter dictionary.

        This is only used if no explicit params were injected into guidance.
        """
        return {
            "mass_kg": float(constants.CAPSULE_MASS_KG),
            "ref_area_m2": float(constants.CAPSULE_REFERENCE_AREA_M2),
            "aero_model": "scheduled_orion_like",
            "CD_const": 1.05,
            "CL_const": 0.32,
        }

    def _rollout_params(self) -> Dict[str, float]:
        """
        Return the parameter dictionary used by predictive rollouts.
        """
        if self.params is not None:
            return dict(self.params)
        return self._default_rollout_params()

    def _state_to_vector(self, state: ReentryState) -> List[float]:
        """
        Convert the control facing state object into the 3DOF translational vector.
        """
        return [
            float(state.r_m),
            float(state.phi_rad),
            float(state.lam_rad),
            float(state.V_mps),
            float(state.gamma_rad),
            float(state.chi_rad),
        ]

    def _bank_sign_from_heading_offset(
        self,
        heading_error_rad: float,
        heading_offset_psi_rad: float,
    ) -> float:
        """
        Determine the bank sign while preserving the current chi convention.

        Important convention note

        Positive sigma in the translational model produces a right turn.
        The codebase heading error is defined as target course minus current chi.
        In the older law the sign was negative sign of heading error.

        The psi based heading offset used here is defined in the observation
        layer as current heading psi minus target azimuth psi.

        With the current chi convention this means:

        heading_offset_psi = current_psi minus target_psi
        heading_error_chi = target_chi minus current_chi

        Because psi = pi over 2 minus chi, the two quantities have the same
        magnitude but opposite sign. Therefore the correct bank sign is the
        negative sign of heading_offset_psi, not the positive sign.

        This preserves the older negative sign of heading_error rule while
        keeping the current chi convention unchanged.
        """
        # Prefer the psi based heading offset rule when it is nonzero.
        # Negative heading_offset_psi means the target lies to the right
        # of the current heading, so the commanded bank must be positive.
        bank_sign = -sign_no_zero(float(heading_offset_psi_rad))

        # Fall back to the old heading error sign rule if needed.
        if bank_sign == 0.0:
            bank_sign = -sign_no_zero(float(heading_error_rad))

        # If both are exactly zero, return zero so a zero bank candidate can win.
        return float(bank_sign)

    def _bank_sign_with_reversal(
        self,
        heading_offset_psi_rad: float,
        heading_error_rad: float,
        V_mps: float,
    ) -> float:
        """
        Apollo / Lu 2014 bank-reversal logic with a velocity-dependent deadband.

        The bank SIGN steers crossrange. Rather than recompute the sign from
        the instantaneous heading error every cycle (which flips at every
        zero-crossing -> chatter -> large actuator reversals -> erratic range
        bleed), we HOLD the current sign and only reverse when the heading
        error crosses the deadband on the opposite side. The deadband is wide
        at high speed (few reversals) and narrows near the end (tight
        tracking), reusing the existing velocity schedule.

            heading_offset_psi >  +deadband  -> target decisively left  -> s = -1
            heading_offset_psi <  -deadband  -> target decisively right -> s = +1
            |heading_offset_psi| <= deadband -> HOLD current sign (hysteresis)
        """
        deadband = self.heading_deadband_rad(float(V_mps))
        psi = float(heading_offset_psi_rad)

        if psi > deadband:
            self._current_bank_sign = -1.0
        elif psi < -deadband:
            self._current_bank_sign = +1.0
        # else: inside the deadband -> hold the current sign (no reversal)

        # First cycle (or if we start inside the deadband): seed from the
        # instantaneous best sign so we don't sit at zero bank.
        if self._current_bank_sign == 0.0:
            self._current_bank_sign = self._bank_sign_from_heading_offset(
                heading_error_rad=heading_error_rad,
                heading_offset_psi_rad=psi,
            )
        return float(self._current_bank_sign)

    def _candidate_magnitude_list_deg(self, heading_error_rad: float, deadband_rad: float) -> List[float]:
        """
        Build the candidate magnitude list in degrees.

        The main predictor candidate set lives between min_bank_deg and max_bank_deg.
        A zero bank candidate can also be included for coasting behavior.
        """
        # Start from the configured candidate list.
        mags_deg: List[float] = []

        # Clamp every configured candidate into the valid bank range.
        for bank_deg in self.candidate_bank_deg:
            bank_deg_clamped = clamp(float(bank_deg), float(self.min_bank_deg), float(self.max_bank_deg))

            # Avoid duplicate entries caused by clamping.
            if all(abs(bank_deg_clamped - existing) > 1.0e-12 for existing in mags_deg):
                mags_deg.append(bank_deg_clamped)

        # If the configured candidate set was empty, fall back to the bank limits.
        if not mags_deg:
            mags_deg = [float(self.min_bank_deg), float(self.max_bank_deg)]

        # Sort the candidate magnitudes so logs are stable and readable.
        mags_deg = sorted(mags_deg)

        # Add a zero bank candidate if requested.
        # This is useful when heading error is already inside the deadband or
        # when a coasting trajectory is geometrically and thermally preferable.
        if self.include_zero_bank_candidate or abs(float(heading_error_rad)) <= float(deadband_rad):
            mags_deg = [0.0] + mags_deg

        return mags_deg

    def _predicted_terminal_errors(self, final_state: List[float]) -> Dict[str, float]:
        """
        Compute terminal range, heading, and cross track values for a rollout.
        """
        # Evaluate the same spherical geometry used by the observation layer.
        spherical = spherical_range_and_course(
            phi_rad=float(final_state[1]),
            lam_rad=float(final_state[2]),
            target_phi_rad=float(self.target_phi_rad),
            target_lam_rad=float(self.target_lam_rad),
            radius_m=float(constants.RADIUS_EARTH),
        )

        # Pull the predicted target course in the codebase chi convention.
        target_course_chi_rad = float(spherical["target_course_chi_rad"])

        # Compare predicted heading against predicted target course.
        heading_error_rad = wrap_to_pi(target_course_chi_rad - float(final_state[5]))

        # Resolve the predicted range into along track and cross track components.
        along_track_error_m = float(spherical["range_to_go_m"]) * math.cos(heading_error_rad)
        cross_track_error_m = float(spherical["range_to_go_m"]) * math.sin(heading_error_rad)

        return {
            "range_to_go_m": float(spherical["range_to_go_m"]),
            "heading_error_rad": float(heading_error_rad),
            "cross_track_error_m": float(cross_track_error_m),
            "along_track_error_m": float(along_track_error_m),
            "target_course_chi_rad": float(target_course_chi_rad),
            "target_azimuth_north_rad": float(spherical["target_azimuth_north_rad"]),
        }

    def _score_rollout_geometry(self, rollout_terminal: Dict[str, float]) -> float:
        """
        Compute the purely geometric portion of the candidate cost.
        """
        # Range cost rewards candidates that reduce predicted range to go.
        range_term = self.weight_range * abs(float(rollout_terminal["range_to_go_m"]))

        # Heading term rewards candidates that reduce terminal heading error.
        heading_term = self.weight_heading * abs(float(rollout_terminal["heading_error_rad"]))

        # Cross track term rewards candidates that reduce terminal cross track error.
        cross_track_term = self.weight_cross_track * abs(float(rollout_terminal["cross_track_error_m"]))

        # Return the total geometric cost.
        return float(range_term + heading_term + cross_track_term)

    # ------------------------------------------------------------------
    # Continuous bank-magnitude solver (Lu 2014, Fixes 1 + 2)
    # ------------------------------------------------------------------
    @staticmethod
    def _golden_section_minimize(f, a: float, b: float, iters: int) -> float:
        """Minimize a unimodal scalar function f on [a, b] by golden section.

        Returns the argmin. Uses ~iters + 2 evaluations of f.
        """
        invphi = (math.sqrt(5.0) - 1.0) / 2.0    # 0.6180339887
        invphi2 = (3.0 - math.sqrt(5.0)) / 2.0   # 0.3819660113
        c = a + invphi2 * (b - a)
        d = a + invphi * (b - a)
        fc = f(c)
        fd = f(d)
        for _ in range(int(iters)):
            if fc < fd:
                b, d, fd = d, c, fc
                c = a + invphi2 * (b - a)
                fc = f(c)
            else:
                a, c, fc = c, d, fd
                d = a + invphi * (b - a)
                fd = f(d)
        return 0.5 * (a + b)

    def _terminal_miss_distance(self, sigma_cmd_rad, x_nominal_old, rollout_params, t_s):
        """Predicted great-circle miss distance to the target if the capsule
        holds this constant bank to drogue-deploy altitude (Lu Eq. 18-20)."""
        roll = rollout_nominal_to_terminal(
            x_nominal_old=x_nominal_old,
            sigma_cmd_rad=float(sigma_cmd_rad),
            params=rollout_params,
            dt_s=float(self.long_prediction_dt_s),
            terminal_altitude_m=float(self.terminal_prediction_altitude_m),
            max_steps=int(self.long_prediction_max_steps),
        )
        xf = roll["x_final"]
        sph = spherical_range_and_course(
            phi_rad=float(xf[1]), lam_rad=float(xf[2]),
            target_phi_rad=float(self.target_phi_rad),
            target_lam_rad=float(self.target_lam_rad),
            radius_m=float(constants.RADIUS_EARTH),
        )
        miss_m = float(sph["range_to_go_m"])
        # If the rollout skipped out / capped without reaching terminal
        # altitude, penalize so the search avoids that bank region. Skip-out
        # happens at LOW bank (too much vertical lift), so this penalty keeps
        # the cost unimodal and pushes the solution toward higher bank.
        if not roll["reached_terminal"]:
            miss_m += 5.0e6
        return miss_m, roll

    def _solve_sigma_magnitude(self, bank_sign, x_nominal_old, rollout_params, t_s):
        """Golden-section search for the bank magnitude minimizing terminal
        miss distance. Returns (best_sigma_mag_rad, best_miss_m, eval_log)."""
        sigma_max = math.radians(float(self.sigma_search_max_deg))
        eval_log: List[Dict[str, float]] = []

        def cost(mag_rad):
            """Return the predicted terminal miss for a bank magnitude, and log it."""
            if mag_rad < 1.0e-9 or bank_sign == 0.0:
                sigma = 0.0
            else:
                sigma = float(bank_sign) * float(mag_rad)
            sigma = clamp(sigma, float(constants.SIGMA_MIN_RAD), float(constants.SIGMA_MAX_RAD))
            miss_m, roll = self._terminal_miss_distance(sigma, x_nominal_old, rollout_params, t_s)
            eval_log.append({
                "sigma_mag_deg": math.degrees(float(mag_rad)),
                "sigma_cmd_deg": math.degrees(float(sigma)),
                "miss_m": float(miss_m),
                "reached_terminal": bool(roll["reached_terminal"]),
            })
            return miss_m

        best_mag = self._golden_section_minimize(cost, 0.0, sigma_max, int(self.sigma_search_iters))
        best_miss = cost(best_mag)
        return float(best_mag), float(best_miss), eval_log

    def _compute_sigma_cmd_continuous(
        self, t_s, state, obs, bank_sign,
        heading_error_rad, deadband_rad,
        current_altitude_m, current_altitude_rate_mps,
        current_range_to_go_m, current_cross_track_error_m,
    ) -> float:
        """Continuous predictor-corrector path (Fixes 1+2): solve bank
        magnitude by minimizing predicted terminal miss, keep sign from the
        heading-offset logic, run one short heat screen for feasibility."""
        x_nominal_old = self._state_to_vector(state)
        rollout_params = self._rollout_params()
        ctx = self._prediction_context

        best_mag_rad, best_miss_m, eval_log = self._solve_sigma_magnitude(
            bank_sign, x_nominal_old, rollout_params, t_s
        )

        if best_mag_rad < 1.0e-9 or bank_sign == 0.0:
            sigma_cmd_rad = 0.0
        else:
            sigma_cmd_rad = float(bank_sign) * float(best_mag_rad)
        sigma_cmd_rad = clamp(
            sigma_cmd_rad,
            float(constants.SIGMA_MIN_RAD),
            float(constants.SIGMA_MAX_RAD),
        )

        # Short interval/heat screen at the chosen bank, for feasibility flags
        # and the heating debug keys the loggers expect.
        screen = rollout_predictor_corrector_candidate(
            x_nominal_old=x_nominal_old,
            sigma_cmd_rad=float(sigma_cmd_rad),
            params=rollout_params,
            dt_s=float(self.prediction_dt_s),
            horizon_steps=int(self.predictor_horizon_steps),
            x_interval_old=ctx.x_interval_old,
            supervisor_cfg=ctx.supervisor_cfg,
            heat_shield=ctx.heat_shield,
            t0_s=float(t_s),
            heat_rate_limit=self.heat_rate_limit,
            heat_load_limit=self.heat_load_limit,
            interval_screen_active=bool(ctx.interval_active),
        )
        heat_feasible = bool(screen.heat_feasible)

        # Build candidate_dicts from the search evaluations so the candidate
        # plots still have data. Mark the closest-to-chosen as selected.
        candidate_dicts: List[Dict[str, Any]] = []
        for i, e in enumerate(eval_log):
            candidate_dicts.append({
                "candidate_index": int(i),
                "sigma_mag_deg": float(e["sigma_mag_deg"]),
                "sigma_cmd_deg": float(e["sigma_cmd_deg"]),
                "bank_sign": float(bank_sign),
                "terminal_range_to_go_m": float(e["miss_m"]),
                "geometry_cost": float(e["miss_m"]),
                "heat_penalty": 0.0,
                "total_cost": float(e["miss_m"]),
                "heat_feasible": bool(heat_feasible),
                "fully_feasible": bool(heat_feasible),
                "reached_terminal": bool(e["reached_terminal"]),
                "failure_reason": "",
            })
        sel_idx = min(range(len(candidate_dicts)),
                      key=lambda i: candidate_dicts[i]["total_cost"]) if candidate_dicts else -1

        self._last_debug = {
            "guidance_cycle_index": int(ctx.guidance_cycle_index),
            "trajectory_id": str(ctx.trajectory_id),
            "current_altitude_m": float(current_altitude_m),
            "current_altitude_rate_mps": float(current_altitude_rate_mps),
            "current_range_to_go_m": float(current_range_to_go_m),
            "current_heading_error_rad": float(heading_error_rad),
            "current_cross_track_error_m": float(current_cross_track_error_m),
            "heading_deadband_deg": math.degrees(deadband_rad),
            "predictor_mode": "continuous",
            "candidate_dicts": list(candidate_dicts),
            "candidate_sigma_deg": [float(c["sigma_cmd_deg"]) for c in candidate_dicts],
            "candidate_total_cost": [float(c["total_cost"]) for c in candidate_dicts],
            "candidate_heat_feasible": [bool(c["fully_feasible"]) for c in candidate_dicts],
            "candidate_failure_reason": [str(c["failure_reason"]) for c in candidate_dicts],
            "selected_candidate_index": int(sel_idx),
            "selected_candidate_heat_feasible": bool(heat_feasible),
            "selected_interval_screen_used": bool(screen.interval_screen_used),
            "selected_interval_failure_kind": str(screen.interval_failure_kind),
            "any_feasible_candidate": bool(heat_feasible),
            "chosen_sigma_cmd_deg": math.degrees(float(sigma_cmd_rad)),
            "chosen_sigma_mag_deg": math.degrees(float(best_mag_rad)),
            "selected_geometry_cost": float(best_miss_m),
            "selected_heat_penalty": float(self.weight_heat * float(screen.heat_penalty)),
            "selected_total_cost": float(best_miss_m),
            "selected_failure_reason": "" if heat_feasible else "heat_infeasible_at_solution",
            "selected_violation_amount": float(screen.violation_amount),
            "selected_first_violation_step": int(screen.first_violation_step),
            "selected_first_violation_time_s": float(screen.first_violation_time_s),
            "selected_max_heating_rate_lo": float(screen.max_heating_rate_interval.lo),
            "selected_max_heating_rate_hi": float(screen.max_heating_rate_interval.hi),
            "selected_max_heat_load_lo": float(screen.max_heat_load_interval.lo),
            "selected_max_heat_load_hi": float(screen.max_heat_load_interval.hi),
            "predicted_terminal_miss_m": float(best_miss_m),
        }
        return float(sigma_cmd_rad)

    def compute_sigma_cmd(
        self,
        t_s: float,
        state: ReentryState,
        obs: Dict[str, float],
    ) -> float:
        """
        Return the ideal commanded bank angle in radians.

        This keeps the same public interface as the older guidance law,
        but the internal implementation is now a short horizon heat aware
        predictor corrector.

        The selected sigma_cmd is still only a command.
        The downstream bank actuator and roll controller remain unchanged.
        """
        # Extract the current guidance features from the observation dictionary.
        heading_error_rad = float(obs.get("heading_error_rad", 0.0))
        heading_offset_psi_rad = float(obs.get("heading_offset_psi_rad", -heading_error_rad))
        current_range_to_go_m = float(obs.get("range_to_go_m", 0.0))
        current_cross_track_error_m = float(obs.get("cross_track_error_m", 0.0))
        current_altitude_m = float(obs.get("altitude_m", max(0.0, state.r_m - constants.RADIUS_EARTH)))

        # Compute the current altitude rate for debug output and for the
        # conceptual link to Lu style altitude rate aware logic.
        current_altitude_rate_mps = float(state.V_mps) * math.sin(float(state.gamma_rad))

        # Compute the current heading deadband from speed.
        deadband_rad = self.heading_deadband_rad(float(state.V_mps))

        # Build a default debug dictionary first so every return path provides
        # a stable set of keys for the notebook logger.
        self._last_debug = {
            "guidance_cycle_index": int(self._prediction_context.guidance_cycle_index),
            "trajectory_id": str(self._prediction_context.trajectory_id),
            "current_altitude_m": float(current_altitude_m),
            "current_altitude_rate_mps": float(current_altitude_rate_mps),
            "current_range_to_go_m": float(current_range_to_go_m),
            "current_heading_error_rad": float(heading_error_rad),
            "current_cross_track_error_m": float(current_cross_track_error_m),
            "heading_deadband_deg": math.degrees(deadband_rad),
            "candidate_dicts": [],
            "candidate_sigma_deg": [],
            "candidate_total_cost": [],
            "candidate_heat_feasible": [],
            "candidate_failure_reason": [],
            "selected_candidate_index": -1,
            "selected_candidate_heat_feasible": True,
            "any_feasible_candidate": True,
            "chosen_sigma_cmd_deg": 0.0,
            "chosen_sigma_mag_deg": 0.0,
            "selected_geometry_cost": 0.0,
            "selected_heat_penalty": 0.0,
            "selected_total_cost": 0.0,
            "selected_failure_reason": "",
            "selected_violation_amount": 0.0,
            "selected_first_violation_step": -1,
            "selected_first_violation_time_s": math.nan,
            "selected_max_heating_rate_lo": 0.0,
            "selected_max_heating_rate_hi": 0.0,
            "selected_max_heat_load_lo": 0.0,
            "selected_max_heat_load_hi": 0.0,
        }

        # Optionally disable guidance above a speed threshold.
        # This keeps compatibility with the old class behavior.
        if float(state.V_mps) > float(self.velocity_enable_mps):
            self._last_debug["selected_failure_reason"] = "guidance_disabled_above_velocity_enable"
            return 0.0

        # Gamma safety override: if the capsule has pitched into a steep dive,
        # abandon cross-range tracking and roll wings-level (sigma = 0) so all
        # available lift goes vertical. Without this, the predictor-corrector
        # keeps commanding ~70 deg bank toward an unreachable target while
        # gamma collapses past -78 deg and trips the cos_gamma safety gate
        # well above the drogue deploy altitude.
        gamma_safety_threshold_rad = math.radians(30.0)
        if abs(float(state.gamma_rad)) > gamma_safety_threshold_rad:
            self._last_debug["selected_failure_reason"] = "gamma_safety_override_wings_level"
            self._last_debug["chosen_sigma_cmd_deg"] = 0.0
            self._last_debug["chosen_sigma_mag_deg"] = 0.0
            return 0.0

        # Determine commanded bank sign. With bank reversal enabled the sign
        # is held through a velocity-dependent deadband and only reversed when
        # the heading error decisively crosses it (Apollo / Lu 2014); this is
        # what makes the classical baseline competent. Otherwise fall back to
        # the naive instantaneous-sign rule.
        if self.use_bank_reversal:
            bank_sign = self._bank_sign_with_reversal(
                heading_offset_psi_rad=heading_offset_psi_rad,
                heading_error_rad=heading_error_rad,
                V_mps=float(state.V_mps),
            )
        else:
            bank_sign = self._bank_sign_from_heading_offset(
                heading_error_rad=heading_error_rad,
                heading_offset_psi_rad=heading_offset_psi_rad,
            )

        # Continuous predictor-corrector path (Lu 2014, Fixes 1 + 2): solve
        # the bank magnitude that minimizes predicted terminal miss distance.
        if self.use_continuous_predictor:
            return self._compute_sigma_cmd_continuous(
                t_s=t_s, state=state, obs=obs, bank_sign=bank_sign,
                heading_error_rad=heading_error_rad, deadband_rad=deadband_rad,
                current_altitude_m=current_altitude_m,
                current_altitude_rate_mps=current_altitude_rate_mps,
                current_range_to_go_m=current_range_to_go_m,
                current_cross_track_error_m=current_cross_track_error_m,
            )

        # --- Legacy discrete-bin selector (use_continuous_predictor=False) ---
        # Build the candidate magnitude list in degrees.
        candidate_mags_deg = self._candidate_magnitude_list_deg(
            heading_error_rad=heading_error_rad,
            deadband_rad=deadband_rad,
        )

        # Convert the current state object into the translational vector used
        # by the rollout utilities in math_3d.py.
        x_nominal_old = self._state_to_vector(state)

        # Pull rollout parameters once for all candidates.
        rollout_params = self._rollout_params()

        # Copy the latest context object locally for readability.
        ctx = self._prediction_context

        # Store full candidate telemetry here.
        candidate_dicts: List[Dict[str, Any]] = []

        # Evaluate each candidate bank magnitude in turn.
        for candidate_index, sigma_mag_deg in enumerate(candidate_mags_deg):
            # Convert candidate magnitude into radians.
            sigma_mag_rad = math.radians(float(sigma_mag_deg))

            # Apply the chosen bank sign.
            # The zero bank candidate remains exactly zero.
            if abs(sigma_mag_rad) <= 1.0e-12 or bank_sign == 0.0:
                sigma_cmd_candidate_rad = 0.0
            else:
                sigma_cmd_candidate_rad = float(bank_sign) * float(sigma_mag_rad)

            # Clamp the candidate into the global sigma limits.
            sigma_cmd_candidate_rad = clamp(
                sigma_cmd_candidate_rad,
                float(constants.SIGMA_MIN_RAD),
                float(constants.SIGMA_MAX_RAD),
            )

            # Roll the candidate forward through the corrected translational
            # model plus the interval heat shield model.
            rollout = rollout_predictor_corrector_candidate(
                x_nominal_old=x_nominal_old,
                sigma_cmd_rad=float(sigma_cmd_candidate_rad),
                params=rollout_params,
                dt_s=float(self.prediction_dt_s),
                horizon_steps=int(self.predictor_horizon_steps),
                x_interval_old=ctx.x_interval_old,
                supervisor_cfg=ctx.supervisor_cfg,
                heat_shield=ctx.heat_shield,
                t0_s=float(t_s),
                heat_rate_limit=self.heat_rate_limit,
                heat_load_limit=self.heat_load_limit,
                interval_screen_active=bool(ctx.interval_active),
            )

            # Compute predicted terminal targeting quantities from the final state.
            terminal_errors = self._predicted_terminal_errors(
                final_state=rollout.final_nominal_state
            )

            # Compute the geometry portion of the cost.
            geometry_cost = self._score_rollout_geometry(terminal_errors)

            # Start the heat penalty from the rollout internal violation amount.
            heat_penalty = self.weight_heat * float(rollout.heat_penalty)

            # Add a very large penalty if the candidate violates heat feasibility.
            if not rollout.heat_feasible:
                heat_penalty += float(self.infeasible_penalty)

            # Add an even larger penalty only when interval screening was active
            # and interval propagation became numerically invalid.
            if bool(rollout.interval_screen_used) and not rollout.interval_valid:
                heat_penalty += float(self.invalid_interval_penalty)

            # Total cost is geometry plus heat related penalties.
            total_cost = float(geometry_cost + heat_penalty)

            # A candidate is fully feasible when heat remained feasible and
            # any active interval screen also remained numerically valid.
            candidate_fully_feasible = bool(
                rollout.heat_feasible
                and (not bool(rollout.interval_screen_used) or bool(rollout.interval_valid))
            )

            # Convert the rollout result into a debug rich dictionary.
            candidate_debug = {
                "candidate_index": int(candidate_index),
                "sigma_mag_deg": float(sigma_mag_deg),
                "sigma_cmd_rad": float(sigma_cmd_candidate_rad),
                "sigma_cmd_deg": math.degrees(float(sigma_cmd_candidate_rad)),
                "bank_sign": float(bank_sign),
                "terminal_range_to_go_m": float(terminal_errors["range_to_go_m"]),
                "terminal_heading_error_rad": float(terminal_errors["heading_error_rad"]),
                "terminal_cross_track_error_m": float(terminal_errors["cross_track_error_m"]),
                "geometry_cost": float(geometry_cost),
                "heat_penalty": float(heat_penalty),
                "total_cost": float(total_cost),
                "heat_feasible": bool(rollout.heat_feasible),
                "interval_valid": bool(rollout.interval_valid),
                "interval_screen_used": bool(rollout.interval_screen_used),
                "interval_failure_kind": str(rollout.interval_failure_kind),
                "fully_feasible": bool(candidate_fully_feasible),
                "violation_amount": float(rollout.violation_amount),
                "first_violation_step": int(rollout.first_violation_step),
                "first_violation_time_s": float(rollout.first_violation_time_s),
                "failure_reason": str(rollout.failure_reason),
                "max_heating_rate_lo": float(rollout.max_heating_rate_interval.lo),
                "max_heating_rate_hi": float(rollout.max_heating_rate_interval.hi),
                "max_heat_load_lo": float(rollout.max_heat_load_interval.lo),
                "max_heat_load_hi": float(rollout.max_heat_load_interval.hi),
            }

            # Append the candidate telemetry to the full list.
            candidate_dicts.append(candidate_debug)

        # Determine whether any candidate is fully feasible.
        feasible_indices = [
            idx for idx, cand in enumerate(candidate_dicts)
            if bool(cand["fully_feasible"])
        ]

        # Prefer the lowest total cost among fully feasible candidates.
        if feasible_indices:
            best_index = min(
                feasible_indices,
                key=lambda idx: float(candidate_dicts[idx]["total_cost"]),
            )
            any_feasible_candidate = True
        else:
            # If every candidate is infeasible, still choose the least bad one.
            best_index = min(
                range(len(candidate_dicts)),
                key=lambda idx: float(candidate_dicts[idx]["total_cost"]),
            )
            any_feasible_candidate = False

        # Pull the chosen candidate dictionary.
        chosen = candidate_dicts[best_index]

        # Record rich debug output for logging and RL ready export.
        self._last_debug = {
            "guidance_cycle_index": int(ctx.guidance_cycle_index),
            "trajectory_id": str(ctx.trajectory_id),
            "current_altitude_m": float(current_altitude_m),
            "current_altitude_rate_mps": float(current_altitude_rate_mps),
            "current_range_to_go_m": float(current_range_to_go_m),
            "current_heading_error_rad": float(heading_error_rad),
            "current_cross_track_error_m": float(current_cross_track_error_m),
            "heading_deadband_deg": math.degrees(deadband_rad),
            "candidate_dicts": list(candidate_dicts),
            "candidate_sigma_deg": [float(c["sigma_cmd_deg"]) for c in candidate_dicts],
            "candidate_total_cost": [float(c["total_cost"]) for c in candidate_dicts],
            "candidate_heat_feasible": [bool(c["fully_feasible"]) for c in candidate_dicts],
            "candidate_failure_reason": [str(c["failure_reason"]) for c in candidate_dicts],
            "selected_candidate_index": int(best_index),
            "selected_candidate_heat_feasible": bool(chosen["fully_feasible"]),
            "selected_interval_screen_used": bool(chosen.get("interval_screen_used", True)),
            "selected_interval_failure_kind": str(chosen.get("interval_failure_kind", "")),
            "any_feasible_candidate": bool(any_feasible_candidate),
            "chosen_sigma_cmd_deg": float(chosen["sigma_cmd_deg"]),
            "chosen_sigma_mag_deg": float(chosen["sigma_mag_deg"]),
            "selected_geometry_cost": float(chosen["geometry_cost"]),
            "selected_heat_penalty": float(chosen["heat_penalty"]),
            "selected_total_cost": float(chosen["total_cost"]),
            "selected_failure_reason": str(chosen["failure_reason"]),
            "selected_violation_amount": float(chosen["violation_amount"]),
            "selected_first_violation_step": int(chosen["first_violation_step"]),
            "selected_first_violation_time_s": float(chosen["first_violation_time_s"]),
            "selected_max_heating_rate_lo": float(chosen["max_heating_rate_lo"]),
            "selected_max_heating_rate_hi": float(chosen["max_heating_rate_hi"]),
            "selected_max_heat_load_lo": float(chosen["max_heat_load_lo"]),
            "selected_max_heat_load_hi": float(chosen["max_heat_load_hi"]),
        }

        # Return the chosen sigma command in radians.
        return float(chosen["sigma_cmd_rad"])


# Guidance scheduler ------------------------------------------------------------------------------
class GuidanceScheduler:
    """
    Gate guidance updates so sigma_cmd refreshes more slowly than the physics loop.
    """

    def __init__(self, period_s: float):
        """Build the scheduler with the guidance update period in seconds."""
        # Store the allowed guidance update period.
        self.period_s = float(period_s)

        # This tracks the last time a guidance update was granted.
        self._last_update_time: Optional[float] = None

    def reset(self) -> None:
        """
        Reset the internal update timer.
        """
        self._last_update_time = None

    def should_update(self, t_s: float) -> bool:
        """
        Return True only when guidance is allowed to update.
        """
        # The first call should always allow an update.
        if self._last_update_time is None:
            self._last_update_time = t_s
            return True

        # Compute elapsed time since the last accepted update.
        elapsed = t_s - self._last_update_time

        # Grant an update once the configured period has passed.
        if elapsed >= self.period_s:
            self._last_update_time = t_s
            return True

        # Otherwise guidance keeps the previous sigma_cmd.
        return False


# Aero authority utility ------------------------------------------------------------------------------
def compute_aero_accel_mps2(
    lift_N: float,
    drag_N: float,
    mass_kg: float,
) -> float:
    """
    Return aerodynamic acceleration magnitude.
    """
    # Mass must stay positive for the force to acceleration conversion.
    if mass_kg <= 0.0:
        raise ValueError("mass_kg must be positive")

    # Combine lift and drag into a total aerodynamic force magnitude.
    total_force = math.sqrt(lift_N * lift_N + drag_N * drag_N)

    # Convert the net force into an acceleration magnitude.
    return total_force / mass_kg


# Observation layer ------------------------------------------------------------------------------
class BasicObservationProvider:
    """
    Build guidance features from the raw reentry state.

    This uses a spherical Earth targeting model and preserves the current
    chi convention rather than silently changing it.

    Returned features include both the native chi based heading error and
    the north referenced heading offset used to express Lu style sign logic.
    """

    def __init__(
        self,
        target_phi_rad: float,
        target_lam_rad: float,
        params: Optional[dict] = None,
    ):
        """
        Build the observation provider.

        Parameters:
            target_phi_rad: target latitude in radians.
            target_lam_rad: target longitude in radians.
            params:         optional parameter dict shared with guidance.
        """
        # Store the target location.
        self.target_phi = float(target_phi_rad)
        self.target_lam = float(target_lam_rad)

        # Keep a params dictionary so it can be shared with guidance if needed.
        self.params = params or {}

    def reset(self) -> None:
        """
        Reset hook for interface symmetry.
        """
        pass

    def observe(self, state: ReentryState, t_s: float) -> Dict[str, float]:
        """
        Return a feature dictionary for guidance.
        """
        # Time is not directly used in the current observation logic.
        del t_s

        # Convert radius into geometric altitude.
        altitude_m = max(0.0, state.r_m - constants.RADIUS_EARTH)

        # Query the nominal atmosphere at the current altitude.
        atm = AtmosphereModel.US_Standard_ATM(
            altitude_m,
            lat_deg=math.degrees(state.phi_rad),
            lon_deg=math.degrees(state.lam_rad),
            date=self.params.get("entry_datetime_utc", None),
            f107=self.params.get("f107", 150.0),
            f107a=self.params.get("f107a", 150.0),
            ap=self.params.get("ap", 7.0),
        )

        # Use zero density outside the supported atmosphere range.
        rho = atm["rho_kgm3"] if atm["rho_kgm3"] is not None else 0.0

        # Dynamic pressure is one half rho V squared.
        dynamic_pressure_pa = 0.5 * rho * state.V_mps * state.V_mps

        # Compute spherical range and course geometry to the target.
        spherical = spherical_range_and_course(
            phi_rad=float(state.phi_rad),
            lam_rad=float(state.lam_rad),
            target_phi_rad=float(self.target_phi),
            target_lam_rad=float(self.target_lam),
            radius_m=float(constants.RADIUS_EARTH),
        )

        # Handle the exact target coincidence case cleanly.
        if spherical["range_to_go_m"] <= 1.0e-9:
            target_course_chi_rad = float(state.chi_rad)
            target_azimuth_north_rad = chi_to_psi_north_rad(float(state.chi_rad))
            heading_error_rad = 0.0
            heading_offset_psi_rad = 0.0
            current_heading_psi_rad = chi_to_psi_north_rad(float(state.chi_rad))
        else:
            # Pull the target course in both chi and north azimuth form.
            target_course_chi_rad = spherical["target_course_chi_rad"]
            target_azimuth_north_rad = spherical["target_azimuth_north_rad"]

            # Compute current heading in the north referenced azimuth convention.
            current_heading_psi_rad = chi_to_psi_north_rad(float(state.chi_rad))

            # Native chi based heading error is target course minus current heading.
            heading_error_rad = wrap_to_pi(target_course_chi_rad - float(state.chi_rad))

            # Lu style heading offset in psi form is current heading minus target azimuth.
            heading_offset_psi_rad = wrap_to_pi(current_heading_psi_rad - target_azimuth_north_rad)

        # Resolve the current range to go into along track and cross track components.
        along_track_error_m = spherical["range_to_go_m"] * math.cos(heading_error_rad)
        cross_track_error_m = spherical["range_to_go_m"] * math.sin(heading_error_rad)

        # Return a rich observation dictionary for guidance and logging.
        return {
            "altitude_m": float(altitude_m),
            "density_kgm3": float(rho),
            "dynamic_pressure_pa": float(dynamic_pressure_pa),
            "speed_mps": float(state.V_mps),
            "flight_path_angle_rad": float(state.gamma_rad),
            "sigma_actual_rad": float(state.sigma_actual_rad),
            "roll_rate_rad_s": float(state.roll_rate_rad_s),
            "range_to_go_m": float(spherical["range_to_go_m"]),
            "great_circle_range_m": float(spherical["range_to_go_m"]),
            "great_circle_angle_rad": float(spherical["central_angle_rad"]),
            "target_course_chi_rad": float(target_course_chi_rad),
            "target_azimuth_north_rad": float(target_azimuth_north_rad),
            "current_heading_psi_rad": float(current_heading_psi_rad),
            "heading_error_rad": float(heading_error_rad),
            "heading_offset_psi_rad": float(heading_offset_psi_rad),
            "north_error_m": float(spherical["north_error_m"]),
            "east_error_m": float(spherical["east_error_m"]),
            "along_track_error_m": float(along_track_error_m),
            "cross_track_error_m": float(cross_track_error_m),
            "crossrange_error_m": float(cross_track_error_m),
            "downrange_error_m": float(along_track_error_m),
            "target_phi_rad": float(self.target_phi),
            "target_lam_rad": float(self.target_lam),
        }


# Bank command shaper ------------------------------------------------------------------------------

@dataclass
class BankActuatorLimits:
    """
    Limits for smoothing sigma_cmd into sigma_target.
    """
    sigma_rate_max_rps: float
    sigma_accel_max_rps2: float


class BankActuator:
    """
    Smooth sigma_cmd into sigma_target using simple second order limits.
    """

    def __init__(self, limits: BankActuatorLimits):
        """Build the actuator with its rate and acceleration limits."""
        # Store the command shaping limits.
        self.limits = limits

        # Initialize the target bank command state.
        self.sigma_target_rad = 0.0
        self.sigma_target_dot_rps = 0.0

    def reset(
        self,
        sigma_target0_rad: float = 0.0,
        sigma_target_dot0_rps: float = 0.0,
    ) -> None:
        """
        Reset the internal target bank state.
        """
        self.sigma_target_rad = float(wrap_to_pi(sigma_target0_rad))
        self.sigma_target_dot_rps = float(sigma_target_dot0_rps)

    def step(
        self,
        dt_s: float,
        sigma_cmd_rad: float,
    ) -> tuple[float, float]:
        """
        Advance the command shaper by one step.
        """
        # A nonpositive time step is invalid for forward propagation.
        if dt_s <= 0.0:
            raise ValueError("dt_s must be positive")

        # Wrap the incoming bank command into the principal angle range.
        sigma_cmd_rad = wrap_to_pi(sigma_cmd_rad)

        # Compute the current target tracking error.
        prev_error = wrap_to_pi(sigma_cmd_rad - self.sigma_target_rad)

        # If the target already matches the command, stop target motion.
        if abs(prev_error) <= 1.0e-12:
            self.sigma_target_rad = sigma_cmd_rad
            self.sigma_target_dot_rps = 0.0
            return self.sigma_target_rad, self.sigma_target_dot_rps

        # Accelerate the target bank state toward the commanded bank sign.
        sigma_ddot_cmd = sign_no_zero(prev_error) * self.limits.sigma_accel_max_rps2

        # Integrate target bank rate.
        self.sigma_target_dot_rps += sigma_ddot_cmd * dt_s

        # Enforce the maximum target bank rate magnitude.
        self.sigma_target_dot_rps = clamp(
            self.sigma_target_dot_rps,
            -self.limits.sigma_rate_max_rps,
            self.limits.sigma_rate_max_rps,
        )

        # Integrate the target bank angle.
        self.sigma_target_rad = wrap_to_pi(
            self.sigma_target_rad + self.sigma_target_dot_rps * dt_s
        )

        # Compute the new residual error after the integration step.
        new_error = wrap_to_pi(sigma_cmd_rad - self.sigma_target_rad)

        # Detect whether the shaped target crossed the command.
        crossed = (
            (prev_error > 0.0 and new_error < 0.0)
            or (prev_error < 0.0 and new_error > 0.0)
            or abs(new_error) <= 1.0e-12
        )

        # Snap exactly to the command after crossing to avoid chatter.
        if crossed:
            self.sigma_target_rad = sigma_cmd_rad
            self.sigma_target_dot_rps = 0.0

        return self.sigma_target_rad, self.sigma_target_dot_rps


# Roll torque controller ------------------------------------------------------------------------------

@dataclass
class RollTorqueController:
    """
    PD controller that converts bank tracking error into a roll torque command.
    """

    kp_Nm_per_rad: float = constants.ROLL_KP_NM_PER_RAD
    kd_Nm_per_rad_s: float = constants.ROLL_KD_NM_PER_RAD_S
    max_torque_Nm: float = constants.ROLL_CMD_MAX_TORQUE_NM

    sigma_deadband_rad: float = constants.ROLL_SIGMA_DEADBAND_RAD
    rate_deadband_rad_s: float = constants.ROLL_RATE_DEADBAND_RAD_S

    def reset(self) -> None:
        """
        Reset hook for interface symmetry.
        """
        pass

    def step(
        self,
        sigma_target_rad: float,
        sigma_actual_rad: float,
        roll_rate_rad_s: float,
    ) -> float:
        """
        Return the commanded body roll torque in newton meters.
        """
        # Compute wrapped sigma tracking error.
        sigma_error = wrap_to_pi(sigma_target_rad - sigma_actual_rad)

        # If both sigma error and roll rate are tiny, command zero torque.
        if (
            abs(sigma_error) < self.sigma_deadband_rad
            and abs(roll_rate_rad_s) < self.rate_deadband_rad_s
        ):
            return 0.0

        # Compute PD roll torque command.
        tau_roll_cmd = (
            self.kp_Nm_per_rad * sigma_error
            - self.kd_Nm_per_rad_s * roll_rate_rad_s
        )

        # Enforce the maximum torque magnitude.
        tau_roll_cmd = clamp(
            tau_roll_cmd,
            -self.max_torque_Nm,
            self.max_torque_Nm,
        )

        return float(tau_roll_cmd)


# Control stack config and debug output ------------------------------------------------------------------------------
@dataclass
class CapsuleControlConfig:
    """
    Configuration shared across the full capsule control stack.
    """
    guidance_period_s: float = constants.GUIDANCE_PERIOD_S


@dataclass
class CapsuleControlOutput:
    """
    Per step control telemetry.
    """
    # Full observation dictionary used by guidance.
    obs: Dict[str, float] = field(default_factory=dict)

    # Whether guidance updated sigma_cmd on this step.
    guidance_updated: bool = False

    # Most recent ideal bank command.
    sigma_cmd_rad: float = 0.0

    # Smoothed target bank command.
    sigma_target_rad: float = 0.0

    # Smoothed target bank rate.
    sigma_target_dot_rps: float = 0.0

    # Actual flown bank angle.
    sigma_actual_rad: float = 0.0

    # Actual roll rate.
    roll_rate_rad_s: float = 0.0

    # Commanded roll torque for the RCS allocator.
    tau_roll_cmd_Nm: float = 0.0

    # Instantaneous aerodynamic acceleration magnitude.
    aero_accel_mps2: float = 0.0

    # Guidance debug packet for logging and RL ready export.
    guidance_debug: Dict[str, Any] = field(default_factory=dict)


# Full stack ------------------------------------------------------------------------------
class CapsuleControlStack:
    """
    Full high level control chain for milestone 1.

    This class does not update sigma_actual itself.
    Actual bank angle is a physical result of the RCS and roll dynamics layer.
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
        """
        Build the full control stack from its component objects.

        Parameters:
            cfg:             control stack configuration.
            scheduler:       gates how often guidance updates.
            guidance:        the guidance law that proposes the bank.
            obs_provider:    builds guidance features from the state.
            bank_actuator:   smooths the commanded bank into a target.
            roll_controller: turns bank error into a roll torque command.
        """
        # Store all stack components.
        self.cfg = cfg
        self.scheduler = scheduler
        self.guidance = guidance
        self.obs_provider = obs_provider
        self.bank_actuator = bank_actuator
        self.roll_controller = roll_controller

        # Hold the most recent sigma command between guidance update times.
        self._sigma_cmd_hold_rad = 0.0

        # Initialize the per step debug object.
        self.debug = CapsuleControlOutput()

        # If guidance has no explicit params but the observation provider does,
        # copy them into guidance so predictor rollouts use the same values.
        if hasattr(self.guidance, "params") and getattr(self.guidance, "params", None) is None:
            if hasattr(self.obs_provider, "params"):
                try:
                    self.guidance.params = dict(self.obs_provider.params)
                except Exception:
                    pass

    def reset(self) -> None:
        """
        Reset all control stack submodules.
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
        Advance the high level control stack by one step.
        """
        # A nonpositive time step is invalid for forward propagation.
        if dt_s <= 0.0:
            raise ValueError("dt_s must be positive")

        # Build the observation dictionary from the current state
        obs = self.obs_provider.observe(state, t_s)

        # Check whether guidance is allowed to update on this step
        guidance_updated = self.scheduler.should_update(t_s)

        # Only recompute sigma_cmd when the scheduler grants an update
        if guidance_updated:
            self._sigma_cmd_hold_rad = float(
                self.guidance.compute_sigma_cmd(t_s, state, obs)
            )

        # Smooth the held sigma command into a target sigma.
        sigma_target_rad, sigma_target_dot_rps = self.bank_actuator.step(
            dt_s=dt_s,
            sigma_cmd_rad=self._sigma_cmd_hold_rad,
        )

        # Pull the currently flown bank and roll rate from the physical state
        sigma_actual_rad = float(state.sigma_actual_rad)
        roll_rate_rad_s = float(state.roll_rate_rad_s)

        # Convert the bank tracking error into a roll torque command
        tau_roll_cmd_Nm = self.roll_controller.step(
            sigma_target_rad=sigma_target_rad,
            sigma_actual_rad=sigma_actual_rad,
            roll_rate_rad_s=roll_rate_rad_s,
        )

        # Compute the aerodynamic acceleration magnitude for logging
        aero_accel_mps2 = compute_aero_accel_mps2(
            lift_N=lift_N,
            drag_N=drag_N,
            mass_kg=mass_kg,
        )

        # Write the latest commanded values back into the state object.
        state.sigma_cmd_rad = float(self._sigma_cmd_hold_rad)
        state.sigma_target_rad = float(sigma_target_rad)

        # Pull the latest guidance debug dictionary if the guidance object
        # exposes one. This makes the notebook logger richer without changing
        # the main control stack interface.
        guidance_debug: Dict[str, Any] = {}
        if hasattr(self.guidance, "get_last_debug"):
            try:
                guidance_debug = dict(self.guidance.get_last_debug())
            except Exception:
                guidance_debug = {}

        # Build the per step control output packet.
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
            guidance_debug=dict(guidance_debug),
        )

        return self.debug


# Backward compatibility aliases
SigmaControlConfig = CapsuleControlConfig
SigmaControlDebug = CapsuleControlOutput
SigmaControlStack = CapsuleControlStack