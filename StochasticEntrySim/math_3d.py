"""
math_3d.py

Interval and nominal translational reentry dynamics for milestone 1.

This file keeps the existing interval supervisor architecture while using a
corrected spherical Earth point mass entry model with bank modulated lift.

State convention kept intact

X = [r, phi, lam, V, gamma, chi]

r      radial distance from Earth center in meters
phi    geocentric latitude in radians
lam    longitude in radians
V      speed magnitude in m per s
gamma  flight path angle in radians
chi    heading angle in radians measured from east toward north

Modeling scope

spherical Earth
non rotating Earth for milestone 1
point mass 3DOF translational dynamics
bank angle modulates lift between vertical plane and lateral components
drag acts only opposite the velocity vector and therefore only enters V_dot

The interval supervisor remains an annotation layer around the nominal path.

Important milestone 1 approximation note

The predictor corrector support added here uses the current corrected 3DOF
translational model plus the current interval heat shield model. That is the
closest clean milestone 1 fit to the current codebase. It is not a literal
line by line reproduction of every equation in Lu.
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import AtmosphereModel
import constants
from interval_math import Interval, box_add, box_scalar_mul, promote


# Keep one canonical ordering for state names.
# This is used by interval width dictionaries and debug output.
STATE_NAMES = ("r", "phi", "lam", "V", "gamma", "chi")

# Configuration and result models ---------------------------------

@dataclass
class IntervalSupervisorConfig:
    """
    Half widths used to inflate a nominal translational state into an interval box.
    """

    # Half width for radius in meters.
    r_half_width_m: float = 0.0

    # Half width for latitude in radians.
    phi_half_width_rad: float = 0.0

    # Half width for longitude in radians.
    lam_half_width_rad: float = 0.0

    # Half width for speed in meters per second.
    V_half_width_mps: float = 0.0

    # Half width for flight path angle in radians.
    gamma_half_width_rad: float = 0.0

    # Half width for heading angle in radians.
    chi_half_width_rad: float = 0.0

    # Optional lower altitude constraint for interval summary checks.
    min_altitude_m: Optional[float] = None

    # Optional upper altitude constraint for interval summary checks.
    max_altitude_m: Optional[float] = None

    # Optional upper speed constraint for interval summary checks.
    max_speed_mps: Optional[float] = None

    # Optional upper dynamic pressure constraint for interval summary checks.
    max_dynamic_pressure_pa: Optional[float] = None

    # When True, interval heating is propagated during annotation.
    include_heating: bool = False


@dataclass
class IntervalAnnotationResult:
    """
    Interval supervisor output for one annotated translational step
    """

    # Interval state before the propagated step
    x_interval_old: List[Interval]

    # Interval state after the propagated step
    x_interval_new: List[Interval]

    # Interval derivative used during propagation
    dx_interval: List[Interval]

    # Bank interval used for this step
    sigma_interval_used: Interval

    # Interval density hull used for this step
    rho_interval: Interval

    # Interval dynamic pressure for this step
    q_interval: Interval

    # Interval geometric altitude for this step
    altitude_interval: Interval

    # Interval temperature hull if available
    temperature_interval: Optional[Interval]

    # Interval pressure hull if available
    pressure_interval: Optional[Interval]

    # Widths of the interval state before stepping
    state_widths_old: Dict[str, float] = field(default_factory=dict)

    # Widths of the interval state after stepping.
    state_widths_new: Dict[str, float] = field(default_factory=dict)

    # Widths of the interval derivative.
    dx_widths: Dict[str, float] = field(default_factory=dict)

    # Maximum cell heating rate interval if heating is enabled
    heating_qdot_max_interval: Optional[Interval] = None

    # Mean heating rate interval if heating is enabled.
    heating_qdot_mean_interval: Optional[Interval] = None

    # Maximum accumulated heat interval if heating is enabled
    heating_Q_max_interval: Optional[Interval] = None

    # Heat shield object after the propagated step
    heat_shield: Optional[Any] = None

    # Summary status from simple bound checks
    safety_status: str = "not_evaluated"

    # Per limit check status values.
    safety_checks: Dict[str, str] = field(default_factory=dict)

    # Atmospheric layer indices intersected by the altitude interval
    layer_indices: List[int] = field(default_factory=list)


@dataclass
class RolloutStepSummary:
    """
    One step summary for a predictor corrector candidate rollout.
    """

    # Step index inside the candidate rollout.
    step_index: int

    # Physical time associated with this step.
    time_s: float

    # Nominal state before the step.
    x_nominal_before: List[float]

    # Nominal state after the step.
    x_nominal_after: List[float]

    # Interval state after the step if propagation was valid.
    x_interval_after: Optional[List[Interval]]

    # Maximum heating rate interval seen at this step.
    max_heating_rate_interval: Interval

    # Maximum accumulated heat interval seen at this step.
    max_heat_load_interval: Interval

    # Whether the candidate is still heat feasible after this step.
    heat_feasible_after_step: bool

    # Whether interval propagation remained valid after this step.
    interval_valid_after_step: bool

    # Failure reason string if the candidate has failed.
    failure_reason: str = ""


@dataclass
class PredictorCorrectorRolloutResult:
    """
    Summary returned to guidance for one candidate rollout.

    Guidance uses this object to combine geometry prediction and heat
    feasibility inside one coherent candidate evaluation.
    """

    # Candidate bank command being evaluated.
    sigma_cmd_rad: float

    # Number of short horizon steps requested.
    horizon_steps: int

    # Sequence of nominal states including the initial state.
    nominal_states: List[List[float]] = field(default_factory=list)

    # Sequence of interval states including the initial state when available.
    interval_states: List[List[Interval]] = field(default_factory=list)

    # Step by step summary packets.
    step_summaries: List[RolloutStepSummary] = field(default_factory=list)

    # Final nominal state.
    final_nominal_state: List[float] = field(default_factory=list)

    # Final interval state if propagation remained valid.
    final_interval_state: Optional[List[Interval]] = None

    # Final heat shield state after the rollout.
    final_heat_shield: Optional[Any] = None

    # Hull of maximum heating rate across the rollout horizon.
    max_heating_rate_interval: Interval = field(default_factory=lambda: Interval(0.0, 0.0))

    # Hull of maximum accumulated heat across the rollout horizon.
    max_heat_load_interval: Interval = field(default_factory=lambda: Interval(0.0, 0.0))

    # Index of the first violating step.
    first_violation_step: int = -1

    # Physical time of the first violating step.
    first_violation_time_s: float = math.nan

    # True only if heat rate and heat load limits stayed satisfied.
    heat_feasible: bool = True

    # True only if interval propagation stayed valid.
    interval_valid: bool = True

    # Scalar violation magnitude used for optimization penalties.
    violation_amount: float = 0.0

    # Scalar heat penalty accumulated during rollout.
    heat_penalty: float = 0.0

    # Human readable reason for failure.
    failure_reason: str = ""

# Small helpers used by both nominal and interval paths

def _get_param(params: Dict[str, Any], *keys: str, default: Optional[float] = None) -> float:
    """
    Fetch the first matching numeric parameter from a dictionary.
    """
    # Try each alias in order so older notebook parameter names keep working.
    for key in keys:
        if key in params:
            return float(params[key])

    # If nothing matched and no default exists, fail loudly.
    if default is None:
        raise KeyError(f"Missing required parameter. Tried keys={keys}")

    # Otherwise return the caller provided default.
    return float(default)


def hull_intervals(intervals: List[Interval]) -> Interval:
    """
    Return the smallest interval containing all input intervals.
    """
    # Empty input is a real usage error because there is no hull to build.
    if not intervals:
        raise ValueError("interval list cannot be empty")

    # Start from the first interval.
    out = intervals[0]

    # Expand the hull one interval at a time.
    for iv in intervals[1:]:
        out = out.hull(iv)

    return out


def nominal_state_to_interval_box(x_nominal: List[float]) -> List[Interval]:
    """
    Convert a nominal float state into a punctual interval box.
    """
    # The translational state must always have six components.
    if len(x_nominal) != 6:
        raise ValueError("x_nominal must have 6 components")

    # Promote every scalar state into a zero width interval.
    return [promote(float(v)) for v in x_nominal]


def inflate_nominal_state_to_interval_box(
    x_nominal: List[float],
    cfg: IntervalSupervisorConfig,
) -> List[Interval]:
    """
    Inflate a nominal float state into an interval box using configured half widths.
    """
    # The translational state must always have six components.
    if len(x_nominal) != 6:
        raise ValueError("x_nominal must have 6 components")

    # Collect half widths in the same order as the state vector.
    half_widths = [
        float(cfg.r_half_width_m),
        float(cfg.phi_half_width_rad),
        float(cfg.lam_half_width_rad),
        float(cfg.V_half_width_mps),
        float(cfg.gamma_half_width_rad),
        float(cfg.chi_half_width_rad),
    ]

    # Build one interval per state component.
    return [
        Interval(float(x_nominal[i]) - half_widths[i], float(x_nominal[i]) + half_widths[i])
        for i in range(6)
    ]


def interval_component_widths(box: List[Interval]) -> Dict[str, float]:
    """
    Compute per component interval widths for a translational box.
    """
    # The interval state box must match the six state convention.
    if len(box) != 6:
        raise ValueError("box must have 6 components")

    # Return a small name to width dictionary for logging and diagnostics.
    return {STATE_NAMES[i]: float(box[i].width()) for i in range(6)}


def interval_dynamic_pressure_bounds(rho_interval: Interval, V_interval: Interval) -> Interval:
    """
    Compute interval dynamic pressure bounds q = 0.5 rho V squared.
    """
    # Dynamic pressure is computed directly in interval form.
    return promote(0.5) * rho_interval * V_interval.pow_int(2)


def atmosphere_interval_hull_from_state_box(x_box: List[Interval]) -> Dict[str, Any]:
    """
    Collect hulls across all atmospheric layers intersected by the altitude box.
    """
    # The translational interval state must have six components.
    if len(x_box) != 6:
        raise ValueError("x_box must have 6 components")

    # Convert radial distance interval into geometric altitude interval.
    z_geometric_altitude = constants.intv_geometric_altitude(x_box[0])

    # Query the interval atmosphere model across all intersected layers.
    atm_by_layer = AtmosphereModel.intv_US_Standard_ATM(z_geometric_altitude)

    # Above the supported atmosphere range, return zero density and no layers.
    if not atm_by_layer:
        return {
            "altitude_m": z_geometric_altitude,
            "T_K": None,
            "p_Pa": None,
            "rho_kgm3": promote(0.0),
            "layer_indices": [],
        }

    # Preserve the specific layers touched for logging.
    layer_indices = list(atm_by_layer.keys())

    # Gather interval values from every touched layer.
    temperature_intervals = [atm_by_layer[k]["T_K"] for k in layer_indices]
    pressure_intervals = [atm_by_layer[k]["p_Pa"] for k in layer_indices]
    density_intervals = [atm_by_layer[k]["rho_kgm3"] for k in layer_indices]

    # Return layerwise hulls so the rest of the code can work with one density,
    # temperature, and pressure enclosure.
    return {
        "altitude_m": z_geometric_altitude,
        "T_K": hull_intervals(temperature_intervals),
        "p_Pa": hull_intervals(pressure_intervals),
        "rho_kgm3": hull_intervals(density_intervals),
        "layer_indices": layer_indices,
    }


def make_interval_heat_shield() -> Any:
    """
    Build a fresh interval compatible heat shield model.
    """
    # This uses the currently defined heat shield geometry and discretization
    # from constants.py, including the updated ring and sector layout.
    return constants.HeatShield(
        radius_m=float(constants.HEAT_SHIELD_RADIUS_M),
        nose_radius_m=float(constants.HEAT_SHIELD_NOSE_RADIUS_M),
        num_rings=int(constants.HEAT_SHIELD_NUM_RINGS),
        num_sectors=int(constants.HEAT_SHIELD_NUM_SECTORS),
        radial_exp=float(constants.HEAT_SHIELD_RADIAL_EXP),
        azimuthal_gain=float(constants.HEAT_SHIELD_AZIMUTHAL_GAIN),
    )


def clone_interval_heat_shield(heat_shield: Optional[Any]) -> Any:
    """
    Return an independent heat shield object for candidate rollouts.

    Candidate evaluations must not contaminate each other.
    """
    # If no shield exists yet, start from a fresh one.
    if heat_shield is None:
        return make_interval_heat_shield()

    # Prefer an explicit clone method if the class provides one.
    if hasattr(heat_shield, "clone") and callable(getattr(heat_shield, "clone")):
        return heat_shield.clone()

    # Deep copy is a safe fallback for the current heat shield structure.
    return copy.deepcopy(heat_shield)


def interval_heating_envelope(
    rho_interval: Interval,
    V_interval: Interval,
    dt_s: float,
    heat_shield: Optional[Any] = None,
) -> Dict[str, Any]:
    """
    Compute a simple interval heating envelope using the existing heat shield model.
    """
    # Build a shield if the caller did not provide one.
    if heat_shield is None:
        heat_shield = make_interval_heat_shield()

    # Advance the shield thermal state using interval density and speed.
    heat_shield.update(rho=rho_interval, V=V_interval, dt=float(dt_s))

    # Return both the updated shield object and useful summary intervals.
    return {
        "heat_shield": heat_shield,
        "qdot_max": heat_shield.qdot_max(),
        "qdot_mean": heat_shield.qdot_mean(),
        "Q_max": heat_shield.Q_max(),
    }


def classify_interval_against_limits(
    x_interval_new: List[Interval],
    q_interval: Interval,
    cfg: Optional[IntervalSupervisorConfig],
) -> Dict[str, Any]:
    """
    Return a simple safety style summary for interval outputs.
    """
    # Without a config there are no limits to check.
    if cfg is None:
        return {"status": "not_evaluated", "checks": {}}

    # Store individual check results here.
    checks: Dict[str, str] = {}

    # Convert radius into altitude because limits are more intuitive in altitude.
    altitude_interval = constants.intv_geometric_altitude(x_interval_new[0])

    # Speed is already one of the native state components.
    speed_interval = x_interval_new[3]

    def classify_upper(iv: Interval, upper: Optional[float]) -> Optional[str]:
        """
        Classify an interval against an upper bound.
        """
        # No bound means no check.
        if upper is None:
            return None

        # Entire interval lies below or on the limit.
        if iv.hi <= upper:
            return "inside"

        # Entire interval lies above the limit.
        if iv.lo > upper:
            return "outside"

        # Otherwise the interval straddles the limit.
        return "mixed"

    def classify_lower(iv: Interval, lower: Optional[float]) -> Optional[str]:
        """
        Classify an interval against a lower bound.
        """
        # No bound means no check.
        if lower is None:
            return None

        # Entire interval lies above or on the limit.
        if iv.lo >= lower:
            return "inside"

        # Entire interval lies below the limit.
        if iv.hi < lower:
            return "outside"

        # Otherwise the interval straddles the limit.
        return "mixed"

    # Check lower altitude bound if requested.
    altitude_min_status = classify_lower(altitude_interval, cfg.min_altitude_m)
    if altitude_min_status is not None:
        checks["min_altitude_m"] = altitude_min_status

    # Check upper altitude bound if requested.
    altitude_max_status = classify_upper(altitude_interval, cfg.max_altitude_m)
    if altitude_max_status is not None:
        checks["max_altitude_m"] = altitude_max_status

    # Check speed upper bound if requested.
    speed_status = classify_upper(speed_interval, cfg.max_speed_mps)
    if speed_status is not None:
        checks["max_speed_mps"] = speed_status

    # Check dynamic pressure upper bound if requested.
    q_status = classify_upper(q_interval, cfg.max_dynamic_pressure_pa)
    if q_status is not None:
        checks["max_dynamic_pressure_pa"] = q_status

    # If no active checks existed, return a neutral status.
    if not checks:
        return {"status": "not_evaluated", "checks": checks}

    # If any check is fully outside, the overall status is outside.
    if any(v == "outside" for v in checks.values()):
        status = "outside"

    # If nothing is outside but at least one check is mixed, mark mixed.
    elif any(v == "mixed" for v in checks.values()):
        status = "mixed"

    # Otherwise all checks were fully inside.
    else:
        status = "inside"

    return {"status": status, "checks": checks}


def _zero_dx_box() -> List[Interval]:
    """
    Build a six component zero derivative box.
    """
    return [Interval(0.0, 0.0) for _ in range(6)]


def make_invalid_interval_annotation(
    x_interval_old: List[Interval],
    sigma_interval_used: Interval,
    heat_shield: Optional[Any],
    failure_reason: str,
) -> IntervalAnnotationResult:
    """
    Build a safe fallback annotation when live interval logging fails.

    This is only used by the nominal path logging wrapper. The predictor
    corrector candidate rollout still handles interval invalidity explicitly
    so guidance can reject or penalize those candidates.
    """
    # Keep the old interval state so the nominal simulator can continue.
    x_interval_new = [Interval(iv.lo, iv.hi) for iv in x_interval_old]

    # Use zero derivatives because a valid propagated derivative is unavailable.
    dx_interval = _zero_dx_box()

    # Build width summaries from the unchanged box.
    widths_old = interval_component_widths(x_interval_old)
    widths_new = interval_component_widths(x_interval_new)
    dx_widths = interval_component_widths(dx_interval)

    # Altitude can still be computed from the old box.
    altitude_interval = constants.intv_geometric_altitude(x_interval_old[0])

    # Return a conservative placeholder annotation.
    return IntervalAnnotationResult(
        x_interval_old=[Interval(iv.lo, iv.hi) for iv in x_interval_old],
        x_interval_new=x_interval_new,
        dx_interval=dx_interval,
        sigma_interval_used=sigma_interval_used,
        rho_interval=Interval(0.0, 0.0),
        q_interval=Interval(0.0, 0.0),
        altitude_interval=altitude_interval,
        temperature_interval=None,
        pressure_interval=None,
        state_widths_old=widths_old,
        state_widths_new=widths_new,
        dx_widths=dx_widths,
        heating_qdot_max_interval=None,
        heating_qdot_mean_interval=None,
        heating_Q_max_interval=None,
        heat_shield=heat_shield,
        safety_status="invalid_interval",
        safety_checks={"interval_valid": str(failure_reason)},
        layer_indices=[],
    )

# Nominal helper functions for mirrored notebook logic

def nominal_aero_forces_from_state(
    x: List[float],
    sigma_rad: float,
    params: Dict[str, Any],
) -> Dict[str, float]:
    """
    Compute aerodynamic force components from a nominal float state.

    The returned components use the local basis [radial, north, east].
    Additional scalar quantities are returned for the vertical plane and
    lateral lift components that directly feed gamma_dot and chi_dot.
    """
    # Enforce the six state convention.
    if len(x) != 6:
        raise ValueError("x must have 6 state components")

    # Unpack the nominal state.
    r_m, phi_rad, lam_rad, V_mps, gamma_rad, chi_rad = x

    # Latitude and longitude are not needed directly inside the local
    # aerodynamic decomposition, so they are explicitly ignored here.
    del phi_rad, lam_rad

    # Convert radius into altitude for the atmosphere lookup.
    altitude_m = max(0.0, constants.altitude(r_m))

    # Query the nominal atmosphere at the current altitude.
    atm = AtmosphereModel.US_Standard_ATM(altitude_m)

    # Use zero density outside the supported atmosphere range.
    rho = atm["rho_kgm3"] if atm["rho_kgm3"] is not None else 0.0

    # Dynamic pressure is one half rho V squared.
    q_pa = 0.5 * rho * V_mps * V_mps

    # Pull aerodynamic parameters while preserving old alias names.
    surface_area = _get_param(params, "ref_area_m2", "S")
    CD = _get_param(params, "CD_const", "CD")
    CL = _get_param(params, "CL_const", "CL")

    # Convert dynamic pressure into drag and lift magnitudes.
    drag_mag_N = q_pa * surface_area * CD
    lift_mag_N = q_pa * surface_area * CL

    # Precompute sines and cosines used repeatedly in the decomposition.
    sin_gamma = math.sin(gamma_rad)
    cos_gamma = math.cos(gamma_rad)
    sin_chi = math.sin(chi_rad)
    cos_chi = math.cos(chi_rad)
    sin_sigma = math.sin(sigma_rad)
    cos_sigma = math.cos(sigma_rad)

    # Drag acts opposite the velocity vector.
    D_r_N = -drag_mag_N * sin_gamma
    D_north_N = -drag_mag_N * cos_gamma * sin_chi
    D_east_N = -drag_mag_N * cos_gamma * cos_chi

    # Bank rotates lift between the vertical plane and lateral direction.
    lift_vertical_plane_N = lift_mag_N * cos_sigma
    lift_lateral_left_N = -lift_mag_N * sin_sigma

    # Resolve lift into radial, north, and east local components.
    L_r_N = lift_vertical_plane_N * cos_gamma
    L_north_N = -lift_mag_N * (cos_sigma * sin_gamma * sin_chi + sin_sigma * cos_chi)
    L_east_N = lift_mag_N * (-cos_sigma * sin_gamma * cos_chi + sin_sigma * sin_chi)

    # Return a rich dictionary for both physics and logging.
    return {
        "rho_kgm3": float(rho),
        "q_pa": float(q_pa),
        "drag_mag_N": float(drag_mag_N),
        "lift_mag_N": float(lift_mag_N),
        "D_r_N": float(D_r_N),
        "D_north_N": float(D_north_N),
        "D_east_N": float(D_east_N),
        "L_r_N": float(L_r_N),
        "L_north_N": float(L_north_N),
        "L_east_N": float(L_east_N),
        "lift_vertical_plane_N": float(lift_vertical_plane_N),
        "lift_lateral_left_N": float(lift_lateral_left_N),
        "Dtheta": float(D_east_N),
        "Dphi": float(D_north_N),
        "Ltheta": float(L_east_N),
        "Lphi": float(L_north_N),
    }


def nominal_eom_step(
    x: List[float],
    sigma_rad: float,
    params: Dict[str, Any],
    dt_s: float,
) -> Dict[str, Any]:
    """
    One Euler step of the corrected nominal translational dynamics.

    This helper exists so notebook code can mirror the same equations used
    by the interval path instead of keeping a stale duplicate.
    """
    # Enforce the six state convention.
    if len(x) != 6:
        raise ValueError("x must have 6 state components")

    # Negative or zero time step is invalid for forward propagation.
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive")

    # Unpack the nominal state.
    r_m, phi_rad, lam_rad, V_mps, gamma_rad, chi_rad = x

    # Compute nominal aerodynamic forces at the current state.
    aero = nominal_aero_forces_from_state(x=x, sigma_rad=sigma_rad, params=params)

    # Pull mass while preserving old alias names.
    mass_kg = _get_param(params, "mass_kg", "m")

    # Compute gravity at the current radius.
    g_mps2 = constants.gravity(r_m)

    # Precompute trigonometric values used by the equations of motion.
    cos_phi = math.cos(phi_rad)
    sin_phi = math.sin(phi_rad)
    cos_gamma = math.cos(gamma_rad)
    sin_gamma = math.sin(gamma_rad)
    sin_chi = math.sin(chi_rad)
    cos_chi = math.cos(chi_rad)

    # Small epsilon used to avoid singular denominators.
    eps = 1.0e-9

    # Radius and speed are clamped away from zero for numerical safety.
    safe_r = max(r_m, constants.RADIUS_EARTH + 1.0)
    safe_V = max(abs(V_mps), eps)

    # Guard cos phi before using it in longitude rate and tan phi.
    if abs(cos_phi) < eps:
        cos_phi = eps if cos_phi >= 0.0 else -eps

    # Guard cos gamma before using it in heading rate.
    if abs(cos_gamma) < eps:
        cos_gamma = eps if cos_gamma >= 0.0 else -eps

    # Build tan phi from the already guarded cosine.
    tan_phi = sin_phi / cos_phi

    # Radial motion follows the vertical component of velocity.
    r_dot = V_mps * sin_gamma

    # Latitude rate depends on the north component of horizontal velocity.
    phi_dot = (V_mps * cos_gamma * sin_chi) / safe_r

    # Longitude rate depends on the east component of horizontal velocity.
    lam_dot = (V_mps * cos_gamma * cos_chi) / (safe_r * cos_phi)

    # Speed decreases from drag and from climbing against gravity.
    V_dot = -aero["drag_mag_N"] / mass_kg - g_mps2 * sin_gamma

    # Flight path angle changes from vertical plane lift and curvature terms.
    gamma_dot = (
        aero["lift_vertical_plane_N"] / (mass_kg * safe_V)
        + (safe_V / safe_r - g_mps2 / safe_V) * cos_gamma
    )

    # Heading changes from lateral lift and spherical geometry coupling.
    chi_dot = (
        aero["lift_lateral_left_N"] / (mass_kg * safe_V * cos_gamma)
        - (safe_V / safe_r) * cos_gamma * cos_chi * tan_phi
    )

    # Advance the nominal state using one explicit Euler step.
    x_next = [
        r_m + r_dot * dt_s,
        phi_rad + phi_dot * dt_s,
        lam_rad + lam_dot * dt_s,
        V_mps + V_dot * dt_s,
        gamma_rad + gamma_dot * dt_s,
        chi_rad + chi_dot * dt_s,
    ]

    # Return both the advanced state and the derivative information for logging.
    return {
        "x_next": x_next,
        "x_dot": [r_dot, phi_dot, lam_dot, V_dot, gamma_dot, chi_dot],
        "aero": aero,
        "gravity_mps2": g_mps2,
    }

# Interval aerodynamic force model
def intv_aero_forces(
    r: Interval,
    V: Interval,
    gamma: Interval,
    chi: Interval,
    sigma: Interval,
    params: Dict[str, Any],
) -> Dict[int, Dict[str, Interval]]:
    """
    Interval version of the aerodynamic force model.

    Returned components use the local basis [radial, north, east].
    Additional scalars are returned for the vertical plane and lateral lift
    components that directly feed gamma_dot and chi_dot.
    """
    # Convert interval radius into interval altitude.
    z_geometric_altitude = constants.intv_geometric_altitude(r)

    # Query the interval atmosphere across all intersected layers.
    atm_intv = AtmosphereModel.intv_US_Standard_ATM(z_geometric_altitude)

    # If the atmosphere model returns nothing, use zero aero forces.
    if not atm_intv:
        zero = promote(0.0)
        return {
            -1: {
                "q": zero,
                "rho": zero,
                "drag_mag": zero,
                "lift_mag": zero,
                "D_r": zero,
                "D_north": zero,
                "D_east": zero,
                "L_r": zero,
                "L_north": zero,
                "L_east": zero,
                "lift_vertical_plane": zero,
                "lift_lateral_left": zero,
                "Dr": zero,
                "Dtheta": zero,
                "Dphi": zero,
                "Lr": zero,
                "Ltheta": zero,
                "Lphi": zero,
            }
        }

    # Create one aerodynamic dictionary per atmospheric layer.
    aero_forces_by_layer: Dict[int, Dict[str, Interval]] = {}

    # Pull aerodynamic parameters while preserving old alias names.
    CD = _get_param(params, "CD_const", "CD")
    CL = _get_param(params, "CL_const", "CL")
    surface_area = _get_param(params, "ref_area_m2", "S")

    # Precompute trigonometric terms in interval form.
    sin_gamma = gamma.sin()
    cos_gamma = gamma.cos()
    sin_chi = chi.sin()
    cos_chi = chi.cos()
    sin_sigma = sigma.sin()
    cos_sigma = sigma.cos()

    # Process each atmospheric layer independently.
    for k in atm_intv.keys():
        # Pull layer specific density.
        rho = atm_intv[k]["rho_kgm3"] if atm_intv[k]["rho_kgm3"] is not None else promote(0.0)

        # Dynamic pressure in interval form.
        q = promote(0.5) * rho * V.pow_int(2)

        # Drag and lift magnitudes in interval form.
        drag_mag = q * surface_area * CD
        lift_mag = q * surface_area * CL

        # Resolve drag opposite the velocity direction.
        D_r = -drag_mag * sin_gamma
        D_north = -drag_mag * cos_gamma * sin_chi
        D_east = -drag_mag * cos_gamma * cos_chi

        # Resolve lift into vertical plane and lateral pieces using bank.
        lift_vertical_plane = lift_mag * cos_sigma
        lift_lateral_left = -lift_mag * sin_sigma

        # Resolve lift into the local radial, north, and east basis.
        L_r = lift_vertical_plane * cos_gamma
        L_north = -lift_mag * (cos_sigma * sin_gamma * sin_chi + sin_sigma * cos_chi)
        L_east = lift_mag * (-cos_sigma * sin_gamma * cos_chi + sin_sigma * sin_chi)

        # Store the full layerwise result.
        aero_forces_by_layer[k] = {
            "q": q,
            "rho": rho,
            "drag_mag": drag_mag,
            "lift_mag": lift_mag,
            "D_r": D_r,
            "D_north": D_north,
            "D_east": D_east,
            "L_r": L_r,
            "L_north": L_north,
            "L_east": L_east,
            "lift_vertical_plane": lift_vertical_plane,
            "lift_lateral_left": lift_lateral_left,
            "Dr": D_r,
            "Dtheta": D_east,
            "Dphi": D_north,
            "Lr": L_r,
            "Ltheta": L_east,
            "Lphi": L_north,
        }

    return aero_forces_by_layer

# Interval equations of motion

def intv_eom_3d(
    t: float,
    X: List[Interval],
    params: Dict[str, Any],
    sigma: Interval,
) -> Dict[str, Interval]:
    """
    Interval valued 3DOF equations of motion for milestone 1 atmospheric entry.

    The underlying model is the non rotating spherical Earth subset of the Lu
    style bank modulated lift equations. The heading convention remains the
    codebase chi convention, not the north referenced psi used in many papers.
    """
    # Time is not used directly in the current translational model.
    del t

    # Enforce the six component state convention.
    if len(X) != 6:
        raise ValueError("X must have 6 interval components")

    # Unpack the interval state.
    r = X[0]
    phi = X[1]
    lam = X[2]
    V = X[3]
    gamma = X[4]
    chi = X[5]

    # Longitude itself is not used directly inside the local EOM terms.
    del lam

    # Pull mass while preserving old alias names.
    mass_kg = _get_param(params, "mass_kg", "m")

    # Compute interval gravity at the current radius interval.
    g = constants.intv_gravity(r)

    # Compute layerwise interval aerodynamic forces.
    aero = intv_aero_forces(
        r=r,
        V=V,
        gamma=gamma,
        chi=chi,
        sigma=sigma,
        params=params,
    )

    # Build layerwise hulls for the force magnitudes used by the EOM.
    drag_mag = hull_intervals([aero[layer]["drag_mag"] for layer in aero])
    lift_vertical_plane = hull_intervals([aero[layer]["lift_vertical_plane"] for layer in aero])
    lift_lateral_left = hull_intervals([aero[layer]["lift_lateral_left"] for layer in aero])

    # Precompute trigonometric interval terms.
    cos_phi = phi.cos()
    tan_phi = phi.sin() / cos_phi
    cos_gamma = gamma.cos()
    sin_gamma = gamma.sin()
    sin_chi = chi.sin()
    cos_chi = chi.cos()

    # Radial motion follows the vertical component of velocity.
    r_dot = V * sin_gamma

    # Latitude rate comes from the north component of horizontal motion.
    phi_dot = (V * cos_gamma * sin_chi) / r

    # Longitude rate comes from the east component of horizontal motion.
    lam_dot = (V * cos_gamma * cos_chi) / (r * cos_phi)

    # Speed rate comes from drag and gravity projection.
    V_dot = -drag_mag / mass_kg - g * sin_gamma

    # Flight path angle rate comes from vertical plane lift and curvature terms.
    gamma_dot = (lift_vertical_plane / (mass_kg * V)) + (V / r - g / V) * cos_gamma

    # Heading rate comes from lateral lift and spherical geometry coupling.
    chi_dot = (lift_lateral_left / (mass_kg * V * cos_gamma)) - (V / r) * cos_gamma * cos_chi * tan_phi

    # Return the derivative dictionary in named form.
    return {
        "r_dot": r_dot,
        "phi_dot": phi_dot,
        "lam_dot": lam_dot,
        "V_dot": V_dot,
        "gamma_dot": gamma_dot,
        "chi_dot": chi_dot,
    }


def f_interval(
    t: float,
    X: List[Interval],
    params: Dict[str, Any],
    sigma_iv: Interval,
) -> List[Interval]:
    """
    Ordered derivative adapter for interval integration.
    """
    # Convert the named derivative dictionary into state order.
    dx = intv_eom_3d(t=t, X=X, params=params, sigma=sigma_iv)

    return [
        dx["r_dot"],
        dx["phi_dot"],
        dx["lam_dot"],
        dx["V_dot"],
        dx["gamma_dot"],
        dx["chi_dot"],
    ]


# Interval supervisor annotation helpers
def build_interval_annotation(
    x_interval_old: List[Interval],
    params: Dict[str, Any],
    sigma_interval_used: Interval,
    dt_s: Optional[float] = None,
    supervisor_cfg: Optional[IntervalSupervisorConfig] = None,
    heat_shield: Optional[Any] = None,
    t_s: float = 0.0,
) -> IntervalAnnotationResult:
    """
    Propagate one interval translational step and return a rich annotation packet.
    """
    # Use the default entry step when a step size is not supplied.
    if dt_s is None:
        dt_s = float(getattr(constants, "ENTRY_DT_S", 0.25))

    # Enforce the six state convention.
    if len(x_interval_old) != 6:
        raise ValueError("x_interval_old must have 6 interval components")

    # Build atmosphere hull information from the current interval state.
    atmospheric_hull = atmosphere_interval_hull_from_state_box(x_interval_old)
    rho_interval = atmospheric_hull["rho_kgm3"]
    altitude_interval = atmospheric_hull["altitude_m"]
    temperature_interval = atmospheric_hull["T_K"]
    pressure_interval = atmospheric_hull["p_Pa"]

    # Compute interval dynamic pressure from density and speed.
    q_interval = interval_dynamic_pressure_bounds(rho_interval, x_interval_old[3])

    # Compute interval derivatives using the corrected interval EOM.
    dx_interval = f_interval(
        t=float(t_s),
        X=x_interval_old,
        params=params,
        sigma_iv=sigma_interval_used,
    )

    # Advance the interval state with one Euler step.
    x_interval_new = box_add(x_interval_old, box_scalar_mul(float(dt_s), dx_interval))

    # Default heating outputs to None unless heating is requested.
    heating_qdot_max_interval = None
    heating_qdot_mean_interval = None
    heating_Q_max_interval = None

    # Propagate heating if the config requests it or if a live heat shield exists.
    if (supervisor_cfg is not None and supervisor_cfg.include_heating) or (heat_shield is not None):
        heating_info = interval_heating_envelope(
            rho_interval=rho_interval,
            V_interval=x_interval_old[3],
            dt_s=float(dt_s),
            heat_shield=heat_shield,
        )

        # Keep the updated shield object for the caller.
        heat_shield = heating_info["heat_shield"]

        # Expose useful heating summaries for logging and predictor screening.
        heating_qdot_max_interval = heating_info["qdot_max"]
        heating_qdot_mean_interval = heating_info["qdot_mean"]
        heating_Q_max_interval = heating_info["Q_max"]

    # Run simple interval bound checks if a config was supplied.
    safety_info = classify_interval_against_limits(
        x_interval_new=x_interval_new,
        q_interval=q_interval,
        cfg=supervisor_cfg,
    )

    # Return the complete annotation packet.
    return IntervalAnnotationResult(
        x_interval_old=list(x_interval_old),
        x_interval_new=x_interval_new,
        dx_interval=dx_interval,
        sigma_interval_used=sigma_interval_used,
        rho_interval=rho_interval,
        q_interval=q_interval,
        altitude_interval=altitude_interval,
        temperature_interval=temperature_interval,
        pressure_interval=pressure_interval,
        state_widths_old=interval_component_widths(x_interval_old),
        state_widths_new=interval_component_widths(x_interval_new),
        dx_widths=interval_component_widths(dx_interval),
        heating_qdot_max_interval=heating_qdot_max_interval,
        heating_qdot_mean_interval=heating_qdot_mean_interval,
        heating_Q_max_interval=heating_Q_max_interval,
        heat_shield=heat_shield,
        safety_status=safety_info["status"],
        safety_checks=safety_info["checks"],
        layer_indices=list(atmospheric_hull["layer_indices"]),
    )


def annotate_nominal_state_with_interval_supervisor(
    x_nominal_old: List[float],
    params: Dict[str, Any],
    sigma_actual_after_rad: float,
    x_interval_old: Optional[List[Interval]] = None,
    supervisor_cfg: Optional[IntervalSupervisorConfig] = None,
    dt_s: Optional[float] = None,
    heat_shield: Optional[Any] = None,
    t_s: float = 0.0,
) -> IntervalAnnotationResult:
    """
    Build an interval annotation for a nominal translational step.

    This wrapper is used by the live nominal simulation path. Unlike the
    predictor corrector candidate rollout, this wrapper catches interval
    propagation failures and returns a safe fallback packet so the nominal
    trajectory can continue running and logging.
    """
    # If the caller did not provide an interval state, create one now.
    if x_interval_old is None:
        if supervisor_cfg is None:
            x_interval_old = nominal_state_to_interval_box(x_nominal_old)
        else:
            x_interval_old = inflate_nominal_state_to_interval_box(x_nominal_old, supervisor_cfg)

    # Promote the flown bank angle into an interval.
    sigma_interval_used = promote(float(sigma_actual_after_rad))

    try:
        # Normal path. Build the real interval annotation.
        return build_interval_annotation(
            x_interval_old=x_interval_old,
            params=params,
            sigma_interval_used=sigma_interval_used,
            dt_s=dt_s,
            supervisor_cfg=supervisor_cfg,
            heat_shield=heat_shield,
            t_s=t_s,
        )
    except Exception as exc:
        # Fallback path. Preserve the nominal run and surface the interval failure.
        return make_invalid_interval_annotation(
            x_interval_old=x_interval_old,
            sigma_interval_used=sigma_interval_used,
            heat_shield=heat_shield,
            failure_reason=f"live_interval_annotation_invalid: {exc}",
        )


# Short horizon rollout utilities for guidance
def _positive_violation_amount(interval_value: Interval, limit_value: Optional[float]) -> float:
    """
    Return normalized positive violation amount using the interval upper bound.
    """
    # Without a limit there is no violation.
    if limit_value is None:
        return 0.0

    # Normalize by the limit magnitude so different limits scale sensibly.
    denom = max(abs(float(limit_value)), 1.0)

    # Only positive excess above the limit counts as a violation.
    return max(0.0, (float(interval_value.hi) - float(limit_value)) / denom)


def _hull_or_default(a: Optional[Interval], b: Optional[Interval]) -> Interval:
    """
    Merge two optional intervals.
    """
    # If neither interval exists yet, return a zero interval.
    if a is None and b is None:
        return Interval(0.0, 0.0)

    # If only one interval exists, return a copy of it.
    if a is None:
        return Interval(float(b.lo), float(b.hi))
    if b is None:
        return Interval(float(a.lo), float(a.hi))

    # Otherwise build the hull.
    return a.hull(b)


def rollout_predictor_corrector_candidate(
    x_nominal_old: List[float],
    sigma_cmd_rad: float,
    params: Dict[str, Any],
    dt_s: float,
    horizon_steps: int,
    x_interval_old: Optional[List[Interval]] = None,
    supervisor_cfg: Optional[IntervalSupervisorConfig] = None,
    heat_shield: Optional[Any] = None,
    t0_s: float = 0.0,
    heat_rate_limit: Optional[float] = None,
    heat_load_limit: Optional[float] = None,
) -> PredictorCorrectorRolloutResult:
    """
    Roll out one candidate sigma command for the guidance predictor corrector.

    This utility is the heat aware prediction engine used by control.py.
    The rollout propagates nominal state, interval state, and the interval
    heat shield forward together so guidance can reject or penalize heat
    infeasible candidates.

    Important milestone 1 note

    The translational propagation is the corrected current 3DOF model.
    The heat feasibility layer uses the current interval heat shield model.
    This is a deliberate approximation rather than a literal implementation
    of every Lu equation.
    """
    # Enforce the six state convention.
    if len(x_nominal_old) != 6:
        raise ValueError("x_nominal_old must have 6 components")

    # Require a positive step size.
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive")

    # Require a positive horizon length.
    if horizon_steps <= 0:
        raise ValueError("horizon_steps must be positive")

    # Build a local supervisor config that definitely enables heating when
    # heat limits are being used during candidate evaluation.
    local_supervisor_cfg = supervisor_cfg

    if local_supervisor_cfg is None and (heat_rate_limit is not None or heat_load_limit is not None):
        local_supervisor_cfg = IntervalSupervisorConfig(include_heating=True)

    elif local_supervisor_cfg is not None and not local_supervisor_cfg.include_heating:
        local_supervisor_cfg = copy.deepcopy(local_supervisor_cfg)
        local_supervisor_cfg.include_heating = True

    # Build or copy the interval state used by this candidate rollout.
    if x_interval_old is None:
        if local_supervisor_cfg is None:
            x_interval = nominal_state_to_interval_box(x_nominal_old)
        else:
            x_interval = inflate_nominal_state_to_interval_box(x_nominal_old, local_supervisor_cfg)
    else:
        x_interval = [Interval(iv.lo, iv.hi) for iv in x_interval_old]

    # Copy the nominal state so this candidate cannot alter the caller state.
    x_nominal = [float(v) for v in x_nominal_old]

    # Promote the candidate sigma command into interval form once.
    sigma_iv = promote(float(sigma_cmd_rad))

    # Clone the heat shield so each candidate evolves independently.
    candidate_heat_shield = clone_interval_heat_shield(heat_shield)

    # Store the full sequence of predicted nominal and interval states.
    nominal_states: List[List[float]] = [list(x_nominal)]
    interval_states: List[List[Interval]] = [list(x_interval)]

    # Store step by step summaries for debug logging.
    step_summaries: List[RolloutStepSummary] = []

    # These track the overall worst heating across the entire horizon.
    max_heating_rate_interval: Optional[Interval] = None
    max_heat_load_interval: Optional[Interval] = None

    # These fields summarize candidate feasibility.
    heat_feasible = True
    interval_valid = True
    violation_amount = 0.0
    heat_penalty = 0.0
    failure_reason = ""
    first_violation_step = -1
    first_violation_time_s = math.nan

    # Propagate the candidate one short horizon step at a time.
    for step_idx in range(int(horizon_steps)):
        # Compute the physical time associated with this rollout step.
        t_step_s = float(t0_s + step_idx * dt_s)

        # Save the nominal state before stepping for debug visibility.
        x_nom_before = list(x_nominal)

        # Nominal prediction always uses the corrected translational model.
        nominal_step = nominal_eom_step(
            x=x_nominal,
            sigma_rad=float(sigma_cmd_rad),
            params=params,
            dt_s=float(dt_s),
        )

        # Advance the nominal state.
        x_nominal = [float(v) for v in nominal_step["x_next"]]

        try:
            # Propagate the interval state and interval heating together.
            annotation = build_interval_annotation(
                x_interval_old=x_interval,
                params=params,
                sigma_interval_used=sigma_iv,
                dt_s=float(dt_s),
                supervisor_cfg=local_supervisor_cfg,
                heat_shield=candidate_heat_shield,
                t_s=float(t_step_s),
            )

            # Update the working interval state for the next step.
            x_interval = annotation.x_interval_new

            # Keep the updated heat shield for the next step.
            candidate_heat_shield = annotation.heat_shield

            # Pull heating summaries if available.
            qdot_interval = (
                annotation.heating_qdot_max_interval
                if annotation.heating_qdot_max_interval is not None
                else Interval(0.0, 0.0)
            )
            Q_interval = (
                annotation.heating_Q_max_interval
                if annotation.heating_Q_max_interval is not None
                else Interval(0.0, 0.0)
            )

            # Update the horizon wide hulls.
            max_heating_rate_interval = _hull_or_default(max_heating_rate_interval, qdot_interval)
            max_heat_load_interval = _hull_or_default(max_heat_load_interval, Q_interval)

            # Track stepwise violation size and reason labels.
            step_violation = 0.0
            step_reason_parts: List[str] = []

            # Check heating rate limit if one was provided.
            if heat_rate_limit is not None:
                rate_violation = _positive_violation_amount(qdot_interval, float(heat_rate_limit))
                if rate_violation > 0.0:
                    step_violation += rate_violation
                    step_reason_parts.append("heat_rate_limit")

            # Check integrated heat load limit if one was provided.
            if heat_load_limit is not None:
                load_violation = _positive_violation_amount(Q_interval, float(heat_load_limit))
                if load_violation > 0.0:
                    step_violation += load_violation
                    step_reason_parts.append("heat_load_limit")

            # If any limit was violated, mark the candidate infeasible.
            if step_violation > 0.0:
                heat_feasible = False

                # Keep the maximum scalar violation across all steps.
                violation_amount = max(violation_amount, float(step_violation))

                # Accumulate heat penalty so guidance can optimize against it.
                heat_penalty += float(step_violation)

                # Record the first violating step only once.
                if first_violation_step < 0:
                    first_violation_step = int(step_idx)
                    first_violation_time_s = float(t_step_s)

                # Preserve the first failure reason for easy logging.
                if not failure_reason:
                    failure_reason = ",".join(step_reason_parts)

            # Build the per step summary packet for this valid step.
            step_summary = RolloutStepSummary(
                step_index=int(step_idx),
                time_s=float(t_step_s),
                x_nominal_before=list(x_nom_before),
                x_nominal_after=list(x_nominal),
                x_interval_after=list(x_interval),
                max_heating_rate_interval=qdot_interval,
                max_heat_load_interval=Q_interval,
                heat_feasible_after_step=bool(heat_feasible),
                interval_valid_after_step=True,
                failure_reason=str(failure_reason),
            )

        except Exception as exc:
            # Any denominator or interval failure makes this candidate invalid.
            interval_valid = False
            heat_feasible = False

            # Store a readable failure reason for logging and RL export.
            failure_reason = f"interval_rollout_invalid: {exc}"

            # Force a nonzero violation so the candidate is strongly penalized.
            violation_amount = max(violation_amount, 1.0)
            heat_penalty += 1.0

            # Record the first invalid step.
            if first_violation_step < 0:
                first_violation_step = int(step_idx)
                first_violation_time_s = float(t_step_s)

            # Build a final step summary that records invalid propagation.
            step_summary = RolloutStepSummary(
                step_index=int(step_idx),
                time_s=float(t_step_s),
                x_nominal_before=list(x_nom_before),
                x_nominal_after=list(x_nominal),
                x_interval_after=None,
                max_heating_rate_interval=_hull_or_default(max_heating_rate_interval, None),
                max_heat_load_interval=_hull_or_default(max_heat_load_interval, None),
                heat_feasible_after_step=False,
                interval_valid_after_step=False,
                failure_reason=str(failure_reason),
            )

            # Record the failing step and stop this candidate rollout.
            step_summaries.append(step_summary)
            break

        # Append successful step outputs after the try block.
        nominal_states.append(list(x_nominal))
        interval_states.append(list(x_interval))
        step_summaries.append(step_summary)

    # If no heating was computed at all, return zero intervals instead of None.
    if max_heating_rate_interval is None:
        max_heating_rate_interval = Interval(0.0, 0.0)
    if max_heat_load_interval is None:
        max_heat_load_interval = Interval(0.0, 0.0)

    # Build the final candidate rollout result for guidance.
    return PredictorCorrectorRolloutResult(
        sigma_cmd_rad=float(sigma_cmd_rad),
        horizon_steps=int(horizon_steps),
        nominal_states=nominal_states,
        interval_states=interval_states,
        step_summaries=step_summaries,
        final_nominal_state=list(x_nominal),
        final_interval_state=list(x_interval) if interval_valid else None,
        final_heat_shield=candidate_heat_shield,
        max_heating_rate_interval=max_heating_rate_interval,
        max_heat_load_interval=max_heat_load_interval,
        first_violation_step=int(first_violation_step),
        first_violation_time_s=float(first_violation_time_s),
        heat_feasible=bool(heat_feasible),
        interval_valid=bool(interval_valid),
        violation_amount=float(violation_amount),
        heat_penalty=float(heat_penalty),
        failure_reason=str(failure_reason),
    )