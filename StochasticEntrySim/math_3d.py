"""
interval and nominal translational reentry dynamics for milestone one

state order

r
phi
lam
v
gamma
chi
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import AtmosphereModel
import constants
from interval_math import Interval, box_add, box_scalar_mul, box_split, promote


STATE_NAMES = ("r", "phi", "lam", "V", "gamma", "chi")


@dataclass
class IntervalSupervisorConfig:
    """
    store the live interval supervisor settings

    this includes the initial half widths heat limits recenter logic
    split logic and denominator safety thresholds

    input
    scalar configuration values

    output
    one config object used by interval propagation
    """

    r_half_width_m: float = 0.0
    phi_half_width_rad: float = 0.0
    lam_half_width_rad: float = 0.0
    V_half_width_mps: float = 0.0
    gamma_half_width_rad: float = 0.0
    chi_half_width_rad: float = 0.0

    min_altitude_m: Optional[float] = None
    max_altitude_m: Optional[float] = None
    max_speed_mps: Optional[float] = None
    max_dynamic_pressure_pa: Optional[float] = None

    include_heating: bool = False
    heat_rate_limit: Optional[float] = constants.HEAT_RATE_LIMIT_DEFAULT
    heat_load_limit: Optional[float] = constants.HEAT_LOAD_LIMIT_DEFAULT

    interval_recenter_enabled: bool = constants.INTERVAL_RECENTER_ENABLED
    interval_recenter_use_cadence: bool = getattr(constants, "INTERVAL_RECENTER_USE_CADENCE", True)
    interval_recenter_cadence_s: float = constants.INTERVAL_RECENTER_CADENCE_S
    interval_recenter_width_thresholds: Dict[str, float] = field(
        default_factory=lambda: dict(constants.INTERVAL_RECENTER_WIDTH_THRESHOLDS)
    )

    interval_box_split_enabled: bool = constants.INTERVAL_BOX_SPLIT_ENABLED
    interval_box_split_max_depth: int = constants.INTERVAL_BOX_SPLIT_MAX_DEPTH
    interval_box_split_width_thresholds: Dict[str, float] = field(
        default_factory=lambda: dict(constants.INTERVAL_BOX_SPLIT_WIDTH_THRESHOLDS)
    )

    interval_denominator_safety_V_mps: float = constants.INTERVAL_DENOMINATOR_SAFETY_V_MPS
    interval_denominator_safety_cos_gamma: float = constants.INTERVAL_DENOMINATOR_SAFETY_COS_GAMMA
    interval_denominator_safety_cos_phi: float = constants.INTERVAL_DENOMINATOR_SAFETY_COS_PHI


@dataclass
class IntervalAnnotationResult:
    """
    store the result of one interval annotation step

    this tracks the new interval box derivative widths heat summary
    safety flags and whether recenter or splitting was used

    input
    values created by one interval propagation step

    output
    one structured result object
    """

    x_interval_old: List[Interval]
    x_interval_new: List[Interval]
    dx_interval: List[Interval]

    sigma_interval_used: Interval
    rho_interval: Interval
    q_interval: Interval
    altitude_interval: Interval
    temperature_interval: Optional[Interval]
    pressure_interval: Optional[Interval]

    state_widths_old: Dict[str, float] = field(default_factory=dict)
    state_widths_new: Dict[str, float] = field(default_factory=dict)
    dx_widths: Dict[str, float] = field(default_factory=dict)

    heating_qdot_max_interval: Optional[Interval] = None
    heating_qdot_mean_interval: Optional[Interval] = None
    heating_Q_max_interval: Optional[Interval] = None
    heat_shield: Optional[Any] = None

    safety_status: str = "not_evaluated"
    safety_checks: Dict[str, str] = field(default_factory=dict)
    layer_indices: List[int] = field(default_factory=list)

    interval_valid: bool = True
    interval_active: bool = True
    interval_failed: bool = False
    interval_failure_kind: str = "none"
    interval_failure_reason: str = ""
    interval_numerical_failure: bool = False
    interval_heat_violation: bool = False

    recentered_this_step: bool = False
    recenter_reason: str = ""
    split_this_step: bool = False
    split_reason: str = ""
    split_depth_used: int = 0
    split_child_count: int = 0

    last_valid_interval_state: Optional[List[Interval]] = None


@dataclass
class RolloutStepSummary:
    """
    store one predictor rollout step summary

    input
    nominal interval and heat data for one rollout step

    output
    one step packet used by guidance logs
    """

    step_index: int
    time_s: float
    x_nominal_before: List[float]
    x_nominal_after: List[float]
    x_interval_after: Optional[List[Interval]]
    max_heating_rate_interval: Interval
    max_heat_load_interval: Interval
    heat_feasible_after_step: bool
    interval_valid_after_step: bool
    failure_reason: str = ""


@dataclass
class PredictorCorrectorRolloutResult:
    """
    store one predictor corrector candidate rollout

    input
    all states and heat results from a short horizon candidate simulation

    output
    one object used by guidance scoring
    """

    sigma_cmd_rad: float
    horizon_steps: int
    nominal_states: List[List[float]] = field(default_factory=list)
    interval_states: List[List[Interval]] = field(default_factory=list)
    step_summaries: List[RolloutStepSummary] = field(default_factory=list)

    final_nominal_state: List[float] = field(default_factory=list)
    final_interval_state: Optional[List[Interval]] = None
    final_heat_shield: Optional[Any] = None

    max_heating_rate_interval: Interval = field(default_factory=lambda: Interval(0.0, 0.0))
    max_heat_load_interval: Interval = field(default_factory=lambda: Interval(0.0, 0.0))

    first_violation_step: int = -1
    first_violation_time_s: float = math.nan

    heat_feasible: bool = True
    interval_valid: bool = True
    interval_screen_used: bool = True
    interval_failure_kind: str = "none"
    violation_amount: float = 0.0
    heat_penalty: float = 0.0
    failure_reason: str = ""


def _get_param(params: Dict[str, Any], *keys: str, default: Optional[float] = None) -> float:
    """
    fetch the first matching numeric parameter

    input
    parameter dictionary and possible key names

    output
    one float value
    """
    for key in keys:
        if key in params:
            return float(params[key])

    if default is None:
        raise KeyError(f"Missing required parameter. Tried keys={keys}")

    return float(default)


def hull_intervals(intervals: List[Interval]) -> Interval:
    """
    build the smallest interval that contains all input intervals

    input
    a nonempty list of intervals

    output
    one hull interval
    """
    if not intervals:
        raise ValueError("interval list cannot be empty")

    out = intervals[0]
    for iv in intervals[1:]:
        out = out.hull(iv)
    return out


def nominal_state_to_interval_box(x_nominal: List[float]) -> List[Interval]:
    """
    convert a nominal state into a zero width interval box

    input
    six component nominal state

    output
    six component interval box
    """
    if len(x_nominal) != 6:
        raise ValueError("x_nominal must have 6 components")

    return [promote(float(v)) for v in x_nominal]


def inflate_nominal_state_to_interval_box(
    x_nominal: List[float],
    cfg: IntervalSupervisorConfig,
) -> List[Interval]:
    """
    inflate a nominal state into an interval box using configured half widths

    input
    six component nominal state and supervisor config

    output
    six component interval box
    """
    if len(x_nominal) != 6:
        raise ValueError("x_nominal must have 6 components")

    half_widths = [
        float(cfg.r_half_width_m),
        float(cfg.phi_half_width_rad),
        float(cfg.lam_half_width_rad),
        float(cfg.V_half_width_mps),
        float(cfg.gamma_half_width_rad),
        float(cfg.chi_half_width_rad),
    ]

    return [
        Interval(float(x_nominal[i]) - half_widths[i], float(x_nominal[i]) + half_widths[i])
        for i in range(6)
    ]


def interval_component_widths(box: List[Interval]) -> Dict[str, float]:
    """
    compute one width value for each state component

    input
    six component interval box

    output
    mapping from state name to width
    """
    if len(box) != 6:
        raise ValueError("box must have 6 components")

    return {STATE_NAMES[i]: float(box[i].width()) for i in range(6)}


def copy_interval(iv: Interval) -> Interval:
    """
    make a shallow copy of one interval

    input
    one interval

    output
    copied interval
    """
    return Interval(float(iv.lo), float(iv.hi))


def copy_interval_box(box: List[Interval]) -> List[Interval]:
    """
    make a shallow copy of one interval box

    input
    one interval box

    output
    copied interval box
    """
    return [copy_interval(iv) for iv in box]


def hull_interval_boxes(boxes: List[List[Interval]]) -> List[Interval]:
    """
    build the componentwise hull across many interval boxes

    input
    nonempty list of interval boxes

    output
    one hull box
    """
    if not boxes:
        raise ValueError("boxes cannot be empty")

    return [
        hull_intervals([box[i] for box in boxes])
        for i in range(len(boxes[0]))
    ]


def _interval_abs_min(iv: Interval) -> float:
    """
    return the smallest absolute value that may appear in an interval

    input
    one interval

    output
    scalar absolute minimum
    """
    if iv.lo <= 0.0 <= iv.hi:
        return 0.0
    return min(abs(float(iv.lo)), abs(float(iv.hi)))


def _crosses_cadence_boundary(t_s: float, dt_s: float, cadence_s: float) -> bool:
    """
    detect whether a step crosses a cadence bucket boundary

    input
    current time step size and cadence size

    output
    true when a new cadence bucket is entered
    """
    if cadence_s <= 0.0:
        return False

    prev_bin = int(max(float(t_s), 0.0) / float(cadence_s))
    next_bin = int(max(float(t_s) + float(dt_s), 0.0) / float(cadence_s))
    return next_bin > prev_bin


def interval_width_exceeds_thresholds(
    box: List[Interval],
    thresholds: Optional[Dict[str, float]],
) -> bool:
    """
    check whether any state width exceeds its configured threshold

    input
    interval box and optional width limits

    output
    true when any width is too large
    """
    if not thresholds:
        return False

    widths = interval_component_widths(box)
    for name, width in widths.items():
        limit = thresholds.get(name)
        if limit is not None and float(width) > float(limit):
            return True

    return False


def dangerous_interval_denominator_flags(
    box: List[Interval],
    cfg: Optional[IntervalSupervisorConfig],
) -> Dict[str, bool]:
    """
    flag denominator regions that are getting close to singular behavior

    the sensitive terms are v cos gamma and cos phi

    input
    interval box and supervisor config

    output
    dictionary of boolean danger flags
    """
    if cfg is None:
        return {
            "V": False,
            "cos_gamma": False,
            "cos_phi": False,
        }

    V_iv = box[3]
    gamma_iv = box[4]
    phi_iv = box[1]

    cos_gamma_iv = gamma_iv.cos()
    cos_phi_iv = phi_iv.cos()

    return {
        "V": float(V_iv.lo) <= float(cfg.interval_denominator_safety_V_mps),
        "cos_gamma": _interval_abs_min(cos_gamma_iv) <= float(cfg.interval_denominator_safety_cos_gamma),
        "cos_phi": _interval_abs_min(cos_phi_iv) <= float(cfg.interval_denominator_safety_cos_phi),
    }


def should_recenter_interval_box(
    x_interval_old: List[Interval],
    x_nominal_center: Optional[List[float]],
    cfg: Optional[IntervalSupervisorConfig],
    t_s: float,
    dt_s: float,
) -> Tuple[bool, str]:
    """
    decide whether the live interval box should be rebuilt around nominal

    this now separates width based health recenter from fixed cadence recenter
    so rl data can ignore maintenance resets when cadence is disabled

    input
    current interval box nominal center config current time and step size

    output
    recenter flag and short reason string
    """
    if cfg is None or not bool(cfg.interval_recenter_enabled):
        return False, ""

    if x_nominal_center is None:
        return False, ""

    if interval_width_exceeds_thresholds(x_interval_old, cfg.interval_recenter_width_thresholds):
        return True, "width_threshold"

    # this block is the main fix
    # cadence recenter only runs when the new config flag allows it
    if bool(cfg.interval_recenter_use_cadence) and _crosses_cadence_boundary(
        t_s=float(t_s),
        dt_s=float(dt_s),
        cadence_s=float(cfg.interval_recenter_cadence_s),
    ):
        return True, "cadence"

    return False, ""


def recenter_interval_box_around_nominal(
    x_nominal_center: List[float],
    cfg: Optional[IntervalSupervisorConfig],
) -> List[Interval]:
    """
    rebuild the live interval box around the current nominal center

    input
    nominal center and optional config

    output
    new interval box centered on nominal
    """
    if cfg is None:
        return nominal_state_to_interval_box(x_nominal_center)

    return inflate_nominal_state_to_interval_box(x_nominal_center, cfg)


def choose_split_dimension(
    box: List[Interval],
    cfg: Optional[IntervalSupervisorConfig],
) -> int:
    """
    choose the dimension used for adaptive splitting

    denominator danger takes priority
    otherwise the widest dimension is split

    input
    interval box and config

    output
    split index
    """
    flags = dangerous_interval_denominator_flags(box, cfg)

    if flags.get("cos_gamma", False):
        return 4
    if flags.get("cos_phi", False):
        return 1
    if flags.get("V", False):
        return 3

    widths = [iv.width() for iv in box]
    return int(max(range(len(widths)), key=lambda i: widths[i]))


def interval_box_needs_split(
    box: List[Interval],
    cfg: Optional[IntervalSupervisorConfig],
) -> bool:
    """
    decide whether adaptive splitting should be attempted

    input
    interval box and config

    output
    true when the box is too wide or near a dangerous denominator
    """
    if cfg is None or not bool(cfg.interval_box_split_enabled):
        return False

    if interval_width_exceeds_thresholds(box, cfg.interval_box_split_width_thresholds):
        return True

    flags = dangerous_interval_denominator_flags(box, cfg)
    return any(bool(v) for v in flags.values())


def split_reason_from_box(
    box: List[Interval],
    cfg: Optional[IntervalSupervisorConfig],
) -> str:
    """
    describe why adaptive splitting is needed

    input
    interval box and config

    output
    short split reason string
    """
    if cfg is None:
        return ""

    if interval_width_exceeds_thresholds(box, cfg.interval_box_split_width_thresholds):
        return "width_threshold"

    flags = dangerous_interval_denominator_flags(box, cfg)
    if bool(flags.get("cos_gamma", False)):
        return "denominator_cos_gamma"
    if bool(flags.get("cos_phi", False)):
        return "denominator_cos_phi"
    if bool(flags.get("V", False)):
        return "denominator_V"

    return ""


def hull_heat_shields(heat_shields: List[Any]) -> Optional[Any]:
    """
    hull many heat shield states cell by cell

    input
    list of heat shield objects

    output
    one hull heat shield
    """
    if not heat_shields:
        return None

    out = clone_interval_heat_shield(heat_shields[0])

    for shield in heat_shields[1:]:
        for i in range(len(out.qdot)):
            out.qdot[i] = out.qdot[i].hull(shield.qdot[i])
            out.Q[i] = out.Q[i].hull(shield.Q[i])

    return out


def classify_heat_violation(
    heating_qdot_max_interval: Optional[Interval],
    heating_Q_max_interval: Optional[Interval],
    cfg: Optional[IntervalSupervisorConfig],
) -> Dict[str, Any]:
    """
    classify interval heating against the configured hard limits

    input
    heat rate interval heat load interval and config

    output
    dictionary with violation flags and a short reason
    """
    rate_violation = False
    load_violation = False
    reason_parts: List[str] = []

    if cfg is not None and heating_qdot_max_interval is not None and cfg.heat_rate_limit is not None:
        if float(heating_qdot_max_interval.hi) > float(cfg.heat_rate_limit):
            rate_violation = True
            reason_parts.append("interval_heat_rate_limit")

    if cfg is not None and heating_Q_max_interval is not None and cfg.heat_load_limit is not None:
        if float(heating_Q_max_interval.hi) > float(cfg.heat_load_limit):
            load_violation = True
            reason_parts.append("interval_heat_load_limit")

    return {
        "heat_violation": bool(rate_violation or load_violation),
        "rate_violation": bool(rate_violation),
        "load_violation": bool(load_violation),
        "reason": ",".join(reason_parts),
    }


def nominal_heating_envelope_from_state(
    x_nominal: List[float],
    dt_s: float,
    heat_shield: Optional[Any] = None,
) -> Dict[str, Any]:
    """
    advance the heat shield using only the nominal state

    input
    nominal state step size and optional shield

    output
    updated heat shield and nominal heat summaries
    """
    if len(x_nominal) != 6:
        raise ValueError("x_nominal must have 6 components")

    altitude_m = max(0.0, constants.altitude(float(x_nominal[0])))
    atm = AtmosphereModel.US_Standard_ATM(altitude_m)
    rho = atm["rho_kgm3"] if atm["rho_kgm3"] is not None else 0.0

    if heat_shield is None:
        heat_shield = make_interval_heat_shield()

    heat_shield.update(rho=float(rho), V=float(x_nominal[3]), dt=float(dt_s))

    return {
        "heat_shield": heat_shield,
        "qdot_max": heat_shield.qdot_max(),
        "qdot_mean": heat_shield.qdot_mean(),
        "Q_max": heat_shield.Q_max(),
        "rho_kgm3": float(rho),
        # Stagnation-point heating component breakdown (convective vs
        # radiative) plus the radiative-equilibrium wall temperature.
        "qdot_conv_stag": heat_shield.qdot_conv_stag,
        "qdot_rad_stag": heat_shield.qdot_rad_stag,
        "qdot_total_stag": heat_shield.qdot_total_stag,
        "T_wall_stag_K": heat_shield.T_wall_stag_K,
    }


def interval_dynamic_pressure_bounds(rho_interval: Interval, V_interval: Interval) -> Interval:
    """
    compute interval dynamic pressure

    equation
    q equals one half rho v squared

    input
    density interval and speed interval

    output
    dynamic pressure interval
    """
    return promote(0.5) * rho_interval * V_interval.pow_int(2)


def atmosphere_interval_hull_from_state_box(x_box: List[Interval]) -> Dict[str, Any]:
    """
    collect atmosphere hulls across all layers touched by the altitude interval

    input
    interval state box

    output
    density pressure temperature and layer indices
    """
    if len(x_box) != 6:
        raise ValueError("x_box must have 6 components")

    z_geometric_altitude = constants.intv_geometric_altitude(x_box[0])
    atm_by_layer = AtmosphereModel.intv_US_Standard_ATM(z_geometric_altitude)

    if not atm_by_layer:
        return {
            "altitude_m": z_geometric_altitude,
            "T_K": None,
            "p_Pa": None,
            "rho_kgm3": promote(0.0),
            "layer_indices": [],
        }

    layer_indices = list(atm_by_layer.keys())
    temperature_intervals = [atm_by_layer[k]["T_K"] for k in layer_indices if atm_by_layer[k]["T_K"] is not None]
    pressure_intervals = [atm_by_layer[k]["p_Pa"] for k in layer_indices if atm_by_layer[k]["p_Pa"] is not None]
    density_intervals = [atm_by_layer[k]["rho_kgm3"] for k in layer_indices]

    return {
        "altitude_m": z_geometric_altitude,
        "T_K": hull_intervals(temperature_intervals) if temperature_intervals else None,
        "p_Pa": hull_intervals(pressure_intervals) if pressure_intervals else None,
        "rho_kgm3": hull_intervals(density_intervals),
        "layer_indices": layer_indices,
    }


def make_interval_heat_shield() -> Any:
    """
    build a fresh interval compatible heat shield

    input
    none

    output
    heat shield object
    """
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
    return an independent heat shield for candidate rollouts

    input
    optional source heat shield

    output
    copied heat shield
    """
    if heat_shield is None:
        return make_interval_heat_shield()

    if hasattr(heat_shield, "clone") and callable(getattr(heat_shield, "clone")):
        return heat_shield.clone()

    return copy.deepcopy(heat_shield)


def interval_heating_envelope(
    rho_interval: Interval,
    V_interval: Interval,
    dt_s: float,
    heat_shield: Optional[Any] = None,
) -> Dict[str, Any]:
    """
    propagate the interval heat shield state for one step

    input
    density interval speed interval step size and optional shield

    output
    updated shield and heat summaries
    """
    if heat_shield is None:
        heat_shield = make_interval_heat_shield()

    heat_shield.update(rho=rho_interval, V=V_interval, dt=float(dt_s))

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
    classify the propagated interval against optional simple bounds

    input
    new interval box dynamic pressure interval and config

    output
    summary status and per limit checks
    """
    if cfg is None:
        return {"status": "not_evaluated", "checks": {}}

    checks: Dict[str, str] = {}
    altitude_interval = constants.intv_geometric_altitude(x_interval_new[0])
    speed_interval = x_interval_new[3]

    def classify_upper(iv: Interval, upper: Optional[float]) -> Optional[str]:
        if upper is None:
            return None
        if iv.hi <= upper:
            return "inside"
        if iv.lo > upper:
            return "outside"
        return "mixed"

    def classify_lower(iv: Interval, lower: Optional[float]) -> Optional[str]:
        if lower is None:
            return None
        if iv.lo >= lower:
            return "inside"
        if iv.hi < lower:
            return "outside"
        return "mixed"

    altitude_min_status = classify_lower(altitude_interval, cfg.min_altitude_m)
    if altitude_min_status is not None:
        checks["min_altitude_m"] = altitude_min_status

    altitude_max_status = classify_upper(altitude_interval, cfg.max_altitude_m)
    if altitude_max_status is not None:
        checks["max_altitude_m"] = altitude_max_status

    speed_status = classify_upper(speed_interval, cfg.max_speed_mps)
    if speed_status is not None:
        checks["max_speed_mps"] = speed_status

    q_status = classify_upper(q_interval, cfg.max_dynamic_pressure_pa)
    if q_status is not None:
        checks["max_dynamic_pressure_pa"] = q_status

    if not checks:
        return {"status": "not_evaluated", "checks": checks}

    if any(v == "outside" for v in checks.values()):
        status = "outside"
    elif any(v == "mixed" for v in checks.values()):
        status = "mixed"
    else:
        status = "inside"

    return {"status": status, "checks": checks}


def _zero_dx_box() -> List[Interval]:
    """
    build a zero derivative box

    input
    none

    output
    six zero intervals
    """
    return [Interval(0.0, 0.0) for _ in range(6)]


def make_invalid_interval_annotation(
    x_interval_old: List[Interval],
    sigma_interval_used: Interval,
    heat_shield: Optional[Any],
    failure_reason: str,
) -> IntervalAnnotationResult:
    """
    build a safe fallback annotation when interval propagation becomes invalid

    input
    old interval box sigma interval heat shield and reason string

    output
    invalid annotation result that preserves the last valid box
    """
    x_interval_new = copy_interval_box(x_interval_old)
    dx_interval = _zero_dx_box()
    widths_old = interval_component_widths(x_interval_old)
    widths_new = interval_component_widths(x_interval_new)
    dx_widths = interval_component_widths(dx_interval)
    altitude_interval = constants.intv_geometric_altitude(x_interval_old[0])

    return IntervalAnnotationResult(
        x_interval_old=copy_interval_box(x_interval_old),
        x_interval_new=x_interval_new,
        dx_interval=dx_interval,
        sigma_interval_used=copy_interval(sigma_interval_used),
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
        safety_status="interval_numerical_failure",
        safety_checks={"interval_failure": str(failure_reason)},
        layer_indices=[],
        interval_valid=False,
        interval_active=False,
        interval_failed=True,
        interval_failure_kind="numerical",
        interval_failure_reason=str(failure_reason),
        interval_numerical_failure=True,
        interval_heat_violation=False,
        last_valid_interval_state=copy_interval_box(x_interval_old),
    )


def _propagate_interval_annotation_core(
    x_interval_old: List[Interval],
    params: Dict[str, Any],
    sigma_interval_used: Interval,
    dt_s: float,
    supervisor_cfg: Optional[IntervalSupervisorConfig],
    heat_shield: Optional[Any],
    t_s: float,
) -> IntervalAnnotationResult:
    """
    propagate one interval step without recovery logic

    input
    old box parameters sigma step size config heat shield and time

    output
    one annotation result
    """
    atmospheric_hull = atmosphere_interval_hull_from_state_box(x_interval_old)
    rho_interval = atmospheric_hull["rho_kgm3"]
    altitude_interval = atmospheric_hull["altitude_m"]
    temperature_interval = atmospheric_hull["T_K"]
    pressure_interval = atmospheric_hull["p_Pa"]

    q_interval = interval_dynamic_pressure_bounds(rho_interval, x_interval_old[3])

    dx_interval = f_interval(
        t=float(t_s),
        X=x_interval_old,
        params=params,
        sigma_iv=sigma_interval_used,
    )

    x_interval_new = box_add(x_interval_old, box_scalar_mul(float(dt_s), dx_interval))

    heating_qdot_max_interval = None
    heating_qdot_mean_interval = None
    heating_Q_max_interval = None

    if (supervisor_cfg is not None and supervisor_cfg.include_heating) or (heat_shield is not None):
        heating_info = interval_heating_envelope(
            rho_interval=rho_interval,
            V_interval=x_interval_old[3],
            dt_s=float(dt_s),
            heat_shield=heat_shield,
        )
        heat_shield = heating_info["heat_shield"]
        heating_qdot_max_interval = heating_info["qdot_max"]
        heating_qdot_mean_interval = heating_info["qdot_mean"]
        heating_Q_max_interval = heating_info["Q_max"]

    safety_info = classify_interval_against_limits(
        x_interval_new=x_interval_new,
        q_interval=q_interval,
        cfg=supervisor_cfg,
    )

    heat_info = classify_heat_violation(
        heating_qdot_max_interval=heating_qdot_max_interval,
        heating_Q_max_interval=heating_Q_max_interval,
        cfg=supervisor_cfg,
    )

    safety_status = safety_info["status"]
    safety_checks = dict(safety_info["checks"])

    if bool(heat_info["rate_violation"]):
        safety_checks["heat_rate_limit"] = "outside"
    if bool(heat_info["load_violation"]):
        safety_checks["heat_load_limit"] = "outside"
    if bool(heat_info["heat_violation"]):
        safety_status = "heat_violation"

    last_valid_interval_state = copy_interval_box(x_interval_new)
    if bool(heat_info["heat_violation"]):
        last_valid_interval_state = copy_interval_box(x_interval_old)

    return IntervalAnnotationResult(
        x_interval_old=copy_interval_box(x_interval_old),
        x_interval_new=copy_interval_box(x_interval_new),
        dx_interval=copy_interval_box(dx_interval),
        sigma_interval_used=copy_interval(sigma_interval_used),
        rho_interval=copy_interval(rho_interval),
        q_interval=copy_interval(q_interval),
        altitude_interval=copy_interval(altitude_interval),
        temperature_interval=copy_interval(temperature_interval) if temperature_interval is not None else None,
        pressure_interval=copy_interval(pressure_interval) if pressure_interval is not None else None,
        state_widths_old=interval_component_widths(x_interval_old),
        state_widths_new=interval_component_widths(x_interval_new),
        dx_widths=interval_component_widths(dx_interval),
        heating_qdot_max_interval=copy_interval(heating_qdot_max_interval) if heating_qdot_max_interval is not None else None,
        heating_qdot_mean_interval=copy_interval(heating_qdot_mean_interval) if heating_qdot_mean_interval is not None else None,
        heating_Q_max_interval=copy_interval(heating_Q_max_interval) if heating_Q_max_interval is not None else None,
        heat_shield=heat_shield,
        safety_status=safety_status,
        safety_checks=safety_checks,
        layer_indices=list(atmospheric_hull["layer_indices"]),
        interval_valid=True,
        interval_active=not bool(heat_info["heat_violation"]),
        interval_failed=bool(heat_info["heat_violation"]),
        interval_failure_kind="heat" if bool(heat_info["heat_violation"]) else "none",
        interval_failure_reason=str(heat_info["reason"]),
        interval_numerical_failure=False,
        interval_heat_violation=bool(heat_info["heat_violation"]),
        last_valid_interval_state=last_valid_interval_state,
    )


def merge_child_interval_annotations(
    child_results: List[IntervalAnnotationResult],
    x_interval_old: List[Interval],
    sigma_interval_used: Interval,
    supervisor_cfg: Optional[IntervalSupervisorConfig],
    recentered_this_step: bool,
    recenter_reason: str,
    split_depth_used: int,
    split_reason: str,
) -> IntervalAnnotationResult:
    """
    hull child annotations back into one parent result

    input
    child results old box sigma config recenter info split depth and split reason

    output
    one merged annotation result
    """
    x_interval_new = hull_interval_boxes([child.x_interval_new for child in child_results])
    dx_interval = hull_interval_boxes([child.dx_interval for child in child_results])
    rho_interval = hull_intervals([child.rho_interval for child in child_results])
    q_interval = hull_intervals([child.q_interval for child in child_results])
    altitude_interval = hull_intervals([child.altitude_interval for child in child_results])

    temperature_candidates = [child.temperature_interval for child in child_results if child.temperature_interval is not None]
    pressure_candidates = [child.pressure_interval for child in child_results if child.pressure_interval is not None]
    qdot_candidates = [child.heating_qdot_max_interval for child in child_results if child.heating_qdot_max_interval is not None]
    qmean_candidates = [child.heating_qdot_mean_interval for child in child_results if child.heating_qdot_mean_interval is not None]
    Q_candidates = [child.heating_Q_max_interval for child in child_results if child.heating_Q_max_interval is not None]

    temperature_interval = hull_intervals(temperature_candidates) if temperature_candidates else None
    pressure_interval = hull_intervals(pressure_candidates) if pressure_candidates else None
    heating_qdot_max_interval = hull_intervals(qdot_candidates) if qdot_candidates else None
    heating_qdot_mean_interval = hull_intervals(qmean_candidates) if qmean_candidates else None
    heating_Q_max_interval = hull_intervals(Q_candidates) if Q_candidates else None
    heat_shield = hull_heat_shields([child.heat_shield for child in child_results if child.heat_shield is not None])

    safety_info = classify_interval_against_limits(
        x_interval_new=x_interval_new,
        q_interval=q_interval,
        cfg=supervisor_cfg,
    )

    heat_info = classify_heat_violation(
        heating_qdot_max_interval=heating_qdot_max_interval,
        heating_Q_max_interval=heating_Q_max_interval,
        cfg=supervisor_cfg,
    )

    safety_status = safety_info["status"]
    safety_checks = dict(safety_info["checks"])
    if bool(heat_info["rate_violation"]):
        safety_checks["heat_rate_limit"] = "outside"
    if bool(heat_info["load_violation"]):
        safety_checks["heat_load_limit"] = "outside"
    if bool(heat_info["heat_violation"]):
        safety_status = "heat_violation"

    failure_reason_parts = []
    for child in child_results:
        if child.interval_failure_reason and child.interval_failure_reason not in failure_reason_parts:
            failure_reason_parts.append(child.interval_failure_reason)

    last_valid_interval_state = copy_interval_box(x_interval_new)
    if bool(heat_info["heat_violation"]):
        last_valid_interval_state = copy_interval_box(x_interval_old)

    layer_indices: List[int] = []
    for child in child_results:
        for layer_index in child.layer_indices:
            if layer_index not in layer_indices:
                layer_indices.append(layer_index)

    return IntervalAnnotationResult(
        x_interval_old=copy_interval_box(x_interval_old),
        x_interval_new=copy_interval_box(x_interval_new),
        dx_interval=copy_interval_box(dx_interval),
        sigma_interval_used=copy_interval(sigma_interval_used),
        rho_interval=copy_interval(rho_interval),
        q_interval=copy_interval(q_interval),
        altitude_interval=copy_interval(altitude_interval),
        temperature_interval=copy_interval(temperature_interval) if temperature_interval is not None else None,
        pressure_interval=copy_interval(pressure_interval) if pressure_interval is not None else None,
        state_widths_old=interval_component_widths(x_interval_old),
        state_widths_new=interval_component_widths(x_interval_new),
        dx_widths=interval_component_widths(dx_interval),
        heating_qdot_max_interval=copy_interval(heating_qdot_max_interval) if heating_qdot_max_interval is not None else None,
        heating_qdot_mean_interval=copy_interval(heating_qdot_mean_interval) if heating_qdot_mean_interval is not None else None,
        heating_Q_max_interval=copy_interval(heating_Q_max_interval) if heating_Q_max_interval is not None else None,
        heat_shield=heat_shield,
        safety_status=safety_status,
        safety_checks=safety_checks,
        layer_indices=layer_indices,
        interval_valid=True,
        interval_active=not bool(heat_info["heat_violation"]),
        interval_failed=bool(heat_info["heat_violation"]),
        interval_failure_kind="heat" if bool(heat_info["heat_violation"]) else "none",
        interval_failure_reason=",".join(failure_reason_parts),
        interval_numerical_failure=False,
        interval_heat_violation=bool(heat_info["heat_violation"]),
        recentered_this_step=bool(recentered_this_step),
        recenter_reason=str(recenter_reason),
        split_this_step=True,
        split_reason=str(split_reason),
        split_depth_used=int(split_depth_used),
        split_child_count=len(child_results),
        last_valid_interval_state=last_valid_interval_state,
    )


def split_and_propagate_interval_annotation(
    x_interval_old: List[Interval],
    params: Dict[str, Any],
    sigma_interval_used: Interval,
    dt_s: float,
    supervisor_cfg: Optional[IntervalSupervisorConfig],
    heat_shield: Optional[Any],
    t_s: float,
    split_depth: int,
    recentered_this_step: bool,
    recenter_reason: str,
    split_reason: str,
) -> IntervalAnnotationResult:
    """
    recursively split and propagate an unhealthy interval box

    input
    old box parameters sigma step size config heat shield time and split state

    output
    one annotation result after recursive splitting
    """
    if supervisor_cfg is None or not bool(supervisor_cfg.interval_box_split_enabled):
        return make_invalid_interval_annotation(
            x_interval_old=x_interval_old,
            sigma_interval_used=sigma_interval_used,
            heat_shield=heat_shield,
            failure_reason="adaptive_box_split_disabled",
        )

    if int(split_depth) >= int(supervisor_cfg.interval_box_split_max_depth):
        return make_invalid_interval_annotation(
            x_interval_old=x_interval_old,
            sigma_interval_used=sigma_interval_used,
            heat_shield=heat_shield,
            failure_reason="adaptive_box_split_depth_exhausted",
        )

    split_index = choose_split_dimension(x_interval_old, supervisor_cfg)
    if float(x_interval_old[split_index].width()) <= 1.0e-12:
        return make_invalid_interval_annotation(
            x_interval_old=x_interval_old,
            sigma_interval_used=sigma_interval_used,
            heat_shield=heat_shield,
            failure_reason="adaptive_box_split_zero_width_dimension",
        )

    left_box, right_box = box_split(copy_interval_box(x_interval_old), idx=int(split_index))
    child_results: List[IntervalAnnotationResult] = []

    for child_box in [left_box, right_box]:
        child_heat_shield = clone_interval_heat_shield(heat_shield)

        try:
            child_result = _propagate_interval_annotation_core(
                x_interval_old=child_box,
                params=params,
                sigma_interval_used=sigma_interval_used,
                dt_s=float(dt_s),
                supervisor_cfg=supervisor_cfg,
                heat_shield=child_heat_shield,
                t_s=float(t_s),
            )

            if (
                child_result.interval_valid
                and not child_result.interval_heat_violation
                and interval_box_needs_split(child_result.x_interval_new, supervisor_cfg)
            ):
                child_result = split_and_propagate_interval_annotation(
                    x_interval_old=child_box,
                    params=params,
                    sigma_interval_used=sigma_interval_used,
                    dt_s=float(dt_s),
                    supervisor_cfg=supervisor_cfg,
                    heat_shield=child_heat_shield,
                    t_s=float(t_s),
                    split_depth=int(split_depth) + 1,
                    recentered_this_step=bool(recentered_this_step),
                    recenter_reason=str(recenter_reason),
                    split_reason=split_reason_from_box(child_result.x_interval_new, supervisor_cfg),
                )

        except Exception as exc:
            if int(split_depth) + 1 < int(supervisor_cfg.interval_box_split_max_depth):
                child_result = split_and_propagate_interval_annotation(
                    x_interval_old=child_box,
                    params=params,
                    sigma_interval_used=sigma_interval_used,
                    dt_s=float(dt_s),
                    supervisor_cfg=supervisor_cfg,
                    heat_shield=child_heat_shield,
                    t_s=float(t_s),
                    split_depth=int(split_depth) + 1,
                    recentered_this_step=bool(recentered_this_step),
                    recenter_reason=str(recenter_reason),
                    split_reason=str(split_reason),
                )
            else:
                child_result = make_invalid_interval_annotation(
                    x_interval_old=child_box,
                    sigma_interval_used=sigma_interval_used,
                    heat_shield=child_heat_shield,
                    failure_reason=f"adaptive_box_split_child_invalid: {exc}",
                )

        if not child_result.interval_valid:
            return child_result

        child_results.append(child_result)

    return merge_child_interval_annotations(
        child_results=child_results,
        x_interval_old=x_interval_old,
        sigma_interval_used=sigma_interval_used,
        supervisor_cfg=supervisor_cfg,
        recentered_this_step=bool(recentered_this_step),
        recenter_reason=str(recenter_reason),
        split_depth_used=int(split_depth) + 1,
        split_reason=str(split_reason),
    )


def build_interval_annotation(
    x_interval_old: List[Interval],
    params: Dict[str, Any],
    sigma_interval_used: Interval,
    dt_s: Optional[float] = None,
    supervisor_cfg: Optional[IntervalSupervisorConfig] = None,
    heat_shield: Optional[Any] = None,
    t_s: float = 0.0,
    x_nominal_center: Optional[List[float]] = None,
) -> IntervalAnnotationResult:
    """
    propagate one interval step with recenter and split recovery logic

    input
    old box parameters sigma optional step size config heat shield time and nominal center

    output
    one annotation result
    """
    if dt_s is None:
        dt_s = float(getattr(constants, "ENTRY_DT_S", 0.25))

    if len(x_interval_old) != 6:
        raise ValueError("x_interval_old must have 6 interval components")

    working_box = copy_interval_box(x_interval_old)
    recentered_this_step = False
    recenter_reason = ""

    # recenter now carries an explicit reason
    # this keeps width based health events separate from fixed cadence maintenance
    recenter_needed, recenter_reason = should_recenter_interval_box(
        x_interval_old=working_box,
        x_nominal_center=x_nominal_center,
        cfg=supervisor_cfg,
        t_s=float(t_s),
        dt_s=float(dt_s),
    )

    if recenter_needed:
        working_box = recenter_interval_box_around_nominal(
            x_nominal_center=list(x_nominal_center),
            cfg=supervisor_cfg,
        )
        recentered_this_step = True

    try:
        if interval_box_needs_split(working_box, supervisor_cfg):
            result = split_and_propagate_interval_annotation(
                x_interval_old=working_box,
                params=params,
                sigma_interval_used=sigma_interval_used,
                dt_s=float(dt_s),
                supervisor_cfg=supervisor_cfg,
                heat_shield=heat_shield,
                t_s=float(t_s),
                split_depth=0,
                recentered_this_step=bool(recentered_this_step),
                recenter_reason=str(recenter_reason),
                split_reason=split_reason_from_box(working_box, supervisor_cfg),
            )
        else:
            result = _propagate_interval_annotation_core(
                x_interval_old=working_box,
                params=params,
                sigma_interval_used=sigma_interval_used,
                dt_s=float(dt_s),
                supervisor_cfg=supervisor_cfg,
                heat_shield=heat_shield,
                t_s=float(t_s),
            )

            if (
                result.interval_valid
                and not result.interval_heat_violation
                and interval_box_needs_split(result.x_interval_new, supervisor_cfg)
            ):
                result = split_and_propagate_interval_annotation(
                    x_interval_old=working_box,
                    params=params,
                    sigma_interval_used=sigma_interval_used,
                    dt_s=float(dt_s),
                    supervisor_cfg=supervisor_cfg,
                    heat_shield=heat_shield,
                    t_s=float(t_s),
                    split_depth=0,
                    recentered_this_step=bool(recentered_this_step),
                    recenter_reason=str(recenter_reason),
                    split_reason=split_reason_from_box(result.x_interval_new, supervisor_cfg),
                )

        result.recentered_this_step = bool(result.recentered_this_step or recentered_this_step)
        if recenter_reason and not result.recenter_reason:
            result.recenter_reason = str(recenter_reason)
        return result

    except Exception as exc:
        return make_invalid_interval_annotation(
            x_interval_old=working_box,
            sigma_interval_used=sigma_interval_used,
            heat_shield=heat_shield,
            failure_reason=f"live_interval_annotation_invalid: {exc}",
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
    annotate one nominal step with an interval tube

    input
    nominal state parameters actual sigma old interval optional config heat shield and time

    output
    one interval annotation result
    """
    if x_interval_old is None:
        if supervisor_cfg is None:
            x_interval_old = nominal_state_to_interval_box(x_nominal_old)
        else:
            x_interval_old = inflate_nominal_state_to_interval_box(x_nominal_old, supervisor_cfg)

    sigma_interval_used = promote(float(sigma_actual_after_rad))

    try:
        return build_interval_annotation(
            x_interval_old=x_interval_old,
            params=params,
            sigma_interval_used=sigma_interval_used,
            dt_s=dt_s,
            supervisor_cfg=supervisor_cfg,
            heat_shield=heat_shield,
            t_s=t_s,
            x_nominal_center=list(x_nominal_old),
        )
    except Exception as exc:
        return make_invalid_interval_annotation(
            x_interval_old=x_interval_old,
            sigma_interval_used=sigma_interval_used,
            heat_shield=heat_shield,
            failure_reason=f"live_interval_annotation_invalid: {exc}",
        )


def aero_coefficients_from_speed_altitude(V_mps: float, altitude_m: float, params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Compute aerodynamic coefficients using a reduced-order speed-based schedule.
    
    This is not the official Orion aerodatabase. It is a scheduled Orion-like
    coefficient model where CD and L/D are scheduled, and CL is computed as CD * L/D.
    It preserves the old constant model as a fallback.
    """
    model_name = params.get("aero_model", "constant")
    V_kms = V_mps / 1000.0
    V_kfps = V_mps * 0.00328084  # proxy for FMV in the Orion CM hypersonic database

    # Properly blended FMV per Bibb et al. 2010 Eq. 1:
    #   FMV = Mach          if U <=  8.8 kfps  (~2680 m/s)
    #   FMV = V_kfps        if U >=  9.8 kfps  (~2987 m/s)
    #   linear blend between
    # Mach needs the local atmospheric temperature. We accept it via params
    # (cheap: the main loop populates params['T_K'] once per step from the
    # atmosphere query it already does); if absent we query atmosphere here
    # so the function still returns the right answer in isolation.
    def _fmv_blended(V_mps_, V_kfps_, T_K_):
        if V_kfps_ >= 9.8:
            return V_kfps_
        if T_K_ is None or not (T_K_ == T_K_) or T_K_ <= 0.0:
            return V_kfps_  # no temperature available -> degrade to old behaviour
        mach_ = V_mps_ / math.sqrt(1.4 * 287.053 * T_K_)
        if V_kfps_ <= 8.8:
            return mach_
        # 8.8 kfps < V_kfps < 9.8 kfps: linear blend
        alpha = (V_kfps_ - 8.8) / (9.8 - 8.8)
        return (1.0 - alpha) * mach_ + alpha * V_kfps_

    _T_K_in_params = params.get("T_K", None)
    if _T_K_in_params is None and model_name == "orion_cm_trim":
        # Pay the atmosphere query so FMV is correct in the subsonic regime
        try:
            import AtmosphereModel as _atm_mod
            _, _T_K_in_params = _atm_mod.get_atmosphere_state(float(altitude_m))
        except Exception:
            _T_K_in_params = None
    FMV_blended = _fmv_blended(V_mps, V_kfps, _T_K_in_params)

    if model_name == "orion_cm_trim":
        # Trim CD and L/D scheduled vs the blended velocity parameter FMV from
        # Bibb, Walker, Brauckmann, Robinson, "Development of the Orion Crew Module
        # Static Aerodynamic Database, Part I: Hypersonic" (NASA Langley / JSC).
        # Paper covers 8 <= FMV <= 40 at trim (Cm_cg = 0); trim alpha varies from
        # ~167 deg at low FMV to ~157 deg at high FMV, giving L/D ~0.27 to 0.40.
        # Below FMV = 8 (i.e. below ~Mach 8, where blended FMV is Mach number)
        # we use the ORION-SPECIFIC subsonic/supersonic static aero database:
        #   Bibb et al., "Development of the Orion Crew Module Static
        #   Aerodynamic Database, Part II: Subsonic/Supersonic," AIAA 2011-3507
        #   / NASA NTRS 20110013645. Covers M = 0.2-4.0 at trim (alpha ~145-175
        #   deg, heat-shield forward).
        # Key Orion-specific facts taken directly from that paper:
        #   - SUBSONIC trim CD is ~0.7-0.9 (frontal-area reference) -- explicitly
        #     LOWER than the historical Apollo flight data, so the previous
        #     Apollo surrogate over-predicted subsonic drag.
        #   - CD rises through the transonic drag rise toward the hypersonic
        #     value ~1.55 (Bibb Part I) above M ~ 5.
        # L/D values below taper smoothly from the Part I hypersonic trim
        # (~0.40 at M=8) down to ~0.10 at very low subsonic Mach where trim
        # lift authority collapses; the paper notes L/D is secondary to drag
        # in this regime.
        # Format: (FMV, CD_trim, LD_trim)
        default_trim_schedule = [
            (0.0,  0.83, 0.10),  # very low subsonic — Orion CD ~0.8, minimal lift
            (0.5,  0.85, 0.15),  # subsonic (Orion Part II: CD 0.7-0.9)
            (0.7,  0.88, 0.20),  # high subsonic (matches M=0.7 fig, CD<Apollo)
            (0.9,  1.02, 0.23),  # transonic onset — drag rise begins
            (1.1,  1.28, 0.25),  # transonic drag rise
            (1.5,  1.38, 0.28),  # low supersonic
            (2.0,  1.43, 0.30),  # supersonic
            (3.0,  1.47, 0.33),  # supersonic
            (5.0,  1.50, 0.36),  # high supersonic, approaching hypersonic
            (8.0,  1.52, 0.40),  # hypersonic low FMV (Bibb Part I) -- blend point
            (10.0, 1.50, 0.39),  # ---- below: Bibb 2010 Part I hypersonic ----
            (15.0, 1.52, 0.33),
            (20.0, 1.54, 0.30),
            (25.0, 1.55, 0.27),
            (30.0, 1.55, 0.27),
            (35.0, 1.55, 0.28),
            (40.0, 1.55, 0.30),
        ]
        schedule = params.get("orion_trim_schedule", default_trim_schedule)
        schedule = sorted(schedule, key=lambda x: x[0])

        # Use the properly-blended FMV from above (Bibb 2010 Eq. 1) so
        # the subsonic/transonic regime is keyed on Mach instead of the
        # numerically-different V_kfps proxy.
        FMV = FMV_blended

        if FMV <= schedule[0][0]:
            CD = schedule[0][1]
            LD = schedule[0][2]
            region_string = "below_table"
        elif FMV >= schedule[-1][0]:
            CD = schedule[-1][1]
            LD = schedule[-1][2]
            region_string = "above_table"
        else:
            CD = schedule[0][1]
            LD = schedule[0][2]
            region_string = "interpolated"
            for i in range(len(schedule) - 1):
                f0, cd0, ld0 = schedule[i]
                f1, cd1, ld1 = schedule[i + 1]
                if f0 <= FMV <= f1:
                    alpha = (FMV - f0) / (f1 - f0) if f1 > f0 else 0.0
                    CD = cd0 + alpha * (cd1 - cd0)
                    LD = ld0 + alpha * (ld1 - ld0)
                    break

        CL = CD * LD

    elif model_name == "scheduled_orion_like":
        # Default simple speed-based schedule with smooth linear interpolation
        # V_kms, CD, LD
        default_schedule = [
            (0.0,  1.50, 0.10),
            (1.0,  1.45, 0.15),
            (2.0,  1.35, 0.22),
            (5.0,  1.25, 0.27),
            (7.5,  1.20, 0.28),
            (11.0, 1.15, 0.24),
        ]
        schedule = params.get("aero_schedule", default_schedule)
        
        # Sort schedule by V_kms just in case
        schedule = sorted(schedule, key=lambda x: x[0])
        
        if V_kms <= schedule[0][0]:
            CD = schedule[0][1]
            LD = schedule[0][2]
        elif V_kms >= schedule[-1][0]:
            CD = schedule[-1][1]
            LD = schedule[-1][2]
        else:
            CD = schedule[0][1]
            LD = schedule[0][2]
            # Linear interpolation
            for i in range(len(schedule) - 1):
                v0, cd0, ld0 = schedule[i]
                v1, cd1, ld1 = schedule[i+1]
                if v0 <= V_kms <= v1:
                    alpha = (V_kms - v0) / (v1 - v0) if v1 > v0 else 0.0
                    CD = cd0 + alpha * (cd1 - cd0)
                    LD = ld0 + alpha * (ld1 - ld0)
                    break
        
        CL = CD * LD
        region_string = "interpolated"
        if V_kms <= schedule[0][0]:
            region_string = "low_speed_clamp"
        elif V_kms >= schedule[-1][0]:
            region_string = "high_speed_clamp"

    else:
        # Fallback to constant model
        CD = _get_param(params, "CD_const", "CD")
        CL = _get_param(params, "CL_const", "CL")
        LD = CL / CD if CD != 0.0 else 0.0
        region_string = "constant"

    # ------------------------------------------------------------------
    # CPAS (parachute) modifiers. The main loop pushes the current chute
    # drag CD*A and lift scale into params via the cpas_* keys below.
    # We fold them into the effective aero coefficients here so the same
    # CD/CL pipeline that feeds the EOMs and the predictor-corrector
    # rollouts automatically sees the chute contribution.
    # ------------------------------------------------------------------
    cpas_drag_cdA_extra_m2 = float(params.get("cpas_drag_cdA_extra_m2", 0.0))
    cpas_lift_scale = float(params.get("cpas_lift_scale", 1.0))
    ref_area_m2 = float(params.get("ref_area_m2", 0.0))
    if cpas_drag_cdA_extra_m2 > 0.0 and ref_area_m2 > 0.0:
        CD = CD + cpas_drag_cdA_extra_m2 / ref_area_m2
    if cpas_lift_scale != 1.0:
        CL = CL * cpas_lift_scale
    LD = (CL / CD) if CD != 0.0 else 0.0

    return {
        "aero_model": model_name,
        "CD": float(CD),
        "CL": float(CL),
        "LD": float(LD),
        "V_kms": float(V_kms),
        "V_kfps": float(V_kfps),
        "FMV": float(FMV_blended),
        "altitude_m": float(altitude_m),
        "schedule_region": region_string,
        "cpas_drag_cdA_extra_m2": float(cpas_drag_cdA_extra_m2),
        "cpas_lift_scale": float(cpas_lift_scale),
    }


def nominal_aero_forces_from_state(
    x: List[float],
    sigma_rad: float,
    params: Dict[str, Any],
) -> Dict[str, float]:
    """
    compute nominal aerodynamic force components

    this resolves drag opposite velocity and rotates lift with bank
    into radial north and east components

    input
    nominal state bank angle and parameter dictionary

    output
    aerodynamic force dictionary
    """
    if len(x) != 6:
        raise ValueError("x must have 6 state components")

    r_m, phi_rad, lam_rad, V_mps, gamma_rad, chi_rad = x

    altitude_m = max(0.0, constants.altitude(r_m))

    atm = AtmosphereModel.US_Standard_ATM(
        altitude_m,
        lat_deg=math.degrees(phi_rad),
        lon_deg=math.degrees(lam_rad),
        date=params.get("entry_datetime_utc", None),
        f107=params.get("f107", 150.0),
        f107a=params.get("f107a", 150.0),
        ap=params.get("ap", 7.0),
    )
    rho = atm["rho_kgm3"] if atm["rho_kgm3"] is not None else 0.0
    q_pa = 0.5 * rho * V_mps * V_mps

    surface_area = _get_param(params, "ref_area_m2", "S")
    coeffs = aero_coefficients_from_speed_altitude(
        V_mps=V_mps,
        altitude_m=altitude_m,
        params=params,
    )
    CD = coeffs["CD"]
    CL = coeffs["CL"]

    drag_mag_N = q_pa * surface_area * CD
    lift_mag_N = q_pa * surface_area * CL

    sin_gamma = math.sin(gamma_rad)
    cos_gamma = math.cos(gamma_rad)
    sin_chi = math.sin(chi_rad)
    cos_chi = math.cos(chi_rad)
    sin_sigma = math.sin(sigma_rad)
    cos_sigma = math.cos(sigma_rad)

    D_r_N = -drag_mag_N * sin_gamma
    D_north_N = -drag_mag_N * cos_gamma * sin_chi
    D_east_N = -drag_mag_N * cos_gamma * cos_chi

    lift_vertical_plane_N = lift_mag_N * cos_sigma
    lift_lateral_left_N = -lift_mag_N * sin_sigma

    L_r_N = lift_vertical_plane_N * cos_gamma
    L_north_N = -lift_mag_N * (cos_sigma * sin_gamma * sin_chi + sin_sigma * cos_chi)
    L_east_N = lift_mag_N * (-cos_sigma * sin_gamma * cos_chi + sin_sigma * sin_chi)

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
        "CD": float(CD),
        "CL": float(CL),
        "LD": float(coeffs["LD"]),
        "FMV": float(coeffs.get("FMV", coeffs.get("V_kfps", 0.0))),
        "aero_model": coeffs["aero_model"],
        "aero_schedule_region": coeffs["schedule_region"],
        "V_kms": coeffs["V_kms"],
    }


def nominal_eom_step(
    x: List[float],
    sigma_rad: float,
    params: Dict[str, Any],
    dt_s: float,
) -> Dict[str, Any]:
    """
    advance the corrected nominal translational equations by one euler step

    equations
    drag acts only in v dot
    vertical plane lift enters gamma dot
    lateral lift enters chi dot

    input
    nominal state bank angle parameters and step size

    output
    next state derivative and aero summary
    """
    if len(x) != 6:
        raise ValueError("x must have 6 state components")
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive")

    r_m, phi_rad, lam_rad, V_mps, gamma_rad, chi_rad = x

    aero = nominal_aero_forces_from_state(x=x, sigma_rad=sigma_rad, params=params)
    mass_kg = _get_param(params, "mass_kg", "m")
    g_mps2 = constants.gravity(r_m)

    cos_phi = math.cos(phi_rad)
    sin_phi = math.sin(phi_rad)
    cos_gamma = math.cos(gamma_rad)
    sin_gamma = math.sin(gamma_rad)
    sin_chi = math.sin(chi_rad)
    cos_chi = math.cos(chi_rad)

    eps = 1.0e-9
    safe_r = max(r_m, constants.RADIUS_EARTH + 1.0)
    safe_V = max(abs(V_mps), eps)

    if abs(cos_phi) < eps:
        cos_phi = eps if cos_phi >= 0.0 else -eps

    if abs(cos_gamma) < eps:
        cos_gamma = eps if cos_gamma >= 0.0 else -eps

    tan_phi = sin_phi / cos_phi

    r_dot = V_mps * sin_gamma
    phi_dot = (V_mps * cos_gamma * sin_chi) / safe_r
    lam_dot = (V_mps * cos_gamma * cos_chi) / (safe_r * cos_phi)
    V_dot = -aero["drag_mag_N"] / mass_kg - g_mps2 * sin_gamma

    gamma_dot = (
        aero["lift_vertical_plane_N"] / (mass_kg * safe_V)
        + (safe_V / safe_r - g_mps2 / safe_V) * cos_gamma
    )

    chi_dot = (
        aero["lift_lateral_left_N"] / (mass_kg * safe_V * cos_gamma)
        - (safe_V / safe_r) * cos_gamma * cos_chi * tan_phi
    )

    x_next = [
        r_m + r_dot * dt_s,
        phi_rad + phi_dot * dt_s,
        lam_rad + lam_dot * dt_s,
        V_mps + V_dot * dt_s,
        gamma_rad + gamma_dot * dt_s,
        chi_rad + chi_dot * dt_s,
    ]

    return {
        "x_next": x_next,
        "x_dot": [r_dot, phi_dot, lam_dot, V_dot, gamma_dot, chi_dot],
        "aero": aero,
        "gravity_mps2": g_mps2,
    }


def intv_aero_forces(
    r: Interval,
    V: Interval,
    gamma: Interval,
    chi: Interval,
    sigma: Interval,
    params: Dict[str, Any],
) -> Dict[int, Dict[str, Interval]]:
    """
    compute interval aerodynamic forces over all touched atmosphere layers

    input
    interval state pieces sigma interval and parameters

    output
    dictionary of force dictionaries by layer
    """
    z_geometric_altitude = constants.intv_geometric_altitude(r)
    atm_intv = AtmosphereModel.intv_US_Standard_ATM(z_geometric_altitude)

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

    aero_forces_by_layer: Dict[int, Dict[str, Interval]] = {}
    surface_area = _get_param(params, "ref_area_m2", "S")

    # Evaluate scheduled coefficient model using interval endpoints
    # This is a reduced-order approximation to compute a conservative interval.
    z_mid = float(z_geometric_altitude.lo + z_geometric_altitude.hi) / 2.0
    v_lo = float(V.lo)
    v_hi = float(V.hi)
    
    coeffs_lo = aero_coefficients_from_speed_altitude(V_mps=v_lo, altitude_m=z_mid, params=params)
    coeffs_hi = aero_coefficients_from_speed_altitude(V_mps=v_hi, altitude_m=z_mid, params=params)
    
    CD_iv = Interval(min(coeffs_lo["CD"], coeffs_hi["CD"]), max(coeffs_lo["CD"], coeffs_hi["CD"]))
    CL_iv = Interval(min(coeffs_lo["CL"], coeffs_hi["CL"]), max(coeffs_lo["CL"], coeffs_hi["CL"]))

    sin_gamma = gamma.sin()
    cos_gamma = gamma.cos()
    sin_chi = chi.sin()
    cos_chi = chi.cos()
    sin_sigma = sigma.sin()
    cos_sigma = sigma.cos()

    for k in atm_intv.keys():
        rho = atm_intv[k]["rho_kgm3"] if atm_intv[k]["rho_kgm3"] is not None else promote(0.0)
        q = promote(0.5) * rho * V.pow_int(2)

        drag_mag = q * surface_area * CD_iv
        lift_mag = q * surface_area * CL_iv

        D_r = -drag_mag * sin_gamma
        D_north = -drag_mag * cos_gamma * sin_chi
        D_east = -drag_mag * cos_gamma * cos_chi

        lift_vertical_plane = lift_mag * cos_sigma
        lift_lateral_left = -lift_mag * sin_sigma

        L_r = lift_vertical_plane * cos_gamma
        L_north = -lift_mag * (cos_sigma * sin_gamma * sin_chi + sin_sigma * cos_chi)
        L_east = lift_mag * (-cos_sigma * sin_gamma * cos_chi + sin_sigma * sin_chi)

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


def intv_eom_3d(
    t: float,
    X: List[Interval],
    params: Dict[str, Any],
    sigma: Interval,
) -> Dict[str, Interval]:
    """
    evaluate the interval valued three dof translational equations

    input
    time interval state parameters and bank interval

    output
    named derivative dictionary
    """
    del t

    if len(X) != 6:
        raise ValueError("X must have 6 interval components")

    r = X[0]
    phi = X[1]
    lam = X[2]
    V = X[3]
    gamma = X[4]
    chi = X[5]
    del lam

    mass_kg = _get_param(params, "mass_kg", "m")
    g = constants.intv_gravity(r)

    aero = intv_aero_forces(
        r=r,
        V=V,
        gamma=gamma,
        chi=chi,
        sigma=sigma,
        params=params,
    )

    drag_mag = hull_intervals([aero[layer]["drag_mag"] for layer in aero])
    lift_vertical_plane = hull_intervals([aero[layer]["lift_vertical_plane"] for layer in aero])
    lift_lateral_left = hull_intervals([aero[layer]["lift_lateral_left"] for layer in aero])

    cos_phi = phi.cos()
    tan_phi = phi.sin() / cos_phi
    cos_gamma = gamma.cos()
    sin_gamma = gamma.sin()
    sin_chi = chi.sin()
    cos_chi = chi.cos()

    r_dot = V * sin_gamma
    phi_dot = (V * cos_gamma * sin_chi) / r
    lam_dot = (V * cos_gamma * cos_chi) / (r * cos_phi)
    V_dot = -drag_mag / mass_kg - g * sin_gamma
    gamma_dot = (lift_vertical_plane / (mass_kg * V)) + (V / r - g / V) * cos_gamma
    chi_dot = (lift_lateral_left / (mass_kg * V * cos_gamma)) - (V / r) * cos_gamma * cos_chi * tan_phi

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
    convert the named interval derivative dictionary into state order

    input
    time interval state parameters and bank interval

    output
    ordered derivative box
    """
    dx = intv_eom_3d(t=t, X=X, params=params, sigma=sigma_iv)
    return [
        dx["r_dot"],
        dx["phi_dot"],
        dx["lam_dot"],
        dx["V_dot"],
        dx["gamma_dot"],
        dx["chi_dot"],
    ]


def _positive_violation_amount(interval_value: Interval, limit_value: Optional[float]) -> float:
    """
    return normalized positive violation amount from the interval upper bound

    input
    interval value and optional limit

    output
    nonnegative scalar violation amount
    """
    if limit_value is None:
        return 0.0

    denom = max(abs(float(limit_value)), 1.0)
    return max(0.0, (float(interval_value.hi) - float(limit_value)) / denom)


def _hull_or_default(a: Optional[Interval], b: Optional[Interval]) -> Interval:
    """
    merge two optional intervals into one interval

    input
    two optional intervals

    output
    one interval
    """
    if a is None and b is None:
        return Interval(0.0, 0.0)
    if a is None:
        return Interval(float(b.lo), float(b.hi))
    if b is None:
        return Interval(float(a.lo), float(a.hi))
    return a.hull(b)


def rollout_nominal_to_terminal(
    x_nominal_old: List[float],
    sigma_cmd_rad: float,
    params: Dict[str, Any],
    dt_s: float,
    terminal_altitude_m: float,
    max_steps: int,
) -> Dict[str, Any]:
    """
    Propagate ONLY the nominal translational state forward under a constant
    bank angle until the predicted altitude drops below terminal_altitude_m
    (typically drogue-deploy altitude, where the entry phase ends and the
    bank angle no longer matters), or until max_steps is reached.

    This is the long-horizon predictor used by the continuous bank-magnitude
    solver (Lu 2014 Eq. 18-24). Unlike rollout_predictor_corrector_candidate,
    it does NOT propagate the interval box or the heat shield — it is a cheap
    nominal-only forward integration whose sole job is to estimate where the
    capsule lands as a function of bank angle. Heat feasibility is still
    screened separately on the short horizon by the interval rollout.

    A coarse dt_s (e.g. 1-2 s) keeps the cost low; the predictor only needs
    the terminal POSITION, not a flight-accurate trajectory.

    Returns a dict:
        x_final              -- final 6-state vector
        reached_terminal     -- True if alt fell below terminal_altitude_m
        steps_used           -- number of Euler steps taken
        terminated_reason    -- "terminal_alt" | "max_steps" | "singular"
    """
    if len(x_nominal_old) != 6:
        raise ValueError("x_nominal_old must have 6 components")
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive")
    if max_steps <= 0:
        raise ValueError("max_steps must be positive")

    # Use a CPAS-free copy of params: during the entry/guidance phase the
    # chutes are not deployed, and we never want a frozen chute-drag term
    # leaking into a long forward prediction.
    pred_params = dict(params)
    pred_params["cpas_drag_cdA_extra_m2"] = 0.0
    pred_params["cpas_lift_scale"] = 1.0

    x_nominal = [float(v) for v in x_nominal_old]
    r_earth = float(constants.RADIUS_EARTH)
    terminated_reason = "max_steps"
    steps_used = 0

    for step_idx in range(int(max_steps)):
        steps_used = step_idx + 1

        # Stop if the flight-path geometry is about to go singular (cos gamma
        # near zero makes the longitude/heading channels blow up). Treat the
        # current state as the terminal prediction.
        if abs(math.cos(x_nominal[4])) < 0.05:
            terminated_reason = "singular"
            break

        nominal_step = nominal_eom_step(
            x=x_nominal,
            sigma_rad=float(sigma_cmd_rad),
            params=pred_params,
            dt_s=float(dt_s),
        )
        x_nominal = [float(v) for v in nominal_step["x_next"]]

        # Guard against a NaN / non-finite blow-up in the coarse integration.
        if not all(math.isfinite(v) for v in x_nominal):
            terminated_reason = "singular"
            break

        alt_m = x_nominal[0] - r_earth
        if alt_m <= float(terminal_altitude_m):
            terminated_reason = "terminal_alt"
            break

    return {
        "x_final": list(x_nominal),
        "reached_terminal": bool(terminated_reason == "terminal_alt"),
        "steps_used": int(steps_used),
        "terminated_reason": str(terminated_reason),
    }


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
    interval_screen_active: bool = True,
) -> PredictorCorrectorRolloutResult:
    """
    roll out one guidance candidate through nominal and interval dynamics

    input
    nominal start state candidate bank parameters horizon settings interval state and heat settings

    output
    one rollout result used by guidance scoring
    """
    if len(x_nominal_old) != 6:
        raise ValueError("x_nominal_old must have 6 components")
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive")
    if horizon_steps <= 0:
        raise ValueError("horizon_steps must be positive")

    local_supervisor_cfg = supervisor_cfg
    if bool(interval_screen_active):
        if local_supervisor_cfg is None and (heat_rate_limit is not None or heat_load_limit is not None):
            local_supervisor_cfg = IntervalSupervisorConfig(include_heating=True)
        elif local_supervisor_cfg is not None and not local_supervisor_cfg.include_heating:
            local_supervisor_cfg = copy.deepcopy(local_supervisor_cfg)
            local_supervisor_cfg.include_heating = True

        if local_supervisor_cfg is not None:
            if heat_rate_limit is not None:
                local_supervisor_cfg.heat_rate_limit = float(heat_rate_limit)
            if heat_load_limit is not None:
                local_supervisor_cfg.heat_load_limit = float(heat_load_limit)

    if bool(interval_screen_active):
        if x_interval_old is None:
            if local_supervisor_cfg is None:
                x_interval = nominal_state_to_interval_box(x_nominal_old)
            else:
                x_interval = inflate_nominal_state_to_interval_box(x_nominal_old, local_supervisor_cfg)
        else:
            x_interval = copy_interval_box(x_interval_old)

        candidate_heat_shield = clone_interval_heat_shield(heat_shield)
        interval_states: List[List[Interval]] = [copy_interval_box(x_interval)]
    else:
        x_interval = None
        candidate_heat_shield = None
        interval_states = []

    x_nominal = [float(v) for v in x_nominal_old]
    sigma_iv = promote(float(sigma_cmd_rad))
    nominal_states: List[List[float]] = [list(x_nominal)]
    step_summaries: List[RolloutStepSummary] = []

    max_heating_rate_interval: Optional[Interval] = None
    max_heat_load_interval: Optional[Interval] = None
    heat_feasible = True
    interval_valid = True
    interval_failure_kind = "none"
    violation_amount = 0.0
    heat_penalty = 0.0
    failure_reason = ""
    first_violation_step = -1
    first_violation_time_s = math.nan

    for step_idx in range(int(horizon_steps)):
        t_step_s = float(t0_s + step_idx * dt_s)
        x_nom_before = list(x_nominal)

        nominal_step = nominal_eom_step(
            x=x_nominal,
            sigma_rad=float(sigma_cmd_rad),
            params=params,
            dt_s=float(dt_s),
        )
        x_nominal = [float(v) for v in nominal_step["x_next"]]

        qdot_interval = Interval(0.0, 0.0)
        Q_interval = Interval(0.0, 0.0)
        x_interval_after: Optional[List[Interval]] = None

        if bool(interval_screen_active):
            annotation = build_interval_annotation(
                x_interval_old=x_interval,
                params=params,
                sigma_interval_used=sigma_iv,
                dt_s=float(dt_s),
                supervisor_cfg=local_supervisor_cfg,
                heat_shield=candidate_heat_shield,
                t_s=float(t_step_s),
                x_nominal_center=list(x_nom_before),
            )

            if not annotation.interval_valid:
                interval_valid = False
                interval_failure_kind = "numerical"
                failure_reason = str(annotation.interval_failure_reason or "interval_rollout_invalid")
                violation_amount = max(violation_amount, 1.0)

                if first_violation_step < 0:
                    first_violation_step = int(step_idx)
                    first_violation_time_s = float(t_step_s)

                step_summary = RolloutStepSummary(
                    step_index=int(step_idx),
                    time_s=float(t_step_s),
                    x_nominal_before=list(x_nom_before),
                    x_nominal_after=list(x_nominal),
                    x_interval_after=None,
                    max_heating_rate_interval=_hull_or_default(max_heating_rate_interval, None),
                    max_heat_load_interval=_hull_or_default(max_heat_load_interval, None),
                    heat_feasible_after_step=bool(heat_feasible),
                    interval_valid_after_step=False,
                    failure_reason=str(failure_reason),
                )
                step_summaries.append(step_summary)
                nominal_states.append(list(x_nominal))
                break

            x_interval = copy_interval_box(annotation.x_interval_new)
            candidate_heat_shield = annotation.heat_shield
            x_interval_after = copy_interval_box(x_interval)
            interval_states.append(copy_interval_box(x_interval))

            qdot_interval = annotation.heating_qdot_max_interval if annotation.heating_qdot_max_interval is not None else Interval(0.0, 0.0)
            Q_interval = annotation.heating_Q_max_interval if annotation.heating_Q_max_interval is not None else Interval(0.0, 0.0)

            max_heating_rate_interval = _hull_or_default(max_heating_rate_interval, qdot_interval)
            max_heat_load_interval = _hull_or_default(max_heat_load_interval, Q_interval)

            step_violation = 0.0
            step_reason_parts: List[str] = []

            if annotation.interval_heat_violation:
                interval_failure_kind = "heat"

                if heat_rate_limit is not None:
                    rate_violation = _positive_violation_amount(qdot_interval, float(heat_rate_limit))
                    if rate_violation > 0.0:
                        step_violation += rate_violation
                        step_reason_parts.append("heat_rate_limit")

                if heat_load_limit is not None:
                    load_violation = _positive_violation_amount(Q_interval, float(heat_load_limit))
                    if load_violation > 0.0:
                        step_violation += load_violation
                        step_reason_parts.append("heat_load_limit")

            if step_violation > 0.0:
                heat_feasible = False
                violation_amount = max(violation_amount, float(step_violation))
                heat_penalty += float(step_violation)

                if first_violation_step < 0:
                    first_violation_step = int(step_idx)
                    first_violation_time_s = float(t_step_s)

                if not failure_reason:
                    failure_reason = ",".join(step_reason_parts) if step_reason_parts else str(annotation.interval_failure_reason)

        step_summary = RolloutStepSummary(
            step_index=int(step_idx),
            time_s=float(t_step_s),
            x_nominal_before=list(x_nom_before),
            x_nominal_after=list(x_nominal),
            x_interval_after=x_interval_after,
            max_heating_rate_interval=qdot_interval,
            max_heat_load_interval=Q_interval,
            heat_feasible_after_step=bool(heat_feasible),
            interval_valid_after_step=bool(interval_valid),
            failure_reason=str(failure_reason),
        )

        nominal_states.append(list(x_nominal))
        step_summaries.append(step_summary)

    if max_heating_rate_interval is None:
        max_heating_rate_interval = Interval(0.0, 0.0)
    if max_heat_load_interval is None:
        max_heat_load_interval = Interval(0.0, 0.0)

    return PredictorCorrectorRolloutResult(
        sigma_cmd_rad=float(sigma_cmd_rad),
        horizon_steps=int(horizon_steps),
        nominal_states=nominal_states,
        interval_states=interval_states,
        step_summaries=step_summaries,
        final_nominal_state=list(x_nominal),
        final_interval_state=copy_interval_box(x_interval) if (bool(interval_screen_active) and interval_valid and x_interval is not None) else None,
        final_heat_shield=candidate_heat_shield,
        max_heating_rate_interval=max_heating_rate_interval,
        max_heat_load_interval=max_heat_load_interval,
        first_violation_step=int(first_violation_step),
        first_violation_time_s=float(first_violation_time_s),
        heat_feasible=bool(heat_feasible),
        interval_valid=bool(interval_valid),
        interval_screen_used=bool(interval_screen_active),
        interval_failure_kind=str(interval_failure_kind),
        violation_amount=float(violation_amount),
        heat_penalty=float(heat_penalty),
        failure_reason=str(failure_reason),
    )
