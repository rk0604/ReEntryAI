import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import AtmosphereModel
import constants
from interval_math import Interval, box_add, box_scalar_mul, promote


"""
math_3d.py

This file is the interval math module for translational reentry dynamics.

The nominal milestone 1 simulator path stays elsewhere.
This file only provides interval equations of motion and interval
supervisor utilities that can annotate the nominal simulator with
uncertainty envelopes.

Translational state convention

X = [r, phi, lam, V, gamma, chi]

r      radial distance from Earth center in meters
phi    geocentric latitude in radians
lam    longitude in radians
V      speed magnitude in m per s
gamma  flight path angle in radians
chi    heading angle in radians
"""


STATE_NAMES = ("r", "phi", "lam", "V", "gamma", "chi")


# ============================================================================
# Configuration and result models
# ============================================================================

@dataclass
class IntervalSupervisorConfig:
    """
    Half widths used to inflate a nominal translational state into an interval box.

    Limits are optional. When provided, a simple interval safety summary can be
    generated for logging.
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


@dataclass
class IntervalAnnotationResult:
    """
    Interval supervisor output for one annotated translational step.
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
    safety_status: str = "not_evaluated"
    safety_checks: Dict[str, str] = field(default_factory=dict)
    layer_indices: List[int] = field(default_factory=list)


# ============================================================================
# Small helpers used by the interval path
# ============================================================================

def _get_param(params: Dict[str, Any], *keys: str, default: Optional[float] = None) -> float:
    """
    Fetch the first matching numeric parameter from a dictionary.
    """
    for key in keys:
        if key in params:
            return float(params[key])

    if default is None:
        raise KeyError(f"Missing required parameter. Tried keys={keys}")

    return float(default)


def hull_intervals(intervals: List[Interval]) -> Interval:
    """
    Return the smallest interval containing all input intervals.
    """
    if not intervals:
        raise ValueError("interval list cannot be empty")

    out = intervals[0]
    for iv in intervals[1:]:
        out = out.hull(iv)
    return out


def nominal_state_to_interval_box(x_nominal: List[float]) -> List[Interval]:
    """
    Convert a nominal float state into a punctual interval box.
    """
    if len(x_nominal) != 6:
        raise ValueError("x_nominal must have 6 components")

    return [promote(float(v)) for v in x_nominal]


def inflate_nominal_state_to_interval_box(
    x_nominal: List[float],
    cfg: IntervalSupervisorConfig,
) -> List[Interval]:
    """
    Inflate a nominal float state into an interval box using configured half widths.
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
    Compute per component interval widths for a translational box.

    These widths are useful for ML logging because they tell the training
    pipeline how large the uncertainty envelope is in each state channel.
    """
    if len(box) != 6:
        raise ValueError("box must have 6 components")

    return {STATE_NAMES[i]: float(box[i].width()) for i in range(6)}


def interval_dynamic_pressure_bounds(rho_interval: Interval, V_interval: Interval) -> Interval:
    """
    Compute interval dynamic pressure bounds q = 0.5 rho V squared.
    """
    return promote(0.5) * rho_interval * V_interval.pow_int(2)


def atmosphere_interval_hull_from_state_box(x_box: List[Interval]) -> Dict[str, Any]:
    """
    Collect hulls across all atmospheric layers intersected by the altitude box.

    Layer hulls are used because an altitude interval can cross layer boundaries.
    The supervisor should return one conservative enclosure that covers all
    admissible layers touched by the current state box.
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
    temperature_intervals = [atm_by_layer[k]["T_K"] for k in layer_indices]
    pressure_intervals = [atm_by_layer[k]["p_Pa"] for k in layer_indices]
    density_intervals = [atm_by_layer[k]["rho_kgm3"] for k in layer_indices]

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
    return constants.HeatShield(
        radius_m=float(constants.HEAT_SHIELD_RADIUS_M),
        nose_radius_m=float(constants.HEAT_SHIELD_NOSE_RADIUS_M),
        num_rings=int(constants.HEAT_SHIELD_NUM_RINGS),
        radial_exp=float(constants.HEAT_SHIELD_RADIAL_EXP),
    )


def interval_heating_envelope(
    rho_interval: Interval,
    V_interval: Interval,
    dt_s: float,
    heat_shield: Optional[Any] = None,
) -> Dict[str, Any]:
    """
    Compute a simple interval heating envelope using the existing heat shield model.

    This stays lightweight on purpose. It produces a conservative heating
    annotation without turning the interval supervisor into the main simulator.
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
    Return a simple safety style summary for interval outputs.

    Each check is marked inside, mixed, or outside.
    The aggregate status is outside if any check is outside.
    It is mixed if no check is outside and at least one is mixed.
    It is inside if every evaluated check is inside.
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


# ============================================================================
# Interval aerodynamic force model
# ============================================================================

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

    The output is keyed by atmospheric layer because an altitude interval can
    intersect more than one layer.
    """
    z_geometric_altitude = constants.intv_geometric_altitude(r)
    atm_intv = AtmosphereModel.intv_US_Standard_ATM(z_geometric_altitude)

    if not atm_intv:
        zero = promote(0.0)
        return {
            -1: {
                "q": zero,
                "Dr": zero,
                "Dtheta": zero,
                "Dphi": zero,
                "Lr": zero,
                "Ltheta": zero,
                "Lphi": zero,
                "rho": zero,
            }
        }

    aero_forces_by_layer: Dict[int, Dict[str, Interval]] = {}

    CD = _get_param(params, "CD_const", "CD")
    CL = _get_param(params, "CL_const", "CL")
    surface_area = _get_param(params, "ref_area_m2", "S")

    cos_gamma = gamma.cos()
    sin_gamma = gamma.sin()
    cos_chi = chi.cos()
    sin_chi = chi.sin()
    cos_sigma = sigma.cos()
    sin_sigma = sigma.sin()

    for k in atm_intv.keys():
        rho = atm_intv[k]["rho_kgm3"] if atm_intv[k]["rho_kgm3"] is not None else promote(0.0)
        q = promote(0.5) * rho * V.pow_int(2)

        drag_mag = q * surface_area * CD
        lift_mag = q * surface_area * CL

        v_r = sin_gamma
        v_theta = cos_gamma * cos_chi
        v_phi = cos_gamma * sin_chi

        D_r = -drag_mag * v_r
        Dtheta = drag_mag * v_theta
        Dphi = -drag_mag * v_phi

        e1_r = cos_gamma
        e1_theta = -sin_gamma * cos_chi
        e1_phi = -sin_gamma * sin_chi

        e2_r = promote(0.0)
        e2_theta = sin_chi
        e2_phi = -cos_chi

        Lr = lift_mag * (e1_r * cos_sigma + e2_r * sin_sigma)
        Ltheta = lift_mag * (e1_theta * cos_sigma + e2_theta * sin_sigma)
        Lphi = lift_mag * (e1_phi * cos_sigma + e2_phi * sin_sigma)

        aero_forces_by_layer[k] = {
            "q": q,
            "Dr": D_r,
            "Dtheta": Dtheta,
            "Dphi": Dphi,
            "Lr": Lr,
            "Ltheta": Ltheta,
            "Lphi": Lphi,
            "rho": rho,
        }

    return aero_forces_by_layer


# ============================================================================
# Interval equations of motion
# ============================================================================

def intv_eom_3d(
    t: float,
    X: List[Interval],
    params: Dict[str, Any],
    sigma: Interval,
) -> Dict[str, Interval]:
    """
    Interval valued 3D equations of motion for atmospheric entry.
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

    m = _get_param(params, "mass_kg", "m")
    g = constants.intv_gravity(r)

    aero = intv_aero_forces(
        r=r,
        V=V,
        gamma=gamma,
        chi=chi,
        sigma=sigma,
        params=params,
    )

    Dr = hull_intervals([aero[layer]["Dr"] for layer in aero])
    Dtheta = hull_intervals([aero[layer]["Dtheta"] for layer in aero])
    Dphi = hull_intervals([aero[layer]["Dphi"] for layer in aero])
    Lr = hull_intervals([aero[layer]["Lr"] for layer in aero])
    Ltheta = hull_intervals([aero[layer]["Ltheta"] for layer in aero])
    Lphi = hull_intervals([aero[layer]["Lphi"] for layer in aero])

    cos_phi = phi.cos()
    cos_gamma = gamma.cos()
    sin_gamma = gamma.sin()
    sin_chi = chi.sin()
    cos_chi = chi.cos()

    r_dot = V * sin_gamma
    phi_dot = (V * cos_gamma * sin_chi) / r
    lam_dot = (V * cos_gamma * cos_chi) / (r * cos_phi)

    V_dot = (Dr + Lr) / m - g * sin_gamma
    gamma_dot = (Ltheta / (m * V)) + (V / r - g / V) * cos_gamma
    phi_tan = phi.sin() / cos_phi
    chi_dot = (Lphi / (m * V * cos_gamma)) + (V / r) * sin_chi * phi_tan

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
    dx = intv_eom_3d(t=t, X=X, params=params, sigma=sigma_iv)

    return [
        dx["r_dot"],
        dx["phi_dot"],
        dx["lam_dot"],
        dx["V_dot"],
        dx["gamma_dot"],
        dx["chi_dot"],
    ]


# ============================================================================
# Interval supervisor annotation helpers
# ============================================================================

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

    This function is meant to be called after the nominal simulator has already
    chosen the actual bank angle for the step.

    The interval path should annotate uncertainty around that flown decision,
    not create a second control mode. For that reason the nominal bank angle is
    typically converted into a punctual interval before calling this function.
    """
    if dt_s is None:
        dt_s = float(getattr(constants, "ENTRY_DT_S", 0.25))

    if len(x_interval_old) != 6:
        raise ValueError("x_interval_old must have 6 interval components")

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

    # The interval derivative is evaluated around the current state box and then
    # pushed one Euler step forward. This mirrors the nominal fixed step update
    # while preserving a deterministic enclosure for uncertainty logging.
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

    If x_interval_old is not provided, the nominal state is either converted into
    a punctual interval box or inflated using the supplied supervisor config.
    """
    if x_interval_old is None:
        if supervisor_cfg is None:
            x_interval_old = nominal_state_to_interval_box(x_nominal_old)
        else:
            x_interval_old = inflate_nominal_state_to_interval_box(x_nominal_old, supervisor_cfg)

    # The nominal simulator owns control and attitude. The interval path should
    # use the actual bank angle already flown by the nominal path, so the bank
    # input is turned into a punctual interval here.
    sigma_interval_used = promote(float(sigma_actual_after_rad))

    return build_interval_annotation(
        x_interval_old=x_interval_old,
        params=params,
        sigma_interval_used=sigma_interval_used,
        dt_s=dt_s,
        supervisor_cfg=supervisor_cfg,
        heat_shield=heat_shield,
        t_s=t_s,
    )