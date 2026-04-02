import math
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple

import AtmosphereModel
import constants
from control import ReentryState, SigmaControlStack
from ReactionControl import CapsuleRCSSystem


"""
milestone1_nominal.py

This file owns the nominal milestone 1 simulation path.

This includes float translational dynamics, attitude helpers,
control bridging, and the closed loop nominal step.

math_3d.py can then stay focused on interval math and interval
supervisor annotations only.

Translational state convention

X = [r, phi, lam, V, gamma, chi]

r      radial distance from Earth center in meters
phi    geocentric latitude in radians
lam    longitude in radians
V      speed magnitude in m per s
gamma  flight path angle in radians
chi    heading angle in radians
"""


# ============================================================================
# Small math helpers
# ============================================================================

def wrap_to_pi(angle_rad: float) -> float:
    """
    Wrap an angle into the range from minus pi to pi.
    """
    while angle_rad > math.pi:
        angle_rad -= 2.0 * math.pi
    while angle_rad < -math.pi:
        angle_rad += 2.0 * math.pi
    return angle_rad


def clamp(value: float, lower: float, upper: float) -> float:
    """
    Clamp a scalar into a closed interval.
    """
    return max(lower, min(upper, value))


def normalize(vec: List[float], eps: Optional[float] = None) -> List[float]:
    """
    Normalize a 3D vector and return a plain Python list.
    """
    if eps is None:
        eps = getattr(constants, "EPS_NORM", 1.0e-12)

    x = float(vec[0])
    y = float(vec[1])
    z = float(vec[2])

    n = math.sqrt(x * x + y * y + z * z)
    if n < eps:
        return [0.0, 0.0, 0.0]

    return [x / n, y / n, z / n]


def dot3(a: List[float], b: List[float]) -> float:
    """
    Return the 3D dot product of two vectors.
    """
    return float(a[0]) * float(b[0]) + float(a[1]) * float(b[1]) + float(a[2]) * float(b[2])


def cross3(a: List[float], b: List[float]) -> List[float]:
    """
    Return the 3D cross product a x b.
    """
    ax, ay, az = float(a[0]), float(a[1]), float(a[2])
    bx, by, bz = float(b[0]), float(b[1]), float(b[2])

    return [
        ay * bz - az * by,
        az * bx - ax * bz,
        ax * by - ay * bx,
    ]


def mat_vec_mul(M: List[List[float]], v: List[float]) -> List[float]:
    """
    Multiply a 3 by 3 matrix by a 3D vector.
    """
    return [
        float(M[0][0]) * float(v[0]) + float(M[0][1]) * float(v[1]) + float(M[0][2]) * float(v[2]),
        float(M[1][0]) * float(v[0]) + float(M[1][1]) * float(v[1]) + float(M[1][2]) * float(v[2]),
        float(M[2][0]) * float(v[0]) + float(M[2][1]) * float(v[1]) + float(M[2][2]) * float(v[2]),
    ]


def mat_mul_3x3(A: List[List[float]], B: List[List[float]]) -> List[List[float]]:
    """
    Multiply two 3 by 3 matrices.
    """
    return [
        [
            float(A[i][0]) * float(B[0][j])
            + float(A[i][1]) * float(B[1][j])
            + float(A[i][2]) * float(B[2][j])
            for j in range(3)
        ]
        for i in range(3)
    ]


def safe_nonzero(value: float, eps: float = 1.0e-12) -> float:
    """
    Prevent division by zero in the float dynamics path.
    """
    if abs(value) >= eps:
        return value
    return eps if value >= 0.0 else -eps


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


def add_scaled_state(x: List[float], dx: List[float], dt_s: float) -> List[float]:
    """
    Return x plus dt times dx for a flat state vector.
    """
    return [float(x[i]) + float(dt_s) * float(dx[i]) for i in range(len(x))]


# ============================================================================
# Milestone 1 capsule attitude state and frame helpers
# ============================================================================

@dataclass
class CapsuleAttitudeState:
    """
    Minimal rigid body attitude state for milestone 1.

    q_ba stores the quaternion that maps body frame vectors into
    aerodynamic frame vectors.

    omega_b_rad_s stores body angular rates relative to the
    aerodynamic frame.

    sigma_rel_rad stores the actual bank angle.
    """
    q_ba: List[float]
    omega_b_rad_s: List[float]
    sigma_rel_rad: float


@dataclass
class Milestone1StepResult:
    """
    Result packet for one nominal closed loop milestone 1 step.
    """
    t_s: float
    x_old: List[float]
    x_new: List[float]
    att_old: CapsuleAttitudeState
    att_new: CapsuleAttitudeState
    sigma_actual_before_rad: float
    sigma_actual_after_rad: float
    control_out: Any
    roll_step: Any
    dx_trans: List[float]


def make_initial_capsule_attitude(sigma0_rad: float = 0.0) -> CapsuleAttitudeState:
    """
    Construct an initial roll only capsule attitude state.
    """
    q_ba = quat_from_axis_angle([0.0, 0.0, 1.0], sigma0_rad)

    return CapsuleAttitudeState(
        q_ba=q_ba,
        omega_b_rad_s=[0.0, 0.0, 0.0],
        sigma_rel_rad=float(wrap_to_pi(sigma0_rad)),
    )


# ============================================================================
# Quaternion helpers
# ============================================================================

def quat_normalize(q: List[float]) -> List[float]:
    """
    Normalize a quaternion [w, x, y, z].
    """
    w, x, y, z = [float(v) for v in q]
    n = math.sqrt(w * w + x * x + y * y + z * z)

    if n < getattr(constants, "EPS_NORM", 1.0e-12):
        return [1.0, 0.0, 0.0, 0.0]

    return [w / n, x / n, y / n, z / n]


def quat_from_axis_angle(axis: List[float], angle_rad: float) -> List[float]:
    """
    Build a quaternion [w, x, y, z] from an axis and angle.
    """
    ax = normalize(axis)
    half = 0.5 * float(angle_rad)
    s = math.sin(half)

    return quat_normalize([
        math.cos(half),
        float(ax[0]) * s,
        float(ax[1]) * s,
        float(ax[2]) * s,
    ])


def quat_to_dcm(q: List[float]) -> List[List[float]]:
    """
    Convert a quaternion [w, x, y, z] into a 3 by 3 direction cosine matrix.
    """
    w, x, y, z = quat_normalize(q)

    return [
        [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w),       2.0 * (x * z + y * w)],
        [2.0 * (x * y + z * w),       1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w)],
        [2.0 * (x * z - y * w),       2.0 * (y * z + x * w),       1.0 - 2.0 * (x * x + y * y)],
    ]


# ============================================================================
# Local and aerodynamic frame helpers
# ============================================================================

def local_level_frame(phi_rad: float, lam_rad: float) -> Tuple[List[float], List[float], List[float]]:
    """
    Return local north, east, and up unit vectors in the world frame.
    """
    north = [
        -math.sin(phi_rad) * math.cos(lam_rad),
        -math.sin(phi_rad) * math.sin(lam_rad),
        math.cos(phi_rad),
    ]

    east = [
        -math.sin(lam_rad),
        math.cos(lam_rad),
        0.0,
    ]

    up = [
        math.cos(phi_rad) * math.cos(lam_rad),
        math.cos(phi_rad) * math.sin(lam_rad),
        math.sin(phi_rad),
    ]

    return normalize(north), normalize(east), normalize(up)


def position_world_from_state(x: List[float]) -> List[float]:
    """
    Convert [r, phi, lam, ...] into a world frame position vector.
    """
    r = float(x[0])
    phi = float(x[1])
    lam = float(x[2])

    _, _, up = local_level_frame(phi, lam)
    return [r * up[0], r * up[1], r * up[2]]


def velocity_hat_from_state(x: List[float]) -> List[float]:
    """
    Build the unit velocity direction in the world frame.
    """
    _, phi, lam, _, gamma, chi = [float(v) for v in x]
    north, east, up = local_level_frame(phi, lam)

    v_hat = [
        math.cos(gamma) * math.sin(chi) * north[0] + math.cos(gamma) * math.cos(chi) * east[0] + math.sin(gamma) * up[0],
        math.cos(gamma) * math.sin(chi) * north[1] + math.cos(gamma) * math.cos(chi) * east[1] + math.sin(gamma) * up[1],
        math.cos(gamma) * math.sin(chi) * north[2] + math.cos(gamma) * math.cos(chi) * east[2] + math.sin(gamma) * up[2],
    ]

    return normalize(v_hat)


def aero_frame_dcm_from_state(x: List[float]) -> List[List[float]]:
    """
    Build R_wa, the DCM mapping aerodynamic frame vectors into world frame vectors.
    """
    _, phi, lam, _, _, _ = [float(v) for v in x]

    north, _, up = local_level_frame(phi, lam)
    v_hat = velocity_hat_from_state(x)

    z_a_w = normalize([-v_hat[0], -v_hat[1], -v_hat[2]])

    x_a_w = [
        up[0] - dot3(up, z_a_w) * z_a_w[0],
        up[1] - dot3(up, z_a_w) * z_a_w[1],
        up[2] - dot3(up, z_a_w) * z_a_w[2],
    ]

    if math.sqrt(dot3(x_a_w, x_a_w)) < getattr(constants, "EPS_NORM", 1.0e-12):
        x_a_w = [
            north[0] - dot3(north, z_a_w) * z_a_w[0],
            north[1] - dot3(north, z_a_w) * z_a_w[1],
            north[2] - dot3(north, z_a_w) * z_a_w[2],
        ]

    x_a_w = normalize(x_a_w)
    y_a_w = normalize(cross3(z_a_w, x_a_w))

    return [
        [x_a_w[0], y_a_w[0], z_a_w[0]],
        [x_a_w[1], y_a_w[1], z_a_w[1]],
        [x_a_w[2], y_a_w[2], z_a_w[2]],
    ]


def sigma_from_q_ba(q_ba: List[float]) -> float:
    """
    Extract actual bank angle from the body to aero attitude quaternion.
    """
    R_ba = quat_to_dcm(q_ba)
    x_b_in_a = [R_ba[0][0], R_ba[1][0], R_ba[2][0]]
    sigma = math.atan2(x_b_in_a[1], x_b_in_a[0])

    sigma_sign = getattr(constants, "SIGMA_SIGN", 1.0)
    return wrap_to_pi(sigma_sign * sigma)


def body_to_world_dcm(x_trans: List[float], q_ba: List[float]) -> List[List[float]]:
    """
    Return R_wb, the DCM mapping body frame vectors into world frame vectors.
    """
    R_wa = aero_frame_dcm_from_state(x_trans)
    R_ba = quat_to_dcm(q_ba)
    return mat_mul_3x3(R_wa, R_ba)


def heat_shield_normal_world(x_trans: List[float], q_ba: List[float]) -> List[float]:
    """
    Return the heat shield normal in the world frame.
    """
    R_wb = body_to_world_dcm(x_trans, q_ba)
    return mat_vec_mul(R_wb, [0.0, 0.0, 1.0])


def thruster_world_pose(
    x_trans: List[float],
    q_ba: List[float],
    thruster_position_b_m: List[float],
    thruster_direction_b_unit: List[float],
) -> Tuple[List[float], List[float]]:
    """
    Convert a body fixed thruster pose into world coordinates.
    """
    R_wb = body_to_world_dcm(x_trans, q_ba)
    r_cm_w = position_world_from_state(x_trans)

    r_thr_w_local = mat_vec_mul(R_wb, thruster_position_b_m)
    d_thr_w = normalize(mat_vec_mul(R_wb, thruster_direction_b_unit))

    r_thr_w = [
        r_cm_w[0] + r_thr_w_local[0],
        r_cm_w[1] + r_thr_w_local[1],
        r_cm_w[2] + r_thr_w_local[2],
    ]

    return r_thr_w, d_thr_w


def step_capsule_roll_state(
    att: CapsuleAttitudeState,
    tau_rcs_b_Nm: List[float],
    Izz_kgm2: Optional[float] = None,
    dt_s: Optional[float] = None,
) -> CapsuleAttitudeState:
    """
    Advance the milestone 1 roll only attitude state by one fixed step.
    """
    if Izz_kgm2 is None:
        Izz_kgm2 = float(getattr(constants, "CAPSULE_IZZ_KGM2", 20000.0))
    if dt_s is None:
        dt_s = float(getattr(constants, "ENTRY_DT_S", 0.25))

    roll_rate_old = float(att.omega_b_rad_s[2])
    tau_roll = float(tau_rcs_b_Nm[2])

    roll_accel = tau_roll / safe_nonzero(Izz_kgm2)
    roll_rate_new = roll_rate_old + dt_s * roll_accel
    sigma_new = wrap_to_pi(float(att.sigma_rel_rad) + dt_s * roll_rate_new)

    q_ba_new = quat_from_axis_angle([0.0, 0.0, 1.0], sigma_new)

    return CapsuleAttitudeState(
        q_ba=q_ba_new,
        omega_b_rad_s=[0.0, 0.0, roll_rate_new],
        sigma_rel_rad=sigma_new,
    )


# ============================================================================
# Float aerodynamic force model
# ============================================================================

def aero_forces(
    r: float,
    V: float,
    gamma: float,
    chi: float,
    params: Dict[str, Any],
    sigma_actual_rad: Optional[float] = None,
    sigma_cmd_rad: Optional[float] = None,
) -> Tuple[float, float, float, float, float, float, float, float]:
    """
    Compute drag and lift components in the spherical 3D frame.
    """
    sigma_used = float(
        sigma_actual_rad
        if sigma_actual_rad is not None
        else (sigma_cmd_rad if sigma_cmd_rad is not None else 0.0)
    )

    altitude_m = float(r) - float(constants.RADIUS_EARTH)
    altitude_m = max(0.0, altitude_m)

    atm = AtmosphereModel.US_Standard_ATM(altitude_m)
    rho = atm["rho_kgm3"] if atm["rho_kgm3"] is not None else 0.0

    q = 0.5 * rho * V * V

    CD = _get_param(params, "CD_const", "CD")
    CL = _get_param(params, "CL_const", "CL")
    surface_area = _get_param(params, "ref_area_m2", "S")

    drag_mag = q * surface_area * CD
    lift_mag = q * surface_area * CL

    v_r = math.sin(gamma)
    v_theta = math.cos(gamma) * math.cos(chi)
    v_phi = math.cos(gamma) * math.sin(chi)

    D_r = -drag_mag * v_r
    Dtheta = drag_mag * v_theta
    Dphi = -drag_mag * v_phi

    e1_r = math.cos(gamma)
    e1_theta = -math.sin(gamma) * math.cos(chi)
    e1_phi = -math.sin(gamma) * math.sin(chi)

    e2_r = 0.0
    e2_theta = math.sin(chi)
    e2_phi = -math.cos(chi)

    cos_sig = math.cos(sigma_used)
    sin_sig = math.sin(sigma_used)

    Lr = lift_mag * (e1_r * cos_sig + e2_r * sin_sig)
    Ltheta = lift_mag * (e1_theta * cos_sig + e2_theta * sin_sig)
    Lphi = lift_mag * (e1_phi * cos_sig + e2_phi * sin_sig)

    return D_r, Dtheta, Dphi, Lr, Ltheta, Lphi, float(rho), float(q)


# ============================================================================
# Float equations of motion
# ============================================================================

def eom_3d(
    t: float,
    X: List[float],
    params: Dict[str, Any],
    sigma_actual_rad: Optional[float] = None,
    sigma_cmd_rad: Optional[float] = None,
) -> Dict[str, float]:
    """
    Float 3D spherical equations of motion for atmospheric entry.
    """
    r, phi, lam, V, gamma, chi = [float(v) for v in X]

    m = _get_param(params, "mass_kg", "m")
    g = float(constants.gravity(r))

    Dr, Dtheta, Dphi, Lr, Ltheta, Lphi, _, _ = aero_forces(
        r=r,
        V=V,
        gamma=gamma,
        chi=chi,
        params=params,
        sigma_actual_rad=sigma_actual_rad,
        sigma_cmd_rad=sigma_cmd_rad,
    )

    cos_phi = safe_nonzero(math.cos(phi))
    cos_gamma = safe_nonzero(math.cos(gamma))
    V_safe = safe_nonzero(V)
    r_safe = safe_nonzero(r)

    r_dot = V * math.sin(gamma)
    phi_dot = (V * math.cos(gamma) * math.sin(chi)) / r_safe
    lam_dot = (V * math.cos(gamma) * math.cos(chi)) / (r_safe * cos_phi)

    V_dot = (Dr + Lr) / m - g * math.sin(gamma)
    gamma_dot = (Ltheta / (m * V_safe)) + (V_safe / r_safe - g / V_safe) * math.cos(gamma)
    phi_tan = math.sin(phi) / cos_phi
    chi_dot = (Lphi / (m * V_safe * cos_gamma)) + (V_safe / r_safe) * math.sin(chi) * phi_tan

    return {
        "r_dot": float(r_dot),
        "phi_dot": float(phi_dot),
        "lam_dot": float(lam_dot),
        "V_dot": float(V_dot),
        "gamma_dot": float(gamma_dot),
        "chi_dot": float(chi_dot),
    }


def f_float(
    t: float,
    X: List[float],
    params: Dict[str, Any],
    sigma_actual_rad: Optional[float] = None,
    sigma_cmd_rad: Optional[float] = None,
) -> List[float]:
    """
    Ordered derivative adapter for the float equations of motion.
    """
    dx = eom_3d(
        t=t,
        X=X,
        params=params,
        sigma_actual_rad=sigma_actual_rad,
        sigma_cmd_rad=sigma_cmd_rad,
    )

    return [
        dx["r_dot"],
        dx["phi_dot"],
        dx["lam_dot"],
        dx["V_dot"],
        dx["gamma_dot"],
        dx["chi_dot"],
    ]


# ============================================================================
# Control bridge helpers
# ============================================================================

def make_control_step_fn(
    control: SigmaControlStack,
    aero_force_fn: Callable[..., Tuple[float, float, float, float, float, float, float, float]],
    params: Dict[str, Any],
    dt_s: Optional[float] = None,
) -> Callable[[float, List[float], float, float], Any]:
    """
    Build a helper with the form control_step_fn(t, x, sigma_actual_rad, roll_rate_rad_s).
    """
    if dt_s is None:
        dt_s = float(getattr(constants, "ENTRY_DT_S", 0.25))

    mass_kg = _get_param(params, "mass_kg", "m")

    def control_step_fn(
        t: float,
        x: List[float],
        sigma_actual_rad: float,
        roll_rate_rad_s: float,
    ):
        r, phi, lam, V, gamma, chi = [float(v) for v in x]

        state = ReentryState(
            r_m=r,
            phi_rad=phi,
            lam_rad=lam,
            V_mps=V,
            gamma_rad=gamma,
            chi_rad=chi,
            sigma_actual_rad=float(sigma_actual_rad),
            roll_rate_rad_s=float(roll_rate_rad_s),
        )

        Dr, Dtheta, Dphi, Lr, Ltheta, Lphi, _, _ = aero_force_fn(
            r=r,
            V=V,
            gamma=gamma,
            chi=chi,
            params=params,
            sigma_actual_rad=sigma_actual_rad,
        )

        lift_mag = math.sqrt(Lr * Lr + Ltheta * Ltheta + Lphi * Lphi)
        drag_mag = math.sqrt(Dr * Dr + Dtheta * Dtheta + Dphi * Dphi)

        return control.step(
            t_s=float(t),
            dt_s=float(dt_s),
            state=state,
            lift_N=float(lift_mag),
            drag_N=float(drag_mag),
            mass_kg=float(mass_kg),
        )

    return control_step_fn


def make_sigma_fn(
    control: SigmaControlStack,
    aero_force_fn: Callable[..., Tuple[float, float, float, float, float, float, float, float]],
    params: Dict[str, Any],
    sigma_actual_provider: Callable[[], Tuple[float, float]],
    dt_s: Optional[float] = None,
) -> Callable[[float, List[float]], float]:
    """
    Backward compatible adapter for older code that still wants sigma_fn(t, x).
    """
    control_step_fn = make_control_step_fn(
        control=control,
        aero_force_fn=aero_force_fn,
        params=params,
        dt_s=dt_s,
    )

    def sigma_fn(t: float, x: List[float]) -> float:
        sigma_actual_rad, roll_rate_rad_s = sigma_actual_provider()

        _ = control_step_fn(
            t=float(t),
            x=list(x),
            sigma_actual_rad=float(sigma_actual_rad),
            roll_rate_rad_s=float(roll_rate_rad_s),
        )

        return float(sigma_actual_rad)

    return sigma_fn


# ============================================================================
# Milestone 1 closed loop stepping
# ============================================================================

def step_closed_loop_milestone1(
    t_s: float,
    x_trans: List[float],
    att: CapsuleAttitudeState,
    params: Dict[str, Any],
    control_step_fn: Callable[[float, List[float], float, float], Any],
    rcs: CapsuleRCSSystem,
    dt_s: Optional[float] = None,
    Izz_kgm2: Optional[float] = None,
) -> Milestone1StepResult:
    """
    Advance one full nominal milestone 1 closed loop step.

    The nominal step runs first because control, thruster firing,
    attitude propagation, and actual bank angle all belong to the
    primary simulator path.
    """
    if dt_s is None:
        dt_s = float(getattr(constants, "ENTRY_DT_S", 0.25))
    if Izz_kgm2 is None:
        Izz_kgm2 = float(getattr(constants, "CAPSULE_IZZ_KGM2", 20000.0))

    x_old = [float(v) for v in x_trans]
    att_old = CapsuleAttitudeState(
        q_ba=list(att.q_ba),
        omega_b_rad_s=list(att.omega_b_rad_s),
        sigma_rel_rad=float(att.sigma_rel_rad),
    )

    sigma_actual_before_rad = float(sigma_from_q_ba(att.q_ba))
    roll_rate_rad_s = float(att.omega_b_rad_s[2])

    control_out = control_step_fn(
        float(t_s),
        list(x_old),
        float(sigma_actual_before_rad),
        float(roll_rate_rad_s),
    )

    roll_step = rcs.step_roll_channel(
        tau_roll_cmd_Nm=float(control_out.tau_roll_cmd_Nm),
        dt_s=float(dt_s),
    )

    att_new = step_capsule_roll_state(
        att=att_old,
        tau_rcs_b_Nm=list(roll_step.wrench.torque_b_Nm),
        Izz_kgm2=float(Izz_kgm2),
        dt_s=float(dt_s),
    )

    sigma_actual_after_rad = float(sigma_from_q_ba(att_new.q_ba))
    att_new.sigma_rel_rad = sigma_actual_after_rad

    dx_trans = f_float(
        t=float(t_s),
        X=list(x_old),
        params=params,
        sigma_actual_rad=float(sigma_actual_after_rad),
    )

    x_new = add_scaled_state(x_old, dx_trans, float(dt_s))

    return Milestone1StepResult(
        t_s=float(t_s),
        x_old=x_old,
        x_new=x_new,
        att_old=att_old,
        att_new=att_new,
        sigma_actual_before_rad=float(sigma_actual_before_rad),
        sigma_actual_after_rad=float(sigma_actual_after_rad),
        control_out=control_out,
        roll_step=roll_step,
        dx_trans=dx_trans,
    )