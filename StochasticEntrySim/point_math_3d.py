import math
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple

import AtmosphereModel
import constants
from control import ReentryState, SigmaControlStack
from ReactionControl import CapsuleRCSSystem


"""
milestone1_nominal.py

This file contains the nominal milestone 1 simulation path.

It includes the float translational dynamics, attitude helpers,
control bridge logic, and the closed loop nominal step.

This keeps math_3d.py focused on interval math and interval
supervisor annotations.

Translational state convention

X = [r, phi, lam, V, gamma, chi]

r      radial distance from Earth center in meters
phi    geocentric latitude in radians
lam    longitude in radians
V      speed magnitude in meters per second
gamma  flight path angle in radians
chi    heading angle in radians
"""


# Small math helpers used throughout the nominal simulation path.

def wrap_to_pi(angle_rad: float) -> float:
    """
    Wrap an angle into the range from negative pi to positive pi.
    """
    # Repeatedly subtract full turns until the angle is inside the target range.
    while angle_rad > math.pi:
        angle_rad -= 2.0 * math.pi

    # Repeatedly add full turns until the angle is inside the target range.
    while angle_rad < -math.pi:
        angle_rad += 2.0 * math.pi

    # Return the wrapped angle.
    return angle_rad


def clamp(value: float, lower: float, upper: float) -> float:
    """
    Clamp a scalar into a closed interval.
    """
    # Limit the value so it cannot go below the lower bound or above the upper bound.
    return max(lower, min(upper, value))


def normalize(vec: List[float], eps: Optional[float] = None) -> List[float]:
    """
    Normalize a three component vector and return a plain Python list.
    """
    # Pull a default norm tolerance from constants if one was not provided.
    if eps is None:
        eps = getattr(constants, "EPS_NORM", 1.0e-12)

    # Convert all components to float so downstream math stays consistent.
    x = float(vec[0])
    y = float(vec[1])
    z = float(vec[2])

    # Compute the Euclidean norm of the vector.
    n = math.sqrt(x * x + y * y + z * z)

    # Return a zero vector if the norm is too small to normalize safely.
    if n < eps:
        return [0.0, 0.0, 0.0]

    # Return the normalized vector components.
    return [x / n, y / n, z / n]


def dot3(a: List[float], b: List[float]) -> float:
    """
    Return the three dimensional dot product of two vectors.
    """
    # Multiply corresponding components and sum them.
    return float(a[0]) * float(b[0]) + float(a[1]) * float(b[1]) + float(a[2]) * float(b[2])


def cross3(a: List[float], b: List[float]) -> List[float]:
    """
    Return the three dimensional cross product a cross b.
    """
    # Read the first vector components as floats.
    ax, ay, az = float(a[0]), float(a[1]), float(a[2])

    # Read the second vector components as floats.
    bx, by, bz = float(b[0]), float(b[1]), float(b[2])

    # Return the standard cross product result.
    return [
        ay * bz - az * by,
        az * bx - ax * bz,
        ax * by - ay * bx,
    ]


def mat_vec_mul(M: List[List[float]], v: List[float]) -> List[float]:
    """
    Multiply a three by three matrix by a three component vector.
    """
    # Compute each output row explicitly for clarity and predictable ordering.
    return [
        float(M[0][0]) * float(v[0]) + float(M[0][1]) * float(v[1]) + float(M[0][2]) * float(v[2]),
        float(M[1][0]) * float(v[0]) + float(M[1][1]) * float(v[1]) + float(M[1][2]) * float(v[2]),
        float(M[2][0]) * float(v[0]) + float(M[2][1]) * float(v[1]) + float(M[2][2]) * float(v[2]),
    ]


def mat_mul_3x3(A: List[List[float]], B: List[List[float]]) -> List[List[float]]:
    """
    Multiply two three by three matrices.
    """
    # Build the output matrix row by row using the standard matrix product rule.
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
    # Return the original value if it is already far enough from zero.
    if abs(value) >= eps:
        return value

    # Preserve the sign when replacing a near zero value with a small safe fallback.
    return eps if value >= 0.0 else -eps


def _get_param(params: Dict[str, Any], *keys: str, default: Optional[float] = None) -> float:
    """
    Fetch the first matching numeric parameter from a dictionary.
    """
    # Check each candidate key in order and return the first one that exists.
    for key in keys:
        if key in params:
            return float(params[key])

    # Raise an error if no key matched and no default was supplied.
    if default is None:
        raise KeyError(f"Missing required parameter. Tried keys={keys}")

    # Return the fallback value when provided.
    return float(default)


def add_scaled_state(x: List[float], dx: List[float], dt_s: float) -> List[float]:
    """
    Return x plus dt times dx for a flat state vector.
    """
    # Apply a forward Euler style state update component by component.
    return [float(x[i]) + float(dt_s) * float(dx[i]) for i in range(len(x))]


# Capsule attitude state and helper result objects used by milestone 1.

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
    # Build a pure roll quaternion about the aerodynamic z axis.
    q_ba = quat_from_axis_angle([0.0, 0.0, 1.0], sigma0_rad)

    # Start with zero body rate and the wrapped initial bank angle.
    return CapsuleAttitudeState(
        q_ba=q_ba,
        omega_b_rad_s=[0.0, 0.0, 0.0],
        sigma_rel_rad=float(wrap_to_pi(sigma0_rad)),
    )


# Quaternion helpers used by the roll attitude model.

def quat_normalize(q: List[float]) -> List[float]:
    """
    Normalize a quaternion [w, x, y, z].
    """
    # Convert all quaternion components to floats.
    w, x, y, z = [float(v) for v in q]

    # Compute the quaternion magnitude.
    n = math.sqrt(w * w + x * x + y * y + z * z)

    # Fall back to identity if the input magnitude is too small.
    if n < getattr(constants, "EPS_NORM", 1.0e-12):
        return [1.0, 0.0, 0.0, 0.0]

    # Return the unit quaternion.
    return [w / n, x / n, y / n, z / n]


def quat_from_axis_angle(axis: List[float], angle_rad: float) -> List[float]:
    """
    Build a quaternion [w, x, y, z] from an axis and angle.
    """
    # Normalize the rotation axis before constructing the quaternion.
    ax = normalize(axis)

    # Use the half angle form of the axis angle quaternion equation.
    half = 0.5 * float(angle_rad)
    s = math.sin(half)

    # Build and normalize the quaternion to avoid drift from numerical noise.
    return quat_normalize([
        math.cos(half),
        float(ax[0]) * s,
        float(ax[1]) * s,
        float(ax[2]) * s,
    ])


def quat_to_dcm(q: List[float]) -> List[List[float]]:
    """
    Convert a quaternion [w, x, y, z] into a three by three direction cosine matrix.
    """
    # Normalize first so the resulting matrix stays orthonormal.
    w, x, y, z = quat_normalize(q)

    # Return the standard quaternion to DCM conversion.
    return [
        [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w),       2.0 * (x * z + y * w)],
        [2.0 * (x * y + z * w),       1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w)],
        [2.0 * (x * z - y * w),       2.0 * (y * z + x * w),       1.0 - 2.0 * (x * x + y * y)],
    ]


# Frame helpers that move between local level, aerodynamic, body, and world coordinates.

def local_level_frame(phi_rad: float, lam_rad: float) -> Tuple[List[float], List[float], List[float]]:
    """
    Return local north, east, and up unit vectors in the world frame.
    """
    # Build the local north direction from latitude and longitude.
    north = [
        -math.sin(phi_rad) * math.cos(lam_rad),
        -math.sin(phi_rad) * math.sin(lam_rad),
        math.cos(phi_rad),
    ]

    # Build the local east direction.
    east = [
        -math.sin(lam_rad),
        math.cos(lam_rad),
        0.0,
    ]

    # Build the local up direction pointing radially away from Earth center.
    up = [
        math.cos(phi_rad) * math.cos(lam_rad),
        math.cos(phi_rad) * math.sin(lam_rad),
        math.sin(phi_rad),
    ]

    # Normalize all three basis directions before returning them.
    return normalize(north), normalize(east), normalize(up)


def position_world_from_state(x: List[float]) -> List[float]:
    """
    Convert [r, phi, lam, ...] into a world frame position vector.
    """
    # Read the spherical position states.
    r = float(x[0])
    phi = float(x[1])
    lam = float(x[2])

    # Reuse the local up direction as the outward radial direction.
    _, _, up = local_level_frame(phi, lam)

    # Scale the radial unit vector by the radial distance.
    return [r * up[0], r * up[1], r * up[2]]


def velocity_hat_from_state(x: List[float]) -> List[float]:
    """
    Build the unit velocity direction in the world frame.
    """
    # Read the translational state as floats.
    _, phi, lam, _, gamma, chi = [float(v) for v in x]

    # Build the local level basis at the current position.
    north, east, up = local_level_frame(phi, lam)

    # Resolve the velocity direction into world coordinates using gamma and chi.
    v_hat = [
        math.cos(gamma) * math.sin(chi) * north[0] + math.cos(gamma) * math.cos(chi) * east[0] + math.sin(gamma) * up[0],
        math.cos(gamma) * math.sin(chi) * north[1] + math.cos(gamma) * math.cos(chi) * east[1] + math.sin(gamma) * up[1],
        math.cos(gamma) * math.sin(chi) * north[2] + math.cos(gamma) * math.cos(chi) * east[2] + math.sin(gamma) * up[2],
    ]

    # Normalize before returning so it behaves as a unit direction vector.
    return normalize(v_hat)


def aero_frame_dcm_from_state(x: List[float]) -> List[List[float]]:
    """
    Build R_wa, the DCM mapping aerodynamic frame vectors into world frame vectors.
    """
    # Read only the state components needed to form the local frame.
    _, phi, lam, _, _, _ = [float(v) for v in x]

    # Get the local north direction and the radial up direction.
    north, _, up = local_level_frame(phi, lam)

    # Build the velocity direction in world coordinates.
    v_hat = velocity_hat_from_state(x)

    # Define the aerodynamic z axis opposite the velocity direction.
    z_a_w = normalize([-v_hat[0], -v_hat[1], -v_hat[2]])

    # Project the radial up direction into the plane normal to z_a_w to form x_a_w.
    x_a_w = [
        up[0] - dot3(up, z_a_w) * z_a_w[0],
        up[1] - dot3(up, z_a_w) * z_a_w[1],
        up[2] - dot3(up, z_a_w) * z_a_w[2],
    ]

    # If the up vector is nearly aligned with z_a_w, use north as a backup reference.
    if math.sqrt(dot3(x_a_w, x_a_w)) < getattr(constants, "EPS_NORM", 1.0e-12):
        x_a_w = [
            north[0] - dot3(north, z_a_w) * z_a_w[0],
            north[1] - dot3(north, z_a_w) * z_a_w[1],
            north[2] - dot3(north, z_a_w) * z_a_w[2],
        ]

    # Normalize the aerodynamic x axis after projection.
    x_a_w = normalize(x_a_w)

    # Form the aerodynamic y axis from the cross product.
    y_a_w = normalize(cross3(z_a_w, x_a_w))

    # Return the world from aerodynamic direction cosine matrix.
    return [
        [x_a_w[0], y_a_w[0], z_a_w[0]],
        [x_a_w[1], y_a_w[1], z_a_w[1]],
        [x_a_w[2], y_a_w[2], z_a_w[2]],
    ]


def sigma_from_q_ba(q_ba: List[float]) -> float:
    """
    Extract actual bank angle from the body to aero attitude quaternion.
    """
    # Convert the quaternion into a body to aerodynamic rotation matrix.
    R_ba = quat_to_dcm(q_ba)

    # Read the body x axis expressed in the aerodynamic frame.
    x_b_in_a = [R_ba[0][0], R_ba[1][0], R_ba[2][0]]

    # Recover the roll angle about the aerodynamic z axis.
    sigma = math.atan2(x_b_in_a[1], x_b_in_a[0])

    # Apply the configured sign convention before wrapping the result.
    sigma_sign = getattr(constants, "SIGMA_SIGN", 1.0)
    return wrap_to_pi(sigma_sign * sigma)


def body_to_world_dcm(x_trans: List[float], q_ba: List[float]) -> List[List[float]]:
    """
    Return R_wb, the DCM mapping body frame vectors into world frame vectors.
    """
    # Build the world from aerodynamic rotation from the translational state.
    R_wa = aero_frame_dcm_from_state(x_trans)

    # Build the aerodynamic from body rotation from the quaternion.
    R_ba = quat_to_dcm(q_ba)

    # Compose both transforms to get world from body.
    return mat_mul_3x3(R_wa, R_ba)


def heat_shield_normal_world(x_trans: List[float], q_ba: List[float]) -> List[float]:
    """
    Return the heat shield normal in the world frame.
    """
    # Convert from body frame to world frame.
    R_wb = body_to_world_dcm(x_trans, q_ba)

    # Map the body z axis into world coordinates.
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
    # Build the body to world rotation.
    R_wb = body_to_world_dcm(x_trans, q_ba)

    # Get the capsule center of mass position in world coordinates.
    r_cm_w = position_world_from_state(x_trans)

    # Rotate the local thruster position offset into the world frame.
    r_thr_w_local = mat_vec_mul(R_wb, thruster_position_b_m)

    # Rotate and normalize the thruster direction into the world frame.
    d_thr_w = normalize(mat_vec_mul(R_wb, thruster_direction_b_unit))

    # Add the rotated local offset to the center of mass position.
    r_thr_w = [
        r_cm_w[0] + r_thr_w_local[0],
        r_cm_w[1] + r_thr_w_local[1],
        r_cm_w[2] + r_thr_w_local[2],
    ]

    # Return the world position and unit direction of the thruster.
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
    # Pull default inertia if one was not provided.
    if Izz_kgm2 is None:
        Izz_kgm2 = float(getattr(constants, "CAPSULE_IZZ_KGM2", 20000.0))

    # Pull the default fixed time step if one was not provided.
    if dt_s is None:
        dt_s = float(getattr(constants, "ENTRY_DT_S", 0.25))

    # Read the current roll rate from the body rate vector.
    roll_rate_old = float(att.omega_b_rad_s[2])

    # Read the body roll torque from the RCS wrench.
    tau_roll = float(tau_rcs_b_Nm[2])

    # Convert torque into roll angular acceleration.
    roll_accel = tau_roll / safe_nonzero(Izz_kgm2)

    # Integrate roll rate forward by one fixed step.
    roll_rate_new = roll_rate_old + dt_s * roll_accel

    # Integrate bank angle forward using the updated roll rate.
    sigma_new = wrap_to_pi(float(att.sigma_rel_rad) + dt_s * roll_rate_new)

    # Rebuild the attitude quaternion from the updated roll angle.
    q_ba_new = quat_from_axis_angle([0.0, 0.0, 1.0], sigma_new)

    # Return the updated roll only attitude state.
    return CapsuleAttitudeState(
        q_ba=q_ba_new,
        omega_b_rad_s=[0.0, 0.0, roll_rate_new],
        sigma_rel_rad=sigma_new,
    )


# Aerodynamic force model used by the float entry dynamics.

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
    Compute drag and lift components in the spherical three dimensional frame.
    """
    # Prefer actual bank when available, otherwise fall back to commanded bank or zero.
    sigma_used = float(
        sigma_actual_rad
        if sigma_actual_rad is not None
        else (sigma_cmd_rad if sigma_cmd_rad is not None else 0.0)
    )

    # Convert radius to altitude above the Earth surface.
    altitude_m = float(r) - float(constants.RADIUS_EARTH)

    # Prevent negative altitude from entering the atmosphere model.
    altitude_m = max(0.0, altitude_m)

    # Query the atmosphere model at the current altitude.
    atm = AtmosphereModel.US_Standard_ATM(altitude_m)

    # Use zero density if the atmosphere model returns a missing value.
    rho = atm["rho_kgm3"] if atm["rho_kgm3"] is not None else 0.0

    # Compute dynamic pressure.
    q = 0.5 * rho * V * V

    # Read aerodynamic coefficients and reference area from params.
    CD = _get_param(params, "CD_const", "CD")
    CL = _get_param(params, "CL_const", "CL")
    surface_area = _get_param(params, "ref_area_m2", "S")

    # Convert dynamic pressure into total drag and lift magnitudes.
    drag_mag = q * surface_area * CD
    lift_mag = q * surface_area * CL

    # Resolve the velocity direction in the local spherical basis.
    v_r = math.sin(gamma)
    v_theta = math.cos(gamma) * math.cos(chi)
    v_phi = math.cos(gamma) * math.sin(chi)

    # Drag points opposite the velocity direction.
    D_r = -drag_mag * v_r
    Dtheta = drag_mag * v_theta
    Dphi = -drag_mag * v_phi

    # e1 defines the lift direction associated with bank equal to zero.
    e1_r = math.cos(gamma)
    e1_theta = -math.sin(gamma) * math.cos(chi)
    e1_phi = -math.sin(gamma) * math.sin(chi)

    # e2 defines the lift direction ninety degrees away in bank.
    e2_r = 0.0
    e2_theta = math.sin(chi)
    e2_phi = -math.cos(chi)

    # Compute the bank rotation weights.
    cos_sig = math.cos(sigma_used)
    sin_sig = math.sin(sigma_used)

    # Resolve lift into spherical components after bank rotation.
    Lr = lift_mag * (e1_r * cos_sig + e2_r * sin_sig)
    Ltheta = lift_mag * (e1_theta * cos_sig + e2_theta * sin_sig)
    Lphi = lift_mag * (e1_phi * cos_sig + e2_phi * sin_sig)

    # Return force components together with density and dynamic pressure.
    return D_r, Dtheta, Dphi, Lr, Ltheta, Lphi, float(rho), float(q)


# Float equations of motion for the nominal translational state.

def eom_3d(
    t: float,
    X: List[float],
    params: Dict[str, Any],
    sigma_actual_rad: Optional[float] = None,
    sigma_cmd_rad: Optional[float] = None,
) -> Dict[str, float]:
    """
    Float spherical three dimensional equations of motion for atmospheric entry.
    """
    # Unpack the translational state.
    r, phi, lam, V, gamma, chi = [float(v) for v in X]

    # Read the vehicle mass and local gravity.
    m = _get_param(params, "mass_kg", "m")
    g = float(constants.gravity(r))

    # Compute aerodynamic force components at the current state.
    Dr, Dtheta, Dphi, Lr, Ltheta, Lphi, _, _ = aero_forces(
        r=r,
        V=V,
        gamma=gamma,
        chi=chi,
        params=params,
        sigma_actual_rad=sigma_actual_rad,
        sigma_cmd_rad=sigma_cmd_rad,
    )

    # Protect denominators that can become numerically small.
    cos_phi = safe_nonzero(math.cos(phi))
    cos_gamma = safe_nonzero(math.cos(gamma))
    V_safe = safe_nonzero(V)
    r_safe = safe_nonzero(r)

    # Radial rate comes from the vertical component of velocity.
    r_dot = V * math.sin(gamma)

    # Latitude rate depends on the northward component of motion.
    phi_dot = (V * math.cos(gamma) * math.sin(chi)) / r_safe

    # Longitude rate depends on the eastward component and the local cosine latitude factor.
    lam_dot = (V * math.cos(gamma) * math.cos(chi)) / (r_safe * cos_phi)

    # Speed rate includes aerodynamic acceleration along the velocity direction and gravity.
    V_dot = (Dr + Lr) / m - g * math.sin(gamma)

    # Flight path angle rate includes normal lift and spherical curvature effects.
    gamma_dot = (Ltheta / (m * V_safe)) + (V_safe / r_safe - g / V_safe) * math.cos(gamma)

    # Tangent of latitude appears in the heading dynamics.
    phi_tan = math.sin(phi) / cos_phi

    # Heading rate includes lateral lift and spherical transport terms.
    chi_dot = (Lphi / (m * V_safe * cos_gamma)) + (V_safe / r_safe) * math.sin(chi) * phi_tan

    # Return the derivatives in a named dictionary.
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
    # Evaluate the named derivative form first.
    dx = eom_3d(
        t=t,
        X=X,
        params=params,
        sigma_actual_rad=sigma_actual_rad,
        sigma_cmd_rad=sigma_cmd_rad,
    )

    # Return the derivatives in the same order as the state vector.
    return [
        dx["r_dot"],
        dx["phi_dot"],
        dx["lam_dot"],
        dx["V_dot"],
        dx["gamma_dot"],
        dx["chi_dot"],
    ]


# Control bridge helpers that connect the guidance stack to the nominal dynamics.

def make_control_step_fn(
    control: SigmaControlStack,
    aero_force_fn: Callable[..., Tuple[float, float, float, float, float, float, float, float]],
    params: Dict[str, Any],
    dt_s: Optional[float] = None,
) -> Callable[[float, List[float], float, float], Any]:
    """
    Build a helper with the form control_step_fn(t, x, sigma_actual_rad, roll_rate_rad_s).
    """
    # Use the nominal entry time step when one is not supplied explicitly.
    if dt_s is None:
        dt_s = float(getattr(constants, "ENTRY_DT_S", 0.25))

    # Read the mass once since it is passed into every control step.
    mass_kg = _get_param(params, "mass_kg", "m")

    def control_step_fn(
        t: float,
        x: List[float],
        sigma_actual_rad: float,
        roll_rate_rad_s: float,
    ):
        # Unpack the translational state for the controller facing state object.
        r, phi, lam, V, gamma, chi = [float(v) for v in x]

        # Build the controller input state using the current actual bank and roll rate.
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

        # Recompute aerodynamic forces so lift and drag magnitudes can be passed into control.
        Dr, Dtheta, Dphi, Lr, Ltheta, Lphi, _, _ = aero_force_fn(
            r=r,
            V=V,
            gamma=gamma,
            chi=chi,
            params=params,
            sigma_actual_rad=sigma_actual_rad,
        )

        # Convert component forces into total lift magnitude.
        lift_mag = math.sqrt(Lr * Lr + Ltheta * Ltheta + Lphi * Lphi)

        # Convert component forces into total drag magnitude.
        drag_mag = math.sqrt(Dr * Dr + Dtheta * Dtheta + Dphi * Dphi)

        # Advance the control stack by one step and return its output object.
        return control.step(
            t_s=float(t),
            dt_s=float(dt_s),
            state=state,
            lift_N=float(lift_mag),
            drag_N=float(drag_mag),
            mass_kg=float(mass_kg),
        )

    # Return the closure that matches the desired control step signature.
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
    # Reuse the richer control step helper internally.
    control_step_fn = make_control_step_fn(
        control=control,
        aero_force_fn=aero_force_fn,
        params=params,
        dt_s=dt_s,
    )

    def sigma_fn(t: float, x: List[float]) -> float:
        # Query the current actual bank and roll rate from the external provider.
        sigma_actual_rad, roll_rate_rad_s = sigma_actual_provider()

        # Step the controller so any internal state continues to evolve normally.
        _ = control_step_fn(
            t=float(t),
            x=list(x),
            sigma_actual_rad=float(sigma_actual_rad),
            roll_rate_rad_s=float(roll_rate_rad_s),
        )

        # Return the current actual bank angle for legacy call sites.
        return float(sigma_actual_rad)

    # Return the backward compatible adapter.
    return sigma_fn


# Closed loop milestone 1 step that advances control, RCS, attitude, and translation together.

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
    # Fill in default time step and roll inertia when they are not provided.
    if dt_s is None:
        dt_s = float(getattr(constants, "ENTRY_DT_S", 0.25))
    if Izz_kgm2 is None:
        Izz_kgm2 = float(getattr(constants, "CAPSULE_IZZ_KGM2", 20000.0))

    # Copy the incoming translational state so the result keeps both old and new values.
    x_old = [float(v) for v in x_trans]

    # Copy the incoming attitude state so the old state is preserved in the result object.
    att_old = CapsuleAttitudeState(
        q_ba=list(att.q_ba),
        omega_b_rad_s=list(att.omega_b_rad_s),
        sigma_rel_rad=float(att.sigma_rel_rad),
    )

    # Recover the actual bank angle represented by the old quaternion.
    sigma_actual_before_rad = float(sigma_from_q_ba(att.q_ba))

    # Read the current roll rate from the old attitude state.
    roll_rate_rad_s = float(att.omega_b_rad_s[2])

    # Run the controller using the current translational and attitude state.
    control_out = control_step_fn(
        float(t_s),
        list(x_old),
        float(sigma_actual_before_rad),
        float(roll_rate_rad_s),
    )

    # Advance the RCS roll channel using the commanded roll torque.
    roll_step = rcs.step_roll_channel(
        tau_roll_cmd_Nm=float(control_out.tau_roll_cmd_Nm),
        dt_s=float(dt_s),
    )

    # Propagate the roll attitude using the torque actually delivered by the RCS system.
    att_new = step_capsule_roll_state(
        att=att_old,
        tau_rcs_b_Nm=list(roll_step.wrench.torque_b_Nm),
        Izz_kgm2=float(Izz_kgm2),
        dt_s=float(dt_s),
    )

    # Recompute the actual bank angle after the roll state update.
    sigma_actual_after_rad = float(sigma_from_q_ba(att_new.q_ba))

    # Keep the scalar bank value in sync with the quaternion.
    att_new.sigma_rel_rad = sigma_actual_after_rad

    # Evaluate the translational state derivative using the updated actual bank angle.
    dx_trans = f_float(
        t=float(t_s),
        X=list(x_old),
        params=params,
        sigma_actual_rad=float(sigma_actual_after_rad),
    )

    # Advance the translational state with a forward Euler step.
    x_new = add_scaled_state(x_old, dx_trans, float(dt_s))

    # Return a complete packet containing both old and new states plus control and RCS details.
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