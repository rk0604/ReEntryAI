"""
Closed loop point mass dynamics and attitude helpers for the 3DOF simulator.

This module ties the translational state to a simple roll only attitude model
and the reaction control system. It provides:

  Small vector and quaternion helpers (normalize, dot, cross, matrix products,
  axis angle to quaternion, quaternion to direction cosine matrix).
  Frame helpers that move between the local level frame, the aerodynamic frame,
  the body frame, and the world frame.
  The roll only attitude state and its one step update, where a roll torque
  drives a roll rate and the bank angle.
  The aerodynamic force adapter and the live translational derivatives, both of
  which defer to the canonical equations in math_3d.
  step_closed_loop_milestone1, which advances guidance, the timed RCS, the roll
  attitude, and the translational state together for one outer step.

The translational state order is [r, phi, lambda, V, gamma, chi]. The bank angle
is recovered from the attitude quaternion, so the flown bank can differ from the
commanded bank while the actuator and RCS catch up.
"""

import math
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

import AtmosphereModel
import constants
from control import ReentryState, SigmaControlStack
from math_3d import nominal_aero_forces_from_state, nominal_eom_step
from ReactionControl import CapsuleRCSSystem, RCSWrench, ThrusterFireCommand


# Small math helpers used throughout the nominal simulation path

def wrap_to_pi(angle_rad: float) -> float:
    """Wrap an angle in radians into the range minus pi to pi."""
    # Keep angles inside the principal interval used everywhere else in the sim
    while angle_rad > math.pi:
        angle_rad -= 2.0 * math.pi

    while angle_rad < -math.pi:
        angle_rad += 2.0 * math.pi

    return angle_rad


def clamp(value: float, lower: float, upper: float) -> float:
    """Return value held inside the closed range [lower, upper]."""
    # Keep a scalar inside a closed interval
    return max(lower, min(upper, value))


def normalize(vec: List[float], eps: Optional[float] = None) -> List[float]:
    """
    Return the unit length version of a three component vector.

    Parameters:
        vec: a three element sequence.
        eps: smallest length treated as nonzero. Defaults to a small constant.

    Returns:
        The unit vector, or a zero vector when the input is too short to
        normalize safely.
    """
    # Pull the default norm threshold from constants when none is supplied
    if eps is None:
        eps = getattr(constants, "EPS_NORM", 1.0e-12)

    # Convert all components to float so downstream math stays consistent
    x = float(vec[0])
    y = float(vec[1])
    z = float(vec[2])

    # Euclidean norm of the vector
    n = math.sqrt(x * x + y * y + z * z)

    # Return a zero vector if the input is too small to normalize safely
    if n < eps:
        return [0.0, 0.0, 0.0]

    return [x / n, y / n, z / n]


def dot3(a: List[float], b: List[float]) -> float:
    """Return the dot product of two three component vectors."""
    # Standard three dimensional dot product
    return float(a[0]) * float(b[0]) + float(a[1]) * float(b[1]) + float(a[2]) * float(b[2])


def cross3(a: List[float], b: List[float]) -> List[float]:
    """Return the cross product of two three component vectors."""
    # Standard three dimensional cross product
    ax, ay, az = float(a[0]), float(a[1]), float(a[2])
    bx, by, bz = float(b[0]), float(b[1]), float(b[2])

    return [
        ay * bz - az * by,
        az * bx - ax * bz,
        ax * by - ay * bx,
    ]


def mat_vec_mul(M: List[List[float]], v: List[float]) -> List[float]:
    """Return the product of a three by three matrix and a three vector."""
    # Three by three matrix times three component vector
    return [
        float(M[0][0]) * float(v[0]) + float(M[0][1]) * float(v[1]) + float(M[0][2]) * float(v[2]),
        float(M[1][0]) * float(v[0]) + float(M[1][1]) * float(v[1]) + float(M[1][2]) * float(v[2]),
        float(M[2][0]) * float(v[0]) + float(M[2][1]) * float(v[1]) + float(M[2][2]) * float(v[2]),
    ]


def mat_mul_3x3(A: List[List[float]], B: List[List[float]]) -> List[List[float]]:
    """Return the product of two three by three matrices."""
    # Three by three matrix product
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
    """Nudge a value away from zero to keep divisions in the float path safe."""
    # Prevent division by zero in the float dynamics path
    if abs(value) >= eps:
        return value

    return eps if value >= 0.0 else -eps


def _get_param(params: Dict[str, Any], *keys: str, default: Optional[float] = None) -> float:
    """
    Return the first matching numeric parameter from a dictionary.

    Parameters:
        params:  the parameter dictionary.
        keys:    candidate key names tried in order.
        default: value returned when no key matches. When None, a missing key
                 raises KeyError.
    """
    # Return the first matching numeric parameter from a dictionary
    for key in keys:
        if key in params:
            return float(params[key])

    if default is None:
        raise KeyError(f"Missing required parameter Tried keys={keys}")

    return float(default)


def add_scaled_state(x: List[float], dx: List[float], dt_s: float) -> List[float]:
    """
    Take one forward Euler step on a flat state vector.

    Parameters:
        x:    the current state vector.
        dx:   the derivative vector.
        dt_s: the time step in seconds.

    Returns:
        The vector x plus dt_s times dx, element by element.
    """
    # Forward Euler update for a flat state vector
    return [float(x[i]) + float(dt_s) * float(dx[i]) for i in range(len(x))]


@dataclass
class CapsuleAttitudeState:
    """
    Roll only attitude of the capsule.

    Fields:
        q_ba:          quaternion that maps body vectors into the aerodynamic
                       frame.
        omega_b_rad_s: body angular rate relative to the aerodynamic frame.
        sigma_rel_rad: the bank angle the attitude currently represents.
    """
    # q_ba maps body frame vectors into aerodynamic frame vectors
    q_ba: List[float]

    # omega_b_rad_s stores body angular rate relative to the aerodynamic frame
    omega_b_rad_s: List[float]

    # sigma_rel_rad stores the actual bank angle represented by the attitude
    sigma_rel_rad: float


@dataclass
class AggregatedRollChannelStepResult:
    """
    Roll channel telemetry aggregated across the internal RCS sub steps of one
    outer step. Holds the commanded torque, the available capacity, the average
    duty, the time averaged wrench, and the firing counts and backlog.
    """
    # Outer step roll torque command from the controller
    tau_roll_cmd_Nm: float

    # Magnitude of available roll torque for the active direction
    tau_roll_capacity_Nm: float

    # Requested average duty for the outer step
    requested_duty: float

    # True if any internal RCS substep fired
    fired_this_step: bool

    # Average duty by thruster across the full outer step
    fire_cmd: ThrusterFireCommand

    # Time averaged body wrench over the full outer step
    wrench: RCSWrench

    # Remaining undelivered positive roll demand in seconds
    roll_pos_backlog_s: float

    # Remaining undelivered negative roll demand in seconds
    roll_neg_backlog_s: float

    # Positive roll channel firing state after the last internal substep
    roll_pos_is_on: bool

    # Negative roll channel firing state after the last internal substep
    roll_neg_is_on: bool

    # Number of internal RCS substeps inside this outer step
    num_internal_steps: int

    # Number of internal substeps that actually fired any roll thruster
    num_fired_internal_steps: int

    # Full list of internal roll packets for debugging if needed
    substeps: List[Any]


@dataclass
class Milestone1StepResult:
    """
    Full record of one closed loop outer step. Holds the time, the old and new
    translational and attitude states, the bank before and after, the control
    output, the aggregated roll telemetry, and the translational derivative.
    """
    # Outer time at the start of the nominal step
    t_s: float

    # Translational state before the step
    x_old: List[float]

    # Translational state after the step
    x_new: List[float]

    # Attitude state before the step
    att_old: CapsuleAttitudeState

    # Attitude state after the step
    att_new: CapsuleAttitudeState

    # Actual bank before the outer step
    sigma_actual_before_rad: float

    # Actual bank after the outer step
    sigma_actual_after_rad: float

    # High level control output computed once for the outer step
    control_out: Any

    # Time averaged roll realization packet across the outer step
    roll_step: Any

    # Translational derivative used for the outer step Euler update
    dx_trans: List[float]


def make_initial_capsule_attitude(sigma0_rad: float = 0.0) -> CapsuleAttitudeState:
    """
    Build the starting attitude as a pure roll about the aerodynamic z axis.

    Parameters:
        sigma0_rad: initial bank angle in radians.

    Returns:
        A CapsuleAttitudeState at the given bank with zero angular rate.
    """
    # Initialize the capsule in a pure roll state relative to the aerodynamic frame
    q_ba = quat_from_axis_angle([0.0, 0.0, 1.0], sigma0_rad)

    return CapsuleAttitudeState(
        q_ba=q_ba,
        omega_b_rad_s=[0.0, 0.0, 0.0],
        sigma_rel_rad=float(wrap_to_pi(sigma0_rad)),
    )


# Quaternion helpers used by the roll attitude model

def quat_normalize(q: List[float]) -> List[float]:
    """Return the unit quaternion, falling back to identity when q is too small."""
    # Convert quaternion components to floats
    w, x, y, z = [float(v) for v in q]

    # Euclidean magnitude of the quaternion
    n = math.sqrt(w * w + x * x + y * y + z * z)

    # Fall back to identity if the quaternion is too small
    if n < getattr(constants, "EPS_NORM", 1.0e-12):
        return [1.0, 0.0, 0.0, 0.0]

    return [w / n, x / n, y / n, z / n]


def quat_from_axis_angle(axis: List[float], angle_rad: float) -> List[float]:
    """
    Build a unit quaternion for a rotation about an axis.

    Parameters:
        axis:      the rotation axis, normalized internally.
        angle_rad: the rotation angle in radians.

    Returns:
        The rotation quaternion in [w, x, y, z] order.
    """
    # Normalize the rotation axis first
    ax = normalize(axis)

    # Axis angle to quaternion uses the half angle form
    half = 0.5 * float(angle_rad)
    s = math.sin(half)

    return quat_normalize([
        math.cos(half),
        float(ax[0]) * s,
        float(ax[1]) * s,
        float(ax[2]) * s,
    ])


def quat_to_dcm(q: List[float]) -> List[List[float]]:
    """Return the direction cosine matrix for a quaternion, normalized first."""
    # Normalize before converting so the DCM stays orthonormal
    w, x, y, z = quat_normalize(q)

    return [
        [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w),       2.0 * (x * z + y * w)],
        [2.0 * (x * y + z * w),       1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w)],
        [2.0 * (x * z - y * w),       2.0 * (y * z + x * w),       1.0 - 2.0 * (x * x + y * y)],
    ]


# Frame helpers that move between local level aerodynamic body and world coordinates

def local_level_frame(phi_rad: float, lam_rad: float) -> Tuple[List[float], List[float], List[float]]:
    """
    Return the local north, east, and up unit vectors in world coordinates.

    Parameters:
        phi_rad: latitude in radians.
        lam_rad: longitude in radians.

    Returns:
        A tuple (north, east, up) of unit vectors.
    """
    # Local north direction in world coordinates
    north = [
        -math.sin(phi_rad) * math.cos(lam_rad),
        -math.sin(phi_rad) * math.sin(lam_rad),
        math.cos(phi_rad),
    ]

    # Local east direction in world coordinates
    east = [
        -math.sin(lam_rad),
        math.cos(lam_rad),
        0.0,
    ]

    # Local up direction points radially away from Earth center
    up = [
        math.cos(phi_rad) * math.cos(lam_rad),
        math.cos(phi_rad) * math.sin(lam_rad),
        math.sin(phi_rad),
    ]

    return normalize(north), normalize(east), normalize(up)


def position_world_from_state(x: List[float]) -> List[float]:
    """Return the world position vector for a translational state, radius times up."""
    # Convert spherical position into a world position vector
    r = float(x[0])
    phi = float(x[1])
    lam = float(x[2])

    _, _, up = local_level_frame(phi, lam)

    return [r * up[0], r * up[1], r * up[2]]


def velocity_hat_from_state(x: List[float]) -> List[float]:
    """Return the unit velocity direction in world coordinates from gamma and chi."""
    # Read the translational state
    _, phi, lam, _, gamma, chi = [float(v) for v in x]

    # Build the local basis at the current location
    north, east, up = local_level_frame(phi, lam)

    # Resolve the velocity direction using gamma and chi
    v_hat = [
        math.cos(gamma) * math.sin(chi) * north[0] + math.cos(gamma) * math.cos(chi) * east[0] + math.sin(gamma) * up[0],
        math.cos(gamma) * math.sin(chi) * north[1] + math.cos(gamma) * math.cos(chi) * east[1] + math.sin(gamma) * up[1],
        math.cos(gamma) * math.sin(chi) * north[2] + math.cos(gamma) * math.cos(chi) * east[2] + math.sin(gamma) * up[2],
    ]

    return normalize(v_hat)


def aero_frame_dcm_from_state(x: List[float]) -> List[List[float]]:
    """
    Return the aerodynamic frame basis in world coordinates as a matrix.

    The aerodynamic z axis points opposite the velocity, the x axis is the
    projection of local up into the plane normal to z, with north used as a
    fallback when up is nearly aligned with z, and y completes a right handed
    frame.

    Parameters:
        x: the translational state.

    Returns:
        A three by three matrix whose columns are the aerodynamic x, y, and z
        axes expressed in world coordinates.
    """
    # Build the aerodynamic frame basis in world coordinates
    _, phi, lam, _, _, _ = [float(v) for v in x]

    north, _, up = local_level_frame(phi, lam)
    v_hat = velocity_hat_from_state(x)

    # Aerodynamic z points opposite the velocity direction
    z_a_w = normalize([-v_hat[0], -v_hat[1], -v_hat[2]])

    # Aerodynamic x is the projection of local up into the plane normal to aerodynamic z
    x_a_w = [
        up[0] - dot3(up, z_a_w) * z_a_w[0],
        up[1] - dot3(up, z_a_w) * z_a_w[1],
        up[2] - dot3(up, z_a_w) * z_a_w[2],
    ]

    # Use north as a fallback if local up is nearly aligned with aerodynamic z
    if math.sqrt(dot3(x_a_w, x_a_w)) < getattr(constants, "EPS_NORM", 1.0e-12):
        x_a_w = [
            north[0] - dot3(north, z_a_w) * z_a_w[0],
            north[1] - dot3(north, z_a_w) * z_a_w[1],
            north[2] - dot3(north, z_a_w) * z_a_w[2],
        ]

    x_a_w = normalize(x_a_w)

    # Aerodynamic y completes the right handed frame
    y_a_w = normalize(cross3(z_a_w, x_a_w))

    return [
        [x_a_w[0], y_a_w[0], z_a_w[0]],
        [x_a_w[1], y_a_w[1], z_a_w[1]],
        [x_a_w[2], y_a_w[2], z_a_w[2]],
    ]


def sigma_from_q_ba(q_ba: List[float]) -> float:
    """Recover the bank angle in radians from the body to aerodynamic quaternion."""
    # Convert body to aerodynamic attitude into the actual bank angle
    R_ba = quat_to_dcm(q_ba)
    x_b_in_a = [R_ba[0][0], R_ba[1][0], R_ba[2][0]]

    sigma = math.atan2(x_b_in_a[1], x_b_in_a[0])

    sigma_sign = getattr(constants, "SIGMA_SIGN", 1.0)
    return wrap_to_pi(sigma_sign * sigma)


def body_to_world_dcm(x_trans: List[float], q_ba: List[float]) -> List[List[float]]:
    """Return the body to world rotation, world from aero composed with aero from body."""
    # Compose world from aerodynamic with aerodynamic from body
    R_wa = aero_frame_dcm_from_state(x_trans)
    R_ba = quat_to_dcm(q_ba)
    return mat_mul_3x3(R_wa, R_ba)


def heat_shield_normal_world(x_trans: List[float], q_ba: List[float]) -> List[float]:
    """Return the heat shield normal, body z, expressed in world coordinates."""
    # Map body z into world coordinates so the current shield normal is available
    R_wb = body_to_world_dcm(x_trans, q_ba)
    return mat_vec_mul(R_wb, [0.0, 0.0, 1.0])


def thruster_world_pose(
    x_trans: List[float],
    q_ba: List[float],
    thruster_position_b_m: List[float],
    thruster_direction_b_unit: List[float],
) -> Tuple[List[float], List[float]]:
    """
    Convert a body fixed thruster pose into the world frame.

    Parameters:
        x_trans:                   the translational state.
        q_ba:                      the body to aerodynamic quaternion.
        thruster_position_b_m:     thruster position in body coordinates.
        thruster_direction_b_unit: thruster direction in body coordinates.

    Returns:
        A tuple (position_world, direction_world). The position adds the world
        offset to the center of mass position, and the direction is a unit
        vector.
    """
    # Convert a body fixed thruster pose into the world frame
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
    Advance the roll only attitude model by one time step.

    The roll torque about body z makes an angular acceleration through the roll
    inertia. The roll rate is integrated first, then the bank angle, then the
    quaternion is rebuilt from the new bank.

    Parameters:
        att:          the current attitude state.
        tau_rcs_b_Nm: the body torque vector. Only the z component is used.
        Izz_kgm2:     roll inertia. Defaults to the capsule value.
        dt_s:         time step in seconds. Defaults to the entry step.

    Returns:
        The next attitude state.
    """
    # Advance the roll only attitude model for one time step
    if Izz_kgm2 is None:
        Izz_kgm2 = float(getattr(constants, "CAPSULE_IZZ_KGM2", 20000.0))

    if dt_s is None:
        dt_s = float(getattr(constants, "ENTRY_DT_S", 0.25))

    # Current roll rate is the body z rate in this milestone model
    roll_rate_old = float(att.omega_b_rad_s[2])

    # Body z torque is the roll control torque
    tau_roll = float(tau_rcs_b_Nm[2])

    # Angular acceleration follows rigid body dynamics about body z
    roll_accel = tau_roll / safe_nonzero(Izz_kgm2)

    # Integrate roll rate first
    roll_rate_new = roll_rate_old + dt_s * roll_accel

    # Then integrate bank angle using the updated roll rate
    sigma_new = wrap_to_pi(float(att.sigma_rel_rad) + dt_s * roll_rate_new)

    # Rebuild the quaternion from the new roll angle
    q_ba_new = quat_from_axis_angle([0.0, 0.0, 1.0], sigma_new)

    return CapsuleAttitudeState(
        q_ba=q_ba_new,
        omega_b_rad_s=[0.0, 0.0, roll_rate_new],
        sigma_rel_rad=sigma_new,
    )


# Aerodynamic force model used by the float entry dynamics

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
    Return nominal aerodynamic force components using the canonical model

    This helper keeps the legacy return order used by the notebook
    The actual force decomposition is delegated to math_3d so the live path
    and the predictor path share one translational model

    Input
    r radius in meters
    V speed in meters per second
    gamma flight path angle in radians
    chi heading angle in radians measured from east toward north
    params aerodynamic parameter dictionary
    sigma_actual_rad realized bank angle in radians when available
    sigma_cmd_rad commanded bank angle used only as a fallback

    Output
    D_r
    Dtheta
    Dphi
    Lr
    Ltheta
    Lphi
    rho
    q
    """
    # Realized bank must drive the flown translational dynamics
    sigma_used = float(
        sigma_actual_rad
        if sigma_actual_rad is not None
        else (sigma_cmd_rad if sigma_cmd_rad is not None else 0.0)
    )

    # The canonical force model only needs radius speed gamma chi and sigma
    x_state = [
        float(r),
        0.0,
        0.0,
        float(V),
        float(gamma),
        float(chi),
    ]

    # Pull the canonical decomposition from math_3d
    aero = nominal_aero_forces_from_state(
        x=x_state,
        sigma_rad=float(sigma_used),
        params=params,
    )

    # Return the legacy component order so existing notebook code stays valid
    return (
        float(aero["D_r_N"]),
        float(aero["Dtheta"]),
        float(aero["Dphi"]),
        float(aero["L_r_N"]),
        float(aero["Ltheta"]),
        float(aero["Lphi"]),
        float(aero["rho_kgm3"]),
        float(aero["q_pa"]),
    )

def eom_3d(
    t: float,
    X: List[float],
    params: Dict[str, Any],
    sigma_actual_rad: Optional[float] = None,
    sigma_cmd_rad: Optional[float] = None,
) -> Dict[str, float]:
    """
    Return live translational derivatives from the canonical nominal model

    This removes the stale duplicate live dynamics
    The live nominal simulation now uses the same equations as the predictor
    rollout in math_3d

    Input
    t time in seconds
    X six component translational state
    params aerodynamic and mass parameter dictionary
    sigma_actual_rad realized bank angle in radians when available
    sigma_cmd_rad commanded bank angle used only as a fallback

    Output
    Dictionary with r_dot phi_dot lam_dot V_dot gamma_dot chi_dot
    """
    # Time is kept only for signature compatibility
    del t

    # Realized bank must drive the flown translational dynamics
    sigma_used = float(
        sigma_actual_rad
        if sigma_actual_rad is not None
        else (sigma_cmd_rad if sigma_cmd_rad is not None else 0.0)
    )

    # Use the canonical nominal equations from math_3d
    nominal_step = nominal_eom_step(
        x=[float(v) for v in X],
        sigma_rad=float(sigma_used),
        params=params,
        dt_s=1.0,
    )

    r_dot, phi_dot, lam_dot, V_dot, gamma_dot, chi_dot = nominal_step["x_dot"]

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
    Return the translational derivative as an ordered list.

    A thin adapter around eom_3d that returns [r_dot, phi_dot, lam_dot, V_dot,
    gamma_dot, chi_dot] in state order. Parameters match eom_3d.
    """
    # Ordered derivative adapter for the translational state
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


# Control bridge helpers that connect the guidance stack to the nominal dynamics

def make_control_step_fn(
    control: SigmaControlStack,
    aero_force_fn: Callable[..., Tuple[float, float, float, float, float, float, float, float]],
    params: Dict[str, Any],
    dt_s: Optional[float] = None,
) -> Callable[[float, List[float], float, float], Any]:
    """
    Build a closure that runs one control stack step.

    Parameters:
        control:       the guidance and control stack.
        aero_force_fn: function that returns the aero force components.
        params:        vehicle and aero parameters.
        dt_s:          control step in seconds. Defaults to the entry step.

    Returns:
        A function control_step_fn(t, x, sigma_actual_rad, roll_rate_rad_s) that
        builds the control facing state, computes lift and drag, and steps the
        control stack.
    """
    # This helper presents the controller with the current translational state actual bank and roll rate
    if dt_s is None:
        dt_s = float(getattr(constants, "ENTRY_DT_S", 0.25))

    mass_kg = _get_param(params, "mass_kg", "m")

    def control_step_fn(
        t: float,
        x: List[float],
        sigma_actual_rad: float,
        roll_rate_rad_s: float,
    ):
        """Build the control facing state, compute lift and drag, and step the stack."""
        # Build the controller facing state object from the translational state
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

        # Compute total lift and drag magnitudes for the control stack
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
    Build a backward compatible bank function for older call sites.

    Parameters:
        control:               the guidance and control stack.
        aero_force_fn:         function that returns the aero force components.
        params:                vehicle and aero parameters.
        sigma_actual_provider: callable returning the current bank and roll rate.
        dt_s:                  control step in seconds.

    Returns:
        A function sigma_fn(t, x) that steps the controller so its internal
        state advances, then returns the realized bank.
    """
    # Backward compatible adapter for older call sites that still expect sigma_fn
    control_step_fn = make_control_step_fn(
        control=control,
        aero_force_fn=aero_force_fn,
        params=params,
        dt_s=dt_s,
    )

    def sigma_fn(t: float, x: List[float]) -> float:
        """Step the controller from the provided bank and rate, then return the bank."""
        # Query the current realized bank and roll rate from the external provider
        sigma_actual_rad, roll_rate_rad_s = sigma_actual_provider()

        # Step the controller so internal control state still evolves
        _ = control_step_fn(
            t=float(t),
            x=list(x),
            sigma_actual_rad=float(sigma_actual_rad),
            roll_rate_rad_s=float(roll_rate_rad_s),
        )

        # Return the realized bank for older code paths
        return float(sigma_actual_rad)

    return sigma_fn


def _aggregate_roll_substeps(roll_substeps: List[Any]) -> AggregatedRollChannelStepResult:
    """
    Combine the internal RCS sub step results into one outer step packet.

    Parameters:
        roll_substeps: the per sub step roll channel results. Must be non empty.

    Returns:
        An AggregatedRollChannelStepResult with the time averaged wrench, the
        average duty per thruster, the fired counts, and the final channel
        timing telemetry.
    """
    # At least one internal substep must exist
    if not roll_substeps:
        raise ValueError("roll_substeps cannot be empty")

    num_internal_steps = len(roll_substeps)
    num_fired_internal_steps = 0

    # Sum wrench across internal substeps so the average outer step wrench can be formed
    total_force_b = np.zeros(3, dtype=float)
    total_torque_b = np.zeros(3, dtype=float)

    # Accumulate average duty by thruster across the outer step
    duty_totals: Dict[str, float] = {}

    for substep in roll_substeps:
        total_force_b += np.asarray(substep.wrench.force_b_N, dtype=float)
        total_torque_b += np.asarray(substep.wrench.torque_b_Nm, dtype=float)

        if bool(substep.fired_this_step):
            num_fired_internal_steps += 1

        for name, duty in substep.fire_cmd.duty_by_name.items():
            duty_totals[name] = duty_totals.get(name, 0.0) + float(duty)

    # Convert summed duty into average duty over the whole outer step
    avg_duty_by_name = {
        name: total / float(num_internal_steps)
        for name, total in duty_totals.items()
        if total > 0.0
    }

    avg_fire_cmd = ThrusterFireCommand(duty_by_name=avg_duty_by_name)

    # Time averaged outer step wrench
    avg_wrench = RCSWrench(
        force_b_N=total_force_b / float(num_internal_steps),
        torque_b_Nm=total_torque_b / float(num_internal_steps),
    )

    # Keep some channel telemetry from the final internal substep because it is the latest true state
    last = roll_substeps[-1]
    first = roll_substeps[0]

    return AggregatedRollChannelStepResult(
        tau_roll_cmd_Nm=float(first.tau_roll_cmd_Nm),
        tau_roll_capacity_Nm=float(first.tau_roll_capacity_Nm),
        requested_duty=float(first.requested_duty),
        fired_this_step=bool(num_fired_internal_steps > 0),
        fire_cmd=avg_fire_cmd,
        wrench=avg_wrench,
        roll_pos_backlog_s=float(last.roll_pos_backlog_s),
        roll_neg_backlog_s=float(last.roll_neg_backlog_s),
        roll_pos_is_on=bool(last.roll_pos_is_on),
        roll_neg_is_on=bool(last.roll_neg_is_on),
        num_internal_steps=int(num_internal_steps),
        num_fired_internal_steps=int(num_fired_internal_steps),
        substeps=list(roll_substeps),
    )


# Closed loop milestone one step that advances control RCS attitude and translation together

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
    Advance guidance, the RCS, the roll attitude, and the translation together.

    One outer step runs the high level controller once, then steps the timed RCS
    and the roll attitude on the smaller internal sub step, and finally advances
    the translational state with one forward Euler step using the realized bank.

    Parameters:
        t_s:             current time in seconds.
        x_trans:         the six element translational state.
        att:             the current attitude state.
        params:          vehicle and aero parameters.
        control_step_fn: closure that runs one control stack step.
        rcs:             the reaction control system.
        dt_s:            outer step in seconds. Defaults to the entry step.
        Izz_kgm2:        roll inertia. Defaults to the capsule value.

    Returns:
        A Milestone1StepResult with the old and new states, the bank before and
        after, the control output, the aggregated roll telemetry, and the
        translational derivative used for the step.
    """
    # Fill in default outer step size and roll inertia when not supplied
    if dt_s is None:
        dt_s = float(getattr(constants, "ENTRY_DT_S", 0.25))

    if Izz_kgm2 is None:
        Izz_kgm2 = float(getattr(constants, "CAPSULE_IZZ_KGM2", 20000.0))

    # Internal RCS and roll propagation runs on a smaller substep than translational physics
    internal_dt_s = float(getattr(constants, "RCS_INTERNAL_DT_S", dt_s))
    internal_steps = int(round(float(dt_s) / internal_dt_s))

    # The outer step must be an exact multiple of the internal actuator step
    if internal_steps <= 0:
        raise ValueError("internal_steps must be positive")

    if abs(internal_steps * internal_dt_s - float(dt_s)) > 1.0e-12:
        raise ValueError("dt_s must be an exact multiple of RCS_INTERNAL_DT_S")

    # Copy translational state so both old and new values are preserved in the result packet
    x_old = [float(v) for v in x_trans]

    # Copy attitude state so the old state is preserved
    att_old = CapsuleAttitudeState(
        q_ba=list(att.q_ba),
        omega_b_rad_s=list(att.omega_b_rad_s),
        sigma_rel_rad=float(att.sigma_rel_rad),
    )

    # Recover the current actual bank represented by the old quaternion
    sigma_actual_before_rad = float(sigma_from_q_ba(att.q_ba))

    # Current roll rate is the body z angular rate in this roll only model
    roll_rate_rad_s = float(att.omega_b_rad_s[2])

    # High level controller is evaluated once per outer step
    control_out = control_step_fn(
        float(t_s),
        list(x_old),
        float(sigma_actual_before_rad),
        float(roll_rate_rad_s),
    )

    # Timed RCS firing and rigid body roll propagation happen on the faster internal substep
    att_running = CapsuleAttitudeState(
        q_ba=list(att_old.q_ba),
        omega_b_rad_s=list(att_old.omega_b_rad_s),
        sigma_rel_rad=float(att_old.sigma_rel_rad),
    )

    roll_substeps: List[Any] = []

    for substep_idx in range(internal_steps):
        # Use the same outer torque command across all internal actuator substeps
        roll_step_i = rcs.step_roll_channel(
            tau_roll_cmd_Nm=float(control_out.tau_roll_cmd_Nm),
            dt_s=float(internal_dt_s),
        )

        roll_substeps.append(roll_step_i)

        # Integrate actual roll motion using the real torque delivered in this internal substep
        att_running = step_capsule_roll_state(
            att=att_running,
            tau_rcs_b_Nm=list(roll_step_i.wrench.torque_b_Nm),
            Izz_kgm2=float(Izz_kgm2),
            dt_s=float(internal_dt_s),
        )

    # Final attitude after all internal substeps is the realized outer step attitude
    att_new = att_running

    # Recompute the actual bank after the internal actuator and attitude update
    sigma_actual_after_rad = float(sigma_from_q_ba(att_new.q_ba))

    # Keep the scalar bank field synchronized with the quaternion
    att_new.sigma_rel_rad = sigma_actual_after_rad

    # Aggregate the internal actuator packets into one outer step packet for logging
    roll_step = _aggregate_roll_substeps(roll_substeps)

    # Translational dynamics still advance on the outer Euler step using the realized final bank
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
