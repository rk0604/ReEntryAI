"""
Version 1 of the simulation: a planar (2D) entry model.

This is the early two dimensional model, kept for reference. It assumes a
spherical, non rotating Earth and tracks four states: altitude, speed, flight
path angle, and downrange distance. The production simulator uses the three
dimensional model in math_3d.py instead. This file is retained because the
planar equations are simple and useful for checking the basic physics.
"""

import constants
import AtmosphereModel
# import StochasticEntrySim.math_3d as math_3d
import math
from typing import Callable

state = {
    "alt": 0.0,     # geometric altitude
    "speed": 0.0,   # inertial speed
    "gamma": 0.0,   # flight path angle
    "s":0.0         # downrange distance
}

def altitude_rate(speed:float, gamma:float) -> float:
    """
    Return the rate of change of altitude.

    Parameters:
        speed: speed in meters per second.
        gamma: flight path angle in radians, positive when climbing.

    Returns:
        The altitude rate in meters per second, which is speed times sin(gamma).
    """
    h_dot = speed * (math.sin(gamma))
    return h_dot

def velocity_rate(mass:float, drag_mag: float, gamma: float, alt:float) -> float:
    """
    Return the rate of change of speed.

    Speed drops from drag and from the component of gravity along the flight
    path.

    Parameters:
        mass:     vehicle mass in kilograms.
        drag_mag: drag force magnitude in newtons.
        gamma:    flight path angle in radians.
        alt:      geometric altitude in meters, used to find local gravity.

    Returns:
        The speed rate in meters per second squared.
    """
    space_craft_radius = alt + constants.RADIUS_EARTH
    drag_mass_ratio = -1*(drag_mag/mass)
    gravity_effect = (constants.gravity(space_craft_radius)) * (math.sin(gamma))

    v_dot = drag_mass_ratio - gravity_effect
    return v_dot

def flight_path_rate(lift_in_plane: float,
                     sigma: float,
                     speed: float,
                     mass: float,
                     alt: float,
                     gamma: float) -> float:
    """
    Return the rate of change of the flight path angle.

    The angle turns from the in plane lift and from the balance of centrifugal
    curvature against gravity.

    Parameters:
        lift_in_plane: lift component in the plane of motion, in newtons.
        sigma:         bank angle in radians (kept for interface symmetry).
        speed:         speed in meters per second.
        mass:          vehicle mass in kilograms.
        alt:           geometric altitude in meters.
        gamma:         flight path angle in radians.

    Returns:
        The flight path angle rate in radians per second.
    """
    space_craft_radius = alt + constants.RADIUS_EARTH
    lift_mass_ratio = lift_in_plane / (mass * speed)

    # Centrifugal term minus the gravity term, the source of skip behavior.
    curvature_minus_grav = (speed / space_craft_radius) - (
        constants.gravity(space_craft_radius) / speed
    )

    gamma_dot = lift_mass_ratio + curvature_minus_grav * math.cos(gamma)
    return gamma_dot

def downrange_distance_rate(speed:float, gamma: float) -> float:
    """
    Return the rate of change of downrange distance.

    Parameters:
        speed: speed in meters per second.
        gamma: flight path angle in radians.

    Returns:
        The downrange rate in meters per second, which is speed times
        cos(gamma).
    """
    s_dot = speed * (math.cos(gamma))
    return s_dot

def dynamic_pressure(rho: float, speed:float) -> float:
    """
    Return the dynamic pressure.

    Parameters:
        rho:   atmospheric density in kilograms per cubic meter.
        speed: speed in meters per second.

    Returns:
        Dynamic pressure in pascals, which is one half rho times speed squared.
    """
    v_2 = speed ** 2
    dynamic_p = 0.5 * (rho * v_2)
    return dynamic_p

def aerodynamic_forces(q: float, ref_area: float, C_D:float, C_L:float) -> dict:
    """
    Return the drag and lift magnitudes from dynamic pressure and coefficients.

    Parameters:
        q:        dynamic pressure in pascals.
        ref_area: reference area in square meters.
        C_D:      drag coefficient.
        C_L:      lift coefficient.

    Returns:
        A dictionary with keys drag_mag and lift_mag, both in newtons.
    """
    # Drag_mag = q * S * C_D
    drag = q * ref_area * C_D
    # Lift_mag = q*S*C_L
    lift = q * ref_area * C_L
    out = {"drag_mag": drag, "lift_mag": lift}
    return out

# how does bank angle affect lift? helper
"""
Bank angle interpretation for in-plane lift component:
    L_in_plane = L * cos(sigma)

Then:
- If sigma = 0 rad:
    cos(sigma) = 1  =>  L_in_plane = +L
    (all lift is used in the 2D plane of motion - going up)

- If sigma = pi/2 rad:
    cos(sigma) = 0  =>  L_in_plane = 0
    (no lift component in the 2D plane; lift is entirely out-of-plane - parallel to the ground)

- If sigma = pi rad:
    cos(sigma) = -1 =>  L_in_plane = -L
    (lift is flipped, acting downward in the 2D plane - downward)
"""

def banked_lift_prj(lift_mag:float, sigma:float) -> float:
    """
    Project the lift magnitude onto the plane of motion for a given bank.

    Parameters:
        lift_mag: total lift magnitude in newtons.
        sigma:    bank angle in radians.

    Returns:
        The in plane lift in newtons, which is lift_mag times cos(sigma). See
        the note above for the meaning at zero, ninety, and one hundred eighty
        degrees of bank.
    """
    Lift_in_plane = lift_mag * (math.cos(sigma))
    return Lift_in_plane

# main planar dynamics function --------------------------------------------------------------------------- b1fead

# Given time, current state, control (bank angle), and parameters, compute the four state derivatives
state_2d = {
    "alt": 120000.0,             # h0 [m], from init_altitude_m
    "speed": 7700.0,             # V0 [m/s], from init_speed_mps
    "gamma": math.radians(-5.0), # rad, from init_fpa_deg
    "s": 0.0,                    # initial downrange [m]
}

vehicle_params = {
    "mass_kg":     5000.0,
    "ref_area_m2": 10.0,
    "CL_const":    0.30,
    "CD_const":    1.00,
    "nose_radius_m": 1.0,     # not needed for dynamics, but available
}

env_params = {
    "mu_m3s2":   3.986004418e14,
    "R_E_m":     6371000.0,
    "rho0_kgm3": 1.225,
    "H_m":       7200.0,
}

bank_profile = "constant"     # from bank_profile
sigma_const_deg = 0.0         # from sigma_const_deg

# Convert degrees to radians before calling eom_2d.
# sigma_rad = math.radians(sigma_const_deg)
sigma = 0.0                   # placeholder here, but in practice: radians

def eom_2d(t:int, state:dict, sigma: float, vehicle_params:dict, env_params:float) -> dict:
    """
    Compute the four planar state derivatives for one instant.

    Pulls the current state, queries the atmosphere, builds the aerodynamic
    forces, projects lift into the plane, and returns the rates for altitude,
    speed, flight path angle, and downrange distance.

    Parameters:
        t:              time in seconds (unused, kept for integrator symmetry).
        state:          dict with keys alt, speed, gamma, s.
        sigma:          bank angle in radians.
        vehicle_params: dict with mass_kg, ref_area_m2, CL_const, CD_const.
        env_params:     environment parameters (unused in this body).

    Returns:
        A dict with keys alt_rate, vel_rate, gamma_rate, and down_range_rate.
    """
    # fetch state vars
    alt = state["alt"]
    speed = state["speed"]
    gamma = state["gamma"]
    down_range_distance = state["s"]

    mass = vehicle_params["mass_kg"]
    ref_area_m2 = vehicle_params["ref_area_m2"]
    CL = vehicle_params["CL_const"]
    CD = vehicle_params["CD_const"]

    # atmosphere calc.
    atmosphere_dict = AtmosphereModel.US_Standard_ATM(alt) # h is the altitude
    rho = atmosphere_dict["rho_kgm3"]
    temp = atmosphere_dict["T_K"]
    pressure = atmosphere_dict["p_Pa"]

    # aero forces
    dynamic_p_q = dynamic_pressure(rho, speed)
    aero_dict = aerodynamic_forces(dynamic_p_q, ref_area_m2, CD, CL)
    drag_mag = aero_dict["drag_mag"]
    lift_mag = aero_dict["lift_mag"]

    lift_in_plane = banked_lift_prj(lift_mag, sigma)
    # compute gravity at current position
    space_craft_radius = alt + constants.RADIUS_EARTH
    gravity_curr = constants.gravity(space_craft_radius)

    # rate functions
    alt_rate = altitude_rate(speed, gamma)
    vel_rate = velocity_rate(mass, drag_mag, gamma, alt)
    gamma_rate = flight_path_rate(lift_in_plane, sigma, speed, mass, alt, gamma)
    down_range_rate = downrange_distance_rate(speed, gamma)

    out = {
        "alt_rate": alt_rate,
        "vel_rate": vel_rate,
        "gamma_rate": gamma_rate,
        "down_range_rate": down_range_rate
    }

    return out

# Helpers for sim -----------------------------------------------
'''
index 0 = altitude 
index 1 = speed
index 2 = flight-path angle
index 3 = downrange
'''
STATE_VECTOR = []
DERIVATIVE_STATE_VECTOR = [] # list(output of eom_2d.values())

# Vector <-> dict helpers
"""
State vector layout:
    index 0 = altitude h [m]
    index 1 = speed V [m/s]
    index 2 = flight-path angle gamma [rad]
    index 3 = downrange s [m]
"""

def build_state_dict(state_vec: list) -> dict:
    """
    Convert a state vector [h, V, gamma, s] into the dict eom_2d expects.

    Parameters:
        state_vec: list ordered as altitude, speed, flight path angle,
                   downrange.

    Returns:
        A dict with keys alt, speed, gamma, s.
    """
    return {
        "alt": state_vec[0],
        "speed": state_vec[1],
        "gamma": state_vec[2],
        "s": state_vec[3],
    }


def build_derivative_vector(deriv_dict: dict) -> list:
    """
    Convert the eom_2d output dict into the ordered derivative vector.

    Parameters:
        deriv_dict: dict with keys alt_rate, vel_rate, gamma_rate,
                    down_range_rate.

    Returns:
        A list ordered as [h_dot, V_dot, gamma_dot, s_dot].
    """
    return [
        deriv_dict["alt_rate"],
        deriv_dict["vel_rate"],
        deriv_dict["gamma_rate"],
        deriv_dict["down_range_rate"],
    ]


def rhs_2d(t: float,
           y: list,
           sigma: float,
           vehicle_params: dict,
           env_params: dict) -> list:
    """
    Right hand side wrapper for numerical integrators.

    Given a time and a state vector y = [h, V, gamma, s], this returns the
    derivative dy/dt in the same order, by converting to the dict form, calling
    eom_2d, and converting the result back to a vector.

    Parameters:
        t:              time in seconds.
        y:              state vector [h, V, gamma, s].
        sigma:          bank angle in radians.
        vehicle_params: vehicle parameter dict.
        env_params:     environment parameter dict.

    Returns:
        The derivative vector [h_dot, V_dot, gamma_dot, s_dot].
    """
    state_dict = build_state_dict(y)
    deriv_dict = eom_2d(t, state_dict, sigma, vehicle_params, env_params)
    return build_derivative_vector(deriv_dict)
