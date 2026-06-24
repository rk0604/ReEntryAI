"""
Shared per timestep telemetry: derived physics quantities and RL reward
ingredients computed the same way by both drivers (run_sim.py and the notebook),
so the generated dataset is consistent regardless of which driver produced it.

Everything here is a pure function of already logged state, with no simulation
side effects, so the values are fully reproducible and traceable.

Quantities (all SI unless noted):
    mach                       V divided by the local speed of sound
    fmv                        Bibb 2010 blended velocity parameter (passthrough)
    specific_kinetic_e_j_kg    0.5 * V^2
    specific_potential_e_j_kg  g0 * altitude, using the energy height convention
    specific_energy_j_kg       kinetic plus potential, the entry guidance energy
    energy_height_m            E / g0, the energy altitude
    load_factor_g              sqrt(L^2 + D^2) / (m * g0), the sensed aero g load
    aero_decel_mps2            sqrt(L^2 + D^2) / m
    drag_decel_mps2            D / m
    dynamic_pressure_pa        q = 0.5 * rho * V^2, echoed for completeness

RL reward (the documented default; weights can be overridden in the RL env):
    reward_default             a transparent weighted sum of the terms below
    rew_range_term             minus w_range  times range_to_go_m
    rew_heat_term              minus w_heat   times max(0, qdot minus qdot_limit)
    rew_gload_term             minus w_gload  times max(0, load_factor_g minus g_limit)
    rew_fuel_term              minus w_fuel   times rcs_fuel_rate_kg_s
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Optional

import constants

# Air properties used for Mach. Air is treated as calorically perfect below
# about 90 km. Above that the ratio is still a useful logged scalar even though
# Mach loses its usual meaning.
_GAMMA_AIR = 1.4
_R_AIR = 287.053
_G0 = float(constants.STANDARD_GRAVITY_MPS2)


def compute_derived_metrics(
    *,
    r_m: float,
    V_mps: float,
    altitude_m: float,
    T_K: Optional[float],
    rho_kgm3: float,
    drag_mag_N: float,
    lift_mag_N: float,
    mass_kg: float,
    fmv: float = float("nan"),
) -> Dict[str, float]:
    """
    Return a dictionary of derived per step physics scalars.

    Parameters:
        r_m:        distance from the center of the Earth in meters.
        V_mps:      speed in meters per second.
        altitude_m: geometric altitude in meters.
        T_K:        static temperature in kelvin, or None when not available.
        rho_kgm3:   atmospheric density in kilograms per cubic meter.
        drag_mag_N: drag force magnitude in newtons.
        lift_mag_N: lift force magnitude in newtons.
        mass_kg:    vehicle mass in kilograms.
        fmv:        the blended velocity parameter, passed through unchanged.

    Returns:
        A dictionary with the quantities listed in the module docstring. Mach
        and the speed of sound are NaN when the temperature is missing or not
        positive.
    """
    V = float(V_mps)
    alt = float(altitude_m)

    # Mach number. Guard a missing or non positive temperature.
    if T_K is not None and T_K == T_K and T_K > 0.0:
        a_sound = math.sqrt(_GAMMA_AIR * _R_AIR * float(T_K))
        mach = V / a_sound if a_sound > 0.0 else float("nan")
    else:
        a_sound = float("nan")
        mach = float("nan")

    # Energy, using the energy height convention common in entry guidance.
    ke = 0.5 * V * V
    pe = _G0 * alt
    E = ke + pe
    energy_height_m = E / _G0 if _G0 > 0.0 else float("nan")

    # Load factor, the sensed aerodynamic acceleration expressed in g.
    aero_force_N = math.sqrt(float(drag_mag_N) ** 2 + float(lift_mag_N) ** 2)
    m = float(mass_kg)
    aero_decel = aero_force_N / m if m > 0.0 else float("nan")
    drag_decel = float(drag_mag_N) / m if m > 0.0 else float("nan")
    load_factor_g = aero_decel / _G0 if _G0 > 0.0 else float("nan")

    q_pa = 0.5 * float(rho_kgm3) * V * V

    return {
        "mach": float(mach),
        "speed_of_sound_mps": float(a_sound),
        "fmv": float(fmv),
        "specific_kinetic_e_j_kg": float(ke),
        "specific_potential_e_j_kg": float(pe),
        "specific_energy_j_kg": float(E),
        "energy_height_m": float(energy_height_m),
        "aero_force_N": float(aero_force_N),
        "aero_decel_mps2": float(aero_decel),
        "drag_decel_mps2": float(drag_decel),
        "load_factor_g": float(load_factor_g),
        "dynamic_pressure_pa": float(q_pa),
    }


# RCS propellant model


def rcs_propellant_kg(
    fired_thruster_seconds: float,
    thrust_per_thruster_N: float = float(constants.ORION_CM_RCS_THRUST_N),
    isp_s: float = float(constants.ORION_CM_RCS_ISP_S),
) -> float:
    """
    Return the propellant mass burned for a given amount of thruster on time.

    The mass flow per thruster is thrust divided by (Isp times g0), and the
    propellant is that flow times the total fired thruster seconds.

    Parameters:
        fired_thruster_seconds: total on time summed across all thrusters this
                                step. For example three thrusters firing for
                                0.1 s each is 0.3 thruster seconds.
        thrust_per_thruster_N:  thrust of one jet in newtons.
        isp_s:                  specific impulse in seconds.

    Returns:
        The propellant mass burned in kilograms. Returns zero when Isp is not
        positive.
    """
    if isp_s <= 0.0:
        return 0.0
    mdot = thrust_per_thruster_N / (isp_s * _G0)
    return float(mdot) * float(fired_thruster_seconds)


# RL reward (the documented default)


@dataclass
class RewardWeights:
    """
    Default reward weights.

    These can be overridden in the RL environment layer. They exist so the
    generated dataset ships with a transparent, reproducible baseline reward
    column instead of baking the choice in permanently.

    Fields:
        w_range:         penalty per meter of range to go (the miss).
        w_heat:          penalty per W/m^2 over the heat rate limit.
        w_gload:         penalty per g over the load factor limit.
        w_fuel:          penalty per kg/s of RCS propellant flow.
        qdot_limit_w_m2: heat rate limit used in the heat penalty.
        gload_limit_g:   load factor limit used in the g penalty.
    """
    w_range: float = 1.0e-5     # per metre of range-to-go (miss)
    w_heat: float = 1.0e-6      # per W/m^2 over the heat-rate limit
    w_gload: float = 1.0e-1     # per g over the load-factor limit
    w_fuel: float = 1.0e1       # per kg/s of RCS propellant flow
    qdot_limit_w_m2: float = float(constants.HEAT_RATE_LIMIT_DEFAULT)
    gload_limit_g: float = 8.0  # crewed entry structural/physiological cap


def composite_reward(
    *,
    range_to_go_m: float,
    qdot_w_m2: float,
    load_factor_g: float,
    rcs_fuel_rate_kg_s: float,
    weights: Optional[RewardWeights] = None,
) -> Dict[str, float]:
    """
    Return the per step reward ingredients and the composite default.

    The composite is a transparent negative cost sum. The RL environment is
    free to re weight or replace it. Terminal bonuses, such as landing on target
    or a soft touchdown, are intentionally left out here, since they belong to
    the episode layer.

    Parameters:
        range_to_go_m:      remaining great circle distance to target in meters.
        qdot_w_m2:          current heat rate in W/m^2.
        load_factor_g:      current load factor in g.
        rcs_fuel_rate_kg_s: current RCS propellant flow in kg/s.
        weights:            a RewardWeights instance, or None to use defaults.

    Returns:
        A dictionary with the four named reward terms and the composite
        reward_default, which is their sum.
    """
    w = weights or RewardWeights()
    rew_range = -w.w_range * max(0.0, float(range_to_go_m))
    rew_heat = -w.w_heat * max(0.0, float(qdot_w_m2) - w.qdot_limit_w_m2)
    rew_gload = -w.w_gload * max(0.0, float(load_factor_g) - w.gload_limit_g)
    rew_fuel = -w.w_fuel * max(0.0, float(rcs_fuel_rate_kg_s))
    return {
        "rew_range_term": float(rew_range),
        "rew_heat_term": float(rew_heat),
        "rew_gload_term": float(rew_gload),
        "rew_fuel_term": float(rew_fuel),
        "reward_default": float(rew_range + rew_heat + rew_gload + rew_fuel),
    }
