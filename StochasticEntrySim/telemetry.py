"""
Shared per-timestep telemetry: derived physics quantities and RL reward
ingredients computed identically by both drivers (run_sim.py and the
notebook) so the generated dataset is consistent regardless of driver.

Everything here is a pure function of already-logged state — no
simulation side effects — so values are fully reproducible and traceable.

Quantities (all SI unless noted):
    mach                       V / a, local speed of sound from atmosphere T
    fmv                        Bibb 2010 blended velocity parameter (passthrough)
    specific_kinetic_e_j_kg    0.5 * V^2
    specific_potential_e_j_kg  g0 * altitude  (energy-height convention)
    specific_energy_j_kg       kinetic + potential (entry-guidance energy E)
    energy_height_m            E / g0  (the "energy altitude")
    load_factor_g              sqrt(L^2 + D^2) / (m * g0)   sensed aero g's
    aero_decel_mps2            sqrt(L^2 + D^2) / m
    drag_decel_mps2            D / m
    dynamic_pressure_pa        q = 0.5 * rho * V^2  (echo, for completeness)

RL reward (documented default — override weights in the RL env layer):
    reward_default             a transparent weighted sum of the terms below
    rew_range_term             - w_range  * range_to_go_m
    rew_heat_term              - w_heat   * max(0, qdot - qdot_limit)
    rew_gload_term             - w_gload  * max(0, load_factor_g - g_limit)
    rew_fuel_term              - w_fuel   * rcs_fuel_rate_kg_s
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Optional

import constants

# Air properties for Mach (calorically perfect air below ~90 km; above that
# Mach loses its usual meaning but the ratio is still a useful logged scalar).
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
    """Return a dict of derived per-step physics scalars."""
    V = float(V_mps)
    alt = float(altitude_m)

    # Mach number (guard missing / non-physical temperature)
    if T_K is not None and T_K == T_K and T_K > 0.0:
        a_sound = math.sqrt(_GAMMA_AIR * _R_AIR * float(T_K))
        mach = V / a_sound if a_sound > 0.0 else float("nan")
    else:
        a_sound = float("nan")
        mach = float("nan")

    # Energy (energy-height convention used in entry guidance)
    ke = 0.5 * V * V
    pe = _G0 * alt
    E = ke + pe
    energy_height_m = E / _G0 if _G0 > 0.0 else float("nan")

    # Load factor (sensed aerodynamic acceleration in g's)
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


# =============================================================================
# RCS propellant model
# =============================================================================

def rcs_propellant_kg(
    fired_thruster_seconds: float,
    thrust_per_thruster_N: float = float(constants.ORION_CM_RCS_THRUST_N),
    isp_s: float = float(constants.ORION_CM_RCS_ISP_S),
) -> float:
    """Propellant mass burned for a given amount of thruster on-time.

        mdot_per_thruster = thrust / (Isp * g0)
        propellant        = mdot * fired_thruster_seconds

    `fired_thruster_seconds` is the SUM over all thrusters of their on-time
    this step (e.g. 3 thrusters firing for 0.1 s = 0.3 thruster-seconds).
    """
    if isp_s <= 0.0:
        return 0.0
    mdot = thrust_per_thruster_N / (isp_s * _G0)
    return float(mdot) * float(fired_thruster_seconds)


# =============================================================================
# RL reward (documented default)
# =============================================================================

@dataclass
class RewardWeights:
    """Default reward weights. Override in the RL environment layer; these
    exist so the generated dataset ships with a transparent, reproducible
    baseline reward column rather than baking the choice in irreversibly."""
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
    """Return the per-step reward ingredients and the composite default.

    The composite is a transparent negative-cost sum; the RL env is free to
    re-weight or replace it. Terminal bonuses (on-target, soft touchdown)
    are intentionally NOT included here — they belong to the episode layer.
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
