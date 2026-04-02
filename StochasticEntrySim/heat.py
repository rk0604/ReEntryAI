'''
{
  "t": 183.2,
  "triangles": [
    {"id": 0, "qdot": 615000, "Q": 2.41e8},
    {"id": 1, "qdot": 602000, "Q": 2.35e8},
    {"id": 2, "qdot": 540000, "Q": 2.00e8}
    // ...
  ]
}  

Better because it's a fixed-length, control-relevant summary of the full triangle mesh (and stochastic atmosphere), 
so model can learn stable policies without drowning in hundreds of weakly-informative per-cell values
RL observation (recommended): [h, v, gamma, s, delta_rho,
                               qdot_max, Q_max, qdot_mean,
                               frac_over_qdot_limit, frac_over_Q_limit,
                               center_ring_qdot_mean, shoulder_ring_qdot_mean, nonuniformity]

'''
import math
import AtmosphereModel
import constants
from interval_math import Interval, promote, box_add, box_scalar_mul, scalar_times_interval
 
# Compute stagnation-point convective heating using Sutton-Graves
"""
Inputs:
    rho : Interval   (local atmospheric density)
    V   : Interval   (vehicle speed)
    R_n : float      (nose radius)

Output:
    qdot_stag : Interval (W/m^2)

This defines the peak heating on the heat shield and sets the scale
for all spatial heating distributions.
"""
def intv_convective_SG_heating(rho: Interval, V: Interval, R_n: Interval) -> Interval:
    pass

# Radial distribution of convective heating across a circular heat shield
"""
Inputs:
    qdot_stag : Interval (stagnation-point heat flux)
    r         : float    (radial position on shield)
    R_shield  : float    (shield radius)

Output:
    qdot_r : Interval (local heat flux)

Models non-uniform heating from center to shoulder.
"""
def intv_radial_q_dist(qdot_stag: Interval, r: float, R_shield: float) -> Interval:
    pass

# Compute mean heat flux over an annular ring of the heat shield
"""
Inputs:
    qdot(r) : callable or discretized profile
    r_inner, r_outer : float (ring bounds)

Output:
    qdot_mean_ring : Interval

Used to track heating on center and shoulder regions separately.
"""
def intv_ring_avg_heat_flux(qdot: Interval, r_inner: float, r_outter: float) -> Interval:
    pass

# Time-integrated heat load for a heat-shield ring
"""
Inputs:
    qdot_mean_ring : Interval (W/m^2)
    dt             : float

Output:
    Q_ring : Interval (J/m^2)

Tracks thermal stress accumulation over the trajectory.
"""
def intv_total_heat_load(qdot_mean_ring: Interval, dt: float) -> Interval:
    pass

"""
Measure spatial non-uniformity of heat flux across the shield.

Output:
    nonuniformity : Interval

High values indicate localized overheating and structural risk,
even if average heating is acceptable.
"""
def intv_non_uniformity_heat_flux(max_flux: Interval, min_flux: Interval) -> Interval:
    pass

# Aggregate heat-shield metrics used for constraints and learning signals.
"""
Outputs:
    qdot_max  : Interval (peak instantaneous heat flux)
    Q_max     : Interval (peak integrated heat load)
    qdot_mean : Interval (area-averaged heat flux)

These summarize overall thermal severity.
"""
def aggregate_heat_metrics(qdot_max: Interval, Q_max: Interval, qdot_mean: Interval) -> dict:
  pass
