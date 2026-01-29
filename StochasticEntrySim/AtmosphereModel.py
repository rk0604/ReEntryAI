import constants
import math
from interval_math import Interval, promote, box_add, box_scalar_mul, scalar_times_interval

"""
Goal is to implement the 1976 US standard atmosphere model since it works great for altitudes less than 86 kilometers
Input: altitude (geometric, meters)
Output: 
 - Temperature T(h)
 - Pressure p(h)
 - Density rho(h)
"""

# Convert geometric altitude z (spacecraft altitude in meters) into geopotential altitude h
def geopotential_altitude(z) -> float:
    numerator = (constants.RADIUS_EARTH * z)
    denominator = (constants.RADIUS_EARTH + z)
    geo_alt = numerator / denominator
    return geo_alt

# determine which atmospheric layer (0-86km) contains the current geopotential_alt
# It returns an integer layer index (0 to 6) that tells you which set of base values and lapse rate to use.
# The rate at which temperature changes with altitude inside an atmospheric layer.
def find_layer(h: float):
    key = ""
    if h < 0.0:
        return -1
    elif h < 11_000.0:
        key = "troposphere"
    elif h < 20_000.0:
        key = "tropopause"
    elif h < 32_000.0:
        key = "stratosphere_1"
    elif h < 47_000.0:
        key = "stratosphere_2"
    elif h < 51_000.0:
        key = "stratopause"
    elif h < 71_000.0:
        key = "mesosphere_1"
    elif h < 86_000.0:
        key = "mesosphere_2"
    else:
        return -1  # above USSA76 range
    
    layer_index = constants.boundaries[key][0]
    lapse_rate = constants.boundaries[key][1]
    return layer_index, lapse_rate

# helper to retrieve base values for that layer
def get_layer_params(layer_index: int) -> dict:
    return constants.base_layers[layer_index]

# compute temperature at (geo-potential) altitude h inside that layer
def temperature_in_layer(h: float, layer_params: dict) -> float:
    lapse_rate = layer_params['L_K_per_m']
    base_temp = layer_params['T_b_K']   # NOTE: key is T_b_K (capital K)
    base_height = layer_params['h_b_m']

    if lapse_rate != 0.0:
        # T = T_b + L(h - h_b)
        lp_calc = h - base_height
        temp = base_temp + (lapse_rate * lp_calc)
        return temp
    else:
        # isothermal layer: T(h) = T_b
        return base_temp

# compute pressure in layer 
def pressure_in_layer(h: float, layer_params: dict, temp: float) -> float:
    lapse_rate = layer_params['L_K_per_m']
    base_temp = layer_params['T_b_K']   # NOTE: key is T_b_K
    base_height = layer_params['h_b_m']
    base_pressure = layer_params['p_b_Pa']

    if lapse_rate != 0.0:
        # p(h) = p_b * (T/T_b)^(-g0 / (L*R))
        exp_calc_1 = (-1.0 * constants.g_0) / (lapse_rate * constants.R_AIR_DRY)
        temp_ratio = temp / base_temp
        exp_calc_f = math.pow(temp_ratio, exp_calc_1)
        layer_pressure = base_pressure * exp_calc_f
        return layer_pressure
    else:
        # p(h) = p_b * exp( -g0 (h - h_b) / (R * T_b) )
        exponent = (-constants.g_0 * (h - base_height)) / (constants.R_AIR_DRY * base_temp)
        layer_pressure = base_pressure * math.exp(exponent)
        return layer_pressure

# compute air density from pressure and temperature (ideal gas law)
def density_from_pressure(pa: float, temp: float) -> float:
    # density = pressure / (R * T)
    density = pa / (constants.R_AIR_DRY * temp)
    return density

# use for altitude <= 86km
def US_Standard_ATM(z: float) -> dict:
    # 1) geometric -> geopotential altitude
    h = geopotential_altitude(z)

    # 2) layer
    layer_info = find_layer(h)
    if layer_info == -1:
        # outside 0–86 km USSA76 range
        return {
            "T_K": None,
            "p_Pa": None,
            "rho_kgm3": None,
            "geopotential_alt_m": h,
            "layer_index": -1
        }
    
    layer_index, _ = layer_info

    # 3) get base params for that layer
    layer_params = get_layer_params(layer_index)

    # 4) compute temperature at h
    T = temperature_in_layer(h, layer_params)

    # 5) compute pressure at h
    p = pressure_in_layer(h, layer_params, T)

    # 6) compute density at h
    rho = density_from_pressure(p, T)

    # 7) return all three plus some useful metadata
    return {
        "T_K": T,
        "p_Pa": p,
        "rho_kgm3": rho,
        "geopotential_alt_m": h,
        "layer_index": layer_index
    }


def atmosphere_from_radius(r: float):
    """Get standard atmosphere values based on altitude."""
    h = constants.altitude(r)
    if h < 0:
        h = 0.0
    return US_Standard_ATM(h)

# ----------------------------------------------------------- function testing
# example:
result = US_Standard_ATM(15000.0)
# print(result)
''' 
{
 'T_K': 216.65, 
  'p_Pa': 12111.618022126993, 
  'rho_kgm3': 0.1947517559484244, 
  'geopotential_alt_m': 14964.805806259874, 
  'layer_index': 1
} 
'''


# ----------------------------------------------------------- BELOW IS THE INTERVAL MATH VERSION OF THE ABOVE ---------------------------------------------------

# converts geometric altitude (z as an Interval) into geopotential altitude h (Interval)
def intv_geopotential_altitude(z: Interval) -> Interval:
    R = constants.RADIUS_EARTH  # float constant

    numerator = R * z          # Interval via __rmul__
    denominator = R + z        # scalar promoted automatically

    return numerator / denominator

'''
 Interval-safe layer handling (USSA76):
 Given h = Interval(h_lo, h_hi):
 1) Determine which layer boundaries are crossed by [h_lo, h_hi].
 2) If no boundary is crossed, evaluate the single layer normally.
 3) If boundaries are crossed:
    - Split h at each crossed boundary into sub-intervals h1, h2, ...
      where each sub-interval lies entirely within one layer.
    - For each sub-interval:
         compute T_i, p_i, rho_i using that layer’s formulas
    - Combine results with hull (component-wise):
         T = hull(T_1, T_2, ...)
         p = hull(p_1, p_2, ...)
         rho = hull(rho_1, rho_2, ...)
 This guarantees the final outputs enclose all possible true values.
'''
def intv_find_layers(h: Interval):
    """
    Given geopotential altitude interval h, return a list of layer indices
    that h intersects.

    If h is fully inside one layer -> list of length 1
    If h crosses boundaries -> list of multiple layers
    """

    # USSA76 geopotential layer boundaries (meters)
    # layer i is valid on [bounds[i], bounds[i+1])
    bounds = [
        0.0,
        11_000.0,
        20_000.0,
        32_000.0,
        47_000.0,
        51_000.0,
        71_000.0,
        86_000.0
    ]

    key = ""
    layers = []
    for i in range(len(bounds) - 1):
        lo = bounds[i]
        hi = bounds[i + 1]

        # interval intersection test
        if h.lo < hi and h.hi > lo:
            layers.append(i)

    lapse_rates_return = {} # (layer_index, lapse_rate)
    # find the layers' respective lapse rates and layer index
    for l in layers:
        l_r = constants.intv_boundaries[l]
        lapse_rates_return[l] = l_r

    return lapse_rates_return

print("\n--- Testing intv_geopotential_altitude ---")

z_intv = Interval(10_000.0, 13_000.0)   # geometric altitude interval (m)
h_intv = intv_geopotential_altitude(z_intv)

print("z interval:", z_intv)
print("h (geopotential) interval:", h_intv)
print("width (m):", h_intv.width())
print("layers: ", intv_find_layers(h_intv))
'''
--- Testing intv_geopotential_altitude ---
z interval: Interval(10000.0, 13000.0)
h (geopotential) interval: Interval(9979.659213593904, 12979.649661088586)
width (m): 2999.9904474946816
layers:  {0: -0.0065, 1: 0.0}
'''

# Compute interval temperature for each atmospheric layer intersected by h
"""
    Parameters:
    layers: dict
        Keys are layer indices (from intv_find_layers)
    h: Interval
        Geopotential altitude interval

    Returns:
    dict
        layer_index -> Interval temperature
"""
def intv_temp_in_layer(layers: dict, h: Interval) -> dict:

    temps = {}

    for k in layers.keys():
        layer_params = get_layer_params(k)

        lapse_rate = layer_params['L_K_per_m']
        base_temp = layer_params['T_b_K']
        base_height = layer_params['h_b_m']

        if lapse_rate != 0.0:
            temp_interval = base_temp + lapse_rate * (h - base_height)
        else:
            # Isothermal layer: temperature is constant
            temp_interval = Interval(base_temp, base_temp)

        temps[k] = temp_interval

    return temps
        
print("\n---------------------- Testing intv_temp_in_layer ----------------------")
print("interval temperature(K) results: ", intv_temp_in_layer(intv_find_layers(h_intv), h_intv))
# interval temperature(K) results:  {0: Interval(203.78227720292418, 223.2822151116396), 1: Interval(216.65, 216.65)}

# interval math variant of the pressure computation helper
"""
    Parameters:
    layers: dict
        Keys are layer indices (from intv_find_layers)
    h: Interval
        Geopotential altitude interval

    Returns:
    dict
        layer_index -> Interval pressure
"""
def intv_pressure_in_layer(layers: dict, h: Interval, temperatures: dict) -> dict:
    pressure = {}

    for k in layers.keys():
        layer_params = get_layer_params(k)

        lapse_rate    = layer_params['L_K_per_m']
        base_temp     = layer_params['T_b_K']
        base_height   = layer_params['h_b_m']
        base_pressure = layer_params['p_b_Pa']

        T = temperatures[k]  # Interval

        if lapse_rate != 0.0:
            # p = p_b * (T / T_b)^(-g0 / (L R))
            exponent = (-constants.g_0) / (lapse_rate * constants.R_AIR_DRY)
            temp_ratio = T / base_temp            # Interval (positive)
            factor = (temp_ratio.log() * exponent).exp()
            pressure[k] = base_pressure * factor

        else:
            # p = p_b * exp( -g0 (h - h_b) / (R T_b) )
            exponent = (
                -constants.g_0 * (h - base_height)
            ) / (constants.R_AIR_DRY * base_temp)
            pressure[k] = base_pressure * exponent.exp()

    return pressure

print("\n---------------------- Testing intv_pressure_in_layer ----------------------")
temperature_computations = intv_temp_in_layer(intv_find_layers(h_intv), h_intv)
pressure_computations = intv_pressure_in_layer(intv_find_layers(h_intv), h_intv, temperature_computations)
print("interval pressure results(Pascals): ", pressure_computations)
# interval pressure results(Pascals):  {0: Interval(16404.31960584069, 26518.687101234515), 1: Interval(16563.498028596357, 26582.821861831824)}

"""
    Parameters:
    pressures: dict
    temperatures: dict
        Keys are layer indices (from intv_find_layers)

    Returns:
    dict
        layer_index -> layer density
"""
# interval version compute air density from pressure and temperature (ideal gas law)
def intv_density_from_pressure(pressures: dict, temperatures:dict) -> dict:
    # density = pressure / (R * T)
    densities = {}
    for k in pressures.keys():
        density_intv = pressures[k] / (constants.R_AIR_DRY * temperatures[k])
        densities[k] = density_intv
    return densities
print("\n---------------------- Testing intv_density_from_pressure ----------------------")
density_computations = intv_density_from_pressure(pressure_computations, temperature_computations)
print("density copmutations", density_computations)
# density copmutations {0: Interval(0.25594225971588885, 0.453339454740784), 1: Interval(0.2663368610060334, 0.42744505533429455)}

"""
    Parameters:
    h: Interval geometric altitude

    Returns:
    Dict{
        "T_K": T,
        "p_Pa": p,
        "rho_kgm3": rho,
        "geopotential_alt_m": h,
        "layer_index": layer_index
        }
"""
# use for altitude <= 86km
def intv_US_Standard_ATM(z: Interval) -> float:
    # 1) Interval(geometric -> geopotential altitude)
    h = intv_geopotential_altitude(z)

    # 2) layer
    layer_info = intv_find_layers(h) # {0: -0.0065, 1: 0.0}
    if -1 in layer_info.values():
        # outside 0–86 km USSA76 range
        return {
            "T_K": None,
            "p_Pa": None,
            "rho_kgm3": None,
            "geopotential_alt_m": h,
            "layer_index": -1
        }
    
    atmosphere_info = {}

    # 3) Compute temperature at various h
    temperatures = intv_temp_in_layer(layer_info, h)

    # 4) compute pressure at various h
    pressures = intv_pressure_in_layer(layer_info, h, temperatures)

    # 5) compute density at various h
    densities = intv_density_from_pressure(pressures, temperatures)

    for k in layer_info.keys():
        atmosphere_info[k] = {
            "T_K": temperatures[k],
            "p_Pa": pressures[k],
            "rho_kgm3": densities[k],
        }

    #6) return all
    return atmosphere_info
    
print("-------------------- intv_US_Standard_ATM testing --------------------")
atmosphere_info = intv_US_Standard_ATM(z_intv)
for k in atmosphere_info.keys():
    print(f"{k}: {atmosphere_info[k]} \n")

# 0: {'T_K': Interval(203.78227720292418, 223.2822151116396), 'p_Pa': Interval(16404.31960584069, 26518.687101234515), 'rho_kgm3': Interval(0.25594225971588885, 0.453339454740784)} 
# 1: {'T_K': Interval(216.65, 216.65), 'p_Pa': Interval(16563.498028596357, 26582.821861831824), 'rho_kgm3': Interval(0.2663368610060334, 0.42744505533429455)} 
