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
# It returns an integer layer index (0 to 6) that is used to get a set of base values and lapse rate to use.
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


# ----------------------------------------------------------- BELOW IS THE INTERVAL MATH VERSION OF THE ABOVE ---------------------------------------------------

# converts geometric altitude (z as an Interval) into geopotential altitude h (Interval)
def intv_geopotential_altitude(z: Interval) -> Interval:
    R = constants.RADIUS_EARTH  # float constant
    numerator = R * z          # Interval via __rmul__
    denominator = R + z        # scalar promoted automatically
    return numerator / denominator


def intv_find_layers(h: Interval):
    """
    Given geopotential altitude interval h, return a list of layer indices
    that h intersects.

    If h is fully inside one layer -> list of length 1
    If h crosses boundaries -> list of multiple layers
    """

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

    layers = []
    for i in range(len(bounds) - 1):
        lo = bounds[i]
        hi = bounds[i + 1]
        if h.lo < hi and h.hi > lo:
            layers.append(i)

    lapse_rates_return = {}
    for l in layers:
        lapse_rates_return[l] = constants.intv_boundaries[l]

    return lapse_rates_return


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
            temp_interval = Interval(base_temp, base_temp)

        temps[k] = temp_interval

    return temps


def intv_pressure_in_layer(layers: dict, h: Interval, temperatures: dict) -> dict:
    """
    Interval pressure computation for each atmospheric layer.

    IMPORTANT MODELING NOTE:
    Interval arithmetic can cause temperature or pressure intervals
    to include non-physical values (such as <= 0) due to overestimation.
    Since pressure and temperature are strictly positive in reality,
    Enforce this known physical invariant via conservative clamping.

    This prevents undefined operations (such as log of non-positive values)
    while preserving a valid outer enclosure.
    """

    pressure = {}

    # Small positive lower bound to enforce physical domain
    # This does NOT assert pressure is this small, only that it is positive
    EPS = 1e-9

    for k in layers.keys():
        layer_params = get_layer_params(k)

        lapse_rate = layer_params["L_K_per_m"]
        base_temp = layer_params["T_b_K"]
        base_height = layer_params["h_b_m"]
        base_pressure = layer_params["p_b_Pa"]

        T = temperatures[k]

        # -------------------------------
        # Enforce temperature positivity
        # -------------------------------
        # Temperature must be strictly positive.
        # If interval overestimation causes T.lo <= 0,
        # we clamp it to a small positive value.
        T_safe = Interval(
            max(T.lo, EPS),
            max(T.hi, EPS)
        )

        if lapse_rate != 0.0:
            # Pressure law with lapse rate
            exponent = (-constants.g_0) / (lapse_rate * constants.R_AIR_DRY)

            # temp_ratio must be positive for log
            temp_ratio = T_safe / base_temp

            # Enforce positivity before log
            temp_ratio_safe = Interval(
                max(temp_ratio.lo, EPS),
                max(temp_ratio.hi, EPS)
            )

            factor = (temp_ratio_safe.log() * exponent).exp()
            pressure[k] = base_pressure * factor

        else:
            # Isothermal layer case
            exponent = (
                -constants.g_0 * (h - base_height)
            ) / (constants.R_AIR_DRY * base_temp)

            pressure[k] = base_pressure * exponent.exp()

    return pressure


def intv_density_from_pressure(pressures: dict, temperatures: dict) -> dict:
    """
    Interval density computation using ideal gas law.

    Density must be strictly positive.
    Enforce positivity of temperature to avoid division by zero
    or negative values due to interval overestimation.
    """

    densities = {}
    EPS = 1e-9

    for k in pressures.keys():
        T = temperatures[k]

        # Enforce temperature positivity
        T_safe = Interval(
            max(T.lo, EPS),
            max(T.hi, EPS)
        )

        densities[k] = pressures[k] / (constants.R_AIR_DRY * T_safe)

    return densities


def intv_US_Standard_ATM(z: Interval) -> float:
    h_geopotential_altitude = intv_geopotential_altitude(z)
    layer_info = intv_find_layers(h_geopotential_altitude)

    atmosphere_info = {}
    temperatures = intv_temp_in_layer(layer_info, h_geopotential_altitude)
    pressures = intv_pressure_in_layer(layer_info, h_geopotential_altitude, temperatures)
    densities = intv_density_from_pressure(pressures, temperatures)

    for k in layer_info.keys():
        atmosphere_info[k] = {
            "T_K": temperatures[k],
            "p_Pa": pressures[k],
            "rho_kgm3": densities[k],
        }

    return atmosphere_info


# ----------------------------------------------------------- TESTS (ONLY RUN WHEN EXECUTED DIRECTLY) ---------------------------------------------------
if __name__ == "__main__":

    print("\n--- Testing intv_geopotential_altitude ---")
    z_intv = Interval(10_000.0, 13_000.0)
    h_intv = intv_geopotential_altitude(z_intv)
    print("z interval:", z_intv)
    print("h (geopotential) interval:", h_intv)
    print("width (m):", h_intv.width())
    print("layers: ", intv_find_layers(h_intv))

    print("\n---------------------- Testing intv_temp_in_layer ----------------------")
    print("interval temperature(K) results: ", intv_temp_in_layer(intv_find_layers(h_intv), h_intv))

    print("\n---------------------- Testing intv_pressure_in_layer ----------------------")
    temps = intv_temp_in_layer(intv_find_layers(h_intv), h_intv)
    pressures = intv_pressure_in_layer(intv_find_layers(h_intv), h_intv, temps)
    print("interval pressure results(Pascals): ", pressures)

    print("\n---------------------- Testing intv_density_from_pressure ----------------------")
    densities = intv_density_from_pressure(pressures, temps)
    print("density computations", densities)

    print("\n-------------------- intv_US_Standard_ATM testing --------------------")
    atmosphere_info = intv_US_Standard_ATM(z_intv)
    for k in atmosphere_info.keys():
        print(f"{k}: {atmosphere_info[k]} \n")
"""
--- Testing intv_geopotential_altitude ---
z interval: Interval(10000.0, 13000.0)
h (geopotential) interval: Interval(9979.659213593904, 12979.649661088586)
width (m): 2999.9904474946816
layers:  {0: -0.0065, 1: 0.0}

---------------------- Testing intv_temp_in_layer ----------------------
interval temperature(K) results:  {0: Interval(203.78227720292418, 223.2822151116396), 1: Interval(216.65, 216.65)}

---------------------- Testing intv_pressure_in_layer ----------------------
interval pressure results(Pascals):  {0: Interval(16404.31960584069, 26518.687101234515), 1: Interval(16563.498028596357, 26582.821861831824)}

---------------------- Testing intv_density_from_pressure ----------------------
density computations {0: Interval(0.25594225971588885, 0.453339454740784), 1: Interval(0.2663368610060334, 0.42744505533429455)}

-------------------- intv_US_Standard_ATM testing --------------------
0: {'T_K': Interval(203.78227720292418, 223.2822151116396), 'p_Pa': Interval(16404.31960584069, 26518.687101234515), 'rho_kgm3': Interval(0.25594225971588885, 0.453339454740784)}

1: {'T_K': Interval(216.65, 216.65), 'p_Pa': Interval(16563.498028596357, 26582.821861831824), 'rho_kgm3': Interval(0.2663368610060334, 0.42744505533429455)}
"""