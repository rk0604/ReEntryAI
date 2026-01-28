
RADIUS_EARTH = 6378100  # meters
MU = 3.986e14           # 3.986*10^14 m^3/s^2
SEA_LEVEL_DENSITY = 1.225  # (kg/m^3)
H_M = 7200.0            # describes how quickly the atmosphere thins with the altitude
R_AIR_DRY = 287.053 # specific gas constanto of air J/(kg*K) 
R_U = 8.314             # J/(mol*k)
g_0 = 9.80665           # m/s^2

# 1976 atmosphere model layer boundaries
boundaries = {
    "troposphere":   (0, -0.0065),   # 0-11 km, L = -6.5 K/km
    "tropopause":    (1,  0.0),      # 11-20 km, L = 0
    "stratosphere_1":(2,  0.0010),   # 20-32 km, L = +1.0 K/km
    "stratosphere_2":(3,  0.0028),   # 32-47 km, L = +2.8 K/km
    "stratopause":   (4,  0.0),      # 47-51 km, L = 0
    "mesosphere_1":  (5, -0.0028),   # 51-71 km, L = -2.8 K/km
    "mesosphere_2":  (6, -0.0020),   # 71-86 km, L = -2.0 K/km
}

intv_boundaries = {
    0: -0.0065,
    1: 0.0,
    2: 0.0010,
    3: 0.0028,
    4: 0.0,
    5: -0.0028,
    6: -0.0020
}

# Basic environment functions
def altitude(r: float) -> float:
    """Convert radial distance from Earth center to altitude above surface."""
    return r - RADIUS_EARTH

def gravity(r: float) -> float:
    """Compute gravity magnitude at radius r."""
    return MU / (r * r)

# table for base values for the 1976 us std. atmosphere model
base_layers = {
    0: {
        "h_b_m":      0.0,          # 0 km
        "L_K_per_m": -0.0065,       # -6.5 K/km
        "T_b_K":     288.15,        # base temperature at 0 km
        "p_b_Pa":    101325.0       # base pressure at 0 km
    },
    1: {
        "h_b_m":      11000.0,      # 11 km
        "L_K_per_m":  0.0,          # isothermal
        "T_b_K":      216.65,       # base temperature at 11 km
        "p_b_Pa":     22632.1       # base pressure at 11 km
    },
    2: {
        "h_b_m":      20000.0,      # 20 km
        "L_K_per_m":  0.0010,       # +1.0 K/km
        "T_b_K":      216.65,       # base temperature at 20 km
        "p_b_Pa":     5474.89       # base pressure at 20 km
    },
    3: {
        "h_b_m":      32000.0,      # 32 km
        "L_K_per_m":  0.0028,       # +2.8 K/km
        "T_b_K":      228.65,       # base temperature at 32 km
        "p_b_Pa":     868.019       # base pressure at 32 km
    },
    4: {
        "h_b_m":      47000.0,      # 47 km
        "L_K_per_m":  0.0,          # isothermal
        "T_b_K":      270.65,       # base temperature at 47 km
        "p_b_Pa":     110.906       # base pressure at 47 km
    },
    5: {
        "h_b_m":      51000.0,      # 51 km
        "L_K_per_m": -0.0028,       # -2.8 K/km
        "T_b_K":      270.65,       # base temperature at 51 km
        "p_b_Pa":     66.9389       # base pressure at 51 km
    },
    6: {
        "h_b_m":      71000.0,      # 71 km
        "L_K_per_m": -0.0020,       # -2.0 K/km
        "T_b_K":      214.65,       # base temperature at 71 km
        "p_b_Pa":     3.95639       # base pressure at 71 km
        # Top of USSA76 definition is 86 km, where:
        # T = 186.95 K, p = 0.3734 Pa (not a new layer, just the top of this one)
    }
}
