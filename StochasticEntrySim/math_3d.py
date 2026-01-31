import math
import AtmosphereModel
import constants
from interval_math import Interval, promote, box_add, box_scalar_mul, scalar_times_interval

"""
State vector x = [r, phi, lam, V, gamma, chi]

Definitions (ASCII only):

r     = radial distance from Earth center (m)
phi   = geocentric latitude (rad)
lam   = longitude (rad)
V     = inertial speed (m/s)
gamma = flight path angle (rad)
chi   = heading angle (rad)

Control input:
sigma = bank angle (rad)
sigma = 0 means lift is up (normal to velocity plane)
sigma = pi means lift is down
sigma = +/- pi/2 gives pure lateral lift for turning

Vehicle and physical parameters:
m     = mass (kg)
S     = reference area (m^2)
CD    = drag coefficient
CL    = lift coefficient
mu    = Earth gravitational parameter (m^3/s^2)
RE    = mean Earth radius (m)
rho   = atmospheric density (kg/m^3)
q     = dynamic pressure (Pa)
"""

# Aerodynamic forces -------------------------------------------------------------
def aero_forces(r, V, gamma, chi, params, sigma):
    """
    Compute drag and lift components in the spherical 3D frame.

    Inputs:
    r      : radial distance (m)
    V      : speed (m/s)
    gamma  : flight path angle (rad)
    chi    : heading angle (rad)
    params : dict containing m, S, CD, CL
    sigma  : bank angle (rad)

    Output:
    (Dr, Dtheta, Dphi, Lr, Ltheta, Lphi, rho, q)
    """

    z = constants.altitude(r) # radius to geometric alt
    atm = AtmosphereModel.US_Standard_ATM(z)
    rho = atm["rho_kgm3"] if atm["rho_kgm3"] is not None else 0.0

    q = 0.5 * rho * V * V

    CD = float(params["CD"])
    CL = float(params["CL"])
    S  = float(params["S"])

    D = q * S * CD
    L = q * S * CL

    # velocity unit vector
    vr = math.sin(gamma)
    vtheta = math.cos(gamma) * math.cos(chi)
    vphi   = math.cos(gamma) * math.sin(chi)

    # drag components
    Dr     = -D * vr
    Dtheta = -D * vtheta
    Dphi   = -D * vphi

    # lift direction basis
    e1_r = math.cos(gamma)
    e1_theta = -math.sin(gamma) * math.cos(chi)
    e1_phi   = -math.sin(gamma) * math.sin(chi)

    e2_r = 0.0
    e2_theta = math.sin(chi)
    e2_phi   = -math.cos(chi)

    # banked lift direction
    cos_sig = math.cos(sigma)
    sin_sig = math.sin(sigma)

    Lr     = L * (e1_r * cos_sig + e2_r * sin_sig)
    Ltheta = L * (e1_theta * cos_sig + e2_theta * sin_sig)
    Lphi   = L * (e1_phi * cos_sig + e2_phi * sin_sig)

    return Dr, Dtheta, Dphi, Lr, Ltheta, Lphi, rho, q

# -------------------------------------------------------------
# 3D Equations of Motion
# -------------------------------------------------------------

def eom_3d(t, x, params, sigma_fn):
    """
    Compute time derivatives of the 6-state reentry model.

    States x = [r, phi, lam, V, gamma, chi]

    Returns dx/dt as a list of floats.
    """

    r, phi, lam, V, gamma, chi = x
    m = float(params["m"])

    # control input
    sigma = float(sigma_fn(t, x))

    # gravity
    g = constants.gravity(r)

    # aerodynamic forces
    Dr, Dtheta, Dphi, Lr, Ltheta, Lphi, rho, q = aero_forces(r, V, gamma, chi, params, sigma)

    # radial position rate
    r_dot = V * math.sin(gamma)

    # latitude rate
    phi_dot = (V * math.cos(gamma) * math.sin(chi)) / r

    # longitude rate
    lam_dot = (V * math.cos(gamma) * math.cos(chi)) / (r * math.cos(phi))

    # speed rate
    V_dot = (Dr + Lr) / m - g * math.sin(gamma)

    # flight path angle rate
    gamma_dot = (Ltheta / (m * V)) + (V / r - g / V) * math.cos(gamma)

    # heading angle rate
    chi_dot = (Lphi / (m * V * math.cos(gamma))) + (V / r) * math.sin(chi) * math.tan(phi)

    # make absolutely sure we return only floats
    return [
        float(r_dot),
        float(phi_dot),
        float(lam_dot),
        float(V_dot),
        float(gamma_dot),
        float(chi_dot),
    ]

# ------------------------------------------------------------- INTERVAL VERSION OF MATH_3D BELOW -------------------------------------------------------
"""
Compute drag and lift components in the spherical 3D frame.

Inputs: all are Interval()

r      : radial distance (m)
V      : speed (m/s)
gamma  : flight path angle (rad)
chi    : heading angle (rad)
params : dict containing m, S, CD, CL
sigma  : bank angle (rad)

Output:
(Dr, Dtheta, Dphi, Lr, Ltheta, Lphi, rho, q)
"""
def intv_aero_forces(r: Interval, V: Interval, gamma: Interval, chi: Interval, sigma: Interval, params: dict) -> dict:

    z_geometric_altitude = constants.intv_geometric_altitude(r)
    atm_intv = AtmosphereModel.intv_US_Standard_ATM(z_geometric_altitude)
    
    aero_forces = {} # key is the layer index; value is the heat value associated

    for k in atm_intv.keys():
        aero_forces[k] = {} # declare empty dict to hold 
        rho = atm_intv[k]['rho_kgm3'] if atm_intv[k]['rho_kgm3'] is not None else 0.0

        # print(atm_intv)
        # print(rho)
        q = rho * 0.5 * V * V
        aero_forces[k]["q"] = q

        CD = float(params["CD_const"])
        CL = float(params["CL_const"])
        surface_area  = float(params["ref_area_m2"])

        Drag_mag = q * surface_area * CD
        Lift_mag = q * surface_area * CL
        
        # velocity components
        v_r = Interval.sin(gamma)                        # velocity component radially (up or down)
        vTheta = Interval.cos(gamma) * Interval.cos(chi) # velocity component (n/s)
        vPhi = Interval.cos(gamma) * Interval.sin(chi)   # velocity component (e/w)

        # drag components
        D_r = -Drag_mag * v_r                           # drag magnitude in the radial direction
        Dtheta = Drag_mag * vTheta                      # drag magnitude (n/s)
        Dphi = -Drag_mag * vPhi                         # drag magnitude (e/w)

        # lift direction basis
        """
        - Lift must always act perpendicular to the velocity direction.
        - e1 is the “upward” lift direction (changes altitude and flight-path angle).
        - e2 is the “sideways” lift direction (causes turning / crossrange)
        - e1 and e2 are unit vectors and both are perpendicular to velocity
        - The bank angle sigma rotates lift between e1 (up) and e2 (sideways)
        - Final lift = L * (cos(sigma) * e1 + sin(sigma) * e2).
        """

        e1_r = Interval.cos(gamma)
        e1_theta = -1 * Interval.sin(gamma) * Interval.cos(chi) # north/south part of upward lift
        e1_phi = -1 * Interval.sin(gamma) * Interval.sin(chi)   # east/west part of upward lift

        e2_r = promote(0.0)
        e2_theta = Interval.sin(chi)                                # north/south part of sideways lift
        e2_phi   = -Interval.cos(chi)                               # east/west part of sideways lift

        # banked lift direction
        cos_sig = Interval.cos(sigma)
        sin_sig = Interval.sin(sigma)

        Lr     = Lift_mag * (e1_r * cos_sig + e2_r * sin_sig)
        Ltheta = Lift_mag * (e1_theta * cos_sig + e2_theta * sin_sig)
        Lphi   = Lift_mag * (e1_phi * cos_sig + e2_phi * sin_sig)

        aero_forces[k]["Dr"] = D_r
        aero_forces[k]["Dtheta"] = Dtheta
        aero_forces[k]["Dphi"] = Dphi
        aero_forces[k]["Lr"] = Lr
        aero_forces[k]["Ltheta"] = Ltheta
        aero_forces[k]["Lphi"] = Lphi
        aero_forces[k]["rho"] = rho

    return aero_forces


# 3d equations of motion interval variants
"""
Interval-valued 3D equations of motion for atmospheric entry.

Inputs:
    t      : time (seconds), scalar, not interval-propagated here
    X      : state vector [r, phi, lam, V, gamma, chi], each an Interval
                r     -> radial distance from Earth center (m)
                phi   -> geocentric latitude (rad)
                lam   -> longitude (rad)
                V     -> speed (m/s)
                gamma -> flight path angle (rad)
                chi   -> heading angle (rad)
    params : vehicle parameters (mass, aero coefficients, reference area)
    sigma  : bank angle (rad), Interval

Output:
    dict mapping state derivative names to Interval values
"""
def intv_eom_3d(t: float, X: list, params: dict, sigma: Interval) -> dict:
    # unpack state intervals
    r     = X[0]  # radial distance interval
    phi   = X[1]  # latitude interval
    lam   = X[2]  # longitude interval
    V     = X[3]  # speed interval
    gamma = X[4]  # flight path angle interval
    chi   = X[5]  # heading angle interval

    # vehicle mass (treated as exact scalar)
    m = float(params["mass_kg"])

    # gravitational acceleration as an interval function of radius
    g = constants.intv_gravity(r)

    # compute interval aerodynamic forces (possibly layer-wise)
    aero = intv_aero_forces(
        r=r,
        V=V,
        gamma=gamma,
        chi=chi,
        sigma=sigma,
        params=params,
    )

    # hull forces across all intersected atmosphere layers
    Dr     = None  # radial drag interval
    Dtheta = None  # north south drag interval
    Dphi   = None  # east west drag interval
    Lr     = None  # radial lift interval
    Ltheta = None  # north south lift interval
    Lphi   = None  # east west lift interval

    for layer in aero.keys():
        if Dr is None:
            Dr     = aero[layer]["Dr"]
            Dtheta = aero[layer]["Dtheta"]
            Dphi   = aero[layer]["Dphi"]
            Lr     = aero[layer]["Lr"]
            Ltheta = aero[layer]["Ltheta"]
            Lphi   = aero[layer]["Lphi"]
        else:
            Dr     = box_add(Dr, aero[layer]["Dr"])
            Dtheta = box_add(Dtheta, aero[layer]["Dtheta"])
            Dphi   = box_add(Dphi, aero[layer]["Dphi"])
            Lr     = box_add(Lr, aero[layer]["Lr"])
            Ltheta = box_add(Ltheta, aero[layer]["Ltheta"])
            Lphi   = box_add(Lphi, aero[layer]["Lphi"])

    # kinematic equations
    r_dot = V * gamma.sin()                         # radial rate
    phi_dot = (V * gamma.cos() * chi.sin()) / r    # latitude rate
    lam_dot = (V * gamma.cos() * chi.cos()) / (r * phi.cos())  # longitude rate

    # dynamic equations
    V_dot = (Dr + Lr) / m - g * gamma.sin()         # speed rate
    gamma_dot = (Ltheta / (m * V)) + (V / r - g / V) * gamma.cos()  # flight path angle rate
    phi_tan = (Interval.sin(phi)) / (Interval.cos(phi))
    chi_dot = (Lphi / (m * V * gamma.cos())) + (V / r) * chi.sin() * phi_tan  # heading rate

    return {
        "r_dot": r_dot,
        "phi_dot": phi_dot,
        "lam_dot": lam_dot,
        "V_dot": V_dot,
        "gamma_dot": gamma_dot,
        "chi_dot": chi_dot,
    }

# ------------------------------------------ testing ---------------------------------------------------------------
r_iv = Interval(
    constants.RADIUS_EARTH + 58_000.0,
    constants.RADIUS_EARTH + 62_000.0
)

V_iv = Interval(7600.0, 7800.0)                 # m/s
gamma_iv = Interval(
    math.radians(-5.2),
    math.radians(-4.8)
)                                               # rad
chi_iv = Interval(
    math.radians(85.0),
    math.radians(95.0)
)                                               # rad
sigma_iv = Interval(
    math.radians(-5.0),
    math.radians(5.0)
)       
print("\n===== Interval Aero Forces Test =====\n")

vehicle_params = {
    "mass_kg":     5000.0,
    "ref_area_m2": 10.0,
    "CL_const":    0.30,
    "CD_const":    1.00,
    "nose_radius_m": 1.0,
}

env_params = {
    "mu_m3s2":   3.986004418e14,
    "R_E_m":     6371000.0,
    "rho0_kgm3": 1.225,
    "H_m":       7200.0,
}

# test
aero = intv_aero_forces(
    r=r_iv,
    V=V_iv,
    gamma=gamma_iv,
    chi=chi_iv,
    sigma=sigma_iv,
    params=vehicle_params,
)

for layer, vals in aero.items():
    print(f"\nLayer {layer}")
    print("-" * 30)
    for k, v in vals.items():
        print(f"{k:8s}: {v}")
'''
Layer 5
------------------------------
q       : Interval(6608.968320232319, 12666.867962587668)
Dr      : Interval(5530.242156886362, 11480.30926473919)
Dtheta  : Interval(-11001.184362471855, 11001.184362471842)
Dphi    : Interval(-126224.4347377493, -65567.22825058083)
Lr      : Interval(19670.168475174247, 37867.33042132479)
Ltheta  : Interval(-3612.143320979881, 3612.143320979881)
Lphi    : Interval(1495.862507770258, 3732.7500593958284)
rho     : Interval(0.00022884239335984482, 0.00041639934130794444)
'''

# --------------------- Interval EOM 3D Test ---------------------

# interval state vector
X_iv = [
    r_iv,        # radial distance interval
    Interval(0.3, 0.31),   # latitude (rad)
    Interval(1.0, 1.01),   # longitude (rad)
    V_iv,        # speed interval
    gamma_iv,    # flight path angle interval
    chi_iv,      # heading angle interval
]

# call interval equations of motion
xdot = intv_eom_3d(
    t=0.0,
    X=X_iv,
    params=vehicle_params,
    sigma=sigma_iv,
)

# pretty print results
print("\n===== Interval EOM 3D Output =====")
for k, v in xdot.items():
    print(f"{k:10s}: {v}")
'''
===== Interval EOM 3D Output =====
r_dot     : Interval(-706.9341255426853, -635.9516093255977)
phi_dot   : Interval(0.0011707768760396147, 0.0012076636941853193)
lam_dot   : Interval(-0.00011052306627050185, 0.00011052306627050171)
V_dot     : Interval(5.84427931542488, 10.741647587210721)
gamma_dot : Interval(-0.00018077833016133234, 7.492100188617948e-05)
chi_dot   : Interval(0.0004021508495562784, 0.0004868456152008609)
'''