import math
import AtmosphereModel
import constants

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

    h = constants.altitude(r)
    atm = AtmosphereModel.US_Standard_ATM(h)
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
