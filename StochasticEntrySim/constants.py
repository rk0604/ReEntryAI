"""
constants.py

Central location for:
- Earth and atmosphere constants
- simulation cadence constants
- capsule geometry and RCS constants
- first-pass rigid-body attitude control tuning
- heat shield ring model and helper functions

Notes
-----
1. This file keeps your existing atmosphere and heat shield structure.
2. The new capsule and RCS values are first-pass engineering constants
   for milestone 1 of the rigid-body refactor.
3. These are not meant to claim exact Orion flight-certified values.
   They are practical simulation values chosen to fit your current codebase.
"""

import math
import numpy as np

from interval_math import Interval, promote


# ============================================================================
# Simulation cadence
# ============================================================================

# Fixed Euler integration step used by the main physics loop.
# Your current sim runs at 4 Hz, so dt = 0.25 s.
ENTRY_DT_S = 0.25

# Guidance updates more slowly than physics.
# Your current architecture uses 1 Hz guidance.
GUIDANCE_PERIOD_S = 1.0

# Convenience ratio: number of physics steps per guidance update.
GUIDANCE_STEPS = int(round(GUIDANCE_PERIOD_S / ENTRY_DT_S))


# ============================================================================
# Earth and environment constants
# ============================================================================

# Mean Earth radius in meters.
RADIUS_EARTH = 6_378_100.0

# Earth standard gravitational parameter, m^3 / s^2.
MU = 3.986e14

# Sea-level density, kg / m^3.
SEA_LEVEL_DENSITY = 1.225

# Simple atmosphere exponential scale height in meters.
# This is still useful for rough fallback models even though you also
# keep a layered US Standard Atmosphere table below.
H_M = 7200.0

# Gas constant for dry air, J / (kg K).
R_AIR_DRY = 287.053

# Universal gas constant, J / (mol K).
R_U = 8.314

# Standard gravity at Earth surface, m / s^2.
g_0 = 9.80665


# ============================================================================
# Atmosphere layer definitions
# ============================================================================
# These tables preserve your current US Standard Atmosphere 1976 structure.

# Layer names mapped to:
# (layer_index, lapse_rate_K_per_m)
boundaries = {
    "troposphere":    (0, -0.0065),  # 0 to 11 km
    "tropopause":     (1,  0.0),     # 11 to 20 km
    "stratosphere_1": (2,  0.0010),  # 20 to 32 km
    "stratosphere_2": (3,  0.0028),  # 32 to 47 km
    "stratopause":    (4,  0.0),     # 47 to 51 km
    "mesosphere_1":   (5, -0.0028),  # 51 to 71 km
    "mesosphere_2":   (6, -0.0020),  # 71 to 86 km
}

# Same lapse-rate information, keyed only by layer index.
intv_boundaries = {
    0: -0.0065,
    1:  0.0,
    2:  0.0010,
    3:  0.0028,
    4:  0.0,
    5: -0.0028,
    6: -0.0020,
}

# Base values for the US Standard Atmosphere 1976.
# Each layer stores its geometric base altitude, lapse rate,
# base temperature, and base pressure.
base_layers = {
    0: {
        "h_b_m": 0.0,
        "L_K_per_m": -0.0065,
        "T_b_K": 288.15,
        "p_b_Pa": 101325.0,
    },
    1: {
        "h_b_m": 11000.0,
        "L_K_per_m": 0.0,
        "T_b_K": 216.65,
        "p_b_Pa": 22632.1,
    },
    2: {
        "h_b_m": 20000.0,
        "L_K_per_m": 0.0010,
        "T_b_K": 216.65,
        "p_b_Pa": 5474.89,
    },
    3: {
        "h_b_m": 32000.0,
        "L_K_per_m": 0.0028,
        "T_b_K": 228.65,
        "p_b_Pa": 868.019,
    },
    4: {
        "h_b_m": 47000.0,
        "L_K_per_m": 0.0,
        "T_b_K": 270.65,
        "p_b_Pa": 110.906,
    },
    5: {
        "h_b_m": 51000.0,
        "L_K_per_m": -0.0028,
        "T_b_K": 270.65,
        "p_b_Pa": 66.9389,
    },
    6: {
        "h_b_m": 71000.0,
        "L_K_per_m": -0.0020,
        "T_b_K": 214.65,
        "p_b_Pa": 3.95639,
        # Top of USSA76 is 86 km.
        # Approx values there are T = 186.95 K and p = 0.3734 Pa.
    },
}


# ============================================================================
# Basic environment helper functions
# ============================================================================

def altitude(r: float) -> float:
    """
    Convert radial distance from Earth center to geometric altitude
    above the Earth surface.
    """
    return r - RADIUS_EARTH


def intv_geometric_altitude(r) -> Interval:
    """
    Interval version of radial distance to geometric altitude.
    """
    r_iv = promote(r)
    return r_iv - RADIUS_EARTH


def gravity(r: float) -> float:
    """
    Gravitational acceleration magnitude at radius r.
    """
    return MU / (r * r)


def intv_gravity(r) -> Interval:
    """
    Interval version of gravitational acceleration magnitude.
    """
    r_iv = promote(r)
    return MU / (r_iv * r_iv)


# ============================================================================
# Capsule geometry and aerodynamic reference values
# ============================================================================
# These values support the new rigid-body capsule model.
# They are practical milestone-1 approximations.

# Crew module outer diameter, meters.
# Orion is roughly 5 m class, so this is a good first-pass value.
ORION_CM_DIAMETER_M = 5.0

# Outer radius, meters.
ORION_CM_RADIUS_M = 0.5 * ORION_CM_DIAMETER_M

# Heat shield radius, meters.
# For milestone 1, treat this as the same characteristic radius.
HEAT_SHIELD_RADIUS_M = ORION_CM_RADIUS_M

# Effective nose radius for Sutton-Graves style stagnation heating.
# This is a tuning input to your current HeatShield model.
HEAT_SHIELD_NOSE_RADIUS_M = 1.0

# Default heat shield ring discretization.
HEAT_SHIELD_NUM_RINGS = 5

# Radial exponent for nonuniform heat-flux shaping over the shield.
# Larger values push more heating concentration toward the center.
HEAT_SHIELD_RADIAL_EXP = 1.0

# Simple projected reference area for drag/lift models.
# This is the circular projected area of the capsule heat shield.
CAPSULE_REFERENCE_AREA_M2 = math.pi * HEAT_SHIELD_RADIUS_M ** 2

# First-pass capsule mass for milestone 1.
# Keep using your existing mass elsewhere if it already exists.
# This value is here so the rigid-body constants live together.
CAPSULE_MASS_KG = 8500.0


# ============================================================================
# RCS layout constants
# ============================================================================
# These constants support the first-pass 12-thruster body-fixed layout.

# Total number of crew module RCS thrusters in the milestone-1 layout.
CAPSULE_RCS_NUM_THRUSTERS = 12

# Approximate single-thruster force in Newtons.
# 160 lbf is about 712 N.
ORION_CM_RCS_THRUST_N = 712.0

# Minimum pulse width used by the simple fixed-step firing logic.
# Since your sim runs at 0.25 s per step, the simplest first milestone
# is to let the minimum pulse equal one full simulation step.
ORION_CM_RCS_MIN_PULSE_S = ENTRY_DT_S

# Radius of the thruster shoulder ring from the capsule body origin.
CAPSULE_RCS_RING_RADIUS_M = 2.0

# Height of the shoulder ring above the body origin along +z_b.
CAPSULE_RCS_RING_Z_M = 1.2

# Four equally spaced pod azimuths around the body, radians.
# These are useful if you later want to generate the layout procedurally.
CAPSULE_RCS_POD_AZIMUTHS_RAD = (
    0.0,
    0.5 * math.pi,
    1.0 * math.pi,
    1.5 * math.pi,
)


# ============================================================================
# Rigid-body attitude and bank-angle control constants
# ============================================================================
# These are milestone-1 tuning constants for roll tracking.
# They are chosen to be readable and easy to tune.

# Sign convention used when converting body attitude into sigma.
# Keep this at +1.0 unless you discover your body-axis convention
# produces the opposite bank-angle sign.
SIGMA_SIGN = 1.0

# Numerical epsilon for normalization and safety checks.
EPS_NORM = 1.0e-12

# Simple rate and angle deadbands for roll control logic.
ROLL_SIGMA_DEADBAND_RAD = math.radians(0.5)
ROLL_RATE_DEADBAND_RAD_S = math.radians(0.25)

# First-pass body inertia values, kg m^2.
# Use a diagonal tensor for milestone 1.
# Izz is the main one used for roll-only control in the first stage.
CAPSULE_IXX_KGM2 = 12_000.0
CAPSULE_IYY_KGM2 = 12_000.0
CAPSULE_IZZ_KGM2 = 20_000.0

# Full inertia matrix in body coordinates.
CAPSULE_I_B_KGM2 = np.diag([
    CAPSULE_IXX_KGM2,
    CAPSULE_IYY_KGM2,
    CAPSULE_IZZ_KGM2,
])

# Roll controller gains for sigma tracking.
# tau_roll_cmd = kp * sigma_error - kd * roll_rate
ROLL_KP_NM_PER_RAD = 8_000.0
ROLL_KD_NM_PER_RAD_S = 4_000.0

# Saturation limit for commanded roll torque.
ROLL_CMD_MAX_TORQUE_NM = 6_000.0


# ============================================================================
# Optional bank-angle limits
# ============================================================================
# These are useful if you want one place to clamp commands.

SIGMA_MIN_RAD = -math.pi
SIGMA_MAX_RAD = math.pi


# ============================================================================
# Heat shield geometry model
# ============================================================================

class HeatShield:
    """
    Ring-discretized circular heat shield model.

    This class tracks per-ring heat flux and integrated heat load
    using the same Interval class throughout.

    Geometry assumptions
    --------------------
    - The shield is circular in the body frame.
    - It is split into concentric rings.
    - Each ring stores:
      * instantaneous heat flux qdot
      * accumulated heat load Q

    Thermal model assumptions
    -------------------------
    - Stagnation heating uses a Sutton-Graves style relation.
    - Off-stagnation heating is scaled by a radial shape factor.
    - This is a practical engineering approximation, not a full TPS solver.
    """

    def __init__(
        self,
        radius_m: float,
        nose_radius_m: float,
        num_rings: int = HEAT_SHIELD_NUM_RINGS,
        radial_exp: float = HEAT_SHIELD_RADIAL_EXP,
    ):
        self.radius = radius_m
        self.nose_radius = nose_radius_m
        self.num_rings = num_rings
        self.radial_exp = radial_exp

        # Ring edge radii from center to outer edge.
        self.r_edges = [
            i * radius_m / num_rings
            for i in range(num_rings + 1)
        ]

        # Ring center radii used for applying the radial heat shape.
        self.r_centers = [
            0.5 * (self.r_edges[i] + self.r_edges[i + 1])
            for i in range(num_rings)
        ]

        # Area of each ring annulus.
        self.areas = [
            math.pi * (self.r_edges[i + 1] ** 2 - self.r_edges[i] ** 2)
            for i in range(num_rings)
        ]

        # Initialize interval-valued thermal state.
        base_iv = promote(0.0)
        self.qdot = [type(base_iv)(0.0, 0.0) for _ in range(num_rings)]
        self.Q = [type(base_iv)(0.0, 0.0) for _ in range(num_rings)]

    def stagnation_qdot(self, rho, V):
        """
        Sutton-Graves style stagnation-point convective heating.

        Expected units
        --------------
        rho : atmospheric density
        V   : speed

        This function is written to stay compatible with your Interval math.
        """
        IV = type(rho)
        k = IV(1.83e-4, 1.83e-4)
        return k * (rho / self.nose_radius).sqrt() * V.pow_int(3)

    def radial_shape_factor(self, r: float) -> float:
        """
        Shape factor that reduces heating away from the centerline.

        Returns a scalar in [0, 1].
        """
        return max(
            0.0,
            (1.0 - (r / self.radius) ** 2) ** self.radial_exp
        )

    def update(self, rho, V, dt: float):
        """
        Update heat flux and integrated heat load for all rings.

        Parameters
        ----------
        rho : density-like interval quantity
        V   : speed-like interval quantity
        dt  : time step in seconds
        """
        qdot_stag = self.stagnation_qdot(rho, V)

        for i, r in enumerate(self.r_centers):
            shape = self.radial_shape_factor(r)

            # Build the shape factor using the same Interval type.
            IV = type(self.qdot[i])
            shape_iv = IV(shape, shape)

            self.qdot[i] = qdot_stag * shape_iv
            self.Q[i] = self.Q[i] + self.qdot[i] * dt

    def qdot_max(self):
        """
        Return an interval spanning the min and max ring heat flux.
        """
        return type(self.qdot[0])(
            min(q.lo for q in self.qdot),
            max(q.hi for q in self.qdot),
        )

    def Q_max(self):
        """
        Return an interval spanning the min and max accumulated ring heat load.
        """
        return type(self.Q[0])(
            min(Q.lo for Q in self.Q),
            max(Q.hi for Q in self.Q),
        )

    def qdot_mean(self):
        """
        Area-weighted mean heat flux over the shield.
        """
        total_area = sum(self.areas)
        weighted = type(self.qdot[0])(0.0, 0.0)

        for q, A in zip(self.qdot, self.areas):
            weighted = weighted + q * A

        return weighted / total_area