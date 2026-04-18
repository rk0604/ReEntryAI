"""
constants.py

Central location for Earth, atmosphere, simulation cadence, capsule,
RCS, heat shield, and predictor corrector tuning constants.

This version keeps the updated ring and sector heat shield model intact.
Only the additional constants needed by the heat aware milestone 1
predictor corrector are added here.
"""

import math
from typing import Optional

import numpy as np

from interval_math import Interval, promote


# Simulation cadence

# Fixed Euler integration step used by the main physics loop.
ENTRY_DT_S = 0.25

# Guidance updates more slowly than physics.
GUIDANCE_PERIOD_S = 1.0

# Number of physics steps per guidance update.
GUIDANCE_STEPS = int(round(GUIDANCE_PERIOD_S / ENTRY_DT_S))


# Earth and environment constants

# Mean Earth radius in meters.
RADIUS_EARTH = 6_378_100.0

# Earth standard gravitational parameter in m^3 / s^2.
MU = 3.986e14

# Sea level density in kg / m^3.
SEA_LEVEL_DENSITY = 1.225

# Simple atmosphere exponential scale height in meters.
H_M = 7200.0

# Gas constant for dry air in J / (kg K).
R_AIR_DRY = 287.053

# Universal gas constant in J / (mol K).
R_U = 8.314

# Standard gravity at Earth surface in m / s^2.
g_0 = 9.80665


# Atmosphere layer definitions

# Layer names mapped to layer index and lapse rate in K / m.
boundaries = {
    "troposphere":    (0, -0.0065),
    "tropopause":     (1,  0.0),
    "stratosphere_1": (2,  0.0010),
    "stratosphere_2": (3,  0.0028),
    "stratopause":    (4,  0.0),
    "mesosphere_1":   (5, -0.0028),
    "mesosphere_2":   (6, -0.0020),
}

# Same lapse rate information keyed by layer index.
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
    },
}


# Basic environment helper functions

def altitude(r: float) -> float:
    """
    Convert radial distance from Earth center to geometric altitude.
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


# Capsule geometry and aerodynamic reference values

# Crew module outer diameter in meters.
ORION_CM_DIAMETER_M = 5.0

# Outer radius in meters.
ORION_CM_RADIUS_M = 0.5 * ORION_CM_DIAMETER_M

# Heat shield radius in meters.
HEAT_SHIELD_RADIUS_M = ORION_CM_RADIUS_M

# Effective nose radius for stagnation heating in meters.
HEAT_SHIELD_NOSE_RADIUS_M = 1.0

# Radial exponent for heating concentration.
HEAT_SHIELD_RADIAL_EXP = 1.0

# Circular projected reference area in square meters.
CAPSULE_REFERENCE_AREA_M2 = math.pi * HEAT_SHIELD_RADIUS_M ** 2

# First pass capsule mass in kg.
CAPSULE_MASS_KG = 8500.0


# RCS layout constants

# Total number of crew module RCS thrusters in the milestone 1 layout.
CAPSULE_RCS_NUM_THRUSTERS = 12

# Approximate single thruster force in Newtons.
ORION_CM_RCS_THRUST_N = 712.0

# Minimum pulse width used by the fixed step firing logic.
ORION_CM_RCS_MIN_PULSE_S = ENTRY_DT_S

# Radius of the thruster shoulder ring from the capsule body origin.
CAPSULE_RCS_RING_RADIUS_M = 2.0

# Height of the shoulder ring above the body origin.
CAPSULE_RCS_RING_Z_M = 1.2

# Four equally spaced pod azimuths around the body in radians.
CAPSULE_RCS_POD_AZIMUTHS_RAD = (
    0.0,
    0.5 * math.pi,
    1.0 * math.pi,
    1.5 * math.pi,
)


# Rigid body attitude and bank angle control constants

# Sign convention used when converting body attitude into sigma.
SIGMA_SIGN = 1.0

# Numerical epsilon for normalization and safety checks.
EPS_NORM = 1.0e-12

# Simple rate and angle deadbands for roll control logic.
ROLL_SIGMA_DEADBAND_RAD = math.radians(0.5)
ROLL_RATE_DEADBAND_RAD_S = math.radians(0.25)

# First pass body inertia values in kg m^2.
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
ROLL_KP_NM_PER_RAD = 8_000.0
ROLL_KD_NM_PER_RAD_S = 4_000.0

# Saturation limit for commanded roll torque.
ROLL_CMD_MAX_TORQUE_NM = 6_000.0


# Optional bank angle limits

SIGMA_MIN_RAD = -math.pi
SIGMA_MAX_RAD = math.pi


# Heat aware predictor corrector constants

# This milestone 1 version uses a discrete candidate search rather than
# a literal Lu scalar corrector solve. The candidate search is easier to
# integrate cleanly into the current architecture while keeping the
# guidance to bank actuator to roll controller to RCS chain intact.

# Number of short horizon rollout steps evaluated per guidance update.
PREDICTOR_CORRECTOR_HORIZON_STEPS = 12

# Candidate bank magnitudes in degrees used by the short horizon search.
PREDICTOR_CANDIDATE_BANK_DEG = (
    40.0,
    45.0,
    50.0,
    55.0,
    60.0,
    65.0,
    70.0,
)

# Cost weights for the geometry part of the short horizon search.
# These defaults are mild and should be tuned after seeing real logs.
PREDICTOR_COST_WEIGHT_RANGE = 1.0e-5
PREDICTOR_COST_WEIGHT_HEADING = 25.0
PREDICTOR_COST_WEIGHT_CROSS_TRACK = 1.0e-5
PREDICTOR_COST_WEIGHT_HEAT = 1.0

# Large penalties used when a candidate is not heat feasible or when
# interval propagation becomes invalid.
PREDICTOR_HEAT_INFEASIBLE_PENALTY = 1.0e6
PREDICTOR_INVALID_INTERVAL_PENALTY = 1.0e7

# Default heat feasibility limits for milestone 1 guidance screening.
# These are conservative simulation defaults and not flight certified values.
HEAT_RATE_LIMIT_DEFAULT = 1.5e7
HEAT_LOAD_LIMIT_DEFAULT = 2.5e9

# Default csv paths used by the notebook logging.
SUCCESS_TRAJECTORY_CSV_DEFAULT = "trajectory_success.csv"
FAILED_HEAT_STEP_CSV_DEFAULT = "trajectory_failed_heat_steps.csv"
FAILED_HEAT_EPISODE_CSV_DEFAULT = "trajectory_failed_heat_episodes.csv"


# Heat shield geometry model

# Default heat shield discretization.
# This keeps the updated sector count already present in the pasted file.
HEAT_SHIELD_NUM_RINGS = 2
HEAT_SHIELD_NUM_SECTORS = 8

# Radial exponent controls how strongly heating drops from the center.
HEAT_SHIELD_RADIAL_EXP = 1.0

# Azimuthal gain controls directional heating contrast.
# This keeps the updated value already present in the pasted file.
HEAT_SHIELD_AZIMUTHAL_GAIN = 0.0


class HeatShield:
    """
    Ring and sector heat shield model.

    The shield is divided into concentric radial rings and equal angular sectors.
    Each cell stores interval valued instantaneous heat flux and accumulated
    heat load. This model is the conservative heat tracking layer used by the
    milestone 1 heat aware predictor corrector.

    This class keeps the current geometry and thermal state behavior intact.
    A small clone helper is added so candidate rollouts can copy the shield
    state cleanly without altering the live trajectory state.
    """

    def __init__(
        self,
        radius_m: float,
        nose_radius_m: float,
        num_rings: int = HEAT_SHIELD_NUM_RINGS,
        num_sectors: int = HEAT_SHIELD_NUM_SECTORS,
        radial_exp: float = HEAT_SHIELD_RADIAL_EXP,
        azimuthal_gain: float = HEAT_SHIELD_AZIMUTHAL_GAIN,
    ):
        if radius_m <= 0.0:
            raise ValueError("radius_m must be positive")
        if nose_radius_m <= 0.0:
            raise ValueError("nose_radius_m must be positive")
        if int(num_rings) <= 0:
            raise ValueError("num_rings must be positive")
        if int(num_sectors) <= 0:
            raise ValueError("num_sectors must be positive")
        if azimuthal_gain < 0.0 or azimuthal_gain > 1.0:
            raise ValueError("azimuthal_gain must be between 0 and 1")

        self.radius = float(radius_m)
        self.nose_radius = float(nose_radius_m)
        self.num_rings = int(num_rings)
        self.num_sectors = int(num_sectors)
        self.radial_exp = float(radial_exp)
        self.azimuthal_gain = float(azimuthal_gain)

        # Radial geometry.
        self.r_edges = [
            i * self.radius / self.num_rings
            for i in range(self.num_rings + 1)
        ]
        self.r_centers = [
            0.5 * (self.r_edges[i] + self.r_edges[i + 1])
            for i in range(self.num_rings)
        ]

        # Angular geometry.
        self.theta_edges = [
            j * 2.0 * math.pi / self.num_sectors
            for j in range(self.num_sectors + 1)
        ]
        self.theta_centers = [
            0.5 * (self.theta_edges[j] + self.theta_edges[j + 1])
            for j in range(self.num_sectors)
        ]

        # Flat cell metadata.
        self.cell_ring_index = []
        self.cell_sector_index = []
        self.cell_r_inner = []
        self.cell_r_outer = []
        self.cell_r_center = []
        self.cell_theta_lo = []
        self.cell_theta_hi = []
        self.cell_theta_center = []
        self.cell_area = []

        for ring_idx in range(self.num_rings):
            r_inner = self.r_edges[ring_idx]
            r_outer = self.r_edges[ring_idx + 1]
            r_center = 0.5 * (r_inner + r_outer)

            for sector_idx in range(self.num_sectors):
                theta_lo = self.theta_edges[sector_idx]
                theta_hi = self.theta_edges[sector_idx + 1]
                theta_center = 0.5 * (theta_lo + theta_hi)

                area = 0.5 * (r_outer ** 2 - r_inner ** 2) * (theta_hi - theta_lo)

                self.cell_ring_index.append(ring_idx)
                self.cell_sector_index.append(sector_idx)
                self.cell_r_inner.append(r_inner)
                self.cell_r_outer.append(r_outer)
                self.cell_r_center.append(r_center)
                self.cell_theta_lo.append(theta_lo)
                self.cell_theta_hi.append(theta_hi)
                self.cell_theta_center.append(theta_center)
                self.cell_area.append(area)

        self.areas = list(self.cell_area)

        self.reset_thermal_state()

    def reset_thermal_state(self) -> None:
        """
        Reset all interval valued heat flux and heat load states to zero.
        """
        base_iv = promote(0.0)
        self.qdot = [type(base_iv)(0.0, 0.0) for _ in range(len(self.cell_area))]
        self.Q = [type(base_iv)(0.0, 0.0) for _ in range(len(self.cell_area))]

    def clone(self) -> "HeatShield":
        """
        Return a deep copy of the heat shield state.

        This is useful for predictor corrector candidate evaluation because
        each candidate rollout should evolve an independent shield state.
        """
        new_obj = HeatShield(
            radius_m=self.radius,
            nose_radius_m=self.nose_radius,
            num_rings=self.num_rings,
            num_sectors=self.num_sectors,
            radial_exp=self.radial_exp,
            azimuthal_gain=self.azimuthal_gain,
        )

        new_obj.qdot = [type(iv)(iv.lo, iv.hi) for iv in self.qdot]
        new_obj.Q = [type(iv)(iv.lo, iv.hi) for iv in self.Q]
        return new_obj

    def stagnation_qdot(self, rho, V):
        """
        Sutton Graves style stagnation point convective heating.

        Float inputs are promoted to interval form so the class works with
        both nominal values and interval values.
        """
        rho_iv = promote(rho)
        V_iv = promote(V)
        iv_type = type(rho_iv)
        k = iv_type(1.83e-4, 1.83e-4)
        return k * (rho_iv / self.nose_radius).sqrt() * V_iv.pow_int(3)

    def radial_shape_factor(self, r: float) -> float:
        """
        Return the radial heating factor for a cell center.
        """
        return max(
            0.0,
            (1.0 - (r / self.radius) ** 2) ** self.radial_exp
        )

    def _wrap_angle_0_2pi(self, theta_rad: float) -> float:
        """
        Wrap an angle into the range from 0 to 2 pi.
        """
        two_pi = 2.0 * math.pi
        return theta_rad % two_pi

    def _angle_difference(self, a_rad: float, b_rad: float) -> float:
        """
        Return the smallest signed angular difference.
        """
        diff = (a_rad - b_rad + math.pi) % (2.0 * math.pi) - math.pi
        return diff

    def azimuthal_shape_factor(
        self,
        theta_rad: float,
        hot_theta_rad: float,
        azimuthal_gain: Optional[float] = None,
    ) -> float:
        """
        Return the angular heating factor for a cell center.

        hot_theta_rad is the direction of peak heating in the shield plane.
        """
        gain = self.azimuthal_gain if azimuthal_gain is None else float(azimuthal_gain)

        if gain < 0.0 or gain > 1.0:
            raise ValueError("azimuthal_gain must be between 0 and 1")

        theta_wrapped = self._wrap_angle_0_2pi(theta_rad)
        hot_wrapped = self._wrap_angle_0_2pi(hot_theta_rad)
        delta = self._angle_difference(theta_wrapped, hot_wrapped)

        return (1.0 - gain) + gain * 0.5 * (1.0 + math.cos(delta))

    def cell_index(self, ring_idx: int, sector_idx: int) -> int:
        """
        Convert a ring index and sector index into the flat storage index.
        """
        if ring_idx < 0 or ring_idx >= self.num_rings:
            raise IndexError("ring_idx out of range")
        if sector_idx < 0 or sector_idx >= self.num_sectors:
            raise IndexError("sector_idx out of range")
        return ring_idx * self.num_sectors + sector_idx

    def qdot_grid(self):
        """
        Return qdot as a ring by sector nested list.
        """
        return [
            [
                self.qdot[self.cell_index(ring_idx, sector_idx)]
                for sector_idx in range(self.num_sectors)
            ]
            for ring_idx in range(self.num_rings)
        ]

    def Q_grid(self):
        """
        Return Q as a ring by sector nested list.
        """
        return [
            [
                self.Q[self.cell_index(ring_idx, sector_idx)]
                for sector_idx in range(self.num_sectors)
            ]
            for ring_idx in range(self.num_rings)
        ]

    def update(
        self,
        rho,
        V,
        dt: float,
        hot_theta_rad: float = 0.0,
        azimuthal_gain: Optional[float] = None,
    ) -> None:
        """
        Update heat flux and accumulated heat load for every ring sector cell.
        """
        if dt <= 0.0:
            raise ValueError("dt must be positive")

        qdot_stag = self.stagnation_qdot(rho, V)

        for cell_idx in range(len(self.qdot)):
            r_center = self.cell_r_center[cell_idx]
            theta_center = self.cell_theta_center[cell_idx]

            radial_factor = self.radial_shape_factor(r_center)
            azimuthal_factor = self.azimuthal_shape_factor(
                theta_rad=theta_center,
                hot_theta_rad=hot_theta_rad,
                azimuthal_gain=azimuthal_gain,
            )

            total_shape = radial_factor * azimuthal_factor
            iv_type = type(self.qdot[cell_idx])
            shape_iv = iv_type(total_shape, total_shape)

            self.qdot[cell_idx] = qdot_stag * shape_iv
            self.Q[cell_idx] = self.Q[cell_idx] + self.qdot[cell_idx] * dt

    def qdot_max(self):
        """
        Return an interval spanning the minimum and maximum cell heat flux.
        """
        return type(self.qdot[0])(
            min(q.lo for q in self.qdot),
            max(q.hi for q in self.qdot),
        )

    def Q_max(self):
        """
        Return an interval spanning the minimum and maximum accumulated heat load.
        """
        return type(self.Q[0])(
            min(q.lo for q in self.Q),
            max(q.hi for q in self.Q),
        )

    def qdot_mean(self):
        """
        Return the area weighted mean heat flux over the shield face.
        """
        total_area = sum(self.areas)
        weighted = type(self.qdot[0])(0.0, 0.0)

        for q, area in zip(self.qdot, self.areas):
            weighted = weighted + q * area

        return weighted / total_area