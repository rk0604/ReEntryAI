"""
Atmosphere model for the entry simulator.

This module returns air density and temperature as a function of altitude. The
values come from the NRLMSIS empirical atmosphere through the pymsis package.
Two query paths are provided:

  1. A direct pymsis call for full fidelity. This path uses latitude,
     longitude, date, and space weather indices.
  2. A precomputed altitude table that interpolates density and temperature.
     The table path is much faster and is intended for reinforcement learning
     training, where a pymsis call on every step would dominate the runtime.

Interval valued helpers return a guaranteed density range over an altitude
interval. The interval supervisor uses that range to bound heating and load
factor.

Naming note: the function US_Standard_ATM is a historical name. It does not
return the 1976 US Standard Atmosphere. It returns NRLMSIS values through the
shared query path defined below.
"""

import constants
import math
import datetime
from interval_math import Interval, promote, box_add, box_scalar_mul, scalar_times_interval

# The pymsis and numpy imports are resolved once at module import so that
# get_atmosphere_density does not resolve them again on every simulation step.
# Two pymsis call styles exist across versions. The newer calculate function
# and the older msis.run function are both probed here, and whichever is found
# is recorded in _PYMSIS_BACKEND for later use.
import numpy as _np  # noqa: E402
try:
    from pymsis import calculate as _pymsis_calculate
    _PYMSIS_BACKEND = "calculate"
except ImportError:
    try:
        from pymsis import msis as _pymsis_msis
        _PYMSIS_BACKEND = "msis"
    except ImportError:
        _PYMSIS_BACKEND = None


# Shared gas constants, defined once so call sites do not repeat them.
SPEED_OF_SOUND_GAMMA = 1.4
R_AIR_DRY_J_KG_K = 287.053


def mach_from_V_T(V_mps, T_K, gamma=SPEED_OF_SOUND_GAMMA, R=R_AIR_DRY_J_KG_K):
    """
    Return the Mach number for a given speed and temperature.

    The speed of sound uses the perfect gas relation for dry air,
    a = sqrt(gamma * R * T), so the Mach number is V / a.

    Parameters:
        V_mps: speed in meters per second.
        T_K:   static temperature in kelvin.
        gamma: ratio of specific heats for air (default 1.4).
        R:     specific gas constant for dry air in J/(kg*K).

    Returns:
        The Mach number as a float. Returns None when the temperature is
        missing, not a number, or not positive (for example above the
        sensible atmosphere, where temperature is not defined).
    """
    # Reject a missing, NaN, or non positive temperature before the sqrt.
    if T_K is None or not (T_K == T_K) or T_K <= 0:
        return None
    return float(V_mps) / math.sqrt(gamma * R * float(T_K))


def _pymsis_state(h_m, lat_deg, lon_deg, date, f107, f107a, ap):
    """
    Query pymsis at one point and return density and temperature.

    Parameters:
        h_m:     geometric altitude in meters. Values below zero are clamped
                 to zero.
        lat_deg: latitude in degrees.
        lon_deg: longitude in degrees.
        date:    a datetime for the query. A fixed default date is used when
                 this is None.
        f107:    daily F10.7 solar flux index.
        f107a:   81 day averaged F10.7 index.
        ap:      geomagnetic activity index.

    Returns:
        A tuple (rho_kgm3, T_K). In the pymsis output, index 0 is total mass
        density and index 10 is temperature. Above 1000 km the air is treated
        as vacuum and the function returns zero density and a NaN temperature.

    The newer pymsis calculate function and the older msis.run function take
    slightly different arguments, so both call styles are handled here.
    """
    # A None date is replaced by a fixed reference date for reproducibility.
    if date is None:
        date = datetime.datetime(2026, 1, 1, 0, 0, 0)
    h_m = max(0.0, float(h_m))
    h_km = h_m / 1000.0
    # Above 1000 km treat the air as vacuum.
    if h_km > 1000.0:
        return 0.0, float("nan")
    if _PYMSIS_BACKEND == "calculate":
        out = _np.squeeze(_pymsis_calculate(
            date, lon_deg, lat_deg, h_km, f107, f107a, [[ap] * 7], version=0))
        return float(out[0]), float(out[10])
    if _PYMSIS_BACKEND == "msis":
        out = _np.squeeze(_pymsis_msis.run(
            date, lon_deg, lat_deg, h_km, f107, f107a, ap, version=0))
        return float(out[0]), float(out[10])
    raise RuntimeError("PyMSIS is missing. Please install it with: pip install pymsis")


# Fast tabulated atmosphere.
# A pymsis call is the per step bottleneck, on the order of 7 ms each and run
# about 12 to 20 times per guidance step. For training, density and temperature
# are precomputed once on a fine altitude grid and then interpolated. The
# interpolation uses the log of density against altitude, which stays accurate
# across the roughly 12 orders of magnitude that density spans. This is on the
# order of 100 to 1000 times faster. The table is built at one fixed latitude,
# longitude, and date, so it drops the latitude and longitude variation of the
# atmosphere. That is acceptable for training. For a final evaluation, disable
# the table and use the full pymsis path.
_ATM_TABLE = None  # dict with keys: alt, log_rho, T, alt_max


def build_atmosphere_table(alt_max_m=160_000.0, dz_m=100.0,
                           lat_deg=0.0, lon_deg=0.0, date=None,
                           f107=150.0, f107a=150.0, ap=7.0):
    """
    Precompute the density and temperature altitude table and switch all
    later lookups to use it.

    Parameters:
        alt_max_m: highest altitude in the table, in meters.
        dz_m:      altitude spacing between table points, in meters.
        lat_deg:   latitude used for every sample, in degrees.
        lon_deg:   longitude used for every sample, in degrees.
        date:      datetime used for every sample. A default is used when None.
        f107:      daily F10.7 solar flux index.
        f107a:     81 day averaged F10.7 index.
        ap:        geomagnetic activity index.

    Returns:
        The number of pymsis samples taken to build the table.

    Side effect: sets the module level table so get_atmosphere_state and
    get_atmosphere_density interpolate the table instead of calling pymsis.
    """
    global _ATM_TABLE
    alts = _np.arange(0.0, float(alt_max_m) + float(dz_m), float(dz_m))
    rho = _np.empty_like(alts)
    T = _np.empty_like(alts)
    for i, a in enumerate(alts):
        r, t = _pymsis_state(float(a), lat_deg, lon_deg, date, f107, f107a, ap)
        # Floor density so the later log is always defined.
        rho[i] = max(float(r), 1.0e-30)
        # Replace a NaN temperature with zero so the table stays numeric.
        T[i] = float(t) if t == t else 0.0
    # Store the log of density so lookups can interpolate it linearly.
    _ATM_TABLE = {
        "alt": alts, "log_rho": _np.log(rho), "T": T, "alt_max": float(alts[-1]),
    }
    return int(len(alts))


def disable_atmosphere_table():
    """
    Clear the precomputed table so later lookups call pymsis directly.

    Use this before a final evaluation to recover full fidelity, including the
    latitude and longitude variation that the table drops.
    """
    global _ATM_TABLE
    _ATM_TABLE = None


def atmosphere_table_active() -> bool:
    """
    Report whether the fast precomputed table is currently in use.

    Returns:
        True when a table has been built and not cleared, otherwise False.
    """
    return _ATM_TABLE is not None


def get_atmosphere_state(
    h_m,
    lat_deg=0.0,
    lon_deg=0.0,
    date=None,
    f107=150.0,
    f107a=150.0,
    ap=7.0,
):
    """
    Return air density and temperature at a geometric altitude.

    When a precomputed table is active, the result is interpolated from it
    using the log of density against altitude, which is far faster than a
    pymsis call. When no table is active, pymsis is queried directly. Above
    1000 km the air is treated as vacuum.

    Parameters:
        h_m:     geometric altitude in meters. Values below zero are clamped.
        lat_deg: latitude in degrees (ignored when the table is active).
        lon_deg: longitude in degrees (ignored when the table is active).
        date:    datetime for the query (ignored when the table is active).
        f107:    daily F10.7 solar flux index.
        f107a:   81 day averaged F10.7 index.
        ap:      geomagnetic activity index.

    Returns:
        A tuple (rho_kgm3, T_K). Density is in kg/m^3 and temperature is in
        kelvin. Temperature is NaN where it is not defined (for example above
        the sensible atmosphere).
    """
    h_m = max(0.0, float(h_m))
    if _ATM_TABLE is not None:
        # Above 1000 km the table does not apply; treat as vacuum.
        if h_m / 1000.0 > 1000.0:
            return 0.0, float("nan")
        tbl = _ATM_TABLE
        # Clamp to the top of the table so interpolation stays in range.
        h = min(h_m, tbl["alt_max"])
        # Interpolate the log of density, then exponentiate back to density.
        rho = math.exp(float(_np.interp(h, tbl["alt"], tbl["log_rho"])))
        T = float(_np.interp(h, tbl["alt"], tbl["T"]))
        return rho, (T if T > 0.0 else float("nan"))
    return _pymsis_state(h_m, lat_deg, lon_deg, date, f107, f107a, ap)


def get_atmosphere_density(
    h_m,
    lat_deg=0.0,
    lon_deg=0.0,
    date=None,
    f107=150.0,
    f107a=150.0,
    ap=7.0,
):
    """
    Return only the air density at a geometric altitude.

    This is a thin wrapper around get_atmosphere_state that keeps an older
    call signature. Many callers across the simulator use this name, so it is
    kept for compatibility.

    Parameters:
        Same as get_atmosphere_state.

    Returns:
        Density in kg/m^3 as a float.
    """
    rho, _T = get_atmosphere_state(h_m, lat_deg, lon_deg, date, f107, f107a, ap)
    return rho

def US_Standard_ATM(
    z: float,
    lat_deg=0.0,
    lon_deg=0.0,
    date=None,
    f107=150.0,
    f107a=150.0,
    ap=7.0,
) -> dict:
    """
    Return a small atmosphere dictionary at one altitude.

    Despite the name, this returns NRLMSIS values through
    get_atmosphere_state, not the 1976 US Standard Atmosphere. The name is
    kept because many call sites already use it.

    Parameters:
        z:       geometric altitude in meters.
        lat_deg: latitude in degrees.
        lon_deg: longitude in degrees.
        date:    datetime for the query.
        f107:    daily F10.7 solar flux index.
        f107a:   81 day averaged F10.7 index.
        ap:      geomagnetic activity index.

    Returns:
        A dictionary with keys T_K (kelvin), p_Pa (always None, since pressure
        is not produced by this path), rho_kgm3 (kg/m^3), geopotential_alt_m
        (always None), and layer_index (always minus one). The None and minus
        one entries exist so the return shape matches older callers.
    """
    rho, T_K = get_atmosphere_state(
        h_m=z,
        lat_deg=lat_deg,
        lon_deg=lon_deg,
        date=date,
        f107=f107,
        f107a=f107a,
        ap=ap,
    )

    return {
        "T_K": T_K,
        "p_Pa": None,
        "rho_kgm3": rho,
        "geopotential_alt_m": None,
        "layer_index": -1,
    }

def atmosphere_from_radius(r: float):
    """
    Return the atmosphere dictionary for a position given by radius.

    Parameters:
        r: distance from the center of the Earth in meters.

    Returns:
        The same dictionary as US_Standard_ATM. The radius is converted to a
        geometric altitude first, and any negative altitude is clamped to zero
        so the query stays valid at or below the surface.
    """
    h = constants.altitude(r)
    if h < 0:
        h = 0.0
    return US_Standard_ATM(h)

def intv_US_Standard_ATM(
    z: Interval,
    lat_deg=0.0,
    lon_deg=0.0,
    date=None,
    f107=150.0,
    f107a=150.0,
    ap=7.0,
) -> dict:
    """
    Return a guaranteed density range for an altitude interval.

    Density decreases as altitude increases, so the density is largest at the
    low end of the altitude interval and smallest at the high end. The two
    endpoints are queried and ordered into a single density interval that is
    guaranteed to contain the true density everywhere inside the altitude box.

    Parameters:
        z:       an altitude Interval in meters, with attributes lo and hi.
        lat_deg: latitude in degrees.
        lon_deg: longitude in degrees.
        date:    datetime for the query.
        f107:    daily F10.7 solar flux index.
        f107a:   81 day averaged F10.7 index.
        ap:      geomagnetic activity index.

    Returns:
        A dictionary keyed by a single layer index 0. Its value is a
        dictionary with T_K (None here), p_Pa (None here), and rho_kgm3 set to
        the density Interval. The layer keyed shape matches the layered
        atmosphere interface used elsewhere.
    """
    # Query both altitude endpoints. The high altitude gives the low density.
    pymsis_lo = get_atmosphere_density(z.hi, lat_deg, lon_deg, date, f107, f107a, ap)
    pymsis_hi = get_atmosphere_density(z.lo, lat_deg, lon_deg, date, f107, f107a, ap)

    # Order the two values into a valid interval regardless of which is larger.
    rho_interval = Interval(min(pymsis_lo, pymsis_hi), max(pymsis_lo, pymsis_hi))

    return {
        0: {
            "T_K": None,
            "p_Pa": None,
            "rho_kgm3": rho_interval,
        }
    }

if __name__ == "__main__":
    pass
