import constants
import math
import datetime
from interval_math import Interval, promote, box_add, box_scalar_mul, scalar_times_interval

# Hoist PyMSIS + numpy imports to module scope so get_atmosphere_density()
# doesn't re-resolve them on every step of the sim. We probe both the new
# `calculate` API and the legacy `msis.run` API at import time and pick
# whichever is available.
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


# Standard atmosphere defaults (avoid repeating in every call site).
SPEED_OF_SOUND_GAMMA = 1.4
R_AIR_DRY_J_KG_K = 287.053


def mach_from_V_T(V_mps, T_K, gamma=SPEED_OF_SOUND_GAMMA, R=R_AIR_DRY_J_KG_K):
    """Mach number from speed and temperature (perfect-gas dry air).
    Returns None when T is missing/non-positive (e.g. above the atmosphere)."""
    if T_K is None or not (T_K == T_K) or T_K <= 0:  # NaN-safe
        return None
    return float(V_mps) / math.sqrt(gamma * R * float(T_K))


def _pymsis_state(h_m, lat_deg, lon_deg, date, f107, f107a, ap):
    """Raw PyMSIS query -> (rho_kgm3, T_K). Index 0 = density, 10 = temperature."""
    if date is None:
        date = datetime.datetime(2026, 1, 1, 0, 0, 0)
    h_m = max(0.0, float(h_m))
    h_km = h_m / 1000.0
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


# --- Fast tabulated atmosphere -----------------------------------------------
# PyMSIS is the per-step bottleneck (~7 ms/query, called ~12-20x per guidance
# step). For RL training we precompute density + temperature on a fine altitude
# grid ONCE and interpolate (log-density vs altitude, accurate for the ~12
# orders-of-magnitude density span). This is ~100-1000x faster. The table is at
# a single nominal lat/lon/date, so it drops lat/lon atmospheric variation --
# acceptable for training; use full PyMSIS (table disabled) for final eval.
_ATM_TABLE = None  # dict: alt_grid_m, log_rho, T_K, alt_max_m


def build_atmosphere_table(alt_max_m=160_000.0, dz_m=100.0,
                           lat_deg=0.0, lon_deg=0.0, date=None,
                           f107=150.0, f107a=150.0, ap=7.0):
    """Precompute the (rho, T) altitude table and switch lookups to it.
    Returns the number of PyMSIS samples taken."""
    global _ATM_TABLE
    alts = _np.arange(0.0, float(alt_max_m) + float(dz_m), float(dz_m))
    rho = _np.empty_like(alts)
    T = _np.empty_like(alts)
    for i, a in enumerate(alts):
        r, t = _pymsis_state(float(a), lat_deg, lon_deg, date, f107, f107a, ap)
        rho[i] = max(float(r), 1.0e-30)
        T[i] = float(t) if t == t else 0.0  # NaN guard
    _ATM_TABLE = {
        "alt": alts, "log_rho": _np.log(rho), "T": T, "alt_max": float(alts[-1]),
    }
    return int(len(alts))


def disable_atmosphere_table():
    """Revert to full PyMSIS queries (high-fidelity, slow)."""
    global _ATM_TABLE
    _ATM_TABLE = None


def atmosphere_table_active() -> bool:
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
    Return (rho_kgm3, T_K) at the given geodetic altitude.

    If a precomputed table is active (build_atmosphere_table), interpolate it
    (log-density vs altitude) -- ~100-1000x faster than PyMSIS. Otherwise query
    PyMSIS directly. Above 1000 km: vacuum (rho=0, T=nan).
    """
    h_m = max(0.0, float(h_m))
    if _ATM_TABLE is not None:
        if h_m / 1000.0 > 1000.0:
            return 0.0, float("nan")
        tbl = _ATM_TABLE
        h = min(h_m, tbl["alt_max"])
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
    Return PyMSIS total mass density in kg/m^3. Backward-compatible
    wrapper around get_atmosphere_state() — kept because callers all over
    the simulator use it.
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
    """Get standard atmosphere values based on altitude."""
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
    pymsis_lo = get_atmosphere_density(z.hi, lat_deg, lon_deg, date, f107, f107a, ap)
    pymsis_hi = get_atmosphere_density(z.lo, lat_deg, lon_deg, date, f107, f107a, ap)
    
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
