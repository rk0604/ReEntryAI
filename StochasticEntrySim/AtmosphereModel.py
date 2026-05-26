import constants
import math
import datetime
from interval_math import Interval, promote, box_add, box_scalar_mul, scalar_times_interval

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
    Return PyMSIS total mass density in kg/m^3.
    h_m is geometric/geodetic altitude in meters.
    """
    if date is None:
        date = datetime.datetime(2026, 1, 1, 0, 0, 0)

    h_m = max(0.0, float(h_m))
    h_km = h_m / 1000.0

    if h_km > 1000.0:
        raise ValueError("Altitude must be <= 1000 km for PyMSIS.")

    try:
        from pymsis import calculate
        import numpy as np
        aps = [[ap] * 7]
        out = calculate(date, lon_deg, lat_deg, h_km, f107, f107a, aps, version=0)
        out = np.squeeze(out)
        return float(out[0])
    except ImportError:
        try:
            from pymsis import msis
            import numpy as np
            try:
                msis_out = msis.run(date, lon_deg, lat_deg, h_km, f107, f107a, ap, version=0)
            except Exception:
                msis_out = msis.run(date, lon_deg, lat_deg, h_km, f107, f107a, ap, version=0)
            return float(np.squeeze(msis_out)[0])
        except ImportError:
            raise RuntimeError("PyMSIS is missing. Please install it with: pip install pymsis")

def US_Standard_ATM(
    z: float,
    lat_deg=0.0,
    lon_deg=0.0,
    date=None,
    f107=150.0,
    f107a=150.0,
    ap=7.0,
) -> dict:
    rho = get_atmosphere_density(
        h_m=z,
        lat_deg=lat_deg,
        lon_deg=lon_deg,
        date=date,
        f107=f107,
        f107a=f107a,
        ap=ap,
    )

    return {
        "T_K": None,
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
