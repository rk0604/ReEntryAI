"""
jepa/testset_v3.py

The v3 closed-loop test set (jepa/v3_plan.md section 1): 16 bridge cases (the
frozen v1/v2 ICs) plus 48 new entries that vary altitude, speed, flight-path
angle, approach direction, downrange and crossrange widely, every one of them
reachable in all four test atmospheres. The target stays fixed (15 N, 15 E).

Reachability comes from physics only, never from a controller's result:
    1. For an entry state (h0, V0, gamma0), fly constant-bank trajectories
       (sigma = 0, +-5, +-10, ... deg) with the simulator's own equations of
       motion (math_3d.nominal_eom_step) and the true perturbed aero
       parameters of each atmosphere. The dynamics are rotationally symmetric
       (non-rotating spherical Earth, altitude-only atmosphere), so one sweep
       from (0 N, 0 E) heading east gives the footprint in (downrange,
       crossrange) for any entry position and heading.
    2. Only banks whose peak load stays within the study's 4.5 g limit count
       (steeper banks shorten the range but pull up to 10 g). The landing
       points of those banks, ordered by bank, bound the reachable region.
       Mixing banks (reversals, late lift-up) reaches at least this region,
       so it is a conservative inner estimate. Intersect the four
       atmospheres' regions and keep >= margin km inside.
    3. Choose the target's position inside that region (Sobol fractions), then
       rotate the whole picture so the target lands on 15 N, 15 E with the
       approach direction from the Sobol sequence. This gives the entry
       latitude, longitude and heading.

    python jepa/testset_v3.py --out configs/v3_testset.csv
"""

from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

import numpy as np

import constants
from jepa import data as D

R_KM = float(constants.RADIUS_EARTH) / 1000.0
G0 = 9.80665
H_DROGUE_M = 7600.0
ATMOSPHERES = {
    "nominal": {"density_scale": 1.00, "lift_scale": 1.00, "mass_scale": 1.0},
    "mismatch": {"density_scale": 0.88, "lift_scale": 0.88, "mass_scale": 1.0},
    "mismatch_high": {"density_scale": 1.12, "lift_scale": 1.12, "mass_scale": 1.0},
    "out_of_range": {"density_scale": 0.78, "lift_scale": 0.80, "mass_scale": 1.0},
}


@dataclass
class TestsetSpec:
    """Every number that defines the v3 test set (saved next to it)."""
    n_bridge: int = 16
    n_new: int = 48
    sobol_seed: int = 20261018
    altitude_m: Tuple[float, float] = (105_000.0, 135_000.0)
    velocity_mps: Tuple[float, float] = (7_200.0, 8_400.0)
    gamma_deg: Tuple[float, float] = (-4.0, -1.5)
    skip_guard_speed_mps: Tuple[float, float] = (7_600.0, 8_600.0)   # v2's guard: shallowest angle steepens with speed
    skip_guard_gamma_deg: float = -4.0
    downrange_frac: Tuple[float, float] = (0.10, 0.90)               # inside the reachable downrange band
    crossrange_frac: Tuple[float, float] = (-0.70, 0.70)             # of the reachable half-width at that downrange
    margin_km: float = 50.0
    max_entry_lat_deg: float = 30.0
    max_liftup_peak_g: float = 4.0
    load_limit_g: float = 4.5                                        # footprint banks must stay within this
    bank_step_deg: float = 5.0
    max_candidates: int = 4096
    target_phi_deg: float = 15.0
    target_lambda_deg: float = 15.0
    footprint_dt_s: float = 1.0


# ============================================================================
# Sphere geometry
# ============================================================================
def unit(lat: float, lon: float) -> np.ndarray:
    """Unit vector of a latitude/longitude (radians)."""
    return np.array([math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)])


def latlon(v: np.ndarray) -> Tuple[float, float]:
    """Latitude/longitude (radians) of a vector."""
    v = v / np.linalg.norm(v)
    return math.asin(max(-1.0, min(1.0, v[2]))), math.atan2(v[1], v[0])


def heading_basis(E: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Local east and north unit vectors at E."""
    east = np.cross([0.0, 0.0, 1.0], E)
    east /= np.linalg.norm(east)
    return east, np.cross(E, east)


def entry_frame(lat: float, lon: float, chi: float) -> np.ndarray:
    """
    Rotation Q whose columns are (entry, heading direction, entry x heading).
    Q maps the local frame (entry at 0 N 0 E, flying east) to the real one.
    chi is the simulator's heading convention: from east toward north.
    """
    E = unit(lat, lon)
    east, north = heading_basis(E)
    h = math.cos(chi) * east + math.sin(chi) * north
    return np.column_stack([E, h, np.cross(E, h)])


def to_local(Q: np.ndarray, lat: float, lon: float) -> Tuple[float, float]:
    """Real latitude/longitude -> (downrange km, crossrange km) in the entry frame."""
    p = Q.T @ unit(lat, lon)
    return R_KM * math.atan2(p[1], p[0]), R_KM * math.asin(max(-1.0, min(1.0, p[2])))


def _rotation_between(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Smallest rotation taking unit vector a to unit vector b (Rodrigues)."""
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-12:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


def _rotation_about(axis: np.ndarray, ang: float) -> np.ndarray:
    """Rotation by ang (radians) about a unit axis."""
    x, y, z = axis
    c, s = math.cos(ang), math.sin(ang)
    C = 1 - c
    return np.array([[c + x * x * C, x * y * C - z * s, x * z * C + y * s],
                     [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
                     [z * x * C - y * s, z * y * C + x * s, c + z * z * C]])


def place_entry(downrange_km: float, crossrange_km: float, approach_rad: float,
                target_lat: float, target_lon: float) -> Tuple[float, float, float]:
    """
    Entry (lat, lon, chi) such that the target sits at (downrange, crossrange)
    in the entry frame. approach_rad rotates the whole geometry about the
    target, so it sets the direction the capsule comes from.
    """
    t_local = unit(crossrange_km / R_KM, downrange_km / R_KM)
    T = unit(target_lat, target_lon)
    Q = _rotation_about(T, approach_rad) @ _rotation_between(t_local, T)
    E, h = Q[:, 0], Q[:, 1]
    lat, lon = latlon(E)
    east, north = heading_basis(E / np.linalg.norm(E))
    return lat, lon, math.atan2(float(h @ north), float(h @ east))


# ============================================================================
# Footprints
# ============================================================================
def truth_params(domain: Dict[str, float]) -> dict:
    """The nominal aero/mass params with one atmosphere's truth applied (as the env does)."""
    from jepa.evaluate import nominal_params
    return D.apply_domain(dict(nominal_params()), domain)


def fly_constant_bank(h0: float, V0: float, gamma0: float, sigma: float, params: dict,
                      dt: float = 1.0, max_s: float = 2000.0) -> Dict[str, float]:
    """
    Fly a constant bank from (0 N, 0 E) heading east down to drogue altitude.
    Returns landing (downrange_km, crossrange_km), time, peak g, and whether
    it landed (False for a skip-out above 160 km or a timeout).
    """
    from math_3d import nominal_aero_forces_from_state, nominal_eom_step
    Re = float(constants.RADIUS_EARTH)
    x = [h0 + Re, 0.0, 0.0, V0, gamma0, 0.0]
    m = float(params["mass_kg"])
    gmax, t = 0.0, 0.0
    while t < max_s:
        a = nominal_aero_forces_from_state(x=x, sigma_rad=sigma, params=params)
        gmax = max(gmax, math.hypot(a["drag_mag_N"], a["lift_mag_N"]) / m / G0)
        xn = nominal_eom_step(x=x, sigma_rad=sigma, params=params, dt_s=dt)["x_next"]
        t += dt
        if xn[0] - Re <= H_DROGUE_M:
            f = (x[0] - Re - H_DROGUE_M) / max(x[0] - xn[0], 1e-6)
            lat, lon = x[1] + f * (xn[1] - x[1]), x[2] + f * (xn[2] - x[2])
            Q = np.eye(3)
            dr, cr = to_local(Q, lat, lon)
            return {"downrange_km": dr, "crossrange_km": cr, "t_s": t, "peak_g": gmax, "landed": True}
        if xn[0] - Re > D.SKIP_OUT_ALT_M:
            break
        x = xn
    return {"downrange_km": float("nan"), "crossrange_km": float("nan"), "t_s": t, "peak_g": gmax, "landed": False}


def footprint(h0: float, V0: float, gamma0: float, params: dict, dt: float = 1.0,
              step_deg: float = 5.0, g_cap: float = 4.5) -> List[dict]:
    """
    Constant-bank landing points ordered by bank. The sweep walks out from
    lift-up (0 deg) on each side and stops after the first bank whose peak g
    exceeds g_cap (that point is kept, flagged within_load False, for plots).
    Peak g grows with |bank|, so the within-load points are contiguous.
    """
    def fly(s):
        return {"sigma_deg": s, **fly_constant_bank(h0, V0, gamma0, math.radians(s), params, dt)}

    def flag(p):
        p["within_load"] = bool(p["landed"] and p["peak_g"] <= g_cap)
        return p

    pts = [flag(fly(0.0))]
    n_side = int(round(90.0 / step_deg))
    for sign in (-1.0, 1.0):
        if not pts[0]["within_load"]:
            break
        for k in range(1, n_side + 1):
            p = flag(fly(sign * k * step_deg))
            pts.append(p)
            if not p["within_load"]:
                break
    return sorted(pts, key=lambda p: p["sigma_deg"])


def _polygon(fp: Sequence[dict]) -> np.ndarray:
    """Closed boundary of the reachable region: within-load landing points in bank order."""
    return np.array([[p["downrange_km"], p["crossrange_km"]] for p in fp if p["within_load"]]).reshape(-1, 2)


def _inside(poly: np.ndarray, pt: np.ndarray) -> bool:
    """Even-odd point-in-polygon test."""
    x, y = pt
    inside = False
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1 + 1e-300) + x1:
            inside = not inside
    return inside


def margin_km(fp: Sequence[dict], pt: Tuple[float, float]) -> float:
    """Signed distance from pt to the reachable region's boundary (km, > 0 inside)."""
    poly = _polygon(fp)
    if len(poly) < 3 or not np.isfinite(poly).all():
        return float("-inf")
    p = np.asarray(pt, dtype=float)
    d = math.inf
    for i in range(len(poly)):
        a, b = poly[i], poly[(i + 1) % len(poly)]
        ab = b - a
        t = float(np.clip(np.dot(p - a, ab) / max(np.dot(ab, ab), 1e-12), 0.0, 1.0))
        d = min(d, float(np.linalg.norm(p - (a + t * ab))))
    return d if _inside(poly, p) else -d


def _crossrange_halfwidth(fp: Sequence[dict], x: float) -> float:
    """Largest symmetric crossrange half-width of the region at downrange x."""
    poly = _polygon(fp)
    if len(poly) < 3:
        return 0.0
    ys = []
    for i in range(len(poly)):
        (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % len(poly)]
        if (x1 - x) * (x2 - x) <= 0 and x1 != x2:
            ys.append(y1 + (x - x1) * (y2 - y1) / (x2 - x1))
    if len(ys) < 2:
        return 0.0
    return max(0.0, min(max(ys), -min(ys)))


def _max_bank(fp: Sequence[dict]) -> float:
    """Largest |bank| (deg) that both sides of the sweep can fly within the load limit."""
    ok = [p["sigma_deg"] for p in fp if p["within_load"]]
    return min(max(ok, default=0.0), -min(ok, default=0.0))


def _shallowest_gamma(V: float, spec: TestsetSpec) -> float:
    """v2's skip-out guard: the shallowest allowed flight-path angle steepens with speed."""
    v0, v1 = spec.skip_guard_speed_mps
    frac = min(max((V - v0) / (v1 - v0), 0.0), 1.0)
    return spec.gamma_deg[1] + frac * (spec.skip_guard_gamma_deg - spec.gamma_deg[1])


# ============================================================================
# The test set
# ============================================================================
def _candidate(u: np.ndarray, spec: TestsetSpec, params: Dict[str, dict]) -> Tuple[Optional[dict], str]:
    """One Sobol point -> an accepted case (with its footprints) or a rejection reason."""
    h0 = spec.altitude_m[0] + u[0] * (spec.altitude_m[1] - spec.altitude_m[0])
    V0 = spec.velocity_mps[0] + u[1] * (spec.velocity_mps[1] - spec.velocity_mps[0])
    g_hi = _shallowest_gamma(V0, spec)
    gam = spec.gamma_deg[0] + u[2] * (g_hi - spec.gamma_deg[0])
    approach = 2.0 * math.pi * u[3]
    fps = {tag: footprint(h0, V0, math.radians(gam), p, spec.footprint_dt_s, spec.bank_step_deg, spec.load_limit_g)
           for tag, p in params.items()}
    up = {tag: next(pt for pt in fp if pt["sigma_deg"] == 0) for tag, fp in fps.items()}
    if not all(p["landed"] for p in up.values()):
        return None, "skip-out or no landing at lift-up"
    if max(p["peak_g"] for p in up.values()) > spec.max_liftup_peak_g:
        return None, "lift-up peak g above limit"
    lo = max(min(pt["downrange_km"] for pt in fp if pt["within_load"]) for fp in fps.values()) + spec.margin_km
    hi = min(up[tag]["downrange_km"] for tag in fps) - spec.margin_km
    if hi <= lo:
        return None, "empty downrange intersection"
    df = spec.downrange_frac[0] + u[4] * (spec.downrange_frac[1] - spec.downrange_frac[0])
    dr = lo + df * (hi - lo)
    w = min(_crossrange_halfwidth(fp, dr) for fp in fps.values()) - spec.margin_km
    cf = spec.crossrange_frac[0] + u[5] * (spec.crossrange_frac[1] - spec.crossrange_frac[0])
    cr = cf * max(w, 0.0)
    margins = {tag: margin_km(fp, (dr, cr)) for tag, fp in fps.items()}
    if min(margins.values()) < spec.margin_km:
        return None, "margin below limit"
    tphi, tlam = math.radians(spec.target_phi_deg), math.radians(spec.target_lambda_deg)
    lat, lon, chi = place_entry(dr, cr, approach, tphi, tlam)
    back = to_local(entry_frame(lat, lon, chi), tphi, tlam)
    assert abs(back[0] - dr) < 0.1 and abs(back[1] - cr) < 0.1, (back, dr, cr)
    if abs(math.degrees(lat)) > spec.max_entry_lat_deg:
        return None, "entry latitude outside training range"
    row = {"altitude_m": h0, "velocity_mps": V0, "gamma_deg": gam, "chi_deg": math.degrees(chi),
           "phi_deg": math.degrees(lat), "lambda_deg": math.degrees(lon), "downrange_km": dr, "crossrange_km": cr,
           "approach_deg": math.degrees(approach)}
    for tag in fps:
        row[f"margin_km_{tag}"] = margins[tag]
        row[f"feasible_{tag}"] = margins[tag] > 0.0
        row[f"liftup_peak_g_{tag}"] = up[tag]["peak_g"]
        row[f"liftup_time_s_{tag}"] = up[tag]["t_s"]
        row[f"max_bank_deg_{tag}"] = _max_bank(fps[tag])
    return {"row": row, "footprints": fps}, "accepted"


def _bridge_row(ic: dict, spec: TestsetSpec, params: Dict[str, dict]) -> dict:
    """Footprint margins of a frozen v1/v2 IC (not filtered: it is kept as is)."""
    Q = entry_frame(math.radians(ic["phi_deg"]), math.radians(ic["lambda_deg"]), math.radians(ic["chi_deg"]))
    dr, cr = to_local(Q, math.radians(spec.target_phi_deg), math.radians(spec.target_lambda_deg))
    fps = {tag: footprint(ic["altitude_m"], ic["velocity_mps"], math.radians(ic["gamma_deg"]), p,
                          spec.footprint_dt_s, spec.bank_step_deg, spec.load_limit_g)
           for tag, p in params.items()}
    row = dict(ic, downrange_km=dr, crossrange_km=cr, approach_deg=float("nan"))
    for tag, fp in fps.items():
        up = next(pt for pt in fp if pt["sigma_deg"] == 0)
        m = margin_km(fp, (dr, cr))
        row[f"margin_km_{tag}"] = m
        row[f"feasible_{tag}"] = m > 0.0
        row[f"liftup_peak_g_{tag}"] = up["peak_g"]
        row[f"liftup_time_s_{tag}"] = up["t_s"]
        row[f"max_bank_deg_{tag}"] = _max_bank(fp)
    return {"row": row, "footprints": fps}


def _candidate_job(job: dict) -> Tuple[Optional[dict], str]:
    """Pool worker: evaluate one Sobol candidate."""
    spec = TestsetSpec(**job["spec"])
    params = {tag: truth_params(dom) for tag, dom in ATMOSPHERES.items()}
    return _candidate(np.asarray(job["u"]), spec, params)


def build_testset(spec: Optional[TestsetSpec] = None, workers: int = 1, bridge_csv: Optional[str] = None,
                  batch: int = 64, log=print) -> dict:
    """
    Build the 64-case test set. Returns {"cases": DataFrame, "footprints":
    DataFrame, "rejections": {reason: count}, "spec": dict}.
    """
    import pandas as pd
    import torch
    spec = spec or TestsetSpec()
    params = {tag: truth_params(dom) for tag, dom in ATMOSPHERES.items()}
    rows, fps_rows = [], []
    bridge_csv = bridge_csv or str(SIM_DIR / "configs" / "heldout_testset.csv")
    for i, ic in enumerate(D.load_test_ics(bridge_csv, n=spec.n_bridge)):
        b = _bridge_row(ic, spec, params)
        rows.append({"case": i, "set": "bridge", "sobol_index": -1, **b["row"]})
        fps_rows += [{"case": i, "atmosphere": t, **pt} for t, fp in b["footprints"].items() for pt in fp]
    sob = torch.quasirandom.SobolEngine(dimension=6, scramble=True, seed=spec.sobol_seed)
    pts = sob.draw(spec.max_candidates).double().numpy()
    rejections: Dict[str, int] = {}
    accepted, k = 0, 0
    pool = mp.get_context("spawn").Pool(int(workers)) if workers > 1 else None
    try:
        while accepted < spec.n_new and k < len(pts):
            chunk = list(range(k, min(k + batch, len(pts))))
            k = chunk[-1] + 1
            jobs = [{"u": pts[j].tolist(), "spec": asdict(spec)} for j in chunk]
            res = pool.map(_candidate_job, jobs) if pool else [_candidate_job(j) for j in jobs]
            for j, (cand, why) in zip(chunk, res):
                if accepted >= spec.n_new:
                    k = j          # candidates after the 48th acceptance are not counted
                    break
                rejections[why] = rejections.get(why, 0) + 1
                if cand is None:
                    continue
                case = spec.n_bridge + accepted
                rows.append({"case": case, "set": "new", "sobol_index": j, **cand["row"]})
                fps_rows += [{"case": case, "atmosphere": t, **pt} for t, fp in cand["footprints"].items() for pt in fp]
                accepted += 1
            log(f"[testset] {accepted}/{spec.n_new} new cases from {k} Sobol candidates")
    finally:
        if pool:
            pool.close()
            pool.join()
    if accepted < spec.n_new:
        raise RuntimeError(f"only {accepted} feasible cases in {len(pts)} candidates; widen the spec")
    cases = pd.DataFrame(rows)
    rejections.pop("accepted", None)
    return {"cases": cases, "footprints": pd.DataFrame(fps_rows), "rejections": rejections, "spec": asdict(spec),
            "n_candidates": k}


def planning_horizon_steps(cases, macro_dt_s: float = 5.0, pad_s: float = 60.0) -> int:
    """One MPPI horizon for every learned controller: the longest lift-up flight + pad (plan 4.2)."""
    t = max(cases[f"liftup_time_s_{tag}"].max() for tag in ATMOSPHERES)
    return int(math.ceil((t + pad_s) / macro_dt_s))


def load_cases(path: str) -> List[dict]:
    """The frozen test set as a list of dicts (ic fields + metadata)."""
    import pandas as pd
    return pd.read_csv(path).to_dict("records")


IC_FIELDS = ("altitude_m", "velocity_mps", "gamma_deg", "chi_deg", "phi_deg", "lambda_deg")


def ic_of(case: dict) -> dict:
    """The ic_override dict of a test-set row."""
    return {k: float(case[k]) for k in IC_FIELDS}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(SIM_DIR / "configs" / "v3_testset.csv"))
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    r = build_testset(workers=a.workers)
    r["cases"].to_csv(a.out, index=False)
    r["footprints"].to_csv(a.out.replace(".csv", "_footprints.csv"), index=False)
    with open(a.out.replace(".csv", "_spec.json"), "w") as f:
        json.dump({"spec": r["spec"], "rejections": r["rejections"], "n_candidates": r["n_candidates"],
                   "horizon_steps": planning_horizon_steps(r["cases"])}, f, indent=1)
    print(r["cases"][["case", "set", "altitude_m", "velocity_mps", "gamma_deg", "chi_deg", "phi_deg", "lambda_deg",
                      "downrange_km", "crossrange_km", "margin_km_out_of_range"]].round(2).to_string())
    print("rejections:", r["rejections"], "horizon steps:", planning_horizon_steps(r["cases"]))
