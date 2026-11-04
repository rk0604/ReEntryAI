"""
jepa/benchmark_v5.py

The frozen v5 benchmark (jepa/v5_plan.md section 3). Three suites live in
configs/benchmark_v5/ and are verified against its SHA256SUMS every time they load:

    A256    configs/v3_testset.csv (64 cases) x testset_v3.ATMOSPHERES = 256 conditions,
            unchanged; A256.json records the file's checksum and the atmospheres
    ORB-1   the standard scenario (Orion-class LEO return): Lu (2014) CEV entry state with
            gamma0 = -1.40 deg, three targets on the entry great circle at 25/50/75 % of the
            nominal 90-to-0 deg constant-bank footprint (measured in the full simulator),
            Bairstow (2006) Table C.1 dispersions (EI speed 15 m/s 3-sigma, decision D4),
            3,000 cases in three blocks of 1,000
    DEV48   12 new testset_v3 cases (Sobol seed 20261104, A256's Sobol indices skipped) x the
            four atmospheres: the only suite any baseline is tuned on

score() turns a flights table into every section-5 metric per arm x suite x (atmosphere or
target range): success at 2/5/10/25 km, log miss, the Orion EM-1 requirement test (one-sided
90 % Clopper-Pearson lower bound on success@10km >= 0.99865), NASA-STD-3000 g-duration
compliance (Bairstow App. D.2 procedure), authority saturation, compute per decision.

    python jepa/benchmark_v5.py --build --workers 11     # writes configs/benchmark_v5/ (then frozen)
    python jepa/benchmark_v5.py --verify
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import multiprocessing as mp
import os
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

import numpy as np

import constants
from jepa import testset_v3 as T3
from jepa.data import NOMINAL_DOMAIN

BENCH_DIR = SIM_DIR / "configs" / "benchmark_v5"
BENCH_REL = "configs/benchmark_v5"
A256_CSV = SIM_DIR / "configs" / "v3_testset.csv"
SUITES = ("A256", "ORB1", "DEV48")
FLOAT_FMT = "%.6f"
R_KM = float(constants.RADIUS_EARTH) / 1000.0
FT = 0.3048
# NASA-STD-3000 duration-based g limits, deconditioned, ill or injured crew (Bairstow 2006 Table D.3).
NASA_STD_3000_DECONDITIONED = ((0.5, 14.0), (8.0, 10.0), (30.0, 7.8), (100.0, 4.0), (1000.0, 4.0))
EM1_REQUIREMENT = 0.99865          # P(miss <= 10 km), McNamara & Braun AAS 14-094, 3-sigma
EM1_CONFIDENCE = 0.90
TRAIN_DOMAIN = {"density_scale": (0.70, 1.30), "lift_scale": (0.70, 1.30), "mass_scale": (0.90, 1.10)}


class ChecksumError(RuntimeError):
    """A benchmark file does not match configs/benchmark_v5/SHA256SUMS."""


# ============================================================================
# ORB-1 specification
# ============================================================================
@dataclass(frozen=True)
class ORB1Spec:
    """Every number that defines ORB-1 (v5_plan.md section 3.2); saved as ORB1_spec.json."""
    n_cases: int = 3000
    n_blocks: int = 3
    seed: int = 20261103
    altitude_m: float = 121_920.0           # 400,000 ft (Lu 2014, Bairstow 2006, Orion)
    velocity_mps: float = 7_623.5           # Lu 2014 Table 1, CEV orbital entry
    gamma_deg: float = -1.40                # D3
    phi_deg: float = 0.0
    lambda_deg: float = 0.0
    chi_deg: float = 45.99                  # the A256 corridor (toward 15 N, 15 E)
    target_fracs: Tuple[float, ...] = (0.25, 0.50, 0.75)
    sigma3_fpa_deg: float = 0.10            # Bairstow Table C.1
    sigma3_azimuth_deg: float = 0.10
    sigma3_speed_mps: float = 15.0          # D4 (Bairstow's 300 m/s is lunar-return geometry)
    sigma3_mass_frac: float = 0.05
    sigma3_density_frac: float = 0.30
    sigma3_cd_frac: float = 0.05            # hypersonic value, applied at every Mach
    sigma3_cl_frac: float = 0.05
    truncate_sigma: float = 3.0
    n_seeds: int = 5                        # learned arms: case i flies training seed i mod 5
    footprint_dt_s: float = 1.0
    footprint_step_deg: float = 5.0
    load_limit_g: float = 4.5
    target_round_km: float = 0.1            # target downranges are stored to this precision

    def nominal_ic(self) -> dict:
        return {"altitude_m": self.altitude_m, "velocity_mps": self.velocity_mps, "gamma_deg": self.gamma_deg,
                "chi_deg": self.chi_deg, "phi_deg": self.phi_deg, "lambda_deg": self.lambda_deg}


def target_index(case: int, spec: ORB1Spec = ORB1Spec()) -> int:
    """Target range of an ORB-1 case: cycles 0, 1, 2 inside each block (334/333/333 per block)."""
    return (int(case) % (spec.n_cases // spec.n_blocks)) % len(spec.target_fracs)


def block_of(case: int, spec: ORB1Spec = ORB1Spec()) -> int:
    return int(case) // (spec.n_cases // spec.n_blocks)


def truncated_normals(rng: np.random.Generator, k: int, cut: float) -> np.ndarray:
    """k independent standard normals truncated at +-cut, by rejection, one variable at a time."""
    out = np.empty(k)
    for i in range(k):
        z = rng.standard_normal()
        while abs(z) > cut:
            z = rng.standard_normal()
        out[i] = z
    return out


ORB1_DRAWS = ("fpa", "azimuth", "speed", "mass", "density", "cd", "cl")


def orb1_dispersions(spec: ORB1Spec, n: Optional[int] = None) -> List[dict]:
    """The dispersed EI state and truth multipliers of the first n cases (all if None), in case order."""
    rng = np.random.default_rng(spec.seed)
    rows = []
    for i in range(spec.n_cases if n is None else int(n)):
        z = dict(zip(ORB1_DRAWS, truncated_normals(rng, len(ORB1_DRAWS), spec.truncate_sigma)))
        rho = 1.0 + spec.sigma3_density_frac / 3.0 * z["density"]
        cd = 1.0 + spec.sigma3_cd_frac / 3.0 * z["cd"]
        cl = 1.0 + spec.sigma3_cl_frac / 3.0 * z["cl"]
        m = 1.0 + spec.sigma3_mass_frac / 3.0 * z["mass"]
        rows.append({
            "case": i, "block": block_of(i, spec), "target_index": target_index(i, spec),
            "train_seed": i % spec.n_seeds,
            "altitude_m": spec.altitude_m, "velocity_mps": spec.velocity_mps + spec.sigma3_speed_mps / 3.0 * z["speed"],
            "gamma_deg": spec.gamma_deg + spec.sigma3_fpa_deg / 3.0 * z["fpa"],
            "chi_deg": spec.chi_deg + spec.sigma3_azimuth_deg / 3.0 * z["azimuth"],
            "phi_deg": spec.phi_deg, "lambda_deg": spec.lambda_deg,
            "rho_mult": rho, "cd_mult": cd, "cl_mult": cl, "mass_mult": m,
            # data.apply_domain knobs: drag ~ rho S CD c_D and lift ~ rho S CL c_L exactly
            "density_scale": rho * cd, "lift_scale": cl / cd, "mass_scale": m,
        })
    return rows


def in_training_domain(row: dict) -> bool:
    return all(lo <= float(row[k]) <= hi for k, (lo, hi) in TRAIN_DOMAIN.items())


def domain_of(row: dict) -> Dict[str, float]:
    return {k: float(row[k]) for k in ("density_scale", "lift_scale", "mass_scale")}


# ============================================================================
# Full-simulator constant-bank flights (gate S0) and target placement
# ============================================================================
def fly_constant_bank_sim(ic: dict, sigma_deg: float, domain: Optional[dict] = None,
                          config_path: str = "configs/default.json", record: bool = False) -> dict:
    """
    One full-simulator flight (rl_env policy mode: actuator, roll control, RCS) commanding a constant
    bank from the entry interface to the drogue. Returns the landing point's along-track downrange and
    crossrange in the entry frame, time, peak g, outcome, and (record=True) data.physical_state rows.
    """
    import rl_env
    from jepa.data import SIGMA_MAX_RAD, apply_domain, physical_state
    env = rl_env.OrionEntryEnv(config_path=config_path, controller_mode="policy", dispersion=False,
                               max_episode_steps=1400, fast_atmosphere=True)
    apply_domain(env.params, dict(domain or NOMINAL_DOMAIN))
    env.reset(seed=0, options={"ic_override": dict(ic)})
    a = np.array([math.radians(sigma_deg) / SIGMA_MAX_RAD], dtype=np.float32)
    ys = [physical_state(env)] if record else []
    done, info = False, {}
    while not done:
        _, _, term, trunc, info = env.step(a)
        if record:
            ys.append(physical_state(env))
        done = bool(term or trunc)
    Q = T3.entry_frame(math.radians(ic["phi_deg"]), math.radians(ic["lambda_deg"]), math.radians(ic["chi_deg"]))
    dr, cr = T3.to_local(Q, float(env.x[1]), float(env.x[2]))
    out = {"sigma_deg": float(sigma_deg), "downrange_km": dr, "crossrange_km": cr, "t_s": float(env.t_s),
           "peak_g": float(info.get("peak_g", float("nan"))), "outcome": str(info.get("outcome")),
           "end_altitude_m": float(env.x[0] - constants.RADIUS_EARTH), "end_gamma_deg": math.degrees(float(env.x[4]))}
    if record:
        out["y"] = np.asarray(ys)
    return out


def great_circle_point(ic: dict, downrange_km: float) -> Tuple[float, float]:
    """(lat, lon) in degrees of the point downrange_km along the entry great circle (zero crossrange)."""
    Q = T3.entry_frame(math.radians(ic["phi_deg"]), math.radians(ic["lambda_deg"]), math.radians(ic["chi_deg"]))
    lat, lon = T3.latlon(Q @ T3.unit(0.0, downrange_km / R_KM))
    return math.degrees(lat), math.degrees(lon)


def orb1_targets(d90_km: float, d0_km: float, spec: ORB1Spec = ORB1Spec()) -> List[dict]:
    """The three target ranges (rounded) and their positions on the nominal entry great circle."""
    out = []
    for k, f in enumerate(spec.target_fracs):
        d = round((d90_km + f * (d0_km - d90_km)) / spec.target_round_km) * spec.target_round_km
        lat, lon = great_circle_point(spec.nominal_ic(), d)
        out.append({"target_index": k, "frac": f, "downrange_km": round(d, 3), "phi_deg": round(lat, 9),
                    "lambda_deg": round(lon, 9), "config_path": f"{BENCH_REL}/orb1_target{k}.json"})
    return out


def target_config(target: dict) -> dict:
    """configs/default.json with only the mission target (and the run labels) changed."""
    cfg = json.loads((SIM_DIR / "configs" / "default.json").read_text())
    cfg["run_id"] = f"orb1_target{target['target_index']}"
    cfg["run_name"] = (f"v5 benchmark ORB-1 target {target['target_index']} "
                       f"({target['frac']:.0%} of footprint, {target['downrange_km']:.1f} km)")
    cfg["mission"]["target_phi_deg"] = float(target["phi_deg"])
    cfg["mission"]["target_lambda_deg"] = float(target["lambda_deg"])
    return cfg


# ============================================================================
# Footprint margins (testset_v3 code) for ORB-1 cases, in a process pool
# ============================================================================
def _orb1_margin_job(job: dict) -> dict:
    """
    Footprint margins of one case under its true parameters (testset_v3 code). The targets are placed on the
    0-90 deg constant-bank footprint, so `margin_km` uses that footprint with no load cap; `margin_km_4p5g`
    keeps only the banks whose peak g stays within the 4.5 g limit (testset_v3's A256 rule), a conservative
    bound on what a load-limited guidance can reach (load relief varies the bank, a constant bank cannot).
    """
    row, tgt, spec = job["row"], job["target"], ORB1Spec(**job["spec"])
    params = T3.truth_params(domain_of(row))
    fp = T3.footprint(row["altitude_m"], row["velocity_mps"], math.radians(row["gamma_deg"]), params,
                      spec.footprint_dt_s, spec.footprint_step_deg, float("inf"))
    capped = [dict(p, within_load=bool(p["landed"] and p["peak_g"] <= spec.load_limit_g)) for p in fp]
    Q = T3.entry_frame(math.radians(row["phi_deg"]), math.radians(row["lambda_deg"]), math.radians(row["chi_deg"]))
    dr, cr = T3.to_local(Q, math.radians(tgt["phi_deg"]), math.radians(tgt["lambda_deg"]))
    up = next(p for p in fp if p["sigma_deg"] == 0)
    return {"case": row["case"], "downrange_km": dr, "crossrange_km": cr, "margin_km": T3.margin_km(fp, (dr, cr)),
            "margin_km_4p5g": T3.margin_km(capped, (dr, cr)), "liftup_time_s": up["t_s"],
            "liftup_peak_g": up["peak_g"], "max_bank_deg_4p5g": T3._max_bank(capped)}


def _finite(x: float, fill: float = -1e9) -> float:
    return round(float(x), 3) if math.isfinite(x) else fill


def _pool_map(fn, jobs: Sequence[dict], workers: int, log=print, label: str = "jobs") -> List[dict]:
    if workers <= 1:
        return [fn(j) for j in jobs]
    out = []
    with mp.get_context("spawn").Pool(int(workers)) as pool:
        for i, r in enumerate(pool.imap(fn, jobs, chunksize=4)):
            out.append(r)
            if (i + 1) % 250 == 0:
                log(f"[benchmark] {label}: {i + 1}/{len(jobs)}")
    return out


def build_orb1(d90_km: float, d0_km: float, spec: ORB1Spec = ORB1Spec(), n: Optional[int] = None,
               workers: int = 1, log=print):
    """The ORB-1 case table (pandas) and its targets. n < n_cases builds only the first n cases."""
    import pandas as pd
    targets = orb1_targets(d90_km, d0_km, spec)
    rows = orb1_dispersions(spec, n)
    jobs = [{"row": r, "target": targets[r["target_index"]], "spec": asdict(spec)} for r in rows]
    margins = {m["case"]: m for m in _pool_map(_orb1_margin_job, jobs, workers, log, "ORB-1 footprints")}
    for r in rows:
        t = targets[r["target_index"]]
        m = margins[r["case"]]
        r.update({"target_frac": t["frac"], "target_downrange_km": t["downrange_km"], "config_path": t["config_path"],
                  "target_phi_deg": t["phi_deg"], "target_lambda_deg": t["lambda_deg"],
                  "downrange_km": m["downrange_km"], "crossrange_km": m["crossrange_km"],
                  "margin_km": _finite(m["margin_km"]), "feasible": bool(m["margin_km"] > 0.0),
                  "margin_km_4p5g": _finite(m["margin_km_4p5g"]), "feasible_4p5g": bool(m["margin_km_4p5g"] > 0.0),
                  "in_training_domain": in_training_domain(r), "liftup_time_s": m["liftup_time_s"],
                  "liftup_peak_g": m["liftup_peak_g"], "max_bank_deg_4p5g": m["max_bank_deg_4p5g"]})
    cols = ["case", "block", "target_index", "target_frac", "target_downrange_km", "target_phi_deg",
            "target_lambda_deg", "config_path", "train_seed", *T3.IC_FIELDS, "rho_mult", "cd_mult", "cl_mult",
            "mass_mult", "density_scale", "lift_scale", "mass_scale", "in_training_domain", "downrange_km",
            "crossrange_km", "margin_km", "feasible", "margin_km_4p5g", "feasible_4p5g", "liftup_time_s",
            "liftup_peak_g", "max_bank_deg_4p5g"]
    return pd.DataFrame(rows)[cols], targets


# ============================================================================
# DEV48: new testset_v3 cases, A256's Sobol indices skipped
# ============================================================================
def build_dev48(workers: int = 1, n_new: int = 12, sobol_seed: int = 20261104, log=print):
    """testset_v3's candidate filter with a new Sobol seed; any Sobol index A256 used is skipped."""
    import pandas as pd
    import torch
    spec = T3.TestsetSpec(n_bridge=0, n_new=int(n_new), sobol_seed=int(sobol_seed))
    a256 = pd.read_csv(A256_CSV)
    skip = {int(i) for i in a256["sobol_index"] if int(i) >= 0}
    pts = torch.quasirandom.SobolEngine(dimension=6, scramble=True, seed=spec.sobol_seed).draw(
        spec.max_candidates).double().numpy()
    order = [j for j in range(len(pts)) if j not in skip]
    rows, rejections, k = [], {}, 0
    pool = mp.get_context("spawn").Pool(int(workers)) if workers > 1 else None
    try:
        while len(rows) < spec.n_new and k < len(order):
            chunk = order[k:k + max(8, int(workers) * 2)]
            k += len(chunk)
            jobs = [{"u": pts[j].tolist(), "spec": asdict(spec)} for j in chunk]
            res = pool.map(T3._candidate_job, jobs) if pool else [T3._candidate_job(j) for j in jobs]
            for j, (cand, why) in zip(chunk, res):
                if len(rows) >= spec.n_new:
                    break
                rejections[why] = rejections.get(why, 0) + 1
                if cand is not None:
                    rows.append({"case": len(rows), "set": "dev", "sobol_index": j, **cand["row"]})
            log(f"[benchmark] DEV48: {len(rows)}/{spec.n_new} cases from {k} Sobol candidates")
    finally:
        if pool:
            pool.close()
            pool.join()
    if len(rows) < spec.n_new:
        raise RuntimeError(f"DEV48: only {len(rows)} feasible cases")
    rejections.pop("accepted", None)
    return pd.DataFrame(rows), {"spec": asdict(spec), "skipped_sobol_indices": sorted(skip),
                                "rejections": rejections, "n_candidates": k}


# ============================================================================
# Writing, checksums, loading
# ============================================================================
def _csv_bytes(df) -> bytes:
    buf = io.StringIO()
    df.to_csv(buf, index=False, float_format=FLOAT_FMT, lineterminator="\n")
    return buf.getvalue().encode()


def _json_bytes(obj) -> bytes:
    return (json.dumps(obj, indent=1, sort_keys=True) + "\n").encode()


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def a256_meta() -> dict:
    return {"suite": "A256", "file": "configs/v3_testset.csv", "sha256": sha256_file(A256_CSV), "n_cases": 64,
            "atmospheres": T3.ATMOSPHERES, "n_conditions": 256,
            "note": "the v3 test set unchanged; every v3 and v4 flight pairs with v5 flights on it"}


def write_files(files: Dict[str, bytes], out_dir: Path = BENCH_DIR) -> None:
    """Write the suite files and SHA256SUMS (every file except SHA256SUMS itself, sorted)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, data in files.items():
        tmp = out_dir / (name + ".tmp")
        tmp.write_bytes(data)
        os.replace(tmp, out_dir / name)
    names = sorted(p.name for p in out_dir.iterdir() if p.is_file() and p.name != "SHA256SUMS"
                   and not p.name.endswith(".tmp"))
    sums = "".join(f"{sha256_file(out_dir / n)}  {n}\n" for n in names)
    (out_dir / "SHA256SUMS").write_bytes(sums.encode())


def verify(bench_dir: Path = BENCH_DIR) -> Dict[str, str]:
    """Check every file listed in SHA256SUMS (and A256's pointer); raise ChecksumError on any mismatch."""
    bench_dir = Path(bench_dir)
    sums_path = bench_dir / "SHA256SUMS"
    if not sums_path.exists():
        raise ChecksumError(f"{sums_path} missing")
    want = {}
    for line in sums_path.read_text().splitlines():
        if line.strip():
            h, name = line.split(None, 1)
            want[name.strip()] = h
    bad = [n for n, h in want.items() if not (bench_dir / n).exists() or sha256_file(bench_dir / n) != h]
    meta = json.loads((bench_dir / "A256.json").read_text())
    if sha256_file(A256_CSV) != meta["sha256"]:
        bad.append("configs/v3_testset.csv (A256.json pointer)")
    if bad:
        raise ChecksumError(f"benchmark files changed since they were frozen: {bad}")
    return want


def load_suite(name: str, bench_dir: Path = BENCH_DIR, check: bool = True) -> List[dict]:
    """
    The conditions of one suite as dicts, checksums verified first. Every condition carries case, suite,
    atmosphere (ORB-1: "orb1"), ic, domain, config_path and the suite's own columns.
    """
    import pandas as pd
    bench_dir = Path(bench_dir)
    if check:
        verify(bench_dir)
    if name == "A256" or name == "DEV48":
        path = A256_CSV if name == "A256" else bench_dir / "DEV48.csv"
        out = []
        for c in pd.read_csv(path).to_dict("records"):
            for atm, dom in T3.ATMOSPHERES.items():
                out.append({**c, "suite": name, "atmosphere": atm, "ic": T3.ic_of(c), "domain": dict(dom),
                            "config_path": "configs/default.json"})
        return out
    if name == "ORB1":
        out = []
        for c in pd.read_csv(bench_dir / "ORB1.csv").to_dict("records"):
            out.append({**c, "suite": "ORB1", "atmosphere": "orb1", "ic": T3.ic_of(c), "domain": domain_of(c)})
        return out
    raise ValueError(f"unknown suite {name!r}")


def condition_key(cond: dict) -> str:
    """File stem of one condition: <atm>_case<NN> (A256, DEV48) or case<NNNN> (ORB-1)."""
    if cond["suite"] == "ORB1":
        return f"case{int(cond['case']):04d}"
    return f"{cond['atmosphere']}_case{int(cond['case']):02d}"


README = """# v5 benchmark (`configs/benchmark_v5/`)

Frozen on {date} by `jepa/benchmark_v5.py --build` (git commit recorded in `build_info.json`).
Every file here is listed in `SHA256SUMS` and verified whenever a suite loads
(`benchmark_v5.verify`); a mismatch raises. Nothing here is a simulator setting: the per-target
mission configs are copies of `configs/default.json` with only the target changed.

## Purpose

A fixed, documented set of test conditions for comparing entry-guidance controllers on the
ReEntryAI 3DOF simulator (Orion-like capsule, non-rotating spherical Earth, 1-D tabulated
atmosphere, bank actuator 20 deg/s, drogue deploy at 7.6 km). Plan and preregistration:
`jepa/v5_plan.md` section 3.

## Suites

| Suite | Conditions | Use |
|---|---|---|
| A256 | `configs/v3_testset.csv` (64 cases) x 4 atmospheres (`A256.json`) | results; continuity with v3 and v4 |
| ORB-1 | `ORB1.csv`, 3,000 dispersed cases, 3 blocks of 1,000 | results; the standard scenario |
| DEV48 | `DEV48.csv`, 12 cases x 4 atmospheres | tuning only, never a result |

### A256
The v3 test set unchanged (sha256 `{a256_sha}`), each case in four atmospheres:
nominal; thin (density x0.88, L/D x0.88); dense (x1.12, x1.12); out of range (x0.78, L/D x0.80).

### ORB-1 (Orion-class LEO return)
Nominal entry interface: 121.92 km (400,000 ft), 7,623.5 m/s Earth-relative (Lu 2014 Table 1,
CEV), flight-path angle -1.40 deg, at (0 N, 0 E) heading chi = 45.99 deg (east toward north).
Final condition: drogue deploy at 7.6 km (Lu's final altitude).

Targets on the entry great circle at 25 %, 50 % and 75 % of the nominal footprint between a
90 deg and a 0 deg constant bank, both flown in the full simulator from the nominal entry
(90 deg: {d90:.1f} km, {t90:.0f} s, peak {g90:.1f} g, ends with |cos gamma| < 0.05 at {end90:.1f} km
altitude, falling almost vertically; 0 deg: {d0:.1f} km, {t0:.0f} s, peak {g0:.1f} g):

| Target | Fraction | Downrange (km) | Latitude (deg) | Longitude (deg) | Mission config |
|---|---|---|---|---|---|
{target_rows}

Dispersions: Gaussian, truncated at +-3 sigma, independent, drawn case by case with
`numpy.random.default_rng({seed})` in the order {draws} (rejection per variable):

| Parameter | Nominal | 3 sigma | Source |
|---|---|---|---|
| EI flight-path angle | -1.40 deg | 0.10 deg | Bairstow 2006 Table C.1 |
| EI azimuth | 45.99 deg | 0.10 deg | Bairstow |
| EI speed | 7,623.5 m/s | 15 m/s | decision D4 (no published LEO value) |
| vehicle mass | 8,500 kg | 5 % | Bairstow |
| density multiplier | 1 | 30 % | Bairstow ("day of flight") |
| CD multiplier | 1 | 5 % | Bairstow (hypersonic) |
| CL multiplier | 1 | 5 % | Bairstow (hypersonic) |

Simulator knobs (`data.apply_domain`): `density_scale = rho_mult x cd_mult`,
`lift_scale = cl_mult / cd_mult`, `mass_scale = mass_mult`, so drag scales as rho S CD c_D and
lift as rho S CL c_L exactly. Case i is in block i // 1000, flies target (i mod 1000) mod 3, and
learned arms fly it with training seed i mod 5.

Columns of `ORB1.csv`: `case, block, target_index, target_frac, target_downrange_km,
target_phi_deg, target_lambda_deg, config_path, train_seed`, the entry state (`altitude_m,
velocity_mps, gamma_deg, chi_deg, phi_deg, lambda_deg`), the multipliers and knobs,
`in_training_domain` (density and lift scale in 0.70-1.30, mass scale in 0.90-1.10: the world
models' training randomization), the target in the case's own entry frame (`downrange_km,
crossrange_km`), `margin_km` (signed distance of the target inside the 0-90 deg constant-bank
footprint under the case's true parameters, `testset_v3.footprint` with no load cap, the footprint
the targets are defined on; > 0 is `feasible`), `margin_km_4p5g` and `feasible_4p5g` (the same
with only the banks whose peak stays within 4.5 g, testset_v3's A256 rule: a conservative bound,
since a load-relieving guidance varies the bank), `max_bank_deg_4p5g`, and the lift-up flight's
time and peak g. {n_feasible} of {n_cases} cases are feasible ({feas_pct:.2f} %); {n_feasible45}
are feasible within 4.5 g at a constant bank ({feas45_by_target}); {n_in_domain} are inside the
training domain.

Sequential stopping: every arm flies block 1; after each block an arm with >= 2 flights that are
not (drogue and miss <= 10 km) has failed the requirement and stops.

### DEV48
`testset_v3` candidate filter (same `TestsetSpec`, Sobol seed 20261104, A256's Sobol indices
skipped), 12 cases, same columns as `configs/v3_testset.csv`, in the four A256 atmospheres.
Every tunable baseline is tuned here and only here.

## Scoring (`benchmark_v5.score`)
Per arm x suite x (atmosphere or target range): success within 2/5/10/25 km (outcome `drogue`),
mean log miss, the Orion EM-1 requirement test (P(miss <= 10 km) >= 0.99865 with 90 %
confidence: one-sided Clopper-Pearson lower bound, at most one miss in 3,000), NASA-STD-3000
g-duration compliance for a deconditioned crew (Bairstow Table D.3, App. D.2 procedure: 0.2 g
levels from the peak down to 1 g, time above each level against the log-log limit curve),
authority saturation, planning time per decision.

## Sources
Lu, P., "Entry Guidance: A Unified Method," JGCD 37(3), 2014. Bairstow, S. H., "Reentry
Guidance with Extended Range Capability for Low L/D Spacecraft," MIT SM thesis, 2006.
McNamara, L. W., and Braun, R. D., AAS 14-094 (Orion EM-1 entry requirements).

## How to cite
"ReEntryAI v5 benchmark (A256, ORB-1, DEV48), configs/benchmark_v5 at commit {commit}."
"""


def build_all(workers: int = 1, out_dir: Path = BENCH_DIR, log=print, n_orb1: Optional[int] = None) -> dict:
    """Generate every benchmark file (S0's full-simulator footprint first) and write them with SHA256SUMS."""
    import datetime
    import subprocess
    spec = ORB1Spec()
    ic = spec.nominal_ic()
    log("[benchmark] S0 footprint: 90 deg and 0 deg constant bank in the full simulator")
    f90, f0 = fly_constant_bank_sim(ic, 90.0), fly_constant_bank_sim(ic, 0.0)
    for f in (f90, f0):
        # A 90 deg bank from -1.40 deg falls almost vertically at the end (|cos gamma| < 0.05 under 1 km above
        # the drogue altitude); its landing point is where it stopped (the rest is < 0.1 km of downrange).
        if not (f["outcome"] == "drogue" or (f["outcome"] == "cos_gamma_fail" and f["end_altitude_m"] < 20_000.0)):
            raise RuntimeError(f"S0 footprint flight at {f['sigma_deg']} deg ended {f['outcome']} "
                               f"at {f['end_altitude_m']:.0f} m")
    orb1, targets = build_orb1(f90["downrange_km"], f0["downrange_km"], spec, n_orb1, workers, log)
    dev48, dev_meta = build_dev48(workers, log=log)
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=SIM_DIR, capture_output=True,
                                text=True).stdout.strip() or "unknown"
    except OSError:
        commit = "unknown"
    n_feas, n_dom = int(orb1["feasible"].sum()), int(orb1["in_training_domain"].sum())
    n_feas45 = int(orb1["feasible_4p5g"].sum())
    feas45_by_target = {int(k): int(v) for k, v in orb1.groupby("target_index")["feasible_4p5g"].sum().items()}
    footprint = {"sim": "rl_env policy mode, constant command", "ic": ic,
                 "sigma90": {k: v for k, v in f90.items() if k != "y"},
                 "sigma0": {k: v for k, v in f0.items() if k != "y"}}
    files = {
        "A256.json": _json_bytes(a256_meta()),
        "ORB1.csv": _csv_bytes(orb1),
        "ORB1_spec.json": _json_bytes({"spec": asdict(spec), "draw_order": list(ORB1_DRAWS), "targets": targets,
                                       "footprint": footprint}),
        "DEV48.csv": _csv_bytes(dev48),
        "DEV48_spec.json": _json_bytes(dev_meta),
    }
    for t in targets:
        files[f"orb1_target{t['target_index']}.json"] = _json_bytes(target_config(t))
    rows = "\n".join(f"| {t['target_index']} | {t['frac']:.0%} | {t['downrange_km']:.1f} | {t['phi_deg']:.4f} | "
                     f"{t['lambda_deg']:.4f} | `{Path(t['config_path']).name}` |" for t in targets)
    files["README.md"] = README.format(
        date=datetime.date.today().isoformat(), a256_sha=a256_meta()["sha256"][:16], d90=f90["downrange_km"],
        t90=f90["t_s"], d0=f0["downrange_km"], t0=f0["t_s"], target_rows=rows, seed=spec.seed,
        draws=", ".join(ORB1_DRAWS), n_feasible=n_feas, n_cases=len(orb1), feas_pct=100.0 * n_feas / len(orb1),
        n_feasible45=n_feas45, feas45_by_target=", ".join(f"target {k}: {v}" for k, v in feas45_by_target.items()),
        n_in_domain=n_dom, end90=f90["end_altitude_m"] / 1000.0, g90=f90["peak_g"], g0=f0["peak_g"],
        commit=commit[:12]).encode()
    files["build_info.json"] = _json_bytes({"git_commit": commit, "n_orb1": int(len(orb1)), "n_dev48": int(len(dev48)),
                                            "orb1_feasible": n_feas, "orb1_feasible_4p5g": n_feas45,
                                            "orb1_feasible_4p5g_by_target": feas45_by_target,
                                            "orb1_in_training_domain": n_dom})
    write_files(files, Path(out_dir))
    return {"targets": targets, "footprint": footprint, "orb1_feasible_frac": n_feas / len(orb1)}


# ============================================================================
# Scoring
# ============================================================================
def clopper_pearson(k: int, n: int, conf: float = 0.95) -> Tuple[float, float]:
    """Two-sided exact interval for a binomial proportion."""
    from scipy.stats import beta
    if n == 0:
        return float("nan"), float("nan")
    a = 1.0 - conf
    lo = 0.0 if k == 0 else float(beta.ppf(a / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - a / 2, k + 1, n - k))
    return lo, hi


def requirement_lcb(k: int, n: int, conf: float = EM1_CONFIDENCE) -> float:
    """One-sided Clopper-Pearson lower confidence bound on P(success) at level conf."""
    from scipy.stats import beta
    if n == 0:
        return float("nan")
    if k == n:
        return float((1.0 - conf) ** (1.0 / n))
    if k == 0:
        return 0.0
    return float(beta.ppf(1.0 - conf, k, n - k + 1))


def g_limit_at(duration_s: float, table=NASA_STD_3000_DECONDITIONED) -> float:
    """The limit curve g(duration), log-log interpolation; flat beyond its ends."""
    d = np.array([t for t, _ in table], dtype=float)
    g = np.array([v for _, v in table], dtype=float)
    if duration_s <= d[0]:
        return float(g[0])
    if duration_s >= d[-1]:
        return float(g[-1])
    return float(np.exp(np.interp(math.log(duration_s), np.log(d), np.log(g))))


def g_duration(g_hist: Sequence[float], dt_s: float, step_g: float = 0.2, floor_g: float = 1.0,
               table=NASA_STD_3000_DECONDITIONED) -> dict:
    """
    Bairstow App. D.2: from the peak g rounded down to a 0.2 g level, down to 1 g, the time spent above
    each level (samples above it x dt). Compliant if every (time, level) point is under the limit curve.
    Returns compliant, the smallest margin (limit g minus level, g) and the curve points.
    """
    g = np.asarray(g_hist, dtype=float)
    g = g[np.isfinite(g)]
    if not len(g):
        return {"compliant": True, "margin_g": float("nan"), "points": []}
    top = math.floor(float(g.max()) / step_g + 1e-9) * step_g
    pts, margin = [], float("inf")
    lev = top
    while lev >= floor_g - 1e-9:
        dur = float((g > lev).sum()) * dt_s
        if dur > 0.0:
            m = g_limit_at(dur, table) - lev
            margin = min(margin, m)
            pts.append((round(lev, 3), dur))
        lev = round(lev - step_g, 6)
    return {"compliant": bool(margin >= 0.0), "margin_g": float(margin) if pts else float("nan"), "points": pts}


SCORE_KEYS = ("suite", "arm", "group")


def score(flights, by_seed: bool = False):
    """
    Section-5 metrics per arm x suite x group (atmosphere for A256/DEV48, target range for ORB-1), plus an
    "all" group per arm x suite. flights: DataFrame with v4's result columns plus suite, arm,
    g_duration_compliant, sat_frac, plan_ms_median/p99 (fly_v5 writes them). Returns a DataFrame.
    """
    import pandas as pd
    df = flights.copy()
    if not len(df):
        return pd.DataFrame()
    df["group"] = np.where(df["suite"] == "ORB1", "target" + df.get("target_index", pd.Series(0, index=df.index))
                           .fillna(-1).astype(int).astype(str), df["atmosphere"].astype(str))
    landed = df["outcome"] == "drogue"
    for r in (2, 5, 10, 25):
        df[f"s{r}"] = landed & (df["miss_km"] <= r if r == 10 else df["miss_km"] < r)
    df["log_miss"] = np.log(np.maximum(df["miss_km"].astype(float), 0.01))
    keys = ["suite", "arm"] + (["seed"] if by_seed else [])
    rows = []
    for gk, g in list(df.groupby(keys + ["group"], dropna=False)) + [
            ((*k, "all") if isinstance(k, tuple) else (k, "all"), g) for k, g in df.groupby(keys, dropna=False)]:
        n = len(g)
        k10 = int(g["s10"].sum())
        # v3/v4 flights have no g-duration column: they count as unknown, not as compliant
        comp = (g["g_duration_compliant"].dropna().map(lambda v: str(v).lower() in ("true", "1", "1.0"))
                if "g_duration_compliant" in g else pd.Series([], dtype=bool))
        lo, hi = clopper_pearson(int(comp.sum()), int(len(comp)))
        row = dict(zip(keys + ["group"], gk))
        row.update({
            "n": n, "success_2km": g["s2"].mean(), "success_5km": g["s5"].mean(), "success_10km": g["s10"].mean(),
            "success_25km": g["s25"].mean(), "mean_log_miss": g["log_miss"].mean(),
            "median_miss_km": g["miss_km"].median(), "p95_miss_km": g["miss_km"].quantile(0.95),
            "landed_frac": landed.loc[g.index].mean(), "peak_g_mean": g["peak_g"].mean(), "peak_g_max": g["peak_g"].max(),
            "heat_load_MJm2_mean": g.get("heat_load_MJm2", pd.Series(np.nan)).mean(),
            "peak_qdot_MWm2_mean": g.get("peak_qdot_MWm2", pd.Series(np.nan)).mean(),
            "rcs_fuel_kg_mean": g.get("rcs_fuel_kg", pd.Series(np.nan)).mean(),
            "n_reversals_mean": g.get("n_reversals", pd.Series(np.nan)).mean(),
            "plan_ms_median": g.get("plan_ms_median", g.get("plan_ms_mean", pd.Series(np.nan))).median(),
            "plan_ms_p99": g.get("plan_ms_p99", g.get("plan_ms_max", pd.Series(np.nan))).quantile(0.99),
            "sat_frac_mean": g.get("sat_frac", pd.Series(np.nan)).mean(),
            "n_miss_10km": n - k10, "requirement_lcb": requirement_lcb(k10, n),
            "requirement_met": bool(requirement_lcb(k10, n) >= EM1_REQUIREMENT),
            "g_compliant_frac": comp.mean() if len(comp) else np.nan, "g_compliant_lo": lo, "g_compliant_hi": hi,
        })
        rows.append(row)
    return pd.DataFrame(rows).sort_values(keys + ["group"]).reset_index(drop=True)


def orb1_eligible(flights, arm: str, block: int) -> bool:
    """Sequential stopping: an arm flies block b only if blocks < b gave fewer than 2 non-successes at 10 km."""
    if block == 0:
        return True
    d = flights[(flights["arm"] == arm) & (flights["suite"] == "ORB1") & (flights["block"] < block)]
    bad = int((~((d["outcome"] == "drogue") & (d["miss_km"] <= 10.0))).sum())
    return bad < 2


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--out", default=str(BENCH_DIR))
    a = ap.parse_args()
    if a.build:
        r = build_all(a.workers, Path(a.out))
        print(json.dumps({k: v for k, v in r.items() if k != "footprint"}, indent=1))
    if a.verify or not a.build:
        print(f"verified {len(verify(Path(a.out)))} files")
