"""
Phase 3 -- baseline characterization.

Run the classical predictor-corrector (OrionEntryEnv in "baseline" mode) over
the frozen held-out test set and record per-case metrics. This produces the
metric *distributions* the RL policy must beat (Phase 5). The baseline is
deterministic for a given IC (the controlled segment ends at drogue deploy,
before any stochastic CPAS effects), so results are reproducible.

The predictor is the slow path (~minutes/run: it does many long nominal
rollouts per guidance cycle), so cases are run in parallel across CPU cores.
Each worker process builds one env (and one atmosphere table) and reuses it for
all the cases it's handed.

Usage:
    python run_baseline.py                       # full test set, all cores
    python run_baseline.py --limit 24            # quick subset first
    python run_baseline.py --workers 8 --out runs/baseline_metrics.csv
    python run_baseline.py --hifi-atmosphere     # disable the fast table

Output: a CSV of per-case metrics (one row per IC) + a printed distribution
summary. Read it later with pandas for the Pareto comparison.
"""

from __future__ import annotations

import argparse
import csv
import os
import time
from pathlib import Path
from typing import Dict, List, Optional

THIS_DIR = Path(__file__).resolve().parent
DEFAULT_OUT = THIS_DIR / "runs" / "baseline_metrics.csv"

METRIC_COLUMNS = (
    "idx", "alt0_m", "V0_mps", "gamma0_deg", "chi0_deg", "phi0_deg", "lam0_deg",
    "outcome", "miss_km", "peak_g", "peak_qdot_MWm2", "rcs_fuel_kg",
    "success", "t_s", "wall_s",
)

# --- per-worker globals (built once per process) -----------------------------
_ENV = None
_CONFIG = None
_FAST_ATM = True


def _worker_init(config_path: str, fast_atmosphere: bool) -> None:
    """Build the baseline env once in each worker process."""
    global _ENV, _CONFIG, _FAST_ATM
    import rl_env  # imported in the child so the atmosphere table builds here
    _CONFIG = config_path
    _FAST_ATM = fast_atmosphere
    _ENV = rl_env.OrionEntryEnv(
        config_path=config_path, controller_mode="baseline",
        fast_atmosphere=fast_atmosphere)


def _run_one(case: Dict) -> Dict:
    """Run the predictor-corrector on one frozen IC; return its metric row."""
    global _ENV
    idx = int(case["idx"])
    ic = case["ic"]
    t0 = time.time()
    _ENV.reset(options={"ic_override": ic})
    done = False
    info = {}
    # action is ignored in baseline mode; sample() is just a placeholder.
    a = _ENV.action_space.sample()
    while not done:
        _, _, term, trunc, info = _ENV.step(a)
        done = term or trunc
    wall = time.time() - t0
    return {
        "idx": idx,
        "alt0_m": ic["altitude_m"], "V0_mps": ic["velocity_mps"],
        "gamma0_deg": ic["gamma_deg"], "chi0_deg": ic["chi_deg"],
        "phi0_deg": ic["phi_deg"], "lam0_deg": ic["lambda_deg"],
        "outcome": info.get("outcome", "?"),
        "miss_km": float(info.get("miss_km", float("nan"))),
        "peak_g": float(info.get("peak_g", float("nan"))),
        "peak_qdot_MWm2": float(info.get("peak_qdot_MWm2", float("nan"))),
        "rcs_fuel_kg": float(info.get("rcs_fuel_kg", float("nan"))),
        "success": int(bool(info.get("success", False))),
        "t_s": float(info.get("t_s", float("nan"))),
        "wall_s": round(wall, 2),
    }


def _percentile(vals: List[float], q: float) -> float:
    s = sorted(v for v in vals if v == v)  # drop NaN
    if not s:
        return float("nan")
    i = min(len(s) - 1, max(0, int(round(q / 100.0 * (len(s) - 1)))))
    return s[i]


def _summarize(rows: List[Dict]) -> str:
    n = len(rows)
    succ = sum(r["success"] for r in rows)
    lines = [f"  cases             : {n}",
             f"  success (<=target): {succ}/{n} ({100.0*succ/max(1,n):.1f}%)"]
    outc: Dict[str, int] = {}
    for r in rows:
        outc[r["outcome"]] = outc.get(r["outcome"], 0) + 1
    lines.append(f"  outcomes          : {dict(sorted(outc.items()))}")
    for col, label in (("miss_km", "miss [km]"), ("peak_g", "peak g"),
                       ("peak_qdot_MWm2", "peak q [MW/m2]"), ("rcs_fuel_kg", "RCS fuel [kg]")):
        vals = [r[col] for r in rows]
        lines.append(f"  {label:16s}: median={_percentile(vals,50):.2f}  "
                     f"p90={_percentile(vals,90):.2f}  max={_percentile(vals,100):.2f}")
    return "\n".join(lines)


def run_batch(config_path: str, test_set_path: str, out_path: Path,
              limit: Optional[int], workers: int, fast_atmosphere: bool) -> List[Dict]:
    from make_test_set import load_test_set
    ics = load_test_set(test_set_path)
    if limit is not None:
        ics = ics[:limit]
    cases = [{"idx": i, "ic": ic} for i, ic in enumerate(ics)]
    n = len(cases)
    workers = max(1, min(workers, n))
    print(f"[run_baseline] {n} cases | {workers} workers | "
          f"fast_atmosphere={fast_atmosphere}")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: List[Dict] = []
    t0 = time.time()
    if workers == 1:
        _worker_init(config_path, fast_atmosphere)
        for k, c in enumerate(cases):
            rows.append(_run_one(c))
            _progress(k + 1, n, t0)
    else:
        import multiprocessing as mp
        ctx = mp.get_context("spawn")
        with ctx.Pool(processes=workers, initializer=_worker_init,
                      initargs=(config_path, fast_atmosphere)) as pool:
            for k, row in enumerate(pool.imap_unordered(_run_one, cases, chunksize=1)):
                rows.append(row)
                _progress(k + 1, n, t0)

    rows.sort(key=lambda r: r["idx"])
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=METRIC_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    print(f"\n[run_baseline] wrote {len(rows)} rows -> {out_path}  "
          f"({time.time()-t0:.0f}s wall)")
    print(_summarize(rows))
    return rows


def _progress(done: int, total: int, t0: float) -> None:
    el = time.time() - t0
    rate = done / el if el > 0 else 0.0
    eta = (total - done) / rate if rate > 0 else float("inf")
    print(f"  [{done}/{total}] {el:.0f}s elapsed, ~{eta:.0f}s ETA", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description="Characterize the classical baseline.")
    ap.add_argument("--config", default="configs/default.json")
    ap.add_argument("--test-set", default="configs/heldout_testset.csv")
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--limit", type=int, default=None, help="run only the first N cases")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--hifi-atmosphere", action="store_true",
                    help="disable the fast atmosphere table (full PyMSIS)")
    args = ap.parse_args()
    run_batch(args.config, args.test_set, Path(args.out), args.limit,
              args.workers, not args.hifi_atmosphere)


if __name__ == "__main__":
    main()
