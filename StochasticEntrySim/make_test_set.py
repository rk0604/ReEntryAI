"""
Phase 2 -- frozen held-out test set generator.

Samples N dispersed initial conditions from a FIXED seed and writes them to a
CSV. This set is the apples-to-apples benchmark: both the RL policy and the
classical predictor-corrector are evaluated on these *exact* entry conditions
(replayed via `OrionEntryEnv.reset(options={"ic_override": row})`), and the set
is NEVER used during training. Freezing the IC values (not just the seed) makes
the comparison reproducible regardless of NumPy/RNG version.

Usage:
    python make_test_set.py                 # 1000 cases -> configs/heldout_testset.csv
    python make_test_set.py 500 --seed 7    # 500 cases, seed 7
    python make_test_set.py --config configs/scenario_a_lunar.json

Load it (in eval code):
    from make_test_set import load_test_set
    for row in load_test_set("configs/heldout_testset.csv"):
        env.reset(options={"ic_override": row})
        ...
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List

import numpy as np

import mission_config

THIS_DIR = Path(__file__).resolve().parent
DEFAULT_OUT = THIS_DIR / "configs" / "heldout_testset.csv"
DEFAULT_SEED = 20260604  # fixed -> the test set is the same every time it's built

# Column order in the CSV (the six dispersable IC variables, config units).
IC_COLUMNS = (
    "altitude_m", "velocity_mps", "gamma_deg", "chi_deg", "phi_deg", "lambda_deg",
)


def generate_test_set(config_path: str = "configs/default.json",
                      n: int = 1000, seed: int = DEFAULT_SEED) -> List[Dict[str, float]]:
    """Sample `n` dispersed ICs from `seed` using the config's dispersion spec."""
    cfg = mission_config.load_config(config_path)
    spec = mission_config.dispersion_spec(cfg)
    if not spec:
        raise ValueError(
            f"{config_path} has no `dispersion` block to sample from. "
            "Add per-variable ranges before building a test set.")
    rng = np.random.default_rng(seed)
    return [mission_config.sample_dispersed_ic(cfg, rng) for _ in range(int(n))]


def write_test_set(rows: List[Dict[str, float]], out_path: Path,
                   config_path: str, seed: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    provenance = {
        "config": str(config_path), "seed": int(seed), "n": len(rows),
        "columns": list(IC_COLUMNS),
    }
    with out_path.open("w", newline="") as f:
        # Provenance as a leading '#' comment so load_test_set can recover it.
        f.write("# heldout_testset " + json.dumps(provenance) + "\n")
        w = csv.writer(f)
        w.writerow(("idx",) + IC_COLUMNS)
        for i, r in enumerate(rows):
            w.writerow((i,) + tuple(f"{float(r[c]):.6f}" for c in IC_COLUMNS))


def load_test_set(path: str | Path) -> List[Dict[str, float]]:
    """Read a frozen test set CSV -> list of IC dicts (config units).

    '#'-prefixed provenance lines are skipped. Each returned dict has the six
    IC_COLUMNS keys and is ready for `reset(options={"ic_override": row})`.
    """
    p = Path(path)
    if not p.is_absolute() and not p.exists():
        p = THIS_DIR / p
    rows: List[Dict[str, float]] = []
    with p.open("r", newline="") as f:
        lines = [ln for ln in f if not ln.lstrip().startswith("#")]
    reader = csv.DictReader(lines)
    for rec in reader:
        rows.append({c: float(rec[c]) for c in IC_COLUMNS})
    return rows


def _summarize(rows: List[Dict[str, float]]) -> str:
    out = []
    for c in IC_COLUMNS:
        vals = np.array([r[c] for r in rows])
        out.append(f"    {c:12s} [{vals.min():.3f}, {vals.max():.3f}]  mean={vals.mean():.3f}")
    return "\n".join(out)


def main() -> None:
    ap = argparse.ArgumentParser(description="Build a frozen held-out IC test set.")
    ap.add_argument("n", nargs="?", type=int, default=1000, help="number of cases")
    ap.add_argument("--config", default="configs/default.json")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()

    rows = generate_test_set(args.config, args.n, args.seed)
    out_path = Path(args.out)
    write_test_set(rows, out_path, args.config, args.seed)

    # Reload and verify round-trip reproducibility.
    reloaded = load_test_set(out_path)
    ok = (len(reloaded) == len(rows) and all(
        abs(reloaded[i][c] - rows[i][c]) < 1e-6
        for i in range(len(rows)) for c in IC_COLUMNS))
    print(f"[make_test_set] wrote {len(rows)} cases from seed {args.seed}")
    print(f"               config = {args.config}")
    print(f"               out    = {out_path}")
    print(f"               reload round-trip ok = {ok}")
    print("  IC ranges:")
    print(_summarize(rows))


if __name__ == "__main__":
    main()
