"""
Phase 2 frozen held out test set generator.

This samples N dispersed initial conditions from a fixed seed and writes them to
a CSV. The set is the like for like benchmark: both the RL policy and the
classical predictor corrector are evaluated on these exact entry conditions,
replayed with OrionEntryEnv.reset(options={"ic_override": row}), and the set is
never used during training. Freezing the initial condition values, not just the
seed, makes the comparison reproducible regardless of the NumPy or RNG version.

Usage:
    python make_test_set.py                 1000 cases to configs/heldout_testset.csv
    python make_test_set.py 500 --seed 7    500 cases, seed 7
    python make_test_set.py --config configs/scenario_a_lunar.json

Loading it in evaluation code:
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
DEFAULT_SEED = 20260604  # fixed so the test set is the same every time it is built

# Column order in the CSV. These are the six dispersable initial condition
# variables, in config units.
IC_COLUMNS = (
    "altitude_m", "velocity_mps", "gamma_deg", "chi_deg", "phi_deg", "lambda_deg",
)


def generate_test_set(config_path: str = "configs/default.json",
                      n: int = 1000, seed: int = DEFAULT_SEED) -> List[Dict[str, float]]:
    """
    Sample dispersed initial conditions using the config dispersion spec.

    Parameters:
        config_path: mission config that defines the dispersion ranges.
        n:           number of cases to sample.
        seed:        RNG seed. The same seed always produces the same set.

    Returns:
        A list of initial condition dictionaries, each in config units.

    Raises:
        ValueError when the config has no dispersion block to sample from.
    """
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
    """
    Write the sampled cases to a CSV with a provenance header.

    Parameters:
        rows:        list of initial condition dictionaries to write.
        out_path:    destination CSV path. Parent folders are created.
        config_path: config that produced the cases, recorded for provenance.
        seed:        seed that produced the cases, recorded for provenance.

    The first line is a comment that starts with # and holds a small JSON
    record of the config, seed, count, and columns, so the source of the file
    can be recovered later. Each value is written with six decimals.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    provenance = {
        "config": str(config_path), "seed": int(seed), "n": len(rows),
        "columns": list(IC_COLUMNS),
    }
    with out_path.open("w", newline="") as f:
        # Provenance is stored as a leading comment so load_test_set can skip it.
        f.write("# heldout_testset " + json.dumps(provenance) + "\n")
        w = csv.writer(f)
        w.writerow(("idx",) + IC_COLUMNS)
        for i, r in enumerate(rows):
            w.writerow((i,) + tuple(f"{float(r[c]):.6f}" for c in IC_COLUMNS))


def load_test_set(path: str | Path) -> List[Dict[str, float]]:
    """
    Read a frozen test set CSV into a list of initial condition dictionaries.

    Parameters:
        path: CSV path. A relative path that does not exist is looked up next
              to this file.

    Returns:
        A list of dictionaries, each with the six IC_COLUMNS keys in config
        units, ready for reset(options={"ic_override": row}). Lines that start
        with # are treated as provenance and skipped.
    """
    p = Path(path)
    if not p.is_absolute() and not p.exists():
        p = THIS_DIR / p
    rows: List[Dict[str, float]] = []
    with p.open("r", newline="") as f:
        # Drop the provenance comment lines before parsing.
        lines = [ln for ln in f if not ln.lstrip().startswith("#")]
    reader = csv.DictReader(lines)
    for rec in reader:
        rows.append({c: float(rec[c]) for c in IC_COLUMNS})
    return rows


def _summarize(rows: List[Dict[str, float]]) -> str:
    """
    Build a short text summary of the per column ranges and means.

    Parameters:
        rows: list of initial condition dictionaries.

    Returns:
        A multi line string with the min, max, and mean of each column.
    """
    out = []
    for c in IC_COLUMNS:
        vals = np.array([r[c] for r in rows])
        out.append(f"    {c:12s} [{vals.min():.3f}, {vals.max():.3f}]  mean={vals.mean():.3f}")
    return "\n".join(out)


def main() -> None:
    """
    Command line entry point. Build a test set, write it, and verify it.

    Parses the case count, config, seed, and output path, generates the cases,
    writes the CSV, then reloads it and checks that every value matches within a
    small tolerance. Prints the counts, paths, the round trip check, and the
    per column ranges.
    """
    ap = argparse.ArgumentParser(description="Build a frozen held-out IC test set.")
    ap.add_argument("n", nargs="?", type=int, default=1000, help="number of cases")
    ap.add_argument("--config", default="configs/default.json")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()

    rows = generate_test_set(args.config, args.n, args.seed)
    out_path = Path(args.out)
    write_test_set(rows, out_path, args.config, args.seed)

    # Reload and confirm the values survive the write and read unchanged.
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
