"""
jepa/analyze_v4.py

Everything the v4 study reports (jepa/v4_plan.md sections 5, 6, 8, 10), written
to 05_analysis/ as numbered tables (T##) and figures (F##) with captions.json.

    open_loop_seed     footprint error and along-track bias of one model and seed
                       on the held-out test and validation flights (03_open_loop/)
    load_flights_v4    v3's and v4's flights in one table (arm, load, seed, case,
                       atmosphere, metrics), v3 read-only
    Boot               the paired hierarchical bootstrap (analyze_v3.hboot's
                       resampling: seeds, then cases within each atmosphere),
                       generalized to any statistic of several arms and cached in
                       05_analysis/cache/ keyed by comparison and a hash of its
                       input values, so a rerun recomputes only what changed
    analyze            T01-T13, F10-F18 and the criteria (P1-P6, L1-L3, every A
                       and B variant, R0-check, M)

Definitions (section 5): success@5km / @2km = outcome drogue and miss < 5 / 2 km;
log miss = ln(max(miss, 0.01 km)), failures keep their miss at termination;
paired log-miss difference = mean over conditions of log miss(A) - log miss(B).
Deterministic arms (classical, classical_oracle) have one run per condition,
broadcast over seeds; physics arms' "seed s" is noise index s (v2 seed s's stream).
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
import zlib
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

import numpy as np
import pandas as pd

from jepa import analyze_v3 as A3
from jepa import study_v3 as S
from jepa import study_v4 as S4
from jepa import testset_v3 as TS

ATM = list(TS.ATMOSPHERES)
ATM_LAB = A3.ATM_LAB
DETERMINISTIC = ("classical", "classical_oracle", "classical_E0a")
G = "4.5g"
CL, OR = ("classical", "own"), ("classical_oracle", "own")
MN, MNm, MT, MP = ("mppi_nominal", G), ("mppi_nominal", "matched"), ("mppi_truth", G), ("mppi_prior", G)
V2, V2m, V2S, V1 = ("v2", G), ("v2", "matched"), ("v2single", G), ("v1", G)
LOAD_ARMS = [("v2", "matched"), ("v2", "3.5g"), ("v2", "4g"), ("v2", G), ("v2", "8g")]
E1_ARMS = [CL, OR, MN, MP, V2, MT]
LAB = {"classical": "classical", "classical_oracle": "classical, knows density & mass", "mppi_nominal": "MPPI + nominal physics",
       "mppi_truth": "MPPI + true physics", "mppi_prior": "MPPI + v2's physics prior", "v2": "v2 ensemble",
       "v2single": "v2 single member", "v1": "v1 (fair budget)", "R0": "R0 (v2 recipe, 1 member)"}
COL = {"classical": "#6250d6", "classical_oracle": "#a3a29c", "mppi_nominal": "#1baf7a", "mppi_truth": "#0b0b0b",
       "mppi_prior": "#8fb8ea", "v2": "#eb6834", "v2single": "#eda100", "v1": "#2a78d6", "R0": "#c2410c"}
NON_TABLE_ARMS = ("diag_case10", "classical_E0a", "v1_v4run", "v2_v4run", "v2single_v4run", "E0c")


def label(arm: str, load: str = G) -> str:
    base = LAB.get(arm, arm)
    return base if load in (G, "own") else f"{base} @{load}"


def color(arm: str) -> str:
    return COL.get(arm, "#52514e")


# ============================================================================
# Captions, figures, tables (guarded writes)
# ============================================================================
class Captions:
    """05_analysis/captions.json: number -> caption and data source."""

    def __init__(self, L: S4.LayoutV4):
        self.L, self.path = L, L.analysis / "captions.json"
        self.d = json.loads(self.path.read_text()) if self.path.exists() else {}

    def add(self, key: str, caption: str, data: str) -> None:
        self.d[key] = {"caption": caption, "data": data}
        self.L.write_json(self.path, self.d)


def _fig(L, caps, plt, fig, name, caption, data) -> Path:
    path = L.figures / f"{name}.png"
    L.save_fig(fig, path)
    plt.close(fig)
    caps.add(name, caption, data)
    return path


def _table(L, caps, df: pd.DataFrame, name, caption, data) -> Path:
    path = L.tables / f"{name}.csv"
    L.write_csv(df, path)
    caps.add(name, caption, data)
    return path


# ============================================================================
# Open loop (section 4.3: per variant and seed, test and validation flights)
# ============================================================================
def splits(cfg, data_path: str) -> Dict[str, list]:
    """The v3 split (seed 0, 80/10/10) of the 6,000 flights, macro-stepped: {"val", "test"}."""
    from jepa import data as D
    from jepa import train_v2 as T2
    eps = T2.load_episodes(data_path)
    if cfg.n_episodes:
        eps = eps[:cfg.n_episodes]
    _, va, te = D.split_episodes([D.to_macro(e, D.MACRO_DT_S) for e in eps], (0.8, 0.1, 0.1), seed=cfg.split_seed)
    return {"val": va, "test": te}


def _ol_models(name: str, ckpt: str, device: str) -> dict:
    from jepa import ablate_v4 as AB
    from jepa import fly_v4 as F4
    from jepa import models_v2 as V
    if name == "v1" or (name in AB.DESCRIPTIONS and AB.family(name) == "v1"):
        m = F4.load_v1_family(ckpt)
        return {name: {"wm": m["wm"], "prober": m["pi"], "stats": m["stats"], "rf": m["wm"].rf}}
    m = V.load_checkpoint_v2(ckpt, device)
    return {name: {**V.as_v1(m), "stats": m["stats"], "rf": m["wm"].rf}}


def open_loop_seed(cfg, L: S4.LayoutV4, L3: S.Layout, eps: Dict[str, list], name: str, seed: int,
                   device: str = "cpu", h_cross: float = 7600.0) -> Optional[pd.DataFrame]:
    """Footprint errors of one model (a variant, or v3's v1) and seed; cached as 03_open_loop/footprint_*.csv."""
    import torch
    from jepa import evaluate_v2 as E2
    path = L.open_loop / f"footprint_{name}_seed{seed}.csv"
    if path.exists():
        return pd.read_csv(path)
    ck = L3.ckpt("v1", seed) if name == "v1" else L.ckpt(name, seed)
    if not Path(ck).exists():
        return None
    with torch.no_grad():
        models = _ol_models(name, str(ck), device)
        rows = []
        for split, ee in eps.items():
            r = pd.DataFrame(E2.footprint_points(models, ee, cfg.ttg_s, h_cross))
            r["split"] = split
            rows.append(r)
    rows = pd.concat(rows, ignore_index=True)
    rows["seed"] = seed
    rows["beta_scale"] = rows.dom_density_scale / rows.dom_mass_scale
    rows["density_bin"] = pd.cut(rows.beta_scale, [0, 0.85, 1.15, 10], labels=["thin", "nominal", "dense"]).astype(str)
    L.write_csv(rows, path)
    along = rows.groupby(["model", "split", "ttg_s", "density_bin"]).agg(
        along_mean_km=("along_km", "mean"), along_median_km=("along_km", "median"),
        err_median_km=("err_km", "median"), n=("err_km", "size")).reset_index()
    along["seed"] = seed
    L.write_csv(along, L.open_loop / f"along_track_{name}_seed{seed}.csv")
    return rows


def open_loop_table(L: S4.LayoutV4) -> pd.DataFrame:
    """Median footprint error per model, seed, split and time-to-go (all cached open-loop files)."""
    files = sorted(L.open_loop.glob("footprint_*_seed*.csv"))
    if not files:
        return pd.DataFrame(columns=["model", "seed", "split", "ttg_s", "err_median_km", "err_p90_km", "n"])
    fp = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    return fp.groupby(["model", "seed", "split", "ttg_s"]).agg(
        err_median_km=("err_km", "median"), err_p90_km=("err_km", lambda s: s.quantile(0.9)),
        n=("err_km", "size")).reset_index()


# ============================================================================
# Flights
# ============================================================================
def _read_results(root: Path) -> pd.DataFrame:
    rows = []
    for p in sorted((root / "04_flights" / "results").rglob("*.json")):
        try:
            rows.append(json.loads(p.read_text()))
        except (OSError, json.JSONDecodeError):
            print(f"  [collect] unreadable {p}, skipped")
    return pd.DataFrame(rows)


def load_flights_v4(cfg, L: S4.LayoutV4, L3: S.Layout) -> pd.DataFrame:
    """
    v3's flights (its 04_flights/flights.csv, else its result files; read-only)
    and v4's (results rebuilt into v4's flights.csv) in one table. If E0b failed,
    v4's re-flown v1, v2 and v2single replace v3's at 4.5 g (section 4.0).
    """
    f3 = L3.flights / "flights.csv"
    d3 = pd.read_csv(f3) if f3.exists() else _read_results(L3.root)
    d3["source"] = "v3"
    d4 = S4.collect_flights(L, write=True)
    if len(d4):
        d4["source"] = "v4"
    df = pd.concat([d3, d4], ignore_index=True) if len(d4) else d3
    df["arm"] = df.controller
    if S4.e0b_failed(L):
        df = df[~((df.source == "v3") & df.arm.isin(["v1", "v2", "v2single"]) & (df.load == G))]
        ren = {"v1_v4run": "v1", "v2_v4run": "v2", "v2single_v4run": "v2single"}
        df["arm"] = df.arm.replace(ren)
    df["seed"] = pd.to_numeric(df.seed, errors="coerce").fillna(-1).astype(int)
    df["mppi_rep"] = pd.to_numeric(df.mppi_rep, errors="coerce").fillna(-1).astype(int)
    df["load"] = df.load.astype(str)
    df["case"] = df.case.astype(int)
    df["success_5km"] = df.success_5km.astype(str).str.lower().isin(["true", "1", "1.0"])
    df["success_2km"] = (df.outcome == "drogue") & (df.miss_km < cfg.success2_km)
    df["log_miss"] = np.log(df.miss_km.clip(lower=cfg.log_miss_floor_km))
    cases = S4.load_testset(cfg)
    df = df[df.case.isin(cases.case)]
    feas = cases.melt(id_vars=["case"], value_vars=[f"feasible_{t}" for t in ATM], var_name="atmosphere",
                      value_name="feasible")
    feas["atmosphere"] = feas.atmosphere.str.replace("feasible_", "", regex=False)
    return df.merge(feas, on=["case", "atmosphere"], how="left")


METRICS: Dict[str, Callable] = {
    **{k: v for k, v in A3.METRICS.items()},
    "success2_%": lambda g: 100 * g.success_2km.mean(),
    "log_miss_mean": lambda g: g.log_miss.mean(),
}


def metric_table(df: pd.DataFrame, by=("arm", "load", "atmosphere")) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """analyze_v3.metric_table with success@2km and mean log miss: per seed, then mean/min/max over seeds + pooled."""
    by = list(by)
    rows = [{**dict(zip(by + ["seed"], key)), "n": len(g), **{k: f(g) for k, f in METRICS.items()}}
            for key, g in df.groupby(by + ["seed"], dropna=False)]
    per_seed = pd.DataFrame(rows)
    out = []
    for key, g in per_seed.groupby(by, dropna=False):
        key = key if isinstance(key, tuple) else (key,)
        mask = np.ones(len(df), dtype=bool)
        for col, v in zip(by, key):
            mask &= (df[col] == v).to_numpy()
        pooled = df[mask]
        r = {**dict(zip(by, key)), "n_seeds": g.seed.nunique(), "n_flights": int(g.n.sum())}
        for k in METRICS:
            r[k], r[f"{k}_min"], r[f"{k}_max"], r[f"{k}_pooled"] = g[k].mean(), g[k].min(), g[k].max(), METRICS[k](pooled)
        out.append(r)
    return pd.DataFrame(out), per_seed


def cube(df: pd.DataFrame, arm_load: tuple, col: str, seeds: Sequence[int], cases: Sequence[int]) -> np.ndarray:
    """Values [n_seeds, n_atm, n_cases] (NaN where a flight is missing); deterministic arms [1, A, C]."""
    arm, load = arm_load
    d = df[(df.arm == arm) & (df.load == load) & (df.mppi_rep <= 0)]
    ss = [-1] if arm in DETERMINISTIC else list(seeds)
    X = np.full((len(ss), len(ATM), len(cases)), np.nan)
    si, ai, ci = {s: i for i, s in enumerate(ss)}, {a: i for i, a in enumerate(ATM)}, {c: i for i, c in enumerate(cases)}
    for r in d[["seed", "atmosphere", "case", col]].itertuples(index=False):
        if r.seed in si and r.atmosphere in ai and r.case in ci:
            X[si[r.seed], ai[r.atmosphere], ci[r.case]] = float(getattr(r, col))
    return X


def present(df: pd.DataFrame, arm_load: tuple) -> bool:
    return bool(((df.arm == arm_load[0]) & (df.load == arm_load[1])).any())


# ============================================================================
# Bootstrap (cached)
# ============================================================================
def _m(x: np.ndarray) -> np.ndarray:
    """Mean over everything but the leading resample axis, NaN-safe."""
    with np.errstate(all="ignore"):
        return np.nanmean(x.reshape(x.shape[0], -1), axis=1)


FNS: Dict[str, Callable] = {
    "diff_mean": lambda a: _m(a[0]) - _m(a[1]),
    "paired_mean": lambda a: _m(a[0] - a[1]),
    "share": lambda a: (_m(a[0]) - _m(a[1])) / (_m(a[2]) - _m(a[3])),
    "gap_closed": lambda a: (_m(a[0]) - _m(a[1])) / (_m(a[0]) - _m(a[2])),
}


class Boot:
    """
    Paired hierarchical bootstrap over [seeds|1, atmospheres, cases] cubes that
    share every index (analyze_v3.hboot's scheme and random stream: resample
    seeds, then cases within each atmosphere; deterministic arms repeated over
    seeds), for any statistic in FNS. Results are cached by comparison key and
    a hash of the inputs; `computed` / `reused` count what happened.
    """

    def __init__(self, L: S4.LayoutV4, n_boot: int):
        self.L, self.n_boot = L, int(n_boot)
        self.computed = self.reused = 0

    def __call__(self, key: str, cubes: Sequence[np.ndarray], fn: str, atm: Optional[int] = None,
                 stream: Optional[str] = None) -> dict:
        """stream: comparisons given the same stream name draw the same resamples (default: their own)."""
        stream = stream or key
        h = hashlib.sha1(json.dumps([key, fn, self.n_boot, atm, stream]).encode())
        for X in cubes:
            h.update(str(X.shape).encode())
            h.update(np.ascontiguousarray(X, dtype=np.float64).tobytes())
        h = h.hexdigest()
        safe = re.sub(r"[^A-Za-z0-9_.@-]+", "_", key) + (f"_atm{atm}" if atm is not None else "")
        path = self.L.cache / f"{safe}.json"
        if path.exists():
            try:
                d = json.loads(path.read_text())
                if d.get("hash") == h:
                    self.reused += 1
                    return d["result"]
            except (OSError, json.JSONDecodeError):
                pass
        res = self.run(cubes, fn, atm, seed=zlib.crc32(f"{stream}|{atm}".encode()) & 0x7FFFFFFF)
        self.L.write_json(path, {"key": key, "fn": fn, "atm": atm, "hash": h, "result": res})
        self.computed += 1
        return res

    def run(self, cubes: Sequence[np.ndarray], fn: str, atm: Optional[int], seed: int) -> dict:
        S_ = max(X.shape[0] for X in cubes)
        Xs = [np.repeat(X, S_, 0) if X.shape[0] == 1 else X for X in cubes]
        if atm is not None:
            Xs = [X[:, atm:atm + 1] for X in Xs]
        A, C = Xs[0].shape[1:]
        f = FNS[fn]
        with np.errstate(all="ignore"):
            point = float(f([X[None] for X in Xs])[0])
        rng = np.random.default_rng(seed)
        vals = []
        for b0 in range(0, self.n_boot, 1000):
            nb = min(1000, self.n_boot - b0)
            si = rng.integers(0, S_, (nb, S_))
            ci = rng.integers(0, C, (nb, A, C))
            ix = (si[:, :, None, None], np.arange(A)[None, None, :, None], ci[:, None, :, :])
            with np.errstate(all="ignore"):
                vals.append(f([X[ix] for X in Xs]))
        d = np.concatenate(vals)
        d = d[np.isfinite(d)]
        if not len(d):
            return {"diff": point, "ci_lo": np.nan, "ci_hi": np.nan, "ci90_lo": np.nan, "ci90_hi": np.nan,
                    "p_boot": np.nan, "n_boot_ok": 0}
        lo, hi = np.percentile(d, [2.5, 97.5])
        lo90, hi90 = np.percentile(d, [5.0, 95.0])
        p = min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean()))
        return {"diff": point, "ci_lo": float(lo), "ci_hi": float(hi), "ci90_lo": float(lo90), "ci90_hi": float(hi90),
                "p_boot": float(max(p, 1.0 / len(d))), "n_boot_ok": int(len(d))}


def per_seed_tests(XA: np.ndarray, XB: np.ndarray, seeds: Sequence[int], kind: str) -> List[dict]:
    """Each seed alone: success (exact McNemar, difference in pp) or log miss (exact sign test, mean paired difference)."""
    rows = []
    n = max(XA.shape[0], XB.shape[0])
    for k in range(n):
        a, b = XA[k if XA.shape[0] > 1 else 0].ravel(), XB[k if XB.shape[0] > 1 else 0].ravel()
        s = seeds[k] if n > 1 else -1
        if kind == "success":
            n10, n01, p = A3.mcnemar_p(a, b)
            with np.errstate(all="ignore"):
                rows.append({"seed": s, "diff": 100 * (np.nanmean(a) - np.nanmean(b)), "a_only": n10, "b_only": n01,
                             "p": p})
        else:
            w, nn, p = A3.sign_test(a, b)
            with np.errstate(all="ignore"):
                rows.append({"seed": s, "diff": float(np.nanmean(a - b)), "a_wins": w, "nonties": nn, "p": p})
    return rows


def _holds(rows: List[dict], sign: float) -> Tuple[int, int]:
    """Seeds whose effect has the given sign (the pooled effect's), and the number of seeds."""
    k = sum(1 for r in rows if np.isfinite(r["diff"]) and np.sign(r["diff"]) == sign and r["diff"] != 0)
    return k, len(rows)


# ============================================================================
# The analysis
# ============================================================================
def _plt():
    return A3._plt()


def analyze(cfg, L: S4.LayoutV4, L3: S.Layout, log=print) -> dict:
    """All tables, statistics, criteria and figures. Returns a summary dict for the report."""
    caps = Captions(L)
    boot = Boot(L, cfg.n_boot)
    df = load_flights_v4(cfg, L, L3)
    cases = sorted(int(c) for c in S4.load_testset(cfg).case)
    seeds = [s for s in cfg.seeds if ((df.arm == "v2") & (df.load == G) & (df.seed == s)).any()] or list(cfg.seeds)
    out: Dict[str, object] = {"seeds": seeds, "n_flights": int(len(df)), "n_v4_flights": int((df.source == "v4").sum()),
                              "boot": boot}
    main = df[~df.arm.isin(NON_TABLE_ARMS) & (df.mppi_rep <= 0)]

    def cb(al, col):
        return cube(df, al, col, seeds, cases)

    # ---- T01 / T02 ----
    t01, per_seed = metric_table(main)
    _table(L, caps, t01.round(4), "T01_metrics_by_arm",
           "Every metric per arm (v3 and v4) and atmosphere: mean over seeds, seed min and max, and pooled over seeds "
           "(deterministic arms: one run).", "04_flights/flights.csv + the v3 folder's flights")
    t02, _ = metric_table(main[main.feasible.astype(bool)])
    _table(L, caps, t02.round(4), "T02_metrics_feasible_only",
           "T01 on conditions whose target is inside the conservative 4.5 g footprint in that atmosphere.",
           "as T01 + configs/v3_testset.csv")
    out["T01"], out["per_seed"] = t01, per_seed

    # ---- RQ1: planner vs model ----
    def comp(name, a, b, metric):
        """One comparison, pooled and per atmosphere, with per-seed tests."""
        if not (present(df, a) and present(df, b)):
            return []
        col, fn, scale = {"success5": ("success_5km", "diff_mean", 100.0), "success2": ("success_2km", "diff_mean", 100.0),
                          "log_miss": ("log_miss", "paired_mean", 1.0)}[metric]
        XA, XB = cb(a, col), cb(b, col)
        ps = per_seed_tests(XA, XB, seeds, "success" if metric.startswith("success") else "logmiss")
        rows = []
        for i, t in [(None, "all")] + list(enumerate(ATM)):
            r = boot(f"{name}|{a}|{b}|{metric}", [XA, XB], fn, i)
            r = {k: (v * scale if k in ("diff", "ci_lo", "ci_hi", "ci90_lo", "ci90_hi") else v) for k, v in r.items()}
            row = {"id": name, "a": label(*a), "b": label(*b), "metric": metric, "atmosphere": t, **r}
            if t == "all":
                k, n = _holds(ps, np.sign(r["diff"]) if np.isfinite(r["diff"]) and r["diff"] != 0 else 1.0)
                row.update({"seeds_holding": k, "n_seeds": n, "per_seed": json.dumps([round(p["diff"], 4) for p in ps])})
            rows.append(row)
        return rows

    rq1 = [("P1", MN, CL), ("P2", V2, MN), ("P3", MT, V2), ("P4", OR, CL), ("P5", V2, MP), ("total", V2, CL),
           ("P1@matched", MNm, CL), ("P2@matched", V2m, MNm), ("total@matched", V2m, CL), ("truth-nominal", MT, MN)]
    t04 = pd.DataFrame([r for nm, a, b in rq1 for m in ("success5", "log_miss") for r in comp(nm, a, b, m)])
    _table(L, caps, t04.round(4), "T04_rq1_paired_differences",
           "RQ1, planner vs model at 4.5 g (and at the matched load): paired differences A − B with 95% and 90% "
           f"hierarchical-bootstrap intervals ({cfg.n_boot} resamples: seeds, then cases within each atmosphere), "
           "bootstrap p, and seeds where the pooled effect holds. success5 in percentage points; log_miss = mean "
           "paired difference of ln(miss).", "04_flights/flights.csv + v3 flights")
    out["T04"] = t04

    def get(t, id_, metric="success5", atm="all"):
        if not len(t):
            return None
        r = t[(t.id == id_) & (t.metric == metric) & (t.atmosphere == atm)]
        return r.iloc[0] if len(r) else None

    # ---- T06: gap decomposition (P6) ----
    rows = []
    for load, (v2_, mn_) in ((G, (V2, MN)), ("matched", (V2m, MNm))):
        if not all(present(df, x) for x in (v2_, mn_, CL)):
            continue
        for metric, col in (("success5", "success_5km"), ("log_miss", "log_miss")):
            Xv, Xm, Xc = cb(v2_, col), cb(mn_, col), cb(CL, col)
            sc = 100.0 if metric == "success5" else 1.0
            st = f"gap|{load}|{metric}"                     # one resample stream for every piece (section 8.1, P6)
            for i, t in [(None, "all")] + list(enumerate(ATM)):
                tot = boot(f"{st}|total", [Xv, Xc], "diff_mean", i, st)
                pl = boot(f"{st}|planner", [Xm, Xc], "diff_mean", i, st)
                mo = boot(f"{st}|model", [Xv, Xm], "diff_mean", i, st)
                psh = boot(f"{st}|planner_share", [Xm, Xc, Xv, Xc], "share", i, st)
                msh = boot(f"{st}|model_share", [Xv, Xm, Xv, Xc], "share", i, st)
                ceil = boot(f"{st}|ceiling", [cb(MT, col), Xv], "diff_mean", i, st) \
                    if load == G and present(df, MT) else None
                rows.append({"load": load, "metric": metric, "atmosphere": t,
                             "classical": float(np.nanmean(Xc[:, i] if i is not None else Xc)) * sc,
                             "total": tot["diff"] * sc, "total_ci_lo": tot["ci_lo"] * sc, "total_ci_hi": tot["ci_hi"] * sc,
                             "planner": pl["diff"] * sc, "planner_ci_lo": pl["ci_lo"] * sc, "planner_ci_hi": pl["ci_hi"] * sc,
                             "model": mo["diff"] * sc, "model_ci_lo": mo["ci_lo"] * sc, "model_ci_hi": mo["ci_hi"] * sc,
                             "planner_share": psh["diff"], "planner_share_ci_lo": psh["ci_lo"],
                             "planner_share_ci_hi": psh["ci_hi"], "model_share": msh["diff"],
                             "model_share_ci_lo": msh["ci_lo"], "model_share_ci_hi": msh["ci_hi"],
                             "ceiling_truth_minus_v2": ceil["diff"] * sc if ceil else np.nan})
    t06 = pd.DataFrame(rows)
    _table(L, caps, t06.round(4), "T06_gap_decomposition",
           "The classical → v2 gap split into the planner's share (MPPI + nominal physics − classical) and the model's "
           "(v2 − MPPI + nominal physics), at 4.5 g and at the matched load, on success@5km (pp) and mean log miss; "
           "shares with bootstrap intervals from the same resamples. ceiling = MPPI + true physics − v2.",
           "05_analysis/tables/T04, flights")
    out["T06"] = t06

    # ---- RQ2: load ----
    g_mean = main.groupby(["arm", "load", "atmosphere"]).peak_g.mean()
    rows = []
    for al in [CL, MN, MNm] + LOAD_ARMS:
        d = main[(main.arm == al[0]) & (main.load == al[1])]
        if not len(d):
            continue
        for t, g in [("all", d)] + [(t, d[d.atmosphere == t]) for t in ATM]:
            if not len(g):
                continue
            lim = pd.to_numeric(g.g_limit, errors="coerce")
            rows.append({"arm": label(*al), "arm_key": al[0], "load": al[1], "atmosphere": t, "n": len(g),
                         "n_seeds": g.seed.nunique(), "g_limit_mean": float(lim.mean()) if lim.notna().any() else np.nan,
                         "peak_g_mean": g.peak_g.mean(), "success5_%": 100 * g.success_5km.mean(),
                         "success2_%": 100 * g.success_2km.mean(), "log_miss_mean": g.log_miss.mean(),
                         "miss_median_km": g.miss_km.median()})
    t07 = pd.DataFrame(rows)
    _table(L, caps, t07.round(4), "T07_load_curve",
           "L3, what load buys: v2 at the matched limit (classical's own peak g per condition), 3.5, 4.0, 4.5 (v3) and "
           "8 g (v3, seed 0 only), with the classical controller and MPPI + nominal physics for reference: mean peak g, "
           "success@5km and @2km, mean log miss.", "flights")
    out["T07"] = t07

    # ---- RQ3: ablations ----
    variants = [v for v in cfg.variants if v != "R0"]
    rows, seed_rows = [], []
    from jepa import ablate_v4 as AB
    for v in variants:
        ref = ("R0", G) if v.startswith("A") else V1
        va = (v, G)
        if not (present(df, va) and present(df, ref)):
            rows.append({"variant": v, "reference": ref[0], "change": AB.DESCRIPTIONS[v], "ran": False})
            continue
        XA, XB = cb(va, "log_miss"), cb(ref, "log_miss")
        r = boot(f"abl|{v}|log_miss", [XA, XB], "paired_mean")
        ps = per_seed_tests(XA, XB, seeds, "logmiss")
        want = 1.0 if v.startswith("A") else -1.0
        k, n = _holds(ps, want)
        row = {"variant": v, "reference": ref[0], "change": AB.DESCRIPTIONS[v], "ran": True,
               "d_log_miss": r["diff"], "ci_lo": r["ci_lo"], "ci_hi": r["ci_hi"], "p_boot": r["p_boot"],
               "seeds_holding": k, "n_seeds": n, "geo_miss_ratio": math.exp(r["diff"]) if np.isfinite(r["diff"]) else np.nan}
        for metric, col in (("success5", "success_5km"), ("success2", "success_2km")):
            rr = boot(f"abl|{v}|{metric}", [cb(va, col), cb(ref, col)], "diff_mean")
            row[f"d_{metric}_pp"], row[f"d_{metric}_ci_lo"], row[f"d_{metric}_ci_hi"] = (
                100 * rr["diff"], 100 * rr["ci_lo"], 100 * rr["ci_hi"])
            row[f"{metric}_%"] = 100 * float(np.nanmean(cb(va, col)))
        for col in ("peak_g", "rcs_fuel_kg"):
            row[f"{col}_mean"] = float(np.nanmean(cb(va, col)))
            row[f"{col}_mean_ref"] = float(np.nanmean(cb(ref, col)))
        if v.startswith("B") and present(df, ("R0", G)):
            gc = boot(f"abl|{v}|gap_closed", [cb(V1, "log_miss"), XA, cb(("R0", G), "log_miss")], "gap_closed")
            row.update({"gap_closed": gc["diff"], "gap_closed_ci_lo": gc["ci_lo"], "gap_closed_ci_hi": gc["ci_hi"]})
        try:
            row["params_per_member"] = AB.param_count(cfg, v)["per_member"]
        except Exception:
            row["params_per_member"] = np.nan
        rows.append(row)
        for p in ps:
            seed_rows.append({"variant": v, "reference": ref[0], "seed": p["seed"], "d_log_miss": p["diff"],
                              "sign_test_p": p["p"], "a_wins": p["a_wins"], "nonties": p["nonties"],
                              "success5_%": 100 * float(np.nanmean(cb(va, "success_5km")[seeds.index(p["seed"])]))
                              if p["seed"] in seeds else np.nan,
                              "success2_%": 100 * float(np.nanmean(cb(va, "success_2km")[seeds.index(p["seed"])]))
                              if p["seed"] in seeds else np.nan})
    t08 = pd.DataFrame(rows)
    ran = [i for i, r in enumerate(rows) if r.get("ran") and np.isfinite(r.get("p_boot", np.nan))]
    t08["p_holm"] = np.nan
    if ran:
        adj = A3.holm([rows[i]["p_boot"] for i in ran])
        for i, pa in zip(ran, adj):
            t08.loc[i, "p_holm"] = pa
    out["T08"] = t08
    t09 = pd.DataFrame(seed_rows)
    _table(L, caps, t09.round(4), "T09_ablations_per_seed",
           "Every ablation on each seed alone: mean paired log-miss difference against its reference (R0 for A, v1 for "
           "B), exact sign test on the paired misses, success@5km and @2km.", "flights")

    # ---- open loop (T10) and the mechanism (M) ----
    ol = open_loop_table(L)
    _table(L, caps, ol.round(4), "T10_open_loop_by_variant",
           "Open loop: median and p90 landing-point error [km] per model, seed, split (held-out test / validation "
           "flights) and time-to-go, with the commands actually flown.", "03_open_loop/footprint_*_seed*.csv")
    out["T10"] = ol
    files = sorted(L.open_loop.glob("along_track_*_seed*.csv"))
    at = pd.concat([pd.read_csv(f) for f in files], ignore_index=True) if files else pd.DataFrame(
        columns=["model", "split", "ttg_s", "density_bin", "along_mean_km", "along_median_km", "err_median_km", "n", "seed"])
    _table(L, caps, at.round(4), "T10b_along_track_bias",
           "Open loop: signed along-track error (+ = predicted landing beyond the true one) by air-density bin "
           "(density/mass scale < 0.85 thin, > 1.15 dense), per model, seed, split and time-to-go.",
           "03_open_loop/along_track_*_seed*.csv")
    mech = mechanism(cfg, df, ol, seeds)
    out["M"] = mech

    # ---- criteria (T05), with T08's verdicts ----
    crit = criteria(cfg, df, main, t04, t06, t07, t08, mech, g_mean, seeds, cases, boot, cb)
    out["T05"] = crit
    t08["met"] = [next((c["met"] for c in crit.to_dict("records") if c["id"] == v), None) for v in t08.variant]
    _table(L, caps, t08.round(4), "T08_ablations",
           "RQ3: every variant against its reference (A: R0, a change removed; B: v1, a change added). Paired log-miss "
           "difference (variant − reference; > 0 = misses more), 95% interval, bootstrap and Holm p (13 variants), "
           "seeds holding, success@5km and @2km, mean peak g and RCS fuel (variant and reference), gap closed (B), "
           "parameters per member.", "flights")
    _table(L, caps, crit, "T05_criteria",
           "Pre-registered criteria (v4_plan.md section 8, fixed before any v4 flight): what was measured and whether "
           "each is met. Holm within RQ1 (P2, P4, P5), RQ2 (L2) and RQ3 (the 13 variants).", "05_analysis/tables")

    # ---- T03: E0 and physics fidelity ----
    out["T03"] = e0_table(cfg, L)
    _table(L, caps, out["T03"], "T03_reproducibility_and_fidelity",
           "E0 (section 4.0): classical determinism (E0a), the v4 runner re-flying v3's learned flights (E0b: outcome "
           "agreement, paired log-miss difference v4 − v3 = the GPU-nondeterminism floor), runner equivalence on the "
           "CPU (E0c); physics-model fidelity (section 4.1.3: the torch port per step, and open-loop landing error).",
           "01_checks/")

    # ---- T11 case 10, T12 thresholds, T13 ensemble re-score ----
    out["T11"] = case10_table(cfg, L, L3, df)
    _table(L, caps, out["T11"], "T11_case10_diagnosis",
           f"RQ4, case {cfg.diag_case} (descriptive): part 'by_arm' = success and median miss of every arm on this case; "
           "part 'diag' = per diagnostic flight (decision log): where it separates from the classical trajectory, "
           "the planner's predicted miss and prediction error over the last 20 decisions, and how often the chosen "
           "plan's predicted peak g sat at the soft limit.", "04_flights/ (diag_case10 telemetry), flights")
    out["T12"] = threshold_table(cfg, main)
    _table(L, caps, out["T12"].round(4), "T12_success_vs_threshold",
           "Success rate against the miss threshold (0.5 to 25 km; outcome drogue and miss < threshold), per arm and "
           "atmosphere, pooled over seeds.", "flights")
    out["T13"] = ensemble_rescore(cfg, df, cb, boot)
    _table(L, caps, out["T13"].round(4), "T13_ensemble_rescore",
           "v3's ensemble question re-scored with the v4 metrics: v2 ensemble − v2 single member, paired log miss and "
           "success@2km (and @5km), pooled and per atmosphere.", "v3 flights")
    out["figures"] = figures(cfg, L, L3, caps, df, main, t01, per_seed, t04, t06, t07, t08, ol, mech, seeds, cases)
    log(f"[analyze] {len(df)} flights ({out['n_v4_flights']} v4), seeds {seeds}; bootstrap: "
        f"{boot.computed} computed, {boot.reused} from cache")
    return out


# ============================================================================
# Criteria
# ============================================================================
def criteria(cfg, df, main, t04, t06, t07, t08, mech, g_mean, seeds, cases, boot, cb) -> pd.DataFrame:
    rows = []

    def row(id_, family_, q, measured, met, value=np.nan, ci=(np.nan, np.nan), p=None, holding=None, valid=None):
        rows.append({"id": id_, "family": family_, "question": q, "measured": measured, "value": value,
                     "ci_lo": ci[0], "ci_hi": ci[1], "p": p, "p_holm": None, "seeds_holding": holding,
                     "met": met, "valid": valid})

    def g(id_, metric="success5", atm="all"):
        if not len(t04):
            return None
        r = t04[(t04.id == id_) & (t04.metric == metric) & (t04.atmosphere == atm)]
        return r.iloc[0] if len(r) else None

    def per_atm(id_, metric="success5"):
        return ", ".join(f"{ATM_LAB[t]} {r['diff']:+.1f} [{r['ci_lo']:+.1f}, {r['ci_hi']:+.1f}]"
                         for t in ATM if (r := g(id_, metric, t)) is not None)

    def seeds_ok(r):
        return int(r["seeds_holding"]) >= max(1, int(r["n_seeds"]) - 1)

    # P1
    r = g("P1")
    row("P1", "RQ1", "How much does the planner add? (MPPI + nominal physics − classical)",
        "not run" if r is None else f"success@5km {r['diff']:+.1f} pp [{r['ci_lo']:+.1f}, {r['ci_hi']:+.1f}] pooled; "
        + per_atm("P1"), None, *(() if r is None else (r["diff"], (r["ci_lo"], r["ci_hi"]))))
    for id_, q in (("P2", "Does the learned model add beyond the planner? (v2 − MPPI + nominal physics)"),
                   ("P5", "Do the learned corrections matter? (v2 − MPPI + v2's physics prior)")):
        r = g(id_)
        if r is None:
            row(id_, "RQ1", q, "not run", False)
            continue
        met = r["diff"] > 0 and r["ci_lo"] > 0 and seeds_ok(r)
        hs = f"{int(r['seeds_holding'])}/{int(r['n_seeds'])}"
        row(id_, "RQ1", q, f"success@5km {r['diff']:+.1f} pp [{r['ci_lo']:+.1f}, {r['ci_hi']:+.1f}], holds in "
            f"{hs} seeds", met, r["diff"], (r["ci_lo"], r["ci_hi"]), r["p_boot"], hs)
    # P3: equivalence per atmosphere
    rs = {t: g("P3", "success5", t) for t in ATM}
    if all(v is not None for v in rs.values()):
        ok = all(-cfg.equiv_pp < v["ci90_lo"] and v["ci90_hi"] < cfg.equiv_pp for v in rs.values())
        row("P3", "RQ1", "Is v2's model close to a perfect one? (MPPI + true physics − v2, equivalent within ±5 pp)",
            ", ".join(f"{ATM_LAB[t]} {v['diff']:+.1f} (90% [{v['ci90_lo']:+.1f}, {v['ci90_hi']:+.1f}])"
                      for t, v in rs.items()), ok)
    else:
        row("P3", "RQ1", "Is v2's model close to a perfect one?", "not run", False)
    r = g("P4")
    if r is None:
        row("P4", "RQ1", "Does the classical controller suffer from its model? (knowing true density & mass − classical; a lower bound: its predictor keeps the nominal L/D)", "not run", False)
    else:
        row("P4", "RQ1", "Does the classical controller suffer from its model? (knowing true density & mass − classical; a lower bound: its predictor keeps the nominal L/D)",
            f"success@5km {r['diff']:+.1f} pp [{r['ci_lo']:+.1f}, {r['ci_hi']:+.1f}]; " + per_atm("P4"),
            r["diff"] > 0 and r["ci_lo"] > 0, r["diff"], (r["ci_lo"], r["ci_hi"]), r["p_boot"])
    # P6
    if len(t06):
        parts = []
        for load in (G, "matched"):
            for metric in ("success5", "log_miss"):
                x = t06[(t06.load == load) & (t06.metric == metric) & (t06.atmosphere == "all")]
                if len(x):
                    x = x.iloc[0]
                    parts.append(f"{metric} @{load}: planner {x.planner_share:.2f} [{x.planner_share_ci_lo:.2f}, "
                                 f"{x.planner_share_ci_hi:.2f}], model {x.model_share:.2f} [{x.model_share_ci_lo:.2f}, "
                                 f"{x.model_share_ci_hi:.2f}]")
        row("P6", "RQ1", "Decomposition: planner and model shares of the classical → v2 gap", "; ".join(parts), None)
    else:
        row("P6", "RQ1", "Decomposition", "not run", None)
    # L1, L2, L3
    l1 = None
    if present(df, V2m) and present(df, CL):
        gaps = {t: g_mean.get(("v2", "matched", t), np.nan) - g_mean.get(("classical", "own", t), np.nan) for t in ATM}
        l1 = all(np.isfinite(v) and v <= cfg.l1_tol_g for v in gaps.values())
        row("L1", "RQ2", "Is the matched arm really matched? (v2@matched mean peak g ≤ classical's + 0.1 g)",
            ", ".join(f"{ATM_LAB[t]} {v:+.2f} g" for t, v in gaps.items()), l1)
    else:
        row("L1", "RQ2", "Is the matched arm really matched?", "not run", False)
    r = g("total@matched")
    if r is None:
        row("L2", "RQ2", "Does v2 beat classical at matched load?", "not run", False)
    else:
        worse = [t for t in ATM if (rr := g("total@matched", "success5", t)) is not None and rr["ci_hi"] < 0]
        met = r["diff"] > 0 and r["ci_lo"] > 0 and seeds_ok(r) and not worse
        hs = f"{int(r['seeds_holding'])}/{int(r['n_seeds'])}"
        row("L2", "RQ2", "Does v2 beat classical at matched load? (v2@matched − classical)",
            f"success@5km {r['diff']:+.1f} pp [{r['ci_lo']:+.1f}, {r['ci_hi']:+.1f}], holds in {hs} seeds; "
            f"classical significantly better in: {', '.join(worse) if worse else 'none'}"
            + ("" if l1 else " (FLAGGED: L1 not met, the comparison is not load-matched)"),
            met, r["diff"], (r["ci_lo"], r["ci_hi"]), r["p_boot"], hs, valid=l1)
    lc = t07[(t07.atmosphere == "all") & (t07.arm_key == "v2")] if len(t07) else t07
    row("L3", "RQ2", "What does load buy? (v2: success@5km and mean log miss against mean peak g)",
        "; ".join(f"{x.load}: {x.peak_g_mean:.2f} g, {x['success5_%']:.1f}%, log miss {x.log_miss_mean:.2f}"
                  for _, x in lc.iterrows()) or "not run", None)
    # A / B
    for _, x in t08.iterrows():
        v = x.variant
        if not x.get("ran", False):
            row(v, "RQ3", f"{'Necessary' if v.startswith('A') else 'Sufficient'}: {x.change}", "not run", False)
            continue
        k, n = int(x.seeds_holding), int(x.n_seeds)
        if v.startswith("A"):
            met = x.d_log_miss > 0 and x.ci_lo > 0 and k >= max(1, n - 1)
            meas = (f"log miss vs R0 {x.d_log_miss:+.3f} [{x.ci_lo:+.3f}, {x.ci_hi:+.3f}] (×{x.geo_miss_ratio:.2f} "
                    f"geometric-mean miss), holds in {k}/{n} seeds")
        else:
            met = x.d_log_miss < 0 and x.ci_hi < 0 and k >= max(1, n - 1)
            gc = x.get("gap_closed", np.nan)
            own = bool(np.isfinite(gc) and gc >= cfg.gap_closed_min and x.get("gap_closed_ci_lo", np.nan) > 0)
            meas = (f"log miss vs v1 {x.d_log_miss:+.3f} [{x.ci_lo:+.3f}, {x.ci_hi:+.3f}], holds in {k}/{n} seeds; "
                    f"gap closed {gc:.2f} [{x.get('gap_closed_ci_lo', np.nan):.2f}, {x.get('gap_closed_ci_hi', np.nan):.2f}]"
                    f" -> {'sufficient on its own' if own else 'not sufficient on its own'}")
        row(v, "RQ3", f"{'Necessary' if v.startswith('A') else 'Sufficient'}: {x.change}", meas, met, x.d_log_miss,
            (x.ci_lo, x.ci_hi), x.p_boot, f"{k}/{n}")
    # R0-check
    if present(df, ("R0", G)) and present(df, V2S):
        d = {t: 100 * (np.nanmean(cube(df, ("R0", G), "success_5km", seeds, cases)[:, i])
                       - np.nanmean(cube(df, V2S, "success_5km", seeds, cases)[:, i])) for i, t in enumerate(ATM)}
        row("R0-check", "check", "Does retraining reproduce v3's v2single? (R0 − v2single within ±3 pp per atmosphere)",
            ", ".join(f"{ATM_LAB[t]} {v:+.1f} pp" for t, v in d.items()),
            all(np.isfinite(v) and abs(v) <= cfg.r0_check_pp for v in d.values()))
    else:
        row("R0-check", "check", "Does retraining reproduce v3's v2single?", "not run", False)
    # M
    if mech.get("rho") is not None and np.isfinite(mech["rho"]):
        row("M", "RQ3", "Mechanism: does open-loop error predict closed-loop log miss across the 15 models?",
            f"Spearman rho {mech['rho']:+.2f} [{mech['ci_lo']:+.2f}, {mech['ci_hi']:+.2f}] over {mech['n_models']} models "
            f"(open loop at {cfg.mechanism_ttg_s} s, test split)", None, mech["rho"], (mech["ci_lo"], mech["ci_hi"]))
    else:
        row("M", "RQ3", "Mechanism: open-loop error vs closed-loop log miss", "not run", None)
    # Holm within families
    for fam, ids in (("RQ1", ("P2", "P4", "P5")), ("RQ2", ("L2",)),
                     ("RQ3", tuple(v for v in cfg.variants if v != "R0"))):
        idx = [i for i, r in enumerate(rows) if r["id"] in ids and r["p"] is not None and np.isfinite(r["p"])]
        if idx:
            for i, pa in zip(idx, A3.holm([rows[i]["p"] for i in idx])):
                rows[i]["p_holm"] = pa
                rows[i]["met"] = bool(rows[i]["met"] and pa < 0.05)
    return pd.DataFrame(rows)


def mechanism(cfg, df, ol: pd.DataFrame, seeds) -> dict:
    """Spearman rho across models (R0, A, B, v1) between open-loop median error and closed-loop mean log miss."""
    models = [m for m in list(cfg.variants) + ["v1"]]
    x = ol[(ol.split == "test") & (ol.ttg_s == cfg.mechanism_ttg_s)] if len(ol) else ol
    rows = []
    for m in models:
        for s in seeds:
            xo = x[(x.model == m) & (x.seed == s)].err_median_km if len(x) else []
            d = df[(df.arm == m) & (df.load == G) & (df.seed == s)]
            if len(xo) and len(d):
                rows.append({"model": m, "seed": s, "ol_err_km": float(xo.iloc[0]), "log_miss": float(d.log_miss.mean())})
    tab = pd.DataFrame(rows)
    if not len(tab) or tab.model.nunique() < 3:
        return {"rho": None, "table": tab, "n_models": int(tab.model.nunique()) if len(tab) else 0}
    piv_x = tab.pivot_table(index="model", columns="seed", values="ol_err_km")
    piv_y = tab.pivot_table(index="model", columns="seed", values="log_miss")
    rho = A3.spearman(piv_x.mean(1), piv_y.mean(1))
    rng = np.random.default_rng(zlib.crc32(b"mechanism"))
    S_ = piv_x.shape[1]
    vals = []
    for _ in range(cfg.n_boot):
        si = rng.integers(0, S_, S_)
        vals.append(A3.spearman(piv_x.iloc[:, si].mean(1), piv_y.iloc[:, si].mean(1)))
    v = np.asarray(vals, dtype=float)
    v = v[np.isfinite(v)]
    lo, hi = np.percentile(v, [2.5, 97.5]) if len(v) else (np.nan, np.nan)
    return {"rho": rho, "ci_lo": float(lo), "ci_hi": float(hi), "n_models": int(piv_x.shape[0]), "table": tab}


# ============================================================================
# T03, T11, T12, T13
# ============================================================================
def e0_table(cfg, L: S4.LayoutV4) -> pd.DataFrame:
    rows = []
    rd = lambda p: pd.read_csv(p) if p.exists() else None
    a = rd(L.checks / "E0a.csv")
    if a is not None and len(a):
        ea = a[a.check == "E0a"]
        rows.append({"check": "E0a", "item": "classical re-flown == v3 (outcome, miss 1e-6 km, fuel)",
                     "value": f"{int(ea.identical.sum())}/{len(ea)} identical, max |d miss| {ea.d_miss_km.max():.2e} km",
                     "pass": bool(ea.identical.all())})
        o = a[a.check == "oracle_vs_classical_nominal"]
        if len(o):
            rows.append({"check": "E0a", "item": "classical_oracle == classical in the nominal atmosphere",
                         "value": f"{int(o.identical.sum())}/{len(o)} identical, max |d miss| {o.d_miss_km.max():.2e} km",
                         "pass": bool(o.identical.all())})
    p = L.checks / "E0b_summary.json"
    if p.exists():
        s = json.loads(p.read_text())
        for q in s.get("pairs", []):
            rows.append({"check": "E0b", "item": f"{q['v4_arm']} vs v3 {q['v3_arm']} (seed {cfg.e0b_seed})",
                         "value": f"outcome agreement {100 * q['outcome_agreement']:.1f}%, mean d log miss "
                                  f"{q['mean_d_log_miss']:+.4f} [{q['ci_lo']:+.4f}, {q['ci_hi']:+.4f}], sd "
                                  f"{q['sd_d_log_miss']:.4f}, 5-95% [{q['p05']:+.3f}, {q['p95']:+.3f}], identical misses "
                                  f"{100 * q['frac_identical_miss']:.0f}%",
                         "pass": bool(q["outcome_agreement"] >= cfg.e0b_agree_min and q["ci_lo"] <= 0 <= q["ci_hi"])})
    c = rd(L.checks / "E0c.csv")
    if c is not None and len(c):
        rows.append({"check": "E0c", "item": "run_lockstep_v4 == run_lockstep_v3 on the CPU (v1, v2 member 0)",
                     "value": f"{int(c.identical.sum())}/{len(c)} flights identical", "pass": bool(c.identical.all())})
    p = L.checks / "E0_physics_port.json"
    if p.exists():
        s = json.loads(p.read_text())
        rows.append({"check": "physics 1", "item": f"torch port vs math_3d.nominal_eom_step ({s['n']} states, float64)",
                     "value": f"max relative error per step {s['max_rel_err']:.1e} (increment {s['max_rel_err_increment']:.1e}, "
                              f"load {s['max_rel_err_gload']:.1e})", "pass": bool(s["max_rel_err"] < 1e-9)})
    f = rd(L.checks / "E0_physics_fidelity.csv")
    if f is not None and len(f):
        for ttg, g in f.groupby("ttg_s"):
            rows.append({"check": "physics 2", "item": f"open-loop landing error, {int(ttg)} s to go ({len(g)} flights)",
                         "value": f"true parameters: median {g.err_truth_km.median():.2f} km, p90 "
                                  f"{g.err_truth_km.quantile(0.9):.2f} km; nominal: median {g.err_nominal_km.median():.2f}"
                                  f" km, p90 {g.err_nominal_km.quantile(0.9):.2f} km", "pass": None})
    return pd.DataFrame(rows, columns=["check", "item", "value", "pass"])


def _track_sep(t1: Dict[str, np.ndarray], t2: Dict[str, np.ndarray], km: float):
    """First time two 1 Hz flights' ground positions are more than km apart: (t, h, V) of the first, or None."""
    from jepa.evaluate import haversine_km
    n = min(len(t1["t_s"]), len(t2["t_s"]))
    d = haversine_km(np.radians(t1["phi_deg"][:n].astype(float)), np.radians(t1["lambda_deg"][:n].astype(float)),
                     np.radians(t2["phi_deg"][:n].astype(float)), np.radians(t2["lambda_deg"][:n].astype(float)))
    k = np.nonzero(d > km)[0]
    if not len(k):
        return None
    i = int(k[0])
    return float(t1["t_s"][i]), float(t1["h_m"][i]), float(t1["V_mps"][i])


def case10_table(cfg, L: S4.LayoutV4, L3: S.Layout, df: pd.DataFrame) -> pd.DataFrame:
    c = cfg.diag_case
    rows = []
    d = df[(df.case == c) & (df.mppi_rep <= 0)]
    for (arm, load, t), g in d.groupby(["arm", "load", "atmosphere"]):
        rows.append({"part": "by_arm", "arm": arm, "load": load, "atmosphere": t, "seed": None, "n": len(g),
                     "success5_%": 100 * g.success_5km.mean(), "miss_median_km": g.miss_km.median()})
    diag = d[d.arm == "diag_case10"]
    for r in diag.itertuples():
        p = L.telemetry_path(r.flight_id)
        if not p.exists():
            continue
        tel = S.load_telemetry(L, r.flight_id)
        cfid = f"classical/{r.atmosphere}_case{c:02d}"
        sep = _track_sep(tel, S.load_telemetry(L3, cfid), cfg.separation_km) if L3.telemetry_path(cfid).exists() else None
        n = len(tel.get("dec_t_s", []))
        last = slice(max(0, n - 20), n)
        from jepa.evaluate import haversine_km
        pred_err = haversine_km(np.radians(tel["dec_phi_pred_deg"][last].astype(float)),
                                np.radians(tel["dec_lambda_pred_deg"][last].astype(float)),
                                math.radians(float(tel["phi_deg"][-1])), math.radians(float(tel["lambda_deg"][-1]))) \
            if n else np.array([np.nan])
        pk = tel["dec_peak_g_pred"].astype(float) if n else np.array([np.nan])
        rows.append({"part": "diag", "arm": r.arm, "load": r.load, "atmosphere": r.atmosphere, "seed": int(r.seed), "n": 1,
                     "success5_%": 100.0 * bool(r.success_5km), "miss_median_km": r.miss_km, "outcome": r.outcome,
                     "g_limit": r.g_limit, "peak_g": r.peak_g,
                     "sep_t_s": sep[0] if sep else np.nan, "sep_h_km": sep[1] / 1000 if sep else np.nan,
                     "sep_V_mps": sep[2] if sep else np.nan,
                     "pred_miss_last20_km": float(np.nanmean(tel["dec_miss_best_km"][last])) if n else np.nan,
                     "pred_err_last20_km": float(np.nanmean(pred_err)),
                     "frac_load_binding": float(np.mean(pk >= float(r.g_limit) - 0.1)) if n else np.nan,
                     "peak_g_pred_max": float(np.nanmax(pk)) if n else np.nan, "decisions": n})
    return pd.DataFrame(rows)


def threshold_table(cfg, main: pd.DataFrame) -> pd.DataFrame:
    rows = []
    th = np.asarray(cfg.thresholds_km, dtype=float)
    for (arm, load), g in main.groupby(["arm", "load"]):
        for t, gg in [("all", g)] + [(t, g[g.atmosphere == t]) for t in ATM]:
            if not len(gg):
                continue
            ok = (gg.outcome == "drogue").to_numpy()[:, None] & (gg.miss_km.to_numpy()[:, None] < th[None])
            for k, v in zip(th, ok.mean(0)):
                rows.append({"arm": arm, "load": load, "atmosphere": t, "threshold_km": float(k), "success_%": 100 * v})
    return pd.DataFrame(rows)


def ensemble_rescore(cfg, df, cb, boot: Boot) -> pd.DataFrame:
    """v3's v2 vs v2single with the v4 metrics (section 4.5)."""
    if not (present(df, V2) and present(df, V2S)):
        return pd.DataFrame(columns=["metric", "atmosphere", "diff", "ci_lo", "ci_hi", "p_boot"])
    rows = []
    for metric, col, fn, sc in (("log_miss", "log_miss", "paired_mean", 1.0), ("success2", "success_2km", "diff_mean", 100.0),
                                ("success5", "success_5km", "diff_mean", 100.0)):
        XA, XB = cb(V2, col), cb(V2S, col)
        for i, t in [(None, "all")] + list(enumerate(ATM)):
            r = boot(f"ens|{metric}", [XA, XB], fn, i)
            rows.append({"metric": metric, "atmosphere": t, "diff": r["diff"] * sc, "ci_lo": r["ci_lo"] * sc,
                         "ci_hi": r["ci_hi"] * sc, "p_boot": r["p_boot"]})
    return pd.DataFrame(rows)


# ============================================================================
# Figures
# ============================================================================
def figures(cfg, L, L3, caps, df, main, t01, per_seed, t04, t06, t07, t08, ol, mech, seeds, cases) -> Dict[str, Path]:
    plt = _plt()
    figs: Dict[str, Path] = {}
    x = np.arange(len(ATM))

    # F10 the E1 grid
    arms = [a for a in E1_ARMS if present(main, a)]
    fig, ax = plt.subplots(figsize=(11, 3.6))
    w = 0.8 / max(1, len(arms))
    for j, al in enumerate(arms):
        r = t01[(t01.arm == al[0]) & (t01.load == al[1])].set_index("atmosphere").reindex(ATM)
        xs = x + (j - (len(arms) - 1) / 2) * w
        ax.bar(xs, r["success5_%"], w * 0.9, color=color(al[0]), label=label(*al))
        if al[0] not in DETERMINISTIC:
            ax.errorbar(xs, r["success5_%"], yerr=[r["success5_%"] - r["success5_%_min"], r["success5_%_max"] - r["success5_%"]],
                        fmt="none", ecolor="#0b0b0b", lw=0.8, capsize=2)
            ps = per_seed[(per_seed.arm == al[0]) & (per_seed.load == al[1])]
            for i, t in enumerate(ATM):
                v = ps[ps.atmosphere == t]["success5_%"].values
                ax.scatter(np.full(len(v), xs[i]), v, s=5, color="#0b0b0b", zorder=3)
    ax.set_xticks(x, [ATM_LAB[t] for t in ATM])
    ax.set(ylim=(0, 105), ylabel="success@5km [%]", title="E1: planner × model at 4.5 g (bars: mean over seeds; "
                                                           "whiskers: seed range; dots: seeds)")
    ax.legend(ncol=6, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.08), frameon=False)
    fig.tight_layout()
    figs["F10"] = _fig(L, caps, plt, fig, "F10_success_by_arm_and_atmosphere",
                       "success@5km of the planner × model grid at 4.5 g: classical (nominal model), classical knowing the true "
                       "density and mass (not L/D), MPPI through the nominal physics, through v2's physics prior, through v2's learned model, "
                       "and through the true physics.", "05_analysis/tables/T01_metrics_by_arm.csv")

    # F11 gap decomposition (waterfall per atmosphere)
    d = t06[(t06.metric == "success5") & (t06.load == G)] if len(t06) else t06
    fig, ax = plt.subplots(figsize=(11, 3.6))
    for i, t in enumerate(["all"] + ATM):
        r = d[d.atmosphere == t]
        if not len(r):
            continue
        r = r.iloc[0]
        base = r.classical
        x0 = i * 5
        ax.bar(x0, base, 0.8, color=color("classical"))
        ax.bar(x0 + 1, r.planner, 0.8, bottom=base, color=color("mppi_nominal"))
        ax.bar(x0 + 2, r.model, 0.8, bottom=base + r.planner, color=color("v2"))
        if np.isfinite(r.ceiling_truth_minus_v2):
            ax.bar(x0 + 3, r.ceiling_truth_minus_v2, 0.8, bottom=base + r.total, color=color("mppi_truth"), alpha=0.6)
        ax.text(x0 + 1.5, -9, "pooled" if t == "all" else ATM_LAB[t], ha="center", fontsize=7)
    for c_, l_, al in (("classical", "classical", 1.0), ("mppi_nominal", "+ planner (MPPI, nominal physics)", 1.0),
                       ("v2", "+ learned model (= v2)", 1.0), ("mppi_truth", "+ remaining model gap (true physics)", 0.6)):
        ax.bar([0], [0], color=color(c_), alpha=al, label=l_)
    ax.set(ylim=(-12, 105), ylabel="success@5km [%]", xticks=[],
           title="Where the classical → v2 gap comes from (4.5 g; a step down is a loss)")
    ax.legend(fontsize=7, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.04), frameon=False)
    fig.tight_layout()
    figs["F11"] = _fig(L, caps, plt, fig, "F11_gap_decomposition",
                       "The classical → v2 gap at 4.5 g, stacked: classical, + what the planner adds with the same nominal "
                       "model, + what the learned model adds with the same planner (= v2), + what a perfect model would "
                       "still add (MPPI + true physics). Pooled and per atmosphere.", "05_analysis/tables/T06_gap_decomposition.csv")

    # F12 success vs mean peak g
    fig, ax = plt.subplots(1, 4, figsize=(13, 3.1), sharey=True)
    for i, t in enumerate(ATM):
        lc = t07[(t07.atmosphere == t)] if len(t07) else t07
        v2c = lc[lc.arm_key == "v2"].sort_values("peak_g_mean") if len(lc) else lc
        if len(v2c):
            ax[i].plot(v2c.peak_g_mean, v2c["success5_%"], marker="o", color=color("v2"), label="v2 (loads)")
            for _, r in v2c.iterrows():
                ax[i].annotate(r.load, (r.peak_g_mean, r["success5_%"]), fontsize=6, xytext=(3, 3), textcoords="offset points")
        for key, mk in (("classical", "*"), ("mppi_nominal", "s")):
            rr = lc[lc.arm_key == key] if len(lc) else lc
            for _, r in rr.iterrows():
                m_ = "D" if r.load == "matched" else mk
                ax[i].scatter([r.peak_g_mean], [r["success5_%"]], marker=m_, s=60 if m_ == "*" else 25, color=color(key),
                              label=label(key, r.load), zorder=3)
        ax[i].set(xlabel="mean peak g", title=ATM_LAB[t])
    ax[0].set_ylabel("success@5km [%]")
    ax[0].legend(fontsize=6)
    fig.tight_layout()
    figs["F12"] = _fig(L, caps, plt, fig, "F12_success_vs_peak_g",
                       "L3, what load buys: success@5km against mean peak g for v2 at the matched limit, 3.5, 4.0, 4.5 and 8 g "
                       "(8 g: seed 0), with the classical controller (star) and MPPI + nominal physics (square: 4.5 g; "
                       "diamond: matched load).",
                       "05_analysis/tables/T07_load_curve.csv")

    # F13 ablation forest
    tt = t08[t08.get("ran", False) == True] if len(t08) else t08  # noqa: E712
    fig, ax = plt.subplots(figsize=(8, 0.9 + 0.32 * max(1, len(tt))))
    for i, (_, r) in enumerate(tt.iterrows()):
        y = -i
        sig = bool(r.get("p_holm", 1.0) < 0.05) if np.isfinite(r.get("p_holm", np.nan)) else False
        col = "#c2410c" if r.variant.startswith("A") else "#2a78d6"
        ax.plot([r.ci_lo, r.ci_hi], [y, y], color=col, lw=1.6)
        ax.scatter([r.d_log_miss], [y], s=30, color=col if sig else "white", edgecolor=col, zorder=3)
    ax.axvline(0, color="#0b0b0b", lw=0.7)
    ax.set_yticks([-i for i in range(len(tt))], [f"{r.variant} (vs {r.reference})" for _, r in tt.iterrows()], fontsize=7)
    ax.set(xlabel="paired log-miss difference, variant − reference (> 0: misses more; 95% CI)",
           title="Ablations: A remove one v2 change from R0, B add one to v1 (filled: Holm p < 0.05)")
    fig.tight_layout()
    figs["F13"] = _fig(L, caps, plt, fig, "F13_ablation_forest",
                       "Every ablation's paired log-miss difference against its reference with 95% intervals; filled when "
                       "significant after Holm over the 13 variants. A variants above 0 = the change is necessary; B "
                       "variants below 0 = the change helps v1.", "05_analysis/tables/T08_ablations.csv")

    # F14 size curve
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.2))
    sizes = [v for v in ("A7b_size_33k", "A7a_size_250k", "R0") if v in cfg.variants]
    from jepa import ablate_v4 as AB
    for k, (what, yl) in enumerate((("log_miss", "closed loop: mean log miss"), ("ol", f"open loop: median error at "
                                                                                     f"{cfg.mechanism_ttg_s} s [km]"))):
        xs, ys, lo, hi = [], [], [], []
        for v in sizes:
            n = AB.param_count(cfg, v)["per_member"]
            if what == "log_miss":
                vals = [main[(main.arm == v) & (main.load == G) & (main.seed == s)].log_miss.mean() for s in seeds]
            else:
                o = ol[(ol.model == v) & (ol.split == "test") & (ol.ttg_s == cfg.mechanism_ttg_s)] if len(ol) else ol
                vals = list(o.err_median_km) if len(o) else []
            vals = [x_ for x_ in vals if np.isfinite(x_)]
            if vals:
                xs.append(n), ys.append(np.mean(vals)), lo.append(np.mean(vals) - min(vals)), hi.append(max(vals) - np.mean(vals))
        if xs:
            ax[k].errorbar(xs, ys, yerr=[lo, hi], marker="o", color="#c2410c", capsize=3)
            for x_, y_, v in zip(xs, ys, [v for v in sizes]):
                ax[k].annotate(v.split("_")[0], (x_, y_), fontsize=7, xytext=(4, 4), textcoords="offset points")
        ax[k].set(xscale="log", xlabel="parameters per member", ylabel=yl)
    ax[0].set_title("Size curve (whiskers: seed range)")
    fig.tight_layout()
    figs["F14"] = _fig(L, caps, plt, fig, "F14_size_curve",
                       "Three model sizes about a decade apart (A7b ≈ 33k, A7a ≈ 250k, R0 ≈ 2.45M parameters per member): "
                       "closed-loop mean log miss and open-loop error against size.", "T08, T10")

    # F15 mechanism
    fig, ax = plt.subplots(figsize=(6, 4))
    tab = mech.get("table")
    if tab is not None and len(tab):
        m_ = tab.groupby("model")[["ol_err_km", "log_miss"]].mean()
        for name, r in m_.iterrows():
            ax.scatter([r.ol_err_km], [r.log_miss], color="#2a78d6" if name.startswith("B") or name == "v1" else "#c2410c")
            ax.annotate(name.split("_")[0], (r.ol_err_km, r.log_miss), fontsize=6.5, xytext=(3, 3), textcoords="offset points")
    rho = mech.get("rho")
    ax.set(xscale="log", xlabel=f"open-loop median landing error at {cfg.mechanism_ttg_s} s (test split) [km]",
           ylabel="closed-loop mean log miss",
           title=f"Mechanism: Spearman rho = {rho:+.2f} [{mech['ci_lo']:+.2f}, {mech['ci_hi']:+.2f}]"
           if rho is not None and np.isfinite(rho) else "Mechanism (not enough models yet)")
    fig.tight_layout()
    figs["F15"] = _fig(L, caps, plt, fig, "F15_open_vs_closed_loop",
                       "Across the 15 models (R0, A1–A8, B1–B3, v1; mean over seeds): open-loop prediction error against "
                       "closed-loop log miss. Do the design changes act through prediction quality?", "T10, flights")

    # F16 success vs threshold
    th = threshold_table(cfg, main)
    pick = [CL, OR, MN, MT, MP, V1, V2, V2S, V2m, ("R0", G)]
    fig, ax = plt.subplots(1, 4, figsize=(13, 3.1), sharey=True)
    for i, t in enumerate(ATM):
        for al in pick:
            d = th[(th.arm == al[0]) & (th.load == al[1]) & (th.atmosphere == t)]
            if len(d):
                ax[i].plot(d.threshold_km, d["success_%"], color=color(al[0]), lw=1.2, label=label(*al),
                           ls="--" if al[1] == "matched" else "-")
        ax[i].axvline(5, color="#0b0b0b", lw=0.6, ls=":")
        ax[i].set(xlabel="miss threshold [km]", title=ATM_LAB[t])
    ax[0].set_ylabel("success [%]")
    ax[0].legend(fontsize=5.5)
    fig.tight_layout()
    figs["F16"] = _fig(L, caps, plt, fig, "F16_success_vs_threshold",
                       "Success rate against the miss threshold (0.5–25 km), per atmosphere, pooled over seeds.",
                       "05_analysis/tables/T12_success_vs_threshold.csv")

    # F17 case 10
    figs["F17"] = case10_figure(cfg, L, L3, caps, plt, df)

    # F18 miss-by-case heatmap of the new arms
    new = [CL, OR, MN, MNm, MT, MP, ("v2", "4g"), ("v2", "3.5g"), V2m, V2]
    hm = main.groupby(["case", "arm", "load", "atmosphere"]).miss_km.median().unstack(["arm", "load", "atmosphere"])
    cols = [(a, l_, t) for a, l_ in new for t in ATM if (a, l_, t) in hm.columns]
    fig, ax = plt.subplots(figsize=(min(14, 1.0 + 0.36 * max(1, len(cols))), 1.6 + 0.12 * max(1, len(hm.index))))
    if cols:
        hm = hm.reindex(columns=pd.MultiIndex.from_tuples(cols))
        im = ax.imshow(np.log10(hm.values.clip(0.05, 2000)), aspect="auto", cmap="RdYlGn_r", vmin=-1, vmax=2)
        ax.set_xticks(range(len(cols)), [f"{label(a, l_)}\n{t[:8]}" for a, l_, t in cols], fontsize=5, rotation=90)
        ax.set_yticks(range(len(hm.index)), hm.index, fontsize=5)
        fig.colorbar(im, ax=ax, shrink=0.6).set_label("log10 median miss [km]")
    ax.set(ylabel="case", title="Median miss over seeds, per case: the E1/E2 arms (v2 at 4.5 g for reference)")
    fig.tight_layout()
    figs["F18"] = _fig(L, caps, plt, fig, "F18_miss_by_case_new_arms",
                       "Median miss over seeds for every case, atmosphere and new arm (as v3's F15). Green ≤ 1 km, red ≥ 100 km.",
                       "flights")
    return figs


def case10_figure(cfg, L, L3, caps, plt, df) -> Path:
    c = cfg.diag_case
    diag = df[(df.arm == "diag_case10") & (df.case == c)].sort_values(["atmosphere", "load", "seed"])
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.6))
    cols = {"v2_4.5g": "#eb6834", "v2_8g": "#f4a582", "v2single_4.5g": "#eda100", "v2single_8g": "#f7d48a"}
    for t in cfg.diag_atmospheres:
        fid = f"classical/{t}_case{c:02d}"
        if L3.telemetry_path(fid).exists():
            tel = S.load_telemetry(L3, fid)
            ax[0].plot(tel["lambda_deg"], tel["phi_deg"], color=color("classical"), lw=1.4,
                       ls="-" if t == cfg.diag_atmospheres[0] else "--", label=f"classical {t[:8]}")
    for r in diag.itertuples():
        if not L.telemetry_path(r.flight_id).exists():
            continue
        tel = S.load_telemetry(L, r.flight_id)
        ls = "-" if r.atmosphere == cfg.diag_atmospheres[0] else "--"
        lab = f"{r.load} {r.atmosphere[:8]}" if r.seed == cfg.seeds[0] else None
        ax[0].plot(tel["lambda_deg"], tel["phi_deg"], color=cols.get(r.load, "#52514e"), lw=0.8, ls=ls, label=lab)
        if len(tel.get("dec_t_s", [])):
            ax[1].plot(tel["dec_t_s"], tel["dec_miss_best_km"], color=cols.get(r.load, "#52514e"), lw=0.8, ls=ls)
            ax[1].axhline(r.miss_km, color=cols.get(r.load, "#52514e"), lw=0.4, ls=":")
            ax[2].plot(tel["dec_t_s"], tel["dec_peak_g_pred"], color=cols.get(r.load, "#52514e"), lw=0.8, ls=ls)
    ax[0].scatter([15], [15], marker="*", s=90, color="#0b0b0b", zorder=3)
    for g in (cfg.g_primary, cfg.g_legacy):
        ax[2].axhline(g, color="#0b0b0b", lw=0.6, ls=":")
    ax[0].set(xlabel="longitude [deg]", ylabel="latitude [deg]", title=f"case {c}: ground tracks (dashed: "
                                                                        f"{cfg.diag_atmospheres[-1]})")
    ax[1].set(yscale="log", xlabel="time [s]", ylabel="predicted miss of the chosen plan [km]",
              title="Predicted miss per decision (dotted: the real miss)")
    ax[2].set(xlabel="time [s]", ylabel="predicted peak g of the chosen plan", title="Load: predicted peak g vs limits")
    ax[0].legend(fontsize=6)
    fig.tight_layout()
    return _fig(L, caps, plt, fig, "F17_case10_diagnosis",
                f"Case {c} diagnostic flights (v2 ensemble and single member, 4.5 and 8 g, every seed; solid: "
                f"{cfg.diag_atmospheres[0]}, dashed: {cfg.diag_atmospheres[-1]}): ground tracks, the planner's predicted "
                "miss at each decision against the real miss, and the chosen plan's predicted peak g against the limits.",
                "04_flights/telemetry/diag_case10/ (decision log)")
