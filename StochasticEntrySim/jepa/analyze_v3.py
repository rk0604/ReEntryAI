"""
jepa/analyze_v3.py

Everything the v3 study reports (jepa/v3_plan.md sections 5-8), written to
05_analysis/ as numbered tables (T##) and figures (F##) with captions.json,
and the PDF report in 06_report/.

    testset_figures       F01 coverage, F02 footprint margins (01_testset/)
    open_loop_seed        per model and seed: footprint errors, along-track
                          bias, ensemble spread (03_open_loop/)
    analyze               flights.csv -> metric tables, hierarchical bootstrap,
                          per-seed McNemar / sign tests, Holm, criteria C1-C7 +
                          Bridge, load-matching cost, planner-noise share,
                          RCS analyses, figures
    build_report          06_report/SpaceJEPA_v3_report.pdf

Definitions (plan section 5): success@5km = outcome drogue and miss < 5 km;
median / p90 miss = 50th / 90th percentile (linear interpolation) of the
miss distance over a group; every failure counts with its miss at termination.
Groups are summarized per seed over the 64 cases, then mean and range over
seeds, and pooled over all seeds.
"""

from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

import numpy as np
import pandas as pd

from jepa import study_v3 as S
from jepa import testset_v3 as TS

ATM = list(TS.ATMOSPHERES)
ATM_LAB = {"nominal": "nominal", "mismatch": "thin ×0.88", "mismatch_high": "dense ×1.12",
           "out_of_range": "×0.78 / L/D ×0.80"}
CTRL = ["classical", "v1", "v2", "v2single", "classical_oracle"]
LAB = {"classical": "classical", "v1": "v1 (fair budget)", "v2": "v2 ensemble", "v2single": "v2 single member",
       "classical_oracle": "classical, knows truth"}
COL = {"classical": "#6250d6", "v1": "#2a78d6", "v2": "#eb6834", "v2single": "#eda100",
       "classical_oracle": "#a3a29c"}
BRIDGE_REF_MEDIAN_KM = {"nominal": 0.90, "mismatch": 1.29, "mismatch_high": 1.89, "out_of_range": 2.42}


def _plt():
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                         "grid.color": "#e4e3df", "grid.linewidth": 0.6, "legend.frameon": False,
                         "axes.titlesize": 8.5, "axes.edgecolor": "#8a8984", "axes.labelcolor": "#52514e",
                         "xtick.color": "#52514e", "ytick.color": "#52514e"})
    return plt


class Captions:
    """05_analysis/captions.json: number -> caption and data source."""

    def __init__(self, layout: S.Layout):
        self.path = layout.analysis / "captions.json"
        self.d = json.loads(self.path.read_text()) if self.path.exists() else {}

    def add(self, key: str, caption: str, data: str) -> None:
        self.d[key] = {"caption": caption, "data": data}
        S.write_json(self.path, self.d)


def _fig(layout: S.Layout, caps: Captions, plt, fig, name: str, caption: str, data: str, folder=None) -> Path:
    path = Path(folder or layout.figures) / f"{name}.png"
    S.save_fig(fig, path)
    plt.close(fig)
    caps.add(name, caption, data)
    return path


def _table(layout: S.Layout, caps: Captions, df: pd.DataFrame, name: str, caption: str, data: str) -> Path:
    path = layout.tables / f"{name}.csv"
    S.write_csv(df, path)
    caps.add(name, caption, data)
    return path


# ============================================================================
# Test set
# ============================================================================
def testset_figures(layout: S.Layout, cases: pd.DataFrame, fps: pd.DataFrame) -> List[Path]:
    """F01 coverage of the entry space (new vs bridge), F02 footprint margins."""
    plt = _plt()
    caps = Captions(layout)
    out = []
    new, br = cases[cases.set == "new"], cases[cases.set == "bridge"]
    fig, ax = plt.subplots(1, 4, figsize=(13, 3.1))
    pairs = [("velocity_mps", "gamma_deg", "entry speed [m/s]", "flight-path angle [deg]"),
             ("altitude_m", "velocity_mps", "entry altitude [m]", "entry speed [m/s]"),
             ("downrange_km", "crossrange_km", "target downrange [km]", "target crossrange [km]")]
    for a, (x, y, xl, yl) in zip(ax[:3], pairs):
        a.scatter(new[x], new[y], s=14, color="#eb6834", label="48 new (Sobol)")
        a.scatter(br[x], br[y], s=14, color="#2a78d6", marker="s", label="16 bridge (v1/v2)")
        a.set(xlabel=xl, ylabel=yl)
    ax[0].legend(fontsize=7)
    a = ax[3]
    a.scatter(new.lambda_deg, new.phi_deg, s=14, color="#eb6834")
    a.scatter(br.lambda_deg, br.phi_deg, s=14, color="#2a78d6", marker="s")
    for r in cases.itertuples():
        a.plot([r.lambda_deg, TS.TestsetSpec.target_lambda_deg], [r.phi_deg, TS.TestsetSpec.target_phi_deg],
               color="#a3a29c", lw=0.4, zorder=0)
    a.scatter([TS.TestsetSpec.target_lambda_deg], [TS.TestsetSpec.target_phi_deg], marker="*", s=120,
              color="#0b0b0b", zorder=3)
    a.set(xlabel="longitude [deg]", ylabel="latitude [deg]", title="entry points → target (15N 15E)")
    fig.tight_layout()
    out.append(_fig(layout, caps, plt, fig, "F01_testset_coverage",
                    "The 64 entry cases: 48 new Sobol cases (orange) against the 16 frozen bridge cases (blue). "
                    "Entry state, target position in the entry frame, and entry points on the map.",
                    "01_testset/v3_testset.csv", layout.testset))
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.4), gridspec_kw={"width_ratios": [2.2, 1]})
    m = cases.set_index("case")[[f"margin_km_{t}" for t in ATM]]
    for i, t in enumerate(ATM):
        ax[0].bar(m.index + (i - 1.5) * 0.2, m[f"margin_km_{t}"].clip(-150, 400), 0.2,
                  color=["#6250d6", "#2a78d6", "#eb6834", "#1baf7a"][i], label=ATM_LAB[t])
    ax[0].axhline(0, color="#0b0b0b", lw=0.7)
    ax[0].axhline(50, color="#a3a29c", lw=0.7, ls="--")
    ax[0].set(xlabel="case", ylabel="margin to the footprint edge [km]",
              title="Target inside the 4.5 g footprint (> 0 = reachable; dashed = 50 km sampling margin)")
    ax[0].legend(ncol=4, fontsize=7)
    c = int(new.case.iloc[0])
    for i, t in enumerate(ATM):
        f = fps[(fps.case == c) & (fps.atmosphere == t)].sort_values("sigma_deg")
        ok = f[f.within_load.astype(bool)]
        ax[1].plot(list(ok.downrange_km) + [ok.downrange_km.iloc[0]], list(ok.crossrange_km) + [ok.crossrange_km.iloc[0]],
                   color=["#6250d6", "#2a78d6", "#eb6834", "#1baf7a"][i], lw=1.2, label=ATM_LAB[t])
    r = cases.set_index("case").loc[c]
    ax[1].scatter([r.downrange_km], [r.crossrange_km], marker="*", s=120, color="#0b0b0b", zorder=3)
    ax[1].set(xlabel="downrange [km]", ylabel="crossrange [km]", title=f"case {c}: reachable regions at ≤ 4.5 g")
    ax[1].legend(fontsize=6.5)
    fig.tight_layout()
    out.append(_fig(layout, caps, plt, fig, "F02_footprint_margins",
                    "Left: distance from the target to the edge of each atmosphere's constant-bank footprint "
                    "(banks whose peak load stays within 4.5 g; a conservative estimate of what is reachable). "
                    "Right: the four footprints of one new case and its target.",
                    "01_testset/v3_testset.csv, 01_testset/footprints.csv", layout.testset))
    return out


# ============================================================================
# Open loop (plan section 7)
# ============================================================================
def test_split(cfg: S.StudyConfig, data_path: str) -> List[dict]:
    """The held-out test split of the 6,000 flights (same split as training)."""
    from jepa import data as D
    from jepa import train_v2 as T2
    eps = T2.load_episodes(data_path)
    if cfg.n_episodes:
        eps = eps[:cfg.n_episodes]
    macro = [D.to_macro(e, D.MACRO_DT_S) for e in eps]
    return D.split_episodes(macro, (0.8, 0.1, 0.1), seed=cfg.split_seed)[2]


def open_loop_seed(cfg: S.StudyConfig, layout: S.Layout, test_eps: Sequence[dict], kind: str, seed: int,
                   h_cross: float = 7600.0) -> Optional[pd.DataFrame]:
    """Footprint errors (and v2 member spread) of one model and seed on the test split; cached as CSV."""
    import torch
    from jepa import evaluate_v2 as E2
    from jepa import models as M1
    from jepa import models_v2 as V
    path = layout.open_loop / f"footprint_{kind}_seed{seed}.csv"
    if path.exists():
        return pd.read_csv(path)
    ck = layout.ckpt(kind, seed)
    if not ck.exists():
        return None
    torch.set_grad_enabled(False)
    if kind == "v1":
        m = M1.load_checkpoint(str(ck))
        models = {"v1": {"wm": m["wm"], "prober": m["pi"], "stats": m["stats"], "rf": m["wm"].rf}}
    else:
        m = V.load_checkpoint_v2(str(ck))
        models = {"v2": {**V.as_v1(m), "stats": m["stats"], "rf": m["wm"].rf},
                  "v2 member 0": {**V.as_v1(m, [0]), "stats": m["stats"], "rf": m["wm"].rf}}
    rows = pd.DataFrame(E2.footprint_points(models, test_eps, cfg.ttg_s, h_cross))
    rows["seed"] = seed
    rows["beta_scale"] = rows.dom_density_scale / rows.dom_mass_scale
    rows["density_bin"] = pd.cut(rows.beta_scale, [0, 0.85, 1.15, 10], labels=["thin", "nominal", "dense"]).astype(str)
    S.write_csv(rows, path)
    along = rows.groupby(["model", "ttg_s", "density_bin"]).agg(
        along_mean_km=("along_km", "mean"), along_median_km=("along_km", "median"), err_median_km=("err_km", "median"),
        n=("err_km", "size")).reset_index()
    along["seed"] = seed
    S.write_csv(along, layout.open_loop / f"along_track_{kind}_seed{seed}.csv")
    if kind == "v2":
        sp = pd.DataFrame(E2.member_spread(m, test_eps, cfg.ttg_s, h_cross))
        sp["seed"] = seed
        S.write_csv(sp, layout.open_loop / f"spread_v2_seed{seed}.csv")
    torch.set_grad_enabled(True)
    return rows


def open_loop_summary(cfg: S.StudyConfig, layout: S.Layout) -> Dict[str, object]:
    """T09 (median error per model, seed, time-to-go), F30-F32."""
    plt = _plt()
    caps = Captions(layout)
    files = sorted(layout.open_loop.glob("footprint_*_seed*.csv"))
    if not files:
        return {}
    fp = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    per = fp.groupby(["model", "seed", "ttg_s"]).err_km.median().rename("err_median_km").reset_index()
    t09 = per.pivot_table(index=["model", "ttg_s"], columns="seed", values="err_median_km")
    t09["mean"], t09["min"], t09["max"] = t09.mean(1), t09.min(1), t09.max(1)
    t09 = t09.reset_index()
    _table(layout, caps, t09.round(3), "T09_open_loop_error_by_seed",
           "Open loop: median landing-point error [km] on the held-out test flights, per model, time-to-go and seed.",
           "03_open_loop/footprint_*_seed*.csv")
    out = {"T09": t09}
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.2))
    for name, col in (("v1", COL["v1"]), ("v2", COL["v2"]), ("v2 member 0", COL["v2single"])):
        d = per[per.model == name]
        if not len(d):
            continue
        for s_, g in d.groupby("seed"):
            ax[0].plot(g.ttg_s, g.err_median_km, color=col, lw=0.6, alpha=0.5)
        mm = d.groupby("ttg_s").err_median_km.mean()
        ax[0].plot(mm.index, mm.values, color=col, lw=2, marker="o", ms=3, label=name)
    ax[0].set(yscale="log", xlabel="time-to-go [s]", ylabel="median landing error [km]",
              title="Open-loop error (thin: each seed; thick: mean)")
    ax[0].invert_xaxis()
    ax[0].legend()
    b = fp[fp.ttg_s == (200 if 200 in set(fp.ttg_s) else fp.ttg_s.max())]
    bias = b.groupby(["model", "density_bin", "seed"]).along_km.mean().reset_index()
    bins = ["thin", "nominal", "dense"]
    for j, name in enumerate(("v1", "v2", "v2 member 0")):
        d = bias[bias.model == name]
        for i, bn in enumerate(bins):
            v = d[d.density_bin == bn].along_km.values
            if len(v):
                ax[1].scatter(np.full(len(v), i + (j - 1) * 0.22), v, s=12,
                              color=[COL["v1"], COL["v2"], COL["v2single"]][j], label=name if i == 0 else None)
    ax[1].axhline(0, color="#0b0b0b", lw=0.7)
    ax[1].set_xticks(range(3), bins)
    ax[1].set(ylabel="mean along-track error [km]", title="Along-track bias by air density (200 s to go; dot = seed)")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    out["F30"] = _fig(layout, caps, plt, fig, "F30_open_loop_error_and_bias",
                      "Left: open-loop median landing-point error vs time-to-go, every seed. Right: mean along-track "
                      "error by air-density bin (+ = predicted landing beyond the true one), one dot per seed.",
                      "03_open_loop/footprint_*_seed*.csv")
    sp_files = sorted(layout.open_loop.glob("spread_v2_seed*.csv"))
    if sp_files:
        sp = pd.concat([pd.read_csv(f) for f in sp_files], ignore_index=True)
        cal = sp.groupby(["seed", "ttg_s"]).apply(
            lambda g: pd.Series({"spearman_rho": spearman(g.spread_km, g.err_mean_km),
                                 "spread_median_km": g.spread_km.median(), "err_median_km": g.err_mean_km.median()}),
            include_groups=False).reset_index()
        _table(layout, caps, cal.round(3), "T10_ensemble_spread_calibration",
               "v2 ensemble spread vs actual error of the ensemble-mean landing point (Spearman rho), per seed and "
               "time-to-go.", "03_open_loop/spread_v2_seed*.csv")
        fig, ax = plt.subplots(figsize=(4.8, 3.2))
        for s_, g in cal.groupby("seed"):
            ax.plot(g.ttg_s, g.spearman_rho, marker="o", ms=3, label=f"seed {s_}")
        ax.set(xlabel="time-to-go [s]", ylabel="Spearman rho (spread vs error)", ylim=(-0.2, 1),
               title="Is the ensemble's disagreement informative?")
        ax.invert_xaxis()
        ax.legend(fontsize=6.5)
        out["F31"] = _fig(layout, caps, plt, fig, "F31_spread_calibration",
                          "Rank correlation between the v2 members' disagreement and the actual landing error, per seed.",
                          "03_open_loop/spread_v2_seed*.csv")
        out["T10"] = cal
    return out


# ============================================================================
# Flights: tables and statistics
# ============================================================================
def load_flights(layout: S.Layout, cfg: S.StudyConfig, df: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """flights.csv joined with the test set's feasibility flags."""
    df = S.collect_flights(layout, write=True) if df is None else df
    if not len(df):
        return df
    cases = S.load_testset(cfg)
    feas = cases.melt(id_vars=["case"], value_vars=[f"feasible_{t}" for t in ATM], var_name="atmosphere",
                      value_name="feasible")
    feas["atmosphere"] = feas.atmosphere.str.replace("feasible_", "", regex=False)
    df = df.merge(feas, on=["case", "atmosphere"], how="left")
    df = df.merge(cases[["case", "set", "altitude_m", "velocity_mps", "gamma_deg", "downrange_km", "crossrange_km"]]
                  .rename(columns={"set": "case_set"}), on="case", how="left")
    df["seed"] = df.seed.fillna(-1).astype(int)
    df["mppi_rep"] = df.mppi_rep.fillna(-1).astype(int)
    df["log10_miss"] = np.log10(df.miss_km.clip(lower=1e-3))
    return df


def primary(df: pd.DataFrame, cfg: S.StudyConfig, legacy: bool = False) -> pd.DataFrame:
    """The load-matched arm (classical + learned at g_primary, MPPI rep 0), or the legacy 8 g arm."""
    load = f"{cfg.g_legacy if legacy else cfg.g_primary:g}g"
    learned = df.controller.isin(S.LEARNED) & (df.load == load) & (df.mppi_rep == 0)
    return df[learned | df.controller.str.startswith("classical")]


def spearman(a, b) -> float:
    """Spearman rank correlation (average ranks for ties), NaN-safe, no SciPy needed."""
    d = pd.DataFrame({"a": np.asarray(a, dtype=float), "b": np.asarray(b, dtype=float)}).dropna()
    if len(d) < 3 or d.a.nunique() < 2 or d.b.nunique() < 2:
        return float("nan")
    return float(np.corrcoef(d.a.rank(), d.b.rank())[0, 1])


def _p90(x):
    x = np.asarray(x, dtype=float)
    return float(np.nanpercentile(x, 90)) if np.isfinite(x).any() else float("nan")


METRICS: Dict[str, Callable] = {
    "success5_%": lambda g: 100 * g.success_5km.mean(),
    "success25_%": lambda g: 100 * g.success_25km.mean(),
    "miss_median_km": lambda g: g.miss_km.median(),
    "miss_p90_km": lambda g: _p90(g.miss_km),
    "peak_g_mean": lambda g: g.peak_g.mean(),
    "peak_qdot_MWm2_mean": lambda g: g.peak_qdot_MWm2.mean(),
    "heat_load_MJm2_mean": lambda g: g.heat_load_MJm2.mean(),
    "rcs_fuel_kg_mean": lambda g: g.rcs_fuel_kg.mean(),
    "jet_s_mean": lambda g: g.jet_s.mean(),
    "reversals_mean": lambda g: g.n_reversals.mean(),
    "flight_s_mean": lambda g: g.t_s.mean(),
    "plan_ms_mean": lambda g: g.plan_ms_mean.mean(),
}


def metric_table(df: pd.DataFrame, by=("controller", "atmosphere")) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Per group and seed: every metric; then mean / min / max over seeds and the pooled value."""
    by = list(by)
    rows = []
    for key, g in df.groupby(by + ["seed"], dropna=False):
        rows.append({**dict(zip(by + ["seed"], key)), "n": len(g), **{k: f(g) for k, f in METRICS.items()}})
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
            r[f"{k}"] = g[k].mean()
            r[f"{k}_min"] = g[k].min()
            r[f"{k}_max"] = g[k].max()
            r[f"{k}_pooled"] = METRICS[k](pooled)
        out.append(r)
    return pd.DataFrame(out), per_seed


def cube(df: pd.DataFrame, controller: str, col: str, seeds: Sequence[int], cases: Sequence[int],
         load: Optional[str] = None, rep: int = 0) -> np.ndarray:
    """Values [n_seeds, n_atm, n_cases] (NaN where a flight is missing). Classical: [1, A, C]."""
    d = df[df.controller == controller]
    if controller in S.LEARNED:
        d = d[(d.load == load) & (d.mppi_rep == rep)]
        ss = list(seeds)
    else:
        ss = [-1]
    X = np.full((len(ss), len(ATM), len(cases)), np.nan)
    si = {s: i for i, s in enumerate(ss)}
    ai = {a: i for i, a in enumerate(ATM)}
    ci = {c: i for i, c in enumerate(cases)}
    for r in d[["seed", "atmosphere", "case", col]].itertuples(index=False):
        if r.seed in si and r.atmosphere in ai and r.case in ci:
            X[si[r.seed], ai[r.atmosphere], ci[r.case]] = float(getattr(r, col))
    return X


STATS: Dict[str, Callable] = {
    "mean": lambda x: np.nanmean(x, axis=-1),
    "median": lambda x: np.nanmedian(x, axis=-1),
    "p90": lambda x: np.nanpercentile(x, 90, axis=-1),
}


def hboot(XA: np.ndarray, XB: np.ndarray, stat: str, n_boot: int, seed: int = 0,
          atm: Optional[int] = None) -> Dict[str, float]:
    """
    Paired hierarchical bootstrap of stat(A) - stat(B): resample seeds, then
    cases within each atmosphere; A and B share every index (common ICs and
    MPPI noise). A deterministic arm ([1, A, C]) is repeated over seeds.
    Returns point estimate, 95% CI, and a two-sided bootstrap p-value.
    """
    S_ = max(XA.shape[0], XB.shape[0])
    XA = np.repeat(XA, S_, 0) if XA.shape[0] == 1 else XA
    XB = np.repeat(XB, S_, 0) if XB.shape[0] == 1 else XB
    if atm is not None:
        XA, XB = XA[:, atm:atm + 1], XB[:, atm:atm + 1]
    A, C = XA.shape[1:]
    f = STATS[stat]
    point = float(f(XA.reshape(1, -1))[0] - f(XB.reshape(1, -1))[0])
    rng = np.random.default_rng(seed)
    diffs = []
    for b0 in range(0, n_boot, 1000):
        nb = min(1000, n_boot - b0)
        si = rng.integers(0, S_, (nb, S_))
        ci = rng.integers(0, C, (nb, A, C))
        ix = (si[:, :, None, None], np.arange(A)[None, None, :, None], ci[:, None, :, :])
        diffs.append(f(XA[ix].reshape(nb, -1)) - f(XB[ix].reshape(nb, -1)))
    d = np.concatenate(diffs)
    d = d[np.isfinite(d)]
    lo, hi = np.percentile(d, [2.5, 97.5]) if len(d) else (np.nan, np.nan)
    p = min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean())) if len(d) else np.nan
    return {"diff": point, "ci_lo": float(lo), "ci_hi": float(hi), "p_boot": float(max(p, 1.0 / max(len(d), 1)))}


def mcnemar_p(a: np.ndarray, b: np.ndarray) -> Tuple[int, int, float]:
    """Exact McNemar on paired booleans: (A-only successes, B-only successes, two-sided p)."""
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok].astype(bool), b[ok].astype(bool)
    n10, n01 = int((a & ~b).sum()), int((~a & b).sum())
    return n10, n01, _binom_two_sided(n10, n10 + n01)


def _binom_two_sided(k: int, n: int) -> float:
    if n == 0:
        return 1.0
    k = min(k, n - k)
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n)


def sign_test(x: np.ndarray, y: np.ndarray) -> Tuple[int, int, float]:
    """Paired sign test of x < y (smaller miss wins): (wins of x, non-ties, two-sided p)."""
    ok = np.isfinite(x) & np.isfinite(y)
    d = x[ok] - y[ok]
    w, n = int((d < 0).sum()), int((d != 0).sum())
    return w, n, _binom_two_sided(w, n)


def holm(pvals: Sequence[float]) -> List[float]:
    """Holm-Bonferroni adjusted p-values (same order as given)."""
    p = np.asarray(pvals, dtype=float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (len(p) - rank) * p[i]))
        adj[i] = run
    return adj.tolist()


COMPARISONS = [("v2", "v1"), ("v2", "classical"), ("v2", "v2single"), ("v1", "classical")]


def analyze(cfg: S.StudyConfig, layout: S.Layout, df: Optional[pd.DataFrame] = None) -> dict:
    """All flight tables, statistics, criteria and figures. Returns a summary dict for the report."""
    plt = _plt()
    caps = Captions(layout)
    df = load_flights(layout, cfg, df)
    if not len(df):
        print("[analyze] no flights yet")
        return {}
    cases = sorted(S.load_testset(cfg).case.tolist())
    load = f"{cfg.g_primary:g}g"
    P = primary(df, cfg)
    seeds = sorted(int(s) for s in P[P.controller == "v2"].seed.unique() if s >= 0)
    out: Dict[str, object] = {"seeds": seeds, "n_flights": len(df), "load": load}

    # ---- T01 / T02: the main metric tables ----
    t01, per_seed = metric_table(P)
    _table(layout, caps, t01.round(3), "T01_metrics_primary",
           f"Every metric per controller and atmosphere at the {cfg.g_primary:g} g load limit: mean over seeds, "
           "min and max over seeds, and pooled over seeds (classical: deterministic, one run).",
           "04_flights/flights.csv")
    _table(layout, caps, per_seed.round(3), "T01b_metrics_per_seed",
           "T01 for every seed separately.", "04_flights/flights.csv")
    t02, _ = metric_table(P[P.feasible.astype(bool)])
    _table(layout, caps, t02.round(3), "T02_metrics_feasible_only",
           "T01 restricted to conditions whose target is inside the conservative 4.5 g footprint in that atmosphere.",
           "04_flights/flights.csv, 01_testset/v3_testset.csv")
    out["T01"], out["T02"], out["per_seed"] = t01, t02, per_seed

    # ---- cubes for paired statistics ----
    def cb(ctrl, col, ld=load, rep=0, sub=None):
        return cube(df if sub is None else sub, ctrl, col, seeds, cases, ld, rep)

    present = set(P.controller)
    rows, fam = [], []
    for a_, b_ in COMPARISONS:
        if a_ not in present or b_ not in present or not seeds:
            continue
        for metric, col, stat in (("success5", "success_5km", "mean"), ("miss_median", "miss_km", "median"),
                                  ("miss_p90", "miss_km", "p90")):
            XA, XB = cb(a_, col), cb(b_, col)
            r = hboot(XA, XB, stat, cfg.n_boot, seed=17)
            rows.append({"a": a_, "b": b_, "metric": metric, "atmosphere": "all", **r})
            for i, t in enumerate(ATM):
                rows.append({"a": a_, "b": b_, "metric": metric, "atmosphere": t,
                             **hboot(XA, XB, stat, cfg.n_boot, seed=17 + i + 1, atm=i)})
    t03 = pd.DataFrame(rows)
    if len(t03):
        t03.loc[t03.metric == "success5", ["diff", "ci_lo", "ci_hi"]] *= 100
    out["T03"] = t03

    # ---- per-seed tests ----
    rows = []
    for a_, b_ in COMPARISONS:
        if a_ not in present or b_ not in present:
            continue
        SA, SB = cb(a_, "success_5km"), cb(b_, "success_5km")
        MA, MB = cb(a_, "miss_km"), cb(b_, "miss_km")
        for k, s_ in enumerate(seeds):
            ka, kb = (k if SA.shape[0] > 1 else 0), (k if SB.shape[0] > 1 else 0)
            n10, n01, pm = mcnemar_p(SA[ka].ravel(), SB[kb].ravel())
            w, n, ps = sign_test(MA[ka].ravel(), MB[kb].ravel())
            rows.append({"a": a_, "b": b_, "seed": s_,
                         "success5_diff_pp": 100 * (np.nanmean(SA[ka]) - np.nanmean(SB[kb])),
                         "a_only_success": n10, "b_only_success": n01, "mcnemar_p": pm,
                         "miss_wins_a": w, "miss_nonties": n, "sign_p": ps})
    t04 = pd.DataFrame(rows)
    _table(layout, caps, t04.round(4), "T04_per_seed_tests",
           "Each comparison on each seed alone: success@5km difference (percentage points), exact McNemar test on "
           "the paired successes, and the exact sign test on the paired miss distances.", "04_flights/flights.csv")
    out["T04"] = t04

    # ---- criteria ----
    crit = criteria(cfg, df, P, t03, t04, seeds, cases)
    # Holm over the hypothesis-test criteria (C1, C2, C7): adjust their bootstrap p-values together.
    tested = [i for i, c in enumerate(crit) if c.get("p") is not None and np.isfinite(c["p"])]
    if tested:
        adj = holm([crit[i]["p"] for i in tested])
        for i, pa in zip(tested, adj):
            crit[i]["p_holm"] = pa
            crit[i]["met"] = bool(crit[i]["met"] and pa < 0.05)
    t05 = pd.DataFrame(crit)
    _table(layout, caps, t05, "T05_criteria",
           "Pre-registered criteria (v3_plan.md section 8.2) with what was measured and whether each is met. "
           "C1, C2, C7 also require a Holm-adjusted bootstrap p < 0.05.", "05_analysis/tables/T03, T04")
    if len(t03):
        _table(layout, caps, t03.round(4), "T03_paired_differences",
               "Paired differences A − B with 95% hierarchical-bootstrap intervals (resample seeds, then cases in each "
               f"atmosphere; {cfg.n_boot} resamples). success5 in percentage points, miss in km.",
               "04_flights/flights.csv")
    out["T05"] = t05

    # ---- load matching, planner noise, RCS, per-case ----
    out.update(load_matching(cfg, layout, caps, df, seeds))
    out.update(noise_share(cfg, layout, caps, df, cases))
    out.update(rcs_tables(cfg, layout, caps, P))
    out.update(case_table(cfg, layout, caps, P))
    out["figures"] = flight_figures(cfg, layout, caps, df, P, t01, per_seed, t03, seeds, cases)
    return out


def criteria(cfg, df, P, t03, t04, seeds, cases) -> List[dict]:
    """C1-C7 + Bridge (v3_plan.md section 8.2)."""
    out = []

    def row(id_, q, measured, met, p=None):
        out.append({"id": id_, "question": q, "measured": measured, "met": bool(met), "p": p, "p_holm": None})

    def get(a, b, metric, atm="all"):
        r = t03[(t03.a == a) & (t03.b == b) & (t03.metric == metric) & (t03.atmosphere == atm)]
        return r.iloc[0] if len(r) else None

    def seeds_holding(a, b):
        d = t04[(t04.a == a) & (t04.b == b)]
        return int((d.success5_diff_pp > 0).sum()), len(d)

    for id_, a, b, q in (("C1", "v2", "v1", "Is v2 better than v1 at an equal budget?"),
                         ("C2", "v2", "classical", "Does v2 beat the classical controller?")):
        r = get(a, b, "success5")
        if r is None:
            row(id_, q, "not run", False)
            continue
        k, n = seeds_holding(a, b)
        met = r["diff"] > 0 and r["ci_lo"] > 0 and k >= max(1, n - 1)
        extra = ""
        if id_ == "C2":
            worse = [t for t in ATM if (rr := get(a, b, "success5", t)) is not None and rr["ci_hi"] < 0]
            met = met and not worse
            extra = f"; classical significantly better in: {', '.join(worse) if worse else 'none'}"
        row(id_, q, f"success@5km {r['diff']:+.1f} pp [{r['ci_lo']:+.1f}, {r['ci_hi']:+.1f}], holds in {k}/{n} seeds"
            + extra, met, r["p_boot"])
    g = P.groupby(["controller", "atmosphere"]).peak_g.mean()
    if "v2" in g.index.get_level_values(0) and "classical" in g.index.get_level_values(0):
        gaps = {t: g["v2"].get(t, np.nan) - g["classical"].get(t, np.nan) for t in ATM}
        row("C3", "Is the comparison really load-matched?",
            "v2 − classical mean peak g: " + ", ".join(f"{ATM_LAB[t]} {v:+.2f}" for t, v in gaps.items()),
            all(np.isfinite(v) and v <= 0.5 for v in gaps.values()))
    fe = P[P.feasible.astype(bool) & (P.controller == "v2")]
    if len(fe):
        s5 = fe.groupby("atmosphere").success_5km.mean() * 100
        row("C4", "Does v2 reach the goal (≥ 90% within 5 km, feasible conditions)?",
            ", ".join(f"{ATM_LAB[t]} {s5.get(t, np.nan):.1f}%" for t in ATM), all(s5.get(t, 0) >= 90 for t in ATM))
    v2 = P[P.controller == "v2"].set_index(["seed", "case", "atmosphere"])
    cl = P[P.controller == "classical"].set_index(["case", "atmosphere"])
    if len(v2) and len(cl):
        j = v2.join(cl[["peak_qdot_MWm2", "heat_load_MJm2"]], on=["case", "atmosphere"], rsuffix="_cl")
        rq = float((j.peak_qdot_MWm2 / j.peak_qdot_MWm2_cl).median())
        rl = float((j.heat_load_MJm2 / j.heat_load_MJm2_cl).median())
        row("C5", "Is v2 no worse on heating (median paired ratio ≤ 1.10)?",
            f"peak heat rate ratio {rq:.3f}, heat load ratio {rl:.3f}", rq <= 1.10 and rl <= 1.10)
    ps = P[P.controller == "v2"].groupby(["atmosphere", "seed"]).success_5km.mean().unstack() * 100
    if len(ps) and ps.shape[1] > 1:
        rng_ = (ps.max(1) - ps.min(1))
        row("C6", "Is v2 stable across seeds (success@5km range ≤ 10 pp)?",
            ", ".join(f"{ATM_LAB[t]} {rng_.get(t, np.nan):.1f} pp" for t in ATM), all(rng_.get(t, 99) <= 10 for t in ATM))
    r = get("v2", "v2single", "success5")
    if r is not None:
        row("C7", "Does the ensemble matter (secondary)?",
            f"v2 − single member success@5km {r['diff']:+.1f} pp [{r['ci_lo']:+.1f}, {r['ci_hi']:+.1f}]",
            r["diff"] > 0 and r["ci_lo"] > 0, r["p_boot"])
    leg = df[(df.controller == "v2") & (df.load == f"{cfg.g_legacy:g}g") & (df.seed == 0) & (df.case < 16)]
    if len(leg):
        med = leg.groupby("atmosphere").miss_km.median()
        dev = {t: abs(med.get(t, np.nan) - BRIDGE_REF_MEDIAN_KM[t]) for t in ATM}
        row("Bridge", "Does v3 reproduce v2's first run (seed 0, 8 g, 16 bridge cases)?",
            ", ".join(f"{ATM_LAB[t]} {med.get(t, np.nan):.2f} vs {BRIDGE_REF_MEDIAN_KM[t]:.2f} km" for t in ATM),
            all(np.isfinite(v) and v <= 1.0 for v in dev.values()))
    return out


def load_matching(cfg, layout, caps, df, seeds) -> dict:
    """T06: what the 4.5 g limit costs (seed 0, v1 and v2 at 8 g vs 4.5 g)."""
    rows = []
    for ctrl in ("v1", "v2"):
        for ld in (f"{cfg.g_legacy:g}g", f"{cfg.g_primary:g}g"):
            d = df[(df.controller == ctrl) & (df.load == ld) & (df.mppi_rep == 0) & df.seed.isin(cfg.legacy_seeds)]
            for t, g in d.groupby("atmosphere"):
                rows.append({"controller": ctrl, "load": ld, "atmosphere": t, "n": len(g),
                             **{k: METRICS[k](g) for k in ("success5_%", "miss_median_km", "miss_p90_km", "peak_g_mean",
                                                           "peak_qdot_MWm2_mean", "rcs_fuel_kg_mean")}})
    t06 = pd.DataFrame(rows)
    if len(t06):
        _table(layout, caps, t06.round(3), "T06_load_matching_cost",
               f"The legacy {cfg.g_legacy:g} g limit against the {cfg.g_primary:g} g limit, same seed(s) "
               f"{list(cfg.legacy_seeds)} and noise: what load matching costs.", "04_flights/flights.csv")
    return {"T06": t06}


def noise_share(cfg, layout, caps, df, cases) -> dict:
    """T07: MPPI's own randomness (extra reps, one seed) against seed-to-seed variation."""
    d = df[(df.controller == "v2") & (df.load == f"{cfg.g_primary:g}g")]
    reps = d[d.seed == cfg.noise_seed]
    seeds = d[d.mppi_rep == 0]
    if reps.mppi_rep.nunique() < 2 or seeds.seed.nunique() < 2:
        return {}
    var_rep = reps.groupby(["case", "atmosphere"]).log10_miss.var().mean()
    var_seed = seeds.groupby(["case", "atmosphere"]).log10_miss.var().mean()
    s_rep = reps.groupby(["atmosphere", "mppi_rep"]).success_5km.mean().unstack() * 100
    s_seed = seeds.groupby(["atmosphere", "seed"]).success_5km.mean().unstack() * 100
    t07 = pd.DataFrame({"atmosphere": ATM,
                        "success5_range_over_mppi_reps_pp": [(s_rep.loc[t].max() - s_rep.loc[t].min()) if t in s_rep.index else np.nan for t in ATM],
                        "success5_range_over_seeds_pp": [(s_seed.loc[t].max() - s_seed.loc[t].min()) if t in s_seed.index else np.nan for t in ATM]})
    t07.loc[len(t07)] = {"atmosphere": "per-condition var of log10 miss: MPPI reps / seeds",
                         "success5_range_over_mppi_reps_pp": var_rep, "success5_range_over_seeds_pp": var_seed}
    _table(layout, caps, t07.round(4), "T07_planner_noise_share",
           f"Planner noise vs training seeds (v2, {cfg.g_primary:g} g): success@5km range over MPPI repeats of seed "
           f"{cfg.noise_seed} against the range over seeds, and the mean per-condition variance of log10(miss) "
           f"(ratio = share of seed-to-seed variation that MPPI noise alone produces: {var_rep / var_seed:.2f}).",
           "04_flights/flights.csv")
    return {"T07": t07, "noise_share": float(var_rep / var_seed)}


def rcs_tables(cfg, layout, caps, P) -> dict:
    """T08: RCS usage per controller and atmosphere, and its rank correlation with miss."""
    rows = []
    for (c, t), g in P.groupby(["controller", "atmosphere"]):
        rows.append({"controller": c, "atmosphere": t, "fuel_kg": g.rcs_fuel_kg.mean(), "jet_s": g.jet_s.mean(),
                     "fuel_above60km": g.fuel_kg_above60km.mean(), "fuel_60_25km": g.fuel_kg_60_25km.mean(),
                     "fuel_below25km": g.fuel_kg_below25km.mean(), "bank_track_rms_deg": g.bank_track_rms_deg.mean(),
                     "roll_rate_max_dps": g.roll_rate_max_dps.mean(), "roll_sat_frac": g.roll_sat_frac.mean(),
                     "rho_fuel_vs_miss": spearman(g.rcs_fuel_kg, g.miss_km),
                     "rho_track_vs_miss": spearman(g.bank_track_rms_deg, g.miss_km)})
    t08 = pd.DataFrame(rows)
    _table(layout, caps, t08.round(3), "T08_rcs_usage",
           "RCS usage (mean per flight) per controller and atmosphere at the primary load limit: fuel, jet-on seconds, "
           "fuel by phase, bank-tracking RMS (commanded − actual), roll rate, fraction of time within 10% of the 20°/s "
           "limit, and Spearman correlations with the miss.", "04_flights/flights.csv")
    return {"T08": t08}


def case_table(cfg, layout, caps, P) -> dict:
    """T11: median miss per case, atmosphere and controller (over seeds) - the hard cases."""
    t = P.groupby(["case", "atmosphere", "controller"]).miss_km.median().unstack("controller").reset_index()
    s5 = P.groupby(["case", "atmosphere", "controller"]).success_5km.mean().unstack("controller")
    s5.columns = [f"success5_{c}" for c in s5.columns]
    t = t.merge(s5.reset_index(), on=["case", "atmosphere"])
    _table(layout, caps, t.round(3), "T11_per_case",
           "Per case and atmosphere: median miss over seeds [km] and success@5km fraction over seeds, per controller.",
           "04_flights/flights.csv")
    return {"T11": t}


# ============================================================================
# Figures
# ============================================================================
def flight_figures(cfg, layout, caps, df, P, t01, per_seed, t03, seeds, cases) -> Dict[str, Path]:
    plt = _plt()
    figs: Dict[str, Path] = {}
    ctrls = [c for c in CTRL if c in set(P.controller)]
    x = np.arange(len(ATM))
    w = 0.8 / max(1, len(ctrls))

    # F10 success@5km by atmosphere, seed range + dots
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.4))
    for j, c in enumerate(ctrls):
        r = t01[t01.controller == c].set_index("atmosphere").reindex(ATM)
        xs = x + (j - (len(ctrls) - 1) / 2) * w
        for a_, col in ((ax[0], "success5_%"), (ax[1], "success25_%")):
            a_.bar(xs, r[col], w * 0.9, color=COL[c], label=LAB[c])
            if c in S.LEARNED:
                a_.errorbar(xs, r[col], yerr=[r[col] - r[f"{col}_min"], r[f"{col}_max"] - r[col]], fmt="none",
                            ecolor="#0b0b0b", lw=0.8, capsize=2)
                ps = per_seed[per_seed.controller == c]
                for i, t in enumerate(ATM):
                    v = ps[ps.atmosphere == t][col].values
                    a_.scatter(np.full(len(v), xs[i]), v, s=5, color="#0b0b0b", zorder=3)
    for a_, tt in ((ax[0], "Success within 5 km (the goal)"), (ax[1], "Success within 25 km (v1/v2 legacy)")):
        a_.set_xticks(x, [ATM_LAB[t] for t in ATM])
        a_.set(ylim=(0, 105), ylabel="% of flights", title=tt)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, ncol=len(ctrls), fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, 1.07))
    fig.tight_layout()
    figs["F10"] = _fig(layout, caps, plt, fig, "F10_success_by_atmosphere",
                       f"Success rate per atmosphere at the {cfg.g_primary:g} g limit: bars = mean over seeds, whiskers = "
                       "seed range, dots = individual seeds (64 cases each). Classical is deterministic.",
                       "05_analysis/tables/T01_metrics_primary.csv")

    # F11 miss ECDF per atmosphere
    fig, ax = plt.subplots(1, 4, figsize=(13, 3.0), sharey=True)
    for i, t in enumerate(ATM):
        for c in ctrls:
            v = np.sort(P[(P.controller == c) & (P.atmosphere == t)].miss_km.values)
            if len(v):
                ax[i].plot(v, np.arange(1, len(v) + 1) / len(v), color=COL[c], label=LAB[c], lw=1.4)
        ax[i].axvline(5, color="#0b0b0b", lw=0.7, ls="--")
        ax[i].set(xscale="log", xlabel="miss [km]", title=ATM_LAB[t])
    ax[0].set_ylabel("fraction of flights ≤ miss")
    ax[0].legend(fontsize=6.5, loc="upper left")
    fig.tight_layout()
    figs["F11"] = _fig(layout, caps, plt, fig, "F11_miss_ecdf_by_atmosphere",
                       "Cumulative distribution of the miss distance, pooled over seeds (dashed line: the 5 km goal). "
                       "Higher and further left is better.", "04_flights/flights.csv")

    # F12 paired differences (forest plot)
    if len(t03):
        d = t03[t03.metric == "success5"]
        pairs = list(dict.fromkeys(zip(d.a, d.b)))
        fig, ax = plt.subplots(figsize=(7.5, 0.8 + 0.26 * len(pairs) * 5.6))
        y = 0
        ticks, labels = [], []
        for a_, b_ in pairs:
            for t in ["all"] + ATM:
                r = d[(d.a == a_) & (d.b == b_) & (d.atmosphere == t)]
                if not len(r):
                    continue
                r = r.iloc[0]
                col = COL[a_] if r.ci_lo > 0 else (COL[b_] if r.ci_hi < 0 else "#a3a29c")
                ax.plot([r.ci_lo, r.ci_hi], [y, y], color=col, lw=2.2 if t == "all" else 1.2)
                ax.scatter([r["diff"]], [y], color=col, s=22 if t == "all" else 12, zorder=3)
                ticks.append(y)
                labels.append(f"{LAB[a_]} − {LAB[b_]}" if t == "all" else f"   {ATM_LAB[t]}")
                y -= 1
            y -= 0.6
        ax.axvline(0, color="#0b0b0b", lw=0.7)
        ax.set_yticks(ticks, labels, fontsize=7)
        ax.set(xlabel="difference in success@5km [percentage points] (95% hierarchical bootstrap)",
               title="Paired differences: right of 0 favors the first controller")
        fig.tight_layout()
        figs["F12"] = _fig(layout, caps, plt, fig, "F12_paired_differences",
                           "Paired success@5km differences with 95% intervals (resampling seeds, then cases). Colored when "
                           "the interval excludes 0.", "05_analysis/tables/T03_paired_differences.csv")

    # F13 load and heating
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.1))
    for k, (col, yl) in enumerate((("peak_g", "peak load [g]"), ("peak_qdot_MWm2", "peak heat rate [MW/m²]"),
                                   ("heat_load_MJm2", "heat load [MJ/m²]"))):
        for j, c in enumerate(ctrls):
            data = [P[(P.controller == c) & (P.atmosphere == t)][col].dropna().values for t in ATM]
            pos = x + (j - (len(ctrls) - 1) / 2) * w
            bp = ax[k].boxplot([dd if len(dd) else [np.nan] for dd in data], positions=pos, widths=w * 0.8,
                               patch_artist=True, showfliers=False, medianprops={"color": "#0b0b0b"})
            for b in bp["boxes"]:
                b.set(facecolor=COL[c], alpha=0.85, linewidth=0.5)
        ax[k].set_xticks(x, [ATM_LAB[t] for t in ATM], fontsize=7)
        ax[k].set(ylabel=yl)
    ax[0].axhline(cfg.g_primary, color="#0b0b0b", lw=0.7, ls="--")
    ax[0].set_title(f"Peak load (dashed: {cfg.g_primary:g} g limit)")
    ax[1].set_title("Peak heat rate")
    ax[2].set_title("Total heat load")
    fig.tight_layout()
    figs["F13"] = _fig(layout, caps, plt, fig, "F13_load_and_heating",
                       "Peak g, peak stagnation heat rate, and heat load per controller and atmosphere (boxes: quartiles; "
                       "whiskers: 1.5 IQR; colors as in F10).", "04_flights/flights.csv")

    # F14 load matching cost
    lm = df[df.controller.isin(["v1", "v2"]) & df.seed.isin(cfg.legacy_seeds) & (df.mppi_rep == 0)]
    if lm.load.nunique() > 1:
        fig, ax = plt.subplots(1, 2, figsize=(10, 3.0))
        for k, (col, yl) in enumerate((("success_5km", "success@5km [%]"), ("peak_g", "mean peak g"))):
            for j, (c, ld) in enumerate([(c, ld) for c in ("v1", "v2") for ld in (f"{cfg.g_legacy:g}g", f"{cfg.g_primary:g}g")]):
                v = lm[(lm.controller == c) & (lm.load == ld)].groupby("atmosphere")[col].mean().reindex(ATM)
                v = v * (100 if col == "success_5km" else 1)
                ax[k].bar(x + (j - 1.5) * 0.2, v.values, 0.18, color=COL[c], alpha=1.0 if ld.startswith(f"{cfg.g_primary:g}") else 0.45,
                          label=f"{c} {ld}", hatch=None if ld.startswith(f"{cfg.g_primary:g}") else "///")
            ax[k].set_xticks(x, [ATM_LAB[t] for t in ATM], fontsize=7)
            ax[k].set(ylabel=yl)
        ax[0].legend(fontsize=7, ncol=2)
        ax[0].set_title(f"What the {cfg.g_primary:g} g limit costs (seed {cfg.legacy_seeds[0]})")
        ax[1].set_title("... and what it buys")
        fig.tight_layout()
        figs["F14"] = _fig(layout, caps, plt, fig, "F14_load_matching_cost",
                           f"Same models and noise, legacy {cfg.g_legacy:g} g limit (hatched) against {cfg.g_primary:g} g.",
                           "05_analysis/tables/T06_load_matching_cost.csv")

    # F15 miss by case heatmap
    hm = P.groupby(["case", "controller", "atmosphere"]).miss_km.median().unstack(["controller", "atmosphere"])
    cols = [(c, t) for c in ctrls for t in ATM if (c, t) in hm.columns]
    if cols:
        hm = hm.reindex(columns=pd.MultiIndex.from_tuples(cols))
        fig, ax = plt.subplots(figsize=(min(13, 1.0 + 0.55 * len(cols)), 1.6 + 0.12 * len(hm.index)))
        im = ax.imshow(np.log10(hm.values.clip(0.05, 2000)), aspect="auto", cmap="RdYlGn_r", vmin=-1, vmax=2)
        ax.set_xticks(range(len(cols)), [f"{c}\n{t[:8]}" for c, t in cols], fontsize=6, rotation=90)
        ax.set_yticks(range(len(hm.index)), hm.index, fontsize=5)
        n_bridge = int((np.asarray(hm.index) < 16).sum())
        if 0 < n_bridge < len(hm.index):
            ax.axhline(n_bridge - 0.5, color="#0b0b0b", lw=0.8)
        for k in range(1, len(ctrls)):
            ax.axvline(k * len(ATM) - 0.5, color="white", lw=1.5)
        cb_ = fig.colorbar(im, ax=ax, shrink=0.6)
        cb_.set_label("log10 median miss [km] (green ≤ 1 km, red ≥ 100 km)")
        ax.set(ylabel="case (above the line: 16 bridge cases)", title="Median miss over seeds, per case")
        fig.tight_layout()
        figs["F15"] = _fig(layout, caps, plt, fig, "F15_miss_by_case",
                           "Median miss over seeds for every case, controller and atmosphere. Rows above the black line are "
                           "the 16 bridge cases.", "05_analysis/tables/T11_per_case.csv")

    # F16 success vs entry conditions
    sv = P.groupby(["controller", "case"]).agg(s5=("success_5km", "mean"), V=("velocity_mps", "first"),
                                               gam=("gamma_deg", "first"), h=("altitude_m", "first")).reset_index()
    fig, ax = plt.subplots(1, len(ctrls), figsize=(3.2 * len(ctrls), 2.9), sharey=True, squeeze=False)
    for k, c in enumerate(ctrls):
        d = sv[sv.controller == c]
        sc = ax[0, k].scatter(d.V, d.gam, c=100 * d.s5, cmap="RdYlGn", vmin=0, vmax=100, s=18, edgecolor="#52514e",
                              linewidth=0.3)
        ax[0, k].set(xlabel="entry speed [m/s]", title=LAB[c])
    ax[0, 0].set_ylabel("flight-path angle [deg]")
    fig.colorbar(sc, ax=ax[0, -1], label="success@5km over seeds and atmospheres [%]")
    fig.tight_layout()
    figs["F16"] = _fig(layout, caps, plt, fig, "F16_success_vs_entry",
                       "Which entries are hard: success@5km per case (all seeds and atmospheres) against entry speed and "
                       "flight-path angle.", "04_flights/flights.csv")

    # F17 seed stability
    learned = [c for c in ctrls if c in S.LEARNED]
    if learned and len(seeds) > 1:
        fig, ax = plt.subplots(1, 4, figsize=(13, 2.8), sharey=True)
        for i, t in enumerate(ATM):
            for c in learned:
                d = per_seed[(per_seed.controller == c) & (per_seed.atmosphere == t)].sort_values("seed")
                ax[i].plot(d.seed, d["success5_%"], marker="o", ms=3, color=COL[c], label=LAB[c])
            cl_ = per_seed[(per_seed.controller == "classical") & (per_seed.atmosphere == t)]["success5_%"]
            if len(cl_):
                ax[i].axhline(float(cl_.iloc[0]), color=COL["classical"], ls="--", lw=1, label="classical")
            ax[i].set(xlabel="seed", title=ATM_LAB[t], xticks=seeds)
        ax[0].set_ylabel("success@5km [%]")
        ax[0].legend(fontsize=6.5)
        fig.tight_layout()
        figs["F17"] = _fig(layout, caps, plt, fig, "F17_seed_stability",
                           "success@5km of every training seed (dashed: the deterministic classical controller).",
                           "05_analysis/tables/T01b_metrics_per_seed.csv")

    # F20-F22 RCS
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.2))
    for c in ctrls:
        d = P[P.controller == c]
        ax[0].scatter(d.rcs_fuel_kg, d.miss_km, s=4, alpha=0.5, color=COL[c], label=LAB[c])
        ax[1].scatter(d.bank_track_rms_deg, d.miss_km, s=4, alpha=0.5, color=COL[c])
    for a_, xl in ((ax[0], "RCS fuel per flight [kg]"), (ax[1], "bank-tracking RMS error [deg]")):
        a_.set(yscale="log", xlabel=xl, ylabel="miss [km]")
        a_.axhline(5, color="#0b0b0b", lw=0.7, ls="--")
    ax[0].legend(fontsize=6.5, markerscale=3)
    ax[0].set_title("Fuel vs miss")
    ax[1].set_title("Tracking error vs miss")
    ph = P.groupby("controller")[["fuel_kg_above60km", "fuel_kg_60_25km", "fuel_kg_below25km"]].mean().reindex(ctrls)
    bottom = np.zeros(len(ctrls))
    for k, (col, lab) in enumerate((("fuel_kg_above60km", "above 60 km"), ("fuel_kg_60_25km", "60–25 km"),
                                    ("fuel_kg_below25km", "below 25 km"))):
        ax[2].bar(range(len(ctrls)), ph[col].values, bottom=bottom, color=["#2a78d6", "#8fb8ea", "#d4e4f7"][k], label=lab)
        bottom += np.nan_to_num(ph[col].values)
    ax[2].set_xticks(range(len(ctrls)), [LAB[c] for c in ctrls], fontsize=6.5, rotation=15)
    ax[2].set(ylabel="mean RCS fuel per flight [kg]", title="Fuel by flight phase")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    figs["F20"] = _fig(layout, caps, plt, fig, "F20_rcs_fuel_tracking_phase",
                       "RCS: fuel per flight against miss, bank-tracking error against miss, and mean fuel by flight phase "
                       f"(all at the {cfg.g_primary:g} g limit, pooled over seeds and atmospheres).",
                       "05_analysis/tables/T08_rcs_usage.csv")

    # F21 example flights from telemetry
    ex = example_flights(cfg, layout, P, seeds)
    if ex:
        fig, ax = plt.subplots(2, len(ex), figsize=(4.2 * len(ex), 5.4), squeeze=False)
        for k, (title, recs) in enumerate(ex):
            for c, tel in recs:
                ax[0, k].plot(tel["lambda_deg"], tel["phi_deg"], color=COL[c], lw=1.2, label=LAB[c])
                ax[1, k].plot(tel["t_s"], tel["sigma_deg"], color=COL[c], lw=1.0)
                ax[1, k].plot(tel["t_s"], tel["sigma_cmd_deg"], color=COL[c], lw=0.6, ls=":")
            ax[0, k].scatter([15], [15], marker="*", s=90, color="#0b0b0b", zorder=3)
            ax[0, k].set(title=title, xlabel="longitude [deg]", ylabel="latitude [deg]")
            ax[1, k].set(xlabel="time [s]", ylabel="bank [deg] (dotted: command)")
        ax[0, 0].legend(fontsize=6.5)
        fig.tight_layout()
        figs["F21"] = _fig(layout, caps, plt, fig, "F21_example_flights",
                           "Example flights from the 1 Hz telemetry (seed 0): ground tracks to the target (star), and "
                           "actual bank (solid) against the command (dotted).", "04_flights/telemetry/")
    return figs


def example_flights(cfg, layout, P, seeds, n: int = 3):
    """Telemetry of a few contrasting conditions: v2's best, its median, and its worst (seed 0)."""
    if not seeds:
        return []
    v2 = P[(P.controller == "v2") & (P.seed == seeds[0])].sort_values("miss_km")
    if not len(v2):
        return []
    picks = [v2.iloc[0], v2.iloc[len(v2) // 2], v2.iloc[-1]][:n]
    out = []
    for r in picks:
        recs = []
        for c in ("classical", "v1", "v2"):
            d = P[(P.controller == c) & (P.case == r.case) & (P.atmosphere == r.atmosphere) &
                  ((P.seed == seeds[0]) | (P.seed == -1))]
            if len(d) and layout.telemetry_path(d.iloc[0].flight_id).exists():
                recs.append((c, S.load_telemetry(layout, d.iloc[0].flight_id)))
        out.append((f"case {int(r.case)}, {ATM_LAB[r.atmosphere]} (v2 miss {r.miss_km:.1f} km)", recs))
    return out


# ============================================================================
# Report
# ============================================================================
def build_report(cfg: S.StudyConfig, layout: S.Layout, A: dict, OL: Optional[dict] = None,
                 out_pdf: Optional[str] = None) -> str:
    """06_report/SpaceJEPA_v3_report.pdf from the tables and figures (all numbers are read, none typed in)."""
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.units import inch
    from reportlab.platypus import PageBreak, Paragraph, SimpleDocTemplate
    from jepa.make_report_v2 import MUT, RULE, TINT_YELLOW, Doc, _fonts

    _fonts()
    d = Doc()
    out_pdf = out_pdf or str(layout.report / "SpaceJEPA_v3_report.pdf")
    t01, t05 = A["T01"], A["T05"]
    seeds = A["seeds"]

    def m(c, t, k, fmt="{:.1f}"):
        r = t01[(t01.controller == c) & (t01.atmosphere == t)]
        return fmt.format(r.iloc[0][k]) if len(r) else "–"

    d.add(Paragraph("SpaceJEPA v3: v1 vs v2 vs classical", d.s["title"]),
          Paragraph(f"{len(seeds)} training seeds · 64 entry cases × 4 atmospheres · load-matched at "
                    f"{cfg.g_primary:g} g · goal: land within {cfg.success_km:g} km · {A['n_flights']} flights",
                    d.s["subtitle"]))
    lines = []
    for t in ATM:
        lines.append(f"<b>{ATM_LAB[t]}</b>: success@5km classical {m('classical', t, 'success5_%')}%, "
                     f"v1 {m('v1', t, 'success5_%')}%, v2 {m('v2', t, 'success5_%')}% "
                     f"(median miss {m('classical', t, 'miss_median_km', '{:.2f}')} / {m('v1', t, 'miss_median_km', '{:.2f}')} / "
                     f"{m('v2', t, 'miss_median_km', '{:.2f}')} km)")
    met = t05[t05.id.str.startswith("C")]
    lines.append(f"<b>Criteria met: {int(met.met.sum())} of {len(met)}</b> pre-registered (C1–C7), "
                 f"Bridge check: {'met' if t05[t05.id == 'Bridge'].met.any() else 'not met or not run'}.")
    d.callout("The numbers (learned controllers: mean over seeds)", lines, TINT_YELLOW, colors.HexColor("#eda100"))
    d.h2("Pre-registered criteria (v3_plan.md §8.2, fixed before any v3 flight)")
    rows = [["#", "Question", "Measured", "p (Holm)", "Met"]]
    for r in t05.itertuples():
        ph = "" if r.p_holm is None or (isinstance(r.p_holm, float) and not np.isfinite(r.p_holm)) else f"{r.p_holm:.3g}"
        rows.append([r.id, r.question, r.measured, ph, "yes" if r.met else "<b>no</b>"])
    d.table(rows, [0.45, 1.9, 3.3, 0.6, 0.45])

    figs = A.get("figures", {})
    d.add(PageBreak())
    d.h1("1. Test set")
    d.p("Sixteen bridge cases are the frozen v1/v2 entries. Forty-eight new entries come from a scrambled Sobol "
        "sequence over entry altitude, speed, flight-path angle, approach heading, and the target's position inside "
        "the reachable region. Every new target is reachable in all four atmospheres, within 4.5 g, with at least "
        "50 km to spare, according to constant-bank sweeps of the simulator's own equations of motion.")
    for name in ("F01_testset_coverage", "F02_footprint_margins"):
        p = layout.testset / f"{name}.png"
        if p.exists():
            d.fig(str(p), Captions(layout).d.get(name, {}).get("caption", name), "real")
    d.h1("2. Closed loop: success and miss")
    for k in ("F10", "F11", "F12"):
        if k in figs:
            d.fig(str(figs[k]), Captions(layout).d.get(Path(figs[k]).stem, {}).get("caption", ""), "real",
                  width=(5.4 if k == "F12" else 7.0) * inch)
    rows = [["Controller", "Atmosphere", "success@5km % (seed range)", "median km", "p90 km", "peak g", "fuel kg"]]
    for _, r in t01.iterrows():
        rng = (f" ({r['success5_%_min']:.0f}–{r['success5_%_max']:.0f})"
               if r.controller in S.LEARNED and r.n_seeds > 1 else "")
        rows.append([LAB.get(r.controller, r.controller), ATM_LAB[r.atmosphere], f"{r['success5_%']:.1f}{rng}",
                     f"{r.miss_median_km:.2f}", f"{r.miss_p90_km:.2f}", f"{r.peak_g_mean:.2f}",
                     f"{r.rcs_fuel_kg_mean:.1f}"])
    d.table(rows, [1.4, 1.2, 1.3, 0.8, 0.8, 0.7, 0.8])
    d.h1("3. Load, heating, and the cost of load matching")
    for k in ("F13", "F14"):
        if k in figs:
            d.fig(str(figs[k]), Captions(layout).d.get(Path(figs[k]).stem, {}).get("caption", ""), "real")
    d.h1("4. Which cases are hard, and seed stability")
    for k in ("F15", "F16", "F17"):
        if k in figs:
            d.fig(str(figs[k]), Captions(layout).d.get(Path(figs[k]).stem, {}).get("caption", ""), "real")
    if "noise_share" in A:
        d.p(f"<b>Planner noise.</b> Re-flying seed {cfg.noise_seed} with {len(cfg.noise_reps)} more MPPI noise streams "
            f"reproduces {100 * A['noise_share']:.0f}% of the seed-to-seed variance of log miss (T07): that share of "
            "the differences between seeds is the planner's own randomness, not training.")
    d.h1("5. RCS usage")
    for k in ("F20", "F21"):
        if k in figs:
            d.fig(str(figs[k]), Captions(layout).d.get(Path(figs[k]).stem, {}).get("caption", ""), "real")
    if OL:
        d.h1("6. Open loop, per seed")
        for k in ("F30", "F31"):
            if k in OL:
                d.fig(str(OL[k]), Captions(layout).d.get(Path(OL[k]).stem, {}).get("caption", ""), "real",
                      width=(4.2 if k == "F31" else 7.0) * inch)
    d.h1("Methods and definitions")
    d.bullets([
        "<b>success@5km</b>: the flight reaches drogue deploy (7.6 km altitude) and its great-circle distance to "
        "15°N 15°E is below 5 km. Any other ending is a failure; its miss is the distance at termination.",
        "<b>median / p90 miss</b>: 50th / 90th percentile (linear interpolation) over the group, failures included.",
        "<b>peak g</b>: maximum aerodynamic load factor over every physics substep. <b>peak heat rate</b>: maximum "
        "stagnation-point heat rate (simulator heating model). <b>heat load</b>: trapezoid integral of the 1 Hz heat rate.",
        "<b>RCS fuel</b>: the simulator's propellant count from jet on-time (linear: kg = jet-seconds × thrust ÷ "
        "(Isp·g0)); jet-seconds are recovered exactly from it.",
        f"<b>Statistics</b>: paired by condition (same IC, same MPPI noise stream for every learned controller); "
        f"hierarchical bootstrap ({cfg.n_boot} resamples: seeds, then cases within each atmosphere); per-seed exact "
        "McNemar and sign tests; Holm correction over C1, C2, C7.",
        f"<b>Planner</b>: MPPI, {cfg.mppi.get('n_samples')} samples, 6 knots, 30° noise, 4 iterations at entry then 1, "
        f"horizon {S.horizon_steps(cfg)} steps of 5 s from the test set, g limit {cfg.g_primary:g} g (legacy "
        f"{cfg.g_legacy:g} g); v2 adds the ensemble disagreement at {cfg.w_unc_km:g} km per km.",
        "<b>Feasibility</b>: the conservative constant-bank footprint within 4.5 g; T02 repeats T01 on feasible "
        "conditions only.",
    ])
    d.p(f"All tables and figures, with their data sources, are in the study folder's 05_analysis/ "
        f"(captions.json). Git commit: {json.loads((layout.study / 'study_config.json').read_text()).get('git_commit', '?')[:10]}.")

    def page(canvas, doc):
        if doc.page == 1:
            return
        canvas.saveState()
        canvas.setFont("DV", 7.5)
        canvas.setFillColor(MUT)
        canvas.drawString(0.75 * inch, 0.5 * inch, "SpaceJEPA v3: full study")
        canvas.drawRightString(letter[0] - 0.75 * inch, 0.5 * inch, str(doc.page))
        canvas.setStrokeColor(RULE)
        canvas.setLineWidth(0.5)
        canvas.line(0.75 * inch, letter[1] - 0.55 * inch, letter[0] - 0.75 * inch, letter[1] - 0.55 * inch)
        canvas.restoreState()

    layout.report.mkdir(parents=True, exist_ok=True)
    tmp = out_pdf + ".tmp.pdf"
    doc = SimpleDocTemplate(tmp, pagesize=letter, leftMargin=0.75 * inch, rightMargin=0.75 * inch,
                            topMargin=0.75 * inch, bottomMargin=0.75 * inch, title="SpaceJEPA v3 report")
    doc.build(d.story, onFirstPage=page, onLaterPages=page)
    os.replace(tmp, out_pdf)
    return out_pdf
