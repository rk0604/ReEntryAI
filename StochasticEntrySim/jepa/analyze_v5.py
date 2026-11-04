"""
jepa/analyze_v5.py

Everything the v5 study reports (jepa/v5_plan.md sections 5, 6, 8, 10), written to 05_analysis/ as numbered
tables (T01-T12, CSV) and figures (F01-F10, PNG and PDF at 200 dpi) with captions.json.

    load_flights_v5   v3's flights (A256: classical, v1, v2, v2single, ...), v4's (R0, v2's load curve) and v5's
                      in one table (arm, load, suite, seed, case, atmosphere / target, metrics); v3 and v4 read-only
    cube, comp        A256 paired comparisons with analyze_v4.Boot (hierarchical bootstrap: seeds, then cases within
                      atmospheres; cached in 05_analysis/cache/ by comparison and input hash)
    orb_paired        ORB-1 paired comparisons: each case one unit, cases resampled within training-seed strata
    criteria          G1-G7, J-each / J-equiv / J-mech, R-each / R-imit / R-gap, Q1-Q3 (T05)
    analyze           the tables, the figures, the criteria; returns a summary for the report

Definitions (section 5): success@5km / @2km = outcome drogue and miss < 5 / 2 km; success@10km = drogue and
miss <= 10 km; log miss = ln(max(miss, 0.01 km)). Deterministic arms (classical, fnpeg, predguid, their
exploratory variants, legacy RL) fly once per condition and are broadcast over seeds; MPPI arms' "seed s" is
model seed s and noise index s; RL arms' seeds are training seeds 0-2 (comparisons with them pair seeds 0-2).
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
import zlib
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

import numpy as np
import pandas as pd

from jepa import analyze_v3 as A3
from jepa import analyze_v4 as A4
from jepa import benchmark_v5 as B
from jepa import models_v5 as M5
from jepa import study_v3 as S
from jepa import study_v4 as S4
from jepa import study_v5 as S5
from jepa import testset_v3 as TS

ATM = list(TS.ATMOSPHERES)
ATM_LAB = {**A3.ATM_LAB, "all": "pooled"}
G = "4.5g"
OWN = "own"
CL, FN, PG = ("classical", OWN), ("fnpeg", OWN), ("predguid", OWN)
MF, V2S, V2Sm = ("mppi_filtered", G), ("v2single", G), ("v2single", "fnpeg_matched")
R0, V1, V2, MN, MT = ("R0", G), ("v1", G), ("v2", G), ("mppi_nominal", G), ("mppi_truth", G)
FN180, PG180, FNO = ("fnpeg_180", OWN), ("predguid_180", OWN), ("fnpeg_oracle", OWN)
E3 = [(a, G) for a in M5.ARMS_E3]
RL1 = [(a, OWN) for a in S5.RL_TIER1]
RL2 = [(a, OWN) for a in S5.RL_TIER2]
NON_RESULT_STAGES = ("gates", "tuning", "timing", "E4_curve")
NON_TABLE_PREFIX = ("gate_", "tune_", "E0_", "diag_case10", "classical_E0a", "E0c")
LOAD_CURVE = [("v2", "3.5g"), ("v2", "4g"), ("v2", G), ("v2", "8g"), ("v2", "matched")]

LAB = {"classical": "classical (v3)", "fnpeg": "FNPEG (Lu 2014)", "predguid": "PredGuid (Bairstow 2006)",
       "mppi_filtered": "MPPI + filtered physics", "v2single": "v2 single member", "v2": "v2 ensemble",
       "R0": "R0 (v2 recipe, 1 member)", "v1": "v1", "mppi_nominal": "MPPI + nominal physics",
       "mppi_truth": "MPPI + true physics", "fnpeg_180": "FNPEG |σ| ≤ 180°", "predguid_180": "PredGuid |σ| ≤ 180°",
       "fnpeg_oracle": "FNPEG, true parameters", "mlp_direct": "MLP direct", "gru_direct": "GRU direct",
       "gru_pi": "GRU + PI prober", "r0_nojepa": "R0 without JEPA", "dagger_fnpeg": "DAgger (FNPEG expert)",
       "bcppo_fnpeg": "BC → PPO (FNPEG)", "ppo_hist": "PPO", "sac_hist": "SAC", "ppo_matched": "PPO, full sim"}
COL = {"classical": "#6250d6", "fnpeg": "#1baf7a", "predguid": "#2a78d6", "mppi_filtered": "#8fb8ea",
       "v2single": "#eda100", "v2": "#eb6834", "R0": "#c2410c", "v1": "#a3a29c", "mppi_truth": "#0b0b0b",
       "mlp_direct": "#b5179e", "gru_direct": "#7209b7", "gru_pi": "#3a0ca3", "r0_nojepa": "#f77f00",
       "dagger_fnpeg": "#2b9348", "bcppo_fnpeg": "#55a630", "ppo_hist": "#d00000", "sac_hist": "#9d0208",
       "ppo_matched": "#6a040f", "fnpeg_180": "#80ed99", "predguid_180": "#90e0ef", "fnpeg_oracle": "#38a3a5"}


def key_of(arm: str, load: str) -> str:
    return arm if load in (G, OWN) else f"{arm}@{load}"


def label(arm: str, load: str = G) -> str:
    base = LAB.get(arm, arm)
    return base if load in (G, OWN) else f"{base} @{load}"


def color(arm: str) -> str:
    return COL.get(arm.split("@")[0], COL.get(arm, "#52514e"))


# ============================================================================
# Writing (every figure as PNG and PDF at 200 dpi, with its caption)
# ============================================================================
def _fig(L, caps, plt, fig, name, caption, data) -> Path:
    png = L.figures / f"{name}.png"
    L.save_fig(fig, png, dpi=200)
    L.save_fig(fig, L.figures / f"{name}.pdf", dpi=200)
    plt.close(fig)
    caps.add(name, caption, data)
    return png


def _table(L, caps, df: pd.DataFrame, name, caption, data) -> Path:
    return A4._table(L, caps, df, name, caption, data)


def _b(v) -> bool:
    return str(v).lower() in ("true", "1", "1.0")


# ============================================================================
# Flights
# ============================================================================
def _read(layout) -> pd.DataFrame:
    if layout is None:
        return pd.DataFrame()
    f = layout.flights / "flights.csv"
    if f.exists():
        return pd.read_csv(f, low_memory=False)
    if (layout.flights / "results").exists():
        return A4._read_results(layout.root)
    return pd.DataFrame()


def load_flights_v5(cfg: S5.StudyConfigV5, L: S5.LayoutV5, L3: S.Layout, L4=None) -> pd.DataFrame:
    """v3's, v4's and v5's flights in one table (v3 and v4 read-only)."""
    parts = []
    d3 = _read(L3)
    if len(d3):
        d3["source"] = "v3"
        parts.append(d3)
    d4 = _read(L4)
    if len(d4):
        d4["source"] = "v4"
        parts.append(d4)
    d5 = S5.collect_flights(L, write=True)
    if len(d5):
        d5["source"] = "v5"
        parts.append(d5)
    if not parts:
        return pd.DataFrame(columns=["arm", "load", "suite", "seed", "case", "atmosphere", "miss_km", "outcome"])
    df = pd.concat(parts, ignore_index=True, sort=False)
    if "arm" not in df:
        df["arm"] = np.nan
    df["arm"] = df["arm"].where(df["arm"].notna(), df["controller"]).astype(str)
    if L4 is not None and len(d4) and S4.e0b_failed(L4):
        # v4's E0b failed: its re-flown v1, v2 and v2single replace v3's at 4.5 g (v4_plan.md section 4.0)
        df = df[~((df.source == "v3") & df.arm.isin(["v1", "v2", "v2single"]) & (df.load == G))]
        df["arm"] = df.arm.replace({"v1_v4run": "v1", "v2_v4run": "v2", "v2single_v4run": "v2single"})
    df["suite"] = df.get("suite", pd.Series(np.nan, index=df.index)).fillna("A256").astype(str)
    df["stage"] = df.get("stage", pd.Series(np.nan, index=df.index)).fillna("").astype(str)
    df.loc[df.source != "v5", "stage"] = df.loc[df.source != "v5", "source"]
    df["load"] = df["load"].fillna(OWN).astype(str).replace({"nan": OWN, "None": OWN})
    df["seed"] = pd.to_numeric(df["seed"], errors="coerce").fillna(-1).astype(int)
    df["mppi_rep"] = pd.to_numeric(df.get("mppi_rep"), errors="coerce").fillna(-1).astype(int)
    df["case"] = pd.to_numeric(df["case"], errors="coerce").astype(int)
    df["miss_km"] = pd.to_numeric(df["miss_km"], errors="coerce")
    df["peak_g"] = pd.to_numeric(df["peak_g"], errors="coerce")
    landed = df.outcome == "drogue"
    df["success_2km"] = landed & (df.miss_km < 2.0)
    df["success_5km"] = landed & (df.miss_km < 5.0)
    df["success_10km"] = landed & (df.miss_km <= 10.0)
    df["success_25km"] = landed & (df.miss_km < 25.0)
    df["log_miss"] = np.log(df.miss_km.clip(lower=cfg.log_miss_floor_km))
    df["block"] = pd.to_numeric(df.get("block"), errors="coerce")
    df["target_index"] = pd.to_numeric(df.get("target_index"), errors="coerce")
    df["key"] = [key_of(a, l) for a, l in zip(df.arm, df.load)]
    # A256: the conditions of this configuration only, with their feasibility under the true parameters
    a256 = {(int(c["case"]), c["atmosphere"]): c.get("feasible") for c in S5.conditions(cfg, "A256")}
    is_a = df.suite == "A256"
    keep = ~is_a | np.array([(c, a) in a256 for c, a in zip(df.case, df.atmosphere)])
    df = df[keep].copy()
    is_a = df.suite == "A256"
    feas = df.get("feasible", pd.Series(np.nan, index=df.index)).map(lambda v: np.nan if pd.isna(v) else _b(v))
    feas[is_a] = [a256.get((c, a)) for c, a in zip(df.case[is_a], df.atmosphere[is_a])]
    df["feasible"] = feas.astype(object)
    return df.reset_index(drop=True)


def main_flights(df: pd.DataFrame) -> pd.DataFrame:
    """Result flights: no gate, tuning, timing or learning-curve flights, one MPPI repetition."""
    m = ~df.stage.isin(NON_RESULT_STAGES) & (df.mppi_rep <= 0)
    m &= ~df.arm.str.startswith(NON_TABLE_PREFIX)
    return df[m]


# ============================================================================
# A256 comparisons
# ============================================================================
def cube(df: pd.DataFrame, al: tuple, col: str, seeds: Sequence[int], cases: Sequence[int]) -> np.ndarray:
    """A256 values [n_seeds | 1, n_atm, n_cases] (NaN where missing); deterministic arms [1, A, C]."""
    arm, load = al
    d = df[(df.arm == arm) & (df.load == load) & (df.suite == "A256")]
    det = len(d) > 0 and bool((d.seed < 0).all())
    ss = [-1] if det else list(seeds)
    X = np.full((len(ss), len(ATM), len(cases)), np.nan)
    si, ai, ci = {s: i for i, s in enumerate(ss)}, {a: i for i, a in enumerate(ATM)}, {c: i for i, c in enumerate(cases)}
    for s, a, c, v in zip(d.seed, d.atmosphere, d.case, d[col]):
        if s in si and a in ai and c in ci:
            X[si[s], ai[a], ci[c]] = float(v)
    return X


def present(df: pd.DataFrame, al: tuple, suite: str = "A256") -> bool:
    return bool(((df.arm == al[0]) & (df.load == al[1]) & (df.suite == suite)).any())


class Comparer:
    """Paired A256 comparisons A - B, pooled and per atmosphere, with per-seed tests (analyze_v4's scheme)."""

    def __init__(self, df, boot: A4.Boot, cases):
        self.df, self.boot, self.cases = df, boot, cases

    def __call__(self, name: str, a: tuple, b: tuple, metric: str, seeds: Sequence[int]) -> List[dict]:
        if not (present(self.df, a) and present(self.df, b)):
            return []
        col, fn, scale = {"success5": ("success_5km", "diff_mean", 100.0), "success2": ("success_2km", "diff_mean", 100.0),
                          "log_miss": ("log_miss", "paired_mean", 1.0)}[metric]
        XA, XB = cube(self.df, a, col, seeds, self.cases), cube(self.df, b, col, seeds, self.cases)
        if XA.shape[0] > 1 and XB.shape[0] > 1 and XA.shape[0] != XB.shape[0]:
            return []
        ps = A4.per_seed_tests(XA, XB, list(seeds), "success" if metric.startswith("success") else "logmiss")
        rows = []
        for i, t in [(None, "all")] + list(enumerate(ATM)):
            r = self.boot(f"v5|{name}|{a}|{b}|{metric}|{list(seeds)}", [XA, XB], fn, i)
            r = {k: (v * scale if k in ("diff", "ci_lo", "ci_hi", "ci90_lo", "ci90_hi") else v) for k, v in r.items()}
            row = {"id": name, "a": label(*a), "b": label(*b), "a_key": key_of(*a), "b_key": key_of(*b),
                   "metric": metric, "atmosphere": t, **r}
            if t == "all":
                k, n = A4._holds(ps, np.sign(r["diff"]) if np.isfinite(r["diff"]) and r["diff"] != 0 else 1.0)
                row.update({"seeds_holding": k, "n_seeds": n,
                            "per_seed": json.dumps([round(float(p["diff"]), 4) for p in ps])})
            rows.append(row)
        return rows


def _get(t: pd.DataFrame, id_: str, metric: str = "success5", atm: str = "all"):
    if t is None or not len(t):
        return None
    r = t[(t.id == id_) & (t.metric == metric) & (t.atmosphere == atm)]
    return r.iloc[0] if len(r) else None


# ============================================================================
# ORB-1 comparisons (each case one unit; cases resampled within training-seed strata)
# ============================================================================
def orb_paired(L, key: str, a: np.ndarray, b: np.ndarray, strata: np.ndarray, n_boot: int) -> dict:
    a, b, strata = np.asarray(a, float), np.asarray(b, float), np.asarray(strata, int)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b, strata = a[ok], b[ok], strata[ok]
    h = hashlib.sha1(json.dumps([key, n_boot]).encode())
    for x in (a, b, strata.astype(float)):
        h.update(np.ascontiguousarray(x).tobytes())
    h = h.hexdigest()
    path = L.cache / (re.sub(r"[^A-Za-z0-9_.@-]+", "_", key) + ".json")
    if path.exists():
        try:
            d = json.loads(path.read_text())
            if d.get("hash") == h:
                return d["result"]
        except (OSError, ValueError):
            pass
    n = len(a)
    if n == 0:
        res = {"diff": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "p_boot": np.nan, "n": 0}
    else:
        d = a - b
        rng = np.random.default_rng(zlib.crc32(key.encode()) & 0x7FFFFFFF)
        tot = np.zeros(n_boot)
        for s in np.unique(strata):
            ds = d[strata == s]
            tot += ds[rng.integers(0, len(ds), (n_boot, len(ds)))].sum(1)
        v = tot / n
        lo, hi = np.percentile(v, [2.5, 97.5])
        p = min(1.0, 2 * min((v <= 0).mean(), (v >= 0).mean()))
        res = {"diff": float(d.mean()), "ci_lo": float(lo), "ci_hi": float(hi),
               "p_boot": float(max(p, 1.0 / n_boot)), "n": int(n)}
    L.write_json(path, {"key": key, "hash": h, "result": res})
    return res


# ============================================================================
# Telemetry-derived ORB-1 quantities (landing errors, g-duration curves), cached per flight
# ============================================================================
def orb1_telemetry(L, df_orb: pd.DataFrame) -> pd.DataFrame:
    """Per ORB-1 flight: landing error (downrange, crossrange, km) in the nominal entry frame, g-duration points."""
    path = L.cache / "orb1_telemetry.csv"
    have = pd.read_csv(path) if path.exists() else pd.DataFrame(columns=["flight_id"])
    known = set(have.flight_id)
    spec = json.loads((B.BENCH_DIR / "ORB1_spec.json").read_text())
    ic = B.ORB1Spec().nominal_ic()
    Q = TS.entry_frame(math.radians(ic["phi_deg"]), math.radians(ic["lambda_deg"]), math.radians(ic["chi_deg"]))
    tgt = {int(t["target_index"]): TS.to_local(Q, math.radians(t["phi_deg"]), math.radians(t["lambda_deg"]))
           for t in spec["targets"]}
    rows = []
    for r in df_orb[~df_orb.flight_id.isin(known)].itertuples(index=False):
        try:
            t = S.load_telemetry(L, r.flight_id)
        except (OSError, ValueError, KeyError):
            continue
        dr, cr = TS.to_local(Q, math.radians(float(t["phi_deg"][-1])), math.radians(float(t["lambda_deg"][-1])))
        td, tc = tgt.get(int(r.target_index) if np.isfinite(r.target_index) else 0, (np.nan, np.nan))
        tt = np.asarray(t["t_s"], float)
        dt = float(np.median(np.diff(tt))) if len(tt) > 1 else 1.0
        gd = B.g_duration(np.asarray(t["g"], float)[1:], dt)
        rows.append({"flight_id": r.flight_id, "dr_err_km": dr - td, "cr_err_km": cr - tc,
                     "g_points": json.dumps(gd["points"])})
    if rows:
        have = pd.concat([have, pd.DataFrame(rows)], ignore_index=True)
        L.write_csv(have, path)
    return have


# ============================================================================
# Criteria (T05)
# ============================================================================
def criteria(cfg, df, t04, t06, t07, mech, orb: dict, gates: dict, seeds, rl_seeds) -> pd.DataFrame:
    rows = []

    def row(id_, fam, q, measured, met, value=np.nan, ci=(np.nan, np.nan), p=None, holding=None, valid=None):
        rows.append({"id": id_, "family": fam, "question": q, "measured": measured, "value": value, "ci_lo": ci[0],
                     "ci_hi": ci[1], "p": p, "p_holm": None, "seeds_holding": holding, "met": met, "valid": valid})

    def per_atm(t, id_, metric="success5", fmt="{:+.1f}"):
        return ", ".join(f"{ATM_LAB[a]} {fmt.format(r['diff'])} [{fmt.format(r['ci_lo'])}, {fmt.format(r['ci_hi'])}]"
                         for a in ATM if (r := _get(t, id_, metric, a)) is not None)

    def s5_text(r, unit="pp"):
        return (f"success@5km {r['diff']:+.1f} {unit} [{r['ci_lo']:+.1f}, {r['ci_hi']:+.1f}], holds in "
                f"{int(r['seeds_holding'])}/{int(r['n_seeds'])} seeds")

    # ---- RQ1 ----
    for id_, q in (("G1", "Does v2 beat FNPEG? (v2single − fnpeg)"), ("G2", "Does v2 beat PredGuid? (v2single − predguid)"),
                   ("G4", "Does v2 beat FNPEG at FNPEG's own load? (v2single@fnpeg_matched − fnpeg)"),
                   ("G5", "Does the learned model beat explicit estimation? (v2single − mppi_filtered)")):
        r = _get(t04, id_)
        if r is None:
            row(id_, "RQ1", q, "not run", False)
            continue
        worse = [a for a in ATM if (x := _get(t04, id_, "success5", a)) is not None and x["ci_hi"] < 0]
        need = 4 if int(r["n_seeds"]) >= 5 else max(1, int(r["n_seeds"]) - 1)
        met = bool(r["diff"] > 0 and r["ci_lo"] > 0 and int(r["seeds_holding"]) >= need and not worse)
        row(id_, "RQ1", q, s5_text(r) + "; " + per_atm(t04, id_) + (f"; B significantly better in {worse}" if worse else ""),
            met, r["diff"], (r["ci_lo"], r["ci_hi"]), r["p_boot"], f"{int(r['seeds_holding'])}/{int(r['n_seeds'])}")
    # G3: validity of the load comparison
    a = df[df.suite == "A256"]
    pk = {k: a[a.key == k].groupby("atmosphere").peak_g.mean() for k in ("fnpeg", "v2single@fnpeg_matched")}
    if len(pk["fnpeg"]) and len(pk["v2single@fnpeg_matched"]):
        f_ok = bool((pk["fnpeg"] <= cfg.tune_peak_g_max + 1e-9).all())
        m_ok = all(pk["v2single@fnpeg_matched"].get(t, np.inf) <= pk["fnpeg"].get(t, -np.inf) + 0.1 for t in pk["fnpeg"].index)
        row("G3", "RQ1", "Is the load comparison valid? (fnpeg mean peak g ≤ 4.6 g; matched arm ≤ fnpeg + 0.1 g, per atmosphere)",
            "; ".join(f"{ATM_LAB[t]}: fnpeg {pk['fnpeg'][t]:.2f} g, matched {pk['v2single@fnpeg_matched'].get(t, np.nan):.2f} g"
                      for t in pk["fnpeg"].index), f_ok and m_ok)
    else:
        row("G3", "RQ1", "Is the load comparison valid?", "not run", False)
    g3 = rows[-1]["met"]
    for r_ in rows:
        if r_["id"] == "G4":
            r_["valid"] = bool(g3)
    # G6: competence
    gp = {k: v.get("pass") for k, v in gates.items()}
    parts = [f"gates {', '.join(f'{k} {_pass(v)}' for k, v in gp.items())}"]
    for id_ in ("G6a", "G6b"):
        r = _get(t04, id_)
        if r is not None:
            parts.append(f"{r['a']} − classical: {per_atm(t04, id_)}")
    row("G6", "RQ1", "Are the baselines competent? (V1-V4 passed; fnpeg and predguid vs classical per atmosphere)",
        "; ".join(parts), all(v is True for k, v in gp.items() if k in ("V1", "V2", "V3", "V4")))
    for id_, q in (("G7a", "What does the 90° handle cost FNPEG? (fnpeg_180 − fnpeg, exploratory)"),
                   ("G7b", "What does the 90° handle cost PredGuid? (predguid_180 − predguid, exploratory)"),
                   ("G7c", "How far is FNPEG from its own ceiling? (fnpeg_oracle − fnpeg, exploratory)")):
        r = _get(t04, id_)
        row(id_, "RQ1", q, "not run" if r is None else f"success@5km {r['diff']:+.1f} pp [{r['ci_lo']:+.1f}, "
            f"{r['ci_hi']:+.1f}]; {per_atm(t04, id_)}", None, *(() if r is None else (r["diff"], (r["ci_lo"], r["ci_hi"]))))
    _holm(rows, ("G1", "G2", "G4", "G5"))

    # ---- RQ2 ----
    jr = []
    for arm in M5.ARMS_E3:
        id_ = f"J-{arm}"
        r = _get(t06, id_, "log_miss")
        q = f"Does R0 beat {arm}? (paired log miss {arm} − R0)"
        if r is None:
            row(id_, "RQ2", q, "not run", False)
            continue
        need = 4 if int(r["n_seeds"]) >= 5 else max(1, int(r["n_seeds"]) - 1)
        met = bool(r["diff"] > 0 and r["ci_lo"] > 0 and int(r["seeds_holding"]) >= need)
        row(id_, "RQ2", q, f"mean paired log miss {r['diff']:+.3f} [{r['ci_lo']:+.3f}, {r['ci_hi']:+.3f}] "
            f"(miss ratio ×{math.exp(r['diff']):.2f}), holds in {int(r['seeds_holding'])}/{int(r['n_seeds'])} seeds; "
            + per_atm(t06, id_, "log_miss", "{:+.3f}"), met, r["diff"], (r["ci_lo"], r["ci_hi"]), r["p_boot"],
            f"{int(r['seeds_holding'])}/{int(r['n_seeds'])}")
        jr.append(id_)
    _holm(rows, jr)
    r = _get(t06, "J-r0_nojepa", "log_miss")
    if r is not None:
        eq = bool(-cfg.equiv_log < r["ci90_lo"] and r["ci90_hi"] < cfg.equiv_log)
        each = next((x["met"] and (x["p_holm"] or 1) < 0.05 for x in rows if x["id"] == "J-r0_nojepa"), False)
        ans = "the JEPA objective matters (R0 better)" if each else ("equivalent: the objective does not matter" if eq
                                                                      else "inconclusive")
        row("J-equiv", "RQ2", f"Is r0_nojepa equivalent to R0? (90% CI inside ±{cfg.equiv_log})",
            f"90% CI [{r['ci90_lo']:+.3f}, {r['ci90_hi']:+.3f}]: {ans}", eq, r["diff"], (r["ci90_lo"], r["ci90_hi"]))
    else:
        row("J-equiv", "RQ2", "Is r0_nojepa equivalent to R0?", "not run", False)
    if mech.get("rho") is not None:
        row("J-mech", "RQ2", f"Does open-loop error at {mech['ttg_s']:.0f} s predict closed-loop log miss? (Spearman ρ "
            f"across {mech['n_models']} models)", f"ρ = {mech['rho']:+.2f} [{mech['ci_lo']:+.2f}, {mech['ci_hi']:+.2f}]",
            None, mech["rho"], (mech["ci_lo"], mech["ci_hi"]))
    else:
        row("J-mech", "RQ2", "Does open-loop error predict closed-loop log miss?", "fewer than 3 models with both", None)

    # ---- RQ3 ----
    rr = []
    for arm in S5.RL_TIER1 + S5.RL_TIER2:
        id_ = f"R-{arm}"
        r = _get(t07, id_)
        if r is None:
            if arm in S5.RL_TIER1:
                row(id_, "RQ3", f"Does v2 beat {arm}?", "not run", False)
            continue
        n = int(r["n_seeds"])
        met = bool(r["diff"] > 0 and r["ci_lo"] > 0 and int(r["seeds_holding"]) >= n)
        row(id_, "RQ3" if arm in S5.RL_TIER1 else "RQ3 (Tier 2, outside the Holm family)",
            f"Does v2 beat {arm}? (v2single − {arm}, seeds {list(rl_seeds)})", s5_text(r), met, r["diff"],
            (r["ci_lo"], r["ci_hi"]), r["p_boot"], f"{int(r['seeds_holding'])}/{n}")
        if arm in S5.RL_TIER1:
            rr.append(id_)
    _holm(rows, rr)
    r = _get(t07, "R-imit")
    row("R-imit", "RQ3", "Imitation gap (dagger_fnpeg − fnpeg, estimate)", "not run" if r is None else
        f"success@5km {r['diff']:+.1f} pp [{r['ci_lo']:+.1f}, {r['ci_hi']:+.1f}]", None,
        *(() if r is None else (r["diff"], (r["ci_lo"], r["ci_hi"]))))
    gaps = t07.attrs.get("gaps", {}) if t07 is not None else {}
    row("R-gap", "RQ3", "Sim-to-sim gap (success@5km fast env − full simulator, same policy and A256 cases)",
        "; ".join(f"{a} {v:+.1f} pp" for a, v in gaps.items()) or "not run", None)

    # ---- RQ4 ----
    sc = orb.get("score")
    if sc is not None and len(sc):
        for r in sc[sc.group == "all"].itertuples():
            row(f"Q1-{r.arm}", "RQ4", f"Does {r.arm} meet the accuracy requirement? (P(miss ≤ 10 km) LCB ≥ "
                f"{B.EM1_REQUIREMENT})", f"{int(r.n)} flights, {int(r.n_miss_10km)} beyond 10 km, LCB {r.requirement_lcb:.5f}",
                bool(r.requirement_met), r.requirement_lcb)
        for r in sc[sc.group == "all"].itertuples():
            row(f"Q2-{r.arm}", "RQ4", f"Does {r.arm} respect the crew g-duration limits? (descriptive)",
                f"{100 * r.g_compliant_frac:.1f}% of flights compliant [{100 * r.g_compliant_lo:.1f}, "
                f"{100 * r.g_compliant_hi:.1f}]" if np.isfinite(r.g_compliant_frac) else "no g history",
                None, r.g_compliant_frac, (r.g_compliant_lo, r.g_compliant_hi))
    q3 = []
    for id_, other in (("Q3a", "fnpeg"), ("Q3b", "predguid")):
        r = orb.get(id_)
        q = f"Does v2 beat {other} on ORB-1 block 1? (v2single − {other})"
        if not r or not r.get("n"):
            row(id_, "RQ4", q, "not run", False)
            continue
        met = bool(r["diff"] > 0 and r["ci_lo"] > 0)
        row(id_, "RQ4", q, f"success@5km {100 * r['diff']:+.1f} pp [{100 * r['ci_lo']:+.1f}, {100 * r['ci_hi']:+.1f}] "
            f"over {r['n']} cases", met, 100 * r["diff"], (100 * r["ci_lo"], 100 * r["ci_hi"]), r["p_boot"])
        q3.append(id_)
    _holm(rows, q3)
    return pd.DataFrame(rows)


def _pass(v) -> str:
    return "pass" if v is True else ("FAIL" if v is False else "pending")


def _holm(rows: List[dict], ids: Sequence[str]) -> None:
    """Holm-adjust p within a family; a pass/fail criterion also needs p_holm < 0.05."""
    fam = [r for r in rows if r["id"] in ids and r["p"] is not None and np.isfinite(r["p"])]
    if not fam:
        return
    adj = A3.holm([r["p"] for r in fam])
    for r, p in zip(fam, adj):
        r["p_holm"] = float(p)
        if r["met"] is not None:
            r["met"] = bool(r["met"] and p < 0.05)


# ============================================================================
# Small tables
# ============================================================================
def gates_table(L) -> Tuple[pd.DataFrame, dict]:
    rows, out = [], {}
    for g in ("S0", "V1", "V2", "V3", "V4", "E0"):
        p = L.checks / f"{g}.json"
        if not p.exists():
            rows.append({"gate": g, "pass": None, "done": False, "summary": "not run"})
            continue
        d = json.loads(p.read_text())
        out[g] = d
        keep = {k: v for k, v in d.items() if not isinstance(v, (list, dict)) and k not in ("gate", "git_commit")}
        rows.append({"gate": g, "pass": d.get("pass"), "done": (L.checks / f"{g}.done").exists(),
                     "summary": json.dumps(keep, default=str)})
    return pd.DataFrame(rows), out


def tuning_table(L) -> pd.DataFrame:
    parts = []
    for m in sorted(p for p in L.tuning.glob("*") if p.is_dir()) if L.tuning.exists() else []:
        f = m / "configs.csv"
        if f.exists():
            d = pd.read_csv(f)
            ch = json.loads((m / "chosen.json").read_text()) if (m / "chosen.json").exists() else {}
            d.insert(0, "method", m.name)
            d["chosen"] = d.config_index == ch.get("config_index", -1)
            parts.append(d)
    vf = L.tuning / "fnpeg" / "vf.json"
    out = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["method", "config_index"])
    if vf.exists():
        out.attrs["vf"] = json.loads(vf.read_text())
    return out


def threshold_table(cfg, main: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (suite, key), g in main.groupby(["suite", "key"]):
        landed = (g.outcome == "drogue").to_numpy()
        miss = g.miss_km.to_numpy()
        for thr in cfg.thresholds_km:
            rows.append({"suite": suite, "arm": key, "threshold_km": thr, "n": len(g),
                         "success_%": 100 * float((landed & (miss < thr)).mean())})
    return pd.DataFrame(rows)


def mechanism(cfg, df, ol: pd.DataFrame, seeds) -> dict:
    """Spearman ρ across R0, v1 and the E3 models: open-loop median error vs closed-loop mean log miss (A256, 4.5 g)."""
    if not len(ol):
        return {"rho": None, "n_models": 0, "table": pd.DataFrame()}
    av = sorted(ol.ttg_s.unique())
    ttg = 300.0 if 300.0 in av else float(max(av))
    x = ol[(ol.split == "test") & (ol.ttg_s == ttg)]
    rows = []
    for m in ["R0", "v1"] + list(M5.ARMS_E3):
        for s in seeds:
            xo = x[(x.model == m) & (x.seed == s)].err_median_km
            d = df[(df.arm == m) & (df.load == G) & (df.seed == s) & (df.suite == "A256")]
            if len(xo) and len(d):
                rows.append({"model": m, "seed": s, "ol_err_km": float(xo.iloc[0]), "log_miss": float(d.log_miss.mean())})
    tab = pd.DataFrame(rows)
    if not len(tab) or tab.model.nunique() < 3:
        return {"rho": None, "n_models": int(tab.model.nunique()) if len(tab) else 0, "table": tab, "ttg_s": ttg}
    px = tab.pivot_table(index="model", columns="seed", values="ol_err_km")
    py = tab.pivot_table(index="model", columns="seed", values="log_miss")
    rho = A3.spearman(px.mean(1), py.mean(1))
    rng = np.random.default_rng(zlib.crc32(b"v5 mechanism"))
    vals = []
    for _ in range(cfg.n_boot):
        si = rng.integers(0, px.shape[1], px.shape[1])
        vals.append(A3.spearman(px.iloc[:, si].mean(1), py.iloc[:, si].mean(1)))
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    lo, hi = np.percentile(v, [2.5, 97.5]) if len(v) else (np.nan, np.nan)
    return {"rho": float(rho), "ci_lo": float(lo), "ci_hi": float(hi), "n_models": int(px.shape[0]), "table": tab,
            "ttg_s": ttg}


def _params_of(L, arm: str, seed: int) -> Optional[int]:
    p = L.model_dir(arm, seed) / "model.json"
    if p.exists():
        return json.loads(p.read_text()).get("params")
    return None


# ============================================================================
# The analysis
# ============================================================================
def analyze(cfg: S5.StudyConfigV5, L: S5.LayoutV5, L3: S.Layout, L4=None, log=print) -> dict:
    """Every table, figure and criterion. Returns a summary dict for the report."""
    caps = A4.Captions(L)
    boot = A4.Boot(L, cfg.n_boot)
    df = load_flights_v5(cfg, L, L3, L4)
    main = main_flights(df)
    cases = sorted({int(c["case"]) for c in S5.conditions(cfg, "A256")})
    seeds = list(cfg.seeds)
    rl_seeds = list(cfg.rl_seeds)
    comp = Comparer(main, boot, cases)
    out: Dict[str, object] = {"n_flights": int(len(df)), "n_v5_flights": int((df.source == "v5").sum()), "boot": boot,
                              "seeds": seeds}
    data = "04_flights/flights.csv (+ the v3 and v4 folders' flights, read only)"

    # ---- T01 / T02 ----
    sc = B.score(main.assign(arm=main.key)) if len(main) else pd.DataFrame()
    if len(sc):
        ps = B.score(main.assign(arm=main.key), by_seed=True)
        ps = ps[ps.seed >= 0]
        spread = ps.groupby(["suite", "arm", "group"]).success_5km.agg(success_5km_seed_min="min",
                                                                      success_5km_seed_max="max",
                                                                      n_seeds="count").reset_index()
        sc = sc.merge(spread, on=["suite", "arm", "group"], how="left")
    _table(L, caps, sc.round(5), "T01_metrics_by_arm",
           "Every section-5 metric per arm × suite × atmosphere (A256, DEV48 excluded) or target range (ORB-1), plus the "
           "arm's pooled row (group all): success at 2/5/10/25 km, mean log miss, median and p95 miss, peak g, heat, RCS, "
           "reversals, planning ms, authority saturation, the ORB-1 requirement LCB and g-duration compliance with its "
           "Clopper–Pearson interval; seed min/max of success@5km for seeded arms. v3 and v4 arms included.", data)
    fe = main[main.feasible.map(lambda v: v is True or _b(v))]
    t02 = B.score(fe.assign(arm=fe.key)) if len(fe) else pd.DataFrame()
    _table(L, caps, t02.round(5), "T02_metrics_feasible_only",
           "T01 on conditions whose target is reachable within the load limit under the true parameters (A256: the "
           "conservative 4.5 g footprint; ORB-1: constant-bank feasibility at 4.5 g).", data + " + configs/benchmark_v5")
    out["T01"] = sc

    # ---- T03 gates ----
    t03, gates = gates_table(L)
    _table(L, caps, t03, "T03_gates", "Gates S0 (benchmark), V1–V4 (FNPEG and PredGuid verification) and E0 (runner "
           "equivalence): pass flag and the measured numbers (01_checks/<gate>.json).", "01_checks/")
    out["gates"] = gates

    # ---- T04 RQ1 ----
    rq1 = [("G1", V2S, FN), ("G2", V2S, PG), ("G4", V2Sm, FN), ("G5", V2S, MF), ("G6a", FN, CL), ("G6b", PG, CL),
           ("G7a", FN180, FN), ("G7b", PG180, PG), ("G7c", FNO, FN)]
    t04 = pd.DataFrame([r for nm, a, b in rq1 for m in ("success5", "log_miss") for r in comp(nm, a, b, m, seeds)])
    _table(L, caps, t04.round(5), "T04_rq1_paired_differences",
           "RQ1 at 4.5 g: paired differences A − B, pooled and per atmosphere, with 95% and 90% hierarchical-bootstrap "
           f"intervals ({cfg.n_boot} resamples: seeds, then cases within each atmosphere), bootstrap p, and the seeds "
           "where the pooled effect holds. success5 in percentage points; log_miss = mean paired difference of ln(miss). "
           "G1/G2/G4/G5 are the Holm family; G6 the competence check; G7 exploratory.", data)
    out["T04"] = t04

    # ---- T06 RQ2 ----
    rows = [r for a in E3 for m in ("log_miss", "success5") for r in comp(f"J-{a[0]}", a, R0, m, seeds)]
    t06 = pd.DataFrame(rows)
    ol = pd.concat([A4.open_loop_table(L)] + ([A4.open_loop_table(L4)] if L4 is not None else []), ignore_index=True)
    ol_s = ol[ol.split == "test"].groupby(["model", "ttg_s"]).err_median_km.mean().unstack("ttg_s") if len(ol) else None
    if len(t06):
        pr = {a: _params_of(L, a, seeds[0]) for a in M5.ARMS_E3}
        t06["params_a"] = [pr.get(k) for k in t06.a_key]
        r0p = next((json.loads((L.model_dir(a, seeds[0]) / "model.json").read_text()).get("r0_params")
                    for a in M5.ARMS_E3 if (L.model_dir(a, seeds[0]) / "model.json").exists()), None)
        t06["params_R0"] = r0p
        if ol_s is not None:
            for ttg in ol_s.columns:
                t06[f"ol_err_{int(ttg)}s_a"] = [ol_s[ttg].get(k, np.nan) for k in t06.a_key]
                t06[f"ol_err_{int(ttg)}s_R0"] = ol_s[ttg].get("R0", np.nan)
    _table(L, caps, t06.round(5), "T06_rq2_vs_r0",
           "RQ2: each non-JEPA model vs R0 (v4's one-member v2 recipe) under the same MPPI at 4.5 g on A256: mean paired "
           "log-miss difference (model − R0; > 0 means R0 misses less) and success@5km, with bootstrap intervals; "
           "parameter counts; open-loop median footprint error on the test split by time to go (mean over seeds).",
           data + " + 03_open_loop/ (v5 and v4)")
    if len(ol):
        _table(L, caps, ol.round(4), "T06b_open_loop", "Open-loop footprint error (median and p90, km) per model, seed, "
               "split and time to go: the E3 models (v5) and R0 / v1 (v4).", "03_open_loop/")
    mech = mechanism(cfg, df, ol, seeds)
    out["T06"], out["mech"] = t06, mech

    # ---- T07 RQ3 ----
    rows = []
    for a in RL1 + RL2:
        rows += comp(f"R-{a[0]}", V2S, a, "success5", rl_seeds)
    rows += comp("R-imit", ("dagger_fnpeg", OWN), FN, "success5", rl_seeds)
    t07 = pd.DataFrame(rows)
    rl_rows, gaps = [], {}
    for arm, _ in RL1 + RL2:
        for s in rl_seeds:
            d = L.policy_dir(arm, s)
            if not (d / "policy.pt").exists():
                continue
            run = json.loads((d / "run.json").read_text()) if (d / "run.json").exists() else {}
            fe_ = json.loads((d / "fast_eval.json").read_text()) if (d / "fast_eval.json").exists() else {}
            full = main[(main.arm == arm) & (main.seed == s) & (main.suite == "A256")]
            fa = (fe_.get("A256") or {}).get("success_5km_fast", np.nan)
            rl_rows.append({"arm": arm, "seed": s, "full_success5_%": 100 * full.success_5km.mean() if len(full) else np.nan,
                            "full_median_miss_km": full.miss_km.median() if len(full) else np.nan,
                            "fast_A256_success5_%": 100 * fa,
                            "fast_train_dist_success5_%": 100 * (fe_.get("train_distribution") or {}).get(
                                "success_5km_fast", np.nan),
                            "decisions": run.get("decisions"), "iterations": run.get("iterations"),
                            "episodes": run.get("episodes"), "train_s": run.get("train_s", run.get("wall_s"))})
    t07b = pd.DataFrame(rl_rows)
    if len(t07b):
        t07b["gap_fast_minus_full_pp"] = t07b["fast_A256_success5_%"] - t07b["full_success5_%"]
        gaps = t07b.groupby("arm").gap_fast_minus_full_pp.mean().round(2).to_dict()
    t07.attrs["gaps"] = gaps
    _table(L, caps, t07.round(5), "T07_rq3_rl", "RQ3: v2single − each RL arm (success@5km, A256, 4.5 g; RL training "
           f"seeds {rl_seeds} paired with v2 seeds {rl_seeds}), and the imitation gap dagger_fnpeg − fnpeg.", data)
    _table(L, caps, t07b.round(4), "T07b_rl_runs", "RL runs: full-simulator success@5km and median miss on A256, the "
           "fast environment's success on the same A256 cases and on its training distribution, the fast-vs-full gap "
           "(pp), and the training budget used (decisions, DAgger iterations, episodes, seconds).",
           "02_policies/<arm>/seed_<k>/{run.json, fast_eval.json} + flights")
    out["T07"] = t07

    # ---- T08 RQ4 ----
    orb = {}
    od = main[main.suite == "ORB1"]
    if len(od):
        s8 = B.score(od.assign(arm=od.key))
        blocks = od.groupby("key").block.agg(lambda b: ",".join(str(int(x) + 1) for x in sorted(b.dropna().unique())))
        s8["blocks_flown"] = s8.arm.map(blocks)
        orb["score"] = s8
        b1 = od[od.block == 0]
        piv = b1.pivot_table(index="case", columns="key", values="success_5km", aggfunc="first")
        strata = (piv.index.to_numpy() % len(cfg.seeds))
        for id_, other in (("Q3a", "fnpeg"), ("Q3b", "predguid")):
            if "v2single" in piv and other in piv:
                orb[id_] = orb_paired(L, f"v5|{id_}|orb1_block1", piv["v2single"].astype(float).to_numpy(),
                                      piv[other].astype(float).to_numpy(), strata, cfg.n_boot)
    else:
        s8 = pd.DataFrame()
    _table(L, caps, s8.round(5), "T08_rq4_orb1",
           "RQ4, ORB-1: per arm and target range (and pooled): blocks flown (sequential stopping rule), success at "
           f"2/5/10 km, the one-sided 90% Clopper–Pearson lower bound on P(miss ≤ 10 km) against {B.EM1_REQUIREMENT}, "
           "g-duration compliance (NASA-STD-3000 deconditioned crew) with its interval.", data)
    if any(k in orb for k in ("Q3a", "Q3b")):
        out["Q3"] = {k: orb[k] for k in ("Q3a", "Q3b") if k in orb}

    # ---- T05 criteria ----
    t05 = criteria(cfg, main, t04, t06, t07, mech, orb, gates, seeds, rl_seeds)
    _table(L, caps, t05, "T05_criteria", "Pre-registered criteria (v5_plan.md section 8): measured value, interval, "
           "bootstrap p, Holm-adjusted p within each family (RQ1: G1, G2, G4, G5; RQ2: the four J-each; RQ3: the Tier 1 "
           "R-each; RQ4: Q3a, Q3b), seeds holding, met (blank = an estimate, not pass/fail), valid (G4 needs G3).",
           "05_analysis/tables/T04, T06, T07, T08")
    out["T05"] = t05

    # ---- T09 compute ----
    cp = L.checks / "compute_per_decision.csv"
    t09 = pd.read_csv(cp) if cp.exists() else pd.DataFrame(columns=["arm", "device", "ms_median", "ms_p99", "params"])
    _table(L, caps, t09, "T09_rq5_compute", "RQ5: compute per guidance decision, median and p99 ms, single CPU thread "
           "and GPU, one A256 flight per arm (nominal air), and parameter counts.", "01_checks/compute_per_decision.csv")

    # ---- T10 tuning ----
    t10 = tuning_table(L)
    vf = t10.attrs.get("vf", {})
    _table(L, caps, t10.round(5), "T10_tuning", "DEV48 tuning (section 4.5): every configuration tried per method with "
           "its sampled settings and DEV48 scores (success@5km, mean log miss, mean peak g); chosen = the frozen choice. "
           + (f"FNPEG final speed Vf = {vf.get('vf_mps')} m/s (median at 7.6 km over {vf.get('n_flights')} classical "
              "DEV48 flights)." if vf else ""), "03_tuning/")

    # ---- T11 legacy ----
    inv = L.study / "legacy_rl_inventory.csv"
    t11 = pd.read_csv(inv) if inv.exists() else pd.DataFrame(columns=["path", "kind"])
    lg = df[df.arm.str.startswith("rl_legacy_")]
    if len(lg):
        ev = lg.groupby(["arm", "suite"]).agg(n=("miss_km", "size"), success5_pct=("success_5km", lambda s: 100 * s.mean()),
                                              success25_pct=("success_25km", lambda s: 100 * s.mean()),
                                              median_miss_km=("miss_km", "median"),
                                              landed_pct=("outcome", lambda s: 100 * (s == "drogue").mean())).reset_index()
        ev["note"] = np.where(ev.suite == "A256", "out of its training distribution", "its own held-out set")
        t11 = pd.concat([t11.assign(part="inventory"), ev.assign(part="evaluation")], ignore_index=True)
    _table(L, caps, t11, "T11_legacy_rl", "Legacy RL (section 4.6): every file found under MyDrive/ReEntryAI/ (path, run, "
           "steps, size, date) and each final model re-evaluated with its own settings on heldout_testset.csv (first "
           f"{cfg.legacy_heldout_n} cases) and on A256 (out of its training distribution).",
           "00_study/legacy_rl_inventory.csv + flights")

    # ---- T12 thresholds ----
    t12 = threshold_table(cfg, main[main.suite.isin(["A256", "ORB1"])])
    _table(L, caps, t12, "T12_success_vs_threshold", "Success (outcome drogue and miss below the threshold) vs threshold, "
           "0.5–25 km, per arm and suite (pooled over seeds and atmospheres or targets).", data)

    # ---- figures ----
    try:
        out["figures"] = figures(cfg, L, L3, caps, df, main, sc, t04, t06, t12, ol, mech, seeds, rl_seeds, cases)
    except Exception as e:                           # a figure must never cost the tables
        log(f"[analysis] figures failed: {type(e).__name__}: {e}")
        raise
    out["captions"] = caps.d
    log(f"[analysis] {len(df)} flights; bootstrap: {boot.computed} computed, {boot.reused} reused from cache")
    return out


# ============================================================================
# Figures
# ============================================================================
def figures(cfg, L, L3, caps, df, main, sc, t04, t06, t12, ol, mech, seeds, rl_seeds, cases) -> Dict[str, Path]:
    plt = A4._plt()
    figs: Dict[str, Path] = {}
    a256 = main[main.suite == "A256"]
    order = ["classical", "fnpeg", "predguid", "mppi_filtered", "v2single", "v2single@fnpeg_matched", "R0"] + \
        list(M5.ARMS_E3) + list(S5.RL_TIER1) + list(S5.RL_TIER2)
    arms = [k for k in order if (a256.key == k).any()]

    # F01: success@5km by arm and atmosphere, seed whiskers
    fig, ax = plt.subplots(figsize=(11, 4.2))
    w = 0.8 / max(len(arms), 1)
    for i, k in enumerate(arms):
        d = a256[a256.key == k]
        for j, t in enumerate(ATM):
            dt = d[d.atmosphere == t]
            if not len(dt):
                continue
            per = dt.groupby("seed").success_5km.mean() * 100
            x = j + (i - (len(arms) - 1) / 2) * w
            ax.bar(x, per.mean(), w, color=color(k), label=label(*_split(k)) if j == 0 else None)
            if len(per) > 1:
                ax.plot([x, x], [per.min(), per.max()], color="#222", lw=0.8)
    ax.set_xticks(range(len(ATM)), [ATM_LAB[t] for t in ATM])
    ax.set_ylabel("success@5km (%)")
    ax.set_ylim(0, 100)
    ax.legend(fontsize=6.5, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.1), frameon=False)
    figs["F01"] = _fig(L, caps, plt, fig, "F01_success_by_arm",
                       "Success@5km on A256 by arm and atmosphere: mean over seeds (bars) with the seed min–max "
                       "(whiskers); deterministic arms fly once.", "T01")

    # F02: success vs mean peak g
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    pts = []
    for arm, load in LOAD_CURVE:
        d = a256[(a256.arm == arm) & (a256.load == load)]
        if len(d):
            pts.append((d.peak_g.mean(), 100 * d.success_5km.mean(), load))
    if pts:
        pts.sort()
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "-o", color=color("v2"), label="v2 ensemble, soft g limits (v3/v4)")
        for x, y, l in pts:
            ax.annotate(l, (x, y), fontsize=6.5, xytext=(3, 3), textcoords="offset points")
    for k in ("classical", "fnpeg", "predguid", "mppi_filtered", "v2single", "v2single@fnpeg_matched"):
        d = a256[a256.key == k]
        if len(d):
            ax.scatter(d.peak_g.mean(), 100 * d.success_5km.mean(), color=color(k), s=40, label=label(*_split(k)), zorder=3)
    ax.axvline(4.5, color="#999", lw=0.8, ls="--")
    ax.set_xlabel("mean peak g")
    ax.set_ylabel("success@5km (%)")
    ax.legend(fontsize=7, frameon=False)
    figs["F02"] = _fig(L, caps, plt, fig, "F02_success_vs_load",
                       "Success@5km vs mean peak g on A256 (all atmospheres): v2's load curve (v3/v4 soft limits), "
                       "FNPEG, PredGuid, MPPI + filtered physics, v2 single member, and v2 at FNPEG's realized peak g.",
                       "T01")

    # F03: RQ2 forest plot
    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    t = t06[(t06.metric == "log_miss")] if len(t06) else t06
    labs = []
    for i, k in enumerate(M5.ARMS_E3):
        for j, at in enumerate(["all"] + ATM):
            r = t[(t.a_key == k) & (t.atmosphere == at)] if len(t) else []
            if not len(r):
                continue
            r = r.iloc[0]
            y = i * 6 + j
            ax.errorbar(r["diff"], -y, xerr=[[r["diff"] - r["ci_lo"]], [r["ci_hi"] - r["diff"]]], fmt="o",
                        color=color(k), ms=4 if at == "all" else 2.5, lw=1.2 if at == "all" else 0.7)
            labs.append((-y, f"{LAB.get(k, k)}" if at == "all" else f"   {ATM_LAB[at]}"))
    ax.axvline(0, color="#555", lw=0.8)
    for s in (-cfg.equiv_log, cfg.equiv_log):
        ax.axvline(s, color="#bbb", lw=0.6, ls=":")
    if labs:
        ax.set_yticks([y for y, _ in labs], [l for _, l in labs], fontsize=6.5)
    ax.set_xlabel("mean paired log-miss difference, model − R0 (> 0: R0 misses less)")
    figs["F03"] = _fig(L, caps, plt, fig, "F03_rq2_forest",
                       "RQ2 forest plot: mean paired log-miss difference of each non-JEPA model vs R0 (same MPPI, 4.5 g, "
                       "A256), pooled (large) and per atmosphere, 95% bootstrap intervals; dotted lines: the ±0.15 "
                       "equivalence margin.", "T06")

    # F04: open loop vs closed loop
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    tab = mech.get("table")
    if tab is not None and len(tab):
        g = tab.groupby("model")[["ol_err_km", "log_miss"]].mean()
        for m, r in g.iterrows():
            ax.scatter(r.ol_err_km, r.log_miss, color=color(m), s=40)
            ax.annotate(LAB.get(m, m), (r.ol_err_km, r.log_miss), fontsize=7, xytext=(3, 3), textcoords="offset points")
        for m, d in tab.groupby("model"):
            ax.scatter(d.ol_err_km, d.log_miss, color=color(m), s=8, alpha=0.4)
    ax.set_xlabel(f"open-loop median footprint error at {mech.get('ttg_s', 300):.0f} s to go (km, test split)")
    ax.set_ylabel("closed-loop mean log miss (A256, 4.5 g)")
    if mech.get("rho") is not None:
        ax.set_title(f"Spearman ρ = {mech['rho']:+.2f} [{mech['ci_lo']:+.2f}, {mech['ci_hi']:+.2f}]", fontsize=9)
    figs["F04"] = _fig(L, caps, plt, fig, "F04_open_vs_closed_loop",
                       "J-mech: open-loop footprint error vs closed-loop mean log miss for R0, v1 and the four E3 models "
                       "(large: mean over seeds; small: each seed), with Spearman's ρ and its bootstrap interval over seeds.",
                       "T06b, flights")

    # F05: RL learning curves
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    curve = main.iloc[0:0]
    ecur = df[df.stage == "E4_curve"]
    for arm in list(S5.RL_TIER1) + list(S5.RL_TIER2):
        for s in rl_seeds:
            d = L.policy_dir(arm, s)
            lc = d / "learning_curve.csv"
            if lc.exists():
                c = pd.read_csv(lc)
                if len(c) and "success_5km_fast" in c:
                    x = c["decisions"] if c["decisions"].notna().any() else c.index
                    axs[0].plot(x, 100 * c["success_5km_fast"], color=color(arm), lw=1, alpha=0.8,
                                label=LAB.get(arm, arm) if s == rl_seeds[0] else None)
            e = ecur[ecur.arm.str.startswith(f"{arm}@") & (ecur.seed == s)]
            if len(e):
                pts = []
                for point, g in e.groupby(e.arm.str.split("@").str[1]):
                    x = _snapshot_x(d, point)
                    pts.append((x, 100 * g.success_5km.mean()))
                pts.sort()
                axs[1].plot([p[0] for p in pts], [p[1] for p in pts], "-o", ms=3, color=color(arm), lw=1,
                            label=LAB.get(arm, arm) if s == rl_seeds[0] else None)
    axs[0].set_title("fast environment (training distribution)", fontsize=9)
    axs[1].set_title(f"full simulator, DEV48 nominal ({cfg.rl_curve_cases} cases), snapshots", fontsize=9)
    for a_ in axs:
        a_.set_xlabel("decisions (DAgger: iteration)")
        a_.set_ylabel("success@5km (%)")
        a_.legend(fontsize=7, frameon=False)
    figs["F05"] = _fig(L, caps, plt, fig, "F05_rl_learning_curves",
                       "RL learning curves: success@5km in the fast environment during training (left; one line per seed) "
                       "and full-simulator success of the saved snapshots on the DEV48 learning-curve cases (right).",
                       "02_policies/<arm>/seed_<k>/learning_curve.csv, E4_curve flights")

    # F06-F08: ORB-1
    od = main[main.suite == "ORB1"]
    tel = orb1_telemetry(L, od) if len(od) else pd.DataFrame(columns=["flight_id"])
    ot = od.merge(tel, on="flight_id", how="left")
    okeys = [k for k in ["classical", "fnpeg", "predguid", "mppi_filtered", "v2single", "mppi_truth"] +
             list(M5.ARMS_E3) + list(S5.RL_TIER1) if (od.key == k).any()]
    n = max(len(okeys), 1)
    cols = min(n, 3)
    fig, axs = plt.subplots(int(math.ceil(n / cols)), cols, figsize=(3.6 * cols, 3.4 * math.ceil(n / cols)),
                            squeeze=False)
    for ax, k in zip(axs.ravel(), okeys):
        d = ot[ot.key == k]
        ax.scatter(d.dr_err_km, d.cr_err_km, s=3, color=color(k), alpha=0.6)
        ax.add_patch(plt.Circle((0, 0), 10.0, fill=False, color="#d00", lw=0.8))
        lim = max(15.0, float(np.nanpercentile(np.hypot(d.dr_err_km, d.cr_err_km), 99)) if d.dr_err_km.notna().any() else 15.0)
        lim = min(lim, 200.0)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_title(f"{label(*_split(k))} (n={len(d)})", fontsize=8)
        ax.set_xlabel("downrange error (km)", fontsize=7)
        ax.set_ylabel("crossrange error (km)", fontsize=7)
    for ax in axs.ravel()[len(okeys):]:
        ax.axis("off")
    figs["F06"] = _fig(L, caps, plt, fig, "F06_orb1_landing_scatter",
                       "ORB-1 landing points relative to the target (downrange and crossrange error in the nominal entry "
                       "frame) per arm, with the 10 km requirement circle; axes cut at the 99th percentile (≤ 200 km).",
                       "flights + telemetry (05_analysis/cache/orb1_telemetry.csv)")
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for k in okeys:
        m = np.sort(od[od.key == k].miss_km.to_numpy())
        if len(m):
            ax.step(np.maximum(m, 0.01), np.arange(1, len(m) + 1) / len(m), where="post", color=color(k),
                    label=f"{label(*_split(k))} (n={len(m)})")
    ax.axvline(10.0, color="#d00", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("miss (km)")
    ax.set_ylabel("fraction of flights")
    ax.legend(fontsize=7, frameon=False)
    figs["F07"] = _fig(L, caps, plt, fig, "F07_orb1_miss_cdf",
                       "ORB-1 miss distance CDFs per arm (all blocks flown), with the 10 km requirement line.", "T08")
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    tb = np.array(B.NASA_STD_3000_DECONDITIONED)
    ax.plot(tb[:, 0], tb[:, 1], color="#d00", lw=1.4, label="NASA-STD-3000, deconditioned crew")
    for k in okeys:
        env: Dict[float, float] = {}
        for gp in ot[ot.key == k].g_points.dropna():
            for lev, dur in json.loads(gp):
                env[lev] = max(env.get(lev, 0.0), dur)
        if env:
            lv = sorted(env)
            ax.plot([env[x] for x in lv], lv, color=color(k), lw=1, label=label(*_split(k)))
    ax.set_xscale("log")
    ax.set_xlabel("time above the level (s)")
    ax.set_ylabel("g level")
    ax.legend(fontsize=7, frameon=False)
    figs["F08"] = _fig(L, caps, plt, fig, "F08_g_duration",
                       "ORB-1 g-duration envelopes (Bairstow App. D procedure, 0.2 g bins): for each g level, the longest "
                       "time any flight of the arm spent above it, against the NASA-STD-3000 deconditioned-crew limit.",
                       "telemetry (05_analysis/cache/orb1_telemetry.csv)")

    # F09: success vs threshold
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    for ax, suite in zip(axs, ("A256", "ORB1")):
        d = t12[t12.suite == suite]
        for k in [k for k in order + ["mppi_truth", "v2", "v1"] if (d.arm == k).any()]:
            x = d[d.arm == k]
            ax.plot(x.threshold_km, x["success_%"], color=color(k), lw=1, label=label(*_split(k)))
        ax.axvline(5.0, color="#999", lw=0.6, ls=":")
        ax.axvline(10.0, color="#999", lw=0.6, ls=":")
        ax.set_title(suite, fontsize=9)
        ax.set_xlabel("threshold (km)")
        ax.set_ylabel("success (%)")
        ax.legend(fontsize=6, frameon=False, ncol=2)
    figs["F09"] = _fig(L, caps, plt, fig, "F09_success_vs_threshold",
                       "Success (landed within the threshold) vs threshold, A256 and ORB-1, per arm.", "T12")

    # F10: example trajectories and bank histories
    nom = a256[(a256.atmosphere == "nominal")]
    ex_c = [c for c in cases if all(((nom.key == k) & (nom.case == c) & (nom.seed.isin([-1, seeds[0]]))).any()
                                    for k in ("fnpeg", "predguid", "v2single"))][:3]
    fig, axs = plt.subplots(max(len(ex_c), 1), 2, figsize=(10, 2.8 * max(len(ex_c), 1)), squeeze=False)
    for i, c in enumerate(ex_c):
        for k in ("fnpeg", "predguid", "v2single"):
            r = nom[(nom.key == k) & (nom.case == c) & (nom.seed.isin([-1, seeds[0]]))].iloc[0]
            lay = L3 if r.source == "v3" else L
            try:
                t = S.load_telemetry(lay, r.flight_id)
            except (OSError, ValueError, KeyError):
                continue
            axs[i, 0].plot(np.asarray(t["V_mps"]) / 1000, np.asarray(t["h_m"]) / 1000, color=color(k), lw=1,
                           label=f"{label(*_split(k))}: {r.miss_km:.1f} km")
            axs[i, 1].plot(t["t_s"], t["sigma_deg"], color=color(k), lw=0.9)
        axs[i, 0].set_title(f"A256 case {c}, nominal air", fontsize=8)
        axs[i, 0].set_xlabel("speed (km/s)", fontsize=7)
        axs[i, 0].set_ylabel("altitude (km)", fontsize=7)
        axs[i, 0].legend(fontsize=6.5, frameon=False)
        axs[i, 1].set_xlabel("time (s)", fontsize=7)
        axs[i, 1].set_ylabel("bank angle (deg)", fontsize=7)
    figs["F10"] = _fig(L, caps, plt, fig, "F10_example_trajectories",
                       "Example flights on the same three A256 conditions (nominal air): altitude vs speed (left) and the "
                       "realized bank angle (right) for FNPEG, PredGuid and v2 single member (seed "
                       f"{seeds[0]}); the legend gives each miss.", "telemetry (v5; v3 for v2single)")
    return figs


def _split(k: str) -> Tuple[str, str]:
    return (k.split("@")[0], k.split("@")[1]) if "@" in k else (k, G)


def _snapshot_x(d: Path, point: str) -> float:
    """Decisions (PPO, SAC) or iteration (DAgger) of a snapshot, from its saved extra."""
    p = d / "snapshots" / f"policy_{point}.pt"
    try:
        import torch
        ex = torch.load(p, map_location="cpu", weights_only=False).get("extra", {})
        return float(ex.get("decisions", ex.get("iteration", ex.get("point", np.nan))))
    except Exception:
        m = re.search(r"(\d+)", point)
        return float(m.group(1)) if m else np.nan
