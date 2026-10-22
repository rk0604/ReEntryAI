"""
jepa/results_pdf_v3.py

Everything the four v3 notebooks saved to the study folder (jepa/v3_plan.md
section 9), in one PDF:

    0  what exists     files per stage and what is finished
    1  settings        00_study/study_config.json
    2  test set        01_testset/                         (notebook 0)
    3  training        02_models/<v1|v2>/seed_<k>/         (notebook 1): logs, curves, checkpoint notes
    4  flights         04_flights/flights.csv              (notebook 2): counts, outcomes, first look
    5  closed loop     05_analysis/                        (notebook 3): criteria, every table and figure
    6  open loop       03_open_loop/, T09, T10, F30, F31   (notebook 3)
    A  appendix        the long tables, anything not shown above, every caption, the folder index

It only reads files and recomputes nothing, so it takes minutes on a CPU
runtime. With a partly finished study it shows what exists and names what is
missing. Captions and data sources come from 05_analysis/captions.json.

    python jepa/results_pdf_v3.py <study folder> [--out <pdf>]

Default output: <study folder>/06_report/SpaceJEPA_v3_full_results.pdf
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional
from xml.sax.saxutils import escape

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

import numpy as np
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import CondPageBreak, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

from jepa import analyze_v3 as A
from jepa import make_report_v2 as R2
from jepa import study_v3 as S

TITLE = "SpaceJEPA v3: full results"
OUT_NAME = "SpaceJEPA_v3_full_results.pdf"
MAX_FIG_H = 7.4 * inch      # the letter frame is 9.5 in tall: room for the caption and a heading
MAX_COLS = 10               # wider tables are split into column slices
# Identifier columns, repeated on every slice of a split table (at most 3, in the table's own order).
KEY_COLS = ("id", "controller", "model", "kind", "arm", "a", "b", "metric", "load", "atmosphere", "set", "case",
            "seed", "stage", "ttg_s")
T01_COLS = ["controller", "atmosphere", "n_seeds", "success5_%", "success5_%_min", "success5_%_max", "miss_median_km",
            "miss_p90_km", "peak_g_mean", "peak_qdot_MWm2_mean", "heat_load_MJm2_mean", "rcs_fuel_kg_mean"]
CASE_COLS = ["case", "set", "altitude_m", "velocity_mps", "gamma_deg", "chi_deg", "phi_deg", "lambda_deg",
             "downrange_km", "crossrange_km", "approach_deg"]


# ============================================================================
# Formatting
# ============================================================================
def _fmt(v) -> str:
    if v is None:
        return "–"
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else "no"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        v = float(v)
        if not np.isfinite(v):
            return "–"
        if v.is_integer() and abs(v) < 1e7:
            return str(int(v))
        a = abs(v)
        if a >= 100:
            return f"{v:,.1f}"
        if a >= 1:
            return f"{v:.2f}"
        if a >= 0.001:
            return f"{v:.3f}"
        return f"{v:.2e}"
    s = str(v)
    return "–" if s in ("nan", "None") else s


def _size(n: float) -> str:
    return f"{n / 1e6:.1f} MB" if n >= 1e5 else f"{n / 1e3:.1f} kB"


def _widths(df: pd.DataFrame, small: bool) -> List[float]:
    """
    Column widths: each column's one-line width (longest cell, or longest header word, since headers wrap at
    underscores), then, if the table is wider than the page, only the widest columns are narrowed (to a common
    cap) and wrap; the short ones keep their width.
    """
    k = (6.9 if small else 8.2) / 72 * inch      # font size in points
    w = []
    for c in df.columns:
        body = max((len(_fmt(v)) for v in df[c].head(400)), default=1)
        head = max(len(p) for p in str(c).split("_"))
        w.append(max(0.62 * k * min(body, 70), 0.72 * k * head) + 9)
    if sum(w) <= R2.TEXT_W:
        return w
    lo, hi = 0.0, max(w)
    for _ in range(40):                          # the cap at which the narrowed widths fill the text width
        cap = (lo + hi) / 2
        lo, hi = (cap, hi) if sum(min(x, cap) for x in w) < R2.TEXT_W else (lo, cap)
    return [min(x, lo) for x in w]


def _csv(path) -> Optional[pd.DataFrame]:
    try:
        return pd.read_csv(path)
    except (OSError, pd.errors.EmptyDataError):
        return None


def _json(path) -> dict:
    try:
        return json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def _truthy(s: pd.Series) -> pd.Series:
    """Booleans read back from CSV (True/False, or strings when the column also has blanks)."""
    return s.map(lambda v: str(v).strip().lower() in ("true", "1", "1.0", "yes"))


def _find(folder: Path, prefix: str, ext: str) -> Optional[Path]:
    """The first file in `folder` named <prefix>_*.<ext> (e.g. F10_*.png)."""
    hits = sorted(folder.glob(f"{prefix}_*.{ext}")) if folder.exists() else []
    return hits[0] if hits else None


class Doc(R2.Doc):
    """make_report_v2.Doc with figures that always fit the page and tables straight from DataFrames."""

    def __init__(self):
        super().__init__()
        self.s["cell_s"] = ParagraphStyle("cell_s", parent=self.s["cell"], fontSize=6.9, leading=8.6)
        self.s["cellb_s"] = ParagraphStyle("cellb_s", parent=self.s["cellb"], fontSize=6.9, leading=8.6)
        self.s["tcap"] = ParagraphStyle("tcap", parent=self.s["caption"], spaceBefore=4, spaceAfter=4)
        self.s["missing"] = ParagraphStyle("missing", parent=self.s["small"], textColor=R2.ORANGE)
        self.table_n = 0
        self.shown: set = set()        # files already in the PDF; the appendix adds every other one

    def fig(self, path, caption, kind="real", width=R2.TEXT_W):
        """A figure no larger than it was drawn (a small plot is not blown up) and never taller than MAX_FIG_H."""
        from PIL import Image as PILImage
        with PILImage.open(path) as im:
            w, h = im.size
            dpi = float((im.info.get("dpi") or (150,))[0]) or 150.0
        self.shown.add(Path(path))
        return super().fig(str(path), caption, kind, width=min(width, w / dpi * inch, MAX_FIG_H * w / h))

    def note(self, text):
        self.add(Paragraph(text, self.s["small"]), Spacer(1, 4))

    def missing(self, what, path, notebook):
        self.add(Paragraph(f"Not found yet: {escape(what)} (<font name='DVM'>{escape(str(path))}</font>), written by "
                           f"notebook {notebook}.", self.s["missing"]), Spacer(1, 4))

    def df(self, df: Optional[pd.DataFrame], caption: str = "", max_cols: int = MAX_COLS) -> None:
        """A DataFrame as a table; a wide one becomes column slices that each repeat its key columns."""
        if df is None or not len(df):
            self.note("(empty table)")
            return
        self.table_n += 1
        df = df.reset_index(drop=True)
        keys = [c for c in df.columns if c in KEY_COLS][:3]
        rest = [c for c in df.columns if c not in keys]
        step = max(1, max_cols - len(keys))
        parts = [rest[i:i + step] for i in range(0, len(rest), step)] or [[]]
        for j, cols in enumerate(parts):
            # A long table may start low on a page and split; it just needs room for its caption and a few rows
            # (a heading right before it moves with it).
            head = []
            if self.story and isinstance(self.story[-1], Paragraph) and self.story[-1].style.name in ("h1", "h2"):
                head = [self.story.pop()]
            part = f" (part {j + 1} of {len(parts)})" if len(parts) > 1 else ""
            cap = f"<font name='DVB'>Table {self.table_n}{part}.</font> {caption if j == 0 else ''}"
            self.add(CondPageBreak(1.8 * inch), *head, Paragraph(cap, self.s["tcap"]))
            self._table(df[keys + cols])

    def _table(self, df: pd.DataFrame) -> None:
        small = len(df.columns) > 6 or len(df) > 40
        sh, sc = (self.s["cellb_s"], self.s["cell_s"]) if small else (self.s["cellb"], self.s["cell"])
        head = [Paragraph(escape(str(c)).replace("_", " "), sh) for c in df.columns]
        body = [[Paragraph(escape(_fmt(v)), sc) for v in row] for row in df.itertuples(index=False)]
        t = Table([head] + body, colWidths=_widths(df, small), repeatRows=1, hAlign="LEFT")
        pad = 2.2 if small else 3.2
        st = [("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 0), (-1, 0), 0.8, R2.INK),
              ("LINEBELOW", (0, -1), (-1, -1), 0.5, R2.RULE), ("TOPPADDING", (0, 0), (-1, -1), pad),
              ("BOTTOMPADDING", (0, 0), (-1, -1), pad), ("LEFTPADDING", (0, 0), (-1, -1), 3),
              ("RIGHTPADDING", (0, 0), (-1, -1), 3)]
        st += [("BACKGROUND", (0, i), (-1, i), colors.HexColor("#f6f5f2")) for i in range(2, len(df) + 1, 2)]
        t.setStyle(TableStyle(st))
        self.add(t, Spacer(1, 8))

    def table_file(self, path: Path, caps: dict, df: Optional[pd.DataFrame] = None) -> None:
        """A saved T## table with its caption and data source from captions.json."""
        self.shown.add(Path(path))
        self.df(df if df is not None else _csv(path), _caption(caps, Path(path).stem, Path(path).name))

    def fig_file(self, path: Path, caps: dict, fallback: str = "") -> None:
        self.fig(path, _caption(caps, Path(path).stem, fallback or Path(path).name))


def _caption(caps: dict, name: str, fallback: str = "") -> str:
    """'<name>: caption. Data: source.' from captions.json, escaped for reportlab."""
    c = caps.get(name, {})
    text = escape(c.get("caption") or fallback)
    data = f" <i>Data: {escape(c['data'])}.</i>" if c.get("data") else ""
    return f"<font name='DVM'>{escape(name)}</font>: {text}{data}"


# ============================================================================
# Sections
# ============================================================================
def _inventory(root: Path) -> pd.DataFrame:
    """Every file in the study folder: path (relative), stage, size."""
    rows = []
    for dirpath, _, names in os.walk(root):
        for n in names:
            if n.startswith("."):
                continue
            p = Path(dirpath) / n
            try:
                size = p.stat().st_size
            except OSError:
                continue
            rel = p.relative_to(root).as_posix()
            rows.append({"path": rel, "stage": rel.split("/")[0], "size": size})
    return pd.DataFrame(rows, columns=["path", "stage", "size"])


def cover(d: Doc, L: S.Layout, info: dict) -> None:
    cfg = info["config"]
    t01, t05 = info["T01"], info["T05"]
    seeds = cfg.get("seeds") or []
    d.add(Spacer(1, 0.8 * inch), Paragraph(TITLE, d.s["title"]),
          Paragraph(f"v1 vs v2 vs the classical predictor-corrector · {len(seeds)} training seeds · "
                    f"{info['n_cases']} entry cases × {len(A.ATM)} atmospheres · load limit {cfg.get('g_primary', '?')} g "
                    f"· goal: land within {cfg.get('success_km', '?')} km", d.s["subtitle"]))
    lines = []
    if t01 is not None and "success5_%" in t01:
        for t in A.ATM:
            parts = []
            for c in A.CTRL:
                r = t01[(t01.controller == c) & (t01.atmosphere == t)]
                if not len(r):
                    continue
                r = r.iloc[0]
                rng = (f" ({r['success5_%_min']:.1f}–{r['success5_%_max']:.1f})"
                       if c in S.LEARNED and r.get("n_seeds", 1) > 1 and "success5_%_min" in r else "")
                parts.append(f"{A.LAB.get(c, c)} {r['success5_%']:.1f}%{rng}")
            lines.append(f"<b>{A.ATM_LAB[t]}</b>: " + ", ".join(parts))
    else:
        lines.append("T01 is not in 05_analysis/tables yet (notebook 3).")
    if t05 is not None and "met" in t05:
        met = _truthy(t05.met)
        is_c = t05.id.astype(str).str.startswith("C")
        unmet = ", ".join(t05.id[is_c & ~met].astype(str)) or "none"
        bridge = met[t05.id == "Bridge"]
        lines.append(f"<b>Pre-registered criteria met: {int((is_c & met).sum())} of {int(is_c.sum())}</b> "
                     f"(not met: {unmet}). Bridge check: "
                     f"{'met' if len(bridge) and bool(bridge.iloc[0]) else 'not met or not run'}.")
    d.callout("Success within 5 km at the load limit (learned controllers: mean over seeds, seed range in brackets)",
              lines, R2.TINT_YELLOW, colors.HexColor("#eda100"))
    d.add(Paragraph("<br/>".join([
        "0. What is in the study folder", "1. Settings", "2. Test set (notebook 0)",
        "3. Training: every model and seed (notebook 1)", "4. Flights (notebook 2)",
        "5. Closed loop: criteria, tables, figures (notebook 3)", "6. Open loop, per seed (notebook 3)",
        "Appendix: long tables, everything else, every caption, folder index"]), d.s["toc"]), Spacer(1, 14))
    d.note(f"Study folder: <font name='DVM'>{escape(str(L.root))}</font><br/>"
           f"Git commit of the study: {escape(str(info['study'].get('git_commit', '?'))[:10])} · study created "
           f"{escape(str(info['study'].get('created', '?')))} · this PDF made {time.strftime('%Y-%m-%d %H:%M')}<br/>"
           "Every number and figure is read from the files in the study folder; nothing is typed in or recomputed.")


def section_status(d: Doc, L: S.Layout, info: dict) -> None:
    d.add(PageBreak())
    d.h1("0. What is in the study folder")
    inv, cfg = info["inventory"], info["config"]
    seeds = cfg.get("seeds") or []
    n_ck = {k: sum(L.ckpt(k, s).exists() for s in seeds) for k in ("v1", "v2")}
    n_fp = len(list(L.open_loop.glob("footprint_*_seed*.csv"))) if L.open_loop.exists() else 0
    flights = info["flights"]
    status = {
        "00_study": "settings frozen" if (L.study / "study_config.json").exists() else "no study_config.json",
        "01_testset": f"{info['n_cases']} cases" if (L.testset / "v3_testset.csv").exists() else "not frozen yet",
        "02_models": f"v2: {n_ck['v2']} of {len(seeds)} seeds, v1: {n_ck['v1']} of {len(seeds)} seeds",
        "03_open_loop": f"{n_fp} of {2 * len(seeds)} model × seed evaluations",
        "04_flights": f"{len(flights)} flights in flights.csv" if flights is not None else "no flights.csv yet",
        "05_analysis": (f"{len(list(L.tables.glob('T*.csv')))} tables, {len(list(L.figures.glob('F*.png')))} figures"
                        if L.analysis.exists() else "not analyzed yet"),
        "06_report": ", ".join(p.name for p in sorted(L.report.glob("*.pdf"))) if L.report.exists() else "",
    }
    rows = []
    for stage, what in S.STAGES.items():
        sub = inv[inv.stage == stage]
        rows.append({"folder": stage, "what": what, "files": len(sub), "size": _size(sub["size"].sum()),
                     "status": status.get(stage, "")})
    d.df(pd.DataFrame(rows), "The study folder, stage by stage (jepa/v3_plan.md §9).")


def section_settings(d: Doc, L: S.Layout, info: dict) -> None:
    d.h1("1. Settings")
    study = info["study"]
    if not study:
        d.missing("the frozen settings", L.study / "study_config.json", "0, 1, 2 or 3 (any of them)")
        return
    d.p(f"Frozen on the first run of any v3 notebook ({escape(str(study.get('created', '?')))}, git "
        f"{escape(str(study.get('git_commit', '?'))[:10])}); every later session must use the same settings.")
    rows = [{"setting": k, "value": v if isinstance(v, (int, float, str, bool)) or v is None
             else json.dumps(v, separators=(", ", ": "))[:400]} for k, v in info["config"].items()]
    d.df(pd.DataFrame(rows), "Every setting of the study (00_study/study_config.json, \"config\").", max_cols=2)
    env = study.get("environment") or {}
    if env:
        d.df(pd.DataFrame([{"item": k, "value": json.dumps(v) if isinstance(v, (dict, list)) else v}
                           for k, v in env.items()]),
             "The environment that froze the settings (00_study/study_config.json, \"environment\").", max_cols=2)


def section_testset(d: Doc, L: S.Layout, info: dict) -> None:
    d.add(PageBreak())
    d.h1("2. Test set (notebook 0)")
    cases, caps = info["cases"], info["caps"]
    if cases is None:
        d.missing("the frozen test set", L.testset / "v3_testset.csv", 0)
        return
    spec = _json(L.testset / "testset_spec.json")
    n_new = int((cases.set == "new").sum()) if "set" in cases else len(cases)
    d.p(f"{len(cases)} entry cases: {len(cases) - n_new} bridge cases (the frozen v1/v2 entries) and {n_new} new "
        f"cases" + (f" from {spec['n_candidates']} Sobol candidates (rejected: {escape(json.dumps(spec.get('rejections')))})"
                    if spec.get("n_candidates") else "") +
        (f". MPPI planning horizon {spec['horizon_steps']} steps ({spec['horizon_steps'] * 5} s)."
         if spec.get("horizon_steps") else "."))
    feas_cols = [c for c in cases.columns if c.startswith("feasible_")]
    if feas_cols and "set" in cases:
        feas = cases.groupby("set")[feas_cols].mean().mul(100).round(1).reset_index()
        feas.columns = ["set"] + [A.ATM_LAB.get(c[len("feasible_"):], c) + " %" for c in feas_cols]
        d.df(feas, "Share of cases whose target is inside the conservative 4.5 g footprint, per atmosphere.")
        bad = [{"atmosphere": A.ATM_LAB.get(c[len("feasible_"):], c),
                "infeasible cases": ", ".join(str(int(x)) for x in cases.case[~cases[c].astype(bool)]) or "none"}
               for c in feas_cols]
        d.df(pd.DataFrame(bad), "Cases flagged infeasible at 4.5 g. They stay in the primary results; T02 leaves them out.")
    ic = [c for c in CASE_COLS[2:] if c in cases]
    if ic:
        new = cases[cases.set == "new"] if "set" in cases else cases
        rng = new[ic].describe().loc[["min", "50%", "max"]].T.reset_index().rename(
            columns={"index": "quantity", "50%": "median"})
        d.df(rng, "Range of the new cases' entry conditions and target placement.")
    for name in ("F01", "F02"):
        p = _find(L.testset, name, "png")
        if p:
            d.fig_file(p, caps)
        else:
            d.missing(f"figure {name}", L.testset, "0 or 3")
    d.shown.update({L.testset / "v3_testset.csv"})


def _ckpt_extra(path: Path) -> dict:
    try:
        import torch
        return torch.load(path, map_location="cpu", weights_only=False).get("extra") or {}
    except Exception as e:      # an unreadable checkpoint should not stop the PDF
        return {"_error": f"{type(e).__name__}: {e}"[:120]}


def section_training(d: Doc, L: S.Layout, info: dict) -> None:
    d.add(PageBreak())
    d.h1("3. Training: every model and seed (notebook 1)")
    seeds = info["config"].get("seeds") or []
    for kind, label in (("v2", "v2 ensemble"), ("v1", "v1 (fair budget)")):
        d.h2(label)
        notes, finals, health = [], [], []
        for s in seeds:
            md = L.model_dir(kind, s)
            row = {"seed": s, "checkpoint": "done" if L.ckpt(kind, s).exists() else "not trained yet"}
            if L.ckpt(kind, s).exists():
                ex = _ckpt_extra(L.ckpt(kind, s))
                row["train_min"] = round(ex["train_s"] / 60, 1) if ex.get("train_s") else None
                row["train_episodes"] = ex.get("n_train_episodes")
                if kind == "v1":
                    ext = ex.get("extended") or {}
                    row["extended (wm, prober)"] = ", ".join(
                        "yes" if (v or {}).get("extend") else "no" for v in (ext.get("wm"), ext.get("prober")))
                if ex.get("_error"):
                    row["note"] = ex["_error"]
            notes.append(row)
            log = _csv(md / "train_log.csv")
            if log is not None and len(log):
                d.shown.add(md / "train_log.csv")
                for stg, g in log.groupby("stage", sort=False):
                    last = g.iloc[-1]
                    finals.append({"seed": s, "stage": stg,
                                   **{c: last[c] for c in g.columns if c != "stage" and pd.notna(last[c])}})
            h = _csv(md / "latent_health.csv")
            if h is not None and len(h):
                d.shown.add(md / "latent_health.csv")
                health.append(h.assign(seed=s))
        d.df(pd.DataFrame(notes), f"{label}: checkpoints (02_models/{kind}/seed_&lt;k&gt;/ckpt.pt).")
        if finals:
            d.df(pd.DataFrame(finals), f"{label}: the last logged values of each training stage "
                                       f"(02_models/{kind}/seed_&lt;k&gt;/train_log.csv).")
        if health:
            hh = pd.concat(health, ignore_index=True)
            d.df(hh[["seed"] + [c for c in hh.columns if c != "seed"]],
                 f"{label}: latent health on the validation split (02_models/{kind}/seed_&lt;k&gt;/latent_health.csv).")
        for s in seeds:
            p = L.model_dir(kind, s) / "F_curves.png"
            if p.exists():
                d.fig(p, f"<font name='DVM'>02_models/{kind}/seed_{s}/F_curves.png</font>: training curves of {label}, "
                         f"seed {s}: every logged loss on a log scale, world model then prober. "
                         f"<i>Data: train_log.csv in the same folder.</i>")
            else:
                d.missing(f"{label} seed {s} training curves", p, 1)


def section_flights(d: Doc, L: S.Layout, info: dict) -> None:
    d.add(PageBreak())
    d.h1("4. Flights (notebook 2)")
    df = info["flights"]
    if df is None or not len(df):
        d.missing("the merged flight table", L.flights / "flights.csv", 2)
        return
    d.shown.add(L.flights / "flights.csv")
    n_cond = info["n_cases"] * len(A.ATM)
    df = df.copy()
    df["load"] = df["load"].fillna("own guidance").astype(str) if "load" in df else "–"
    df["mppi_rep"] = df["mppi_rep"].fillna(0).astype(int) if "mppi_rep" in df else 0
    arms = []
    for (c, ld, rep), g in df.groupby(["controller", "load", "mppi_rep"], sort=False):
        seeds = sorted(int(s) for s in g.seed.dropna().unique()) if "seed" in g else []
        n_seed_runs = max(1, len([s for s in seeds if s >= 0]))
        arms.append({"controller": A.LAB.get(c, c), "load": ld, "mppi_rep": rep,
                     "seeds": ", ".join(str(s) for s in seeds if s >= 0) or "deterministic",
                     "flights": len(g), "of": n_seed_runs * n_cond,
                     "success5_%": 100 * g.success_5km.mean(), "miss_median_km": g.miss_km.median()})
    d.p(f"{len(df)} flights. Each is one JSON result and one 1 Hz telemetry file in 04_flights/; flights.csv "
        "merges the results. 'of' is the number of flights the arm has when every seed it has run is complete "
        f"({info['n_cases']} cases × {len(A.ATM)} atmospheres per seed).")
    d.df(pd.DataFrame(arms), "Flights per controller arm (all atmospheres together).")
    if "outcome" in df:
        out = pd.crosstab([df.controller.map(lambda c: A.LAB.get(c, c)), df["load"]], df.outcome).reset_index()
        d.df(out, "How the flights ended (drogue = reached drogue deploy at 7.6 km; anything else is a failure).")
    first = df[df.mppi_rep == 0].assign(arm=lambda x: x.controller.map(lambda c: A.LAB.get(c, c)) + " · " + x["load"])
    agg = {"n": ("miss_km", "size"), "success5_%": ("success_5km", lambda s: 100 * s.mean()),
           "median_km": ("miss_km", "median"), "p90_km": ("miss_km", lambda s: s.quantile(0.9))}
    for col, name in (("peak_g", "peak_g"), ("rcs_fuel_kg", "fuel_kg")):
        if col in first:
            agg[name] = (col, "mean")
    look = first.groupby(["arm", "atmosphere"], sort=False).agg(**agg).reset_index()
    look["atmosphere"] = look.atmosphere.map(lambda t: A.ATM_LAB.get(t, t))
    d.df(look, "First look per arm and atmosphere, pooled over seeds (the table notebook 2 prints). "
               "The statistics are in section 5.")
    tcols = {"t_s": "flight_s", "n_reversals": "bank_reversals", "plan_ms_mean": "plan_ms_mean",
             "plan_ms_max": "plan_ms_max", "wall_s": "wall_s"}
    tcols = {k: v for k, v in tcols.items() if k in df}
    if tcols:
        tim = df.groupby("controller")[list(tcols)].mean().rename(columns=tcols).reset_index()
        tim["controller"] = tim.controller.map(lambda c: A.LAB.get(c, c))
        d.df(tim, "Flight time, bank reversals and planning time, mean per flight. The classical controller's "
                  "planning time is its whole 1 s guidance step including the physics, so it is not comparable "
                  "with MPPI's per-decision time.")


def section_closed_loop(d: Doc, L: S.Layout, info: dict) -> None:
    d.add(PageBreak())
    d.h1("5. Closed loop: criteria, tables, figures (notebook 3)")
    caps, t05, t01 = info["caps"], info["T05"], info["T01"]
    if t05 is None and t01 is None:
        d.missing("the closed-loop analysis", L.tables, 3)
        return
    if t05 is not None:
        d.h2("Pre-registered criteria (jepa/v3_plan.md §8.2, fixed before any v3 flight)")
        d.table_file(L.tables / "T05_criteria.csv", caps,
                     t05[[c for c in ("id", "question", "measured", "p_holm", "met") if c in t05]])
    d.h2("The main metrics at the load limit")
    for name in ("T01_metrics_primary", "T02_metrics_feasible_only"):     # ~60 columns each: the full ones are in A1
        t = _csv(L.tables / f"{name}.csv")
        if t is not None:
            d.df(t[[c for c in T01_COLS if c in t]], _caption(caps, name) +
                 " Selected columns; the full table is in the appendix.", max_cols=len(T01_COLS))
        else:
            d.missing(f"table {name[:3]}", L.tables, 3)
    groups = [("Success and miss", ["F10", "F11", "F12"], ["T03", "T04"]),
              ("Load, heating, and the cost of load matching", ["F13", "F14"], ["T06"]),
              ("Which cases are hard, and seed stability", ["F15", "F16", "F17"], ["T07"]),
              ("RCS usage", ["F20", "F21"], ["T08"])]
    for title, figs, tabs in groups:
        d.h2(title)
        for k in figs:
            p = _find(L.figures, k, "png")
            d.fig_file(p, caps) if p else d.missing(f"figure {k}", L.figures, 3)
        for k in tabs:
            p = _find(L.tables, k, "csv")
            d.table_file(p, caps) if p else d.missing(f"table {k}", L.tables, 3)


def section_open_loop(d: Doc, L: S.Layout, info: dict) -> None:
    d.add(PageBreak())
    d.h1("6. Open loop, per seed (notebook 3)")
    caps = info["caps"]
    d.p("How well each trained model predicts the landing point on the 600 held-out test flights, with no "
        "controller involved, from 400, 300, 200, 100 and 50 s before landing.")
    for k in ("F30", "F31"):
        p = _find(L.figures, k, "png")
        d.fig_file(p, caps) if p else d.missing(f"figure {k}", L.figures, 3)
    for k in ("T09", "T10"):
        p = _find(L.tables, k, "csv")
        d.table_file(p, caps) if p else d.missing(f"table {k}", L.tables, 3)
    along = info["along"]
    if along is not None and "along_mean_km" in along:
        piv = along.pivot_table(index=["model", "ttg_s"], columns="density_bin", values="along_mean_km").reset_index()
        piv.columns = [str(c) if c in ("model", "ttg_s") else f"{c} km" for c in piv.columns]
        d.df(piv, "Mean along-track error [km] by air density (+ = predicted landing beyond the true one), "
                  "averaged over seeds. Per seed: appendix. <i>Data: 03_open_loop/along_track_*_seed*.csv.</i>")


def appendix(d: Doc, L: S.Layout, info: dict) -> None:
    d.add(PageBreak())
    d.h1("Appendix")
    caps = info["caps"]
    d.h2("A1. The long tables")
    for k in ("T01", "T01b", "T02", "T11"):
        p = _find(L.tables, k, "csv")
        if p:
            d.table_file(p, caps)
    cases = info["cases"]
    if cases is not None:
        d.h2("A2. The 64 test cases")
        cols = [c for c in CASE_COLS if c in cases] + [c for c in cases.columns if c.startswith("feasible_")]
        d.df(cases[cols], "The frozen test set (01_testset/v3_testset.csv; footprint margins per atmosphere are in "
                          "01_testset/footprints.csv).")
    along = info["along"]
    if along is not None:
        d.h2("A3. Open-loop along-track bias, every seed")
        d.df(along, "03_open_loop/along_track_&lt;model&gt;_seed&lt;k&gt;.csv, all files together.")
        d.shown.update(L.open_loop.glob("along_track_*_seed*.csv"))
    inv = info["inventory"]
    rest_t = [p for p in sorted(L.tables.glob("*.csv")) if p not in d.shown] if L.tables.exists() else []
    rest_f = [p for p in (L.root / r for r in sorted(inv.path) if r.endswith(".png")) if p not in d.shown]
    if rest_t or rest_f:
        d.h2("A4. Every other table and figure in the study folder")
        for p in rest_t:
            d.table_file(p, caps)
        for p in rest_f:
            d.fig(p, _caption(caps, p.stem, str(p.relative_to(L.root).as_posix())))
    if caps:
        d.h2("A5. Every caption and data source (05_analysis/captions.json)")
        d.df(pd.DataFrame([{"id": k, "caption": v.get("caption", ""), "data": v.get("data", "")}
                           for k, v in sorted(caps.items())]), "captions.json", max_cols=3)
    if len(inv):
        d.h2("A6. Folder index")
        parts = inv.path.str.split("/")
        inv = inv.assign(group=[("/".join(p[:3]) + "/…") if p[0] == "04_flights" and len(p) > 3 else "/".join(p)
                                for p in parts])
        idx = inv.groupby("group", sort=True).agg(files=("path", "size"), size=("size", lambda s: _size(s.sum())))
        d.df(idx.reset_index().rename(columns={"group": "path"}),
             "Every file in the study folder (the per-flight results and telemetry grouped by controller).")


# ============================================================================
def _page(canvas, doc):
    if doc.page == 1:
        return
    canvas.saveState()
    canvas.setFont("DV", 7.5)
    canvas.setFillColor(R2.MUT)
    canvas.drawString(0.75 * inch, 0.5 * inch, TITLE)
    canvas.drawRightString(letter[0] - 0.75 * inch, 0.5 * inch, str(doc.page))
    canvas.setStrokeColor(R2.RULE)
    canvas.setLineWidth(0.5)
    canvas.line(0.75 * inch, letter[1] - 0.55 * inch, letter[0] - 0.75 * inch, letter[1] - 0.55 * inch)
    canvas.restoreState()


def build(root, out_pdf: Optional[str] = None, log=print) -> str:
    """Write the full-results PDF of the study folder `root`; returns its path."""
    L = S.Layout(root)
    if not L.root.exists():
        raise FileNotFoundError(f"no study folder at {L.root}")
    t0 = time.time()
    R2._fonts()
    log("reading the study folder ...")
    study = _json(L.study / "study_config.json")
    cases = _csv(L.testset / "v3_testset.csv")
    cfg = study.get("config") or {}
    along_files = sorted(L.open_loop.glob("along_track_*_seed*.csv")) if L.open_loop.exists() else []
    info: Dict[str, object] = {
        "study": study, "config": cfg, "cases": cases,
        # cases flown per seed: the config's subset (smoke mode), else the whole frozen set
        "n_cases": len(cfg["cases"]) if cfg.get("cases") else (len(cases) if cases is not None else 64),
        "caps": _json(L.analysis / "captions.json"),
        "T01": _csv(L.tables / "T01_metrics_primary.csv"), "T05": _csv(L.tables / "T05_criteria.csv"),
        "flights": _csv(L.flights / "flights.csv"),
        "along": pd.concat([_csv(f) for f in along_files], ignore_index=True) if along_files else None,
    }
    log(f"  indexing files ({time.time() - t0:.0f}s) ...")
    info["inventory"] = _inventory(L.root)
    log(f"  {len(info['inventory'])} files ({time.time() - t0:.0f}s)")
    d = Doc()
    for name, fn in (("cover", cover), ("status", section_status), ("settings", section_settings),
                     ("test set", section_testset), ("training", section_training), ("flights", section_flights),
                     ("closed loop", section_closed_loop), ("open loop", section_open_loop), ("appendix", appendix)):
        fn(d, L, info)
        log(f"  {name} ({time.time() - t0:.0f}s)")
    out_pdf = str(out_pdf or L.report / OUT_NAME)
    Path(out_pdf).parent.mkdir(parents=True, exist_ok=True)
    tmp = out_pdf + ".tmp.pdf"
    doc = SimpleDocTemplate(tmp, pagesize=letter, leftMargin=0.75 * inch, rightMargin=0.75 * inch,
                            topMargin=0.75 * inch, bottomMargin=0.75 * inch, title=TITLE)
    doc.build(d.story, onFirstPage=_page, onLaterPages=_page)
    os.replace(tmp, out_pdf)
    log(f"wrote {out_pdf}: {d.fig_n} figures, {d.table_n} tables, {Path(out_pdf).stat().st_size / 1e6:.1f} MB "
        f"({time.time() - t0:.0f}s)")
    return out_pdf


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("root", help="the study folder (e.g. MyDrive/ReEntryAI_jepa_v3)")
    ap.add_argument("--out", default=None, help=f"output PDF (default: <root>/06_report/{OUT_NAME})")
    a = ap.parse_args()
    build(a.root, a.out)
