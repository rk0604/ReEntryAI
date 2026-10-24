"""
jepa/report_v4.py

06_report/SpaceJEPA_v4_report.pdf (jepa/v4_plan.md section 10). It reads only
the files analyze_v4 wrote (05_analysis/tables, figures, captions.json) and the
plan's section 15 (deviations); no number is typed in. Building blocks:
results_pdf_v3.Doc and make_report_v2._fonts (DejaVu from the system, else
matplotlib's copies, so it builds on Colab and Windows).

Sections: summary callout (headline numbers, criteria met), RQ1 to RQ4,
methods and definitions, deviations, and an appendix with every table and
any figure not shown above.

    python jepa/report_v4.py <v4 study folder>
"""

from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path
from typing import Optional, Tuple
from xml.sax.saxutils import escape

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

import numpy as np
import pandas as pd

TITLE = "SpaceJEPA v4: why is v2 better?"
PLAN = SIM_DIR / "jepa" / "v4_plan.md"
ATM_LAB = {"nominal": "nominal", "mismatch": "thin ×0.88", "mismatch_high": "dense ×1.12",
           "out_of_range": "×0.78 / L/D ×0.80", "all": "pooled"}


def _csv(path: Path) -> Optional[pd.DataFrame]:
    try:
        return pd.read_csv(path)
    except (OSError, pd.errors.EmptyDataError):
        return None


def _f(v, fmt="{:+.1f}") -> str:
    try:
        v = float(v)
    except (TypeError, ValueError):
        return "–"
    return fmt.format(v) if np.isfinite(v) else "–"


MAX_ROWS = 300


def compact(key: str, df: Optional[pd.DataFrame]) -> Tuple[Optional[pd.DataFrame], str]:
    """
    A long table as the report shows it (the CSV keeps every row): T12 at 1, 2, 5, 10 and 25 km, the open-loop
    tables as means over seeds, anything else cut at MAX_ROWS rows. Returns (table, note).
    """
    if df is None or len(df) <= MAX_ROWS:
        return df, ""
    if key == "T12":
        d = df[df.threshold_km.isin([1.0, 2.0, 5.0, 10.0, 25.0])]
        d = d.pivot_table(index=["arm", "load", "atmosphere"], columns="threshold_km", values="success_%").reset_index()
        d.columns = [c if isinstance(c, str) else f"success_%_{c:g}km" for c in d.columns]
        return d, f"Shown at 1, 2, 5, 10, 25 km; the CSV has all {len(df)} rows (every 0.5 km)."
    if key == "T10":
        d = df.groupby(["model", "split", "ttg_s"]).err_median_km.mean().unstack("ttg_s").reset_index()
        d.columns = [c if isinstance(c, str) else f"err_median_km_{int(c)}s" for c in d.columns]
        return d, f"Mean over seeds; the CSV has every seed ({len(df)} rows)."
    if key == "T10b":
        d = df.groupby(["model", "split", "ttg_s", "density_bin"]).along_mean_km.mean().unstack("density_bin").reset_index()
        return d, f"Mean along-track error over seeds; the CSV has every seed ({len(df)} rows)."
    return df.head(MAX_ROWS), f"First {MAX_ROWS} of {len(df)} rows; the CSV has all of them."


def deviations_text() -> str:
    """Section 15 of the plan, as written (the record of any change after the full run started)."""
    try:
        txt = PLAN.read_text(encoding="utf-8")
    except OSError:
        return "(v4_plan.md not found)"
    m = re.search(r"^## 15\. Deviations\s*$(.*)", txt, flags=re.M | re.S)
    return m.group(1).strip() if m else "(no section 15 in v4_plan.md)"


def build_report(L, out_pdf: Optional[str] = None) -> str:
    """Write the report from the files in L (a study_v4.LayoutV4); returns its path."""
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.units import inch
    from reportlab.platypus import PageBreak, Paragraph, SimpleDocTemplate
    from jepa import make_report_v2 as R2
    from jepa import results_pdf_v3 as RP

    R2._fonts()
    d = RP.Doc()
    caps = json.loads((L.analysis / "captions.json").read_text()) if (L.analysis / "captions.json").exists() else {}
    cfgdoc = json.loads((L.study / "study_config.json").read_text()) if (L.study / "study_config.json").exists() else {}
    cfg = cfgdoc.get("config", {})
    T = {p.stem.split("_")[0]: p for p in sorted(L.tables.glob("T*.csv"))}
    Fg = {p.stem.split("_")[0]: p for p in sorted(L.figures.glob("F*.png"))}
    t = {k: _csv(p) for k, p in T.items()}
    t05, t01, t04, t06, t08 = (t.get(k) for k in ("T05", "T01", "T04", "T06", "T08"))

    d.add(Paragraph(TITLE, d.s["title"]),
          Paragraph(f"Planner vs model, load matching, ablations of the v1 → v2 changes, case {cfg.get('diag_case', 10)} · "
                    f"{len(cfg.get('seeds', []))} seeds · mode {cfg.get('mode', '?')} · git "
                    f"{str(cfgdoc.get('git_commit', '?'))[:10]}", d.s["subtitle"]))

    # ---- summary callout ----
    lines = []
    if t01 is not None and len(t01):
        def s5(arm, load):
            r = t01[(t01.arm == arm) & (t01.load == load)]
            return ", ".join(f"{ATM_LAB.get(x.atmosphere, x.atmosphere)} {x['success5_%']:.1f}%" for _, x in r.iterrows()) or "–"
        for arm, load, lab in (("classical", "own", "classical"), ("mppi_nominal", "4.5g", "MPPI + nominal physics"),
                               ("v2", "4.5g", "v2 (learned model)"), ("mppi_truth", "4.5g", "MPPI + true physics"),
                               ("v2", "matched", "v2 at classical's load")):
            lines.append(f"<b>{lab}</b>, success@5km: {s5(arm, load)}")
    if t06 is not None and len(t06):
        r = t06[(t06.load == "4.5g") & (t06.metric == "success5") & (t06.atmosphere == "all")]
        if len(r):
            r = r.iloc[0]
            lines.append(f"<b>Gap classical → v2 at 4.5 g</b>: {_f(r.total)} pp = planner {_f(r.planner)} + model "
                         f"{_f(r.model)} (shares {_f(r.planner_share, '{:.2f}')} / {_f(r.model_share, '{:.2f}')})")
    if t05 is not None and len(t05):
        test = t05[t05.met.notna() & (t05.met.astype(str) != "")]
        met = test[test.met.astype(str).str.lower() == "true"]
        lines.append(f"<b>Criteria met: {len(met)} of {len(test)}</b> pass/fail criteria (estimates P1, P6, L3, M are not "
                     f"pass/fail): {', '.join(met.id.astype(str)) or 'none'}.")
    d.callout("The numbers (read from 05_analysis/tables)", lines or ["No analysis tables yet."],
              R2.TINT_YELLOW, colors.HexColor("#eda100"))
    if t05 is not None and len(t05):
        d.h2("Pre-registered criteria (v4_plan.md section 8, fixed before any v4 flight)")
        rows = [["#", "Question", "Measured", "p (Holm)", "Met"]]
        for r in t05.itertuples():
            ph = _f(r.p_holm, "{:.3g}")
            met = "estimate" if str(r.met) in ("None", "nan", "") else ("yes" if str(r.met).lower() == "true" else "<b>no</b>")
            rows.append([escape(str(r.id)), escape(str(r.question)), escape(str(r.measured)), ph, met])
        d.table(rows, [0.95, 1.7, 3.1, 0.6, 0.55])

    def fig(k):
        if k in Fg:
            d.fig_file(Fg[k], caps)

    def tab(k, df=None):
        if k in T:
            df, note = compact(k, t.get(k) if df is None else df)
            d.table_file(T[k], caps, df)
            if note:
                d.note(note)

    d.add(PageBreak())
    d.h1("RQ1. Planner or model?")
    d.p("Each arm is a planner flying with a model. Comparing arms that share one of the two isolates the other: "
        "MPPI + nominal physics − classical is the planner's contribution, v2 − MPPI + nominal physics the learned "
        "model's, MPPI + true physics − v2 what a perfect model would still add, classical knowing the true density and mass − classical "
        "the classical controller's own model error, and v2 − MPPI + v2's physics prior what the learned corrections add.")
    fig("F10"), fig("F11")
    if t04 is not None and len(t04):
        tab("T04", t04[t04.atmosphere == "all"].drop(columns=["per_seed"], errors="ignore"))
    tab("T06")
    fig("F18")
    d.h1("RQ2. Load")
    d.p("v2 flown with lower soft load limits, down to the classical controller's own peak g on each condition "
        "(matched). L1 checks the realized load; L2 asks whether v2 still beats classical there.")
    fig("F12"), tab("T07")
    d.h1("RQ3. Which v1 → v2 changes matter?")
    d.p("A variants remove one v2 change from the one-member reference R0 (is it necessary?); B variants add one to v1 "
        "(is it sufficient?). The primary metric is the paired log-miss difference; its mean is the log of the "
        "geometric-mean miss ratio.")
    fig("F13"), tab("T08"), fig("F14"), fig("F15"), tab("T10")
    d.h1(f"RQ4. Case {cfg.get('diag_case', 10)} (descriptive)")
    d.p("Diagnostic flights with the decision log: where the trajectory separates from the classical one, whether the "
        "planner's predicted miss tracked the real outcome (model error) or the chosen plan was predicted well and still "
        "missed (planner error), whether the predicted peak g sat at the soft limit, and which arms land the case. "
        "No claim of cause beyond these comparisons.")
    fig("F17")
    t11 = t.get("T11")
    if t11 is not None and len(t11):
        tab("T11", t11[t11.part == "diag"].dropna(axis=1, how="all"))
    d.h1("Success against the threshold, and the ensemble re-scored")
    fig("F16"), tab("T13")
    d.h1("Methods and definitions")
    d.bullets([
        "<b>success@5km / @2km</b>: outcome drogue (7.6 km) and great-circle miss below 5 / 2 km. <b>log miss</b> = "
        "ln(max(miss, 0.01 km)); failures keep their miss at termination.",
        "<b>Paired comparisons</b>: the unit is one (case, atmosphere) condition, paired across arms by entry state and "
        "MPPI noise stream (common random numbers; physics arms use noise index s = v2 seed s's stream). Seeds are "
        "repeated measures; deterministic arms are broadcast over seeds.",
        f"<b>Uncertainty</b>: hierarchical bootstrap, {cfg.get('n_boot')} resamples (seeds, then cases within each "
        "atmosphere), 95% percentile intervals (90% for the P3 equivalence test); bootstrap p-values; per-seed exact "
        "McNemar and sign tests; Holm within RQ1 (P2, P4, P5), RQ2 (L2), RQ3 (13 variants).",
        "<b>Physics models</b>: mppi_nominal / mppi_truth fly MPPI through a float64 torch port of the classical "
        "predictor's own model (math_3d.nominal_eom_step, the tabulated atmosphere, Orion trim aero) with its 2 s Euler "
        "step, the realized bank rate-limited at 20°/s; mppi_prior through v2's physics backbone with zero corrections.",
        "<b>Planner</b>: MPPI as in v3 (256 samples, 6 knots, 30° noise, 4 iterations at entry then 1, horizon 188 steps, "
        "soft g limit 50 km per g² above the limit).",
    ])
    d.h1("Deviations (v4_plan.md section 15)")
    for para in [p for p in deviations_text().split("\n\n") if p.strip()]:
        d.p(escape(para.strip()))
    d.add(PageBreak())
    d.h1("Appendix: every table")
    for k, p in T.items():
        if p not in d.shown:
            tab(k)
    rest = [p for k, p in Fg.items() if p not in d.shown]
    if rest:
        d.h1("Appendix: other figures")
        for p in rest:
            d.fig_file(p, caps)

    def page(canvas, doc):
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

    out_pdf = str(out_pdf or (L.report / "SpaceJEPA_v4_report.pdf"))
    tmp = out_pdf + ".tmp.pdf"
    L.guard(tmp), L.guard(out_pdf)
    L.report.mkdir(parents=True, exist_ok=True)
    doc = SimpleDocTemplate(tmp, pagesize=letter, leftMargin=0.75 * inch, rightMargin=0.75 * inch,
                            topMargin=0.75 * inch, bottomMargin=0.75 * inch, title=TITLE)
    doc.build(d.story, onFirstPage=page, onLaterPages=page)
    os.replace(tmp, out_pdf)
    return out_pdf


if __name__ == "__main__":
    from jepa import study_v4 as S4
    print(build_report(S4.LayoutV4(sys.argv[1])))
