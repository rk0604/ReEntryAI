"""
jepa/report_v5.py

06_report/SpaceJEPA_v5_report.pdf (jepa/v5_plan.md section 10). It reads only the files analyze_v5 wrote
(05_analysis/tables, figures, captions.json), the frozen settings, the tuned configurations (03_tuning/*/chosen.json)
and the plan's section 15 (deviations); no number is typed in. Building blocks as v4 (results_pdf_v3.Doc,
make_report_v2._fonts).

Sections: summary callout, criteria, RQ1-RQ5, benchmark datasheet, methods (with the exact FNPEG and PredGuid
settings), deviations, appendix with every table and any figure not shown above.

    python jepa/report_v5.py <v5 study folder>
"""

from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path
from typing import Optional
from xml.sax.saxutils import escape

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

import numpy as np
import pandas as pd

from jepa import report_v4 as R4

TITLE = "SpaceJEPA v5: strong baselines and a standard benchmark"
PLAN = SIM_DIR / "jepa" / "v5_plan.md"
ATM_LAB = R4.ATM_LAB
_csv, _f = R4._csv, R4._f


def deviations_text() -> str:
    try:
        txt = PLAN.read_text(encoding="utf-8")
    except OSError:
        return "(v5_plan.md not found)"
    m = re.search(r"^## 15\. Deviations\s*$(.*)", txt, flags=re.M | re.S)
    return m.group(1).strip() if m else "(no section 15 in v5_plan.md)"


def compact(key: str, df: Optional[pd.DataFrame]):
    if df is None or len(df) <= R4.MAX_ROWS:
        return df, ""
    if key == "T12":
        d = df[df.threshold_km.isin([1.0, 2.0, 5.0, 10.0, 25.0])]
        d = d.pivot_table(index=["suite", "arm"], columns="threshold_km", values="success_%").reset_index()
        d.columns = [c if isinstance(c, str) else f"success_%_{c:g}km" for c in d.columns]
        return d, f"Shown at 1, 2, 5, 10, 25 km; the CSV has all {len(df)} rows."
    return df.head(R4.MAX_ROWS), f"First {R4.MAX_ROWS} of {len(df)} rows; the CSV has all of them."


def _inputs(L) -> list:
    from jepa import benchmark_v5 as B
    files = [L.study / "study_config.json", PLAN, Path(__file__), B.BENCH_DIR / "ORB1_spec.json"]
    files += [p for p in L.analysis.rglob("*") if p.is_file() and "cache" not in p.parts]
    files += list(L.tuning.glob("*/chosen.json"))
    return files


def is_current(L, pdf: Path) -> bool:
    """The report is newer than everything it reads (analysis outputs, config, plan, ORB-1 spec, tuning choices)."""
    try:
        t = pdf.stat().st_mtime
    except OSError:
        return False
    return all(p.stat().st_mtime <= t for p in _inputs(L) if p.exists())


def build_report(L, out_pdf: Optional[str] = None, force: bool = False) -> str:
    """Write the report from the files in L (a study_v5.LayoutV5); returns its path. The study's report is not
    rebuilt when it is newer than every file it reads (unless force)."""
    if out_pdf is None and not force and is_current(L, L.report / "SpaceJEPA_v5_report.pdf"):
        return str(L.report / "SpaceJEPA_v5_report.pdf")
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
    t01, t05, t08, t09 = (t.get(k) for k in ("T01", "T05", "T08", "T09"))

    d.add(Paragraph(TITLE, d.s["title"]),
          Paragraph(f"Published guidance laws, explicit estimation, non-JEPA models, RL, and ORB-1 · tier "
                    f"{cfg.get('tier', '?')} · mode {cfg.get('mode', '?')} · git {str(cfgdoc.get('git_commit', '?'))[:10]}",
                    d.s["subtitle"]))

    # ---- summary callout ----
    lines = []
    if t01 is not None and len(t01):
        a = t01[(t01.suite == "A256") & (t01.group == "all")]
        for k in ("classical", "fnpeg", "predguid", "mppi_filtered", "v2single", "v2single@fnpeg_matched"):
            r = a[a.arm == k]
            if len(r):
                r = r.iloc[0]
                lines.append(f"<b>{escape(k)}</b> (A256, pooled): success@5km {100 * r.success_5km:.1f}%, mean peak g "
                             f"{r.peak_g_mean:.2f}, median miss {r.median_miss_km:.2f} km")
    if t08 is not None and len(t08):
        for r in t08[t08.group == "all"].itertuples():
            lines.append(f"<b>ORB-1 {escape(str(r.arm))}</b>: {int(r.n)} flights (blocks {r.blocks_flown}), success@10km "
                         f"{100 * r.success_10km:.2f}%, requirement LCB {r.requirement_lcb:.4f} "
                         f"({'met' if str(r.requirement_met).lower() == 'true' else 'not met'})")
    if t05 is not None and len(t05):
        test = t05[t05.met.notna() & (t05.met.astype(str) != "")]
        met = test[test.met.astype(str).str.lower() == "true"]
        lines.append(f"<b>Criteria met: {len(met)} of {len(test)}</b> pass/fail criteria: "
                     f"{', '.join(met.id.astype(str)) or 'none'}.")
    d.callout("The numbers (read from 05_analysis/tables)", lines or ["No analysis tables yet."],
              R2.TINT_YELLOW, colors.HexColor("#eda100"))
    if t05 is not None and len(t05):
        d.h2("Pre-registered criteria (v5_plan.md section 8, fixed before any v5 flight)")
        rows = [["#", "Question", "Measured", "p (Holm)", "Met"]]
        for r in t05.itertuples():
            met = "estimate" if str(r.met) in ("None", "nan", "") else ("yes" if str(r.met).lower() == "true" else "<b>no</b>")
            if str(getattr(r, "valid", "")).lower() == "false":
                met += " (invalid: G3)"
            rows.append([escape(str(r.id)), escape(str(r.question)), escape(str(r.measured)), _f(r.p_holm, "{:.3g}"), met])
        d.table(rows, [0.8, 1.8, 3.1, 0.6, 0.6])

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
    d.h1("RQ1. Published guidance laws and explicit estimation")
    d.p("FNPEG (Lu 2014, with predictive load relief and fading-memory filters) and PredGuid (Bairstow 2006 App. B) are "
        "the published, load-limited, model-adapting classical side; MPPI + filtered physics is v2's planner with the "
        "physics model's density and lift scales set every decision from FNPEG's filters (explicit estimation). v2 "
        "single member flies v3's checkpoints; at FNPEG's own realized peak g it is the load-matched comparison.")
    fig("F01"), fig("F02")
    t04 = t.get("T04")
    if t04 is not None and len(t04):
        tab("T04", t04[t04.atmosphere == "all"].drop(columns=["per_seed"], errors="ignore"))
    tab("T03")
    d.h1("RQ2. Is the JEPA objective needed?")
    d.p("Four non-JEPA dynamics models trained on the same 6,000 flights with free-running multi-step rollouts, sized "
        "to R0's parameter count, flown under the same MPPI at 4.5 g. The primary metric is the paired log-miss "
        "difference against R0 (v4's one-member v2 recipe).")
    fig("F03"), tab("T06"), fig("F04")
    d.h1("RQ3. Reinforcement-learning baselines")
    d.p("Policies trained in a batched fast environment (the physics port, 5 s decisions, the history observation) and "
        "evaluated in the full simulator like every other arm; the fast-vs-full gap is reported per arm.")
    fig("F05"), tab("T07"), tab("T07b"), tab("T11")
    d.h1("RQ4. ORB-1: an Orion-class LEO return")
    d.p("Published dispersions (Bairstow Table C.1, entry speed per D4), three target ranges inside the nominal "
        "footprint, a sequential stopping rule over three blocks of 1,000 cases. Q1 is the accuracy requirement "
        "(one-sided 90% lower bound on P(miss ≤ 10 km) ≥ 0.99865), Q2 crew g-duration limits, Q3 v2 against the "
        "published laws on block 1.")
    fig("F06"), fig("F07"), fig("F08"), tab("T08")
    d.h1("RQ5. Compute per decision")
    tab("T09")
    d.h1("Success against the threshold")
    fig("F09"), fig("F10")
    d.h1("Benchmark datasheet")
    try:
        from jepa import benchmark_v5 as B
        spec = json.loads((B.BENCH_DIR / "ORB1_spec.json").read_text())
        fp = spec.get("footprint", {})
        d.bullets([
            "<b>A256</b>: v3's 64 cases × 4 atmospheres (nominal, thin ×0.88, dense ×1.12, out of range ×0.78 / L/D "
            "×0.80), unchanged (configs/v3_testset.csv, checksum pointer in A256.json).",
            "<b>ORB-1</b>: " + escape(", ".join(f"{k} {v}" for k, v in (spec.get("spec") or {}).items()
                                               if not isinstance(v, (list, dict)))[:900]),
            "<b>Footprint</b> (constant bank, nominal entry): " + escape(json.dumps(fp)[:600]),
            "<b>Targets</b>: " + escape("; ".join(f"{x['target_index']}: {x['downrange_km']} km" for x in spec.get("targets", []))),
            "<b>DEV48</b>: tuning only, never a result (configs/benchmark_v5/DEV48.csv).",
        ])
    except Exception as e:                                   # the datasheet must not cost the report
        d.p(escape(f"(benchmark datasheet unavailable: {e})"))
    d.h1("Methods and definitions")
    tuned = {}
    for m in ("fnpeg", "predguid"):
        p = L.tuning / m / "chosen.json"
        if p.exists():
            tuned[m] = json.loads(p.read_text())
    d.bullets([
        "<b>success@5km / @2km</b>: outcome drogue (7.6 km) and great-circle miss below 5 / 2 km; <b>success@10km</b>: "
        "miss ≤ 10 km. <b>log miss</b> = ln(max(miss, 0.01 km)); failures keep their miss at termination.",
        "<b>Paired comparisons</b>: the unit is one condition, paired across arms by entry state and MPPI noise stream "
        "(common random numbers). Deterministic arms are broadcast over seeds; RL comparisons pair seeds 0–2.",
        f"<b>Uncertainty</b>: hierarchical bootstrap, {cfg.get('n_boot')} resamples (seeds, then cases within each "
        "atmosphere), 95% percentile intervals (90% for TOST); ORB-1: cases resampled within training-seed strata; "
        "per-seed exact McNemar and sign tests; Holm within each family.",
        "<b>FNPEG settings (chosen on DEV48)</b>: " + escape(json.dumps(tuned.get("fnpeg", {}).get("config", "defaults"),
                                                                         default=str)),
        "<b>PredGuid settings (chosen on DEV48)</b>: " + escape(json.dumps(tuned.get("predguid", {}).get("config", "defaults"),
                                                                            default=str)),
        "<b>Planner</b>: MPPI as in v3 (256 samples, 6 knots, 30° noise, 4 iterations at entry then 1, horizon 188 "
        "steps of 5 s, soft g limit 50 km per g² above the limit).",
    ])
    d.h1("Deviations (v5_plan.md section 15)")
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

    out_pdf = Path(out_pdf or (L.report / "SpaceJEPA_v5_report.pdf"))
    L.report.mkdir(parents=True, exist_ok=True)
    import tempfile
    tmp = Path(tempfile.gettempdir()) / f"jepa_v5_report_{os.getpid()}.pdf"   # built locally, then one guarded write
    doc = SimpleDocTemplate(str(tmp), pagesize=letter, leftMargin=0.75 * inch, rightMargin=0.75 * inch,
                            topMargin=0.75 * inch, bottomMargin=0.75 * inch, title=TITLE)
    doc.build(d.story, onFirstPage=page, onLaterPages=page)
    L.save_large(tmp.read_bytes(), out_pdf)
    tmp.unlink(missing_ok=True)
    return str(out_pdf)


if __name__ == "__main__":
    from jepa import study_v5 as S5
    print(build_report(S5.LayoutV5(sys.argv[1])))
