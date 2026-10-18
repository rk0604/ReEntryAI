"""
jepa/make_results_report_v2.py

Five-page PDF of the first full SpaceJEPA v2 run (Colab GPU, one seed):
closed loop, open loop, dense air, load, uncertainty, verdict, caveats, next steps.
Inputs are the notebook's CSVs from MyDrive/ReEntryAI_jepa_v2/full/ (plus the
along-track bias figure the notebook draws, saved as a PNG).

    python jepa/make_results_report_v2.py --dir <folder with the CSVs> [--bias-png <png>] [--out <pdf>]

Expected files in --dir: closed_loop.csv, footprint_wide_test.csv,
footprint_v1_home.csv, decision_criteria.csv (names may carry a prefix).
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import sys
import tempfile
from pathlib import Path

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer

from jepa.make_report_v2 import MUT, RULE, TINT_YELLOW, Doc, _fonts
from reportlab.lib import colors

COL = {"classical": "#6250d6", "v1": "#2a78d6", "v1-wide": "#1baf7a", "v2": "#eb6834", "v2 (no unc.)": "#eda100",
       "physics-nominal": "#6250d6", "Physics prior only": "#a3a29c"}
CTRL = ["classical", "v1", "v1-wide", "v2", "v2 (no unc.)"]
TRUTHS = ["nominal", "mismatch", "mismatch_high", "out_of_range"]
TLAB = {"nominal": "nominal", "mismatch": "thin ×0.88", "mismatch_high": "dense ×1.12", "out_of_range": "×0.78/0.80"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": "#e4e3df", "legend.frameon": False,
                     "axes.titlesize": 8.2, "savefig.dpi": 220})


def _find(d, name):
    hits = sorted(glob.glob(os.path.join(d, f"*{name}")))
    if not hits:
        raise FileNotFoundError(f"{name} not in {d}")
    return hits[0]


def sign_p(wins, n):
    """Exact two-sided sign test."""
    k = min(wins, n - wins)
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n)


def analyse(d):
    cl = pd.read_csv(_find(d, "closed_loop.csv"))
    fw = pd.read_csv(_find(d, "footprint_wide_test.csv"))
    fh = pd.read_csv(_find(d, "footprint_v1_home.csv"))
    miss = cl.pivot_table(index=["tag", "case"], columns="controller", values="miss_km")
    g = cl.pivot_table(index=["tag", "case"], columns="controller", values="peak_g")
    summ = cl.groupby(["tag", "controller"]).agg(succ=("success", "mean"), med=("miss_km", "median"),
                                                 p90=("miss_km", lambda s: s.quantile(0.9)),
                                                 g=("peak_g", "mean"), fuel=("rcs_fuel_kg", "mean"))
    paired = {}
    for t in TRUTHS:
        for o in ("classical", "v1", "v1-wide", "v2 (no unc.)"):
            dlt = miss.loc[t, "v2"] - miss.loc[t, o]
            w, n = int((dlt < 0).sum()), int((dlt != 0).sum())
            paired[(t, o)] = (float(dlt.median()), w, n, sign_p(w, n))
    fails = {t: {c: sorted(miss.loc[t].index[miss.loc[t, c] >= 25].tolist()) for c in CTRL} for t in TRUTHS}
    return dict(cl=cl, fw=fw, fh=fh, miss=miss, g=g, summ=summ, paired=paired, fails=fails,
                crit=pd.read_csv(_find(d, "decision_criteria.csv")))


def fig_closed(A, out):
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.45), gridspec_kw={"width_ratios": [1.35, 1]})
    x = np.arange(4)
    w = 0.16
    for j, c in enumerate(CTRL):
        v = [A["summ"].loc[(t, c), "med"] for t in TRUTHS]
        axes[0].bar(x + (j - 2) * w, v, w * 0.92, color=COL[c], label=c,
                    hatch="///" if c == "v2 (no unc.)" else None, edgecolor="white", linewidth=0.3)
    axes[0].set_xticks(x, [TLAB[t] for t in TRUTHS])
    axes[0].set(yscale="log", ylabel="median miss [km]", title="Median miss per atmosphere (16 cases each)")
    axes[0].axhline(25, color="#52514e", lw=0.7, ls=":")
    axes[0].legend(ncol=5, fontsize=6.2, loc="upper center", bbox_to_anchor=(0.5, -0.13))
    for i, t in enumerate(TRUTHS):
        d = A["miss"].loc[t, "v2"] - A["miss"].loc[t, "classical"]
        jit = (np.arange(len(d)) % 5 - 2) * 0.05
        axes[1].scatter(np.full(len(d), i) + jit, d.clip(-30, 10), s=9, color=["#eb6834" if v < 0 else "#6250d6" for v in d])
    axes[1].axhline(0, color="#52514e", lw=0.8)
    axes[1].set_xticks(range(4), [TLAB[t].split(" ")[0] if t != "out_of_range" else "×0.78" for t in TRUTHS])
    axes[1].set(ylabel="v2 − classical miss [km] (clipped)", title="Per case: below 0 = v2 better")
    fig.tight_layout(w_pad=1.2)
    p = os.path.join(out, "r_closed.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def fig_open(A, out):
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.35), sharey=True)
    for ax, df, title in ((axes[0], A["fw"], "Wide test set (~520 flights)"),
                          (axes[1], A["fh"], "v1's own test flights (46)")):
        m = df.groupby(["model", "ttg_s"]).err_km.median().unstack(0)
        for name in ("physics-nominal", "Physics prior only", "v1", "v1-wide", "v2"):
            ax.plot(m.index, m[name], marker="o", ms=3.5, color=COL[name], label=name,
                    ls="--" if name == "physics-nominal" else ":" if name == "Physics prior only" else "-")
        ax.set(yscale="log", xlabel="time-to-go when predicting [s]", title=title)
        ax.invert_xaxis()
    axes[0].set_ylabel("median landing-point error [km]")
    axes[0].legend(fontsize=6.2, loc="lower left")
    fig.tight_layout(w_pad=1.0)
    p = os.path.join(out, "r_open.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def fig_load(A, out):
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.2), gridspec_kw={"width_ratios": [1.3, 1]})
    x = np.arange(4)
    w = 0.16
    for j, c in enumerate(CTRL):
        axes[0].bar(x + (j - 2) * w, [A["summ"].loc[(t, c), "g"] for t in TRUTHS], w * 0.92, color=COL[c], label=c,
                    hatch="///" if c == "v2 (no unc.)" else None, edgecolor="white", linewidth=0.3)
    axes[0].set_xticks(x, [TLAB[t] for t in TRUTHS])
    axes[0].set(ylabel="mean peak load [g]", title="Peak g-load (soft limit for learned planners: 8 g)")
    axes[0].legend(ncol=3, fontsize=6.0, loc="upper right")
    axes[0].set_ylim(0, 8.5)
    cl = A["cl"]
    v = cl[cl.controller.isin(["v2", "v2 (no unc.)"])]
    sp = v.pivot_table(index=["tag", "case"], columns="controller", values="spread_km_mean")
    ms = v.pivot_table(index=["tag", "case"], columns="controller", values="miss_km")
    axes[1].scatter(sp["v2"], ms["v2"], s=9, color=COL["v2"], label="v2")
    axes[1].set(xlabel="mean ensemble spread at decisions [km]", ylabel="miss [km]", yscale="log",
                title="Closed loop: spread vs final miss")
    fig.tight_layout(w_pad=1.2)
    p = os.path.join(out, "r_load.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def _page(canvas, doc):
    """Footer with the title and page number."""
    canvas.saveState()
    canvas.setFont("DV", 7.5)
    canvas.setFillColor(MUT)
    canvas.drawString(0.75 * inch, 0.45 * inch, "SpaceJEPA v2: first full results")
    canvas.drawRightString(letter[0] - 0.75 * inch, 0.45 * inch, str(doc.page))
    canvas.restoreState()


def build(A, figs, bias_png, out_pdf):
    _fonts()
    d = Doc()
    S, P = A["summ"], A["paired"]
    s = lambda t, c, k: S.loc[(t, c), k]
    fw = A["fw"].groupby(["model", "ttg_s"]).err_km.median().unstack(0)
    fh = A["fh"].groupby(["model", "ttg_s"]).err_km.median().unstack(0)

    d.add(Paragraph("SpaceJEPA v2: first full results", d.s["title"]),
          Paragraph("Colab GPU run, one seed · 16 frozen test cases × 4 atmospheres × 5 controllers · "
                    "open-loop tests on ~520 held-out flights", d.s["subtitle"]))
    d.callout("The short version",
              [f"<b>v2 is a large step up from v1 and beats the classical controller in 3 of 4 atmospheres.</b> In nominal "
               f"air its median miss is {s('nominal','v2','med'):.2f} km (classical {s('nominal','classical','med'):.2f}, "
               f"v1 {s('nominal','v1','med'):.2f}). In thin air it is {s('mismatch','v2','med'):.2f} km (classical "
               f"{s('mismatch','classical','med'):.2f}). In the extreme ×0.78 atmosphere it is "
               f"{s('out_of_range','v2','med'):.2f} km, while the classical controller misses by "
               f"{s('out_of_range','classical','med'):.1f} km and fails 10 of 16 cases.",
               f"<b>v1's dense-air problem is fixed:</b> {s('mismatch_high','v1','med'):.1f} km → "
               f"{s('mismatch_high','v2','med'):.2f} km. The classical controller is still slightly better there "
               f"({s('mismatch_high','classical','med'):.2f} km); the gap is small and not robust (below).",
               "<b>Open loop, v2 predicts the landing point 10× better than the classical controller's own physics "
               f"model</b> 5 minutes out ({fw.loc[300,'v2']:.1f} vs {fw.loc[300,'physics-nominal']:.1f} km).",
               "<b>Two honest catches:</b> v2 pulls more g than the classical controller (about 5–6 g against 3.5–4), "
               "so the comparison is not load-matched. And the ensemble's uncertainty penalty made no measurable difference to "
               "miss distance. 5 of the 6 pre-registered criteria passed."], TINT_YELLOW, colors.HexColor("#eda100"))
    d.h2("Pre-registered criteria (fixed before the run)")
    rows = [["#", "Criterion", "Measured", "Met"]]
    names = ["No collapse of the latent codes", "Gain comes from the model, not just the data",
             "Beats the classical controller's physics model", "Dense-air bias gone", "Robust in closed loop",
             "Uncertainty is informative"]
    for i, r in A["crit"].iterrows():
        rows.append([str(i + 1), names[i], r["measured"], "yes" if str(r["met"]) == "True" else "<b>no</b>"])
    d.table(rows, [0.25, 2.25, 3.95, 0.55])
    d.p("Criterion 5 failed on two counts. The first is v2's 81% success in the ×0.78 atmosphere, but its 3 failures "
        "(cases 4, 5, 14) are cases every controller fails, and v2 misses each of them by the least. The second is "
        "a 0.82 km gap to the classical controller in dense air, below the plan's own ~1 km noise level. Criterion "
        "6 passed only narrowly, and in closed loop the uncertainty term did nothing measurable (Section 4).")

    d.h1("1. Closed loop: flying the frozen test cases")
    d.fig(figs["closed"], "Left: median miss per atmosphere (log scale; dotted line = 25 km success radius). Right: "
                          "v2 minus classical, per test case. Orange = v2 better, violet = classical better.", "real")
    rows = [["Atmosphere", "classical", "v1", "v1-wide", "v2", "v2 (no unc.)"]]
    for t in TRUTHS:
        rows.append([TLAB[t]] + [f"{s(t,c,'med'):.2f} km · {100*s(t,c,'succ'):.0f}%" for c in CTRL])
    d.table(rows, [1.3, 1.14, 1.14, 1.14, 1.14, 1.14])
    pr = lambda t, o: P[(t, o)]
    d.p("<b>Paired, case by case</b> (the same 16 starting conditions for every controller; exact sign test). Each "
        "entry is the median of v2 − other, how many cases v2 wins, and the p-value:")
    rows = [["Atmosphere", "v2 vs classical", "v2 vs v1", "v2 vs v1-wide"]]
    for t in TRUTHS:
        rows.append([TLAB[t]] + [f"{pr(t,o)[0]:+.2f} km, wins {pr(t,o)[1]}/{pr(t,o)[2]}, p={pr(t,o)[3]:.3f}"
                                 for o in ("classical", "v1", "v1-wide")])
    d.table(rows, [1.3, 1.9, 1.9, 1.9])
    d.bullets(["<b>Clear wins:</b> v2 beats the classical controller in nominal air (14/16 cases), thin air (15/16) and "
               "the ×0.78 atmosphere (16/16). It beats v1 in thin, dense and ×0.78 air (15–16/16).",
               "<b>Not established:</b> v2 against v1 in nominal air (11/16, p = 0.21). Both are already around 1–2 km there.",
               f"<b>Dense air, the one loss:</b> the classical controller is better in 13 of 16 cases, but by a median of only "
               f"{abs(pr('mismatch_high','classical')[0]):.2f} km. The result is significant before correcting for the 16 comparisons "
               "(p = 0.02) and not after (Holm). It also does not repeat with the near-identical no-uncertainty run. "
               "Treat it as \"about equal, classical slightly ahead\".",
               "<b>Failures overlap:</b> every v2 failure is also a classical and v1 failure. Under ×0.78 air the "
               "classical controller misses 6 cases by between 23 and 27 km, so its 37.5% success depends on where "
               "the 25 km line is drawn.",
               "<b>The cleanest comparison is the ×0.78 atmosphere,</b> because there the load is matched (v2 4.04 g, "
               "classical 4.06 g): v2 wins all 16 cases."])

    d.h1("2. Open loop: predicting where the capsule lands")
    d.fig(figs["open"], "Median error of the predicted parachute point, given the state at a time-to-go and the "
                        "commands actually flown. Left: held-out wide-envelope flights. Right: v1's own narrow-corridor "
                        "test flights, which v2 never saw.", "real")
    rows = [["Time-to-go", "physics-nominal", "v1", "v1-wide", "v2", "v2 gain vs physics-nominal"]]
    for t in (400.0, 300.0, 200.0, 100.0, 50.0):
        rows.append([f"{t:.0f} s", f"{fw.loc[t,'physics-nominal']:.2f}", f"{fw.loc[t,'v1']:.2f}",
                     f"{fw.loc[t,'v1-wide']:.2f}", f"<b>{fw.loc[t,'v2']:.2f}</b>",
                     f"{fw.loc[t,'physics-nominal']/fw.loc[t,'v2']:.0f}×"])
    d.table(rows, [1.0, 1.25, 1.0, 1.0, 1.0, 1.75])
    bins = A["fw"].assign(b=pd.cut(A["fw"].dom_density_scale / A["fw"].dom_mass_scale, [0, .85, 1.15, 9],
                                   labels=["thin", "near", "dense"]))
    bm = bins[bins.ttg_s == 200].pivot_table(index="model", columns="b", values="err_km", aggfunc="median",
                                             observed=True)
    d.bullets(["<b>v2 wins at every time-to-go</b>, and the margin grows with the horizon: 2× better than the classical "
               "controller's model at 50 s, 13× at 400 s. The single members are almost as good as the ensemble mean, so "
               "the gain comes from the architecture and data, not from averaging.",
               f"<b>Its accuracy barely depends on the atmosphere</b>, which is the sign of real system identification. At "
               f"200 s: thin {bm.loc['v2','thin']:.1f}, near-nominal {bm.loc['v2','near']:.1f}, dense "
               f"{bm.loc['v2','dense']:.1f} km. The classical model gets {bm.loc['physics-nominal','thin']:.1f} / "
               f"{bm.loc['physics-nominal','near']:.1f} / {bm.loc['physics-nominal','dense']:.1f} km.",
               f"<b>It does not lose on v1's home turf:</b> at 300 s, v2 {fh.loc[300,'v2']:.1f} km against v1 "
               f"{fh.loc[300,'v1']:.1f} km and physics-nominal {fh.loc[300,'physics-nominal']:.1f} km.",
               "<b>The data alone was not enough:</b> v1's architecture retrained on the same 6,000 flights (v1-wide) is no "
               "better than the physics model in open loop. In closed loop it is worse than v1 in nominal and dense air. "
               "It did get a smaller training budget (6,000 steps on CPU)."])

    d.h1("3. The dense-air question")
    if bias_png and os.path.exists(bias_png):
        d.fig(bias_png, "Mean signed along-track error of the predicted landing point against the true density/mass "
                        "factor (+ = predicted beyond the truth). v1 and v1-wide swing by tens of km at the edges; v2 "
                        "stays near zero at every density and time-to-go.", "real", width=6.6 * inch)
    d.p(f"The bias that drove v1's 12–14 km dense-air misses is gone in prediction (criterion 4: mean along-track "
        f"error +0.2 km for v2 against +1.1 km for v1 at 200 s). In closed loop, v2's dense-air median dropped from "
        f"{s('mismatch_high','v1','med'):.1f} to {s('mismatch_high','v2','med'):.2f} km, a win in 15 of 16 cases. The "
        f"remaining gap to the classical controller ({s('mismatch_high','classical','med'):.2f} km) is small. In "
        f"dense air the nominal model the classical controller relies on is only 12% off, so it is the classical "
        f"controller's best mismatched case.")

    d.h1("4. Load, uncertainty, and cost")
    d.fig(figs["load"], "Left: mean peak g per atmosphere. Right: v2's mean ensemble spread at decision time against "
                        "its final miss, one dot per flight.", "real")
    d.bullets([f"<b>v2 is not load-matched.</b> Its mean peak is {s('nominal','v2','g'):.1f} g against the classical controller's "
               f"{s('nominal','classical','g'):.1f} g in nominal air, and {s('mismatch_high','v2','g'):.1f} vs "
               f"{s('mismatch_high','classical','g'):.1f} g in dense air. That is v1's pattern too: the 8 g soft limit lets "
               "the planner trade load for accuracy. Within v2, higher g does not go with lower miss, but a "
               "g-matched rerun is needed before claiming a win over the classical controller.",
               "<b>The uncertainty penalty does what it optimizes but not what we hoped.</b> It lowers the planner's "
               "spread in 57 of 64 flights (by ~5%). It changes the miss in no atmosphere (sign tests p = 0.8–1.0), "
               "and the spread did not flag the large-miss cases. In open loop the spread does track error (Spearman "
               "0.44 at 200 s).",
               f"<b>Fuel:</b> the learned planners use about half the classical controller's RCS fuel "
               f"({s('nominal','v2','fuel'):.0f} vs {s('nominal','classical','fuel'):.0f} kg in nominal air). The "
               f"saving comes from MPPI's smoother commands and is shared by v1."])

    d.h1("5. Caveats and next steps")
    d.p("<b>What limits these conclusions:</b> one training seed; 16 cases per atmosphere, the same 16 starting "
        "conditions reused across atmospheres; no load matching; a smaller training budget for v1-wide; v1 evaluated "
        "mostly outside its training range in the wide open-loop test. v2 also flew through the new GPU runner and "
        "v1 through the original CPU runner. They were checked to make the same decisions for a single model, not "
        "for the 5-model ensemble. The code audit found no data leakage: splits are by flight, statistics and the λ "
        "choice use training and validation data only, and the encoder cannot see the future.")
    d.numbered(["<b>Load-matched rerun:</b> set the MPPI g-limit to ~4.5 g (or add a g-preference term) and fly the same "
                "cases. This is the one result needed before claiming v2 beats the classical controller. Cost: ~15 min of closed loop.",
                "<b>More cases:</b> raise n_cases to 48 (classical runs are cached, so only new cases cost time) to settle "
                "the dense-air comparison and the success rates.",
                "<b>Train longer:</b> the stage-1 loss was still falling at 24,000 steps with train ≈ val, so there is "
                "no overfitting yet. Try 48,000 steps, and a second seed to measure run-to-run spread.",
                "<b>Rethink the uncertainty term:</b> it has no closed-loop effect at weight 1. Try a larger weight "
                "only late in flight, or drop it and keep plain ensemble averaging.",
                "<b>Infeasible cases:</b> a full lift-up flight in the true ×0.78 atmosphere lands 74 km (case 5) and "
                "82 km (case 14) short of the target, so no controller can reach it. Score success only on reachable "
                "cases, and check case 4 the same way."])

    doc = SimpleDocTemplate(out_pdf, pagesize=letter, leftMargin=0.75 * inch, rightMargin=0.75 * inch,
                            topMargin=0.6 * inch, bottomMargin=0.7 * inch, title="SpaceJEPA v2: first full results")
    doc.build(d.story, onFirstPage=_page, onLaterPages=_page)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--bias-png", default="")
    ap.add_argument("--out", default=str(SIM_DIR / "jepa" / "SpaceJEPA_v2_results.pdf"))
    a = ap.parse_args()
    A = analyse(a.dir)
    tmp = tempfile.mkdtemp(prefix="jepa_results_")
    figs = {"closed": fig_closed(A, tmp), "open": fig_open(A, tmp), "load": fig_load(A, tmp)}
    build(A, figs, a.bias_png, a.out)
    from pypdf import PdfReader
    print(f"wrote {a.out} ({len(PdfReader(a.out).pages)} pages)")
