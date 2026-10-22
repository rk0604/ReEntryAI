"""
jepa/make_report_v2.py

Build "SpaceJEPA v2: a beginner's guide" (PDF): part 1 explains the pipeline
and the experiments, part 2 the models and the architecture. Figures come from
jepa/report_v2_figures.py; every caption says whether it shows real data,
numbers from the code, an illustration, or an estimate.

    python jepa/make_report_v2.py --wide <episodes pickle or shard dir> \
        [--dense <dense_runs.pkl>] [--out jepa/SpaceJEPA_v2_report.pdf]

--dense is a pickle of jepa.mppi.run_episode results (classical and v1 under
density/L-D x1.12 on the first 8 frozen test cases). With --run-dense and a
--dense path that does not exist yet, those 16 flights are flown first (~2 min
on 4 cores) and saved there. Without either, that figure and its paragraph are
left out.
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (Image, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table,
                                TableStyle)

FONT_DIR = "/usr/share/fonts/truetype/dejavu"
INK, SEC, MUT = colors.HexColor("#1b1b1a"), colors.HexColor("#52514e"), colors.HexColor("#8a8984")
BLUE, ORANGE, VIOLET = colors.HexColor("#2a78d6"), colors.HexColor("#eb6834"), colors.HexColor("#4a3aa7")
TINT_BLUE, TINT_ORANGE, TINT_YELLOW = colors.HexColor("#eef4fc"), colors.HexColor("#fdf0ea"), colors.HexColor("#fdf6e3")
RULE = colors.HexColor("#d9d8d3")
TEXT_W = 7.0 * inch


def _font_dir() -> str:
    """The system DejaVu fonts, else the copies matplotlib ships (Colab and Windows have no system DejaVu)."""
    if os.path.exists(os.path.join(FONT_DIR, "DejaVuSans.ttf")):
        return FONT_DIR
    import matplotlib
    return os.path.join(matplotlib.get_data_path(), "fonts", "ttf")


def _fonts():
    """DejaVu Sans covers Greek and math symbols (the built-in PDF fonts do not)."""
    fd = _font_dir()
    pdfmetrics.registerFont(TTFont("DV", os.path.join(fd, "DejaVuSans.ttf")))
    pdfmetrics.registerFont(TTFont("DVB", os.path.join(fd, "DejaVuSans-Bold.ttf")))
    pdfmetrics.registerFont(TTFont("DVM", os.path.join(fd, "DejaVuSansMono.ttf")))
    from reportlab.pdfbase.pdfmetrics import registerFontFamily
    registerFontFamily("DV", normal="DV", bold="DVB", italic="DV", boldItalic="DVB")


def _styles():
    s = {}
    s["title"] = ParagraphStyle("title", fontName="DVB", fontSize=26, leading=31, textColor=INK, spaceAfter=10)
    s["subtitle"] = ParagraphStyle("subtitle", fontName="DV", fontSize=13, leading=18, textColor=SEC, spaceAfter=14)
    s["part"] = ParagraphStyle("part", fontName="DVB", fontSize=20, leading=25, textColor=BLUE, spaceAfter=6)
    s["h1"] = ParagraphStyle("h1", fontName="DVB", fontSize=14, leading=18, textColor=INK, spaceBefore=12, spaceAfter=6,
                             keepWithNext=1)
    s["h2"] = ParagraphStyle("h2", fontName="DVB", fontSize=11, leading=14, textColor=INK, spaceBefore=8, spaceAfter=4,
                             keepWithNext=1)
    s["body"] = ParagraphStyle("body", fontName="DV", fontSize=9.6, leading=14, textColor=INK, spaceAfter=6)
    s["bullet"] = ParagraphStyle("bullet", parent=s["body"], leftIndent=14, bulletIndent=3, spaceAfter=3)
    s["caption"] = ParagraphStyle("caption", fontName="DV", fontSize=8, leading=11, textColor=SEC, spaceBefore=3,
                                  spaceAfter=10)
    s["callout_t"] = ParagraphStyle("callout_t", fontName="DVB", fontSize=9.4, leading=13, textColor=INK, spaceAfter=3)
    s["callout"] = ParagraphStyle("callout", fontName="DV", fontSize=9.2, leading=13, textColor=INK)
    s["cell"] = ParagraphStyle("cell", fontName="DV", fontSize=8.2, leading=10.8, textColor=INK)
    s["cellb"] = ParagraphStyle("cellb", fontName="DVB", fontSize=8.2, leading=10.8, textColor=INK)
    s["toc"] = ParagraphStyle("toc", fontName="DV", fontSize=10, leading=16, textColor=INK)
    s["small"] = ParagraphStyle("small", fontName="DV", fontSize=8.4, leading=12, textColor=SEC)
    return s


KIND = {"real": ("REAL DATA", colors.HexColor("#1baf7a")), "code": ("FROM THE CODE", BLUE),
        "illus": ("ILLUSTRATION", MUT), "v1": ("v1 RESULTS (PLAN.md §9)", VIOLET), "est": ("ESTIMATE", ORANGE)}


class Doc:
    """Small helper that turns calls into a reportlab story."""

    def __init__(self):
        self.s = _styles()
        self.story = []
        self.fig_n = 0

    def add(self, *flow):
        self.story.extend(flow)

    def part(self, text, sub):
        self.add(PageBreak(), Spacer(1, 1.2 * inch), Paragraph(text, self.s["part"]),
                 Paragraph(sub, self.s["subtitle"]), Spacer(1, 6))

    def h1(self, text):
        self.add(Paragraph(text, self.s["h1"]))

    def h2(self, text):
        self.add(Paragraph(text, self.s["h2"]))

    def p(self, text):
        self.add(Paragraph(text, self.s["body"]))

    def bullets(self, items):
        for it in items:
            self.add(Paragraph(it, self.s["bullet"], bulletText="•"))
        self.add(Spacer(1, 4))

    def numbered(self, items):
        for i, it in enumerate(items, 1):
            self.add(Paragraph(it, self.s["bullet"], bulletText=f"{i}."))
        self.add(Spacer(1, 4))

    def fig(self, path, caption, kind, width=TEXT_W):
        from PIL import Image as PILImage
        self.fig_n += 1
        w, h = PILImage.open(path).size
        img = Image(path, width=width, height=width * h / w)
        tag, col = KIND[kind]
        cap = Paragraph(f'<font name="DVB" color="{col.hexval()}">{tag}</font>  '
                        f'<font name="DVB">Figure {self.fig_n}.</font> {caption}', self.s["caption"])
        head = []
        # A heading directly before a figure moves with it (keepWithNext does not reach into KeepTogether).
        if self.story and isinstance(self.story[-1], Paragraph) and self.story[-1].style.name in ("h1", "h2"):
            head = [self.story.pop()]
        self.add(KeepTogether(head + [img, cap]))
        return self.fig_n

    def callout(self, title, text, tint=TINT_BLUE, edge=BLUE):
        body = [Paragraph(title, self.s["callout_t"])]
        for para in (text if isinstance(text, list) else [text]):
            body.append(Paragraph(para, self.s["callout"]))
        t = Table([[body]], colWidths=[TEXT_W])
        t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), tint), ("LINEBEFORE", (0, 0), (0, -1), 3, edge),
                               ("LEFTPADDING", (0, 0), (-1, -1), 10), ("RIGHTPADDING", (0, 0), (-1, -1), 10),
                               ("TOPPADDING", (0, 0), (-1, -1), 7), ("BOTTOMPADDING", (0, 0), (-1, -1), 8)]))
        self.add(Spacer(1, 2), t, Spacer(1, 8))

    def table(self, rows, widths, header=True, zebra=True):
        data = [[Paragraph(str(c), self.s["cellb"] if (header and i == 0) else self.s["cell"]) for c in r]
                for i, r in enumerate(rows)]
        t = Table(data, colWidths=[w * inch for w in widths], repeatRows=1 if header else 0)
        st = [("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 0), (-1, 0), 0.8, INK),
              ("LINEBELOW", (0, -1), (-1, -1), 0.5, RULE), ("TOPPADDING", (0, 0), (-1, -1), 3.5),
              ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5), ("LEFTPADDING", (0, 0), (-1, -1), 4),
              ("RIGHTPADDING", (0, 0), (-1, -1), 4)]
        if zebra:
            for i in range(1, len(rows)):
                if i % 2 == 0:
                    st.append(("BACKGROUND", (0, i), (-1, i), colors.HexColor("#f6f5f2")))
        t.setStyle(TableStyle(st))
        self.add(t, Spacer(1, 10))


def _page(canvas, doc):
    """Header rule + footer with the title and page number (not on the cover)."""
    if doc.page == 1:
        return
    canvas.saveState()
    canvas.setFont("DV", 7.5)
    canvas.setFillColor(MUT)
    canvas.drawString(0.75 * inch, 0.5 * inch, "SpaceJEPA v2: a beginner's guide")
    canvas.drawRightString(letter[0] - 0.75 * inch, 0.5 * inch, str(doc.page))
    canvas.setStrokeColor(RULE)
    canvas.setLineWidth(0.5)
    canvas.line(0.75 * inch, letter[1] - 0.55 * inch, letter[0] - 0.75 * inch, letter[1] - 0.55 * inch)
    canvas.restoreState()


# ============================================================================
def build(figs: dict, info: dict, out_pdf: str) -> None:
    _fonts()
    d = Doc()
    F, I = figs, info
    dense = I.get("dense")
    mp = I["mppi"]
    feat = I["features"]
    sg = I["sigreg"]
    rep = I["params"]["rep"]
    order = ["basics", "jepa_idea", "pipeline", "data_examples", "coverage", "mppi_snapshot", "mppi_loop",
             "v1_results", "dense", "compute", "architecture", "features", "masks", "action_tokens", "parallel",
             "sigreg", "schedule", "masked_windows", "prober", "ensemble", "params"]
    N = {k: i + 1 for i, k in enumerate(k for k in order if k in F)}

    # ---------------------------------------------------------------- cover
    d.add(Spacer(1, 1.3 * inch),
          Paragraph("SpaceJEPA v2", d.s["title"]),
          Paragraph("A beginner's guide to the learned world model that steers the Orion capsule through reentry, "
                    "and to the bigger, stronger model built next", d.s["subtitle"]),
          Spacer(1, 4),
          Paragraph("Part 1: the pipeline and the experiments<br/>Part 2: the models and the architecture", d.s["toc"]),
          Spacer(1, 18))
    d.callout("Read this first: what exists today",
              ["<b>v1</b>, the first JEPA world model (26 thousand parameters), has been trained and tested. Its results "
               "are in Part 1. <b>v2</b>, the much larger model this report is about (12.3 million parameters), is fully "
               "built and has been checked for correctness on a CPU. It has <b>not</b> been trained at full size yet. That "
               "happens when you run the notebook <font name='DVM'>ReEntryAI_JEPA_v2.ipynb</font> on a Colab GPU.",
               "So this report explains how v2 works and how it will be tested. It contains no v2 accuracy numbers. "
               "Every number in it comes from v1, from real simulator runs, or from the v2 code itself."],
              TINT_YELLOW, colors.HexColor("#eda100"))
    d.add(Paragraph("<b>How figures are labeled.</b> Each caption starts with a tag:", d.s["body"]))
    d.table([["Tag", "Meaning"],
             ["REAL DATA", "Flights from the simulator: the training data, or controllers flying the frozen test cases."],
             ["v1 RESULTS", "Numbers measured in v1's full run (jepa/PLAN.md §9)."],
             ["FROM THE CODE", "Computed by the v2 code itself: attention masks, schedules, parameter counts."],
             ["ILLUSTRATION", "A drawing, or a made-up example that explains an idea."],
             ["ESTIMATE", "A planning number, not a measurement."]], [1.4, 5.6])
    d.add(Spacer(1, 10), Paragraph("ReEntryAI project · branch <font name='DVM'>jepa-v2</font> · "
                                   "pull request rk0604/ReEntryAI#1 · October 2026", d.s["small"]))

    # ================================================================ PART 1
    d.part("Part 1: the pipeline and the experiments",
           "What problem we solve, how a world model is trained and used to fly, what v1 showed, and how v2 will be tested")

    d.h1("1. The problem in one page")
    d.p("When the Orion capsule comes home, it hits the top of the atmosphere at about 120 km altitude and almost "
        "8 km per second, around 25 times the speed of sound. Over the next 6 to 8 minutes, air drag slows it down. "
        "At 7.6 km altitude the small <b>drogue parachute</b> opens. The goal is to make that happen as close as "
        "possible to a chosen target point. Landing within 25 km of it counts as a success.")
    d.p("The capsule has no wings and no steering engine. It flies slightly tilted, so its blunt shape makes a "
        "little <b>lift</b>, about 0.3 times its drag. The only thing guidance controls is the <b>bank angle</b> σ: "
        "how far the capsule rolls around its direction of flight. Rolling points the lift in a different direction "
        "(Figure 1, right). Lift pointed up makes the flight flatter and longer. Lift pointed sideways makes the "
        "capsule turn and drop faster. Real guidance flies at some bank angle and now and then flips its sign "
        "(a <b>bank reversal</b>) to keep the capsule headed for the target.")
    d.fig(F["basics"], "Left: altitude over time for 14 flights from the v2 training data, with one highlighted. "
                       "Each flight starts at the entry interface and ends when the parachute opens. Right: the bank "
                       "angle σ tilts the lift. (Left: real data; right: illustration.)", "real")
    d.p("<b>Why it is hard.</b> Nobody knows the exact air density or the capsule's exact aerodynamics on the day. "
        "If the air is 12% denser than expected, the capsule slows down sooner and falls short, unless guidance "
        "notices and corrects.")
    d.p("<b>The bar to beat</b> is the <b>classical predictor-corrector</b>. Every second it simulates the rest of the "
        "flight with a fixed, nominal physics model for several bank angles and picks the one that lands on target. "
        "On 1,000 frozen test cases it reaches a 2.05 km median miss and 99.9% success. Its weak spot is that its "
        "model never changes: it does not learn what the air is really like on this flight.")

    d.h1("2. The idea: learn a model of the world, then plan with it")
    d.p("Three earlier attempts in this project learned a <b>policy</b>, a function that maps the current state "
        "straight to a bank command. All three fell far short:")
    d.table([["Approach", "What went wrong", "Result"],
             ["Reinforcement learning (PPO) from scratch", "never discovered the bank-reversal pattern by trial and error", "~1% success"],
             ["Behavior cloning (copy the classical controller)", "small mistakes take it to states it never saw, where it is lost", "91 km median miss"],
             ["Cloning, then PPO", "the RL step destroyed the cloned skill", "438 km median miss"]],
            [2.4, 3.3, 1.3])
    d.p("The new approach is adapted from the SkyJEPA paper (Rao et al., 2026, quadrotors). It learns a <b>world "
        "model</b> instead: a function that answers <i>\"if I fly this sequence of bank angles from here, what "
        "happens?\"</i> A <b>planner</b> then uses it like a flight simulator in its head. Every 5 seconds it tries 256 "
        "candidate plans inside the model, and flies the one predicted to land closest to the target.")
    d.callout("In plain words: why a world model helps",
              ["It can be learned from <b>any</b> flights, even random ones, because it only has to predict what "
               "happens next. It needs no reward and no trial-and-error exploration.",
               "Mistakes do not snowball: the planner looks at the real state every 5 seconds and plans again.",
               "The target only appears in the planner's score, so the same model could fly to any target."])
    d.h2("What \"JEPA\" means")
    d.p("JEPA stands for <b>Joint-Embedding Predictive Architecture</b>. An <b>encoder</b> translates the recent "
        "history of the flight into a short list of numbers, the <b>latent code</b> (or just \"code\"). It is like a "
        "summary written in the model's own language. A <b>predictor</b> then learns to predict future codes from the "
        "current code and the planned bank commands. Training compares predicted codes with the codes the encoder "
        "produces for the real future (Figure 2). Because it never has to reproduce every raw number, the code is "
        "free to keep only what matters for the dynamics.")
    d.fig(F["jepa_idea"], "Two ways to build a model that predicts the future. Top: predict the raw state and feed "
                          "the prediction back in, so errors pile up. Bottom: JEPA predicts in the space of codes. A "
                          "separate prober later reads physical quantities out of the codes.", "illus")
    d.p("v1 tested the top option as a baseline (the \"predictive model\"). In closed loop it reached a 13 km median "
        "miss against 1.7 km for the JEPA model (Section 7).")

    d.h1("3. The pipeline, step by step")
    d.fig(F["pipeline"], "The four steps of the v2 pipeline. Each step's output feeds the next. The real simulator "
                         "appears twice: it makes the training data, and it flies the final tests.", "illus")
    d.numbered(["<b>Data.</b> Fly 6,000 random flights through the full simulator, each with its own randomized "
                "atmosphere and starting point.",
                "<b>Stage 1, the world model.</b> Train the encoder and predictor so that predicted codes match real "
                "future codes. No physical units are involved yet.",
                "<b>Stage 2, the prober.</b> Freeze the world model, then train a small \"prober\" that turns predicted "
                "codes into a physical trajectory (altitude, position, speed) by correcting a simple physics model.",
                "<b>Fly and test.</b> Plan with the model (MPPI), fly the frozen test cases in the real simulator, and "
                "compare with v1 and the classical controller."])

    d.h1("4. Step 1: the training data")
    d.p("Each training flight follows a <b>random bank profile</b>, of one of two kinds. <i>Smooth random</i> (\"GP\") "
        "profiles wander over all bank angles. <i>Reversal</i> profiles hold a slowly changing angle and flip its "
        "sign at random times, the pattern real guidance produces. Each profile is flown through the full "
        "simulator, so the data contains everything real flight has: the bank actuator, which can only roll 20° per "
        "second; the roll controller; 12 thruster jets; 3D equations of motion; and a realistic atmosphere table. "
        "One bank command is held for 5 seconds, a <b>macro step</b>, exactly as the planner will issue them.")
    d.fig(F["data_examples"], "Two training flights. Gray steps: the commanded bank angle, held for 5 s each. Blue: "
                              "what the capsule actually did. It lags and overshoots, because the actuator can only "
                              "roll 20° per second. The world model has to learn this lag.", "real")
    d.p("Every flight also gets its own <b>randomized truth</b>, so the model cannot memorize one atmosphere. It has to "
        "work out the conditions of each flight from what it observes:")
    d.bullets(["Entry point: altitude 100–140 km, speed 6.8–8.6 km/s, flight-path angle −1° to −6°, any heading.",
               "Air density ×0.7–1.3, lift-to-drag ratio (L/D) ×0.7–1.3, mass ×0.9–1.1. v1's data used only "
               "×0.85–1.15 and a single entry corridor.",
               "The flights do not aim at the target, and most land hundreds of km from it. That is fine: the model "
               "learns physics, and the target only matters to the planner."])
    d.fig(F["coverage"], f"The v2 data (orange, {I['n_wide']} flights shown) against v1's data (blue, "
                         f"{I['n_v1']} flights). Left: v1's flights all started in one small box. Middle: v2 visits "
                         "far more of the altitude–speed space. Right: v2's atmospheres cover all four test "
                         "atmospheres (dotted) with room to spare.", "real")
    d.p("In total: 6,000 flights, a median flight of about 6 minutes (75 macro steps), roughly 450,000 macro steps. "
        "The flights are split by flight: 80% for training, 10% for validation (used to choose settings), and 10% "
        "kept for testing.")

    d.h1("5. Steps 2 and 3: training in two stages")
    d.p("<b>Stage 1</b> trains the world model alone. Its only job is to make predicted codes match the real future "
        "codes, while a regularizer called <b>SIGReg</b> stops the codes from collapsing (Part 2, Section 14). "
        "<b>Stage 2</b> freezes the world model and trains the prober on top of it.")
    d.callout("In plain words: why two stages?",
              ["Splitting the work keeps the code general. Stage 1 learns \"how flights evolve\" without being told "
               "which quantities we care about. Stage 2 is then a small, cheap translation step. It also lets us "
               "test whether the code really contains information: if a prober can read the true air density out "
               "of the code, the world model has learned it."])

    d.h1("6. Step 4: flying with the model (MPPI)")
    d.p("<b>MPPI</b> (Model Predictive Path Integral control) is a simple, sturdy way to plan with a model. Every "
        "5 seconds:")
    d.numbered(["Take the last 2 minutes of measured history and give it to the world model.",
                "Make 256 candidate bank plans. Each is the previous best plan plus a smooth random change.",
                "Fly all 256 plans inside the model down to parachute altitude, in one batch.",
                "Score each plan. The score (cost, in km) is the predicted miss distance, plus small penalties for "
                "jerky commands, too many g, and diving too steeply. In v2 it also includes how much the ensemble's "
                "models disagree.",
                "Average the plans, giving cheaper plans much more weight. Fly the first 5 seconds of the result, "
                "then repeat."])
    d.fig(F["mppi_snapshot"], f"One real planning step of the shipped v1 model, 150 s into frozen test case 0 (normal "
                              f"air). Left: the 256 candidate bank plans, darker = lower cost, with the weighted "
                              f"average dashed. Right: where the model predicts each plan lands. At this moment the "
                              f"best plan is predicted to miss by {mp['miss_best']:.0f} km, but the planner keeps "
                              f"replanning: the same flight, continued, opened its parachute "
                              f"{mp['final_miss']:.2f} km from the target at t = {mp['final_t']:.0f} s.", "real")
    d.fig(F["mppi_loop"], "The v2 planning loop. The difference from v1 is the ensemble: every plan is flown by 5 "
                          "independently trained world models, and their disagreement is part of the cost.", "illus")
    d.callout("In plain words: why there are no hand-written rules",
              ["The classical controller has rules such as \"reverse the bank when the heading error is too large\". MPPI "
               "has none. Bank reversals appear by themselves, because the plans that reverse at the right time "
               "are the ones predicted to land on target."])

    d.h1("7. What v1 showed")
    d.p("v1 (26 thousand parameters, trained on 320 narrow-corridor flights) was the first test of the whole idea. "
        f"Figure {N['v1_results']} shows its main results.")
    d.fig(F["v1_results"], "Left: open-loop prediction of the landing point from different times before the "
                           "parachute. The JEPA model plus prober beats the predictive baseline and the physics prior "
                           "alone, and is close to the classical controller's own high-fidelity model. Right: closed "
                           "loop on the first 16 frozen test cases; bar labels are success rates.", "v1")
    d.p("<b>The good news:</b> with 26 thousand parameters, v1 matched the classical controller in normal air "
        "(1.69 km against 2.24 km median miss, both 100% success). It was far better than every earlier learned "
        "controller, and better than planning with the predictive baseline.")
    if dense:
        cl, v1m = dense["summary"]["classical"]["median"], dense["summary"]["jepa"]["median"]
        v1c, clc = dense["per_case"]["jepa"], dense["per_case"]["classical"]
        big = [m for m in v1c if m > 5.0]
        small = [m for m in v1c if m <= 5.0]
        rest = ("" if not small else f" (in the remaining case, {small[0]:.1f} km)" if len(small) == 1 else
                f" (in the other {len(small)} cases, {', '.join(f'{m:.1f}' for m in small)} km)")
        d.p(f"<b>The bad news: dense air.</b> For this report, the shipped v1 model and the classical controller were "
            f"flown again on the first 8 frozen test cases with density and L/D both 12% above nominal (Figure {N['dense']}). "
            f"The classical controller's median miss was {cl:.1f} km (every case between {min(clc):.1f} and "
            f"{max(clc):.1f} km); v1's was {v1m:.1f} km. Every v1 flight still landed inside the 25 km circle, but in "
            f"{len(big)} of 8 cases it landed {min(big):.0f}–{max(big):.0f} km off{rest}. Missing by about the same "
            f"amount again and again is a <b>bias</b>: the model gets something about dense air systematically wrong, "
            f"most likely how much drag there is.")
        d.fig(F["dense"], "The shipped v1 model against the classical controller, density and L/D ×1.12, 8 frozen test "
                          f"cases, v1's own settings (256 MPPI samples). Right: the end of the ground track for case "
                          f"{dense['worst']}, v1's worst.", "real")
    d.p("Other things v1 left open:")
    d.bullets(["<b>No clear sign of system identification.</b> On flights with unusual atmospheres, v1's predictions were "
               "no better than the nominal physics model's. So there was no evidence that it inferred the air density "
               "from its history, which is the whole point of a learned model.",
               "<b>Short training horizon.</b> v1 learned to predict 24 steps (2 minutes) ahead, but plans 110 steps "
               "(9 minutes) ahead: 4.6 times further than it was ever trained on.",
               "<b>Narrow data.</b> All 320 flights started in one corridor.",
               "<b>Not g-matched.</b> v1 pulled about 5.6 g peak against the classical controller's 3.7 g. Part of its "
               "accuracy may have been bought with load."])

    d.h1("8. The v2 experiment")
    d.h2("8.1 What changes")
    d.p("Part 2 explains every change in detail. In one table:")
    d.table([["", "v1", "v2", "Aimed at"],
             ["Data", "320 flights, one corridor", "6,000 flights, wide envelope", "narrow data, dense air"],
             ["Inputs", "10 measured quantities", "+ the drag and L/D corrections measured directly", "dense-air bias"],
             ["History", "35 s (TCN)", "120 s (attention)", "system identification"],
             ["Predictor", "GRU (one step at a time)", "transformer (all steps at once)", "long horizons"],
             ["Training horizon", "24 steps", "24 → 48 → 96 steps", "extrapolation"],
             ["Prober", "sees the code", "sees the code + the rollout's altitude and speed", "dense-air bias"],
             ["Models", "1", "ensemble of 5, disagreement used by the planner", "risky plans"],
             ["Parameters", "32,675", f"{rep['ensemble (x5)']:,}", ""]],
            [1.15, 1.7, 2.45, 1.7])
    d.h2("8.2 What the notebook runs")
    d.table([["Notebook section", "What happens", "Est. time (L4 GPU)"],
             ["0. Setup", "clone the jepa-v2 branch, install pymsis, mount Google Drive, check the GPU", "2 min"],
             ["1. Data", "load the 6,000 flights from Drive, split 80/10/10, plot coverage", "1–2 min"],
             ["2. Model", "build the 5-member ensemble, print sizes and a time estimate", "seconds"],
             ["3. Background", "start the classical-controller flights and the v1-wide retraining on the CPUs", "runs in parallel"],
             ["4. λ sweep", "try 5 strengths of the SIGReg regularizer; keep the best one", "~5 min"],
             ["5. Stage 1", "train the ensemble world model (24,000 steps)", "~30 min"],
             ["6. Stage 2", "train the 5 probers (6,000 steps) and save the final model", "~10 min"],
             ["7. Open-loop tests", "landing-point error, dense-air bias, calibration, system ID", "~5 min"],
             ["8. Closed loop", "fly 16 test cases × 4 atmospheres × 5 controllers", "~15 min"],
             ["9. Verdict", "score the 6 pass/fail criteria", "seconds"]],
            [1.35, 4.25, 1.4])
    d.h2("8.3 Controls: how we will know what helped")
    d.p("A bigger model trained on more data could be better for many reasons. The notebook includes controls that "
        "separate them:")
    d.bullets(["<b>v1-wide</b>: v1's exact small architecture, retrained on the same 6,000 flights. If v2 beats v1-wide, "
               "the gain comes from the new model, not only from the new data.",
               "<b>v2 (no unc.)</b>: the same ensemble flown without the disagreement penalty. It shows what the "
               "uncertainty term buys.",
               "<b>physics-nominal</b>: the classical controller's own model, used only for prediction. It is the "
               "reference for open-loop accuracy.",
               "<b>Physics prior only</b>: the prober's crude physics with all learned corrections set to zero. It shows "
               "what the learning adds.",
               "<b>λ sweep</b>: the regularizer strength is chosen on validation flights, never on test flights."])
    d.h2("8.4 The tests")
    d.p("<b>Open loop</b> (no controller; the true future commands are fed in):")
    d.bullets(["<b>Landing-point error vs. time-to-go</b>: from 400, 300, 200, 100 and 50 s before the parachute, how "
               "far is the predicted parachute point from the real one? This is the number guidance depends on.",
               "<b>Dense-air bias</b>: the same error split into along-track and cross-track parts and binned by true "
               "density. A bias shows up as a line that rises with density.",
               "<b>v1's home turf</b>: v1's own test flights, to check that the wider model has not become worse on "
               "v1's narrow corridor.",
               "<b>Calibration</b>: does the ensemble's disagreement grow when its error grows?",
               "<b>System identification probe</b>: can a straight-line (linear) readout of the code recover the true "
               "density and L/D, and how quickly after entry?"])
    d.p("<b>Closed loop</b>: every controller flies the same 16 frozen test cases under four atmospheres:")
    d.table([["Atmosphere", "Density", "L/D", "Why"],
             ["nominal", "×1.00", "×1.00", "the classical controller's own model is exactly right"],
             ["mismatch", "×0.88", "×0.88", "thinner air"],
             ["mismatch_high", "×1.12", "×1.12", "denser air: v1's weak spot"],
             ["out_of_range", "×0.78", "×0.80", "outside v1's training range (inside v2's)"]],
            [1.3, 0.8, 0.8, 4.1])
    d.h2("8.5 Pass/fail criteria, fixed before the run")
    d.p("Writing down what counts as success <i>before</i> seeing the results stops us from moving the goalposts "
        "afterwards. The notebook's last cell scores these automatically:")
    d.table([["#", "Question", "Passes if"],
             ["1", "Did the codes stay healthy (no collapse)?", "every member's effective rank > 16 of 64"],
             ["2", "Is it the model, not just the data?", "v2's median landing error < v1-wide's at 200 s and 300 s"],
             ["3", "Does it beat the classical controller's model?", "v2 < physics-nominal at 200 s and 300 s"],
             ["4", "Is the dense-air bias gone?", "mean along-track error at 200 s < 2 km, and smaller than v1's"],
             ["5", "Does it fly robustly?", "≥ 95% success in all 4 atmospheres, and median miss in dense air ≤ classical's"],
             ["6", "Is the uncertainty useful?", "spread–error rank correlation > 0.4, and the p90 miss is no worse "
                                                 "than without the penalty in ≥ 3 of 4 atmospheres"]],
            [0.3, 2.7, 4.0])
    d.p("Sixteen cases is still a small sample. Differences under about 1 km in the median should be treated as noise. "
        "Before claiming a win over the classical controller, rerun with 48 cases.")
    d.h2("8.6 Time, cost, and which GPU")
    d.fig(F["compute"], "Estimated wall time and cost of the full notebook. The estimates come from the exact "
                        "training compute (33 PFLOP for stage 1) at an assumed ~15% of each GPU's peak speed, "
                        "plus measured CPU costs for the simulator. Colab compute-unit rates are approximate and "
                        "change; the runtime picker shows the live rate.", "est")
    d.p("<b>Recommendation: the L4 GPU.</b> It costs about the same as a T4 but finishes in half the time. It "
        "supports the faster bf16 number format and has enough memory to skip the memory-saving mode. Its larger "
        "runtime also has more CPUs for the simulator work. Use a T4 only on the free tier. Everything is "
        "resumable, so a disconnect costs at most 10 minutes of training. The notebook prints the measured speed and "
        "the time remaining after its first few hundred training steps.")
    d.h2("8.7 Where the results are stored")
    d.table([["Location on Google Drive", "What it holds"],
             ["MyDrive/ReEntryAI_jepa_data/wide_6000/", "the input: 6,000 flights (written by the data notebook)"],
             ["MyDrive/ReEntryAI_jepa_v2/full/v2_ckpt.pt", "the final v2 model (all 5 members and their probers)"],
             ["…/full/v1_wide_ckpt.pt", "the v1-wide control model"],
             ["…/full/wm_stage.pt, pi_stage.pt, sweep/", "training checkpoints, used to resume after a disconnect"],
             ["…/full/*.csv", "results: λ sweep, landing-point errors, closed-loop flights and summary, the verdict"],
             ["MyDrive/ReEntryAI_jepa_v2/closed_loop_cache_full.pkl", "every closed-loop flight, saved one by one"]],
            [3.3, 3.7])
    d.p("The figures are shown only inside the notebook. Use <i>File → Save a copy in Drive</i> after the run to "
        "keep them.")
    d.h2("8.8 What has already been checked")
    d.p("This machine has no GPU, so v2 was checked on a CPU with a tiny version of the model and the first 160 "
        "real flights of the dataset. The numbers from that run say nothing about v2's real accuracy, but the "
        "checks show the code is right:")
    d.bullets(["The encoder gives the same code whether it reads a long sequence or a single 2-minute window "
               "(difference 10<super>−11</super>).",
               "Running the predictor one step at a time (as v1's code expects) gives exactly the same results as "
               "running it in one parallel pass.",
               "v2's regularizer gives the same value as v1's on the same data.",
               "Both training stages lower their validation loss; the prober beats the physics prior alone.",
               "Checkpoints save, load, and resume correctly.",
               "The new GPU flight runner makes the same bank decisions as v1's unchanged runner on the same model "
               "(within 0.001°).",
               "The whole notebook runs from start to finish in its small \"smoke\" mode (under 3 minutes on 4 CPUs)."])

    # ================================================================ PART 2
    d.part("Part 2: the models and the architecture",
           "Every block of v2: what it does, why it was chosen, and how it differs from v1")

    d.h1("9. The big picture")
    d.p(f"v2 is an <b>ensemble</b> of 5 identical-looking but independently trained models. Each one has three "
        f"parts: a state encoder, an action encoder, and a predictor (Figure {N['architecture']}). A fourth part, the prober, is added "
        f"in stage 2. Each member has {rep['per member']:,} parameters; the ensemble has {rep['ensemble (x5)']:,}.")
    d.fig(F["architecture"], "One v2 ensemble member during stage 1. The encoder turns 2 minutes of history into a "
                             "64-number code. The predictor reads that code plus one token per planned bank "
                             "command, and outputs a predicted code for every future step at once. The loss "
                             "compares them with the real future codes. SIGReg keeps the codes spread out.", "illus")

    d.h1("10. What the model sees: inputs")
    d.p("Every 5 seconds the model receives 15 numbers about the capsule. Ten are the same as v1's; five are new:")
    d.table([["Input", "Meaning", "New in v2?"],
             ["altitude, speed, flight-path angle γ", "where the capsule is and how it moves", ""],
             ["sin and cos of heading, latitude", "direction of flight and position on Earth", ""],
             ["realized bank angle, roll rate", "what the capsule is actually doing (lags the command)", ""],
             ["drag and lift acceleration (in g)", "what the accelerometer measures", ""],
             ["sin and cos of the bank angle", "the bank's effect on lift is through sin and cos", "yes"],
             ["log of drag acceleration", "drag grows more than 10,000× between 120 km and 30 km; the log tames that range", "yes"],
             ["<b>measured drag correction</b> c<sub>D,meas</sub>", "how much more or less drag than the simple physics prior predicts, right now", "<b>yes</b>"],
             ["<b>measured L/D correction</b> c<sub>L,meas</sub>", "how much more or less lift per drag than the prior's 0.30", "<b>yes</b>"]],
            [2.3, 3.8, 0.9])
    d.h2("Physics-residual features: reading the air density off the accelerometer")
    d.p("Drag depends on air density, speed squared, and the capsule's shape and mass. The prober (Section 16) "
        "starts from a deliberately simple physics formula and learns corrections to it. v2 computes those "
        "corrections directly from the measurements: c<sub>D,meas</sub> = log(measured drag ÷ drag the simple "
        "formula predicts at this altitude and speed). If the air is 20% denser than the formula assumes, this "
        "number moves by about log 1.2 ≈ 0.18.")
    d.fig(F["features"], f"Left: the measured drag correction along {I['n_wide']} real training flights, colored by "
                         f"their true density (÷ mass) factor. The color bands separate cleanly: the feature shifts "
                         f"with density. The bends with altitude come from where the real atmosphere differs from "
                         f"the simple formula. Right: the median of the feature between 35 and 60 km against the "
                         f"true factor. The slope is {feat['slope']:.2f} and the correlation r = {feat['r']:.2f}, "
                         f"so the feature reads the truth almost exactly.", "real")
    d.callout("In plain words: why this targets the dense-air bias",
              ["v1 had to work out \"how much denser is the air than expected\" by combining altitude, speed, and "
               "drag inside a small network, a multiply-and-divide over a 10,000× range. A small, systematic slip in "
               "that calculation is exactly what a dense-air bias looks like. v2 hands the network the answer to "
               "that calculation. The network still has to decide how to use it, but it no longer has to rediscover "
               "the physics."])

    d.h1("11. The state encoder: attention over 2 minutes")
    d.p("v1's encoder was a small <b>TCN</b> (temporal convolution network) that looked at the last 7 steps (35 s). v2 "
        "uses a <b>transformer</b>, the building block behind modern language models. Its key operation is "
        "<b>attention</b>: each time step looks at other time steps and decides how much each one matters. That is "
        "useful for system identification. To judge how dense the air is, the model can compare drag a minute ago "
        "with drag now, at different altitudes and speeds.")
    d.p("v2's attention is <b>banded and causal</b>. Each step may only look at itself and the 8 steps before it "
        f"(causal: never the future). Stacking 3 such layers gives a window of exactly 25 steps = 120 s (Figure {N['masks']}). "
        "Why insist on an exact window? The same encoder makes the predictor's starting code and its training targets. "
        "With an exact window, the code for \"now\" is the same whether the encoder reads a long training sequence "
        "or only the last 2 minutes in flight, so training and flying match.")
    d.fig(F["masks"], "Attention patterns computed by the v2 code (colored = allowed to look). Left: one encoder "
                      "layer looks back 8 steps. Middle: after 3 layers, each step depends on exactly the last 25 "
                      "steps. Right: the predictor's causal pattern: the output for future step k sees the starting "
                      "code and the actions up to step k, never later ones.", "code")
    d.p("Positions are encoded with <b>rotary embeddings</b> (RoPE). They tell attention how far apart two steps "
        "are, not their absolute position, so the encoder treats every 2-minute window the same way. The encoder "
        "has 3 layers, each 128 numbers wide with 4 attention heads, and outputs a code of 64 numbers (v1: 32).")

    d.h1("12. Action tokens: packing the recent commands")
    d.p("The predictor needs to know the planned bank commands. v2 turns each 5-second step into one <b>action "
        "token</b>, built from that step's command and the 3 before it (20 seconds). Each command is described by its "
        "value and its sine and cosine. Why 20 seconds? Reversing the bank from +90° to −90° at the actuator's "
        "20°/s limit takes 9 seconds, so 20 s of commands contains everything the actuator lag depends on.")
    d.fig(F["action_tokens"], "Each action token packs the current command and the 3 before it. Tokens overlap, the "
                              "way a sliding window does.", "illus")

    d.h1("13. The predictor: a causal transformer")
    d.p("v1's predictor was a <b>GRU</b>, a recurrent network that steps forward one 5-second step at a time, so "
        "training on 96 steps means 96 steps in a row. v2's predictor is a <b>causal transformer</b>: it reads the "
        f"starting code and all action tokens together and outputs every future code in one pass (Figure {N['parallel']}). "
        "\"Causal\" means each output may only look at earlier tokens. So the prediction for step k depends only on "
        "the starting code and the commands up to step k, the same rule the GRU followed by construction.")
    d.fig(F["parallel"], "v1 predicts one step after another. v2 predicts all steps in one pass. That makes long "
                         "training horizons affordable on a GPU, and lets the planner score all 256 plans × 110 "
                         "steps × 5 models at once.", "illus")
    d.bullets(["Size: 4 layers, 192 numbers wide, 6 attention heads, about 1.8 million parameters.",
               "The output is the starting code <i>plus</i> a learned change, and the change starts at zero. So before "
               "training, the model's guess is \"nothing changes\", a calm starting point.",
               "For compatibility with v1's evaluation and planning code, the predictor can also run one step at a "
               "time, remembering earlier steps (a \"KV cache\"). Both modes give identical results."])

    d.h1("14. The training objective: prediction + SIGReg")
    d.p("The main loss is simple: the average squared difference between each predicted code and the real future "
        "code. But that loss alone has a cheating solution. If the encoder outputs <i>the same code for every "
        "state</i>, prediction is trivially perfect and the model has learned nothing. This is called <b>collapse</b>.")
    d.p("<b>SIGReg</b> (Sketched Isotropic Gaussian Regularization, from LeJEPA) prevents it. It pushes the codes "
        "across many flights to spread out like a round bell curve (a standard normal distribution) in every "
        "direction. To check this cheaply, it projects the codes onto 32 random directions and measures how far each "
        "1-D projection is from a bell curve. That uses the Epps-Pulley test, which compares the data's "
        "\"characteristic function\" with the bell curve's.")
    d.fig(F["sigreg"], f"Made-up 2-D codes scored by the real SIGReg code. Collapsed codes are charged heavily "
                       f"(point: {sg['collapsed to a point']:.0f}, line: {sg['collapsed to a line']:.0f}). Healthy, "
                       f"bell-shaped codes cost almost nothing ({sg['healthy: round bell curve']:.2f}). The scale "
                       f"is logarithmic.", "illus")
    d.bullets(["<b>Per step:</b> v2 applies SIGReg separately at each prediction step, as the SkyJEPA paper does. v1 "
               "pooled all steps together.",
               "<b>How strong?</b> The total loss is prediction + λ × SIGReg. Too small a λ allows partial collapse; too "
               "large a λ forces spread at the expense of prediction. The notebook trains 5 members with λ = 0.01, "
               "0.02, 0.04, 0.08 and 0.16 side by side. It keeps the one whose prober predicts landing points best on "
               "validation flights.",
               "Like v1 and LeJEPA, v2 needs no other tricks against collapse: no stop-gradient, no slowly-updated "
               "target network, no decoder."])

    d.h1("15. Long horizons: curriculum and masked windows")
    d.p("The planner looks up to 110 steps (9 minutes) ahead. v1 trained on 24-step windows, so it had to extrapolate "
        "4.6× beyond its training. v2 trains on longer and longer windows: 24 steps first, then 48, then 96 (a "
        f"<b>curriculum</b>, Figure {N['schedule']}). Starting short makes early training easy and fast; ending long matches how "
        "the model is used. The planner then extrapolates only 1.15× beyond training.")
    d.fig(F["schedule"], "The training schedules used in the full notebook. Left: prediction horizon during "
                         "training, for the world model (stage 1) and the prober (stage 2). Right: learning rate: "
                         "a short warm-up, then a smooth cosine decay.", "code")
    d.p(f"A long window raises a practical problem: most wide-envelope flights are shorter than 96 steps (Figure {N['masked_windows']}). "
        "Normally such flights would simply be left out of long windows, and the long-horizon training would then "
        "only see long, gliding flights, a biased sample. v2 instead lets windows run past the landing and switches "
        "the missing steps off with a <b>mask</b>.")
    d.fig(F["masked_windows"], "Left: flight lengths in the v2 data (first 160 real flights) and v1's data, with the "
                               "three training horizons. Right: a window that runs past the landing. Only the real "
                               "steps count in the loss. (Left: real data; right: illustration.)", "real")

    d.h1("16. The physics-inspired prober")
    d.p("The prober turns predicted codes into a physical trajectory. It does not predict altitude or position "
        "directly. Instead it corrects a deliberately <b>crude physics prior</b>: one exponential atmosphere, one "
        "ballistic coefficient, L/D = 0.30, and lift pointed along the commanded bank. At each step a small network "
        "outputs three corrections: a drag factor, an L/D factor, and a bank offset (for the actuator lag). "
        f"The corrected forces then go through the <b>exact 3D equations of motion</b> (Figure {N['prober']}). This was v1's "
        "single biggest win: geometry and gravity come for free, and the network only learns what is hard to model.")
    d.fig(F["prober"], "The stage-2 prober. New in v2: besides the predicted code, the small network sees where "
                       "the rollout actually is (altitude h and speed V). Each step's output feeds the next step.",
          "illus")
    d.bullets(["<b>Why add altitude and speed?</b> The difference between the real atmosphere and the exponential "
               "formula depends on altitude, and drag and lift depend on Mach number (speed). Those are properties "
               "of where the <i>predicted trajectory</i> is. In dense air the capsule slows down higher up; if the "
               "prediction drifts, a correction that ignores the rollout's own altitude is applied in the wrong "
               f"place. Figure {N['features']} (left) shows how strongly the needed correction varies with altitude.",
               "<b>Starts as pure physics:</b> the network's last layer starts at zero, so before training the prober "
               "<i>is</i> the crude prior. That gives the free baseline \"physics prior only\".",
               "<b>Robust loss:</b> errors are measured in friendly units (1 km of altitude, 5 km of position, 50 m/s, "
               "half a degree) through a pseudo-Huber function. It behaves like a squared error for small errors and "
               "grows only linearly for large ones, so a few bad long-horizon windows cannot dominate training.",
               "<b>Always full precision:</b> the physics runs in 32-bit floating point even when the networks use "
               "faster 16-bit math on the GPU."])

    d.h1("17. The ensemble and uncertainty")
    d.p("A single model can be confidently wrong, especially in situations it rarely saw. v2 trains 5 members that "
        "differ only in their random starting weights and the random order of their training data. Where the data "
        "was plentiful they tend to agree; where it was thin they tend to disagree. Their disagreement is a "
        "practical estimate of <b>epistemic uncertainty</b>, the uncertainty that comes from limited knowledge.")
    d.fig(F["ensemble"], "How MPPI uses the ensemble (made-up numbers). Each plan is flown by all 5 models. Plan B's "
                         "average landing point is fine, but the models disagree widely about it. The cost adds the "
                         "average miss and the spread, so MPPI prefers plan A, the plan the models agree on.", "illus")
    d.bullets(["Cost of a plan = average over the 5 models of v1's cost, + 1 km per km of spread in the predicted "
               "landing points.",
               "The 5 members are stored as one network whose weights have an extra \"member\" dimension. Training and "
               "planning therefore run all 5 at once, as one set of large matrix multiplications.",
               "Whether the spread is trustworthy is tested, not assumed: the calibration test (Section 8.4) and the "
               "\"v2 (no unc.)\" control (Section 8.3)."])

    d.h1("18. Size and compute")
    d.fig(F["params"], "Left: parameters per component, v1 against one v2 member (log scale). Right: training "
                       "compute per step as the prediction horizon grows. Longer horizons cost more, which is one "
                       "reason for the curriculum.", "code")
    d.bullets([f"<b>{rep['per member']:,} parameters per member, {rep['ensemble (x5)']:,} in total</b>, about 375× v1. "
               "Most of them are in the predictor.",
               "<b>Mixed precision:</b> the networks compute in 16-bit (bf16 on L4/A100, fp16 on T4) for speed. Codes, "
               "losses, and physics stay in 32-bit.",
               "<b>Memory-saving mode</b> (activation checkpointing) switches on automatically on GPUs under 20 GB, "
               "such as the T4. It recomputes some values instead of storing them: about 30% slower, about 10× less memory.",
               "<b>Resumable:</b> every 10 minutes the model, the optimizer, and the training history are saved to Drive."])

    d.h1("19. Glossary")
    d.table([["Term", "Meaning"],
             ["bank angle σ", "how far the capsule is rolled around its direction of flight; the only steering control"],
             ["drogue parachute", "the first, small parachute; it opens at 7.6 km altitude and marks the end of guidance"],
             ["L/D", "lift-to-drag ratio; about 0.3 for Orion"],
             ["classical predictor-corrector", "the hand-built guidance baseline; it re-simulates the flight with a fixed nominal model every second"],
             ["world model", "a model that predicts what happens if a given sequence of actions is applied"],
             ["policy", "a model that maps the current state straight to an action (what RL and cloning learned)"],
             ["JEPA", "Joint-Embedding Predictive Architecture: learns to predict future codes rather than raw data"],
             ["latent code", "the short list of numbers (64 in v2) the encoder uses to summarize the recent flight"],
             ["encoder / predictor / prober", "history → code / code + plan → future codes / codes → physical trajectory"],
             ["collapse", "the failure where every input maps to the same code; prediction becomes trivially perfect and useless"],
             ["SIGReg", "a regularizer that keeps codes spread like a round bell curve, which prevents collapse"],
             ["λ (lambda)", "how strongly SIGReg is weighted against the prediction loss"],
             ["attention / transformer", "a network layer in which every position weighs information from other positions"],
             ["causal", "may only use the past, never the future"],
             ["curriculum", "training on easy (short) problems first, then harder (longer) ones"],
             ["mask", "a 1/0 flag that switches individual loss terms on or off"],
             ["ensemble", "several independently trained models used together; their disagreement signals uncertainty"],
             ["MPPI", "Model Predictive Path Integral control: sample many plans, score them in the model, average the good ones, repeat"],
             ["macro step", "5 seconds; one bank command is held for one macro step"],
             ["open loop / closed loop", "prediction with the true commands fed in / a controller actually flying the simulator"],
             ["system identification", "inferring the unknown conditions (density, L/D) from the observed flight"],
             ["effective rank", "roughly how many independent directions the codes really use (64 = all of them)"]],
            [1.9, 5.1])
    d.p("Code: <font name='DVM'>jepa/models_v2.py</font> (architecture), <font name='DVM'>jepa/train_v2.py</font> "
        "(training), <font name='DVM'>jepa/mppi_v2.py</font> (planner), <font name='DVM'>jepa/evaluate_v2.py</font> "
        "(diagnostics), <font name='DVM'>ReEntryAI_JEPA_v2.ipynb</font> (the experiment). Design notes: "
        "<font name='DVM'>jepa/PLAN.md</font> §11. This report: <font name='DVM'>jepa/make_report_v2.py</font>.")

    doc = SimpleDocTemplate(out_pdf, pagesize=letter, leftMargin=0.75 * inch, rightMargin=0.75 * inch,
                            topMargin=0.75 * inch, bottomMargin=0.75 * inch, title="SpaceJEPA v2: a beginner's guide",
                            author="ReEntryAI", subject="JEPA world model and MPPI for Orion entry guidance")
    doc.build(d.story, onFirstPage=_page, onLaterPages=_page)


if __name__ == "__main__":
    from jepa.report_v2_figures import make_all
    ap = argparse.ArgumentParser()
    ap.add_argument("--wide", required=True, help="wide-envelope episodes: pickle or shard directory")
    ap.add_argument("--dense", default="", help="pickle of dense-air closed-loop results (optional)")
    ap.add_argument("--run-dense", action="store_true", help="fly the dense-air cases if --dense does not exist")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--figs", default="", help="figure directory (default: a temp dir)")
    ap.add_argument("--out", default=str(SIM_DIR / "jepa" / "SpaceJEPA_v2_report.pdf"))
    a = ap.parse_args()
    if a.run_dense and a.dense and not os.path.exists(a.dense):
        import pickle
        from jepa.report_v2_figures import run_dense_cases
        with open(a.dense, "wb") as f:
            pickle.dump(run_dense_cases(workers=a.workers), f)
    r = make_all(a.wide, a.dense, a.figs or tempfile.mkdtemp(prefix="jepa_report_"))
    build(r["figs"], r["info"], a.out)
    print(f"wrote {a.out}")
