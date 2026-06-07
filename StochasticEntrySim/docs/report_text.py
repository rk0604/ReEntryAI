"""Text and layout for the ReEntryAI technical overview PDF (reportlab)."""
from reportlab.lib.pagesizes import LETTER
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.platypus import (
    BaseDocTemplate, PageTemplate, Frame, Paragraph, Spacer, Image, Table,
    TableStyle, PageBreak, NextPageTemplate, KeepTogether)

NAVY = colors.HexColor("#16294d")
BLUE = colors.HexColor("#2d6cdf")
ORANGE = colors.HexColor("#e8833a")
GREEN = colors.HexColor("#2e9e5b")
PURPLE = colors.HexColor("#7b5cd6")
MUTED = colors.HexColor("#5b6b7c")
LGRAY = colors.HexColor("#f4f6fa")
TEXT = colors.HexColor("#16202b")

S = getSampleStyleSheet()
H1 = ParagraphStyle("H1", parent=S["Heading1"], textColor=NAVY, fontSize=20,
                    spaceBefore=10, spaceAfter=6, leading=24)
H2 = ParagraphStyle("H2", parent=S["Heading2"], textColor=BLUE, fontSize=13.5,
                    spaceBefore=12, spaceAfter=3, leading=17)
BODY = ParagraphStyle("Body", parent=S["BodyText"], textColor=TEXT, fontSize=10.7,
                      leading=15.5, spaceAfter=6, alignment=TA_LEFT)
CAP = ParagraphStyle("Cap", parent=S["BodyText"], textColor=MUTED, fontSize=8.8,
                     leading=11, alignment=TA_CENTER, spaceBefore=2, spaceAfter=12)
PART = ParagraphStyle("Part", parent=S["Heading1"], textColor=colors.white,
                      fontSize=26, leading=30, alignment=TA_LEFT)
COVER_T = ParagraphStyle("CoverT", parent=S["Title"], textColor=NAVY, fontSize=30,
                         leading=35, alignment=TA_CENTER)
COVER_S = ParagraphStyle("CoverS", parent=S["Title"], textColor=MUTED, fontSize=14,
                         leading=20, alignment=TA_CENTER)
BULLET = ParagraphStyle("Bul", parent=BODY, leftIndent=14, bulletIndent=2, spaceAfter=3)


def P(t):
    return Paragraph(t, BODY)


def bul(t):
    return Paragraph(t, BULLET, bulletText="•")


def fig(path, w_in, cap):
    from PIL import Image as PILImage
    iw, ih = PILImage.open(path).size
    w = w_in * inch
    h = w * ih / iw
    return [Image(path, width=w, height=h), Paragraph(cap, CAP)]


def kv_table(rows, col_w):
    t = Table(rows, colWidths=[c * inch for c in col_w])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), NAVY),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 9.5),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, LGRAY]),
        ("TEXTCOLOR", (0, 1), (-1, -1), TEXT),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#cdd6e4")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    return t


def _cover(c, d):
    c.saveState()
    c.setFillColor(NAVY); c.rect(0, 0, LETTER[0], LETTER[1], fill=1, stroke=0)
    c.setFillColor(BLUE); c.rect(0, LETTER[1] - 0.35 * inch, LETTER[0], 0.35 * inch, fill=1, stroke=0)
    c.setFillColor(ORANGE); c.rect(0, 0, LETTER[0], 0.2 * inch, fill=1, stroke=0)
    c.restoreState()


def _part(c, d, label):
    c.saveState()
    c.setFillColor(NAVY); c.rect(0, 0, LETTER[0], LETTER[1], fill=1, stroke=0)
    c.setFillColor(ORANGE); c.rect(0.9 * inch, LETTER[1] - 3.7 * inch, 1.6 * inch, 0.12 * inch, fill=1, stroke=0)
    c.restoreState()


def _body(c, d):
    c.saveState()
    c.setStrokeColor(colors.HexColor("#cdd6e4")); c.setLineWidth(0.6)
    c.line(0.9 * inch, 0.72 * inch, LETTER[0] - 0.9 * inch, 0.72 * inch)
    c.setFont("Helvetica", 8); c.setFillColor(MUTED)
    c.drawString(0.9 * inch, 0.55 * inch, "ReEntryAI  Technical Overview")
    c.drawRightString(LETTER[0] - 0.9 * inch, 0.55 * inch, "Page %d" % d.page)
    c.restoreState()


def build_pdf(figs, out_path):
    doc = BaseDocTemplate(out_path, pagesize=LETTER,
                          leftMargin=0.9 * inch, rightMargin=0.9 * inch,
                          topMargin=0.85 * inch, bottomMargin=0.9 * inch)
    fw = LETTER[0] - 1.8 * inch
    fh = LETTER[1] - 1.75 * inch
    frame = Frame(0.9 * inch, 0.9 * inch, fw, fh, id="body")
    pframe = Frame(0.9 * inch, 0.9 * inch, fw, fh, id="part")
    doc.addPageTemplates([
        PageTemplate(id="cover", frames=[Frame(0.9 * inch, 1.2 * inch, fw, fh - 0.3 * inch)], onPage=_cover),
        PageTemplate(id="part1", frames=[pframe], onPage=lambda c, d: _part(c, d, "1")),
        PageTemplate(id="part2", frames=[pframe], onPage=lambda c, d: _part(c, d, "2")),
        PageTemplate(id="body", frames=[frame], onPage=_body),
    ])
    s = []

    # ---------- COVER ----------
    s += [Spacer(1, 1.6 * inch)]
    s += [Paragraph("ReEntryAI", ParagraphStyle("ct", parent=COVER_T, textColor=colors.white, fontSize=40))]
    s += [Spacer(1, 0.2 * inch)]
    s += [Paragraph("A Simulator and Reinforcement Learning Framework<br/>for Orion Crew Module Atmospheric Entry Guidance",
                    ParagraphStyle("cs", parent=COVER_S, textColor=colors.HexColor("#bcd0f5"), fontSize=15))]
    s += [Spacer(1, 0.5 * inch)]
    s += [Paragraph("Technical Overview", ParagraphStyle("cs2", parent=COVER_S, textColor=ORANGE, fontSize=13))]
    s += [NextPageTemplate("body"), PageBreak()]

    # ---------- CONTENTS ----------
    s += [Paragraph("What this document covers", H1)]
    s += [P("This document explains how the project works in two parts. The first part describes the "
            "simulator, the physics model that flies a capsule back through the atmosphere. The second "
            "part describes the reinforcement learning system, which trains a neural network to steer the "
            "capsule and adds guaranteed safety on top.")]
    s += [Spacer(1, 6)]
    s += [kv_table([
        ["Part", "Topic"],
        ["Part 1", "The Simulator: the physics, the heating, the controls, and the parachutes"],
        ["Part 2", "The Reinforcement Learning System: the learned controller and certified safety"],
    ], [0.9, 5.0])]
    s += [Spacer(1, 10)]
    s += [P("Throughout, the goal is the same. A capsule arrives very fast and very high. It must slow "
            "down, survive intense heat and strong forces, and arrive at a precise landing site. The "
            "simulator models all of this. The learning system tries to fly it better than traditional methods.")]

    # =========================================================================
    # PART 1
    # =========================================================================
    s += [NextPageTemplate("part1"), PageBreak()]
    s += [Spacer(1, 2.4 * inch)]
    s += [Paragraph("Part 1", ParagraphStyle("p1", parent=PART, textColor=ORANGE))]
    s += [Paragraph("The Simulator", ParagraphStyle("p1b", parent=PART, fontSize=34))]
    s += [Spacer(1, 0.15 * inch)]
    s += [Paragraph("How a capsule is flown back through the atmosphere, modeled in software.",
                    ParagraphStyle("p1c", parent=COVER_S, textColor=colors.HexColor("#bcd0f5"),
                                   fontSize=13, alignment=TA_LEFT))]
    s += [NextPageTemplate("body"), PageBreak()]

    s += [Paragraph("1.1  What the simulator is", H1)]
    s += [P("The simulator is a detailed software model of an Orion style crew capsule returning from space. "
            "It tracks where the capsule is, how fast it is going, how hot it gets, how hard the deceleration "
            "pushes on the crew, and where it will land. It runs the whole flight in small time steps and "
            "writes out a complete record of everything that happened. It is the world that everything else "
            "in the project lives in.")]

    s += [Paragraph("1.2  The entry problem and the state", H2)]
    s += [P("At entry the capsule is roughly 120 km up and moving at many kilometers per second. The air is "
            "thin at first, then thickens quickly. The capsule rides on its blunt heat shield, using the small "
            "amount of lift the shield produces to steer. At every instant the capsule is described by six "
            "numbers: its distance from the center of the Earth, its latitude and longitude, its speed, its "
            "flight path angle (how steeply it is descending), and its heading (which way it points across the "
            "ground).")]
    s += fig(figs["geometry"], 4.6, "Figure 1. The entry geometry. Six numbers describe the vehicle at any moment.")

    s += [Paragraph("1.3  Equations of motion", H2)]
    s += [P("Those six numbers change over time according to the laws of motion for flight around a round "
            "planet. In plain terms:")]
    s += [bul("Altitude changes with the vertical part of the speed."),
          bul("Position over the ground changes with the horizontal part of the speed."),
          bul("Speed drops because of drag from the air and the pull of gravity."),
          bul("The descent angle bends from the upward part of lift fighting gravity."),
          bul("The heading turns from the sideways part of lift.")]
    s += [P("A single control, the <b>bank angle</b>, decides how the lift is split between pulling the nose up "
            "and turning it sideways. Rolling the capsule to the left or right is the only steering it has. The "
            "model advances these equations forward in steps of a quarter of a second.")]

    s += [Paragraph("1.4  Atmosphere", H2)]
    s += [P("To know the drag and heating, the model needs the air density and temperature at each altitude. "
            "These come from a standard empirical atmosphere model. For speed during long training runs, the "
            "values are precomputed once into a lookup table and read back instantly.")]

    s += [Paragraph("1.5  Aerodynamics", H2)]
    s += [P("The lift and drag come from a published Orion crew module aerodynamic database. The drag "
            "coefficient and the lift to drag ratio change with speed: blunt and draggy when very fast, "
            "different again through the sound barrier and below. The capsule always flies heat shield first. "
            "This is what gives it a real, vehicle accurate flight shape rather than a generic guess.")]

    s += [Paragraph("1.6  Heating", H2)]
    s += [P("Surviving the heat is the hardest part of entry. The model adds up two sources. The first is "
            "convective heating, the hot gas scraping the surface. The second is radiative heating, the glow "
            "of the shock layer, which only becomes important at the very high speeds of a return from the "
            "Moon. Together they give the heat rate at the nose, and from the balance of heat coming in versus "
            "heat radiated away the model computes the surface temperature.")]
    s += fig(figs["heating"], 6.0, "Figure 2. The heating model. Two sources add to a total, which sets the surface temperature.")
    s += [P("The heat shield itself is modeled as a grid of small patches, each tracking its own heat and "
            "temperature, so the total heat load over the whole flight can be measured. Two limits must be "
            "respected at all times.")]
    s += [kv_table([
        ["Quantity", "Value"],
        ["Vehicle mass", "8500 kg"],
        ["Reference area", "about 19.6 square meters"],
        ["Heat shield nose radius", "6 m"],
        ["Heat rate limit", "15 MW per square meter"],
        ["Total heat load limit", "2500 MJ per square meter"],
        ["Integration step", "0.25 s"],
        ["Guidance update rate", "once per second"],
    ], [2.6, 2.4])]

    s += [Paragraph("1.7  Controls: how the capsule actually rolls", H2)]
    s += [P("The bank angle is not changed instantly. A realistic chain of parts carries it out. Guidance "
            "asks for a bank angle. A bank actuator limits how fast it can change. A roll controller computes "
            "the torque needed. A set of twelve small thrusters, the reaction control system, fires in short "
            "pulses to produce that torque and roll the capsule. The fuel used by those thrusters is tracked "
            "throughout. This means the steering has real lag and real cost, just like the actual vehicle.")]

    s += [Paragraph("1.8  Parachutes", H2)]
    s += [P("Below about 7.6 km the parachute system takes over. Drogue chutes deploy first to stabilize and "
            "slow the capsule, then pilot chutes pull out the main chutes, which open in stages to avoid a "
            "violent jolt. The model includes the swinging of the capsule under the chutes, partial collapses, "
            "and a small random chance of a chute failing, so the landing phase is realistic rather than ideal.")]

    s += [Paragraph("1.9  Guaranteed bounds: the interval supervisor", H2)]
    s += [P("Alongside the normal best estimate trajectory, the simulator can carry a guaranteed band. Using "
            "interval arithmetic, it propagates not a single state but a small box of possible states, and it "
            "computes the worst case heat and worst case force that box could ever produce. This gives a hard, "
            "provable bound, not just a best guess. It is the key tool that later makes the learned controller "
            "safe.")]
    s += fig(figs["interval"], 5.6, "Figure 3. The interval band. The worst case edge reaches a limit before the best estimate does.")

    s += [Paragraph("1.10  The traditional controller", H2)]
    s += [P("For comparison, the simulator includes a classical guidance method called a predictor corrector. "
            "Once per second it quickly simulates the rest of the flight for several candidate bank angles, "
            "finds the bank that would land closest to the target, and flips the bank left or right to steer "
            "the heading. This is the established, hand designed approach, and it is the standard that the "
            "learned controller is measured against.")]

    s += [Paragraph("1.11  Putting it together", H2)]
    s += [P("In one cycle the controls choose and carry out a bank angle while the environment computes air, "
            "forces, and heat. The equations of motion then step the state forward, the guaranteed bounds "
            "update, and at low altitude the parachutes engage. The new state feeds back in and the cycle "
            "repeats until landing.")]
    s += fig(figs["architecture"], 7.0, "Figure 4. One cycle of the closed loop simulation. Controls on top, physics in the middle, integration and safety below.")

    s += [Paragraph("1.12  An example flight", H2)]
    s += [P("A return from the Moon arrives at about 11 km per second. To avoid crushing forces, the capsule "
            "uses a skip: it dips into the air, slows, skips back up out of the dense atmosphere, coasts, then "
            "re enters and descends to its parachutes. The simulator reproduces this naturally when given lunar "
            "entry conditions and a far target.")]
    s += fig(figs["skip"], 6.0, "Figure 5. A simulated lunar return. The altitude dips, skips back up, then descends to landing.")

    # =========================================================================
    # PART 2
    # =========================================================================
    s += [NextPageTemplate("part2"), PageBreak()]
    s += [Spacer(1, 2.4 * inch)]
    s += [Paragraph("Part 2", ParagraphStyle("p2", parent=PART, textColor=ORANGE))]
    s += [Paragraph("The Learning System", ParagraphStyle("p2b", parent=PART, fontSize=34))]
    s += [Spacer(1, 0.15 * inch)]
    s += [Paragraph("Training a neural network to steer, and proving it stays safe.",
                    ParagraphStyle("p2c", parent=COVER_S, textColor=colors.HexColor("#bcd0f5"),
                                   fontSize=13, alignment=TA_LEFT))]
    s += [NextPageTemplate("body"), PageBreak()]

    s += [Paragraph("2.1  The idea", H1)]
    s += [P("The traditional controller is carefully hand designed. The central question of this project is "
            "whether a controller that <b>learns</b> from experience can do better, especially when conditions "
            "vary or things go wrong, and whether such a learned controller can be made provably safe enough "
            "for a crewed vehicle. The learned controller is called a surrogate, because it stands in for the "
            "traditional guidance.")]

    s += [Paragraph("2.2  Steering as a learning problem", H2)]
    s += [P("Reinforcement learning frames steering as a loop. At each step a policy, which is a neural "
            "network, looks at the current situation and chooses a bank angle. The environment, the simulator, "
            "flies one second forward and returns a new situation and a score. Over millions of seconds of "
            "simulated flight, the network gradually learns which bank angles lead to high scores.")]
    s += fig(figs["rlloop"], 6.0, "Figure 6. The learning loop. The policy acts, the environment responds with a new situation and a score.")

    s += [Paragraph("2.3  What the policy sees, does, and is scored on", H2)]
    s += [kv_table([
        ["Element", "Description"],
        ["Observation", "About 13 numbers: altitude, speed, descent angle, heading, range to target, "
                        "tracking errors, energy, current bank, and how much heat and force margin remain."],
        ["Action", "One number: the bank angle, from one side to the other. The same control the "
                   "traditional method uses, so the comparison is fair."],
        ["Reward", "A score each step. Getting closer to the target earns points. Approaching the heat "
                   "or force limits loses points. Landing on target and staying survivable earns a large "
                   "bonus at the end. Losing control is heavily penalized."],
    ], [1.4, 4.5])]
    s += [Spacer(1, 6)]
    s += [P("The reward is the heart of the design. A dense progress signal, a small reward every step for "
            "closing the distance to the target, is what lets the network learn to steer at all. Without it, "
            "the only feedback would come at the very end of a long flight, far too sparse to learn from.")]

    s += [Paragraph("2.4  How the policy improves: PPO", H2)]
    s += [P("Training uses an algorithm called Proximal Policy Optimization. It keeps two networks. The actor "
            "chooses bank angles. The critic estimates how good a situation is. The trainer lets the actor try "
            "many flights, compares how each action turned out against the critic's expectation, and nudges the "
            "actor toward the actions that did better than expected. The nudges are deliberately small, so the "
            "policy improves steadily without suddenly forgetting what already worked.")]

    s += [Paragraph("2.5  The environment and the fair baseline", H2)]
    s += [P("The simulator is wrapped in a standard learning interface so any modern algorithm can drive it. "
            "The same wrapper can run in two modes. In policy mode the learned network flies. In baseline mode "
            "the traditional predictor corrector flies. Because both run through the exact same physics, the "
            "two can be compared with no hidden advantages.")]

    s += [Paragraph("2.6  Random conditions and a frozen test set", H2)]
    s += [P("If the policy only ever saw one entry, it would memorize that one case. Instead, each training "
            "episode randomizes the entry conditions, the speed, the angle, the heading, and the starting "
            "point, so the policy must learn to handle a whole range of situations. Separately, a fixed set of "
            "one thousand cases is set aside and never used for training. Both the learned and the traditional "
            "controller are scored on these identical cases, which makes the final comparison trustworthy.")]

    s += [Paragraph("2.7  The training and evaluation pipeline", H2)]
    s += [P("Heavy training runs in the cloud on many copies of the environment at once, because the slow part "
            "is the physics, not the network. Progress is logged continuously and checkpoints are saved so a "
            "long run can survive interruptions. Once trained, the policy and the traditional controller are "
            "both run over the frozen test set, and their results are compared as full distributions rather "
            "than single numbers.")]
    s += fig(figs["pipeline"], 7.0, "Figure 7. Training on random conditions, then a fair head to head evaluation on identical frozen cases.")

    s += [Paragraph("2.8  Bringing guaranteed safety into learning", H2)]
    s += [P("A learned controller usually has no safety guarantee, which is a serious problem for a crewed "
            "vehicle. The guaranteed bounds from Part 1 are brought into the learning system in three "
            "independent ways:")]
    s += [bul("<b>Awareness.</b> The worst case heat and force margins are added to what the policy sees, so it "
              "knows how much uncertainty it is carrying."),
          bul("<b>Robust reward.</b> The heat and force penalties are based on the guaranteed worst case rather "
              "than the best estimate, so the policy learns to keep a safe buffer."),
          bul("<b>Safety shield.</b> Before any bank angle is flown, a quick guaranteed check asks whether the "
              "worst case would break a limit in the next step. If so, the bank is replaced with the nearest "
              "one that stays safe.")]
    s += fig(figs["shield"], 6.0, "Figure 8. The safety shield. The learned choice is checked by guaranteed math before it flies.")

    s += [Paragraph("2.9  Safety and performance together", H2)]
    s += [P("This combination is the main contribution. The learned policy supplies performance and the "
            "ability to adapt. The guaranteed math supplies certified safety. The learned controller proposes, "
            "and the provable bounds dispose. Evaluation therefore looks beyond miss distance, force, and heat "
            "to two more measures: how often a safety limit is broken, which should fall to zero, and how often "
            "the shield has to step in, which should shrink as the policy learns the safe envelope on its own.")]
    s += [P("The honest limits are stated plainly. The shield checks a short horizon, so it is rigorous for "
            "the instant by instant limits rather than the whole future, and the guarantees hold under a stated "
            "model of the uncertainty. These are the right caveats for an approach that aims to be both useful "
            "and trustworthy.")]

    s += [Paragraph("2.10  Where the project stands", H2)]
    s += [kv_table([
        ["Stage", "Status"],
        ["Simulator with heating, controls, parachutes, and guaranteed bounds", "built"],
        ["Strong traditional controller as the baseline", "built"],
        ["Learning environment and fair baseline mode", "built"],
        ["Random conditions and frozen test set", "built"],
        ["Interval awareness, robust reward, and safety shield", "built"],
        ["Cloud training and head to head evaluation", "in progress"],
    ], [4.1, 1.4])]
    s += [Spacer(1, 8)]
    s += [P("In short, the simulator provides a faithful world, the traditional controller provides a strong "
            "bar to clear, and the learning system aims to clear it while staying provably safe. That pairing, "
            "a learned controller checked by guaranteed bounds, is what makes the approach suitable for the "
            "demanding setting of returning a crew from space.")]

    doc.build(s)
