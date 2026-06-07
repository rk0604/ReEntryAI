"""
Generate the ReEntryAI technical overview PDF (diagrams + text).
Run: python docs/make_report.py
Output: docs/ReEntryAI_Technical_Overview.pdf
"""
from __future__ import annotations
import os, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Wedge, Polygon
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "_figs")
os.makedirs(FIG, exist_ok=True)

# ---- palette ----------------------------------------------------------------
NAVY = "#16294d"; BLUE = "#2d6cdf"; ORANGE = "#e8833a"; GREEN = "#2e9e5b"
GRAY = "#e9edf3"; LGRAY = "#f4f6fa"; RED = "#d2433a"; PURPLE = "#7b5cd6"
TEXT = "#16202b"; MUTED = "#5b6b7c"

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})


def _box(ax, x, y, w, h, text, fc=GRAY, ec=NAVY, tc=TEXT, fs=10.5, bold=False, lw=1.4):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle="round,pad=0.015,rounding_size=0.05",
                 fc=fc, ec=ec, lw=lw, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs, color=tc, zorder=3,
            fontweight="bold" if bold else "normal")
    return (x + w / 2, y + h / 2, x, y, w, h)


def _arrow(ax, p1, p2, color=NAVY, lw=1.7, rad=0.0, style="-|>"):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle=style, mutation_scale=15,
                 color=color, lw=lw, connectionstyle=f"arc3,rad={rad}", zorder=1))


def _fig(w=10, h=6):
    f, ax = plt.subplots(figsize=(w, h))
    ax.set_xlim(0, 10); ax.set_ylim(0, h / w * 10)
    ax.axis("off")
    return f, ax


def save(f, name):
    p = os.path.join(FIG, name)
    f.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(f)
    return p


# ============================================================================
# D1 : entry geometry / coordinate frame
# ============================================================================
def d_geometry():
    f, ax = plt.subplots(figsize=(8, 6)); ax.set_aspect("equal"); ax.axis("off")
    R = 3.0
    earth = Circle((0, -R + 0.2), R, fc="#0a1e3a", ec=BLUE, lw=1.5, zorder=1)
    ax.add_patch(earth)
    ax.add_patch(Wedge((0, -R + 0.2), R, 60, 120, width=0.05, fc=GREEN, zorder=2))
    # capsule position
    ang = math.radians(78)
    cx, cy = R * math.cos(ang), -R + 0.2 + R * math.sin(ang)
    cx += 0.0
    px, py = cx * 1.55, (cy + R - 0.2) * 1.55 - R + 0.2
    ax.plot([0, px], [-R + 0.2, py], color=MUTED, lw=1.0, ls=":", zorder=2)
    ax.add_patch(Circle((px, py), 0.16, fc=ORANGE, ec=NAVY, lw=1.2, zorder=4))
    # velocity vector
    vx, vy = 1.6, -0.9
    _arrow(ax, (px, py), (px + vx, py + vy), color=RED, lw=2.2)
    ax.text(px + vx + 0.05, py + vy - 0.1, "velocity V", color=RED, fontsize=11)
    # local horizon
    hx, hy = 1.7, 0.0
    nx, ny = px / math.hypot(px, py + R - 0.2), (py + R - 0.2) / math.hypot(px, py + R - 0.2)
    tx, ty = -ny, nx
    ax.plot([px - tx * 1.4, px + tx * 1.4], [py - ty * 1.4, py + ty * 1.4],
            color=MUTED, lw=1.2, zorder=2)
    ax.text(px - tx * 1.5, py - ty * 1.5 - 0.12, "local horizon",
            color=MUTED, fontsize=10, ha="right")
    ax.text(px, py + 0.46, "capsule", color=NAVY, fontsize=11, ha="center", va="bottom")
    # gamma label
    ax.text(px + 0.65, py - 0.30, "flight path\nangle gamma", color=NAVY, fontsize=9.5)
    # r and altitude
    ax.annotate("", xy=(px, py), xytext=(0, -R + 0.2),
                arrowprops=dict(arrowstyle="<->", color=NAVY, lw=1.0))
    ax.text(px / 2 - 0.55, (py - R + 0.2) / 2 + 0.1, "radius r", color=NAVY, fontsize=10, rotation=58)
    ax.text(0, -R + 0.2 - 0.05, "Earth center", color="white", fontsize=9, ha="center", va="top")
    ax.text(-2.8, 1.7, "State of the vehicle (six numbers):\n"
            "radius r, latitude, longitude,\nspeed V, flight path angle,\nheading angle",
            fontsize=10.5, color=TEXT, va="top",
            bbox=dict(boxstyle="round,pad=0.5", fc=LGRAY, ec=NAVY, lw=1.0))
    ax.set_xlim(-3.2, 3.2); ax.set_ylim(-3.2, 2.2)
    return save(f, "d_geometry.png")


# ============================================================================
# D2 : simulator architecture / data flow
# ============================================================================
def d_architecture():
    f, ax = _fig(11, 7.6)
    H = ax.get_ylim()[1]
    # guidance / control chain (top)
    g = _box(ax, 0.2, H - 1.3, 1.9, 0.9, "Guidance\n(chooses bank angle)", fc="#dce8ff", ec=BLUE, bold=True)
    ba = _box(ax, 2.5, H - 1.3, 1.7, 0.9, "Bank actuator\n(rate limits)", fc="#dce8ff", ec=BLUE)
    pd = _box(ax, 4.6, H - 1.3, 1.7, 0.9, "Roll controller\n(PD law)", fc="#dce8ff", ec=BLUE)
    rcs = _box(ax, 6.7, H - 1.3, 1.8, 0.9, "Reaction control\n(12 thrusters)", fc="#dce8ff", ec=BLUE)
    tq = _box(ax, 8.9, H - 1.3, 1.0, 0.9, "torque", fc="#dce8ff", ec=BLUE)
    for a, b in [(g, ba), (ba, pd), (pd, rcs), (rcs, tq)]:
        _arrow(ax, (a[0] + a[4] / 2, a[1] + 0.45), (b[0] - b[4] / 2, b[1] + 0.45), color=BLUE)
    # environment chain (middle)
    st = _box(ax, 0.2, H - 3.1, 1.9, 0.9, "Vehicle state\n(six numbers)", fc="#fff0e0", ec=ORANGE, bold=True)
    atm = _box(ax, 2.5, H - 3.1, 1.7, 0.9, "Atmosphere\n(density, temp)", fc=GRAY)
    aero = _box(ax, 4.6, H - 3.1, 1.7, 0.9, "Aerodynamics\n(lift, drag)", fc=GRAY)
    heat = _box(ax, 6.7, H - 3.1, 1.8, 0.9, "Heating\n(heat shield)", fc=GRAY)
    for a, b in [(st, atm), (atm, aero), (aero, heat)]:
        _arrow(ax, (a[0] + a[4] / 2, a[1] + 0.45), (b[0] - b[4] / 2, b[1] + 0.45), color=ORANGE)
    # integrator
    eom = _box(ax, 3.4, H - 5.0, 3.2, 1.0,
               "Equations of motion integrator\n(steps the state forward 0.25 s)",
               fc="#e3f6ea", ec=GREEN, bold=True)
    _arrow(ax, (tq[0], tq[1] - 0.45), (eom[0] + 1.3, eom[1] + 0.5), color=BLUE, rad=-0.25)
    _arrow(ax, (heat[0], heat[1] - 0.45), (eom[0] + 0.9, eom[1] + 0.5), color=ORANGE, rad=0.15)
    _arrow(ax, (aero[0], aero[1] - 0.45), (eom[0], eom[1] + 0.5), color=ORANGE, rad=0.1)
    # feedback to state
    _arrow(ax, (eom[2], eom[1] + 0.5), (st[0], st[1] - 0.45), color=GREEN, rad=0.35)
    ax.text(1.0, H - 4.2, "new state\nfeeds back", color=GREEN, fontsize=9, ha="center")
    # CPAS + supervisor side
    cpas = _box(ax, 7.1, H - 5.0, 2.6, 1.0,
                "Parachutes (CPAS)\ndrogue, pilot, main\nbelow 7.6 km", fc="#f0e8ff", ec=PURPLE)
    _arrow(ax, (eom[0] + 1.6, eom[1]), (cpas[0] - 1.3, eom[1]), color=PURPLE)
    sup = _box(ax, 0.3, H - 5.0, 2.6, 1.0,
               "Interval supervisor\nguaranteed worst case\nbounds on heat and g",
               fc="#fdeaea", ec=RED)
    _arrow(ax, (eom[2], eom[1]), (sup[0] + 1.3, eom[1]), color=RED, style="<|-|>")
    ax.text(5.0, H - 0.25, "Closed loop simulation, one cycle", fontsize=12.5,
            color=NAVY, ha="center", fontweight="bold")
    return save(f, "d_architecture.png")


# ============================================================================
# D3 : heating model
# ============================================================================
def d_heating():
    f, ax = _fig(10, 4.6)
    H = ax.get_ylim()[1]
    conv = _box(ax, 0.3, H - 1.6, 2.6, 1.0,
                "Convective heating\n(hot gas touching\nthe surface)", fc="#fff0e0", ec=ORANGE)
    rad = _box(ax, 0.3, H - 3.1, 2.6, 1.0,
               "Radiative heating\n(glowing shock layer,\nmatters above 9 km/s)", fc="#fdeaea", ec=RED)
    tot = _box(ax, 3.8, H - 2.35, 2.4, 1.0, "Total heat rate\non the nose", fc=GRAY, ec=NAVY, bold=True)
    wall = _box(ax, 7.0, H - 2.35, 2.6, 1.0,
                "Surface temperature\n(balance of heat in\nand heat radiated out)",
                fc="#e3f6ea", ec=GREEN)
    _arrow(ax, (conv[0] + 1.3, conv[1]), (tot[0] - 1.2, tot[1] + 0.2), color=ORANGE, rad=-0.15)
    _arrow(ax, (rad[0] + 1.3, rad[1]), (tot[0] - 1.2, tot[1] - 0.2), color=RED, rad=0.15)
    _arrow(ax, (tot[0] + 1.2, tot[1]), (wall[0] - 1.3, wall[1]), color=NAVY)
    ax.text(5.0, H - 0.3, "How the heat shield load is computed", fontsize=12.5,
            color=NAVY, ha="center", fontweight="bold")
    ax.text(5.0, 0.15, "Limits the vehicle must respect: heat rate under 15 MW/m^2,"
            " total heat load under 2500 MJ/m^2",
            fontsize=9.5, color=MUTED, ha="center")
    return save(f, "d_heating.png")


# ============================================================================
# D4 : interval tube concept
# ============================================================================
def d_interval():
    f, ax = plt.subplots(figsize=(9, 4.4)); ax.axis("off")
    t = np.linspace(0, 10, 200)
    nom = 6 + 3.2 * np.exp(-((t - 5) ** 2) / 4)
    width = 0.25 + 0.06 * t + 0.5 * np.exp(-((t - 5) ** 2) / 6)
    ax.fill_between(t, nom - width, nom + width, color=BLUE, alpha=0.18, label="guaranteed worst case band")
    ax.plot(t, nom, color=BLUE, lw=2.0, label="best estimate (nominal)")
    ax.plot(t, nom + width, color=BLUE, lw=0.8, ls=":")
    ax.plot(t, nom - width, color=BLUE, lw=0.8, ls=":")
    ax.axhline(9.4, color=RED, lw=1.6)
    ax.text(0.2, 9.55, "safety limit", color=RED, fontsize=10)
    i = np.argmax(nom + width > 9.4)
    if nom[i] + width[i] > 9.4:
        ax.scatter([t[i]], [9.4], color=RED, zorder=5, s=40)
        ax.annotate("worst case edge reaches the\nlimit before the best estimate does",
                    xy=(t[i], 9.4), xytext=(6.6, 4.2), color=RED, fontsize=9.5,
                    ha="left", arrowprops=dict(arrowstyle="->", color=RED, lw=1.2))
    ax.set_xlim(0, 10); ax.set_ylim(2, 10.4)
    ax.set_xlabel("time"); ax.set_ylabel("heat rate or g load")
    ax.set_title("Interval arithmetic: a guaranteed band around the trajectory",
                 color=NAVY, fontsize=12.5, fontweight="bold")
    ax.legend(loc="lower right", fontsize=9, frameon=True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return save(f, "d_interval.png")


# ============================================================================
# D5 : Artemis skip altitude profile (real data if available)
# ============================================================================
def d_skip():
    f, ax = plt.subplots(figsize=(9, 4.3))
    try:
        import pandas as pd
        df = pd.read_csv(os.path.join(HERE, "..", "revision_v1", "trajectory.csv"))
        t = df["t_s"].values; alt = df["alt_m"].values / 1000.0
        ax.plot(t, alt, color=GREEN, lw=2.2)
        ax.set_xlim(0, t.max())
    except Exception:
        t = np.linspace(0, 720, 400)
        alt = 122 * np.exp(-t / 220) + 30 * np.exp(-((t - 150) ** 2) / 1200)
        ax.plot(t, alt, color=GREEN, lw=2.2)
    ax.set_ylim(0, 130)
    ax.set_xlabel("time after entry (s)"); ax.set_ylabel("altitude (km)")
    ax.set_title("Example lunar return: a skip entry at 11 km/s",
                 color=NAVY, fontsize=12.5, fontweight="bold")
    ax.annotate("dips into the air,\nslows down", xy=(95, 57), xytext=(150, 90),
                fontsize=9.5, color=NAVY,
                arrowprops=dict(arrowstyle="->", color=NAVY))
    ax.annotate("skips back up,\nthen re enters", xy=(180, 86), xytext=(250, 105),
                fontsize=9.5, color=NAVY,
                arrowprops=dict(arrowstyle="->", color=NAVY))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return save(f, "d_skip.png")


# ============================================================================
# D6 : RL closed loop (MDP)
# ============================================================================
def d_rlloop():
    f, ax = _fig(9.5, 4.8)
    H = ax.get_ylim()[1]
    env = _box(ax, 5.6, H - 3.1, 3.4, 1.6,
               "Environment\n(the full simulator)\nflies the vehicle one\nsecond forward",
               fc="#fff0e0", ec=ORANGE, bold=True, fs=11)
    pol = _box(ax, 0.5, H - 3.1, 3.4, 1.6,
               "Policy\n(a neural network)\nlooks at the situation,\npicks a bank angle",
               fc="#dce8ff", ec=BLUE, bold=True, fs=11)
    _arrow(ax, (pol[0] + 1.7, pol[1] + 0.45), (env[0] - 1.7, env[1] + 0.45), color=BLUE)
    ax.text(4.75, env[1] + 0.95, "action:\nbank angle", color=BLUE, fontsize=9.5, ha="center")
    _arrow(ax, (env[0] - 1.7, env[1] - 0.45), (pol[0] + 1.7, pol[1] - 0.45), color=ORANGE)
    ax.text(4.75, env[1] - 1.0, "observation (situation)\nand reward (score)",
            color=ORANGE, fontsize=9.5, ha="center")
    ax.text(4.75, H - 0.3, "The control loop as a learning problem",
            color=NAVY, fontsize=12.5, ha="center", fontweight="bold")
    ax.text(4.75, 0.15, "Repeat for millions of seconds of flight; the network slowly"
            " learns bank angles that score well.",
            color=MUTED, fontsize=9.5, ha="center")
    return save(f, "d_rlloop.png")


# ============================================================================
# D7 : training + evaluation pipeline
# ============================================================================
def d_pipeline():
    f, ax = _fig(11, 5.6)
    H = ax.get_ylim()[1]
    disp = _box(ax, 0.2, H - 1.6, 2.3, 1.0, "Random entry\nconditions\n(every episode)", fc=GRAY)
    env = _box(ax, 2.9, H - 1.6, 2.0, 1.0, "Environment\n(simulator)", fc="#fff0e0", ec=ORANGE)
    ppo = _box(ax, 5.3, H - 1.6, 2.2, 1.0, "PPO trainer\n(improves the\npolicy)", fc="#dce8ff", ec=BLUE, bold=True)
    pol = _box(ax, 7.9, H - 1.6, 2.0, 1.0, "Trained\npolicy", fc="#dce8ff", ec=BLUE, bold=True)
    for a, b in [(disp, env), (env, ppo), (ppo, pol)]:
        _arrow(ax, (a[0] + a[4] / 2, a[1] + 0.5), (b[0] - b[4] / 2, b[1] + 0.5), color=NAVY)
    _arrow(ax, (ppo[0], ppo[1] - 0.5), (env[0], env[1] - 0.5), color=BLUE, rad=0.4)
    ax.text(4.1, H - 2.4, "policy keeps trying", color=BLUE, fontsize=8.5, ha="center")
    # evaluation
    test = _box(ax, 0.2, H - 3.9, 2.3, 1.0, "Frozen test set\n1000 fixed cases\n(never trained on)", fc="#f0e8ff", ec=PURPLE)
    both = _box(ax, 2.9, H - 3.9, 2.6, 1.0, "Run BOTH controllers\nclassical and learned\non the same cases", fc=LGRAY, ec=NAVY, bold=True)
    metrics = _box(ax, 6.0, H - 3.9, 3.7, 1.0,
                   "Compare distributions\nmiss, peak g, heat, fuel,\nsafety violations", fc="#e3f6ea", ec=GREEN)
    _arrow(ax, (test[0] + 1.15, test[1]), (both[0] - 1.3, both[1]), color=PURPLE)
    _arrow(ax, (both[0] + 1.3, both[1]), (metrics[0] - 1.85, metrics[1]), color=NAVY)
    _arrow(ax, (pol[0], pol[1] - 0.5), (both[0] + 0.6, both[1] + 0.5), color=BLUE, rad=-0.3)
    ax.text(5.5, H - 0.3, "Training (top) and fair evaluation (bottom)",
            color=NAVY, fontsize=12.5, ha="center", fontweight="bold")
    ax.text(5.5, 0.15, "Heavy training runs in the cloud with full logging;"
            " both controllers face identical conditions.",
            color=MUTED, fontsize=9.5, ha="center")
    return save(f, "d_pipeline.png")


# ============================================================================
# D8 : interval safety shield decision
# ============================================================================
def d_shield():
    f, ax = _fig(9.5, 5.2)
    H = ax.get_ylim()[1]
    pol = _box(ax, 0.3, H - 1.4, 2.4, 0.9, "Policy proposes\na bank angle", fc="#dce8ff", ec=BLUE, bold=True)
    chk = _box(ax, 3.2, H - 1.55, 3.0, 1.2,
               "Interval check:\nwould the worst case\nbreak a heat or g limit\nnext step?",
               fc="#fdeaea", ec=RED, bold=True, fs=10)
    _arrow(ax, (pol[0] + 1.2, pol[1] + 0.0), (chk[0] - 1.5, chk[1] + 0.05), color=BLUE)
    safe = _box(ax, 7.0, H - 1.45, 2.3, 0.9, "No: use the\npolicy bank", fc="#e3f6ea", ec=GREEN, bold=True)
    _arrow(ax, (chk[0] + 1.5, chk[1] + 0.05), (safe[0] - 1.2, safe[1] + 0.0), color=GREEN)
    ax.text(6.72, chk[1] + 0.42, "safe", color=GREEN, fontsize=9.5)
    fix = _box(ax, 3.0, H - 3.5, 3.4, 1.0,
               "Yes: replace with the\nnearest bank that stays\nwithin the limits", fc="#fff0e0", ec=ORANGE, bold=True)
    _arrow(ax, (chk[0], chk[1] - 0.6), (fix[0], fix[1] + 0.5), color=RED)
    ax.text(chk[0] + 0.25, chk[1] - 0.95, "unsafe", color=RED, fontsize=9.5)
    fly = _box(ax, 3.4, H - 5.0, 2.6, 0.9, "Fly the safe bank", fc=GRAY, ec=NAVY, bold=True)
    _arrow(ax, (fix[0], fix[1] - 0.5), (fly[0], fly[1] + 0.45), color=ORANGE)
    _arrow(ax, (safe[0], safe[1] - 0.45), (fly[0] + 1.3, fly[1] + 0.45), color=GREEN, rad=0.3)
    ax.text(4.75, H - 0.25, "The safety shield: learned choice, checked by guaranteed math",
            color=NAVY, fontsize=12, ha="center", fontweight="bold")
    return save(f, "d_shield.png")


def build_all():
    return {
        "geometry": d_geometry(), "architecture": d_architecture(),
        "heating": d_heating(), "interval": d_interval(), "skip": d_skip(),
        "rlloop": d_rlloop(), "pipeline": d_pipeline(), "shield": d_shield(),
    }


if __name__ == "__main__":
    figs = build_all()
    from report_text import build_pdf
    build_pdf(figs, os.path.join(HERE, "ReEntryAI_Technical_Overview.pdf"))
    print("done")
