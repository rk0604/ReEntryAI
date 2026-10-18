"""
jepa/report_v2_figures.py

Figures for the SpaceJEPA v2 beginner's report (jepa/make_report_v2.py).

Every figure is one of three kinds, and its caption in the report says which:
    real data      flights from the simulator (wide-envelope episodes, v1's
                   episodes, closed-loop runs of the shipped v1 checkpoint and
                   the classical controller)
    from the code  numbers computed by the v2 code itself (attention masks,
                   schedules, parameter counts, FLOPs)
    illustration   hand-drawn schematics and synthetic examples

    python jepa/report_v2_figures.py --wide <pickle or shard dir> --dense <dense_runs.pkl> --out <dir>
    (jepa/make_report_v2.py calls make_all and builds the PDF.)
"""

from __future__ import annotations

import math
import os
import pickle
import sys
from pathlib import Path

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle, Wedge

from jepa import data as D
from jepa import models_v2 as V
from jepa.make_figures import AQUA, BLUE, GRAY, INK, MUTED, ORANGE, VIOLET, YELLOW, arrow, box, canvas, panel

# Entity colors for charts (same slots as the v2 notebook).
C_V1, C_V2, C_V1W, C_CL = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
C_MUTED, C_SEC = "#a3a29c", "#52514e"
DIVERGE = LinearSegmentedColormap.from_list("thin_dense", ["#2a78d6", "#e8e7e3", "#eb6834"])
W = 7.0                                # figure width, inches (= the report's text width)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans", "font.size": 8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e4e3df",
    "grid.linewidth": 0.6, "axes.edgecolor": "#8a8984", "axes.labelcolor": C_SEC, "axes.titlecolor": "#0b0b0b",
    "axes.titlesize": 8.5, "xtick.color": C_SEC, "ytick.color": C_SEC, "lines.linewidth": 1.5,
    "legend.frameon": False, "legend.fontsize": 7, "savefig.dpi": 220,
})


def save(fig, out, name):
    """Write a PNG and return its path."""
    path = os.path.join(out, name)
    fig.savefig(path, bbox_inches="tight", facecolor="white", pad_inches=0.04)
    plt.close(fig)
    return path


def diagram(h, xmax=100.0):
    """A blank schematic canvas W x h inches with x in [0, xmax], y in [0, xmax * h / W]."""
    return canvas(W, h, xmax, xmax * h / W)


# ============================================================================
# Part 1: the problem, the data, the pipeline
# ============================================================================
def fig_basics(eps, out):
    """Real altitude profiles (left) and the bank-angle idea (right)."""
    fig = plt.figure(figsize=(W, 2.6))
    ax = fig.add_axes([0.07, 0.17, 0.48, 0.72])
    sel = [e for e in eps if e["outcome"] == "drogue"][:14]
    for i, e in enumerate(sel):
        t = np.arange(len(e["y"]))
        ax.plot(t, e["y"][:, 0] / 1000, color=C_V2 if i == 0 else C_MUTED, lw=1.8 if i == 0 else 0.8,
                zorder=3 if i == 0 else 2)
    e = sel[0]
    ax.axhline(7.6, color=C_SEC, lw=0.8, ls=":")
    ax.text(330, 11, "parachute opens at 7.6 km", fontsize=7, color=C_SEC)
    ax.annotate(f"entry: {e['y'][0, 0] / 1000:.0f} km, {e['y'][0, 3] / 1000:.1f} km/s",
                xy=(0, e["y"][0, 0] / 1000), xytext=(60, e["y"][0, 0] / 1000 + 6), fontsize=7, color=INK,
                arrowprops=dict(arrowstyle="-", color=C_SEC, lw=0.7))
    ax.set(xlabel="time since entry [s]", ylabel="altitude [km]", ylim=(0, 150),
           title="14 real flights from the training data")
    # Bank-angle sketch.
    ax2 = fig.add_axes([0.60, 0.02, 0.40, 0.95])
    ax2.set_xlim(-1.4, 1.6)
    ax2.set_ylim(-1.15, 1.35)
    ax2.set_aspect("equal")
    ax2.axis("off")
    ax2.add_patch(Circle((0, 0), 0.42, fc=GRAY[1], ec=GRAY[0], lw=1.2))
    ax2.text(0, 0, "capsule\n(seen from\nbehind)", ha="center", va="center", fontsize=6.5, color=C_SEC)
    s = math.radians(35)
    L = 1.05
    ax2.add_patch(FancyArrowPatch((0, 0.42), (0, L), arrowstyle="-|>", mutation_scale=9, color=C_MUTED, lw=1, ls="--"))
    ax2.add_patch(FancyArrowPatch((0.42 * math.sin(s), 0.42 * math.cos(s)), (L * math.sin(s), L * math.cos(s)),
                                  arrowstyle="-|>", mutation_scale=11, color=C_V2, lw=2))
    ax2.add_patch(Wedge((0, 0), 0.62, 90 - 35, 90, fc="none", ec=INK, lw=0.8))
    ax2.text(0.18, 0.72, "σ", fontsize=10, color=INK)
    ax2.text(L * math.sin(s) + 0.05, L * math.cos(s) + 0.02, "lift", fontsize=8, color=C_V2, fontweight="bold")
    ax2.text(-0.08, L + 0.05, "up", fontsize=7, color=C_SEC, ha="right")
    ax2.text(-1.35, -0.62, "Rolling to bank angle σ tilts the lift:\n"
             "• lift up (σ = 0°): flatter, longer flight\n"
             "• lift sideways (σ = ±90°): turns, drops faster\n"
             "Bank angle is the only steering control.", fontsize=7, color=INK, va="top", linespacing=1.45)
    return save(fig, out, "p1_basics.png")


def fig_pipeline(out):
    """The v2 pipeline in four stages (illustration)."""
    fig, ax = diagram(2.75)
    fs = 5.9
    cols = [(1, "(1) Data", GRAY), (26, "(2) World model", BLUE), (51, "(3) Prober", VIOLET),
            (76, "(4) Fly and test", ORANGE)]
    for x0, title, col in cols:
        panel(ax, x0, 1, x0 + 23, 38.5, title, col)
    box(ax, 12.5, 29.5, 21.5, 6.5, "6,000 flights,\nrandom bank profiles", GRAY, fs)
    box(ax, 12.5, 19.5, 21.5, 6.5, "full simulator:\nactuator, RCS, 3D physics", GRAY, fs)
    box(ax, 12.5, 9.0, 21.5, 7.5, "random truth per flight:\ndensity, L/D ×0.7–1.3,\nany entry condition", GRAY, fs)
    arrow(ax, (12.5, 26.2), (12.5, 22.8))
    arrow(ax, (12.5, 12.8), (12.5, 16.2))
    box(ax, 37.5, 29.5, 21.5, 6.5, "encoder: history →\nlatent code s", BLUE, fs)
    box(ax, 37.5, 19.5, 21.5, 6.5, "predictor: s + bank plan\n→ future codes", ORANGE, fs)
    box(ax, 37.5, 9.0, 21.5, 7.5, "loss: predicted codes ≈\nreal future codes\n+ SIGReg (no collapse)", YELLOW, fs)
    arrow(ax, (37.5, 26.2), (37.5, 22.8))
    arrow(ax, (37.5, 16.2), (37.5, 12.8))
    box(ax, 62.5, 29.5, 21.5, 6.5, "frozen world model\n(5 ensemble members)", BLUE, fs, dashed=True)
    box(ax, 62.5, 19.5, 21.5, 6.5, "prober: code → drag,\nlift and bank corrections", VIOLET, fs)
    box(ax, 62.5, 9.0, 21.5, 7.5, "exact 3D equations of\nmotion → altitude,\nposition, speed …", VIOLET, fs)
    arrow(ax, (62.5, 26.2), (62.5, 22.8))
    arrow(ax, (62.5, 16.2), (62.5, 12.8))
    box(ax, 87.5, 29.5, 21.5, 6.5, "MPPI: try 256 bank plans\nin the model every 5 s", ORANGE, fs)
    box(ax, 87.5, 19.5, 21.5, 6.5, "fly the best one in the\nreal simulator, repeat", GRAY, fs)
    box(ax, 87.5, 9.0, 21.5, 7.5, "compare with v1 and the\nclassical controller\nunder 4 atmospheres", YELLOW, fs)
    arrow(ax, (87.5, 26.2), (87.5, 22.8))
    arrow(ax, (87.5, 16.2), (87.5, 12.8))
    for xa, xb in ((24.2, 26.8), (49.2, 51.8), (74.2, 76.8)):
        arrow(ax, (xa, 19.5), (xb, 19.5), lw=1.6)
    return save(fig, out, "p1_pipeline.png")


def fig_jepa_idea(out):
    """Predict raw states vs. predict codes (illustration)."""
    fig, ax = diagram(2.9)
    fs = 6.8
    panel(ax, 1, 21.5, 99, 40.0, "Option A: predict the raw future directly (a \"predictive\" model)", GRAY)
    box(ax, 12, 29, 17, 6, "state now\n(altitude, speed, …)", GRAY, fs)
    box(ax, 35, 29, 16, 6, "model", ORANGE, fs, bold=True)
    box(ax, 58, 29, 17, 6, "predicted next\nraw state", GRAY, fs)
    box(ax, 86, 29, 20, 6, "feed it back in,\nerrors pile up", YELLOW, fs)
    arrow(ax, (20.5, 29), (27, 29))
    arrow(ax, (43, 29), (49.5, 29))
    arrow(ax, (66.5, 29), (76, 29))
    arrow(ax, (86, 25.8), (35, 25.8), ls="--", rad=-0.0, col=C_MUTED)
    panel(ax, 1, 1, 99, 20.0, "JEPA: translate into a compact code, and predict the future code", BLUE)
    box(ax, 10, 11.5, 14, 6, "history now\n(2 min)", GRAY, fs)
    box(ax, 26, 11.5, 11, 6, "encoder", BLUE, fs, bold=True)
    box(ax, 39, 11.5, 7, 5, "code\ns", BLUE, fs)
    box(ax, 53.5, 11.5, 14, 6, "predictor\n+ bank plan", ORANGE, fs, bold=True)
    box(ax, 68.5, 11.5, 9, 5, "future\ncode", ORANGE, fs)
    box(ax, 82, 11.5, 8, 5, "real\ncode", BLUE, fs)
    box(ax, 93, 11.5, 9, 6, "encoder\n(future)", BLUE, fs)
    arrow(ax, (17, 11.5), (20.5, 11.5))
    arrow(ax, (31.5, 11.5), (35.5, 11.5))
    arrow(ax, (42.5, 11.5), (46.5, 11.5))
    arrow(ax, (60.5, 11.5), (64, 11.5))
    arrow(ax, (88.5, 11.5), (86, 11.5))
    ax.annotate("", xy=(73.2, 11.5), xytext=(78, 11.5), arrowprops=dict(arrowstyle="<->", color=YELLOW[0], lw=1.4))
    ax.text(75.6, 15.3, "match", ha="center", fontsize=6.6, color=YELLOW[0], fontweight="bold")
    ax.text(50, 3.0, "The code can drop details that don't matter for predicting the flight.\nA later \"prober\" "
            "turns codes back into altitude, position and speed.", ha="center", fontsize=6.3, color=C_SEC)
    return save(fig, out, "p1_jepa_idea.png")


def fig_data_examples(eps, out):
    """Commanded vs. realized bank angle on real flights: the actuator lag."""
    gp = next(e for e in eps if e["kind"] == "gp" and e["outcome"] == "drogue")
    rv = next(e for e in eps if e["kind"] == "reversal" and e["outcome"] == "drogue")
    fig, axes = plt.subplots(2, 1, figsize=(W, 3.0), sharex=True)
    for ax, e, name in ((axes[0], gp, "smooth random (\"GP\") profile"), (axes[1], rv, "reversal profile")):
        t = np.arange(len(e["cmd"]))
        ax.step(t, np.degrees(e["cmd"]), where="post", color=C_MUTED, lw=1.1, label="commanded bank")
        ax.plot(t, np.degrees(e["y"][1:, 6]), color=C_V1, lw=1.5, label="realized bank (actuator lag)")
        d = e["domain"]
        ax.set_title(f"{name}: density ×{d['density_scale']:.2f}, L/D ×{d['lift_scale']:.2f}, "
                     f"flight {e['t_final_s']:.0f} s", loc="left")
        ax.set_ylabel("bank [deg]")
        ax.set_ylim(-115, 115)
        ax.set_yticks([-90, 0, 90])
    axes[0].legend(loc="lower right", ncol=2, bbox_to_anchor=(1.0, 1.12))
    axes[1].set_xlabel("time since entry [s]")
    fig.tight_layout(h_pad=0.6)
    return save(fig, out, "p1_data_examples.png")


def fig_coverage(eps, v1_eps, out):
    """Wide-envelope data vs. v1's data (real)."""
    fig, axes = plt.subplots(1, 3, figsize=(W, 2.35))
    icw = np.array([[e["ic"]["velocity_mps"], e["ic"]["gamma_deg"]] for e in eps])
    ico = np.array([[e["ic"]["velocity_mps"], e["ic"]["gamma_deg"]] for e in v1_eps])
    axes[0].scatter(icw[:, 0] / 1000, icw[:, 1], s=5, color=C_V2, alpha=0.7, label="v2 data")
    axes[0].scatter(ico[:, 0] / 1000, ico[:, 1], s=5, color=C_V1, alpha=0.7, label="v1 data")
    axes[0].set(xlabel="entry speed [km/s]", ylabel="entry flight-path angle [deg]", title="Where flights start")
    axes[0].legend(loc="lower left", markerscale=2)
    Yw = np.concatenate([e["y"][::10] for e in eps])
    Yo = np.concatenate([e["y"][::10] for e in v1_eps])
    axes[1].scatter(Yw[:, 3] / 1000, Yw[:, 0] / 1000, s=1, color=C_V2, alpha=0.3)
    axes[1].scatter(Yo[:, 3] / 1000, Yo[:, 0] / 1000, s=1, color=C_V1, alpha=0.3)
    axes[1].set(xlabel="speed [km/s]", ylabel="altitude [km]", title="States visited")
    bw = [e["domain"]["density_scale"] / e["domain"]["mass_scale"] for e in eps]
    bo = [e["domain"]["density_scale"] / e["domain"]["mass_scale"] for e in v1_eps]
    bins = np.linspace(0.6, 1.45, 30)
    axes[2].hist(bw, bins=bins, color=C_V2, alpha=0.75, density=True, label="v2 data")
    axes[2].hist(bo, bins=bins, color=C_V1, alpha=0.75, density=True, label="v1 data")
    for x in (0.78, 0.88, 1.12):
        axes[2].axvline(x, color=INK, lw=0.8, ls=":")
    axes[2].text(1.13, axes[2].get_ylim()[1] * 0.92, "test\natmospheres", fontsize=6.5, color=C_SEC, va="top")
    axes[2].set(xlabel="true density ÷ mass factor", ylabel="share of flights", title="Atmosphere randomization")
    fig.tight_layout(w_pad=1.2)
    return save(fig, out, "p1_coverage.png")


def fig_mppi_snapshot(out, case=0, t_decide=150):
    """
    One real MPPI decision of the shipped v1 model on frozen test case 0
    (nominal air): fly 150 s with v1-MPPI, then redo the next decision by hand
    to keep every one of its 256 samples.
    """
    import rl_env
    from jepa import models as M1, mppi as P
    torch.manual_seed(0)
    models = M1.load_checkpoint(str(SIM_DIR / "jepa_runs" / "full" / "jepa_ckpt.pt"))
    ic = D.load_test_ics(str(SIM_DIR / "configs" / "heldout_testset.csv"), n=case + 1)[case]
    env = rl_env.OrionEntryEnv(config_path=str(SIM_DIR / "configs" / "default.json"), controller_mode="policy",
                               dispersion=False, max_episode_steps=1400, fast_atmosphere=True)
    D.apply_domain(env.params, dict(D.NOMINAL_DOMAIN))
    env.reset(seed=0, options={"ic_override": ic})
    cfg = P.MPPIConfig(n_samples=256, n_sub=2)
    ctrl = P.LearnedModelController("jepa", models, cfg, env.target_phi_rad, env.target_lam_rad, env.drogue_alt_m)
    ctrl.reset(D.physical_state(env))
    sigma, past, cmds = 0.0, [D.physical_state(env)], []
    for k in range(t_decide):
        if k % 5 == 0:
            sigma = ctrl.act(D.physical_state(env))
        cmds.append(sigma)
        env.step(np.array([sigma / D.SIGMA_MAX_RAD], dtype=np.float32))
        past.append(D.physical_state(env))
    y_now = D.physical_state(env)
    y_hist_saved = ctrl.y_hist
    ctrl.y_hist = torch.cat([ctrl.y_hist[1:], torch.as_tensor(y_now, dtype=torch.float32)[None]])
    mp = ctrl.mppi
    basis = mp._noise_basis(cfg.horizon)
    eps_ = (torch.randn(cfg.n_samples, cfg.n_knots, generator=mp.gen) * math.radians(cfg.noise_std_deg)) @ basis
    eps_[0] = 0.0
    u = (mp.u_nom + eps_).clamp(-D.SIGMA_MAX_RAD, D.SIGMA_MAX_RAD)
    with torch.no_grad():
        o = ctrl._rollout(u)
        J, kc, miss = mp._cost(o, u, ctrl.u_hist[-1:])
        temp = max(cfg.temp_kappa * float(torch.quantile(J, cfg.temp_quantile) - J.min()), 1e-6)
        w = torch.softmax(-(J - J.min()) / temp, 0)
        u_new = (w[:, None] * u).sum(0)
        ob = ctrl._rollout(u_new[None])
        _, kcb, missb = mp._cost(ob, u_new[None], ctrl.u_hist[-1:])
    # Fly the same flight on to the parachute with the normal planner, to report how it really ends.
    ctrl.y_hist = y_hist_saved
    k, done, fin = t_decide, False, {}
    while not done:
        if k % 5 == 0:
            sigma = ctrl.act(D.physical_state(env))
        _, _, term, trunc, fin = env.step(np.array([sigma / D.SIGMA_MAX_RAD], dtype=np.float32))
        k += 1
        done = bool(term or trunc)
    order = torch.argsort(J, descending=True)
    rank = torch.empty_like(J)
    rank[order] = torch.linspace(0, 1, len(J))
    cmap = LinearSegmentedColormap.from_list("cost", ["#f9dccd", "#eb6834", "#7a2d0c"])
    fig, axes = plt.subplots(1, 2, figsize=(W, 2.75), gridspec_kw={"width_ratios": [1.15, 1]})
    past = np.array(past)
    tp = np.arange(len(cmds))
    axes[0].step(tp, np.degrees(cmds), where="post", color=INK, lw=1.2, label="flown so far")
    tf = t_decide + np.arange(cfg.horizon) * 5.0
    for i in order.tolist():
        k = int(kc[i])
        axes[0].plot(tf[:k], np.degrees(u[i, :k]), color=cmap(float(rank[i])), lw=0.5, alpha=0.8)
    axes[0].plot(tf[:int(kcb[0])], np.degrees(u_new[:int(kcb[0])]), color=INK, lw=2.0, ls="--",
                 label="weighted average plan")
    axes[0].axvline(t_decide, color=C_SEC, lw=0.8, ls=":")
    axes[0].set(xlabel="time since entry [s]", ylabel="bank command [deg]", ylim=(-95, 95),
                title=f"256 candidate bank plans at t = {t_decide} s")
    axes[0].legend(loc="lower left", fontsize=6.5)
    y6 = o["y6"].numpy()
    for i in order.tolist():
        k = int(kc[i]) + 1
        axes[1].plot(np.degrees(y6[i, :k, 2]), np.degrees(y6[i, :k, 1]), color=cmap(float(rank[i])), lw=0.5, alpha=0.8)
    axes[1].plot(np.degrees(past[:, 2]), np.degrees(past[:, 1]), color=INK, lw=1.4, label="flown so far")
    yb = ob["y6"][0].numpy()
    axes[1].plot(np.degrees(yb[:int(kcb[0]) + 1, 2]), np.degrees(yb[:int(kcb[0]) + 1, 1]), color=INK, lw=2.0,
                 ls="--", label=f"weighted plan: predicted miss {float(missb[0]):.1f} km")
    tphi, tlam = np.degrees(env.target_phi_rad), np.degrees(env.target_lam_rad)
    axes[1].plot([tlam], [tphi], marker="*", ms=12, color="#0b0b0b", ls="none", label="target")
    axes[1].set(xlabel="longitude [deg]", ylabel="latitude [deg]", title="Where the model says each plan lands (zoom)",
                xlim=(tlam - 3.0, tlam + 1.2), ylim=(tphi - 3.0, tphi + 1.2))
    axes[1].legend(loc="lower right", fontsize=6.0)
    sm = plt.cm.ScalarMappable(cmap=cmap)
    cb = fig.colorbar(sm, ax=axes[1], fraction=0.05, pad=0.02, ticks=[0, 1])
    cb.ax.set_yticklabels(["high\ncost", "low\ncost"], fontsize=6.3)
    fig.tight_layout(w_pad=1.0)
    info = {"t_decide": t_decide, "miss_best": float(missb[0]), "n": cfg.n_samples,
            "miss_med": float(miss.median()), "ess": float(1.0 / (w ** 2).sum()),
            "final_miss": float(fin.get("miss_km", float("nan"))), "final_t": float(fin.get("t_s", float("nan"))),
            "final_outcome": fin.get("outcome")}
    return save(fig, out, "p1_mppi_snapshot.png"), info


def fig_mppi_loop(out):
    """The v2 planning loop (illustration)."""
    fig, ax = diagram(3.3)
    fs = 6.0
    panel(ax, 1, 1, 70, 46, "Planner (every 5 s, one batched GPU pass)", ORANGE)
    panel(ax, 72, 1, 99, 46, "Real simulator", GRAY)
    box(ax, 14, 35, 22, 8, "history: last 25 states\n+ 24 bank commands\n(2 minutes)", GRAY, fs)
    box(ax, 14, 19, 22, 9, "256 bank plans =\nlast best plan +\nsmooth random changes", ORANGE, fs)
    box(ax, 37.5, 27, 18, 15, "5 world models\n(the ensemble)\n\neach one flies\nevery plan down to\nparachute altitude", BLUE, fs)
    box(ax, 59, 34.5, 19.5, 15, "cost of a plan =\naverage predicted miss\n+ disagreement of the 5\n+ smoothness\n+ too many g\n+ too steep a dive",
        YELLOW, fs)
    box(ax, 59, 14.5, 19.5, 11, "weight plans by\nexp(−cost), average\nthem → new best plan\n(kept for next time)", ORANGE, fs)
    box(ax, 85.5, 37, 23, 7, "fly the first 5 s\nof the best plan", ORANGE, fs)
    box(ax, 85.5, 24.5, 23, 9, "actuator → roll control\n→ RCS jets\n→ 3D physics", GRAY, fs)
    box(ax, 85.5, 12, 23, 7, "new state + measured\ndrag and lift", GRAY, fs)
    arrow(ax, (25, 35), (28.5, 31))
    arrow(ax, (25, 19), (28.5, 23))
    arrow(ax, (46.5, 31), (49.25, 33))
    arrow(ax, (59, 27), (59, 20))
    ax.plot([68.75, 71, 71], [14.5, 14.5, 37], color=INK, lw=1.1, zorder=1)
    arrow(ax, (71, 37), (74, 37))
    arrow(ax, (85.5, 33.5), (85.5, 29))
    arrow(ax, (85.5, 20), (85.5, 15.5))
    ax.plot([74, 72, 72, 2.5, 2.5], [12, 12, 4.5, 4.5, 35], color=INK, lw=1.1, zorder=1)
    arrow(ax, (2.5, 35), (3, 35))
    ax.text(36, 6.3, "replan from the new history every 5 s", fontsize=6.0, color=C_SEC, ha="center")
    return save(fig, out, "p1_mppi_loop.png")


def fig_v1_results(out):
    """v1's published results (PLAN.md section 9, 16 cases, one seed)."""
    fig, axes = plt.subplots(1, 2, figsize=(W, 2.6))
    ttg = [400, 300, 200, 100, 50]
    series = {"JEPA+PI (v1)": ([40.5, 19.7, 6.0, 1.1, 0.31], C_V1, "-"),
              "physics prior only": ([57.2, 36.0, 11.3, 0.9, 0.47], C_MUTED, ":"),
              "predictive model": ([53.8, 27.1, 9.3, 4.4, 2.58], C_V1W, "-"),
              "classical's physics model": ([32.3, 22.6, 5.7, 0.6, 0.12], C_CL, "--")}
    for name, (v, c, ls) in series.items():
        axes[0].plot(ttg, v, marker="o", ms=4, color=c, ls=ls, label=name)
    axes[0].set(xscale="linear", yscale="log", xlabel="time left before the parachute [s]",
                ylabel="median landing-point error [km]", title="Open loop: predicting where it lands")
    axes[0].invert_xaxis()
    axes[0].legend(fontsize=6.3, loc="lower left")
    ctrl = ["classical", "JEPA-MPPI (v1)", "predictive-MPPI"]
    cols = [C_CL, C_V1, C_V1W]
    med = {"nominal air": [2.24, 1.69, 13.1], "thin air (×0.88)": [3.85, 5.07, 10.2]}
    succ = {"nominal air": [100, 100, 87.5], "thin air (×0.88)": [87.5, 93.8, 81.3]}
    x = np.arange(2)
    bw = 0.26
    for j, (c, col) in enumerate(zip(ctrl, cols)):
        vals = [med[k][j] for k in med]
        bars = axes[1].bar(x + (j - 1) * bw, vals, bw - 0.03, color=col, label=c)
        for b, k in zip(bars, med):
            axes[1].text(b.get_x() + b.get_width() / 2, b.get_height() + 0.25, f"{succ[k][j]:.0f}%",
                         ha="center", fontsize=6.2, color=C_SEC)
    axes[1].set_xticks(x, list(med))
    axes[1].set(ylabel="median miss [km]", title="Closed loop: 16 test cases (labels: success rate)")
    axes[1].legend(fontsize=6.3, loc="upper left")
    fig.tight_layout(w_pad=1.5)
    return save(fig, out, "p1_v1_results.png")


def fig_dense(runs, out):
    """Real: the shipped v1 model vs. the classical controller in 12% denser air, 8 frozen test cases."""
    import pandas as pd
    df = pd.DataFrame([{k: v for k, v in r.items() if k not in ("traj", "target")} for r in runs])
    fig, axes = plt.subplots(1, 2, figsize=(W, 2.7), gridspec_kw={"width_ratios": [1, 1.05]})
    for c, col, mk, lab in (("classical", C_CL, "s", "classical"), ("jepa", C_V1, "o", "v1 JEPA-MPPI")):
        d = df[df.controller == c].sort_values("case")
        axes[0].plot(d.case, d.miss_km, marker=mk, ms=5, ls="none", color=col,
                     label=f"{lab}: median {d.miss_km.median():.1f} km")
    axes[0].axhline(25, color=C_SEC, lw=0.8, ls=":")
    axes[0].text(7.4, 19, "success radius 25 km", fontsize=6.3, color=C_SEC, ha="right")
    axes[0].set(yscale="log", ylim=(0.4, 120), xlabel="test case", ylabel="miss at parachute [km]",
                title="Density and L/D ×1.12: miss per case")
    axes[0].set_xticks(range(8))
    axes[0].legend(fontsize=6.5, loc="upper left")
    worst = int(df[df.controller == "jepa"].sort_values("miss_km").iloc[-1].case)
    for c, col, lab in (("classical", C_CL, "classical"), ("jepa", C_V1, "v1")):
        r = next(r for r in runs if r["controller"] == c and r["case"] == worst)
        T = r["traj"]
        axes[1].plot(np.degrees(T[:, 2]), np.degrees(T[:, 1]), color=col, label=f"{lab}: {r['miss_km']:.1f} km")
    tphi, tlam = np.degrees(r["target"])
    axes[1].plot([tlam], [tphi], marker="*", ms=12, color="#0b0b0b", ls="none", label="target")
    axes[1].set(xlim=(tlam - 1.3, tlam + 0.35), ylim=(tphi - 1.3, tphi + 0.35), xlabel="longitude [deg]",
                ylabel="latitude [deg]", title=f"Case {worst}: last part of the ground track (1° ≈ 111 km)")
    axes[1].legend(fontsize=6.5, loc="lower right")
    fig.tight_layout(w_pad=1.5)
    summ = df.groupby("controller").miss_km.agg(["median", "max"]).to_dict("index")
    succ = df.groupby("controller").success.mean().to_dict()
    per_case = {c: df[df.controller == c].sort_values("case").miss_km.tolist() for c in ("classical", "jepa")}
    return save(fig, out, "p1_dense.png"), {"summary": summ, "success": succ, "worst": worst, "per_case": per_case}


def fig_compute(out):
    """Estimated notebook time and cost per GPU (estimates from PLAN.md 11.6, not measured)."""
    stages = ["setup + data", "λ sweep", "stage 1", "stage 2", "open-loop tests", "closed loop"]
    mins = {"T4": [3, 9, 75, 12, 6, 32], "L4": [3, 5, 30, 10, 5, 15], "A100": [3, 3, 12, 8, 4, 12]}
    rate = {"T4": 2.0, "L4": 4.8, "A100": 13.0}         # approx. Colab compute units per hour
    cols = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
    fig, axes = plt.subplots(1, 2, figsize=(W, 2.2), gridspec_kw={"width_ratios": [1.6, 1]})
    gpus = list(mins)
    for i, g in enumerate(gpus):
        left = 0
        for j, (s, m) in enumerate(zip(stages, mins[g])):
            axes[0].barh(i, m, left=left, color=cols[j], edgecolor="white", linewidth=1.0,
                         label=s if i == 0 else None, height=0.6)
            left += m
        axes[0].text(left + 2, i, f"{left / 60:.1f} h", va="center", fontsize=7, color=INK)
    axes[0].set_yticks(range(3), gpus)
    axes[0].set(xlabel="minutes (estimate)", title="Estimated wall time of the full notebook", xlim=(0, 175))
    axes[0].legend(ncol=3, fontsize=6.2, loc="upper center", bbox_to_anchor=(0.5, -0.38))
    axes[0].grid(axis="y", visible=False)
    cost = [sum(mins[g]) / 60 * rate[g] for g in gpus]
    bars = axes[1].bar(gpus, cost, color=[C_MUTED, C_V2, C_MUTED], width=0.55)
    for b, c in zip(bars, cost):
        axes[1].text(b.get_x() + b.get_width() / 2, c + 0.6, f"~{c:.0f}", ha="center", fontsize=7, color=INK)
    axes[1].set(ylabel="Colab compute units", title="Estimated cost (orange = recommended)", ylim=(0, 16))
    fig.tight_layout(w_pad=1.5)
    return save(fig, out, "p1_compute.png")


# ============================================================================
# Part 2: the architecture
# ============================================================================
def fig_architecture(out):
    """Stage 1 of v2: encoder, action tokens, causal transformer predictor, losses (illustration)."""
    fig, ax = diagram(4.0)
    fs = 6.0
    panel(ax, 1, 1, 99, 56, "One ensemble member (v2 trains 5 of these side by side)", BLUE)
    box(ax, 12, 45, 20, 8, "state history\n25 steps × 15 features\n(120 s)", GRAY, fs)
    box(ax, 12, 31.5, 20, 10, "state encoder:\nbanded causal\ntransformer\n(3 layers, width 128)", BLUE, fs)
    box(ax, 12, 19.5, 14, 6, "code $s_t$\n(64 numbers)", BLUE, fs)
    box(ax, 12, 7.5, 20, 8, "bank commands:\nlast 4 per step\n→ action tokens", AQUA, fs)
    arrow(ax, (12, 41), (12, 36.5))
    arrow(ax, (12, 26.5), (12, 22.5))
    xs = [32, 42, 52, 62, 77]
    labels = ["$s_t$", "$a_t$", "$a_{t+1}$", "$a_{t+2}$", "$a_{t+K-1}$"]
    for x, lab, col in zip(xs, labels, [BLUE, AQUA, AQUA, AQUA, AQUA]):
        box(ax, x, 14, 8, 5, lab, col, fs + 0.6)
    ax.text(69.5, 14, "…", fontsize=10, ha="center", va="center", color=C_SEC)
    box(ax, 54.5, 24.5, 54, 6.5, "predictor: causal transformer (4 layers, width 192)\n"
                                 "each output only looks left, at the past", ORANGE, fs)
    for x in xs:
        arrow(ax, (x, 16.6), (x, 21.1), lw=0.8)
    outs = ["$\\tilde{s}_{t+1}$", "$\\tilde{s}_{t+2}$", "$\\tilde{s}_{t+3}$", "$\\tilde{s}_{t+K}$"]
    reals = ["$s_{t+1}$", "$s_{t+2}$", "$s_{t+3}$", "$s_{t+K}$"]
    for x, lab, real in zip(xs[1:], outs, reals):
        arrow(ax, (x, 27.9), (x, 32.4), lw=0.8)
        box(ax, x, 35, 8, 5, lab, ORANGE, fs + 0.6)
        box(ax, x, 48, 8, 5, real, BLUE, fs + 0.6, dashed=True)
        ax.annotate("", xy=(x, 37.6), xytext=(x, 45.4), arrowprops=dict(arrowstyle="<->", color=YELLOW[0], lw=1.1))
    ax.text(69.5, 35, "…", fontsize=10, ha="center", va="center", color=C_SEC)
    ax.text(69.5, 48, "…", fontsize=10, ha="center", va="center", color=C_SEC)
    ax.text(83, 48, "real future codes\n(same encoder,\nlater histories)", fontsize=5.8, color=BLUE[0], va="center")
    ax.text(83, 41.5, "predicted ≈ real\n= the loss", fontsize=6.0, color=YELLOW[0], va="center", fontweight="bold")
    arrow(ax, (19, 19.5), (27.5, 15.5), rad=-0.15, lw=0.9)
    arrow(ax, (22, 7.5), (37.8, 12.2), rad=0.1, col=AQUA[0], lw=0.9)
    box(ax, 76, 6.0, 40, 7.5, "+ SIGReg: keep the codes spread out like a\nbell curve in every direction (no collapse)",
        YELLOW, fs)
    return save(fig, out, "p2_architecture.png")


def fig_features(eps, out):
    """Real: the measured drag correction reads the true air density."""
    fig, axes = plt.subplots(1, 2, figsize=(W, 2.6), gridspec_kw={"width_ratios": [1.25, 1]})
    norm = plt.Normalize(0.7, 1.3)
    pts = []
    for e in eps:
        y = torch.from_numpy(e["y"]).float()
        f = V.raw_features_v2(y).numpy()
        h = e["y"][:, 0] / 1000
        cd = f[:, V.FEAT_NAMES_V2.index("cD_meas")]
        beta = e["domain"]["density_scale"] / e["domain"]["mass_scale"]
        ok = (h > 12) & (h < 100) & (e["y"][:, 8] > 0.05)
        axes[0].plot(cd[ok], h[ok], color=DIVERGE(norm(beta)), lw=0.6, alpha=0.8)
        band = ok & (h > 35) & (h < 60)
        if band.sum() > 5:
            pts.append((beta, float(np.median(cd[band]))))
    axes[0].set(xlabel="measured drag correction  cD_meas", ylabel="altitude [km]",
                title=f"{len(eps)} real flights, colored by true density")
    sm = plt.cm.ScalarMappable(cmap=DIVERGE, norm=norm)
    cb = fig.colorbar(sm, ax=axes[0], fraction=0.05, pad=0.02)
    cb.set_label("true density ÷ mass factor", fontsize=6.5)
    pts = np.array(pts)
    axes[1].scatter(np.log(pts[:, 0]), pts[:, 1], s=9, color=C_V2, alpha=0.8)
    a, b = np.polyfit(np.log(pts[:, 0]), pts[:, 1], 1)
    xx = np.linspace(np.log(0.65), np.log(1.45), 10)
    axes[1].plot(xx, a * xx + b, color=INK, lw=1, ls="--", label=f"fit: slope {a:.2f}")
    r = np.corrcoef(np.log(pts[:, 0]), pts[:, 1])[0, 1]
    axes[1].set(xlabel="log(true density ÷ mass factor)", ylabel="median cD_meas, 35–60 km",
                title=f"The feature reads the truth (r = {r:.2f})")
    axes[1].legend(loc="upper left")
    fig.tight_layout(w_pad=1.2)
    return save(fig, out, "p2_features.png"), {"slope": float(a), "r": float(r), "n": len(pts)}


def fig_masks(out):
    """From the code: encoder band mask, its 3-layer receptive field, the predictor's causal mask."""
    L = 40
    m = V.band_mask(L, 9, "cpu").float().numpy()
    reach = m.copy()
    for _ in range(2):
        reach = (reach @ m > 0).astype(float)
    fig, axes = plt.subplots(1, 3, figsize=(W, 2.45))
    cm = LinearSegmentedColormap.from_list("on", ["#fbfbfa", "#2a78d6"])
    axes[0].imshow(m, cmap=cm, interpolation="nearest")
    axes[0].set(title="one encoder layer:\nlooks back 8 steps", xlabel="looks at step", ylabel="step")
    axes[1].imshow(reach, cmap=cm, interpolation="nearest")
    axes[1].set(title=f"after 3 layers: exactly\n{int(reach[-1].sum())} steps (120 s)", xlabel="looks at step")
    n = 13
    cmask = np.tril(np.ones((n, n)))
    axes[2].imshow(cmask, cmap=LinearSegmentedColormap.from_list("on2", ["#fbfbfa", "#eb6834"]), interpolation="nearest")
    axes[2].set(title="predictor: each output looks\nat the code + earlier actions", xlabel="input token")
    axes[2].set_xticks([0, 4, 8, 12], ["s_t", "a+3", "a+7", "a+11"])
    axes[2].set_yticks([0, 4, 8, 12], ["", "s̃+4", "s̃+8", "s̃+12"])
    for ax in axes:
        ax.grid(False)
    fig.tight_layout(w_pad=1.0)
    return save(fig, out, "p2_masks.png")


def fig_action_tokens(out):
    """Illustration: each action token packs the current and 3 previous bank commands."""
    fig, ax = plt.subplots(figsize=(W, 1.9))
    rng = np.random.default_rng(3)
    u = np.clip(np.cumsum(rng.normal(0, 25, 14)) + 20, -85, 85)
    t = np.arange(14)
    ax.step(np.r_[t, 14] * 5, np.r_[u, u[-1]], where="post", color=AQUA[0], lw=1.6)
    for j, (k, col) in enumerate(((8, ORANGE), (11, BLUE))):
        x0, x1 = (k - 3) * 5, (k + 1) * 5
        ax.add_patch(Rectangle((x0, -98), x1 - x0, 196, fc=col[1], ec=col[0], lw=1, alpha=0.6, zorder=0))
        ax.text((x0 + x1) / 2, 103, f"token for step {k}:\ncommands {k - 3}…{k}", ha="center", fontsize=6.5,
                color=col[0], va="bottom")
    ax.set(xlabel="time [s] (one command is held for 5 s)", ylabel="bank cmd [deg]", ylim=(-100, 140), xlim=(0, 70))
    ax.set_yticks([-90, 0, 90])
    return save(fig, out, "p2_action_tokens.png")


def fig_parallel(out):
    """Illustration: GRU (one step at a time) vs. causal transformer (all steps in one pass)."""
    fig, ax = diagram(2.1)
    fs = 6.5
    panel(ax, 1, 1, 49, 28.5, "v1: GRU, one step after another", GRAY)
    panel(ax, 51, 1, 99, 28.5, "v2: transformer, all steps in one pass", ORANGE)
    for i in range(4):
        x = 8 + i * 11
        box(ax, x, 13, 8, 6, "GRU" if i < 3 else "…", ORANGE if i < 3 else GRAY, fs)
        if i < 3:
            arrow(ax, (x + 4, 13), (x + 7, 13))
            box(ax, x, 22, 7, 4, f"s̃+{i + 1}", ORANGE, fs)
            arrow(ax, (x, 16), (x, 20))
            box(ax, x, 5, 7, 4, f"a+{i}", AQUA, fs)
            arrow(ax, (x, 7), (x, 10))
    box(ax, 75, 13, 44, 6, "causal transformer: each output sees only the left", ORANGE, fs)
    for i in range(5):
        x = 57 + i * 9
        box(ax, x, 22, 7, 4, f"s̃+{i + 1}" if i < 4 else "…", ORANGE, fs)
        arrow(ax, (x, 16), (x, 20))
        box(ax, x, 5, 7, 4, "s_t" if i == 0 else f"a+{i - 1}", BLUE if i == 0 else AQUA, fs)
        arrow(ax, (x, 7), (x, 10))
    return save(fig, out, "p2_parallel.png")


def fig_sigreg(out):
    """Illustration with synthetic numbers, scored by the real SIGReg code."""
    rng = np.random.default_rng(0)
    n = 1024
    sets = {"collapsed to a point": rng.normal(0, 0.03, (n, 2)) + 0.4,
            "collapsed to a line": np.c_[rng.normal(0, 1, n), 0.05 * rng.normal(0, 1, n)] @ np.array([[0.8, 0.6], [-0.6, 0.8]]),
            "healthy: round bell curve": rng.normal(0, 1, (n, 2))}
    sig = V.SIGRegV2(n_proj=64)
    fig = plt.figure(figsize=(W, 2.3))
    scores = {}
    for i, (name, x) in enumerate(sets.items()):
        ax = fig.add_axes([0.02 + i * 0.21, 0.12, 0.19, 0.7])
        col = C_V2 if i < 2 else C_V1
        ax.scatter(x[:, 0], x[:, 1], s=1.5, color=col, alpha=0.5)
        ax.set(xlim=(-3.5, 3.5), ylim=(-3.5, 3.5), title=name)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(True)
            sp.set_color("#d9d8d3")
        torch.manual_seed(1)
        scores[name] = float(sig(torch.tensor(x, dtype=torch.float32)[None])[0][0])
    ax = fig.add_axes([0.70, 0.18, 0.29, 0.62])
    names = list(scores)
    ax.barh(range(3), [scores[k] for k in names], color=[C_V2, C_V2, C_V1], height=0.55)
    for i, k in enumerate(names):
        ax.text(scores[k] * 1.15, i, f"{scores[k]:.3g}", va="center", fontsize=6.5)
    ax.set_yticks(range(3), ["point", "line", "healthy"])
    ax.set(xscale="log", xlabel="SIGReg penalty (log)", title="What SIGReg charges")
    ax.grid(axis="y", visible=False)
    return save(fig, out, "p2_sigreg.png"), scores


def fig_schedule(out):
    """From the code: the horizon curriculum and the learning-rate schedule (full mode)."""
    wm = V.TrainCfgV2(steps=24_000, warmup=1_000, lr=6e-4, curriculum=((0.0, 24), (0.2, 48), (0.5, 96)))
    pi = V.TrainCfgV2(steps=6_000, warmup=300, lr=1e-3, curriculum=((0.0, 16), (0.3, 32), (0.6, 48)))
    fig, axes = plt.subplots(1, 2, figsize=(W, 2.1))
    for cfg, col, lab in ((wm, C_V2, "stage 1: world model"), (pi, VIOLET[0], "stage 2: prober")):
        s = np.arange(cfg.steps)
        frac = s / cfg.steps
        axes[0].step(frac * 100, [V.horizon_at(i, cfg) for i in s], where="post", color=col, label=lab)
        axes[1].plot(frac * 100, [V._lr(i, cfg) * 1e4 for i in s], color=col, label=lab)
    axes[0].set(xlabel="training progress [%]", ylabel="prediction horizon T [steps]",
                title="Curriculum: short first, then long (×5 s)")
    axes[0].set_yticks([16, 24, 32, 48, 96])
    axes[0].legend(loc="upper left")
    axes[1].set(xlabel="training progress [%]", ylabel="learning rate [×1e-4]", title="Warm-up, then cosine decay")
    fig.tight_layout(w_pad=1.5)
    return save(fig, out, "p2_schedule.png")


def fig_masked_windows(eps, v1_eps, out):
    """Real flight lengths vs. the training horizons, and a masked window (illustration)."""
    fig, axes = plt.subplots(1, 2, figsize=(W, 2.3), gridspec_kw={"width_ratios": [1.1, 1]})
    Lw = np.array([len(e["cmd"]) // 5 for e in eps])
    Lo = np.array([len(e["cmd"]) // 5 for e in v1_eps])
    bins = np.arange(30, 250, 8)
    axes[0].hist(Lw, bins=bins, color=C_V2, alpha=0.75, label="v2 data")
    axes[0].hist(Lo, bins=bins, color=C_V1, alpha=0.6, label="v1 data")
    for T in (24, 48, 96):
        axes[0].axvline(T, color=INK, lw=0.8, ls=":")
        axes[0].text(T + 1.5, axes[0].get_ylim()[1] * 0.95, f"T={T}", fontsize=6.5, va="top")
    share = (Lw < 96).mean() * 100
    axes[0].set(xlabel="flight length [5 s steps]", ylabel="flights",
                title=f"{share:.0f}% of these v2 flights are shorter than 96 steps")
    axes[0].legend(loc="upper right")
    ax = axes[1]
    ax.set_xlim(0, 30)
    ax.set_ylim(0, 7)
    ax.axis("off")
    ax.set_title("A training window that runs past the landing", fontsize=8.5)
    x0 = 3.0
    for i in range(22):
        real = i < 14
        ax.add_patch(Rectangle((x0 + i * 1.2, 3.4), 1.05, 1.3, fc=ORANGE[1] if real else "#f2f1ed",
                               ec=ORANGE[0] if real else C_MUTED, lw=0.8, hatch=None if real else "////"))
        ax.text(x0 + 0.52 + i * 1.2, 2.5, "1" if real else "0", fontsize=6, ha="center", color=INK if real else C_MUTED)
    ax.text(x0, 5.3, "real future steps (in the loss)", fontsize=6.3, color=ORANGE[0])
    ax.text(x0 + 14 * 1.2, 5.3, "after landing (masked)", fontsize=6.3, color=C_SEC)
    ax.text(0.2, 2.5, "mask", fontsize=6.3, color=C_SEC)
    ax.text(0.2, 0.8, "Short flights still teach long horizons\nfor the part of the window they cover.",
            fontsize=6.3, color=C_SEC, va="bottom")
    fig.tight_layout(w_pad=1.0)
    return save(fig, out, "p2_masked_windows.png")


def fig_prober(out):
    """Stage 2: the state-aware physics-inspired prober (illustration)."""
    fig, ax = diagram(3.0)
    fs = 6.0
    panel(ax, 1, 1, 99, 41.5, "Stage 2: from predicted codes to a physical trajectory", VIOLET)
    box(ax, 11, 30, 17, 7, "predicted code $\\tilde{s}$\n(stage 1, frozen)", ORANGE, fs, dashed=True)
    box(ax, 11, 15, 17, 7, "where the rollout is:\naltitude h, speed V", GRAY, fs)
    box(ax, 30, 22.5, 12, 7, "small MLP\n(3 layers)", VIOLET, fs, bold=True)
    box(ax, 49, 22.5, 18, 13, "3 corrections:\nlog drag factor $c_D$\nlog L/D factor $c_L$\nbank offset Δσ", VIOLET, fs)
    box(ax, 44, 6.5, 41, 6.5, "crude prior: one exponential atmosphere, one\nballistic coefficient, L/D = 0.30, lift along the command",
        GRAY, fs)
    box(ax, 72, 22.5, 14, 11, "exact 3D\nequations of\nmotion", VIOLET, fs)
    box(ax, 90, 22.5, 17, 11, "next altitude,\nlatitude, longitude,\nspeed, γ, heading", GRAY, fs)
    arrow(ax, (19.5, 30), (24, 25))
    arrow(ax, (19.5, 15), (24, 20))
    arrow(ax, (36, 22.5), (40, 22.5))
    arrow(ax, (58, 22.5), (65, 22.5))
    arrow(ax, (65, 7.5), (70, 17), rad=0.25)
    arrow(ax, (79, 22.5), (81.4, 22.5))
    ax.annotate("", xy=(19.6, 12.6), xytext=(90, 17),
                arrowprops=dict(arrowstyle="-|>", color=C_MUTED, lw=1, ls="--",
                                connectionstyle="angle,angleA=-90,angleB=0,rad=0"))
    ax.text(88.5, 10.4, "feeds the next step", fontsize=5.8, color=C_SEC, ha="right")
    ax.text(50, 36.5, "New in v2: the correction also sees where the rollout actually is (h, V).",
            fontsize=6.0, color=C_SEC, ha="center")
    return save(fig, out, "p2_prober.png")


def fig_ensemble(out):
    """Illustration with synthetic points: how the ensemble's disagreement enters the plan's cost."""
    rng = np.random.default_rng(4)
    fig, axes = plt.subplots(1, 2, figsize=(W, 2.5), gridspec_kw={"width_ratios": [1.3, 1]})
    ax = axes[0]
    target = np.array([0.0, 0.0])
    A = np.array([3.0, 2.0]) + rng.normal(0, 0.6, (5, 2))
    B = np.array([-1.5, 0.8]) + rng.normal(0, 4.5, (5, 2))
    for P_, col, name in ((A, C_V1, "plan A"), (B, C_V2, "plan B")):
        mean = P_.mean(0)
        spread = np.sqrt(((P_ - mean) ** 2).sum(1).mean())
        ax.add_patch(Circle(mean, spread, fc="none", ec=col, lw=1, ls="--"))
        ax.scatter(P_[:, 0], P_[:, 1], s=22, color=col, label=f"{name}: 5 models' landing points", zorder=3)
    ax.plot(*target, marker="*", ms=13, color="#0b0b0b", ls="none", label="target")
    ax.set(xlim=(-11, 11), ylim=(-9, 9), xlabel="east [km]", ylabel="north [km]", title="Each plan, flown by all 5 models")
    ax.set_aspect("equal")
    ax.legend(fontsize=6.3, loc="lower right")
    def cost(P_):
        miss = np.mean(np.hypot(P_[:, 0], P_[:, 1]))
        mean = P_.mean(0)
        return miss, np.sqrt(((P_ - mean) ** 2).sum(1).mean())
    (ma, sa), (mb, sb) = cost(A), cost(B)
    ax = axes[1]
    ax.barh([1, 0], [ma, mb], color=[C_V1, C_V2], height=0.5, label="average predicted miss")
    ax.barh([1, 0], [sa, sb], left=[ma, mb], color=[C_V1, C_V2], alpha=0.35, height=0.5, hatch="///",
            label="disagreement (spread)")
    for y, m, s in ((1, ma, sa), (0, mb, sb)):
        ax.text(m + s + 0.3, y, f"{m + s:.1f}", va="center", fontsize=7)
    ax.set_yticks([1, 0], ["plan A", "plan B"])
    ax.set(xlabel="cost [km]", title="MPPI prefers the plan the models agree on", xlim=(0, 22))
    ax.legend(fontsize=6.3, loc="upper right")
    ax.grid(axis="y", visible=False)
    fig.tight_layout(w_pad=1.5)
    return save(fig, out, "p2_ensemble.png")


def fig_params(out):
    """From the code: parameter counts (v1 vs v2) and training compute per step vs horizon."""
    from jepa import models as M1
    v1 = M1.JEPAWorldModel()
    groups_v1 = {"state encoder": sum(p.numel() for m in (v1.state_tcn, v1.state_head) for p in m.parameters()),
                 "action encoder": sum(p.numel() for m in (v1.act_tcn, v1.act_head) for p in m.parameters()),
                 "predictor": sum(p.numel() for m in (v1.h_init, v1.gru, v1.readout) for p in m.parameters()),
                 "prober": M1.n_params(M1.PIProber(32))}
    wm = V.WorldModelV2()
    rep = V.param_report(wm, V.PIProberV2(5, 64))
    groups_v2 = {"state encoder": rep["state encoder"], "action encoder": rep["action tokens"],
                 "predictor": rep["predictor"], "prober": rep["PI prober"]}
    fig, axes = plt.subplots(1, 2, figsize=(W, 2.35), gridspec_kw={"width_ratios": [1.3, 1]})
    names = list(groups_v1)
    y = np.arange(len(names))
    axes[0].barh(y + 0.18, [groups_v1[k] for k in names], 0.34, color=C_V1, label="v1 (one model)")
    axes[0].barh(y - 0.18, [groups_v2[k] for k in names], 0.34, color=C_V2, label="v2 (per ensemble member)")
    for i, k in enumerate(names):
        axes[0].text(groups_v2[k] * 1.2, i - 0.18, f"{groups_v2[k]:,}", va="center", fontsize=6.2)
        axes[0].text(groups_v1[k] * 1.2, i + 0.18, f"{groups_v1[k]:,}", va="center", fontsize=6.2, color=C_SEC)
    axes[0].set_yticks(y, names)
    axes[0].invert_yaxis()
    axes[0].set(xscale="log", xlabel="parameters (log scale)", xlim=(1e3, 3e7),
                title=f"Size: v1 {sum(groups_v1.values()):,} vs v2 {rep['ensemble (x5)']:,} (5 members)")
    axes[0].legend(loc="upper center", bbox_to_anchor=(0.45, -0.25), ncol=2, fontsize=6.3)
    axes[0].grid(axis="y", visible=False)
    Ts = np.arange(8, 100, 4)
    axes[1].plot(Ts, [V.train_flops_per_step(wm, 256, int(T)) / 1e12 for T in Ts], color=C_V2)
    for T in (24, 48, 96):
        axes[1].plot([T], [V.train_flops_per_step(wm, 256, T) / 1e12], marker="o", ms=4, color=C_V2)
    axes[1].set(xlabel="prediction horizon T [steps]", ylabel="TFLOP per training step",
                title="Training cost per step (5 members, 256 windows each)")
    fig.tight_layout(w_pad=1.5)
    return save(fig, out, "p2_params.png"), {"v1": groups_v1, "v2": groups_v2, "rep": rep}


def run_dense_cases(n_cases: int = 8, workers: int = 4) -> list:
    """Fly the classical controller and the shipped v1 model under density and L/D x1.12 (frozen test ICs)."""
    from jepa import mppi as P
    ics = D.load_test_ics(str(SIM_DIR / "configs" / "heldout_testset.csv"), n=n_cases)
    dom = {"density_scale": 1.12, "lift_scale": 1.12, "mass_scale": 1.0}
    cfg = str(SIM_DIR / "configs" / "default.json")
    jobs = []
    for i, ic in enumerate(ics):
        jobs.append({"controller": "classical", "ic": ic, "domain": dom, "tag": "mismatch_high", "case": i,
                     "record": True, "config_path": cfg})
        jobs.append({"controller": "jepa", "ic": ic, "domain": dom, "ckpt": str(SIM_DIR / "jepa_runs" / "full" / "jepa_ckpt.pt"),
                     "mppi": {"n_samples": 256, "n_sub": 2}, "tag": "mismatch_high", "case": i, "record": True,
                     "config_path": cfg})
    jobs.sort(key=lambda j: j["controller"] != "classical")
    return D.run_parallel(P.run_episode, jobs, workers, label="dense air")


# ============================================================================
def make_all(wide_path: str, dense_path: str, out: str, n_wide: int = 600) -> dict:
    """Render every figure; returns {name: path} plus a few numbers quoted in the text."""
    from jepa.train_v2 import load_episodes
    os.makedirs(out, exist_ok=True)
    eps = load_episodes(wide_path)[:n_wide]
    v1_eps = pickle.load(open(SIM_DIR / "jepa_runs" / "episodes_900.pkl", "rb"))
    figs, info = {}, {"n_wide": len(eps), "n_v1": len(v1_eps)}
    figs["basics"] = fig_basics(eps, out)
    figs["pipeline"] = fig_pipeline(out)
    figs["jepa_idea"] = fig_jepa_idea(out)
    figs["data_examples"] = fig_data_examples(eps, out)
    figs["coverage"] = fig_coverage(eps, v1_eps, out)
    figs["mppi_snapshot"], info["mppi"] = fig_mppi_snapshot(out)
    figs["mppi_loop"] = fig_mppi_loop(out)
    figs["v1_results"] = fig_v1_results(out)
    if dense_path and os.path.exists(dense_path):
        figs["dense"], info["dense"] = fig_dense(pickle.load(open(dense_path, "rb")), out)
    figs["compute"] = fig_compute(out)
    figs["architecture"] = fig_architecture(out)
    figs["features"], info["features"] = fig_features(eps, out)
    figs["masks"] = fig_masks(out)
    figs["action_tokens"] = fig_action_tokens(out)
    figs["parallel"] = fig_parallel(out)
    figs["sigreg"], info["sigreg"] = fig_sigreg(out)
    figs["schedule"] = fig_schedule(out)
    figs["masked_windows"] = fig_masked_windows(eps, v1_eps, out)
    figs["prober"] = fig_prober(out)
    figs["ensemble"] = fig_ensemble(out)
    figs["params"], info["params"] = fig_params(out)
    return {"figs": figs, "info": info}


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--wide", required=True)
    ap.add_argument("--dense", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    r = make_all(a.wide, a.dense, a.out)
    print(r["info"])
