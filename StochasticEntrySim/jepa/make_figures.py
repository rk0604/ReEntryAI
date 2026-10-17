"""
Paper-style schematic figures for the JEPA world-model experiment.

    python jepa/make_figures.py        (from StochasticEntrySim/)

Writes PNGs (300 dpi) to jepa/figures/:
    fig1_pipeline.png           the end-to-end experiment
    fig2_jepa_architecture.png  stage 1 (JEPA training) and stage 2 (PI prober)
    fig3_mppi_loop.png          planning with the world model in closed loop
"""

import os

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")

# Categorical slots (fixed meaning across all figures): edge color, light fill.
BLUE = ("#2a78d6", "#e3edfa")     # state encoder / latents
AQUA = ("#1baf7a", "#dcf3ea")     # action encoder
ORANGE = ("#eb6834", "#fde6dc")   # predictor / planner
VIOLET = ("#4a3aa7", "#e6e3f5")   # physics / prober
YELLOW = ("#eda100", "#fdf1d3")   # losses / cost
GRAY = ("#52514e", "#f2f1ed")     # data / simulator
INK, MUTED = "#1b1b1a", "#7a7974"

plt.rcParams.update({"font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
                     "font.size": 9, "savefig.dpi": 300})


def canvas(w, h, xmax, ymax):
    """Blank figure with a data coordinate system [0, xmax] x [0, ymax]."""
    fig, ax = plt.subplots(figsize=(w, h))
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.axis("off")
    return fig, ax


def box(ax, x, y, w, h, text, col=GRAY, fs=8.5, bold=False, dashed=False, tc=INK, r=0.6):
    """Rounded box centered at (x, y). Returns its anchor points."""
    p = FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                       fc=col[1], ec=col[0], lw=1.2, ls=(0, (3, 2)) if dashed else "-", zorder=2)
    ax.add_patch(p)
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=tc,
            fontweight="bold" if bold else "normal", zorder=3, linespacing=1.35)
    return {"l": (x - w / 2, y), "r": (x + w / 2, y), "t": (x, y + h / 2), "b": (x, y - h / 2), "c": (x, y)}


def arrow(ax, a, b, col=INK, lw=1.1, ls="-", rad=0.0, head=8):
    """Arrow from point a to point b, optionally curved."""
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=head, lw=lw, color=col,
                                 ls=ls, connectionstyle=f"arc3,rad={rad}", zorder=1,
                                 shrinkA=1, shrinkB=1))


def panel(ax, x0, y0, x1, y1, title, col=GRAY):
    """Light background panel with a bold title, to group a stage."""
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0,rounding_size=1.0",
                                fc="#fbfbfa", ec="#d9d8d3", lw=0.9, zorder=0))
    ax.text(x0 + 1.0, y1 - 1.2, title, ha="left", va="top", fontsize=9.5, fontweight="bold", color=col[0])


def save(fig, name):
    """Write a figure to jepa/figures/."""
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    fig.savefig(path, bbox_inches="tight", facecolor="white", pad_inches=0.05)
    plt.close(fig)
    print("wrote", path)


# ============================================================================
# Figure 1: the experiment pipeline
# ============================================================================
def fig1():
    """Data synthesis -> stage 1 -> stage 2 -> closed-loop MPPI -> evaluation."""
    fig, ax = canvas(11.5, 3.9, 115, 39)
    fs = 7.6
    cols = [(14, "(a) Data synthesis", GRAY), (43, "(b) Stage 1: JEPA", BLUE),
            (72, "(c) Stage 2: PI prober", VIOLET), (101, "(d) Closed loop + eval", ORANGE)]
    rows = [
        [("Random bank profiles\nGP + reversal-structured", GRAY, False),
         ("Real closed loop (OrionEntryEnv)\nactuator → roll PD → RCS → 3DOF", GRAY, False),
         ("Domain randomization\nρ, L/D ×[0.85, 1.15],  m ±5%", GRAY, False)],
        [("320 episodes → 5 s steps\n35 s history windows", GRAY, False),
         ("Enc$_\\theta$, Enc$_\\phi$ (causal TCN)\n+ GRU predictor (26k params)", BLUE, False),
         ("$\\mathcal{L}=\\sum_k\\|\\tilde s_{t+k}-s_{t+k}\\|^2+\\lambda\\,$SIGReg\n(no decoder)", YELLOW, False)],
        [("Frozen latent rollout\n$\\tilde s_{t+1},\\ldots,\\tilde s_{t+T}$", BLUE, True),
         ("$\\psi(\\tilde s)\\rightarrow(c_D,\\,c_L,\\,\\Delta\\sigma)$\ncorrections to nominal aero", VIOLET, False),
         ("Exact 3DOF EOM\n→ $h,\\phi,\\lambda,V,\\gamma,\\chi$", VIOLET, False)],
        [("MPPI every 5 s\n256 bank plans → drogue", ORANGE, False),
         ("Frozen test set, 2 truths\nnominal / ρ, L/D ×0.88", GRAY, False),
         ("vs. classical pred.-corrector\nand predictive-model MPPI", YELLOW, False)],
    ]
    for (x, title, col), items in zip(cols, rows):
        panel(ax, x - 13, 3, x + 13, 37, title, col)
        bs = [box(ax, x, y, 23, 6.5, t, c, fs=fs, dashed=d) for y, (t, c, d) in zip((27.5, 18, 8.5), items)]
        arrow(ax, bs[0]["b"], bs[1]["t"])
        if x == 14:
            arrow(ax, bs[2]["t"], bs[1]["b"])
        else:
            arrow(ax, bs[1]["b"], bs[2]["t"])
    for x in (27, 56, 85):
        arrow(ax, (x - 0.2, 18), (x + 3.2, 18), lw=1.6, head=11)
    save(fig, "fig1_pipeline.png")


# ============================================================================
# Figure 2: JEPA architecture (stage 1) and the PI prober (stage 2)
# ============================================================================
def fig2():
    """Encoders, latent predictor unroll, SIGReg, and the physics prober."""
    fig, ax = canvas(12, 7.2, 120, 72)

    # ---------------- Stage 1 ----------------
    panel(ax, 1, 33, 119, 71, "Stage 1 — self-supervised latent dynamics (no reconstruction)", BLUE)
    xs = box(ax, 11, 56, 17, 7, "State history $X_t$\n7 × 10  (35 s)", GRAY, fs=8)
    xa = box(ax, 11, 44.5, 17, 6.5, "Bank commands\n$\\sigma_{t-6},\\ldots,\\sigma_{t+T-1}$", GRAY, fs=8)
    es = box(ax, 30, 56, 13, 7, "Enc$_\\theta$\ncausal TCN", BLUE, bold=True)
    ea = box(ax, 30, 44.5, 13, 6.5, "Enc$_\\phi$\ncausal TCN", AQUA, bold=True)
    arrow(ax, xs["r"], es["l"])
    arrow(ax, xa["r"], ea["l"])
    s0 = box(ax, 44, 56, 7, 5.5, "$s_t$", BLUE, fs=10)
    arrow(ax, es["r"], s0["l"])

    xpos = [56, 71, 86, 104]
    zl = ["$z_t$", "$z_{t+1}$", None, "$z_{t+T-1}$"]
    sl = ["$\\tilde s_{t+1}$", "$\\tilde s_{t+2}$", None, "$\\tilde s_{t+T}$"]
    preds = []
    for x, z, st in zip(xpos, zl, sl):
        if z is None:
            ax.text(x, 49, "$\\cdots$", ha="center", va="center", fontsize=14, color=ORANGE[0])
            preds.append(None)
            continue
        p = box(ax, x, 49, 9, 6, "GRU\nPred", ORANGE, fs=8, bold=True)
        preds.append(p)
        ax.text(x, 41.2, z, ha="center", va="center", fontsize=9, color=AQUA[0])
        arrow(ax, (x, 42.6), p["b"], col=AQUA[0])
        b = box(ax, x, 58, 8, 4.6, st, BLUE, fs=9.5)
        arrow(ax, p["t"], b["b"])
        ax.annotate("", xy=(x, 60.5), xytext=(x, 64.3),
                    arrowprops=dict(arrowstyle="<->", color=YELLOW[0], lw=1.1))
    arrow(ax, ea["r"], (50, 42.5), col=AQUA[0])
    arrow(ax, s0["b"], preds[0]["l"], rad=0.25)
    arrow(ax, preds[0]["r"], preds[1]["l"])
    arrow(ax, preds[1]["r"], (81, 49))
    arrow(ax, (91, 49), preds[3]["l"])
    box(ax, 93, 66.6, 46, 4.4, "targets $s_{t+k}=\\mathrm{Enc}_\\theta(X_{t+k})$ — same encoder, gradients flow",
        BLUE, fs=7.8, dashed=True)
    ax.text(63.5, 62.4, "$\\|\\tilde s-s\\|^2$", ha="center", va="center", fontsize=8.5, color=YELLOW[0])
    box(ax, 22, 36.3, 38, 4.6,
        "SIGReg on all $s,\\tilde s$: random projections $\\xi_m$,\neach 1-D marginal → $\\mathcal{N}(0,1)$ (anti-collapse)",
        YELLOW, fs=7.3)
    ax.text(80, 36.2, "$\\mathcal{L}_{1}=\\frac{1}{T}\\sum_{k=1}^{T}\\|\\tilde s_{t+k}-s_{t+k}\\|^2+\\lambda\\,\\mathrm{SIGReg}$,"
            "   $\\lambda=0.02$,  $T=24$ (120 s)", ha="center", va="center", fontsize=9, color=INK)

    # ---------------- Stage 2 ----------------
    panel(ax, 1, 1, 119, 31, "Stage 2 — physics-inspired prober on frozen latents (stop-gradient)", VIOLET)
    y = 19
    fl = box(ax, 12, y, 18, 7, "$\\tilde s_{t+k+1}$\n(frozen JEPA rollout)", BLUE, dashed=True, fs=8.2)
    mlp = box(ax, 31, y, 11, 7, "$\\psi$  (MLP)", VIOLET, bold=True)
    cr = box(ax, 50, y, 19, 10, "$c_D$  log drag factor\n$c_L$  log L/D factor\n$\\Delta\\sigma$  bank offset (lag)",
             VIOLET, fs=8)
    eom = box(ax, 80, y, 22, 10, "3DOF EOM\n$\\dot V=-a_D-g\\sin\\gamma$\n$\\dot\\gamma=a_{Lv}/V+\\ldots$\n"
              "$\\dot\\chi=-a_{Lh}/(V\\cos\\gamma)+\\ldots$", VIOLET, fs=7.6)
    out = box(ax, 107, y, 16, 8, "$\\hat x_{t+k+1}$\n$h,\\phi,\\lambda,V,\\gamma,\\chi$", GRAY, fs=8.5)
    arrow(ax, fl["r"], mlp["l"])
    arrow(ax, mlp["r"], cr["l"])
    arrow(ax, cr["r"], eom["l"])
    arrow(ax, eom["r"], out["l"])
    nom = box(ax, 60, 6.8, 44, 6.0,
              "nominal aero:  $a_D=\\frac{1}{2}\\rho_0 e^{-h/H}V^2/\\beta\\cdot e^{c_D}$,  $a_L=a_D\\,0.3\\,e^{c_L}$\n"
              "lift along $\\sigma_{cmd}+\\Delta\\sigma$, evaluated at the integrated $h, V$", GRAY, fs=7.3)
    arrow(ax, cr["b"], (cr["b"][0], 9.8), col=VIOLET[0])
    arrow(ax, (eom["b"][0] - 4, 9.8), (eom["b"][0] - 4, eom["b"][1]), col=MUTED)
    ax.text(18, 6.8, "$\\sigma_{cmd}$", fontsize=10, color=AQUA[0], ha="center", va="center")
    arrow(ax, (21, 6.8), nom["l"], col=AQUA[0])
    arrow(ax, out["b"], (nom["r"][0], 6.0), col=MUTED, rad=-0.3, ls=(0, (3, 2)))
    ax.text(110.5, 4.2, "integrated $h, V$ feed\nthe next step's forces", fontsize=7, color=MUTED, ha="center")
    ax.text(98, 27.5, "$\\mathcal{L}_2=\\sum_k\\|\\hat x_{t+k}-x_{t+k}\\|^2$ + force match", fontsize=7.8,
            color=YELLOW[0], ha="center", va="center")
    save(fig, "fig2_jepa_architecture.png")


# ============================================================================
# Figure 3: MPPI with the world model in the loop
# ============================================================================
def fig3():
    """Receding-horizon MPPI on the learned model, in the real closed loop."""
    fig, ax = canvas(11.5, 5.6, 115, 56)
    ax.set_xlim(-3, 115)

    panel(ax, 1, 12, 82, 55, "Planner (every 5 s, one batched call)", ORANGE)
    hist = box(ax, 13, 44, 18, 7.5, "history: 7 states\n+ 6 past commands", GRAY, fs=8)
    samp = box(ax, 13, 27, 18, 9, "sample $S{=}256$ bank\nsequences: warm start\n+ noise on 6 knots", ORANGE, fs=8)
    wm = box(ax, 37, 35.5, 18, 15, "JEPA world model\n+ PI prober\n\nroll all $S$ samples\nto drogue altitude\n"
             "(≤110 steps, 550 s)", BLUE, fs=8)
    arrow(ax, hist["r"], (wm["l"][0], 40))
    arrow(ax, samp["r"], (wm["l"][0], 31))
    cost = box(ax, 64, 41, 30, 13, "cost per sample  $J^{(s)}$ =\nmiss at drogue crossing [km]\n"
               "+ smoothness  + g over 8 g\n+ $|\\gamma|$ over 82° before drogue\n+ penalty if never reaches drogue",
               YELLOW, fs=7.8)
    arrow(ax, (wm["r"][0], 40), (cost["l"][0], 41))
    wts = box(ax, 64, 22, 30, 8.5, "$w^{(s)}\\propto\\exp(-(J^{(s)}-J_{min})/\\lambda)$\n"
              "$u\\leftarrow\\sum_s w^{(s)}u^{(s)}$,  shift warm start", ORANGE, fs=8)
    arrow(ax, cost["b"], wts["t"])
    arrow(ax, (wts["l"][0], 19.5), samp["b"], rad=-0.2, col=MUTED, ls=(0, (3, 2)))
    ax.text(33, 15.0, "warm start for next decision", fontsize=7.2, color=MUTED, ha="center")

    panel(ax, 85, 2, 114, 55, "Real simulator", GRAY)
    cmd = box(ax, 99.5, 41, 25, 7, "first command $\\sigma_{cmd}$\nheld for 5 s", ORANGE, fs=8)
    env = box(ax, 99.5, 27, 25, 11, "OrionEntryEnv\nbank actuator (20°/s)\n→ roll PD → 12-jet RCS\n"
              "→ 3DOF + NRLMSISE", GRAY, fs=7.8)
    st = box(ax, 99.5, 8, 25, 7, "new state + sensed\ndrag/lift (IMU)", GRAY, fs=8)
    arrow(ax, (wts["r"][0], 22), (cmd["l"][0], 39), rad=0.2, lw=1.5, head=10)
    arrow(ax, cmd["b"], env["t"])
    arrow(ax, env["b"], st["t"])
    # Feedback path routed under the planner panel and up its left edge.
    ax.plot([st["l"][0], -0.5, -0.5], [8, 8, 44], color=INK, lw=1.5, zorder=1)
    arrow(ax, (-0.5, 44), hist["l"], lw=1.5, head=10)
    ax.text(43, 9.6, "feedback: replan from the new history every 5 s (no policy, no sign logic)",
            fontsize=7.8, color=INK, ha="center")
    save(fig, "fig3_mppi_loop.png")


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
