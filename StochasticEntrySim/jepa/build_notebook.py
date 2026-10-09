"""
Build ReEntryAI_JEPA.ipynb, the end-to-end JEPA world-model experiment.

All logic lives in the jepa/ package; the notebook only configures, calls, and
plots. Edit this builder and re-run it rather than hand-editing the .ipynb:

    python jepa/build_notebook.py        (from StochasticEntrySim/)
"""
import os
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "ReEntryAI_JEPA.ipynb")

nb = nbf.v4.new_notebook()
cells = []


def md(text):
    """Append a markdown cell."""
    cells.append(nbf.v4.new_markdown_cell(text.strip("\n")))


def code(text):
    """Append a code cell."""
    cells.append(nbf.v4.new_code_cell(text.strip("\n")))


# ---------------------------------------------------------------------------
md(r"""
# ReEntryAI × SkyJEPA — a latent world model + MPPI for entry guidance

**Question.** Can a JEPA-style latent world model, trained only on randomized
closed-loop sim data, steer Orion to the drogue-deploy target by planning
(MPPI) rather than by a learned policy?

**Why this, after the RL results.** PPO from scratch never found the reversal
schedule (~1% success), and the behavior clone stalls at a ~91 km median miss
because of compounding errors. A planner on a learned *model* avoids both: there
is no exploration problem (the model is learned from any data, not from reward),
and feedback comes from re-planning every 5 s, not from a policy that drifts
off the expert's state distribution. The bar is the classical
predictor-corrector: **2.05 km median miss, 99.9% success** on the frozen test set.

**Recipe** (Rao et al., *SkyJEPA*, arXiv 2606.23444, mapped to this sim — full
design and rationale in [`jepa/PLAN.md`](jepa/PLAN.md)):

| SkyJEPA (quadrotor) | Here (Orion entry) |
|---|---|
| GP position references tracked by NMPC/MPPI | GP and reversal-structured **bank profiles** flown through the real closed loop (actuator → roll PD → RCS → 3DOF) |
| 500 randomized airframes (mass, inertia, drag, motors) | per-episode **density, L/D, mass** randomization + entry-interface dispersion |
| TCN encoders + GRU predictor, multi-step latent loss + SIGReg | same, ~26k params, 5 s macro steps, 35 s history |
| physics-inspired prober: nominal thrust + learned residual, SO(3) integrator | **nominal aero (exponential atmosphere, one β, lift along bank) + 3 latent corrections**, exact 3DOF EOM |
| MPPI tracking cost, 512 samples × 15 steps | MPPI on **predicted drogue-point miss**, 256 samples × up to 110 steps (to the ground) |

**Run modes.** `MODE="quick"` is a ~5 min smoke test of every stage (numbers
are meaningless). `MODE="full"` is the real experiment, ~20–25 min on a
12-core laptop CPU, no GPU. Collected data, checkpoints, and classical-baseline
results are cached under `jepa_runs/`, so re-runs skip what is already done.
""")

# ---------------------------------------------------------------------------
md("## 0. Setup")
code(r"""
import os, sys, math, time, pickle
from pathlib import Path

SIM_DIR = Path.cwd()                     # run from StochasticEntrySim/
assert (SIM_DIR / "rl_env.py").exists(), "start Jupyter in StochasticEntrySim/"
sys.path.insert(0, str(SIM_DIR))

import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt
%matplotlib inline

from jepa import data as D, models as M, evaluate as E, mppi as P

MODE = "quick"        # "quick" (smoke test, ~5 min) | "full" (the experiment, ~20-25 min)
CFG = {
    "quick": dict(n_episodes=24, horizon=12, prober_horizon=12, batch=256,
                  steps_jepa=300, steps_pi=200, steps_mlp=150, steps_pred=300,
                  k_eval=30, ttg_s=(300, 200, 100), n_cases=2, mppi_samples=128),
    "full":  dict(n_episodes=320, horizon=24, prober_horizon=16, batch=256,
                  steps_jepa=3000, steps_pi=1500, steps_mlp=800, steps_pred=2000,
                  k_eval=60, ttg_s=(400, 300, 200, 100, 50), n_cases=16, mppi_samples=256),
}[MODE]
MAX_STEPS = 900          # cap on 1 Hz steps per data episode (skip-outs time out here)
MACRO_DT = D.MACRO_DT_S  # 5 s: one bank command held per world-model step
N_SUB = 2                # Euler substeps per macro step inside the prober
LAMBDA_SIG = 0.02        # SIGReg weight (paper value)
H_DROGUE = 7600.0        # end of the controlled segment (configs/default.json)
WORKERS = max(1, (os.cpu_count() or 2) - 1)
# One thread is measurably fastest for these ~25k-parameter models on CPU:
# the step time is dominated by per-op dispatch, which threads only add to.
torch.set_num_threads(1)
torch.manual_seed(0)

RUN_DIR = SIM_DIR / "jepa_runs" / MODE
RUN_DIR.mkdir(parents=True, exist_ok=True)
CKPT = RUN_DIR / "jepa_ckpt.pt"
print(f"MODE={MODE}  workers={WORKERS}  run dir={RUN_DIR}")
print(CFG)
""")

code(r"""
# Plot style: one axis per chart, recessive grid, fixed color per model/controller
# (categorical slots in fixed order; the color follows the entity everywhere).
plt.rcParams.update({
    "figure.dpi": 110, "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6,
    "axes.edgecolor": "#8a8984", "axes.labelcolor": "#52514e", "axes.titlecolor": "#0b0b0b",
    "xtick.color": "#52514e", "ytick.color": "#52514e", "lines.linewidth": 1.6,
    "legend.frameon": False,
})
COL = {"JEPA+PI": "#2a78d6", "Predictive": "#eb6834", "JEPA+MLP": "#1baf7a",
       "Physics prior only": "#eda100", "physics-nominal": "#4a3aa7",
       "jepa": "#2a78d6", "pred": "#eb6834", "classical": "#4a3aa7"}
MUTED = "#a3a29c"
""")

# ---------------------------------------------------------------------------
md(r"""
## 1. Data synthesis (paper §IV-E)

The paper's three data requirements, and how each is met:

* **Diverse** — half the episodes fly a squashed multi-length-scale GP bank
  profile (20 s / 50 s / 150 s length scales: broad coverage of magnitude and
  sign), half fly a *reversal-structured* profile (slowly varying magnitude, sign
  flips every ~60–200 s), the shape real entry guidance produces. Commands are
  held per 5 s macro step, exactly as the planner will issue them.
* **Dynamically feasible** — every profile is flown through `OrionEntryEnv` in
  policy mode, so the bank actuator's 20 °/s rate limit and the roll/RCS loop
  are in the data. The world model must learn that lag.
* **Robust** — each episode draws a truth with density ×U(0.85, 1.15), L/D
  ×U(0.85, 1.15), mass ×U(0.95, 1.05), on top of the config's IC dispersion.
  (Implemented through `ref_area_m2`, `cpas_lift_scale`, `mass_kg` in the env
  params: no physics module is touched.)

Episodes are random, so most land hundreds of km from the target. That is
intended: the world model is task-agnostic; the target only enters in the MPPI cost.
""")
code(r"""
t0 = time.time()
episodes = D.collect_dataset(CFG["n_episodes"], workers=WORKERS, max_steps=MAX_STEPS,
                             cache_path=str(SIM_DIR / "jepa_runs" / f"episodes_{MAX_STEPS}.pkl"))
print(f"{len(episodes)} episodes in {time.time() - t0:.0f}s")
ep_df = pd.DataFrame([{"kind": e["kind"], "outcome": e["outcome"], "flight_s": e["t_final_s"],
                       "miss_km": e["miss_km"], **e["domain"]} for e in episodes])
display(ep_df.groupby(["kind", "outcome"]).agg(n=("flight_s", "size"), flight_s=("flight_s", "median"),
                                               miss_km=("miss_km", "median")).round(1))
display(ep_df[["density_scale", "lift_scale", "mass_scale"]].describe().loc[["min", "mean", "max"]].round(3))
""")
code(r"""
fig, axes = plt.subplots(3, 1, figsize=(8, 5.4), sharex=True)
for ax, e in zip(axes, episodes[:3]):
    t = np.arange(len(e["cmd"]))
    ax.plot(t, np.degrees(e["cmd"]), color=MUTED, lw=1.2, label="commanded")
    ax.plot(t, np.degrees(e["y"][1:, 6]), color=COL["JEPA+PI"], label="realized (actuator lag)")
    ax.set_ylabel("bank [deg]")
    ax.set_title(f'{e["kind"]} profile, outcome={e["outcome"]}, '
                 f'density x{e["domain"]["density_scale"]:.2f}, L/D x{e["domain"]["lift_scale"]:.2f}',
                 loc="left", fontsize=9)
axes[0].legend(loc="upper right", ncol=2)
axes[-1].set_xlabel("time since entry interface [s]")
fig.tight_layout()

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
Y = np.concatenate([e["y"][::5] for e in episodes])
dens = np.concatenate([np.full(len(e["y"][::5]), e["domain"]["density_scale"]) for e in episodes])
sc = axes[0].scatter(Y[:, 3] / 1000, Y[:, 0] / 1000, c=dens, cmap="Blues", s=3, vmin=0.8, vmax=1.15)
fig.colorbar(sc, ax=axes[0], label="truth density scale")
axes[0].set(xlabel="speed [km/s]", ylabel="altitude [km]", title="State coverage (every 5 s)")
for e in episodes[:60]:
    axes[1].plot(np.degrees(e["y"][:, 2]), np.degrees(e["y"][:, 1]), color=MUTED, lw=0.8)
axes[1].plot([15], [15], marker="*", ms=14, color="#e34948", ls="none", label="target")
axes[1].set(xlabel="longitude [deg]", ylabel="latitude [deg]", title="Ground tracks (random bank profiles)")
axes[1].legend(loc="lower right")
fig.tight_layout()
""")

# ---------------------------------------------------------------------------
md(r"""
## 2. Macro-step dataset

Each 1 Hz episode is subsampled to 5 s macro steps. A training window holds the
history the encoders need (receptive field rf = 7 macro steps = 35 s, front-padded
with the entry state and wings-level commands so windows can start at entry
interface, as the planner must), then *T* future steps. Splits are by episode.
""")
code(r"""
macro = [D.to_macro(e, MACRO_DT) for e in episodes]
train_eps, val_eps, test_eps = D.split_episodes(macro, (0.8, 0.1, 0.1), seed=0)
stats = D.compute_stats(train_eps)
RF = M.JEPAWorldModel().rf
ds_tr = D.MacroDataset(train_eps, RF, CFG["horizon"], stats)
ds_va = D.MacroDataset(val_eps, RF, CFG["horizon"], stats)
ds_tr_p = D.MacroDataset(train_eps, RF, CFG["prober_horizon"], stats)
ds_va_p = D.MacroDataset(val_eps, RF, CFG["prober_horizon"], stats)
print(f"episodes train/val/test = {len(train_eps)}/{len(val_eps)}/{len(test_eps)}")
print(f"windows (T={CFG['horizon']}) train/val = {len(ds_tr)}/{len(ds_va)};  rf = {RF} macro steps "
      f"= {RF * MACRO_DT:.0f} s of history")
""")

# ---------------------------------------------------------------------------
md(r"""
## 3. Stage 1 — JEPA latent dynamics (paper §IV-A, IV-B)

* `Enc_θ`: causal TCN over the last 7 macro states (altitude, speed, γ, heading,
  latitude, realized bank, roll rate, sensed drag and lift) → latent *s* ∈ ℝ³².
  Longitude is left out: the dynamics don't depend on it.
* `Enc_φ`: causal TCN over the last 7 bank commands → action code *z* ∈ ℝ¹⁶
  (this is where actuator lag can be represented).
* `Pred`: GRU unrolled *T* steps on the encoded actions; trained so the rollout
  matches the *encoded* future, `‖s̃ₜ₊ₖ − sₜ₊ₖ‖²`, plus `λ·SIGReg` to keep the
  latents an isotropic Gaussian. There is no decoder, no EMA target, no stop-gradient.

What to watch: `pred` (train) and `val` should fall well below ~2, which is what
two independent N(0, I) latents would score. `sigreg` should settle near a
small constant. A collapsed encoder shows up as a tiny `pred` together with a large
`sigreg` and a low effective rank below.
""")
code(r"""
wm = M.JEPAWorldModel()
sigreg = M.SIGReg()
print(f"world model: {M.n_params(wm):,} parameters  (paper: ~9k)")
g = torch.Generator().manual_seed(1)
vb = ds_va.sample(512, g)

def jepa_val():
    with torch.no_grad():
        return M.jepa_loss(wm, sigreg, vb, 0.0)[1]["pred"]

t0 = time.time()
hist_jepa = M.fit(wm.parameters(), lambda: M.jepa_loss(wm, sigreg, ds_tr.sample(CFG["batch"]), LAMBDA_SIG),
                  CFG["steps_jepa"], val_fn=jepa_val, log_every=max(50, CFG["steps_jepa"] // 10), label="jepa")
wm.eval()
print(f"stage 1: {time.time() - t0:.0f}s")
""")
code(r"""
fig, axes = plt.subplots(1, 3, figsize=(12, 3.3))
axes[0].plot(hist_jepa["step"], hist_jepa["pred"], color=COL["JEPA+PI"], label="train")
axes[0].plot(hist_jepa["step"], hist_jepa["val"], color=COL["JEPA+PI"], ls="--", label="val")
axes[0].set(xlabel="step", ylabel="latent MSE", title="Multi-step latent prediction", yscale="log")
axes[0].legend()
axes[1].plot(hist_jepa["step"], hist_jepa["sigreg"], color=COL["JEPA+PI"])
axes[1].set(xlabel="step", ylabel="Epps-Pulley statistic", title="SIGReg", yscale="log")

health = E.latent_health(wm, ds_va)
axes[2].bar(np.arange(len(health["std"])), health["std"], color=COL["JEPA+PI"], width=0.7)
axes[2].axhline(1.0, color=MUTED, lw=1)
axes[2].set(xlabel="latent dim", ylabel="std", title=f'Latent std (eff. rank {health["eff_rank"]:.1f} / {wm.d_latent})')
fig.tight_layout()
print(f'effective rank {health["eff_rank"]:.1f} of {wm.d_latent}, mean |corr| {health["mean_abs_corr"]:.2f}')

# Does the latent organize by flight regime? PCA of encoded states, colored by altitude.
b = ds_va.batch(torch.arange(min(len(ds_va), 3000)))
with torch.no_grad():
    s = wm.encode_states(b["feats"])[:, 0]
U, S_, V_ = torch.pca_lowrank(s, q=2)
pc = (s - s.mean(0)) @ V_
fig, ax = plt.subplots(figsize=(4.6, 3.6))
sc = ax.scatter(pc[:, 0], pc[:, 1], c=b["y"][:, 0, 0] / 1000, cmap="Blues", s=4)
fig.colorbar(sc, ax=ax, label="altitude [km]")
ax.set(xlabel="PC1", ylabel="PC2", title="Encoded latents (val)")
fig.tight_layout()
""")

# ---------------------------------------------------------------------------
md(r"""
## 4. Stage 2 — probers on frozen latents (paper §IV-C)

The world model is frozen; only the prober trains, on the rollout loss
`Σₖ ‖x̃ₜ₊ₖ − xₜ₊ₖ‖²` (altitude, position, speed, γ, heading, made dimensionless).

* **PI prober** (the method). The latent outputs three corrections: a log drag
  factor `c_D`, a log L/D factor `c_L`, and a bank offset. These correct a deliberately crude
  nominal model: drag = ½ρ_exp(h)V²/β_nom·e^{c_D} with a single exponential
  atmosphere, lift = drag·0.30·e^{c_L}, lift directed along the *commanded* bank
  plus the offset. The forces are evaluated at the **integrated** altitude and speed and
  pushed through the exact 3DOF EOM. The latent has to supply the real
  density profile, the Mach-scheduled CD and L/D, the actuator lag, and the
  randomized truth (inferred from the 35 s history). An auxiliary loss matches the
  forces to the sensed ones, so the g-load the planner reads off is physical.
* **MLP prober** (ablation, the paper's "Ours + Prober"). The latent maps to a free-form
  state increment per step, with no physics.
* **Physics prior only** (free ablation). The PI prober with all corrections at
  zero, i.e. the crude nominal model with instantaneous bank and no learning.
""")
code(r"""
gp = torch.Generator().manual_seed(2)
vbp = ds_va_p.sample(512, gp)
pi = M.PIProber(wm.d_latent)
mlp = M.MLPProber(wm.d_latent, stats.dy_std[:6])
prior_only = M.PIProber(wm.d_latent)          # untrained: zero corrections

def prober_val(p):
    with torch.no_grad():
        return M.prober_loss(wm, p, vbp, MACRO_DT, N_SUB)[1]["state"]

t0 = time.time()
hist_pi = M.fit(pi.parameters(), lambda: M.prober_loss(wm, pi, ds_tr_p.sample(CFG["batch"]), MACRO_DT, N_SUB, w_acc=1.0),
                CFG["steps_pi"], val_fn=lambda: prober_val(pi), log_every=max(50, CFG["steps_pi"] // 8), label="PI prober")
hist_mlp = M.fit(mlp.parameters(), lambda: M.prober_loss(wm, mlp, ds_tr_p.sample(CFG["batch"]), MACRO_DT, N_SUB),
                 CFG["steps_mlp"], val_fn=lambda: prober_val(mlp), log_every=max(50, CFG["steps_mlp"] // 4), label="MLP prober")
pi.eval(); mlp.eval()
print(f"stage 2: {time.time() - t0:.0f}s   |  val state loss: physics prior only = {prober_val(prior_only):.3g}")
""")

# ---------------------------------------------------------------------------
md(r"""
## 5. Baseline — autoregressive predictive model (the paper's "Predictive")

Same encoders and a similar size, but it regresses the next macro step's physical
state directly (one-step, teacher-forced, paper Eq. 4) and is rolled out on its
own predictions (Eq. 5–6). This is the model class SkyJEPA argues compounds error.
""")
code(r"""
pm = M.PredictiveModel(stats.dy_std)
print(f"predictive model: {M.n_params(pm):,} parameters")
vb1 = ds_va.sample(512, torch.Generator().manual_seed(3))
t0 = time.time()
hist_pred = M.fit(pm.parameters(), lambda: M.predictive_loss(pm, ds_tr.sample(CFG["batch"])), CFG["steps_pred"],
                  val_fn=lambda: M.predictive_loss(pm, vb1)[1]["onestep"], log_every=max(50, CFG["steps_pred"] // 5),
                  label="predictive")
pm.eval()
print(f"predictive: {time.time() - t0:.0f}s")
M.save_checkpoint(str(CKPT), wm, pi, mlp, pm, stats, extra={"mode": MODE, "cfg": CFG, "macro_dt": MACRO_DT, "n_sub": N_SUB})
print(f"saved -> {CKPT}")
""")

# ---------------------------------------------------------------------------
md(r"""
## 6. Open-loop evaluation (paper §VI-A, VI-B)

Held-out episodes, with the true future commands fed in. Two views:

1. **Rollout error vs. unroll length** (paper Fig. 8 / Table III), out to
   `k_eval` steps. That is 2.5× the JEPA training horizon in full mode, which
   tests the paper's long-horizon claim.
2. **Footprint error vs. time-to-go**, the number guidance actually depends on: from
   the state at time-to-go τ, where does each model say the drogue deploy will
   be, and how far is that from where it really was? `physics-nominal` is the
   classical predictor's view, with the full NRLMSISE-table + Orion aero schedule at
   *nominal* parameters and instantaneous bank. It is a strong reference: it
   only errs through the randomized truth and the actuator lag, and those are
   exactly the two things the learned model is supposed to pick up from history.
""")
code(r"""
models = {"JEPA+PI": {"kind": "jepa", "wm": wm, "prober": pi},
          "JEPA+MLP": {"kind": "jepa", "wm": wm, "prober": mlp},
          "Physics prior only": {"kind": "jepa", "wm": wm, "prober": prior_only},
          "Predictive": {"kind": "pred", "pm": pm}}
K = CFG["k_eval"]
curves = E.rollout_errors(models, test_eps + val_eps, stats, RF, K=K, n_windows=512, dt=MACRO_DT, n_sub=N_SUB)
kk = np.arange(1, K + 1) * MACRO_DT
fig, axes = plt.subplots(1, 4, figsize=(14, 3.2))
for ax, (q, lab) in zip(axes, [("pos_km", "position RMSE [km]"), ("h_km", "altitude RMSE [km]"),
                               ("V_mps", "speed RMSE [m/s]"), ("gamma_deg", "γ RMSE [deg]")]):
    for name, c in curves.items():
        ax.plot(kk, c[q], color=COL[name], label=name)
    ax.axvline(CFG["horizon"] * MACRO_DT, color=MUTED, lw=1, ls=":")
    ax.set(xlabel="prediction horizon [s]", ylabel=lab, yscale="log")
axes[0].legend(fontsize=8)
axes[0].set_title("dotted = JEPA training horizon", loc="left", fontsize=8)
fig.tight_layout()
display(pd.DataFrame({n: {f"pos_km@{int(kk[i])}s": c["pos_km"][i] for i in (min(11, K - 1), K // 2, K - 1)}
                      for n, c in curves.items()}).T.round(2))
""")
code(r"""
nominal = E.nominal_params()
t0 = time.time()
fp = pd.DataFrame(E.footprint_errors(models, test_eps + val_eps, stats, RF, ttg_s=CFG["ttg_s"],
                                     h_cross=H_DROGUE, dt=MACRO_DT, n_sub=N_SUB, nominal=nominal))
print(f"footprint eval: {time.time() - t0:.0f}s on {fp[fp.model == 'JEPA+PI'].shape[0] // len(CFG['ttg_s'])} episodes")
med = fp.groupby(["model", "ttg_s"]).err_km.median().unstack(0)
fig, ax = plt.subplots(figsize=(6.2, 3.6))
for name in med.columns:
    ax.plot(med.index, med[name], marker="o", ms=5, color=COL[name], label=name)
ax.set(xlabel="time-to-go when predicting [s]", ylabel="median footprint error [km]", yscale="log",
       title="Predicted vs. true drogue-deploy point")
ax.invert_xaxis()
ax.legend(fontsize=8)
fig.tight_layout()
display(med.round(2))

# Does history conditioning buy robustness to the randomized truth? Split by how far the
# truth was from nominal (|ln density scale| + |ln L/D scale|).
fp["mismatch"] = np.abs(np.log(fp.dom_density_scale)) + np.abs(np.log(fp.dom_lift_scale))
fp["truth"] = np.where(fp.mismatch > fp.mismatch.median(), "far from nominal", "near nominal")
display(fp.pivot_table(index="model", columns=["ttg_s", "truth"], values="err_km", aggfunc="median").round(1))
""")
code(r"""
# Temporal straightening (paper Sec. VI-A): mean cosine between consecutive latent steps.
dsK = D.MacroDataset(test_eps + val_eps, RF, K, stats)
bK = dsK.batch(torch.randperm(len(dsK))[:256])
P_ = RF - 1
with torch.no_grad():
    lat_pred = M.jepa_predict(wm, pi, bK["feats"][:, :P_ + 1], bK["y"][:, 0, :6], bK["u"][:, :P_],
                              bK["u"][:, P_:], MACRO_DT, N_SUB)["latents"]
    lat_true = wm.encode_states(bK["feats"])
display(pd.Series({"JEPA rollout latents": E.straightness(lat_pred).mean(),
                   "JEPA encoded latents (truth)": E.straightness(lat_true).mean(),
                   "normalized state features (truth)": E.straightness(bK["feats"][:, P_:]).mean()},
                  name="straightness (1 = straight)").round(3))
""")

# ---------------------------------------------------------------------------
md(r"""
## 7. Closed loop — MPPI on the world model (paper §IV-D)

Every 5 s the planner samples bank sequences around its warm start. The noise lives on 6
knots spread over the predicted time-to-go. It rolls every sample through the
world model + PI prober to the drogue altitude in one batch and scores the
predicted miss, plus smoothness, g-load over 8 g, and flight-path angle steeper than
82° before drogue. That last one is a soft constraint on the prober's interpretable γ. The env
ends the flight at |γ| > 87°, and the classical controller avoids that with a hard
wings-level rule. MPPI then softmax-averages the samples and flies the first command. There is no
sign logic and no deadband: reversals come from the cost.

Head-to-head on the first `n_cases` frozen test-set ICs under two truths:

* **nominal** — the truth matches the classical predictor's model.
* **mismatch** — the truth has density ×0.88 and L/D ×0.88. That is inside the world model's
  training range, but the classical predictor still assumes nominal. Can 35 s of
  history do the implicit system ID the classical controller lacks?

`pred` is MPPI on the predictive baseline (the paper's closed-loop comparison).
Classical runs (~70 s each) are cached in `jepa_runs/` because they don't depend
on the learned model.
""")
code(r"""
CONDITIONS = {"nominal": dict(D.NOMINAL_DOMAIN),
              "mismatch": {"density_scale": 0.88, "lift_scale": 0.88, "mass_scale": 1.0}}
MPPI_CFG = {"n_samples": CFG["mppi_samples"], "n_sub": N_SUB}
ics = D.load_test_ics("configs/heldout_testset.csv", n=CFG["n_cases"])
cache_path = SIM_DIR / "jepa_runs" / "classical_cache.pkl"
classical_cache = pickle.load(open(cache_path, "rb")) if cache_path.exists() else {}

jobs = []
for tag, dom in CONDITIONS.items():
    for i, ic in enumerate(ics):
        if (tag, i) not in classical_cache:
            jobs.append({"controller": "classical", "ic": ic, "domain": dom, "tag": tag, "case": i, "record": True})
        for ctrl in ("jepa", "pred"):
            jobs.append({"controller": ctrl, "ic": ic, "domain": dom, "ckpt": str(CKPT), "mppi": MPPI_CFG,
                         "tag": tag, "case": i, "record": True})
jobs.sort(key=lambda j: j["controller"] != "classical")   # slow classical jobs first, for load balance
t0 = time.time()
new = D.run_parallel(P.run_episode, jobs, WORKERS, label="closed loop")
print(f"closed loop: {len(jobs)} episodes in {time.time() - t0:.0f}s")
for r in new:
    if r["controller"] == "classical":
        classical_cache[(r["tag"], r["case"])] = r
pickle.dump(classical_cache, open(cache_path, "wb"))
results = [r for r in new if r["controller"] != "classical"] + \
          [classical_cache[(t, i)] for t in CONDITIONS for i in range(len(ics))]
display(P.summarize(results))
""")
code(r"""
# One case, both truths: where each controller put the capsule, and how it banked.
def traj(tag, ctrl, case=0):
    return next(r for r in results if r["tag"] == tag and r["controller"] == ctrl and r["case"] == case)

fig, axes = plt.subplots(2, 2, figsize=(12, 6.6))
for row, tag in enumerate(CONDITIONS):
    for ctrl in ("classical", "jepa", "pred"):
        r = traj(tag, ctrl)
        T_ = r["traj"]
        t = np.arange(len(T_)) * MACRO_DT
        lab = f'{ctrl}: {r["miss_km"]:.1f} km'
        axes[row, 0].plot(np.degrees(T_[:, 2]), np.degrees(T_[:, 1]), color=COL[ctrl], label=lab)
        axes[row, 1].step(t, np.degrees(T_[:, 10]), where="post", color=COL[ctrl], label=ctrl)
    tphi, tlam = np.degrees(r["target"])
    axes[row, 0].plot([tlam], [tphi], marker="*", ms=14, color="#e34948", ls="none")
    # Zoom on the terminal approach: at full scale the tracks overlap and a km-level miss is invisible.
    axes[row, 0].set(xlim=(tlam - 1.2, tlam + 0.4), ylim=(tphi - 1.2, tphi + 0.4), xlabel="longitude [deg]",
                     ylabel="latitude [deg]", title=f"{tag}: terminal approach (case 0), 1 deg ≈ 111 km")
    axes[row, 0].legend(fontsize=8, loc="lower right")
    axes[row, 1].set(xlabel="time [s]", ylabel="bank command [deg]", title=f"{tag}: bank command")
fig.tight_layout()

fig, axes = plt.subplots(1, 2, figsize=(11, 3.4), sharey=True)
df = pd.DataFrame([{k: v for k, v in r.items() if k not in ("traj", "target")} for r in results])
for ax, tag in zip(axes, CONDITIONS):
    for ctrl in ("classical", "jepa", "pred"):
        m = np.sort(df[(df.tag == tag) & (df.controller == ctrl)].miss_km.values)
        ax.step(m, np.arange(1, len(m) + 1) / len(m), where="post", color=COL[ctrl], label=ctrl)
    ax.axvline(25.0, color=MUTED, lw=1, ls=":")
    ax.set(xscale="log", xlabel="miss at drogue deploy [km]", title=f"{tag} truth  (dotted: 25 km success radius)")
axes[0].set_ylabel("fraction of cases")
axes[0].legend(fontsize=8)
fig.tight_layout()
""")

# ---------------------------------------------------------------------------
md(r"""
## 8. Reading the results

The decision criteria are fixed in advance in [`jepa/PLAN.md`](jepa/PLAN.md) §6. In short:

* **The world model is sound** if the latent std is ≈1 with an effective rank well above 1,
  the JEPA+PI rollout error grows smoothly past the training horizon, and
  JEPA+PI beats both *Physics prior only* (the latent adds information) and
  *Predictive* (the paper's central claim) on the footprint at τ ≥ 200 s.
* **The method is interesting for guidance** if JEPA-MPPI succeeds (miss < 25 km) on most cases.
  That would already beat every learned controller so far (the BC clone is at 91 km).
  It is also interesting if it degrades *less* than the classical controller from nominal → mismatch.
* **Matching the 2 km classical bar is not expected** from this first pass. The
  plan lists the next levers: DAgger-style closed-loop data, a longer training
  horizon, and terminal-phase replanning at 1 Hz.

Planning cost per decision is in the `plan_ms` column. Compare it with ~124 ms per 1 s
guidance step for the classical predictor-corrector.
""")

nb["cells"] = cells
nb["metadata"] = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                  "language_info": {"name": "python"}}
nbf.write(nb, OUT)
print(f"wrote {OUT} ({len(cells)} cells)")
