"""
Build ReEntryAI_JEPA_v2.ipynb: train SpaceJEPA v2 on a Colab GPU and evaluate
it against v1 and the classical predictor-corrector.

All logic lives in the jepa/ package (models_v2, train_v2, evaluate_v2,
mppi_v2, plus v1's evaluate and mppi); the notebook configures, calls, and
plots. Edit this builder and re-run it rather than hand-editing the .ipynb:

    python jepa/build_notebook_v2.py        (from StochasticEntrySim/)
"""
import os
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "ReEntryAI_JEPA_v2.ipynb")

cells = []


def md(text):
    """Append a markdown cell."""
    cells.append(nbf.v4.new_markdown_cell(text.strip("\n")))


def code(text):
    """Append a code cell."""
    cells.append(nbf.v4.new_code_cell(text.strip("\n")))


# ---------------------------------------------------------------------------
md(r"""
# SpaceJEPA v2 — an attention-based world-model ensemble for entry guidance

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rk0604/ReEntryAI/blob/jepa-v2/StochasticEntrySim/ReEntryAI_JEPA_v2.ipynb)

v1 (`ReEntryAI_JEPA.ipynb`, 26k parameters, 320 narrow-corridor flights) matched the
classical predictor-corrector on nominal truth. But it showed a systematic bias when the air
is denser than nominal: 13 km median miss against 1 km for the classical controller.
v2 trains a much larger model on the **6,000-flight wide-envelope dataset** and targets that weakness directly.
The full design rationale is in [`jepa/PLAN.md`](jepa/PLAN.md) §11.

| | v1 | **v2** |
|---|---|---|
| Data | 320 flights, one entry corridor, density/L-D ×0.85–1.15 | **6,000 flights**, 100–140 km, 6.8–8.6 km/s, −1 to −6°, any heading, density/L-D ×0.7–1.3 |
| Inputs | 10 state features | + **physics residuals**: the drag and L/D corrections *measured* from the IMU at the current state |
| History encoder | causal TCN, 35 s | **banded causal transformer**, exact 120 s receptive field |
| Actions | TCN code | one token per 5 s step from the last 4 commands (multi-step action encoding) |
| Predictor | GRU | **causal transformer** over [context, action tokens]; parallel in training and planning |
| Objective | latent MSE + pooled SIGReg, T = 24 | masked latent MSE + **per-step SIGReg**, curriculum **T = 24 → 48 → 96**, λ sweep |
| Prober | physics-inspired (corrections to a crude aero model, exact 3DOF EOM) | same, plus the rollout's own altitude and speed as inputs; prober horizon up to 48 |
| Model count | 1 | **ensemble of 5**; disagreement on the landing point is an uncertainty penalty in MPPI |
| Parameters | 26k (+5k prober) | **2.45M per member, 12.3M for the ensemble** |

**What this notebook does.** (1) Loads the 6,000 flights (from your Drive, where
`ReEntryAI_JEPA_data.ipynb` saved them). (2) Sweeps the SIGReg weight. (3) Trains the
ensemble and its probers on the GPU, with checkpoints on Drive. (4) Retrains v1's
architecture on the same data (CPU, in the background), which separates "more data" from
"better model". (5) Evaluates open loop with the unchanged `evaluate.py`. (6) Flies closed loop
on the frozen test set under four truths with the unchanged `mppi.py` (classical, v1) and its
v2 extension `mppi_v2.py`.

**Runtime.** Use **Runtime → Change runtime type → GPU**. Rough totals for `MODE = "full"`:
**A100 ≈ 45–60 min, L4 ≈ 1.2–1.5 h, T4 ≈ 2.3–3 h** (cell 2 prints a refined estimate, and every
training stage prints its measured ms/step and ETA). **Everything is resumable.** Training
checkpoints every 10 min, and closed-loop flights are cached one by one on Drive. After a
disconnect, run all cells again: finished work is skipped.
""")

# ---------------------------------------------------------------------------
md(r"""
## 0. Setup

On Colab this clones the `jepa-v2` branch, installs `pymsis`, and mounts Google Drive.
The data is read from `MyDrive/ReEntryAI_jepa_data/wide_6000/`, where the data notebook
wrote it: the merged `episodes_wide_6000.pkl`, or its `shard_*.pkl` files if the merge step
was skipped. Outputs go to `MyDrive/ReEntryAI_jepa_v2/`. Locally, start Jupyter in
`StochasticEntrySim/` and point `JEPA_WIDE_DATA` at the pickle or shard folder.
""")
code(r"""
import os, sys
IN_COLAB = "google.colab" in sys.modules
BRANCH = "jepa-v2"
if IN_COLAB:
    REPO = "/content/ReEntryAI"
    if not os.path.exists(REPO):
        !git clone -q -b {BRANCH} https://github.com/rk0604/ReEntryAI.git {REPO}
    else:
        !git -C {REPO} pull -q
    %pip install -q pymsis
    os.chdir(f"{REPO}/StochasticEntrySim")
    from google.colab import drive
    drive.mount("/content/drive")
    DRIVE = "/content/drive/MyDrive"
    RUN_ROOT = f"{DRIVE}/ReEntryAI_jepa_v2"                                   # checkpoints, caches, results
    DATA_CANDIDATES = [f"{DRIVE}/ReEntryAI_jepa_data/wide_6000/episodes_wide_6000.pkl",
                       f"{DRIVE}/ReEntryAI_jepa_data/wide_6000"]               # shards, if the merge was skipped
    LOCAL_DATA = "/content/data"                                              # local-disk copy (faster rereads)
    !nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
else:
    RUN_ROOT = "jepa_runs/v2"
    DATA_CANDIDATES = [os.environ.get("JEPA_WIDE_DATA", ""), "jepa_runs/episodes_wide_6000.pkl",
                       "jepa_runs/wide_6000"]
    LOCAL_DATA = "jepa_runs"
print("working dir:", os.getcwd())
""")

code(r"""
import math, time, pickle, multiprocessing as mp
from pathlib import Path

SIM_DIR = Path.cwd()                                   # StochasticEntrySim/
assert (SIM_DIR / "rl_env.py").exists(), "start Jupyter in StochasticEntrySim/"
sys.path.insert(0, str(SIM_DIR))

import numpy as np
import pandas as pd
import torch
import matplotlib.pyplot as plt
%matplotlib inline

from jepa import data as D, models as M, evaluate as E, mppi as P                 # v1, unchanged
from jepa import models_v2 as V, train_v2 as T2, evaluate_v2 as E2, mppi_v2 as P2  # v2

MODE = os.environ.get("JEPA_V2_MODE", "full")   # "full": the experiment (GPU) | "smoke": every stage, ~10 min, numbers meaningless
RUN_SWEEP = True          # SIGReg lambda sweep before the main run (~5-10 min on a GPU)
TRAIN_V1_WIDE = True      # v1's architecture retrained on the wide data (CPU, background): data vs. model ablation
CFG = {
    "smoke": dict(
        n_episodes=160, members=3,
        arch=dict(d_latent=16, enc_dim=32, enc_layers=2, enc_heads=2, enc_band=5, pred_dim=48, pred_layers=2,
                  pred_heads=2, act_chunk=4, mlp_ratio=2),
        wm=dict(steps=300, batch=32, warmup=30, curriculum=((0.0, 8), (0.3, 16), (0.6, 24)), log_every=100, amp=True),
        pi=dict(steps=120, batch=32, warmup=10, lr=2e-3, curriculum=((0.0, 8), (0.5, 16)), log_every=40),
        sweep=dict(lams=(0.02, 0.08),
                   wm=dict(steps=100, batch=32, warmup=10, curriculum=((0.0, 8), (0.5, 16)), log_every=50),
                   pi=dict(steps=40, batch=32, warmup=5, lr=2e-3, curriculum=((0.0, 8),), log_every=20)),
        v1_wide=dict(steps_jepa=150, steps_pi=60),
        ttg_s=(300, 200, 100), k_eval=20, n_cases=1, mppi_samples=48, eval_eps=60, probe_eps=40),
    "full": dict(
        n_episodes=None, members=5,
        arch=dict(),        # models_v2 defaults: latent 64; encoder 128 wide x 3 layers, band 9 (rf 25 = 120 s); predictor 192 x 4
        wm=dict(steps=24_000, batch=256, warmup=1_000, lr=6e-4, curriculum=((0.0, 24), (0.2, 48), (0.5, 96)),
                log_every=500),
        pi=dict(steps=6_000, batch=512, warmup=300, lr=1e-3, weight_decay=1e-4,
                curriculum=((0.0, 16), (0.3, 32), (0.6, 48)), log_every=250),
        sweep=dict(lams=(0.01, 0.02, 0.04, 0.08, 0.16),
                   wm=dict(steps=4_000, batch=256, warmup=300, lr=6e-4, curriculum=((0.0, 24), (0.4, 48)), log_every=500),
                   pi=dict(steps=1_500, batch=512, warmup=100, lr=1e-3, weight_decay=1e-4,
                           curriculum=((0.0, 16), (0.5, 32)), log_every=250)),
        v1_wide=dict(steps_jepa=6_000, steps_pi=3_000),
        ttg_s=(400, 300, 200, 100, 50), k_eval=60, n_cases=16, mppi_samples=256, eval_eps=None, probe_eps=400),
}[MODE]

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
GPU_MEM_GB = torch.cuda.get_device_properties(0).total_memory / 1e9 if DEVICE.type == "cuda" else 0.0
if DEVICE.type == "cuda":
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True
    print(f"GPU: {torch.cuda.get_device_name(0)} ({GPU_MEM_GB:.0f} GB), "
          f"AMP dtype {V.amp_dtype(DEVICE)}")
elif MODE == "full":
    print("WARNING: no GPU. Runtime -> Change runtime type -> GPU (T4 works; L4/A100 are 2-5x faster).")
GRAD_CKPT = DEVICE.type == "cuda" and GPU_MEM_GB < 20   # T4: recompute activations (~1.3x compute) to fit 15 GB at T = 96
CPUS = os.cpu_count() or 2
BG_WORKERS = max(1, CPUS - 1)       # CPU pool: classical controller, v1, v1-wide (this process keeps a core for the GPU)
torch.set_num_threads(1 if DEVICE.type == "cuda" else max(1, CPUS // 2))
torch.manual_seed(0)

WM_CFG, PI_CFG = V.TrainCfgV2(**CFG["wm"]), V.TrainCfgV2(**CFG["pi"])
SW_WM_CFG, SW_PI_CFG = V.TrainCfgV2(**CFG["sweep"]["wm"]), V.TrainCfgV2(**CFG["sweep"]["pi"])
MACRO_DT, N_SUB, H_DROGUE = D.MACRO_DT_S, 2, 7600.0
LAMBDA_DEFAULT = 0.04     # used if RUN_SWEEP is off: per-step SIGReg at N = 256 rows ~ v1's 0.02 pooled over 1024
RUN_DIR = Path(RUN_ROOT) / MODE
RUN_DIR.mkdir(parents=True, exist_ok=True)
CKPT_V1 = SIM_DIR / "jepa_runs" / "full" / "jepa_ckpt.pt"   # the shipped v1 model (320 narrow-corridor flights)
CKPT_V1W = RUN_DIR / "v1_wide_ckpt.pt"
CKPT_V2 = RUN_DIR / "v2_ckpt.pt"
print(f"MODE={MODE}  device={DEVICE}  cpus={CPUS}  grad_ckpt={GRAD_CKPT}  run dir={RUN_DIR}")
""")

code(r"""
# Plot style: recessive grid, one color per entity everywhere (fixed categorical order).
plt.rcParams.update({
    "figure.dpi": 110, "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6,
    "axes.edgecolor": "#8a8984", "axes.labelcolor": "#52514e", "axes.titlecolor": "#0b0b0b",
    "xtick.color": "#52514e", "ytick.color": "#52514e", "lines.linewidth": 1.6, "legend.frameon": False,
})
COL = {"v1": "#2a78d6", "v2": "#eb6834", "v1-wide": "#1baf7a", "v2 (no unc.)": "#eda100",
       "classical": "#6250d6", "physics-nominal": "#6250d6", "Physics prior only": "#a3a29c"}
LS = {"v2 (no unc.)": "--", "physics-nominal": "--", "Physics prior only": ":"}   # second cue besides color
MUTED, MEMBER = "#a3a29c", "#f4b79a"
""")

# ---------------------------------------------------------------------------
md(r"""
## 1. Data: 6,000 wide-envelope flights

Same episode format as `jepa/data.py collect_episode`. The flights are random bank profiles (GP
and reversal-structured), flown through the full simulator under per-flight density, L/D, and
mass randomization. Each flight is split into 5 s macro steps, and the episodes are split
80/10/10 into train, validation, and test.
""")
code(r"""
t0 = time.time()
src = T2.find_wide_dataset(DATA_CANDIDATES, download_dir=LOCAL_DATA)
DATA_PATH = T2.local_copy(src, LOCAL_DATA)
episodes = T2.load_episodes(DATA_PATH)
if CFG["n_episodes"]:
    episodes = episodes[:CFG["n_episodes"]]
print(f"{len(episodes)} flights from {src} in {time.time() - t0:.0f}s")
ep_df = pd.DataFrame([{**e["ic"], **e["domain"], "kind": e["kind"], "outcome": e["outcome"],
                       "flight_s": e["t_final_s"]} for e in episodes])
display(ep_df.groupby(["kind", "outcome"]).size().unstack(fill_value=0))
display(ep_df[["altitude_m", "velocity_mps", "gamma_deg", "chi_deg", "density_scale", "lift_scale",
               "mass_scale", "flight_s"]].describe().loc[["min", "50%", "max"]].round(2))
""")
code(r"""
ARCH = dict(CFG["arch"])
RF = V.WorldModelV2(n_members=1, **ARCH).rf
MAX_T = max(T for _, T in WM_CFG.curriculum + PI_CFG.curriculum + SW_WM_CFG.curriculum + SW_PI_CFG.curriculum)
t0 = time.time()
data = T2.prepare(episodes, MAX_T, DEVICE, RF, split_seed=0)
test_eps = data["test"][:CFG["eval_eps"]] if CFG["eval_eps"] else data["test"]
print(f"episodes train/val/test = {len(data['train'])}/{len(data['val'])}/{len(data['test'])}; "
      f"{len(data['ds_tr']):,} training windows (every start, every horizon up to {MAX_T}); "
      f"history rf = {RF} steps = {RF * MACRO_DT:.0f} s; {time.time() - t0:.0f}s")

# v1's own data (shipped), for the coverage comparison and v1's home-turf test below.
v1_eps = pickle.load(open(SIM_DIR / "jepa_runs" / "episodes_900.pkl", "rb"))
_, v1_va, v1_te = D.split_episodes([D.to_macro(e, MACRO_DT) for e in v1_eps], (0.8, 0.1, 0.1), seed=0)
v1_home = (v1_te + v1_va)[:CFG["eval_eps"]] if CFG["eval_eps"] else v1_te + v1_va

fig, axes = plt.subplots(1, 3, figsize=(13, 3.4))
Yw = np.concatenate([e["y"][::10] for e in episodes[::max(1, len(episodes) // 600)]])
Yo = np.concatenate([e["y"][::10] for e in v1_eps])
axes[0].scatter(Yw[:, 3] / 1000, Yw[:, 0] / 1000, s=1, color=COL["v2"], alpha=0.25, label="wide (v2)")
axes[0].scatter(Yo[:, 3] / 1000, Yo[:, 0] / 1000, s=1, color=COL["v1"], alpha=0.4, label="v1 data")
axes[0].set(xlabel="speed [km/s]", ylabel="altitude [km]", title="State coverage")
axes[0].legend(markerscale=8, loc="lower right")
dens_w = ep_df.density_scale / ep_df.mass_scale
dens_o = np.array([e["domain"]["density_scale"] / e["domain"]["mass_scale"] for e in v1_eps])
bins = np.linspace(0.6, 1.45, 35)
axes[1].hist(dens_w, bins=bins, color=COL["v2"], alpha=0.7, label="wide (v2)", density=True)
axes[1].hist(dens_o, bins=bins, color=COL["v1"], alpha=0.7, label="v1 data", density=True)
for x, lab in ((0.88, "x0.88"), (1.12, "x1.12"), (0.78, "x0.78")):
    axes[1].axvline(x, color=MUTED, lw=1, ls=":")
axes[1].set(xlabel="truth density / mass scale", ylabel="density", title="Truth randomization (dotted: closed-loop truths)")
axes[1].legend()
axes[2].scatter(ep_df.velocity_mps / 1000, ep_df.gamma_deg, s=2, color=COL["v2"], alpha=0.4, label="wide (v2)")
ic_o = pd.DataFrame([e["ic"] for e in v1_eps])
axes[2].scatter(ic_o.velocity_mps / 1000, ic_o.gamma_deg, s=4, color=COL["v1"], label="v1 data")
axes[2].set(xlabel="entry speed [km/s]", ylabel="entry flight-path angle [deg]", title="Entry interface")
axes[2].legend(markerscale=3, loc="lower left")
fig.tight_layout()
""")

# ---------------------------------------------------------------------------
md(r"""
## 2. The v2 world model

Each component, and the reason for it (the full argument is in PLAN.md §11):

* **Physics-residual features.** `cD_meas = log(a_D / a_D,nominal(h, V))` and `cL_meas = log((L/D)/0.30)` are the
  corrections the prober needs, read straight off the IMU at the current state. v1 had to infer
  them as a nonlinear product of normalized features, and its dense-air bias suggests it did not do that well.
* **State encoder.** A causal transformer whose attention is banded (each layer sees ±40 s; three layers stack
  to an exact 120 s receptive field) with relative (rotary) positions. One pass over a window therefore still gives
  every sliding encoding, as v1's TCN did, and a longer memory helps system ID.
* **Action tokens.** One token per 5 s step, built from the last 4 commands (20 s), which is where the
  actuator lag lives.
* **Predictor.** A causal transformer over `[context latent, action tokens]`. The latent at step *k* sees only
  `s_t` and commands up to *k*, the same causal structure as v1's GRU. Training is one parallel pass, so
  a 96-step (480 s) horizon is affordable. MPPI also rolls all samples in one pass.
* **Ensemble.** 5 members, each with its own encoder, predictor, and prober, initialized and batched
  independently. They train as one set of batched matmuls (every weight has a leading member dimension).
""")
code(r"""
wm = V.WorldModelV2(n_members=CFG["members"], **ARCH).to(DEVICE)
wm.grad_ckpt = GRAD_CKPT
pi = V.PIProberV2(wm.M, wm.d_latent).to(DEVICE)
v1 = M.load_checkpoint(str(CKPT_V1))
rep = V.param_report(wm, pi)
display(pd.DataFrame({"v2": rep}).T.style.format("{:,}"))
print(f"v1: world model {M.n_params(v1['wm']):,} + PI prober {M.n_params(v1['pi']):,} parameters "
      f"-> v2 is {rep[f'ensemble (x{wm.M})'] / (M.n_params(v1['wm']) + M.n_params(v1['pi'])):.0f}x larger")

flops = sum(V.train_flops_per_step(wm, WM_CFG.batch, V.horizon_at(s, WM_CFG)) for s in range(WM_CFG.steps))
print(f"stage 1: {flops / 1e15:.1f} PFLOP over {WM_CFG.steps:,} steps "
      f"({V.train_flops_per_step(wm, WM_CFG.batch, max(T for _, T in WM_CFG.curriculum)) / 1e12:.2f} TFLOP/step at the longest horizon)")
# Effective rates for these narrow (128-192 wide) layers including the memory-bound ops; measured ms/step is
# printed while training, and that, not this, is the real number.
for gpu, tflops, ck in (("T4", 10, 1.33), ("L4", 18, 1.0), ("A100", 45, 1.0)):
    print(f"  ~{flops / (tflops * 1e12) * ck / 60:.0f} min on {gpu} (assumes ~{tflops} effective TFLOP/s"
          f"{', activation checkpointing' if ck > 1 else ''})")
""")
code(r"""
# Under the hood: the pieces that are new in v2, straight from jepa/models_v2.py and jepa/mppi_v2.py.
import inspect
for obj in (V.raw_features_v2, V.WorldModelV2.encode_states, V.WorldModelV2.action_tokens, V.WorldModelV2.predict,
            V.per_step_sigreg, V.wm_loss_v2, V.PIProberV2.corrections, P2.EnsembleMPPI._cost):
    print(inspect.getsource(obj))
    print("#" + "-" * 78)
""")

# ---------------------------------------------------------------------------
md(r"""
## 3. Start the CPU work in the background

Three things need CPUs but not the GPU, so they run in a background process pool while the GPU trains:

* **v1-wide**: v1's exact architecture and recipe (26k parameters, T = 24, models.fit), retrained on the
  wide training split. If v2 beats v1-wide, the gain is the model and not only the data.
* **Classical predictor-corrector** on the frozen test-set ICs under the four truths. It costs about 100 s per
  flight and does not depend on any learned model. Results are cached on Drive one flight at a time,
  seeded with the nominal/mismatch runs shipped with v1.

The four truths (the same tags as v1's notebook and `size_sweep.py`):

| tag | density | L/D | note |
|---|---|---|---|
| `nominal` | ×1.00 | ×1.00 | the classical predictor's own model |
| `mismatch` | ×0.88 | ×0.88 | thinner air |
| `mismatch_high` | ×1.12 | ×1.12 | **denser air: v1's weak spot** |
| `out_of_range` | ×0.78 | ×0.80 | outside v1's training range, inside v2's |
""")
code(r"""
CONDITIONS = {"nominal": dict(D.NOMINAL_DOMAIN),
              "mismatch": {"density_scale": 0.88, "lift_scale": 0.88, "mass_scale": 1.0},
              "mismatch_high": {"density_scale": 1.12, "lift_scale": 1.12, "mass_scale": 1.0},
              "out_of_range": {"density_scale": 0.78, "lift_scale": 0.80, "mass_scale": 1.0}}
ics = D.load_test_ics("configs/heldout_testset.csv", n=CFG["n_cases"])
cache = P2.ResultCache(Path(RUN_ROOT) / f"closed_loop_cache_{MODE}.pkl",
                       seed_classical=str(SIM_DIR / "jepa_runs" / "classical_cache.pkl"))
bg = mp.get_context("spawn").Pool(BG_WORKERS)
bg_tasks = {}
if TRAIN_V1_WIDE:
    bg_tasks["v1-wide"] = bg.apply_async(T2.train_v1_baseline, ({
        "episodes_path": DATA_PATH, "n_episodes": CFG["n_episodes"], "ckpt": str(CKPT_V1W), "split_seed": 0,
        **CFG["v1_wide"]},), error_callback=lambda e: print("v1-wide failed:", repr(e)))
n_new = 0
for tag, dom in CONDITIONS.items():
    for i, ic in enumerate(ics):
        key = ("classical", tag, i, "classical")
        if key not in cache:
            job = {"controller": "classical", "ic": ic, "domain": dom, "tag": tag, "case": i, "record": True,
                   "label": "classical", "key": key}
            bg_tasks[key] = bg.apply_async(P2.run_episode_keyed, (job,), callback=cache.put_result,
                                           error_callback=lambda e: print("classical flight failed:", repr(e)))
            n_new += 1
print(f"background: {'v1-wide training + ' if TRAIN_V1_WIDE else ''}{n_new} classical flights "
      f"({4 * len(ics) - n_new} cached) on {BG_WORKERS} worker(s); ~{n_new * 110 / BG_WORKERS / 60:.0f} min")
""")

# ---------------------------------------------------------------------------
md(r"""
## 4. SIGReg weight sweep

One ensemble member per λ is trained in a single batched run, then each member gets its prober and is
scored on **validation footprint error** (the geometric mean of the medians at τ = 300, 200, 100 s). The latent MSE
cannot be the selection metric: a more isotropic latent is harder to predict, so a lower MSE can just mean a
weaker regularizer. SIGReg here is per rollout step (the paper's form), with the statistic × N convention
of v1. λ = 0.04 at N = 256 rows applies about the same pressure as v1's 0.02 pooled over 1,024 rows.
""")
code(r"""
if RUN_SWEEP:
    t0 = time.time()
    sweep = T2.lambda_sweep(data, CFG["sweep"]["lams"], ARCH, SW_WM_CFG, SW_PI_CFG, DEVICE, str(RUN_DIR / "sweep"),
                            ttg_s=(300, 200, 100), h_cross=H_DROGUE, n_eval=150)
    sweep["table"].to_csv(RUN_DIR / "lambda_sweep.csv")
    display(sweep["table"].round(3))
    LAMBDA = sweep["best"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.0))
    t = sweep["table"]
    axes[0].plot(t.index, t.score_km, marker="o", color=COL["v2"])
    axes[0].set(xscale="log", xlabel="SIGReg weight λ", ylabel="val footprint error [km]", title="Selection metric (lower is better)")
    axes[1].plot(t.index, t.eff_rank, marker="o", color=COL["v2"])
    axes[1].axhline(wm.d_latent, color=MUTED, lw=1, ls=":")
    axes[1].set(xscale="log", xlabel="SIGReg weight λ", ylabel="effective rank", title=f"Latent spread (dotted: {wm.d_latent} = full)")
    fig.tight_layout()
    del sweep
    if DEVICE.type == "cuda":
        torch.cuda.empty_cache()
    print(f"sweep: {(time.time() - t0) / 60:.1f} min")
else:
    LAMBDA = LAMBDA_DEFAULT
print(f"SIGReg weight for the ensemble: lambda = {LAMBDA:g}")
""")

# ---------------------------------------------------------------------------
md(r"""
## 5. Stage 1: the JEPA ensemble

The loss is the masked multi-step latent MSE plus λ·SIGReg, averaged over rollout steps, on both the encoded
targets and the predictions. The horizon curriculum is 24 → 48 → 96 steps (120 → 480 s). Planning
rolls out up to 110 steps, so the planner only extrapolates 1.15× past training, against v1's 4.6×.
Windows run past the end of short flights under a mask. The loop checkpoints to Drive every 10 minutes
and resumes on rerun.

What to watch: `val_pred` falls at each horizon (it steps up when the horizon doubles, which is expected), and
`sigreg` settles. Collapse would show as a tiny `pred` together with a large `sigreg` and a low effective rank below.
""")
code(r"""
t0 = time.time()
hist_wm = T2.train_world_model(wm, data, WM_CFG, LAMBDA, DEVICE, str(RUN_DIR / "wm_stage.pt"))
print(f"stage 1: {(time.time() - t0) / 60:.1f} min this session, {hist_wm['sec'][-1] / 60:.1f} min of training in total")
""")
code(r"""
fig, axes = plt.subplots(1, 2, figsize=(11, 3.3))
st = np.array(hist_wm["step"])
axes[0].plot(st, hist_wm["pred"], color=COL["v2"], label="train")
axes[0].plot(st, [v.get("pred", np.nan) for v in hist_wm["val"]], color=COL["v2"], ls="--", label="val")
axes[1].plot(st, hist_wm["sigreg"], color=COL["v2"])
for ax in axes:
    for frac, T in WM_CFG.curriculum[1:]:
        ax.axvline(frac * WM_CFG.steps, color=MUTED, lw=1, ls=":")
        ax.text(frac * WM_CFG.steps, 1.0, f" T={T}", transform=ax.get_xaxis_transform(), va="top", fontsize=8, color="#52514e")
axes[0].set(xlabel="step", ylabel="latent MSE (mean over members)", yscale="log", title="Multi-step latent prediction")
axes[0].legend(loc="lower left")
axes[1].set(xlabel="step", ylabel="per-step Epps-Pulley statistic", yscale="log", title="SIGReg")
fig.tight_layout()
health = pd.DataFrame(T2.latent_health_v2(wm, data["ds_va"])).set_index("member")
display(health.round(2))
""")

# ---------------------------------------------------------------------------
md(r"""
## 6. Stage 2: state-aware physics-inspired probers

v1's prober is kept: the latent emits a log drag factor, a log L/D factor, and a bank offset that correct a
deliberately crude nominal aero model (one exponential atmosphere, one ballistic coefficient), pushed
through the exact 3DOF EOM. That was v1's biggest single win. v2 adds the **rollout's own altitude and
speed** as prober inputs. The density-profile and Mach parts of the correction depend on where the
integrated trajectory is, which the latent cannot know once the rollout drifts. The loss is pseudo-Huber on
dimensionless state errors, so 48-step (240 s) prober horizons stay stable, plus v1's auxiliary force match.
Each member trains its own prober on its own frozen latents.
""")
code(r"""
prior_only = V.PIProberV2(wm.M, wm.d_latent).to(DEVICE)          # zero corrections: the crude physics prior
gv = torch.Generator(device=data["ds_va"].device).manual_seed(7)
vbp = {T: data["ds_va"].sample((wm.M, 256), T, gv) for _, T in PI_CFG.curriculum}

def prober_val(p):
    with torch.no_grad():
        return {T: V.prober_loss_v2(wm, p, b, MACRO_DT, N_SUB, 0.0)[1]["state"] for T, b in vbp.items()}

t0 = time.time()
hist_pi = T2.train_prober(wm, pi, data, PI_CFG, DEVICE, str(RUN_DIR / "pi_stage.pt"))
V.save_checkpoint_v2(str(CKPT_V2), wm, pi, data["stats"],
                     extra={"mode": MODE, "cfg": CFG, "lambda": LAMBDA, "macro_dt": MACRO_DT, "n_sub": N_SUB})
print(f"stage 2: {(time.time() - t0) / 60:.1f} min this session; saved -> {CKPT_V2}")
display(pd.DataFrame({"Physics prior only": prober_val(prior_only), "v2 PI prober": prober_val(pi)})
        .rename_axis("prober horizon T").round(3))
fig, ax = plt.subplots(figsize=(6, 3.2))
ax.plot(hist_pi["step"], hist_pi["state"], color=COL["v2"], label="train")
ax.plot(hist_pi["step"], [v.get("state", np.nan) for v in hist_pi["val"]], color=COL["v2"], ls="--", label="val")
for frac, T in PI_CFG.curriculum[1:]:
    ax.axvline(frac * PI_CFG.steps, color=MUTED, lw=1, ls=":")
ax.set(xlabel="step", ylabel="pseudo-Huber state loss", yscale="log", title="Prober rollout loss (dotted: horizon steps)")
ax.legend()
fig.tight_layout()
""")

# ---------------------------------------------------------------------------
md(r"""
## 7. Open loop

All on held-out flights, with the true future commands. The headline metric uses the **unchanged
`evaluate.footprint_errors`**: from the state at time-to-go τ, where does each model say the drogue
deploy will be, and how far is that from where it really was?

* `v2` is the ensemble mean: the members' prober corrections are averaged and integrated once (one
  deterministic model, through v1's `jepa_predict`). Each `v2 member i` is shown too.
* `physics-nominal` is the classical predictor's model (high fidelity, nominal parameters, instant bank).
  `Physics prior only` is the crude prior with zero corrections.
* `v1` is the shipped v1 model, trained on 320 narrow-corridor flights, so most of this test set is new
  to it. `v1-wide` is the same architecture retrained on this data.
""")
code(r"""
if TRAIN_V1_WIDE:
    t0 = time.time()
    print("v1-wide:", bg_tasks["v1-wide"].get(), f"(waited {time.time() - t0:.0f}s)")
    v1w = M.load_checkpoint(str(CKPT_V1W))
models_v2 = {"wm": wm.eval(), "pi": pi.eval(), "stats": data["stats"]}
FAM_V1 = {"v1": {"kind": "jepa", "wm": v1["wm"], "prober": v1["pi"]},
          "Physics prior only": {"kind": "jepa", "wm": v1["wm"], "prober": M.PIProber(v1["wm"].d_latent)}}
FAM_V1W = {"v1-wide": {"kind": "jepa", "wm": v1w["wm"], "prober": v1w["pi"]}} if TRAIN_V1_WIDE else {}
FAM_V2 = {"v2": V.as_v1(models_v2), **{f"v2 member {m}": V.as_v1(models_v2, [m]) for m in range(wm.M)}}
MAIN = ["physics-nominal", "Physics prior only", "v1", "v1-wide", "v2"]

def footprint(eps, ttg_s, nominal=True):
    # Each family has its own features and history length, so one evaluate.footprint_errors call each;
    # all calls select the same (episode, time-to-go) pairs.
    rows = E.footprint_errors(FAM_V1, eps, v1["stats"], v1["wm"].rf, ttg_s=ttg_s, h_cross=H_DROGUE, dt=MACRO_DT,
                              n_sub=N_SUB, nominal=E.nominal_params() if nominal else None)
    if FAM_V1W:
        rows += E.footprint_errors(FAM_V1W, eps, v1w["stats"], v1w["wm"].rf, ttg_s=ttg_s, h_cross=H_DROGUE,
                                   dt=MACRO_DT, n_sub=N_SUB)
    rows += E.footprint_errors(FAM_V2, eps, data["stats"], wm.rf, ttg_s=ttg_s, h_cross=H_DROGUE, dt=MACRO_DT, n_sub=N_SUB)
    return pd.DataFrame(rows)

t0 = time.time()
fp = footprint(test_eps, CFG["ttg_s"])
fp.to_csv(RUN_DIR / "footprint_wide_test.csv", index=False)
med = fp.groupby(["model", "ttg_s"]).err_km.median().unstack(0)
print(f"footprint eval: {time.time() - t0:.0f}s on {fp[fp.model == 'v2'].shape[0] // len(CFG['ttg_s'])} flights")
fig, ax = plt.subplots(figsize=(6.6, 3.8))
for i, c in enumerate([c for c in med.columns if c.startswith("v2 member")]):
    ax.plot(med.index, med[c], color=MEMBER, lw=1, label="v2 members" if i == 0 else None)
for name in [m for m in MAIN if m in med.columns]:
    ax.plot(med.index, med[name], marker="o", ms=5, color=COL[name], ls=LS.get(name, "-"), label=name)
ax.set(xlabel="time-to-go when predicting [s]", ylabel="median footprint error [km]", yscale="log",
       title="Predicted vs. true drogue-deploy point (wide test set)")
ax.invert_xaxis()
ax.legend(fontsize=8)
fig.tight_layout()
display(med[[m for m in MAIN if m in med.columns] + [c for c in med.columns if c.startswith("v2 member")]].round(2))
""")

md(r"""
### The dense-air question

v1's weakness was a *systematic* bias in denser-than-nominal air. A systematic error shows up as a nonzero
**mean along-track error** that tracks the truth density (positive means the model predicts the landing beyond the
true point, i.e. it under-predicts drag). The plot below bins the test flights by their truth density/mass
scale (the ballistic coefficient, which is what the forces reveal) and plots that signed error.
""")
code(r"""
fp["truth"] = pd.cut(fp.dom_density_scale / fp.dom_mass_scale, [0, 0.85, 1.15, 10],
                     labels=["thin (<0.85)", "near nominal", "dense (>1.15)"])
display(fp[fp.model.isin(MAIN)].pivot_table(index="model", columns=["ttg_s", "truth"], values="err_km",
                                             aggfunc="median", observed=False).round(1))
SIGNED = {"v1": {**FAM_V1["v1"], "stats": v1["stats"], "rf": v1["wm"].rf},
          **({"v1-wide": {**FAM_V1W["v1-wide"], "stats": v1w["stats"], "rf": v1w["wm"].rf}} if FAM_V1W else {}),
          "v2": {**FAM_V2["v2"], "rf": wm.rf}}
pts = pd.DataFrame(E2.footprint_points(SIGNED, test_eps, (300, 200, 100), H_DROGUE, dt=MACRO_DT, n_sub=N_SUB))
pts["beta_scale"] = pts.dom_density_scale / pts.dom_mass_scale
edges = np.linspace(0.6, 1.45, 9)
pts["bin"] = pd.cut(pts.beta_scale, edges)
fig, axes = plt.subplots(1, 3, figsize=(13, 3.4), sharey=True)
for ax, ttg in zip(axes, (300, 200, 100)):
    for name in SIGNED:
        g = pts[(pts.model == name) & (pts.ttg_s == ttg)].groupby("bin", observed=True).along_km
        mid = [b.mid for b in g.mean().index]
        ax.plot(mid, g.mean().values, marker="o", ms=5, color=COL[name], label=name)
    ax.axhline(0.0, color=MUTED, lw=1)
    ax.axvline(1.0, color=MUTED, lw=1, ls=":")
    ax.set(xlabel="truth density / mass scale", title=f"τ = {ttg} s")
axes[0].set_ylabel("mean along-track error [km]\n(+ = predicts beyond the truth)")
axes[0].legend(fontsize=8)
fig.tight_layout()
display(pts[pts.beta_scale > 1.15].groupby(["model", "ttg_s"])[["along_km", "cross_km", "err_km"]].mean().round(2)
        .rename(columns=lambda c: f"mean {c} (dense > 1.15)"))
""")

md(r"""
### v1's home turf, ensemble calibration, implicit system ID, rollout error

* **v1's own test set** (narrow corridor, never seen by v2): does the wide model give up accuracy on v1's distribution?
* **Calibration**: MPPI penalizes member disagreement, which only helps if the spread tracks the actual error.
* **System ID**: a linear (ridge) readout of the latent, fit on training flights, predicting each test flight's truth
  factors, by time since entry. The reference is the current features alone (no memory). v2's features already
  contain the measured drag ratio, so the memory has to beat that.
""")
code(r"""
fp_home = footprint(v1_home, CFG["ttg_s"])
fp_home.to_csv(RUN_DIR / "footprint_v1_home.csv", index=False)
med_home = fp_home.groupby(["model", "ttg_s"]).err_km.median().unstack(0)
display(med_home[[m for m in MAIN if m in med_home.columns]].round(2)
        .style.set_caption("median footprint error [km] on v1's own held-out flights"))

sp = pd.DataFrame(E2.member_spread(models_v2, test_eps, CFG["ttg_s"], H_DROGUE, dt=MACRO_DT, n_sub=N_SUB))
calib = sp.groupby("ttg_s").apply(lambda g: pd.Series({
    "spread_med_km": g.spread_km.median(), "err_med_km": g.err_mean_km.median(),
    "spearman_rho": g.spread_km.rank().corr(g.err_mean_km.rank()),
    "err<=2*spread": (g.err_mean_km <= 2 * g.spread_km).mean()}), include_groups=False)
display(calib.round(2))
fig, ax = plt.subplots(figsize=(4.8, 4.0))
s2 = sp[sp.ttg_s == 200]
ax.scatter(s2.spread_km.clip(lower=1e-2), s2.err_mean_km.clip(lower=1e-2), s=8, color=COL["v2"], alpha=0.6)
lim = [1e-2, max(1.0, float(s2[["spread_km", "err_mean_km"]].max().max()) * 1.5)]
ax.plot(lim, lim, color=MUTED, lw=1, ls=":")
ax.set(xscale="log", yscale="log", xlabel="member spread [km]", ylabel="ensemble-mean error [km]",
       title="Spread vs. error at τ = 200 s (dotted: y = x)")
fig.tight_layout()
""")
code(r"""
probe = {"v1": {"wm": v1["wm"], "stats": v1["stats"], "rf": v1["wm"].rf},
         **({"v1-wide": {"wm": v1w["wm"], "stats": v1w["stats"], "rf": v1w["wm"].rf}} if TRAIN_V1_WIDE else {}),
         "v2 (member 0)": {"wm": V.as_v1(models_v2, [0])["wm"], "stats": data["stats"], "rf": wm.rf},
         "current features only": {"features": data["stats"]}}
sid = pd.DataFrame(E2.system_id_probe(probe, data["train"], test_eps, max_eps=CFG["probe_eps"]))
display(sid.pivot(index="model", columns="t_since_entry", values=["R2_density/mass", "R2_L/D"]).round(2))

K = CFG["k_eval"]
curves = E.rollout_errors(FAM_V1, test_eps, v1["stats"], v1["wm"].rf, K=K, n_windows=512, dt=MACRO_DT, n_sub=N_SUB)
if FAM_V1W:
    curves.update(E.rollout_errors(FAM_V1W, test_eps, v1w["stats"], v1w["wm"].rf, K=K, n_windows=512, dt=MACRO_DT, n_sub=N_SUB))
curves.update(E.rollout_errors({"v2": FAM_V2["v2"]}, test_eps, data["stats"], wm.rf, K=K, n_windows=512, dt=MACRO_DT, n_sub=N_SUB))
kk = np.arange(1, K + 1) * MACRO_DT
fig, axes = plt.subplots(1, 3, figsize=(13, 3.2))
for ax, (q, lab) in zip(axes, [("pos_km", "position RMSE [km]"), ("h_km", "altitude RMSE [km]"), ("V_mps", "speed RMSE [m/s]")]):
    for name, c in curves.items():
        ax.plot(kk, c[q], color=COL[name], ls=LS.get(name, "-"), label=name)
    ax.axvline(24 * MACRO_DT, color=COL["v1"], lw=1, ls=":")
    ax.set(xlabel="prediction horizon [s]", ylabel=lab, yscale="log")
axes[0].set_title("dotted: v1 training horizon (v2 trains to 96 steps = 480 s)", loc="left", fontsize=8)
axes[0].legend(fontsize=8)
fig.tight_layout()
""")

# ---------------------------------------------------------------------------
md(r"""
## 8. Closed loop on the frozen test set

Every controller flies the same frozen test-set ICs under the same four truths, through the same actuator chain.

* **classical**: `OrionEntryEnv` baseline mode; its predictor keeps the nominal model (cached, from section 3).
* **v1**, **v1-wide**: the unchanged `mppi.run_episode` (MPPI on the world model + PI prober), CPU pool.
* **v2**: `mppi_v2.run_lockstep`. This is v1's MPPI (same samples, knots, temperature, and cost terms, also run
  with 256 samples), with the cost averaged over the 5 members plus 1 km per km of member disagreement on the
  landing point. All flights of a truth run in lockstep, so each 5 s decision is one batched GPU rollout.
  With one member and no uncertainty term it makes the same decisions as `mppi.run_episode`; `smoke_v2.py`
  checks this.
* **v2 (no unc.)**: the same ensemble with the disagreement penalty off. It isolates what the uncertainty term buys.

The g-limit is v1's (8 g soft limit), so v1 and v2 are compared like for like. PLAN.md §9 notes that v1 was not
g-matched to the classical controller. Watch `peak_g_mean`.
""")
code(r"""
MPPI_V1 = {"n_samples": CFG["mppi_samples"], "n_sub": N_SUB}
V2_VARIANTS = {"v2": {"w_unc_km": 1.0}, "v2 (no unc.)": {"w_unc_km": 0.0}}
IDS = {"v1": P2.file_id(str(CKPT_V1)), "v2": P2.file_id(str(CKPT_V2))}
if TRAIN_V1_WIDE:
    IDS["v1-wide"] = P2.file_id(str(CKPT_V1W))

def key(label, tag, i):
    if label == "classical":
        return ("classical", tag, i, "classical")
    return (label, tag, i, f"{IDS['v2' if label.startswith('v2') else label]}/{CFG['mppi_samples']}")

# v1 and v1-wide on the CPU pool, through the unchanged mppi.run_episode.
for label, ck in (("v1", CKPT_V1), ("v1-wide", CKPT_V1W)):
    if label not in IDS:
        continue
    for tag, dom in CONDITIONS.items():
        for i, ic in enumerate(ics):
            k = key(label, tag, i)
            if k not in cache:
                job = {"controller": "jepa", "ic": ic, "domain": dom, "ckpt": str(ck), "mppi": MPPI_V1, "tag": tag,
                       "case": i, "record": True, "label": label, "key": k}
                bg_tasks[k] = bg.apply_async(P2.run_episode_keyed, (job,), callback=cache.put_result,
                                             error_callback=lambda e: print("flight failed:", repr(e)))
# v2 on the GPU in this process, one truth at a time (cached after each).
models_gpu = V.load_checkpoint_v2(str(CKPT_V2), DEVICE)
t_all = time.time()
for tag, dom in CONDITIONS.items():
    jobs = [{"ic": ic, "domain": dom, "tag": tag, "case": i, "record": True, "label": label, "key": key(label, tag, i),
             "mppi": over} for label, over in V2_VARIANTS.items() for i, ic in enumerate(ics)
            if key(label, tag, i) not in cache]
    if jobs:
        t0 = time.time()
        res = P2.run_lockstep(jobs, models_gpu, MPPI_V1, n_sub=N_SUB, verbose=False)
        for j, r in zip(jobs, res):
            cache.put(j["key"], r)
        print(f"v2 {tag}: {len(jobs)} flights in {time.time() - t0:.0f}s")
print(f"v2 closed loop: {time.time() - t_all:.0f}s; waiting for the CPU pool ...")
t0 = time.time()
for t in bg_tasks.values():
    t.wait()
bg.close()
bg.join()
print(f"CPU pool done (+{time.time() - t0:.0f}s)")

LABELS = ["classical", "v1", "v1-wide", "v2", "v2 (no unc.)"]
results = [cache.get(key(l, t, i)) for l in LABELS for t in CONDITIONS for i in range(len(ics))
           if (l != "v1-wide" or TRAIN_V1_WIDE) and cache.get(key(l, t, i)) is not None]
summary = P.summarize(results)
summary.to_csv(RUN_DIR / "closed_loop_summary.csv")
pd.DataFrame([{k_: v for k_, v in r.items() if k_ not in ("traj", "target", "key")} for r in results]) \
    .to_csv(RUN_DIR / "closed_loop.csv", index=False)
display(summary)
""")
code(r"""
df = pd.DataFrame([{k_: v for k_, v in r.items() if k_ not in ("traj", "target", "key")} for r in results])
present = [l for l in LABELS if l in set(df.controller)]
display(df.pivot_table(index="controller", columns="tag", values="miss_km", aggfunc="median").loc[present, list(CONDITIONS)]
        .round(2).style.set_caption("median miss [km] by truth"))
fig, axes = plt.subplots(1, len(CONDITIONS), figsize=(16, 3.4), sharey=True)
for ax, tag in zip(axes, CONDITIONS):
    for c in present:
        m = np.sort(df[(df.tag == tag) & (df.controller == c)].miss_km.values)
        ax.step(m, np.arange(1, len(m) + 1) / len(m), where="post", color=COL[c], ls=LS.get(c, "-"), label=c)
    ax.axvline(25.0, color=MUTED, lw=1, ls=":")
    ax.set(xscale="log", xlabel="miss at drogue deploy [km]", title=tag)
axes[0].set_ylabel("fraction of cases")
axes[0].legend(fontsize=8)
fig.tight_layout()

# The dense-air case 0: ground track near the target and the bank command, per controller.
def flight(label, tag, case=0):
    return cache.get(key(label, tag, case))

fig, axes = plt.subplots(1, 2, figsize=(12, 3.8))
for c in [l for l in ("classical", "v1", "v2") if flight(l, "mismatch_high") is not None]:
    r = flight(c, "mismatch_high")
    T_ = r["traj"]
    axes[0].plot(np.degrees(T_[:, 2]), np.degrees(T_[:, 1]), color=COL[c], label=f'{c}: {r["miss_km"]:.1f} km')
    axes[1].step(np.arange(len(T_)) * MACRO_DT, np.degrees(T_[:, 10]), where="post", color=COL[c], label=c)
tphi, tlam = np.degrees(r["target"])
axes[0].plot([tlam], [tphi], marker="*", ms=14, color="#0b0b0b", ls="none")
axes[0].set(xlim=(tlam - 1.2, tlam + 0.4), ylim=(tphi - 1.2, tphi + 0.4), xlabel="longitude [deg]",
            ylabel="latitude [deg]", title="mismatch_high (density x1.12), case 0: terminal approach")
axes[0].legend(fontsize=8, loc="lower right")
axes[1].set(xlabel="time [s]", ylabel="bank command [deg]", title="bank command")
fig.tight_layout()
""")

# ---------------------------------------------------------------------------
md(r"""
## 9. Decision criteria (fixed before the full run; PLAN.md §11.5)

| # | Criterion | Met if |
|---|---|---|
| 1 | No collapse | every member's effective rank > d_latent / 4 (16 of 64) |
| 2 | The model, not just the data | v2 median footprint error < v1-wide's at τ = 200 s and 300 s (wide test set) |
| 3 | Beats the classical predictor's model | v2 < physics-nominal at τ = 200 s and 300 s |
| 4 | Dense-air bias gone | for truth density/mass > 1.15, \|mean along-track error\| at τ = 200 s < 2 km and < v1's |
| 5 | Closed loop, robust | v2 success ≥ 95% under all four truths, and median miss ≤ classical's under `mismatch_high` |
| 6 | Uncertainty is informative | spread–error Spearman ρ > 0.4 at τ = 200 s, and v2's p90 miss ≤ v2 (no unc.)'s under ≥ 3 of 4 truths |

Criteria 1–3 decide whether v2 replaces v1 as the world model. Criterion 4 is the question that motivated v2.
Criteria 5–6 decide whether the ensemble planner is worth its cost. n = 16 cases is still small: treat
differences under ~1 km in the median as noise, and rerun with `n_cases = 48` on a high-CPU runtime before
claiming anything about the classical controller.
""")
code(r"""
def m_at(frame, model, ttg):
    s = frame[(frame.model == model) & (frame.ttg_s == ttg)].err_km
    return float(s.median()) if len(s) else np.nan

def p(df_, ctrl, tag, fn):
    s = df_[(df_.controller == ctrl) & (df_.tag == tag)].miss_km
    return float(fn(s)) if len(s) else np.nan

crit = []
er = health.eff_rank.min()
crit.append(("1 no collapse", f"min effective rank {er:.1f} / {wm.d_latent}", er > wm.d_latent / 4))
v2a, v2b = m_at(fp, "v2", 200), m_at(fp, "v2", 300)
if TRAIN_V1_WIDE:
    wa, wb = m_at(fp, "v1-wide", 200), m_at(fp, "v1-wide", 300)
    crit.append(("2 model, not data", f"v2 {v2a:.1f}/{v2b:.1f} km vs v1-wide {wa:.1f}/{wb:.1f} km (200/300 s)", v2a < wa and v2b < wb))
na, nb = m_at(fp, "physics-nominal", 200), m_at(fp, "physics-nominal", 300)
crit.append(("3 beats physics-nominal", f"v2 {v2a:.1f}/{v2b:.1f} km vs {na:.1f}/{nb:.1f} km", v2a < na and v2b < nb))
dense = pts[(pts.beta_scale > 1.15) & (pts.ttg_s == 200)].groupby("model").along_km.mean()
b2, b1 = dense.get("v2", np.nan), dense.get("v1", np.nan)
crit.append(("4 dense-air bias", f"mean along-track v2 {b2:+.1f} km, v1 {b1:+.1f} km", abs(b2) < 2 and abs(b2) < abs(b1)))
succ = [100 * df[(df.controller == "v2") & (df.tag == t)].success.mean() for t in CONDITIONS]
mh_v2, mh_cl = p(df, "v2", "mismatch_high", np.median), p(df, "classical", "mismatch_high", np.median)
crit.append(("5 closed loop", f"v2 success {min(succ):.0f}% (worst truth); mismatch_high median v2 {mh_v2:.2f} vs "
             f"classical {mh_cl:.2f} km", min(succ) >= 95 and mh_v2 <= mh_cl))
rho = float(calib.loc[200.0, "spearman_rho"]) if 200.0 in calib.index else np.nan
wins = sum(p(df, "v2", t, lambda s: s.quantile(0.9)) <= p(df, "v2 (no unc.)", t, lambda s: s.quantile(0.9)) for t in CONDITIONS)
crit.append(("6 uncertainty", f"Spearman rho {rho:.2f}; p90 no worse than no-unc. in {wins}/4 truths", rho > 0.4 and wins >= 3))
verdict = pd.DataFrame(crit, columns=["criterion", "measured", "met"])
verdict.to_csv(RUN_DIR / "decision_criteria.csv", index=False)
display(verdict)
print(f"all outputs in {RUN_DIR}")
""")

nb = nbf.v4.new_notebook()
nb["cells"] = cells
nb["metadata"] = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                  "language_info": {"name": "python"}, "accelerator": "GPU",
                  "colab": {"gpuType": "T4", "provenance": []}}
nbf.write(nb, OUT)
print(f"wrote {OUT} ({len(cells)} cells)")
