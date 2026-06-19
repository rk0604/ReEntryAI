"""
Generator for the ReEntryAI RL Colab notebook.

Writing .ipynb JSON by hand is error-prone, so we build it programmatically:
each cell is a (kind, source) pair and we emit a valid nbformat-4 notebook.
Re-run this whenever the notebook needs to change, then commit both files.

    python colab/build_colab_notebook.py
"""

from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "ReEntryAI_RL_Colab.ipynb"

REPO_URL = "https://github.com/rk0604/ReEntryAI.git"

cells = []


def md(src: str) -> None:
    cells.append(("markdown", src.strip("\n")))


def code(src: str) -> None:
    cells.append(("code", src.strip("\n")))


# ---------------------------------------------------------------------------
md(r"""
# ReEntryAI — RL surrogate controller (Colab)

Train a reinforcement-learning bank-angle policy for Orion crew-module entry
guidance and benchmark it against the classical predictor-corrector on a
**frozen, shared** held-out test set.

**All heavy compute lives here (Colab), not on the local machine.** This
notebook covers:

| Phase | What this notebook does |
|---|---|
| 3 | Characterize the classical baseline on the 1,000-case held-out set |
| 4 | Train a PPO policy (`OrionEntryEnv`, dispersed ICs) with W&B logging |
| 5 | Evaluate the trained policy on the **same** held-out set; compare |

The simulator, Gym env (`rl_env.OrionEntryEnv`), dispersion, and test set were
all built in Phases 0–2 and live in the repo. See `docs/rl_mdp_spec.md`.

> **Runtime tip:** set `Runtime → Change runtime type → GPU` and (if you have
> Colab Pro) a **High-RAM / more-CPU** machine. The env step is CPU-bound, so
> the number of vCPUs drives training throughput more than the GPU.
""")

# ---------------------------------------------------------------------------
md("## 0 · Environment setup")

code(r"""
# GPU / machine check
import multiprocessing, os
print("CPUs:", multiprocessing.cpu_count())
!nvidia-smi -L || echo "no GPU (CPU-only is fine; training is mostly CPU-bound)"
""")

code(
    "# Clone the repo, or pull the latest if it's already here (so code fixes\n"
    "# pushed to main show up without deleting the runtime).\n"
    "import os\n"
    f'REPO_URL = "{REPO_URL}"\n'
    'if not os.path.isdir("ReEntryAI"):\n'
    "    !git clone --depth 1 {REPO_URL}\n"
    "else:\n"
    "    !git -C ReEntryAI pull --ff-only\n"
    '    print("pulled latest")'
)

code(r"""
# Install dependencies
!pip -q install pymsis gymnasium "stable-baselines3>=2.3.0" wandb pandas pyarrow tensorboard
print("deps installed")
""")

code(r"""
# Put the simulator package on the path and import the env
import sys, os
SIM_DIR = os.path.abspath("ReEntryAI/StochasticEntrySim")
assert os.path.isdir(SIM_DIR), SIM_DIR
os.chdir(SIM_DIR)                 # configs/ paths resolve relative to here
if SIM_DIR not in sys.path:
    sys.path.insert(0, SIM_DIR)

import numpy as np
import rl_env, mission_config, telemetry
from make_test_set import load_test_set
print("env module loaded:", rl_env.__file__)
print("observation features:", rl_env.OBS_NAMES)
""")

code(r"""
# Persist everything we generate to Google Drive under MyDrive/ReEntryAI so it
# survives Colab disconnects (and so the cache/resume cells below can reuse it).
# Set USE_DRIVE = False to keep artifacts only in the ephemeral Colab VM.
USE_DRIVE = True
if USE_DRIVE:
    from google.colab import drive
    drive.mount("/content/drive")
    ARTIFACT_DIR = "/content/drive/MyDrive/ReEntryAI"
else:
    ARTIFACT_DIR = "/content/ReEntryAI_artifacts"

# Subfolders for the different artifact kinds.
DATA_DIR  = f"{ARTIFACT_DIR}/data"          # metrics CSVs (baseline, policy)
CKPT_DIR  = f"{ARTIFACT_DIR}/checkpoints"   # periodic SB3 checkpoints (resume)
MODEL_DIR = f"{ARTIFACT_DIR}/models"        # final saved models
TB_DIR    = f"{ARTIFACT_DIR}/tb"            # tensorboard logs
for d in (ARTIFACT_DIR, DATA_DIR, CKPT_DIR, MODEL_DIR, TB_DIR):
    os.makedirs(d, exist_ok=True)
print("artifacts ->", ARTIFACT_DIR)
!ls -la "{ARTIFACT_DIR}"
""")

# ---------------------------------------------------------------------------
md(r"""
## 1 · Sanity check the env

Confirm the Gym contract (SB3's `check_env`) and run one quick episode so we
know imports, spaces, and physics all work before committing compute.
""")

code(r"""
from stable_baselines3.common.env_checker import check_env

env = rl_env.OrionEntryEnv(controller_mode="policy", dispersion=True)
check_env(env, warn=True)        # raises if the Gym API is violated

obs, info = env.reset(seed=0)
print("obs shape:", obs.shape, "| action space:", env.action_space)
done = 0
for _ in range(2000):
    obs, r, term, trunc, info = env.step(env.action_space.sample())
    if term or trunc:
        done = 1; break
print("one random episode ->", info["outcome"],
      "| miss %.1f km | peak_g %.1f | dispersed IC alt=%.0f m"
      % (info["miss_km"], info["peak_g"], info["ic"]["altitude_m"]))
""")

# ---------------------------------------------------------------------------
md(r"""
## 2 · Phase 3 — classical baseline on the held-out set

Run the predictor-corrector over the frozen test set to get the metric
distributions RL must beat. This is embarrassingly parallel; set `WORKERS` to
the number of vCPUs. Each predictor run is ~40 s, so 1,000 cases ≈
`1000 * 40 / WORKERS` seconds.

Start with a subset (`LIMIT`) to estimate wall-time, then run the full set.
""")

code(r"""
import run_baseline, multiprocessing
from pathlib import Path

LIMIT         = 50                 # set to None for the full 1,000-case set
WORKERS       = multiprocessing.cpu_count()
FORCE_BASELINE = False             # True = recompute even if the cache exists
BASE_OUT = Path(DATA_DIR) / "baseline_metrics.csv"

# Cache: skip the (slow) baseline if it's already on Drive.
if BASE_OUT.exists() and not FORCE_BASELINE:
    print(f"[cache] baseline already exists -> {BASE_OUT} (set FORCE_BASELINE=True to redo)")
else:
    run_baseline.run_batch(
        config_path="configs/default.json",
        test_set_path="configs/heldout_testset.csv",
        out_path=BASE_OUT, limit=LIMIT, workers=WORKERS, fast_atmosphere=True)
    print("baseline metrics ->", BASE_OUT)
""")

code(r"""
# Inspect the baseline distribution
import pandas as pd
base_df = pd.read_csv(BASE_OUT)
print("baseline cases:", len(base_df), "| success rate %.1f%%"
      % (100.0 * base_df["success"].mean()))
base_df[["miss_km", "peak_g", "peak_qdot_MWm2", "rcs_fuel_kg"]].describe(
    percentiles=[0.5, 0.9, 0.99])
""")

# ---------------------------------------------------------------------------
md(r"""
## 3 · Phase 4 — train the RL policy (PPO)

`OrionEntryEnv` in **policy** mode with **dispersed ICs**. Observations are
already normalized in the env, so no `VecNormalize` is needed. We run several
envs in parallel (`SubprocVecEnv`) because the bottleneck is CPU-side physics.

**Reward (rebuilt).** The previous run's steering signal was ~20x too weak, so
PPO collapsed into a lazy ~200 km-miss policy that barely used the thrusters
(only ~1% success, barely better than a random policy). This version makes
steering the dominant signal:
- a strong **dense progress reward** (`W_PROGRESS_PER_KM`, per km of
  range-to-target closed each step),
- a **heading-shaping** term (`W_HEADING_PER_RAD`, reward for aiming at the
  target even while range-to-go is still flat),
- a small **entropy bonus** (`ENT_COEF`) so the policy actually explores banking,
- and a **curriculum**: train with `DISPERSION=False` first so it can master one
  trajectory before facing scattered entry conditions.

**Run it in two stages:**
1. **Diagnostic** (the defaults below): `TOTAL_STEPS=150_000`, `DISPERSION=False`.
   Watch `ep/miss_km` in W&B — it must fall well below 100 km within ~150k steps
   (~12 min). If it does, the reward works; go to stage 2.
2. **Full run:** set `TOTAL_STEPS=1_000_000`, `DISPERSION=True`, and rename
   `RUN_NAME` to `..._v5_full` (new name = fresh checkpoints + curve).

> ⚠️ **A new reward needs a new `RUN_NAME`**, or the resume logic reloads the old
> policy and skips training.

**Interval ablation ladder.** The env exposes three independent interval knobs
(`INTERVAL_OBS`, `INTERVAL_REWARD`, `INTERVAL_SHIELD`, set in the cell below).
Run them as an ablation, each with its own `RUN_NAME`:
1. all off — plain RL baseline
2. `INTERVAL_OBS` — uncertainty-aware policy
3. `+ INTERVAL_REWARD` — penalize the worst-case heat/g
4. `+ INTERVAL_SHIELD` — provably-safe RL (1-step reachability veto)

Watch `ep/interval_violations` (should fall toward zero) and
`ep/shield_intervention_rate` (the policy learning the safe envelope). The
shield makes constraint satisfaction certified rather than hoped-for, which is
the core safety contribution. Note the shield is ~10x slower per step, so use it
once the reward is tuned.
""")

code(r"""
import wandb
wandb.login()        # paste your API key when prompted (or set WANDB_API_KEY)
""")

code(r"""
from stable_baselines3 import PPO
from stable_baselines3.common.vec_env import SubprocVecEnv
from stable_baselines3.common.monitor import Monitor
from stable_baselines3.common.callbacks import CheckpointCallback, BaseCallback
from wandb.integration.sb3 import WandbCallback

N_ENVS   = max(1, multiprocessing.cpu_count())

# === run knobs ============================================================
# STAGE 1 (diagnostic, the defaults here): short run, NO dispersion (curriculum),
# strong steering reward, exploration on. Watch ep/miss_km -- it should drop well
# below 100 km within ~150k steps (~12 min). If it does, the reward works.
# STAGE 2 (full run): set TOTAL_STEPS=1_000_000, DISPERSION=True, and bump
# RUN_NAME to ..._v5_full (a new name = fresh curve + checkpoints).
TOTAL_STEPS = 150_000
DISPERSION  = False
RUN_NAME    = "ppo_orion_entry_v4_diag"

# Reward shaping (the real fix). The v3 run's steering reward was ~20x too weak,
# so PPO sat in a lazy ~200 km-miss local optimum (barely better than random).
# Make steering the dominant signal, add a heading "aim at the target" term, and
# turn on exploration so the policy actually tries banking.
W_PROGRESS_PER_KM = 1.0      # was 0.05 -- reward per km of range-to-target closed/step
W_HEADING_PER_RAD = 10.0     # reward per rad of heading-to-target error removed/step
ENT_COEF          = 0.01     # was 0.0 -- entropy bonus so PPO explores banking
REWARD_WEIGHTS = telemetry.RewardWeights(w_range=1.0e-3, w_fuel=0.05)
# ==========================================================================

# --- Interval-arithmetic ablation switches (the 3 ways intervals enter RL) ---
# Flip these to run the ablation ladder. Use a DISTINCT RUN_NAME per setting.
#   INTERVAL_OBS    : append worst-case interval features to the observation
#   INTERVAL_REWARD : penalize the interval upper bound on heat / g
#   INTERVAL_SHIELD : veto unsafe banks via 1-step interval reachability
INTERVAL_OBS    = False
INTERVAL_REWARD = False
INTERVAL_SHIELD = False

def make_env(rank: int, seed: int = 0):
    def _init():
        import rl_env, telemetry
        e = rl_env.OrionEntryEnv(
            controller_mode="policy", dispersion=DISPERSION,
            max_episode_steps=1400,
            reward_weights=telemetry.RewardWeights(w_range=1.0e-3, w_fuel=0.05),
            w_progress_per_km=W_PROGRESS_PER_KM, w_heading_per_rad=W_HEADING_PER_RAD,
            interval_obs=INTERVAL_OBS, interval_reward=INTERVAL_REWARD,
            interval_shield=INTERVAL_SHIELD)
        e = Monitor(e)
        e.reset(seed=seed + rank)
        return e
    return _init
""")

code(r"""
# Custom callback: log per-episode mission metrics (miss, peak-g, success...)
# to W&B, beyond SB3's default reward/length curves -> "log everything".
class MissionMetrics(BaseCallback):
    def _on_step(self) -> bool:
        for info in self.locals.get("infos", []):
            if "episode" in info:                 # episode just ended
                wandb.log({
                    "ep/miss_km":        info.get("miss_km", float("nan")),
                    "ep/peak_g":         info.get("peak_g", float("nan")),
                    "ep/peak_qdot_MWm2": info.get("peak_qdot_MWm2", float("nan")),
                    "ep/rcs_fuel_kg":    info.get("rcs_fuel_kg", float("nan")),
                    "ep/success":        float(info.get("success", False)),
                    "ep/return":         info["episode"]["r"],
                    "ep/length":         info["episode"]["l"],
                    # interval diagnostics (zero unless an interval mode is on)
                    "ep/peak_g_hi":            info.get("peak_g_hi", float("nan")),
                    "ep/interval_violations":  info.get("interval_violations", 0),
                    "ep/shield_intervention_rate": info.get("shield_intervention_rate", 0.0),
                    "global_step":       self.num_timesteps,
                })   # no explicit step= : tensorboard syncing owns the step
        return True
""")

code(r"""
import glob, re

# Resume from W&B (same run id) so curves continue across Colab disconnects.
run = wandb.init(project="ReEntryAI", id=RUN_NAME, name=RUN_NAME, resume="allow",
                 sync_tensorboard=True,
                 config={"total_steps": TOTAL_STEPS, "n_envs": N_ENVS,
                         "algo": "PPO", "w_range": 1.0e-3, "w_fuel": 0.05,
                         "w_progress": W_PROGRESS_PER_KM, "w_heading": W_HEADING_PER_RAD,
                         "dispersion": DISPERSION, "ent_coef": ENT_COEF})

vec = SubprocVecEnv([make_env(i) for i in range(N_ENVS)])

# Cache/resume: pick up the latest checkpoint on Drive if one exists.
def _latest_checkpoint(d):
    cks = glob.glob(f"{d}/ppo_*_steps.zip")
    if not cks:
        return None, 0
    def steps(p):
        m = re.search(r"ppo_(\d+)_steps", p); return int(m.group(1)) if m else 0
    latest = max(cks, key=steps)
    return latest, steps(latest)

ckpt_run_dir = f"{CKPT_DIR}/{RUN_NAME}"
os.makedirs(ckpt_run_dir, exist_ok=True)
resume_path, done_steps = _latest_checkpoint(ckpt_run_dir)

if resume_path:
    print(f"[resume] loading checkpoint {resume_path} ({done_steps} steps done)")
    model = PPO.load(resume_path, env=vec, device="auto",
                     tensorboard_log=f"{TB_DIR}/{RUN_NAME}")
else:
    print("[fresh] no checkpoint found; starting a new run")
    # batch_size divides the rollout buffer (n_steps * N_ENVS) for any vCPU
    # count, so SB3 never silently clamps it; net_arch [128,128] gives the
    # 13-16 dim obs a bit more capacity than the [64,64] default.
    model = PPO("MlpPolicy", vec, verbose=1, device="auto",
                n_steps=1024, batch_size=1024, n_epochs=10,
                gamma=0.999, gae_lambda=0.95, clip_range=0.2,
                ent_coef=ENT_COEF, learning_rate=3e-4,
                policy_kwargs=dict(net_arch=[128, 128]),
                tensorboard_log=f"{TB_DIR}/{RUN_NAME}")

ckpt = CheckpointCallback(save_freq=max(1, 50_000 // N_ENVS),
                          save_path=ckpt_run_dir, name_prefix="ppo")
cbs = [WandbCallback(gradient_save_freq=0, verbose=1), MissionMetrics(), ckpt]

remaining = max(0, TOTAL_STEPS - done_steps)
model.learn(total_timesteps=remaining, callback=cbs, progress_bar=True,
            reset_num_timesteps=(resume_path is None))
model.save(f"{MODEL_DIR}/{RUN_NAME}_final")
vec.close()
print("training done ->", f"{MODEL_DIR}/{RUN_NAME}_final.zip")
""")

# ---------------------------------------------------------------------------
md(r"""
## 4 · Phase 5 — evaluate the policy on the held-out set

Run the trained policy on the **same** frozen ICs the baseline used (via
`reset(options={"ic_override": ...})`), collect the identical metrics, and
compare the two distributions head-to-head. This is the apples-to-apples result.
""")

code(r"""
from stable_baselines3 import PPO
import pandas as pd

EVAL_LIMIT = None                      # None = all 1,000 held-out cases
model = PPO.load(f"{MODEL_DIR}/{RUN_NAME}_final")

# Deterministic env (ICs come from the test set). MUST use the SAME interval
# flags as training so the obs dimension matches and the shield is applied.
eval_env = rl_env.OrionEntryEnv(
    controller_mode="policy", dispersion=False,
    interval_obs=INTERVAL_OBS, interval_reward=INTERVAL_REWARD,
    interval_shield=INTERVAL_SHIELD)
cases = load_test_set("configs/heldout_testset.csv")
if EVAL_LIMIT:
    cases = cases[:EVAL_LIMIT]

rows = []
for i, ic in enumerate(cases):
    obs, info = eval_env.reset(options={"ic_override": ic})
    done = False
    while not done:
        action, _ = model.predict(obs, deterministic=True)
        obs, r, term, trunc, info = eval_env.step(action)
        done = term or trunc
    rows.append({"idx": i, "outcome": info["outcome"], "miss_km": info["miss_km"],
                 "peak_g": info["peak_g"], "peak_qdot_MWm2": info["peak_qdot_MWm2"],
                 "rcs_fuel_kg": info["rcs_fuel_kg"], "success": int(info["success"]),
                 "interval_violations": info.get("interval_violations", 0),
                 "shield_intervention_rate": info.get("shield_intervention_rate", 0.0)})
    if (i + 1) % 50 == 0:
        print(f"  [{i+1}/{len(cases)}]")

policy_df = pd.DataFrame(rows)
policy_df.to_csv(f"{DATA_DIR}/policy_metrics.csv", index=False)
print("policy metrics -> %s/policy_metrics.csv" % DATA_DIR)
print("policy success rate: %.1f%%" % (100.0 * policy_df["success"].mean()))
""")

code(r"""
# Head-to-head comparison (RL vs classical) on identical conditions
import pandas as pd
base_df   = pd.read_csv(f"{DATA_DIR}/baseline_metrics.csv")
policy_df = pd.read_csv(f"{DATA_DIR}/policy_metrics.csv")

def summary(df):
    return pd.Series({
        "success_%":   100.0 * df["success"].mean(),
        "miss_median": df["miss_km"].median(),
        "miss_p90":    df["miss_km"].quantile(0.90),
        "peak_g_med":  df["peak_g"].median(),
        "fuel_med":    df["rcs_fuel_kg"].median(),
    })

cmp = pd.DataFrame({"classical": summary(base_df), "RL": summary(policy_df)})
print(cmp.round(2))
""")

code(r"""
# Distribution plots: RL vs classical
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
for col, a, title in zip(
        ["miss_km", "peak_g", "rcs_fuel_kg"], ax,
        ["Miss distance [km]", "Peak g", "RCS fuel [kg]"]):
    a.hist(base_df[col], bins=40, alpha=0.6, label="classical")
    a.hist(policy_df[col], bins=40, alpha=0.6, label="RL")
    a.set_title(title); a.legend()
plt.tight_layout(); plt.show()
""")

md(r"""
---
### Next steps / tuning knobs
- **Reward shaping** — `w_range`, on-target/survive bonuses, optional progress
  shaping. Watch `ep/miss_km` in W&B; if it plateaus high, raise `w_range`.
- **Throughput** — more vCPUs / `N_ENVS`; `LIMIT`/`EVAL_LIMIT` for quick passes.
- **Harder scenario** — swap `config_path="configs/scenario_a_lunar.json"`
  (lunar return) once Scenario B works.
- **Fairness** — both controllers already share ICs, observations, action
  authority, physics, and metrics (see `docs/rl_mdp_spec.md` §5).
""")

# ---------------------------------------------------------------------------
nb = {
    "cells": [
        {"cell_type": k, "metadata": {},
         **({"outputs": [], "execution_count": None} if k == "code" else {}),
         "source": s}
        for (k, s) in cells
    ],
    "metadata": {
        "colab": {"provenance": [], "toc_visible": True},
        "kernelspec": {"display_name": "Python 3", "name": "python3"},
        "language_info": {"name": "python"},
        "accelerator": "GPU",
    },
    "nbformat": 4,
    "nbformat_minor": 0,
}

OUT.write_text(json.dumps(nb, indent=1), encoding="utf-8")
print(f"wrote {OUT}  ({len(cells)} cells)")
