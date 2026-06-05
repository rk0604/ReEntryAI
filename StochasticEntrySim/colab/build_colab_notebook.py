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
    "# Clone the repo (re-clone fresh each session)\n"
    "import os\n"
    f'REPO_URL = "{REPO_URL}"\n'
    'if not os.path.isdir("ReEntryAI"):\n'
    "    !git clone --depth 1 {REPO_URL}\n"
    "else:\n"
    '    print("repo already cloned")'
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
# (Optional) Mount Google Drive for persistent checkpoints + logs.
# Skip this cell if you don't want Drive; artifacts then live in the Colab VM.
USE_DRIVE = True
if USE_DRIVE:
    from google.colab import drive
    drive.mount("/content/drive")
    ARTIFACT_DIR = "/content/drive/MyDrive/ReEntryAI"
else:
    ARTIFACT_DIR = "/content/ReEntryAI_artifacts"
os.makedirs(ARTIFACT_DIR, exist_ok=True)
print("artifacts ->", ARTIFACT_DIR)
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
import run_baseline, multiprocessing, shutil
from pathlib import Path

LIMIT   = 50                       # set to None for the full 1,000-case set
WORKERS = multiprocessing.cpu_count()
BASE_OUT = Path(ARTIFACT_DIR) / "baseline_metrics.csv"

rows = run_baseline.run_batch(
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

**Reward note:** the env's default `w_range` (1e-5 / m) makes the terminal miss
penalty small next to the +50 on-target bonus. We pass a stronger `w_range`
here for a denser gradient toward the target — tune this.
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
TOTAL_STEPS = 1_000_000
RUN_NAME = "ppo_orion_entry"

# Stronger miss penalty than the env default for a trainable gradient.
REWARD_WEIGHTS = telemetry.RewardWeights(w_range=1.0e-3)

def make_env(rank: int, seed: int = 0):
    def _init():
        import rl_env, telemetry
        e = rl_env.OrionEntryEnv(
            controller_mode="policy", dispersion=True,
            max_episode_steps=1400,
            reward_weights=telemetry.RewardWeights(w_range=1.0e-3))
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
                }, step=self.num_timesteps)
        return True
""")

code(r"""
run = wandb.init(project="ReEntryAI", name=RUN_NAME, sync_tensorboard=True,
                 config={"total_steps": TOTAL_STEPS, "n_envs": N_ENVS,
                         "algo": "PPO", "w_range": 1.0e-3})

vec = SubprocVecEnv([make_env(i) for i in range(N_ENVS)])

model = PPO("MlpPolicy", vec, verbose=1, device="auto",
            n_steps=1024, batch_size=4096, gamma=0.999, gae_lambda=0.95,
            ent_coef=0.0, learning_rate=3e-4,
            tensorboard_log=f"{ARTIFACT_DIR}/tb/{run.id}")

ckpt = CheckpointCallback(save_freq=max(1, 50_000 // N_ENVS),
                          save_path=f"{ARTIFACT_DIR}/checkpoints/{run.id}",
                          name_prefix="ppo")
cbs = [WandbCallback(gradient_save_freq=0, verbose=1), MissionMetrics(), ckpt]

model.learn(total_timesteps=TOTAL_STEPS, callback=cbs, progress_bar=True)
model.save(f"{ARTIFACT_DIR}/{RUN_NAME}_final")
vec.close()
print("training done ->", f"{ARTIFACT_DIR}/{RUN_NAME}_final.zip")
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
model = PPO.load(f"{ARTIFACT_DIR}/{RUN_NAME}_final")

# Single deterministic env (no dispersion sampling; ICs come from the test set).
eval_env = rl_env.OrionEntryEnv(controller_mode="policy", dispersion=False)
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
                 "rcs_fuel_kg": info["rcs_fuel_kg"], "success": int(info["success"])})
    if (i + 1) % 50 == 0:
        print(f"  [{i+1}/{len(cases)}]")

policy_df = pd.DataFrame(rows)
policy_df.to_csv(f"{ARTIFACT_DIR}/policy_metrics.csv", index=False)
print("policy success rate: %.1f%%" % (100.0 * policy_df["success"].mean()))
""")

code(r"""
# Head-to-head comparison (RL vs classical) on identical conditions
import pandas as pd
base_df   = pd.read_csv(f"{ARTIFACT_DIR}/baseline_metrics.csv")
policy_df = pd.read_csv(f"{ARTIFACT_DIR}/policy_metrics.csv")

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
