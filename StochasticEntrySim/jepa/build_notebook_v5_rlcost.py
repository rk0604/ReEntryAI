"""
Build ReEntryAI_JEPA_v5_rl_cost.ipynb: measure what v5's RL work would cost on a Colab T4 (jepa/rlcost_v5.py).
Edit this builder, never the .ipynb:

    python jepa/build_notebook_v5_rlcost.py          (from StochasticEntrySim/)

Colab rules (CLAUDE.md): opens on a T4 runtime (the cheapest GPU, and the GPU the study's RL would use), saves its
result to Drive, then releases the runtime by itself (also when a cell fails); AUTO_DISCONNECT = False keeps it.
Outside Colab (the smoke test) it reads JEPA_V5_MODE=smoke and JEPA_V5_RLCOST_ROOT.
"""
import os
import sys

import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from build_notebooks_v3 import AUTO_DISCONNECT_CODE, DONE, DONE_MD  # noqa: E402

NAME = "ReEntryAI_JEPA_v5_rl_cost.ipynb"
BADGE = ("[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)]"
         "(https://colab.research.google.com/github/rk0604/ReEntryAI/blob/space-jepa5/StochasticEntrySim/" + NAME + ")")

INTRO = f"""
# SpaceJEPA v5 · What would RL cost on Colab?

{BADGE}

A short measurement, not part of the study: it writes nothing into `MyDrive/ReEntryAI_jepa_v5/`. It times the real
RL training code on this runtime with the study's frozen settings and projects the hours and Colab units for all of
v5's Tier 1 RL work (36 tuning configurations and 15 training runs), to decide between training RL on Colab and on
Kaggle's free T4 (`jepa/kaggle_rl.md`).

**How to run.** It opens on a **T4 GPU** runtime (keep it: T4 is the GPU the study's RL would use). **Runtime → Run
all.** About 30–45 minutes (≈ 1–1.5 units). The result is printed at the end and saved to
`MyDrive/ReEntryAI_v5_rl_cost/`. **The notebook disconnects itself** when it finishes or a cell fails, so it never
idles on paid units; set `AUTO_DISCONNECT = False` in the Setup cell to keep the runtime.

What it times: PPO's decisions per second (each PPO run stops at 200 M decisions or 3 h, whichever first), one chunk
of FNPEG-labelled episodes for DAgger and BC (not time-capped; it runs on the CPU), the supervised fits, the DEV48
policy flights each tuning configuration needs, and each run's fast evaluation. FNPEG flies its default settings
(tuning may change its cost somewhat), so treat the result as ±20 %.
"""

SETTINGS = r"""
import os, sys
IN_COLAB = "google.colab" in sys.modules
MODE = os.environ.get("JEPA_V5_MODE", "full")        # "smoke": tiny sizes (the smoke test; numbers meaningless)
MINUTES = 5.0 if MODE == "full" else 0.05            # minutes of PPO per measurement
"""

SETUP = r"""
BRANCH = "space-jepa5"
if IN_COLAB:
    REPO = "/content/ReEntryAI_v5"
    if not os.path.exists(REPO):
        !git clone -q -b {BRANCH} https://github.com/rk0604/ReEntryAI.git {REPO}
    else:
        !git -C {REPO} pull -q origin {BRANCH}
    %pip install -q pymsis
    os.chdir(f"{REPO}/StochasticEntrySim")
    from google.colab import drive
    drive.mount("/content/drive")
    OUT = "/content/drive/MyDrive/ReEntryAI_v5_rl_cost"
    SCRATCH = "/content/rl_cost_scratch"
else:
    OUT = os.environ.get("JEPA_V5_RLCOST_ROOT", "jepa_runs/v5_rl_cost")
    SCRATCH = os.path.join(OUT, "scratch")

import json, time
from pathlib import Path
sys.path.insert(0, os.getcwd())
import torch
from jepa import study_v5 as S5, rlcost_v5 as RC

CFG = S5.StudyConfigV5.smoke() if MODE == "smoke" else S5.StudyConfigV5()
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
GPU = torch.cuda.get_device_name(0) if DEVICE == "cuda" else "cpu"
CPUS = os.cpu_count() or 2
L = S5.LayoutV5(SCRATCH).make()                       # scratch only: never the study folder
os.makedirs(OUT, exist_ok=True)
print(f"MODE={MODE} device={DEVICE} ({GPU}) cpus={CPUS}; git {S5._GIT.get()[:10]}; result -> {OUT}")
if MODE == "full" and "T4" not in GPU:
    print(f"NOTE: this runtime's GPU is {GPU!r}, not a T4: the projection applies to this GPU, not to a T4.")
""" + AUTO_DISCONNECT_CODE

MEASURE = r"""
SMOKE = {} if MODE == "full" else dict(chunk=4, tune_chunk=2, n_flights=2, fast_eval_episodes=16,
                                       ppo_kw={"n_envs": 64, "n_steps": 8, "minibatch": 256})
t0 = time.time()
M = RC.measure(CFG, L, DEVICE, CPUS, MINUTES, **SMOKE)
print(f"measured in {(time.time() - t0) / 60:.1f} min")
"""

PROJECT = r"""
P = RC.project(CFG, M)
TEXT = RC.summary_text(M, P)
print(TEXT)
stamp = time.strftime("%Y%m%d_%H%M")
name = f"rl_cost_{GPU.replace(' ', '_')}_{stamp}"
Path(OUT, name + ".json").write_text(json.dumps({"measured": M, "projection": P, "git": S5._GIT.get()}, indent=1,
                                                default=str))
Path(OUT, name + ".txt").write_text(TEXT)
print(f"saved {OUT}/{name}.json and .txt")
"""

CELLS = [
    ("md", INTRO),
    ("md", "## 0. Settings"),
    ("code", SETTINGS),
    ("md", "## 1. Setup\n\nClone the `space-jepa5` code, mount Drive, a scratch folder on local disk, and the "
           "auto-disconnect hook. A minute."),
    ("code", SETUP),
    ("md", "## 2. Measure\n\nTwo 5-minute PPO runs (`ppo_hist`, `bcppo_fnpeg`'s PPO), one FNPEG-labelled chunk at the "
           "training size (200 episodes) and the tuning size (20), one supervised fit, a few DEV48 policy flights, one "
           "fast evaluation. 30–45 minutes."),
    ("code", MEASURE),
    ("md", "## 3. Projection\n\nEvery Tier 1 RL tuning configuration and training run with the study's settings: "
           "hours, Colab units on this runtime, and weeks of Kaggle's free quota. Saved to "
           "`MyDrive/ReEntryAI_v5_rl_cost/`."),
    ("code", PROJECT),
    ("md", DONE_MD),
    ("code", DONE),
]


def build() -> str:
    nb = nbf.v4.new_notebook()
    nb["cells"] = [nbf.v4.new_markdown_cell(t.strip("\n")) if k == "md" else nbf.v4.new_code_cell(t.strip("\n"))
                   for k, t in CELLS]
    for i, c in enumerate(nb["cells"]):
        c["id"] = f"v5c-cell-{i:02d}"
    nb["metadata"] = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                      "language_info": {"name": "python"}, "accelerator": "GPU",
                      "colab": {"gpuType": "T4", "provenance": []}}
    path = os.path.join(SIM, NAME)
    nbf.validate(nb)
    nbf.write(nb, path)
    print(f"wrote {path} ({len(CELLS)} cells)")
    return path


if __name__ == "__main__":
    build()
