"""
Build ReEntryAI_JEPA_v5_kaggle_rl.ipynb: every v5 RL run (tuning and training) on Kaggle's free GPU, one platform
and one GPU type for all of them (jepa/kaggle_rl.md). Edit this builder, never the .ipynb:

    python jepa/build_notebook_v5_kaggle.py          (from StochasticEntrySim/)

Outside Kaggle (the smoke test) the notebook reads JEPA_V5_MODE, JEPA_V5_KAGGLE_LIVE (the session's study folder),
JEPA_V5_KAGGLE_OUT (stands in for /kaggle/working), JEPA_V5_KAGGLE_INPUTS (stands in for /kaggle/input, os.pathsep
separated), JEPA_V5_KAGGLE_HOURS and JEPA_V5_KAGGLE_MAX_RUNS.

Not a Colab notebook, so CLAUDE.md's Colab auto-disconnect does not apply: a Kaggle "Save & Run All" version ends by
itself, and this notebook stops starting runs before Kaggle's 12-hour limit. A failed cell still saves the state.
"""
import os

import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
NAME = "ReEntryAI_JEPA_v5_kaggle_rl.ipynb"
BRANCH = "space-jepa5"

INTRO = r"""
# SpaceJEPA v5 · RL on Kaggle

Every v5 RL run, tuning and training, trains here, on **one GPU type (T4)** and one pinned Kaggle environment, so
all RL arms are comparable (`jepa/kaggle_rl.md`). Each run stops at its frozen budget (decisions or 3 hours of wall
time, whichever first); the GPU and the decisions reached are recorded per run. Everything else in v5 runs in
`ReEntryAI_JEPA_v5.ipynb` on Colab.

**One-time setup**
1. On Colab, run `ReEntryAI_JEPA_v5.ipynb` until FNPEG is tuned; its section 11 then writes
   `MyDrive/ReEntryAI_jepa_v5/00_study/kaggle/to_kaggle/v5_rl_inputs.json`. Download that file.
2. Kaggle → Datasets → New dataset `reentryai-v5-rl-inputs` with that one file.
3. Kaggle → Create → New notebook → File → Import notebook → GitHub: `rk0604/ReEntryAI`, branch `space-jepa5`,
   `StochasticEntrySim/ReEntryAI_JEPA_v5_kaggle_rl.ipynb`.
4. Notebook settings: **Accelerator: GPU T4 x2** (one GPU is used), **Internet: on**,
   **Environment: pin to original environment**. Add input: the dataset of step 2.

**Each session**
5. **Save Version → Save & Run All (Commit).** It runs in the background (≤ 12 h), starts no run it cannot finish
   in time, and saves `v5_rl_state.tar` and `v5_rl_manifest.json` as its output after every run (and if a cell fails).
6. When it has finished: **Add input → Your Work → this notebook's output**, the newest version (an older attached
   version is ignored when a newer one is attached). Then Save Version again. Repeat until the last cell prints
   **RL COMPLETE** (≈ 55 T4-hours at most for Tier 1, about two weeks of the free 30 h/week).
7. Any time: download the newest version's two output files into `MyDrive/ReEntryAI_jepa_v5/00_study/kaggle/from_kaggle/`
   and run section 11 of the Colab notebook; it imports the finished runs (checked by sha256, nothing overwritten).

Interactive sessions use quota while idle and lose their files when stopped: use Save & Run All.
"""

SETTINGS = r"""
import os, sys, time
T_START = time.time()
IN_KAGGLE = os.path.isdir("/kaggle/working")
MODE = os.environ.get("JEPA_V5_MODE", "full")      # "full" | "smoke" (the smoke test; numbers meaningless)
RL_GPU = "T4"                                     # the one GPU type of every RL run (refuses any other in full mode)
SESSION_HOURS = float(os.environ.get("JEPA_V5_KAGGLE_HOURS", "11.5"))   # Kaggle stops a version at 12 h
MAX_RUNS = int(os.environ["JEPA_V5_KAGGLE_MAX_RUNS"]) if os.environ.get("JEPA_V5_KAGGLE_MAX_RUNS") else None
CPU_WORKERS = None                                # DEV48 tuning flights (None = every core)
"""

SETUP = r"""
BRANCH = "space-jepa5"
if IN_KAGGLE:
    REPO = "/tmp/ReEntryAI_v5"                    # not under /kaggle/working: only the state files are output
    if not os.path.exists(REPO):
        !git clone -q -b {BRANCH} https://github.com/rk0604/ReEntryAI.git {REPO}
    %pip install -q pymsis
    os.chdir(f"{REPO}/StochasticEntrySim")
    LIVE, OUT, INPUT_ROOTS = "/tmp/ReEntryAI_jepa_v5", "/kaggle/working", ["/kaggle/input"]
    os.environ["JEPA_V5_RL_PLATFORM"] = "kaggle"
else:
    LIVE, OUT = os.environ["JEPA_V5_KAGGLE_LIVE"], os.environ["JEPA_V5_KAGGLE_OUT"]
    INPUT_ROOTS = os.environ["JEPA_V5_KAGGLE_INPUTS"].split(os.pathsep)
    os.environ.setdefault("JEPA_V5_RL_PLATFORM", "local-smoke")
sys.path.insert(0, os.getcwd())
import json
import pandas as pd, torch
from jepa import study_v5 as S5, rl_v5 as RL, kaggle_v5 as K, benchmark_v5 as B

CFG = S5.StudyConfigV5.smoke() if MODE == "smoke" else S5.StudyConfigV5()
CFG.tier = 1
RL_ARMS = list(S5.RL_TIER1)
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
GPU = torch.cuda.get_device_name(0) if DEVICE == "cuda" else "cpu"
if MODE == "full" and RL_GPU not in GPU:
    raise RuntimeError(f"this session has {GPU!r}; every v5 RL run trains on a {RL_GPU} "
                       "(notebook settings -> Accelerator -> GPU T4 x2)")
CPUS = os.cpu_count() or 2
CPU_WORKERS = CPU_WORKERS or (CPUS if MODE == "full" else min(2, CPUS))
torch.set_num_threads(1 if DEVICE == "cuda" else max(1, CPUS // 2))     # as the Colab notebook

L = S5.LayoutV5(LIVE).make()
LOG = L.logger("kaggle")
SEQ = K.restore(L, INPUT_ROOTS, LOG) + 1          # this session's number
S5.repair(L, LOG)
K.copy_inputs(L, INPUT_ROOTS, LOG)
B.verify()
FROZEN = S5.freeze_config_v5(CFG, L)              # the Colab study's settings (raises if this code's differ)
GIT = S5._GIT.get()
HW = K.check_platform(L, RL.hardware_info(DEVICE), LOG)
LIMIT = S5.SessionLimit(deadline=T_START + SESSION_HOURS * 3600, max_runs=MAX_RUNS,
                        after_run=lambda: K.pack(L, OUT, SEQ, LOG))

def _save_on_error(result):
    # a failed cell (Run all stops there) still leaves this session's state in the output
    err = result.error_in_exec
    if err is not None and not isinstance(err, KeyboardInterrupt):
        try:
            K.pack(L, OUT, SEQ, LOG)
        except Exception as e:
            print("could not save the state:", e)

_ip = get_ipython()
for _f in list(_ip.events.callbacks["post_run_cell"]):
    if getattr(_f, "__name__", "") == "_save_on_error":
        _ip.events.unregister("post_run_cell", _f)
_ip.events.register("post_run_cell", _save_on_error)
K.pack(L, OUT, SEQ, LOG)
LOG(f"session {SEQ}: MODE={MODE} device={DEVICE} ({GPU}) cpus={CPUS} git={GIT[:10]} (study frozen at "
    f"{str(FROZEN.get('git_commit'))[:10]}); RL platform {HW}; stops starting runs after "
    f"{SESSION_HOURS} h")
"""

TUNE = r"""
S5.tune_rl(CFG, L, RL_ARMS, DEVICE, CPU_WORKERS, LOG, limit=LIMIT)
for a in RL_ARMS:
    p = L.tuning / f"rl_{a}" / "chosen.json"
    print(a, "->", json.loads(p.read_text())["score"] if p.exists() else "not chosen yet")
"""

TRAIN = r"""
S5.train_rl_stage(CFG, L, RL_ARMS, DEVICE, LOG, limit=LIMIT)
"""

SUMMARY = r"""
MAN = K.pack(L, OUT, SEQ, LOG)
rows = []
for a in RL_ARMS:
    for s in CFG.rl_seeds:
        d = L.policy_dir(a, s)
        r = json.loads((d / "run.json").read_text()) if (d / "run.json").exists() else {}
        rows.append({"arm": a, "seed": s, "done": (d / "policy.pt").exists(),
                     "decisions": r.get("decisions"), "train_h": round((r.get("train_s") or 0) / 3600, 2),
                     "gpu": (r.get("hardware") or {}).get("gpu")})
display(pd.DataFrame(rows))
n_tuned = sum((L.tuning / f"rl_{a}" / "chosen.json").exists() for a in RL_ARMS)
n_done = sum(r["done"] for r in rows)
if n_tuned == len(RL_ARMS) and n_done == len(rows):
    print("RL COMPLETE: download v5_rl_state.tar and v5_rl_manifest.json (Output) into "
          "MyDrive/ReEntryAI_jepa_v5/00_study/kaggle/from_kaggle/ and run section 11 of the Colab notebook.")
else:
    print(f"session {SEQ} saved: {n_tuned}/{len(RL_ARMS)} methods tuned, {n_done}/{len(rows)} runs trained. "
          "NEXT: attach this version's output as an input and Save & Run All again.")
"""

CELLS = [
    ("md", INTRO),
    ("md", "## 0. Settings"),
    ("code", SETTINGS),
    ("md", "## 1. Setup\n\nClone, restore the newest earlier session, the inputs, the frozen settings and the one RL "
           "platform (a session on other hardware refuses to run)."),
    ("code", SETUP),
    ("md", "## 2. RL tuning\n\n12 configurations per method at 10 % of the budget, seed 0, scored on DEV48 in the full "
           "simulator (`v5_plan.md` §4.5)."),
    ("code", TUNE),
    ("md", "## 3. RL training\n\n`dagger_fnpeg`, `bcppo_fnpeg`, `ppo_hist`, seeds 0–4, each to its frozen budget; "
           "then each policy's fast-environment success."),
    ("code", TRAIN),
    ("md", "## 4. Save and summary"),
    ("code", SUMMARY),
]


def build() -> str:
    nb = nbf.v4.new_notebook()
    nb["cells"] = [nbf.v4.new_markdown_cell(t.strip("\n")) if k == "md" else nbf.v4.new_code_cell(t.strip("\n"))
                   for k, t in CELLS]
    for i, c in enumerate(nb["cells"]):
        c["id"] = f"v5k-cell-{i:02d}"
    nb["metadata"] = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                      "language_info": {"name": "python"},
                      "kaggle": {"accelerator": "nvidiaTeslaT4", "isInternetEnabled": True, "isGpuEnabled": True}}
    path = os.path.join(SIM, NAME)
    nbf.validate(nb)
    nbf.write(nb, path)
    print(f"wrote {path} ({len(CELLS)} cells)")
    return path


if __name__ == "__main__":
    build()
