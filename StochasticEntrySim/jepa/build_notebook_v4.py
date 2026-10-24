"""
Build ReEntryAI_JEPA_v4.ipynb, the one notebook of the v4 study (jepa/v4_plan.md section 9).

All logic lives in jepa/ (study_v4, fly_v4, physics_v4, ablate_v4, analyze_v4, report_v4); the
notebook configures, calls, and shows. Edit this builder, never the .ipynb:

    python jepa/build_notebook_v4.py          (from StochasticEntrySim/)

The auto-disconnect code is the v3 builder's (repo rule, CLAUDE.md), imported, not copied.
Smoke mode (every stage, tiny, numbers meaningless): JEPA_V4_MODE=smoke, with JEPA_V4_ROOT,
JEPA_V3_ROOT and JEPA_WIDE_DATA set (jepa/smoke_v4.py does this).
"""
import os
import sys

import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from build_notebooks_v3 import AUTO_DISCONNECT_CODE, DONE, DONE_MD  # noqa: E402

NAME = "ReEntryAI_JEPA_v4.ipynb"
BRANCH = "space-jepa4"
BADGE = ("[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)]"
         f"(https://colab.research.google.com/github/rk0604/ReEntryAI/blob/{BRANCH}/StochasticEntrySim/{NAME})")

# The Settings cell owns AUTO_DISCONNECT (so it can be switched off before anything runs); the v3 code sets it
# itself, so that one line is made to keep the setting.
_V3_LINE = "AUTO_DISCONNECT = IN_COLAB      # set to False to keep the runtime (e.g. to explore results by hand)"
assert _V3_LINE in AUTO_DISCONNECT_CODE, "build_notebooks_v3.AUTO_DISCONNECT_CODE changed: update this builder"
AUTO_DISCONNECT_V4 = AUTO_DISCONNECT_CODE.replace(
    _V3_LINE, "AUTO_DISCONNECT = AUTO_DISCONNECT   # set in the Settings cell (section 0)")

INTRO = f"""
# SpaceJEPA v4 · Why is v2 better?

{BADGE}

**How to run.** Open in Colab (badge above), choose an **L4 GPU** runtime (it opens that way), and
**Runtime → Run all**. That is all. Colab can disconnect when idle or when its session limit ends;
**after any disconnect, Run all again**: every stage reads Drive, skips finished work in seconds, and
resumes what is half done. No cell ever needs editing to resume. The notebook **disconnects itself**
when the last cell finishes or a cell fails (so it never idles on paid compute); set
`AUTO_DISCONNECT = False` in the Settings cell to keep the runtime. For the analysis and report only,
you can switch to a cheap **CPU runtime**: the GPU stages then print what remains and skip.

v2 beats the classical predictor-corrector by 54 points of success within 5 km (v3), but it differs
from it in two things at once (the MPPI planner and the learned model), and from v1 in about eight
design changes. This study attributes the gains (`jepa/v4_plan.md`):

| | Question | Experiment |
|---|---|---|
| RQ1 | How much of v2's advantage is the planner, how much the model? Is v2's model close to a perfect one? | E1: MPPI through the classical controller's own physics (nominal and true parameters) and through v2's physics prior; the classical controller knowing the true density and mass (not L/D: its predictor resets the L/D scale, see v4_plan.md §14) |
| RQ2 | Does v2 still beat classical when it pulls no more load? | E2: v2 at 4.0 g, 3.5 g and at classical's own peak g per condition |
| RQ3 | Which v1 → v2 changes are necessary, which sufficient? | E3: 13 single-change variants of a one-member v2 (R0) and of v1, 5 seeds each |
| RQ4 | Why does v2 fail case 10? | E4: 40 diagnostic flights with the planner's decision log (descriptive) |

Everything is paired with v3's 5,120 flights: same 64 cases, 4 atmospheres, seeds, MPPI settings and
noise streams. Inputs (read-only): `MyDrive/ReEntryAI_jepa_v3/` and
`MyDrive/ReEntryAI_jepa_data/wide_6000/`. Output: **`MyDrive/ReEntryAI_jepa_v4/`** (layout in
`v4_plan.md` §3 and the auto-written `00_study/README.md`). v4 refuses to write anywhere else.

**Rough cost on an L4** (from v3's measured rates; the status cell measures the real ones): E0 ≈ 1 h,
E1 ≈ 4–7 h (the classical oracle is ≈ 6 CPU-hours), E2 ≈ 3–5 h, E4 ≈ 5 min, E3 training ≈ 17–20 h (70
runs), E3 flights ≈ 15–25 h, open loop ≈ 1–2 h, analysis ≈ 15 min. Flights are limited by CPU cores
(the simulator), not the GPU. It spans many sessions; that is what the resume rule is for.
"""

SETTINGS = r"""
import os, sys
IN_COLAB = "google.colab" in sys.modules
MODE = os.environ.get("JEPA_V4_MODE", "full")     # "full": the study | "smoke": every stage, tiny, numbers meaningless
RUN = {"E0": True, "E1": True, "E2": True, "E4": True, "E3_train": True, "E3_fly": True,   # stages to run this session
       "open_loop": True, "analysis": True, "report": True}
AUTO_DISCONNECT = IN_COLAB      # release the runtime when the notebook finishes or a cell fails; False keeps it
CPU_WORKERS = None              # classical flights, one simulator per core (None = every core)
GPU_WORKERS = None              # learned and physics flights sharing the GPU (None = 6 GB each, 2 GB spare, at most 4)
                                # (on a CPU, i.e. smoke mode: half the cores, at most 4, as v3's smoke run)
"""

SETUP = r"""
BRANCH = "space-jepa4"
if IN_COLAB:
    REPO = "/content/ReEntryAI"
    if not os.path.exists(REPO):
        !git clone -q -b {BRANCH} https://github.com/rk0604/ReEntryAI.git {REPO}
    else:
        !git -C {REPO} pull -q
    %pip install -q pymsis reportlab
    os.chdir(f"{REPO}/StochasticEntrySim")
    from google.colab import drive
    drive.mount("/content/drive")
    DRIVE = "/content/drive/MyDrive"
    V4_ROOT = f"{DRIVE}/ReEntryAI_jepa_v4" + ("" if MODE == "full" else f"_{MODE}")
    V3_ROOT = f"{DRIVE}/ReEntryAI_jepa_v3" + ("" if MODE == "full" else f"_{MODE}")
    DATA_CANDIDATES = [f"{DRIVE}/ReEntryAI_jepa_data/wide_6000/episodes_wide_6000.pkl",
                       f"{DRIVE}/ReEntryAI_jepa_data/wide_6000"]
    LOCAL_DATA = "/content/data"
else:
    V4_ROOT = os.environ.get("JEPA_V4_ROOT", f"jepa_runs/v4_{MODE}")
    V3_ROOT = os.environ.get("JEPA_V3_ROOT", f"jepa_runs/v3_{MODE}")
    DATA_CANDIDATES = [os.environ.get("JEPA_WIDE_DATA", ""), "jepa_runs/episodes_wide_6000.pkl", "jepa_runs/wide_6000"]
    LOCAL_DATA = "jepa_runs"

import json, time
from pathlib import Path
SIM_DIR = Path.cwd()
assert (SIM_DIR / "rl_env.py").exists(), "start Jupyter in StochasticEntrySim/"
sys.path.insert(0, str(SIM_DIR))
import numpy as np, pandas as pd, torch
import matplotlib.pyplot as plt
from IPython.display import Image, display
%matplotlib inline
from jepa import study_v4 as S4, fly_v4 as F4, ablate_v4 as AB, analyze_v4 as A4, report_v4 as R4, train_v2 as T2

CFG = S4.StudyConfigV4.smoke() if MODE == "smoke" else S4.StudyConfigV4()
L3 = S4.v3_layout(V3_ROOT)                                            # read-only v3 folder (raises if missing)
DATA_DIRS = [c if os.path.isdir(c) else os.path.dirname(c) for c in DATA_CANDIDATES if c and os.path.exists(c)]
L = S4.LayoutV4(V4_ROOT, read_only=[V3_ROOT] + DATA_DIRS).make()      # every v4 write is checked against these
FROZEN = S4.freeze_config_v4(CFG, L, L3)                               # raises if v3's or v4's settings differ
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
CPUS = os.cpu_count() or 2
GPU_MEM_GB = torch.cuda.get_device_properties(0).total_memory / 1e9 if DEVICE == "cuda" else 0.0
CPU_WORKERS = CPU_WORKERS or CPUS
GPU_WORKERS = GPU_WORKERS or (max(1, min(4, CPUS, int((GPU_MEM_GB - 2) // 6))) if DEVICE == "cuda"
                              else max(1, min(4, CPUS // 2)))          # CPU (smoke): as v3's smoke run
CAN_GPU = DEVICE == "cuda" or MODE == "smoke"      # GPU stages; in smoke mode the CPU stands in for the GPU
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")   # inherited by the GPU workers
if DEVICE == "cuda":
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True
GRAD_CKPT = DEVICE == "cuda" and GPU_MEM_GB < 20   # T4: recompute activations to fit T = 96 (as v3)
torch.set_num_threads(1 if DEVICE == "cuda" else max(1, CPUS // 2))
DATA = None                                        # the 6,000 flights, loaded only if a stage needs them

def need_data():
    # Load the 6,000 flights once (copied to local disk first, as v3 does).
    global DATA, DATA_PATH
    if DATA is None:
        t0 = time.time()
        DATA_PATH = T2.local_copy(T2.find_wide_dataset(DATA_CANDIDATES, download_dir=LOCAL_DATA), LOCAL_DATA)
        DATA = AB.load_data_v4(CFG, DATA_PATH, DEVICE)
        print(f"{DATA['n_episodes']} flights: train/val/test = {len(DATA['train'])}/{len(DATA['val'])}/"
              f"{len(DATA['test'])} episodes ({time.time() - t0:.0f}s)")
    return DATA

print(f"MODE={MODE}  v4 folder={V4_ROOT}  v3 folder={V3_ROOT}  device={DEVICE}"
      + (f" ({torch.cuda.get_device_name(0)}, {GPU_MEM_GB:.0f} GB)" if DEVICE == "cuda" else "")
      + f"  cpus={CPUS}  workers: {CPU_WORKERS} CPU, {GPU_WORKERS} GPU  git={FROZEN['git_commit'][:10]}")
if MODE == "full" and DEVICE != "cuda":
    print("No GPU: the flight, training and open-loop stages will only report what remains. "
          "Runtime -> Change runtime type -> L4 for those.")
"""

STATUS = r"""
JOBS = S4.build_jobs(CFG, L, L3)
ST = S4.status(CFG, L, JOBS)
display(pd.DataFrame(ST["stages"]))
arms = pd.DataFrame(ST["arms"])
display(arms[arms.done < arms.total] if (arms.done < arms.total).any() else arms.groupby("stage")[["done", "total"]].sum())
tr = pd.DataFrame(ST["training"])
print(f"{ST['flights_done']}/{ST['flights_total']} flights, {ST['models_done']}/{ST['models_total']} models trained"
      + (f", ~{ST['training_hours_left']} h of training left" if ST["training_hours_left"] else "")
      + ("  (E0b failed: v3's learned arms are re-flown, see 01_checks/)" if ST["e0b_failed"] else ""))
"""

E0 = r"""
if RUN["E0"]:
    e0_jobs = [j for j in JOBS if j["stage"] == "E0"]
    S4.run_classical_v4(e0_jobs, L, workers=CPU_WORKERS)                    # E0a: 32 classical re-flights
    E0C_OK = S4.e0c_run(CFG, L, L3)                                         # E0c: runner equivalence, CPU
    if CAN_GPU:
        S4.run_flights_v4(CFG, e0_jobs, L, device=DEVICE, workers=GPU_WORKERS)   # E0b: v3 flights re-flown
    else:
        print("E0b waits for a GPU runtime")
    PHYS = S4.physics_checks(CFG, L, L3)                                     # section 4.1.3 checks 1 and 2
    E0 = S4.e0_checks(CFG, L, L3)
    if E0["E0b"] is False:                                                   # contingency: re-fly v3's learned arms
        JOBS = S4.build_jobs(CFG, L, L3)
        if CAN_GPU:
            S4.run_flights_v4(CFG, [j for j in JOBS if j["stage"] == "E0"], L, device=DEVICE, workers=GPU_WORKERS)
    for f in ("E0a.csv", "E0c.csv"):
        if (L.checks / f).exists():
            display(pd.read_csv(L.checks / f).head(12))
    if (L.checks / "E0b_summary.json").exists():
        display(pd.DataFrame(json.loads((L.checks / "E0b_summary.json").read_text())["pairs"]))
    print("physics:", json.dumps(PHYS, indent=1)[:1500])
    print("E0:", E0)
    if E0["E0a"] is False:
        raise RuntimeError("E0a FAILED: the classical controller no longer reproduces v3's flights (01_checks/E0a.csv). "
                           "Stop and report: every v4-v3 comparison assumes the same simulator and guidance.")
    if E0C_OK is False:
        raise RuntimeError("E0c FAILED: run_lockstep_v4 differs from run_lockstep_v3 on the CPU (01_checks/E0c.csv).")
"""

FLY_E124 = r"""
todo = [j for j in JOBS if j["stage"] in ("E1", "E2", "E4") and RUN[j["stage"]]]
S4.run_classical_v4(todo, L, workers=CPU_WORKERS)                           # classical_oracle (CPU pool)
if CAN_GPU:
    S4.run_flights_v4(CFG, todo, L, device=DEVICE, workers=GPU_WORKERS)     # physics and v2 arms, seed-major
else:
    print(f"{sum(not S4.flight_done(L, j) for j in todo if j['model']['kind'] != 'classical')} GPU flights wait "
          "for a GPU runtime")
"""

TRAIN = r"""
if RUN["E3_train"]:
    todo = [(v, s) for s in CFG.seeds for v in CFG.variants if not L.ckpt(v, s).exists()]   # seed-major
    print(f"{len(todo)} of {len(CFG.seeds) * len(CFG.variants)} trainings to do")
    if todo and CAN_GPU:
        need_data()
        for v, s in todo:
            t0 = time.time()
            r = AB.train_variant_seed(CFG, DATA, L, v, s, DEVICE, grad_ckpt=GRAD_CKPT)
            print(f"== {v} seed {s}: {'cached' if r['cached'] else f'trained in {(time.time() - t0) / 60:.1f} min'}")
    elif todo:
        print("training waits for a GPU runtime")
    p = L.model_dir("R0", CFG.seeds[0]) / "F_curves.png"
    if p.exists():
        display(Image(str(p)))
"""

FLY_E3 = r"""
if RUN["E3_fly"]:
    todo = [j for j in JOBS if j["stage"] == "E3"]
    if CAN_GPU:
        S4.run_flights_v4(CFG, todo, L, device=DEVICE, workers=GPU_WORKERS)
    else:
        print(f"{sum(not S4.flight_done(L, j) for j in todo)} E3 flights wait for a GPU runtime")
"""

OPEN_LOOP = r"""
if RUN["open_loop"]:
    models = list(CFG.variants) + ["v1"]
    ckpt = lambda m, s: L3.ckpt("v1", s) if m == "v1" else L.ckpt(m, s)
    todo = [(m, s) for s in CFG.seeds for m in models
            if ckpt(m, s).exists() and not (L.open_loop / f"footprint_{m}_seed{s}.csv").exists()]
    print(f"{len(todo)} model-seeds to evaluate open loop")
    if todo and CAN_GPU:
        eps = {"val": need_data()["val"], "test": DATA["test"]}
        for m, s in todo:
            t0 = time.time()
            A4.open_loop_seed(CFG, L, L3, eps, m, s, device=DEVICE)
            print(f"  {m} seed {s}: {time.time() - t0:.0f}s")
    elif todo:
        print("open loop waits for a GPU runtime")
    ol = A4.open_loop_table(L)
    if len(ol):
        display(ol[(ol.split == "test")].pivot_table(index="model", columns="ttg_s", values="err_median_km").round(2))
"""

ANALYSIS = r"""
if RUN["analysis"]:
    t0 = time.time()
    R = A4.analyze(CFG, L, L3)
    pd.set_option("display.max_colwidth", 140)
    display(R["T05"][["id", "question", "measured", "p_holm", "met"]])
    for k in ("F10", "F11", "F12", "F13", "F17"):
        if k in R["figures"]:
            display(Image(str(R["figures"][k])))
    print(f"analysis in {time.time() - t0:.0f}s")
"""

REPORT = r"""
S4.status(CFG, L, S4.build_jobs(CFG, L, L3))
if RUN["report"]:
    print("report:", R4.build_report(L))
S4.update_manifest(L)
print(open(L.study / "README.md").read()[:3000])
"""

CELLS = [
    ("md", INTRO),
    ("md", r"""
## 0. Settings

`MODE`, the stages to run this session (`RUN`), auto-disconnect, and worker counts (automatic unless you
set them). Nothing else needs editing, ever: resuming after a disconnect is Run all.
"""),
    ("code", SETTINGS),
    ("md", r"""
## 1. Setup

Clones (or pulls) the `space-jepa4` branch, installs `pymsis` and `reportlab`, mounts Drive, checks that
the v3 folder exists and was run with the settings v4 depends on (frozen seeds, test set, MPPI, budgets),
and freezes `00_study/study_config.json` (a later session with different settings raises instead of mixing
results). Registers the auto-disconnect hook (CLAUDE.md). Reads `00_study/` of v3; writes `00_study/` of v4.
Takes about a minute.
"""),
    ("code", SETUP + AUTO_DISCONNECT_V4),
    ("md", r"""
## 2. Status

Every stage and arm: flights done / total, the measured rate, hours left, and each model's training state.
Counts result files directory by directory (fast, no tree walk). Writes `00_study/status.json`. Seconds.
"""),
    ("code", STATUS),
    ("md", r"""
## 3. E0: reproducibility checks (gate everything else)

* **E0a**: v3's classical flights re-flown on 32 conditions (CPU pool): must be identical (outcome, miss
  within 1e-6 km, fuel), or the notebook **stops**. Also compares the classical oracle with classical in
  the nominal atmosphere once E1 has flown it.
* **E0b**: v3's `v2single` and `v1` (seed 0, 4.5 g, 256 conditions each) re-flown by the v4 runner from
  v3's checkpoints: outcome agreement ≥ 98% and the paired log-miss difference's 95% CI must include 0.
  Its spread is the GPU-nondeterminism floor of every v4-v3 comparison. If it fails, v3's learned arms
  are re-flown with the v4 runner (3,840 flights) and used instead; record that in `v4_plan.md` §15.
* **E0c**: the v4 runner equals v3's bit for bit on the CPU (3 flights, v1 and v2 member 0).
* **Physics fidelity**: the torch port of the classical model against `math_3d` (1,000 states), and its
  open-loop landing error on 64 recorded v3 flights (300, 200, 100 s to go).

Reads v3's flights and checkpoints; writes `01_checks/` and `04_flights/` of v4. About 1 h (E0b ≈ 30–50
min). Finished checks and flights are skipped.
"""),
    ("code", E0),
    ("md", r"""
## 4. E1, E2, E4: flights that need no training

* **E1** planner × model: `classical_oracle` (CPU pool, 256), and MPPI through the classical controller's
  own physics model with nominal (`mppi_nominal`, at 4.5 g and at the matched load) and true parameters
  (`mppi_truth`), and through v2's physics prior with no learning (`mppi_prior`): 5,120 flights.
* **E2** load: v2 (v3's checkpoints) at 4.0 g, 3.5 g, and the matched load (classical's measured peak g on
  each condition): 3,840 flights.
* **E4** case 10: v2 and v2single, 4.5 g and 8 g, nominal and dense air, 5 seeds, with the decision log: 40.

Run order E1 → E2 → E4, seed-major. Each flight's JSON and telemetry (`.npz`, with the decision log) is
written to `04_flights/` the moment it lands; finished flights are skipped and flights cut by a disconnect
are flown again. About 7–12 h on an L4.
"""),
    ("code", FLY_E124),
    ("md", r"""
## 5. E3 training: R0 and 13 variants × 5 seeds

R0 is v2's recipe with one member; A1–A8 remove one v2 change from it; B1–B3 add one v2 change to v1 (v1's
recipe and v3's fair budget). Seed-major. Each run writes `02_models/<variant>/seed_<k>/`: `variant.json`
(the exact configuration and parameter count), the train state every 5 min (model, optimizer, scaler, step,
RNG state, loss history; a rerun resumes there), then `ckpt.pt` (marks it done; the train state is deleted),
`train_log.csv`, `latent_health.csv`, `F_curves.png`. Loads the 6,000 flights only if something is left to
train. About 17–20 h on an L4 for all 70 (v2-family ≈ 15–20 min each, B ≈ 5–10 min).
"""),
    ("code", TRAIN),
    ("md", r"""
## 6. E3 flights: the 14 models × 5 seeds at 4.5 g

17,920 flights through the same lockstep runner, CRN and settings as v3. A model whose checkpoint does not
exist yet is skipped with a note (train it in section 5 first). About 15–25 h on an L4.
"""),
    ("code", FLY_E3),
    ("md", r"""
## 7. Open loop

Per model (R0, the 13 variants, and v3's v1) and seed: predicted vs true drogue point on the 600 held-out
test and 600 validation flights, from 400 to 50 s to go, with the commands actually flown
(`03_open_loop/footprint_<model>_seed<k>.csv`, along-track bias by density bin). Cached per file. Needs a
GPU (≈ 1 min per model-seed on an L4).
"""),
    ("code", OPEN_LOOP),
    ("md", r"""
## 8. Analysis

v3 and v4 flights together: every table `05_analysis/tables/T01–T13`, figure `figures/F10–F18`,
`captions.json`, and the pre-registered criteria (§8: P1–P6, L1–L3, each A and B variant, R0-check, M),
with the hierarchical bootstrap (10,000 resamples) and Holm within each family. Bootstrap results are
cached in `05_analysis/cache/` by comparison and a hash of their input flights, so a rerun after new
flights recomputes only what changed. Runs on any runtime at any point; with a partly finished study it
reports what exists. ≈ 10–20 min the first time.
"""),
    ("code", ANALYSIS),
    ("md", r"""
## 9. Report

`06_report/SpaceJEPA_v4_report.pdf` (summary, RQ1–RQ4, methods, deviations from `v4_plan.md` §15, every
table), the status file, and a refreshed `00_study/manifest.csv` and `README.md`. Rebuilt every run (cheap).
"""),
    ("code", REPORT),
    ("md", DONE_MD),
    ("code", DONE),
]


def build() -> str:
    nb = nbf.v4.new_notebook()
    nb["cells"] = [nbf.v4.new_markdown_cell(t.strip("\n")) if k == "md" else nbf.v4.new_code_cell(t.strip("\n"))
                   for k, t in CELLS]
    for i, c in enumerate(nb["cells"]):
        c["id"] = f"v4-cell-{i:02d}"                 # stable ids: rebuilding gives the same file
    nb["metadata"] = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                      "language_info": {"name": "python"}, "accelerator": "GPU",
                      "colab": {"gpuType": "L4", "provenance": []}}
    path = os.path.join(SIM, NAME)
    nbf.validate(nb)
    nbf.write(nb, path)
    print(f"wrote {path} ({len(CELLS)} cells)")
    return path


if __name__ == "__main__":
    build()
