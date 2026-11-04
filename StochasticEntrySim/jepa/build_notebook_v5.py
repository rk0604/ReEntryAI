"""
Build ReEntryAI_JEPA_v5.ipynb, the one notebook of the v5 study (jepa/v5_plan.md section 9).

All logic lives in jepa/ (benchmark_v5, guidance_fnpeg, guidance_predguid, fly_v5, models_v5, rlenv_fast_v5, rl_v5,
study_v5, analyze_v5, report_v5); the notebook configures, calls, and shows. Edit this builder, never the .ipynb:

    python jepa/build_notebook_v5.py          (from StochasticEntrySim/)

The auto-disconnect code is the v3 builder's (repo rule, CLAUDE.md), imported, not copied. Smoke mode (every stage,
tiny, numbers meaningless): JEPA_V5_MODE=smoke with JEPA_V5_ROOT, JEPA_V3_ROOT, JEPA_V4_ROOT, JEPA_WIDE_DATA and
JEPA_V5_LEGACY set (jepa/smoke_v5.py does this). JEPA_V5_FAIL_AT=<stage>:<n> simulates a disconnect (smoke check 14).
"""
import os
import sys

import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from build_notebooks_v3 import AUTO_DISCONNECT_CODE, DONE, DONE_MD  # noqa: E402

NAME = "ReEntryAI_JEPA_v5.ipynb"
BRANCH = "space-jepa5"
BADGE = ("[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)]"
         f"(https://colab.research.google.com/github/rk0604/ReEntryAI/blob/{BRANCH}/StochasticEntrySim/{NAME})")

_V3_LINE = "AUTO_DISCONNECT = IN_COLAB      # set to False to keep the runtime (e.g. to explore results by hand)"
assert _V3_LINE in AUTO_DISCONNECT_CODE, "build_notebooks_v3.AUTO_DISCONNECT_CODE changed: update this builder"
AUTO_DISCONNECT_V5 = AUTO_DISCONNECT_CODE.replace(
    _V3_LINE, "AUTO_DISCONNECT = AUTO_DISCONNECT   # set in the Settings cell (section 0)")

INTRO = f"""
# SpaceJEPA v5 · Strong baselines and a standard benchmark

{BADGE}

**How to run.** Open in Colab (badge above), choose an **L4 GPU** runtime (it opens that way), and
**Runtime → Run all**. **After any disconnect, crash or runtime release, Run all again**: every stage
reads Drive, skips finished work in seconds and resumes half-done work from its last checkpoint (at most
5 minutes of any training is lost). No cell ever needs editing to resume; Colab may only ask you to
re-authorize the Drive mount. The notebook **disconnects itself** when the last cell finishes or a cell
fails (so it never idles on paid compute units); set `AUTO_DISCONNECT = False` in the Settings cell to
keep the runtime. CPU-only stages (gates, tuning, the guidance laws, classical ORB-1 flights) are limited
by CPU cores: for long CPU stretches a **high-CPU CPU runtime** is cheaper; GPU stages then print what
remains and skip.

**Tiers (decision D1).** `TIER = 1` (default) runs the gates, legacy RL re-evaluation, tuning, the A256
arms (FNPEG, PredGuid, MPPI + filtered physics, v2 at FNPEG's load), ORB-1 block 1, the four non-JEPA
models and three RL baselines. `TIER = 2` (after reading Tier 1's results) adds ORB-1 blocks 2–3 for the
arms still eligible, `mppi_truth` and the best non-JEPA and RL arms on ORB-1, SAC, PPO in the full
simulator, and the exploratory arms. Everything is built and smoke-tested either way.

| | Question | Experiment |
|---|---|---|
| RQ1 | Does v2 still beat classical guidance when the classical side is a published, tuned, load-limited, model-adapting law, and does the learned model beat explicit online estimation under the same planner? | E1: FNPEG (Lu 2014) and PredGuid (Bairstow 2006); E2: MPPI + filtered physics, v2 at FNPEG's own peak g |
| RQ2 | Is the JEPA objective needed? | E3: MLP, GRU, GRU + PI prober, R0 without JEPA, same data and planner |
| RQ3 | How close do properly built RL baselines get? | E4: DAgger and BC → PPO from FNPEG, PPO (and SAC, full-sim PPO in Tier 2) |
| RQ4 | On a standard Orion-class scenario, which controllers meet an Orion-style accuracy requirement and crew g limits? | ORB-1: 3 × 1,000 dispersed LEO returns, sequential stopping |
| RQ5 | What does each controller cost per decision? | timing table |

Inputs (read-only, guarded): `MyDrive/ReEntryAI_jepa_v3/`, `MyDrive/ReEntryAI_jepa_v4/`,
`MyDrive/ReEntryAI_jepa_data/wide_6000/`, `MyDrive/ReEntryAI/` (legacy RL). Output:
**`MyDrive/ReEntryAI_jepa_v5/`** (layout in `jepa/v5_plan.md` §9 and the auto-written
`00_study/README.md`); stage logs in `00_study/logs/`, progress in `00_study/status.json` and
`status_history.csv`. v5 refuses to write anywhere else.
"""

SETTINGS = r"""
import os, sys
IN_COLAB = "google.colab" in sys.modules
MODE = os.environ.get("JEPA_V5_MODE", "full")     # "full": the study | "smoke": every stage, tiny, numbers meaningless
TIER = int(os.environ.get("JEPA_V5_TIER", "1"))   # D1: 1 = Tier 1 (default); 2 = after reading Tier 1's results
RUN_OVERRIDE = {}               # e.g. {"ORB1_block1": False} to skip a stage this session (RUN is derived from TIER)
AUTO_DISCONNECT = IN_COLAB      # release the runtime when the notebook finishes or a cell fails; False keeps it
CPU_WORKERS = None              # guidance-law, policy and classical flights, one simulator per core (None = every core)
GPU_WORKERS = None              # MPPI flights sharing the GPU (None = 6 GB each, 2 GB spare, at most 4)
"""

SETUP = r"""
BRANCH = "space-jepa5"
if IN_COLAB:
    REPO = "/content/ReEntryAI_v5"
    if not os.path.exists(REPO):
        !git clone -q -b {BRANCH} https://github.com/rk0604/ReEntryAI.git {REPO}
    else:
        !git -C {REPO} pull -q origin {BRANCH}
    %pip install -q pymsis reportlab
    os.chdir(f"{REPO}/StochasticEntrySim")
    from google.colab import drive
    drive.mount("/content/drive")
    DRIVE = "/content/drive/MyDrive"
    SFX = "" if MODE == "full" else f"_{MODE}"
    V5_ROOT, V4_ROOT, V3_ROOT = (f"{DRIVE}/ReEntryAI_jepa_v{k}{SFX}" for k in (5, 4, 3))
    LEGACY_ROOT = f"{DRIVE}/ReEntryAI"
    DATA_CANDIDATES = [f"{DRIVE}/ReEntryAI_jepa_data/wide_6000/episodes_wide_6000.pkl",
                       f"{DRIVE}/ReEntryAI_jepa_data/wide_6000"]
    LOCAL_DATA = "/content/data"
else:
    V5_ROOT = os.environ.get("JEPA_V5_ROOT", f"jepa_runs/v5_{MODE}")
    V4_ROOT = os.environ.get("JEPA_V4_ROOT", f"jepa_runs/v4_{MODE}")
    V3_ROOT = os.environ.get("JEPA_V3_ROOT", f"jepa_runs/v3_{MODE}")
    LEGACY_ROOT = os.environ.get("JEPA_V5_LEGACY", "")
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
from jepa import study_v4 as S4, study_v5 as S5, analyze_v5 as A5, report_v5 as R5, train_v2 as T2

CFG = S5.StudyConfigV5.smoke() if MODE == "smoke" else S5.StudyConfigV5()
CFG.tier = TIER
RUN = {**CFG.stages(), **RUN_OVERRIDE}
L3 = S4.v3_layout(V3_ROOT)                                             # read-only v3 folder (raises if missing)
L4 = S5.v4_layout(V4_ROOT)                                             # read-only v4 folder (None: RQ2's R0 waits)
DATA_DIRS = [c if os.path.isdir(c) else os.path.dirname(c) for c in DATA_CANDIDATES if c and os.path.exists(c)]
L = S5.LayoutV5(V5_ROOT, read_only=[V3_ROOT, V4_ROOT, LEGACY_ROOT] + DATA_DIRS).make()   # every write is checked
LOG = L.logger("setup")
REPAIR = S5.repair(L, LOG)                                             # stale .tmp, unreadable flights, broken .done
FROZEN = S5.freeze_config_v5(CFG, L, inputs={"v3": V3_ROOT, "v4": V4_ROOT if L4 else None,
                                             "legacy": LEGACY_ROOT or None})
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
RUNTIME = torch.cuda.get_device_name(0) if DEVICE == "cuda" else "cpu"
CPUS = os.cpu_count() or 2
GPU_MEM_GB = torch.cuda.get_device_properties(0).total_memory / 1e9 if DEVICE == "cuda" else 0.0
CPU_WORKERS = CPU_WORKERS or (CPUS if MODE == "full" else min(2, CPUS))
GPU_WORKERS = GPU_WORKERS or (max(1, min(4, CPUS, int((GPU_MEM_GB - 2) // 6))) if DEVICE == "cuda" else 1)
CAN_GPU = DEVICE == "cuda" or MODE == "smoke"       # GPU stages; in smoke mode the CPU stands in for the GPU
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
if DEVICE == "cuda":
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True
torch.set_num_threads(1 if DEVICE == "cuda" else max(1, CPUS // 2))
DATA = None

def need_data():
    # The 6,000 flights, loaded once (copied to local disk first).
    global DATA, DATA_PATH
    if DATA is None:
        t0 = time.time()
        DATA_PATH = T2.local_copy(T2.find_wide_dataset(DATA_CANDIDATES, download_dir=LOCAL_DATA), LOCAL_DATA)
        DATA = S5.load_data_v5(CFG, DATA_PATH, DEVICE)
        LOG(f"{DATA['n_episodes']} flights: train/val/test = {len(DATA['train'])}/{len(DATA['val'])}/"
            f"{len(DATA['test'])} episodes ({time.time() - t0:.0f}s)")
    return DATA

LOG(f"MODE={MODE} TIER={TIER} v5 folder={V5_ROOT} device={DEVICE} ({RUNTIME}) cpus={CPUS} workers: {CPU_WORKERS} "
    f"CPU, {GPU_WORKERS} GPU; git={FROZEN['git_commit'][:10]}; v4 folder {'found' if L4 else 'NOT found'}; "
    f"repair: {REPAIR}")
LOG("RUN this session: " + ", ".join(k for k, v in RUN.items() if v))
if MODE == "full" and DEVICE != "cuda":
    LOG("No GPU: training, MPPI flights and open loop only report what remains (Runtime -> Change runtime type -> L4).")
"""

STATUS = r"""
ST = S5.status(CFG, L, S5.all_stage_jobs(CFG, L, L3), RUNTIME)
display(pd.DataFrame(ST["stages"]))
arms = pd.DataFrame(ST["arms"])
if len(arms):
    display(arms[arms.done < arms.total] if (arms.done < arms.total).any() else arms.groupby("stage")[["done", "total"]].sum())
print(f"{ST['flights_done']}/{ST['flights_total']} flights, {ST['trainings_done']}/{ST['trainings_total']} trainings; "
      f"~{ST['flight_hours_left']} flight-hours left (~{ST['flight_units_left']} units at {ST['units_per_hour']}/h); "
      f"rates (flights/h) by arm type: {ST['rates_by_type']}")
"""

BENCH = r"""
if RUN["benchmark"]:
    log = L.logger("benchmark")
    S0 = S5.gate_s0(CFG, L, None if (L.checks / "S0.done").exists() else need_data(), log)
    log(f"[benchmark] gate S0 pass={S0.get('pass')}, checksums={S0.get('checksums')}, "
        f"targets inside the footprint={S0.get('targets_inside')}, "
        f"feasible at 4.5 g={S0.get('feasible_4p5g_frac')} (01_checks/S0.json)")
    print(json.dumps({k: v for k, v in S0.items() if k not in ("checksum_error",)}, indent=1, default=str)[:2500])
    if S0.get("pass") is False and CFG.gates_stop:
        raise RuntimeError("gate S0 FAILED (01_checks/S0.json): the benchmark is not what the plan froze. Stop and report.")
"""

LEGACY = r"""
if RUN["legacy"]:
    log = L.logger("legacy")
    FINALS = S5.legacy_inventory(CFG, L, LEGACY_ROOT if LEGACY_ROOT and os.path.isdir(LEGACY_ROOT) else None, log)
    if FINALS:
        try:
            import stable_baselines3
        except ImportError:
            if IN_COLAB:
                %pip install -q stable-baselines3
            else:
                raise
        S5.run_cpu(CFG, S5.legacy_jobs(CFG, L, FINALS), L, CPU_WORKERS, log, "legacy RL")
    display(pd.read_csv(L.study / "legacy_rl_inventory.csv").head(30))
"""

GATES = r"""
if RUN["gates"]:
    log = L.logger("gates")
    GATES = S5.gates(CFG, L, L3, CPU_WORKERS, log)
    display(pd.DataFrame([{"gate": k, "pass": v.get("pass"), **{x: y for x, y in v.items()
                           if isinstance(y, (int, float, bool)) and x != "pass"}} for k, v in GATES.items()]))
    failed = [k for k, v in GATES.items() if v.get("pass") is False]
    if failed and CFG.gates_stop:
        raise RuntimeError(f"gate(s) {failed} FAILED (01_checks/): no test flight may run on a failed gate (v5_plan.md "
                           "section 4.0). Stop and report.")
"""

TUNING = r"""
if RUN["tuning"]:
    log = L.logger("tuning")
    TUNED = S5.tune_guidance(CFG, L, CPU_WORKERS, log)
    RL_ARMS = list(S5.RL_TIER1) + (list(S5.RL_TIER2) if CFG.tier >= 2 else [])
    if CAN_GPU:
        S5.tune_rl(CFG, L, RL_ARMS, DEVICE, CPU_WORKERS, log)
    else:
        log("RL tuning waits for a GPU runtime")
    for m in ("fnpeg", "predguid") + tuple(f"rl_{a}" for a in RL_ARMS):
        p = L.tuning / m / "chosen.json"
        print(m, "->", json.loads(p.read_text())["score"] if p.exists() else "not chosen yet")
"""

A256 = r"""
if RUN["A256"]:
    log = L.logger("A256")
    S5.run_jobs(CFG, S5.jobs_a256(CFG, L, L3, matched=False), L, DEVICE, CPU_WORKERS, GPU_WORKERS, CAN_GPU, log, "A256")
    S5.run_jobs(CFG, S5.jobs_matched(CFG, L, L3), L, DEVICE, CPU_WORKERS, GPU_WORKERS, CAN_GPU, log, "A256 matched")
"""

ORB1 = r"""
if RUN["ORB1_block1"]:
    log = L.logger("ORB1")
    S5.run_jobs(CFG, S5.orb1_jobs(CFG, L, L3, S5.ORB1_ARMS, [0], "ORB1", log), L, DEVICE, CPU_WORKERS, GPU_WORKERS,
                CAN_GPU, log, "ORB-1 block 1")
"""

E3_TRAIN = r"""
if RUN["E3_train"]:
    log = L.logger("E3_train")
    if CAN_GPU:
        S5.train_e3_stage(CFG, L, need_data, DEVICE, log)
    else:
        log("E3 training waits for a GPU runtime")
    p = L.model_dir(CFG.e3_arms[0], CFG.seeds[0]) / "F_curves.png"
    if p.exists():
        display(Image(str(p)))
"""

E3_FLY = r"""
if RUN["E3_fly"]:
    log = L.logger("E3_fly")
    S5.run_jobs(CFG, S5.jobs_e3(CFG, L, L3), L, DEVICE, CPU_WORKERS, GPU_WORKERS, CAN_GPU, log, "E3")
"""

E4_TRAIN = r"""
RL_ARMS = list(S5.RL_TIER1) + (list(S5.RL_TIER2) if RUN["tier2_rl"] else [])
if RUN["E4_train"]:
    log = L.logger("E4_train")
    if CAN_GPU:
        S5.train_rl_stage(CFG, L, RL_ARMS, DEVICE, log)
    else:
        log("RL training waits for a GPU runtime")
    p = L.policy_dir(RL_ARMS[0], CFG.rl_seeds[0]) / "F_curves.png"
    if p.exists():
        display(Image(str(p)))
"""

E4_FLY = r"""
if RUN["E4_fly"]:
    log = L.logger("E4_fly")
    S5.run_jobs(CFG, S5.jobs_e4(CFG, L, L3, RL_ARMS), L, DEVICE, CPU_WORKERS, GPU_WORKERS, CAN_GPU, log, "E4")
"""

ORB1_MORE = r"""
if RUN["ORB1_more"]:
    log = L.logger("ORB1")
    for b in (1, 2):                    # block 3 becomes eligible only once block 2 is complete
        S5.run_jobs(CFG, S5.orb1_jobs(CFG, L, L3, S5.ORB1_ARMS, [b], "ORB1", log), L, DEVICE, CPU_WORKERS,
                    GPU_WORKERS, CAN_GPU, log, f"ORB-1 block {b + 1}")
    extra = ["mppi_truth"] + S5.best_arms(CFG, L, L3)
    S5.run_jobs(CFG, S5.orb1_jobs(CFG, L, L3, extra, [0], "ORB1", log), L, DEVICE, CPU_WORKERS, GPU_WORKERS, CAN_GPU,
                log, "ORB-1 block 1 extra")
if RUN["explore"]:
    log = L.logger("explore")
    S5.run_jobs(CFG, S5.jobs_explore(CFG, L, L3), L, DEVICE, CPU_WORKERS, GPU_WORKERS, CAN_GPU, log, "exploratory")
"""

OPEN_LOOP = r"""
if RUN["open_loop"]:
    log = L.logger("open_loop")
    if CAN_GPU:
        S5.open_loop_stage(CFG, L, need_data, DEVICE, log)
    else:
        log("open loop waits for a GPU runtime")
"""

ANALYSIS = r"""
if RUN["timing"] and (CAN_GPU or MODE == "smoke"):
    S5.timing_stage(CFG, L, L3, DEVICE, L.logger("timing"))
if RUN["analysis"]:
    t0 = time.time()
    R = A5.analyze(CFG, L, L3, L4, log=L.logger("analysis"))
    pd.set_option("display.max_colwidth", 140)
    display(R["T05"][["id", "question", "measured", "p_holm", "met"]])
    for k in ("F01", "F02", "F03", "F06"):
        if k in R.get("figures", {}):
            display(Image(str(R["figures"][k])))
    print(f"analysis in {time.time() - t0:.0f}s")
"""

REPORT = r"""
S5.status(CFG, L, S5.all_stage_jobs(CFG, L, L3), RUNTIME)
if RUN["report"]:
    print("report:", R5.build_report(L))
S5.update_manifest(L)
print(open(L.study / "README.md", encoding="utf-8").read()[:3000])
"""

CELLS = [
    ("md", INTRO),
    ("md", r"""
## 0. Settings

`MODE`, `TIER` (decision D1; the stages each tier runs are listed in the intro), `RUN_OVERRIDE` to skip a stage
this session, auto-disconnect, and worker counts (automatic unless you set them). Nothing else needs editing:
resuming after a disconnect is Run all. **The notebook disconnects itself** when it finishes or a cell fails;
`AUTO_DISCONNECT = False` turns that off.
"""),
    ("code", SETTINGS),
    ("md", r"""
## 1. Setup

Clones (or pulls) the `space-jepa5` branch into its own folder, installs `pymsis` and `reportlab`, mounts Drive,
opens the v3 folder (required) and the v4 folder (R0's flights and open loop; RQ2's comparison waits without it)
read-only, repairs the v5 folder (stale `.tmp` files, unreadable flight files, `.done` markers whose artifacts are
missing; v5_plan.md §9.1 rule 8), and freezes `00_study/study_config.json` (a later session with different
non-runtime settings raises instead of mixing results). Registers the auto-disconnect hook (CLAUDE.md). A minute.
"""),
    ("code", SETUP + AUTO_DISCONNECT_V5),
    ("md", r"""
## 2. Status

Done / total per stage and arm, flights per hour per arm type (measured on each arm's first 20 flights), projected
hours and Colab units, and each training's state. Writes `00_study/status.json` and appends to
`status_history.csv`. Seconds.
"""),
    ("code", STATUS),
    ("md", r"""
## 3. Benchmark and gate S0

Verifies the frozen suites (`configs/benchmark_v5/`, SHA256SUMS), that the three ORB-1 target ranges lie inside the
nominal footprint, that ≥ 97% of ORB-1 cases are feasible, that the longest nominal flight is shorter than MPPI's
940 s horizon, and that nominal ORB-1 trajectories stay inside the training split's feature ranges (loads the
6,000 flights once). Copies the suites to `01_benchmark/`. Skips when `01_checks/S0.done` exists.
"""),
    ("code", BENCH),
    ("md", r"""
## 4. Legacy RL inventory and re-evaluation

Lists `MyDrive/ReEntryAI/{models, checkpoints, data}` (`00_study/legacy_rl_inventory.csv`: path, run, steps, size,
date) and flies every final model with its own settings on the first 256 cases of `heldout_testset.csv` (the
README's numbers) and on A256 (out of its training distribution). A missing folder is recorded as "not found".
"""),
    ("code", LEGACY),
    ("md", r"""
## 5. Gates V1–V4 and E0

V1: FNPEG with the true parameters lands within 2 km on ≥ 15 of 16 cases, Newton converging in ≤ 6 iterations in
≥ 99% of cycles. V2: load relief at 3 g. V3: the fading-memory filters reach the true drag scale within 3% by 150 s
after 0.05 g. V4: PredGuid within 3.5 km on ≥ 15 of 16. E0: the v5 runner equals `run_lockstep_v4` bit for bit.
**A failed gate stops the notebook** (no test flight may run on it). Passed gates are skipped.
"""),
    ("code", GATES),
    ("md", r"""
## 6. Tuning (DEV48 only, frozen before any test flight)

Vf from the DEV48 classical flights, then random search: ≤ 24 configurations each for FNPEG and PredGuid
(configuration 0 = the defaults), ≤ 12 per RL method at a reduced budget. Objective: success@5km, then mean log
miss, subject to mean peak g ≤ 4.6 g. Each configuration's DEV48 score is saved as `03_tuning/<method>/config_<j>.json`
when its flights are done; `chosen.json` freezes the choice. Test flights of a method wait until it is chosen.
"""),
    ("code", TUNING),
    ("md", r"""
## 7. A256: FNPEG, PredGuid, MPPI + filtered physics, v2 at FNPEG's load

The guidance laws on a CPU pool (one lockstep batch per task), `mppi_filtered` for 5 noise indices on the GPU, then
`v2single@fnpeg_matched` (soft g limit = FNPEG's realized peak g on each condition). One JSON + `.npz` per flight
the moment it lands; finished flights are skipped.
"""),
    ("code", A256),
    ("md", r"""
## 8. ORB-1 block 1

`classical`, `fnpeg`, `predguid`, `mppi_filtered`, `v2single` on cases 0–999 (case i flies training seed i mod 5).
Classical costs ~90 CPU-seconds per flight: a high-CPU CPU runtime is the cheap way to fly it.
"""),
    ("code", ORB1),
    ("md", r"""
## 9. E3 training: four non-JEPA dynamics models × 5 seeds

`mlp_direct`, `gru_direct`, `gru_pi`, `r0_nojepa`, each sized to R0's parameter count, trained with free-running
rollouts; learning rate chosen from {×1, ×0.3, ×3} on the validation split's open-loop loss. Train state to Drive
every 5 minutes (and `.prev`); `ckpt.pt` marks a model done.
"""),
    ("code", E3_TRAIN),
    ("md", r"""
## 10. E3 flights

Each model × 5 seeds under v2's MPPI at 4.5 g on A256. A model not trained yet is skipped with a note.
"""),
    ("code", E3_FLY),
    ("md", r"""
## 11. E4 RL training

`dagger_fnpeg`, `bcppo_fnpeg`, `ppo_hist` (Tier 2: `sac_hist`, `ppo_matched`), 3 seeds each, in the batched fast
environment (ppo_matched: the full simulator). Checkpoints every 5 minutes with optimizer, normalizer and RNG
state; learning curves redrawn at each; `policy.pt` marks a run done; then each policy's fast-environment success
(the fast side of the fast-vs-full gap).
"""),
    ("code", E4_TRAIN),
    ("md", r"""
## 12. E4 flights

Every RL policy in the full simulator on A256 (deterministic, one flight per condition and training seed) and its
snapshots on the DEV48 learning-curve cases.
"""),
    ("code", E4_FLY),
    ("md", r"""
## 13. Tier 2: ORB-1 blocks 2–3, extra ORB-1 arms, exploratory arms

Blocks 2 and 3 for the arms still eligible (fewer than 2 misses beyond 10 km in the blocks before; recomputed from
Drive), `mppi_truth` and the best non-JEPA and RL arms on block 1, and `fnpeg_180`, `predguid_180`,
`fnpeg_oracle` on A256. Off in Tier 1.
"""),
    ("code", ORB1_MORE),
    ("md", r"""
## 14. Open loop

Footprint error of each E3 model and seed on the test split from 400 to 50 s to go
(`03_open_loop/footprint_<arm>_seed<k>.csv`). Cached per file.
"""),
    ("code", OPEN_LOOP),
    ("md", r"""
## 15. Compute per decision, then the analysis

RQ5's timing table (one flight per arm, single CPU thread and GPU), then every table T01–T12, figure F01–F10 (PNG
and PDF, 200 dpi), `captions.json` and the criteria (G1–G7, J, R, Q). Bootstrap results are cached by comparison
and input hash. Runs on any runtime and reports whatever exists.
"""),
    ("code", ANALYSIS),
    ("md", r"""
## 16. Report

`06_report/SpaceJEPA_v5_report.pdf` (summary, criteria, RQ1–RQ5, benchmark datasheet, methods with the chosen
FNPEG and PredGuid settings, deviations from `v5_plan.md` §15, every table), the status file, and
`00_study/README.md`. Rebuilt every run.
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
        c["id"] = f"v5-cell-{i:02d}"                 # stable ids: rebuilding gives the same file
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
