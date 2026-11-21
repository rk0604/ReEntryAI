"""
Build ReEntryAI_JEPA_v5_needed.ipynb, notebook 1 of the v5 study: "Is JEPA needed?" (jepa/plan_jepa_needed_v5.md).

It replaces ReEntryAI_JEPA_v5.ipynb for the rest of the study (that notebook, revision 2, now stops at Setup with
"superseded"): the same Drive folder (MyDrive/ReEntryAI_jepa_v5), the same study code (study_v5's stage functions),
and only what the plan adds: plan 2's order (E2 and robustness flights, R0 on AS1/DS1, then the E3 arms r0_nojepa,
gru_direct, mlp_direct, gru_pi, each trained, flown on A256, then on AS1/DS1), plan 3's learning-rate fallback
(study_v5.train_e3_stage), plan 4b's stop rules (runguard_v5.stop_run; 00_study/STOPPED.json and ACK_STOP), and no
E4 / RL (notebook 2). Edit this builder, never the .ipynb:

    python jepa/build_notebook_v5_needed.py          (from StochasticEntrySim/)

The auto-disconnect code is the v3 builder's (repo rule, CLAUDE.md), imported, not copied. Smoke mode:
JEPA_V5_MODE=smoke with JEPA_V5_ROOT, JEPA_V3_ROOT, JEPA_V4_ROOT, JEPA_WIDE_DATA set (jepa/smoke_v5.py check 18,
on a copy of the smoke study the old notebook made). JEPA_V5_ACK_STOP sets ACK_STOP (smoke test only).
"""
import os
import sys

import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from build_notebooks_v3 import AUTO_DISCONNECT_CODE, DONE, DONE_MD  # noqa: E402

NAME = "ReEntryAI_JEPA_v5_needed.ipynb"
BRANCH = "space-jepa5-needed"   # until merged into space-jepa5 (then change back)
BADGE = ("[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)]"
         f"(https://colab.research.google.com/github/rk0604/ReEntryAI/blob/{BRANCH}/StochasticEntrySim/{NAME})")
E3_ORDER = ("r0_nojepa", "gru_direct", "mlp_direct", "gru_pi")      # study_v5.NEEDED_E3_ORDER (checked in smoke 18)

_V3_LINE = "AUTO_DISCONNECT = IN_COLAB      # set to False to keep the runtime (e.g. to explore results by hand)"
assert _V3_LINE in AUTO_DISCONNECT_CODE, "build_notebooks_v3.AUTO_DISCONNECT_CODE changed: update this builder"
AUTO_DISCONNECT_V5 = AUTO_DISCONNECT_CODE.replace(
    _V3_LINE, "AUTO_DISCONNECT = AUTO_DISCONNECT   # set in the Settings cell (section 0)")

INTRO = f"""
# SpaceJEPA v5 · Notebook 1: Is JEPA needed?

{BADGE}

**The question.** Is the JEPA world model (v2 / R0 + MPPI) needed, or does something simpler do as well? Three
simpler alternatives, each its own experiment (`jepa/plan_jepa_needed_v5.md` §0): a published, tuned guidance law
(FNPEG, PredGuid; already flown), the same MPPI planner with a physics model and explicit online estimation
(`mppi_filtered`, E2), and a learned model of the same size, data and planner trained without the JEPA objective
(E3). How each outcome reads was fixed before any of this notebook's flights; the analysis fills it in on the
report's "Is JEPA needed?" page (table T15).

**How to run.** Open in Colab (badge above). It opens on an **L4 GPU** runtime: **Runtime → Run all**. Nothing needs
editing. **After any disconnect, crash or release, Run all again**: everything resumes from Drive (finished flights
and models are skipped, training resumes from its last checkpoint). The work takes about 110–120 L4-hours (~240
units), so at least 5 sessions (Colab ends a session after 24 h at the latest); the status cell prints what is left
per block at the start of every session. When the GPU work is done for a session the notebook releases the runtime;
then **switch to a CPU runtime** (Runtime → Change runtime type → CPU) and Run all for the analysis and the report
(a GPU runtime never runs them).

**The notebook disconnects itself** (so it never idles on paid compute units): when the last cell finishes, when a
cell fails, when the GPU work of a pass is done, when the watchdog sees no progress for `STALL_MIN` (45) minutes,
when the session's unit budget (`MAX_SESSION_UNITS`) is used up, and when a **stop rule** fires. Set
`AUTO_DISCONNECT = False` in the Settings cell to keep the runtime.

**Stop rules** (plan §4b). The notebook stops and releases the runtime whenever a run is failing in a way more
compute will not fix, never on a bad result: an E3 arm that diverged at every learning-rate multiplier; an E3 run
whose validation loss is non-finite, diverging, not falling, or at 25% of its steps above 3× its own pilot's final
validation loss; an MPPI arm with 5 or more of its first 20 flights failing (an exception or a non-finite state);
3 failure records in one session; any stage that raises. It writes **`00_study/STOPPED.json`** (the trigger and the
evidence) and `00_study/RUN_SUMMARY.txt`. **The next Run all stops at Setup and prints it** until you set
`ACK_STOP = "<the trigger id>"` in the Settings cell, so a stopped problem is never re-run blindly. Watchdog, budget
and 24 h releases are not stop rules and resume normally.

**Learning rates** (plan §3). The pilots and the choice stay as preregistered. A full run whose training loss goes
non-finite keeps its `failed.json`; the seed is retrained from scratch at the next-best pilot multiplier in its own
folder (`seed_<k>/lr_x<m>/`), and the arm's remaining seeds start there. The report's T10b says which seeds used
which multiplier.

**Order** (plan §2, §5): E2 and robustness flights (`v2single` and `mppi_filtered` on AS1 and DS1, `mppi_filtered`
and `v2single@fnpeg_matched` on A256, both on ORB-1 block 1) → R0 on AS1/DS1 → E3 `r0_nojepa` → `gru_direct` →
`mlp_direct` → `gru_pi`, each: pilots if needed, 5 seeds, A256 flights, AS1/DS1 flights → open loop → compute per
decision → (GPU pass ends) → analysis and report on a CPU runtime. **Not in this notebook:** E4 / RL of any kind
(notebook 2); Tier 2.

**One writer.** This notebook writes `MyDrive/ReEntryAI_jepa_v5/`, the folder `ReEntryAI_JEPA_v5.ipynb` wrote. That
notebook now stops at Setup ("superseded"); stop any session still running it before the first Run all here.
Inputs (read-only, guarded): `MyDrive/ReEntryAI_jepa_v3/`, `MyDrive/ReEntryAI_jepa_v4/` (R0's checkpoints and
flights), `MyDrive/ReEntryAI_jepa_data/wide_6000/`.
"""

SETTINGS = r"""
import os, sys
IN_COLAB = "google.colab" in sys.modules
MODE = os.environ.get("JEPA_V5_MODE", "full")     # "full": the study | "smoke": tiny, numbers meaningless
AUTO_DISCONNECT = IN_COLAB      # release the runtime when the notebook finishes or a cell fails; False keeps it
CPU_WORKERS = None              # CPU-pool workers (None = every core); this notebook flies only MPPI arms
GPU_WORKERS = None              # MPPI flights sharing the GPU (None = 6 GB each, 2 GB spare, at most 4)
GPU_TYPE = "L4" if MODE == "full" else None   # the study's one GPU type (00_study/gpu_platform.json)
# Unattended runs (sections 15 of 2026-11-18 and 2026-11-19)
STALL_MIN = 45                  # release the runtime after this many minutes without progress (None: watchdog off)
MAX_SESSION_UNITS = None        # e.g. 40: start no model or flight group that would not fit in this many units
UNITS_PER_HOUR = None           # this runtime's rate (None: L4 1.98, High-RAM CPU 0.26, CPU 0.08)
STOP_RULES = True               # plan 4b: stop and release when a run is failing (00_study/STOPPED.json)
ACK_STOP = os.environ.get("JEPA_V5_ACK_STOP") or None   # after a stop rule: its trigger id, e.g. "e3_did_not_train:mlp_direct"
NOTEBOOK_REV = 3                # this notebook's revision (an old copy stops at Setup: open it from GitHub)
"""

SETUP = r"""
BRANCH = "space-jepa5-needed"   # until merged into space-jepa5 (then change back)
if IN_COLAB:
    REPO = "/content/ReEntryAI_v5"
    if not os.path.exists(REPO):
        !git clone -q -b {BRANCH} https://github.com/rk0604/ReEntryAI.git {REPO}
    else:
        !git -C {REPO} fetch -q origin {BRANCH} && git -C {REPO} checkout -q -B {BRANCH} FETCH_HEAD
    %pip install -q pymsis reportlab
    os.chdir(f"{REPO}/StochasticEntrySim")
    from google.colab import drive
    drive.mount("/content/drive")
    DRIVE = "/content/drive/MyDrive"
    SFX = "" if MODE == "full" else f"_{MODE}"
    V5_ROOT, V4_ROOT, V3_ROOT = (f"{DRIVE}/ReEntryAI_jepa_v{k}{SFX}" for k in (5, 4, 3))
    DATA_CANDIDATES = [f"{DRIVE}/ReEntryAI_jepa_data/wide_6000/episodes_wide_6000.pkl",
                       f"{DRIVE}/ReEntryAI_jepa_data/wide_6000"]
    LOCAL_DATA = "/content/data"
else:
    V5_ROOT = os.environ.get("JEPA_V5_ROOT", f"jepa_runs/v5_{MODE}")
    V4_ROOT = os.environ.get("JEPA_V4_ROOT", f"jepa_runs/v4_{MODE}")
    V3_ROOT = os.environ.get("JEPA_V3_ROOT", f"jepa_runs/v3_{MODE}")
    DATA_CANDIDATES = [os.environ.get("JEPA_WIDE_DATA", ""), "jepa_runs/episodes_wide_6000.pkl", "jepa_runs/wide_6000"]
    LOCAL_DATA = "jepa_runs"
""" + AUTO_DISCONNECT_V5 + r"""
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
from jepa import runguard_v5 as RG
S5.notebook_check(NOTEBOOK_REV)           # an old copy of this notebook stops here (and the hook releases the runtime)

CFG = S5.StudyConfigV5.smoke() if MODE == "smoke" else S5.StudyConfigV5()
L3 = S4.v3_layout(V3_ROOT)                                             # read-only v3 folder (raises if missing)
L4 = S5.v4_layout(V4_ROOT)                                             # read-only v4 folder (R0's checkpoints)
DATA_DIRS = [c if os.path.isdir(c) else os.path.dirname(c) for c in DATA_CANDIDATES if c and os.path.exists(c)]
L = S5.LayoutV5(V5_ROOT, read_only=[V3_ROOT, V4_ROOT] + DATA_DIRS).make()   # every write is checked
LOG = L.logger("setup")

# ---- Unattended runs: summary, stop rules, session budget, watchdog (section 15) ----
_release_plain = release_runtime

def release_runtime(reason):
    # RUN_SUMMARY.txt says why the runtime was released, then the auto-disconnect release (CLAUDE.md).
    if str(reason).startswith("a cell failed"):
        if RG.STATE.get("stopped"):           # a stop rule released the runtime already (00_study/STOPPED.json)
            print(f"[release] already released by stop rule {RG.STATE['stopped']}")
            return
        RG.SUMMARY.failure("cell", "notebook", reason)
    RG.SUMMARY.released(reason)
    RG.SUMMARY.write()
    _release_plain(reason)

def _watchdog_off_on_stop(result):
    if isinstance(result.error_in_exec, KeyboardInterrupt):       # you pressed stop: you are here
        RG.start_watchdog(None, None, log=print)
        print("watchdog off because you pressed stop; run Setup again to turn it back on")

for _cb in list(_ip.events.callbacks["post_run_cell"]):
    if getattr(_cb, "__name__", "") == "_watchdog_off_on_stop":
        _ip.events.unregister("post_run_cell", _cb)
_ip.events.register("post_run_cell", _watchdog_off_on_stop)
ACKED = RG.check_stopped(L, ACK_STOP, LOG)   # an earlier stop rule: refuse until ACK_STOP names its trigger
_RUNTIME0 = torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu"
RATE = RG.units_per_hour(_RUNTIME0, os.cpu_count() or 2, S5.UNITS_PER_HOUR, UNITS_PER_HOUR)
LIMIT = S5.SessionLimit(deadline=RG.deadline(MAX_SESSION_UNITS, RATE))
S5.set_session_limit(LIMIT)
RG.configure(L, {"notebook": S5.NOTEBOOK_NAME, "runtime": _RUNTIME0, "mode": MODE, "git": S5._GIT.get()[:10],
                 "notebook revision": NOTEBOOK_REV, "stop rules": STOP_RULES,
                 **({"acknowledged stop": ACKED.get("trigger")} if ACKED else {})}, RATE, release_runtime, LIMIT,
             MAX_SESSION_UNITS, False, next_note=lambda: S5.next_pass_note(GPU_TYPE), stop_rules=STOP_RULES)
WATCHDOG = RG.start_watchdog(STALL_MIN, RG.stall_action(L, release_runtime), log=LOG)
LOG(f"session: {RATE:g} units/h" + (f", budget {MAX_SESSION_UNITS:g} units (deadline "
    f"{time.strftime('%H:%M', time.localtime(LIMIT.deadline))})" if LIMIT.deadline else ", no unit budget")
    + f"; STOP_RULES={STOP_RULES}; summary in 00_study/RUN_SUMMARY.txt")
S5.mark_deprecated(L, LOG)               # DAgger's first tuning: kept, marked, ignored (section 15, 2026-11-18)
REPAIR = S5.repair(L, LOG)                                             # stale .tmp, unreadable flights, broken .done
FROZEN = S5.freeze_config_v5(CFG, L, inputs={"v3": V3_ROOT, "v4": V4_ROOT if L4 else None})   # raises on a change
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
RUNTIME = torch.cuda.get_device_name(0) if DEVICE == "cuda" else "cpu"
CPUS = os.cpu_count() or 2
GPU_MEM_GB = torch.cuda.get_device_properties(0).total_memory / 1e9 if DEVICE == "cuda" else 0.0
CPU_WORKERS = CPU_WORKERS or (CPUS if MODE == "full" else min(2, CPUS))
GPU_WORKERS = GPU_WORKERS or (max(1, min(4, CPUS, int((GPU_MEM_GB - 2) // 6))) if DEVICE == "cuda" else 1)
if DEVICE == "cuda":                                                   # one GPU type for every model and MPPI flight
    GPU_PLATFORM = S5.check_gpu_platform(L, RUNTIME, GPU_TYPE, {"torch": torch.__version__, "cuda": torch.version.cuda},
                                         LOG)
CAN_GPU = DEVICE == "cuda" or MODE == "smoke"       # GPU stages; in smoke mode the CPU stands in for the GPU
CAN_CPU = DEVICE != "cuda"                          # the analysis and the report: a CPU runtime only
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
        with RG.allow("copying and loading the 6,000 flights", 60):     # long and silent: the watchdog waits
            DATA_PATH = T2.local_copy(T2.find_wide_dataset(DATA_CANDIDATES, download_dir=LOCAL_DATA), LOCAL_DATA)
            DATA = S5.load_data_v5(CFG, DATA_PATH, DEVICE)
        LOG(f"{DATA['n_episodes']} flights: train/val/test = {len(DATA['train'])}/{len(DATA['val'])}/"
            f"{len(DATA['test'])} episodes ({time.time() - t0:.0f}s)")
    return DATA

def fly(name, jobs, log):
    # One flight block: MPPI flights on the GPU (a CPU runtime only counts them), inside a recorded stage.
    with RG.stage(name, log):
        S5.run_jobs(CFG, jobs, L, DEVICE, CPU_WORKERS, GPU_WORKERS, CAN_GPU, log, name, can_cpu=CAN_CPU)

def e3_arm(arm):
    # One E3 arm (plan 2.2, 3): pilots if needed, its seeds (with the learning-rate fallback), then its A256 and
    # AS1/DS1 flights. A stop rule raises here (StopRun) after writing STOPPED.json and releasing the runtime.
    log = L.logger(f"E3_{arm}")
    if arm not in CFG.e3_arms:
        log(f"{arm} is not an E3 arm of this configuration: skipped")
        return
    with RG.stage(f"{arm} training", log):
        if CAN_GPU:
            S5.train_e3_stage(CFG, L, need_data, DEVICE, log, arms=[arm])
        else:
            log(f"{arm} training waits for a GPU runtime")
    fly(f"{arm} A256", S5.jobs_e3_arm(CFG, L, L3, arm), log)
    fly(f"{arm} AS1/DS1", S5.jobs_robust(CFG, L, L3, L4, arm), log)
    p = S5.e3_model_dir(L, arm, CFG.seeds[0]) / "F_curves.png"
    if p.exists():
        display(Image(str(p)))
    plan = S5.e3_arm_plan(L, CFG, arm)
    print(f"{arm}: seeds {plan['seeds']}; remaining seeds at x{plan['mult']}; diverged "
          f"{ {f'x{m:g}': [x[0] for x in v] for m, v in plan['diverged'].items()} }")

missing = S5.gates_missing(L) + [m for m in ("fnpeg", "predguid") if not S5.is_tuned(L, m)]
LOG(f"MODE={MODE} v5 folder={V5_ROOT} device={DEVICE} ({RUNTIME}) cpus={CPUS} workers: {CPU_WORKERS} CPU, "
    f"{GPU_WORKERS} GPU; git={FROZEN['git_commit'][:10]}; v4 folder {'found' if L4 else 'NOT found (R0 waits)'}; "
    f"repair: {REPAIR}" + (f"; NOT DONE YET (from ReEntryAI_JEPA_v5.ipynb's CPU passes): {missing}" if missing else ""))
if MODE == "full" and DEVICE != "cuda":
    LOG(f"CPU runtime: training and flights only report what remains (Runtime -> Change runtime type -> {GPU_TYPE}); "
        "the analysis and the report run.")
"""

STATUS = r"""
STAGE_JOBS = S5.needed_stage_jobs(CFG, L, L3, L4)
ST = S5.status(CFG, L, STAGE_JOBS, RUNTIME, rl=False, default_rates=S5.NEEDED_RATES)
PROJ = pd.DataFrame(S5.needed_projection(CFG, L, ST, S5.UNITS_PER_HOUR["L4"]))
pd.set_option("display.max_rows", 60)
display(PROJ)
display(pd.DataFrame(ST["training"]))
print(f"{ST['flights_done']}/{ST['flights_total']} flights, {ST['trainings_done']}/{ST['trainings_total']} E3 models; "
      f"~{PROJ.hours_left.sum():.1f} L4-hours (~{PROJ.units_left.sum():.0f} units at "
      f"{S5.UNITS_PER_HOUR['L4']}/h) left; measured rates (flights/h): {ST['rates_by_type']}")
"""

E2 = r"""
log = L.logger("E2_needed")
for name, jobs in S5.needed_e2_blocks(CFG, L, L3, log):
    fly(name, jobs, log)
"""

R0 = r"""
log = L.logger("R0_robust")
if L4 is None:
    log("v4's folder was not found: R0's AS1/DS1 flights wait for it")
fly("R0 AS1/DS1", S5.jobs_robust(CFG, L, L3, L4, S5.R0_ARM), log)
"""

OPEN_LOOP = r"""
log = L.logger("open_loop")
with RG.stage("open loop", log):
    if CAN_GPU:
        S5.open_loop_stage(CFG, L, need_data, DEVICE, log)
    else:
        log("open loop waits for a GPU runtime")
"""

TIMING = r"""
log = L.logger("timing")
with RG.stage("timing", log):
    if CAN_GPU:
        S5.timing_stage(CFG, L, L3, DEVICE, log)
    else:
        log(f"the timing table waits for a {GPU_TYPE} runtime (its MPPI arms are timed on the GPU)")
"""

GPU_END = r"""
if not CAN_CPU:                         # a GPU pass: its GPU work is done; the analysis and report are CPU work
    S5.status(CFG, L, S5.needed_stage_jobs(CFG, L, L3, L4), RUNTIME, rl=False, default_rates=S5.NEEDED_RATES)
    S5.update_manifest(L)
    print("The analysis and the report run on a CPU runtime (Runtime -> Change runtime type -> CPU, then Run all).")
    print(RG.SUMMARY.text())
    release_runtime("GPU pass finished: the analysis and the report wait for a CPU runtime")
else:
    print("CPU runtime: continuing with the analysis and the report")
"""

ANALYSIS = r"""
if not CAN_CPU:
    print("the analysis waits for a CPU runtime")
else:
    log = L.logger("analysis")
    with RG.stage("analysis", log):
        t0 = time.time()
        R = A5.analyze(CFG, L, L3, L4, log=log, exclude_e4=True)
        pd.set_option("display.max_colwidth", 160)
        T15 = R.get("T15")
        if T15 is not None:
            display(T15[T15.part == "outcome"][["question", "conclusion", "evidence", "status"]])
            print("IS JEPA NEEDED?", A5.needed_answer(T15))
        display(R["T05"][["id", "question", "measured", "p_holm", "status"]])
        PENDING = R["T05"][R["T05"].status == "pending"].id.tolist()
        if PENDING:
            print(f"PROVISIONAL: {len(PENDING)} criteria are pending (flights still missing): {', '.join(PENDING)}")
        print(f"analysis in {time.time() - t0:.0f}s")
"""

REPORT = r"""
if CAN_CPU:
    S5.status(CFG, L, S5.needed_stage_jobs(CFG, L, L3, L4), RUNTIME, rl=False, default_rates=S5.NEEDED_RATES)
    with RG.stage("report", L.logger("report")):
        print("report:", R5.build_report(L))
    S5.update_manifest(L)
else:
    print("the report waits for a CPU runtime")
print(RG.SUMMARY.text())
"""


def _e3_md(arm: str, n: int) -> str:
    what = {"r0_nojepa": "R0's architecture trained without the JEPA objective. It alone answers Q-C: J-equiv is "
                         "defined for it only, so it runs first (plan §2.2).",
            "gru_direct": "black-box recurrent dynamics. Seed 0 resumes from its `train_state.pt`.",
            "mlp_direct": "black-box feedforward dynamics. Its preregistered ×3 diverged on all 5 seeds (kept); the "
                          "seeds retrain at the next-best pilot multiplier, ×1, in `seed_<k>/lr_x1/` (plan §3).",
            "gru_pi": "the physics-informed decoder without JEPA."}[arm]
    return f"""
## {n}. E3 `{arm}`

{what} Learning-rate pilots on seed 0 if they have not chosen yet, then seeds 0–4 (the fallback of plan §3 when a
run's loss goes non-finite), then its A256 flights (5 seeds × 256 under v2's MPPI at 4.5 g, paired with R0's v4
flights), then its AS1/DS1 flights (case mod 5, paired with R0's; estimates only, T16). A stop rule stops the
notebook here (`00_study/STOPPED.json`).
"""


CELLS = [
    ("md", INTRO),
    ("md", r"""
## 0. Settings

Nothing needs editing: resuming after a disconnect is Run all. **The notebook disconnects itself** when it finishes,
a cell fails, a GPU pass's work is done, the watchdog sees no progress for `STALL_MIN` minutes, the budget
(`MAX_SESSION_UNITS`) is used up, or a stop rule fires; **`AUTO_DISCONNECT = False` turns that off** (and
`STALL_MIN = None` the watchdog alone). After a stop rule, read `00_study/STOPPED.json`, then set `ACK_STOP` to its
trigger id to continue.
"""),
    ("code", SETTINGS),
    ("md", r"""
## 1. Setup

Clones (or updates) `space-jepa5-needed`, installs `pymsis` and `reportlab`, mounts Drive, opens the v3 and v4 folders
read-only. Registers the auto-disconnect hook first (CLAUDE.md), so a failure anywhere in Setup releases the runtime
too; checks the notebook's revision (an old copy stops here); **stops if an earlier session ended on a stop rule**
and `ACK_STOP` does not name it; starts the watchdog, the session budget and the stop rules; repairs the v5 folder
and checks the frozen `00_study/study_config.json` and the GPU type (L4). A minute.
"""),
    ("code", SETUP),
    ("md", r"""
## 2. Status: what is left, per block

Done / total per block, and projected L4-hours and units per block (flights at each arm type's measured rate, else
the plan's; training at a finished seed's own time, else the arm's learning-rate pilots scaled to the full steps,
else the plan's estimate, named in the `source` column; pilots if they have not chosen), and each E3 model's state
(done, partial, diverged, fallback todo). Writes `00_study/status.json` and `status_history.csv`. **Mid-session:**
as soon as an arm's pilots finish, the E3 cell prints the measured training hours left and writes them to
`00_study/e3_projection.json` (readable on Drive while the notebook runs).
"""),
    ("code", STATUS),
    ("md", r"""
## 3. E2 and the robustness suites (plan §2.1)

In this order: `v2single` AS1 (128), `mppi_filtered` AS1 (128), `v2single` DS1 (256), `mppi_filtered` DS1 (256),
`mppi_filtered` A256 (5 noise indices × 256), `v2single@fnpeg_matched` A256 (5 × 256), `v2single` and
`mppi_filtered` on ORB-1 block 1 (1,000 each). Every flight is an MPPI flight on the GPU; one JSON + `.npz` per flight
the moment it lands; finished flights are skipped.
"""),
    ("code", E2),
    ("md", r"""
## 4. R0 on AS1 and DS1 (plan §2.2b)

v4's R0 checkpoints (seeds 0–4, v4's adapter and settings), one flight per condition with training seed case mod 5
(128 + 256 flights), so each E3 arm's AS1/DS1 flights pair with R0's by condition. Exploratory estimates only.
"""),
    ("code", R0),
]
for _i, _a in enumerate(E3_ORDER):
    CELLS += [("md", _e3_md(_a, 5 + _i)), ("code", f'e3_arm("{_a}")')]
CELLS += [
    ("md", r"""
## 9. Open loop

Footprint error of each finished E3 model (the fallback's where a seed retrained) on the test split from 400 to 50 s
to go (`03_open_loop/footprint_<arm>_seed<k>.csv`). Cached per file.
"""),
    ("code", OPEN_LOOP),
    ("md", r"""
## 10. Compute per decision (GPU runtime)

RQ5's timing table: one flight per arm, single CPU thread, and on the GPU for the MPPI arms. Cached per (arm, device).
"""),
    ("code", TIMING),
    ("md", r"""
## 11. A GPU pass ends here

On the L4 every GPU stage is above: this cell records the status and **releases the runtime**, so the analysis and
the report never run on (and pay for) the GPU. Switch to a CPU runtime for them. On a CPU runtime it does nothing.
"""),
    ("code", GPU_END),
    ("md", r"""
## 12. Analysis (CPU runtime)

Every table and figure of `analyze_v5` and the criteria, plus T15 **"Is JEPA needed?"** (plan §0's outcome table
filled from G1, G2, G4, G5, D1–D4, A1–A4, Q3, J-each, J-equiv and J-mech, each met / not met / pending), T16 (E3 and
R0 on AS1/DS1, estimates) and T10b (the learning rates each seed used). E4 / RL criteria say "not in this notebook".
A criterion whose flights are not all there is **pending**, never met or not met.
"""),
    ("code", ANALYSIS),
    ("md", r"""
## 13. Report

`06_report/SpaceJEPA_v5_report.pdf`, with the "Is JEPA needed?" page after the summary, the status file and
`00_study/README.md`.
"""),
    ("code", REPORT),
    ("md", DONE_MD.replace("in the setup cell", "in the Settings cell")),
    ("code", DONE),
]


def build() -> str:
    nb = nbf.v4.new_notebook()
    nb["cells"] = [nbf.v4.new_markdown_cell(t.strip("\n")) if k == "md" else nbf.v4.new_code_cell(t.strip("\n"))
                   for k, t in CELLS]
    for i, c in enumerate(nb["cells"]):
        c["id"] = f"v5n-cell-{i:02d}"                # stable ids: rebuilding gives the same file
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
