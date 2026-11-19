"""
Build ReEntryAI_RL_v5.ipynb, notebook 2 of SpaceJEPA v5 (jepa/plan_rl_v5.md): RQ3 with RL baselines that never
decide the bank-reversal sign (residual RL on FNPEG, Liu-style reference tracking) and the preregistered PPO.

All logic lives in jepa/rl2_v5.py (arms, environments, gates, flight laws) and jepa/study_rl_v5.py (layout, stages,
routing, stop rules, analysis); the notebook configures, calls and shows. Edit this builder, never the .ipynb:

    python jepa/build_notebook_v5_rl.py          (from StochasticEntrySim/)

The auto-disconnect code is the v3 builder's (repo rule, CLAUDE.md), imported, not copied. Smoke mode (every stage,
tiny, numbers meaningless): JEPA_RL2_MODE=smoke with JEPA_RL2_ROOT and JEPA_V5_ROOT (notebook 1's folder, read
only) set; jepa/smoke_rl_v5.py does this.
"""
import os
import sys

import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from build_notebooks_v3 import AUTO_DISCONNECT_CODE, DONE, DONE_MD  # noqa: E402

NAME = "ReEntryAI_RL_v5.ipynb"
BRANCH = "space-jepa5"
REV = 1                                   # must equal study_rl_v5.NOTEBOOK_REV_RL (smoke check 11)
BADGE = ("[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)]"
         f"(https://colab.research.google.com/github/rk0604/ReEntryAI/blob/{BRANCH}/StochasticEntrySim/{NAME})")

_V3_LINE = "AUTO_DISCONNECT = IN_COLAB      # set to False to keep the runtime (e.g. to explore results by hand)"
assert _V3_LINE in AUTO_DISCONNECT_CODE, "build_notebooks_v3.AUTO_DISCONNECT_CODE changed: update this builder"
AUTO_DISCONNECT_RL = AUTO_DISCONNECT_CODE.replace(
    _V3_LINE, "AUTO_DISCONNECT = AUTO_DISCONNECT   # set in the Settings cell")

INTRO = f"""
# SpaceJEPA v5, notebook 2 · RL baselines that can fly (RQ3)

{BADGE}

**What it answers.** RQ3: how close do properly built RL baselines get to v2 and to FNPEG? Three arms, none of which
lets the policy decide the bank-reversal sign except the preregistered PPO (`jepa/plan_rl_v5.md`):

| Arm | What it is |
|---|---|
| `resid_fnpeg` | residual RL on the tuned FNPEG: the policy adds a bounded correction (±10–30°) to FNPEG's bank magnitude; FNPEG keeps the sign and every reversal; at initialization the arm **is** FNPEG |
| `liu_track` | Liu-style reference tracking (Liu et al., Adv. Space Res. 77 (2026) 1031–1049): PPO flies the bank magnitude to track a reference velocity vector designed per case; reversals by the heading-error deadband (FNPEG's tuned one) |
| `ppo_hist` | PPO from scratch with the world models' memory (preregistered, unchanged) |

`dagger_fnpeg` and `bcppo_fnpeg` failed their competence gates in the laptop pilots and are reported as negative
results with their gate numbers; they do not run here.

**How to run.** Open in Colab (badge above). It opens on a **CPU runtime** (choose the High-RAM CPU runtime:
Runtime → Change runtime type), then **Runtime → Run all**. Nothing needs editing. **After any disconnect, crash or
release, Run all again**: every stage reads Drive, skips finished work and resumes half-done training from its last
checkpoint (at most 5 minutes lost). The notebook **disconnects itself** when the last cell finishes or a cell fails
(so it never idles on paid compute units); set `AUTO_DISCONNECT = False` in the Settings cell to keep the runtime.

**Passes (switch the runtime type yourself).** (1) **CPU**: FNPEG on the 512 fixed episodes, gate RI, the three
pilots (gate RP, ~1 h each) and their CPU rates; no L4 time is spent before the gates. (2) **L4**
(`GPU_TYPE = "L4"`): a 5-minute timing run per arm that passed, then tuning and training of the arms routed to the
L4 (the cheaper runtime per million decisions, frozen in `01_checks/route_<arm>.json`); it releases itself when its
GPU work is done. (3) **CPU**: tuning and training of CPU-routed arms (`resid_fnpeg` is expected here: FNPEG's labels
are CPU work), the full-simulator flights (A256, DS1, AS1, ORB-1 block 1) and the analysis. Status cell and
`00_study/RUN_SUMMARY.txt` say what the next pass must do. Notebook 1 (`ReEntryAI_JEPA_v5`) has priority on the L4
and must have tuned FNPEG first.

**Stop rules (stop on a technical failure, never on a bad result).** A run whose weights or return go non-finite
(or that its training check stops), a failed gate RI, 5 of an arm's first 20 flights raising, 3 failure records in
one session, or any stage raising: the notebook writes `00_study/STOPPED.json` (trigger, evidence) and releases the
runtime. The next Run all **refuses to continue** until you set `ACK_STOP = "<trigger>"` in the Settings cell (the
record then moves to `00_study/stops/`). A failed gate RP or poor accuracy is a result and never stops anything. A
**watchdog** releases the runtime after `STALL_MIN` minutes without progress; `MAX_SESSION_UNITS` caps a session.

**Drive.** Output: **`MyDrive/ReEntryAI_rl_v5/`**. Read only (a write there raises): `MyDrive/ReEntryAI_jepa_v5/`
(FNPEG's tuned configuration, copied with its SHA-256 into `00_study/inputs.json`; its flights for pairing) and the
v3/v4 folders (`v2single`'s A256 flights). Notebook 2 never writes into notebook 1's folder; the analysis joins them.
"""

SETTINGS = rf"""
import os, sys
IN_COLAB = "google.colab" in sys.modules
MODE = os.environ.get("JEPA_RL2_MODE", "full")    # "full": the study | "smoke": every stage, tiny, numbers meaningless
AUTO_DISCONNECT = IN_COLAB      # release the runtime when the notebook finishes or a cell fails; False keeps it
CPU_WORKERS = None              # simulator / FNPEG-label processes (None = every core)
GPU_TYPE = "L4" if MODE == "full" else None    # the one GPU type of notebook 2's GPU work (00_study/gpu_platform.json)
STALL_MIN = 45                  # release the runtime after this many minutes without progress (None: watchdog off)
MAX_SESSION_UNITS = None        # e.g. 20: start no run or flight group that would not fit in this many units
UNITS_PER_HOUR = None           # this runtime's rate (None: L4 1.98, High-RAM CPU 0.26, CPU 0.08)
ACK_STOP = None                 # after a stop: the trigger id in 00_study/STOPPED.json, e.g. "flight_errors"
NOTEBOOK_REV_RL = {REV}             # this notebook's revision (an old copy stops at Setup: open it from GitHub)
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
    DRIVE = "/content/drive/MyDrive"
    SFX = "" if MODE == "full" else f"_{MODE}"
    RL_ROOT = f"{DRIVE}/ReEntryAI_rl_v5{SFX}"
    JEPA_ROOT, V4_ROOT, V3_ROOT = (f"{DRIVE}/ReEntryAI_jepa_v{k}{SFX}" for k in (5, 4, 3))
else:
    RL_ROOT = os.environ.get("JEPA_RL2_ROOT", f"jepa_runs/rl_v5_{MODE}")
    JEPA_ROOT = os.environ.get("JEPA_V5_ROOT", f"jepa_runs/v5_{MODE}")
    V4_ROOT = os.environ.get("JEPA_V4_ROOT", "")
    V3_ROOT = os.environ.get("JEPA_V3_ROOT", "")
""" + AUTO_DISCONNECT_RL + r"""
import json, time
from pathlib import Path
SIM_DIR = Path.cwd()
assert (SIM_DIR / "rl_env.py").exists(), "start Jupyter in StochasticEntrySim/"
sys.path.insert(0, str(SIM_DIR))
import numpy as np, pandas as pd, torch
from IPython.display import display
from jepa import study_v4 as S4, study_v5 as S5, study_rl_v5 as SR, runguard_v5 as RG
SR.notebook_check(NOTEBOOK_REV_RL)       # an old copy of this notebook stops here (and the hook releases the runtime)

CFG = SR.RLStudyConfig.smoke() if MODE == "smoke" else SR.RLStudyConfig()
V3_ROOT = V3_ROOT if V3_ROOT and os.path.isdir(V3_ROOT) else None
V4_ROOT = V4_ROOT if V4_ROOT and os.path.isdir(V4_ROOT) else None
L, J = SR.layout(RL_ROOT, JEPA_ROOT, V3_ROOT, V4_ROOT)             # notebook 1's, v3's and v4's folders read-only
L3 = S4.v3_layout(V3_ROOT) if V3_ROOT else None
L4 = S5.v4_layout(V4_ROOT) if V4_ROOT else None
LOG = L.logger("setup")

_release_plain = release_runtime

def release_runtime(reason):
    # RUN_SUMMARY.txt says why the runtime was released, then the auto-disconnect release (CLAUDE.md).
    if str(reason).startswith("a cell failed"):
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
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
RUNTIME = torch.cuda.get_device_name(0) if DEVICE == "cuda" else "cpu"
CPUS = os.cpu_count() or 2
RATE = RG.units_per_hour(RUNTIME, CPUS, S5.UNITS_PER_HOUR, UNITS_PER_HOUR)
LIMIT = S5.SessionLimit(deadline=RG.deadline(MAX_SESSION_UNITS, RATE))
S5.set_session_limit(LIMIT)
FN = None
RG.configure(L, {"runtime": RUNTIME, "mode": MODE, "git": S5._GIT.get()[:10], "notebook": "RL (notebook 2)",
                 "notebook revision": NOTEBOOK_REV_RL}, RATE, release_runtime, LIMIT, MAX_SESSION_UNITS, True,
             next_note=lambda: SR.next_pass(CFG, L, FN))
WATCHDOG = RG.start_watchdog(STALL_MIN, RG.stall_action(L, release_runtime), log=LOG)
SR.check_stopped(L, ACK_STOP, LOG)        # a stopped study refuses to continue until ACK_STOP names its trigger
REPAIR = S5.repair(L, LOG)
FN = SR.fnpeg_input(L, J, LOG)            # FNPEG's tuned configuration (notebook 1), copied with its SHA-256
FROZEN = SR.freeze(CFG, L, inputs={"jepa_v5": JEPA_ROOT, "v3": V3_ROOT, "v4": V4_ROOT})
CPU_WORKERS = CPU_WORKERS or (CPUS if MODE == "full" else min(2, CPUS))
KIND = SR.runtime_kind(DEVICE, RUNTIME)  # "cpu" | "L4"
if DEVICE == "cuda":
    S5.check_gpu_platform(L, RUNTIME, GPU_TYPE, {"torch": torch.__version__, "cuda": torch.version.cuda}, LOG)
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True
torch.set_num_threads(1 if DEVICE == "cuda" else max(1, CPUS // 2))
ON_CPU = KIND == "cpu"
ON_L4 = KIND == "L4" or MODE == "smoke"   # smoke: the CPU stands in for the L4

def here(arm):
    # Does this runtime tune and train the arm? (its route; smoke: always)
    r = SR.route(L, arm)
    return MODE == "smoke" or (r is not None and r["runtime"] == KIND)

def device_for(arm):
    return "cuda" if (KIND == "L4" and DEVICE == "cuda") else "cpu"

LOG(f"MODE={MODE} folder={RL_ROOT} runtime={RUNTIME} ({KIND}) cpus={CPUS} workers={CPU_WORKERS} "
    f"git={FROZEN['git_commit'][:10]}; notebook 1 folder {JEPA_ROOT} (read only); repair: {REPAIR}")
"""

STATUS = r"""
ST = SR.status(CFG, L, FN)
display(pd.DataFrame(ST["arms"]))
print("gate RI:", ST["RI"], "| NEXT:", ST["next_pass"])
"""

GATES = r"""
log = L.logger("gates")
if ON_CPU or MODE == "smoke":
    with SR.stage(L, CFG, "gates", log):
        SR.fnpeg_fixed(CFG, L, FN, CPU_WORKERS, log)
        if "resid_fnpeg" in CFG.arms:
            SR.gate_ri(CFG, L, FN, CPU_WORKERS, log)
        for arm in CFG.arms:
            SR.gate_rp(CFG, L, arm, FN, "cpu", CPU_WORKERS, "cpu", log)
elif any(SR.gate_passed(CFG, L, a) is None for a in CFG.arms):
    print("the gates run on a CPU runtime first; this L4 session does only the arms whose gates are done")
G = [json.loads(p.read_text()) for p in sorted(L.checks.glob("RP_*.json"))]
if G:
    display(pd.DataFrame([{"arm": g["arm"], "passed": g["passed"], "reason": g.get("reason")} for g in G]))
"""

ROUTING = r"""
log = L.logger("routing")
if ON_L4:
    with SR.stage(L, CFG, "timing runs", log):
        for arm in CFG.arms:
            SR.timing_run(CFG, L, arm, FN, device_for(arm), CPU_WORKERS, "L4", log)
R = []
for arm in CFG.arms:
    r = SR.route(L, arm)
    R.append({"arm": arm, "route": (r or {}).get("runtime", "waits for both rates"),
              **{f"{k}_units_per_million": v for k, v in ((r or {}).get("units_per_million") or {}).items()}})
display(pd.DataFrame(R))
"""

TRAIN = r"""
log = L.logger("train")
for arm in CFG.arms:
    if not SR.gate_passed(CFG, L, arm) or not here(arm):
        continue
    with SR.stage(L, CFG, f"tuning {arm}", log):
        SR.tune(CFG, L, arm, FN, device_for(arm), CPU_WORKERS, log)
    with SR.stage(L, CFG, f"training {arm}", log):
        SR.train(CFG, L, arm, FN, device_for(arm), CPU_WORKERS, log)
    if RG.budget_stopped():
        break
"""

FLIGHTS = r"""
log = L.logger("flights")
if ON_CPU:
    for arm in CFG.arms:
        for suite in SR.SUITES:
            js = SR.flight_jobs(CFG, L, arm, suite, FN)
            if js:
                with SR.stage(L, CFG, f"flights {arm} {suite}", log):
                    SR.run_flights(CFG, L, js, CPU_WORKERS, log, f"{arm} {suite}")
else:
    print("the full-simulator flights wait for a CPU runtime")
"""

GPU_END = r"""
if not ON_CPU:                          # an L4 pass: its GPU work is done; the rest is CPU work
    ST = SR.status(CFG, L, FN)
    print("NEXT:", ST["next_pass"])
    print(RG.SUMMARY.text())
    release_runtime("L4 pass finished: the flights and the analysis wait for a CPU runtime")
else:
    print("CPU runtime: continuing with the analysis")
"""

ANALYSIS = r"""
if ON_CPU:
    log = L.logger("analysis")
    with SR.stage(L, CFG, "analysis", log):
        A = SR.analyze(CFG, L, J, L3, L4, log)
        pd.set_option("display.max_colwidth", 140)
        display(A["criteria"])
        for k in ("gates", "ORB1", "gap", "cost"):
            print(k)
            display(A[k])
    ST = SR.status(CFG, L, FN)
    print("NEXT:", ST["next_pass"])
"""

CELLS = [
    ("md", INTRO),
    ("md", "## 0. Settings\n\nNothing needs editing. `ACK_STOP` only after a stop (see the introduction)."),
    ("code", SETTINGS),
    ("md", r"""
## 1. Setup

Clones the branch, mounts Drive, checks the notebook revision, refuses a stopped study without `ACK_STOP`, copies
FNPEG's tuned configuration from notebook 1 (once, with its SHA-256), freezes the settings, starts the watchdog and
installs the auto-disconnect hook (the runtime is released when a cell fails).
"""),
    ("code", SETUP),
    ("md", "## 2. Status\n\nPer arm: gate, rates, route, tuning, seeds trained, flights per suite; and what the next "
           "pass must do."),
    ("code", STATUS),
    ("md", r"""
## 3. Gates (CPU)

FNPEG alone on the 512 fixed fast-environment episodes (the reference); gate **RI** (`resid_fnpeg` at init equals
FNPEG bit for bit in the fast environment and on DEV48's full-simulator flights: a failure stops the study); gate
**RP** for each arm (an untuned 2·10⁷-decision or 1-hour pilot; `resid_fnpeg` needs ≥ FNPEG − 2 pp and a rising
return, `liu_track` ≥ 50 % of FNPEG, `ppo_hist` ≥ 50 % of FNPEG plus the old RV3). A failed RP is a result: the arm
gets no tuning, training or flights. Each pilot also records its decisions per second on the CPU.
"""),
    ("code", GATES),
    ("md", r"""
## 4. Routing (L4 timing runs)

On the L4: a 5-minute timing run for each arm that passed its gate. Units per million decisions = 10⁶ / rate / 3,600
× units per hour (L4 1.98, High-RAM CPU 0.26); each arm tunes and trains on the cheaper runtime, frozen once both
rates exist.
"""),
    ("code", ROUTING),
    ("md", r"""
## 5. Tuning and training (on each arm's routed runtime)

DEV48 random search (≤ 12 configurations, seed 0, 10 % of the budget, each scored by its 48 full-simulator flights;
an arm whose first 3 configurations all score 0 % stops, a result), then seeds 0–4 at the full budget (decisions
bind; the wall clock is a safety cap of max(3 h, 1.5 × the projected time)), checkpoints every 5 minutes, then each
policy's fast-environment success on the fixed episodes and on A256's cases (the fast side of R-gap).
"""),
    ("code", TRAIN),
    ("md", r"""
## 6. Flights (CPU)

Every trained policy in the full simulator, deterministic: A256 (5 seeds × 256), DS1 and AS1 (seed = case mod 5),
ORB-1 block 1 (1,000 cases, seed = case mod 5), paired by condition with notebook 1's `fnpeg` and v3/v4's
`v2single`.
"""),
    ("code", FLIGHTS),
    ("md", r"""
## 7. An L4 pass ends here

On the L4 this cell prints what the next pass has to do and **releases the runtime**, so the flights and the
analysis never run on (and pay for) the GPU. On a CPU runtime it does nothing.
"""),
    ("code", GPU_END),
    ("md", r"""
## 8. Analysis (CPU)

RQ3 with notebook 1's and v3/v4's flights (read only): R-each (`v2single − arm` on A256, success@5km, 95 % CI
excludes 0, ≥ 4 of 5 seeds, Holm over the three arms; an arm without a comparison enters with p = 1), `arm − fnpeg`
estimates (A256, DS1, AS1, ORB-1 block 1 with Clopper–Pearson bounds), the sim-to-sim gap, gates (with the dropped
arms' numbers), tuning and cost. Tables in `05_analysis/tables/RQ3_*.csv`, summary in `05_analysis/RQ3_summary.md`.
A criterion whose flights are not all there is **pending**.
"""),
    ("code", ANALYSIS),
    ("md", DONE_MD),
    ("code", DONE),
]


def build() -> str:
    nb = nbf.v4.new_notebook()
    nb["cells"] = [nbf.v4.new_markdown_cell(t.strip("\n")) if k == "md" else nbf.v4.new_code_cell(t.strip("\n"))
                   for k, t in CELLS]
    for i, c in enumerate(nb["cells"]):
        c["id"] = f"rl5-cell-{i:02d}"                # stable ids: rebuilding gives the same file
    # opens on a CPU runtime: the first pass (gates) is CPU only; the L4 pass is chosen by the user (CLAUDE.md)
    nb["metadata"] = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                      "language_info": {"name": "python"}, "colab": {"provenance": []}}
    path = os.path.join(SIM, NAME)
    nbf.validate(nb)
    nbf.write(nb, path)
    print(f"wrote {path} ({len(CELLS)} cells)")
    return path


if __name__ == "__main__":
    build()
