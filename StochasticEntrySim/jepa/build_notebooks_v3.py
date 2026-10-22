"""
Build the four SpaceJEPA v3 study notebooks (jepa/v3_plan.md section 10):

    ReEntryAI_JEPA_v3_0_testset.ipynb   freeze the 64 cases to Drive, coverage figures  (any runtime)
    ReEntryAI_JEPA_v3_1_train.ipynb     v2 and v1 for seeds 0-4, resumable             (GPU)
    ReEntryAI_JEPA_v3_2_fly.ipynb       the flight queue: classical + learned            (GPU + CPUs)
    ReEntryAI_JEPA_v3_3_analyze.ipynb   open loop, statistics, figures, tables, report   (CPU, High-RAM)

All logic lives in jepa/ (testset_v3, study_v3, fly_v3, analyze_v3); the
notebooks configure, call, and show. Edit this builder, not the .ipynb files:

    python jepa/build_notebooks_v3.py        (from StochasticEntrySim/)

Smoke mode (every stage, tiny, numbers meaningless): set JEPA_V3_MODE=smoke.
Locally, JEPA_V3_ROOT sets the study folder and JEPA_WIDE_DATA the data.
"""
import os
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
BRANCH = "jepa-v3"
BADGE = ("[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)]"
         f"(https://colab.research.google.com/github/rk0604/ReEntryAI/blob/{BRANCH}/StochasticEntrySim/{{name}})")

SETUP_MD = r"""
## 0. Setup

On Colab this clones the `jepa-v3` branch, installs `pymsis`, and mounts Google Drive. Every output
goes to **`MyDrive/ReEntryAI_jepa_v3/`** (`ReEntryAI_jepa_v3_smoke/` in smoke mode); the layout is in
`jepa/v3_plan.md` §9 and in the auto-written `00_study/README.md`. The 6,000 flights are read from
`MyDrive/ReEntryAI_jepa_data/wide_6000/`. Every stage is resumable: after a disconnect, run all
cells again and finished work is skipped.
"""

# Release the Colab runtime when the notebook finishes or a cell fails, so it never idles on paid compute
# (repo rule, see CLAUDE.md). Locally (not Colab) it only prints what it would do.
AUTO_DISCONNECT_CODE = r"""
# ---- Auto-disconnect: stop paying for the runtime once the work is done or a cell fails ----
AUTO_DISCONNECT = IN_COLAB      # set to False to keep the runtime (e.g. to explore results by hand)

def release_runtime(reason):
    # Flush Drive, then release the Colab runtime so it stops using compute units.
    if not (IN_COLAB and AUTO_DISCONNECT):
        print(f"[auto-disconnect off] {reason}")
        return
    print(f"[auto-disconnect] {reason}. Flushing Drive and releasing the runtime; outputs stay visible.", flush=True)
    try:
        drive.flush_and_unmount()
    except Exception as e:
        print("  Drive flush:", e)
    from google.colab import runtime
    runtime.unassign()

def _release_on_error(result):
    err = result.error_in_exec or result.error_before_exec
    if err is not None and not isinstance(err, KeyboardInterrupt):     # pressing stop is not a failure
        release_runtime(f"a cell failed ({type(err).__name__}: {str(err)[:200]})")

_ip = get_ipython()
for _cb in list(_ip.events.callbacks["post_run_cell"]):
    if getattr(_cb, "__name__", "") == "_release_on_error":
        _ip.events.unregister("post_run_cell", _cb)
_ip.events.register("post_run_cell", _release_on_error)
print("auto-disconnect on: the runtime is released when the last cell finishes or a cell fails" if AUTO_DISCONNECT
      else "auto-disconnect off (not on Colab, or switched off): the runtime stays up until you disconnect it")
"""

DONE_MD = r"""
## Done: release the runtime

Everything this notebook made is already on Drive. This cell flushes Drive and releases the Colab
runtime, so it stops using compute units (a failed cell does the same). To keep the runtime, set
`AUTO_DISCONNECT = False` in the setup cell before running.
"""
DONE = r"""
release_runtime("all cells finished")
"""

SETUP = r"""
import os, sys
IN_COLAB = "google.colab" in sys.modules
BRANCH = "jepa-v3"
MODE = os.environ.get("JEPA_V3_MODE", "full")     # "full": the study | "smoke": every stage, tiny, numbers meaningless
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
    STUDY_ROOT = f"{DRIVE}/ReEntryAI_jepa_v3" + ("" if MODE == "full" else f"_{MODE}")
    DATA_CANDIDATES = [f"{DRIVE}/ReEntryAI_jepa_data/wide_6000/episodes_wide_6000.pkl",
                       f"{DRIVE}/ReEntryAI_jepa_data/wide_6000"]
    LOCAL_DATA = "/content/data"
else:
    STUDY_ROOT = os.environ.get("JEPA_V3_ROOT", f"jepa_runs/v3_{MODE}")
    DATA_CANDIDATES = [os.environ.get("JEPA_WIDE_DATA", ""), "jepa_runs/episodes_wide_6000.pkl", "jepa_runs/wide_6000"]
    LOCAL_DATA = "jepa_runs"

import time
from pathlib import Path
SIM_DIR = Path.cwd()
assert (SIM_DIR / "rl_env.py").exists(), "start Jupyter in StochasticEntrySim/"
sys.path.insert(0, str(SIM_DIR))
import numpy as np, pandas as pd, torch
import matplotlib.pyplot as plt
from IPython.display import Image, display
%matplotlib inline
from jepa import testset_v3 as TS, study_v3 as S, fly_v3 as F, analyze_v3 as A

CFG = S.StudyConfig.smoke() if MODE == "smoke" else S.StudyConfig()
L = S.Layout(STUDY_ROOT).make()
FROZEN = S.freeze_config(CFG, L)
CPUS = os.cpu_count() or 2
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
print(f"MODE={MODE}  study folder={STUDY_ROOT}  device={DEVICE}  cpus={CPUS}  git={FROZEN['git_commit'][:10]}")
""" + AUTO_DISCONNECT_CODE


def notebook(name, cells, gpu=True):
    """Write a notebook; every one ends by releasing its Colab runtime. gpu=False opens on a High-RAM CPU runtime."""
    cells = list(cells) + [("md", DONE_MD), ("code", DONE)]
    nb = nbf.v4.new_notebook()
    nb["cells"] = [nbf.v4.new_markdown_cell(t.strip("\n")) if k == "md" else nbf.v4.new_code_cell(t.strip("\n"))
                   for k, t in cells]
    nb["metadata"] = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                      "language_info": {"name": "python"},
                      **({"accelerator": "GPU", "colab": {"gpuType": "L4", "provenance": []}} if gpu else
                         {"colab": {"machine_shape": "hm", "provenance": []}})}
    path = os.path.join(SIM, name)
    nbf.write(nb, path)
    print(f"wrote {path} ({len(cells)} cells)")


# ============================================================================
# 0. Test set
# ============================================================================
notebook("ReEntryAI_JEPA_v3_0_testset.ipynb", gpu=False, cells=[
    ("md", f"""
# SpaceJEPA v3 · 0. The test set

{BADGE.format(name="ReEntryAI_JEPA_v3_0_testset.ipynb")}

v3 is a full study of v1 vs v2 vs the classical predictor-corrector (`jepa/v3_plan.md`). This notebook
freezes its **64 entry cases** to Drive and draws their coverage. It needs no GPU (~2 min).

* **16 bridge cases**: the frozen v1/v2 test entries, to connect with earlier results.
* **48 new cases** from a scrambled Sobol sequence (seed 20261018), varying entry altitude (105–135 km),
  speed (7.2–8.4 km/s), flight-path angle (−1.5° to −4°), approach heading (0–360°), and where the
  target sits inside the reachable region. The **target stays at 15°N 15°E**; each entry is placed
  backwards from it.
* Every new target is reachable in all **4 atmospheres** at **≤ 4.5 g** with ≥ 50 km to spare. The
  check uses constant-bank sweeps of the simulator's own equations of motion (no controller is
  involved), restricted to banks whose peak load stays within the 4.5 g study limit.

The set was built once by `python jepa/testset_v3.py` and committed as `configs/v3_testset.csv`.
Here it is copied to Drive and spot-checked against the simulator.
"""),
    ("md", SETUP_MD), ("code", SETUP),
    ("md", "## 1. Freeze the committed test set to Drive"),
    ("code", r"""
for src, dst in (("configs/v3_testset.csv", "v3_testset.csv"), ("configs/v3_testset_footprints.csv", "footprints.csv"),
                 ("configs/v3_testset_spec.json", "testset_spec.json")):
    data = (SIM_DIR / src).read_bytes()
    if not (L.testset / dst).exists() or (L.testset / dst).read_bytes() != data:
        S.atomic_bytes(L.testset / dst, data)
cases = pd.read_csv(L.testset / "v3_testset.csv")
fps = pd.read_csv(L.testset / "footprints.csv")
spec = __import__("json").loads((L.testset / "testset_spec.json").read_text())
print(f"{len(cases)} cases ({(cases.set == 'new').sum()} new from {spec['n_candidates']} Sobol candidates; "
      f"rejected: {spec['rejections']}); MPPI horizon {spec['horizon_steps']} steps "
      f"({spec['horizon_steps'] * 5} s)")
display(cases[cases.set == "new"][["altitude_m", "velocity_mps", "gamma_deg", "chi_deg", "phi_deg", "lambda_deg",
                                   "downrange_km", "crossrange_km", "approach_deg"]].describe()
        .loc[["min", "50%", "max"]].round(2))
"""),
    ("md", r"""
## 2. Feasibility per atmosphere

`feasible_<atmosphere>` = the target is inside that atmosphere's constant-bank footprint within 4.5 g.
Mixing banks reaches at least that region, so this is conservative. The new cases are feasible
everywhere by construction. Some bridge cases are not, e.g. 5 and 14 in the ×0.78 atmosphere,
which v2 also could not reach. Those conditions stay in the primary results; table T02 repeats every
metric without them.
"""),
    ("code", r"""
feas = cases.groupby("set")[[f"feasible_{t}" for t in TS.ATMOSPHERES]].mean().mul(100).round(1)
display(feas)
for t in TS.ATMOSPHERES:
    bad = cases.case[~cases[f"feasible_{t}"]].tolist()
    print(f"{t:>14}: infeasible cases {bad}  max lift-up peak g {cases[f'liftup_peak_g_{t}'].max():.2f}")
"""),
    ("md", "## 3. Spot check: the footprint sweep against real simulator flights"),
    ("code", r"""
import math, rl_env
from jepa import data as D
params = {t: TS.truth_params(d) for t, d in TS.ATMOSPHERES.items()}
r = cases[cases.set == "new"].iloc[0]
Q = TS.entry_frame(*np.radians([r.phi_deg, r.lambda_deg, r.chi_deg]))
fp = TS.footprint(r.altitude_m, r.velocity_mps, math.radians(r.gamma_deg), params["nominal"])
for sdeg in (0.0, 30.0, -30.0):
    env = rl_env.OrionEntryEnv(controller_mode="policy", dispersion=False, fast_atmosphere=True)
    env.reset(seed=0, options={"ic_override": TS.ic_of(r)})
    done = False
    while not done:
        *_, te, tr, info = env.step(np.array([sdeg / 90.0], dtype=np.float32)); done = te or tr
    dr, cr = TS.to_local(Q, env.x[1], env.x[2])
    p = next(q for q in fp if q["sigma_deg"] == sdeg)
    print(f"bank {sdeg:+5.0f} deg: env lands at {dr:7.1f} / {cr:+6.1f} km, sweep says {p['downrange_km']:7.1f} / "
          f"{p['crossrange_km']:+6.1f} km (peak g {info['peak_g']:.2f} vs {p['peak_g']:.2f})")
"""),
    ("md", "## 4. Coverage figures (F01, F02)"),
    ("code", r"""
for p in A.testset_figures(L, cases, fps):
    display(Image(str(p)))
S.update_manifest(L)
print("saved to", L.testset)
"""),
])

# ============================================================================
# 1. Training
# ============================================================================
notebook("ReEntryAI_JEPA_v3_1_train.ipynb", [
    ("md", f"""
# SpaceJEPA v3 · 1. Train v2 and v1, five seeds each

{BADGE.format(name="ReEntryAI_JEPA_v3_1_train.ipynb")}

Same 6,000 flights, one fixed split (episode seed 0, 80/10/10), normalization from the training split.

| | v2 | v1 (fair budget) |
|---|---|---|
| model | 5-member ensemble, 12.3M parameters | 26k-parameter TCN + GRU (+ 6k prober), one model |
| stage 1 | 24,000 steps × 256 windows/member, horizon 24 → 48 → 96 | **24,000 steps × 256**, T = 24 (v1-wide in v2 had 6,000) |
| stage 2 | 6,000 × 512, horizon 16 → 32 → 48 | 6,000 × 512, T = 16 |
| SIGReg | per-step, **λ = 0.08 fixed** (v2's validation sweep) | pooled, λ = 0.02 |

**Seed s** sets the weight init, the data order and window sampling, and the SIGReg projections. v1 also
follows the plan's **convergence rule**: if its validation loss still improved by more than 2% over the
last 20% of training, it trains once more for 50% of the steps, with a warm restart at 30% of the peak
learning rate. The extension is recorded in the checkpoint.

**Runtime** (measured in v2's run): v2 ≈ 50 min per seed on an **L4**, ≈ 1.6 h on a **T4**; v1 ≈ 5–10
min per seed. All five seeds: **L4 ≈ 4.5 h, T4 ≈ 8.5 h**. Checkpoints are written to Drive every 10 min,
so the work can span several sessions: re-run all cells and it resumes. Each finished seed is
`02_models/<v1|v2>/seed_<k>/ckpt.pt`, with its training log and curves.
"""),
    ("md", SETUP_MD), ("code", SETUP),
    ("code", r"""
from jepa import models_v2 as V, models as M1, train_v2 as T2
GPU_MEM_GB = torch.cuda.get_device_properties(0).total_memory / 1e9 if DEVICE == "cuda" else 0.0
if DEVICE == "cuda":
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True
    print(f"GPU {torch.cuda.get_device_name(0)} ({GPU_MEM_GB:.0f} GB), AMP {V.amp_dtype(torch.device('cuda'))}")
elif MODE == "full":
    print("WARNING: no GPU. Runtime -> Change runtime type -> L4 or T4.")
GRAD_CKPT = DEVICE == "cuda" and GPU_MEM_GB < 20      # T4 (15 GB): recompute activations to fit T = 96
torch.set_num_threads(1 if DEVICE == "cuda" else max(1, CPUS // 2))
SEEDS_THIS_SESSION = list(CFG.seeds)                  # e.g. [0, 1] to split the work over sessions
"""),
    ("md", "## 1. Data"),
    ("code", r"""
t0 = time.time()
src = T2.find_wide_dataset(DATA_CANDIDATES, download_dir=LOCAL_DATA)
DATA_PATH = T2.local_copy(src, LOCAL_DATA)
data = S.load_data(CFG, DATA_PATH, DEVICE)
print(f"{data['n_episodes']} flights from {src}: train/val/test = {len(data['train'])}/{len(data['val'])}/"
      f"{len(data['test'])} episodes, {len(data['ds_tr']):,} v2 training windows ({time.time() - t0:.0f}s)")
wm_ = V.WorldModelV2(n_members=CFG.v2_members, **CFG.v2_arch)
print("v2 parameters:", V.param_report(wm_, V.PIProberV2(wm_.M, wm_.d_latent)))
print("v1 parameters:", M1.n_params(M1.JEPAWorldModel()), "+", M1.n_params(M1.PIProber(32)), "prober")
del wm_
"""),
    ("md", r"""
## 2. Train, seed by seed

v2 first, then v1. A finished seed is skipped, and an interrupted one resumes from its last 10-minute
checkpoint.
"""),
    ("code", r"""
log = []
for seed in SEEDS_THIS_SESSION:
    t0 = time.time()
    r2 = S.train_v2_seed(CFG, data, L, seed, DEVICE, grad_ckpt=GRAD_CKPT)
    r1 = S.train_v1_seed(CFG, data, L, seed, DEVICE)
    log.append({"seed": seed, "v2_cached": r2["cached"], "v1_cached": r1["cached"],
                "minutes": round((time.time() - t0) / 60, 1), "v1_extended": r1.get("extended")})
    print(f"seed {seed}: done in {(time.time() - t0) / 60:.1f} min")
    S.update_manifest(L)
display(pd.DataFrame(log))
"""),
    ("md", "## 3. What was trained"),
    ("code", r"""
rows = []
for kind in ("v2", "v1"):
    for seed in CFG.seeds:
        ck = L.ckpt(kind, seed)
        if not ck.exists():
            rows.append({"model": kind, "seed": seed, "status": "not trained yet"}); continue
        ex = torch.load(ck, map_location="cpu", weights_only=False)["extra"]
        rows.append({"model": kind, "seed": seed, "status": "done", "train_min": round(ex.get("train_s", 0) / 60, 1),
                     "extended": None if kind == "v2" else {k: (v or {}).get("extend") for k, v in ex["extended"].items()},
                     "min_eff_rank": min(h["eff_rank"] for h in ex["health"]) if kind == "v2" else None})
display(pd.DataFrame(rows))
for kind in ("v2", "v1"):
    p = L.model_dir(kind, CFG.seeds[0]) / "F_curves.png"
    if p.exists():
        display(Image(str(p)))
"""),
])

# ============================================================================
# 2. Flights
# ============================================================================
notebook("ReEntryAI_JEPA_v3_2_fly.ipynb", [
    ("md", f"""
# SpaceJEPA v3 · 2. Fly the study

{BADGE.format(name="ReEntryAI_JEPA_v3_2_fly.ipynb")}

Every flight is one job keyed by *(controller, load limit, seed, MPPI repeat, case, atmosphere)*. Its
result (`04_flights/results/…json`) and its 1 Hz telemetry with RCS fuel and jet-on time
(`04_flights/telemetry/…npz`) are written to Drive the moment it lands, so a re-run skips finished flights.

| Controller | Seeds | Load limit | Flights |
|---|---|---|---|
| classical predictor-corrector | deterministic | its own guidance | 256 |
| v1 (fair budget), v2 ensemble, v2 single member | 0–4 | **4.5 g** | 3 × 1,280 |
| v1, v2 legacy | 0 | 8 g | 512 |
| v2, two extra MPPI noise streams | 0 | 4.5 g | 512 |

The learned controllers all fly through one lockstep MPPI runner, sharing env setup, planner code, cost
terms, and timing. At each 5 s decision, every flight's 256 candidate plans go through the model in one
GPU batch. A given (seed, case, atmosphere) gets **the same MPPI noise stream for every learned
controller** (common random numbers), so paired comparisons differ only in the model.

**What limits the speed is CPU cores, not the GPU.** The simulator steps on the CPU. The classical
controller is about 6 CPU-hours for 256 flights, and the learned flights about 8 CPU-hours of
simulation. Standard Colab GPU runtimes often have only 2 vCPUs. Each stage prints its measured rate
and an ETA, and you can stop and resume at any time. The classical stage can also run on any machine
with more cores, writing into the same Drive folder.
"""),
    ("md", SETUP_MD), ("code", SETUP),
    ("code", r"""
JOBS = S.build_jobs(CFG, L)
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")   # inherited by the GPU workers
GPU_MEM_GB = torch.cuda.get_device_properties(0).total_memory / 1e9 if DEVICE == "cuda" else 0.0
CPU_WORKERS = CPUS                         # classical: one simulator per core
# learned: processes sharing the GPU, each stepping its own simulators on one core. Budget 6 GB of GPU
# memory per worker (measured peak of one rollout batch: <= 5 GB), keeping 2 GB spare: L4 -> 3, T4 -> 2.
# If a worker still runs out of GPU memory, lower GPU_WORKERS and re-run: finished flights are kept.
GPU_WORKERS = max(1, min(4, CPUS, int((GPU_MEM_GB - 2) // 6) if DEVICE == "cuda" else CPUS))
RUN_CLASSICAL, RUN_LEARNED = True, True
tab = pd.DataFrame([{"controller": j["controller"], "load": j["load"], "seed": j["seed"], "mppi_rep": j["mppi_rep"],
                     "done": S.flight_done(L, j)} for j in JOBS])
display(tab.groupby(["controller", "load", "seed", "mppi_rep"], dropna=False).done.agg(["sum", "size"])
        .rename(columns={"sum": "done", "size": "flights"}))
print(f"{len(JOBS)} flights, {int(tab.done.sum())} done. MPPI horizon {S.horizon_steps(CFG)} steps. "
      f"CPU workers {CPU_WORKERS}, GPU workers {GPU_WORKERS} on {DEVICE}.")
"""),
    ("md", r"""
## 1. Classical predictor-corrector (CPU)

All 256 conditions run again, including the 64 bridge flights that v1/v2 cached, for two reasons:
every classical flight then has telemetry (the cached results had none), and the cache came from an
earlier simulator version. The guidance is deterministic: re-running gives identical flights.
"""),
    ("code", r"""
if RUN_CLASSICAL:
    S.run_classical(JOBS, L, workers=CPU_WORKERS)
"""),
    ("md", r"""
## 2. Learned controllers (GPU + CPUs), seed by seed

Seed 0 runs completely (all controllers and arms), then seed 1, and so on. If the study stops early,
the finished seeds are complete paired comparisons. A model run whose checkpoint is not trained yet is
skipped with a note.
"""),
    ("code", r"""
if RUN_LEARNED:
    S.run_learned(CFG, JOBS, L, device=DEVICE, workers=GPU_WORKERS)
"""),
    ("md", "## 3. Merge and a first look"),
    ("code", r"""
df = S.collect_flights(L)
print(f"{len(df)} flights in {L.flights / 'flights.csv'}")
if len(df):
    df["arm"] = df.controller + " " + df.load.astype(str)
    display(df[df.mppi_rep.fillna(0) == 0].groupby(["arm", "atmosphere"]).agg(
        n=("miss_km", "size"), success5=("success_5km", "mean"), median_km=("miss_km", "median"),
        p90_km=("miss_km", lambda s: s.quantile(0.9)), peak_g=("peak_g", "mean"), fuel_kg=("rcs_fuel_kg", "mean"))
        .round(3).unstack("atmosphere"))
S.update_manifest(L);
"""),
])

# ============================================================================
# 3. Analysis
# ============================================================================
notebook("ReEntryAI_JEPA_v3_3_analyze.ipynb", gpu=False, cells=[
    ("md", f"""
# SpaceJEPA v3 · 3. Analysis and report

{BADGE.format(name="ReEntryAI_JEPA_v3_3_analyze.ipynb")}

**Runtime: CPU with High-RAM** (it opens that way; nothing here uses a GPU, so a GPU only costs credits).
It takes roughly 30–60 min and releases the runtime by itself when it finishes.

Reads the flights and checkpoints from Drive and writes:
* `03_open_loop/`: footprint error, along-track bias, and ensemble spread per model and seed, on the
  held-out test flights;
* `05_analysis/tables/T##_*.csv` and `figures/F##_*.png`, with `captions.json` saying what each shows
  and which file its data came from;
* `06_report/SpaceJEPA_v3_report.pdf`.

The statistics follow `v3_plan.md` §8, fixed before any v3 flight:
* paired by condition (same IC and the same MPPI noise stream);
* a hierarchical bootstrap with 10,000 resamples (seeds, then cases within each atmosphere);
* per-seed exact McNemar and sign tests;
* Holm correction over the hypothesis-test criteria;
* criteria C1–C7 and the Bridge check.

It can run at any point. With only some seeds finished, it reports those.
"""),
    ("md", SETUP_MD), ("code", SETUP),
    ("md", "## 1. Open loop, per model and seed (cached as CSV)"),
    ("code", r"""
from jepa import train_v2 as T2
try:
    DATA_PATH = T2.local_copy(T2.find_wide_dataset(DATA_CANDIDATES, download_dir=LOCAL_DATA), LOCAL_DATA)
    test_eps = A.test_split(CFG, DATA_PATH)
    print(f"{len(test_eps)} held-out test flights")
    for seed in CFG.seeds:
        for kind in ("v1", "v2"):
            t0 = time.time()
            if A.open_loop_seed(CFG, L, test_eps, kind, seed) is not None:
                print(f"  {kind} seed {seed}: {time.time() - t0:.0f}s")
except FileNotFoundError as e:
    print("open loop skipped:", e)
OL = A.open_loop_summary(CFG, L)
for k in ("F30", "F31"):
    if k in OL:
        display(Image(str(OL[k])))
"""),
    ("md", "## 2. Closed loop: tables, statistics, criteria, figures"),
    ("code", r"""
if not (L.testset / "F01_testset_coverage.png").exists():
    A.testset_figures(L, pd.read_csv(SIM_DIR / CFG.testset_csv), pd.read_csv(SIM_DIR / "configs/v3_testset_footprints.csv"))
t0 = time.time()
R = A.analyze(CFG, L)
print(f"analysis of {R.get('n_flights', 0)} flights, seeds {R.get('seeds')}, in {time.time() - t0:.0f}s")
if R:
    pd.set_option("display.max_colwidth", 120)
    display(R["T05"][["id", "question", "measured", "p_holm", "met"]])
"""),
    ("code", r"""
if R:
    cols = ["controller", "atmosphere", "n_seeds", "success5_%", "success5_%_min", "success5_%_max", "miss_median_km",
            "miss_p90_km", "peak_g_mean", "peak_qdot_MWm2_mean", "heat_load_MJm2_mean", "rcs_fuel_kg_mean"]
    display(R["T01"][cols].round(2))
    for k, p in R["figures"].items():
        display(Image(str(p)))
"""),
    ("md", "## 3. Report"),
    ("code", r"""
if R:
    pdf = A.build_report(CFG, L, R, OL)
    print("report:", pdf)
S.update_manifest(L)
print(open(L.study / "README.md").read())
"""),
])
