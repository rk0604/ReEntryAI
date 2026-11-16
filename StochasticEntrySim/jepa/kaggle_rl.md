# v5 RL training on Kaggle (ground truth for the build; DRAFT 2026-11-16, not implemented)

Decision (the user, 2026-11-16): **every v5 RL run, tuning included, trains on Kaggle's free GPU**, so all RL
arms share one hardware and one wall-clock budget. Everything else stays where `v5_plan.md` puts it. This file
is the spec; when it is built, `v5_plan.md` §9 and §15/§16 get the summary and this file stays the detail.

## 1. Why, and the rule that follows

* Colab units are the constraint (≈ 50 left; 100 units ≈ $10). RL was the largest GPU item: 15 Tier 1 runs
  capped at 3 h each plus 36 tuning runs at 10 % of that, up to ≈ 55 GPU-hours.
* Every RL run stops at `budget_decisions` **or** `budget_hours` of **wall time**, whichever first
  (`rl_v5.py:349`, `:496`, `:968`; `RLCfg.budget_hours = 3.0`, `budget_decisions = 200_000_000`). A slower
  machine reaches fewer decisions in the same hours, so the hardware is part of the result.
* **Rule: one platform and one GPU type for all RL work.** Every RL tuning run and every RL training run,
  all arms, all seeds, start to finish on Kaggle with the same GPU type. A run never resumes on another GPU
  type. Mixing would make arms (or seeds of one arm) differ by hardware, not by method.
* The cap stays as frozen (3 h); the decisions each run reached are reported per run (T07 "training budget
  used"), with the GPU named. The limitation is written in §16: "RL trained on a Kaggle <GPU>, wall-clock cap
  3 h per seed; the learned controllers' compute was not matched to the RL arms' hardware".

## 2. What runs where

| Work | Platform | Notes |
|---|---|---|
| Setup, benchmark S0, legacy RL inventory, gates V1–V4/E0, guidance tuning (FNPEG, PredGuid), CPU-pool test flights | Colab CPU runtime | as now |
| **RL tuning** (`tune_rl`: 12 configurations per method at 10 % budget + their DEV48 policy flights) | **Kaggle** | needs FNPEG's `chosen.json` first (its expert and clone source) |
| **RL training** (`train_rl_stage`: `dagger_fnpeg`, `bcppo_fnpeg`, `ppo_hist` × seeds 0–4) and each run's `fast_eval.json` | **Kaggle** | needs the RL `chosen.json` of its method |
| E3 training, MPPI flights, open loop, GPU timing | Colab L4 | as now (cost not yet measured) |
| E4 policy flights (RL policies in the full simulator), legacy RL flights, analysis, report | Colab CPU runtime | as now; they read the policies imported from Kaggle |
| Tier 2 RL (`sac_hist`, `ppo_matched`) | Kaggle, same GPU type | only if `TIER = 2` later |

Order: Colab CPU pass (until FNPEG and PredGuid are tuned) → Kaggle RL tuning → Kaggle RL training →
import to Drive → Colab L4 pass (can run in parallel with the Kaggle weeks) → Colab CPU pass (policy flights,
analysis, report).

## 3. Kaggle facts the design relies on (check on first use; Kaggle changes these)

* Free GPU quota ≈ 30 h/week (resets weekly), sessions ≤ 12 h; GPUs "T4 ×2" or "P100". The code uses one
  GPU (`cuda:0`). **Choose T4** unless the first run shows P100 is offered and faster; whichever is chosen
  is fixed for the whole RL study.
* ≈ 4 CPU cores per GPU session. DAgger's FNPEG expert labelling and RL tuning's DEV48 flights run on CPU
  inside the same wall-clock cap, so the CPU matters too (another reason for one platform).
* No Google Drive. Inputs come from Kaggle Datasets (`/kaggle/input/...`, read-only). Only `/kaggle/working`
  (≈ 20 GB) persists, and only as a **saved version's output**.
* Internet (for `git clone`) needs a phone-verified Kaggle account.
* Use **Save Version → Save & Run All (commit)**: it runs headless to the end (≤ 12 h) and keeps
  `/kaggle/working` as that version's output. Interactive sessions burn quota while idle and lose
  `/kaggle/working` when stopped, so they are only for a first look.

## 4. Persistence and resume on Kaggle

* Root: `/kaggle/working/ReEntryAI_jepa_v5/`, **the same layout as the Drive folder** (`00_study/`,
  `02_policies/<arm>/seed_<k>/`, `03_tuning/rl_<arm>/...`, `04_flights/` for the DEV48 tuning flights), so
  every file maps 1:1 to Drive. Written through `LayoutV5` (guarded, atomic, fsync'd) as on Colab.
* Each version starts by copying the **previous version's output** (attached as an input) into
  `/kaggle/working`, then `repair`, then continues: finished runs (`policy.pt`) are skipped, a run in
  progress resumes from `train_state.pt` (saved every 5 min, `RLCfg.save_every_s = 300`). Loss per
  interruption: ≤ 5 min of training; a version that hits the 12 h limit loses ≤ 5 min.
* A version stops itself before Kaggle's 12 h limit: it starts no new run when < 3.3 h remain (a run's cap is
  3 h) and saves; the next version continues.
* `train_state.pt` and `run.json` record the GPU name. **Resume refuses another GPU type** (raises with the
  two names).

## 5. Inputs (Kaggle dataset `reentryai-v5-rl-inputs`, uploaded once by the user from Drive)

* `00_study/study_config.json`: the frozen v5 settings. The Kaggle notebook freezes against it
  (`freeze_config_v5` must find equal settings; it raises otherwise).
* `03_tuning/fnpeg/chosen.json` (FNPEG expert and clone source; `tuned(L, "fnpeg")`).
* Nothing else: the benchmark suites and code come from the `space-jepa5` clone; `verify()` runs as on Colab.
  RL Tier 1 does not read the 6,000-episode data (Tier 2's `ppo_matched` does not either: it draws episodes
  from the fast env).
* The notebook prints the commit it cloned; it must equal the commit in `study_config.json` or the user is
  told (a code change between Colab and Kaggle would be a different study).

## 6. Getting results back to Drive

* Kaggle side writes `manifest.json` in its output: every file under `02_policies/` and `03_tuning/rl_*`
  (and the DEV48 tuning flights) with size and sha256, the GPU name, the commit, the Kaggle version number.
* Colab side, new cell **"Import RL from Kaggle"** (CPU, cheap): reads the output either by the Kaggle API
  (`kaggle kernels output <user>/<kernel>`, token in a Colab secret; the user's own account, nothing sent
  elsewhere) or from a zip the user downloaded and put in Drive (`MyDrive/ReEntryAI_v5_kaggle/`).
  It checks every sha256, refuses to overwrite a different file already on Drive, copies with `LayoutV5`
  writes, and logs what it imported. Re-running it is a no-op.
* Only finished artifacts are imported (a run is imported when its `policy.pt` exists; partial `train_state`
  stays on Kaggle).

## 7. Changes to the code (when built)

1. `build_notebook_v5.py`: setting `RL_ON_KAGGLE = True`. When True the Colab notebook skips `tune_rl` and
   `train_rl_stage` with a message ("RL trains on Kaggle; run the Import cell"); E4 policy flights already wait
   while a `policy.pt` is missing (`run_cpu`, "wait for tuning or training").
2. New builder output `ReEntryAI_JEPA_v5_kaggle_rl.ipynb` (Kaggle notebook): clone, inputs, restore previous
   output, setup (freeze check, `verify`, GPU name), `tune_rl` then `train_rl_stage` with the time guard,
   `manifest.json`, end. No Colab auto-disconnect (not Colab; a committed version ends by itself), but the
   same failure behaviour: a failed cell ends the version with its output saved.
3. `rl_v5`: GPU name into `train_state.pt` and `run.json`; resume check.
4. `study_v5`: `import_kaggle_rl(L, src)` with the manifest checks; `rl_platform` and the GPU name recorded in
   `00_study/`, shown in T07.
5. `smoke_v5`: Kaggle mode on CPU with a temp root standing in for `/kaggle/working` (restore from a previous
   "output", resume, refuse a GPU-name change, manifest, import into a Drive-like folder, refusal on a
   mismatched file, import twice is a no-op); `RL_ON_KAGGLE` skips RL on the Colab side.
6. `v5_plan.md`: §9 (where RL runs), §15 entry, §16 limitation (§1 above).

## 8. Open points

* Is the study already frozen on Drive (did setup complete after commit f4503a7)? If yes and nothing was
  flown, the frozen settings are unchanged by this design (`RL_ON_KAGGLE` is a notebook setting, not a study
  setting), so no new folder is needed. If anything platform-related is added to the frozen settings, it
  must be decided before the first freeze.
* Calendar: ≤ ≈ 55 T4-hours of capped RL work ≈ 2 weeks of quota; T4 vs L4 speed is unmeasured.
* E3 training and MPPI flights still need Colab L4 units; measure their rate in the first L4 pass before
  deciding anything else.
