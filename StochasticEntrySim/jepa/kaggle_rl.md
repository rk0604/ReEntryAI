# v5 RL training on Kaggle (ground truth; built 2026-11-16)

Built as: `jepa/kaggle_v5.py` (export, inputs, platform, restore, pack, import), `study_v5.SessionLimit` (no run
starts that cannot finish before the session ends), `rl_v5.check_hardware` (a run never resumes on other hardware),
`jepa/build_notebook_v5_kaggle.py` → `ReEntryAI_JEPA_v5_kaggle_rl.ipynb`, the Colab notebook's `RL_ON_KAGGLE`
(section 11 exports the inputs and imports the results), smoke check 16. Where the build differs from the draft
below, the build is right and the draft has been corrected: the state travels as **one tar + manifest per
session** (Kaggle keeps only a limited number of files per output), and the inputs are **one JSON file**
(`v5_rl_inputs.json`, written by the Colab notebook into `00_study/kaggle/to_kaggle/`).

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
| E3 training, MPPI flights, open loop, GPU timing | Colab T4 (L4 until 2026-11-17) | one GPU type, locked by `00_study/gpu_platform.json` (`v5_plan.md` §15) |
| E4 policy flights (RL policies in the full simulator), legacy RL flights, analysis, report | Colab CPU runtime | as now; they read the policies imported from Kaggle |
| Tier 2 RL (`sac_hist`, `ppo_matched`) | Kaggle, same GPU type | only if `TIER = 2` later |

Order: Colab CPU pass (until FNPEG and PredGuid are tuned) → Kaggle RL tuning → Kaggle RL training →
import to Drive → Colab T4 pass (can run in parallel with the Kaggle weeks) → Colab CPU pass (policy flights,
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

* The session's study folder is `/tmp/ReEntryAI_jepa_v5/` (the repo is cloned to `/tmp/ReEntryAI_v5`), **the
  same layout as the Drive folder** (`00_study/`, `02_policies/<arm>/seed_<k>/`, `03_tuning/rl_<arm>/...`,
  `04_flights/` for the DEV48 tuning flights), written through `LayoutV5` as on Colab. The output
  (`/kaggle/working`) holds only `v5_rl_state.tar` (the whole folder) and `v5_rl_manifest.json` (session number,
  the tar's sha256, every file's size and sha256, the finished runs), rewritten after every run, at setup,
  and by a `post_run_cell` hook when a cell fails.
* Each version starts by unpacking the **newest attached earlier output** (highest session number among
  `/kaggle/input`), checked against its manifest, then `repair`, then continues: finished runs (`policy.pt`)
  are skipped, a run in progress resumes from `train_state.pt` (saved every 5 min, `RLCfg.save_every_s =
  300`). Loss per interruption: ≤ 5 min of training.
* A version stops itself before Kaggle's 12 h limit (`SESSION_HOURS = 11.5`): `SessionLimit` starts a run
  only if 1.1 × its wall-clock cap + the work after it (DEV48 flights, fast eval) + 15 min fit; the next
  version continues.
* Each run writes `hardware.json` (GPU, CPU, torch, CUDA, platform) at its first start and into `run.json`;
  **a resume on another GPU type, torch, CUDA or platform raises.** The study-level
  `00_study/rl_platform.json` is written by session 1; **every later session must match it** or the setup
  cell raises, so no RL run of the study trains elsewhere. The full-mode setup also refuses any GPU that is
  not a T4.

## 5. Inputs (Kaggle dataset `reentryai-v5-rl-inputs`, uploaded once by the user from Drive)

One file, `v5_rl_inputs.json`, written by the Colab notebook's section 11 once FNPEG is tuned
(`kaggle_v5.export_inputs`), holding the text and sha256 of:

* `00_study/study_config.json`: the frozen v5 settings. The Kaggle notebook freezes against it
  (`freeze_config_v5` must find equal settings; it raises otherwise).
* `03_tuning/fnpeg/chosen.json` (FNPEG expert and clone source; `tuned(L, "fnpeg")`) and `vf.json`.
* Nothing else: the benchmark suites and code come from the `space-jepa5` clone; `verify()` runs as on Colab.
  RL Tier 1 does not read the 6,000-episode data (Tier 2's `ppo_matched` does not either: it draws episodes
  from the fast env).
* The notebook prints the commit it cloned; it must equal the commit in `study_config.json` or the user is
  told (a code change between Colab and Kaggle would be a different study).

## 6. Getting results back to Drive

* The user downloads the newest version's `v5_rl_state.tar` and `v5_rl_manifest.json` (Output tab) into
  `MyDrive/ReEntryAI_jepa_v5/00_study/kaggle/from_kaggle/` and runs the Colab notebook's section 11
  (`kaggle_v5.import_rl`, any runtime). The tar is cumulative, so the newest session is enough.
* Import checks the tar's and every file's sha256, that `rl_platform.json` matches any platform already
  imported and that every finished run's `hardware.json` matches it; then, if any target file on Drive
  differs, it raises **before writing anything**; otherwise it writes the new files with `LayoutV5`.
  Importing the same session again writes nothing.
* Imported: finished runs only (`02_policies/<arm>/seed_<k>/` and `03_tuning/rl_<arm>/config_<j>/` whose
  `policy.pt` exists, without `train_state*`), `03_tuning/rl_<arm>/` result files, the DEV48 tuning flights
  (`04_flights/{results,telemetry}/DEV48/tune_<rl arm>_c*`), `00_study/rl_platform.json`, and the Kaggle logs
  (into `00_study/kaggle/logs/`, the newest session's copy replacing the older).

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

* The study was frozen on Drive on 2026-11-16 (the user started the Colab CPU pass at f4503a7). Nothing in
  this design is a study setting (`RL_ON_KAGGLE`, the platform files and the session limit live outside
  `study_config.json`), so that folder continues as is. The Colab CPU pass at f4503a7 never trains RL
  (RL needs a GPU runtime there), so no RL work exists outside Kaggle.
* **The Colab GPU (T4) pass must run the code with `RL_ON_KAGGLE`** (it pulls the newest `space-jepa5` on a fresh
  runtime). A GPU session started on f4503a7's code would tune and train RL on Colab.
* Kaggle details to confirm on first use: whether a failed "Save & Run All" version keeps `/kaggle/working`
  (the hook saves before it ends either way), and how an attached notebook output is updated to its newest
  version (the setup cell prints which session it restored).
* Measured 2026-11-17 on a Colab T4 (`ReEntryAI_JEPA_v5_rl_cost.ipynb`, `jepa/rlcost_v5.py`; FNPEG at its
  defaults): PPO 48,500 decisions/s (`ppo_hist`), 43,600 (`bcppo_fnpeg`), so every PPO run reaches its 200 M
  decisions in ≈ 1.2–1.3 h, **before** the 3 h cap: the cap does not bind on a T4, and the decisions each run
  gets do not depend on the GPU. FNPEG expert 212 s per 200-episode chunk (2 CPU cores), 89 s per 20-episode
  chunk; DEV48 policy flights 6.7 s each. All Tier 1 RL: ≈ 30 T4-hours ≈ 54 Colab units, or ≈ 1 week of
  Kaggle quota (per run: `dagger_fnpeg` 0.61 h, `bcppo_fnpeg` 1.87 h, `ppo_hist` 1.15 h).
* E3 training and MPPI flights still need Colab T4 units; measure their rate in the first T4 pass before
  deciding anything else.
