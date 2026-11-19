# SpaceJEPA v5, notebook 1: "Is JEPA needed?" (plan for a building session)

Written 2026-11-19 after an audit of `MyDrive/ReEntryAI_jepa_v5`. This file is the brief. The study's full
specification stays `jepa/v5_plan.md`; section numbers below (§4.2, §8.1, ...) point there. Where this file and
`v5_plan.md` disagree, this file wins only on the points marked **change**. Each of those also gets a
`v5_plan.md` §15 entry, written before any affected flight or training.

---

## 0. The question and how the answer is decided

**Objective.** Decide whether the JEPA world model (v2 / R0 + MPPI) is needed, or whether something simpler does
as well. "Simpler" means three things, each its own experiment:

| # | Simpler alternative | Experiment | Criterion |
|---|---|---|---|
| Q-A | a published, tuned guidance law (FNPEG, PredGuid), no learning at all | E1 (§4.1), already flown on A256, ORB-1, DS1, AS1 | G1, G2 (A256), Q3 (ORB-1), D1–D2 (DS1), A1–A2 (AS1) |
| Q-B | the same MPPI planner with a physics model plus explicit online estimation (`mppi_filtered`) | E2 (§4.2) | G5 (A256), D3 (DS1), A3 (AS1) |
| Q-C | a learned dynamics model of the same size, data and planner, trained without the JEPA objective | E3 (§4.3) | J-each, J-equiv (A256, log miss vs R0) |

**How the answer reads** (fixed now, before any new flight):

| Outcome | Conclusion written in the report |
|---|---|
| FNPEG ≥ v2 on every suite, and `mppi_filtered` ≈ v2 | Not needed: a classical law or a physics planner with filters does as well |
| v2 > FNPEG on a suite where the model is wrong (AS1, DS1), or v2 > `mppi_filtered` (G5) | A learned model is needed there; then Q-C decides whether the JEPA objective is |
| Q-C: J-each passes for `r0_nojepa` (R0 better) | The JEPA objective is needed |
| Q-C: J-equiv passes for `r0_nojepa` | The JEPA objective is not needed (the architecture and data do the work) |
| Q-C: neither | Inconclusive (as §8.2 says) |
| Q-C: `mlp_direct` / `gru_direct` / `gru_pi` ≈ R0 | Even the physics-informed decoder is not needed |

---

## 1. What already exists (reuse; nothing below is re-flown or retrained)

The notebook uses the **same Drive folder**, `MyDrive/ReEntryAI_jepa_v5`, so every item is picked up by the
existing resume logic (flight ids, `.done` files, `train_state.pt`). No copying.

| Item | State on Drive (2026-11-19) |
|---|---|
| Benchmark (A256, ORB-1, DEV48, DS1, AS1) and gate S0 | done, checksums pass |
| Gates V1–V4, E0 | done, all pass |
| FNPEG tuning (36 configs, chose config 28, 100% on DEV48) | done |
| PredGuid tuning (36 configs, chose config 29, 81% on DEV48) | done |
| A256 flights: `fnpeg`, `predguid`, `fnpeg_180`, `predguid_180`, `fnpeg_oracle`, `fnpeg_norelief` | 6 × 256, done |
| ORB-1 block 1: `classical`, `fnpeg`, `predguid` | 3 × 1,000, done |
| DS1: `classical`, `fnpeg`, `predguid` | 3 × 256, done |
| AS1: `classical`, `fnpeg`, `predguid` | 3 × 128, done |
| v3/v4 A256 flights: `v2single` (5 seeds), `v2`, R0 (5 seeds), `mppi_truth`, `mppi_nominal`, v4 ablations | read from the v3/v4 folders, as now |
| E3 `gru_direct`: learning-rate pilots (chose ×3), seed 0 training | in progress (step ~18,500 of 36,000 at 17:19 UTC) |
| E3 `mlp_direct`: learning-rate pilots (chose ×3) | done; **all 5 seeds diverged** (NaN at steps 3,070–3,714), `failed.json` written |
| `03_tuning/rl_dagger_fnpeg/` | deprecated (supervised-fit bug), ignored; RL is out of this notebook |
| Legacy RL flights (2,560) | done; reported as they are, not part of this notebook's question |

All 8,248 v5 flights were checked against the code: the commits since they flew change only the PredGuid
design cache key (same result for one config), MPPI adapter caching (no MPPI arm has flown), notebook plumbing
and RL. Every flown flight stands.

### What the existing data already says (A256, complete pairs, from `05_analysis/tables`, 2026-11-18)

| Arm | success@5km | success@2km | median miss |
|---|---|---|---|
| `fnpeg` | 96.9% | 96.9% | 0.53 km |
| `v2single` (JEPA + MPPI, 5 seeds) | 96.7% | 54.8% | 1.88 km |
| R0 (v4) | 96.7% | 62.3% | 1.72 km |
| `mppi_truth` (physics MPPI, true parameters) | 96.3% | 57.7% | 1.84 km |
| `predguid` | 75.8% | 48.4% | 2.03 km |
| `mppi_nominal` (physics MPPI, no estimation) | 36.7% | 17.5% | 9.58 km |

* **G1 (v2 beats FNPEG): decided, does not pass.** success@5km −0.2 pp [−1.3, +0.9]. On log miss FNPEG is
  better by 1.16 [1.03, 1.29], about 3× smaller misses.
* **G2 (v2 beats PredGuid): passes**, +20.9 pp [+17.0, +24.8], 5 of 5 seeds.
* success@5km saturates near 97% on A256 (8 of 256 conditions appear unreachable for every arm). That is why
  RQ2 uses paired log miss, not success.
* On the suites where the model is wrong, FNPEG holds on DS1 (98.4%) and collapses on AS1 (6.3%; PredGuid
  14.8%, classical 0%). **AS1 is where a learned model can still show it is needed**, and it costs about
  0.3 L4-hours to find out.

---

## 2. Experiments, arms, seeds and flight counts

All flights use the v5 runner, common random numbers and the 4.5 g limit unless stated. "Seeded" arms on
ORB-1, DS1 and AS1 fly training seed = case mod 5 (§3.2, §3.4). On A256 each seeded arm flies all 5 seeds × 256.

### 2.1 E2 and the robustness suites: the learned model against estimation and model error (GPU, MPPI)

| Order | Arm | Suite | Seeds / noise indices | Flights | L4 h (measured rate) |
|---|---|---|---|---|---|
| 1 | `v2single` | AS1 | case mod 5 | 128 | ~0.2 |
| 2 | `mppi_filtered` | AS1 | MPPI noise index case mod 5 | 128 | ~0.1 |
| 3 | `v2single` | DS1 | case mod 5 | 256 | ~0.4 |
| 4 | `mppi_filtered` | DS1 | case mod 5 | 256 | ~0.2 |
| 5 | `mppi_filtered` | A256 | 5 noise indices × 256 | 1,280 | ~1.2 |
| 6 | `v2single@fnpeg_matched` | A256 | 5 × 256 | 1,280 | ~2.2 |
| 7 | `v2single` | ORB-1 block 1 | case mod 5 | 1,000 | ~1.9 |
| 8 | `mppi_filtered` | ORB-1 block 1 | case mod 5 | 1,000 | ~1.0 |
| | **total** | | | **5,328** | **~7.2** |

`v2single` = v3's v2 checkpoints, member 0, seeds 0–4 (unchanged). `mppi_filtered` as §4.2, with FNPEG's tuned
β. Rates are from `MyDrive/ReEntryAI_v5_gpu_cost/gpu_cost_NVIDIA_L4_20261117_2356.txt`.

**change (order only):** AS1 and DS1 first, then A256, then ORB-1. `v5_plan.md` §4.8 put A256 and ORB-1 first
and DS1/AS1 after E3/E4. The current notebook goes further and runs all of E3 training (~95 L4-h) before any
MPPI flight. Reordering changes no flight. It means the cheapest, most decisive flights exist before the
budget is spent.

Tier 2 (`mppi_truth` on ORB-1, ORB-1 blocks 2–3) stays off.

### 2.2 E3: is the JEPA objective needed (GPU, training then MPPI flights)

Four models, parameter count within ±10% of R0's 2,454,531 (as built in `jepa/models_v5.py`). Each is trained
on R0's data, split, normalization and recipe (36,000 steps × 256, AdamW 6e-4 → 2e-5 cosine, 1,000 warm-up,
weight decay 0.05, clip 1, curriculum 24 → 48 → 96 at 20% / 50%, AMP, free-running rollouts). Each is flown
under v2's MPPI on A256 (5 seeds × 256 = 1,280 flights per arm) and paired with R0's v4 flights.

| Priority | Arm | Isolates | Seeds | Training L4 h | Flight L4 h |
|---|---|---|---|---|---|
| 1 | `r0_nojepa` | the JEPA objective, R0's architecture held fixed | 0–4 | ~37.7 (incl. pilots) | ~2.1 |
| 2 | `gru_direct` | black-box recurrent dynamics | 0–4 (seed 0 running) | ~3.5 left | ~0.9 |
| 3 | `mlp_direct` | black-box feedforward dynamics | 0–4 retrained (see §3) | ~15.2 | ~1.3 |
| 4 | `gru_pi` | physics-informed decoder without JEPA | 0–4 | ~38.7 (incl. pilots) | ~1.2 |

**change (order only):** `r0_nojepa` first. It alone answers Q-C (J-equiv is defined for it only). The current
order is `mlp_direct`, `gru_direct`, `gru_pi`, `r0_nojepa`.

Open loop (§4.7) for every finished model: ~0.01 h, runs after training.

### 2.3 Cheap work on a CPU runtime

* Analysis and report (§10), after each GPU pass, on the CPU runtime.
* Compute per decision (RQ5, §5): the existing stage, on the L4 at the end of a pass (minutes).

### 2.4 Not in this notebook

E4 / RL of any kind (`dagger_fnpeg`, `bcppo_fnpeg`, `ppo_hist`, the RL gates RV1–RV3, RL tuning, the Kaggle
path) moves to notebook 2. The legacy RL flights stay on Drive and in the report as they are. The analysis must
treat missing E4 arms as "not in this notebook", not as pending failures.

---

## 3. E3 learning-rate fix (decision for the user before building)

**What happened.** The pilot trains each multiplier for 3,600 steps with its own short schedule, so its rate
has decayed by the end. ×3 won for `mlp_direct` and `gru_direct`. In the 36,000-step run, ×3 is still near its
peak (1.8e-3) at step ~3,400, and `mlp_direct` diverged there on all 5 seeds. `gru_direct` at ×3 passed that
point and is training normally. `gru_pi` and `r0_nojepa` have not chosen yet and could hit the same problem.

**Proposed rule** (a §15 entry; applies to all four arms, prospectively):
1. A full run that goes NaN is recorded as now (`failed.json`, kept, never deleted).
2. That seed is then retrained once at the next-best pilot multiplier (for `mlp_direct`: ×1, pilot val 0.861
   vs ×3's 0.753) in its own folder (`seed_k/lr_x1/` or similar). The arm's results use that run.
3. If an arm's chosen multiplier diverges on seed 0, the remaining seeds start directly at the fallback.
4. Every multiplier used, and why, goes in the report's tuning table.
5. `gru_direct`'s ×3 runs that finish without NaN stand.

The alternative is to rerun the pilots with the full run's schedule (pilot = first 3,600 steps of the real
36,000-step schedule). That is better science but re-pilots every arm. **The user picks one.**

Also check, with no behaviour change: the early stop (Part A) stops a NaN run within minutes, so a diverged seed
costs ~5 minutes, not hours.

---

## 4. Budget

Units on hand ~126 (129.77 at the screenshot, minus the current session). L4 = 1.98 units/h.

| Block | L4 h | Units | Running total (units) |
|---|---|---|---|
| §2.1 E2 + robustness flights | ~7.2 | ~14 | 14 |
| `gru_direct` (finish) + flights | ~4.4 | ~9 | 23 |
| `r0_nojepa` 5 seeds + flights | ~39.8 | ~79 | 102 |
| `mlp_direct` 5 seeds + flights | ~16.5 | ~33 | 135 |
| `gru_pi` 5 seeds + flights | ~39.9 | ~79 | 214 |

So the current balance covers §2.1, `gru_direct` and `r0_nojepa`, and part of `mlp_direct`. **Decision for the
user:**
- (a) buy ~90 more units for the full E3;
- (b) `mlp_direct` at 3 seeds (~21 units) and drop `gru_pi`;
- (c) `gru_pi` and `mlp_direct` at fewer seeds.

Fewer seeds is a **change** to §4.3 and §8.2's "4 of 5 seeds" rule. It must be a §15 entry with the seed rule
restated (for example "3 of 3") before the training starts.

**Spending guards** (reuse Part A's): watchdog `STALL_MIN = 45`; `MAX_SESSION_UNITS` set by the user per
session; failure records and skip; release after the GPU stages and in the last cell. Each Colab session ends
by its 24 h limit at the latest; Run all resumes.

---

## 5. Notebook design

* **Name:** `StochasticEntrySim/ReEntryAI_JEPA_v5_needed.ipynb`, generated by `jepa/build_notebook_v5_needed.py`.
  Reuse `study_v5.py`'s stage functions; no new study logic beyond §2's order, §3's rule and §2.4's exclusion.
* **Runtime:** L4 for the GPU pass (the study's GPU type is locked to the L4 in `00_study/gpu_platform.json`);
  the analysis pass on a CPU runtime. The notebook's metadata requests an L4 and its text says to switch to CPU
  for the analysis pass. All CPU flights are already done, so there is no CPU flight pass.
* **CLAUDE.md:** `AUTO_DISCONNECT = IN_COLAB`, `release_runtime()` with `drive.flush_and_unmount()` first, the
  `post_run_cell` failure hook (ignores `KeyboardInterrupt`, unregisters an old copy), release in the last
  cell, and text saying it disconnects itself and that `AUTO_DISCONNECT = False` turns it off.
* **Cells, in order:** Settings → Setup (clone/pull `space-jepa5`, mount Drive, repair, platform check) →
  Status (what is done, what is left, projected hours and units per block) → §2.1 flights in the order given →
  E3 `r0_nojepa` → `gru_direct` → `mlp_direct` → `gru_pi` (each: pilot if needed, seeds, then its A256 flights
  right after its seeds) → open loop → compute per decision → release (GPU pass ends) → analysis (CPU) → report
  → release.
* **One writer at a time.** The old `ReEntryAI_JEPA_v5.ipynb` writes the same folder. Bump `NOTEBOOK_REV` so an
  old copy stops at Setup with "superseded by ReEntryAI_JEPA_v5_needed.ipynb". The running old session must be
  stopped before the new notebook first runs. `gru_direct` resumes from `train_state.pt`, so stopping loses
  nothing.
* **Analysis:** the existing `analyze_v5` tables and criteria. Additionally, a one-page "Is JEPA needed?"
  summary that fills §0's outcome table from G1, G2, G4, G5, D1–D4, A1–A4, Q3, J-each, J-equiv and J-mech,
  each with its status (met / not met / pending).

---

## 6. Rules for the building session

* Do not change the simulator, the guidance laws, MPPI, v2/R0, or any arm that has flown. Nothing may change a
  flown flight's id or result.
* Never push to `space-jepa4` or `jepa-v3`. Push only with `git push origin HEAD:space-jepa5`, after the user
  says so.
* Drive is read-only for the building session; nothing outside the study code writes to Drive without asking
  the user.
* Do not run the full study. The user runs it on Colab.
* Smoke test (`python jepa/smoke_v5.py`, CPU) must pass, plus new checks:
  * the stage order of §2;
  * the E3 fallback rule (a forced NaN at the chosen multiplier retrains at the next one and keeps
    `failed.json`);
  * no E4 / RL stage runs;
  * the notebook has the CLAUDE.md release code;
  * an old `NOTEBOOK_REV` stops at Setup.
* Never `git stash` bare; long commands in the background; commit messages end with the repo's
  `Co-Authored-By` line.

---

## 7. Decisions needed from the user before building

1. E3 learning-rate fix: the fallback rule (§3, recommended) or re-piloting with the full schedule.
2. Budget: buy units for full E3, or reduce `mlp_direct` / `gru_pi` seeds (§4).
3. Optional, exploratory (estimates only, no criterion): fly the finished E3 models on AS1 and DS1 too (~0.5
   L4-h). It shows whether a non-JEPA model also handles the model error FNPEG cannot. This is new relative to
   `v5_plan.md` and would need its own §15 entry before any of it flies.
4. Stop the current Colab run of `ReEntryAI_JEPA_v5.ipynb` (now or after `gru_direct` seed 0) before the new
   notebook is pushed.
