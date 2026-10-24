# SpaceJEPA v4: why is v2 better? (plan and implementation brief)

Branch: `space-jepa4` (from `jepa-v3`). Status: **planned, not implemented.** This file is the
preregistration of experiment 4 and the brief for the session that builds it. Read all of it before
writing code. The v1/v2 design is in `jepa/PLAN.md`; v3 is in `jepa/v3_plan.md` (its §13 lists the
implementation decisions v3 made, most of which v4 inherits).

---

## 0. How to use this file (for the implementing session)

1. Work only on the branch `space-jepa4`. Push to `origin space-jepa4`. Never push to `jepa-v3`.
2. Read `CLAUDE.md` (repo rules: every Colab notebook releases its runtime by itself), `jepa/v3_plan.md`,
   and the v3 modules listed in §11 before writing code.
3. Build exactly what §4 to §10 specify. Where the plan leaves a choice open, make it, and record it in
   §14 ("Implementation notes") of this file, as v3 did in its §13.
4. Do **not** change the simulator (`rl_env.py`, `math_3d.py`, physics, RCS, guidance, `configs/`) and do
   not change the behavior of any v3 code path. New behavior goes in new `*_v4.py` modules or behind new
   arguments whose defaults reproduce v3 exactly (§11.3).
5. Do not run the full study. Build, smoke-test on CPU (§12), generate the notebook, commit, push. The
   user runs the notebook on Colab.
6. Do not change the criteria in §8 after the full run has started. If something must change, log it in
   §15 ("Deviations") with the reason and the date, before looking at the affected results.

---

## 1. Where the project stands (context)

ReEntryAI guides an Orion-like capsule through atmospheric entry to a parachute (drogue) deploy point at
15°N 15°E. A learned JEPA world model predicts the trajectory for any bank-command sequence, and MPPI
plans through it (256 samples, replanning every 5 s). Earlier RL attempts failed (PPO ~1% success,
behavior cloning 91 km median miss, clone then PPO 438 km).

* **v1** (`PLAN.md` §9): 26k-parameter TCN + GRU world model with a physics-inspired prober. On 16 test
  cases, nominal air: 100% success within 25 km, 1.69 km median miss (classical 2.24 km), but at 5.6 g
  peak against classical's 3.7 g.
* **v2** (`PLAN.md` §11, `SpaceJEPA_v2_results.pdf`): 12.3M-parameter, 5-member transformer ensemble
  trained on 6,000 wide-envelope flights. Open loop it predicts the landing point 9x better than the
  classical controller's physics model 300 s out (4.1 vs 37.7 km). Closed loop it was not load matched
  (5.8 vs 3.7 g).
* **v3** (`v3_plan.md`, executed `ReEntryAI_JEPA_v3_3_analyze.ipynb`): 5 training seeds per model, 64 entry
  cases × 4 atmospheres, learned controllers held to a 4.5 g soft limit, 5,120 flights. Success within
  5 km (mean over seeds):

  | Atmosphere | classical | v1 | v2 ensemble | v2 single member |
  |---|---|---|---|---|
  | nominal | 51.6% | 85.9% | 97.5% | 97.5% |
  | thin ×0.88 | 21.9% | 29.7% | 98.4% | 98.4% |
  | dense ×1.12 | 90.6% | 26.3% | 98.4% | 98.8% |
  | ×0.78, L/D ×0.80 | 6.3% | 6.6% | 92.2% | 92.2% |

  Criteria: C1 (v2 > v1, +59.5 pp) and C2 (v2 > classical, +54.1 pp, Holm p = 0.0003) met in 5/5 seeds;
  C4, C5, C6 and Bridge met. **C3 (load matched) not met**: v2's mean peak g exceeded classical's by
  +0.50 (nominal), +0.37, +0.56 (dense), +0.18 g. **C7 (ensemble matters) not met**: −0.1 pp.
  Seed 0 at 8 g vs 4.5 g (F14): v2 loses only about 1–2 pp at 4.5 g. Case 10 (a bridge case, γ = −2.16°,
  102.5 km footprint margin in dense air, lift-up peak 2.5 g) is v2's one persistent failure: in the
  v3 per-case heatmap (F15), v2 and v2single miss it in nominal and dense air, while classical lands it
  in every atmosphere and v1 lands it in nominal and dense air.

**The open question v4 answers:** v2 beats classical by 54 points, but v2 differs from the classical
controller in two ways at once: its **planner** (MPPI over a 188-step horizon instead of the predictor
corrector's constant-bank search) and its **model** (learned, history-conditioned, instead of a fixed
nominal model). And v2 differs from v1 in about eight design changes at once. v4 attributes the gains.

---

## 2. Scope

v4 is **only** about understanding why v2 is better. No improvements, no fixes, no stress tests (those
are v5).

| In scope | Out of scope (v5 or never) |
|---|---|
| Planner vs model decomposition (§4.1) | Stress tests: navigation noise, RCS strength, decoupled atmospheres |
| Load-matched comparison, settling C3 (§4.2) | Any simulator change |
| Ablations of the v1 to v2 design changes (§4.3) | Retuning the classical controller |
| Case 10 diagnosis, descriptive (§4.4) | Fixing case 10 or any other failure |
| Success-vs-threshold curves, ensemble re-score with new metrics (§4.5) | New ensemble arms (C7 settled it), new test cases |

**Frozen from v3 (identical, so every v4 flight pairs with every v3 flight):**
* the 6,000 training flights (`MyDrive/ReEntryAI_jepa_data/wide_6000/`), the episode split (seed 0,
  80/10/10), normalization fit on the training split;
* the 64 test cases (`configs/v3_testset.csv`) and the 4 atmospheres (`testset_v3.ATMOSPHERES`);
* seeds 0–4 and what a seed controls (`v3_plan.md` §2.5);
* MPPI settings (256 samples, 6 knots, 30° noise, 4 iterations at entry then 1, horizon 188 steps),
  cost terms and weights, the 4.5 g primary soft limit, `w_unc_km` = 1.0 for the ensemble;
* common random numbers: flight noise seed `fly_v3.crn_seed(seed, case, atmosphere, rep)`;
* success definitions: success@5km = outcome `drogue` and miss < 5 km (`v3_plan.md` §5).

**Why the same test set is correct here:** v4 makes no design decision from test results (no tuning, no
variant selection: every variant is flown; no fixes). Reusing the cases keeps every comparison paired
with v3's 5,120 flights. Fresh cases are v5's job, because v5 changes the controller.

---

## 3. Inputs, outputs, and the Drive layout

**Read-only inputs** (v4 code must never write under these; enforce with an assertion in v4's write
helpers that the target path is inside the v4 study folder):
* `MyDrive/ReEntryAI_jepa_v3/`: `00_study/study_config.json` (the frozen v3 settings),
  `01_testset/`, `02_models/{v1,v2}/seed_<k>/ckpt.pt`, `04_flights/flights.csv` and
  `04_flights/{results,telemetry}/` (5,120 flights).
* `MyDrive/ReEntryAI_jepa_data/wide_6000/` (the 6,000 flights; copied to local disk once per session
  with `train_v2.local_copy`, as v3 does).

**Output: `MyDrive/ReEntryAI_jepa_v4/`** (`ReEntryAI_jepa_v4_smoke/` in smoke mode):

```
MyDrive/ReEntryAI_jepa_v4/
  00_study/
    study_config.json       frozen v4 settings + git commit + the v3 config it was checked against
    manifest.csv, README.md auto-written index (study_v3.update_manifest works on any Layout)
    status.json             what is done / remaining per stage (rewritten by the status cell)
  01_checks/
    E0_*.csv                reproducibility checks (§4.0), physics-model fidelity checks (§4.1.3)
  02_models/
    <variant>/seed_<k>/     ckpt.pt (final), train_state_wm.pt / train_state_pi.pt (while training),
                            train_log.csv, latent_health.csv, F_curves.png, variant.json
  03_open_loop/
    footprint_<variant>_seed<k>.csv, along_track_<variant>_seed<k>.csv
  04_flights/
    results/<arm>/<load>/seed<k>/mppi<m>/<atmosphere>_case<NN>.json    one file per flight
    telemetry/<same path>.npz                                           1 Hz telemetry + decision log
    flights.csv             merged v4 flights (rebuilt from results/)
  05_analysis/
    tables/T##_*.csv, figures/F##_*.png, captions.json (number -> caption, data source)
    cache/                  bootstrap results keyed by comparison + input hash (resumable analysis)
  06_report/
    SpaceJEPA_v4_report.pdf
```

Writes are atomic (temp file + rename, `study_v3.atomic_bytes`), so a disconnect never leaves a
half-written file. Arm and variant names are fixed strings (§4), used in paths, tables and figures.

---

## 4. Experiments

Every learned and MPPI arm flies all 256 conditions (64 cases × 4 atmospheres) for each of the 5 seeds
(1,280 flights per arm), with MPPI rep 0 and the v3 CRN stream, unless the table says otherwise. For
arms with no training seed (physics models), "seed s" means "noise index s": it uses
`crn_seed(s, case, atmosphere, 0)`, the same stream as v2 seed s, so the pairing is exact.

### 4.0 E0: reproducibility checks (run first; they gate everything else)

| Check | Flights | Pass condition |
|---|---|---|
| E0a classical determinism | v3 `classical` re-flown on 32 conditions (8 cases × 4 atmospheres, CPU) | identical to v3's classical flights: same outcome, miss within 1e-6 km, same fuel. If not, stop and report. Also report whether `classical_oracle` in the nominal atmosphere equals `classical` there (it should, unless `apply_domain`'s nominal values differ from the guidance's default parameters; if they differ, say by how much). |
| E0b learned re-fly | v3 `v2single` seed 0 and `v1` seed 0 at 4.5 g, all 256 conditions each, flown by the v4 runner from the v3 checkpoints | outcome agreement ≥ 98%; the 95% CI of the mean paired log-miss difference (v4 − v3) includes 0. Report the spread of the paired differences: it is the GPU-nondeterminism floor for every cross-study comparison. |
| E0c v4 runner equivalence (smoke, CPU) | 3 flights | `run_lockstep_v4` reproduces `run_lockstep_v3` bit-for-bit on CPU for the v1 and v2 adapters (§11.3) |

GPU float32 is not bit-reproducible, and the closed loop amplifies round-off (`v3_plan.md` §13), so E0b
is statistical, not exact. If E0b fails, v4 re-flies the v3 learned arms it compares against (v1, v2,
v2single at 4.5 g, 3,840 flights) with the v4 runner and uses those instead of v3's; record it in §15.

### 4.1 E1: planner vs model (RQ1)

A 2 × 4 grid. Each cell is a planner flying with a model.

| Planner \ model | nominal physics | true physics | physics prior only | learned v2 |
|---|---|---|---|---|
| classical predictor corrector | `classical` (v3) | **`classical_oracle`** (new) | n/a | n/a |
| MPPI (v3 settings) | **`mppi_nominal`** (new) | **`mppi_truth`** (new) | **`mppi_prior`** (new) | `v2` (v3) |

New arms:

| Arm | What it is | Loads | Flights |
|---|---|---|---|
| `classical_oracle` | `fly_v3.fly_classical` with `classical_knows_truth=True` (the guidance predictor gets the true density, L/D and mass). Deterministic, one run per condition. | own guidance | 256 |
| `mppi_nominal` | MPPI planning through the classical controller's own physics model (§4.1.1) with **nominal** parameters | 4.5 g and matched (§4.2) | 2 × 1,280 |
| `mppi_truth` | the same physics model with the flight's **true** parameters (density, lift and mass scales from the job's domain): the ceiling for a perfect model | 4.5 g | 1,280 |
| `mppi_prior` | MPPI through the learned models' own physics backbone with **zero learned corrections** (`PIProberV2`'s crude nominal model: exponential atmosphere, one ballistic coefficient), i.e. v2 with the learning removed | 4.5 g | 1,280 |

What each comparison means (all paired by condition and noise stream):
* `mppi_nominal − classical`: the **planner's** contribution (same nominal model, different planner).
* `v2 − mppi_nominal`: the **learned model's** contribution (same planner, different model).
* `mppi_truth − v2`: what a perfect model would still add (the remaining model gap).
* `classical_oracle − classical`: the classical controller's own model error.
* `v2 − mppi_prior`: what the learned corrections add on top of v2's physics backbone.

#### 4.1.1 The physics model for MPPI (`jepa/physics_v4.py`)

A batched torch port of the translational model the classical guidance predicts with:
`math_3d.nominal_eom_step` and `math_3d.nominal_aero_forces_from_state`. Requirements:
* **Same equations, same parameters, same integration** as the classical guidance's predictor. Read
  the guidance code (the classical branch of `rl_env.OrionEntryEnv`, `env.guidance`) to find its
  integration step and method, and use those. Parameters: `evaluate.nominal_params()` for nominal; for
  truth, apply the flight's domain exactly as `data.apply_domain` does (`ref_area_m2` × density scale,
  `cpas_lift_scale`, `mass_kg` × mass scale).
* **Bank dynamics:** the realized bank starts at the flight's current realized bank (`physical_state`
  index 6) and moves toward each 5 s command at no more than the env's roll-rate limit
  (`fly_v3.SIGMA_RATE_MAX_DPS` = 20°/s), per integration step. Instantaneous bank would make the
  oracle optimistic.
* **State convention:** the runner passes `y0 = [h, φ, λ, V, γ, χ]` (altitude, not radius); convert
  with `r = h + R_EARTH` and back. Run in float64 on the GPU.
* **Outputs in the adapter format** (§4.1.2): `y6 [B, 1, K'+1, 6]` at the 5 s macro steps and
  `gload [B, 1, K']` = aerodynamic acceleration magnitude / g0, taking the **maximum over the
  integration substeps** inside each macro step (record this choice). Stop when every row is below
  drogue altitude (as `V1Adapter.rollout` does).
* `mppi_prior` uses `models.integrate_eom` with zero corrections (the "physics prior only" model of
  `PLAN.md` §3.4) and the same bank-rate handling as v2's prober (whatever v2's prober does: it is
  v2 minus learning, so match v2, not the oracle).

#### 4.1.2 Plugging it into the lockstep runner

`fly_v3.run_lockstep_v3` flies any adapter with `context(y_hists) -> [M, E, D]` and
`rollout(s0, y0, uh, u, n_sub, h_stop) -> (y6, gload)`. A physics adapter needs per-flight
information (realized bank, roll rate, and, for the oracle, the true parameters). Do this:
* `fly_v4.run_lockstep_v4` = a copy of `run_lockstep_v3` with two changes only: (1) it calls
  `adapter.context(y_hists, jobs)` with the deciding flights' jobs (the v1/v2 adapters accept and
  ignore `jobs`); (2) every decision appends `mppi.last` (`J_min`, `miss_best_km`, `spread_km`,
  `ess`, `k_rem`) plus the chosen command and the predicted peak g of the best sample to a per-flight
  **decision log**, saved in the telemetry `.npz` (needed for §4.4).
* `PhysicsAdapter.context` returns `[1, E, 5]` = (realized bank, roll rate, density scale, lift
  scale, mass scale) per flight (scales = 1 for `mppi_nominal`). The runner already repeats the
  context per sample row, so `rollout` reads each row's parameters from `s0[0]`.
* MPPI cost, horizon (188), temperature, iterations, and the g-limit term are unchanged.

#### 4.1.3 Physics-model fidelity checks (stored in `01_checks/`)

1. The torch port matches `math_3d.nominal_eom_step` (float64, 1,000 random states and banks):
   max relative error per step < 1e-9.
2. Open loop: with true parameters and the commands actually flown, the physics model's predicted
   drogue point against the real landing point of 64 recorded v3 flights (any arm, nominal and dense
   air), from 300, 200, 100 s to go. Report the error. This is the oracle's own model error (attitude
   dynamics and integrator differences); `v3_plan.md` §13 found constant-bank sweeps land within
   about 4 km of real flights. The oracle is "the best non-learned model available", not literal truth.

### 4.2 E2: load matching (RQ2, settles C3)

| Arm | Model | g limit | Flights |
|---|---|---|---|
| `v2` at `4.0g` | v3 v2 ensemble checkpoints (no retraining) | 4.0 g | 1,280 |
| `v2` at `3.5g` | same | 3.5 g | 1,280 |
| `v2` at `matched` | same | per condition: classical's measured peak g on that (case, atmosphere), from v3's classical flights | 1,280 |
| `mppi_nominal` at `matched` | §4.1 | same per-condition limit | (counted in §4.1) |

The limit stays a soft cost term (50 km per g² above the limit, as in v3), so the realized load is
checked (criterion L1). Together with v3's 4.5 g arm (5 seeds) and 8 g arm (seed 0), the v2 arms give
a success-vs-peak-g curve.

### 4.3 E3: ablations of the v1 to v2 changes (RQ3)

All ablations use **single-member** models (C7: the ensemble adds nothing), retrained for seeds 0–4.
The reference is retrained through the same pipeline, so every variant differs from it in one thing.

**Reference and anchors**

| Name | What | Trained in v4? |
|---|---|---|
| `R0` | v2's recipe with `n_members = 1`: `WorldModelV2` and `PIProberV2` defaults, λ = 0.08 per-step SIGReg, stage 1 24,000 steps × 256 (curriculum 24 → 48 → 96), stage 2 6,000 × 512 (16 → 32 → 48) | yes, 5 seeds |
| `v1` | v3's v1 with the fair budget | no (v3 checkpoints and flights) |

**A: remove one v2 change from R0 (necessity)**

| Variant | Change from R0 | Implementation |
|---|---|---|
| `A1_inputs_v1` | encoder inputs = v1's 10 features (drop sin σ, cos σ, log a_D, cD_meas, cL_meas) | feature mask over `FEAT_NAMES_V2` |
| `A1b_no_measured_aero` | drop only `cD_meas`, `cL_meas` | feature mask |
| `A2_memory_35s` | encoder receptive field 120 s → 35 s | `enc_band = 3` (rf = 1 + 3·2 = 7) |
| `A3_blind_prober` | prober sees only the latent, not the rollout's h and V | `PIProberV2(state_aware=False)` |
| `A4_train_T24` | training horizon v1's: stage 1 T = 24, stage 2 T = 16, same steps and batches | curricula `((0.0, 24),)` and `((0.0, 16),)` |
| `A5_sigreg_pooled` | v1's regularizer: SIGReg pooled over time, λ = 0.02, v1's projection settings | new `sigreg_mode="pooled"` option in `wm_loss_v2` |
| `A6_latent_32` | latent 64 → 32 | `d_latent = 32` |
| `A7a_size_250k` | about 250k parameters per member (world model + prober), about 10x smaller | scale widths (`enc_dim`, `pred_dim`, prober `hidden`); keep depths, band, latent 64, chunk; heads must divide widths |
| `A7b_size_33k` | about 33k parameters (v1's 26,208 + 6,467 = 32,675) | same, smaller |
| `A8_actions_v1` | v1's command window and input inside v2's tokenizer: `act_chunk = 7`, action feature = u/σmax only (no sin, cos) | new `action_mode="v1"` option |

For A7a and A7b, choose the widths with `models_v2.param_report` before any training, within ±10% of
the target, and record the exact configuration in `variant.json` and §14. A7a, A7b and R0 give a
three-point size curve (about one decade apart), which is what a paper needs to talk about scale.

**B: add one v2 change to v1 (sufficiency)**, each with v1's recipe and v3's fair budget (24,000 × 256,
6,000 × 512, lr 3e-3, the convergence-extension rule of `v3_plan.md` §2.3):

| Variant | Change from v1 | Implementation |
|---|---|---|
| `B1_v1_plus_inputs` | v1 + the five v2 input features (10 → 15) | `n_feat` and the feature function in v1's windows (`study_v3.V1Windows`) |
| `B2_v1_plus_state_prober` | v1's prober also sees the rollout's h and V (v2's normalization `H_NORM`, `V_NORM`) | new `state_aware` option on `models.PIProber` |
| `B3_v1_plus_memory_120s` | v1's TCNs with a 120 s receptive field | kernel 3, dilations (1, 2, 3, 6), channels (32, 32, 32, 32) / (16, 16, 16, 16): rf = 1 + 2·12 = 25; verify `wm.rf == 25` |

Left out on purpose: the zero-initialized residual readout (small change, low prior), and swapping the
whole transformer backbone for TCN + GRU on its own (A7b plus B1–B3 bracket it, and a backbone swap
mixes architecture with size).

**Totals:** 14 trained models (R0, 10 A, 3 B) × 5 seeds = 70 trainings. Each flies 1,280 flights at
4.5 g: 17,920 flights. Each also gets the open-loop evaluation (`analyze_v3.open_loop_seed` style) on the
600 held-out test flights and the 600 validation flights.

### 4.4 E4: case 10 diagnosis (RQ4, descriptive only, no fix)

* **Existing data:** v3 telemetry for case 10 (all arms, all atmospheres) and seed 0 at 8 g vs 4.5 g.
* **New diagnostic flights** with the decision log (§4.1.2): case 10, nominal and dense air, v2 ensemble
  and v2single, seeds 0–4, at 4.5 g and 8 g: 2 × 2 × 5 × 2 = 40 flights, under arm `diag_case10`.
* Also read case 10 off every E1, E2 and E3 arm (free).
* Questions to answer (pre-specified): where in the flight (time, altitude, velocity) v2's trajectory
  separates from classical's; whether the planner's predicted miss at each decision tracked the real
  outcome (model error) or the planner chose a bad plan its model predicted well (planner error);
  whether the predicted peak g of the chosen plan sat at the soft limit (load binding); whether
  `mppi_truth` and `mppi_nominal` land case 10 (does a perfect model fix it?); which ablations change it.
* Output: one figure and one table, with a short written description in the report. No claims of cause
  beyond what these comparisons show.

### 4.5 Extra analyses (no new flights)

* **Success vs threshold:** success rate against the miss threshold from 0.5 to 25 km, every arm, per
  atmosphere, v3 and v4 arms together.
* **Ensemble re-score:** v3's v2 vs v2single with the v4 metrics (paired log miss, success@2km).
* **Gap decomposition:** the classical → v2 gap split into planner and model shares (§8.1), on
  success@5km and on log miss, at 4.5 g and at the matched load.

### 4.6 Flight totals and run order

| Stage | Arms | Flights | Needs training? |
|---|---|---|---|
| E0 | checks | 512 learned + 32 classical | no |
| E1 | classical_oracle, mppi_nominal (4.5 g, matched), mppi_truth, mppi_prior | 256 + 5,120 | no |
| E2 | v2 at 4.0, 3.5, matched | 3,840 | no |
| E3 | R0 + 13 variants | 17,920 | yes (70 runs) |
| E4 | diag_case10 | 40 | no |
| **total** | | **about 27,700** | |

**Run order** (so the most important results exist first if the budget runs out): E0 → E1 → E2 → E4
→ E3 training → E3 flights → open loop → analysis → report. Within each arm, seed-major (seed 0 for
every arm of the stage, then seed 1, ...), as v3 §4.3.

---

## 5. Metrics (definitions; everything else as `v3_plan.md` §5)

| Metric | Definition |
|---|---|
| success@5km (**headline**) | outcome `drogue` and miss < 5 km |
| success@2km | outcome `drogue` and miss < 2 km |
| log miss (**primary for ablations**) | ln(max(miss_km, 0.01)). Failures keep their miss at termination. |
| paired log-miss difference | log miss(A) − log miss(B) on the same (seed, case, atmosphere, noise stream); its mean over conditions is the log of the geometric-mean miss ratio. Below 0 = A misses less. |
| gap closed (B variants) | g = (LM(v1) − LM(B)) / (LM(v1) − LM(R0)), with LM = mean log miss over all conditions and seeds |
| peak g, heat rate, heat load, RCS fuel, bank reversals, planning ms | as v3 |

---

## 6. Statistics (fixed before any v4 flight)

* **Unit:** one (case, atmosphere) condition, paired across arms by IC and noise stream. Seeds (noise
  indices for physics arms) are repeated measures. Deterministic arms (`classical`, `classical_oracle`)
  have one run per condition, broadcast across seeds (as `analyze_v3.cube` does).
* **Uncertainty:** hierarchical bootstrap, 10,000 resamples (seeds, then cases within each atmosphere),
  `analyze_v3.hboot`, extended with a `"mean"` statistic over paired log-miss differences.
  95% percentile intervals; bootstrap p-values as v3.
* **Per-seed checks:** exact McNemar on success, exact sign test on miss, per seed
  (`analyze_v3.mcnemar_p`, `sign_test`). "Holds in a seed" = same sign as the pooled effect.
* **Multiple comparisons:** Holm within each family: RQ1 (P2, P4, P5 tests), RQ2 (L2), RQ3 (the 13
  variant tests). Threshold checks are not Holm-adjusted.
* **Equivalence (P3):** two one-sided tests via the 90% bootstrap interval against ±5 pp.
* Report every estimate with its interval, including the ones that "fail". Secondary metrics are
  reported without correction and labelled exploratory.

---

## 7. Research questions

* **RQ1** How much of v2's advantage over classical comes from the MPPI planner, and how much from the
  learned model? Is v2's model close to a perfect one?
* **RQ2** Does v2 still beat classical when it pulls no more load?
* **RQ3** Which of the v1 to v2 design changes are necessary, and which are sufficient, for v2's
  closed-loop performance? Do they act through prediction quality?
* **RQ4** Why does v2 fail case 10 (descriptive)?

---

## 8. Preregistered criteria

### 8.1 RQ1: planner vs model (4.5 g unless stated)

| # | Question | Measured | Passes if |
|---|---|---|---|
| P1 | How much does the planner add? | success@5km `mppi_nominal − classical` | estimate only (reported with CI, per atmosphere and pooled) |
| P2 | Does the learned model add beyond the planner? | success@5km `v2 − mppi_nominal`, pooled | > 0, 95% CI excludes 0, holds in ≥ 4 of 5 seeds, Holm p < 0.05 |
| P3 | Is v2's model close to a perfect one? | success@5km `mppi_truth − v2`, per atmosphere | equivalent within ±5 pp in every atmosphere (90% CI inside ±5 pp) |
| P4 | Does the classical controller suffer from its model? | success@5km `classical_oracle − classical`, pooled | > 0, 95% CI excludes 0, Holm p < 0.05 |
| P5 | Do the learned corrections matter (v2 vs its own physics backbone)? | success@5km `v2 − mppi_prior`, pooled | > 0, 95% CI excludes 0, ≥ 4 of 5 seeds, Holm p < 0.05 |
| P6 | Decomposition | planner share = P1 / (v2 − classical), model share = P2 / (v2 − classical), with bootstrap CIs (same resamples); also at the matched load | estimate only |

### 8.2 RQ2: load

| # | Question | Passes if |
|---|---|---|
| L1 | Validity: is the matched arm really matched? | `v2@matched` mean peak g ≤ classical's + 0.1 g in every atmosphere (if not, L2 is reported but flagged invalid) |
| L2 | Does v2 beat classical at matched load? | success@5km `v2@matched − classical` > 0 pooled, 95% CI excludes 0, ≥ 4 of 5 seeds, Holm p < 0.05, and no atmosphere where classical is significantly better |
| L3 | What does load buy? | success@5km and log miss against mean peak g for v2 at matched, 3.5, 4.0, 4.5 (v3), 8 g (v3, seed 0): estimate only |

### 8.3 RQ3: ablations (primary metric: paired log miss; Holm over the 13 variants)

| # | Question | Passes if |
|---|---|---|
| A-each | Is the change necessary? (each A variant vs R0) | mean paired log-miss difference (variant − R0) > 0, 95% CI excludes 0, Holm p < 0.05, holds in ≥ 4 of 5 seeds |
| B-each | Is the change sufficient? (each B variant vs v1) | mean paired log-miss difference (B − v1) < 0, 95% CI excludes 0, Holm p < 0.05, holds in ≥ 4 of 5 seeds; and gap closed g ≥ 0.5 with its CI lower bound > 0 to call it "sufficient on its own" |
| R0-check | Retraining reproduces v3 | R0 vs v3 `v2single`: success@5km within ±3 pp in every atmosphere (a check, not a hypothesis) |
| M | Mechanism | Spearman ρ across the 15 models (R0, A, B, v1) between open-loop median error at 300 s (test split) and closed-loop mean log miss: estimate with a bootstrap CI over seeds |

Also reported for every variant: success@5km (headline), success@2km, peak g, fuel, open-loop error
at 400 to 50 s, along-track bias by density bin.

### 8.4 RQ4

Descriptive (§4.4). No pass/fail.

---

## 9. The notebook: `ReEntryAI_JEPA_v4.ipynb` (one notebook, generated)

Built by `jepa/build_notebook_v4.py` (edit the builder, never the `.ipynb`; reuse
`AUTO_DISCONNECT_CODE`, `DONE_MD`, `DONE` from `build_notebooks_v3.py` by importing them; that module
writes all five v3 notebooks at import time (its `ONLY` filter reads `sys.argv`, so it does not protect
an import), so first move its `notebook(...)` calls under `if __name__ == "__main__":` and confirm its
outputs are unchanged). Colab metadata: GPU, L4. Every section starts with a markdown cell saying what it does,
what it reads and writes on Drive, and roughly how long it takes; code cells are commented.

**Resume rule: after any disconnect, Runtime → Run all continues from where it stopped.** Each stage
reads Drive to find what is finished, skips it in seconds, and resumes what is half done. No cell ever
redoes finished work, and no cell needs editing to resume.

| # | Section | What it does | Resume behavior |
|---|---|---|---|
| 0 | Settings | `MODE` ("full" or "smoke"), `RUN = {"E0": True, "E1": True, ...}` per stage, `AUTO_DISCONNECT = IN_COLAB`, worker counts (auto, overridable) | n/a |
| 1 | Setup | clone or pull `space-jepa4`, `pip install pymsis reportlab`, mount Drive, check the v3 folder exists and its `study_config.json` matches v4's frozen v3 settings, freeze `00_study/study_config.json` (raise if settings differ), register the auto-disconnect hook (`CLAUDE.md`) | idempotent |
| 2 | Status | a table of every stage and arm: done / total, measured rate, estimated hours left; writes `00_study/status.json` | always fast (counts result files per arm directory with `os.scandir`, no full tree walk) |
| 3 | E0 checks | §4.0; stops the notebook with a clear message if E0a fails | skips finished checks |
| 4 | E1, E2, E4 flights | classical oracle on a CPU process pool; physics and learned arms through `run_lockstep_v4` on GPU workers (v3's memory budget: 6 GB per worker, 2 GB spare) | each flight's JSON + `.npz` written the moment it lands; finished flights skipped; in-flight ones re-flown |
| 5 | E3 training | for each variant and seed: train, checkpoint, curves | `train_state_*.pt` saved to Drive every 5 min (`save_every_s = 300`) with model, optimizer, scaler, step, RNG state, and loss history; resume from the last state; `ckpt.pt` marks done; train state deleted after |
| 6 | E3 flights | the 14 models × 5 seeds at 4.5 g | as section 4 |
| 7 | Open loop | §4.3, per variant and seed | per-file CSV cache |
| 8 | Analysis | all tables, statistics, criteria, figures (§8, §10) | bootstrap results cached in `05_analysis/cache/` keyed by comparison + a hash of its input flights, so a rerun after new flights recomputes only what changed |
| 9 | Report | `06_report/SpaceJEPA_v4_report.pdf` and a refreshed `00_study/README.md` | rebuilt each run (cheap) |
| 10 | Done | `release_runtime("all cells finished")` | n/a |

Further requirements:
* **No GPU** (CPU runtime): GPU stages print what remains and skip, so the analysis and report can run
  on a cheap CPU runtime with the same notebook.
* Long stages print progress every few minutes and flush Drive state often enough that a disconnect
  loses at most one checkpoint interval of training or the flights in flight.
* Load the 6,000 flights only if a stage that needs them has work left (training, open loop).
* Use `torch.backends.cuda.matmul.allow_tf32` and AMP exactly as v3 does; gradient checkpointing on GPUs
  under 20 GB (T4), as v3.
* Colab can disconnect when idle or after its session limit; the resume rule covers it. Say this in
  the first markdown cell, with the `AUTO_DISCONNECT = False` switch (`CLAUDE.md`).

---

## 10. Tables, figures, and the report

Naming as v3: `T##_<what>.csv`, `F##_<what>_<by-what>.png`, each with a `captions.json` entry (caption
and data source). Minimum set:

| # | Content |
|---|---|
| T01 | every metric per arm × atmosphere (v3 + v4 arms), mean over seeds, seed range, pooled |
| T02 | T01 on feasible conditions only |
| T03 | E0 results (determinism, re-fly agreement, GPU noise floor) and physics-model fidelity |
| T04 | RQ1 paired differences with CIs (P1–P6), per atmosphere and pooled |
| T05 | criteria table (P1–P6, L1–L3, every A and B variant, R0-check, M) with measured value, p (Holm), met |
| T06 | gap decomposition (planner vs model share) at 4.5 g and matched load |
| T07 | load curve (L3) |
| T08 | ablations: per variant, paired log miss vs reference (CI, p, seeds holding), success@5km, success@2km, gap closed, parameter count |
| T09 | ablations per seed |
| T10 | open loop per variant, seed and time-to-go |
| T11 | case 10 diagnosis |
| T12 | success vs threshold (per arm, atmosphere, threshold) |
| T13 | ensemble re-score (v3 data) |
| F10 | success@5km by arm and atmosphere (the E1 grid as bars, with seed whiskers) |
| F11 | gap decomposition (stacked: classical → planner → model → v2 → oracle ceiling) |
| F12 | success vs mean peak g (L3), with classical's point |
| F13 | ablation forest plot: paired log miss vs R0 (A) and vs v1 (B), CIs, Holm-significant marked |
| F14 | size curve (A7b, A7a, R0): log miss and open-loop error against parameters |
| F15 | open-loop error vs closed-loop log miss across the 15 models (mechanism, M) |
| F16 | success vs threshold curves per atmosphere |
| F17 | case 10: trajectories, predicted vs real miss per decision, predicted peak g vs limit |
| F18 | miss-by-case heatmap for the new arms (as v3 F15) |

The report (`jepa/report_v4.py`) reads only these files (no numbers typed in), like
`analyze_v3.build_report` and `results_pdf_v3.py`; reuse `results_pdf_v3.Doc` and
`make_report_v2._fonts` (which falls back to matplotlib's DejaVu fonts on Colab). Sections: summary
callout (headline numbers, criteria met), RQ1 to RQ4 in order, methods and definitions, deviations
(§15), and an appendix with every table.

---

## 11. Code map

### 11.1 v3 modules to reuse (read these first)

| Module | Use in v4 |
|---|---|
| `jepa/study_v3.py` | `StudyConfig` (v3 numbers), `Layout`, atomic writes, `freeze_config`, `update_manifest`, `load_data`, `fit_resumable` (resumable v1 training), `train_v2_seed` / `train_v1_seed` (patterns to generalize), `build_jobs`, `run_classical`, `run_learned` (worker pool pattern), `collect_flights`, `horizon_steps` |
| `jepa/fly_v3.py` | `crn_seed`, `flight_id`, `Recorder`, `flight_summary`, `finish`, `fly_classical` (incl. `classical_knows_truth`), `V1Adapter`, `V2Adapter`, `load_adapter`, `run_lockstep_v3` |
| `jepa/mppi_v2.py` | `EnsembleMPPI`, `MPPIConfigV2`, `_Flight` |
| `jepa/models_v2.py` | `WorldModelV2`, `PIProberV2` (`state_aware`), `raw_features_v2`, `FEAT_NAMES_V2`, `wm_loss_v2`, `per_step_sigreg`, `TrainCfgV2`, `train_stage`, `param_report`, checkpoints |
| `jepa/train_v2.py` | `train_world_model`, `train_prober` (resumable via state paths), `find_wide_dataset`, `local_copy`, `latent_health_v2` |
| `jepa/models.py` | v1 `JEPAWorldModel`, `PIProber`, `SIGReg` (pooled), `integrate_eom`, losses |
| `jepa/evaluate.py`, `jepa/evaluate_v2.py` | `nominal_params`, `physics_nominal_crossing`, `footprint_points`, `member_spread` |
| `jepa/analyze_v3.py` | `hboot`, `mcnemar_p`, `sign_test`, `holm`, `cube`, `metric_table`, labels and colors, `open_loop_seed` |
| `jepa/results_pdf_v3.py`, `jepa/make_report_v2.py` | PDF building blocks |
| `jepa/build_notebooks_v3.py` | notebook builder pattern, `AUTO_DISCONNECT_CODE`, `DONE_MD`, `DONE` |
| `math_3d.py` | `nominal_eom_step`, `nominal_aero_forces_from_state` (read only; port, do not edit) |

### 11.2 New modules

| Module | Contents |
|---|---|
| `jepa/study_v4.py` | `StudyConfigV4` (every number of this plan, plus `smoke()`), `LayoutV4` (+ read-only v3 layout), arm and variant registries, job building per stage, the run-order queue, status counting |
| `jepa/physics_v4.py` | batched torch port of the nominal EOM and aero (§4.1.1), `PhysicsAdapter` (nominal, truth), `PriorAdapter` |
| `jepa/fly_v4.py` | `run_lockstep_v4` (§4.1.2), decision log, matched-g limits from v3 classical flights, the E0 checks |
| `jepa/ablate_v4.py` | variant definitions → model and training configs; `train_variant_seed` (v2-family and v1-family), resumable |
| `jepa/analyze_v4.py` | statistics, criteria, tables, figures (§8, §10), cached bootstrap |
| `jepa/report_v4.py` | the PDF |
| `jepa/smoke_v4.py` | the CPU smoke test (§12) |
| `jepa/build_notebook_v4.py` | writes `ReEntryAI_JEPA_v4.ipynb` |

### 11.3 Changes to existing modules (additive only)

Allowed: new keyword arguments whose defaults reproduce v3 exactly. Needed:
* `models_v2`: an input feature mask (A1, A1b), `sigreg_mode` ("per_step" default, "pooled") in the
  world-model loss (A5), `action_mode` ("v2" default, "v1") in `action_features` / the action tokenizer
  (A8). `n_members = 1` should already work (v3 smoke mode uses 2); verify it trains and saves.
* `models` (v1): `n_feat` already exists; add the v2 extra features to v1's window features (B1), a
  `state_aware` option on `PIProber` (B2); TCN channels and dilations are already arguments (B3; the
  channel tuples must have one entry per dilation).
* `study_v3.V1Windows` or a v4 subclass for the B1 features.
* `build_notebooks_v3.py`: only the import guard (§9).

Proof required (in `smoke_v4.py`): with default arguments, every changed function returns bit-identical
outputs to the `jepa-v3` code on fixed CPU inputs (compare against values computed from a pristine
checkout of `jepa-v3`, or recompute with the defaults before and after the change).

---

## 12. Smoke test (`python jepa/smoke_v4.py`, CPU, numbers meaningless)

The implementing session may not have the 6,000 flights. Generate about 200 wide-envelope flights with
`jepa.data.collect_dataset` / `collect_shards` and `WIDE_DOMAIN` (as `ReEntryAI_JEPA_data.ipynb` does),
then build a smoke v3 study folder with the v3 code in smoke mode (`StudyConfig.smoke()`: tiny v1 and
v2, 2 seeds, 2 cases, a few flights), then run v4 on top of it in smoke mode. Checks:

1. Default-argument equivalence of every changed v3 function (§11.3).
2. `run_lockstep_v4` = `run_lockstep_v3` on CPU for v1 and v2 adapters (E0c).
3. Physics port vs `math_3d.nominal_eom_step` (§4.1.3 check 1); bank-rate limit honored; gload shape and
   units; `mppi_truth` with nominal parameters equals `mppi_nominal` exactly.
4. CRN: physics arm index s and v2 seed s get the same noise for the same condition.
5. Every variant builds, trains a few steps, saves, resumes after an interrupted save (kill mid-stage,
   rerun, step count continues), and has the parameter count it claims (A7a, A7b within ±10% of target
   at full size, computed without training).
6. Job queue: finished flights skipped, a deleted result re-flown, matched-g limits read correctly.
7. Analysis on smoke data: every table and figure is produced, criteria table complete, bootstrap cache
   reused on a second run (no recomputation), then recomputed after one result file changes.
8. Report PDF builds (no system DejaVu fonts needed).
9. The notebook, executed locally in smoke mode (`jupyter nbconvert --execute`, `JEPA_V4_MODE=smoke`,
   `JEPA_V4_ROOT`, `JEPA_V3_ROOT`, `JEPA_WIDE_DATA` set), runs end to end, then a second execution
   finds everything done and finishes in under a minute.
10. The v4 code refuses to write under the v3 folder.

---

## 13. Definition of done (for the implementing session)

* All modules in §11.2, the changes in §11.3, the notebook and its builder, committed on `space-jepa4`
  and pushed.
* `smoke_v4.py` passes; its output is pasted into §14.
* §14 lists every choice the session made (integration step of the physics model, A7 widths, gload
  convention, worker counts, anything else) and the measured physics-model fidelity on smoke data.
* `v3_plan.md` and the v3 notebooks unchanged; v3 notebooks still build identically
  (`python jepa/build_notebooks_v3.py` output byte-identical, except the import guard).
* A short "how to run" at the top of the notebook: open in Colab (badge), Run all, rerun after any
  disconnect, switch to a CPU runtime for analysis only if desired.

---

## 14. Implementation notes

(Filled in by the implementing session.)

---

## 15. Deviations

(Any change to §4 to §8 after the full run starts: what, why, date, and which results it affects.)
