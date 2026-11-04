# SpaceJEPA v5: strong baselines and a standard benchmark (plan and implementation brief)

Branch: `space-jepa5` (from `space-jepa4`). Status: **specified, not built; waiting for the user's review
(§0.2).** Written 2026-11-03, before any v4 result was analyzed, so nothing here was chosen after
seeing v4 numbers. This file is the preregistration of experiment 5 and the brief for whoever builds
it. The v1/v2 design is in `jepa/PLAN.md`, v3 in `jepa/v3_plan.md`, v4 in `jepa/v4_plan.md`; v5 reuses
their code and conventions (v4 §11) unless this file says otherwise.

---

## 0. How to use this file

### 0.1 Rules for the building session

1. Work only on `space-jepa5`. Push to `origin space-jepa5`. **Never push to `space-jepa4` while the
   v4 Colab run is going**: the v4 notebook clones that branch on every restart. Never push to
   `jepa-v3`.
2. Read `CLAUDE.md` (every Colab notebook releases its runtime by itself), then `v4_plan.md` §11 and the
   modules listed in §11 below, before writing code.
3. Build exactly what §3 to §12 specify. Where this plan leaves a choice open, make it and record it in
   §14, as v4 did.
4. **No simulator changes** (`rl_env.py`, `math_3d.py`, `control.py`, physics, RCS, existing files in
   `configs/`; the new `configs/benchmark_v5/` folder is data, not simulator). Every dispersion in §3
   maps onto existing knobs (`ic_override`, `density_scale`,
   `lift_scale`, `mass_scale`). The one optional suite that would need a simulator hook (§2, "shear")
   stays out unless the user approves it.
5. Do not change v3 or v4 behavior: new code goes in new `*_v5.py` modules, or behind new arguments
   whose defaults reproduce v4 exactly (proved in the smoke test, §12).
6. Do not run the full study. Build, smoke-test on CPU, generate the notebook, commit, push. The user
   runs it on Colab.
7. The test suites (§3) are frozen with checksums **before** any baseline is tuned. Tuning uses only
   the dev suite (§3.3). Criteria (§8) do not change after the first test flight; any change goes into
   §15 with the date and reason, before the affected results are looked at.
8. Persistence and resume (§9.1) are requirements, not polish: regular checkpoints, everything
   (results, figures, logs) on Google Drive as it is produced, and "Run all after a disconnect picks up
   where it stopped". The smoke test proves it (§12, checks 12–15).

### 0.2 Decisions for the user before building

| # | Decision | This plan's default | Alternative |
|---|---|---|---|
| D1 | Compute tier to build and run (§4.9) | **decided by the user (2026-11-03): Tier 1 (core) first, Tier 2 only after the user has read Tier 1's results.** The notebook's default settings run Tier 1 only | (rejected: a leaner Tier 1) |
| D2 | Bank authority for the published guidance laws | 90° (the same handle every learned controller has); 180° versions as exploratory arms | 180° as primary (unfair to the learned arms, whose training data never exceeds 90°) |
| D3 | Standard scenario ORB-1 (§3.2) | Lu (2014) CEV orbital entry state, adapted (γ0 = −1.40°, targets at 25/50/75% of footprint), Bairstow (2006) dispersions, Orion EM-1 accuracy requirement | Lu's exact 3,769 km range (infeasible in our simulator at a 90° bank limit, see §1.3) |
| D4 | EI speed dispersion for ORB-1 | 15 m/s (3σ), our choice: no published LEO value was found | another value the user has a source for |
| D5 | Density-shear suite (needs a simulator hook) | out | in, as Tier 3, with an explicit `configs/`-free hook |

---

## 1. Where the project stands, and what v5's research found

### 1.1 Context

v3 (5 seeds, 64 cases × 4 atmospheres, 4.5 g soft limit) found v2 lands within 5 km in 92–98% of
conditions against the classical controller's 6–91%. v4 (running) decomposes planner versus model,
matches load, and ablates the v1 → v2 changes. A publication-readiness review (2026-11-03, this
session) found the method strong and the **benchmarking weak**:

1. The classical baseline is our own constant-bank search (`control.SimpleBankGuidance`; its own
   docstring says "not a literal full reproduction of Lu"). It keeps the nominal model under mismatch:
   **no online density or L/D estimation**, which every modern predictor-corrector has. 51.6% success
   within 5 km in nominal air is low for this class of method.
2. No published guidance law is compared.
3. No learned dynamics model *without* the JEPA objective is compared (v1's "Predictive" baseline was
   one-step teacher-forced, a weak recipe; none of v4's 13 ablations removes the JEPA losses).
4. The RL baselines are strawmen: PPO from scratch ran as 150k-step diagnostics (1M suggested for a
   full run; ~1% success), behavior cloning used 60 expert episodes (91 km median), clone → PPO
   collapsed (438 km). `ReadMe.md` already names
   DAgger and a warm-critic handoff as the fixes.
5. The test cases are ours. No published scenario, dispersion model or requirement is used.

v5 fixes all five. It is **benchmarking only**: no change to v2, no new JEPA variants, no fixes.

### 1.2 Sources read for this plan

| Source | What v5 takes from it |
|---|---|
| Lu, P., "Entry Guidance: A Unified Method," *JGCD*, 2014, doi:10.2514/1.62605 (`re_entry_ai_research/lu-2014-entry-guidance-a-unified-method.pdf`) | FNPEG: bank magnitude linear in energy, σ0 by step-size-controlled Newton/secant on the range error, energy as the independent variable, integrate range-to-go (Eq. 17) rather than recompute it, heading-offset deadband reversals, guidance at 1 Hz only while aero acceleration ≥ 1.52 m/s², predictive load relief (§VI-B, Eqs. 35–42; for the CEV δ = 0.02 nondimensional = 16 s and k0 = 100, demonstrated at amax = 4 g and 3 g; ~0.1 g margin), the CEV case (Table 1: 121.9 km, 7,623.5 m/s, 2,035 nmi, final 7.6 km and 150 m/s, L/D ≈ 0.28) |
| Brunner, C., and Lu, P., "Skip Entry Trajectory Planning and Guidance," *JGCD* 31(5), 2008 (cited by Lu 2014 as ref. 4) | first-order fading-memory filters on L/L* and D/D* (measured over modeled lift and drag), used inside the predictor |
| Bairstow, S. H., *Reentry Guidance with Extended Range Capability for Low L/D Spacecraft*, MIT SM thesis, 2006 (`research/Reentry Guidance with Extended Range Capability for Low L-D.pdf`) | PredGuid: the complete algorithm in Appendix B (executive, targeting, initial roll, Huntest/constant drag, down/upcontrol, PredGuid predictor and corrector, ballistic, final phase), the final-phase reference-trajectory and gain design procedure (§5.5), the Monte Carlo dispersion table (Table C.1), several target ranges to exercise every logic branch (§3.2.1), duration-based crew g limits (NASA-STD-3000, Tables D.1–D.3), margin metrics (vertical L/D fraction, final-phase energy bucket) |
| McNamara, L. W., and Braun, R. D., AAS 14-094, Orion EM-1 entry design trajectory (NTRS 20140004224) | Orion flies PredGuid (numeric predictor-corrector + Apollo final phase) with "gain set and reference trajectories retuned or at least recalculated for varying L/D and entry range"; the **landing accuracy requirement: within 10 km at 99.865% (3σ) with 90% confidence, i.e. at most 1 miss in a 3,000-case Monte Carlo**; GRAM-2010 atmosphere dispersions (not available to us) |
| Putnam, Z., Neave, M., and Barton, G., "PredGuid entry guidance for Orion return from low Earth orbit" (not available offline) | for LEO return the phase sequence is Initial Roll → Final → Terminal; follow it if the paper can be obtained, otherwise Bairstow App. B's logic, which selects the same phases for short ranges |
| `StochasticEntrySim/colab/build_colab_notebook.py`, `ReadMe.md`, `docs/rl_mdp_spec.md` | the legacy RL setup and its artifacts (§1.4) |

### 1.3 Facts checked in our simulator for this plan

* Trim L/D of our capsule is 0.27–0.30 at Mach 20–40 (`physics_v4.TRIM_SCHEDULE`), mass 8,500 kg,
  reference area 19.63 m²: the same vehicle class as Lu's CEV (L/D ≈ 0.28, 8,383 kg, 5.5 m).
* Constant-bank sweeps with `physics_v4.rollout_physics` (nominal parameters) from Lu's entry state
  (121.9 km, 7,623.5 m/s), range at drogue altitude (7.6 km, which is exactly Lu's final altitude):

  | γ0 | bank 90° | bank 70° | bank 45° | bank 0° |
  |---|---|---|---|---|
  | −1.00° | 2,664 km | 2,896 km | 3,336 km | 3,817 km |
  | −1.25° | 2,418 | 2,640 | 3,068 | 3,538 |
  | −1.50° | 2,210 | 2,421 | 2,836 | 3,295 |
  | −2.00° | 1,878 | 2,068 | 2,457 | 2,893 |

  Lu's 3,769 km is reachable only at γ0 ≈ −1.0° with the lift vector almost fully up. That leaves no
  margin, and it is outside the world models' training envelope (at 7,623.5 m/s the shallowest trained
  γ is −1.07°). Our simulator has a non-rotating Earth and caps bank at 90°, Lu's does not. Hence D3:
  keep Lu's entry state, choose γ0 and target ranges that fit.
* World-model training envelope (`jepa/data.py`): entry altitude 100–140 km, speed 6,800–8,600 m/s,
  γ −6° to −1° (shallowest steepens linearly from −1° at 7,600 m/s to −4° at 8,600 m/s), any heading,
  latitude ±30°, longitude ±20°; truth randomization density ×0.70–1.30, lift ×0.70–1.30,
  mass ×0.90–1.10 (uniform). Bairstow's 3σ density dispersion (30%) is exactly this range.
* Every learned controller and its training data use |σ| ≤ 90° (`data.SIGMA_MAX_RAD = π/2`, "same
  handle as the RL env"). Lu and Apollo/PredGuid normally use up to 180° (lift down). Hence D2.
* The simulator's atmosphere is a 1-D table in altitude, and the Earth does not rotate, so the dynamics
  depend only on the geometry relative to the target, not on where the great circle lies.

### 1.4 Legacy RL artifacts (from the code; Google Drive is not mounted on this machine)

The RL notebook (`colab/ReEntryAI_RL_Colab.ipynb`, built by `colab/build_colab_notebook.py`) wrote
everything to **`MyDrive/ReEntryAI/`**: `data/` (`baseline_metrics.csv`, `policy_metrics.csv`,
`bc_dataset.npz` with 60 classical episodes, `bc_policy_metrics.csv`), `checkpoints/<RUN_NAME>/`
(`ppo_*_steps.zip`), `models/` (`<RUN_NAME>_final.zip`, `bc_orion_entry_v1.zip`,
`ppo_ft_from_bc_v1_final.zip`), `tb/` (TensorBoard), and Weights & Biases project `ReEntryAI`.
Run names in the code: `ppo_orion_entry_v4_diag` (an earlier "v3" reward version also ran; the code
suggests `..._v5_full` for a full run), `bc_orion_entry_v1`, `ppo_ft_from_bc_v1`. Policies are SB3
`MlpPolicy` [128, 128] on `rl_env.OrionEntryEnv` (13 normalized features, action a ∈ [−1, 1] → σ =
90°·a at 1 Hz). Evaluation set: `configs/heldout_testset.csv` (1,000 cases from the narrow default
dispersion of `configs/default.json`), where classical lands 99.9% (2.05 km median). The v5 notebook's
first stage inventories this folder (§4.6).

---

## 2. Scope

| In scope (v5) | Out of scope |
|---|---|
| Published guidance laws: FNPEG (Lu 2014) and PredGuid (Bairstow 2006 / Orion), tuned on a dev suite | changing v2 or any JEPA model |
| Online density and L/D estimation for physics-model planners (`mppi_filtered`) | navigation errors, sensor noise, RCS failures (v6) |
| Non-JEPA learned dynamics models (MLP, GRU, GRU + physics prober, R0 without JEPA losses) under the same MPPI | lunar-return or skip-entry scenarios (the world models are not trained there) |
| RL baselines done properly (DAgger, warm-critic BC → PPO, PPO and SAC with history), plus re-evaluation of the legacy policies | retuning MPPI or the world models |
| A frozen, checksummed benchmark (suites A256, ORB-1, DEV48) with a scoring module | density-shear atmospheres (D5, optional Tier 3, needs a simulator hook) |
| Requirement-style metrics: Orion EM-1 accuracy test, NASA-STD-3000 g-duration compliance, compute per decision | any claim about flight hardware |

---

## 3. The benchmark (test-case set)

All suites live in the repo under `StochasticEntrySim/configs/benchmark_v5/`, are generated by
`jepa/benchmark_v5.py` with fixed seeds, and are verified against `SHA256SUMS` every time they load.
A `README.md` there is the datasheet: purpose, generation code and seeds, column definitions, sources,
scoring, how to cite. The suites are copied to `MyDrive/ReEntryAI_jepa_v5/01_benchmark/` at setup.

### 3.1 A256 (continuity)

The v3 test set unchanged: `configs/v3_testset.csv` (64 cases) × `testset_v3.ATMOSPHERES` (nominal;
thin ×0.88; dense ×1.12; ×0.78 with L/D ×0.80) = 256 conditions. Every v3 and v4 flight pairs with v5
flights on it. `A256.json` records the file's checksum and the atmosphere definitions.

### 3.2 ORB-1 (standard scenario: Orion-class LEO return)

**Nominal mission.** Entry interface 121.92 km (400,000 ft: Lu 2014, Bairstow 2006, Orion), Earth-relative
speed 7,623.5 m/s (Lu 2014 Table 1, CEV orbital entry), γ0 = −1.40° (D3; Lu does not publish γ0; even
the shallowest dispersed case, −1.30° at 7,638.5 m/s, stays 0.18° inside the training envelope's
shallow edge), final condition drogue deploy at 7.6 km (Lu's final altitude). Entry point (0°, 0°),
heading χ0 = 45.99° (the A256 corridor, toward 15°N 15°E), targets on that great circle at **three
downranges: 25%, 50% and 75%** of the nominal footprint between 90° and 0° constant bank, measured in
the full simulator by the S0 gate (§4.0) (about 2,570, 2,840 and 3,120 km, interpolated from the
physics-model sweep in §1.3). Several ranges exercise every branch of the guidance logic (Bairstow
§3.2.1).

**Dispersions** (Gaussian, truncated at ±3σ, independent; Bairstow 2006 Table C.1 unless stated):

| Parameter | Nominal | 3σ | Note |
|---|---|---|---|
| EI flight-path angle | −1.40° | 0.10° | Bairstow (0.5 nmi vacuum periapsis error) |
| EI azimuth | toward the target | 0.10° | Bairstow (orbital inclination error) |
| EI speed | 7,623.5 m/s | 15 m/s | **D4**: Bairstow's 300 m/s is lunar-return geometry, not error; no published LEO value found |
| EI altitude, latitude, longitude | 121.92 km, 0°, 0° | 0 | EI defined by altitude |
| vehicle mass | 8,500 kg | 5% | Bairstow |
| density multiplier (constant with altitude) | 1 | 30% | Bairstow ("day of flight") |
| CD multiplier | 1 | 5% | Bairstow's hypersonic value applied at every Mach (the simulator has no Mach-dependent multiplier) |
| CL multiplier | 1 | 5% | same; also stands in for Bairstow's CG-offset dispersion (no CG in a 3DOF simulator) |

Mapping to simulator knobs (`data.apply_domain`): `density_scale = ρ_mult × c_D`,
`lift_scale = c_L / c_D`, `mass_scale = m_mult`, so lift ∝ ρ·S·CL·c_L and drag ∝ ρ·S·CD·c_D exactly.
Each case records `in_training_domain` (density_scale and lift_scale inside 0.70–1.30, mass_scale
inside 0.90–1.10; about 0.3% of cases fall just outside) and its footprint margin under its true
parameters (`testset_v3` footprint code); analyses report all cases and feasible-only.

**Size and blocks.** 3,000 cases (the EM-1 Monte Carlo size), seed 20261103, in **three blocks of
1,000** (each block 1/3 per target range). Learned arms: case *i* is flown by training seed *i* mod 5
(200 per seed per block), so the requirement test is for the model family (a randomly chosen training
run), and per-seed counts are reported. Deterministic arms fly every case once.

**Sequential stopping (saves compute, cannot change a pass).** Every ORB-1 arm flies block 1. After
each block, an arm that already has ≥ 2 flights that are not (outcome `drogue` and miss ≤ 10 km) has
failed Q1 (§8.4) and stops; the rest fly the next block. Comparisons between arms (Q3) use block 1
only, which every arm has.

### 3.3 DEV48 (tuning only, never a result)

12 new cases from the v3 test-set generator (`testset_v3`, same spec, new Sobol seed 20261104; any
Sobol index already in A256 is skipped) × the 4 A256 atmospheres = 48 conditions. Every tunable baseline
(FNPEG, PredGuid, RL) is tuned here and only here. Reported only in the tuning appendix.

### 3.4 Scoring (`benchmark_v5.score`)

Input: a flights table (one row per flight, v4's columns). Output: every metric of §5, the Q1
requirement test, g-duration compliance, per arm × suite × (atmosphere or target range).

---

## 4. Experiments

Arms are fixed strings, used in paths, tables and figures. All flights use the v4 lockstep runner and
v3's common random numbers (`fly_v3.crn_seed`), so every arm pairs with every other on the same
(seed or noise index, case, atmosphere, rep). Deterministic arms (no MPPI noise, no training seed) fly
once per condition and are broadcast across seeds (as `analyze_v3.cube`).

### 4.0 Gates (run first, gate everything else)

| Gate | What | Pass condition |
|---|---|---|
| S0 benchmark | build/verify the suites; ORB-1 footprint in the full simulator (90° and 0° constant bank from the nominal EI) defines the three target ranges; then the suite is frozen | checksums match; every target range strictly inside the nominal footprint; ≥ 97% of ORB-1 cases feasible under their true parameters; the longest nominal ORB-1 flight is shorter than MPPI's horizon (188 × 5 s = 940 s); nominal ORB-1 trajectories stay inside the training split's feature ranges |
| V1 FNPEG verification | on 16 A256 cases, nominal air, predictor = true parameters, no load limit | miss < 2 km in ≥ 15 of 16 (Lu reports < 1 km in his own simulator; ours adds roll dynamics the predictor does not model), Newton iteration converges in ≤ 6 iterations per cycle in ≥ 99% of cycles |
| V2 load relief | same 16 cases, amax = 3.0 g, cases whose unconstrained peak exceeds 3.5 g | realized peak ≤ amax + 0.15 g and miss < 2 km where feasible; as amax drops below feasibility the miss grows gradually (Lu Fig. 17), checked on 2 cases at 2.7/2.5/2.3 g |
| V3 filters | thin (×0.88) and dense (×1.12) air, 8 cases each | the filter's drag ratio K_D is within 3% of the true scale by 150 s after the 0.05 g point |
| V4 PredGuid verification | 16 A256 cases, nominal air, nominal model | miss < 3.5 km in ≥ 15 of 16 (Bairstow's Monte Carlo bound); the phase sequence is logged and sensible (initial roll → final for these ranges) |
| E0 runner | v5 runner vs `run_lockstep_v4` on CPU, v2 and physics adapters, 3 flights | bit-identical |

v4's E0b (statistical re-fly of v3) already covers GPU nondeterminism; v5 cites v4's floor.

### 4.1 E1: published guidance laws (RQ1)

| Arm | What | Suites | Flights |
|---|---|---|---|
| `fnpeg` | Lu (2014) baseline + predictive load relief at amax = 4.5 g (0.1 g margin inside) + Brunner–Lu fading-memory L and D filters in the predictor; \|σ\| ≤ 90° | A256, ORB-1 | 256 + ≤ 3,000 |
| `predguid` | Bairstow (2006) App. B with the final-phase reference and gains designed preflight, per case, from **nominal** parameters (§5.5 procedure; Orion recomputes them per mission); g limiter at 4.5 g; \|σ\| ≤ 90° (a lift-down L/D-fraction command saturates at 90°; record how often) | A256, ORB-1 | 256 + ≤ 3,000 |
| `classical` | v3's controller | ORB-1 (A256 exists in v3) | ≤ 3,000 |
| `fnpeg_180`, `predguid_180` | exploratory: the same laws with \|σ\| ≤ 180° (how much does the 90° handle cost?) | A256 | 2 × 256 |
| `fnpeg_oracle` | exploratory: FNPEG whose predictor knows the true parameters, filters off (ceiling of the law) | A256 | 256 |

**FNPEG specification** (Lu 2014 §III and §VI-B, Brunner–Lu 2008; anything not fixed here is tuned, §4.5):
* Predictor model: `physics_v4`'s port of the classical predictor's equations (non-rotating Earth,
  same atmosphere table and trim schedule), float64, batched on GPU, with the bank rate limit.
  Integrate in time (RK4, step tuned in {1, 2} s) and carry range-to-go s by Lu's Eq. (17),
  ṡ = −V cos γ / r, signed (undershoot vs overshoot).
* Bank profile: |σ(e)| = σ0 + (e − e0)/(ef − e0)·(σf − σ0), e = μ/r − V²/2. ef from hf = 7.6 km and
  Vf = the median speed at 7.6 km over DEV48 classical flights (record it; Lu's CEV used 150 m/s).
  Stop the prediction at h ≤ 7.6 km or e ≥ ef, whichever comes first; error z = s_end.
* Corrector: σ0 by Lu's step-size-controlled Newton step (Eq. 24) with a finite-difference
  derivative, secant (Eq. 25) after the first iteration, λ = 1/2^i halving until |z| decreases, stop
  on Eq. (23); warm start from the previous cycle; σ0 clipped to [0°, 90°] (or 180°).
* Lateral: Lu's heading-offset deadband, linear in V (two parameters, tuned); instantaneous sign flip
  in the prediction, the simulator's rate limit in flight.
* Load relief: Eqs. (36)–(42), applied **inside** every predicted trajectory and once more on the
  current navigation data (Lu's point 3), δ and k0 tuned around Lu's CEV values (0.02, i.e. 16 s, and
  100).
* Filters: K_D ← β·K_D + (1 − β)·D_meas/D_model, same for L, using the measured aerodynamic
  accelerations the v2 encoder also sees (`cD_meas`, `cL_meas` inputs); active from 0.05 g; the
  predictor multiplies modeled D and L by K_D and K_L along the whole prediction. β tuned.
* Cadence: 1 Hz while aero acceleration ≥ 1.52 m/s² (Lu), else hold the current bank. If 1 Hz is too
  slow on Colab even batched, use 2 s and record it in §15 before any test flight.

**PredGuid specification:** implement Bairstow App. B literally, including the guidance cycle it
specifies; design the final-phase reference trajectory and gain tables for our vehicle with the §5.5
procedure from a nominal trajectory to each target (the "preflight" step, nominal parameters only);
the PredGuid predictor uses the same `physics_v4` model as FNPEG and Bairstow's L/D estimator if the
appendix defines one (record). Log the active phase per guidance cycle.

### 4.2 E2: learned model vs explicit estimation, and load fairness (RQ1)

| Arm | What | Suites | Flights |
|---|---|---|---|
| `mppi_filtered` | v4's `mppi_nominal` with the physics model's density and lift scales set every decision from FNPEG's filters (density_scale = K_D, lift_scale = K_L / K_D, mass 1). Same planner and settings as v2. The direct competitor to "implicit system identification from history." | A256 (5 noise indices), ORB-1 | 1,280 + ≤ 3,000 |
| `v2single` | v3's v2 checkpoints, member 0 (the v3 `v2single` arm), 4.5 g | ORB-1 (A256 exists in v3) | ≤ 3,000 |
| `v2single@fnpeg_matched` | same model, per-condition soft g limit = `fnpeg`'s realized peak g on that condition (as v4's matched arm) | A256 | 1,280 |
| `mppi_truth` | exploratory ceiling (v4 arm) | ORB-1 block 1 | 1,000 |

### 4.3 E3: non-JEPA dynamics models under the same MPPI (RQ2)

Each trained for seeds 0–4 on the same 6,000 flights, split, normalization and inputs as R0 (v2's
features over a 120 s history, v2's action features), flown with v2's MPPI settings at 4.5 g on A256
(5 × 256 each). Parameter count within ±10% of R0's 2,454,531.

| Arm | Model | Isolates |
|---|---|---|
| `mlp_direct` | MLP over the flattened 120 s history window + next action features → normalized 5 s state change and aerodynamic g-load; rolled out on its own outputs | black-box feedforward dynamics |
| `gru_direct` | GRU encoder over the history, GRU cell driven by action features in the rollout, linear head → state change and g-load | black-box recurrent dynamics |
| `gru_pi` | `gru_direct`'s recurrent core feeding v2's `PIProberV2` (physics backbone + learned force corrections) | physics-informed decoder without the JEPA objective |
| `r0_nojepa` | R0's architecture exactly (`WorldModelV2` + `PIProberV2`), trained end to end on the prober's state loss only: no latent prediction loss, no SIGReg | the JEPA objective, architecture held fixed |

**Training (all four):** free-running multi-step rollouts (backpropagation through the model's own
predictions, never teacher forcing: v1's "Predictive" lost partly to teacher forcing, so this removes
that excuse), the prober's dimensionless state loss (`PLAN.md` §3.4: 1 km, 5 km, 50 m/s, 0.5°, 0.5°)
plus v2's auxiliary g-load term at its weight, R0's optimizer, schedule and AMP, 36,000 steps × 256
(the same number of samples as R0's two stages, 24,000 × 256 + 6,000 × 512) with the horizon
curriculum 24 → 48 → 96 at 20% and 50% (R0's stage 1). Learning rate from
{R0's, ×0.3, ×3} chosen on the **validation split's open-loop loss only**, never on flights. Resumable
as v4 (`train_state_*.pt`, RNG state). Open-loop error at 400–50 s to go on the test split, as v4 T10.

### 4.4 E4: RL baselines (RQ3)

**Common to every new RL arm.** Decision every 5 s (MPPI's cadence), action = bank command in
[−90°, 90°], held for 5 s through the simulator's actuator. Observation = `rl_env`'s 13 features +
the last 24 decisions (120 s, the world models' memory) of [h, V, γ, measured a_D, cD_meas, cL_meas,
σ], normalized. Reward = the negative of MPPI's own cost, realized: terminal −miss_km at the drogue
crossing (or at termination), −2 km/rad × |Δσ| per decision, −50 km/g² × max(0, g − 4.5)² per
decision, plus potential-based shaping Φ = −range-to-go (does not change the optimal policy). Training
initial conditions and truth parameters from the world models' training distribution (`EnvelopeSpec`,
`WIDE_DOMAIN`). Hyperparameters tuned on DEV48 only (≤ 12 configurations per method). Seeds 0–2.

**Fast training environment** (`jepa/rlenv_fast_v5.py`): `physics_v4`'s truth model batched on GPU
(thousands of environments), 5 s steps with (2, 2, 1) s substeps and the bank rate limit, domain
randomization as above, termination at drogue altitude, |cos γ| < 0.05 or skip-out. Every policy is
**evaluated in the full simulator** like every other arm; the fast-vs-full gap is reported per arm.

| Arm | Method | Tier |
|---|---|---|
| `rl_legacy_*` | every policy found in `MyDrive/ReEntryAI/models/`, re-evaluated with its own settings on `heldout_testset.csv` (first 256 cases, to reproduce the README numbers) and on A256 (out of its training distribution; labelled) | 1 |
| `dagger_fnpeg` | DAgger with `fnpeg` as the expert (queried on the learner's own states, including its filter state), 10 iterations, β_i = 0.5^i, 200 rollouts per iteration in the fast env; policy outputs \|σ\| (regression) and sign (classification), so reversals are not averaged away (the README's diagnosis of BC) | 1 |
| `bcppo_fnpeg` | BC from 2,000 `fnpeg` episodes in the fast env, then critic warm-up on frozen-actor rollouts, then PPO with a KL penalty to the clone and a small learning rate (the README's "proper clone → RL handoff") | 1 |
| `ppo_hist` | PPO from scratch (in-house, CleanRL-style, torch), fast env, up to 200M decisions or 3 GPU-hours per seed | 1 |
| `sac_hist` | SAC, fast env, up to 20M decisions or 3 GPU-hours per seed | 2 |
| `ppo_matched` | PPO in the full simulator with the world models' data budget (6,000 episodes): the sample-efficiency point | 2 |

### 4.5 Tuning protocol (DEV48 only, frozen before any test flight)

Objective, in order: success@5km, then mean log miss, subject to mean peak g ≤ 4.6 g. Budget ≤ 24
configurations for FNPEG (σf, deadband slope and intercept, δ, k0, β, RK4 step) and for PredGuid
(deadband, initial-roll logic constants, g limiter, anything App. B leaves to the designer), ≤ 12 per
RL method. Random search from documented ranges, seed fixed. The chosen configuration, all tried
configurations and their DEV48 scores go to `03_tuning/` and §14. MPPI and the world models are **not**
tuned (they were designed in v1–v3).

### 4.6 Legacy RL inventory (first stage of the notebook)

List `MyDrive/ReEntryAI/{models,checkpoints,data}` and the W&B project if a key is set, write
`00_study/legacy_rl_inventory.csv` (path, run name, steps, size, date), then evaluate every final
model as `rl_legacy_<name>` (§4.4). Missing folder: record "not found" and continue.

### 4.7 Open loop for the new models

As v4 §4.3 open loop: footprint error at 400, 300, 200, 100, 50 s to go on the test split, for every
E3 model and seed, against R0 (v4) and v1 (v3).

### 4.8 Run order (most important results first)

Gates (§4.0) → legacy RL inventory → tuning (§4.5) → **A256**: `fnpeg`, `predguid`, `mppi_filtered`,
`v2single@fnpeg_matched` → **ORB-1 block 1**: `classical`, `fnpeg`, `predguid`, `mppi_filtered`,
`v2single` → E3 training and flights → E4 training and flights → ORB-1 blocks 2–3 (arms still
eligible) → exploratory arms → open loop → analysis → report. Within an arm, seed-major.

### 4.9 Compute tiers (D1)

Rates must be measured, not assumed: the status cell measures each arm type's flights per hour on the
first 20 flights and prints projected hours and Colab units (at the runtime's rate). Rough sizes:

| Tier | Contents | Flights (new) | Training |
|---|---|---|---|
| 1 | gates, legacy RL eval, tuning (2 × 24 × 48 guidance flights + RL tuning), the A256 arms of §4.1–4.3 except the exploratory ones, `dagger_fnpeg`, `bcppo_fnpeg`, `ppo_hist` (3 seeds each), ORB-1 block 1 for `classical`, `fnpeg`, `predguid`, `mppi_filtered`, `v2single` | ≈ 21,000 (≈ 6,000 of them cheap policy flights) | 20 dynamics models + 9 RL runs |
| 2 | ORB-1 blocks 2–3 (eligible arms only), ORB-1 block 1 for the non-JEPA and RL arm with the highest A256 success@5km (descriptive), `sac_hist`, `ppo_matched`, `mppi_truth` on ORB-1, the exploratory arms | ≈ 5,000–9,000 | 6 RL runs |
| 3 | density-shear suite (D5) | decided with the user | none |

For scale: v4 is ~27,700 flights plus 70 trainings. `classical` costs ~90 CPU-seconds per flight (v3:
6.4 CPU-hours for 256), so its ORB-1 block is ~25 CPU-hours: run CPU-only stages on a CPU runtime with
many cores, not on a GPU runtime. Per D1 the notebook runs Tier 1 by default (`TIER = 1`); Tier 2
stages exist, are built and smoke-tested, but stay switched off until the user sets `TIER = 2` after
reading Tier 1's results.

---

## 5. Metrics

As v4 §5 (success@5km headline, success@2km, log miss, paired log-miss difference, peak g, heat rate,
heat load, RCS fuel, reversals, planning ms), plus:

| Metric | Definition |
|---|---|
| success@10km | outcome `drogue` and miss ≤ 10 km (the EM-1 requirement radius) |
| requirement LCB | one-sided 90% Clopper–Pearson lower bound on P(success@10km); passes at ≥ 0.99865 (≈ ≤ 1 miss in 3,000) |
| g-duration compliance | the flight's g history stays under the NASA-STD-3000 deconditioned-crew curve (Bairstow Table D.3: 14 g for 0.5 s, 10 g for 8 s, 7.8 g for 30 s, 4 g for 100 s, 4 g for 1,000 s; log-log interpolation, Bairstow App. D.2 procedure with 0.2 g bins) |
| authority saturation | fraction of guided time with \|σ\| at 0° or at the bank limit (Bairstow's vertical L/D fraction margin) |
| compute per decision | median and p99 ms per guidance decision, CPU single thread and GPU, per arm, on the same Colab runtime; parameter count |
| fast-vs-full gap (RL) | success@5km in the fast env minus in the full simulator, same policy and cases |

---

## 6. Statistics (fixed before any v5 flight)

* A256: exactly v4 §6 (hierarchical bootstrap 10,000 resamples over seeds then cases within
  atmospheres, 95% percentile intervals, per-seed McNemar and sign tests, Holm within each family, TOST
  with 90% intervals).
* ORB-1: each case is one unit; learned arms have one seed per case, so bootstrap cases within seed
  blocks; Q1 uses exact Clopper–Pearson; Q3 uses block 1 paired by case.
* RL arms: seeds 0–2; "holds in a seed" needs 3 of 3.
* Every estimate is reported with its interval, including failures. Secondary metrics are exploratory.

---

## 7. Research questions

* **RQ1** Does v2 still beat classical guidance when the classical side is a published, tuned,
  load-limited, model-adapting law (FNPEG, PredGuid), and does the learned model beat explicit online
  estimation under the same planner?
* **RQ2** Is the JEPA objective needed, or does any learned dynamics model trained on the same data,
  with the same planner, do as well?
* **RQ3** How close do properly built RL baselines get to v2 and to FNPEG?
* **RQ4** On a standard Orion-class scenario with published dispersions, which controllers meet an
  Orion-style accuracy requirement and crew g-duration limits?
* **RQ5** What does each controller cost per decision? (descriptive)

---

## 8. Preregistered criteria

### 8.1 RQ1 (A256, 4.5 g; Holm family {G1, G2, G4, G5})

| # | Question | Passes if |
|---|---|---|
| G1 | v2 beats FNPEG | success@5km `v2single − fnpeg` > 0 pooled, 95% CI excludes 0, Holm p < 0.05, holds in ≥ 4 of 5 seeds, and no atmosphere where `fnpeg` is significantly better |
| G2 | v2 beats PredGuid | same with `predguid` |
| G3 | validity of the load comparison | `fnpeg` mean peak g ≤ 4.6 g in every atmosphere, and `v2single@fnpeg_matched` mean peak g ≤ `fnpeg`'s + 0.1 g in every atmosphere (if not, G4 is reported but flagged invalid) |
| G4 | v2 beats FNPEG at FNPEG's own load | success@5km `v2single@fnpeg_matched − fnpeg`, as G1 |
| G5 | learned model beats explicit estimation | success@5km `v2single − mppi_filtered`, as G1 |
| G6 | baseline competence (a check, not a hypothesis) | V1–V4 passed before test flights; reported: `fnpeg` and `predguid` vs `classical` per atmosphere |
| G7 | what the 90° handle costs (exploratory) | `fnpeg_180 − fnpeg`, `predguid_180 − predguid`: estimates only |

### 8.2 RQ2 (A256; paired log miss vs R0, v4 flights; Holm over the 4 arms)

| # | Question | Passes if |
|---|---|---|
| J-each | R0 beats this model | mean paired log-miss difference (arm − R0) > 0, 95% CI excludes 0, Holm p < 0.05, ≥ 4 of 5 seeds |
| J-equiv | the JEPA objective does not matter (`r0_nojepa` only) | 90% CI of the mean paired log-miss difference inside ±0.15 (TOST). If neither J-each nor J-equiv passes for `r0_nojepa`, the answer is "inconclusive" |
| J-mech | mechanism | Spearman ρ across R0, v1 and the 4 models between open-loop error at 300 s and closed-loop mean log miss, with a bootstrap CI over seeds (estimate) |

### 8.3 RQ3 (A256; Holm over the Tier 1 RL arms)

| # | Question | Passes if |
|---|---|---|
| R-each | v2 beats this RL arm | success@5km `v2single − arm` > 0 pooled, 95% CI excludes 0, Holm p < 0.05, holds in 3 of 3 RL seeds |
| R-imit | imitation gap | `dagger_fnpeg − fnpeg` success@5km: estimate |
| R-gap | sim-to-sim gap | fast-vs-full gap per arm: estimate |

### 8.4 RQ4 (ORB-1)

| # | Question | Passes if |
|---|---|---|
| Q1 | accuracy requirement, per arm | requirement LCB ≥ 0.99865 over the flights the arm flew (descriptive pass/fail, no p-value) |
| Q2 | crew g limits, per arm | fraction of g-duration-compliant flights with a Clopper–Pearson interval (descriptive) |
| Q3 | v2 vs the published laws on the standard scenario (Holm over 2) | success@5km `v2single − fnpeg` and `v2single − predguid` on block 1, > 0 with 95% CI excluding 0 and Holm p < 0.05 |

### 8.5 RQ5

Descriptive table (§5 compute per decision, parameter counts).

---

## 9. The notebook: `ReEntryAI_JEPA_v5.ipynb` (one notebook, generated)

Built by `jepa/build_notebook_v5.py` (reuse `AUTO_DISCONNECT_CODE`, `DONE_MD`, `DONE` from
`build_notebooks_v3.py` and v4's Settings-cell pattern: `MODE`, `RUN` per stage, `AUTO_DISCONNECT`,
`CPU_WORKERS`, `GPU_WORKERS`). Colab metadata: GPU, L4. **Resume rule as v4: after any disconnect, Run
all continues where it stopped**; nothing finished is redone. CPU-only runtimes skip GPU stages with a
message, so analysis and report can run cheaply.

| # | Section | Resume behavior |
|---|---|---|
| 0 | Settings: `MODE`, `TIER = 1` (D1), `RUN` per stage (derived from `TIER`, overridable), `AUTO_DISCONNECT`, worker counts | n/a |
| 1 | Setup: clone/pull `space-jepa5`, `pip install pymsis reportlab` (+ `stable-baselines3` only for legacy RL), mount Drive, check v3/v4 folders and the frozen settings, freeze `00_study/study_config.json`, auto-disconnect hook | idempotent |
| 2 | Status: done/total per stage and arm, measured rates, projected hours and units | fast |
| 3 | Benchmark: verify suites and checksums, S0 gate | skips if frozen |
| 4 | Legacy RL inventory and evaluation | per-file cache |
| 5 | Gates V1–V4, E0 | skips passed gates; stops the notebook on a failure |
| 6 | Tuning (DEV48) | per-configuration results cached |
| 7 | Tier 1 flights: A256 baselines, `mppi_filtered`, matched arm | one JSON + `.npz` per flight, as v4 |
| 8 | ORB-1 blocks with the stopping rule | per flight; block eligibility recomputed from Drive |
| 9 | E3 training | `train_state_*.pt` every 5 min, as v4 |
| 10 | E3 flights | as 7 |
| 11 | E4 RL training | checkpoints every 5 min with optimizer and RNG state |
| 12 | E4 RL flights | as 7 |
| 13 | Exploratory arms | as 7 |
| 14 | Open loop | per-file cache |
| 15 | Analysis | bootstrap cache keyed by comparison + input hash (v4) |
| 16 | Report: `06_report/SpaceJEPA_v5_report.pdf` | rebuilt each run |
| 17 | Done: `release_runtime("all cells finished")` | n/a |

Drive layout `MyDrive/ReEntryAI_jepa_v5/` (smoke: `ReEntryAI_jepa_v5_smoke/`): `00_study/`,
`01_benchmark/`, `01_checks/` (gates), `02_models/<arm>/seed_<k>/`, `02_policies/<arm>/seed_<k>/`,
`03_tuning/`, `03_open_loop/`, `04_flights/results/<suite>/<arm>/<load>/seed<k>/mppi<m>/<cond>.json`
(+ `telemetry/`, `flights.csv`), `05_analysis/`, `06_report/`. Read only, guarded as v4:
`ReEntryAI_jepa_v3/`, `ReEntryAI_jepa_v4/`, `ReEntryAI_jepa_data/`, `ReEntryAI/` (legacy RL).

### 9.1 Persistence, checkpointing and resume (requirements; user request 2026-11-03)

**The contract.** After any disconnect, crash, idle timeout, runtime release or kernel restart, the
user opens the notebook and presses **Runtime → Run all**. Finished work is skipped in seconds,
half-done work continues from its last checkpoint, nothing finished is redone, and no cell needs
editing. The only manual step Colab may ask for is re-authorizing the Drive mount. The notebook's first
markdown cell says this in one sentence, next to the `AUTO_DISCONNECT = False` switch (`CLAUDE.md`).

**Rules (every stage, every module):**

1. **Drive is the only source of truth.** Every result, checkpoint, log, table and figure goes to
   `MyDrive/ReEntryAI_jepa_v5/` as soon as it exists. The VM's disk holds only caches that can be
   rebuilt (the copied training data, the atmosphere table, a checkpoint staged for upload).
2. **Atomic writes.** Write `<name>.tmp` next to the target, flush and `fsync`, check the size, then
   `os.replace` (v4's `atomic_bytes`). Large files (checkpoints, buffer shards) are first written to
   local disk, then copied to the Drive `.tmp` name, size-checked, then renamed: Drive's FUSE layer is
   slow with many small writes to one big file and can stall on a large one. A disconnect therefore
   never leaves a half-written file under a final name.
3. **Done markers.** A unit of work is done only when its final artifact exists: `ckpt.pt` for a
   training run, the result `.json` for a flight, `config_<j>.json` for a tuning configuration,
   `iter_<i>.npz` for a DAgger iteration, a `.done` file for a gate or an analysis comparison. Every
   stage starts by listing what is done (with `os.scandir`, never a full tree walk) and does the rest.
4. **Checkpoint cadence and contents:**

   | Work | Saved to Drive | When | Contents | On resume |
   |---|---|---|---|---|
   | E3 model training | `02_models/<arm>/seed_<k>/train_state.pt` and `train_state.prev.pt` | every 5 min of wall time, at every curriculum switch, and at the end | model, optimizer, AMP scaler, LR schedule, step, curriculum stage, RNG states (Python, NumPy, torch CPU and CUDA), data-sampler state, loss history | load the newest state that loads cleanly (else `.prev`), continue at the same step; `ckpt.pt` is written last and marks done; then the train states are deleted (as v4) |
   | RL training (PPO, BC → PPO, SAC, DAgger's learner) | `02_policies/<arm>/seed_<k>/train_state.pt` and `.prev` | every 5 min and at the end | actor, critic(s), optimizers, observation-normalizer statistics, global step, update count, environment RNG and episode state, learning-curve CSV | as above |
   | SAC replay buffer | `02_policies/sac_hist/seed_<k>/buffer/shard_<n>.npz` | with each checkpoint: only the transitions added since the previous one (append-only shards) | transitions | reload the newest shards up to the buffer's capacity |
   | DAgger data | `02_policies/dagger_fnpeg/seed_<k>/data/iter_<i>.npz` | after each iteration's rollouts are labeled | states, expert labels, iteration | skip finished iterations |
   | Tuning | `03_tuning/<method>/config_<j>.json` | after a configuration's DEV48 flights | configuration, scores | skip scored configurations; the flights inside are per-flight files |
   | Flights (every suite) | one `.json` + `.npz` per flight | the moment the flight lands | as v4 | finished flights skipped; flights in progress re-flown on the same random-number stream |
   | ORB-1 blocks | per flight | — | — | block eligibility recomputed from the files |
   | Open loop | one CSV per (model, seed) | when done | errors by time to go | skip |
   | Analysis | `05_analysis/cache/` keyed by comparison + input hash; each table and figure | as soon as each is computed | — | recompute only what changed |
   | Report | `06_report/SpaceJEPA_v5_report.pdf` | every run | — | rebuilt (cheap) |

5. **Outputs, not just checkpoints.** Every figure is saved as PNG (200 dpi) and PDF, with its
   `captions.json` entry, the moment it is made. Training curves (`F_curves.png`) and RL learning
   curves are redrawn at every checkpoint, so progress can be inspected on Drive mid-run. Every stage
   appends its progress lines (the same lines the notebook prints) to `00_study/logs/<stage>.log`:
   Colab's cell output is lost on a disconnect, the log is not. The status cell rewrites
   `00_study/status.json` and appends a row to `00_study/status_history.csv` (done/total, rate, hours
   left per stage).
6. **Cells do not depend on each other's memory.** Only Settings and Setup define shared state (paths,
   frozen config). Every stage cell reloads what it needs from Drive, so Run all after a disconnect and
   running one section by itself after Setup both work.
7. **Frozen config and code version.** As v4: `00_study/study_config.json` frozen at first use;
   Setup raises if a non-runtime setting differs. Every artifact's JSON records the git commit that
   made it. A code change mid-study that could change results is logged in §15 before continuing.
8. **Self-repair at Setup.** Delete stale `*.tmp` files; a checkpoint that fails to load falls back to
   `.prev` (logged); a flight `.json` that fails to parse is deleted and re-flown; a `.done` marker
   whose artifact is missing is removed.
9. **Bounded loss.** A disconnect loses at most 5 minutes of any training, the flights that were in
   the air, one tuning configuration's unfinished flights, and one analysis comparison.
10. **Release only after flushing.** `release_runtime` (from `CLAUDE.md`) calls
    `drive.flush_and_unmount()` before `runtime.unassign()`, on success and on failure.

---

## 10. Tables, figures, and the report

Naming as v3/v4 (`T##_*.csv`, `F##_*.png`, `captions.json`). Minimum set:

| # | Content |
|---|---|
| T01 | every metric per arm × suite × atmosphere/target range (v3, v4 and v5 arms) |
| T02 | T01 on feasible conditions only |
| T03 | gates S0, V1–V4, E0 |
| T04 | RQ1 paired differences (G1–G7) per atmosphere and pooled |
| T05 | criteria table (every G, J, R, Q) with value, CI, Holm p, met |
| T06 | RQ2: paired log miss vs R0, open-loop errors, parameter counts |
| T07 | RQ3: RL arms, success, miss, fast-vs-full gap, training budget used |
| T08 | RQ4: ORB-1 per arm: blocks flown, success@2/5/10km, requirement LCB, g compliance, per target range |
| T09 | RQ5: compute per decision and parameters |
| T10 | tuning appendix (DEV48) |
| T11 | legacy RL inventory and re-evaluation |
| F01 | success@5km by arm and atmosphere (A256) with seed whiskers |
| F02 | success vs mean peak g: v2 load curve (v4), `fnpeg`, `predguid`, matched arm |
| F03 | RQ2 forest plot (paired log miss vs R0) |
| F04 | open-loop error vs closed-loop log miss (J-mech) |
| F05 | RL learning curves (fast env) with full-simulator checkpoints |
| F06 | ORB-1 landing scatter (downrange vs crossrange error) per arm, 10 km circle |
| F07 | ORB-1 miss CDFs per arm with the 10 km line |
| F08 | g-duration plot (Bairstow App. D style) with the NASA-STD-3000 curve, per arm |
| F09 | success vs threshold curves (A256 and ORB-1) |
| F10 | example trajectories and bank histories: `fnpeg`, `predguid`, `v2single` on the same 3 conditions |

Report as v4 (`report_v4.py` pattern: reads only these files). Sections: summary callout, RQ1–RQ5,
benchmark datasheet, methods (incl. the exact FNPEG and PredGuid settings), deviations (§15),
appendix with every table.

---

## 11. Code map

### 11.1 Reuse (read first)

v4 §11.1 modules, plus `jepa/study_v4.py` (`LayoutV4` guard, atomic writes, job queue, status),
`jepa/fly_v4.py` (`run_lockstep_v4`, adapters, decision log), `jepa/physics_v4.py` (EOM port,
`PhysicsAdapter`), `jepa/analyze_v4.py` (`Boot`, criteria helpers, figures), `jepa/report_v4.py`,
`jepa/ablate_v4.py` (resumable training pattern), `jepa/testset_v3.py` (case generation, footprint
margin), `rl_env.py` (read only), `colab/build_colab_notebook.py` (legacy env settings).

### 11.2 New modules

| Module | Contents |
|---|---|
| `jepa/benchmark_v5.py` | suite generation (A256 pointer, ORB-1, DEV48), checksums, loading, `score` (§3.4, §5) |
| `jepa/guidance_fnpeg.py` | batched FNPEG (predictor, corrector, lateral logic, load relief), fading-memory filters |
| `jepa/guidance_predguid.py` | Bairstow App. B phases, PredGuid predictor/corrector, final-phase reference and gain design |
| `jepa/fly_v5.py` | 1 Hz lockstep runner for guidance laws (all flights in a batch step 1 s at a time), policy runner for RL arms, `FilteredPhysicsAdapter` (`mppi_filtered`), adapters for the E3 models, matched-g limits from `fnpeg` flights |
| `jepa/models_v5.py` | `mlp_direct`, `gru_direct`, `gru_pi`, `r0_nojepa`; training (resumable); open loop |
| `jepa/rlenv_fast_v5.py` | the batched fast environment |
| `jepa/rl_v5.py` | PPO, SAC (in-house, torch), DAgger, BC → critic warm-up → PPO with KL; legacy SB3 loading |
| `jepa/study_v5.py`, `jepa/analyze_v5.py`, `jepa/report_v5.py`, `jepa/smoke_v5.py`, `jepa/build_notebook_v5.py` | as their v4 counterparts |

### 11.3 Changes to existing modules (additive only)

None expected. If one is needed (for example a hook in `run_lockstep_v4` for a 1 s cadence), it gets a
new keyword argument whose default reproduces v4 bit for bit, proved in the smoke test.

---

## 12. Smoke test (`python jepa/smoke_v5.py`, CPU, numbers meaningless)

On the v4 smoke fixtures (`smoke_v4.py --work`, default `<tempdir>/jepa_v4_smoke/`: the ~200-flight
wide dataset, `v3_smoke/`, `v4_smoke/`; build them with `smoke_v4.py` first if they are absent):

1. Suites: generated twice give identical bytes; checksum mismatch raises; ORB-1 mapping to simulator
   knobs reproduces drag and lift scaling exactly on 100 random draws.
2. FNPEG: on a toy case with predictor = truth, converges and lands < 1 km; batched equals unbatched;
   load relief keeps a 3 g case under 3.15 g; filters converge to a known density scale.
3. PredGuid: phase sequence on one short and one long case; reference/gain design reproducible.
4. `mppi_filtered` with filters frozen at 1 equals `mppi_nominal` exactly; v5 runner equals
   `run_lockstep_v4` bit for bit for v2 and physics adapters (E0).
5. Every E3 model builds at full size within ±10% of R0's parameters (no training), trains a few
   steps, saves, resumes after an interrupted save with the step count continuing.
6. Fast env: matches `physics_v4.rollout_physics` exactly for the same commands; a random policy runs
   10,000 batched steps; PPO, SAC, DAgger and BC → PPO each run a few iterations and checkpoint/resume.
7. Legacy RL: a dummy SB3 zip is found by the inventory and evaluated.
8. Job queue: finished flights skipped, a deleted result re-flown, ORB-1 stopping rule applied from
   files on disk.
9. Analysis on smoke data: every table and figure produced, criteria table complete, bootstrap cache
   reused; report PDF builds.
10. The notebook executed locally in smoke mode runs end to end; a second run finds everything done in
    under a minute. v5 refuses to write under the v3, v4, data and legacy RL folders.
11. v3 and v4 notebooks still build identically.
12. **Resume equivalence (CPU, bit-exact).** For one E3 model, PPO, SAC (with its buffer shards) and
    DAgger: train N steps uninterrupted; separately train, kill the process right after a checkpoint
    at about N/2, restart, finish. Final weights, optimizer state and learning curves are identical.
13. **Interrupted writes.** Leave a stale `.tmp`, truncate the newest `train_state.pt`, corrupt one
    flight `.json`, delete one artifact behind a `.done` marker: Setup repairs each case (§9.1 rule 8)
    and the run finishes with the same outputs as an untouched run.
14. **Simulated disconnect in the notebook.** Execute the smoke notebook with a forced failure in the
    middle of each long stage (an environment variable such as `JEPA_V5_FAIL_AT=<stage>:<n>`), then
    execute it again unchanged: it completes, finished flights keep their file modification times, and
    the outputs equal an uninterrupted smoke run's (flights and training bit-exact on CPU).
15. **Outputs on Drive.** After the smoke run, every table and figure of §10 (PNG and PDF), the
    curves of every trained model and policy, `captions.json`, the stage logs, `status.json`,
    `status_history.csv` and the report exist in the smoke study folder.

---

## 13. Definition of done (for the building session)

* §11.2 modules, the notebook and its builder, `configs/benchmark_v5/` with `README.md` and
  `SHA256SUMS`, committed on `space-jepa5` and pushed to `origin space-jepa5`.
* `smoke_v5.py` passes, including the resume checks 12–15; its output is pasted into §14.
* §14 lists every choice made (FNPEG and PredGuid details and tuning ranges, Vf, model widths, RL
  hyperparameter ranges, cadence, anything else) and the measured smoke throughput per arm type.
* v3 and v4 plans, code paths and notebooks unchanged.
* A short "how to run" at the top of the notebook (badge, Run all, rerun after a disconnect, tiers).

---

## 14. Implementation notes

Built 2026-11-03/04 on branch `space-jepa5` (commits 981c7c1 and the one that adds this section).
Everything below is a choice the plan left open; forced changes are in §15.

**Benchmark (`configs/benchmark_v5/`, `jepa/benchmark_v5.py`).**

* ORB-1 footprint (S0, nominal EI, full simulator, constant command): 0° bank lands at 3,379.2 km
  (peak 2.46 g, 749 s); 90° bank ends at 2,273.3 km by the simulator's |cos γ| < 0.05 gate (peak 9.4 g,
  438 s), used as the near edge. Targets at 25/50/75 %: **2,549.8, 2,826.2, 3,102.7 km** (rounded to
  0.1 km; `orb1_target{0,1,2}.json`).
* Feasibility (constant bank, true parameters): 1,002/1,002, 999/999 and 999/999 per target (S0's
  ≥ 97 % holds). At ≤ 4.5 g with constant bank (`feasible_4p5g`; the largest constant bank under
  4.5 g is about 60°): **only 27/1,002 for target 0**, 999/999 for targets 1 and 2. Whether the arms reach target 0
  under 4.5 g with their modulated bank histories is part of what ORB-1 measures; the 4.5 g
  constant-bank flag is stored per case and not used as a filter. **Open for the user** (see the
  building session's summary): keep target 0 as is, or move it.
* ORB-1 truth knobs: drag ∝ ρ·S·C_D and lift ∝ ρ·S·C_L map exactly onto `data.apply_domain`
  (`density_scale = ρ·c_D`, `lift_scale = c_L/c_D`, `mass_scale = m`); smoke check 1 proves it on the
  simulator's own aero (100 draws, relative error < 1e-12).
* Draw order per case: FPA, azimuth, speed, mass, density, C_D, C_L, each a 3σ-truncated normal by
  rejection, one `default_rng(20261103)` stream; case i uses training seed i mod 5.
* DEV48: `testset_v3`'s candidate filter, Sobol seed 20261104, A256's Sobol indices skipped, 12 cases ×
  4 atmospheres. Files are written with `%.6f` floats and LF line endings (`.gitattributes -text`), so
  the SHA-256 sums hold on Windows and Linux checkouts.

**FNPEG (`jepa/guidance_fnpeg.py`).** Lu (2014) as in the module docstring. Choices:

* Predictor: `physics_v4`'s longitudinal port of the simulator (non-rotating Earth, the simulator's
  atmosphere table and trim schedule), RK4 at `dt_s = 2 s` (tuned in {1, 2}), float64 numpy,
  vectorized over every prediction of a corrector round across flights; the bank follows the profile
  through the 20°/s rate limit inside the prediction.
* σ profile: constant (σ_f = σ0) by default; linear in energy when tuning draws a σ_f.
* `Vf = 150 m/s` until `03_tuning/fnpeg/vf.json` (median drogue speed of DEV48 classical flights)
  replaces it; tuning reads it before drawing configurations.
* Lateral deadband `1° + 5.2e-4 °/(m/s) · V` (1° at rest, 5° at 7,700 m/s), tuned.
* Load relief: Lu's Eqs. 36–41 inside every prediction and once on navigation data; δ = 0.02, k0 = 100
  (nondimensional), amax 4.5 g with a 0.1 g margin.
* Filters: fading memory, β = 0.9 at 1 Hz, on from 0.05 g; `fnpeg_oracle` gets the true K_D, K_L.
* Guidance cycle 1 Hz while the aerodynamic acceleration ≥ 1.52 m/s² (Lu), command held below it and
  for the last 10 s of predicted time to go.
* Corrector: Newton with a finite-difference first step (1°), secant after, step halving (≤ 5), at most 8
  iterations, stop on |z| ≤ 50 m or Eq. 23 with ε = 5 km²/rad.

**PredGuid (`jepa/guidance_predguid.py`).** Bairstow (2006) App. B, every Table B.4/B.7 constant except:
`GMAX = 4.5 g`; vehicle constants from our capsule (LAD = 0.28, C_D estimate 1.55, mass and area);
non-rotating Earth (WIE = 0). Guidance cycle 2 s (Apollo). Final-phase design (§5.5) per flight from
its own entry state with nominal parameters (A256/DEV48) or from ORB-1's nominal EI (ORB-1, the
mission's preflight design), cached per (IC, target, config). LOD by root finding (Bairstow's fixed
0.6·LAD available as `lod_fixed_frac`). Up-control: low-loft variant (B.50–B.53).

**Tuning ranges (§4.5; random search, seed 20261104, configuration 0 = the defaults).**

| Method | Ranges |
|---|---|
| FNPEG (24) | σ_f: none (25 %) or U(20°, 80°); deadband c0 U(0.5°, 3°), c1 U(2e-4, 1e-3) °/(m/s); δ U(0.01, 0.04); k0 logU(50, 200); β U(0.8, 0.97); RK4 step {1, 2} s |
| PredGuid (24) | LATBIAS U(0.2, 0.8); KLAT logU(0.006, 0.025); KA U(45, 85); Q7F U(3, 9); Q19 U(0.1, 0.4); overcontrol U(2, 6); K_ρ filter gain U(0.02, 0.10); GMAX U(4.2, 4.5) g |
| ppo_hist (6) | lr logU(1e-4, 1e-3); n_steps {32, 64, 128}; entropy {0, 1e-3, 1e-2}; GAE λ U(0.9, 0.98); clip {0.1, 0.2, 0.3} |
| bcppo_fnpeg (4) | lr logU(1e-5, 1e-4); KL coefficient logU(0.01, 1); critic warm-up {10, 20, 50} updates |
| dagger_fnpeg (4) | lr logU(1e-4, 3e-3); epochs {20, 30, 50}; hidden {256², 512²} |
| sac_hist (4, Tier 2) | lr logU(1e-4, 1e-3); τ {0.005, 0.01}; batch {512, 1024, 2048} |

Objective as §4.5 (success@5km, then mean log miss, mean peak g ≤ 4.6 g; if none qualifies, the
lowest mean peak g). RL tuning runs at 10 % of the budget on the fast env and is scored by DEV48
full-simulator flights; RL tuning starts only once FNPEG is chosen (the expert of two methods).
Untuned guidance flights are never written as test flights (`untuned` jobs wait).

**E3 models (`jepa/models_v5.py`).** Each arm's width is the one whose parameter count is closest to
R0's (bisection): R0 has 2,454,531 parameters; `mlp_direct` 4 hidden layers of width 831 over the flattened 120 s window (2,454,781, +0.01 %), `gru_direct` 2-layer GRU, hidden 449 (2,454,690, +0.01 %), `gru_pi` 2-layer GRU, hidden 441, with v2's state-aware PI prober (hidden 128) (2,455,767, +0.05 %), `r0_nojepa` = R0 (2,454,531). Training: R0's recipe (36,000 steps, batch 256, AdamW lr 6e-4 → 2e-5
cosine, 1,000 warm-up, weight decay 0.05, clip 1, curriculum 24/48/96 steps at 0/20/50 %, AMP), free-running
rollouts (no teacher forcing), loss = the prober's pseudo-Huber state loss + v2's force term at weight 1 (PI
models: the force vector; direct models: its magnitude, the g-load the planner uses); learning rate ×
{1, 0.3, 3} chosen on seed 0 by a 10 % pilot's validation open-loop state loss (horizon 96). `r0_nojepa` =
R0's architecture exactly (WorldModelV2 + PIProberV2, one member) trained on the prober's state loss only:
no latent prediction loss, no SIGReg; saved as a v2 checkpoint, flown by v4's V2 adapter.

**RL (`jepa/rl_v5.py`, `jepa/rlenv_fast_v5.py`).**

* Fast env: `physics_v4`'s truth model, float64, 5 s decisions in (2, 2, 1) s Euler substeps, bit-exact
  to `rollout_physics` (smoke check 6). Targets on each episode's own footprint (15–85 % of the
  90°-to-0° range, ±60 km crossrange). The observation's heat-rate margin uses Sutton–Graves only (the
  full simulator's interval bound adds a radiative term; the fast-vs-full gap is reported per arm).
* Reward: MPPI's full realized cost including its |γ| > 82° penalty (2 km/deg²) and −3,000 km for a
  flight that never reaches the drogue altitude, × 0.01, plus Φ = −range-to-go shaping.
* Policy: MLP 256-256 with a running observation normalizer. DAgger and BC → PPO: "magsign" head
  (|σ| Gaussian on a sigmoid mean, sign Bernoulli), so reversals are not averaged away; PPO: Gaussian on
  σ (tanh MLP); SAC: tanh-squashed Gaussian (ReLU MLP).
* DAgger expert = FNPEG on the fast env's flights at 5 s decisions, filters at β^5.
* `ppo_matched`: PPO in the full simulator, complete episodes per update (50), the world models' 6,000
  episodes.
* Fast-env evaluation (`fast_eval.json`): the training distribution (2,048 episodes) and the A256 cases
  in the fast env, for R-gap.

**Analysis (`jepa/analyze_v5.py`).** v4's hierarchical bootstrap (seeds, then cases) and cache; R0's
flights are v4's (E0b-corrected); R-each pairs RL seeds 0–2 with R0 seeds 0–2; Holm as §8; `met`
requires p_Holm < 0.05 and the seed rule; J-mech uses 300 s to go (the longest available if 300 s is
missing). ORB-1 Q1/Q2 per arm with Clopper–Pearson (90 %), Q3 by a stratified paired bootstrap over
cases. Landing errors are drawn in the nominal entry frame (downrange/crossrange of the nominal EI).

**Throughput (smoke).** Measured on this laptop's CPU (12 logical cores) in two smoke studies: the
first (check 10's study) while four smoke processes ran in parallel, the second (check 14's study,
same jobs) with two. Both are pessimistic, and smoke models are tiny, so the E3/v2 rows say nothing
about full-size GPU cost. Colab's real rates are measured by the status cell on the first 20 flights
of each arm type.

| Arm type | Smoke flights (1st / 2nd study) | Wall s per flight, median (1st / 2nd) | ms per decision in the flight, median (1st / 2nd) |
|---|---|---|---|
| classical | 5 / 5 | 166 / 153 | 325 / 298 |
| FNPEG (batches of 4) | 23 / 17 | 300 / 255 | 163 / 115 |
| PredGuid (batches of 4) | 19 / 17 | 96 / 85 | 1.0 / 0.56 |
| MPPI + filtered physics (32 samples, smoke) | 17 / 17 | 350 / 390 | 371 / 315 |
| MPPI + v2 (tiny, smoke) | 33 / 33 | 270 / 185 | 162 / 146 |
| MPPI + E3 models (tiny, smoke) | 48 / 48 | 257 / 91 | 101 / 52 |
| RL policy (new arms) | 45 / 45 | 96 / 55 | 1.6 / 0.69 |
| legacy RL (SB3, 1 Hz) | 10 / 10 | 40 / 41 | 1.6 / 1.2 |

Single flight, single thread (`01_checks/compute_per_decision.csv`, smoke): FNPEG 630 ms median per
decision (p99 1,132), PredGuid 1.5 ms, classical 393 ms, RL policies 0.2–0.3 ms.

FNPEG's cost is Python/numpy overhead in the RK4 predictor (thousands of small array operations per
corrector round), so it falls per flight when flights are batched (`guidance_batch = 16` in the study);
the fast env runs about 26,000 environment-steps/s on the CPU with 1,000 environments. DAgger and BC
with the FNPEG expert dominate the RL smoke time: the expert labels every visited state.

**Smoke test output** (`python jepa/smoke_v5.py`, 2026-11-04). The checks ran as several
`--only` invocations in parallel on the laptop (checks 1–4, 7, 8, 11 in one; 5, 6 in one; the notebook
study with 9, 10, 13, 15 in one; 14 in one), and checks 10, 11 and 15 again after their fixes. Each
line below is the latest result of its check; the 12 checks (bit-exact resume) are inside 5, 6 and 14.

```
1. suites
  [PASS] frozen suites verify against SHA256SUMS (10 files)
  [PASS] ORB-1 regenerated (first 12 cases, footprints flown again) == the frozen ORB1.csv bytes (35s)
  [PASS] DEV48 regenerated == the frozen DEV48.csv bytes (102s)
  [PASS] a changed benchmark file fails the checksum
  [PASS] ORB-1 knobs reproduce the dispersed drag and lift (100 draws, simulator aero) (max relative error 4.4e-16)

2. FNPEG
  [PASS] predictor = truth: lands < 1 km (misses [0.505, 0.759] km, 824s for 4 flights)
  [PASS] batched == unbatched (results and telemetry bit-identical)
  [PASS] load relief at amax 3 g keeps the peak under 3.15 g (peak 3.13 g, miss 0.43 km)
  [PASS] filters converge to the true drag scale (within 3%, 150 s after 0.05 g; gate V3's rule) (K_D 1.1200 vs true 1.1200)

3. PredGuid
  [PASS] phase sequences logged and starting with the initial roll (case 30: initial_roll-ballistic-final; case 51: initial_roll-ballistic-final; misses [1.95, 0.91] km, 107s)
  [PASS] final-phase reference and gain design reproducible (1575.5 km ref)

4. mppi_filtered frozen == mppi_nominal; v5 runner == run_lockstep_v4 (E0)
  [PASS] mppi_filtered (filters frozen at 1) == mppi_nominal (2/2 flights identical)
  [PASS] v5 runner == run_lockstep_v4 (physics_nominal) (2/2 flights identical)
  [PASS] v5 runner == run_lockstep_v4 (v2) (2/2 flights identical)
  [PASS] E0 gate passes (997s)

5. E3 models: full size, a few steps, save, resume (12: bit-exact; 13: truncated train state)
  [PASS] every E3 model at full size within ±10% of R0's parameters (R0 2,454,531; mlp_direct 2,454,781 (+0.01%), gru_direct 2,454,690 (+0.01%), gru_pi 2,455,767 (+0.05%), r0_nojepa 2,454,531 (+0.00%))
  [PASS] 12. mlp_direct: killed after the step-3 checkpoint, resumed == uninterrupted (weights True, history True; 25s per run)
  [PASS] 12. gru_direct: killed after the step-3 checkpoint, resumed == uninterrupted (weights True, history True; 22s per run)
  [PASS] 12. gru_pi: killed after the step-3 checkpoint, resumed == uninterrupted (weights True, history True; 27s per run)
  [PASS] 12. r0_nojepa: killed after the step-3 checkpoint, resumed == uninterrupted (weights True, history True; 12s per run)
  [PASS] 13. truncated newest train state: resumes from .prev (logged) and ends identical (weights True, history True)

6. fast env and RL methods (12: checkpoint/resume bit-exact)
  [PASS] fast env == physics_v4.rollout_physics for the same commands (max |diff| 0.0e+00)
  [PASS] a random policy runs 10,000 batched steps (0.4s, 26033 steps/s on the CPU, 0 episodes ended)
  [PASS] 12. ppo_hist: runs, checkpoints, killed and resumed == uninterrupted (weights True, optimizer True, curve True; 10s per run)
  [PASS] 12. sac_hist: runs, checkpoints, killed and resumed == uninterrupted (weights True, optimizer True, curve True; 26s per run)
  [PASS] 12. SAC replay buffer saved as append-only shards (4 shards)
  [PASS] 12. dagger_fnpeg: runs, checkpoints, killed and resumed == uninterrupted (weights True, optimizer True, curve True; 170s per run)
  [PASS] 12. bcppo_fnpeg: runs, checkpoints, killed and resumed == uninterrupted (weights True, optimizer True, curve True; 138s per run)
  [PASS] ppo_matched: runs, checkpoints, killed and resumed == uninterrupted (weights True, optimizer True, curve True; 72s per run)

7. legacy RL
  [PASS] the dummy SB3 zip is found by the inventory (and the checkpoint listed) (3 rows: {'final model': 1, 'checkpoints': 1, 'not found': 1})
  [PASS] the legacy policy is flown with its own settings (1 Hz, rl_env's observation) (['HELD256', 'HELD256'], misses [357.3, 592.4] km)

8. job queue
  [PASS] finished flight skipped, deleted flight re-flown
  [PASS] ORB-1 stopping rule from files: 1 miss > 10 km keeps the arm, 2 stop it (eligible; stopped: 2 flights of earlier blocks missed 10 km (failed Q1))
  [PASS] ORB-1 block 2 waits until block 1 is complete (waits for 3 flights of earlier blocks)

9. analysis and report on the smoke data
  [PASS] every table T01-T12 (missing [])
  [PASS] every figure F01-F10 as PNG and PDF (missing [])
  [PASS] every table and figure has a caption and a data source
  [PASS] criteria table complete (G1-G7, J-each/equiv/mech, R-each/imit/gap, Q1-Q3) (32 rows; missing []; not run: ['G7a', 'G7b', 'G7c'])
  [PASS] bootstrap cache reused on a forced rerun (no recomputation) (120 reused)
  [PASS] an unforced rerun with unchanged inputs skips the analysis
  [PASS] report written by the notebook (1448 kB)
  [PASS] report builds with matplotlib's DejaVu (no system fonts)

10. v5 refuses to write under the v3, v4, data and legacy RL folders
  [PASS] write into the data folder refused
  [PASS] write into the v3 folder refused
  [PASS] write into the v4 folder refused
  [PASS] write into the legacy folder refused
  [PASS] write to outside every folder refused
  [PASS] write to a path escaping the v5 folder refused
  [PASS] write to a Windows extended-length path inside v3 refused
  [PASS] a flight id escaping the results folder refused
  [PASS] an extended-length path inside the v5 folder is accepted (Windows guard fix)

10. the notebook in smoke mode (first run = the whole smoke study)
  [PASS] notebook run 1 ends without error (109.3 min)
  [PASS] notebook run 2 ends without error (36s)
  [PASS] run 2 finds everything done in under a minute (36s)
  [PASS] run 2: no flight flown, no model trained (9 queues checked)
  [PASS] v3, v4 and legacy folders untouched by the whole smoke study

10. the notebook in smoke mode, rerun on the same study after the check 15 fixes (--keep)
  [PASS] notebook rerun on the kept smoke study ends without error (56s)
  [PASS] notebook run 2 ends without error (57s)
  [PASS] run 2 finds everything done in under a minute (57s)
  [PASS] run 2: no flight flown, no model trained (9 queues checked)
  [PASS] v3, v4 and legacy folders untouched by the whole smoke study

11. v3 and v4 code, plans and notebooks unchanged
  [PASS] v3/v4 modules, plans, simulator and the v4 notebook identical to the v4 commit (changed: [])
  [PASS] the v4 notebook builder still writes the committed notebook byte for byte
  [PASS] the v3 builder reproduces the cells of the 5 committed v3 notebooks (differ: []; saved from Colab by the user, checked by git diff only: ['ReEntryAI_JEPA_v3_3_analyze.ipynb'])

13. interrupted writes: Setup repairs each case and the rerun ends with the same outputs
  [PASS] repair: stale .tmp deleted, corrupt flight JSON deleted, .done with a missing artifact removed ({'tmp_deleted': 1, 'json_deleted': 1, 'done_removed': 1, 'files_seen': 777, 'sec': 0.2})
  [PASS] the rerun re-flies the corrupt flight, recomputes E0, and every flight equals the untouched run's (200/200 identical, 620s; differing [])

14. simulated disconnects in the notebook's long stages, then an unchanged rerun
  [FAIL] each forced failure stopped the run (simulated disconnect) ([('cpu_flights:3', True, 186), ('rl_train:1', False, 1139) <- hung launch, killed by hand, ('mppi_flights:2', True, 4652), ('e3_train:4', True, 1221), ('mppi_flights:1', True, 383), ('rl_train:1', True, 1674), ('cpu_flights:2', True, 434)])
  [PASS] the unchanged rerun completes (16.9 min)
  [PASS] finished flights keep their modification times (187 flights before the rerun)
  [PASS] every flight equals the uninterrupted run's (bit-exact on CPU) (200/200 identical; differing []; only in one run [])
  [PASS] every trained model and policy equals the uninterrupted run's (bit-exact on CPU) (19/19 identical; differing [])

15. outputs in the smoke study folder
  [PASS] status.json, status_history.csv, README, frozen config, captions, report, inventory, timing (missing [])
  [PASS] status history appended by every run (8 rows)
  [PASS] a log per stage in 00_study/logs/ (missing [])
  [PASS] curves of every trained model and policy (missing [])
  [PASS] tables and figures (PNG and PDF) on disk (14 tables, (10, 10) figures)

75 PASS, 1 FAIL
```

The one FAIL is the second of check 14's seven launches (`rl_train:1`, right after the `cpu_flights:3`
launch): it never started a kernel and sat for 19 min until it was killed by hand, so its forced failure
never fired. The other six forced failures stopped the notebook where intended, and the unchanged rerun
after them completed and matched the uninterrupted study bit for bit (200/200 flights, 19/19 models),
so resume is unaffected. The cause was not found: later kernel kills, including one with a live CPU pool
(`cpu_flights:2`), did not hang. The harness now gives the kernel a 180 s start timeout and kills a
launch that writes nothing for 25 min (recorded in `notes.stalled_runs`). A rerun
of the same two launches on a fresh study with the final harness (2026-11-04) did not hang: both forced
failures fired (49 s and 1,280 s), with no stall.

---

## 15. Deviations

* **2026-11-03, FNPEG and PredGuid: wings level when |γ| > 30°.** The simulator ends a flight when
  |cos γ| < 0.05 (γ steeper than about 87°). Late in a flight both published laws can still command
  large banks while the path steepens, and the flight then ends at that gate short of the drogue
  altitude. The classical controller already has this override (`gamma_safe`), so both laws get the same
  rule: σ = 0 once |γ| > 30° (`gamma_safe_deg = 30`, in flight and inside both predictors and
  PredGuid's reference trajectory). It acts only in the steep terminal descent; `None` turns it off.
* **2026-11-03, PredGuid: `VFINAL` = 29,000 fps (Bairstow: 25,000 fps).** At 0.05 g Apollo's initial-roll
  logic (B.18) sends an entry still faster than `VFINAL` into Huntest, the lunar-return lofting logic, and
  a slower one straight to the final phase. Our entries start at 24,400-27,600 fps (A256) and at
  25,011 fps (ORB-1's 7,623.5 m/s), so with Bairstow's value the branch would flip from case to case
  across the middle of both suites. With 29,000 fps every LEO entry flies initial roll → final, the
  branch Apollo's guidance takes for entries below lunar-return speed (Putnam).
* **2026-11-04, smoke check 11 (not a study change).** The committed `ReEntryAI_JEPA_v3_3_analyze.ipynb`
  is the copy the user saved from Colab (commit c747b92, "Created using Colab", an extra `%pip install`
  cell), so the v3 builder cannot reproduce it cell for cell. The check compares it by `git diff`
  against the v4 commit (unchanged) and rebuilds the other four v3 notebooks exactly.
