# SpaceJEPA v5: strong baselines and a standard benchmark (plan and implementation brief)

Branch: `space-jepa5` (from `space-jepa4`). Status: **built and smoke-tested (§14); revised
2026-11-11 to 14 after a pre-run review of what a reviewer would object to (§15 lists every change, §16 the
limitations that remain); the user runs it on Colab.** Written 2026-11-03, before any v4 result was analyzed, so nothing here was chosen after
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
| D3 | Standard scenario ORB-1 (§3.2) | Lu (2014) CEV orbital entry state, adapted (γ0 = −1.40°, targets at 25/50/75% of footprint), Bairstow (2006) dispersions, Orion EM-1 accuracy requirement. **Target 0 kept by the user (2026-11-04)** although only 27/1,002 cases reach it at constant bank under 4.5 g (§14) | Lu's exact 3,769 km range (infeasible in our simulator at a 90° bank limit, see §1.3) |
| D4 | EI speed dispersion for ORB-1 | 15 m/s (3σ), our choice: no published LEO value was found | another value the user has a source for |
| D5 | Density-shear suite (needs a simulator hook) | out | in, as Tier 3, with an explicit `configs/`-free hook. **Decided by the user (2026-11-15): in, as DS1 in Tier 1 (§3.4), with a default-off simulator hook (§15)** |

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
| A frozen, checksummed benchmark (suites A256, ORB-1, DEV48; DS1 since 2026-11-15) with a scoring module | winds, Mach-dependent aerodynamic errors (density-shear atmospheres moved in scope on 2026-11-15 as DS1, §3.4) |
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

### 3.4 DS1 (density errors that change during flight; added 2026-11-15, D5)

A256's and ORB-1's errors are constant multipliers for the whole flight, which is exactly what FNPEG's
and PredGuid's filters, `mppi_filtered`'s scales and v2's training randomization assume. DS1 tests the
case where that assumption is wrong: the density error changes with altitude, so a ratio learned high
up is wrong lower down.

**Conditions.** A256's 64 entry cases (`configs/v3_testset.csv`) × 4 atmospheres = 256 conditions, all
with lift and mass scales 1 except where stated:

| Atmosphere | Density factor f(h) on the NRLMSIS table |
|---|---|
| `shear_dense_high` | 1.15 above 70 km, 0.85 below 40 km, linear between |
| `shear_dense_low` | 0.85 above 70 km, 1.15 below 40 km, linear between |
| `gram_like` | exp(σ(h) z(h)): z a unit Gaussian process in altitude (squared-exponential, 10 km correlation, truncated at ±3), σ from 4% at 0 km to 12% at 100 km and above; its own draw per case |
| `gram_like_lift` | another draw of the same process, and L/D × 0.92 |

The GRAM-like profiles imitate the shape of Earth-GRAM's density perturbations (a few percent near the
ground, growing with altitude, correlated over kilometres); they are not drawn from GRAM. f(h) is stored
on a 2 km grid from 0 to 160 km and interpolated linearly (`configs/benchmark_v5/DS1.json`, seed
20261115, generated by `benchmark_v5.build_ds1`; checksummed with the other suites; the A256, ORB-1 and
DEV48 files are unchanged).

**Simulator hook.** `math_3d.nominal_aero_forces_from_state` multiplies the table density by f(h) when
the parameters carry a `density_profile`, which `data.apply_domain` sets from the domain. DS1's domains
name the profile by a number (`density_profile_id` = 1000 × (atmosphere index + 1) + case, resolved from
`DS1.json`), because every flight record copies the domain's values as floats (`fly_v3.result_record`).
It is off when absent, so every v3, v4 and earlier v5 flight is bit-identical (smoke check 1 flies a
whole flight with a profile of ones and without one and compares them bit for bit; check 11 verifies
that `math_3d.py` and `data.py` differ from the v4 commit by added hook lines only). Because that one function feeds the
equations of motion, the sensed accelerations and the g-load, every controller sees the shear only
through its effect on the flight; the heating model does not see it (heat rate and heat load on DS1
use the table density).

**Arms.** `classical`, `fnpeg`, `predguid`, `mppi_filtered`, `v2single` and the Tier 1 RL arms, one
flight per condition (2,048 flights). A seeded arm flies seed case mod 5, as on ORB-1. No feasibility
flag is computed (no footprint under an altitude-varying density).

### 3.4b AS1 (the published Orion aerodynamics; added 2026-11-15, user decision)

The simulator's trim schedule gives the capsule more lift than the published Orion database below about
6 km/s (`jepa/simulator_fidelity_audit.md` §4.1). AS1 asks whether RQ1's answer depends on that.

**Conditions.** A256's 64 entry cases and targets, nominal air, × 2 true vehicles = 128 conditions:

| Vehicle | True trim schedule |
|---|---|
| `bibb_ld` | L/D from Bibb Part I Fig. 27 at FMV ≥ 8 (0.280 at 8, 0.263 at 15, 0.250 at 25, 0.245 at 40), capped at 0.28 below FMV 8; the simulator's CD |
| `bibb_ld_cd` | the same L/D, and CD from Bibb Part II Fig. 40(b) and Part I Fig. 26 (1.37 at FMV 2, 1.34 at 3, 1.38 at 8 and above) |

Values are read off the papers' figures (about ±0.03) and frozen in `configs/benchmark_v5/AS1.json`
(`benchmark_v5.build_as1`; checksummed; the other suite files are unchanged). No controller is told:
FNPEG's and PredGuid's predictors, `mppi_filtered`'s physics model, the classical predictor and the
learned models keep what they had, so the truth differs from every model by an aerodynamic error that
changes with Mach number (L/D × 0.70 at FMV 8 to × 0.93 at FMV 25; inside the 0.70–1.30 lift range the
learned models were trained on, but never constant in training).

**Reachability.** The targets were placed for the simulator's vehicle. Each condition carries
`feasible`: its target is inside the true vehicle's 4.5 g constant-bank footprint, by the code and rule
that made A256's own flags (`testset_v3.footprint`, margin > 0). 60 of 64 for either vehicle. Results
are reported on all conditions and on the reachable ones.

**Hook.** `data.apply_domain` sets the simulator's existing `orion_trim_schedule` parameter when a
domain carries `aero_schedule_id` (1 or 2) and restores the parameters otherwise. `math_3d.py` is not
touched; `data.py`'s addition sits inside the density hook's block. Off when absent (smoke check 1: a
later domain without the id gives parameters equal to fresh ones).

**Arms.** DS1's eight (`classical`, `fnpeg`, `predguid`, `mppi_filtered`, `v2single`, the Tier 1 RL
arms), one flight per condition (1,024 flights), DS1's seed rule.

### 3.5 Scoring (`benchmark_v5.score`)

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
| V1 FNPEG verification | on the 12 DEV48 cases, nominal air, predictor = true parameters, no load limit (changed 2026-11-11, §15: was 16 A256 cases; a gate must not look at test cases) | miss < 2 km in ≥ 11 of 12 (was 15 of 16; Lu reports < 1 km in his own simulator; ours adds roll dynamics the predictor does not model), Newton iteration converges in ≤ 6 iterations per cycle in ≥ 99% of cycles |
| V2 load relief | same 12 cases, amax = 3.0 g, cases whose unconstrained peak exceeds 3.5 g | realized peak ≤ amax + 0.15 g and miss < 2 km where feasible; as amax drops below feasibility the miss grows gradually (Lu Fig. 17), checked on 2 cases at 2.7/2.5/2.3 g **Changed 2026-11-14 (§15):** a flight also passes when the limit was infeasible for it: relief on before the load crossed amax, full lift-up commanded for the 20 s before the peak, and a lower peak than the unrelieved flight |
| V3 filters | thin (×0.88) and dense (×1.12) air, the first 8 DEV48 cases each | the filter's drag ratio K_D is within 3% of the true scale by 150 s after the 0.05 g point |
| V4 PredGuid verification | the 12 DEV48 cases, nominal air, nominal model | miss < 3.5 km in ≥ 10 of 12 (amended 2026-11-17 from 11 of 12, §15, before any test flight; was 15 of 16; Bairstow's Monte Carlo bound); the phase sequence is logged and sensible (initial roll → final for these ranges) |
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
| `fnpeg_norelief` | exploratory (added 2026-11-15): the tuned `fnpeg` with Lu's load relief switched off, nothing else changed and no separate tuning. Shows whether the relief costs FNPEG accuracy (`jepa/guidance_audit.md` §3) | A256 | 256 |

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
  predictor multiplies modeled D and L by K_D and K_L along the whole prediction. β tuned; the tuned β is
  also `mppi_filtered`'s filter and the DAgger/BC expert's (§15, 2026-11-17).
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
| `mppi_filtered` | v4's `mppi_nominal` with the physics model's density and lift scales set every decision from FNPEG's filters (density_scale = K_D, lift_scale = K_L / K_D, mass 1), at FNPEG's tuned β (its flights wait for FNPEG's tuning; §15, 2026-11-17). Same planner and settings as v2. The direct competitor to "implicit system identification from history." | A256 (5 noise indices), ORB-1 | 1,280 + ≤ 3,000 |
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
`WIDE_DOMAIN`). Hyperparameters tuned on DEV48 only (≤ 12 configurations per method). Seeds 0–4
(since 2026-11-15, §15; was 0–2).

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

**Changed 2026-11-11 (§15):** each guidance law gets a second round of 12 configurations drawn around its
round-1 best (36 in total, same objective, the choice is made over both rounds), and every RL method
uses the full 12 configurations. More tuning for the baselines only; nothing else is tuned.

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
`v2single` → E3 training and flights → E4 training and flights → DS1 (§3.4) → AS1 (§3.4b) → ORB-1 blocks 2–3 (arms
still eligible) → exploratory arms → open loop → analysis → report. Within an arm, seed-major.

### 4.9 Compute tiers (D1)

Rates must be measured, not assumed: the status cell measures each arm type's flights per hour on the
first 20 flights and prints projected hours and Colab units (at the runtime's rate). Rough sizes:

| Tier | Contents | Flights (new) | Training |
|---|---|---|---|
| 1 | gates, legacy RL eval, tuning (2 × 24 × 48 guidance flights + RL tuning), the A256 arms of §4.1–4.3 except the exploratory ones, `dagger_fnpeg`, `bcppo_fnpeg`, `ppo_hist` (5 seeds each since 2026-11-15; were 3), ORB-1 block 1 for `classical`, `fnpeg`, `predguid`, `mppi_filtered`, `v2single`; since 2026-11-11 (§15) also `fnpeg_180`, `predguid_180`, `fnpeg_oracle` on A256 and the larger tuning budget (2 × 36 × 48 guidance flights, 12 RL configurations per method); since 2026-11-15 (§15) DS1 (§3.4: 8 arms × 256 conditions), AS1 (§3.4b: 8 arms × 128 conditions) and `fnpeg_norelief` (256) | ≈ 27,900 (+1,280 for AS1 and `fnpeg_norelief`; ≈ 9,700 of them cheap policy flights; was ≈ 23,000 before the 2026-11-15 additions: +1,536 E4 policy flights for seeds 3–4, +2,048 DS1 flights) | 20 dynamics models + 15 RL runs (was 9) (+ 36 RL tuning runs at 10% budget) |
| 2 | ORB-1 blocks 2–3 (eligible arms only), ORB-1 block 1 for the non-JEPA and RL arm with the highest A256 success@5km (descriptive), `sac_hist`, `ppo_matched`, `mppi_truth` on ORB-1 | ≈ 4,000–8,000 | 6 RL runs |
| 3 | density-shear suite (D5): moved to Tier 1 as DS1 on 2026-11-15 (§3.4) | — | none |

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
* ORB-1: each case is one unit; learned arms have one seed per case, so the bootstrap draws the seed
  blocks with replacement and then cases within each drawn block (two-stage; since 2026-11-17, §15; before,
  cases within fixed blocks); Q1 uses exact Clopper–Pearson; Q3 uses block 1 paired by case.
* RL arms: seeds 0–4; "holds in a seed" needs 4 of 5, as for the learned models (since 2026-11-15, §15;
  was seeds 0–2 and 3 of 3).
* DS1 (§3.4, added 2026-11-15): each case is one unit (its four atmospheres averaged); learned arms have
  one seed per case (case mod 5), so the seeds and then the cases within each seed are resampled, as on
  ORB-1 (AS1 likewise).
* Seed-level t-test (added 2026-11-17, §15, before any comparison had data): every paired comparison also
  reports the t-test on its per-seed differences (the 5 training seeds as the units, 4 df), next to its
  bootstrap. Verdicts follow the rules below; a verdict the t-interval does not support is marked in T05
  (`boot_t_disagree`) and in the report.
* Every estimate is reported with its interval, including failures. Secondary metrics are exploratory.
* Added 2026-11-11 (§15), no criterion changed:
  * **Verdicts need complete data.** The analysis runs on every "Run all". A pass/fail criterion is decided
    only when every paired flight of both arms exists; before that its status is "pending" and its numbers
    are marked provisional. Q1: an arm passes when its bound is reached, fails at 2 misses (the stopping
    rule), and is "pending" otherwise (too few flights for the bound, as in Tier 1).
  * **Clustered sensitivity.** A256's 256 conditions come from 64 entry cases. Next to every pooled
    interval the report gives the interval with each entry case resampled together with its four
    atmospheres, and flags a "met" verdict that this interval does not support.

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
* **DS1** (added 2026-11-15) Does RQ1's answer hold when the density error changes with altitude instead
  of staying constant for the whole flight?
* **AS1** (added 2026-11-15) Does RQ1's answer hold when the true vehicle has the published Orion
  aerodynamics, with less lift than the simulator's own vehicle?

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
| G7d, G7e | what the load relief does to FNPEG (exploratory, added 2026-11-15) | `fnpeg_norelief − fnpeg` and `v2single − fnpeg_norelief`, success@5km, with both arms' peak g: estimates only |

### 8.2 RQ2 (A256; paired log miss vs R0, v4 flights; Holm over the 4 arms)

| # | Question | Passes if |
|---|---|---|
| J-each | R0 beats this model | mean paired log-miss difference (arm − R0) > 0, 95% CI excludes 0, Holm p < 0.05, ≥ 4 of 5 seeds |
| J-equiv | the JEPA objective does not matter (`r0_nojepa` only) | 90% CI of the mean paired log-miss difference inside ±0.15 (TOST). If neither J-each nor J-equiv passes for `r0_nojepa`, the answer is "inconclusive" |
| J-mech | mechanism | Spearman ρ across R0, v1 and the 4 models between open-loop error at 300 s and closed-loop mean log miss, with a bootstrap CI over seeds (estimate) |

### 8.3 RQ3 (A256; Holm over the Tier 1 RL arms)

| # | Question | Passes if |
|---|---|---|
| R-each | v2 beats this RL arm | success@5km `v2single − arm` > 0 pooled, 95% CI excludes 0, Holm p < 0.05, holds in ≥ 4 of 5 RL seeds (since 2026-11-15; was 3 of 3 seeds) |
| R-imit | imitation gap | `dagger_fnpeg − fnpeg` success@5km: estimate |
| R-gap | sim-to-sim gap | fast-vs-full gap per arm: estimate |

### 8.4 RQ4 (ORB-1)

| # | Question | Passes if |
|---|---|---|
| Q1 | accuracy requirement, per arm | requirement LCB ≥ 0.99865 over the flights the arm flew (descriptive pass/fail, no p-value) |
| Q2 | crew g limits, per arm | fraction of g-duration-compliant flights with a two-sided 95% Clopper–Pearson interval (descriptive) |
| Q3 | v2 vs the published laws on the standard scenario (Holm over 2) | success@5km `v2single − fnpeg` and `v2single − predguid` on block 1, > 0 with 95% CI excluding 0 and Holm p < 0.05 |

Added 2026-11-11 (§15), estimates only (table T08b): Q3 again (a) on the load-feasible cases (the frozen
`feasible_4p5g` flag: all of targets 1 and 2, 27 of target 0's 1,002), (b) with success counted only when
the flight also kept peak g ≤ 4.6 g, and (c) per target. Target 0 is mostly unreachable within 4.5 g, so
there a law that enforces the limit gives up range and a soft limit can be exceeded; these rows show
whether the pooled Q3 is driven by that. Every table also reports success@5km within the load limit.

### 8.5 RQ5

Descriptive table (§5 compute per decision, parameter counts).

### 8.6 DS1 (added 2026-11-15, before any test flight; own Holm family {D1, D2})

| # | Question | Passes if |
|---|---|---|
| D1 | v2 beats FNPEG under density shear | success@5km `v2single − fnpeg` on DS1 > 0 (cases as units, atmospheres averaged), 95% CI excludes 0, Holm p < 0.05 |
| D2 | v2 beats PredGuid under density shear | same with `predguid` |
| D3, D4, D-*rl arm* | v2 vs `mppi_filtered`, `classical` and each Tier 1 RL arm on DS1 | estimates only |
| D-shift | what each arm loses from A256 to DS1 on the same 64 cases | descriptive (T13b) |

A seed rule does not apply (each case is flown by one seed). DS1 is a separate family so that adding it
cannot change any RQ1–RQ4 verdict.

### 8.7 AS1 (added 2026-11-15, before any test flight; own Holm family {A1, A2})

| # | Question | Passes if |
|---|---|---|
| A1 | v2 beats FNPEG on the published Orion aerodynamics | success@5km `v2single − fnpeg` on AS1 > 0 (cases as units, the two vehicles averaged), 95% CI excludes 0, Holm p < 0.05 |
| A2 | v2 beats PredGuid on the published Orion aerodynamics | same with `predguid` |
| A1-reach, A2-reach | the same on the reachable conditions only | estimates |
| A3, A4, A-*rl arm* | v2 vs `mppi_filtered`, `classical` and each Tier 1 RL arm on AS1 | estimates only |
| A-shift | what each arm loses from A256's nominal atmosphere to AS1 on the same 64 cases | descriptive (T14b) |

Statistics as DS1 (§6). AS1 is a separate family, so it cannot change any RQ1–RQ4 or DS1 verdict.

---

## 9. The notebook: `ReEntryAI_JEPA_v5.ipynb` (one notebook, generated)

Built by `jepa/build_notebook_v5.py` (reuse `AUTO_DISCONNECT_CODE`, `DONE_MD`, `DONE` from
`build_notebooks_v3.py` and v4's Settings-cell pattern: `MODE`, `RUN` per stage, `AUTO_DISCONNECT`,
`CPU_WORKERS`, `GPU_WORKERS`). Colab metadata: GPU, L4 (a T4 for part of 2026-11-17, before any GPU work;
§15). **Resume rule as v4: after any disconnect, Run
all continues where it stopped**; nothing finished is redone. CPU-only runtimes skip GPU stages with a
message. Since 2026-11-16 (§15) a GPU runtime skips the CPU-pool test flights (`GPU_PASS_SKIPS_CPU = True` in
the Settings cell), and since 2026-11-17 every CPU stage: the gates and S0 (a GPU pass refuses to start if they
are not done), FNPEG/PredGuid tuning, the analysis and the report; section 15b releases a GPU pass's runtime once
its GPU work is done. RL tuning's DEV48 scoring flights are the exception and fly in the GPU pass. The run is done
in passes the user switches by hand, CPU → L4 → CPU, each releasing the runtime at its end and printing what the
next pass must do.

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
| T05 | criteria table (every G, J, R, Q) with value, CI, Holm p, met, status (met / not met / pending / estimate / not run), paired flights present and expected, clustered interval |
| T06 | RQ2: paired log miss vs R0, open-loop errors, parameter counts |
| T07 | RQ3: RL arms, success, miss, fast-vs-full gap, training budget used |
| T08 | RQ4: ORB-1 per arm: blocks flown, success@2/5/10km, requirement LCB, g compliance, per target range |
| T08b | RQ4: Q3 and its sensitivity analyses (load-feasible cases, within the load limit, per target) |
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
benchmark datasheet, methods (incl. the exact FNPEG and PredGuid settings), limitations (§16),
deviations (§15), appendix with every table. While any criterion is pending the summary says
"provisional report".

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
| `jepa/benchmark_v5.py` | suite generation (A256 pointer, ORB-1, DEV48, DS1), checksums, loading, `score` (§3.5, §5) |
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
  constant-bank flag is stored per case and not used as a filter. **User decision 2026-11-04: keep
  target 0 as frozen** (the benchmark files and checksums are unchanged).
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

**PredGuid's tuning is a 4-D search on DEV48** (added 2026-11-17, audit C3, §15). Every DEV48 flight goes initial
roll → ballistic → final, so four of the eight ranges never change a command there: with Q7F, Q19, KA and K_ρ's
filter gain all at the low end, or all at the high end, each of the 48 DEV48 conditions flies bit-identically to
the defaults (commanded bank, trajectory and miss). Only the logged density estimate differs: K_ρ's gain sets it,
and only the up-control predictor reads it. Q7F and Q19 act in HUNTEST and up-control only. The initial roll's
D > KA test comes out the same for every KA in the range. LATBIAS at its high end, the positive control, changes
39 of the 48 (miss −2.4 to +4.4 km). So the tuning searches LATBIAS, KLAT, overcontrol and GMAX. The ranges and the
24 configurations stay as frozen. The chosen values of the other four are random draws that DEV48 could not
rank, and they would matter only on a test entry that enters HUNTEST or up-control. Each flight records its
phase sequence (`predguid_phases`).

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
* DAgger expert = FNPEG on the fast env's flights at 5 s decisions, filters at β^5 with β the tuned FNPEG's
  (`fnpeg_cfg["filter_beta"]`; `RLCfg.expert_beta_1hz = 0.9` only when the chosen configuration has none,
  i.e. configuration 0; §15, 2026-11-17).
* `ppo_matched`: PPO in the full simulator, complete episodes per update (50), the world models' 6,000
  episodes.
* Fast-env evaluation (`fast_eval.json`): the training distribution (2,048 episodes) and the A256 cases
  in the fast env, for R-gap.

**Analysis (`jepa/analyze_v5.py`).** v4's hierarchical bootstrap (seeds, then cases) and cache; R0's
flights are v4's (E0b-corrected); R-each pairs RL seed k with `v2single` seed k (seeds 0–4 since
2026-11-15; were 0–2); Holm as §8; `met`
requires p_Holm < 0.05 and the seed rule; J-mech uses 300 s to go (the longest available if 300 s is
missing). ORB-1 Q1 per arm with the one-sided 90% Clopper–Pearson lower bound, Q2 with a two-sided 95%
Clopper–Pearson interval (the code's; this line said 90% until 2026-11-17, §15), Q3 by a two-stage paired
bootstrap (training seeds, then cases within each seed) with the seed-level t-test beside it. Landing errors are drawn in the nominal entry frame (downrange/crossrange of the nominal EI).

**Throughput (smoke, 2026-11-14).** Measured on this laptop's CPU (12 logical cores) in the smoke study of
check 10 while check 14's study ran beside it, so every number is pessimistic; smoke models are tiny (the
E3 and v2 rows say nothing about full-size GPU cost). Colab's real rates are measured by the status cell
on the first 20 flights of each arm type.

| Arm type | Smoke flights | Wall s per flight (median) | ms per decision in the flight (median) |
|---|---|---|---|
| classical | 5 | 160 | 314 |
| FNPEG, including the 180° and oracle arms (batches of 4) | 37 | 320 | 163 |
| PredGuid, including the 180° arm (batches of 4) | 29 | 124 | 1.2 |
| MPPI + filtered physics (32 samples, smoke) | 17 | 475 | 361 |
| MPPI + v2 (tiny, smoke) | 33 | 173 | 116 |
| MPPI + E3 models (tiny, smoke) | 48 | 101 | 39 |
| RL policy (new arms) | 45 | 41 | 0.55 |
| legacy RL (SB3, 1 Hz) | 10 | 48 | 1.8 |

Single flight, single thread (`01_checks/compute_per_decision.csv`, smoke of 2026-11-04): FNPEG 630 ms
median per decision (p99 1,132), PredGuid 1.5 ms, classical 393 ms, RL policies 0.2–0.3 ms. With the full
study's batch of 16 on 11 workers, the 55 gate flights on DEV48 (43 FNPEG, 12 PredGuid) took 10 minutes.

FNPEG's cost is Python/numpy overhead in the RK4 predictor (thousands of small array operations per
corrector round), so it falls per flight when flights are batched (`guidance_batch = 16` in the study);
the fast env runs about 26,000 environment-steps/s on the CPU with 1,000 environments. DAgger and BC
with the FNPEG expert dominate the RL smoke time: the expert labels every visited state.

**Smoke test output** (`python jepa/smoke_v5.py`, 2026-11-14, the code of this commit). The checks ran as
separate `--only` invocations: 1–8 in one (7 once more after its fixture file, removed by a temp-folder
cleanup, was restored); a fresh smoke study with 10, 9, 15, 13 in one (run 1 took
189.9 min of wall time, which includes about 45 min of the laptop asleep); 14 beside it; then 10, 9, 15
and 11 once more on the same study after two fixes the first pass prompted (the analysis fingerprint and
gate V2's rule, §15). Each line below is the latest result of its check; the 12 checks (bit-exact resume)
are inside 5, 6 and 14.

```
1. suites
  [PASS] frozen suites verify against SHA256SUMS (10 files)
  [PASS] ORB-1 regenerated (first 12 cases, footprints flown again) == the frozen ORB1.csv bytes (15s)
  [PASS] DEV48 regenerated == the frozen DEV48.csv bytes (34s)
  [PASS] a changed benchmark file fails the checksum
  [PASS] ORB-1 knobs reproduce the dispersed drag and lift (100 draws, simulator aero) (max relative error 4.4e-16)

2. FNPEG
  [PASS] predictor = truth: lands < 1 km (misses [0.505, 0.759] km, 290s for 4 flights)
  [PASS] batched == unbatched (results and telemetry bit-identical)
  [PASS] load relief at amax 3 g keeps the peak under 3.15 g (peak 3.13 g, miss 0.43 km)
  [PASS] filters converge to the true drag scale (within 3%, 150 s after 0.05 g; gate V3's rule) (K_D 1.1200 vs true 1.1200)

3. PredGuid
  [PASS] phase sequences logged and starting with the initial roll (case 30: initial_roll-ballistic-final; case 51: initial_roll-ballistic-final; misses [1.95, 0.91] km, 52s)
  [PASS] final-phase reference and gain design reproducible (1575.5 km ref)

4. mppi_filtered frozen == mppi_nominal; v5 runner == run_lockstep_v4 (E0)
  [PASS] mppi_filtered (filters frozen at 1) == mppi_nominal (2/2 flights identical)
  [PASS] v5 runner == run_lockstep_v4 (physics_nominal) (2/2 flights identical)
  [PASS] v5 runner == run_lockstep_v4 (v2) (2/2 flights identical)
  [PASS] E0 gate passes (546s)

5. E3 models: full size, a few steps, save, resume (12: bit-exact; 13: truncated train state)
  [PASS] every E3 model at full size within ±10% of R0's parameters (R0 2,454,531; mlp_direct 2,454,781 (+0.01%), gru_direct 2,454,690 (+0.01%), gru_pi 2,455,767 (+0.05%), r0_nojepa 2,454,531 (+0.00%))
  [PASS] 12. mlp_direct: killed after the step-3 checkpoint, resumed == uninterrupted (weights True, history True; 11s per run)
  [PASS] 12. gru_direct: killed after the step-3 checkpoint, resumed == uninterrupted (weights True, history True; 10s per run)
  [PASS] 12. gru_pi: killed after the step-3 checkpoint, resumed == uninterrupted (weights True, history True; 10s per run)
  [PASS] 12. r0_nojepa: killed after the step-3 checkpoint, resumed == uninterrupted (weights True, history True; 11s per run)
  [PASS] 13. truncated newest train state: resumes from .prev (logged) and ends identical (weights True, history True)

6. fast env and RL methods (12: checkpoint/resume bit-exact)
  [PASS] fast env == physics_v4.rollout_physics for the same commands (max |diff| 0.0e+00)
  [PASS] a random policy runs 10,000 batched steps (0.3s, 37751 steps/s on the CPU, 0 episodes ended)
  [PASS] 12. ppo_hist: runs, checkpoints, killed and resumed == uninterrupted (weights True, optimizer True, curve True; 14s per run)
  [PASS] 12. sac_hist: runs, checkpoints, killed and resumed == uninterrupted (weights True, optimizer True, curve True; 15s per run)
  [PASS] 12. SAC replay buffer saved as append-only shards (4 shards)
  [PASS] 12. dagger_fnpeg: runs, checkpoints, killed and resumed == uninterrupted (weights True, optimizer True, curve True; 151s per run)
  [PASS] 12. bcppo_fnpeg: runs, checkpoints, killed and resumed == uninterrupted (weights True, optimizer True, curve True; 129s per run)
  [PASS] ppo_matched: runs, checkpoints, killed and resumed == uninterrupted (weights True, optimizer True, curve True; 72s per run)

7. legacy RL
  [PASS] the dummy SB3 zip is found by the inventory (and the checkpoint listed) (3 rows: {'final model': 1, 'checkpoints': 1, 'not found': 1})
  [PASS] the legacy policy is flown with its own settings (1 Hz, rl_env's observation) (['HELD256', 'HELD256'], misses [357.3, 592.4] km)

8. job queue
  [PASS] finished flight skipped, deleted flight re-flown
  [PASS] ORB-1 stopping rule from files: 1 miss > 10 km keeps the arm, 2 stop it (eligible; stopped: 2 flights of earlier blocks missed 10 km (failed Q1))
  [PASS] ORB-1 block 2 waits until block 1 is complete (waits for 3 flights of earlier blocks)
  [PASS] tuning round 2: drawn once round 1 is scored, around its best, recorded, identical on a rerun (24 then 36 configurations; center = configuration 7)
  [PASS] Tier 1 flies the 180° and oracle arms; gates V1-V4 fly DEV48 only; RL tuning uses 12 configurations (12 V1 flights on ['DEV48'])

9. analysis and report on the smoke data
  [PASS] every table T01-T12 (missing [])
  [PASS] every figure F01-F10 as PNG and PDF (missing [])
  [PASS] every table and figure has a caption and a data source
  [PASS] criteria table complete (G1-G7, J-each/equiv/mech, R-each/imit/gap, Q1-Q3) (34 rows; missing []; not run: [])
  [PASS] every criterion has a status; on the complete smoke study only Q1 can be pending (too few cases) (pending: ['Q1-classical', 'Q1-fnpeg', 'Q1-mppi_filtered', 'Q1-predguid', 'Q1-v2single'])
  [PASS] paired-flight counts and clustered-case intervals in T04; Q3 sensitivity table T08b (18 pooled comparisons, 18 complete)
  [PASS] bootstrap cache reused on a forced rerun (no recomputation) (150 + 30 clustered reused)
  [PASS] an unforced rerun with unchanged inputs skips the analysis
  [PASS] report written by the notebook (1525 kB)
  [PASS] report builds with matplotlib's DejaVu (no system fonts)
  [PASS] a criterion with a missing flight is pending (no verdict, numbers marked provisional); others decided (G1: 14 of 16 paired flights, status pending; G2: not met)

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
  [PASS] notebook run 1 ends without error (189.9 min)
  [PASS] v3, v4 and legacy folders untouched by the whole smoke study
  (run 2 of this first pass took 108 s: it redid the analysis; fixed, see the rerun below)

10. the notebook in smoke mode, rerun on the same study with the final code (--keep)
  [PASS] notebook rerun on the kept smoke study ends without error (60s)
  [PASS] notebook run 2 ends without error (21s)
  [PASS] run 2 finds everything done in under a minute (21s)
  [PASS] run 2: no flight flown, no model trained (9 queues checked)
  [PASS] v3, v4 and legacy folders untouched by the whole smoke study

11. v3 and v4 code, plans and notebooks unchanged
  [PASS] v3/v4 modules, plans, simulator and the v4 notebook identical to the v4 commit (changed: [])
  [PASS] the v4 notebook builder still writes the committed notebook byte for byte
  [PASS] the v3 builder reproduces the cells of the 5 committed v3 notebooks (differ: []; saved from Colab by the user, checked by git diff only: ['ReEntryAI_JEPA_v3_3_analyze.ipynb'])

13. interrupted writes: Setup repairs each case and the rerun ends with the same outputs
  [PASS] repair: stale .tmp deleted, corrupt flight JSON deleted, .done with a missing artifact removed ({'tmp_deleted': 1, 'json_deleted': 1, 'done_removed': 1, 'files_seen': 912, 'sec': 0.3})
  [PASS] the rerun re-flies the corrupt flight, recomputes E0, and every flight equals the untouched run's (230/230 identical, 655s; differing [])

14. simulated disconnects in the notebook's long stages, then an unchanged rerun
  [PASS] each forced failure stopped the run (simulated disconnect) ([('cpu_flights:3', True, 139), ('rl_train:1', True, 3382), ('mppi_flights:2', True, 2063), ('e3_train:4', True, 3343), ('mppi_flights:1', True, 166), ('rl_train:1', True, 928), ('cpu_flights:2', True, 166)])
  [PASS] the unchanged rerun completes (21.9 min)
  [PASS] finished flights keep their modification times (185 flights before the rerun)
  [PASS] every flight equals the uninterrupted run's (bit-exact on CPU) (230/230 identical; differing []; only in one run [])
  [PASS] every trained model and policy equals the uninterrupted run's (bit-exact on CPU) (19/19 identical; differing [])

15. outputs in the smoke study folder
  [PASS] status.json, status_history.csv, README, frozen config, captions, report, inventory, timing (missing [])
  [PASS] status history appended by every run (12 rows)
  [PASS] a log per stage in 00_study/logs/ (missing [])
  [PASS] curves of every trained model and policy (missing [])
  [PASS] gate V1-V4 flights are on the dev suite, none on a test suite (DEV48: ['gate_V1_fnpeg', 'gate_V2_fnpeg_2.7g', 'gate_V2_fnpeg_3g', 'gate_V3_fnpeg', 'gate_V4_predguid']; A256: [])
  [PASS] tuning: round 2 recorded (refine.json) and both rounds scored before the choice ({'fnpeg': (True, 3, [1, 2]), 'predguid': (True, 3, [1, 2])})
  [PASS] the 180° and oracle arms flew in Tier 1 (['fnpeg_180', 'fnpeg_oracle', 'predguid_180'])
  [PASS] tables and figures (PNG and PDF) on disk (15 tables, (10, 10) figures)

81 PASS, 0 FAIL
```

**2026-11-15, after DS1, the density hook and five RL seeds: 85 PASS, 0 FAIL (127 min).** The listing
above is the 2026-11-14 run; the new checks are DS1's (the regenerated file, 256 conditions and
their profiles, the hook's drag ratio and its clearing by the next domain, a flight with a profile of
ones against one without), the hook-only diff of `math_3d.py` and `data.py` against the v4 commit, and the DS1 tables
T13 and T13b with criteria D1, D2 and D-shift.

**2026-11-16, after AS1, `fnpeg_norelief`, the interval fields and the classical timing: 92 PASS,
0 FAIL (490 min on a busy laptop; the smoke study is larger and the machine was shared).** New checks:
AS1's regenerated file and its 128 conditions, the aero-schedule hook and its exact restore, a flight on
the AS1 vehicle landing shorter, T14 and T14b with A1 and A2 in their own family, the interval fields in
every v5 flight record (no shield step, no intervention), and `classical_noscreen` in T09.

History. On 2026-11-04 (the first build) the same checks gave 75 PASS and 1 FAIL: one of check 14's
seven launches never started a kernel and was killed by hand; the cause was not found and it did not
recur (2026-11-04 rerun of two launches, 2026-11-14 all seven). On 2026-11-13 a rerun was lost when the
app that had started it closed: the test's console went with it, the notebook's worker processes could
no longer start and the run hung. The smoke test now says so in `run_nb`, and its stall watchdog ignores
a sleep of the machine. Run it from a console that stays open.

**Revision of 2026-11-11 (pre-run review; the reasons are in §15).**

* Gates: `gate_suite = "DEV48"`, every DEV48 case (12); V1 and V4 need 11 of 12. E0 (runner equivalence,
  v2 and physics adapters, no baseline) still uses three A256 conditions.
* Tier 1 now flies `fnpeg_180`, `predguid_180` and `fnpeg_oracle` on A256 (768 flights, stage `explore`).
* Guidance tuning round 2 (`tune_refine_n = 12`, `tune_refine_shrink = 0.25`): once round 1's 24
  configurations are scored, 12 more are drawn around its best: a numeric setting uniformly (in log space
  for a log range) within ±25% of its range's width around the best value, clipped to the range; a choice
  keeps the best value with probability 1/2; σ_f keeps "none" with probability 1/2 if the best has none.
  Seed `tune_seed + 100` (+1 for PredGuid). The center and the draws are written to
  `03_tuning/<law>/refine.json` before any round-2 flight, so a resumed run flies the same configurations;
  `configs.csv` has a `round` column. RL: `rl_tune_n = 12` for every method.
* "Within the load limit" = peak g ≤ 4.6 g (`benchmark_v5.LOAD_OK_G`: the 4.5 g limit plus the 0.1 g margin
  G3 and the tuning objective already use); column `success_5km_load` in T01, T02 and T08.
* T02 (feasible only) on ORB-1 now uses the suite's `feasible_4p5g` flag, as its caption always said (the
  flight files carry the no-load-cap flag, which is true for every case).
* Clustered bootstrap: `analyze_v5.ClusterBoot`, v4's bootstrap with one case draw shared by the four
  atmospheres; cached like the others. On synthetic data with a case effect shared across atmospheres its
  interval is about twice as wide; without one the two agree.
* G6 is "pending" until V1–V4 have all been evaluated (before: an empty gate list counted as passed).
* The analysis fingerprint hashes the content of `study_config.json` and `legacy_rl_inventory.csv` instead of
  their timestamps: the inventory is rewritten by every run, which made every "Run all" redo the analysis
  and the report (found by the smoke test: run 2 took 108 s; 24 s after the fix).
* The smoke test's stall watchdog ignores a pause of the whole machine (a laptop's sleep).

**Gates V1–V4 flown on the 12 DEV48 cases with the full study's settings (2026-11-14, this laptop, 55
dev-suite flights in 10 min; E0 needs the real v3 checkpoints and runs on Colab).** V1 passes: 12 of 12
within 2 km (largest miss 1.15 km), every Newton cycle within 6 iterations. V3 passes: K_D equals the true
drag scale in all 16 flights. V4 passes at the threshold: 11 of 12 within 3.5 km (the twelfth 3.57 km). (On Colab, 2026-11-17, the same
code scored 10 of 12: the amendment of §15, 2026-11-17.)
V2, as first written, fails, and would have stopped the notebook:

| DEV48 case | unrelieved peak (V1) | peak at amax = 3 g | miss at 3 g | relief on before the load crossed 3 g | full lift-up for the 20 s before the peak |
|---|---|---|---|---|---|
| 0 | 4.04 g | 3.81 g | 1.05 km | 17 s | yes |
| 2 | 4.12 | 3.90 | 0.46 | 17 s | yes |
| 3 | 3.99 | 3.76 | 0.28 | 17 s | yes |
| 5 | 3.75 | 3.55 | 0.32 | 18 s | yes |
| 7 | 4.22 | 3.98 | 1.09 | 16 s | yes |
| 8 | 4.06 | 3.82 | 0.64 | 17 s | yes |
| 9 | 3.82 | 3.65 | 0.47 | 18 s | yes |
| 10 | 3.73 | 3.48 | 2.84 | 18 s | yes |
| 11 | 3.71 | 3.47 | 0.28 | 17 s | yes |

The relief does what Lu's law prescribes (it comes on one prediction horizon, 16 s, ahead and saturates
at full lift-up), and it lowers the peak by 0.17–0.25 g, but on these entries (γ0 from −1.6° to −3.8°) the
load still passes 3 g. Lowering amax further does not lower the peak (case 0: 3.81, 3.83, 3.87, 3.94 g at
3.0, 2.7, 2.5, 2.3 g) and the miss grows (1.05 → 2.58 km), as in Lu's Fig. 17. Where the limit can be held
it is: smoke check 2's case stays at 3.13 g, and case 10 holds 3.10–3.14 g for about 160 s before its
final rise. With the rule of §15 (2026-11-14) V2 passes: 9 cases, 0 held, 9 infeasible with the law
saturated.

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
* **2026-11-11, pre-run review.** The session that built v5 reviewed it as a referee would and the user
  approved seven changes. No v5 test flight had been reported when they were made. A study folder started
  with the earlier settings is refused by the frozen-config check (`00_study/study_config.json`), so
  flights under old and new settings cannot mix. The changes:
* **2026-11-11, gates V1–V4 fly DEV48, not A256 (§4.0).** A failed gate stops the notebook and leads to a
  fix of the baseline; with 16 of the 64 A256 cases in the gates, that fix would have been informed by
  test cases. Now 12 DEV48 cases in nominal air (V3: 8 in thin and 8 in dense air), 11 of 12 to pass.
* **2026-11-11, `fnpeg_180`, `predguid_180`, `fnpeg_oracle` moved from Tier 2 to Tier 1 (D1's tier
  contents, §4.9).** The laws were published with a 180° bank range; with only the 90° versions flown by
  default, the baselines could be called handicapped. They stay exploratory (G7: estimates, no criterion).
* **2026-11-11, more tuning for the baselines (§4.5).** Guidance laws: 36 configurations instead of at
  most 24 (a second round of 12 around the best of the first). RL: 12 per method (the plan's cap; the
  build had used 4 to 6). v2 was designed over three studies, so the baselines' budget was the weak side.
* **2026-11-11, analysis additions; no criterion changed (§6, §8.4).** (1) A verdict is given only on
  complete data ("pending" before); for Q1 an arm with fewer than 2 misses and too few flights for the
  bound is "pending" instead of "not met". (2) Clustered-case intervals next to every pooled A256
  interval. (3) Q3 also on the load-feasible cases, within the load limit, and per target (T08b), because
  only 27 of 1,002 target-0 cases are reachable within 4.5 g at a constant bank. (4) A limitations
  section (§16) printed in the report.
* **2026-11-14, gate V2 accepts a flight on which the 3 g limit was infeasible (§4.0).** Found by flying
  the gates on DEV48 before any test flight (§14): at amax = 3 g none of the 9 high-load cases stays under
  3.15 g, although the relief engages a full prediction horizon early and commands full lift-up for more
  than 20 s before the peak. The first rule ("peak ≤ amax + 0.15 g" for every flight) tested the scenario,
  not the implementation, and would have stopped the notebook at the gates. New rule: a flight passes if
  the limit is held (then it must also land within 2 km), or if the relief was on before the load crossed
  amax, the command was full lift-up (within 1°) for the 20 s before the peak, and the peak is at least
  0.05 g below the unrelieved V1 flight's. The gate's JSON records which of the two applied to each flight.
  The study's own load check is unchanged (G3: `fnpeg` mean peak g ≤ 4.6 g per atmosphere).
* **2026-11-15, second pre-run review (the user's decisions).** No v5 test flight existed. The changes
  add to the study; none alters a flight, comparison or criterion of RQ1, RQ2 or RQ4. RQ3 gains two
  seeds per RL arm, and R-each's seed rule changes with them.
  A study folder frozen with the earlier settings is refused by the frozen-config check.
* **2026-11-15, DS1 and the simulator's density hook (D5 moved in, §3.4, §8.6).** The user asked that
  "a win is a claim about constant errors only" be dealt with without invalidating experiments 1–4.
  DS1 flies A256's 64 cases under four altitude-varying density profiles, 8 arms, 2,048 flights, its
  own Holm family {D1, D2}. It needs the first v5 change to the simulator: an optional density factor
  in `math_3d.nominal_aero_forces_from_state`, set through `data.apply_domain`, absent by default. Smoke
  check 1 flies a whole flight with a profile of ones and without one and finds them bit-identical;
  check 11 confirms that `math_3d.py` and `data.py` differ from the v4 commit by added hook lines only.
  So every v1–v4 flight and the 6,000 training flights stand as they are.
* **2026-11-15, five RL seeds (§4.4, §6, §8.3).** RQ3 compared v2 (5 seeds) with RL arms of 3 seeds and
  required 3 of 3. Every Tier 1 RL arm now trains seeds 0–4 and R-each needs ≥ 4 of 5 seeds, the rule of
  G1 and J-each. Cost: 6 more RL runs and 1,536 more policy flights (Tier 1 ≈ 26,600 flights, §4.9).
* **2026-11-16, every RL run trains on Kaggle (user decision; `jepa/kaggle_rl.md` is the detail).** Colab
  units ran short, and RL tuning and training were the largest GPU item. All of it (12 tuning configurations
  per method and seeds 0–4 of `dagger_fnpeg`, `bcppo_fnpeg`, `ppo_hist`) now trains in
  `ReEntryAI_JEPA_v5_kaggle_rl.ipynb` on Kaggle's free T4, **one platform and one GPU type for every RL run**,
  because each run's budget is wall-clock capped (3 h) and the hardware therefore shapes the result. Enforced in
  code: the study's `00_study/rl_platform.json` (GPU, torch, CUDA, platform) is fixed by the first Kaggle session
  and every later session must match it; a run never resumes on other hardware (`rl_v5.check_hardware`); full
  mode refuses any GPU but a T4; Colab's `RL_ON_KAGGLE` (full mode) never tunes or trains RL and instead exports
  the inputs (`v5_rl_inputs.json`) and imports finished runs, checked by sha256 and refusing any overwrite. No
  study setting changed, so the study folder frozen on 2026-11-16 continues. RL policy flights stay on Colab's
  CPU runtime like every other CPU-pool flight. Limitation (§16): RL trained on a Kaggle T4 with a 3 h cap per
  run; the learned models trained on Colab GPUs; each RL run's decisions reached are reported (T07). Smoke
  check 16.
* **2026-11-16, GPU passes skip CPU work (user request; notebook plumbing only, no flight, setting or
  result changes).** `run_jobs(..., can_cpu=False)` leaves the CPU-pool flights (guidance laws, classical,
  RL policies) for a CPU runtime and counts them; the notebook passes `can_cpu = not (GPU_PASS_SKIPS_CPU and
  a GPU runtime)` to every flight stage and to the legacy RL flights. Gates and tuning are prerequisites of
  the GPU stages and still run on any runtime. The last cell prints the flights left for each runtime type.
  The setting is a notebook variable, not a study setting, so the frozen `study_config.json` is unaffected.
  Smoke check 1 tests the routing.
* **2026-11-15, AS1, `fnpeg_norelief`, interval fields and the classical timing (the user's decisions
  after the audits; no v5 test flight existed, no study folder had been started).**
  (1) AS1 (§3.4b, §8.7): 128 conditions, 8 arms, 1,024 flights, own Holm family {A1, A2}. It answers the
  fidelity audit's main finding (more lift than Orion's database) with flights instead of a caveat, and
  is v5's only Mach-dependent aerodynamic error. It uses the simulator's existing schedule override;
  `math_3d.py` is unchanged and `data.py` gains default-off lines (smoke checks 1 and 11).
  (2) `fnpeg_norelief` on A256 (256 flights, G7d and G7e, estimates): the guidance audit found that our
  FNPEG's load relief trades accuracy for load without holding the limit on a Lu-analog mission, so a
  poor `fnpeg` result could be blamed on the relief; this arm shows FNPEG without it.
  (3) Every v5 flight record now stores `use_intervals`, `shield_interventions`, `shield_steps` and
  `interval_violations` from the simulator's own report (`fly_v5.interval_fields`), so "the shield never
  acted" is in the data. (4) The compute table T09 times the classical controller twice: as flown, and
  with its diagnostic interval heat screen off (`classical_noscreen`; the screen cannot change a
  command, `jepa/interval_layer_audit.md`). The flights themselves are flown as before, screen on, so
  classical results stay comparable with v3's.
  RQ1–RQ4 and DS1 are unchanged: no flight, comparison or criterion of theirs is altered.
* **2026-11-15, audits; no code behaviour changed by them.** (1) FNPEG against Lu (2014), equation by
  equation, with a Lu-analog test of the load relief, and PredGuid's constants against Bairstow (2006)
  Tables B.4 and B.7, all equal except the departures already recorded (GMAX, VFINAL, LAD, LOD,
  CD_EST_INITIAL, the 90° cap, no Earth rate); on our entries only the initial roll, ballistic, final
  and lateral logic run, so the paper should call the arm Bairstow's Apollo-derived final phase
  (`jepa/guidance_audit.md`). The scripts and results are in `jepa/audit_runs/`. (2) The simulator
  against its sources (`jepa/simulator_fidelity_audit.md`): its trim L/D exceeds Bibb's Orion database
  below about 6 km/s; kept, because changing it would invalidate v1–v4, and stated in §16. Lu's CEV
  mission itself cannot be flown here (a rotating Earth in Lu; §16). On a Lu-analog (2,450 km, 6.3 g
  unrelieved) FNPEG's relief lowered the peak by 1.0–2.7 g but held its limit in none of 17 flights,
  whether σf = σ0, 30° or 60° and whether capped at 90° or 180°; Lu's held 3 g and 4 g. The relief is
  not changed: that would be a new design chosen after seeing its behaviour. G3 and G4 already keep a
  comparison from being called load-matched where `fnpeg` exceeds the limit (§16).
  (3) FNPEG's filters equal Brunner & Lu (2008) Eq. 22; the paper gives no gain, only "slightly less
  than one". (4) Bibb Part II: our drag matches below Mach 1.5 and is 4–10% high from Mach 2 to 8; the
  L/D excess starts near Mach 2–3. (5) The interval layer (`jepa/interval_layer_audit.md`): the
  shield, interval reward and interval observation were off in every v1–v5 flight (legacy 16-feature
  SB3 policies get the interval observation); the classical baseline alone runs a diagnostic interval
  heat screen every cycle, which cannot change its command but is inside its measured step time, so
  T09's classical row is not comparable with the other arms' and is to be labelled so.
* **2026-11-17, all Colab GPU work on one GPU type, a T4 (user decision; budget, not a study setting).** The
  Colab GPU pass (E3 training, MPPI flights, open loop, GPU timing) moves from L4 to T4: about 2.7× fewer units
  per hour, and this work is a fixed amount (E3: `E3TrainCfg.steps = 36_000`; MPPI: fixed sample counts, the
  planner has no time budget), so a slower GPU changes its duration, not its result beyond rounding. Only the
  timings (T09) depend on the GPU, and they name it. Enforced in code: `GPU_TYPE = "T4"` (Settings cell, full
  mode) refuses any other GPU at setup, and the first GPU session writes `00_study/gpu_platform.json`; a later
  session on another GPU type raises (`study_v5.check_gpu_platform`), so every E3 model and MPPI arm shares one
  GPU type. torch/CUDA versions are logged per session, not locked (Colab updates them). No GPU work had run
  when this was made (only CPU passes), so the study folder frozen on 2026-11-16 continues. Smoke check 1.
* **2026-11-17, gate V4 amended to 10 of 12 (user decision; before any test flight; `study_v5.AMENDMENTS`).**
  On Colab V4 failed: 10 of 12 DEV48 PredGuid flights within 3.5 km (misses 0.41–3.95 km, all landed, phase
  sequence initial roll → ballistic → final in all 12), where this laptop had 11 of 12 (2026-11-14 and again
  2026-11-17). The cause is not a bug: PredGuid's lateral logic (App. B.61–B.67) reverses the bank on the first
  2 s guidance cycle at which the crossrange crosses its deadband, so a late reversal can move by one cycle and
  the landing by kilometres. Evidence: the same commit flown in a Linux container (Python 3.13, numpy 2.5)
  reproduced Colab's 12 misses to the metre, Windows (numpy 2.3) gave the laptop's; the two flights of case 4
  agree to float precision until `K2ROLL` flips sign at t ≈ 418 s on one and not the other; on Windows a
  relative change of 1e-8 in the entry speed (0.08 mm/s) moved case 4 from 0.48 to 1.22 km and case 7 from 0.51
  to 3.54 km, 1e-7 case 4 to 3.27 km, while 1e-12 changed nothing. pymsis 0.12 and 0.13 give identical densities;
  flying other gate flights first in the same process changes nothing. So "11 of 12 within 3.5 km" decides
  PredGuid by rounding. The alternative, dropping PredGuid (the law Orion flies) from G2, Q3, D2, A2 and G7, was
  rejected by the user. Changed: `v4_min_frac` 11/12 → 10/12; a study folder frozen at 11/12 continues and
  records the amendment in `00_study/study_config.json` (`amendments`, with this reason), any other difference
  still raises; a V4 verdict cached under the old rule is re-scored from the same flights (V4.json keeps
  `pass_before_amendment`). No flight is flown again. Consequences stated in §16 (PredGuid's per-case misses are
  rounding-sensitive; its DEV48 tuning choice is partly chance). Smoke check 1.
* **2026-11-17, all GPU work on an L4, RL included; Kaggle dropped (user decision; replaces the T4 entry
  above and the 2026-11-16 Kaggle entry).** `GPU_TYPE = "L4"` (Settings cell, full mode). E3 training, RL tuning
  and training, MPPI flights, open loop and GPU timing all run on Colab's L4 (`00_study/gpu_platform.json`, fixed
  by the first GPU session; `study_v5.check_gpu_platform` matches the type as a whole word, so an L40S is
  refused). Reasons: (1) RQ2's reference R0 trained and flew on an L4 in v4 (native bf16, TF32), so the E3 arms,
  `r0_nojepa` above all, now share its GPU type and numeric path; on a T4 they would have run emulated bf16
  without TF32 (audit B2), a confound of the arm meant to isolate the JEPA objective. (2) Every RL run now trains
  on the same GPU type as the learned models, in the notebook (`RL_ON_KAGGLE` and the notebook's Kaggle steps are
  removed; `kaggle_v5.py`, its notebook and smoke check 16 stay in the repository, retired). The RL resume lock
  (`rl_v5.HW_SAME`) now locks the GPU type and platform tag only; torch, CUDA and CPU changes are logged in the
  run's `hardware.json` (`resumes`), the policy `check_gpu_platform` already had, because Colab updates torch
  between sessions and a lock would strand a half-trained run. Cost: about 2.7× the T4's units per hour; budget,
  not a study setting. No GPU session and no RL run had started (`00_study` has no `gpu_platform.json`), so the
  study folder frozen on 2026-11-16 continues. Smoke check 1.
* **2026-11-17, a GPU pass does only GPU work and then releases the runtime (user request; notebook only).**
  On a GPU runtime the notebook skips every CPU stage: FNPEG/PredGuid tuning, the CPU-pool flights, the analysis
  and the report; it refuses to start (and so releases the runtime) if gate S0 or V1–V4/E0 is not done. A new
  cell after the last GPU stage (section 15b) writes the status, prints what the next pass must do and releases
  the runtime, so no CPU stage runs on the GPU. The one exception is RL tuning's DEV48 scoring: each RL
  configuration's 48 full-simulator flights fly in the GPU pass, between trainings. GPU work that needs a CPU
  result first waits and is counted (`mppi_filtered` and RL: FNPEG's tuning; `v2single@fnpeg_matched`: FNPEG's
  A256 flights). Smoke check 1 checks the order.
* **2026-11-17, `mppi_filtered` and the RL experts use FNPEG's tuned filter β (audit A2; before any of these
  arms ran).** §4.1 says "β tuned" and `mppi_filtered` is "FNPEG's filters" (§4.2), but the code used β = 0.9 for
  both `mppi_filtered` (`FilteredPhysicsAdapter`) and the DAgger/BC expert (`FastExpert` overwrote the tuned
  value with `RLCfg.expert_beta_1hz`). Now `mppi_filtered`'s model carries `fnpeg_filter_beta` (from
  `03_tuning/fnpeg/chosen.json`, FNPEG's default 0.9 if the chosen configuration has none), its flights wait for
  FNPEG's tuning like FNPEG's own, the adapter cache is keyed on β, and the expert keeps the tuned β^5. So G5
  compares v2 with the estimator FNPEG actually uses, and R-imit measures imitation, not an expert mismatch.
  Flight ids are unchanged (β is not part of them); no `mppi_filtered` or RL flight existed (E0 flies frozen
  filters only).
* **2026-11-17, ORB-1/DS1/AS1 bootstrap resamples the training seeds; a seed-level t-test beside every paired
  verdict (audit A1 and B1; analysis only, before any v5 comparison had data: FNPEG and PredGuid are still
  being tuned and no MPPI flight has run).** On ORB-1, DS1 and AS1 each case is flown by one seed of a learned
  arm (case mod 5). `orb_paired` resampled cases within those seed blocks but kept the 5 seeds fixed, so its
  interval ignored the variation between the 5 trained models. In a simulated null (the recipe no better than the
  comparator, `v5_code_audit.md` A1) it rejected 5.2%, 15.5%, 35.8% and 58.5% of the time on 1,000 cases at a
  per-seed spread (SD) of 0, 2.8, 5.7 and 11.3 pp. Now it is two-stage, A256's scheme: draw the seeds with
  replacement, then the cases within each drawn seed. Same null: 1.3%, 3.3%, 9.7%, 14.3% on 1,000 cases and
  0.8–3.0% on 64. With 5 seeds a bootstrap is still too narrow when the seeds differ a lot, so every paired
  comparison (A256's G, J, R families and J-equiv included, B1) also reports the t-test on its per-seed
  differences (5 units, 4 df; 3.7–5.3% in the same simulation). The verdict rules of §8 are unchanged (bootstrap,
  seed rule, Holm); T05 carries `t_lo`, `t_hi`, `t_p`, `boot_t_disagree`, and the report marks a verdict the
  t-interval does not support "(seed-level t disagrees)". Why report rather than gate: §8's rules were fixed
  before the study, the t-test is a sensitivity analysis in the style of the clustered intervals (2026-11-11),
  and a marked disagreement tells the reader exactly which verdicts rest on the 5-seed bootstrap alone. The
  bootstrap cache is keyed on the new scheme, so nothing cached under the old one is reused. Also removed: dead
  code after `as1_analysis`'s return.
* **2026-11-17, ORB-1 Q2's interval is 95% (audit C1; text corrected to the code).** §14 said "Q1/Q2 per arm with
  Clopper–Pearson (90 %)". Q2's interval (`benchmark_v5.score`) has always been two-sided 95%; Q1's bound is
  one-sided 90% (`requirement_lcb`), as §8.4 says. Q2 is descriptive; nothing else changes.
* **2026-11-17, RL learning-curve snapshots at fractions of either budget (audit C8; before any RL run).** The 5
  snapshots were placed at fractions of `budget_decisions` only, so a run that stopped at its wall-clock cap first
  would lose its later points. Now snapshot i is saved once the run has used (i + 1)/5 of its decisions or of its
  wall time, whichever comes first, and any point not reached when training ends is the final policy (the one
  flown). Each snapshot records its trigger (`decisions`, `wall_clock` or `final`), its decisions and its elapsed
  time; F05 plots decisions. Snapshots are measurements and do not affect training.
* **2026-11-17, three defensive fixes that change no result (audit C4, C6, C9).** (C4) PredGuid's final-phase
  design cache is keyed on the whole configuration, not six of its fields: the design's dry run also reads gmax_g,
  overcontrol, Q7F, Q19 and K_ρ, which tuning varies. The old key was safe only because every tuning configuration
  draws a different KA; flights before and after are bit-identical (6 DEV48 conditions × 6 configurations in one
  process). (C6) `legacy_jobs` raises if two legacy models would get the same arm name (the study's 5 have
  distinct names). (C9) `models_v2.amp_dtype`'s docstring says fp16 on a T4; recent torch gives a T4 emulated bf16 (the L4 native bf16). The correction is a comment at its use in `models_v5.train_e3`, because `models_v2.py` stays byte-identical to the v4 commit (smoke check 11).
* **2026-11-17, PredGuid's tuning is effectively 4-D (audit C3; documentation only).** On the tuning's own ground
  truth (all 48 DEV48 conditions, 12 cases × 4 atmospheres, flown as `tune_guidance` flies them), Q7F, Q19, KA
  and K_RHO_FILTER_GAIN, all four at the low or all at the high end of their ranges, leave every flight
  bit-identical to the defaults, except the logged density estimate (`dec_Krho_est`, `predguid_Krho_final`), which
  no phase flown on DEV48 reads. LATBIAS at its high end, the positive control, changes 39 of 48. Every DEV48 flight
  goes initial roll → ballistic → final. `PREDGUID_SPACE` is not changed, because tuning had started (§14 says
  what this means for T10).
* **2026-11-18, unattended runs: watchdog, failure records, session budget, early stops, RL last and gated
  (Part A; no flown flight, tuned law or trained model changes).** The user leaves GPU passes running unattended,
  so a pass must never idle on paid units and one broken arm must not stop the rest. (1) A **stall watchdog**
  (`jepa/runguard_v5.py`, started at Setup): every study log line, Drive write and MPPI round counts as progress;
  after `STALL_MIN` (45) minutes without any it writes the reason to `00_study/logs/watchdog.log` and
  `RUN_SUMMARY.txt`, then releases the runtime (Drive flush, then `runtime.unassign()`, each bounded in time). The
  known long silent step, copying and loading the 6,000 flights, runs under a 60-minute allowance.
  (2) **`00_study/RUN_SUMMARY.txt`**, rewritten after every stage: stages, failures and why, units used, why the
  runtime was released, what the next pass must do. (3) **Failures are recorded and skipped**: each stage cell runs
  in `RG.stage`; an exception is logged (`logs/failures.log`, the summary) and the notebook goes on;
  `STOP_ON_FAILURE = True` releases at the first one. A flight that raises (or lands with a non-finite miss) is
  retried once (a lockstep batch flight by flight; the laws are batch-independent, §16) and otherwise recorded in
  `04_flights/errors/<flight>.json`; a stage call whose first 50 flights have more than 10 % errors stops and writes
  `01_checks/stage_<label>_failed.json`. Poor accuracy is never an error. (4) **E3 early stop**
  (`models_v5.stop_rule`): a non-finite training or validation loss, a validation loss above 3× its best at the same
  curriculum horizon, or, from 20 % of the steps on, no validation loss at the first horizon below the untrained
  model's, stops that model; it gets `failed.json` (status `diverged`), the status table shows it, later sessions
  skip it, and the other models go on. A learning-rate pilot that stops loses the choice. Being worse than R0 never
  stops a model. Models that train normally are bit-identical to before (the extra untrained validation loss is
  computed without gradients or random draws). (5) **`MAX_SESSION_UNITS`** (None by default) becomes a deadline at
  the session's rate (`UNITS_PER_HOUR`, default L4 1.98, the rate Colab shows; the code's table had 4.8), and
  `SessionLimit` now also gates every E3 model, pilot set and flight group: nothing starts that would not fit; when
  the budget is used up the summary is written and the runtime released. (6) **Order of a GPU pass**: E3 training,
  E3 flights, the MPPI flights (A256, matched, ORB-1, DS1, AS1), open loop, timing, then RL last (cell 15), so a
  slow or failing RL arm never holds up the rest; DS1's and AS1's RL flights, which wait in cells 11 and 11b for the
  policies, fly in cell 15a. (7) **RL tunes and trains only after its competence gates pass**
  (`01_checks/RV1_<arm>.json`, `RV2`, `RV3`; DAgger RV1 + RV2, BC → PPO RV1–RV3, PPO RV3); until then each arm
  waits with a note. The gates and the supervised-fit fix are the next entry's subject (Part B). (8) **DAgger's
  first tuning is deprecated, not deleted.** Configurations 0–9 scored 0–2.1 % success@5km against FNPEG's 100 %;
  `fit_supervised` underfits (sup_batch 4096 × 20–50 epochs is about 10–250 gradient steps, and BC used PPO's
  learning rate); an offline refit (300 epochs, lr 10⁻³) gave 0.86 sign accuracy and 7° bank error. So those
  results measure the bug, not DAgger. `03_tuning/rl_dagger_fnpeg/` and `04_flights/{results,telemetry}/DEV48/
  tune_dagger_fnpeg_c00…c10/` stay exactly as they are, each with a `DEPRECATED.md` (written once by Setup);
  repair and the analysis skip them and the report lists them here. The refit tunes under new keys:
  `03_tuning/rl_dagger_fnpeg_fit2/` and flights `tune_dagger_fnpeg_fit2_cNN`. No test flight of any RL arm had
  flown. (9) **Notebook revision** (`NOTEBOOK_REV = 2`): a notebook whose Setup does not check this code's revision
  (an old downloaded copy that pulls the new code) stops at its first study cell with "open the notebook from
  GitHub", and its crash hook releases the runtime. The auto-disconnect hook is now registered at the top of Setup,
  so a failure anywhere in Setup releases the runtime too. Tested by smoke check 17 and by
  `ReEntryAI_watchdog_test.ipynb` (5 minutes on a CPU runtime: a cell hangs, and the watchdog must release the
  runtime from its background thread).
* **2026-11-18, RL competence gates, the supervised-fit fix and CPU routing (Part B; no flown flight, tuned law or
  trained model changes; no RL policy had been trained).** (1) **The supervised fit counts steps, not epochs**
  (`rl_v5.fit_supervised`): `sup_steps` (4,000) Adam steps at its own rate `sup_lr` (10⁻³; BC used PPO's 3·10⁻⁴),
  minibatches of `sup_batch`; whole episodes (`sup_val_frac` 20 %) are held out, the held-out loss is read every
  `sup_eval_every` (100) steps, the fit stops after `sup_patience` (10) reads without improvement and keeps its best
  weights; a non-finite loss raises. `dagger_epochs` and `bc_epochs` stay in `RLCfg` (unused) so old run files still
  load. DAgger's refit tunes `sup_lr` (log-uniform 3·10⁻⁴–3·10⁻³), `sup_steps` {2,000, 4,000, 8,000} and the width,
  in `03_tuning/rl_dagger_fnpeg_fit2/`. (2) **More expert flights**: a DAgger tuning iteration labels at least 100
  episodes (`rl_tune_frac` alone gave 20, about 2,500 states), and the FNPEG labels run on every core
  (`ParallelExpert`, `label_workers`; one FNPEG flight per env row, so the labels equal the single-process ones,
  smoke check 18). (3) **Competence gates** (`01_checks/RV*_<arm>.json`, pilots in `01_checks/rl_pilots/`, seed 0,
  the arm's settings with a smaller count): **RV1** held-out sign accuracy at least 5 points above the better of the
  two trivial baselines (always positive; the current bank's sign) and a held-out |bank| error of at most 10°;
  **RV2** the policy alone (no expert) reaches at least 50 % of FNPEG's fast-env success@5km (pilots: DAgger 4
  iterations × 100 episodes, BC 400 episodes); **RV3** in a PPO pilot (0.5 h or 2·10⁷ decisions) the return rises
  (last third above first third), success@5km > 0 at some point, no NaN. DAgger needs RV1 and RV2, BC → PPO RV1–RV3
  (its RV3 pilot starts from its clone), ppo_hist and sac_hist RV3, ppo_matched ppo_hist's. A failed gate is final
  for the arm (it is reported, not tuned or trained); deleting the gate file reruns its pilot. (4) **An arm whose
  first 3 tuning configurations all score 0 % success@5km on DEV48 stops** (`03_tuning/<method>/stopped.json`): no
  choice, no training, reported. (5) **Training checks** at 25 % and 50 % of a PPO run's budget (decisions or
  hours, whichever is further along): non-finite weights or returns stop it; at 50 %, no flight within 5 km in any
  learning-curve row and no rise in the return stops it. A stopped run gets `failed.json` (status `stopped`), later
  sessions skip it, the other seeds go on. Poor accuracy never stops a run. (6) **CPU routing**: DAgger (labelling,
  fits, DEV48 scoring, tuning and its 5 seeds) and BC's expert dataset and clone (`clone.pt`, `clone_fit.json`,
  `clone_hardware.json`) run on a CPU pass (cell 6b); BC → PPO's PPO stage, ppo_hist (and sac_hist) on the L4
  (cell 15); Tier 2's ppo_matched (full-simulator PPO, CPU-bound) on the CPU pass after the L4. Each run's
  `hardware.json` records where it trained. BC → PPO's tuning configurations start from seed 0's clone (the clone
  does not depend on the tuned PPO settings). The GPU arms still share the one GPU type (§16). (7) **Cost
  projection**: the pilots' wall times give `01_checks/rl_cost_projection.json` (hours per tuning configuration and
  per seed, on the runtime each part runs on). The notebook's revision is now 3 (`NOTEBOOK_REV`), so a copy
  opened before this change stops at Setup. (8) **Laptop pilots (2026-11-19, a preview: the study's own
  gate files come from Colab's CPU pass, cell 6b; full study settings, FNPEG's tuned configuration, Intel 12-core
  CPU, 151 min).** *DAgger fails RV1 and RV2*: after 4 iterations × 100 episodes (41,300 labelled states), the
  held-out sign accuracy was 0.909 against the current bank's sign at 0.894, the |bank| error 10.2°, and the policy
  alone reached 0.6 % success@5km in the fast env against FNPEG's 87.1 %. *BC → PPO fails RV1 and RV2*: on 400 expert
  episodes (40,100 states), the clone's held-out sign accuracy was 0.968 against 0.928 and its |bank| error 5.1°, but
  alone it reached 2.7 % against FNPEG's 80.6 %. In both fits the held-out loss is lowest early (step 200 or 500
  of 1,200–1,500) while the training loss keeps falling, so the fit does not carry over to unseen dispersions. The
  likely cause (not yet checked) is that the wrong signs fall at the bank reversals, which set the crossrange,
  so a 3–9 % sign error is enough to fly the vehicle off course.
  Neither arm tunes or trains; both are reported as failed gates. *ppo_hist passes RV3*: in 30.5 min (3,722
  decisions/s on the CPU), the return rose from −30.3 to −3.9 with no NaN, and its best success@5km was 0.2 %
  (above 0 as RV3 requires, but far from FNPEG). Cost projection (`rl_cost_projection.json`): DAgger, had it
  passed, ~2.8 h per tuning configuration and ~5.6 h per seed on the CPU (~62 h, ~16 units); BC's clone ~4.5 h
  per seed on the CPU; ppo_hist ~0.35 h per configuration and 3 h per seed (its wall cap) on the L4 (~19 h,
  ~38 units).
* **2026-11-19, RQ3 moves to notebook 2 with new arms (user decision; written before any notebook-2 tuning,
  training or flight).** The brief and the build record are `jepa/plan_rl_v5.md`; the notebook is
  `ReEntryAI_RL_v5.ipynb` (`jepa/build_notebook_v5_rl.py`, logic in `jepa/rl2_v5.py` and `jepa/study_rl_v5.py`),
  writing to `MyDrive/ReEntryAI_rl_v5` and reading this study's folder read-only. (1) **Arms.** `dagger_fnpeg` and
  `bcppo_fnpeg` are dropped after the Part B laptop pilots with the fix of commit 5471ceb (DAgger RV1 fail: 0.919
  vs 0.909, 9.7°; RV2 fail: 0.014 vs FNPEG 0.805; BC → PPO RV1 pass: 0.957 vs 0.818, 5.8°; RV2 fail: 0.031 vs
  0.805) and are reported as failed gates. They are replaced by `resid_fnpeg` (residual RL on the tuned FNPEG:
  σ = s_F clip(|σ_F| + b a, 0, 90°), FNPEG keeps the sign and the reversals, zero-initialized last layer so the arm
  is FNPEG at init; Johannink et al. 2019, Silver et al. 2018) and `liu_track` (Liu-style reference tracking after
  Liu et al., Adv. Space Res. 77 (2026): bank magnitude only, reversals by FNPEG's tuned heading deadband, Liu's
  reward with a 4.5 g limit and a 5 km terminal bonus; every adaptation listed in plan_rl_v5.md §1.3). `ppo_hist`
  stays as preregistered. (2) **RQ3 Holm family** {R-resid_fnpeg, R-liu_track, R-ppo_hist} (was {R-dagger_fnpeg,
  R-bcppo_fnpeg, R-ppo_hist}); the size stays 3 and an arm that fails its gate enters with p = 1. R-each as §8.3
  (`v2single − arm` on A256, ≥ 4 of 5 seeds). R-imit is replaced by the estimates `arm − fnpeg` (A256, DS1, AS1,
  ORB-1 block 1). (3) **Gates.** RI (`resid_fnpeg` at init equals FNPEG bit for bit in the fast env and on DEV48 in
  the full simulator; a failure is a stop) and RP (an untuned pilot of 2·10⁷ decisions or 1 h on the CPU, judged on
  512 fixed fast-env episodes against FNPEG alone on the same episodes: `resid_fnpeg` ≥ FNPEG − 2 pp with a rising
  return, `liu_track` ≥ 50 % of FNPEG, `ppo_hist` the old RV3 **and** ≥ 50 % of FNPEG: RV3 is raised because a pilot
  passed it at 0.2 %). (4) **Budget rule.** A run's budget is its decision count (`ppo_hist`, `liu_track` 2·10⁸;
  `resid_fnpeg` 4·10⁶, bound by FNPEG's CPU labels); the wall clock is a safety cap of max(3 h, 1.5 × the time
  projected at the routed runtime's measured rate), so `ppo_hist`'s preregistered 3 L4-hours are unchanged. Each arm
  trains on the runtime with fewer units per million decisions (pilot rate on the CPU, a 5-minute timing run on the
  L4; `01_checks/route_<arm>.json`). (5) **Stop rules** (plan_rl_v5.md §6): technical failures write
  `00_study/STOPPED.json` and release the runtime; the next session needs `ACK_STOP`. A failed RP is a result.

---

## 16. Limitations

Stated on 2026-11-11, before any v5 result. The report prints this section.

**The model errors are simple.** A256's and ORB-1's conditions change the truth by constant
multipliers on density, lift and mass. FNPEG's two filters estimate exactly that, and its predictor
uses the simulator's own atmosphere table and aerodynamic model, so the classical laws are tested in
their best case there. Since 2026-11-15, DS1 (§3.4) adds density errors that change with altitude,
in their own Holm family so that RQ1–RQ4 are unchanged, and AS1 (§3.4b) one Mach-dependent aerodynamic
error (the published Orion lift and drag). v5 still does not test winds, navigation errors or sensor
noise, and its results do not carry over to them.

**RL runs are capped in wall time** (added 2026-11-16, rewritten 2026-11-17 and 2026-11-18, §15). Every GPU RL
run (BC → PPO's PPO stage, ppo_hist, sac_hist) tuned and trained on the study's one GPU type, an L4
(`00_study/gpu_platform.json`), the GPU type of the E3 models and of v4's R0; DAgger (no wall-time cap: a fixed number
of iterations and fit steps), BC's clone and ppo_matched trained on Colab's CPU runtime. Each PPO run stops at its
decision budget or its wall-time cap (3 h; `ppo_matched` 24 h), whichever comes first. The RL arms are comparable with each other; the cap depends on the hardware, so each run's decisions
reached are reported (T07), and a faster GPU could have given the RL arms more training within the same cap.

**Single landings depend on rounding** (added 2026-11-17, §15; extended the same day, `v5_code_audit.md` B3).
PredGuid's bank reversals fire on the first guidance cycle the crossrange crosses a deadband, so near the end of a
flight a microscopic difference (a platform's last-digit rounding, or 1e-8 of the entry speed) can move a reversal
by one 2 s cycle and the landing by up to about 3 km. The other arms are less sensitive but not immune: FNPEG's
landings moved by up to 0.43 km on one platform from 1e-8 of the entry speed (0.72 km across Windows and Linux),
an MPPI flight with a learned model by about 0.5 km (measured with a smoke-test model), and the classical
controller by under 0.01 km. An MPPI flight also depends on which flights share its lockstep batch (0.55 km for
the smoke-test model, 0.25 m for `mppi_filtered`), so MPPI flights re-flown after a disconnect in the middle of a
model run differ from an uninterrupted run. Resume is bit-exact for the guidance laws and the classical
controller (a flight alone and in a batch are identical); an RL policy's commands can differ in float32's last
digit with its batch (4e-7° in a check, the landing unchanged). Per-case misses are therefore partly chance. Paired comparisons over 256
or more conditions average this out (the noise has no direction), but single-case results, the DEV48 tuning
choices of PredGuid and FNPEG among near-tied configurations (12 cases per configuration), and success near the
5 km threshold carry it. Gate V4 was amended for this reason.

**Five seeds** (added 2026-11-17, §15). Each learned arm has 5 training seeds, and a verdict about a training
recipe rests on those 5. The pre-registered bootstraps resample them, but with 5 clusters a bootstrap interval is
too narrow when the seeds differ: in a simulated null with a per-seed spread of 5.7 or 11.3 pp, ORB-1's
two-stage bootstrap rejected 9.7% or 14.3% of the time instead of 5%, and A256's about 11% at 11 pp
(`v5_code_audit.md` A1, B1). The seed-level t-test beside each verdict (T05) holds its 5%; where it disagrees
with the bootstrap, the report says so.

**Heating is evaluated at the nominal density** (added 2026-11-17, `v5_code_audit.md` C7). The density knobs act
on the aerodynamic forces only: `density_scale` through the reference area (since v1, so v1–v4 are the same) and
DS1's altitude profile through the density in the equations of motion. The heating model
(`math_3d.nominal_heating_envelope_from_state`) and the RL heat-margin feature read the standard atmosphere's
density at the altitude flown. So peak heat rate and heat load in every table are the heating of the trajectory
each arm flew, evaluated at nominal density, not the heating in the perturbed atmosphere: convective heating goes
as √ρ, so a ±12% density error (A256's mismatch atmospheres) is about ±6% in heat rate, DS1's ±15% shears about
±7%. On ORB-1 `density_scale` also carries the drag-coefficient dispersion, which should not change heating at
all. No criterion uses heating. It was kept for continuity with v1–v4 (the simulator is not changed in v5).

**PredGuid is tuned in 4 of its 8 settings** (added 2026-11-17, `v5_code_audit.md` C3, §14). DEV48's entries
never reach PredGuid's HUNTEST or up-control phases, so Q7F, Q19, KA and K_ρ's filter gain cannot change a DEV48
score. The tuning chose LATBIAS, KLAT, overcontrol and GMAX, and the other four are random draws from their
ranges. A test entry that does reach those phases flies them with untuned values. Each flight records its phase
sequence (`predguid_phases`).

**One simulator, no flight data** (added 2026-11-15; `jepa/simulator_fidelity_audit.md`; AS1, §3.4b,
tests the lift excess named here, with values read off published figures). Every
result comes from a 3DOF model of an Orion-class capsule over a non-rotating spherical Earth, with
perfect navigation, one NRLMSIS density profile, and a trim aerodynamic schedule whose L/D exceeds
the published Orion database (Bibb et al. Part I) below about 6 km/s (0.40 against 0.28 at FMV 8;
on ORB-1's nominal entry, 23% more range authority between 0° and 60° bank). Comparisons between
controllers are paired within this simulator; absolute success rates, loads and heating are not
predictions for Orion. The learned models were trained on the same simulator. The schedule was kept,
because changing it would invalidate v1–v4.

**ORB-1 is an adaptation, not a published benchmark.** The entry state is Lu's and the dispersions are
Bairstow's, but the Earth does not rotate, the bank is capped at 90°, and the flight-path angle, the
three target ranges and the speed dispersion are ours. Its numbers cannot be compared with published
Monte Carlo results, and Q1 applies Orion's requirement to our simulator, not to Orion. Lu's own CEV
mission (3,769 km) is out of our vehicle's reach at any γ0 from −1.0° to −2.0° (checked 2026-11-15:
full lift-up lands 16–880 km short), mainly because Lu's Earth rotates.

**Target 0 mixes accuracy and load.** Only 27 of its 1,002 cases are reachable within 4.5 g at a constant
bank. A law that enforces the limit gives up range there; MPPI's limit is soft and can be exceeded. The
pooled ORB-1 results include target 0; T08b and the "within the load limit" columns separate the two.

**The laws enforce the load limit differently.** FNPEG's load relief looks 16 s ahead and can only add lift;
on the dev suite it lowered the peak load by about 0.2 g and could not hold 3 g on any of 9 high-load
entries (§14). On a Lu-analog mission it lowered the peak by 1.0–2.7 g but held 3 g or 4 g in none of
17 flights, where Lu's (2014) held both (`jepa/guidance_audit.md` §3); the audit verifies FNPEG's
equations, not Lu's closed-loop load performance. MPPI plans the whole trajectory with a soft penalty. G3 checks whether FNPEG stays near
4.5 g on A256; where it does not, the load-matched comparison (G4) is flagged invalid.

**Half the bank range.** Every controller is limited to 90° because the learned controllers were trained
that way. FNPEG and PredGuid were designed for 180°. The 180° versions fly as exploratory arms (G7), not
as the primary comparison.

**Unequal tuning history.** v2 and its planner were designed over three studies, and A256 is the test
set v3 and v4 already reported on. The baselines are tuned on DEV48 only (12 entry cases, 36
configurations per guidance law, 12 per RL method). ORB-1 is new to every arm.

**A256 has 64 entry cases.** Its 256 conditions reuse them in four atmospheres. The pre-registered
bootstrap treats the conditions as independent; the clustered intervals show how much that matters.

**RQ2's key arm changes two things.** `r0_nojepa` removes the JEPA losses and also trains end to end
instead of in two stages, and each non-JEPA model's learning rate is picked from three values on one
seed. A difference between it and R0 cannot be assigned to the objective alone. Both train and fly on an L4
(R0 in v4, the E3 arms in v5, §15 2026-11-17), with the torch and CUDA versions of their sessions.

**RL trains in a simplified simulator.** The policies learn in the fast environment (the physics port)
and are tested in the full simulator, with five seeds per method (three before 2026-11-15); the gap
is reported per arm. PPO in
the full simulator is Tier 2 and limited to 6,000 episodes.

**Compute numbers describe these implementations.** FNPEG's predictor is Python and numpy, the learned
models run on a GPU. The table compares what was run here, not flight software.
