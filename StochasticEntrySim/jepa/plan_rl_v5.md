# SpaceJEPA v5, notebook 2: RL baselines that can fly (plan and build record)

Written 2026-11-19 from the design agreed with the user, after the laptop pilots of the imitation arms (commit
5471ceb). This file is the brief for notebook 2 and the record of how it was built. The study's full specification
stays `jepa/v5_plan.md`; section numbers below (§4.4, §8.3, ...) point there. Where this file and `v5_plan.md`
disagree, this file wins only on the points marked **change**; each of them is a `v5_plan.md` §15 entry, written
before any RL flight or training of notebook 2 exists. Notebook 1 (`plan_jepa_needed_v5.md`) has priority on the
L4; this notebook runs after it.

---

## 0. The question and how the answer is decided

**RQ3 (§7):** how close do properly built RL baselines get to v2 and to FNPEG? v5's first RL design could not answer
it: the imitation arms of FNPEG failed their competence gates (§2), and PPO from scratch passed the old gate RV3
with 0.2 % success@5km in the fast environment. Notebook 2 replaces the failed arms with two arms that **never let
the policy decide the bank-reversal sign** (the likely cause of the imitation failures: small errors in reversal
timing compound, §2), and keeps PPO from scratch as preregistered.

| # | Arm | What it asks | Criterion |
|---|---|---|---|
| 1 | `resid_fnpeg` | Can RL improve a tuned published law by correcting its bank magnitude? | R-each (vs v2) and the estimate `resid_fnpeg − fnpeg` |
| 2 | `liu_track` | Can RL track a reference trajectory designed for this capsule (Liu et al. 2026's framework)? | R-each and `liu_track − fnpeg` |
| 3 | `ppo_hist` | Can RL from scratch, with the world models' memory, fly the mission (preregistered)? | R-each and `ppo_hist − fnpeg` |

**How the answer reads** (fixed now):

| Outcome | Written in the report |
|---|---|
| An arm fails its pilot gate (§3) | "did not reach the gate"; its gate numbers are the result; no L4 time spent on it |
| R-each passes for an arm (v2 > arm) | v2 beats that RL baseline |
| R-each does not pass and `arm − fnpeg` ≥ 0 | the RL baseline is as good as the best classical law: RL competes with v2 there |
| `resid_fnpeg − fnpeg` > 0 on AS1 or DS1 | learning corrects a model-based law where its model is wrong |

Every outcome is reported with the same tables, including nulls, losses and failed gates.

---

## 1. Arms

Common to all three (§4.4, unchanged): 5 s decisions, the bank command held through the simulator's actuator
(20°/s, |σ| ≤ 90°), observation = `rl_env`'s 13 features + the last 24 decisions of [h, V, γ, a_D, cD_meas, cL_meas,
σ] (181 numbers, normalized by a running normalizer saved with the policy), training in the batched fast
environment (`jepa/rlenv_fast_v5.py`) on the world models' training distribution (`EnvelopeSpec`, `WIDE_DOMAIN`),
evaluation in the full simulator. In-house PPO (CleanRL-style: GAE 0.95, γ 0.99, clip 0.2, 4 epochs).

### 1.1 `resid_fnpeg`: residual RL on FNPEG

References: Johannink et al. 2019, "Residual Reinforcement Learning for Robot Control" (ICRA); Silver et al. 2018,
"Residual Policy Learning" (arXiv:1812.06298). The user may supply the papers; nothing below depends on details the
abstracts do not state.

* **FNPEG flies as usual**: the tuned configuration (`MyDrive/ReEntryAI_jepa_v5/03_tuning/fnpeg/chosen.json`,
  config 28: σf 34.0°, deadband 1.32° + 3.97·10⁻⁴ °/(m/s) V, δ 0.0235, k0 50.7, β 0.882, RK4 step 1 s, Vf 147.7 m/s),
  its own filters, sign, reversals, hold and load relief. Its memory is never synced to the command flown (the
  flown sign is always its own).
* **Command:** σ = s_F · clip(|σ_F| + b·a, 0, 90°), with σ_F FNPEG's command, s_F its sign (+ when σ_F is 0), a ∈
  [−1, 1] the policy's output (Gaussian PPO; deterministic mean in every evaluation) and **b the bound, tuned in
  {10, 20, 30}°**. FNPEG keeps the sign and every reversal.
* **Observation:** the 181 common features + |σ_F|/90°, s_F, K_D, K_L (FNPEG's filter state): 185.
* **Initialization:** the actor's last layer has zero weights and bias, so its mean is exactly 0 and **at init the
  arm is FNPEG, bit for bit** (gate RI, §3). The exploration noise starts small (std tuned in {0.1, 0.3} of b), and
  the critic trains alone for the first updates (critic burn-in, as in residual policy learning) so the first actor
  updates do not follow a random critic.
* **Cadence:** in the full simulator FNPEG decides every second, exactly as the `fnpeg` arm; the correction b·a is
  recomputed every 5 s and held. In the fast environment FNPEG decides at the 5 s decisions with its filter at β⁵
  (DAgger's expert, `rl_v5.FastExpert`). The difference is part of the fast-vs-full gap (R-gap).
* **Reward:** the fast environment's (minus MPPI's realized cost plus range shaping, as `ppo_hist`).

### 1.2 `liu_track`: Liu-style reference tracking

Reference: Liu, Lei, Xian, Ren, Li, Zhang, "A deep reinforcement learning-based reentry guidance framework with no
fly zone constraints", Adv. Space Res. 77 (2026) 1031–1049 (read from the extracted text; equation numbers are the
paper's). The arm is called "Liu-style": every adaptation is a deviation (§1.3).

* **Reference velocity vector** (Eqs. 16–21), per case, built at the start of each flight:
  v_ref = v0 − (v0 − vf) exp(−s/δs), with v0 the entry speed, s the range to go, s0 its value at the start,
  **δs = k·s0, k tuned in {0.2, 0.3, 0.4}** (Liu: values in this range "worked well"), vf = FNPEG's tuned Vf at
  hf = 7.6 km (the drogue altitude); γ_ref = arcsin((hf − h)/s) (clamped to ±90°); ψ_ref = ψ if |ψ − ψ_LOS| ≤ δψ,
  else ψ_LOS (Eq. 19); V_ref = v_ref [cos γ_ref cos ψ_ref, sin γ_ref, cos γ_ref sin ψ_ref]. For a lifting capsule in
  equilibrium glide, V² ≈ gR (1 − e^(−2s/(R L/D))), which Eq. 16 approximates within ~10 % for k ≈ 0.3 at our L/D, so
  the reference is flyable.
* **Action:** bank magnitude only, |σ| = 45° (a + 1), a ∈ [−1, 1].
* **Reversals:** not the policy's. The sign follows Liu's heading-error deadband (Eq. 19): it turns toward the line
  of sight once |ψ − ψ_LOS| exceeds δψ and is kept inside it. **δψ is FNPEG's tuned deadband** (linear in V), so the
  lateral logic is a published, tuned one and the policy only flies the longitudinal channel.
* **Reward** (Eqs. 31–37, Liu's coefficients from Table 2): r = r_shaping + r_alt + r_ctrl + r_path + r_rev + r_bonus,
  * r_shaping = 0.25 exp(−‖V − V_ref‖² / κ_v²), κ_v = 0.1 v_ref;
  * r_alt = 0.05 exp(−(sin γ − sin γ_ref)² / κ_a²), κ_a = 0.1 max(|sin γ|, sin 0.5°);
  * r_ctrl = −0.02 (Δ|σ_cmd| / 90°)²;
  * r_path = −0.05 exp((n − n_max)/n_max), **n_max = 4.5 g**;
  * r_rev = −0.1 × (reversals in this decision);
  * r_bonus = 10 when the flight crosses the drogue altitude **within 5 km** of the target.
* **Training terminations** (Liu §3.2.2): the fast environment's own (drogue altitude, near-vertical, skip-out,
  1,400 s), plus load above 4.5 g and |ψ − ψ_LOS| > 90°. Evaluations (gates, fast eval, full simulator) never stop a
  flight on these two: they only end where physics ends it.
* **Observation:** the 181 common features + Liu's Eq. 22 terms not already in them: (V − V_ref)/v_ref (3), s/s0,
  v/v_ref − 1, (ψ − ψ_LOS)/π, (h − hf)/100 km, n/4.5 g, γ_ref, the current sign, reversals/10: 192.

### 1.3 Deviations of `liu_track` from Liu et al. (each one stated in the report)

| Liu et al. | Here | Why |
|---|---|---|
| CAV-H glider, commands angle of attack and bank | Orion-class capsule at fixed trim; bank magnitude only | the capsule has no angle-of-attack channel |
| policy outputs bank rate, sign included | policy outputs \|σ\|; sign from the Eq. 19 deadband | design rule of notebook 2: the policy never decides reversals |
| constant δψ | FNPEG's tuned V-linear deadband | the published constant is for a glider; ours was tuned for this vehicle |
| GRU policy and value networks | MLP on the 24-decision history | the same memory as every other v5 RL arm |
| PPO with adaptive learning rate and clip (Eqs. 38–39) | v5's PPO (fixed clip, tuned lr) | one PPO for every arm; the adaptive scheme is listed as untested |
| path penalty on heat rate, dynamic pressure and load | load only (4.5 g) | v5 constrains load only |
| terminal set: range-to-go and energy (Eq. 7) | drogue crossing within 5 km | v5's success criterion |
| LOS perturbation during training (Eqs. 24–28) | not used | it trains no-fly-zone evasion, which v5 does not test |
| reading of κ_alt (garbled in the extracted text) | 0.1 max(\|sin γ\|, sin 0.5°) | the extraction reads "0.1 sin γ"; the floor avoids a zero scale |
| Adams integration at 0.1 s | the fast environment's (2, 2, 1) s Euler substeps | v5's training environment |

### 1.4 `ppo_hist`

Unchanged from §4.4 and the existing code (`rl_v5.ppo_train`): PPO from scratch, Gaussian over the signed bank
command, tuning space of §4.5. Only its gate changes (§3).

---

## 2. Dropped arms: negative results (reported with their gate numbers)

Local pilots with the Part B fix (commit 5471ceb: DAgger's FNPEG follows the command flown, DART noise for BC, RV2
against FNPEG on the same episodes), full study settings, FNPEG's tuned configuration, Intel 12-core laptop.
Pilot folders: `pilot_dagger2` and `pilot_bc2` (session scratchpad of the Part B session).

| Arm | RV1 (held-out sign accuracy vs best trivial baseline; \|bank\| error; need +0.05 and ≤ 10°) | RV2 (policy alone vs FNPEG, success@5km in the fast env; need ≥ 50 %) |
|---|---|---|
| `dagger_fnpeg` | **fail**: 0.919 vs 0.909, 9.7° | **fail**: 0.014 vs 0.805 |
| `bcppo_fnpeg` | pass: 0.957 vs 0.818, 5.8° | **fail**: 0.031 vs 0.805 |

DAgger's learner-plus-expert success fell from 0.871 (iteration 0, expert only) to 0.149 by iteration 3 while its
held-out sign accuracy stayed above 0.92: a few wrong signs at the reversals are enough to fly off course. Neither
arm tunes, trains or flies in notebook 2. The report lists both as "failed the competence gate" with these numbers
(and the study's own RV files if the Colab CPU pass of notebook 1's predecessor wrote them).

---

## 3. Gates (cheap CPU pilots before any L4 time)

Every gate runs on a **CPU runtime** in the fast environment, seed 0, the arm's untuned settings (tuning comes after
the gates). The pilot budget is **2·10⁷ decisions or 1 h, whichever comes first** (the same for every arm).
Competence is measured on **512 fixed fast-environment episodes** (seed 12345; each environment row flies exactly
one episode, so every arm and FNPEG fly the same 512 entries, parameters and targets). FNPEG alone on them is the
reference (`01_checks/rl2_pilots/fnpeg_fixed.json`, computed once).

| Gate | Arm | Passes if |
|---|---|---|
| **RI** (init identity) | `resid_fnpeg` | (1) the initialized actor's mean is exactly 0 on the fast environment's observations, and its command equals FNPEG's bit for bit; (2) on DEV48 (48 conditions), `resid_fnpeg@init` and `fnpeg`, flown in the same session, have identical telemetry arrays and identical miss, outcome and peak g |
| **RP** (pilot) | `resid_fnpeg` | pilot success@5km ≥ FNPEG's − 2 pp on the 512 episodes, the return rises (last third of the learning curve above the first third), no NaN |
| **RP** | `liu_track` | pilot success@5km ≥ 50 % of FNPEG's on the 512 episodes, no NaN |
| **RP** | `ppo_hist` | old RV3 (return rises, success > 0, no NaN) **and** ≥ 50 % of FNPEG's success (**change**: RV3 raised; a laptop pilot passed the old RV3 with 0.2 %) |

* A failed RP is a **result**: the arm is reported with its gate numbers and gets no tuning, training or flights.
  It never stops the notebook.
* A failed RI is a **technical failure** (the code does not do what it claims): stop rule (§6).
* Gate files: `01_checks/RI_resid_fnpeg.json`, `01_checks/RP_<arm>.json`; deleting one reruns its pilot.

**Cost per million decisions.** Each pilot records decisions per second on its CPU runtime. On the first L4 pass, each
arm that passed its gate runs a 5-minute timing run on the L4 (`01_checks/rl2_rates/<arm>_L4.json`). Units per
million decisions = 10⁶ / (decisions/s) / 3,600 × the runtime's units per hour (L4 1.98; High-RAM CPU 0.26). The
arm's tuning and training then run on the cheaper runtime, frozen in `01_checks/route_<arm>.json` once both rates
exist. `resid_fnpeg` is bound by FNPEG's CPU labels, so the CPU is expected to win for it; `ppo_hist` and
`liu_track` are GPU-batched, so the L4 is expected to win for them.

---

## 4. Protocol (as `v5_plan.md` §4.4–§4.5 and §8.3)

### 4.1 Budgets

| Arm | Decisions per seed | Wall-clock safety cap | Environments × steps per update |
|---|---|---|---|
| `ppo_hist` | 2·10⁸ (preregistered) | 3 h on the L4 (preregistered; does not bind at the measured 48,500 decisions/s on a T4) | 4,096 × 64 |
| `liu_track` | 2·10⁸ (as `ppo_hist`) | as `ppo_hist` | 4,096 × 64 |
| `resid_fnpeg` | 4·10⁶ (see below) | 1.5 × the projected time on its routed runtime | 2,048 × 16, FNPEG labels on every core |

**change (budget rule):** a run's budget is its decision count; the wall-clock cap is a safety stop that is
max(3 h, 1.5 × the projected time at the rate measured on the routed runtime), so it does not bind on a healthy run
and the number of decisions does not depend on the runtime. `ppo_hist`'s preregistered "200 M decisions or 3
GPU-hours" is unchanged on the L4.

`resid_fnpeg`'s budget: every decision needs one FNPEG solution (a numerically predicted trajectory, CPU only), so
it is 50× smaller than the PPO arms'. Residual RL starts from a competent policy and needs far fewer samples than
PPO from scratch (Johannink et al.; Silver et al.). The measured label rate (§4.4) sets the cost. Measured on the
laptop (one process, after the guided phase starts): ~5, 34, 67 and 161 FNPEG labels/s at 32, 256, 1,024 and 2,048
rows per process (FNPEG's solver is vectorized over flights, so the rate grows with the batch). With 2,048
environments sharded over 8 processes (256 rows each) that is ~270 decisions/s, so ~4 CPU-hours per seed; the pilot
measures the real rate on Colab's CPU.

### 4.2 Tuning (DEV48 only, ≤ 12 configurations per arm, frozen before any test flight)

Random search, seed fixed, seed 0, 10 % of the budget per configuration, each configuration scored by 48
full-simulator DEV48 flights; objective as §4.5 (success@5km, then mean log miss, subject to mean peak g ≤ 4.6 g).
An arm whose first 3 configurations all score 0 % stops (`stopped.json`, as Part B).

| Arm | Space |
|---|---|
| `resid_fnpeg` | bound b {10, 20, 30}°; lr log-uniform 3·10⁻⁵–3·10⁻⁴; initial exploration std {0.1, 0.3}; critic burn-in {5, 20} updates |
| `liu_track` | lr log-uniform 10⁻⁴–10⁻³; k (δs = k s0) {0.2, 0.3, 0.4}; entropy {0, 10⁻³, 10⁻²}; clip {0.1, 0.2, 0.3} |
| `ppo_hist` | §4.5's, unchanged (lr, n_steps, entropy, GAE λ, clip) |

### 4.3 Seeds and flights

Seeds 0–4 at the full budget; R-each needs ≥ 4 of 5 seeds (§6). Each trained policy then flies (CPU, full
simulator), paired by condition with the `fnpeg` and `v2single` flights that already exist:

| Suite | Flights per arm | Seed rule | Paired with |
|---|---|---|---|
| A256 | 5 × 256 = 1,280 | every seed | `v2single` (5 seeds, v3/v4 folders), `fnpeg` |
| DS1 | 256 | seed = case mod 5 | `fnpeg`, `v2single` (notebook 1) |
| AS1 | 128 | seed = case mod 5 | `fnpeg`, `v2single` (notebook 1) |
| ORB-1 block 1 | 1,000 | seed = case mod 5 | `fnpeg`, `v2single` (notebook 1) |
| **total** | **2,664** (7,992 for three arms) | | |

Plus DEV48 tuning (48 per configuration) and gate RI's 96 flights.

### 4.4 Cost (projection; the status cell prints the measured one)

| Block | Runtime | Hours | Units |
|---|---|---|---|
| Gates: three 1 h pilots, FNPEG's 512 fixed episodes, RI's 96 DEV48 flights | High-RAM CPU | ~4 | ~1 |
| L4 timing runs (5 min × ≤ 3 arms) | L4 | ~0.3 | ~0.6 |
| `ppo_hist`, `liu_track`: 5 seeds × ~1.2 h + 12 × 0.12 h + DEV48 scoring, each | L4 | ~7.5 each | ~15 each |
| `resid_fnpeg`: 5 seeds + 12 configurations at the measured label rate | High-RAM CPU | (measured) | (measured) |
| Flights, 7,992 (`resid_fnpeg`'s cost FNPEG's 1 Hz solutions; the others seconds each) | High-RAM CPU | ~20 | ~5 |
| **total, if every arm passes its gate** | | | **~40** |

### 4.5 Statistics and the RQ3 family

* **change (RQ3 Holm family):** {`R-resid_fnpeg`, `R-liu_track`, `R-ppo_hist`} (was {`R-dagger_fnpeg`,
  `R-bcppo_fnpeg`, `R-ppo_hist`}). R-each as §8.3: success@5km `v2single − arm` > 0 pooled, 95 % CI excludes 0, Holm
  p < 0.05, holds in ≥ 4 of 5 seeds. **The family size stays 3**: an arm that failed its gate enters Holm with p = 1
  (conservative). Written as a `v5_plan.md` §15 entry before any of these arms runs.
* Estimates (no pass/fail): `arm − fnpeg` on A256 (success@5km, paired log miss); on DS1 and AS1 `v2single − arm`
  and `arm − fnpeg` (cases resampled within seed strata, as D/A comparisons); on ORB-1 block 1 `arm − fnpeg` paired
  by case and each arm's success@10km with its Clopper–Pearson bound; R-gap (fast vs full on the same A256 cases).
* Bootstraps, seed-level t-tests and "complete data before a verdict" as §6.

---

## 5. Drive

* **Output:** a new folder, `MyDrive/ReEntryAI_rl_v5` (smoke: `..._rl_v5_smoke`), the v5 layout (`00_study`,
  `01_checks`, `02_policies`, `03_tuning`, `04_flights`, `05_analysis`).
* **Inputs, read-only and guarded** (a write there raises): `MyDrive/ReEntryAI_jepa_v5` (FNPEG's `chosen.json` and
  `vf.json`, its flights and `flights.csv` for pairing), `MyDrive/ReEntryAI_jepa_v3` and `_v4` (`v2single`'s A256
  flights). Notebook 2 never writes into notebook 1's folder; the analysis joins both.
* The FNPEG configuration is copied with its SHA-256 into `00_study/inputs.json` at first use; a later change in
  notebook 1's `chosen.json` stops Setup.

---

## 6. Notebook

`StochasticEntrySim/ReEntryAI_RL_v5.ipynb`, generated by `jepa/build_notebook_v5_rl.py`; logic in
`jepa/rl2_v5.py` (arms, environments, flights, gates) and `jepa/study_rl_v5.py` (layout, stages, routing, stop
rules, analysis). Reuses `rl_v5.py` (Actor, PPO update, FastExpert/ParallelExpert's worker pattern, `_sync_flight`
for DAgger only, `fnpeg_fast_success`), `rlenv_fast_v5.py`, `fly_v5.py`'s runner and `runguard_v5.py`.

* **Runtime:** opens on a **CPU runtime** (the first pass is CPU only: gates and pilots, so no L4 time is spent
  before the gates pass). The text says when to switch to the L4 (`GPU_TYPE = "L4"`, locked in
  `00_study/gpu_platform.json`).
* **Passes:** CPU (gates, pilots, routing for CPU-routed arms, tuning and training of CPU-routed arms, flights,
  analysis) → L4 (timing runs, tuning and training of L4-routed arms; releases when its GPU work is done) → CPU
  (flights of the new policies, analysis). Each pass prints what the next one has to do.
* **CLAUDE.md:** `AUTO_DISCONNECT = IN_COLAB`, `release_runtime()` (Drive `flush_and_unmount()` first), the
  `post_run_cell` failure hook (ignores `KeyboardInterrupt`, unregisters an old copy), release in the last cell, and
  text saying it disconnects itself and that `AUTO_DISCONNECT = False` turns it off.
* **Watchdog and budget:** `STALL_MIN = 45`, `MAX_SESSION_UNITS`, `RUN_SUMMARY.txt`, as notebook 1.
* **Stop-and-disconnect rules** (`plan_jepa_needed_v5.md` §4b; stop on technical failure, never on a bad result):

| Trigger (id) | Action |
|---|---|
| no progress for 45 min (watchdog) | release (not a stop rule) |
| `rl_nonfinite`: a run's weights or return go non-finite, or its training check stops it at 50 % (no flight within 5 km and no rise in the return) | `failed.json`, **STOPPED**, release |
| `ri_failed`: gate RI fails | **STOPPED**, release |
| `flight_errors`: ≥ 5 of an arm's first 20 flights raise or end non-finite (not misses) | stop the arm, **STOPPED**, release |
| `failures`: ≥ 3 failure records of any kind in one session | **STOPPED**, release |
| `stage_exception`: a stage raises | **STOPPED**, release |
| the session reaches `MAX_SESSION_UNITS` | release (not a stop rule) |

  A stop writes `00_study/STOPPED.json` (trigger, evidence, time, commit) and `RUN_SUMMARY.txt`. The next Run all
  prints it at Setup and **refuses to continue** until `ACK_STOP = "<trigger id>"` is set in Settings; then the
  record moves to `00_study/stops/` and the notebook goes on (a stopped run stays skipped until its `failed.json` is
  deleted). Failed gates RP and poor accuracy are results: they never stop.
* **Notebook revision:** `NOTEBOOK_REV_RL`; an old copy stops at Setup.

---

## 7. Smoke checks (CPU, `python jepa/smoke_rl_v5.py`, numbers meaningless)

1. `resid_fnpeg` at init: actor mean exactly 0; the residual environment's command equals FNPEG's bit for bit over
   whole episodes; two DEV48 flights of `resid_fnpeg@init` and `fnpeg` are identical (telemetry and result).
2. The residual environment's FNPEG labels equal `rl_v5.FastExpert`'s (float64), in-process and in workers; its
   state round-trips (`get_state`/`set_state`) so a resumed run continues exactly.
3. `liu_track`: reference at the start equals v0 within e^(−1/k); sign changes only outside the deadband; reward
   terms finite; training terminations on, evaluation terminations off.
4. PPO on the new environments trains a tiny run, resumes from `train_state.pt` to the same weights as an
   uninterrupted run.
5. The 512 fixed episodes are the same entries for every arm.
6. Each gate's verdict on a forced pass and a forced fail; RP failure is recorded and does not stop.
7. Routing chooses the cheaper runtime and waits while a rate is missing.
8. Jobs: flight ids, seed rule (case mod 5), pairing keys with `fnpeg`/`v2single`; the jobs of a dropped arm are
   empty.
9. Stop rules: each trigger, forced, writes `STOPPED.json`; the next Setup refuses without `ACK_STOP`, continues with
   it.
10. Read-only: a write into the notebook 1 folder raises.
11. The notebook has the CLAUDE.md release code, opens on a CPU runtime, and an old `NOTEBOOK_REV_RL` stops.
12. The analysis builds the RQ3 table with Holm over 3 (a missing arm counted as p = 1) on smoke flights.

---

## 8. Rules for building

Do not change the simulator, the guidance laws, MPPI, v2/R0, or anything that alters a flown flight. New logic in
new modules; small additive hooks only in `study_v5.py`, `models_v5.py`, `analyze_v5.py` (owned by notebook 1's
session). Never push to `space-jepa4` or `jepa-v3`; push only when the user says so, with
`git push origin HEAD:space-jepa5`. Drive is read-only for the building session. Never run the full study.

---

## 9. Build record

Built 2026-11-19 on branch `claude/confident-bartik-cbee34` (not pushed).

**Files.**

| File | What |
|---|---|
| `jepa/rl2_v5.py` | `ResidExpert` (FNPEG on every fast-env row, float64, in-process or sharded over processes, state round-trips), `ResidEnv`, `LiuEnv` (+ `LiuCfg`, `liu_terms`, `liu_sign`), `make_actor`/`make_env`, `ppo_train_env` (rl_v5's PPO loop for any of the environments, snapshots carrying the arm), `train_dir`, `eval_fixed` (fixed episodes or given conditions), gates `ri_fast`, `rp_verdict`, `flights_identical`, `route`, `units_per_million`, full-simulator laws `ResidLaw`, `LiuLaw`, `ZeroResidLaw`, `make_law` |
| `jepa/study_rl_v5.py` | `RLStudyConfig` (+ `smoke()`), layout with notebook 1/v3/v4 read-only, `fnpeg_input` (SHA-256 copy), `freeze`, stop rules (`stop`, `check_stopped`, `check_failures`, `stage`), gates (`fnpeg_fixed`, `gate_ri`, `gate_rp`, `timing_run`), routing, `tune`, `train`, `fast_eval`, `flight_jobs`, `run_flights` (flight-error stop), `status`, `next_pass`, `analyze` (RQ3 with Holm m = 3, estimates vs fnpeg, DS1/AS1/ORB-1, R-gap, gates, tuning, cost; `05_analysis/RQ3_summary.md`) |
| `jepa/build_notebook_v5_rl.py` → `ReEntryAI_RL_v5.ipynb` | 21 cells; opens on a CPU runtime (no accelerator in the metadata); CLAUDE.md release code imported from `build_notebooks_v3.py` |
| `jepa/smoke_rl_v5.py` | the 12 checks of §7 plus a smoke execution of the notebook (twice: the second run only skips and reports) |

**Hooks in files owned by notebook 1's session.** One: `analyze_v5.load_flights_v5(..., extra_v5=())`, which
reads more v5 folders from their `flights.csv` and tags rows with `v5_folder`. With the default it behaves exactly
as before. `study_v5.py` and `models_v5.py` are not changed (only imported).

**Measured.** FNPEG labels on the laptop: ~5 / 34 / 67 / 161 per second per process at 32 / 256 / 1,024 / 2,048
rows (§4.1). In-process and 2-worker labels are identical over 25 decisions of 32 rows (410 s on the laptop).

**Differs from the text above.**

* `resid_fnpeg`'s fast-env budget is 2,048 environments × 16 steps (4·10⁶ decisions, lr 10⁻⁴, minibatch 4,096,
  critic burn-in 5 updates before tuning), not 512 × 32: larger batches raise FNPEG's label rate per process.
* The observation normalizer also updates during the critic burn-in in `ppo_train_env` (`rl_v5.ppo_train` freezes
  it then); `ppo_hist` still trains with `rl_v5.ppo_train`, unchanged.
* A tuning configuration that the training check stops for "no learning" scores 0 (a result); only non-finite
  weights or returns there are the stop rule `rl_nonfinite`. Training seeds stop on both, as §6 says.
* The pilots' own training check (rl_v5's at 50 %) makes a non-learning pilot a failed RP, not a stop.
* CPU-routed arms wait for the L4 timing run (routing needs both rates), so `resid_fnpeg`'s tuning and training
  happen on the third pass (CPU), after the L4 pass.
* Notebook 1's own flights of an arm named `resid_fnpeg`, `liu_track` or `ppo_hist` (if any) are dropped from the
  join: RQ3's arms are notebook 2's.
* No PDF report: the analysis writes CSV tables and a Markdown summary; notebook 1's report can read them later.
