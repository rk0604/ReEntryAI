# SpaceJEPA v3: full study plan (v1 vs v2 vs classical)

Branch: `jepa-v3` (from `jepa-v2`). Status: **implemented, ready to run on Colab.** The code is in
`jepa/` (`testset_v3`, `study_v3`, `fly_v3`, `analyze_v3`, `smoke_v3`, `build_notebooks_v3`), and the four
notebooks are `ReEntryAI_JEPA_v3_{0_testset,1_train,2_fly,3_analyze}.ipynb`. §13 lists where the code
had to make a choice the plan left open. This file is the reference for building and running v3. The v1/v2 design is in `PLAN.md` §1–11; v2's first results are in
`SpaceJEPA_v2_results.pdf`.

---

## 0. What v3 is

v3 is a **study**, not a new model. It answers, with enough data to trust the answer:

> On a much more varied set of entry conditions, how do the v1 and v2 world models (each trained 5
> times on the same 6,000 flights, v1 with a fair training budget) compare with the classical
> predictor-corrector, at matched load, against a **5 km** landing goal? And what do their RCS usage,
> heating and g-load look like?

Decisions already made (from the planning discussion):

| # | Decision |
|---|---|
| 1 | v3 is a full study of v1 and v2; no new architecture |
| 2 | Only the **entry** varies; the **target stays fixed** at 15°N, 15°E (from `configs/default.json`) |
| 3 | The same 4 atmospheres as v2 |
| 4 | Report every metric (miss within 5 km, median, p90, g, heat, fuel), each defined below |
| 5 | Colab **L4 or T4** GPU |
| 6 | RCS is **logged only**; varying RCS strength is deferred to v4 |
| 7 | **No sensor noise** (perfect navigation and IMU, as in v1/v2) |
| 8 | v1 gets a **fair training budget** on the **same 6,000 flights** as v2 |
| 9 | Sampling is done properly everywhere: test-set design, training seeds, MPPI noise |
| 10 | **5 seeds**, **64 cases**, success = **miss < 5 km** |
| 11 | Everything (checkpoints, logs, flights, telemetry, figures, tables, report) goes to **Google Drive**, neatly organized and labeled |
| 12 | No simulator changes (`rl_env.py`, physics, RCS, guidance stay untouched) |

---

## 1. Test set: 64 entry cases

### 1.1 Composition
* **16 bridge cases:** the same first 16 frozen ICs used in v1 and v2 (`configs/heldout_testset.csv`,
  rows 0–15). They connect v3 to the earlier results. Their classical runs are flown again so they
  have telemetry (§13).
* **48 new cases:** entry altitude, speed, flight-path angle, heading and entry position varied far
  more widely (below). They are all feasible by construction.

### 1.2 How the 48 new cases are sampled
With the target fixed, an entry can only be varied freely if its **position** moves with it, so the
target stays reachable. Each new case is therefore built backwards from the target:

1. **Draw the entry state** from a scrambled Sobol sequence (fixed seed `20261018`, so the set is
   reproducible), over 6 dimensions:

   | Dimension | Range | v1/v2 test set had |
   |---|---|---|
   | entry altitude | 105 – 135 km | 119 – 121 km |
   | entry speed | 7.2 – 8.4 km/s | 7.70 – 7.89 km/s |
   | flight-path angle | −1.5° to −4.0°, with v2's skip-out guard (steeper required at high speed) | −2.2° to −2.9° |
   | approach heading (direction the capsule flies toward the target) | 0 – 360° | 43.6° – 47.7° |
   | downrange position inside the reachable band | 10 – 90% of the band | fixed ~2,300 km |
   | crossrange offset inside the reachable band | −70% to +70% of the half-width | ~0 |

2. **Compute the reachable footprint** of that entry state in **each of the 4 atmospheres**. Fly
   constant-bank trajectories (σ = 0°, ±15°, ±30°, …, ±90°) with the simulator's own equations of motion
   (`math_3d.nominal_eom_step`) and the true perturbed parameters, down to drogue altitude. The landing
   points span the footprint; their downrange interval and crossrange half-width define the reachable band.

3. **Intersect the 4 footprints** and keep only the part reachable in all four atmospheres, with
   **≥ 50 km margin** (to absorb actuator lag, which the sweep ignores).

4. **Place the entry:** choose downrange D and crossrange c from the Sobol fractions inside that
   intersection. Then compute the entry latitude, longitude and heading by spherical geometry,
   backwards from the target along the approach heading.

5. **Reject** a candidate if:
   * the entry latitude is outside ±30° (the training data's range: latitude is a model input);
   * the lift-up (σ = 0°) trajectory exceeds **4.0 g** in any atmosphere (so the 4.5 g load limit in
     §4.3 is always feasible);
   * the intersection is empty.

   Keep drawing Sobol points until 48 cases are accepted. Log the number of rejections and the reasons.

6. **Never filter by controller performance.** Feasibility comes from physics alone. No controller's
   result is used to choose cases.

7. **Freeze** the result as `configs/v3_testset.csv`, committed to git and copied to Drive. Each row
   records the IC, the downrange/crossrange to the target, the footprint margins in each atmosphere,
   the lift-up peak g, and the Sobol index.

### 1.3 Bridge-case feasibility
The bridge cases are not filtered. The same footprint computation flags them, and v2 already showed
cases 5 and 14 cannot reach the target in the ×0.78 atmosphere. The primary results use all 64 cases.
A secondary table excludes condition pairs flagged infeasible.

### 1.4 Atmospheres (unchanged from v2)
| Tag | Density | L/D | Mass |
|---|---|---|---|
| `nominal` | ×1.00 | ×1.00 | ×1.00 |
| `mismatch` | ×0.88 | ×0.88 | ×1.00 |
| `mismatch_high` | ×1.12 | ×1.12 | ×1.00 |
| `out_of_range` | ×0.78 | ×0.80 | ×1.00 |

64 cases × 4 atmospheres = **256 conditions** per controller and seed.

---

## 2. Models and seeds

### 2.1 Data (shared by every model and seed)
* The same 6,000 wide-envelope flights (`MyDrive/ReEntryAI_jepa_data/wide_6000/`).
* **One fixed split** (episode split seed 0, 80/10/10), identical for every model and seed. The
  normalization statistics are fit on the training split only. Validation is used for monitoring;
  the test split is used only for the open-loop evaluation.

### 2.2 v2 (as in the first full run)
5-member ensemble, 12.3M parameters, the same architecture and recipe. Stage 1: 24,000 steps, batch
256 per member, horizon curriculum 24 → 48 → 96. Stage 2: 6,000 steps, batch 512, curriculum 16 → 32 → 48.
**SIGReg λ is fixed at 0.08**, the value v2's validation sweep chose. There is no re-sweep, so λ is
never tuned on the test set and the cost is saved.

### 2.3 v1 with a fair budget
v1's own architecture and recipe (26k-parameter TCN + GRU, T = 24, pooled SIGReg λ = 0.02, PI prober
with T = 16, lr 3e-3, clip 0.5) on the **same 6,000 flights and split**, with the **same training
budget as one v2 member**:
* stage 1: 24,000 optimizer steps × 256 windows (v1-wide in v2 got 6,000);
* stage 2: 6,000 steps × 512 windows;
* trained on the GPU, with the same data-on-GPU pipeline as v2 (only the device changes, not the recipe).
* **Convergence check:** if v1's validation loss improves by more than 2% over the last 20% of
  training, extend it once by 50%, and record that.

v1 is a single model per seed (no ensemble), as it was designed.

### 2.4 Secondary arm: v2 single member
Member 0 of each seed's ensemble, flown alone (no averaging, no disagreement term). It separates "a better
model" from "5 models". It is cheap: GPU rollouts, already trained.

### 2.5 What a "seed" controls
Seeds are `s ∈ {0, 1, 2, 3, 4}`. Seed s sets:

| Randomness | Seeded by | Notes |
|---|---|---|
| weight initialization | `s` | v2: all 5 members drawn from seed s (25 members in total) |
| training data order and window sampling | `s` | |
| SIGReg random projections | `s` | |
| MPPI sampling noise | `hash(s, case, atmosphere)` | see §4.2 |
| data split, normalization, test set | **fixed** | identical across seeds |

GPU training is not bit-for-bit deterministic. Seeds reproduce statistically, not exactly.

---

## 3. Controllers flown

| Controller | Seeds | Load limit | Conditions | Flights |
|---|---|---|---|---|
| classical predictor-corrector | 1 (deterministic) | its own guidance | 256 | 256 (all re-flown for telemetry) |
| **v1** (fair budget) | 5 | **4.5 g** (primary) | 256 | 1,280 |
| **v2** (ensemble) | 5 | **4.5 g** (primary) | 256 | 1,280 |
| v2 single member | 5 | 4.5 g | 256 | 1,280 |
| v1, v2 legacy | seed 0 | 8 g (as in v1/v2) | 256 each | 512 |
| v2 planner-noise repeat | seed 0, 2 extra MPPI seeds | 4.5 g | 256 × 2 | 512 |
| *optional:* classical with true atmosphere (oracle upper bound) | 1 | — | 192 (3 mismatched atmospheres) | 192 |

* **Primary comparison: v1 vs v2 vs classical at the 4.5 g load limit.** v2's first run pulled about
  5.8 g against the classical controller's 3.6 g. The 4.5 g arm removes that advantage.
* **The legacy 8 g arm (seed 0)** connects to the v2 results and measures what load matching costs.
* **The planner-noise repeat** separates MPPI's own randomness from training-seed variation.

---

## 4. Closed-loop protocol

### 4.1 One runner for all learned controllers
v1 and v2 fly through **the same lockstep MPPI runner** (`mppi_v2.run_lockstep`, extended with a
v1 model adapter). Every learned controller therefore shares identical env setup, MPPI code, cost terms
and timing. The classical controller runs through the unchanged `mppi.run_episode` (baseline mode).
* Equivalence test before the study: for v1 and for single-member v2, the lockstep runner must match
  the unchanged `mppi.run_episode` decisions to < 0.5° over the first 20 decisions on 3 cases. It
  is a CPU smoke test, like `smoke_v2.py`.

### 4.2 MPPI settings and sampling
* 256 samples, 6 noise knots, 30° noise, adaptive temperature, 4 iterations at entry interface and
  1 afterwards. These are v1/v2's settings, unchanged.
* **Common random numbers.** For a given (seed, case, atmosphere), every learned controller uses the
  **same** MPPI noise stream, `hash(s, case, atmosphere)`. Paired comparisons then differ only in
  the model, not in sampling luck. Different seeds and conditions get independent streams.
* **Planning horizon** is set once, from the test set: H = ceil((longest lift-up flight time over
  the 64 cases and 4 atmospheres + 60 s) / 5 s). It is the same for every learned controller and is
  recorded in the study config. If H exceeds v2's 96-step training horizon, the extrapolation
  factor is reported.
* Cost terms are unchanged: predicted miss, smoothness, g over the limit, γ steeper than 82°, and the
  never-reaches-drogue penalty. v2 adds the ensemble disagreement at 1 km/km. The g weight stays at
  50 km per g²; only the limit changes (4.5 g primary, 8 g legacy).

### 4.3 Job queue and ordering
Every flight is one job keyed by `(controller, load_limit, seed, mppi_seed, case, atmosphere)`. Its
result is written to Drive as soon as it finishes, so a rerun skips completed jobs. Jobs run **seed
by seed**: all controllers for seed 0 first, then seed 1, and so on. If the study stops early, the
seeds it did finish are complete, paired comparisons.

---

## 5. Metrics (exact definitions)

Every flight records these. Unless stated, they come from the simulator's own `info` at termination.

| Metric | Definition |
|---|---|
| **miss distance** [km] | Great-circle distance, on a sphere of radius 6,371 km, between the capsule's latitude/longitude when its altitude first drops to 7.6 km (drogue deploy) and the target (15°N, 15°E). If the flight ends any other way (`cos_gamma_fail`, skip-out, timeout), miss is the distance at termination and the flight counts as a failure. |
| **success @ 5 km** (primary) | outcome = `drogue` **and** miss < 5 km |
| success @ 25 km (legacy) | outcome = `drogue` and miss < 25 km, as in v1/v2 |
| **median miss** [km] | median of miss distance over the flights in a group |
| **p90 miss** [km] | 90th percentile of miss distance over the group (linear interpolation) |
| **peak g** [g] | maximum over the flight of the aerodynamic load factor, total aerodynamic force ÷ (mass × 9.80665), sampled at every physics substep (`info["peak_g"]`) |
| **peak heat rate** [MW/m²] | maximum stagnation-point heat rate from the simulator's heating model, sampled at every substep (`info["peak_qdot_MWm2"]`) |
| **heat load** [MJ/m²] | the integral of the stagnation-point heat rate over the flight, computed with the trapezoid rule from the 1 Hz telemetry (`env._last_qdot`) |
| **RCS fuel** [kg] | propellant used by the 12-jet RCS over the flight, accumulated by the simulator from jet on-time (`info["rcs_fuel_kg"]`) |
| flight time [s] | time from entry interface to termination |
| bank reversals | sign changes of the bank command, counting only commands with \|σ\| > 5° |
| planning time [ms] | per 5 s decision (amortized for batched runs) |

**How groups are summarized.** For each (controller, atmosphere) and each seed, compute the metric over
the 64 cases. Report the **mean over the 5 seeds** and the **range across seeds**, and also the
value pooled over all seeds. The classical controller is deterministic, so it has a single value.

---

## 6. RCS logging (log only; v4 will vary RCS strength)

Each flight saves a 1 Hz telemetry file (`.npz`, compressed) with:
* time, altitude, speed, flight-path angle, latitude, longitude;
* commanded bank, actual bank, roll rate;
* RCS fuel increment and cumulative fuel, and **jet-on seconds** (fuel increment ÷ the simulator's
  propellant rate per jet-second, exact because the simulator's fuel model is linear in on-time);
* g-load and heat rate.

All of it is read from the env after each step, with no simulator change.

Per-flight RCS summaries go into the flight table: total jet-on time, fuel by flight phase (above 60 km,
60–25 km, below 25 km), max roll rate, RMS bank-tracking error (commanded − actual), and the
fraction of time the roll rate is within 10% of its limit.

**RCS analyses (descriptive, no pass/fail):**
* fuel vs. miss, per controller;
* bank-tracking error vs. miss (does a smoother, cheaper plan lag more?);
* fuel by phase per controller;
* whether learned controllers' fuel savings change at the 4.5 g limit.

---

## 7. Open-loop evaluation (per seed)

On the held-out wide test split, for each seed: landing-point error vs. time-to-go (400, 300, 200,
100, 50 s), split by density bin, and the signed along-track bias by density. This time the bias
data is saved to CSV so it can be recomputed. Ensemble spread calibration is computed for v2.
Report per-seed values and mean ± range. This shows whether v2's 10× open-loop advantage is stable
across seeds.

---

## 8. Statistics and pre-registered criteria (fixed before any v3 flight)

### 8.1 Method
* **Unit:** one (case, atmosphere) condition, paired across controllers (common ICs and common MPPI
  noise). Seeds are repeated measures.
* **Uncertainty:** a hierarchical bootstrap with 10,000 resamples. Resample the 5 seeds, then the 64
  cases within each atmosphere; this gives 95% CIs for differences in success rate, median and p90.
* **Per-seed check:** each headline claim is also tested on each seed separately (paired sign test
  or McNemar).
* **Multiple comparisons:** Holm correction over the family of criteria below.

### 8.2 Criteria (load-matched 4.5 g arm unless stated)

| # | Question | Passes if |
|---|---|---|
| C1 | Is v2 better than v1 at an equal budget? | pooled success@5km difference v2 − v1 > 0 with CI excluding 0, holding in ≥ 4 of 5 seeds |
| C2 | Does v2 beat the classical controller? | pooled success@5km difference v2 − classical > 0 with CI excluding 0, holding in ≥ 4 of 5 seeds, and no atmosphere where the classical controller is significantly better |
| C3 | Is the comparison really load-matched? (validity check for C2) | v2's mean peak g ≤ classical's + 0.5 g in every atmosphere |
| C4 | Does v2 reach the goal? | v2 success@5km ≥ 90% in every atmosphere (on feasible conditions) |
| C5 | Is v2 no worse on heating? | median paired ratio v2 ÷ classical ≤ 1.10 for peak heat rate and for heat load |
| C6 | Is v2 stable across seeds? | success@5km range across the 5 seeds ≤ 10 percentage points in every atmosphere |
| C7 | Does the ensemble matter? (secondary) | v2 success@5km − v2-single > 0 with CI excluding 0 |
| Bridge | Does v3 reproduce v2's run? | v2, seed 0, 8 g, on the 16 bridge cases: median miss per atmosphere within 1 km of the v2 results |

Also reported (no pass/fail): median and p90 miss, success@25km, fuel, the RCS analyses, what load
matching costs (8 g vs 4.5 g, seed 0), and the planner-noise share of variance.

---

## 9. Google Drive layout

```
MyDrive/ReEntryAI_jepa_v3/
  00_study/
    study_config.json        frozen settings: git commit, package versions, all numbers in this plan
    manifest.csv             every artifact: path, stage, created, size, sha1
    README.md                auto-written index of what is where
  01_testset/
    v3_testset.csv           the 64 cases + footprint margins (also committed to git)
    footprints.csv           per case × atmosphere: reachable band, margins, lift-up peak g
    F01_testset_coverage.png, F02_footprint_margins.png
  02_models/
    v1/seed_0 … seed_4/      ckpt.pt, train_state.pt (resume), train_log.csv, F_curves.png
    v2/seed_0 … seed_4/      same
  03_open_loop/
    footprint_<model>_seed<k>.csv, along_track_<model>_seed<k>.csv, spread_v2_seed<k>.csv
  04_flights/
    results/<controller>/<load>/seed<k>/mppi<m>/<atmosphere>_case<NN>.json   one file per flight
    telemetry/<same path>.npz                                                 1 Hz telemetry incl. RCS
    flights.csv              merged table, one row per flight (rebuilt from results/)
  05_analysis/
    tables/T01_…csv …        every table in the report, numbered
    figures/F10_…png …       every figure, numbered, plus captions.json (number → caption, data source)
  06_report/
    SpaceJEPA_v3_report.pdf  auto-generated from 05_analysis
```

Naming rules: figures are `F<NN>_<what>_<by-what>.png` (for example
`F14_success5km_by_atmosphere.png`), tables `T<NN>_<what>.csv`. Every figure has a caption
entry saying what it shows and which file its data came from. Writes are atomic (temp file + rename),
so a disconnect never leaves a half-written file.

---

## 10. Notebooks (stages that fit in Colab sessions)

| Notebook | Runtime | Does | Est. time (L4) | Est. time (T4) |
|---|---|---|---|---|
| `v3_0_testset` | any (CPU) | build and freeze the 64 cases, footprints, coverage figures | ~10 min | ~10 min |
| `v3_1_train` | GPU | train v2 and v1 for seeds 0–4 (resumable, one seed at a time) | ~4.5 h | ~8.5 h |
| `v3_2_fly` | GPU + CPUs | the job queue: classical (CPU pool), v1/v2 lockstep (GPU), telemetry | ~4–8 h* | ~7–9 h* |
| `v3_3_analyze` | GPU optional | open-loop per seed, statistics, figures, tables, the PDF report | ~30 min | ~45 min |
| `v3_resultspdf` | CPU | everything the four notebooks saved, in one PDF (`06_report/SpaceJEPA_v3_full_results.pdf`, `jepa/results_pdf_v3.py`); reads only, recomputes nothing | ~5–10 min | ~5–10 min |

\*The flight stage is limited by CPU cores, not the GPU:
* classical: ~6.4 CPU-hours (+6.4 h if the optional oracle runs);
* learned controllers: ~8 CPU-hours of simulator stepping.

Standard Colab GPU runtimes often have only 2 vCPUs, and the notebook prints the measured count and
an ETA. Every stage is resumable, so the study can span several sessions. The classical stage can also
run on any machine with more cores (for example a local Jupyter on the 12-core laptop) and upload its
results to the same Drive folder.

---

## 11. Build order (after this plan is approved)

1. **Code** (all inside `jepa/`, nothing else):
   - `testset_v3.py` (footprints, Sobol sampling, freeze);
   - a v1 model adapter for the lockstep runner;
   - common-random-number MPPI seeding and telemetry capture in `mppi_v2`;
   - `study_v3.py` (config, Drive layout, manifest, job queue);
   - `analyze_v3.py` (statistics, figures, tables, report);
   - notebook builders for the 4 notebooks.
2. **CPU smoke test** (`smoke_v3.py`):
   - equivalence of the lockstep runner with v1 and v2;
   - the CRN check (same noise across controllers);
   - telemetry fields and the fuel → jet-seconds identity;
   - resume and skip of the job queue;
   - a tiny end-to-end run of all 4 notebooks in smoke mode.
3. Commit, push to `jepa-v3`, and run on Colab.

---

## 12. Out of scope (deferred)

* **v4:** vary RCS strength (needs simulator parameter changes), sensor and navigation noise, moving
  targets, atmospheres that decouple density from L/D or vary with altitude, mass dispersion in the tests.
* No new world-model architecture in v3. Model changes suggested by v2 (on-policy data, terminal-phase
  replanning, a redesigned uncertainty term) wait until v3 has measured the current models properly.

---

## 13. Implementation notes (decisions made while building)

* **Footprints respect the load limit.** A constant-bank sweep from ±0° to ±90° reaches short ranges only
  by pulling up to 10 g. A target there would be infeasible at the 4.5 g limit (§3), which would break the
  intent of §1.2 step 5. So the reachable region uses only the banks whose peak load stays ≤ 4.5 g:
  * The sweep runs in 5° steps outward from lift-up, stopping at the first bank over the limit.
  * Mixing banks (reversals, late lift-up) reaches at least this region, so it is a conservative
    inner estimate.
  * The sweep was checked against real simulator flights and lands within about 4 km of them
    (`smoke_v3.py` check 1).
* **The frozen test set** (`configs/v3_testset.csv`, built by `python jepa/testset_v3.py`):
  * 148 Sobol candidates gave 48 accepted cases. Rejections: 84 for lift-up peak g > 4.0, 15 for entry
    latitude outside ±30°, and 1 for margin.
  * Coverage: altitude 105–135 km, speed 7.43–8.40 km/s, γ −1.69° to −3.97°, approach heading 8°–355°,
    downrange 1,980–3,100 km.
  * Bridge cases flagged infeasible at 4.5 g: nominal {5}, ×0.88 {4, 5}, ×0.78 {3, 4, 5, 6, 13, 14}.
    They stay in the primary results, and T02 leaves them out.
* **Planning horizon** = 188 steps (940 s), from the longest lift-up flight (876 s in the ×1.12
  atmosphere) + 60 s. That is 1.96× v2's longest training horizon (96 steps) and 7.8× v1's (24). The
  extrapolation is the same for every seed, and is recorded in the study config.
* **Classical:** all 256 conditions are flown again instead of reusing the 64 cached bridge flights, for
  two reasons:
  * The cached results have no telemetry, and §6 needs RCS data for every flight.
  * The cache (`jepa_runs/classical_cache.pkl`, which v2's notebook also reused for nominal and ×0.88)
    comes from an earlier simulator version. Example: bridge case 0 at ×0.88 is 3.43 km cached vs 3.48 km
    now, with flight time 559.25 vs 559.0 s. `fly_v3.fly_classical` reproduces the current
    `mppi.run_episode` exactly.

  The cost is about 1.6 extra CPU-hours.
* **v1 convergence extension (§2.3):**
  * The rule compares the validation loss at the last log at or before 80% of the steps with the final
    one.
  * An extension trains 50% more steps on a fresh cosine cycle, warm-restarted at 30% of the peak
    learning rate.
  * It applies to both v1 stages and is recorded in the checkpoint.
* **Holm (§8.1)** adjusts the bootstrap p-values of the hypothesis-test criteria (C1, C2, C7). Each of them
  then also needs a Holm-adjusted p < 0.05. C3–C6 and Bridge are threshold checks.
* **Classical planning time** is the wall time of the whole 1 s guidance step, which includes the
  physics. It is not comparable with MPPI's amortized per-decision time.
* **Lockstep equivalence (§4.1)**, `smoke_v3.py` check 2:
  * **v1:** the runner reproduces `mppi.run_episode` to within 1e-5° per decision and to the same miss
    and fuel.
  * **Single-member v2:** the runner is identical to v2's own `mppi_v2.run_lockstep` over whole flights.
    Its parallel predictor rollout matches the KV-cached step-by-step rollout that `run_episode` uses to
    within about a metre in landing point.
  * **Why v2 can't match `run_episode` for 20 decisions:** the closed loop amplifies that float32
    round-off. The decisions agree to about 1e-3° for the first several steps and then drift apart (up
    to tens of degrees by decision 20 with a smoke-trained model). So the plan's "< 0.5° over 20
    decisions" is applied to v1, and v2 is checked by the three deterministic tests above.
* **GPU memory:** the first L4 run ran out of memory. With 4,096 samples per rollout and 2,048 rows per
  predictor pass, each worker needed 6.7 GB, against the 5 GB I had estimated, and 4 workers did not fit
  in 22 GB. The fix:
  * The defaults are now 2,048 samples per rollout and 512 rows per pass, whose measured peak is ≤ 5 GB.
  * The notebook budgets 6 GB per worker with 2 GB spare (L4 → 3 workers, T4 → 2) and sets
    `expandable_segments`.
  * These batching knobs change speed and memory, never results, so the frozen-config check ignores them.
