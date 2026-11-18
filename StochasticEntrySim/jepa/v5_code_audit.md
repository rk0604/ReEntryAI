# SpaceJEPA v5: pre-run code audit

2026-11-17. Audited `space-jepa5` at `20daea2` (the frozen study folder: 2026-11-16 at `f4503a7`, amended
2026-11-17), in the worktree `audit/v5-code-audit`. **No code was changed.** Method:

* I read the code paths that `ReEntryAI_JEPA_v5.ipynb` and `ReEntryAI_JEPA_v5_kaggle_rl.ipynb` execute.
* I ran targeted CPU experiments on two platforms: Windows 11 (Python 3.13.14, torch 2.14.0+cpu, numpy 2.3.5), and
  Docker `python:3.13-slim` on Linux (Python 3.13.18, torch 2.15.0+cpu from the CPU index, numpy 2.5.5), which stands
  in for Colab. The two platforms differ in library versions as well as OS, so a cross-platform difference below
  could come from either.
* I ran synthetic statistics tests where the true answer is known.

The scripts are in the appendix. Every number below is from those runs, not from the study's own Drive folder,
which I had no access to.

Impact code for every finding: (a) the fix changes a frozen setting in `00_study/study_config.json`; (b) the fix
invalidates flights already flown (S0, the 2,560 legacy-RL flights, gates V1–V4, E0); (c) neither.

## 0. Resolution (2026-11-17, branch `v5-fixes`, user decisions)

Every fix below is impact (c): no frozen setting changed and no flown flight is invalidated. Each one is logged in
`v5_plan.md` §15 under 2026-11-17. The GPU type is now an **L4** for all GPU work, RL included; Kaggle is dropped.

| # | Decision | Done | Checked by |
|---|---|---|---|
| A1 | two-stage bootstrap + seed-level t-test | `orb_paired` draws the seeds, then the cases within each drawn seed; it returns `seed_t` (t with 4 df). The cache is keyed on the scheme. The dead code after `as1_analysis` is removed | null simulation `orb_test2.py`: 1.3 / 3.3 / 9.7 / 14.3% on 1,000 cases at a per-seed SD of 0 / 2.8 / 5.7 / 11.3 pp (was 5.2 / 15.5 / 35.8 / 58.5%), and 0.8–3.0% on 64 cases. The t-test gives 3.7–5.3% |
| A2 | tuned β from FNPEG's `chosen.json` | `mppi_filtered` carries `fnpeg_filter_beta` and waits for FNPEG's tuning. The adapter cache is keyed on β. `FastExpert` keeps the tuned β^5 | `a2_check.py`: 8/8 checks (untuned wait, β from `chosen.json`, config 0 gets the default, flight id unchanged, cache, expert) |
| B1 | t-interval beside every verdict, flag disagreements | the Comparer's pooled rows carry `seed_t`. T05 gains `t_lo`, `t_hi`, `t_p`, `t_n`, `t_excludes_0` and `boot_t_disagree` (J-equiv uses the 90% t-interval against ±0.15). The report marks "(seed-level t disagrees)". **Verdict rules are unchanged** | smoke checks 9 and 15 |
| B2 | L4 | `GPU_TYPE = "L4"`: R0 (v4, L4) and the E3 arms share a GPU type and numeric path (native bf16, TF32) | smoke check 1 |
| B3 | extend §16 | §16 "Single landings depend on rounding" names FNPEG, learned-model MPPI, the batch dependence, and that resume is bit-exact for the guidance laws and classical only (an RL policy's commands differ by 4e-7° with the batch, landing unchanged) | measured, see the §16 text |
| C1 | plan to 95% | §8.4 and §14: Q2 is two-sided 95%, Q1 one-sided 90% | text |
| C2 | no change | | |
| C3 | test on the tuning's ground truth (all 48 DEV48 conditions) | confirmed on all 48: Q7F, Q19, KA and K_ρ's gain, all four at the low or all at the high end, leave every flight bit-identical except the logged density estimate. Documented in §14 (beside the tuning ranges), §15, §16 and T10's caption. `PREDGUID_SPACE` is unchanged (tuning is running) | `pg_inert48.py`, `pg_inert48b.py`, `pg_telkeys.py`: LATBIAS high (control) changes 39/48 |
| C4 | key on the whole config | `design_for` keys on `asdict(cfg)` | 36 PredGuid flights (6 DEV48 conditions × 6 configurations in one process): bit-identical before and after |
| C5 | moot (Kaggle dropped) | `kaggle_rl.md` marked retired | |
| C6 | assert unique legacy names | `legacy_jobs` raises on a collision | the study's 5 legacy models have distinct names (`status.json`: 256/256 per suite each) |
| C7 | document | §16 "Heating is evaluated at the nominal density": since v1, and in every suite. Simulator unchanged | `math_3d.py:637-668`, `data.py:134` |
| C8 | snapshot by decisions or wall clock | `_maybe_snapshot` fires at fraction (i+1)/5 of either budget. Points not reached when training ends are the final policy. The trigger is recorded | unit check |
| C9 | fix the docstring | corrected in a comment at its use (`models_v5.py`, E3 training), because smoke check 11 requires `models_v2.py` to stay byte-identical to the v4 commit | smoke check 11 |

Notebook: no Kaggle step, an L4 by default, and every CPU stage skipped on a GPU runtime. A GPU pass refuses to
start if a gate is not done, and section 15b releases the runtime after the last GPU stage. RL tuning's DEV48
scoring flights are the one exception: they fly in the GPU pass. The RL resume lock covers the GPU type and
platform only, and torch/CUDA changes are logged. Smoke check 1 checks the cell order.

---

## 1. Summary verdict

**The expensive runs can proceed.** No finding invalidates a flight already flown, and no proposed fix needs a
frozen setting changed. Two findings need a decision, each at a specific point:

| When | Decision |
|---|---|
| **Before the first `mppi_filtered` flight (T4 pass) and before Kaggle RL tuning starts** | **A2**: does `mppi_filtered` (G5) and the DAgger/BC expert use FNPEG's *tuned* filter gain β, as the plan says, or the fixed 0.9 that the code uses? Fixing it later means re-flying those arms. |
| **Before any ORB-1, DS1 or AS1 result is read** | **A1**: the ORB-1 Q3, DS1 D1/D2 and AS1 A1/A2 tests treat the five training seeds as fixed. In a synthetic test that makes them reject a true null 15–59% of the time on ORB-1 when the v2 seeds differ by 3–11 pp. The fix is analysis-only; flights are unaffected. |

The rest weakens claims or is documentation (B1–B3, C1–C9). Everything listed in §3 checked out. The biggest gaps
are GPU behaviour and anything that needs the Drive data (§4).

| # | Rank | Finding | Fix impact |
|---|---|---|---|
| A1 | A | ORB-1/DS1/AS1 confirmatory tests ignore seed-to-seed variance (false-positive rate up to 59%) | (c) |
| A2 | A | `mppi_filtered` and the DAgger/BC expert hard-code β = 0.9; the plan says FNPEG's tuned β | (c) if decided before those arms run |
| B1 | B | A256 hierarchical bootstrap with 5 seeds under-covers when seeds differ; the 4-of-5 rule does not fix it; J-equiv is affected too | (c) |
| B2 | B | RQ2's key arm crosses GPU platforms (R0: v4 on L4; E3: v5 on T4), and on a T4 `amp_dtype` gives emulated bf16, not the fp16 the code and plan say | (c) |
| B3 | B | Rounding and batch sensitivity goes beyond PredGuid: FNPEG up to 0.43 km on one platform (0.72 km across two), learned-model MPPI 0.56 km from 1e-8 of entry speed and 0.55 km from batch composition. Consequences for the tuning tie-break, resume equivalence and §16 | (c) |
| C1–C9 | C | Documentation and defensive checks (§2.3) | (c) |

---

## 2. Findings

### A1. ORB-1 Q3, DS1 D1/D2 and AS1 A1/A2 treat the five training seeds as fixed

**Where.** `jepa/analyze_v5.py:324-357` (`orb_paired`). The callers are `:1093` (ORB-1 Q3, strata
`case % 5`) and `:390` (DS1 and AS1 via `shift_analysis`, strata `case % len(cfg.seeds)`). A third call at `:450`
sits after `as1_analysis`'s `return` (a leftover copy of the DS1 loop) and never runs.

**What.** On these suites each case is flown by one learned-arm seed (`case mod 5`). `orb_paired` resamples
cases *within* each seed stratum and never resamples the seeds. Its interval therefore answers "do these five v2
checkpoints beat FNPEG on these cases?" It does not answer the question RQ1's A256 analysis asks of the same
comparison: "does the v2 recipe beat FNPEG?" `analyze_v4.Boot` resamples seeds for that (`analyze_v4.py:329`), and
v5_plan §3.2 (`v5_plan.md:188`) says the seed-per-case design makes ORB-1 a test "for the model family (a randomly
chosen training run)". Under that reading, the variance between seeds is missing from three confirmatory Holm
families. They would reach "met" (or "significantly worse") far too easily.

**Evidence.** Synthetic null over seeds (`orb_test.py`, appendix A.3). Five seeds, each with a success-rate offset
drawn around a population mean equal to the comparator's, so the true recipe-level difference is 0. The two arms
share a per-case difficulty. 400 replications, `n_boot` 1,000, two-sided p < 0.05.

```
python orb_test.py <sim> 1000          python orb_test.py <sim> 64
SD of per-seed success   orb_paired rejects   seed-level t-test   |   orb_paired (n=64)   t-test
0 pp                      5.2%                 4.5%               |    4.5%               4.2%
2.8 pp                   15.5%                 5.5%               |    5.0%               3.0%
5.7 pp                   35.8%                 5.0%               |    6.8%               3.8%
11.3 pp                  58.5%                 4.5%               |   12.0%               3.5%
```

`orb_paired` is calibrated when the seeds are identical. The inflation grows with the number of cases: ORB-1 block 1
has 1,000, DS1/AS1 have 64 cases with four or two atmospheres averaged. A seed-level test (the five per-seed
differences, t with 4 df) stays at about 5%.

**Consequence.** RQ4 Q3a/Q3b, DS1 D1/D2, AS1 A1/A2, and every estimate row that uses `orb_paired`: the v2-vs-RL
rows (D-*, A-*) and the Q3 sensitivities in T08b.

**How big in practice** depends on how much v3's five `v2single` checkpoints differ. The study already has that
number: per-seed success@5km of `v2single` on A256, a one-liner on the v3 flights. If the seed SD is ≲ 1 pp,
this is cosmetic; at ≳ 3 pp it decides ORB-1 verdicts.

**Proposed fix** (pick one; all are analysis-only):

1. **Recommended:** make `orb_paired` hierarchical. In each replicate draw 5 seeds with replacement, then
   resample cases within each drawn seed's stratum. This is the A256 scheme. With five clusters it is still a
   little liberal (B1), so also report the seed-level t-test next to it.
2. Keep the pre-registered test and state in the report and §16 that ORB-1/DS1/AS1 verdicts are conditional on
   the five trained checkpoints, not on the training recipe.

Either way, log it in §15 before any ORB-1/DS1/AS1 flight is analysed. **Impact (c):** no flight is touched and
`study_config.json` has no field for the bootstrap scheme.

### A2. `mppi_filtered` and the DAgger/BC expert do not use FNPEG's tuned filter gain β

**Where.**
* `jepa/fly_v5.py:358,363`: `FilteredPhysicsAdapter(beta=0.9)` raised to the power `dt_s` = 5.
* `jepa/fly_v5.py:392-403`: `load_adapter_v5` forwards `**kw` to the adapter, but no caller passes a β
  (`study_v5.py:819, 1278, 1280, 2012`). Its adapter cache key is `(kind, ckpt, device)`, so a fix that passes β must
  add β to that key, or a cached adapter with the old β is reused.
* `jepa/study_v5.py:628-629`: the `mppi_filtered` job carries no β, and has no `untuned` flag, so it can fly
  before FNPEG is tuned.
* `jepa/rl_v5.py:220,606`: `FastExpert` overwrites the tuned config's `filter_beta` with `RLCfg.expert_beta_1hz = 0.9`
  to the power 5.

**What the plan says.**
* `v5_plan.md:342`: `mppi_filtered` sets density and lift "every decision from FNPEG's filters".
* `v5_plan.md:328`: the filters are "K ← β·K + (1 − β)·…", and "β tuned".
* `v5_plan.md:848`: β ~ U(0.8, 0.97) is in FNPEG's tuning space (`FNPEG_SPACE`, `study_v5.py:1315`).
* `v5_plan.md:880`: "DAgger expert = FNPEG …, filters at β^5".

The code plumbs the tuned configuration into FNPEG (`tuned(L, "fnpeg")`) and into the expert's other settings, but
not β.

**Evidence.** The code lines above. `grep -n "filter_beta" jepa/*.py` finds no place where a tuned β reaches either
consumer.

**Consequence.** If tuning chooses β ≠ 0.9 (of the 36 tuning configurations, only config 0, the defaults, has
β = 0.9; `study_v5.py:126`):
* **G5** (Holm family, "learned model beats explicit estimation") pits v2 against an estimator that is *not*
  FNPEG's tuned filter.
* `dagger_fnpeg` and `bcppo_fnpeg` imitate an expert that is not the tuned FNPEG, so **R-imit** (`dagger_fnpeg −
  fnpeg`) mixes the imitation gap with an expert mismatch.

A reviewer can say, correctly, that G5 does not test the arm the paper describes. If tuning keeps β = 0.9 the
point is moot.

**Proposed fix.** Pass `tuned(L, "fnpeg", cfg)["filter_beta"]` (default 0.9) into the `mppi_filtered` model dict and
on to `FilteredPhysicsAdapter`, and give `mppi_filtered` jobs the `untuned` wait that the matched arm already has.
`FastExpert` should keep `fnpeg_cfg["filter_beta"]` raised to the power `dt_s` instead of overwriting it. The
alternative is to decide that β is fixed at 0.9 for these arms ("MPPI is not tuned") and say so in §4.2, §14 and §16.

**Impact (c)**, provided it is decided before the first `mppi_filtered` A256/ORB-1 flight and before Kaggle RL
tuning. E0 flew `physics_filtered_frozen` only (`study_v5.py:1284`), so no flown flight depends on β, and β is not a
frozen field. **If those arms fly first, the fix becomes (b) for them.**

### B1. The A256 bootstrap with five seeds under-covers when the seeds differ

**Where.** `jepa/analyze_v4.py:316-343` (`Boot.run`) and `jepa/analyze_v5.py:228-264` (`ClusterBoot`), used by every
A256 comparison: G1–G7, J-each and J-equiv, R-each.

**Evidence.** Synthetic A256 (`stats_test.py`, appendix A.2): 64 cases × 4 correlated atmospheres × 5 seeds,
deterministic comparator, true difference ≈ 0, 300 replications, `n_boot` 1,000. "met" means p < 0.05, the 95%
CI excludes 0, and ≥ 4 of 5 seeds hold, in either direction.

```
seed component of success SD   Boot p<.05   ClusterBoot p<.05   "met" rule   95% CI covers truth
0 pp                            5/300        6/300                5/300       294/300
~4 pp                          10/300       13/300               10/300       289/300
~11 pp                         34/300       36/300               32/300       264/300 (88%)
known effect (+10.7 pp):       183/300 power; coverage 265/300
```

The seed resampling uses 5 clusters, so its variance is too small by about (S−1)/S and its percentile tails are
too short. The 4-of-5 seed rule removes almost nothing (32 vs 34). The same narrow 90% CI makes **J-equiv's TOST**
(±0.15 log miss) declare equivalence too easily.

**Consequence.** Mild at a realistic seed spread (≤ 4 pp gives ≤ 5% false positives). Worth stating for RQ2/RQ3,
where the arms have learned seeds on both sides or RL seeds that differ widely.

**Proposed fix.** Report a seed-level t-interval (5 per-seed differences, 4 df) next to every A256 "met" and every
J-equiv, as a sensitivity in the style of the clustered interval. Flag a verdict that it does not support.
**Impact (c).**

### B2. RQ2's key comparison crosses GPU platforms, and on a T4 `amp_dtype` is emulated bf16

**Where.**
* `jepa/models_v2.py:686-690`: `amp_dtype` returns bf16 if `torch.cuda.is_bf16_supported()`, else fp16. The
  docstring says "fp16 on older ones (T4)".
* `jepa/models_v5.py:393-394`: E3 training uses `amp_dtype` and a GradScaler only for fp16.
* `jepa/models_v5.py:525-535` and `jepa/fly_v3.py:209-215`: the E3 and V2 adapters autocast to `amp_dtype` at flight
  time.

**Evidence.** The installed torch (2.14) defines `is_bf16_supported(including_emulation=True)`:

```
if torch.cuda.get_device_properties(device).major >= 8: return True
if not including_emulation: return False
return _check_bf16_tensor_supported(device)     # torch.tensor([1.0], dtype=torch.bfloat16, device=device)
```

A T4 is sm_75, and creating a bf16 tensor on it succeeds, so **`amp_dtype` returns bf16 (emulated) on the T4**, not
fp16. I could not run this on a T4; see §4.

`R0` (RQ2's reference, the v4 flights) trained and flew on an **L4**: native bf16, and TF32 active through
`allow_tf32` (`build_notebook_v4.py:40,131-132`). The E3 arms, including `r0_nojepa`, train and fly on the **T4**: emulated
bf16, no TF32. v5_plan §15 (`v5_plan.md:1237`) says a slower GPU "changes its duration, not its result beyond
rounding". That holds within v5, but J-each and J-equiv pair v5-T4 flights with v4-L4 flights. B3 shows that
"rounding" moves learned-MPPI landings by about 0.5 km.

**Consequence.** `r0_nojepa − R0` is meant to isolate "the JEPA objective, architecture held fixed" (§4.3). It also
changes the training and flight hardware and the numeric path. That is most likely noise plus a small unknown
systematic, but it is an unacknowledged confound in RQ2's key arm. §16 names two other confounds of this arm and
not this one.

Separately, emulated bf16 has no tensor-core path on a T4, so E3 training may be much slower than the unit
estimate. That is a budget risk, not a validity one.

**Proposed fix.**
1. Log `amp_dtype` and ms/step in the first T4 session. The notebook already logs torch/CUDA per session.
2. Fix the docstring.
3. Re-fly R0's five v4 checkpoints on the T4 (1,280 MPPI flights, an exploratory control). This separates the
   flight platform from the objective. Optionally retrain one R0 seed on the T4.
4. In any case, add the platform difference to §16's RQ2 limitation.

**Impact (c).** This adds an exploratory arm; nothing frozen changes.

### B3. Rounding and batch sensitivity is wider than §16 says

**Evidence 1: entry speed × (1 ± 1e-8)**, about 7.5e-5 m/s, on DEV48 cases 0–3, nominal atmosphere (`sens.py`,
appendix A.1). Miss in km.

| arm (case) | Win −1e-8 | Win 0 | Win +1e-8 | Linux −1e-8 | Linux 0 | Linux +1e-8 | max spread |
|---|---|---|---|---|---|---|---|
| predguid (0) | 2.330 | 2.050 | 2.329 | 2.914 | 0.723 | 0.458 | **2.46** |
| predguid (1) | 0.302 | 2.199 | 2.089 | 0.549 | 2.028 | 2.175 | **1.90** |
| predguid (2) | 3.484 | 2.791 | 1.501 | 0.638 | 0.414 | 1.757 | **3.07** |
| predguid (3) | 3.007 | 3.034 | 3.083 | 3.105 | 3.061 | 3.078 | 0.10 |
| fnpeg (0) | 0.5237 | 0.5285 | 0.5237 | 0.5141 | 0.5281 | 0.5281 | 0.014 |
| fnpeg (1) | 0.3661 | 0.3643 | 0.3719 | 0.3782 | 0.4382 | 0.3664 | 0.074 |
| fnpeg (2) | 0.7303 | 0.9394 | 0.9350 | 0.6554 | 0.6553 | 0.2229 | **0.72** |
| fnpeg (3) | 0.3079 | 0.3097 | 0.3012 | 0.3012 | 0.2990 | 0.2989 | 0.011 |
| classical (0–3) | 3.8468, 19.6133, 8.1533, 2.1510 | 3.8468, 19.6080, 8.1544, 2.1529 | (not finished) | 3.8471, 19.6069, 8.1542, 2.1513 | 3.8472, 19.6069, 8.1546, 2.1548 | 3.8472, 19.6112, 8.1539, 2.1518 | ≤ 0.007 |
| mppi_filtered, CPU (0) | 1.2938 | 1.2694 | 1.3288 | | | | 0.059 |
| v2single, smoke model, CPU (0) | 1.6337 | 2.1343 | 2.1986 | | | | **0.56** |
| RL smoke policies (dagger, bcppo, ppo; 0–1) | | | | | | | < 0.001 |

The smoke RL policies miss by 100–450 km and are not representative of trained ones.

**Evidence 2: batch composition** (`sens.py` with `--cases 0` vs `--cases 0,1`, same eps 0, same code, same
machine):

```
v2single case 0  alone 2.1342644733740403   in a lockstep batch with case 1  1.5848881010345008   (0.55 km)
mppi_filtered    alone 1.2694464657342934   in a lockstep batch with case 1  1.2691973168330821   (0.25 m)
fnpeg case 0     alone 0.5284865562475233   in a batch of 4                   0.5284865562475233   (identical)
classical        run 1 3.8468194229542587   run 2                             3.8468194229542587   (identical)
```

An MPPI flight's result depends on which flights share its lockstep shard. The likely cause is CPU kernels whose
vector and scalar paths round differently with an element's position, plus batch-size-dependent GEMM.
`run_mppi` shards a model run as `todo[i::workers]` (`study_v5.py:910`). After a disconnect *mid-run*, the remaining
flights are re-sharded, so they differ from what an uninterrupted run would give. That can be km-level for a learned
model. Smoke check 14 passes because its fail points sit between model runs.

**Consequences.**

1. **Paired comparisons.** This is zero-mean noise, not bias: success flips only near 5 km, and log miss gets
   per-flight jitter. It matters for single-case results.
2. **Tuning.** FNPEG's DEV48 misses are sub-km, and its tie-break after success@5km is mean log miss (`_choose`,
   `study_v5.py:1422`). A 0.2 → 0.9 km
   rounding swing moves one flight's log miss by 1.4, which is about 0.03 on a 48-flight mean. The choice among
   near-tied FNPEG configurations is therefore partly chance, as §16 already says for PredGuid. This does not bias
   any test flight.
3. **§16's wording.** "All study flights ran on one platform, Colab, so they are reproducible there"
   (`v5_plan.md:1286`) is true for the CPU-pool arms only. It is not true for MPPI flights flown after a mid-run
   resume, nor if Colab hands out CPU runtimes with different microarchitectures (AVX-512 vs AVX2 vector paths;
   not testable here).

**Proposed fix.** Documentation: extend §16 to name FNPEG (sub-km) and learned-model MPPI (≈ 0.5 km, batch- and
resume-dependent), and state that resume equivalence is bit-exact for the CPU-pool arms only. Optionally record the
CPU model in every flight record; `hardware_info` already reads it for RL. **Impact (c).**

### C. Documentation and defensive checks (all impact (c))

* **C1. Q2's interval is 95%, the plan says 90%.** `benchmark_v5.py:958` calls `clopper_pearson(k, n)` with the
  default `conf=0.95`. `v5_plan.md:890` says "ORB-1 Q1/Q2 per arm with Clopper–Pearson (90 %)". Q1 correctly uses
  90% one-sided (`requirement_lcb`). Fix either the call or the plan.
* **C2. The status total (22,408) is not the plan's Tier 1 total (≈ 27,900), and it will grow.** Reproduced from
  `all_stage_jobs` on an empty layout:

  | stage | flights |
  |---|---|
  | A256 | 1,792 |
  | ORB1_block1 | 5,000 |
  | E3 | 5,120 |
  | E4 | 3,840 |
  | explore | 1,024 |
  | DS1 | 2,048 |
  | AS1 | 1,024 |
  | **sum** | **19,848** (all flight IDs unique) |
  | + legacy | 2,560 |
  | **total** | **22,408** |

  - It omits gates and tuning (2 × 36 × 48 guidance flights plus RL tuning), which the plan counts.
  - It will jump by 1,280 once fnpeg's A256 flights exist (`jobs_matched` returns [] until then).
  - It will jump by about 1,140 E4 learning-curve flights once snapshots exist: DAgger 9 + BC→PPO 5 + PPO 5
    snapshots × 12 cases × 5 seeds.

  Say so next to the status table, so a growing total is not read as a bug.
* **C3. PredGuid: 4 of its 8 tuned parameters are inert on these entries.** Flying DEV48 case 0 (nominal) with Q7F
  8.5, Q19 0.35, KA 50, K_ρ 0.09, and all four changed together gives a commanded bank **bit-identical** to the
  default (miss 2.049964370617562 km, same hash of `sigma_cmd_deg`). LATBIAS 0.7, the positive control, moves the
  miss to 2.818 km (`pg_inert.py`).
  - **Confirmed on all 48 DEV48 conditions (2026-11-17, resolution).** All four at the low end, then all at the
    high end. In `pg_inert48.py`, every miss is exactly equal to the default's, but the hash of every telemetry channel
    differs. `pg_telkeys.py` names the one channel that differs: `dec_Krho_est`, the density estimate that K_ρ's gain
    sets, and the result field `predguid_Krho_final`. With that channel left out of the hash (`pg_inert48b.py`),
    48/48 flights are bit-identical. LATBIAS at its high end changes 39/48 (miss −2.4 to +4.4 km).
  - These parameters act in HUNTEST / up-control / the predictor (`_huntest` sets Q7, the Krho estimator feeds
    `_pg_predict` only).
  - These entries go initial roll → ballistic → final (§15's V4 note), so they never enter those phases.

  The tuning is a 4-D search (LATBIAS, KLAT, overcontrol, gmax_g), not 8-D. Say so in §14 and T10. I tested one
  case; the code path argues for all A256/DEV48 entries.
* **C4. The PredGuid design cache key omits design inputs.** `guidance_predguid.py:723-729` keys on IC, target,
  sigma_max, VFINAL, LAD, lod_fixed_frac, K44, KA.
  - `design_final_phase` and `_dry_run_to_final` also read HS, design_dt_s, gamma_safe_deg, hf_m and cycle_s
    directly. None of these varies in the study.
  - `_dry_run_to_final` also flies the longitudinal `executive` (`:408`), which reads the *tuned* gmax_g (`:429`),
    overcontrol (`:595`), Q7F, Q19 and K_RHO_FILTER_GAIN. None of these is in the key.

  This is harmless today, including for the tuning run. KA is in the key, and every tuning configuration draws its
  own continuous KA: round 2 draws inside the range and never clips onto a bound (`study_v5.py:1368`). So two
  configurations never share a cached design, and test flights use one chosen configuration. It becomes a
  stale-design bug if KA is ever held fixed while the others are tuned. Key on the whole configuration (e.g.
  `dataclasses.astuple(cfg)`).
* **C5. The Kaggle import does not check that the RL runs trained on the Colab study's inputs.**
  - `kaggle_v5.import_rl` checks the tar and file SHA-256s, the platform, and conflicts.
  - It does not compare the state's `00_study/study_config.json` and `03_tuning/fnpeg/chosen.json` with Drive's
    (`_import_target` does not even import them).
  - If FNPEG were re-tuned after RL started, policies trained against the old expert would import silently.

  Compare the state's `INPUT_FILES` with Drive's and refuse on a mismatch.
* **C6. Legacy arm names are not checked for uniqueness.**
  - `legacy_jobs` names an arm `rl_legacy_<file stem>` (`study_v5.py:1728`).
  - Two final models with the same stem in different subfolders of `models/` would share flight IDs, and the
    second would never fly.
  - Status would then show that arm below its total forever. You report 2,560/2,560, so this likely did not
    happen; I could not see the inventory.

  Assert unique names.
* **C7. DS1: heating telemetry and the RL heat-margin feature ignore the density profile.**
  `math_3d.nominal_heating_envelope_from_state` (`math_3d.py:637-668`, used by `rl_env._update_aux_metrics`) calls
  the atmosphere directly. The DS1 multiplier enters only `nominal_aero_forces_from_state` (the dynamics). This is
  the same treatment `density_scale` has always had (it scales area), so it is consistent, but DS1's heat columns
  are for the un-sheared atmosphere. Label them.
* **C8. The learning-curve snapshot marks assume the decision budget binds.** `_snapshot_marks` places
  `curve_points` at fractions of `budget_decisions` (200M; `rl_v5.py:428-441`). If the 3 h wall cap binds first on
  the Kaggle T4, `ppo_hist` and `bcppo_fnpeg` get fewer than 5 curve points, and the E4 curve is shorter than
  planned. Place the marks by elapsed fraction too, or report decisions reached (T07 does).
* **C9. A docstring is wrong.** `models_v2.amp_dtype` says "fp16 on … T4" (see B2).

---

## 3. Verified correct, and how

| # | Item | How checked | Result |
|---|---|---|---|
| 1 | All arms fly the same simulator settings | read `fly_v3.py:173-174` (classical), `mppi_v2.py:112` (MPPI), `fly_v5.py:137` (laws, policies) | `rl_env.OrionEntryEnv`, dispersion off, `ic_override`, `max_episode_steps` 1400, `fast_atmosphere` True everywhere |
| 1 | No truth leakage to classical, FNPEG, PredGuid, `mppi_filtered` | `rl_env.py:291-301` gives the guidance `dict(self.params)`, a shallow copy at construction; `data.apply_domain` (`data.py:127-163`) writes top-level keys of `env.params` only, including DS1's `density_profile` and AS1's `orion_trim_schedule`; FNPEG gets the domain only for `fnpeg_oracle`; PredGuid designs on the nominal model (ORB-1: the mission's nominal EI); `mppi_filtered` sees only measured accelerations | none |
| 1 | v3/v4 flights paired with v5 flights fly the same physics | `git log` on rl_env, control, point_math_3d, AtmosphereModel, ReactionControl, constants, default.json: last change 2026-06-24, before v3's runner (2026-10-18); `git diff 2565808 HEAD` of the v3/v4 runners: empty | same physics |
| 1, 11 | DS1/AS1 hooks are default-off and bit-identical | `git diff 2565808 HEAD`: only `math_3d.py` (+35 lines, a `density_profile` branch skipped when absent) and `data.py` (+27 lines in `apply_domain`) changed. One classical flight (dispersed domain 1.07/0.93/1.0) on a `git archive 2565808` tree vs the current tree: miss `8.821464817332688` both, peak g `4.088467251035855` both, telemetry SHA-256 `dedbfb02…1726` both. The AS1 path only sets a key `math_3d` already read (`math_3d.py:1521`) | bit-identical |
| 1 | MPPI settings and noise identical across MPPI arms | `mppi_job` (`study_v5.py:555-562`): `cfg.mppi()` (256 samples, 6 knots, 30°, n_sub 2), the same horizon, 4.5 g, `w_unc_km` 0; `crn_seed` as v3/v4; noise from a per-flight CPU generator | identical |
| 1 | Action scaling | `rl_env.step` maps a ∈ [−1, 1] to a·π/2 (`rl_env.py:638-641`); `PolicyLaw` maps `clip(a)·90°`; legacy and v5 policies agree with the training env; the 180° arms bypass the clip through `set_action` (`fly_v5.py:145-149`) | consistent |
| 7 | RL observation, full simulator vs fast env | `obs_check.py` (A.4): the same entry, the same command schedule, `FullObs` vs `FastEntryEnv.obs()` | identical at decision 0 in all 13 base features and the 168 history entries (diff 0.0). Later only physics drift (feature differences of about 1e-5 at 5 s, growing to about 1e-2 by 300 s) and the bank channel: the full simulator's roll loop overshoots (60° command reads 80° at 5 s, 45° at 10 s) and the fast env has no roll dynamics, as documented. Heading sign, cross-track, energy (G0), g margin (8 g) and history order all match |
| 7 | Legacy SB3 evaluation | `grep -rn VecNormalize`: used nowhere (`colab/build_colab_notebook.py:197`: "no VecNormalize is needed"); `LegacyPolicy` is deterministic at 1 Hz on `env._build_obs()`, with the 16-feature interval observation when the input is wider than 13 | correct |
| 2 | Dev/test separation | gates (`gate_suite` DEV48), `tune_guidance`, `tune_vf` and `tune_rl` all use `conditions(cfg, "DEV48")`; the E3 learning rate comes from `val_open_loop` on `ds_va` (`models_v5.py:474-481, 609-629`); `fast_eval` on A256 is used only for R-gap. Entry ICs: 0 shared between DEV48 and A256, 0 between ORB-1 and A256/DEV48 (computed). E0 uses A256 cases 0/16/40, but only as a bit-equality check of two runners, not a choice | separated |
| 3 | Metrics, hand-computed on 4 real flight records (smoke study folder): A256 fnpeg mismatch c0, A256 predguid nominal c16, ORB-1 fnpeg c0, ORB-1 classical c0 (`hand_metrics.py`, A.5) | miss: haversine from the final telemetry point to the config target; peak g; success 2/5/10 = drogue ∧ miss < 2 / < 5 / ≤ 10; reversals; g-duration compliance and margin (`B.g_duration` on g[1:]); heat load Σq̇·dt; ORB-1 downrange/crossrange in the nominal entry frame | all match. Miss to ≤ 0.08 m (telemetry stores degrees), peak g to ≤ 0.002 (substep vs 1 Hz), heat to ≤ 0.01%, ORB-1 ‖(dr, cr)‖ = miss |
| 3 | The success definition is the same everywhere | `fly_v3.py:144`, `analyze_v5.py:169-174`, `benchmark_v5.py:943-945`, `study_v5.py:1414-1416` | landed ∧ miss < 5 (≤ 10 at 10 km) everywhere |
| 4 | Holm, McNemar, Clopper–Pearson LCB | `stats_test.py` (A.2): `holm([.01,.04,.03,.2])` = [.04,.09,.09,.2]; exact McNemar on 12/3 discordant pairs = 0.035156 (by hand); `requirement_lcb` = a direct binomial-tail inversion for (1000,1000), (3000,3000), (2999,3000), (2998,3000), (998,1000) | exact. Q1 needs ≥ 1,705 flights with no miss to pass, so with block 1's 1,000 it is "pending" or "fail" in Tier 1, as §8 says |
| 4 | Bootstrap structure | `Boot`/`ClusterBoot` resample seeds then cases (per atmosphere / shared); a deterministic arm is broadcast over seeds; with no seed effect the null rejection is 5/300 and 6/300 and coverage 294/300; `_holds` counts per-seed signs correctly (planted 3-of-5 pattern gives (3, 5)); "met" needs p_Holm < 0.05; pending until n_pairs = n_expected | correct (seed caveat: A1, B1) |
| 5 | Classical = v3/v4's | the same function (`fly_classical_v4 = F.fly_classical`), bit-identical on the old and new trees; `classical_knows_truth` False | identical |
| 5 | FNPEG and PredGuid vs the papers | read `guidance_fnpeg.py` (RK4 predictor, Newton/secant with halving, deadband reversal, load relief, Brunner–Lu filters on measured D and L) and `guidance_predguid.py` (App. B phases, final-phase design, lateral logic) against Lu 2014 and Bairstow App. B and the earlier `guidance_audit.md` | no new discrepancy beyond C3 and C4 |
| 6 | Learned models | E3 normalization = v2's stats and `raw_features_v2` (`models_v5.py:130-141` vs `models_v2.py:100-102`); `r0_nojepa` saved as a v2 checkpoint and flown by the V2 adapter; R0 and every E3 arm are 1 member; E3 budget 36,000 × 256 = R0's 24,000 × 256 + 6,000 × 512 samples; parameter counts per §14 | consistent |
| 7 | Kaggle import safety | `kaggle_v5.import_rl`: tar SHA-256 vs manifest, per-file SHA-256, finished runs only (`policy.pt`), the platform equals the study's and every run's `hardware.json`, any differing file aborts before the first write, tar members checked for absolute paths and `..` | correct (gap: C5) |
| 8 | Bookkeeping | 19,848 Tier 1 flight IDs, 0 duplicates; telemetry written before the result JSON, both atomic (`study_v4.py:234-240`); `repair` removes stale `.tmp` and unparseable fresh JSON; ORB-1 `train_seed` = case mod 5 for all 3,000 cases (the strata are the seeds that fly) | correct |
| 8 | CPU-pool determinism | FNPEG case 0 alone vs in a batch of 4: `0.5284865562475233` both; classical rerun: `3.8468194229542587` both | bit-exact (MPPI is not: B3) |
| 10 | Completeness | every T01–T14b and F01–F10 the plan names is written by `analyze_v5` and printed by `report_v5` (T02 and T10 through the "every table" appendix, `report_v5.py:304`); `criteria()` has G1–G7e, G3, G6, J-each/equiv/mech, R-each/imit/gap, Q1–Q3 (+ T08b sensitivities), D1/D2 (+ estimates), A1/A2 (+ reach); §15/§16 printed from the plan | complete |
| 9 | Timing-only change to classical | `edbea3d`'s `classical_noscreen` is used only in `timing_stage` | flights unaffected |

---

## 4. What could not be verified, and why

* **Anything on a GPU** (no CUDA here):
  - which dtype `amp_dtype` really returns on Colab's T4 (inferred from torch's source: bf16 via emulation, B2);
  - E3 training speed in emulated bf16;
  - cuBLAS/cuDNN nondeterminism;
  - GPU-vs-CPU MPPI equality;
  - resume equivalence of E3 training on GPU (smoke check 12 is CPU-only).

  The notebook sets `allow_tf32` (a no-op on a T4) and no `deterministic` flags. MPPI noise comes from CPU
  generators, so the sample draws are reproducible; the model arithmetic is not guaranteed to be.
* **The Drive data**, which sizes several findings:
  - v3's per-seed `v2single` success spread (decides whether A1 and B1 matter);
  - the legacy inventory names (C6);
  - which CPU types Colab gave each CPU pass (B3).
* **Trained policies and models.** I had only smoke RL policies, which miss by 100–450 km, and a smoke v2 model.
  The rounding numbers for trained RL policies and full-size MPPI models may differ.
* **The full smoke suite.** Not rerun: check 10 takes hours, and check 11 rewrites `ReEntryAI_JEPA_v4.ipynb` inside
  the tree. I did items 11 and 8 by `git diff`, an old-tree flight and job counting instead.
* **The Kaggle side end to end** (smoke check 16 covers it); I reviewed the code only.

---

## Appendix: scripts and commands

All scripts ran from `StochasticEntrySim/` on the audit worktree. Linux runs used
`docker run --cpus 1 -v <sim>:/sim:ro -v <scratch>:/out v5audit-py313 python /out/sens.py /sim …`. The image is
`python:3.13-slim` plus `pip install numpy pandas pymsis gymnasium` and torch from the CPU index.

### A.1 `sens.py`: rounding sensitivity

```python
# usage: python sens.py <sim> <out.json> --arms fnpeg,predguid,classical --cases 0,1,2,3 --eps -1e-8,0,1e-8
#        [--policies name=path,...] [--v2ckpt ckpt.pt] [--threads N]
cfg = S5.StudyConfigV5()
conds = [c for c in S5.conditions(cfg, "DEV48") if c["case"] in cases and c["atmosphere"] in atms]
for arm in arms:
    for eps in eps_list:
        jobs = []
        for c in conds:
            c = dict(c); c["ic"] = dict(c["ic"]); c["ic"]["velocity_mps"] *= (1.0 + eps)
            jobs.append(S5.guidance_job(c, arm, "sens", arm, {}) if arm in ("fnpeg", "predguid") else
                        S5.classical_job(c, "sens") if arm == "classical" else
                        S5.policy_job(c, arm, "sens", pols[arm], 0) if arm in pols else
                        S5.mppi_job(cfg, c, arm, "G", 4.5, 0, model, "sens"))   # physics_filtered or v2 member 0
        # fly: run_guidance_lockstep(jobs, make_law(...)) | fly_classical_v4(j) |
        #      run_mppi_v5(jobs, load_adapter_v5(kind, ckpt, "cpu", members, 512, keep_only=True), cfg.mppi(),
        #                  n_sub=cfg.mppi()["n_sub"], max_rows=16384 if physics else 2048)
```

### A.2 `stats_test.py`: A256 bootstrap, seed rule, Holm, McNemar, LCB

```python
def sim_success(rng, effect, rho_case=0.8, seed_sd=0.3):          # [5, 4, 64] learned, [1, 4, 64] deterministic
    u = rng.normal(size=64); atm = np.array([0.0, -0.3, 0.2, -0.6])
    lat = sqrt(rho_case) * u[None, :] + atm[:, None]
    det = (lat + sqrt(1 - rho_case) * rng.normal(size=(4, 64)) > -0.2).astype(float)[None]
    s_off = rng.normal(scale=seed_sd, size=5)
    lrn = (lat[None] + s_off[:, None, None] + sqrt(1 - rho_case) * rng.normal(size=(5, 4, 64)) > -0.2 - effect)
    return lrn.astype(float), det
# per replication: A4.Boot(L, 1000).run([a, b], "diff_mean", None, seed=r), A5.ClusterBoot(...).run(...),
# A4.per_seed_tests + A4._holds for the 4-of-5 rule; truth = mean difference over 2,000 fresh draws
```

### A.3 `orb_test.py`: `orb_paired` under a null that holds over seeds

```python
for r in range(400):
    strata = np.arange(n) % 5; u = rng.normal(size=n); thr = norm.ppf(0.8)
    so = rng.normal(scale=seed_pp / 100 / norm.pdf(thr) * sqrt(2), size=5)      # seeds' offsets, mean 0
    a = (sqrt(.5) * u + sqrt(.5) * rng.normal(size=n) + so[strata] < thr).astype(float)   # learned arm
    b = (sqrt(.5) * u + sqrt(.5) * rng.normal(size=n) < thr).astype(float)                # deterministic arm
    res = A5.orb_paired(L, key, a, b, strata, 1000)                       # rejects if res["p_boot"] < 0.05
    ds = [(a - b)[strata == s].mean() for s in range(5)]                  # seed-level t-test, 4 df
```

The log labels the rows by `seed_pp` (0, 2, 4, 8). Because of the √2 in the latent scale, the SD of the true
per-seed success probability is √2 × `seed_pp` = 0, 2.8, 5.7, 11.3 pp, which is what A1's table reports.

### A.4 `obs_check.py`: RL observation, full simulator vs fast env

`GuidanceFlight` on DEV48 case 0 nominal, and a one-row `FastEntryEnv` whose pool is set from the same IC, domain and
target (as `rl_v5.fast_success_on` does). Both fly the command schedule 60° (150 s), −60° (150 s), then 30°. Every
5 s it compares `FullObs(env, physical_state(env))` with `env.obs()`.

### A.5 `hand_metrics.py`

It reads `04_flights/{results,telemetry}/<id>` of the local v5 smoke study and recomputes:
* miss (haversine, R = `constants.RADIUS_EARTH`);
* peak g, success flags, reversals, g-duration (`B.g_duration(g[1:], dt)`), heat load;
* ORB-1 (dr, cr) through `testset_v3.entry_frame` / `to_local`.

### A.6 Other commands

```
python pg_inert.py <sim> out.json 0                # PredGuid: default, Q7F, Q19, KA, K_RHO, all four, LATBIAS
python pg_inert48.py <sim> out.json 5              # C3 on all 48 DEV48 conditions (all four low / high, LATBIAS high)
python pg_telkeys.py <sim> 2                       # which telemetry channels differ (dec_Krho_est only)
python pg_inert48b.py <sim> out.json 6             # pg_inert48 with dec_Krho_est left out of the hash
python bitid.py <tree> out.json                     # one classical flight; run on `git archive 2565808` and on HEAD
python count_jobs.py <sim>                          # all_stage_jobs on an empty layout: counts and unique IDs
git diff --stat 2565808 HEAD -- jepa/fly_v3.py jepa/fly_v4.py jepa/mppi_v2.py jepa/models_v2.py jepa/physics_v4.py \
    rl_env.py control.py jepa/data.py math_3d.py   # data.py +27, math_3d.py +35, nothing else
```
