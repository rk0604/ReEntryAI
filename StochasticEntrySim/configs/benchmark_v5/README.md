# v5 benchmark (`configs/benchmark_v5/`)

Frozen on 2026-11-03 by `jepa/benchmark_v5.py --build` (git commit recorded in `build_info.json`).
Every file here is listed in `SHA256SUMS` and verified whenever a suite loads
(`benchmark_v5.verify`); a mismatch raises. Nothing here is a simulator setting: the per-target
mission configs are copies of `configs/default.json` with only the target changed.

## Purpose

A fixed, documented set of test conditions for comparing entry-guidance controllers on the
ReEntryAI 3DOF simulator (Orion-like capsule, non-rotating spherical Earth, 1-D tabulated
atmosphere, bank actuator 20 deg/s, drogue deploy at 7.6 km). Plan and preregistration:
`jepa/v5_plan.md` section 3.

## Suites

| Suite | Conditions | Use |
|---|---|---|
| A256 | `configs/v3_testset.csv` (64 cases) x 4 atmospheres (`A256.json`) | results; continuity with v3 and v4 |
| ORB-1 | `ORB1.csv`, 3,000 dispersed cases, 3 blocks of 1,000 | results; the standard scenario |
| DEV48 | `DEV48.csv`, 12 cases x 4 atmospheres | tuning only, never a result |

### A256
The v3 test set unchanged (sha256 `cbf1b279854de260`), each case in four atmospheres:
nominal; thin (density x0.88, L/D x0.88); dense (x1.12, x1.12); out of range (x0.78, L/D x0.80).

### ORB-1 (Orion-class LEO return)
Nominal entry interface: 121.92 km (400,000 ft), 7,623.5 m/s Earth-relative (Lu 2014 Table 1,
CEV), flight-path angle -1.40 deg, at (0 N, 0 E) heading chi = 45.99 deg (east toward north).
Final condition: drogue deploy at 7.6 km (Lu's final altitude).

Targets on the entry great circle at 25 %, 50 % and 75 % of the nominal footprint between a
90 deg and a 0 deg constant bank, both flown in the full simulator from the nominal entry
(90 deg: 2273.3 km, 438 s, peak 9.4 g, ends with |cos gamma| < 0.05 at 8.5 km
altitude, falling almost vertically; 0 deg: 3379.2 km, 749 s, peak 2.5 g):

| Target | Fraction | Downrange (km) | Latitude (deg) | Longitude (deg) | Mission config |
|---|---|---|---|---|---|
| 0 | 25% | 2549.8 | 16.2559 | 16.3604 | `orb1_target0.json` |
| 1 | 50% | 2826.2 | 17.9608 | 18.2492 | `orb1_target1.json` |
| 2 | 75% | 3102.7 | 19.6477 | 20.1752 | `orb1_target2.json` |

Dispersions: Gaussian, truncated at +-3 sigma, independent, drawn case by case with
`numpy.random.default_rng(20261103)` in the order fpa, azimuth, speed, mass, density, cd, cl (rejection per variable):

| Parameter | Nominal | 3 sigma | Source |
|---|---|---|---|
| EI flight-path angle | -1.40 deg | 0.10 deg | Bairstow 2006 Table C.1 |
| EI azimuth | 45.99 deg | 0.10 deg | Bairstow |
| EI speed | 7,623.5 m/s | 15 m/s | decision D4 (no published LEO value) |
| vehicle mass | 8,500 kg | 5 % | Bairstow |
| density multiplier | 1 | 30 % | Bairstow ("day of flight") |
| CD multiplier | 1 | 5 % | Bairstow (hypersonic) |
| CL multiplier | 1 | 5 % | Bairstow (hypersonic) |

Simulator knobs (`data.apply_domain`): `density_scale = rho_mult x cd_mult`,
`lift_scale = cl_mult / cd_mult`, `mass_scale = mass_mult`, so drag scales as rho S CD c_D and
lift as rho S CL c_L exactly. Case i is in block i // 1000, flies target (i mod 1000) mod 3, and
learned arms fly it with training seed i mod 5.

Columns of `ORB1.csv`: `case, block, target_index, target_frac, target_downrange_km,
target_phi_deg, target_lambda_deg, config_path, train_seed`, the entry state (`altitude_m,
velocity_mps, gamma_deg, chi_deg, phi_deg, lambda_deg`), the multipliers and knobs,
`in_training_domain` (density and lift scale in 0.70-1.30, mass scale in 0.90-1.10: the world
models' training randomization), the target in the case's own entry frame (`downrange_km,
crossrange_km`), `margin_km` (signed distance of the target inside the 0-90 deg constant-bank
footprint under the case's true parameters, `testset_v3.footprint` with no load cap, the footprint
the targets are defined on; > 0 is `feasible`), `margin_km_4p5g` and `feasible_4p5g` (the same
with only the banks whose peak stays within 4.5 g, testset_v3's A256 rule: a conservative bound,
since a load-relieving guidance varies the bank), `max_bank_deg_4p5g`, and the lift-up flight's
time and peak g. 3000 of 3000 cases are feasible (100.00 %); 2025
are feasible within 4.5 g at a constant bank (target 0: 27, target 1: 999, target 2: 999); 2998 are inside the
training domain.

Sequential stopping: every arm flies block 1; after each block an arm with >= 2 flights that are
not (drogue and miss <= 10 km) has failed the requirement and stops.

### DEV48
`testset_v3` candidate filter (same `TestsetSpec`, Sobol seed 20261104, A256's Sobol indices
skipped), 12 cases, same columns as `configs/v3_testset.csv`, in the four A256 atmospheres.
Every tunable baseline is tuned here and only here.

## Scoring (`benchmark_v5.score`)
Per arm x suite x (atmosphere or target range): success within 2/5/10/25 km (outcome `drogue`),
mean log miss, the Orion EM-1 requirement test (P(miss <= 10 km) >= 0.99865 with 90 %
confidence: one-sided Clopper-Pearson lower bound, at most one miss in 3,000), NASA-STD-3000
g-duration compliance for a deconditioned crew (Bairstow Table D.3, App. D.2 procedure: 0.2 g
levels from the peak down to 1 g, time above each level against the log-log limit curve),
authority saturation, planning time per decision.

## Sources
Lu, P., "Entry Guidance: A Unified Method," JGCD 37(3), 2014. Bairstow, S. H., "Reentry
Guidance with Extended Range Capability for Low L/D Spacecraft," MIT SM thesis, 2006.
McNamara, L. W., and Braun, R. D., AAS 14-094 (Orion EM-1 entry requirements).

## How to cite
"ReEntryAI v5 benchmark (A256, ORB-1, DEV48), configs/benchmark_v5 at commit af99feb41c7a."

<!-- DS1 begin -->
### DS1 (density shear; added 2026-11-15)

A256's 64 entry cases and targets, each flown in four atmospheres whose density error changes with
altitude (`DS1.json`, built by `python jepa/benchmark_v5.py --build-ds1`; the other files are unchanged).
Every other suite's error is a constant multiplier, which FNPEG's and PredGuid's filters, mppi_filtered's
scales and v2's training randomization can all represent exactly; here none of them can. The factor
multiplies the simulator's density (`math_3d.density_profile_factor`, piecewise linear on a 2 km grid
from 0 to 160 km); `density_scale` and `mass_scale` stay 1.

| Atmosphere | Density factor | L/D |
|---|---|---|
| shear_dense_high | x1.15 above 70 km, x0.85 below 40 km, linear between | x1 |
| shear_dense_low | x0.85 above 70 km, x1.15 below 40 km | x1 |
| gram_like | exp(s(h) z(h)), z a unit Gaussian process (squared-exponential, 10 km correlation), s from 4% at 0 km to 12% at 100 km, one draw per case, truncated at 3 sigma | x1 |
| gram_like_lift | as gram_like, a second draw per case | x0.92 |

The random profiles are of the size Earth-GRAM reports for density dispersions (a few percent low down,
about 10 % at 80-100 km); they are not drawn from Earth-GRAM. Seed 20261115.
<!-- DS1 end -->

<!-- AS1 begin -->
### AS1 (aerodynamic sensitivity; added 2026-11-15)

A256's 64 entry cases and targets in nominal air, flown by a true vehicle whose trim aerodynamics follow
the published Orion database instead of the simulator's own schedule (`AS1.json`, built by
`python jepa/benchmark_v5.py --build-as1`; the other files are unchanged). The simulator's schedule has
more lift than the database below about 6 km/s (`jepa/simulator_fidelity_audit.md` section 4.1). Every
controller's internal model keeps the simulator's schedule, so this is an aerodynamic error that changes
with Mach number.

| Vehicle | L/D | CD | A256 targets still inside the 4.5 g footprint |
|---|---|---|---|
| bibb_ld | Bibb Part I Fig. 27 at FMV >= 8 (0.28 at FMV 8, 0.25 at FMV 25); capped at 0.28 below FMV 8 | the simulator's | 60 of 64 |
| bibb_ld_cd | as bibb_ld | Bibb Part II Fig. 40(b) and Part I Fig. 26: 1.34 at FMV 3, 1.38 at FMV 8 and above | 60 of 64 |

The values are read off the papers' figures (about +-0.03). `feasible` in `AS1.json` uses the rule and
code that made A256's own flags (`testset_v3.footprint`, 4.5 g, margin > 0).
<!-- AS1 end -->
