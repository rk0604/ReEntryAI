# Simulator fidelity audit (2026-11-15)

Every number in v1–v5 comes from one simulator (`rl_env.OrionEntryEnv` → `point_math_3d` /
`math_3d`). There is no flight data. This file records what that simulator is, what was measured
about it, where it departs from the real vehicle and the published sources, and what the paper can
and cannot claim because of it. It replaces the March 2026 fidelity memo, which is out of date (it
describes the 1976 US Standard Atmosphere, constant aerodynamic coefficients and a 1.0 m nose radius;
the code now uses an NRLMSIS table, a Mach-scheduled trim database and a 6.0 m nose radius).

## 1. In one paragraph

The simulator is a careful 3DOF point-mass model with Orion-like mass, size and thrusters, an
empirical atmosphere and a published aerodynamic database. That makes it good for **comparing
controllers against each other**: every controller flies the same vehicle, the same errors and the
same entry states, paired case by case. It is not good for **absolute claims about Orion**. The Earth
does not rotate, navigation is perfect, there is no wind, the atmosphere is one fixed profile, and
(found in this audit) the vehicle has more lift in the second half of entry than Orion's own
database gives it. So "v2 beats FNPEG by X points in this simulator" is a defensible claim, while
"v2 would land Orion within 5 km Y% of the time" is not.

## 2. What the simulator models

| Component | In the simulator | Real vehicle / source | Effect on results |
|---|---|---|---|
| Dynamics | 3DOF point mass, spherical Earth (R = 6,378.1 km), inverse-square gravity (μ = 3.986e14), **no Earth rotation, no J2** | Orion flies over a rotating, oblate Earth; Lu (2014) simulates rotation | Range for a given trajectory differs by up to hundreds of km on eastward flights (§4.3); all arms share it |
| Integration | forward Euler at 0.25 s; guidance at 1 s (FNPEG), 2 s (PredGuid), 5 s (MPPI, RL) | — | ≤ 1.36 km difference in total range against 0.01 s steps (0.05 s: ≤ 0.23 km), over 1,400–3,700 km flights; closed-loop guidance corrects most of it (§3) |
| Attitude | roll only: bank actuator limited to 20°/s and 40°/s², PD roll controller, 12 RCS jets of 712 N, Izz 20,000 kg·m²; angle of attack always at trim | Orion: 12 × 160 lbf (712 N) RCS jets (Orion reference guide); Lu's capsule model: 20°/s and 10°/s² | Thrusters match Orion; the roll acceleration limit is 4× Lu's, so reversals are quicker than in Lu's runs; no pitch/yaw dynamics or trim transients |
| Vehicle | 8,500 kg, 5.0 m diameter (19.6 m²) | Orion CM 5.03 m (16.5 ft); landing weight 8,255 kg (Artemis I) and 9,299 kg (Artemis II) per the reference guide; Lu's CEV 8,383 kg, 5.5 m | Close to an uncrewed Orion; a crewed one is ~10% heavier (higher ballistic coefficient) |
| Aerodynamics | `orion_cm_trim`: CD and L/D vs Bibb's velocity parameter FMV; hypersonic L/D 0.27–0.40, CD 1.50–1.55 | Bibb et al. Part I, Figs. 26–28: nominal trim L/D 0.28 at FMV 8 falling to ≈ 0.245 at FMV 30–38, trim α 161–164°, trim CD ≈ 1.35–1.42; uncertainty ±0.035 in L/D | **L/D too high below FMV ≈ 25 (§4.1)**: ~23% more downrange authority than Orion's database gives |
| Atmosphere | NRLMSIS (pymsis) density table, 0–160 km in 100 m steps, built once at 0° N 0° E, 2026-01-01, F10.7 = 150, Ap = 7; used for every flight | Density varies with latitude, season, time of day and solar activity; GRAM adds perturbations | One mean profile; errors enter only as test multipliers (A256, ORB-1) and DS1's altitude-varying factors |
| Model errors (tests) | constant multipliers on density, lift and mass (A256, ORB-1); since 2026-11-15 also altitude-varying density (DS1, §3.4 of the plan) | Real errors vary with altitude, Mach and time; Bairstow's dispersions are of this constant kind | Constant errors are the best case for FNPEG's and PredGuid's filters and for v2's training randomization; DS1 tests the next step |
| Navigation | perfect: every controller reads the exact state and the exact sensed accelerations | IMU and GPS errors, filter lag | Every arm's accuracy is optimistic; none is helped more than another by construction, but noise-sensitive controllers would lose more |
| Winds | none during guided flight (the parachute model has a 5 m/s wind under canopy only) | Winds matter mostly below the drogue | Small effect on the drogue-deploy miss we score |
| Heating | Sutton–Graves convective + Tauber–Sutton radiative at the stagnation point, nose radius 6.0 m, on the **table density only** | Orion's heat-shield radius ≈ 6.0 m (1.2 D) | Heat metrics ignore the test's density multiplier (A256's density scale is applied through the reference area) and DS1's profile: they compare trajectories, not heating under the true air |
| Termination and miss | drogue deploy at 7.6 km; miss = great-circle distance at that point | Lu's final altitude 7.6 km (25,000 ft) | Parachute drift is not scored |

## 3. What was measured

| Check | Result | Meaning |
|---|---|---|
| FNPEG's predictor vs the simulator's own equations, 36 flights | worst 0.29 km in range | FNPEG's internal model is the simulator's (same trim table, same atmosphere): FNPEG is tested in its best case |
| Euler 0.25 s vs 0.01 s, constant-bank flights | ≤ 1.36 km in total range (0.05 s: ≤ 0.23 km) | Integration error is small next to the 5 km threshold and shared by all arms |
| Density hook (DS1): drag ratio with/without a profile vs the profile's factor | 4e-16 | The hook does exactly what it says; with no profile, a whole flight is bit-identical (smoke check 1) |
| Our aero vs Bibb's nominal L/D, ORB-1 nominal entry at constant bank | full lift-up range 3,379 → 3,144 km; 60° bank 2,647 → 2,578 km; peak g 2.46 → 2.89 and 4.00 → 4.34 | §4.1 |
| Lu's CEV mission (121.9 km, 7,623.5 m/s, 3,769 km) flown by our FNPEG | out of reach at every γ0 from −1.0° to −2.0°: the law saturates at full lift-up and lands 16–880 km short, peak 2.4–2.9 g (Lu: 5.4 g unconstrained) | §4.3 |

## 4. Findings

### 4.1 The aerodynamic schedule has more lift than Orion's database

`math_3d.aero_coefficients_from_speed_altitude`'s `orion_cm_trim` table cites Bibb et al. Part I, but
its hypersonic rows do not match that paper's nominal trim curve (read from Figs. 27 and 28b):

| FMV (≈ speed in kft/s) | 8 | 10 | 15 | 20 | 25 | 30 | 35 | 40 |
|---|---|---|---|---|---|---|---|---|
| Our L/D | 0.40 | 0.39 | 0.33 | 0.30 | 0.27 | 0.27 | 0.28 | 0.30 |
| Bibb nominal L/D | 0.28 | 0.28 | 0.26 | 0.26 | 0.25 | 0.245 | 0.245 | 0.245 |
| Difference | +0.12 | +0.11 | +0.07 | +0.04 | +0.02 | +0.025 | +0.035 | +0.055 |

Bibb's database uncertainty is ±0.035 in L/D (±0.05 with cg uncertainty), so our values are outside
it below FMV ≈ 20 (below about 6 km/s, the second half of the entry, where the final phase steers).
The table's comment also has the trim-α trend backwards: Bibb's nominal trim α rises from about
161.5° at FMV 8 to 163.7° at FMV 30; a trim α closer to 180° gives less L/D, not more. Our hypersonic
CD (1.50–1.55) is the drag near α ≈ 170–175°; at Bibb's trim α the database gives about 1.35–1.42.
Changing the CD to 1.38 alongside the L/D moves the constant-bank ranges by only 16 km; the L/D is
what matters.

Effect, measured on ORB-1's nominal entry: with Bibb's L/D the full-lift-up range is 235 km shorter
and the range span between 0° and 60° bank shrinks from 733 km to 567 km (−23%); peak loads rise by
0.35–0.45 g. Our vehicle can therefore correct larger errors late in the flight and stays further
below the 4.5 g limit than Orion would. Absolute success rates are optimistic for every controller;
ORB-1's farthest target (75% of the footprint, 3,120 km) would sit near the edge of an Orion-aero
footprint.

The schedule was not changed: v1–v4's 6,000 training flights, v2's world models and every earlier
result were produced with it, and changing it would invalidate experiments 1–4. A sensitivity run
with Bibb's curve is possible without touching the simulator's code (the `orion_trim_schedule`
parameter already overrides the table; this audit used it) and is offered as a follow-up.

**Below Mach 8 (Bibb Part II, Figs. 38–40; read off the figures, ±0.03).** Part II plots the nominal
trim CD against FMV for FMV ≤ 4 and the coefficients at α = 166° up to Mach 8; it has no trim L/D curve.

| FMV | 0.3 | 0.5 | 1.1 | 1.5 | 2 | 3 | 4 | 8 |
|---|---|---|---|---|---|---|---|---|
| Our CD | 0.84 | 0.85 | 1.28 | 1.38 | 1.43 | 1.47 | 1.49 | 1.52 |
| Bibb II nominal trim CD (Fig. 40b; Mach 8 from Fig. 38a at 166°) | 0.82 | 0.82 | 1.22 | 1.37 | 1.37 | 1.34 | 1.35 | 1.38 |
| Our L/D | 0.13 | 0.15 | 0.25 | 0.28 | 0.30 | 0.33 | 0.35 | 0.40 |
| Bibb II CL/CD at α = 166° (Figs. 38–39) | 0.25 | 0.25 | 0.25 | 0.23 | 0.23 | 0.22 | 0.22 | 0.21 |

Up to Mach 1.5 our drag is Bibb's (within 0.06). From Mach 2 to 8 ours is 4–10% higher. Our L/D
rises from 0.28 at Mach 1.5 to 0.40 at Mach 8, while Bibb's lift and drag are nearly flat there
(CL/CD 0.21–0.23 at 166°; the nominal trim α above Mach 1.5 is about 160°, Fig. 40a, where Part I
gives L/D 0.28 at Mach 8). So the L/D excess of this section starts near Mach 2–3 and is largest at
Mach 8. Subsonic, Part II's nominal trim α is about 170° (lower L/D than at 166°), so our 0.13–0.15
is plausible; that phase lasts seconds before the drogue and does not matter for range.

Two more sources agree with Bibb on the hypersonic value: Brunner & Lu (2008) fly an Orion-like
capsule with CL 0.3892, CD 1.3479, L/D 0.289 at high Mach; Strahan, Loe and Seiler (*Orion Entry
Flight Control Stability and Performance*, 2007 draft) write that the CM's L/D "is likely to drop
below the 0.35 preferred value", with guidance saturating for 30+ s in several Monte Carlo cases at
that L/D. Our 0.33–0.40 between Mach 3 and 15 is more lift than any of the three gives Orion.

**Roll control and actuator against Strahan et al.** The CM flies a PD attitude controller feeding a
jet trigger (Schmitt trigger or phase plane) and jet selection, RCS only; ours is the same structure
(`RollTorqueController` with deadbands → 12-jet RCS). They quote 20°/s attitude maneuvers (ours:
20°/s bank rate) and an 80 ms transport latency, which we do not model. Their ISS-return atmospheric
fuel budget is 100 lbm (45 kg, worst of 3,000 Monte Carlo runs); our flights' RCS fuel can be set
against it as a plausibility check, not as a validation (3DOF, roll only, no pitch/yaw damping fuel).

### 4.2 The atmosphere is one profile

The NRLMSIS table is built once at the equator, prime meridian, on 2026-01-01, for moderate solar
activity. Every flight in every suite uses it, and every controller's internal model uses the same
table. Day-of-flight density departs from a mean model and the departure changes with altitude;
Bairstow's dispersion for it is 30% (3σ), which ORB-1 uses as one constant multiplier and A256's
thin and dense atmospheres (×0.88, ×1.12) bracket. DS1 now adds the altitude variation. Latitude and
season are not modelled.

### 4.3 Lu's CEV case cannot be reproduced here

Lu (2014) flies the CEV from 121.9 km at 7,623.5 m/s to Edwards AFB, 2,035 nmi (3,769 km) away, over
a rotating Earth with the US Standard Atmosphere, and gets an unconstrained peak of 5.4 g. Our FNPEG,
in our simulator, cannot reach 3,769 km from that state at any γ0 between −1.0° and −2.0°: it flies
full lift-up and still lands 16–880 km short (lu_cev_search, 2026-11-15). The main reason is Earth
rotation: on an eastward flight it adds an upward acceleration of order 2ΩV cos φ ≈ 1 m/s² (about
0.1 g), which buys a lot of range; Lu's lighter, wider capsule (lower ballistic coefficient) helps too.
This confirms the plan's D3 note. A Lu-analog test of FNPEG's load relief is in
`guidance_audit.md` §3: on a 2,450 km analog (6.3 g unrelieved) the relief lowered the peak by
1.0–2.7 g but held its limit in none of 17 flights, where Lu's held 3 g and 4 g. Whether the vehicle,
the mission or the predictor's simplifications cause the difference cannot be settled without Lu's
own case.

### 4.4 Perfect navigation and no winds

Every controller is handed the exact state and exact sensed accelerations. Real guidance runs on a
navigation filter. This removes a source of error that would hit every arm, though not equally:
FNPEG's filters, `mppi_filtered` and the learned models all estimate the air from measured
accelerations, and noise would degrade each of those estimates differently.

### 4.4b The interval shield (rl_env) does not touch the results

`rl_env.OrionEntryEnv` can propagate an interval box of states (±10 m, ±20 m/s, ±0.15° γ, ±0.2° heading)
through the simulator's own equations and atmosphere table, and its shield vetoes a bank whose
one-step worst-case heat rate or g exceeds a limit. It never changes the physics: at most it changes
the command. Every v3–v5 flight builds the env with the shield, interval reward and interval
observation off (`fly_v5.GuidanceFlight`); only legacy SB3 policies trained with the 16-feature
observation get the interval observation, as they were trained. So it neither adds nor removes
fidelity here. The classical baseline separately runs a diagnostic interval heat screen that cannot
change its command; the arm-by-arm table is in `interval_layer_audit.md`. If the legacy "provably safe" shield is cited, qualify it: its bounds hold only inside
this simulator's models (so the L/D excess and single atmosphere pass straight into them), cover state
uncertainty but not the density, lift and mass errors v5 tests, look one step ahead, use an
approximate g bound (nominal g scaled by the dynamic-pressure bound), and are reset by recentering.

### 4.5 Heating does not see the test atmosphere

The heating model reads the table density. A256's density scale is applied through the reference area
(so lift and drag are exact), which the heating model never sees, and DS1's profile is applied inside
the aero function only. Heat rate and heat load in the tables are therefore the heating of each
trajectory in nominal air. They are secondary metrics and are labelled as such.

## 5. What the paper can and cannot claim

**Supported (relative, paired, in this simulator):**
* RQ1–RQ3 and DS1 compare controllers on the same simulated vehicle, the same errors and the same
  entry states, paired case by case with bootstrap intervals. A shared modelling error moves every
  arm; it does not by itself favour one. The integration error (≤ 1.4 km) is far below the 5 km
  threshold.
* The published baselines are implemented from the papers and checked against them where the papers
  allow (FNPEG's equations against Lu; PredGuid's constants against Bairstow), and FNPEG's predictor
  uses the simulator's own aero and atmosphere, so the baseline is not handicapped by model error.
  Not verified: Lu's closed-loop load holding, which a load-relief analog of his Figs. 14 and 17 did
  not reproduce (`guidance_audit.md` §3); G3 and G4 keep the comparison honest about it.

**Caveats a reviewer will raise, and how to state them:**
* Shared errors can still change the ranking. More late lift helps controllers that steer late; perfect
  navigation helps controllers that react to measured accelerations every step. State the ranking as
  conditional on this simulator, and report the DS1 and aero sensitivity results next to it.
* The learned models were trained on this simulator's own physics. Their advantage includes knowing
  its quirks (the excess L/D among them), which no flight vehicle would offer. A sim-to-real gap is
  untested.
* ORB-1 is an Orion-*class* scenario in this simulator, not Orion. Q1 applies Orion's accuracy
  requirement to the simulated vehicle; it cannot say that a controller "meets Orion's requirement".

**Not supported:** absolute accuracy, load or heating predictions for Orion or any flight vehicle;
any statement that a controller is flight-ready.

Suggested wording for the paper's limitations: *"All results come from a 3DOF simulation of an
Orion-class capsule over a non-rotating spherical Earth with perfect navigation, a single NRLMSIS
density profile and a trim aerodynamic schedule whose lift-to-drag ratio exceeds the published
Orion database below about 6 km/s. Comparisons between controllers are paired within this
simulator; absolute performance numbers are not predictions for a flight vehicle."*

## 6. Raising fidelity without flight data (in order of cost)

1. **Aero sensitivity** (cheap, no simulator code change): re-fly a subset of A256 and ORB-1 with
   Bibb's nominal L/D curve and its ±0.035 band. If the ranking holds, the excess-L/D caveat shrinks
   to a footnote. The learned models were trained on the old curve, so this is also a model-form-error
   test for them.
2. **Published flight reconstructions** as a validation target: NASA has published post-flight
   trajectory reconstructions of Orion's EFT-1 (2014) and Artemis I (2022) entries. With an entry
   state and the flown bank history from one of them, our simulator can be flown open loop and its
   altitude, speed and g-load histories compared with the real ones. This is the closest thing to
   flight data available.
3. **Earth rotation and J2** in the equations of motion (a simulator change; it would make every
   earlier result non-comparable, so it belongs to a new study, as would 4–5).
4. **Navigation errors and sensor noise**: an IMU error model on the sensed accelerations.
5. **Atmosphere variability**: NRLMSIS by latitude and date per case, or Earth-GRAM profiles.

## 7. Sources used

Bibb, Walker, Brauckmann, Robinson, *Development of the Orion Crew Module Static Aerodynamic Database,
Part I: Hypersonic* (Figs. 26–29) and *Part II: Subsonic–Supersonic* (Figs. 38–40); Strahan, Loe,
Seiler, *Orion Entry Flight Control Stability and Performance* (2007 draft); Brunner and Lu, *Skip
Entry Trajectory Planning and Guidance*, JGCD 31(5), 2008; NASA *Orion Reference Guide* (2022); Lu, *Entry Guidance: A
Unified Method*, JGCD 37(3), 2014; the March 2026 ReEntryAI fidelity memo; the simulator's code
(`math_3d.py`, `point_math_3d.py`, `rl_env.py`, `AtmosphereModel.py`, `heat.py`, `constants.py`).
