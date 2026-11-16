# Guidance-law audit: FNPEG against Lu (2014), PredGuid against Bairstow (2006) (2026-11-15)

`jepa/guidance_fnpeg.py` checked against Ping Lu, *Entry Guidance: A Unified Method*, Journal of
Guidance, Control, and Dynamics 37(3), 2014, pp. 713–728, doi:10.2514/1.62605 (§1–4).
`jepa/guidance_predguid.py` checked against S. H. Bairstow, *Reentry Guidance with Extended Range
Capability for Low L/D Spacecraft*, S.M. thesis, MIT, 2006 (`research/Reentry Guidance with Extended
Range Capability for Low L-D.pdf`) (§5).

## 1. Equation by equation

| Lu (2014) | `guidance_fnpeg.py` | Verdict |
|---|---|---|
| Eq. 17: range to go integrated, ṡ = −V cos γ / r (the great-circle formula cannot tell undershoot from overshoot) | `_deriv_from`, 4th state | same |
| Eq. 18: 3DOF dynamics with Earth rotation in the predictor | longitudinal only (r, V, γ, s), no rotation: the formulation of Lu's ref. [30] (Lu 2008) | **differs**; Lu: "for short-range missions … makes little difference". Our simulator does not rotate either, so only the heading-offset term is left out; predictor vs simulator: worst 0.29 km over 36 flights |
| Eq. 19: \|σ\| linear in energy e = μ/r − V²/2 from σ0 to σf; σf prescribed or = σ0 | `predict`: `prof = S0 + (e − E0)/(ef − E0) (SF − S0)`; `sigma_f_deg=None` means σf = σ0 | same; σf is tuned (20–80°, or = σ0 with probability 1/4) |
| Eqs. 20–21: z(σ0) = s(ef) − sf, minimize z²/2 | corrector | same |
| Eqs. 22–25: step-size-controlled Newton (finite difference) then secant; λ = 1/2^i, smallest i with f decreasing | corrector (`fd_deg` = 1° first step) | same |
| Eq. 23: stop when \|z ∂z/∂σ0\| ≤ ε | corrector | same |
| Eq. 26: heading offset δψ = ψ − Ψ; reverse when \|δψ\| exceeds a velocity-dependent deadband | `decide`: deadband c0 + c1·V (tuned) | same (Lu leaves the deadband to the designer) |
| Guidance every 1 s while √(L² + D²) ≥ 1.52 m/s²; otherwise hold the current bank | `a_guid_mps2 = 1.52`; `cmd = Y[6]` below it | same |
| Fading-memory filters on L/L̂ and D/D̂, from Lu's ref. [4] (Brunner & Lu 2008) | `filter_beta = 0.9`, on above 0.05 g, tuned 0.80–0.97 | **unverified**: Lu gives no gain; Brunner & Lu (2008) is needed |
| Eqs. 36–40: ȧ = −2aD/V + aVβ_r sin γ; U_a = [amax − a(1 − 2Dδ/V)] / (aVβ_rδ); sin γ_ref = max(sin γ, U_a); ḣ_ref = V sin γ_ref | `relief_cos` | same (signs, the max, and SI units checked: U_a is dimensionless) |
| Eq. 41: L cos σ_cmd = L cos σ_base − k0 (ḣ − ḣ_ref), **inside every predicted trajectory** ("critically important") | `predict` applies `relief_cos` at every step | same |
| Point 3 / Eq. 42: once more on navigation data for the output command | `decide`: sensed L, D, a, V, γ | same |
| δ = 0.02 (16 s), k0 = 100 (nondimensional), amax margin 0.1 g | `lu_delta`: 0.02 × √(R0/g0) = 16.1 s; `lu_k0`: 100 × √(g0/R0) s⁻¹; `margin_g = 0.1` | same; δ and k0 tuned around Lu's values |
| CEV final condition: h_f = 7.6 km (25,000 ft), V_f = 150 m/s (Table 1) | `hf_m = 7600`; `vf_mps` set from DEV48 classical flights (the plan's choice) | same altitude; V_f is ours |
| Termination at e_f | stops at e_f or at h_f, whichever first | **differs** slightly (both define the same end on a nominal trajectory) |
| Bank limits: the capsule may use more than 90° (Fig. 15 reaches about 150°) | 90° for the primary arm (decision D2), 180° for `fnpeg_180` | **differs by design**; `fnpeg_180` is the published version (G7) |
| — | wings level when \|γ\| > 30° | ours (plan §15, 2026-11-03): only in the steep terminal descent |
| Bank rate and acceleration limits of the capsule: 20°/s and 10°/s² | the simulator's actuator: 20°/s and 40°/s²; the predictor models the 20°/s rate | the simulator reverses faster than Lu's model |
| 1976 US Standard Atmosphere, CEV aero (L/D ≈ 0.28), rotating Earth | the simulator's NRLMSIS table and trim schedule (`simulator_fidelity_audit.md`) | different vehicle and planet; FNPEG uses the simulator's own models, as Lu's uses his |

## 2. Lu's own mission cannot be flown here

Lu's CEV case (Table 1: 121.9 km, 7,623.5 m/s, 2,035 nmi = 3,769 km to Edwards AFB, unconstrained
peak 5.4 g) needs a γ0, which Lu does not give. Flying our FNPEG with a 180° cap and no load limit
from that state, γ0 = −1.0° to −2.0° in eight steps: every flight saturates at full lift-up (largest
bank command 0°, except 12° at −1.0°) and lands short, by 16 km at −1.0° to 883 km at −2.0°, with
peaks of 2.4–2.9 g. Our vehicle on a non-rotating Earth does not have the range; Lu's flies over a
rotating Earth (an eastward flight gains roughly 2ΩV cos φ ≈ 1 m/s² of effective lift) with a lighter,
wider capsule. A trajectory-level comparison with Lu's figures is therefore not possible.

## 3. Lu-analog test of the load relief (Figs. 14 and 17)

Same entry state with ORB-1's γ0 = −1.40°, Vf = 150 m/s, no errors, 180° cap unless noted, target
range chosen so that the unconstrained flight peaks near Lu's 5.4 g.

Stage 1 (constant profile, σf = σ0):

| Range | amax | Peak g | Miss | Largest bank |
|---|---|---|---|---|
| 2,450 km | none | 6.31 | 0.91 km | 105° |
| 2,450 km | 4 g | 5.31 | 0.22 km | 120° |
| 2,450 km | 3 g | 4.65 | 41.0 km | 180° |
| 2,200 km | none | 15.2 | 0.93 km | 119° |
| 2,200 km | 4 g | 6.05 | 190 km | 180° |
| 2,200 km | 3 g | 4.65 | 365 km | 180° |
| 1,950 and 1,700 km | any | 4.65–24.1 | 201–1,080 km | 180° |

Lu (Figs. 14 and 17), on a mission with a 5.4 g unconstrained peak: 4 g and 3 g held with misses
≤ 1.5 km; lowest feasible about 2.8 g; 2.7 g → 0.83 km, 2.5 g → 2.71 km, 2.3 g → 9.45 km, with the
load exceeded only in parts of the flight. Our 2,450 km case (6.3 g unconstrained) is the closest
analog: with the constant profile, 4 g becomes 5.3 g and 3 g becomes 4.65 g. The relief lowers the
peak by about 1 g but does not hold the limit, and it keeps the miss small where Lu's keeps the load.
Shorter ranges are not reachable within any limit (the target needs a deep dive).

The 3 g flights at three different ranges reach the identical peak (4.6544 g), all with the corrector
at the 180° cap: every one dives at full lift-down until the relief takes over, so they fly the same
trajectory up to the peak. That is the mechanism Lu describes ("the initial magnitude of the bank
angle … becomes larger as a ≤ amax is enforced"), driven to its limit.

Stage 2 (2,450 km): decreasing profiles (σf = 30° and 60°, as Lu allows), and the 90° cap of the
primary `fnpeg` arm:

| Cap | Profile | amax | Peak g | Miss | Largest bank |
|---|---|---|---|---|---|
| 180° | σf 30° | none | 7.63 | 7.6 km | 109° |
| 180° | σf 30° | 4 g | 5.65 | 20.4 km | 180° |
| 180° | σf 30° | 3 g | 4.65 | 41.2 km | 180° |
| 180° | σf 30° | 2.7 g | 4.31 | 129 km | 180° |
| 180° | σf 60° | 4 g | 5.44 | 1.3 km | 170° |
| 180° | σf 60° | 3 g | 4.65 | 41.2 km | 180° |
| 90° | constant | none | 6.31 | 2.9 km | 90° |
| 90° | constant | 4 g | 4.32 | 104 km | 90° |
| 90° | constant | 3 g | 3.65 | 292 km | 90° |
| 90° | σf 30° | 4 g | 4.32 | 104 km | 90° |
| 90° | σf 30° | 3 g | 3.65 | 285 km | 90° |

* **A decreasing profile does not help.** σf = 60° flies about like the constant profile (4 g: 5.44 g
  and 1.3 km against 5.31 g and 0.22 km); σf = 30° is worse in every row. The 3 g flights reach the
  same 4.65 g peak with every profile, because all of them dive at full lift-down first.
* **The 90° cap changes the trade, not the outcome.** Without full lift-down the corrector cannot dive
  to make up for the lift the relief adds, so the relief lowers the peak far more (6.31 → 4.32 g at
  4 g, → 3.65 g at 3 g) and the flight misses by 100–290 km instead. The limit is still exceeded, by
  0.3–0.65 g.
* **Across all 17 relieved flights of stages 1 and 2, none holds its limit.** At 2,450 km and the 180°
  cap the limit is exceeded by 1.3–1.7 g, and the miss is small only at 4 g with σf = σ0 or 60°
  (0.2–1.3 km; 20–130 km otherwise); at 90° the peak comes within 0.3–0.65 g of the limit and the
  miss is 100–290 km.

**What could not be separated.** Lu's flights hold 3–4 g with misses of about 1 km on a mission whose
unconstrained peak is 5.4 g; ours, on an analog at 6.3 g, do not. The equations match (§1), so the
difference is in the vehicle, the mission, the predictor's simplifications (§1: longitudinal, no
rotation), or in something the equation check cannot see (for example how fast the relief's lift
arrives after a full lift-down dive, with our faster actuator). Telling these apart needs Lu's own
case, which this simulator cannot fly (§2). So the audit verifies FNPEG's equations, not Lu's
closed-loop load-limit performance, and the paper should not claim that performance for our `fnpeg`.

## 4. What this means for v5

* **Describe `fnpeg` precisely.** "FNPEG (Lu 2014) with the load relief of Eqs. 36–42, implemented
  equation by equation; on a Lu-analog mission in our simulator the relief lowered the peak load by
  1.0–2.7 g but did not hold the limit." Not "Lu's load-limited guidance, which holds 3 g" (§3).
* **G3 and G4 already guard the comparison.** G3 checks that `fnpeg`'s mean peak on A256 is ≤ 4.6 g per
  atmosphere; where it is not, G4 (the load-matched comparison) is flagged invalid, so no result is
  called load-matched on the strength of a relief that did not hold. This test makes it likelier that
  G3 fails at some atmospheres, which is the honest outcome, not a defect of the plan.
* **Tuning cannot rescue it, and does the fair thing.** DEV48 tuning draws σf from 20–80° (or σf = σ0
  with probability 1/4), δ and k0 around Lu's values, and picks the best success@5km among settings
  with mean peak ≤ 4.6 g, else the lowest-peak setting (`study_v5._choose`). Stage 2 shows that σf
  does not buy load authority, so if no setting meets 4.6 g, the most load-conservative one is chosen,
  trading accuracy for load as Lu's design intends. `03_tuning/fnpeg/refine.json` records which.
* **The 90° cap matters most.** The primary arm flies at 90° (decision D2), where the relief brings
  the load closest to the limit but costs range. `fnpeg_180` (G7) shows the other side of the trade
  (accurate, overloaded). Report both, with peak-g columns, next to v2.
* **The dev entries are milder than this analog.** DEV48's unrelieved peaks are at most 4.2 g, below
  the study's 4.5 g limit, and at amax = 3 g the relief lowered them by 0.17–0.25 g (plan §14). The analog
  (6.3 g unrelieved) is a stress test of the relief, not a forecast of A256 numbers; G3 measures A256
  itself.
* **No code change follows from this audit.** Changing the relief now would be a design change made
  after looking at its behaviour; it would need its own pre-registration and would not be Lu's law.

## 5. PredGuid against Bairstow (2006)

**Constants.** Every Table B.4 constant in `PredGuidConfig` equals Bairstow's value: C1 1.25, C16
0.01, C17 0.001, C18 500 fps, C20 175 fpss, CHOOK 0.25, CH1 0.75, KA 64.4 fpss, KDMIN 0.5, KLAT
0.0125, K44 44,389,312, LATBIAS 0.4 nm, VLMIN 18,000 fps, VRCONTRL 700 fps, VCORLIM 1,000 fps, TOL
25 nm, VQUIT 1,000 fps, HS 28,500 ft, LEWD1 0.1, LEWD2 0.2, D0 130 fpss, D2 175 fpss, Q7F 6 fpss, Q19
0.2, ASENSIBLE 1.6 fpss. So does every Table B.7 constant the code uses (Bairstow's "new" column):
APOGEE_EPSILON_TIGHTER 40 nm, APOGEE_EPSILON_TIGHT 25 nm, ATMOS_ALT_MAX 600,000 ft, DELTA_T_PRED 10 s,
K_RHO_FILTER_GAIN 0.05, LIFT_INC_CAPTURE −3°, MAX_NUMBER_RUNS 10, and the first-call override of
Eq. B.59 (bank 90°, LIFT_INC_CAPTURE 10°, MAX_NUMBER_RUNS 20) is applied as written.

**Which of it our flights use.** With VFINAL = 29,000 fps, every LEO entry (24,400–27,600 fps) goes
from the initial roll to the final phase, sometimes through a short ballistic phase (smoke: "initial
roll → ballistic → final"). The Huntest, upcontrol and PredGuid
predictor-corrector branches (B.56–B.59, Table B.7, C19, Q7MIN) are implemented and checked but never
run. What the `predguid` arm flies is Apollo's final-phase guidance, with Bairstow's §5.5 reference
design and his constants for the initial roll, final phase and lateral logic. Say so in the paper:
on these entries "PredGuid" means Bairstow's Apollo-derived final phase, not his skip
predictor-corrector.

**Departures, all deliberate and recorded:**

| Bairstow | Ours | Why |
|---|---|---|
| GMAX 322 fpss (10 g, Apollo) | 4.5 g | the v5 load limit, the same for every arm |
| VFINAL 25,000 fps | 29,000 fps | plan §15 (2026-11-03): every LEO entry flies initial roll → final, the branch Apollo used below lunar-return speed (Putnam) |
| LAD 0.35 (Bairstow's CEV), LDCMINR = LAD cos 15° | 0.28 and LAD cos 15° | our capsule's hypersonic trim L/D |
| LOD 0.21 (60% of LAD), Q2, Q3, Q5, Q6 from Bairstow's reference | LOD by root finding to each case's target (option: a fixed fraction), Q2–Q6 from our §5.5 design | Bairstow's §5.5 procedure, applied to our vehicle and mission |
| CD_EST_INITIAL 1.25 | 1.55 | our trim CD |
| BANK_MAX 180° | 90° (`predguid_180`: 180°) | decision D2 |
| rotating, oblate Earth (EARTH_RATE, J2, flattening) | WIE = KWE = 0, spherical | our simulator does not rotate |

**Trajectory-level comparison: not possible.** Bairstow's nominal case is a lunar return (11,032 m/s
inertial, γ0 = −5.9°, L/D 0.35, Apollo-derived aerodynamics, rotating Earth) flown through the skip
branches. That is another speed regime, another vehicle and another planet model, and it is outside
the world models' training envelope; our LEO entries use only the initial-roll and final-phase
logic (above). What can be checked is what was checked: constants, and the gates (V4: PredGuid within 3.5 km
on at least 11 of the 12 DEV48 cases with the true parameters).

## 6. The filters against Brunner & Lu (2008)

Brunner, C., and Lu, P., *Skip Entry Trajectory Planning and Guidance*, JGCD 31(5), 2008,
pp. 1210–1219, doi:10.2514/1.35055 (Lu's ref. [4]; supplied by the user). Eq. 22:
X ← X + (1 − β)(X* − X), X* = measured / modeled acceleration, X0 = 1, applied to lift and drag, its
output multiplying the nominal L and D in the predictor, once per guidance cycle. `update_filters` is
the same: K ← βK + (1 − β)X*, K starts at 1, both K_D and K_L, used in `predict`. **Same.** The paper
gives no number for β, only "a value slightly less than one"; our 0.9 fits, and the tuning range
0.80–0.97 includes values below what the paper describes. Our 0.05 g switch-on is ours.
Brunner & Lu also note that the filter captures constant biases, and they added altitude-dependent
GRAM density errors on purpose so the algorithm "would not benefit unrealistically": the reason DS1
exists (plan §3.4, §16). Their vehicle (5 m, 8,382 kg, hypersonic L/D 0.289, CD 1.348, bank 20°/s and
10°/s²) is a second source close to ours above Mach 20 (simulator_fidelity_audit.md §4.1).

## 7. Still to verify (papers needed)

* **Lu, P., *Predictor-Corrector Entry Guidance for Low Lifting Vehicles*, JGCD 31(4), 2008,
  pp. 1067–1075, doi:10.2514/1.32055** (Lu's ref. [30]): the longitudinal predictor our FNPEG uses,
  and its guidance on choosing σf.
