# Phase 0 — RL MDP Specification

**Project:** Can a reinforcement-learning *surrogate controller* outperform a
traditional predictor-corrector for Orion atmospheric entry guidance?

This document is the **design contract** for the RL experiment. No training code
lives here — it defines the Markov Decision Process (MDP), the baseline, the
evaluation protocol, and the Colab/logging plan that Phases 1–5 implement.

---

## 1. Thesis and what "outperform" means

**Hypothesis:** an RL policy can match a *competent* predictor-corrector on
nominal entries and **beat it under dispersion, failures, and multi-objective
constraints**, where the classical controller's constant-bank / single-objective
assumptions degrade.

"Outperform" is **not** a single nominal-case number. It is a distribution-level
claim over a held-out dispersed test set, on these metrics:

| Metric | Direction | Source column |
|---|---|---|
| Great-circle miss at drogue (km) | minimize | computed from `phi_rad`, `lam_rad` |
| Peak load factor (g) | minimize / constrain | `load_factor_g` |
| Peak heat rate (MW/m²) | minimize / constrain | `qdot_total_stag_hi` |
| Integrated heat load (MJ/m²) | minimize | `nominal_heat_Q_max_hi` |
| RCS propellant (kg) | minimize | `rcs_fuel_used_kg` |
| Success rate (% within target radius) | maximize | terminal flag |

A defensible result is a **Pareto comparison**: RL and classical each run on the
*same* held-out conditions; RL wins if it dominates on the joint objective or on
the dispersed/failure tail.

---

## 2. Scenario

The controlled segment is **entry interface → drogue deploy (7.6 km)** — exactly
what the predictor-corrector controls. CPAS handles everything below identically
for both controllers, so the comparison is apples-to-apples.

| | Scenario B (PRIMARY) | Scenario A (hard case) |
|---|---|---|
| Config | `configs/default.json` | `configs/scenario_a_lunar.json` |
| Regime | Orbital / LEO return | Lunar return |
| Entry altitude | 120 km | 300 km |
| Entry velocity | 7,800 m/s | 10,500 m/s |
| Entry γ | −2.5° | −12° |
| Capture | always (no skip-out) | narrow corridor |
| Peak g | survivable (~3–5 g target) | ~28–30 g (unsurvivable) |
| Target | (15°, 15°) ≈ 2,348 km, at ~45% of footprint | (15°, 15°), at ~35% of footprint |

Scenario B is the clean first result; Scenario A is the stress case for later.

---

## 3. The MDP

### 3.1 Observation  (≈13 normalized features)
A subset of the 186 logged columns, normalized to ~[−1, 1]:

| # | Feature | Normalization |
|---|---|---|
| 1 | altitude | / 120 km |
| 2 | velocity | / 11,000 m/s |
| 3 | γ (flight-path angle) | / (π/2) |
| 4–5 | sin χ, cos χ (heading; no wraparound) | — |
| 6 | range-to-go | / 2,500 km |
| 7 | heading error | / π |
| 8 | cross-track error | / 500 km |
| 9 | specific energy | / 6×10⁷ J/kg |
| 10 | current bank σ | / (π/2) |
| 11 | heat margin (q_lim − q)/q_lim | — |
| 12 | g margin (n_lim − n)/n_lim | — |
| 13 | altitude rate ḣ | / 300 m/s |

### 3.2 Action  (1-D continuous)
```
a ∈ Box(-1, 1, shape=(1,))   ->   σ_cmd = a[0] · (π/2)
```
Bank angle, ±90°. **Same handle and limits as the classical controller** — this
is what makes the comparison fair. Applied at the 1 Hz guidance cadence; the
existing bank actuator → roll PD → RCS → roll dynamics chain executes it.

> Note: the classical baseline decomposes σ into a predictor-solved *magnitude*
> and a deadband-reversed *sign*. The RL policy outputs σ directly (magnitude and
> sign jointly) — a strictly larger action space, so any RL win is not from
> having more authority.

### 3.3 Reward  (terminal-dominated, dense shaping)
Ingredients already computed per step in `telemetry.composite_reward`:
```
Per step (dense, small):
  r_t = − w_heat · max(0, q̇ − q̇_lim)
        − w_g    · max(0, n − n_lim)
        − w_fuel · ṁ_rcs
        + w_prog · Δ(range_to_go)          # optional progress shaping

Terminal at drogue deploy (large):
  r_T = − k_miss · miss_km
        + B_success · 1[miss < R_target]
        + B_survive · 1[peak_g < n_lim]

Hard failure (heat/g breach or NaN), terminated:
  r_fail = − P_violation
```
Start `k_miss`-dominant; tune weights empirically. `reward_default` already
implements the shaping sum as a baseline.

### 3.4 Episode lifecycle
- `reset(seed)` → sample dispersed initial state (Phase 2)
- `step(a)` → one guidance interval through the existing physics
- **terminated**: drogue altitude reached (success) **or** heat/g hard breach
- **truncated**: max steps / timeout
- `info`: full metric dict (miss, peak_g, peak_heat, fuel, success) for logging

### 3.5 Transition
The existing simulator, unchanged: bank actuator → roll PD → 12-jet RCS → roll
dynamics → 3DOF translational EOM → interval supervisor → heat shield → CPAS.
The **only** substitution is the guidance decision: the policy's `a` replaces
`SimpleBankGuidance.compute_sigma_cmd()` (via a `PolicyGuidance` subclass).

---

## 4. The baseline (the thing to beat)

The **continuous predictor-corrector** (`control.py`), now competent after:
- Fix 1 — golden-section bank-magnitude solver
- Fix 2 — long nominal rollout to drogue altitude
- Fix 4 — initial heading to target
- **Fix 3 — Apollo/Lu bank reversal with velocity deadband**

Current baseline performance on the default config: **6.8 km miss** (down from
1,385 km), in the neighborhood of Orion's ~5 km (2.7 nmi) targeting spec. This
strong baseline is what makes any "RL outperforms" claim meaningful rather than a
strawman win.

The baseline is run on the **same** held-out test set to produce its metric
distributions (Phase 3).

---

## 5. Evaluation protocol (fairness rules)

1. **Held-out test set** — a frozen set of 1,000 dispersed ICs
   (`configs/heldout_testset.csv`, built by `make_test_set.py` from fixed seed
   `20260604`), never used in training. The IC *values* are frozen (not just the
   seed), and both controllers replay the exact same conditions via
   `env.reset(options={"ic_override": row})`.
2. **Same observations / same action authority / same physics / same metrics.**
3. **Strong baseline** — the predictor-corrector must be the tuned, bank-reversal
   version (done), not a crippled one.
4. **Statistical comparison** — compare metric *distributions* (median, tail,
   success rate), not single runs. Report Pareto fronts for multi-objective.
5. **No reward hacking** — inspect agent behavior, not just return.

---

## 6. Colab + "log everything" architecture

**Packaging.** `git clone` the repo in Colab; `pip install pymsis numpy pandas
gymnasium stable-baselines3`. The sim modules import as-is; the env wraps them.

**Three logging layers:**
1. **Training curves** → Weights & Biases (survives Colab disconnects, shareable):
   reward, episode length, losses, explained variance.
2. **Per-episode metrics** → growing parquet: miss, peak-g, peak-heat, fuel,
   success, all ICs + seed. This is the RL-vs-classical comparison dataset.
3. **Full trajectories** → for eval/selected episodes, dump the full 186-column
   CSV (same schema `run_sim.py` produces) so any episode replays in the existing
   plots / dashboard.

**Persistence.** Mount Google Drive; checkpoint model + logs every N steps so
training resumes across the ~12 h Colab session limit.

**Speed.** The RL env step is *fast* — replacing the predictor drops its ~18 long
rollouts per cycle, leaving only physics propagation. PPO with 8–16 parallel envs
(`AsyncVectorEnv`) trains in hours. The slow predictor-corrector is run once for
the baseline (Phase 3), which can be done locally and uploaded.

---

## 7. Phase roadmap

| Phase | Deliverable | Status |
|---|---|---|
| **0** | This MDP spec | ✅ done |
| 1 | `rl_env.py` — `OrionEntryEnv(gymnasium.Env)`; validate by running the existing predictor-corrector through it | ✅ done |
| 2 | IC dispersion in `mission_config` + frozen held-out test set (`make_test_set.py`, `configs/heldout_testset.csv`) | ✅ done |
| 3 | Baseline characterization (predictor-corrector on test set) — runner `run_baseline.py` ✅; full run deferred to Colab | ⬜ Colab |
| 4 | RL training (PPO/SAC, Colab, W&B) — notebook `colab/ReEntryAI_RL_Colab.ipynb` ✅; run on Colab | ⬜ Colab |
| 5 | Fair evaluation + Pareto comparison — in the same notebook | ⬜ Colab |

**All heavy compute (Phases 3–5) runs on Colab**, not the local machine. See
`colab/ReEntryAI_RL_Colab.ipynb` (generated by `colab/build_colab_notebook.py`).

---

## 8. Open parameters to fix during Phase 1–2

- Target radius `R_target` for the success bonus (start 25 km)
- Reward weights `w_heat`, `w_g`, `w_fuel`, `w_prog`, `k_miss`, `B_*`, `P_violation`
- Load-factor limit `n_lim` (crewed: ~4 g nominal, ~8–10 g abort)
- Dispersion ranges: γ, V, χ, lat/lon, wind, chute-failure (Phase 2)
- Guidance cadence for the policy (start 1 Hz, matching the classical controller)

---

## 9. Interval-arithmetic integration (safety + performance)

The simulator already carries an **interval supervisor** that propagates a
guaranteed box of the state forward and bounds the worst-case heat rate, g, and
corridor. Originally only the *classical* predictor-corrector used it; the RL
policy was interval-blind. We now feed intervals into RL three independent,
ablatable ways (all in `OrionEntryEnv`):

1. **Interval observation** (`interval_obs`) — append worst-case features
   (`heat_margin_hi`, `g_margin_hi`, `box_width_norm`) so the policy is
   *uncertainty-aware*.
2. **Interval reward** (`interval_reward`) — penalize the interval **upper
   bound** on heat/g instead of the nominal point, so robustness to bounded
   uncertainty is learned into the weights.
3. **Interval safety shield** (`interval_shield`) — before a bank flies, a
   **1-step interval reachability check** vetoes it if the worst-case breaches a
   limit, projecting to the nearest safe candidate bank. The learned policy
   proposes, the sound supervisor disposes.

New `info` diagnostics: `peak_g_hi`, `interval_violations` (steps whose
worst-case breached a limit), `shield_intervention_rate`.

**Honest limits.** Interval bounds blow up over long horizons (wrapping), so the
shield is a *short-horizon* (default 1-step) filter, rigorous for the
*instantaneous* heat-rate and g limits, not a full-trajectory guarantee. The g
bound scales the nominal load factor by the dynamic-pressure upper bound
(assumes ~constant aero coefficient over the small box). Guarantees are
conditional on the disturbance model (the supervisor half-widths).

### Reframing the thesis
From a bake-off ("can RL match the predictor-corrector?") to **safety +
performance**: a *provably-safe RL surrogate* where the learned policy supplies
adaptivity under dispersion/failures and interval reachability supplies
**certified constraint satisfaction**. The contribution is the **hybrid**, which
directly answers the main barrier to learned controllers in crewed aerospace
GNC ("you can't certify a neural net"). Working title:
*Interval-Shielded Reinforcement Learning for Safe Atmospheric Entry Guidance.*

New evaluation axes beyond miss/g/heat/fuel: **constraint-violation rate**
(certified-zero for the shielded policy) and **shield-intervention rate**
(falls over training as the policy learns the safe envelope). Headline Pareto:
performance vs interventions vs violations, with shielded RL ideally achieving
classical-level safety and better adaptivity. The unshielded RL becomes the
ablation that shows why the shield is needed.
