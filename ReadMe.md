# ReEntryAI

A high-fidelity **3-DOF atmospheric entry simulator** for the Orion Crew Module,
built to test a research question: *can a reinforcement-learning controller match
or beat a classical predictor-corrector for entry guidance, while a provable
interval-arithmetic safety layer guarantees heat-rate and g-load limits are never
violated?*

The **simulator, classical guidance, parachute model, uncertainty (interval)
layer, mission-control dashboard, and RL training pipeline are all built and
working.** The **RL-beats-classical result does not yet exist** — see
[Where the RL stands](#where-the-rl-stands) for an honest account.

---

## Table of contents
- [What's in the box](#whats-in-the-box)
- [Repository layout](#repository-layout)
- [Quick start](#quick-start)
  - [1. Run the classical simulator](#1-run-the-classical-simulator)
  - [2. Train / evaluate the RL policy (Colab)](#2-train--evaluate-the-rl-policy-colab)
  - [3. Launch the mission-control dashboard](#3-launch-the-mission-control-dashboard)
- [The simulator in detail](#the-simulator-in-detail)
- [Mission configs](#mission-configs)
- [Outputs](#outputs)
- [The RL pipeline](#the-rl-pipeline)
- [Where the RL stands](#where-the-rl-stands)
- [Known limitations](#known-limitations)
- [Requirements](#requirements)

---

## What's in the box

| Capability | Status | Where |
|---|---|---|
| 3-DOF entry dynamics (explicit Euler, 0.25 s) | ✅ | `math_3d.py`, `point_math_3d.py` |
| NRLMSISE atmosphere (+ fast lookup table) | ✅ | `AtmosphereModel.py` |
| Reduced-order Orion aero (FMV-scheduled CD / L·D) | ✅ | `math_3d.py` |
| Convective + radiative heating, wall temperature | ✅ | `constants.py` (`HeatShield`) |
| Classical predictor-corrector guidance + bank reversals | ✅ | `control.py` |
| 12-jet RCS + bank actuator + roll control | ✅ | `ReactionControl.py`, `control.py` |
| CPAS parachutes (reefing, squidding, pendulum, failures) | ✅ | `cpas.py` |
| Interval-arithmetic uncertainty supervisor | ✅ | `interval_math.py`, `math_3d.py` |
| Gymnasium RL environment | ✅ | `rl_env.py` |
| IC dispersion + frozen held-out test set | ✅ | `mission_config.py`, `make_test_set.py` |
| Classical baseline characterization (parallel) | ✅ | `run_baseline.py` |
| Behavior-cloning warmstart pipeline | ✅ | `bc_dataset.py`, Colab |
| React mission-control dashboard | ✅ | `mission_control_ui/` |
| **A trained RL policy that beats / matches classical** | ❌ **open** | see below |

---

## Repository layout

```
ReEntryAI/
├── ReadMe.md                     ← you are here
├── StochasticEntrySim/           ← the simulator + RL (the heart of the project)
│   ├── run_sim.py                ← canonical CLI driver (config-driven)
│   ├── run_baseline.py           ← parallel classical-baseline runner over a test set
│   ├── make_test_set.py          ← generates the frozen 1,000-case held-out set
│   ├── rl_env.py                 ← OrionEntryEnv: the Gymnasium RL environment
│   ├── bc_dataset.py             ← collects (obs → classical bank) pairs for imitation
│   │
│   ├── point_math_3d.py          ← nominal closed-loop translational step (the physics core)
│   ├── math_3d.py                ← 3-DOF EOM, aero schedule, heating envelope, interval supervisor
│   ├── interval_math.py          ← Interval class + interval arithmetic / box helpers
│   ├── control.py                ← guidance stack (predictor-corrector + policy), actuator, roll ctrl
│   ├── ReactionControl.py        ← 12-thruster RCS model
│   ├── cpas.py                   ← parachute assembly (drogue/pilot/main, reefing, squid, pendulum)
│   ├── constants.py              ← physical constants, HeatShield (Sutton-Graves / Tauber-Sutton)
│   ├── AtmosphereModel.py        ← NRLMSISE atmosphere + tabulation + interval atmosphere
│   ├── telemetry.py              ← derived metrics, RCS propellant, RL reward weights/terms
│   ├── mission_config.py         ← JSON config loader, IC dispersion, targets, aero application
│   ├── plotting.py               ← shared plotting + plot-category filtering
│   │
│   ├── configs/                  ← mission definitions (JSON) + the frozen test set (CSV)
│   ├── colab/                    ← the RL training notebook (run heavy compute here)
│   ├── docs/                     ← technical overview PDF, RL MDP spec
│   ├── revision_v1/  revision_v1_py/   ← example outputs (notebook vs run_sim.py) + figures
│   ├── runs/                     ← example CPAS failure-mode run
│   ├── data/  data_intv/         ← sample logged data
│   ├── helpers/                  ← one-off notebook-surgery scripts (dev tooling, not runtime)
│   ├── old_code/                 ← archived earlier versions
│   ├── heat.py, math_2d.py, RL_ready_sim/  ← legacy / superseded modules
│
├── mission_control_ui/           ← React + Vite "mission control" dashboard
│   └── src/components/            ← telemetry plots, 3D trajectory, ground track, CPAS panels
│
├── RL_Model/                     ← early 2-D RL prototype notebooks (historical)
└── testSim/                      ← early simulator prototype (historical)
```

> **Note on `_py` folders.** `run_sim.py` writes to `<output_dir>_py/` so its
> outputs never clobber the notebook's `<output_dir>/`. That's why you see paired
> folders like `revision_v1/` (notebook) and `revision_v1_py/` (`run_sim.py`).

---

## Quick start

```bash
git clone https://github.com/rk0604/ReEntryAI.git
cd ReEntryAI/StochasticEntrySim
pip install numpy pandas matplotlib pymsis
```

### 1. Run the classical simulator

`run_sim.py` flies one mission end-to-end (entry → parachutes → splashdown),
writing a per-step trajectory CSV, a summary JSON, and a folder of diagnostic
figures. It is driven by a JSON config selected with the `SIM_CONFIG` env var:

```bash
# default mission
python run_sim.py

# a specific scenario
SIM_CONFIG=configs/scenario_artemis_skip.json python run_sim.py
```

Sanity-check the RL environment's physics against the classical controller:

```bash
python rl_env.py configs/default.json     # runs OrionEntryEnv in baseline mode
```

### 2. Train / evaluate the RL policy (Colab)

All heavy RL compute runs on **Google Colab**, not locally (the env needs
`gymnasium` + `stable-baselines3`, and training wants a GPU/many CPUs).

1. Open `StochasticEntrySim/colab/ReEntryAI_RL_Colab.ipynb` in Colab
   (File → Open notebook → GitHub → `rk0604/ReEntryAI`).
2. Run **Section 0** (clones the repo, installs deps, mounts Drive).
3. Then either:
   - **Section 6 — Imitation warmstart (recommended):** clone the classical
     controller, evaluate it, optionally RL-fine-tune. *(This is the current
     frontier — see status below.)*
   - **Sections 2–5 — from-scratch baseline + PPO:** the original path (the
     from-scratch PPO does **not** converge; kept for reference).

The notebook is generated from `colab/build_colab_notebook.py` — **edit the
builder and re-run it**, don't hand-edit the `.ipynb`.

### 3. Launch the mission-control dashboard

```bash
cd mission_control_ui
npm install
npm run dev          # opens a Vite dev server; serves trajectory data live
```

The dashboard renders a flown trajectory (3-D path, ground track, telemetry,
CPAS sequencing, RCS activity). Data is read from the simulator's output folder
(`revision_v1/`); static copies live in `public/data/` for production builds.

---

## The simulator in detail

The state vector is `[r, φ, λ, V, γ, χ]` (radius, latitude, longitude, speed,
flight-path angle, heading). One guidance decision is made each second; the
physics integrates at `dt = 0.25 s` with explicit forward Euler.

| Module | Responsibility |
|---|---|
| **`point_math_3d.py`** | The nominal closed-loop step: synchronize control state, run guidance → bank actuator → roll torque → RCS firing → realized bank, then advance the translational state. The physics core both `run_sim.py` and `rl_env.py` call. |
| **`math_3d.py`** | The equations of motion; the FMV-scheduled Orion-like aero (CD and L/D vs a blended Mach/velocity parameter, per Bibb et al.); the interval **heating envelope**; and the **interval supervisor** that propagates a guaranteed state box and bounds heat/g. Also the predictor-corrector rollout used to score candidate banks. |
| **`control.py`** | `SimpleBankGuidance` (the classical predictor-corrector: golden-section bank-magnitude solver + long rollout + Apollo/Lu-style bank reversals), `PolicyGuidance` (lets an RL policy supply the bank), the observation provider (range-to-go, heading & cross-track error), the bank actuator (rate/accel limited), and the roll controller. |
| **`ReactionControl.py`** | 12-jet body-frame RCS, bang-bang firing, propellant accounting. |
| **`constants.py`** | The `HeatShield` class — Sutton-Graves convective + Tauber-Sutton radiative + Tauber-Wakefield coupling + radiative-equilibrium wall temperature on a discretized shield grid — plus all physical constants and interval-supervisor defaults. |
| **`AtmosphereModel.py`** | NRLMSISE (via `pymsis`), an optional baked log-density lookup table for RL throughput, and the interval-atmosphere hull used by the uncertainty layer. |
| **`cpas.py`** | The Capsule Parachute Assembly: drogue/pilot/main sequencing, staged reefing, squidding (stochastic partial collapse), pendulum limit-cycle (on 2-of-3 main failure), forward-bay-cover jettison, and stochastic chute failure. |
| **`interval_math.py`** | The `Interval` dataclass and conservative interval arithmetic (sin/cos/sqrt/exp/pow, division guards) plus box/vector helpers — the basis of the guaranteed-bounds layer. |
| **`telemetry.py`** | Derived physics (load factor, specific energy, etc.), the RCS propellant model, and the RL reward (`RewardWeights`, `composite_reward`). |
| **`mission_config.py`** | Loads a JSON mission, applies aero, builds the entry-interface state, samples **IC dispersions**, and resolves mission targets. |

---

## Mission configs

`configs/*.json` define a mission (entry state, aero, CPAS, guidance, interval
limits, target, and an optional `dispersion` block). Pick one with `SIM_CONFIG`.

| Config | What it is |
|---|---|
| `default.json` | Nominal orbital entry from ~120 km. The reference mission. |
| `scenario_a_lunar.json` | Faster lunar-return entry. |
| `scenario_artemis_skip.json` | Artemis-style lunar **skip** entry (dip → skip-out → land). Generates the dashboard's data. |
| `example_dispersed_entry.json` | Entry with IC dispersion enabled. |
| `example_failure_2of3.json` | Forces a 2-of-3 main-chute failure (triggers the pendulum). |
| `heldout_testset.csv` | The **frozen 1,000-case** dispersed test set (seed-locked) used for every RL-vs-classical comparison. Built by `make_test_set.py`. |

The `dispersion` block scatters the six entry-interface variables around nominal
(e.g. ±1 km altitude, ±100 m/s speed, ±0.4° flight-path angle, ±3° heading,
±0.3° lat/lon) — realistic delivery-error Monte-Carlo scatter around one planned
entry.

---

## Outputs

Each run produces, in `<output_dir>[_py]/`:

- **`trajectory_success.csv`** — ~190 columns of per-step state, control, RCS,
  interval bounds, heating, and guidance diagnostics.
- **`trajectory_failed_heat_steps.csv`** — guidance cycles flagged heat-infeasible.
- **`run_summary.json` / `landing_summary.json`** — episode-level summary
  (miss distance, peak g, peak heat rate, RCS fuel, outcome).
- **`figures/`** — ground track, 3-D trajectory, density/dynamic-pressure
  envelopes, heating, heat-shield maps, interval-width growth, RCS activity,
  and CPAS diagnostics. Filter with `PLOT_CATEGORIES`.

---

## The RL pipeline

The reinforcement-learning layer wraps the simulator as a standard Gym task.

- **`rl_env.py` — `OrionEntryEnv`:** observation = 13 normalized features
  (+3 optional interval features); action = 1-D continuous bank angle (±90°);
  episode runs from entry interface to drogue deploy. Reward = dense progress
  shaping (range-to-target closed) + optional heading shaping + hinge penalties
  on heat/g + an active RCS-fuel penalty + terminal miss/on-target/survivability
  bonuses. Three **interval-RL** modes can be toggled independently:
  - `interval_obs` — append worst-case heat/g headroom to the observation,
  - `interval_reward` — penalize the *guaranteed worst-case* heat/g, not nominal,
  - `interval_shield` — veto a bank whose 1-step interval reachability breaches a
    limit (a hard safety constraint, not a soft penalty).
- **`make_test_set.py`** — generates the seed-locked dispersed test set.
- **`run_baseline.py`** — runs the classical predictor-corrector over that test
  set in parallel to produce the distribution RL must beat.
- **`bc_dataset.py`** — runs the env in baseline mode and records
  `(observation → classical bank command)` pairs for behavior cloning.
- **`colab/ReEntryAI_RL_Colab.ipynb`** — the end-to-end training/eval notebook
  (baseline characterization, PPO training, imitation warmstart, head-to-head
  evaluation, W&B logging).

The intended paper contribution is the **interval-shielded RL** angle: a learned
controller made *provably safe* by the interval reachability shield.

---

## Where the RL stands

**Honest status: the RL does not yet match the classical controller. The
infrastructure is complete; the headline result is not.**

The classical predictor-corrector lands a **median 2.05 km** from target with
**99.9% success** across the 1,000-case test set. That is the bar. Every RL
attempt so far falls well short:

| Approach | Success | Median miss | Verdict |
|---|---|---|---|
| Classical predictor-corrector (the bar) | 99.9 % | **2.05 km** | — |
| PPO from scratch, v3 (weak reward) | 1.3 % | ~198 km | lazy local optimum |
| PPO from scratch, v4 (strong reward + exploration) | ~1 % | ~91–198 km | flails, never converges |
| **Behavior clone of the classical controller** | 14.8 % | **91 km** | best learned policy so far |
| Clone → PPO fine-tune | 0 % | 438 km | **regressed** (cold-critic collapse) |

**What we learned:**
- **From-scratch PPO cannot crack this problem** — the success region (a precise
  bank schedule with reversals) is too hard to discover by exploration, so the
  policy collapses to either doing nothing (v3) or flailing (v4).
- **Behavior cloning helps but isn't enough.** Cloning the classical controller
  gets to 91 km, but suffers the classic imitation failure: *compounding errors*
  (the policy drifts into states the expert never demonstrated) and the expert's
  bang-bang bank reversals get "averaged out" by a smooth policy.
- **Naive clone→RL fine-tuning makes it worse.** Behavior cloning trains only the
  policy network, leaving the value (critic) network random. PPO then reinforces
  actions using garbage value estimates and destroys the good clone before the
  critic warms up — a well-known BC→RL failure mode. Result: 438 km.

**Where we're leaving off — open directions (none yet chosen):**
1. **DAgger** — iteratively query the classical controller on the *states the
   clone actually visits* and retrain. Directly attacks the 91 km compounding-
   error problem, with no RL instability. Most promising for a strong learned
   controller.
2. **Proper clone→RL handoff** — warm up the critic on frozen-clone rollouts,
   fine-tune with a tiny learning rate and a KL leash to the clone. Salvages the
   "RL beats classical" framing, more engineering, uncertain payoff.
3. **Reframe around safety** — lead the paper with the simulator + interval
   propagation + the provable safety shield, using the cloned controller as the
   learned policy the shield protects, and drop the "beats classical accuracy"
   claim. Most defensible given current results.

Best artifact today: the behavior clone (`bc_orion_entry_v1`, ~91 km). All runs,
models, and metric CSVs are logged to Weights & Biases and Google Drive.

---

## Known limitations

These are deliberate scope choices / unvalidated assumptions, stated plainly:

- **First-order integration.** Explicit forward Euler at 0.25 s — no RK4, no
  convergence study.
- **Non-rotating Earth.** No Coriolis / Earth-rotation terms.
- **Reduced-order aero.** A scheduled Orion-*like* CD/L·D model, **not** the
  official Orion aerodatabase.
- **Calibrated, not validated, radiative heating.** The Tauber-Sutton velocity
  function is anchored to an Apollo-4-style calibration, not validated against a
  reference trajectory. Heating magnitudes are plausible but unverified.
- **Approximate interval g-bound.** The live shield's g bound scales nominal load
  factor by the dynamic-pressure ratio (a rigorous interval-g variant exists as a
  local experiment but is not in the training loop).
- **No automated test suite.**
- **RL training is Colab-only** (gym/SB3 are not installed for local runs).

---

## Requirements

- **Simulator:** Python 3.10+, `numpy`, `pandas`, `matplotlib`, `pymsis`.
- **RL (Colab):** `gymnasium`, `stable-baselines3>=2.3`, `wandb`, `torch`
  (installed by the notebook).
- **Dashboard:** Node 18+, `npm` (React + Vite; see `mission_control_ui/`).

See `StochasticEntrySim/docs/ReEntryAI_Technical_Overview.pdf` and
`docs/rl_mdp_spec.md` for deeper technical write-ups.
