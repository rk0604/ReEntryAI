# Where the interval layer was on or off (audit, 2026-11-15)

Read-only audit of the code and the generated notebooks; no study code was changed. "Interval layer"
means three separate things in this repository:

1. **`rl_env`'s interval options** (`rl_env.py:172-174`): `interval_obs`, `interval_reward`,
   `interval_shield`. Any of them sets `use_intervals` (`rl_env.py:211`), builds `supervisor_cfg`
   (`rl_env.py:265`) and propagates an interval box each substep (`rl_env.py:521-539`). Only
   `interval_shield` can change a command (`rl_env.py:643-644`), and only in policy mode.
2. **The classical guidance's interval heat screen** (`interval_screen_active`,
   `math_3d.py:2178`): `SimpleBankGuidance` rolls its chosen bank through interval dynamics and the
   interval heat-shield model (`control.py:843-856`). It is driven by the guidance's context
   (`GuidancePredictionContext.interval_active`, default **True**, `control.py:307`).
3. **`run_sim.py`'s interval supervisor** (`run_sim.py:321`, `936`, `1086-1092`): the original
   single-trajectory script and `ReEntryAI_run2.ipynb`. No JEPA study uses it.

## 1. Arm × stage × flags

Every environment the JEPA studies build is listed; a flag not passed keeps `rl_env`'s default (False).

| Arm | Stage | Env construction | interval_obs / reward / shield → use_intervals, supervisor_cfg | Classical heat screen |
|---|---|---|---|---|
| all learned models (v1–v5) | (a) training data | `jepa/data.py:224-226`, policy mode, random bank profiles | off / off / off → False, None | not present (policy mode has no `SimpleBankGuidance`) |
| v1 MPPI + world model | (b) closed loop | `jepa/mppi.py:225-227`, policy mode | off → False, None | not present |
| v2 MPPI (v2–v5), MPPI with the physics model (`mppi_nominal`, `mppi_filtered`), JEPA-free models (`mlp_direct`, `gru_direct`, `gru_pi`, `r0_nojepa`) | (b) | `jepa/mppi_v2.py:110-112`, policy mode (v3–v5 fly through this runner) | off → False, None | not present |
| classical `SimpleBankGuidance` | (b) v1/v2 | `jepa/mppi.py:225-227`, baseline mode | off → False, None | **ON** (context default; `rl_env` never calls `set_prediction_context`) |
| classical | (b) v3–v5 | `jepa/fly_v3.py:172-174`, baseline mode | off → False, None | **ON** |
| classical | compute timing (v5 T09) | `jepa/study_v5.py:1830-1831`, baseline mode | off → False, None | **ON** |
| FNPEG, PredGuid (+ `_180`, `_oracle`) | (b) v5 | `jepa/fly_v5.py:107-109`, policy mode | off → False, None | not present |
| DAgger, BC → PPO, PPO, SAC | training | `jepa/rlenv_fast_v5.py` (torch physics port; no interval code) | not applicable | not present |
| DAgger, BC → PPO, PPO, SAC | (b) v5 | `jepa/fly_v5.py:107-109`, policy mode | off → False, None | not present |
| PPO in the full simulator (Tier 2) | training | `jepa/rl_v5.py:887-889`, policy mode | off → False, None | not present |
| legacy SB3 policies (`rl_legacy_*`) | (b) v5 | `jepa/fly_v5.py:107-109` + `env_kwargs` from `jepa/study_v5.py:1583-1587` | **interval_obs = True if the policy's input is wider than 13**, reward off, shield off → use_intervals True, supervisor_cfg built | not present |
| legacy RL (the old Colab notebook) | training and its own evaluation | `colab/build_colab_notebook.py:276-289`, `393-396` | committed values: all three **False**; the cell tells the user to flip them per run, so what each saved policy was trained with is not recorded in the repository | not present |
| legacy BC dataset | data | `bc_dataset.py:58-59`, baseline mode | off | **ON** |
| all | (c) open-loop evaluation | no environment: the models are scored on stored episodes from stage (a) | not applicable | not applicable |
| test sets and footprints | (d) | `jepa/build_notebooks_v3.py:256` and `jepa/benchmark_v5.py:312-313`, policy mode, constant bank | off → False, None | not present |

Notebooks: of the JEPA notebooks only `ReEntryAI_JEPA_v3_0_testset.ipynb` builds an environment
itself (`OrionEntryEnv(controller_mode="policy", dispersion=False, fast_atmosphere=True)`, flags off);
the others call the modules above. No JEPA notebook (v1, data, v2, v3_0–v3_3, v4, v5) contains
`interval_shield`, `interval_obs`, `interval_reward`, `use_intervals`, `supervisor_cfg`,
`set_prediction_context` or `interval_screen`. `colab/ReEntryAI_RL_Colab.ipynb` has
`INTERVAL_OBS = INTERVAL_REWARD = INTERVAL_SHIELD = False`. `ReEntryAI_run2.ipynb` (not a JEPA study)
runs the interval supervisor with `interval_active = True`.

## 2. The three questions

**1. Was `interval_screen_active` ever True in a saved classical run? Yes: in every one.** The
classical baseline of v1, v2, v3, v4 and v5 ran the screen at every guidance cycle, because nothing
turns it off. In the continuous-predictor mode these studies use (`rl_env.py:296`), the screen is
called *after* the bank is chosen (`control.py:827-839` choose, `843-856` screen) and its result only
fills debug fields (`control.py:857-915` returns the same `sigma_cmd_rad`). It cannot veto.
Checked by flying (`jepa/audit_runs/classical_interval_screen.py`, three held-out entries, nominal
air): with the screen on, as every study flew it, the guidance reports it used in 83–84% of the
steps; with it forced off, every state at every step is bit-identical (largest difference 0.0),
with the same miss (1.42, 1.87, 2.45 km) and peak g.

**2. Do the saved flight logs record shield vetoes? No.** `rl_env`'s `info` carries `use_intervals`,
`interval_active`, `interval_violations`, `shield_interventions` and `shield_steps`
(`rl_env.py:806-815`), but the flight records keep none of them: a v5 result JSON has 58 fields and
the telemetry 14 arrays, neither with a shield or interval field; the `flights.csv` of the v3, v4 and
v5 smoke studies have no such column either. So a veto count per arm cannot be read from the logs.
From the code the count is zero for every arm: `interval_shield` is never passed as True anywhere
under `jepa/`, so `_shield` is never called. The saved Drive studies were not inspected (Drive is not
mounted on this machine); they were written by the same code.

**3. Arms not run under identical flags.**

* **Classical vs everyone else.** The classical arm is the only one whose guidance runs an interval
  heat screen. It does not change its commands, but it is inside the classical arm's per-step wall
  time: `plan_ms` for classical is "the whole 1 s guidance step (predictor + physics)"
  (`jepa/fly_v3.py:194`), and v5's compute table T09 times the same steps (`jepa/study_v5.py:1827-1839`).
  In the three check flights the screen added 0%, 23% and 45% to the flight's wall time (one run
  each on a busy laptop, so a rough figure). The classical compute numbers therefore include work the
  other arms do not do, and overstate what that guidance law needs.
* **Legacy SB3 policies with a 16-feature input** get `interval_obs=True`; every other arm flies
  with the interval propagation off. The physics is the same; these policies see three extra
  observation features and pay for the interval propagation in wall time.
* **Training vs test for the RL arms.** They train in `rlenv_fast_v5` (no interval code, heat margin
  from Sutton–Graves only) and fly in `rl_env`, whose heat-margin feature is the upper bound of the
  interval heat-shield model's stagnation heat rate including the radiative term
  (`rl_env.py:449-452`), with or without `use_intervals`. Plan §14 already records this; the gap is
  reported per arm.
* **The heat-rate numbers of every arm** come from that same interval heat-shield object fed the
  nominal state (`math_3d.py:637-675`), so "peak heat rate" is its `.hi` value for all arms alike.

## 3. What this means

* No arm's trajectory was altered by a shield or a screen. The "provably safe" shield of the legacy
  RL work was not part of any v1–v5 result and must not be described as if it were.
* Fix for the paper, without touching results: say that the classical baseline carries a diagnostic
  interval heat screen, and report classical compute separately (or re-time it with the screen off:
  `env.guidance._prediction_context.interval_active = False` after `reset`, as
  `jepa/audit_runs/classical_interval_screen.py` does).
* The flight records should carry `use_intervals` and `shield_interventions` if a later study turns
  the shield on; today they would not show it.
