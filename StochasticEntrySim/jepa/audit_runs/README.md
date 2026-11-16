# Audit runs (2026-11-15)

Scripts and results behind `jepa/guidance_audit.md` and `jepa/simulator_fidelity_audit.md`. They read
the simulator and the v5 guidance code and change nothing; no study step imports them. Run each from
`StochasticEntrySim/`.

| Script | Result | Used in |
|---|---|---|
| `lu_cev_repro.py search` | `lu_cev_search.json` | guidance audit §2, fidelity audit §4.3: Lu's CEV mission (3,769 km) is out of reach |
| `lu_analog.py stage1` | `lu_analog_stage1.json` | guidance audit §3, stage 1: constant profile, 180° cap |
| `lu_analog.py stage2 2450` | `lu_analog_stage2.json` | guidance audit §3, stage 2: σf 30°/60°, 90° cap |
| `aero_sensitivity.py` | `aero_sensitivity.json` | fidelity audit §4.1: range with our trim L/D vs Bibb's |
| `check_fnpeg_predictor.py` | printed | guidance audit §1: FNPEG's predictor vs the simulator's equations (worst 0.29 km) |
| `classical_interval_screen.py` | `classical_interval_screen.json` | interval-layer audit §2: the classical baseline with its interval heat screen on and off |
| `check_euler_dt.py` | printed | fidelity audit: the 0.25 s Euler step's range error |

Each flight uses the full simulator and takes minutes; the `lu_*` scripts fly in parallel on all cores but
one (stage 2: about 10 minutes on a 12-thread laptop).
