"""
One-shot updater: wires the high-fidelity CPAS module into
ReEntryAI_run2.ipynb. Idempotent -- safe to run multiple times.

Modifies:
  Cell 2  (imports)   -- adds `import cpas`
  Cell 13 (main loop) -- instantiates CPAS, calls step() each iteration,
                         applies all the new effects to the state and
                         params, logs the expanded cpas_* columns
  Cell 16 (save data) -- writes cpas_events.csv and adds the CPAS
                         summary block to run_summary.json
"""
import json

NB_PATH = r"C:\Users\Risha\Desktop\ReEntryAI\StochasticEntrySim\ReEntryAI_run2.ipynb"

with open(NB_PATH, "r", encoding="utf-8") as f:
    nb = json.load(f)


def cell_src(i):
    return "".join(nb["cells"][i]["source"])


def set_cell(i, new_src):
    lines = new_src.split("\n")
    nb["cells"][i]["source"] = [ln + "\n" for ln in lines[:-1]] + ([lines[-1]] if lines[-1] else [])


# -----------------------------------------------------------------------------
# Cell 2 - imports : add `import cpas`
# -----------------------------------------------------------------------------
src2 = cell_src(2)
if "import cpas" not in src2:
    src2 = src2.rstrip() + "\nimport cpas\n"
    set_cell(2, src2)
    print("cell 2: added `import cpas`")


# -----------------------------------------------------------------------------
# Cell 13 - main loop
# -----------------------------------------------------------------------------
src13 = cell_src(13)

# 1) CPAS instantiation right before the for-loop
init_inject_marker = "rows = []"
init_inject_block_v1 = """rows = []
failed_heat_rows = []
thruster_fire_events = []  # one row per (step, thruster_name) fire

# --- CPAS (parachute) integration ---
cpas_inst = cpas.CPAS()
cpas_events_log = []   # one row per phase transition
params['cpas_drag_cdA_extra_m2'] = 0.0
params['cpas_lift_scale'] = 1.0"""

init_inject_block_v2 = """rows = []
failed_heat_rows = []
thruster_fire_events = []  # one row per (step, thruster_name) fire

# --- CPAS (parachute) integration ---
import random as _rng_mod
cpas_inst = cpas.CPAS(rng=_rng_mod.Random(42))
cpas_events_log = []   # one row per phase transition
params['cpas_drag_cdA_extra_m2'] = 0.0
params['cpas_lift_scale'] = 1.0
_R_EARTH = float(constants.RADIUS_EARTH)"""

# Replace whichever variant is present
if init_inject_block_v1 in src13:
    src13 = src13.replace(init_inject_block_v1, init_inject_block_v2)
elif init_inject_block_v2 not in src13:
    # Fresh notebook -- inject at the marker
    old_block = """rows = []
failed_heat_rows = []
thruster_fire_events = []  # one row per (step, thruster_name) fire"""
    assert old_block in src13, "cell 13: init block not found"
    src13 = src13.replace(old_block, init_inject_block_v2)

# 2) Replace the CPAS step + termination block at the top of the loop.
#
# We accept either: (a) the freshly-extracted plain loop with no CPAS,
# or (b) a previously-injected v1 CPAS block. Both are rewritten to the
# v2 high-fidelity block.

# Variant A: plain (no CPAS yet)
old_gate_A = """for k in range(num_steps):
    t_s = k * dt_s
    sigma_actual_rad = float(att_nom.sigma_rel_rad)
    roll_rate_rad_s = float(att_nom.omega_b_rad_s[2])

    if abs(math.cos(x_nom[4])) < 0.2:
        termination_reason = 'nominal_cos_gamma_too_small'
        print(f'Stopping at t = {t_s:.2f} s: cos(gamma) too small')
        break
    if x_nom[3] <= 1.0:
        termination_reason = 'nominal_speed_too_low'
        print(f'Stopping at t = {t_s:.2f} s: speed too low')
        break
    if alt_from_r(x_nom[0]) <= 0.0:
        termination_reason = 'nominal_ground_reached'
        print(f'Stopping at t = {t_s:.2f} s: ground reached')
        break"""

# Variant B: v1 CPAS block already present (from previous updater run)
old_gate_B_marker = "cpas_out = cpas_inst.step(t_s=t_s"

new_gate = """for k in range(num_steps):
    t_s = k * dt_s
    sigma_actual_rad = float(att_nom.sigma_rel_rad)
    roll_rate_rad_s = float(att_nom.omega_b_rad_s[2])

    # --- High-fidelity CPAS step ----------------------------------------
    alt_now_m = alt_from_r(x_nom[0])
    _atm_now = atmosphere_from_state_vector(x_nom)
    _T_now = _atm_now.get('T_K', float('nan'))
    if _T_now and _T_now == _T_now and _T_now > 0:
        mach_now = float(x_nom[3]) / math.sqrt(1.4 * 287.0 * float(_T_now))
    else:
        mach_now = None
    cpas_out = cpas_inst.step(t_s=t_s, alt_m=alt_now_m, V_mps=float(x_nom[3]),
                              mach=mach_now, dt_s=dt_s)
    params['cpas_drag_cdA_extra_m2'] = float(cpas_out.drag_area_cdA_m2)
    params['cpas_lift_scale']        = float(cpas_out.lift_scale)
    if cpas_out.mass_shed_delta_kg > 0.0:
        params['mass_kg'] = max(1.0, float(params['mass_kg']) - cpas_out.mass_shed_delta_kg)
    for ev in cpas_out.events:
        cpas_events_log.append({
            't_s': t_s, 'step': k, 'alt_m': alt_now_m, 'V_mps': float(x_nom[3]),
            'mach': mach_now if mach_now is not None else float('nan'),
            'phase': cpas_out.phase, 'event': ev,
            'reefing_stage': cpas_out.reefing_stage,
            'num_mains_active': cpas_out.num_mains_active,
        })
        print(f'  CPAS {ev:>18s} at t={t_s:6.1f}s  alt={alt_now_m:7.1f}m  V={float(x_nom[3]):6.1f}m/s')

    # cos_gamma gate only when truly stowed (chutes -> vertical descent OK).
    if cpas_out.phase == 'stowed' and abs(math.cos(x_nom[4])) < 0.2:
        termination_reason = 'nominal_cos_gamma_too_small'
        print(f'Stopping at t = {t_s:.2f} s: cos(gamma) too small')
        break
    if x_nom[3] <= 0.1:
        termination_reason = 'nominal_speed_too_low'
        print(f'Stopping at t = {t_s:.2f} s: speed too low')
        break
    if alt_now_m <= 0.0:
        termination_reason = 'nominal_ground_reached'
        print(f'Stopping at t = {t_s:.2f} s: ground reached')
        break"""

if old_gate_A in src13:
    src13 = src13.replace(old_gate_A, new_gate)
elif old_gate_B_marker in src13:
    # v1 CPAS already in place. Replace from "for k in range" up through
    # the alt termination break (inclusive) with the new block.
    import re
    pattern = re.compile(
        r"for k in range\(num_steps\):.*?if alt_now_m <= 0\.0:.*?break",
        flags=re.DOTALL,
    )
    src13 = pattern.sub(new_gate, src13, count=1)
else:
    print("[WARN] could not find a known loop top to replace in cell 13")

# 3) Inject wind drift + pendulum lateral velocity AFTER the step result
#    is consumed and BEFORE the next iteration. We'll do it right before
#    "x_nom = list(x_nom_next)" near the bottom of the loop.
wind_inject = """    # Apply wind drift and pendulum lateral velocity under chutes (post-step).
    if cpas_out.phase != 'stowed':
        _phi = float(x_nom_next[1])
        _v_lat_pend = float(cpas_out.pendulum_lateral_v_mps)
        _az = float(cpas_out.pendulum_azimuth_rad)
        _v_east  = float(cpas_out.wind_east_mps)  + _v_lat_pend * math.sin(_az)
        _v_north = float(cpas_out.wind_north_mps) + _v_lat_pend * math.cos(_az)
        x_nom_next[1] += dt_s * _v_north / _R_EARTH
        x_nom_next[2] += dt_s * _v_east  / (_R_EARTH * max(1e-6, math.cos(_phi)))

    x_nom = list(x_nom_next)"""

if "Apply wind drift and pendulum lateral velocity" not in src13:
    # Find the exact "x_nom = list(x_nom_next)" line and inject before it.
    target = "x_nom = list(x_nom_next)"
    assert target in src13, "cell 13: x_nom assignment line not found"
    src13 = src13.replace(target, wind_inject, 1)

# 4) Expand the row dict CPAS columns. Replace the v1 6-column block with v2.
old_row_v1 = """'force_z_from_rcs_N': float(roll_step.wrench.force_b_N[2]),
        # CPAS (parachute) state and contribution
        'cpas_phase': str(cpas_out.phase),
        'cpas_drag_cdA_m2': float(cpas_out.drag_area_cdA_m2),
        'cpas_lift_scale': float(cpas_out.lift_scale),
        'cpas_open_fraction': float(cpas_out.open_fraction),
        'cpas_force_vertical': int(bool(cpas_out.force_vertical)),
        'cpas_events': ','.join(cpas_out.events) if cpas_out.events else '',"""

old_row_plain = "'force_z_from_rcs_N': float(roll_step.wrench.force_b_N[2]),"

new_row_v2 = """'force_z_from_rcs_N': float(roll_step.wrench.force_b_N[2]),
        # CPAS (parachute) state and contribution
        'cpas_phase': str(cpas_out.phase),
        'cpas_drag_cdA_m2': float(cpas_out.drag_area_cdA_m2),
        'cpas_lift_scale': float(cpas_out.lift_scale),
        'cpas_open_fraction': float(cpas_out.open_fraction),
        'cpas_force_vertical': int(bool(cpas_out.force_vertical)),
        'cpas_events': ','.join(cpas_out.events) if cpas_out.events else '',
        'cpas_fbc_jettisoned': int(bool(cpas_out.fbc_jettisoned)),
        'cpas_mass_shed_kg': float(cpas_out.mass_shed_kg),
        'cpas_reefing_stage': int(cpas_out.reefing_stage),
        'cpas_num_mains_active': int(cpas_out.num_mains_active),
        'cpas_num_mains_squidding': int(cpas_out.num_mains_squidding),
        'cpas_pendulum_angle_deg': float(math.degrees(cpas_out.pendulum_angle_rad)),
        'cpas_pendulum_rate_deg_s': float(math.degrees(cpas_out.pendulum_rate_rad_s)),
        'cpas_pendulum_lateral_v_mps': float(cpas_out.pendulum_lateral_v_mps),
        'cpas_wind_east_mps': float(cpas_out.wind_east_mps),
        'cpas_wind_north_mps': float(cpas_out.wind_north_mps),
        'cpas_mass_kg_now': float(params['mass_kg']),"""

if old_row_v1 in src13:
    src13 = src13.replace(old_row_v1, new_row_v2)
elif "cpas_fbc_jettisoned" not in src13:
    assert old_row_plain in src13, "cell 13: row marker not found"
    src13 = src13.replace(old_row_plain, new_row_v2)

set_cell(13, src13)
print("cell 13: high-fidelity CPAS wiring applied")


# -----------------------------------------------------------------------------
# Cell 16 - save outputs
# -----------------------------------------------------------------------------
src16 = cell_src(16)

cpas_save_block = """

# --- CPAS outputs ---
cpas_events_df = pd.DataFrame(cpas_events_log) if cpas_events_log else pd.DataFrame(
    columns=['t_s', 'step', 'alt_m', 'V_mps', 'mach', 'phase', 'event',
             'reefing_stage', 'num_mains_active']
)
cpas_events_path = OUTPUT_DIR / 'cpas_events.csv'
cpas_events_df.to_csv(cpas_events_path, index=False)
print(f'  {cpas_events_path}')

cpas_summary = cpas_inst.summary()
summary['cpas_summary'] = cpas_summary
summary['cpas_events_csv'] = str(cpas_events_path)
summary['cpas_fbc_jettison_t_s'] = cpas_summary['fbc_jettison_t_s']
summary['cpas_drogue_deploy_t_s'] = cpas_summary['drogue_deploy_t_s']
summary['cpas_pilot_deploy_t_s']  = cpas_summary['pilot_deploy_t_s']
summary['cpas_main_deploy_t_s']   = cpas_summary['main_deploy_t_s']
summary['cpas_landed_t_s']        = cpas_summary['landed_t_s']
summary['cpas_mass_shed_kg_total'] = cpas_summary['mass_shed_kg_total']
summary['cpas_mains_operational_mask'] = cpas_summary['mains_operational_mask']
summary['cpas_pendulum_azimuth_deg'] = cpas_summary['pendulum_azimuth_deg']

with open(OUTPUT_DIR / 'run_summary.json', 'w', encoding='utf-8') as f:
    json.dump(summary, f, indent=2, default=str)
print(f'  {OUTPUT_DIR / "run_summary.json"} (with CPAS summary)')
"""

# Remove any previously-appended CPAS block from earlier updater runs
# (everything from the first '--- CPAS outputs ---' marker onward).
if "--- CPAS outputs ---" in src16:
    src16 = src16.split("# --- CPAS outputs ---")[0].rstrip()

src16 = src16.rstrip() + cpas_save_block
set_cell(16, src16)
print("cell 16: CPAS save block refreshed")


# -----------------------------------------------------------------------------
# Save and syntax-check
# -----------------------------------------------------------------------------
with open(NB_PATH, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1)

for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    try:
        compile("".join(c["source"]), f"cell_{i}", "exec")
    except SyntaxError as e:
        print(f"!! cell {i} syntax error: {e}")
        raise

print("ReEntryAI_run2.ipynb updated with high-fidelity CPAS wiring.")
