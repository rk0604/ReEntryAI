"""
One-shot updater: wires the CPAS parachute module into ReEntryAI_run2.ipynb.

Modifies:
  Cell 2  (imports)   -- adds `import cpas`
  Cell 13 (main loop) -- instantiates CPAS, calls step() each iteration,
                         updates params with chute drag/lift effects,
                         bypasses cos_gamma gate when chutes are deployed,
                         logs cpas_* columns
  Cell 16 (save data) -- writes cpas_events.csv and adds cpas summary
                         fields to run_summary.json
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

# 1) Instantiate CPAS + extra trackers right before the `for k in range(num_steps):`
init_marker = "rows = []"
init_inject = """rows = []
failed_heat_rows = []
thruster_fire_events = []  # one row per (step, thruster_name) fire

# --- CPAS (parachute) integration ---
cpas_inst = cpas.CPAS()
cpas_events_log = []   # one row per phase transition
params['cpas_drag_cdA_extra_m2'] = 0.0
params['cpas_lift_scale'] = 1.0
"""
old_init_block = """rows = []
failed_heat_rows = []
thruster_fire_events = []  # one row per (step, thruster_name) fire"""
assert old_init_block in src13, "cell 13: init block not found"
src13 = src13.replace(old_init_block, init_inject.rstrip())

# 2) Inject the CPAS step + termination-gate bypass at the top of the loop body.
old_gate = """for k in range(num_steps):
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

new_gate = """for k in range(num_steps):
    t_s = k * dt_s
    sigma_actual_rad = float(att_nom.sigma_rel_rad)
    roll_rate_rad_s = float(att_nom.omega_b_rad_s[2])

    # --- CPAS state machine: drogue -> pilot -> main triggers ---
    alt_now_m = alt_from_r(x_nom[0])
    T_now = atmosphere_from_state_vector(x_nom).get('T_K', float('nan'))
    mach_now = float(x_nom[3]) / math.sqrt(1.4 * 287.0 * float(T_now)) if T_now and T_now == T_now and T_now > 0 else None
    cpas_out = cpas_inst.step(t_s=t_s, alt_m=alt_now_m, V_mps=float(x_nom[3]), mach=mach_now)
    params['cpas_drag_cdA_extra_m2'] = float(cpas_out.drag_area_cdA_m2)
    params['cpas_lift_scale'] = float(cpas_out.lift_scale)
    for ev in cpas_out.events:
        cpas_events_log.append({'t_s': t_s, 'step': k, 'alt_m': alt_now_m, 'V_mps': float(x_nom[3]),
                                'mach': mach_now if mach_now is not None else float('nan'),
                                'phase': cpas_out.phase, 'event': ev})
        print(f'  CPAS {ev} at t={t_s:.1f}s  alt={alt_now_m:.1f}m  V={float(x_nom[3]):.1f}m/s')

    # cos_gamma safety gate only applies pre-deployment. Once chutes are out
    # the bank/heading channels are dead and a vertical fall is expected.
    if cpas_out.phase == 'stowed':
        if abs(math.cos(x_nom[4])) < 0.2:
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

assert old_gate in src13, "cell 13: gate block not found"
src13 = src13.replace(old_gate, new_gate)

# 3) Add CPAS fields to the logged row. Append them right after the
#    last RCS-related field (force_z_from_rcs_N) for readability.
old_row_marker = "'force_z_from_rcs_N': float(roll_step.wrench.force_b_N[2]),"
new_row_marker = """'force_z_from_rcs_N': float(roll_step.wrench.force_b_N[2]),
        # CPAS (parachute) state and contribution
        'cpas_phase': str(cpas_out.phase),
        'cpas_drag_cdA_m2': float(cpas_out.drag_area_cdA_m2),
        'cpas_lift_scale': float(cpas_out.lift_scale),
        'cpas_open_fraction': float(cpas_out.open_fraction),
        'cpas_force_vertical': int(bool(cpas_out.force_vertical)),
        'cpas_events': ','.join(cpas_out.events) if cpas_out.events else '',"""
assert old_row_marker in src13, "cell 13: row marker not found"
src13 = src13.replace(old_row_marker, new_row_marker)

set_cell(13, src13)
print("cell 13: CPAS state machine + per-step update + row columns wired in")


# -----------------------------------------------------------------------------
# Cell 16 - save outputs : add cpas_events.csv + cpas summary fields
# -----------------------------------------------------------------------------
src16 = cell_src(16)

# Append the CPAS events CSV write + extended summary to the end of cell 16.
if "cpas_events.csv" not in src16:
    cpas_save = """

# --- CPAS outputs ---
cpas_events_df = pd.DataFrame(cpas_events_log) if cpas_events_log else pd.DataFrame(
    columns=['t_s', 'step', 'alt_m', 'V_mps', 'mach', 'phase', 'event']
)
cpas_events_path = OUTPUT_DIR / 'cpas_events.csv'
cpas_events_df.to_csv(cpas_events_path, index=False)
print(f'  {cpas_events_path}')

cpas_summary = cpas_inst.summary()
summary['cpas_summary'] = cpas_summary
summary['cpas_events_csv'] = str(cpas_events_path)
summary['cpas_total_fires'] = 0  # placeholder, no chute "fires" beyond deploys
summary['cpas_drogue_deploy_t_s'] = cpas_summary['drogue_deploy_t_s']
summary['cpas_pilot_deploy_t_s'] = cpas_summary['pilot_deploy_t_s']
summary['cpas_main_deploy_t_s'] = cpas_summary['main_deploy_t_s']
summary['cpas_landed_t_s'] = cpas_summary['landed_t_s']

# Re-save run_summary.json with the CPAS fields included
with open(OUTPUT_DIR / 'run_summary.json', 'w', encoding='utf-8') as f:
    json.dump(summary, f, indent=2, default=str)
print(f'  {OUTPUT_DIR / "run_summary.json"} (with CPAS summary)')
"""
    src16 = src16.rstrip() + cpas_save
    set_cell(16, src16)
    print("cell 16: cpas_events.csv save + run_summary CPAS fields added")


# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
with open(NB_PATH, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1)

print("ReEntryAI_run2.ipynb updated with CPAS wiring.")
