"""
One-shot, idempotent updater: adds the shared telemetry.py derived metrics,
RCS fuel model, and RL reward ingredients to ReEntryAI_run2.ipynb so the
notebook driver logs the same columns as run_sim.py.
"""
import json

NB = r"C:\Users\Risha\Desktop\ReEntryAI\StochasticEntrySim\ReEntryAI_run2.ipynb"
with open(NB, encoding="utf-8") as f:
    nb = json.load(f)


def find(pred):
    for i, c in enumerate(nb["cells"]):
        if c["cell_type"] == "code" and pred("".join(c["source"])):
            return i
    return None


def setcell(i, s):
    lines = s.split("\n")
    nb["cells"][i]["source"] = [ln + "\n" for ln in lines[:-1]] + ([lines[-1]] if lines[-1] else [])


# 1) import telemetry in the imports cell
imp = find(lambda s: "import plotting" in s or "import cpas" in s)
if imp is not None:
    s = "".join(nb["cells"][imp]["source"])
    if "import telemetry" not in s:
        s = s.rstrip() + "\nimport telemetry\n"
        setcell(imp, s)
        print(f"cell {imp}: import telemetry")

# 2) main loop cell
m = find(lambda s: "for k in range(num_steps)" in s)
s = "".join(nb["cells"][m]["source"])

# 2a) cumulative fuel init right after rows = []
if "rcs_fuel_used_kg_total" not in s:
    s = s.replace("rows = []",
                  "rows = []\nrcs_fuel_used_kg_total = 0.0", 1)

# 2b) derived/fuel/reward block right after the atm assignment
anchor = "atm = atmosphere_from_state_vector(x_nom)"
block = anchor + """

    # --- Shared derived telemetry (g-load, energy, Mach, FMV) ---
    derived = telemetry.compute_derived_metrics(
        r_m=float(x_nom[0]), V_mps=float(x_nom[3]), altitude_m=float(atm['alt_m']),
        T_K=atm.get('T_K', None), rho_kgm3=float(atm['rho_kgm3']),
        drag_mag_N=float(nominal_aero['drag_mag_N']), lift_mag_N=float(nominal_aero['lift_mag_N']),
        mass_kg=float(params['mass_kg']), fmv=float(nominal_aero.get('FMV', float('nan'))),
    )
    fired_on_time_s = float(roll_step.num_fired_internal_steps) * float(constants.RCS_INTERNAL_DT_S)
    fired_thruster_seconds = float(len(roll_step.fire_cmd.active_names())) * fired_on_time_s
    rcs_fuel_step_kg = telemetry.rcs_propellant_kg(fired_thruster_seconds)
    rcs_fuel_rate_kg_s = rcs_fuel_step_kg / float(dt_s) if dt_s > 0 else 0.0
    rcs_fuel_used_kg_total += rcs_fuel_step_kg
    reward = telemetry.composite_reward(
        range_to_go_m=float(ctrl_out.obs.get('range_to_go_m', np.nan)),
        qdot_w_m2=float(nominal_heat_qdot.hi), load_factor_g=float(derived['load_factor_g']),
        rcs_fuel_rate_kg_s=float(rcs_fuel_rate_kg_s),
    )"""
if "Shared derived telemetry" not in s:
    assert anchor in s, "atm anchor not found in notebook loop"
    s = s.replace(anchor, block, 1)

# 2c) row columns right after the LD line
ld_line = "'LD': float(nominal_aero['LD']),"
cols = ld_line + """
        'mach': float(derived['mach']),
        'speed_of_sound_mps': float(derived['speed_of_sound_mps']),
        'fmv': float(derived['fmv']),
        'specific_kinetic_e_j_kg': float(derived['specific_kinetic_e_j_kg']),
        'specific_potential_e_j_kg': float(derived['specific_potential_e_j_kg']),
        'specific_energy_j_kg': float(derived['specific_energy_j_kg']),
        'energy_height_m': float(derived['energy_height_m']),
        'aero_force_N': float(derived['aero_force_N']),
        'aero_decel_mps2': float(derived['aero_decel_mps2']),
        'drag_decel_mps2': float(derived['drag_decel_mps2']),
        'load_factor_g': float(derived['load_factor_g']),
        'rcs_fuel_step_kg': float(rcs_fuel_step_kg),
        'rcs_fuel_rate_kg_s': float(rcs_fuel_rate_kg_s),
        'rcs_fuel_used_kg': float(rcs_fuel_used_kg_total),
        'rew_range_term': float(reward['rew_range_term']),
        'rew_heat_term': float(reward['rew_heat_term']),
        'rew_gload_term': float(reward['rew_gload_term']),
        'rew_fuel_term': float(reward['rew_fuel_term']),
        'reward_default': float(reward['reward_default']),"""
if "'load_factor_g':" not in s:
    assert ld_line in s, "LD row line not found"
    s = s.replace(ld_line, cols, 1)

setcell(m, s)
print(f"cell {m}: telemetry computation + {19} new row columns")

with open(NB, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1)

# syntax check
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    try:
        compile("".join(c["source"]), f"cell_{i}", "exec")
    except SyntaxError as e:
        print(f"!! cell {i} syntax error: {e}")
        raise
print("Notebook telemetry wiring complete.")
