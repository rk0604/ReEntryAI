"""
Rewire the notebook so every plot cell calls plotting.py instead of
holding its own copy of the figure code. Adds `import plotting` to the
imports cell and replaces the 20 plot cells with single-line function
calls.
"""
import json

NB_PATH = r"C:\Users\Risha\Desktop\ReEntryAI\StochasticEntrySim\ReEntryAI_run2.ipynb"

with open(NB_PATH, "r", encoding="utf-8") as f:
    nb = json.load(f)


def set_cell(i, new_src):
    lines = new_src.split("\n")
    nb["cells"][i]["source"] = [ln + "\n" for ln in lines[:-1]] + ([lines[-1]] if lines[-1] else [])


def cell_src(i):
    return "".join(nb["cells"][i]["source"])


# 1) imports cell: add `import plotting`
src2 = cell_src(2)
if "import plotting" not in src2:
    src2 = src2.rstrip() + "\nimport plotting\n"
    set_cell(2, src2)
    print("cell 2: added `import plotting`")

# 2) Build a context dict cell right before the plot section so every plot
#    function has access to the companions it needs. We'll piggyback on the
#    existing markdown cell `## 5. Trajectory diagnostics` (cell 17) by
#    appending a code cell after it.

# Plot calls keyed by the figure name they used to produce.
# We replace by name, not by index, so future cell shuffles don't break things.
plot_replacements = {
    "05a_ground_track":             "plotting.plot_ground_track(df, save_fig, target_phi_rad=target_phi_rad, target_lam_rad=target_lam_rad)",
    "05b_3d_descent_tube":          "plotting.plot_3d_descent_tube(df, save_fig)",
    "06a_altitude":                 "plotting.plot_altitude(df, save_fig)",
    "06b_speed":                    "plotting.plot_speed(df, save_fig)",
    "06c_gamma":                    "plotting.plot_gamma(df, save_fig)",
    "06d_chi":                      "plotting.plot_chi(df, save_fig)",
    "06e_lat_lon":                  "plotting.plot_lat_lon(df, save_fig)",
    "07a_state_widths":             "plotting.plot_state_widths(df, save_fig)",
    "07b_density_q":                "plotting.plot_density_q(df, save_fig)",
    "08a_heating_envelope":         "plotting.plot_heating_envelope(df, save_fig)",
    "08b_heat_shield_map":          "plotting.plot_heat_shield_map(df, save_fig, nominal_heat_shield=nominal_heat_shield)",
    "09a_guidance":                 "plotting.plot_guidance(df, save_fig)",
    "09b_candidate_distribution":   "plotting.plot_candidate_distribution(df, save_fig)",
    "10a_bank_error":               "plotting.plot_bank_error(df, save_fig)",
    "10b_roll_rate_accel":          "plotting.plot_roll_rate_accel(df, save_fig)",
    "10c_torque":                   "plotting.plot_torque(df, save_fig)",
    "10d_duty_vs_fired":            "plotting.plot_duty_vs_fired(df, save_fig)",
    "10e_thruster_raster":          "plotting.plot_thruster_raster(df, save_fig, thruster_fires_df=thruster_fires_df, rcs_system=rcs_system)",
    "10f_firing_rate":              "plotting.plot_firing_rate(df, save_fig, dt_s=dt_s)",
    "10g_backlog":                  "plotting.plot_backlog(df, save_fig)",
}

# Map fig name -> cell index
import re
fig_to_cell = {}
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    src = "".join(c["source"])
    m = re.search(r"save_fig\(['\"]([^'\"]+)['\"]", src)
    if m and m.group(1) in plot_replacements:
        fig_to_cell[m.group(1)] = i

replaced = 0
for fig_name, new_src in plot_replacements.items():
    if fig_name not in fig_to_cell:
        print(f"  [WARN] no cell found for {fig_name}; skipping")
        continue
    set_cell(fig_to_cell[fig_name], new_src + "\n")
    replaced += 1
print(f"Replaced {replaced} plot cells with plotting.* calls")

# Insert the CPAS plot cells right after the altitude cell. Use the name to
# find altitude, and insert directly afterwards. This restores the lost cells.
altitude_idx = fig_to_cell.get("06a_altitude")
if altitude_idx is not None:
    # don't insert duplicates if they already exist
    has_cpas_alt = any(
        c["cell_type"] == "code" and "plot_cpas_altitude_phases" in "".join(c["source"])
        for c in nb["cells"]
    )
    if not has_cpas_alt:
        def code_cell(src):
            lines = src.split("\n")
            return {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [ln + "\n" for ln in lines[:-1]] + ([lines[-1]] if lines[-1] else []),
            }
        def md_cell(src):
            lines = src.split("\n")
            return {
                "cell_type": "markdown",
                "metadata": {},
                "source": [ln + "\n" for ln in lines[:-1]] + ([lines[-1]] if lines[-1] else []),
            }
        new_cells = [
            md_cell("### Parachute (CPAS) deployment visualization\n"
                    "Altitude with phase-colored shading and deploy event markers; "
                    "and a terminal-descent zoom showing how chute drag area drives "
                    "the speed collapse."),
            code_cell("plotting.plot_cpas_altitude_phases(df, save_fig)"),
            code_cell("plotting.plot_cpas_speed_dragarea(df, save_fig)"),
        ]
        nb["cells"] = nb["cells"][: altitude_idx + 1] + new_cells + nb["cells"][altitude_idx + 1 :]
        print(f"Inserted CPAS plot cells after cell {altitude_idx} (06a_altitude)")

with open(NB_PATH, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1)

# Verify syntax
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    try:
        compile("".join(c["source"]), f"cell_{i}", "exec")
    except SyntaxError as e:
        print(f"!! cell {i} syntax error: {e}")
        raise
print("All cells compile OK")
