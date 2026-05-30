"""
One-shot updater: inserts a PLOT_CATEGORIES config cell into
ReEntryAI_run2.ipynb and wraps each plot cell with a category guard so
the user can render only the categories they care about for the current
run. Idempotent.

Usage in the notebook:
    PLOT_CATEGORIES = None             # default -- render everything
    PLOT_CATEGORIES = {'cpas'}         # only CPAS plots
    PLOT_CATEGORIES = {'cpas','heating'}  # CPAS + heating
"""
import json
import re

NB_PATH = r"C:\Users\Risha\Desktop\ReEntryAI\StochasticEntrySim\ReEntryAI_run2.ipynb"

# Mirror of the CATEGORIES dict in plotting.py. Kept here so this script
# doesn't import plotting (which pulls in matplotlib etc.).
CATEGORY_FOR_FN = {
    "plot_ground_track":              "trajectory",
    "plot_3d_descent_tube":           "trajectory",
    "plot_altitude":                  "state",
    "plot_speed":                     "state",
    "plot_gamma":                     "state",
    "plot_chi":                       "state",
    "plot_lat_lon":                   "state",
    "plot_cpas_altitude_phases":      "cpas",
    "plot_cpas_speed_dragarea":       "cpas",
    "plot_cpas_reefing_stages":       "cpas",
    "plot_cpas_pendulum":             "cpas",
    "plot_cpas_chute_descent_track":  "cpas",
    "plot_cpas_squidding":            "cpas",
    "plot_state_widths":              "interval",
    "plot_density_q":                 "interval",
    "plot_heating_envelope":          "heating",
    "plot_heat_shield_map":           "heating",
    "plot_guidance":                  "guidance",
    "plot_candidate_distribution":    "guidance",
    "plot_bank_error":                "rcs",
    "plot_roll_rate_accel":           "rcs",
    "plot_torque":                    "rcs",
    "plot_duty_vs_fired":             "rcs",
    "plot_thruster_raster":           "rcs",
    "plot_firing_rate":               "rcs",
    "plot_backlog":                   "rcs",
}

CONFIG_CELL_MARKER = "# === PLOT_CATEGORIES config ==="
CONFIG_CELL_SRC = """# === PLOT_CATEGORIES config ===
# Controls which plots run when "Run All" executes the notebook.
# Set to None (default) to render every category, or pass a set of
# category names to filter. Available categories:
#   trajectory  -- ground track, 3D descent tube
#   state       -- altitude / speed / gamma / chi / lat-lon histories
#   cpas        -- all parachute plots
#   interval    -- supervisor diagnostics (state widths, density/q)
#   heating     -- heating envelope + heat shield map
#   guidance    -- predictor-corrector candidate analysis
#   rcs         -- bank tracking, roll, torque, raster, firing rate
#
# Examples:
#   PLOT_CATEGORIES = None                 # render everything
#   PLOT_CATEGORIES = {'cpas'}             # just parachute plots
#   PLOT_CATEGORIES = {'cpas','heating'}   # parachute + heating
PLOT_CATEGORIES = None
"""


with open(NB_PATH, "r", encoding="utf-8") as f:
    nb = json.load(f)


def code_cell(src):
    lines = src.split("\n")
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [ln + "\n" for ln in lines[:-1]] + ([lines[-1]] if lines[-1] else []),
    }


# ----------------------------------------------------------------------
# 1. Insert / refresh the PLOT_CATEGORIES config cell right after the
#    markdown header for the diagnostics section.
# ----------------------------------------------------------------------
existing_config_idx = next(
    (i for i, c in enumerate(nb["cells"])
     if c["cell_type"] == "code" and CONFIG_CELL_MARKER in "".join(c["source"])),
    None,
)
if existing_config_idx is not None:
    nb["cells"][existing_config_idx] = code_cell(CONFIG_CELL_SRC)
    print(f"cell {existing_config_idx}: refreshed PLOT_CATEGORIES config cell")
else:
    header_idx = next(
        (i for i, c in enumerate(nb["cells"])
         if c["cell_type"] == "markdown" and "Trajectory diagnostics" in "".join(c["source"])),
        None,
    )
    if header_idx is None:
        raise RuntimeError("Could not find the 'Trajectory diagnostics' header to insert after.")
    nb["cells"].insert(header_idx + 1, code_cell(CONFIG_CELL_SRC))
    print(f"Inserted PLOT_CATEGORIES config cell after markdown cell {header_idx}")

# ----------------------------------------------------------------------
# 2. Wrap each plotting.plot_*(...) cell with a category guard.
#    Idempotent: cells already wrapped are left alone.
# ----------------------------------------------------------------------
CALL_RE = re.compile(r"plotting\.(plot_[A-Za-z0-9_]+)\s*\(")

wrapped_count = 0
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    src = "".join(c["source"])
    if "plotting.want(" in src:
        continue  # already wrapped
    m = CALL_RE.search(src)
    if not m:
        continue
    fn_name = m.group(1)
    category = CATEGORY_FOR_FN.get(fn_name)
    if category is None:
        print(f"  [WARN] cell {i}: unknown plot fn {fn_name!r}, leaving untouched")
        continue

    # Build a guarded version of the cell. Indent every existing line by
    # four spaces and prepend the if-statement.
    indented = "\n".join("    " + ln if ln else "" for ln in src.splitlines())
    new_src = f"if plotting.want({category!r}, PLOT_CATEGORIES):\n{indented}\n"

    lines = new_src.split("\n")
    c["source"] = [ln + "\n" for ln in lines[:-1]] + ([lines[-1]] if lines[-1] else [])
    wrapped_count += 1

print(f"Wrapped {wrapped_count} plot cells with category guard")

# ----------------------------------------------------------------------
# Save + syntax-check
# ----------------------------------------------------------------------
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
print("Notebook ready with modular plot filtering.")
