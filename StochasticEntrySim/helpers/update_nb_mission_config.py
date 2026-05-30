"""
One-shot updater: wires the notebook to use mission_config.

Modifies:
  Cell 3  (imports)        -- adds `import mission_config`
                              and `MISSION_CFG = mission_config.load_config(CONFIG_PATH)`
  Cell 4  (output dir)     -- OUTPUT_DIR now comes from the config
  Cell 7  (params dict)    -- aero overrides applied from config
  Cell 11 (initial state)  -- x_nominal from config
                              target_phi_rad / target_lam_rad from config
  Cell 14 (main loop)      -- CPAS instance built from config (replacing
                              the hardcoded CPAS() call)

Idempotent.
"""
import json
import re

NB_PATH = r"C:\Users\Risha\Desktop\ReEntryAI\StochasticEntrySim\ReEntryAI_run2.ipynb"

with open(NB_PATH, "r", encoding="utf-8") as f:
    nb = json.load(f)


def cell_src(i):
    return "".join(nb["cells"][i]["source"])


def set_cell(i, new_src):
    lines = new_src.split("\n")
    nb["cells"][i]["source"] = [ln + "\n" for ln in lines[:-1]] + ([lines[-1]] if lines[-1] else [])


def find_code(predicate):
    for i, c in enumerate(nb["cells"]):
        if c["cell_type"] == "code" and predicate("".join(c["source"])):
            return i
    return None


# ----------------------------------------------------------------------
# Cell 3 -- imports + load config
# ----------------------------------------------------------------------
imp_idx = find_code(lambda s: s.startswith("# Imports") or "import constants" in s.split("\n")[:6])
if imp_idx is None:
    imp_idx = 3  # fallback
src = cell_src(imp_idx)
if "import mission_config" not in src:
    src = src.rstrip() + "\nimport mission_config\n"
if "MISSION_CFG" not in src:
    src = src.rstrip() + (
        "\n\n# Mission config -- single JSON file controls all routinely-varied\n"
        "# parameters. Change CONFIG_PATH below (or the SIM_CONFIG env var) to\n"
        "# switch runs.\n"
        "CONFIG_PATH = 'configs/default.json'\n"
        "MISSION_CFG = mission_config.load_config(CONFIG_PATH)\n"
        "print(mission_config.summarize(MISSION_CFG))\n"
    )
set_cell(imp_idx, src)
print(f"cell {imp_idx}: imports + MISSION_CFG load")


# ----------------------------------------------------------------------
# Cell 4 -- output directory
# ----------------------------------------------------------------------
out_idx = find_code(lambda s: "OUTPUT_DIR" in s and "Path" in s)
if out_idx is not None:
    new_src = (
        "# Output folder comes from the mission config so each config writes\n"
        "# to its own folder (typically revision_v1/ or runs/<run_id>/).\n"
        "from pathlib import Path\n"
        "OUTPUT_DIR = Path(str(MISSION_CFG.output_dir))\n"
        "OUTPUT_DIR.mkdir(parents=True, exist_ok=True)\n"
        "FIG_DIR = OUTPUT_DIR / 'figures'\n"
        "FIG_DIR.mkdir(parents=True, exist_ok=True)\n"
        "print('Output dir :', OUTPUT_DIR)\n"
        "print('Figure dir :', FIG_DIR)\n"
        "\n"
        "def save_fig(name):\n"
        '    """Save the current figure into the figures folder under the run\'s output_dir."""\n'
        "    path = FIG_DIR / f'{name}.png'\n"
        "    plt.savefig(path, bbox_inches='tight')\n"
        "    try:\n"
        "        rel = path.relative_to(OUTPUT_DIR.parent)\n"
        "    except ValueError:\n"
        "        rel = path\n"
        "    print(f'  saved {rel}')\n"
    )
    set_cell(out_idx, new_src)
    print(f"cell {out_idx}: output dir from MISSION_CFG.output_dir")


# ----------------------------------------------------------------------
# Cell 7 -- aero params: apply config overrides
# ----------------------------------------------------------------------
params_idx = find_code(lambda s: "'mass_kg'" in s and "'aero_model'" in s and "ref_area" in s)
if params_idx is not None:
    src = cell_src(params_idx)
    if "apply_aero_to_params" not in src:
        # Append the apply call at the end
        src = src.rstrip() + "\nmission_config.apply_aero_to_params(MISSION_CFG, params)\n"
        set_cell(params_idx, src)
        print(f"cell {params_idx}: aero overrides applied from MISSION_CFG")


# ----------------------------------------------------------------------
# Cell 11 -- initial state + targets
# ----------------------------------------------------------------------
init_idx = find_code(lambda s: "x_nominal" in s and "RADIUS_EARTH" in s and "target_phi_rad" in s)
if init_idx is not None:
    src = cell_src(init_idx)
    new_src = (
        "# Entry interface conditions and landing target -- all from the mission config.\n"
        "x_nominal = mission_config.initial_state_vector(MISSION_CFG, float(constants.RADIUS_EARTH))\n"
        "_targets = mission_config.mission_targets(MISSION_CFG)\n"
        "target_phi_rad = float(_targets['target_phi_rad'])\n"
        "target_lam_rad = float(_targets['target_lam_rad'])\n"
        "COS_GAMMA_TERMINATION_GATE = float(_targets['cos_gamma_termination_gate'])\n"
        "trajectory_id = MISSION_CFG.run_id\n"
        "\n"
        "# Build the closed-loop control state from the initial vector\n"
        "sigma_actual_0_rad = 0.0\n"
        "roll_rate_0_rad_s  = 0.0\n"
        "state_ctrl = control.ReentryState(\n"
        "    r_m=float(x_nominal[0]),\n"
        "    phi_rad=float(x_nominal[1]),\n"
        "    lam_rad=float(x_nominal[2]),\n"
        "    V_mps=float(x_nominal[3]),\n"
        "    gamma_rad=float(x_nominal[4]),\n"
        "    chi_rad=float(x_nominal[5]),\n"
        "    sigma_actual_rad=float(sigma_actual_0_rad),\n"
        "    roll_rate_rad_s=float(roll_rate_0_rad_s),\n"
        "    sigma_cmd_rad=0.0,\n"
        "    sigma_target_rad=0.0,\n"
        ")\n"
        "\n"
        "print(f'Initial: alt = {(x_nominal[0]-constants.RADIUS_EARTH)/1000.0:.1f} km, "
        "V = {x_nominal[3]:.0f} m/s, gamma = {math.degrees(x_nominal[4]):.2f} deg')\n"
        "print(f'Target : phi = {math.degrees(target_phi_rad):.2f} deg, "
        "lam = {math.degrees(target_lam_rad):.2f} deg')\n"
    )
    set_cell(init_idx, new_src)
    print(f"cell {init_idx}: initial state, targets, state_ctrl from MISSION_CFG")


# ----------------------------------------------------------------------
# Cell 14 -- main loop: replace cpas.CPAS() construction
# ----------------------------------------------------------------------
main_idx = find_code(lambda s: "cpas_inst = cpas.CPAS(" in s)
if main_idx is not None:
    src = cell_src(main_idx)
    # Ensure the random module alias is imported at the top of the cell
    if "import random as _rng_mod" not in src:
        src = "import random as _rng_mod\n" + src
    # Replace the construction call
    old_patterns = [
        "cpas_inst = cpas.CPAS(rng=_rng_mod.Random(42))",
        "cpas_inst = cpas.CPAS()",
    ]
    replaced = False
    for old in old_patterns:
        if old in src:
            src = src.replace(
                old,
                "_cpas_config = mission_config.build_cpas_config(MISSION_CFG)\n"
                "cpas_inst = cpas.CPAS(config=_cpas_config, rng=_rng_mod.Random(int(MISSION_CFG.seed)))",
            )
            replaced = True
            break
    # Also patch the cos_gamma threshold check
    src = re.sub(
        r"abs\(math\.cos\(x_nom\[4\]\)\) < 0\.2",
        "abs(math.cos(x_nom[4])) < COS_GAMMA_TERMINATION_GATE",
        src,
    )
    if replaced:
        set_cell(main_idx, src)
        print(f"cell {main_idx}: CPAS instance + cos_gamma gate from MISSION_CFG")


# ----------------------------------------------------------------------
# Save + syntax check
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
print("Notebook is wired to mission_config.")
