import json
import re

nb_path = r"C:\Users\Risha\desktop\ReEntryAI\StochasticEntrySim\interval_ReEntry4.ipynb"

with open(nb_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

for cell in nb["cells"]:
    if cell["cell_type"] != "code":
        continue
        
    source = "".join(cell["source"])

    # Update params
    if "params = {" in source and '"CD_const": 1.15' in source:
        if '"aero_model"' not in source:
            source = source.replace('"CD_const": 1.15,', '"aero_model": "scheduled_orion_like",\n    "CD_const": 1.15,')

    # Add initial print block
    if "x_nominal =" in source and "START_CHI_DEG" in source:
        if "coeffs0 =" not in source:
            addition = """
# Print the current aero coefficients at the starting state
import math_3d
coeffs0 = math_3d.aero_coefficients_from_speed_altitude(
    V_mps=x_nominal[3],
    altitude_m=max(0.0, constants.altitude(x_nominal[0])),
    params=params,
)
print(f"Initial Aero Coefficients: {coeffs0}")
"""
            source = source + addition

    # Add logging block
    if "trajectory_log.append(" in source and "nom_r_m" in source:
        # Just add a comment for where to add logging if it's complex
        if "# TODO: Add CD, CL, LD, aero_model, aero_schedule_region, V_kms to trajectory_log" not in source:
            source = source.replace("trajectory_log.append(", "# TODO: Add CD, CL, LD, aero_model, aero_schedule_region, V_kms to trajectory_log using nominal_step['aero'] dictionary\ntrajectory_log.append(")

    # Write back cell source
    lines = source.split("\n")
    cell["source"] = [line + "\n" for line in lines[:-1]] + ([lines[-1]] if lines[-1] else [])

with open(nb_path, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1)

print("Notebook updated.")
