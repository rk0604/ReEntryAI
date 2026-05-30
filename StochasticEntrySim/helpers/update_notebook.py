import json
import re

nb_path = r"C:\Users\Risha\desktop\ReEntryAI\StochasticEntrySim\interval_ReEntry3 (1).ipynb"

with open(nb_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

for cell in nb["cells"]:
    if cell["cell_type"] != "code":
        continue

    new_source = []
    i = 0
    while i < len(cell["source"]):
        line = cell["source"][i]

        # t_final_s
        if "t_final_s =" in line or "t_final_s=" in line:
            if "num_steps = int(t_final_s" not in line:
                line = re.sub(r"t_final_s\s*=\s*[\d\.]+", "t_final_s = 15000.0", line)
        if "'t_final_s':" in line:
            line = re.sub(r"'t_final_s':\s*[\d\.]+,", "'t_final_s': 15000.0,", line)

        # max_altitude_m
        if "max_altitude_m=" in line:
            line = re.sub(r"max_altitude_m=[\d\_]+(\.0)?", "max_altitude_m=1_050_000.0", line)
        if "'max_altitude_m':" in line:
            line = re.sub(r"'max_altitude_m':\s*[\d\.]+,", "'max_altitude_m': 1050000.0,", line)

        # params dictionary setup
        if "params = {" in line:
            new_source.append(line)
            # Add params right after
            new_source.append("    \"entry_datetime_utc\": datetime.datetime(2026, 1, 1, 0, 0, 0),\n")
            new_source.append("    \"f107\": 150.0,\n")
            new_source.append("    \"f107a\": 150.0,\n")
            new_source.append("    \"ap\": 7.0,\n")
            i += 1
            # Add import datetime at the top of the cell if not present
            if not any("import datetime" in s for s in cell["source"]):
                new_source.insert(0, "import datetime\n")
            continue

        # x_nominal
        if "x_nominal =" in line and "[" in line:
            new_source.append("START_ALTITUDE_M = 900_000.0\n")
            new_source.append("START_SPEED_MPS = 7_400.0\n")
            new_source.append("START_GAMMA_DEG = -0.5\n")
            new_source.append("START_LAT_DEG = 0.0\n")
            new_source.append("START_LON_DEG = 0.0\n")
            new_source.append("START_CHI_DEG = 90.0\n\n")
            new_source.append("x_nominal = [\n")
            new_source.append("    constants.RADIUS_EARTH + START_ALTITUDE_M,\n")
            new_source.append("    math.radians(START_LAT_DEG),\n")
            new_source.append("    math.radians(START_LON_DEG),\n")
            new_source.append("    START_SPEED_MPS,\n")
            new_source.append("    math.radians(START_GAMMA_DEG),\n")
            new_source.append("    math.radians(START_CHI_DEG),\n")
            new_source.append("]\n")
            # skip until ]
            while i < len(cell["source"]) and "]" not in cell["source"][i]:
                i += 1
            i += 1
            continue

        # x_nom_diag
        if "x_nom_diag =" in line and "[" in line:
            new_source.append("x_nom_diag = list(x_nominal)\n")
            while i < len(cell["source"]) and "]" not in cell["source"][i]:
                i += 1
            i += 1
            continue

        new_source.append(line)
        i += 1

    cell["source"] = new_source

with open(nb_path, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1)
    f.write('\n')

print("Notebook updated.")
