import json

nb_path = r"C:\Users\Risha\desktop\ReEntryAI\StochasticEntrySim\interval_ReEntry4.ipynb"

with open(nb_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

for idx, cell in enumerate(nb["cells"]):
    if cell["cell_type"] == "code":
        source = "".join(cell["source"])
        if "500.0" in source and "speed is too low" in source:
            print(f"Found speed check in cell {idx}")
        if "x_nominal =" in source or "START_ALTITUDE_M =" in source or "Initial altitude" in source:
            print(f"Found starting conditions in cell {idx}")
