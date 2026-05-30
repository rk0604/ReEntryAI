import json

nb_path = r"C:\Users\Risha\desktop\ReEntryAI\StochasticEntrySim\old_code\interval_ReEntry4.ipynb"

with open(nb_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

code_lines = []
for cell in nb["cells"]:
    if cell["cell_type"] == "code":
        source = cell["source"]
        if source:
            code_lines.extend(source)
            if not source[-1].endswith("\n"):
                code_lines.append("\n")
            code_lines.append("\n") # Ensure space between cells to avoid indentation issues

with open(r"C:\Users\Risha\desktop\ReEntryAI\StochasticEntrySim\run_sim.py", "w", encoding="utf-8") as f:
    f.writelines(code_lines)

print("Extraction complete.")
