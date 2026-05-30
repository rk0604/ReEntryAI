"""
One-shot updater: adds three high-fidelity CPAS visualization cells
to ReEntryAI_run2.ipynb right after the existing CPAS plot cells
(06a2 and 06a3). Idempotent.

New cells:
  plot_cpas_reefing_stages       -> 06a4
  plot_cpas_pendulum             -> 06a5
  plot_cpas_chute_descent_track  -> 06a6
"""
import json

NB_PATH = r"C:\Users\Risha\Desktop\ReEntryAI\StochasticEntrySim\ReEntryAI_run2.ipynb"

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


need_to_add = [
    ("plot_cpas_reefing_stages",      "plotting.plot_cpas_reefing_stages(df, save_fig)"),
    ("plot_cpas_pendulum",            "plotting.plot_cpas_pendulum(df, save_fig)"),
    ("plot_cpas_chute_descent_track", "plotting.plot_cpas_chute_descent_track(df, save_fig)"),
]

# Filter to those not already present
existing_calls = {
    name for name in (n for n, _ in need_to_add)
    if any(name in "".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code")
}
to_add = [(n, s) for n, s in need_to_add if n not in existing_calls]
print(f"Already present: {sorted(existing_calls) or '(none)'}")
print(f"Adding: {[n for n, _ in to_add] or '(nothing -- all present)'}")

if to_add:
    # Find the index of the cell containing plot_cpas_speed_dragarea -- the
    # last existing CPAS plot. New cells go directly after it.
    target_idx = None
    for i, c in enumerate(nb["cells"]):
        if c["cell_type"] == "code" and "plot_cpas_speed_dragarea" in "".join(c["source"]):
            target_idx = i
    if target_idx is None:
        raise RuntimeError("Could not find the plot_cpas_speed_dragarea cell to insert after.")

    new_cells = [code_cell(src) for _, src in to_add]
    nb["cells"] = nb["cells"][: target_idx + 1] + new_cells + nb["cells"][target_idx + 1 :]
    print(f"Inserted {len(new_cells)} cells after cell {target_idx}")

with open(NB_PATH, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1)

# Syntax check
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    try:
        compile("".join(c["source"]), f"cell_{i}", "exec")
    except SyntaxError as e:
        print(f"!! cell {i} syntax error: {e}")
        raise

print("Notebook updated with the new CPAS plot cells.")
