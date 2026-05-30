"""
One-shot updater: inserts two parachute (CPAS) visualization cells
into ReEntryAI_run2.ipynb right after the altitude-with-interval-band
plot (cell index 22).

Plot 1 -- altitude with CPAS phase shading + deploy event markers
Plot 2 -- speed and chute drag area (CD*A) on a twin axis, zoomed
          to the terminal descent
"""
import json

NB_PATH = r"C:\Users\Risha\Desktop\ReEntryAI\StochasticEntrySim\ReEntryAI_run2.ipynb"

with open(NB_PATH, "r", encoding="utf-8") as f:
    nb = json.load(f)


def code_cell(src: str) -> dict:
    lines = src.split("\n")
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [ln + "\n" for ln in lines[:-1]] + ([lines[-1]] if lines[-1] else []),
    }


def md_cell(src: str) -> dict:
    lines = src.split("\n")
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": [ln + "\n" for ln in lines[:-1]] + ([lines[-1]] if lines[-1] else []),
    }


md_intro = """### Parachute (CPAS) deployment visualization
Two views of the parachute sequence: altitude with phase-colored shading and
deploy event markers, and a zoomed terminal-descent view showing how chute
drag area drives the speed collapse."""

plot1 = """# Plot 1: Altitude with CPAS phase shading and deploy event markers
from matplotlib.patches import Patch as _Patch
import matplotlib.lines as _mlines

fig, ax = plt.subplots(figsize=(11, 4.5))

phase_colors = {
    'stowed': '#cccccc',
    'drogue': '#3a7bd5',
    'pilot':  '#9b59b6',
    'main':   '#27ae60',
    'landed': '#7f8c8d',
}

# Shade each contiguous phase window
if 'cpas_phase' in df.columns and len(df) > 0:
    phases = df['cpas_phase'].astype(str).values
    times  = df['t_s'].values
    seg_start = times[0]
    seg_phase = phases[0]
    for i in range(1, len(df)):
        if phases[i] != seg_phase:
            ax.axvspan(seg_start, times[i],
                       facecolor=phase_colors.get(seg_phase, '#dddddd'), alpha=0.22)
            seg_start = times[i]
            seg_phase = phases[i]
    ax.axvspan(seg_start, times[-1],
               facecolor=phase_colors.get(seg_phase, '#dddddd'), alpha=0.22)

# Altitude curve
ax.plot(df['t_s'], df['alt_m']/1000.0, color='black', linewidth=1.8, label='Altitude')

# Event markers
event_rows = df[df['cpas_events'].astype(str) != '']
for _, r in event_rows.iterrows():
    ax.axvline(float(r['t_s']), color='red', linestyle='--', alpha=0.65, linewidth=1.0)
    ax.annotate(str(r['cpas_events']),
                xy=(float(r['t_s']), float(r['alt_m'])/1000.0),
                xytext=(6, 10), textcoords='offset points', fontsize=8,
                color='red', rotation=12)

# Custom legend
present_phases = [p for p in ['stowed','drogue','pilot','main','landed']
                  if p in set(df['cpas_phase'].astype(str).unique())]
handles = [_Patch(facecolor=phase_colors[p], alpha=0.5, label=f'phase: {p}') for p in present_phases]
handles += [_mlines.Line2D([0],[0], color='black', lw=1.8, label='Altitude'),
            _mlines.Line2D([0],[0], color='red', ls='--', label='CPAS event')]
ax.legend(handles=handles, loc='upper right', framealpha=0.9, fontsize=9)

ax.set_xlabel('Time s')
ax.set_ylabel('Altitude km')
ax.set_title('CPAS deployment timeline overlaid on altitude')
ax.grid(True, alpha=0.4)
save_fig('06a2_cpas_altitude_phases')
plt.show()
"""

plot2 = """# Plot 2: Terminal descent -- speed (black, left) and chute CD*A (green, right)
fig, ax1 = plt.subplots(figsize=(11, 4.5))

# Zoom to a window starting a bit before the first chute deployment
deployed_mask = df['cpas_phase'].astype(str) != 'stowed'
if deployed_mask.any():
    t_first_deploy = float(df.loc[deployed_mask, 't_s'].min())
    t_zoom_start = max(0.0, t_first_deploy - 20.0)
else:
    t_zoom_start = 0.0
mask_zoom = df['t_s'] >= t_zoom_start

ax1.plot(df.loc[mask_zoom, 't_s'], df.loc[mask_zoom, 'V_mps'],
         color='black', linewidth=1.8, label='Speed')
ax1.set_xlabel('Time s')
ax1.set_ylabel('Speed m/s', color='black')
ax1.grid(True, alpha=0.4)

ax2 = ax1.twinx()
ax2.plot(df.loc[mask_zoom, 't_s'], df.loc[mask_zoom, 'cpas_drag_cdA_m2'],
         color='#27ae60', linewidth=1.8, label='Chute CD*A')
ax2.set_ylabel('Chute CD*A  m^2', color='#27ae60')
ax2.tick_params(axis='y', labelcolor='#27ae60')

# Event markers + labels on the speed axis
event_rows = df[(df['cpas_events'].astype(str) != '') & mask_zoom]
for _, r in event_rows.iterrows():
    ax1.axvline(float(r['t_s']), color='red', linestyle='--', alpha=0.65, linewidth=1.0)
    ax1.annotate(str(r['cpas_events']),
                 xy=(float(r['t_s']), float(r['V_mps'])),
                 xytext=(6, 10), textcoords='offset points', fontsize=8, color='red')

ax1.set_title('Terminal descent: speed (black) vs chute drag area (green)')
fig.tight_layout()
save_fig('06a3_cpas_speed_dragarea')
plt.show()
"""

# Find cell 22 (altitude with interval band) and insert after it.
# Be defensive: locate by source content in case indices ever shift.
insert_after = None
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    src = "".join(c["source"])
    if "06a_altitude" in src and "Altitude with interval band" in src:
        insert_after = i
        break

if insert_after is None:
    raise RuntimeError("Could not locate the altitude plot cell to insert after.")

new_cells = [md_cell(md_intro), code_cell(plot1), code_cell(plot2)]
nb["cells"] = nb["cells"][: insert_after + 1] + new_cells + nb["cells"][insert_after + 1 :]

with open(NB_PATH, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1)

print(f"Inserted {len(new_cells)} cells after cell {insert_after} (altitude plot).")
print("New cells:")
print(f"  cell {insert_after+1} [markdown] CPAS visualization intro")
print(f"  cell {insert_after+2} [code]     06a2_cpas_altitude_phases.png")
print(f"  cell {insert_after+3} [code]     06a3_cpas_speed_dragarea.png")
