"""
Build ReEntryAI_JEPA_data.ipynb: collect the 6,000-flight wide-envelope dataset on Colab.

All logic lives in jepa/data.py; the notebook configures, collects in
resumable shards (saved to Google Drive), and summarizes. Edit this builder
and re-run it rather than hand-editing the .ipynb:

    python jepa/build_data_notebook.py        (from StochasticEntrySim/)
"""
import os
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "ReEntryAI_JEPA_data.ipynb")

cells = []


def md(text):
    """Append a markdown cell."""
    cells.append(nbf.v4.new_markdown_cell(text.strip("\n")))


def code(text):
    """Append a code cell."""
    cells.append(nbf.v4.new_code_cell(text.strip("\n")))


md(r"""
# SpaceJEPA — wide-envelope training data (6,000 flights)

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/rk0604/ReEntryAI/blob/jepa-variant/StochasticEntrySim/ReEntryAI_JEPA_data.ipynb)

The first dataset (320 flights) all started in one corridor: entry within ~1 km,
~100 m/s, ~0.4° and ~3° of heading of a single nominal state, with only the
random steering making them differ. This notebook collects **6,000 flights**
across a wide entry envelope and a wider atmosphere, so the world model sees
many kinds of reentry, not one.

| | First dataset (320) | This dataset (6,000) |
|---|---|---|
| Entry altitude | 120 km ± 1 km | **100 – 140 km** |
| Entry speed | 7.8 km/s ± 0.1 | **6.8 – 8.6 km/s** |
| Flight-path angle | −2.5° ± 0.4° | **−1° to −6°**, steeper required at high speed (no skip-outs) |
| Heading | 46° ± 3° | **any (0 – 360°)** |
| Start latitude | 0° ± 0.3° | **−30° to +30°** |
| Air density / L/D | ×0.85 – 1.15 | **×0.70 – 1.30** |
| Mass | ×0.95 – 1.05 | **×0.90 – 1.10** |
| Steering | random GP + reversal profiles | same |

Every flight still goes through the full simulator (bank actuator, roll PD,
12-jet RCS, 3DOF, NRLMSISE table).

**Cost.** About 8 s of single-core simulation per flight, ~18 CPU-hours in
total. Use **Runtime → Change runtime type → High-RAM** (Colab Pro), which
gives more cores, for **~3.5–4.5 hours**; the standard 2-core runtime takes 12+ hours.

**Resumable.** Flights are saved to Google Drive in shards of 250. If Colab
disconnects, re-run the cells top to bottom: finished shards are skipped, and
flight *i* always uses the same seed, so the dataset is identical either way.
""")

md("## 1. Setup (clones the repo and mounts Drive on Colab)")
code(r"""
import os, sys
IN_COLAB = "google.colab" in sys.modules
if IN_COLAB:
    REPO = "/content/ReEntryAI"
    if not os.path.exists(REPO):
        !git clone -q -b jepa-variant https://github.com/rk0604/ReEntryAI.git {REPO}
    else:
        !git -C {REPO} pull -q
    %pip install -q pymsis
    os.chdir(f"{REPO}/StochasticEntrySim")
    from google.colab import drive
    drive.mount("/content/drive")
    OUT_DIR = "/content/drive/MyDrive/ReEntryAI_jepa_data/wide_6000"
else:
    OUT_DIR = "jepa_runs/wide_6000"          # local run: next to the other JEPA data
sys.path.insert(0, os.getcwd())
os.makedirs(OUT_DIR, exist_ok=True)
print("working dir:", os.getcwd())
print("shards go to:", OUT_DIR)
""")

md("## 2. Configuration")
code(r"""
import time, pickle, collections
from dataclasses import asdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
%matplotlib inline
from jepa import data as D
from jepa.plot_trajectories_3d import plot_3d, spread_report

N_TOTAL = 6000
SHARD_SIZE = 250            # flights per saved file; a disconnect loses at most one shard
MAX_STEPS = 1200            # cap on flight length in seconds (long glides need ~900+)
WORKERS = os.cpu_count() or 1
ENVELOPE = D.EnvelopeSpec() # entry-condition ranges (edit fields here to change them)
DOMAIN = D.WIDE_DOMAIN      # atmosphere / mass randomization ranges

print(f"cores: {WORKERS}")
display(pd.Series({**{f"entry {k}": v for k, v in asdict(ENVELOPE).items()},
                   **{f"truth {k}": v for k, v in asdict(DOMAIN).items()}}, name="range"))
est_h = N_TOTAL * 11.0 / (0.7 * WORKERS) / 3600   # ~11 s per flight on a Colab core, ~70% parallel efficiency
print(f"rough estimate: {est_h:.1f} h for {N_TOTAL} flights on {WORKERS} cores (refined after the smoke test)")
""")

md(r"""
## 3. Preview the entry envelope (no simulation)

2,000 draws of the starting condition. The empty corner in the speed–angle plot
is the skip-out guard: fast, shallow entries bounce back to space, so the
sampler requires steeper angles as speed rises.
""")
code(r"""
rng = np.random.default_rng(0)
ics = pd.DataFrame([D.sample_entry_ic(rng, ENVELOPE) for _ in range(2000)])
fig, axes = plt.subplots(1, 3, figsize=(13, 3.4))
axes[0].scatter(ics.velocity_mps / 1000, ics.gamma_deg, s=3, color="#2a78d6")
axes[0].set(xlabel="entry speed [km/s]", ylabel="flight-path angle [deg]", title="speed vs angle")
axes[1].scatter(ics.velocity_mps / 1000, ics.altitude_m / 1000, s=3, color="#2a78d6")
axes[1].set(xlabel="entry speed [km/s]", ylabel="entry altitude [km]", title="speed vs altitude")
axes[2].hist(ics.chi_deg, bins=36, color="#2a78d6")
axes[2].set(xlabel="heading [deg]", ylabel="count", title="heading")
for ax in axes: ax.grid(alpha=0.3)
fig.tight_layout()
""")

md(r"""
## 4. Smoke test (about a minute)

Flies a handful of flights into a separate folder, to confirm the runtime works
and measure the real cost per flight before committing hours.
""")
code(r"""
smoke_dir = os.path.join(OUT_DIR, "_smoke")
n_smoke = max(8, WORKERS)
t0 = time.time()
D.collect_shards(n_smoke, smoke_dir, shard_size=n_smoke, seed0=900_000, workers=WORKERS,
                 max_steps=MAX_STEPS, spec=DOMAIN, envelope=ENVELOPE)
wall = time.time() - t0
smoke = D.load_shards(smoke_dir)
print(collections.Counter(e["outcome"] for e in smoke))
per_round = wall / max(1, -(-n_smoke // WORKERS))
print(f"{len(smoke)} flights in {wall:.0f}s  ->  estimate for {N_TOTAL}: "
      f"~{N_TOTAL / WORKERS * per_round / 3600:.1f} h")
""")

md(r"""
## 5. Collect all 6,000 flights

This is the long cell. Progress prints per shard with a time-remaining estimate.
**If Colab disconnects:** reconnect, run cells 1, 2, and this one again. Shards
already on Drive are skipped.
""")
code(r"""
D.collect_shards(N_TOTAL, OUT_DIR, shard_size=SHARD_SIZE, seed0=100_000, workers=WORKERS,
                 max_steps=MAX_STEPS, spec=DOMAIN, envelope=ENVELOPE)
""")

md("## 6. What did we get?")
code(r"""
eps = D.load_shards(OUT_DIR)
df = pd.DataFrame([{**e["ic"], **e["domain"], "kind": e["kind"], "outcome": e["outcome"],
                    "flight_s": e["t_final_s"],
                    "peak_g": float(np.max(np.hypot(e["y"][:, 8], e["y"][:, 9])) / 9.80665)} for e in eps])
n_steps = sum(len(e["cmd"]) for e in eps)
print(f"{len(eps)} flights | {n_steps:,} one-second steps = {n_steps / 3600:.0f} h of flight | "
      f"{n_steps // 5:,} five-second training steps")
display(df.groupby(["kind", "outcome"]).size().unstack(fill_value=0))
display(df[["altitude_m", "velocity_mps", "gamma_deg", "chi_deg", "phi_deg", "density_scale",
            "lift_scale", "mass_scale", "flight_s", "peak_g"]].describe().loc[["min", "50%", "max"]].round(2))
spread_report(eps)
""")
code(r"""
# State coverage: the new dataset against the first one (if it is available).
fig, ax = plt.subplots(figsize=(7, 4.2))
old_path = "jepa_runs/episodes_900.pkl"
if os.path.exists(old_path):
    old = pickle.load(open(old_path, "rb"))
    Yo = np.concatenate([e["y"][::5] for e in old])
    ax.scatter(Yo[:, 3] / 1000, Yo[:, 0] / 1000, s=1, color="#eb6834", alpha=0.5, label=f"first dataset ({len(old)})")
Y = np.concatenate([e["y"][::10] for e in eps[::3]])
ax.scatter(Y[:, 3] / 1000, Y[:, 0] / 1000, s=1, color="#2a78d6", alpha=0.25, label=f"wide dataset ({len(eps)}, subsampled)")
ax.set(xlabel="speed [km/s]", ylabel="altitude [km]", title="State coverage")
ax.grid(alpha=0.3)
ax.legend(markerscale=8, frameon=False)
fig.tight_layout()
""")
code(r"""
plot_3d(eps, seed=7, out_path=os.path.join(OUT_DIR, "trajectories_3d.png"))
plt.show()
""")

md(r"""
## 7. Merge into one file for training

Writes `episodes_wide_6000.pkl` (~100 MB) next to the shards on Drive.

To train on it in `ReEntryAI_JEPA.ipynb`, replace the `episodes = D.collect_dataset(...)`
line in section 1 with:

```python
episodes = pickle.load(open("/content/drive/MyDrive/ReEntryAI_jepa_data/wide_6000/episodes_wide_6000.pkl", "rb"))
```

and raise `steps_jepa` / `steps_pi` (about 20x more data than before: ~2x the
steps is a sensible first try). The frozen test set and the classical
comparison are unchanged, so the results stay comparable with the first runs.
""")
code(r"""
merged = os.path.join(OUT_DIR, f"episodes_wide_{len(eps)}.pkl")
with open(merged, "wb") as f:
    pickle.dump(eps, f)
print(f"wrote {merged} ({os.path.getsize(merged) / 1e6:.0f} MB)")
""")

nb = nbf.v4.new_notebook()
nb["cells"] = cells
nb["metadata"] = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                  "language_info": {"name": "python"}}
nbf.write(nb, OUT)
print(f"wrote {OUT} ({len(cells)} cells)")
