"""
3D view of randomly chosen training flights, to eyeball dataset variation.

    python jepa/plot_trajectories_3d.py [episodes.pkl] [seed]    (from StochasticEntrySim/)

Writes jepa/figures/trajectories_3d.png (two panels of 10 flights each) and
prints how spread out the whole dataset is.
"""
import os
import pickle
import sys

import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)


def plot_3d(eps, seed=7, out_path=None, n_per_panel=10):
    """Two 3D panels (longitude, latitude, altitude) of randomly sampled flights."""
    n_per_panel = max(1, min(n_per_panel, len(eps) // 2))
    pick =np.random.default_rng(seed).choice(len(eps), 2 * n_per_panel, replace=False)
    cmap = plt.get_cmap("tab10")
    fig = plt.figure(figsize=(15, 6.5))
    for p in range(2):
        ax = fig.add_subplot(1, 2, p + 1, projection="3d")
        for j, i in enumerate(pick[p * n_per_panel:(p + 1) * n_per_panel]):
            e, y = eps[i], eps[i]["y"]
            ic, d = e["ic"], e["domain"]
            lab = (f'{ic["altitude_m"] / 1e3:.0f} km, {ic["velocity_mps"] / 1e3:.2f} km/s, '
                   f'{ic["gamma_deg"]:.1f}°, ρ×{d["density_scale"]:.2f}'
                   f'{"" if e["outcome"] == "drogue" else " (" + e["outcome"] + ")"}')
            ax.plot(np.degrees(y[:, 2]), np.degrees(y[:, 1]), y[:, 0] / 1000, color=cmap(j % 10), lw=1.4, label=lab)
            ax.scatter(np.degrees(y[-1, 2]), np.degrees(y[-1, 1]), y[-1, 0] / 1000, color=cmap(j % 10), s=14)
        ax.set(xlabel="longitude [deg]", ylabel="latitude [deg]", zlabel="altitude [km]",
               title=f"Flights {p * n_per_panel + 1}–{(p + 1) * n_per_panel} (random sample)")
        ax.view_init(elev=22, azim=-60)
        ax.legend(fontsize=6.5, loc="upper left", frameon=False, title="entry alt, speed, γ, density",
                  title_fontsize=7)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
    return fig


def spread_report(eps):
    """Print how varied the whole dataset is."""
    start = np.array([[np.degrees(e["y"][0, 2]), np.degrees(e["y"][0, 1])] for e in eps])
    end = np.array([[np.degrees(e["y"][-1, 2]), np.degrees(e["y"][-1, 1])] for e in eps])
    rng_km = np.hypot(*(end - start).T) * 111
    t = np.array([e["t_final_s"] for e in eps])
    print(f"{len(eps)} flights | flight time {t.min():.0f}..{t.max():.0f} s (mean {t.mean():.0f})")
    print(f"downrange flown: {rng_km.min():.0f}..{rng_km.max():.0f} km (median {np.median(rng_km):.0f})")
    print(f"endpoint lon {end[:, 0].min():.1f}..{end[:, 0].max():.1f} deg, lat {end[:, 1].min():.1f}..{end[:, 1].max():.1f} deg")


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(SIM, "jepa_runs", "episodes_900.pkl")
    eps = pickle.load(open(path, "rb"))
    out = os.path.join(HERE, "figures", "trajectories_3d.png")
    plot_3d(eps, seed=int(sys.argv[2]) if len(sys.argv) > 2 else 7, out_path=out)
    print("wrote", out)
    spread_report(eps)
