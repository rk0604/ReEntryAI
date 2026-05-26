import json
import re

nb_path = r"C:\Users\Risha\desktop\ReEntryAI\StochasticEntrySim\interval_ReEntry4.ipynb"
new_nb_path = r"C:\Users\Risha\desktop\ReEntryAI\StochasticEntrySim\ReEntryAI_run.ipynb"

with open(nb_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

for cell in nb["cells"]:
    if cell["cell_type"] == "code":
        source = "".join(cell["source"])
        
        # Lower terminal speed check
        source = source.replace("x_nom[3] <= 500.0", "x_nom[3] <= 1.0")
        source = source.replace("500 m/s", "1 m/s")
        
        # Modify starting conditions in x_nominal
        # The notebook might have something like constants.RADIUS_EARTH + 115000.0
        source = re.sub(r"constants\.RADIUS_EARTH\s*\+\s*[\d\_\.]+,", "constants.RADIUS_EARTH + 300_000.0,", source)
        
        # The notebook might have 10_500.0 or 6_500.0 or 7_400.0 for speed
        source = re.sub(r"\s+10\_500\.0,", "\n    7500.0,", source)
        source = re.sub(r"\s+6\_500\.0,", "\n    7500.0,", source)
        source = re.sub(r"\s+7\_400\.0,", "\n    7500.0,", source)
        source = re.sub(r"START_SPEED_MPS\s*=\s*[\d\_\.]+", "START_SPEED_MPS = 7500.0", source)
        source = re.sub(r"START_ALTITUDE_M\s*=\s*[\d\_\.]+", "START_ALTITUDE_M = 300000.0", source)

        # Re-write the cell source back
        lines = source.split("\n")
        cell["source"] = [line + "\n" for line in lines[:-1]] + ([lines[-1]] if lines[-1] else [])

# Add new cells for 3D trajectory and map trajectory
nb["cells"].append({
    "cell_type": "markdown",
    "metadata": {},
    "source": [
        "## Trajectory Visualization\n",
        "Plotting the trajectory in 3D and 2D ground track."
    ]
})

nb["cells"].append({
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": [
        "import matplotlib.pyplot as plt\n",
        "from mpl_toolkits.mplot3d import Axes3D\n",
        "import numpy as np\n",
        "\n",
        "# 3D Capsule Motion\n",
        "fig = plt.figure(figsize=(12, 10))\n",
        "ax = fig.add_subplot(111, projection='3d')\n",
        "\n",
        "r = df['nom_r_m']\n",
        "phi = df['nom_phi_rad']\n",
        "lam = df['nom_lam_rad']\n",
        "\n",
        "x = r * np.cos(phi) * np.cos(lam)\n",
        "y = r * np.cos(phi) * np.sin(lam)\n",
        "z = r * np.sin(phi)\n",
        "\n",
        "ax.plot(x, y, z, label='Capsule Trajectory', color='red', linewidth=2)\n",
        "\n",
        "# Draw Earth sphere for reference\n",
        "u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]\n",
        "ex = constants.RADIUS_EARTH * np.cos(u) * np.sin(v)\n",
        "ey = constants.RADIUS_EARTH * np.sin(u) * np.sin(v)\n",
        "ez = constants.RADIUS_EARTH * np.cos(v)\n",
        "ax.plot_surface(ex, ey, ez, color='blue', alpha=0.1, edgecolor='none')\n",
        "\n",
        "ax.set_title('Capsule Motion in 3D')\n",
        "ax.set_xlabel('X (m)')\n",
        "ax.set_ylabel('Y (m)')\n",
        "ax.set_zlabel('Z (m)')\n",
        "ax.legend()\n",
        "plt.show()\n"
    ]
})

nb["cells"].append({
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": [
        "# Ground Track (Lat/Lon)\n",
        "fig, ax = plt.subplots(figsize=(10, 6))\n",
        "ax.plot(np.degrees(df['nom_lam_rad']), np.degrees(df['nom_phi_rad']), label='Ground Track', color='blue')\n",
        "ax.scatter(np.degrees(df['nom_lam_rad'].iloc[0]), np.degrees(df['nom_phi_rad'].iloc[0]), color='green', label='Start', zorder=5)\n",
        "ax.scatter(np.degrees(df['nom_lam_rad'].iloc[-1]), np.degrees(df['nom_phi_rad'].iloc[-1]), color='red', label='End', zorder=5)\n",
        "ax.scatter(np.degrees(target_lam_rad), np.degrees(target_phi_rad), color='orange', marker='X', s=100, label='Target', zorder=5)\n",
        "\n",
        "ax.set_xlabel('Longitude (deg)')\n",
        "ax.set_ylabel('Latitude (deg)')\n",
        "ax.set_title('Trajectory Ground Track')\n",
        "ax.grid(True)\n",
        "ax.legend()\n",
        "plt.show()\n"
    ]
})

with open(new_nb_path, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1)

print("Notebook ReEntryAI_run.ipynb created successfully.")
