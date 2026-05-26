"""
Centralized plotting for the ReEntryAI simulation.

Both ReEntryAI_run2.ipynb and run_sim.py call into this module so the
two parallel drivers always produce the same set of named PNG figures.

Public entry point:
    render_all_figures(df, save_fig, **ctx)
        -- df is the per-step trajectory DataFrame
        -- save_fig(name) saves the current matplotlib figure
        -- ctx kwargs supply optional companions:
             target_phi_rad, target_lam_rad
             thruster_fires_df, rcs_system, nominal_heat_shield, dt_s

Each individual plot is also exposed as plot_<name>(df, save_fig, **ctx)
so a notebook cell can drive a single figure if desired.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  needed for projection='3d'

import constants


# =============================================================================
# helpers
# =============================================================================

def _normalize_columns(df):
    """
    Add aliased columns so the rest of this module can use a single
    naming convention regardless of which driver built the dataframe.

    Notebook columns:     V_lo, gamma_lo, chi_lo, V_width, ...
                          tau_roll_realized_Nm, force_x_from_rcs_N
    run_sim.py columns:   V_mps_lo, gamma_rad_lo, chi_rad_lo, V_mps_width, ...
                          torque_z_from_rcs, force_x_from_rcs
    """
    if "tau_roll_realized_Nm" not in df.columns and "torque_z_from_rcs" in df.columns:
        df["tau_roll_realized_Nm"] = df["torque_z_from_rcs"]
    for axis in ("x", "y", "z"):
        col_new = f"force_{axis}_from_rcs_N"
        col_old = f"force_{axis}_from_rcs"
        if col_new not in df.columns and col_old in df.columns:
            df[col_new] = df[col_old]

    # State interval bands: notebook bare-name -> run_sim.py suffixed-name
    suffixed = {
        "r":     "r_m",
        "phi":   "phi_rad",
        "lam":   "lam_rad",
        "V":     "V_mps",
        "gamma": "gamma_rad",
        "chi":   "chi_rad",
    }
    for short, long_ in suffixed.items():
        for suffix in ("lo", "hi", "width"):
            short_col = f"{short}_{suffix}"
            long_col = f"{long_}_{suffix}"
            if short_col not in df.columns and long_col in df.columns:
                df[short_col] = df[long_col]
    return df


def plot_with_band(ax, t, nom, lo, hi, title, ylabel,
                   nom_label="Nominal", lo_label="Interval low", hi_label="Interval high"):
    """Shared helper used by the band-style state-history plots."""
    ax.plot(t, nom, label=nom_label)
    ax.plot(t, lo, label=lo_label)
    ax.plot(t, hi, label=hi_label)
    ax.set_xlabel("Time s")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True)
    ax.legend()


# =============================================================================
# individual plots
# =============================================================================

def plot_ground_track(df, save_fig, target_phi_rad=0.0, target_lam_rad=0.0, **_):
    track = df.copy()
    track["phi_deg"] = np.degrees(track["phi_rad"])
    track["lam_deg"] = np.degrees(track["lam_rad"])

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(track["lam_deg"], track["phi_deg"], lw=2, label="Nominal ground track", color="C0")
    ax.scatter(track["lam_deg"].iloc[0], track["phi_deg"].iloc[0], s=120, marker="o", color="C0", label="Start")
    ax.scatter(track["lam_deg"].iloc[-1], track["phi_deg"].iloc[-1], s=140, marker="s", color="C1", label="End")
    ax.scatter(math.degrees(target_lam_rad), math.degrees(target_phi_rad),
               s=200, marker="x", color="C2", linewidths=3, label="Target")
    ax.set_xlabel("Longitude deg")
    ax.set_ylabel("Latitude deg")
    ax.set_title("Ground track with start, end, and target")
    ax.legend()
    ax.set_aspect("equal", "datalim")
    save_fig("05a_ground_track")
    plt.show()


def plot_3d_descent_tube(df, save_fig, **_):
    track = df.copy()
    phi0 = float(track["phi_rad"].iloc[0])
    lam0 = float(track["lam_rad"].iloc[0])
    R = float(constants.RADIUS_EARTH)
    track["north_m"] = R * (track["phi_rad"] - phi0)
    track["east_m"] = R * (track["lam_rad"] - lam0) * math.cos(phi0)

    has_interval = track["interval_alt_lo"].notna() if "interval_alt_lo" in track.columns else None
    tr_iv = track[has_interval] if has_interval is not None else track.iloc[0:0]

    fig = plt.figure(figsize=(16, 6))
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    for ax, elev, azim, title in [(ax1, 25, -60, "3D descent tube view 1"),
                                  (ax2, 25, 120, "3D descent tube view 2")]:
        ax.plot(track["north_m"], track["east_m"], track["alt_m"],
                color="C0", lw=2, label="Nominal trajectory")
        if len(tr_iv):
            ax.plot(tr_iv["north_m"], tr_iv["east_m"], tr_iv["interval_alt_lo"],
                    color="C1", ls="--", lw=1, label="Interval low corner")
            ax.plot(tr_iv["north_m"], tr_iv["east_m"], tr_iv["interval_alt_hi"],
                    color="C2", ls="--", lw=1, label="Interval high corner")
        ax.scatter(track["north_m"].iloc[0], track["east_m"].iloc[0], track["alt_m"].iloc[0],
                   s=80, marker="o", color="C0", label="Start")
        ax.scatter(track["north_m"].iloc[-1], track["east_m"].iloc[-1], track["alt_m"].iloc[-1],
                   s=80, marker="s", color="C1", label="End")
        ax.set_xlabel("North m"); ax.set_ylabel("East m"); ax.set_zlabel("Altitude m")
        ax.set_title(title)
        ax.view_init(elev=elev, azim=azim)
        ax.legend(loc="upper left", fontsize=8)
    save_fig("05b_3d_descent_tube")
    plt.show()


def plot_altitude(df, save_fig, **_):
    fig, ax = plt.subplots(figsize=(11, 4))
    plot_with_band(ax, df["t_s"], df["alt_m"] / 1000.0,
                   df["interval_alt_lo"] / 1000.0, df["interval_alt_hi"] / 1000.0,
                   "Altitude with interval band", "Altitude km")
    save_fig("06a_altitude")
    plt.show()


def plot_cpas_altitude_phases(df, save_fig, **_):
    """Altitude curve with CPAS phase shading and deploy event markers."""
    if "cpas_phase" not in df.columns:
        return
    import matplotlib.lines as _mlines

    fig, ax = plt.subplots(figsize=(11, 4.5))
    phase_colors = {
        "stowed": "#cccccc",
        "drogue": "#3a7bd5",
        "pilot":  "#9b59b6",
        "main":   "#27ae60",
        "landed": "#7f8c8d",
    }
    if len(df) > 0:
        phases = df["cpas_phase"].astype(str).values
        times = df["t_s"].values
        seg_start = times[0]
        seg_phase = phases[0]
        for i in range(1, len(df)):
            if phases[i] != seg_phase:
                ax.axvspan(seg_start, times[i],
                           facecolor=phase_colors.get(seg_phase, "#dddddd"), alpha=0.22)
                seg_start = times[i]
                seg_phase = phases[i]
        ax.axvspan(seg_start, times[-1],
                   facecolor=phase_colors.get(seg_phase, "#dddddd"), alpha=0.22)

    ax.plot(df["t_s"], df["alt_m"] / 1000.0, color="black", linewidth=1.8, label="Altitude")

    event_rows = df[df["cpas_events"].astype(str) != ""] if "cpas_events" in df.columns else df.iloc[0:0]
    for _, r in event_rows.iterrows():
        ax.axvline(float(r["t_s"]), color="red", linestyle="--", alpha=0.65, linewidth=1.0)
        ax.annotate(str(r["cpas_events"]),
                    xy=(float(r["t_s"]), float(r["alt_m"]) / 1000.0),
                    xytext=(6, 10), textcoords="offset points", fontsize=8,
                    color="red", rotation=12)

    present_phases = [p for p in ["stowed", "drogue", "pilot", "main", "landed"]
                      if p in set(df["cpas_phase"].astype(str).unique())]
    handles = [Patch(facecolor=phase_colors[p], alpha=0.5, label=f"phase: {p}") for p in present_phases]
    handles += [_mlines.Line2D([0], [0], color="black", lw=1.8, label="Altitude"),
                _mlines.Line2D([0], [0], color="red", ls="--", label="CPAS event")]
    ax.legend(handles=handles, loc="upper right", framealpha=0.9, fontsize=9)
    ax.set_xlabel("Time s"); ax.set_ylabel("Altitude km")
    ax.set_title("CPAS deployment timeline overlaid on altitude")
    ax.grid(True, alpha=0.4)
    save_fig("06a2_cpas_altitude_phases")
    plt.show()


def plot_cpas_speed_dragarea(df, save_fig, **_):
    """Speed (left) and chute CD*A (right) over the terminal descent."""
    if "cpas_phase" not in df.columns or "cpas_drag_cdA_m2" not in df.columns:
        return

    fig, ax1 = plt.subplots(figsize=(11, 4.5))
    deployed_mask = df["cpas_phase"].astype(str) != "stowed"
    if deployed_mask.any():
        t_first_deploy = float(df.loc[deployed_mask, "t_s"].min())
        t_zoom_start = max(0.0, t_first_deploy - 20.0)
    else:
        t_zoom_start = 0.0
    mask_zoom = df["t_s"] >= t_zoom_start

    ax1.plot(df.loc[mask_zoom, "t_s"], df.loc[mask_zoom, "V_mps"],
             color="black", linewidth=1.8, label="Speed")
    ax1.set_xlabel("Time s"); ax1.set_ylabel("Speed m/s", color="black")
    ax1.grid(True, alpha=0.4)

    ax2 = ax1.twinx()
    ax2.plot(df.loc[mask_zoom, "t_s"], df.loc[mask_zoom, "cpas_drag_cdA_m2"],
             color="#27ae60", linewidth=1.8, label="Chute CD*A")
    ax2.set_ylabel("Chute CD*A  m^2", color="#27ae60")
    ax2.tick_params(axis="y", labelcolor="#27ae60")

    event_rows = df[(df["cpas_events"].astype(str) != "") & mask_zoom] \
        if "cpas_events" in df.columns else df.iloc[0:0]
    for _, r in event_rows.iterrows():
        ax1.axvline(float(r["t_s"]), color="red", linestyle="--", alpha=0.65, linewidth=1.0)
        ax1.annotate(str(r["cpas_events"]),
                     xy=(float(r["t_s"]), float(r["V_mps"])),
                     xytext=(6, 10), textcoords="offset points", fontsize=8, color="red")

    ax1.set_title("Terminal descent: speed (black) vs chute drag area (green)")
    fig.tight_layout()
    save_fig("06a3_cpas_speed_dragarea")
    plt.show()


def plot_speed(df, save_fig, **_):
    fig, ax = plt.subplots(figsize=(11, 4))
    plot_with_band(ax, df["t_s"], df["V_mps"], df["V_lo"], df["V_hi"],
                   "Speed with interval band", "V m/s")
    save_fig("06b_speed")
    plt.show()


def plot_gamma(df, save_fig, **_):
    fig, ax = plt.subplots(figsize=(11, 4))
    plot_with_band(ax, df["t_s"],
                   np.degrees(df["gamma_rad"]),
                   np.degrees(df["gamma_lo"]),
                   np.degrees(df["gamma_hi"]),
                   "Conditional uncertainty propagation for flight path angle", "Gamma deg",
                   nom_label="Nominal gamma", lo_label="Interval gamma low", hi_label="Interval gamma high")
    save_fig("06c_gamma")
    plt.show()


def plot_chi(df, save_fig, **_):
    fig, ax = plt.subplots(figsize=(11, 4))
    plot_with_band(ax, df["t_s"],
                   np.degrees(df["chi_rad"]),
                   np.degrees(df["chi_lo"]),
                   np.degrees(df["chi_hi"]),
                   "Conditional uncertainty propagation for heading angle", "Chi deg",
                   nom_label="Nominal chi", lo_label="Interval chi low", hi_label="Interval chi high")
    save_fig("06d_chi")
    plt.show()


def plot_lat_lon(df, save_fig, **_):
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    plot_with_band(axes[0], df["t_s"], np.degrees(df["phi_rad"]),
                   np.degrees(df["phi_lo"]), np.degrees(df["phi_hi"]),
                   "Latitude with interval band", "phi deg")
    plot_with_band(axes[1], df["t_s"], np.degrees(df["lam_rad"]),
                   np.degrees(df["lam_lo"]), np.degrees(df["lam_hi"]),
                   "Longitude with interval band", "lambda deg")
    save_fig("06e_lat_lon")
    plt.show()


def plot_state_widths(df, save_fig, **_):
    width_cols = [("r_width", "r_m"), ("V_width", "V_mps"),
                  ("gamma_width", "gamma_rad"), ("chi_width", "chi_rad")]
    fig, axes = plt.subplots(2, 2, figsize=(13, 7), sharex=True)
    recenter_steps = df[df["interval_step_recentered"] == 1]["t_s"] if "interval_step_recentered" in df.columns else []
    split_steps = df[df["interval_step_split"] == 1]["t_s"] if "interval_step_split" in df.columns else []
    for ax, (col, label) in zip(axes.ravel(), width_cols):
        if col not in df.columns:
            ax.set_title(f"{col} not logged"); continue
        ax.plot(df["t_s"], df[col], lw=1.5)
        for ts in recenter_steps:
            ax.axvline(ts, color="orange", alpha=0.25, lw=0.7)
        for ts in split_steps:
            ax.axvline(ts, color="red", alpha=0.35, lw=0.7)
        ax.set_title(f"Width of {label}")
        ax.set_xlabel("Time s"); ax.set_ylabel("width")
    fig.suptitle("Interval state-component widths (orange = recenter, red = split)")
    fig.tight_layout()
    save_fig("07a_state_widths")
    plt.show()


def plot_density_q(df, save_fig, **_):
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    plot_with_band(axes[0], df["t_s"], df["rho_kgm3"], df["interval_rho_lo"], df["interval_rho_hi"],
                   "Density with interval band", "rho kg/m^3")
    axes[0].set_yscale("log")
    plot_with_band(axes[1], df["t_s"], df["q_pa"], df["interval_q_lo"], df["interval_q_hi"],
                   "Dynamic pressure with interval band", "q Pa")
    save_fig("07b_density_q")
    plt.show()


def plot_heating_envelope(df, save_fig, **_):
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    axes[0].plot(df["t_s"], df["nominal_heat_qdot_max_hi"] / 1e6, color="C0", lw=2, label="Nominal qdot hi")
    if "heating_qdot_max_hi" in df.columns and df["heating_qdot_max_hi"].notna().any():
        axes[0].plot(df["t_s"], df["heating_qdot_max_hi"] / 1e6, color="C2", lw=1, ls="--", label="Interval qdot hi")
    axes[0].axhline(constants.HEAT_RATE_LIMIT_DEFAULT / 1e6, color="red", ls=":", label="heat rate limit")
    axes[0].set_xlabel("Time s"); axes[0].set_ylabel("qdot MW/m^2"); axes[0].set_title("Peak heating rate")
    axes[0].legend()

    axes[1].plot(df["t_s"], df["nominal_heat_Q_max_hi"] / 1e6, color="C0", lw=2, label="Nominal Q hi")
    if "heating_Q_max_hi" in df.columns and df["heating_Q_max_hi"].notna().any():
        axes[1].plot(df["t_s"], df["heating_Q_max_hi"] / 1e6, color="C2", lw=1, ls="--", label="Interval Q hi")
    axes[1].axhline(constants.HEAT_LOAD_LIMIT_DEFAULT / 1e6, color="red", ls=":", label="heat load limit")
    axes[1].set_xlabel("Time s"); axes[1].set_ylabel("Q MJ/m^2"); axes[1].set_title("Accumulated heat load")
    axes[1].legend()
    save_fig("08a_heating_envelope")
    plt.show()


def plot_heat_shield_map(df, save_fig, nominal_heat_shield=None, **_):
    if nominal_heat_shield is None:
        print("No nominal heat shield available")
        return
    shield = nominal_heat_shield
    n_rings = shield.num_rings
    n_sectors = shield.num_sectors
    Q_grid = np.zeros((n_rings, n_sectors))
    for ring in range(n_rings):
        for sector in range(n_sectors):
            idx = shield.cell_index(ring, sector)
            iv = shield.Q[idx]
            Q_grid[ring, sector] = 0.5 * (float(iv.lo) + float(iv.hi))
    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(Q_grid / 1e6, aspect="auto", cmap="hot", origin="lower")
    ax.set_xlabel("Sector index"); ax.set_ylabel("Ring index (0 = inner)")
    ax.set_title("Heat shield accumulated heat load Q (MJ/m^2)")
    plt.colorbar(im, ax=ax, label="Q MJ/m^2")
    save_fig("08b_heat_shield_map")
    plt.show()


def plot_guidance(df, save_fig, **_):
    g = df[df["guidance_updated"] == 1].copy()
    print("Guidance updates:", len(g))
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    axes[0].plot(g["t_s"], g["guidance_chosen_sigma_cmd_deg"], marker="o", ms=3, lw=1, label="chosen sigma_cmd")
    axes[0].plot(df["t_s"], np.degrees(df["sigma_actual_rad"]), color="C1", alpha=0.7, label="sigma_actual")
    axes[0].set_ylabel("Bank deg"); axes[0].set_title("Commanded vs actual bank")
    axes[0].legend()
    axes[1].plot(g["t_s"], g["guidance_selected_total_cost"], color="C3", label="total cost")
    axes[1].plot(g["t_s"], g["guidance_selected_geometry_cost"], color="C0", alpha=0.6, label="geometry cost")
    axes[1].plot(g["t_s"], g["guidance_selected_heat_penalty"], color="C2", alpha=0.6, label="heat penalty")
    axes[1].set_yscale("symlog"); axes[1].set_xlabel("Time s"); axes[1].set_ylabel("Cost")
    axes[1].set_title("Predictor-corrector cost breakdown"); axes[1].legend()
    save_fig("09a_guidance")
    plt.show()


def plot_candidate_distribution(df, save_fig, **_):
    g = df[df["guidance_updated"] == 1].copy()
    if not len(g):
        return
    counts = g["guidance_chosen_sigma_mag_deg"].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(8, 4))
    counts.plot(kind="bar", ax=ax, color="C0")
    ax.set_xlabel("Chosen bank magnitude deg"); ax.set_ylabel("Count")
    ax.set_title("Distribution of chosen predictor-corrector candidates")
    save_fig("09b_candidate_distribution")
    plt.show()


def plot_bank_error(df, save_fig, **_):
    err_deg = np.degrees(df["sigma_cmd_rad"] - df["sigma_actual_rad"])
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(df["t_s"], err_deg, color="C3", lw=1, label="sigma_cmd - sigma_actual")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("Time s"); ax.set_ylabel("Bank tracking error deg")
    ax.set_title("Bank tracking error")
    ax.legend()
    save_fig("10a_bank_error")
    plt.show()


def plot_roll_rate_accel(df, save_fig, **_):
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    axes[0].plot(df["t_s"], np.degrees(df["roll_rate_rad_s"]), color="C0")
    axes[0].set_xlabel("Time s"); axes[0].set_ylabel("Roll rate deg/s"); axes[0].set_title("Roll rate")
    axes[1].plot(df["t_s"], np.degrees(df["roll_accel_rad_s2"]), color="C1")
    axes[1].set_xlabel("Time s"); axes[1].set_ylabel("Roll accel deg/s^2"); axes[1].set_title("Roll acceleration")
    save_fig("10b_roll_rate_accel")
    plt.show()


def plot_torque(df, save_fig, **_):
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(df["t_s"], df["tau_roll_cmd_Nm"], color="C0", lw=1.2, label="Commanded torque")
    ax.plot(df["t_s"], df["tau_roll_realized_Nm"], color="C1", lw=1.0, alpha=0.8, label="Realized RCS torque (Mz)")
    ax.plot(df["t_s"],  df["tau_roll_capacity_Nm"], color="C2", ls="--", lw=0.8, label="+ capacity")
    ax.plot(df["t_s"], -df["tau_roll_capacity_Nm"], color="C2", ls="--", lw=0.8, label="- capacity")
    ax.set_xlabel("Time s"); ax.set_ylabel("Torque Nm")
    ax.set_title("Roll torque: commanded, realized, capacity envelope")
    ax.legend(ncol=2)
    save_fig("10c_torque")
    plt.show()


def plot_duty_vs_fired(df, save_fig, **_):
    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
    axes[0].plot(df["t_s"], df["requested_duty"], color="C0", lw=1, label="Requested duty")
    axes[0].set_ylabel("Duty (0-1)"); axes[0].set_title("Requested RCS duty cycle")
    axes[0].set_ylim(-0.05, 1.05)
    fired_idx = df[df["fired_this_step"] == 1]
    axes[1].scatter(fired_idx["t_s"], np.ones(len(fired_idx)), s=8, color="C3",
                    label=f"fired_this_step (n={len(fired_idx)})")
    axes[1].set_xlabel("Time s"); axes[1].set_ylabel("Fired")
    axes[1].set_yticks([0, 1]); axes[1].set_ylim(-0.2, 1.2)
    axes[1].legend()
    save_fig("10d_duty_vs_fired")
    plt.show()


def plot_thruster_raster(df, save_fig, thruster_fires_df=None, rcs_system=None, **_):
    if thruster_fires_df is None or len(thruster_fires_df) == 0:
        print("No thruster firings recorded.")
        return
    thruster_names = sorted(thruster_fires_df["thruster"].unique())
    y_index = {name: i for i, name in enumerate(thruster_names)}
    fires = thruster_fires_df.copy()
    fires["y"] = fires["thruster"].map(y_index)
    pos_names = set(getattr(rcs_system, "roll_pos_names", [])) if rcs_system is not None else set()
    neg_names = set(getattr(rcs_system, "roll_neg_names", [])) if rcs_system is not None else set()
    colors = ["C0" if n in pos_names else ("C3" if n in neg_names else "0.5") for n in thruster_names]
    fig, ax = plt.subplots(figsize=(13, max(3, 0.3 * len(thruster_names))))
    for name, c in zip(thruster_names, colors):
        sub = fires[fires["thruster"] == name]
        ax.scatter(sub["t_s"], np.full(len(sub), y_index[name]), s=6, color=c)
    ax.set_yticks(list(y_index.values())); ax.set_yticklabels(thruster_names)
    ax.set_xlabel("Time s")
    ax.set_title("RCS thruster firing raster (blue = +roll pod, red = -roll pod)")
    legend_handles = [Patch(color="C0", label="+roll pod"), Patch(color="C3", label="-roll pod"),
                      Patch(color="0.5", label="other")]
    ax.legend(handles=legend_handles, loc="upper right")
    save_fig("10e_thruster_raster")
    plt.show()


def plot_firing_rate(df, save_fig, dt_s=0.25, **_):
    win = max(1, int(round(5.0 / float(dt_s))))
    fire_rate_hz = df["fired_this_step"].rolling(win, min_periods=1).mean() / float(dt_s)
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(df["t_s"], fire_rate_hz, color="C3", lw=1.2)
    ax.set_xlabel("Time s"); ax.set_ylabel("Firing rate Hz (5 s window)")
    ax.set_title("Rolling RCS firing frequency")
    save_fig("10f_firing_rate")
    plt.show()


def plot_backlog(df, save_fig, **_):
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(df["t_s"], df["roll_pos_backlog_s"], color="C0", label="+roll pod backlog s")
    ax.plot(df["t_s"], df["roll_neg_backlog_s"], color="C3", label="-roll pod backlog s")
    ax.set_xlabel("Time s"); ax.set_ylabel("Backlog s")
    ax.set_title("RCS pulse backlog by pod")
    ax.legend()
    save_fig("10g_backlog")
    plt.show()


# =============================================================================
# master entry point
# =============================================================================

ALL_PLOTS = [
    plot_ground_track,
    plot_3d_descent_tube,
    plot_altitude,
    plot_cpas_altitude_phases,
    plot_cpas_speed_dragarea,
    plot_speed,
    plot_gamma,
    plot_chi,
    plot_lat_lon,
    plot_state_widths,
    plot_density_q,
    plot_heating_envelope,
    plot_heat_shield_map,
    plot_guidance,
    plot_candidate_distribution,
    plot_bank_error,
    plot_roll_rate_accel,
    plot_torque,
    plot_duty_vs_fired,
    plot_thruster_raster,
    plot_firing_rate,
    plot_backlog,
]


def render_all_figures(df, save_fig, **ctx: Any) -> list[str]:
    """
    Render every plot in ALL_PLOTS. Any single plot that raises is logged
    but does NOT abort the rest -- this protects the data-save step from
    a single bad plot.

    Returns the list of plots that ran successfully.
    """
    df = _normalize_columns(df)
    ok: list[str] = []
    for fn in ALL_PLOTS:
        try:
            fn(df, save_fig, **ctx)
            ok.append(fn.__name__)
        except Exception as e:  # noqa: BLE001
            print(f"  [WARN] {fn.__name__} failed: {type(e).__name__}: {e}")
        finally:
            plt.close("all")
    return ok
