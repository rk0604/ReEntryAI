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
# Plot categories (for modular display)
# =============================================================================
# Each plot function is tagged with one or more category strings via the
# @_in_category(...) decorator below. render_all_figures() reads
# `fn._categories` directly, so there is no parallel registry to keep in
# sync. The seven categories are:
#
#   trajectory  -- ground track + 3D tube
#   state       -- per-axis state histories with interval bands
#   cpas        -- everything chute-related
#   interval    -- supervisor diagnostics (state widths, density/q)
#   heating     -- heating envelope + heat shield map
#   guidance    -- predictor-corrector candidate analysis
#   rcs         -- bank tracking, roll, torque, raster, firing rate
#
# `want(category, selection)` is the cell-level helper used in the notebook.
# `render_all_figures(..., only=..., exclude=...)` is the bulk entry point.

def _in_category(*categories):
    """Decorator: tag a plot function with one or more category names."""
    def deco(fn):
        fn._categories = tuple(categories)
        return fn
    return deco


def want(category: str, selection) -> bool:
    """Return True if `category` should be rendered given `selection`.

    selection is one of:
      None         -> render everything
      str          -> render only this single category
      iterable[str]-> render any plot in the listed categories
    """
    if selection is None:
        return True
    if isinstance(selection, str):
        return category == selection
    return category in set(selection)


# Events that get text labels on the dense CPAS plots. Squid start/stop
# events still appear as red vertical lines but are not annotated, so the
# plots aren't drowned in overlapping text.
MAJOR_CPAS_EVENTS = {
    "fbc_jettison",
    "drogue_deploy", "drogue_disreef",
    "pilot_deploy",
    "main_deploy", "main_disreef_1", "main_disreef_2",
    "landed",
}


def _is_major_event(event_str: str) -> bool:
    if not event_str:
        return False
    return any(tok.strip() in MAJOR_CPAS_EVENTS for tok in str(event_str).split(","))


# =============================================================================
# helpers
# =============================================================================

def _normalize_columns(df):
    """
    Add aliased + memoized columns so the rest of this module can use a
    single naming convention regardless of which driver built the
    dataframe, and the dense CPAS plots don't each rerun .astype(str) /
    _is_major_event over the whole trajectory.

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

    # Cache the str-cast phase column and the major-event boolean mask
    # once. ~4 plots used to redo this per-row; now they share these.
    if "cpas_phase" in df.columns and "_cpas_phase_str" not in df.columns:
        df["_cpas_phase_str"] = df["cpas_phase"].astype(str)
    if "cpas_events" in df.columns and "_cpas_event_str" not in df.columns:
        df["_cpas_event_str"] = df["cpas_events"].astype(str)
        df["_cpas_has_event"] = df["_cpas_event_str"] != ""
        df["_cpas_has_major_event"] = df["_cpas_event_str"].apply(_is_major_event)
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


def _overlay_cpas_events(ax, df, lanes_y, *,
                          major_style=("red", "-", 0.85, 1.1),
                          minor_style=("red", "--", 0.30, 0.8)):
    """Render CPAS event vertical lines + staggered text labels for the
    major events. Used by plot_cpas_altitude_phases and
    plot_cpas_reefing_stages, which both used to copy-paste this block.

    `df` should already have been through _normalize_columns so the
    cached _cpas_has_event / _cpas_has_major_event columns are present.
    `lanes_y` is the list of y-positions to cycle the major-event labels
    through (units = whatever the y-axis is).
    """
    if "_cpas_has_event" not in df.columns:
        return
    mc, mls, ma, mw = minor_style
    Mc, Mls, Ma, Mw = major_style
    minor = df[df["_cpas_has_event"]]
    for ts in minor["t_s"].values:
        ax.axvline(float(ts), color=mc, linestyle=mls, alpha=ma, linewidth=mw)
    major = df[df["_cpas_has_major_event"]]
    for i, (_, r) in enumerate(major.iterrows()):
        t = float(r["t_s"])
        y = lanes_y[i % len(lanes_y)]
        ax.axvline(t, color=Mc, linestyle=Mls, alpha=Ma, linewidth=Mw)
        ax.annotate(str(r["_cpas_event_str"]),
                    xy=(t, y), xytext=(4, 0), textcoords="offset points",
                    fontsize=8, color=Mc, va="center", ha="left")


def _lat_lon_to_local_enu(phi_rad, lam_rad, phi0_rad, lam0_rad,
                          R_earth_km=None):
    """Flat-Earth local east/north (in metres) from a (phi, lam) array
    relative to a reference point. Used by both 3D-tube and chute-descent
    plots, which used to inline this conversion.
    """
    R = float(constants.RADIUS_EARTH) if R_earth_km is None else R_earth_km * 1000.0
    north_m = R * (phi_rad - phi0_rad)
    east_m = R * (lam_rad - lam0_rad) * math.cos(phi0_rad)
    return north_m, east_m


# =============================================================================
# individual plots
# =============================================================================

@_in_category('trajectory')
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


@_in_category("trajectory")
def plot_3d_descent_tube(df, save_fig, **_):
    track = df.copy()
    phi0 = float(track["phi_rad"].iloc[0])
    lam0 = float(track["lam_rad"].iloc[0])
    track["north_m"], track["east_m"] = _lat_lon_to_local_enu(
        track["phi_rad"], track["lam_rad"], phi0, lam0
    )

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


@_in_category('state')
def plot_altitude(df, save_fig, **_):
    fig, ax = plt.subplots(figsize=(11, 4))
    plot_with_band(ax, df["t_s"], df["alt_m"] / 1000.0,
                   df["interval_alt_lo"] / 1000.0, df["interval_alt_hi"] / 1000.0,
                   "Altitude with interval band", "Altitude km")
    save_fig("06a_altitude")
    plt.show()


@_in_category("cpas")
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
        phases = df["_cpas_phase_str"].values
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

    # Major/squid event overlays — shared with plot_cpas_reefing_stages.
    ylim_top = float(df["alt_m"].max()) / 1000.0
    _overlay_cpas_events(ax, df,
                          lanes_y=[ylim_top * 0.92, ylim_top * 0.82, ylim_top * 0.72])

    present_phases = [p for p in ["stowed", "drogue", "pilot", "main", "landed"]
                      if p in set(df["_cpas_phase_str"].unique())]
    handles = [Patch(facecolor=phase_colors[p], alpha=0.5, label=f"phase: {p}") for p in present_phases]
    handles += [_mlines.Line2D([0], [0], color="black", lw=1.8, label="Altitude"),
                _mlines.Line2D([0], [0], color="red", ls="-", label="major event"),
                _mlines.Line2D([0], [0], color="red", ls="--", alpha=0.3, label="squid event")]
    ax.legend(handles=handles, loc="upper right", framealpha=0.9, fontsize=9)
    ax.set_xlabel("Time s"); ax.set_ylabel("Altitude km")
    ax.set_title("CPAS deployment timeline overlaid on altitude")
    ax.grid(True, alpha=0.4)
    save_fig("06a2_cpas_altitude_phases")
    plt.show()


@_in_category('cpas')
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


@_in_category("cpas")
def plot_cpas_reefing_stages(df, save_fig, **_):
    """Open fraction + integer reefing stage during chute descent."""
    if "cpas_reefing_stage" not in df.columns or "cpas_open_fraction" not in df.columns:
        return
    deployed_mask = df["_cpas_phase_str"] != "stowed"
    if not deployed_mask.any():
        return
    t0 = float(df.loc[deployed_mask, "t_s"].min()) - 5.0
    sub = df[df["t_s"] >= t0]

    fig, ax1 = plt.subplots(figsize=(11, 4.5))
    ax1.step(sub["t_s"], sub["cpas_open_fraction"], where="post",
             color="#27ae60", linewidth=1.8, label="Open fraction")
    ax1.set_xlabel("Time s")
    ax1.set_ylabel("Open fraction", color="#27ae60")
    ax1.set_ylim(-0.1, 3.0)
    ax1.tick_params(axis="y", labelcolor="#27ae60")
    ax1.grid(True, alpha=0.4)

    ax2 = ax1.twinx()
    ax2.step(sub["t_s"], sub["cpas_reefing_stage"], where="post",
             color="#c0392b", linewidth=1.4, linestyle="--",
             label="Reefing stage (0..3)")
    ax2.set_ylabel("Reefing stage", color="#c0392b")
    ax2.set_ylim(-0.5, 3.5)
    ax2.set_yticks([0, 1, 2, 3])
    ax2.tick_params(axis="y", labelcolor="#c0392b")

    # Squid / major event overlays — shared helper.
    _overlay_cpas_events(ax1, sub, lanes_y=[2.7, 2.4, 2.1],
                          minor_style=("red", ":", 0.25, 0.8))

    ax1.set_title("CPAS opening profile: stair-stepped reefing stages + transient shocks")
    fig.tight_layout()
    save_fig("06a4_cpas_reefing_stages")
    plt.show()


@_in_category('cpas')
def plot_cpas_pendulum(df, save_fig, **_):
    """Pendulum swing angle and rate during the chute descent."""
    if "cpas_pendulum_angle_deg" not in df.columns:
        return
    deployed_mask = df["cpas_phase"].astype(str) != "stowed"
    if not deployed_mask.any():
        return
    t0 = float(df.loc[deployed_mask, "t_s"].min()) - 2.0
    sub = df[df["t_s"] >= t0]

    fig, axes = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
    axes[0].plot(sub["t_s"], sub["cpas_pendulum_angle_deg"],
                 color="#2980b9", lw=1.5, label="Pendulum angle")
    axes[0].axhline(0, color="black", lw=0.5)
    axes[0].set_ylabel("Pendulum angle deg")
    axes[0].set_title("Capsule pendulum oscillation under chutes")
    axes[0].grid(True, alpha=0.4)
    axes[0].legend()

    axes[1].plot(sub["t_s"], sub["cpas_pendulum_lateral_v_mps"],
                 color="#8e44ad", lw=1.2, label="Lateral velocity from pendulum")
    axes[1].axhline(0, color="black", lw=0.5)
    axes[1].set_xlabel("Time s")
    axes[1].set_ylabel("Lateral v m/s")
    axes[1].grid(True, alpha=0.4)
    axes[1].legend()
    save_fig("06a5_cpas_pendulum")
    plt.show()


@_in_category('cpas')
def plot_cpas_squidding(df, save_fig, **_):
    """Two-panel squidding diagnostic:
       top    -- number of mains currently squidding over time
       bottom -- effective chute drag area (shows the dips when mains squid)
    """
    if "cpas_num_mains_squidding" not in df.columns or "cpas_drag_cdA_m2" not in df.columns:
        return
    deployed = df["cpas_phase"].astype(str) == "main"
    if not deployed.any():
        return
    t0 = float(df.loc[deployed, "t_s"].min()) - 5.0
    sub = df[df["t_s"] >= t0]

    fig, axes = plt.subplots(2, 1, figsize=(11, 5.5), sharex=True)
    axes[0].step(sub["t_s"], sub["cpas_num_mains_squidding"], where="post",
                 color="#c0392b", lw=1.5)
    axes[0].set_ylabel("# mains squidding")
    axes[0].set_title("Main-chute squidding (asymmetric inflation) timeline")
    axes[0].set_yticks([0, 1, 2, 3])
    axes[0].set_ylim(-0.2, 3.5)
    axes[0].grid(True, alpha=0.4)

    axes[1].plot(sub["t_s"], sub["cpas_drag_cdA_m2"],
                 color="#27ae60", lw=1.4, label="Effective chute CD*A")
    axes[1].set_xlabel("Time s")
    axes[1].set_ylabel("Chute CD*A m^2")
    axes[1].grid(True, alpha=0.4)
    axes[1].legend(loc="upper right")
    fig.tight_layout()
    save_fig("06a7_cpas_squidding")
    plt.show()


@_in_category("cpas", "trajectory")
def plot_cpas_chute_descent_track(df, save_fig, **_):
    """Local-frame ground track from drogue deploy to landing, showing
    wind drift + pendulum wandering."""
    if "cpas_phase" not in df.columns:
        return
    deployed_mask = df["_cpas_phase_str"] != "stowed"
    if not deployed_mask.any():
        return
    sub = df.loc[deployed_mask]   # no .copy() — we only read from sub below
    if len(sub) < 2:
        return

    phi_arr = sub["phi_rad"].values
    lam_arr = sub["lam_rad"].values
    phi0 = float(phi_arr[0])
    lam0 = float(lam_arr[0])
    north_m, east_m = _lat_lon_to_local_enu(phi_arr, lam_arr, phi0, lam0)

    # Color the track by phase
    phase_colors = {"drogue": "#3a7bd5", "pilot": "#9b59b6",
                    "main": "#27ae60", "landed": "#7f8c8d"}
    phase_str = sub["_cpas_phase_str"].values
    fig, ax = plt.subplots(figsize=(9, 7))
    for phase, color in phase_colors.items():
        mask = phase_str == phase
        if mask.any():
            ax.plot(east_m[mask], north_m[mask],
                    color=color, lw=1.8, label=f"under {phase}")
    ax.scatter(east_m[0], north_m[0],
               s=100, marker="o", color="#3a7bd5", label="Drogue deploy", zorder=5)
    ax.scatter(east_m[-1], north_m[-1],
               s=120, marker="s", color="#c0392b", label="Touchdown", zorder=5)

    ax.set_xlabel("East m  (from drogue deploy)")
    ax.set_ylabel("North m  (from drogue deploy)")
    ax.set_title("Chute descent ground track (wind drift + pendulum wandering)")
    ax.grid(True, alpha=0.4)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_aspect("equal", "datalim")
    save_fig("06a6_cpas_chute_descent_track")
    plt.show()


@_in_category('state')
def plot_speed(df, save_fig, **_):
    fig, ax = plt.subplots(figsize=(11, 4))
    plot_with_band(ax, df["t_s"], df["V_mps"], df["V_lo"], df["V_hi"],
                   "Speed with interval band", "V m/s")
    save_fig("06b_speed")
    plt.show()


@_in_category('state')
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


@_in_category('state')
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


@_in_category('aero')
def plot_aero_coefficients(df, save_fig, **_):
    """Trim CD, CL, and L/D from the orion_cm_trim schedule over the run.

    These are the values fed into the EOMs every step. CPAS hooks
    inflate the CD line (and scale CL to zero) once chutes deploy, so
    the rapid swing in CD around T+ ~290–390 s in a default run is the
    drogue + main inflating, not an aero-database artefact.
    """
    if not {"CD", "CL", "LD"}.issubset(df.columns):
        return
    fig, axes = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

    # Top panel: CD and CL on twin y-axes (CD ~1.5, CL ~0.5 — very
    # different scales, so a twin axis avoids one trace dwarfing the other)
    ax_cd = axes[0]
    ax_cl = ax_cd.twinx()
    ax_cd.plot(df["t_s"], df["CD"], color="#5ddbe0", lw=1.4, label="C_D")
    ax_cl.plot(df["t_s"], df["CL"], color="#f6c043", lw=1.4, label="C_L")
    ax_cd.set_ylabel("C_D", color="#5ddbe0")
    ax_cl.set_ylabel("C_L", color="#f6c043")
    ax_cd.tick_params(axis="y", labelcolor="#5ddbe0")
    ax_cl.tick_params(axis="y", labelcolor="#f6c043")
    ax_cd.grid(True, alpha=0.4)
    ax_cd.set_title("Aero coefficients vs time (Bibb 2010 orion_cm_trim schedule)")
    # Combined legend
    h1, l1 = ax_cd.get_legend_handles_labels()
    h2, l2 = ax_cl.get_legend_handles_labels()
    ax_cd.legend(h1 + h2, l1 + l2, loc="upper right")

    # Bottom panel: L/D
    ax_ld = axes[1]
    ax_ld.plot(df["t_s"], df["LD"], color="#46d160", lw=1.4, label="L/D")
    ax_ld.set_xlabel("Time s")
    ax_ld.set_ylabel("L/D")
    ax_ld.grid(True, alpha=0.4)
    ax_ld.legend(loc="upper right")
    ax_ld.set_title("Lift-to-drag ratio")

    # Overlay CPAS events on both panels so you can see the chute kick
    if "_cpas_has_event" in df.columns:
        ylim_top_cd = float(df["CD"].max())
        _overlay_cpas_events(ax_cd, df,
                              lanes_y=[ylim_top_cd * 0.95,
                                       ylim_top_cd * 0.85,
                                       ylim_top_cd * 0.75])
        ylim_top_ld = float(df["LD"].max())
        _overlay_cpas_events(ax_ld, df,
                              lanes_y=[ylim_top_ld * 0.95,
                                       ylim_top_ld * 0.80,
                                       ylim_top_ld * 0.65])

    fig.tight_layout()
    save_fig("06f_aero_coefficients")
    plt.show()


@_in_category('state')
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


@_in_category('interval')
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


@_in_category('interval')
def plot_density_q(df, save_fig, **_):
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    plot_with_band(axes[0], df["t_s"], df["rho_kgm3"], df["interval_rho_lo"], df["interval_rho_hi"],
                   "Density with interval band", "rho kg/m^3")
    axes[0].set_yscale("log")
    plot_with_band(axes[1], df["t_s"], df["q_pa"], df["interval_q_lo"], df["interval_q_hi"],
                   "Dynamic pressure with interval band", "q Pa")
    save_fig("07b_density_q")
    plt.show()


@_in_category('heating')
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


@_in_category('heating')
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


@_in_category('guidance')
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


@_in_category('guidance')
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


@_in_category('rcs')
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


@_in_category('rcs')
def plot_roll_rate_accel(df, save_fig, **_):
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))
    axes[0].plot(df["t_s"], np.degrees(df["roll_rate_rad_s"]), color="C0")
    axes[0].set_xlabel("Time s"); axes[0].set_ylabel("Roll rate deg/s"); axes[0].set_title("Roll rate")
    axes[1].plot(df["t_s"], np.degrees(df["roll_accel_rad_s2"]), color="C1")
    axes[1].set_xlabel("Time s"); axes[1].set_ylabel("Roll accel deg/s^2"); axes[1].set_title("Roll acceleration")
    save_fig("10b_roll_rate_accel")
    plt.show()


@_in_category('rcs')
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


@_in_category('rcs')
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


@_in_category('rcs')
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


@_in_category('rcs')
def plot_firing_rate(df, save_fig, dt_s=0.25, **_):
    win = max(1, int(round(5.0 / float(dt_s))))
    fire_rate_hz = df["fired_this_step"].rolling(win, min_periods=1).mean() / float(dt_s)
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(df["t_s"], fire_rate_hz, color="C3", lw=1.2)
    ax.set_xlabel("Time s"); ax.set_ylabel("Firing rate Hz (5 s window)")
    ax.set_title("Rolling RCS firing frequency")
    save_fig("10f_firing_rate")
    plt.show()


@_in_category('rcs')
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
    plot_cpas_reefing_stages,
    plot_cpas_pendulum,
    plot_cpas_squidding,
    plot_cpas_chute_descent_track,
    plot_speed,
    plot_gamma,
    plot_chi,
    plot_lat_lon,
    plot_aero_coefficients,
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


def render_all_figures(df, save_fig, only=None, exclude=None, **ctx: Any) -> list[str]:
    """
    Render plots from ALL_PLOTS, optionally filtered by category.

    Args:
        df         : trajectory dataframe
        save_fig   : save_fig(name) callback
        only       : None, str, or iterable of category names. If given,
                     only plots whose categories intersect this set are
                     rendered. Default None = render everything.
        exclude    : None, str, or iterable of category names to skip.
        **ctx      : extra kwargs passed to each plot function
                     (target_phi_rad, thruster_fires_df, rcs_system,
                     nominal_heat_shield, dt_s, ...)

    A single plot that raises is logged as [WARN] but does not abort the
    rest. Returns the list of plot function names that ran successfully.
    """
    df = _normalize_columns(df)

    # Normalize only/exclude into sets
    def _to_set(x):
        if x is None:
            return None
        if isinstance(x, str):
            return {x}
        return set(x)
    only_set = _to_set(only)
    exclude_set = _to_set(exclude)

    ok: list[str] = []
    skipped: list[str] = []
    for fn in ALL_PLOTS:
        cats = set(getattr(fn, "_categories", ()))
        if only_set is not None and not (cats & only_set):
            skipped.append(fn.__name__)
            continue
        if exclude_set is not None and (cats & exclude_set):
            skipped.append(fn.__name__)
            continue
        try:
            fn(df, save_fig, **ctx)
            ok.append(fn.__name__)
        except Exception as e:  # noqa: BLE001
            print(f"  [WARN] {fn.__name__} failed: {type(e).__name__}: {e}")
        finally:
            plt.close("all")

    if skipped:
        print(f"  [skipped {len(skipped)} plots not in selected categories: "
              f"{', '.join(sorted(only_set)) if only_set else ''}]")
    return ok
