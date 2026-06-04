# Imports for the interval propagation and physical bank control notebook

import json
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import List

import constants
import AtmosphereModel
import control
import ReactionControl
import cpas
import plotting
import mission_config
import telemetry

from pathlib import Path

from interval_math import Interval, promote


# --- Load mission config -------------------------------------------------
# Single JSON file controls every routinely-varied parameter. Override the
# path via the SIM_CONFIG env var:
#   SIM_CONFIG=configs/example_failure_2of3.json python run_sim.py
MISSION_CFG = mission_config.load_config()
print(mission_config.summarize(MISSION_CFG))

# This file was originally a Jupyter notebook; a few cells still call
# IPython's display(). Provide a no-op shim so the script runs headless.
try:
    from IPython.display import display  # type: ignore
except Exception:  # noqa: BLE001
    def display(*args, **kwargs):  # type: ignore[no-redef]
        for a in args:
            print(a)
from math_3d import (
    IntervalSupervisorConfig,
    annotate_nominal_state_with_interval_supervisor,
    f_interval,
    nominal_state_to_interval_box,
    nominal_aero_forces_from_state,
    nominal_eom_step,
    nominal_heating_envelope_from_state,
)

from control import (
    ReentryState,
    GuidanceScheduler,
    SimpleBankGuidance,
    BasicObservationProvider,
    BankActuatorLimits,
    BankActuator,
    RollTorqueController,
    CapsuleControlConfig,
    CapsuleControlStack,
    wrap_to_pi,
)

from ReactionControl import (
    build_orion_cm_rcs_12,
)

# Helper functions for readable logging, plotting, and milestone 1 roll updates

STATE_NAMES = ["r_m", "phi_rad", "lam_rad", "V_mps", "gamma_rad", "chi_rad"]


def interval_mid(iv):
    # Returns the midpoint of an interval
    return 0.5 * (iv.lo + iv.hi)


def interval_width(iv):
    # returns the width of an interval
    return iv.hi - iv.lo


def box_mid(box):
    # Return midpoints for every interval in the box
    return [interval_mid(iv) for iv in box]


def box_widths(box):
    # Return the widths for every interval in the box
    return [interval_width(iv) for iv in box]


def deg(x_rad):
    # Convert from radians to degrees
    return math.degrees(x_rad)


def alt_from_r(r_m):
    # Convert radius from Earth center into geometric altitude
    return r_m - constants.RADIUS_EARTH

def wrap_angle_rad(x_rad):
    # Wrap an angle into the principal range
    return wrap_to_pi(float(x_rad))


def atmosphere_from_state_vector(x):
    # Read radius and speed from the translational state vector
    r_m, _, _, V_mps, _, _ = x

    # Convert radius into altitude for the atmosphere model
    alt_m = max(0.0, alt_from_r(r_m))

    # Query the nominal atmosphere at this altitude
    atm = AtmosphereModel.US_Standard_ATM(alt_m)

    # Use zero density if the atmosphere model returns None above its range
    rho = atm["rho_kgm3"] if atm["rho_kgm3"] is not None else 0.0
    T_K = atm["T_K"] if atm["T_K"] is not None else np.nan
    p_Pa = atm["p_Pa"] if atm["p_Pa"] is not None else np.nan

    # Compute dynamic pressure from density and speed
    q_pa = 0.5 * rho * V_mps * V_mps

    # Return a compact telemetry dictionary
    return {
        "alt_m": alt_m,
        "rho_kgm3": rho,
        "T_K": T_K,
        "p_Pa": p_Pa,
        "q_pa": q_pa,
    }


def milestone1_roll_update(
    sigma_actual_rad,
    roll_rate_rad_s,
    tau_rcs_z_Nm,
    dt_s,
    Izz_kgm2,
):
    # Convert applied roll torque into roll angular acceleration
    roll_accel_rad_s2 = float(tau_rcs_z_Nm) / float(Izz_kgm2)

    # Integrate roll rate forward one time step
    roll_rate_next_rad_s = float(roll_rate_rad_s) + roll_accel_rad_s2 * float(dt_s)

    # Integrate actual bank angle using the updated roll rate
    sigma_actual_next_rad = wrap_angle_rad(
        float(sigma_actual_rad) + roll_rate_next_rad_s * float(dt_s)
    )

    # Return the new roll state for this step
    return {
        "roll_accel_rad_s2": roll_accel_rad_s2,
        "roll_rate_next_rad_s": roll_rate_next_rad_s,
        "sigma_actual_next_rad": sigma_actual_next_rad,
    }

import constants
import importlib

print(constants.__file__)
print(hasattr(constants, "PREDICTOR_CORRECTOR_HORIZON_STEPS"))

importlib.reload(constants)

print(constants.__file__)
print(hasattr(constants, "PREDICTOR_CORRECTOR_HORIZON_STEPS"))
print(constants.PREDICTOR_CORRECTOR_HORIZON_STEPS)

# Vehicle aerodynamic control and reaction control system setup used by the notebook

# Define the physical and aerodynamic properties of the capsule
# Aero now follows the Orion CM hypersonic trim database (Bibb et al.).
# CD_const / CL_const are retained only as a fallback if aero_model is switched.
params = {
    "mass_kg": float(constants.CAPSULE_MASS_KG),
    "ref_area_m2": float(constants.CAPSULE_REFERENCE_AREA_M2),
    "aero_model": "orion_cm_trim",
    "CD_const": 1.15,
    "CL_const": 0.28,
}
mission_config.apply_aero_to_params(MISSION_CFG, params)

# Define the simulation time step and total duration
dt_s = float(constants.ENTRY_DT_S)
t_final_s = 1000.0
num_steps = int(t_final_s / dt_s)

# Define the moment of inertia for the roll axis
Izz_kgm2 = float(constants.CAPSULE_IZZ_KGM2)

# Define the landing target latitude and longitude (in radians).
# Pulled from the mission config so a single JSON change retargets the run.
_targets = mission_config.mission_targets(MISSION_CFG)
target_phi_rad = float(_targets["target_phi_rad"])
target_lam_rad = float(_targets["target_lam_rad"])
_cos_gamma_term_gate = float(_targets["cos_gamma_termination_gate"])

# Configure the guidance loop and scheduler
control_cfg = CapsuleControlConfig(
    guidance_period_s=float(constants.GUIDANCE_PERIOD_S),
)

guidance_scheduler = GuidanceScheduler(
    period_s=float(constants.GUIDANCE_PERIOD_S),
)

# Setup the heat aware short horizon predictor corrector guidance
# This still returns sigma_cmd only
# The downstream bank actuator and RCS roll chain remain unchanged
guidance_law = SimpleBankGuidance(
    target_phi_rad=float(target_phi_rad),
    target_lam_rad=float(target_lam_rad),
    params=dict(params),
    max_bank_deg=70.0,
    min_bank_deg=40.0,
    downrange_gain=1e-6,
    velocity_enable_mps=40_000.0,
    predictor_horizon_steps=int(constants.PREDICTOR_CORRECTOR_HORIZON_STEPS),
    prediction_dt_s=float(constants.ENTRY_DT_S),
    candidate_bank_deg=[float(v) for v in constants.PREDICTOR_CANDIDATE_BANK_DEG],
    weight_range=float(constants.PREDICTOR_COST_WEIGHT_RANGE),
    weight_heading=float(constants.PREDICTOR_COST_WEIGHT_HEADING),
    weight_cross_track=float(constants.PREDICTOR_COST_WEIGHT_CROSS_TRACK),
    weight_heat=float(constants.PREDICTOR_COST_WEIGHT_HEAT),
    heat_rate_limit=float(constants.HEAT_RATE_LIMIT_DEFAULT),
    heat_load_limit=float(constants.HEAT_LOAD_LIMIT_DEFAULT),
    infeasible_penalty=float(constants.PREDICTOR_HEAT_INFEASIBLE_PENALTY),
    invalid_interval_penalty=float(constants.PREDICTOR_INVALID_INTERVAL_PENALTY),
    include_zero_bank_candidate=True,
)

# Provide the observation model with the landing targets and capsule parameters
observation_provider = BasicObservationProvider(
    target_phi_rad=target_phi_rad,
    target_lam_rad=target_lam_rad,
    params=params,
)

# Define the physical limits for the bank actuator
bank_actuator = BankActuator(
    BankActuatorLimits(
        sigma_rate_max_rps=math.radians(20.0),
        sigma_accel_max_rps2=math.radians(40.0),
    )
)

# Tune the proportional derivative roll torque controller
roll_controller = RollTorqueController(
    kp_Nm_per_rad=float(constants.ROLL_KP_NM_PER_RAD),
    kd_Nm_per_rad_s=float(constants.ROLL_KD_NM_PER_RAD_S),
    max_torque_Nm=float(constants.ROLL_CMD_MAX_TORQUE_NM),
    sigma_deadband_rad=float(constants.ROLL_SIGMA_DEADBAND_RAD),
    rate_deadband_rad_s=float(constants.ROLL_RATE_DEADBAND_RAD_S),
)

# Assemble the full capsule control stack combining guidance and physical controllers
capsule_control = CapsuleControlStack(
    cfg=control_cfg,
    scheduler=guidance_scheduler,
    guidance=guidance_law,
    obs_provider=observation_provider,
    bank_actuator=bank_actuator,
    roll_controller=roll_controller,
)

# Initialize the reaction control system for the capsule
rcs_system = build_orion_cm_rcs_12()

# Configure the interval supervisor with uncertainty widths and fallback tuning
supervisor_cfg = IntervalSupervisorConfig(
    r_half_width_m=10.0,
    phi_half_width_rad=math.radians(0.01),
    lam_half_width_rad=math.radians(0.01),
    V_half_width_mps=20.0,
    gamma_half_width_rad=math.radians(0.15),
    chi_half_width_rad=math.radians(0.20),
    min_altitude_m=0.0,
    max_altitude_m=130_000.0,
    max_speed_mps=80_000.0,
    max_dynamic_pressure_pa=5.0e7,
    include_heating=True,
    heat_rate_limit=float(constants.HEAT_RATE_LIMIT_DEFAULT),
    heat_load_limit=float(constants.HEAT_LOAD_LIMIT_DEFAULT),
    interval_recenter_enabled=bool(constants.INTERVAL_RECENTER_ENABLED),
    interval_recenter_cadence_s=float(constants.INTERVAL_RECENTER_CADENCE_S),
    interval_recenter_width_thresholds=dict(constants.INTERVAL_RECENTER_WIDTH_THRESHOLDS),
    interval_box_split_enabled=bool(constants.INTERVAL_BOX_SPLIT_ENABLED),
    interval_box_split_max_depth=int(constants.INTERVAL_BOX_SPLIT_MAX_DEPTH),
    interval_box_split_width_thresholds=dict(constants.INTERVAL_BOX_SPLIT_WIDTH_THRESHOLDS),
    interval_denominator_safety_V_mps=float(constants.INTERVAL_DENOMINATOR_SAFETY_V_MPS),
    interval_denominator_safety_cos_gamma=float(constants.INTERVAL_DENOMINATOR_SAFETY_COS_GAMMA),
    interval_denominator_safety_cos_phi=float(constants.INTERVAL_DENOMINATOR_SAFETY_COS_PHI),
)

# Reset the control stack and reaction control system to their initial states
capsule_control.reset()
rcs_system.reset()

# Output a summary dictionary of the configuration
{
    "params": params,
    "dt_s": dt_s,
    "t_final_s": t_final_s,
    "num_steps": num_steps,
    "Izz_kgm2": Izz_kgm2,
    "guidance_period_s": constants.GUIDANCE_PERIOD_S,
    "predictor_horizon_steps": int(constants.PREDICTOR_CORRECTOR_HORIZON_STEPS),
    "predictor_candidate_bank_deg": list(constants.PREDICTOR_CANDIDATE_BANK_DEG),
    "heat_rate_limit": float(constants.HEAT_RATE_LIMIT_DEFAULT),
    "heat_load_limit": float(constants.HEAT_LOAD_LIMIT_DEFAULT),
    "interval_recenter_enabled": bool(constants.INTERVAL_RECENTER_ENABLED),
    "interval_recenter_cadence_s": float(constants.INTERVAL_RECENTER_CADENCE_S),
    "interval_box_split_enabled": bool(constants.INTERVAL_BOX_SPLIT_ENABLED),
    "interval_box_split_max_depth": int(constants.INTERVAL_BOX_SPLIT_MAX_DEPTH),
    "roll_pos_thrusters": rcs_system.roll_pos_names,
    "roll_neg_thrusters": rcs_system.roll_neg_names,
}

# Initial nominal translational and roll state

# State order
# r
# phi
# lam
# V
# gamma
# chi

# Define the starting conditions for the spacecraft trajectory
# (sourced from the loaded mission config)
x_nominal = mission_config.initial_state_vector(
    MISSION_CFG, float(constants.RADIUS_EARTH)
)

# Initialize the starting bank angle and roll rate to zero
sigma_actual_0_rad = 0.0
roll_rate_0_rad_s = 0.0

# Create the control state object using the established initial conditions
state_ctrl = ReentryState(
    r_m=float(x_nominal[0]),
    phi_rad=float(x_nominal[1]),
    lam_rad=float(x_nominal[2]),
    V_mps=float(x_nominal[3]),
    gamma_rad=float(x_nominal[4]),
    chi_rad=float(x_nominal[5]),
    sigma_actual_rad=float(sigma_actual_0_rad),
    roll_rate_rad_s=float(roll_rate_0_rad_s),
    sigma_cmd_rad=0.0,
    sigma_target_rad=0.0,
)

# Clear any previous states in the control and reaction control systems
capsule_control.reset()
rcs_system.reset()

# Output the starting conditions for verification
print("Initial altitude m", alt_from_r(x_nominal[0]))
print("Initial speed m per s", x_nominal[3])
print("Initial gamma deg", deg(x_nominal[4]))
print("Initial chi deg", deg(x_nominal[5]))
print("Initial sigma_actual deg", deg(state_ctrl.sigma_actual_rad))
print("Initial roll_rate deg per s", deg(state_ctrl.roll_rate_rad_s))

# interval supervisor settings

# these widths define translational state uncertainty only
# the bank angle used by the nominal path remains deterministic
# the interval layer annotates uncertainty around that flown nominal motion

# turn cadence recenter off for rl style data generation
# this keeps recenter tied to interval health rather than fixed time buckets
interval_recenter_use_cadence = False

supervisor_cfg = IntervalSupervisorConfig(
    r_half_width_m=10.0,
    phi_half_width_rad=math.radians(0.01),
    lam_half_width_rad=math.radians(0.01),
    V_half_width_mps=20.0,
    gamma_half_width_rad=math.radians(0.15),
    chi_half_width_rad=math.radians(0.20),
    min_altitude_m=0.0,
    max_altitude_m=130_000.0,
    max_speed_mps=80_000.0,
    max_dynamic_pressure_pa=5.0e7,
    include_heating=True,
    heat_rate_limit=float(constants.HEAT_RATE_LIMIT_DEFAULT),
    heat_load_limit=float(constants.HEAT_LOAD_LIMIT_DEFAULT),
    interval_recenter_enabled=bool(constants.INTERVAL_RECENTER_ENABLED),
    interval_recenter_use_cadence=bool(interval_recenter_use_cadence),
    interval_recenter_cadence_s=float(constants.INTERVAL_RECENTER_CADENCE_S),
    interval_recenter_width_thresholds=dict(constants.INTERVAL_RECENTER_WIDTH_THRESHOLDS),
    interval_box_split_enabled=bool(constants.INTERVAL_BOX_SPLIT_ENABLED),
    interval_box_split_max_depth=int(constants.INTERVAL_BOX_SPLIT_MAX_DEPTH),
    interval_box_split_width_thresholds=dict(constants.INTERVAL_BOX_SPLIT_WIDTH_THRESHOLDS),
    interval_denominator_safety_V_mps=float(constants.INTERVAL_DENOMINATOR_SAFETY_V_MPS),
    interval_denominator_safety_cos_gamma=float(constants.INTERVAL_DENOMINATOR_SAFETY_COS_GAMMA),
    interval_denominator_safety_cos_phi=float(constants.INTERVAL_DENOMINATOR_SAFETY_COS_PHI),
)

x_interval = None

supervisor_cfg

# Setup objects for the physical closed loop bank control path

# This notebook no longer uses a hand written sigma schedule.
# The nominal path now follows:
# guidance produces sigma_cmd
# bank actuator produces sigma_target
# roll controller produces tau_roll_cmd_Nm
# RCS produces body torque
# roll dynamics produce sigma_actual
# translational dynamics fly sigma_actual

capsule_control.reset()
rcs_system.reset()

trajectory_id = "traj_000"

# This object is carried forward on the live trajectory
# Candidate rollouts clone it so they do not contaminate each other
interval_heat_shield = None

# Store failed heat guidance cycles for later RL teaching
failed_heat_rows = []

# Output folder for the run_sim.py driver. Comes from the mission config
# (typically runs/<run_id>) so each config gets its own folder. We append
# "_py" so the run_sim and notebook drivers can co-exist without trampling
# each other when both run the same config.
PY_OUTPUT_DIR = Path(str(MISSION_CFG.output_dir) + "_py")
PY_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"[run_sim] writing outputs to {PY_OUTPUT_DIR}")

# Define output file names once here so later save cells can reuse them
success_trajectory_csv_path = str(PY_OUTPUT_DIR / "trajectory_success.csv")
failed_heat_step_csv_path = str(PY_OUTPUT_DIR / "trajectory_failed_heat_steps.csv")
failed_heat_episode_csv_path = str(PY_OUTPUT_DIR / "trajectory_failed_heat_episodes.csv")
cpas_events_csv_path = str(PY_OUTPUT_DIR / "cpas_events.csv")
run_summary_json_path = str(PY_OUTPUT_DIR / "run_summary.json")
landing_summary_json_path = str(PY_OUTPUT_DIR / "landing_summary.json")

# Define the failed heat log schema explicitly
# This makes the exported csv stable even if there are zero failed rows
failed_heat_columns = [
    "trajectory_id",
    "step",
    "guidance_cycle_index",
    "time_s",
    "r_m",
    "phi_rad",
    "lam_rad",
    "V_mps",
    "gamma_rad",
    "chi_rad",
    "target_phi_rad",
    "target_lam_rad",
    "sigma_cmd_candidate_deg_json",
    "sigma_cmd_candidate_cost_json",
    "sigma_cmd_candidate_heat_flag_json",
    "sigma_cmd_candidate_failure_reason_json",
    "selected_sigma_cmd_rad",
    "selected_sigma_cmd_deg",
    "sigma_actual_rad",
    "range_to_go_m",
    "heading_error_rad",
    "cross_track_error_m",
    "dynamic_pressure_pa",
    "density_kgm3",
    "predicted_max_heating_rate_lo",
    "predicted_max_heating_rate_hi",
    "predicted_max_heat_load_lo",
    "predicted_max_heat_load_hi",
    "heat_feasible_flag",
    "violation_amount",
    "first_violation_step",
    "first_violation_time_s",
    "failure_reason",
    "interval_propagation_invalid",
    "geometry_cost",
    "heat_penalty",
    "total_cost",
]

# Keep the main log column list readable for inspection
# The DataFrame is still built from row dictionaries in the main loop
log_columns = [
    "step",
    "t_s",
    "guidance_updated",
    "trajectory_id",
    "guidance_cycle_index",
    "sigma_cmd_rad",
    "sigma_target_rad",
    "sigma_actual_rad",
    "roll_rate_rad_s",
    "roll_accel_rad_s2",
    "tau_roll_cmd_Nm",
    "tau_roll_capacity_Nm",
    "requested_duty",
    "fired_this_step",
    "active_thrusters",
    "torque_z_from_rcs",
    "force_x_from_rcs",
    "force_y_from_rcs",
    "force_z_from_rcs",
    "r_m",
    "phi_rad",
    "lam_rad",
    "V_mps",
    "gamma_rad",
    "chi_rad",
    "alt_m",
    "rho_kgm3",
    "q_pa",
    "T_K",
    "p_Pa",
    "drag_mag_N",
    "lift_mag_N",
    "gravity_mps2",
    "range_to_go_m",
    "heading_error_rad",
    "cross_track_error_m",
    "interval_sigma_lo",
    "interval_sigma_hi",
    "interval_rho_lo",
    "interval_rho_hi",
    "interval_q_lo",
    "interval_q_hi",
    "interval_alt_lo",
    "interval_alt_hi",
    "width_r_m",
    "width_phi_rad",
    "width_lam_rad",
    "width_V_mps",
    "width_gamma_rad",
    "width_chi_rad",
    "dx_width_r",
    "dx_width_phi",
    "dx_width_lam",
    "dx_width_V",
    "dx_width_gamma",
    "dx_width_chi",
    "heating_qdot_max_lo",
    "heating_qdot_max_hi",
    "heating_qdot_mean_lo",
    "heating_qdot_mean_hi",
    "heating_Q_max_lo",
    "heating_Q_max_hi",
    "guidance_selected_candidate_index",
    "guidance_any_feasible_candidate",
    "guidance_selected_candidate_heat_feasible",
    "guidance_chosen_sigma_cmd_deg",
    "guidance_chosen_sigma_mag_deg",
    "guidance_heading_deadband_deg",
    "guidance_current_altitude_rate_mps",
    "guidance_selected_geometry_cost",
    "guidance_selected_heat_penalty",
    "guidance_selected_total_cost",
    "guidance_selected_failure_reason",
    "guidance_selected_violation_amount",
    "guidance_selected_first_violation_step",
    "guidance_selected_first_violation_time_s",
    "guidance_selected_max_heating_rate_lo",
    "guidance_selected_max_heating_rate_hi",
    "guidance_selected_max_heat_load_lo",
    "guidance_selected_max_heat_load_hi",
    "guidance_candidate_count",
    "guidance_candidate_sigma_deg_json",
    "guidance_candidate_cost_json",
    "guidance_candidate_heat_flag_json",
    "guidance_candidate_failure_reason_json",
    "safety_status",
    "safety_checks",
    "layer_indices",
]

print("Closed loop setup ready")
print("Guidance period s", float(constants.GUIDANCE_PERIOD_S))
print("Physics dt s", dt_s)
print("Predictor horizon steps", int(constants.PREDICTOR_CORRECTOR_HORIZON_STEPS))
print("Trajectory id", trajectory_id)
print("Positive roll thrusters", rcs_system.roll_pos_names)
print("Negative roll thrusters", rcs_system.roll_neg_names)

# Nominal translational step driven by the physically flown sigma_actual

# The nominal translational dynamics are advanced using the actual bank angle
# produced by the guidance chain, roll controller, RCS allocation, and the
# milestone 1 roll state update.
#
# This wrapper stays intentionally small.
# The corrected translational physics live inside nominal_eom_step.

def nominal_step_closed_loop(x_nominal_old, sigma_actual_rad, params, dt_s):
    # Run one nominal Euler step using the corrected shared implementation
    step_info = nominal_eom_step(
        x=x_nominal_old,
        sigma_rad=float(sigma_actual_rad),
        params=params,
        dt_s=float(dt_s),
    )

    # Extract the advanced state and the derivative vector
    x_nominal_new = [float(v) for v in step_info["x_next"]]
    dx_nominal = [float(v) for v in step_info["x_dot"]]

    # Return both the compact outputs and the full step info
    return x_nominal_new, dx_nominal, step_info

# One step check for the interval supervisor under deterministic flown sigma

# This confirms the interval layer uses the actual bank angle from the nominal
# path as a punctual interval input rather than creating a separate control mode.

t0 = 0.0

result0 = annotate_nominal_state_with_interval_supervisor(
    x_nominal_old=x_nominal,
    params=params,
    sigma_actual_after_rad=float(state_ctrl.sigma_actual_rad),
    x_interval_old=None,
    supervisor_cfg=supervisor_cfg,
    dt_s=dt_s,
    heat_shield=None,
    t_s=t0,
)

print("Sigma actual used deg", deg(state_ctrl.sigma_actual_rad))
print("Altitude interval m", result0.altitude_interval)
print("Density interval", result0.rho_interval)
print("Dynamic pressure interval Pa", result0.q_interval)
print("State widths new", result0.state_widths_new)
print("Safety status", result0.safety_status)

# ensure the heat shield is active
print(result0.heat_shield)
print(result0.heating_qdot_max_interval)
print(result0.heating_Q_max_interval)

# Diagnostic replay for the first 20 physical steps
#
# This version keeps the same high level replay goal
# The main change is that the nominal closed loop step now handles
# the internal RCS substeps automatically
# That avoids calling the timed RCS model with the outer 025 second step

from point_math_3d import make_initial_capsule_attitude, step_closed_loop_milestone1

diag_steps = 20

# Rebuild fresh control and RCS objects so the replay starts from a clean state
guidance_scheduler_diag = GuidanceScheduler(
    period_s=float(constants.GUIDANCE_PERIOD_S),
)

guidance_law_diag = SimpleBankGuidance(
    target_phi_rad=float(target_phi_rad),
    target_lam_rad=float(target_lam_rad),
    params=dict(params),
    max_bank_deg=70.0,
    min_bank_deg=40.0,
    downrange_gain=1e-6,
    velocity_enable_mps=30_000.0,
    predictor_horizon_steps=int(constants.PREDICTOR_CORRECTOR_HORIZON_STEPS),
    prediction_dt_s=float(constants.ENTRY_DT_S),
    candidate_bank_deg=[float(v) for v in constants.PREDICTOR_CANDIDATE_BANK_DEG],
    weight_range=float(constants.PREDICTOR_COST_WEIGHT_RANGE),
    weight_heading=float(constants.PREDICTOR_COST_WEIGHT_HEADING),
    weight_cross_track=float(constants.PREDICTOR_COST_WEIGHT_CROSS_TRACK),
    weight_heat=float(constants.PREDICTOR_COST_WEIGHT_HEAT),
    heat_rate_limit=float(constants.HEAT_RATE_LIMIT_DEFAULT),
    heat_load_limit=float(constants.HEAT_LOAD_LIMIT_DEFAULT),
    infeasible_penalty=float(constants.PREDICTOR_HEAT_INFEASIBLE_PENALTY),
    invalid_interval_penalty=float(constants.PREDICTOR_INVALID_INTERVAL_PENALTY),
    include_zero_bank_candidate=True,
)

observation_provider_diag = BasicObservationProvider(
    target_phi_rad=float(target_phi_rad),
    target_lam_rad=float(target_lam_rad),
    params=dict(params),
)

bank_actuator_diag = BankActuator(
    BankActuatorLimits(
        sigma_rate_max_rps=math.radians(20.0),
        sigma_accel_max_rps2=math.radians(40.0),
    )
)

roll_controller_diag = RollTorqueController(
    kp_Nm_per_rad=float(constants.ROLL_KP_NM_PER_RAD),
    kd_Nm_per_rad_s=float(constants.ROLL_KD_NM_PER_RAD_S),
    max_torque_Nm=float(constants.ROLL_CMD_MAX_TORQUE_NM),
    sigma_deadband_rad=float(constants.ROLL_SIGMA_DEADBAND_RAD),
    rate_deadband_rad_s=float(constants.ROLL_RATE_DEADBAND_RAD_S),
)

capsule_control_diag = CapsuleControlStack(
    cfg=CapsuleControlConfig(
        guidance_period_s=float(constants.GUIDANCE_PERIOD_S),
    ),
    scheduler=guidance_scheduler_diag,
    guidance=guidance_law_diag,
    obs_provider=observation_provider_diag,
    bank_actuator=bank_actuator_diag,
    roll_controller=roll_controller_diag,
)

rcs_system_diag = build_orion_cm_rcs_12()

capsule_control_diag.reset()
rcs_system_diag.reset()

# Reset the initial translational and roll state
x_nom_diag = [
    constants.RADIUS_EARTH + 80_000.0,
    math.radians(0.0),
    math.radians(0.0),
    6_500.0,
    math.radians(-5.5),
    math.radians(90.0),
]

# Keep the realized attitude as the single source of truth for bank and roll rate
att_diag = make_initial_capsule_attitude(0.0)

# Keep interval context available for guidance candidate rollout hooks
x_interval_diag = None
interval_heat_shield_diag = None

diag_rows = []

print("Target latitude deg", math.degrees(target_phi_rad))
print("Target longitude deg", math.degrees(target_lam_rad))
print("Initial chi deg", math.degrees(x_nom_diag[5]))

# Build a controller bridge that matches the nominal closed loop step interface
def control_step_fn_diag(
    t_s: float,
    x_state: List[float],
    sigma_actual_rad: float,
    roll_rate_rad_s: float,
):
    # Convert the translational state into the controller facing state object
    r_m, phi_rad, lam_rad, V_mps, gamma_rad, chi_rad = [float(v) for v in x_state]

    state_ctrl = ReentryState(
        r_m=float(r_m),
        phi_rad=float(phi_rad),
        lam_rad=float(lam_rad),
        V_mps=float(V_mps),
        gamma_rad=float(gamma_rad),
        chi_rad=float(chi_rad),
        sigma_actual_rad=float(sigma_actual_rad),
        roll_rate_rad_s=float(roll_rate_rad_s),
        sigma_cmd_rad=0.0,
        sigma_target_rad=0.0,
    )

    # Use the currently flown bank to compute physical lift and drag magnitudes
    aero_before = nominal_aero_forces_from_state(
        x=x_state,
        sigma_rad=float(sigma_actual_rad),
        params=params,
    )

    return capsule_control_diag.step(
        t_s=float(t_s),
        dt_s=float(dt_s),
        state=state_ctrl,
        lift_N=float(aero_before["lift_mag_N"]),
        drag_N=float(aero_before["drag_mag_N"]),
        mass_kg=float(params["mass_kg"]),
    )

for k in range(diag_steps):
    t_s = k * dt_s

    # Pass the live interval and heat context into guidance
    guidance_law_diag.set_prediction_context(
        x_interval_old=x_interval_diag,
        heat_shield=interval_heat_shield_diag,
        supervisor_cfg=supervisor_cfg,
        trajectory_id="diag_traj",
        guidance_cycle_index=k,
    )

    # Run one full nominal closed loop step
    # This now performs controller update internal RCS substeps
    # roll attitude propagation and outer translational propagation
    step_result = step_closed_loop_milestone1(
        t_s=float(t_s),
        x_trans=list(x_nom_diag),
        att=att_diag,
        params=params,
        control_step_fn=control_step_fn_diag,
        rcs=rcs_system_diag,
        dt_s=float(dt_s),
        Izz_kgm2=float(Izz_kgm2),
    )

    ctrl_out = step_result.control_out
    roll_step = step_result.roll_step
    x_nom_next = list(step_result.x_new)
    dx_nom = list(step_result.dx_trans)

    sigma_actual_next = float(step_result.sigma_actual_after_rad)
    roll_rate_next = float(step_result.att_new.omega_b_rad_s[2])

    # Pull guidance debug fields
    guidance_debug = dict(ctrl_out.guidance_debug)

    diag_rows.append(
        {
            "time_s": float(t_s),
            "guidance_updated": int(bool(ctrl_out.guidance_updated)),
            "heading_error_rad": float(ctrl_out.obs.get("heading_error_rad", np.nan)),
            "heading_offset_psi_rad": float(ctrl_out.obs.get("heading_offset_psi_rad", np.nan)),
            "sigma_cmd_rad": float(ctrl_out.sigma_cmd_rad),
            "sigma_target_rad": float(ctrl_out.sigma_target_rad),
            "sigma_actual_rad": float(sigma_actual_next),
            "tau_roll_cmd_Nm": float(ctrl_out.tau_roll_cmd_Nm),
            "requested_duty": float(roll_step.requested_duty),
            "fired_this_step": int(bool(roll_step.fired_this_step)),
            "num_internal_steps": int(roll_step.num_internal_steps),
            "num_fired_internal_steps": int(roll_step.num_fired_internal_steps),
            "torque_z_from_rcs": float(roll_step.wrench.torque_b_Nm[2]),
            "roll_pos_backlog_s": float(roll_step.roll_pos_backlog_s),
            "roll_neg_backlog_s": float(roll_step.roll_neg_backlog_s),
            "roll_pos_is_on": int(bool(roll_step.roll_pos_is_on)),
            "roll_neg_is_on": int(bool(roll_step.roll_neg_is_on)),
            "chi_rad": float(x_nom_diag[5]),
            "chi_dot": float(dx_nom[5]),
            "selected_candidate_index": int(guidance_debug.get("selected_candidate_index", -1)),
            "candidate_sigma_values": str(guidance_debug.get("candidate_sigma_deg", [])),
            "candidate_total_costs": str(guidance_debug.get("candidate_total_cost", [])),
            "candidate_heat_flags": str(guidance_debug.get("candidate_heat_feasible", [])),
            "selected_failure_reason": str(guidance_debug.get("selected_failure_reason", "")),
        }
    )

    # Carry the state forward
    x_nom_diag = list(x_nom_next)
    att_diag = step_result.att_new

diag_df = pd.DataFrame(diag_rows)

pd.set_option("display.max_colwidth", 200)
diag_df

import json

from control import ReentryState
from point_math_3d import make_initial_capsule_attitude, step_closed_loop_milestone1, aero_forces

# Start from the declared initial nominal state
x_nom = list(x_nominal)

# Let the interval supervisor initialize itself on the first step
x_interval = None

# Carry the live interval heat shield state forward only while interval is active
interval_heat_shield = None

# Carry the live nominal heat shield state through the full run
nominal_heat_shield = None

# Track whether interval propagation is still active
interval_active = True
interval_failed = False
interval_failure_time_s = np.nan
interval_failure_reason = ""
interval_failure_kind = ""
last_valid_interval_state = None
last_valid_interval_time_s = np.nan

# Track nominal heat stop status separately from interval failure
nominal_heat_violation = False
nominal_heat_violation_reason = ""
sim_terminated_by_nominal_heat = False
termination_reason = ""

# Keep the realized attitude as the single source of truth for bank and roll rate
sigma_actual_init = 0.0
if "state_ctrl" in globals():
    sigma_actual_init = float(getattr(state_ctrl, "sigma_actual_rad", 0.0))

att_nom = make_initial_capsule_attitude(sigma_actual_init)

# Main row storage
rows = []

# Build a controller bridge that matches the nominal closed loop step interface
def control_step_fn_main(t_s, x_state, sigma_actual_rad, roll_rate_rad_s):
    # Convert the translational state into the controller facing state object
    r_m, phi_rad, lam_rad, V_mps, gamma_rad, chi_rad = [float(v) for v in x_state]

    state_ctrl_main = ReentryState(
        r_m=float(r_m),
        phi_rad=float(phi_rad),
        lam_rad=float(lam_rad),
        V_mps=float(V_mps),
        gamma_rad=float(gamma_rad),
        chi_rad=float(chi_rad),
        sigma_actual_rad=float(sigma_actual_rad),
        roll_rate_rad_s=float(roll_rate_rad_s),
        sigma_cmd_rad=0.0,
        sigma_target_rad=0.0,
    )

    # Compute aerodynamic force components using the currently flown bank angle
    Dr, Dtheta, Dphi, Lr, Ltheta, Lphi, _, _ = aero_forces(
        r=float(r_m),
        V=float(V_mps),
        gamma=float(gamma_rad),
        chi=float(chi_rad),
        params=params,
        sigma_actual_rad=float(sigma_actual_rad),
    )

    drag_mag_N = math.sqrt(Dr * Dr + Dtheta * Dtheta + Dphi * Dphi)
    lift_mag_N = math.sqrt(Lr * Lr + Ltheta * Ltheta + Lphi * Lphi)

    return capsule_control.step(
        t_s=float(t_s),
        dt_s=float(dt_s),
        state=state_ctrl_main,
        lift_N=float(lift_mag_N),
        drag_N=float(drag_mag_N),
        mass_kg=float(params["mass_kg"]),
    )

# --- CPAS (parachute) integration ---
import random as _cpas_rng_mod
import os as _os_cpas

# Build CPASConfig from the mission config; fall back to defaults for any
# field the JSON omits. CPAS_NUM_MAINS env var still overrides if set, for
# quick deterministic failure-mode testing without editing the JSON.
_cpas_config = mission_config.build_cpas_config(MISSION_CFG)
_cpas_num_env = _os_cpas.environ.get("CPAS_NUM_MAINS", "").strip()
if _cpas_num_env:
    _n_mains = max(0, min(3, int(_cpas_num_env)))
    print(f"[CPAS] CPAS_NUM_MAINS={_n_mains} env override -> deterministic failure mode")
    _cpas_config.enable_stochastic_failure = False
    _cpas_config.num_mains_operational = _n_mains
cpas_inst = cpas.CPAS(config=_cpas_config, rng=_cpas_rng_mod.Random(int(MISSION_CFG.seed)))
print(f"[CPAS] mains operational: {cpas_inst.state.mains_operational_mask}")

cpas_events_log = []  # one row per phase transition / disreef / squid / FBC
params["cpas_drag_cdA_extra_m2"] = 0.0
params["cpas_lift_scale"] = 1.0
_R_EARTH = float(constants.RADIUS_EARTH)

# Cumulative RCS propellant burned (kg), integrated over the whole run.
rcs_fuel_used_kg_total = 0.0

for k in range(num_steps):
    t_s = k * dt_s

    sigma_actual_rad = float(att_nom.sigma_rel_rad)
    roll_rate_rad_s = float(att_nom.omega_b_rad_s[2])

    # CPAS state machine update (before termination checks so it can deploy
    # in the same step the capsule first drops below the trigger altitude).
    alt_now_m = alt_from_r(x_nom[0])
    atm_now = atmosphere_from_state_vector(x_nom)
    T_now = atm_now.get("T_K", float("nan"))
    if T_now is not None and T_now == T_now and T_now > 0.0:
        mach_now = float(x_nom[3]) / math.sqrt(1.4 * 287.0 * float(T_now))
    else:
        mach_now = None
    cpas_out = cpas_inst.step(t_s=t_s, alt_m=alt_now_m, V_mps=float(x_nom[3]),
                              mach=mach_now, dt_s=float(dt_s))
    params["cpas_drag_cdA_extra_m2"] = float(cpas_out.drag_area_cdA_m2)
    params["cpas_lift_scale"] = float(cpas_out.lift_scale)
    if cpas_out.mass_shed_delta_kg > 0.0:
        params["mass_kg"] = max(1.0, float(params["mass_kg"]) - float(cpas_out.mass_shed_delta_kg))
    for ev in cpas_out.events:
        cpas_events_log.append({
            "t_s": float(t_s),
            "step": int(k),
            "alt_m": float(alt_now_m),
            "V_mps": float(x_nom[3]),
            "mach": float(mach_now) if mach_now is not None else float("nan"),
            "phase": str(cpas_out.phase),
            "event": str(ev),
            "reefing_stage": int(cpas_out.reefing_stage),
            "num_mains_active": int(cpas_out.num_mains_active),
        })
        print(f"  CPAS {ev:>18s} at t={t_s:6.1f}s alt={alt_now_m:7.1f}m V={float(x_nom[3]):6.1f}m/s")

    nom_cos_gamma = math.cos(x_nom[4])

    # cos_gamma safety gate only applies pre-deployment. Once chutes are out
    # the bank/heading channels are dead and a near-vertical descent is the
    # expected physical mode.
    if cpas_out.phase == "stowed" and abs(nom_cos_gamma) < _cos_gamma_term_gate:
        termination_reason = "nominal_cos_gamma_too_small"
        print(f"Stopping at t = {t_s:.2f} s because nominal cos gamma is too small for stable heading propagation")
        break

    if x_nom[3] <= 0.1:
        termination_reason = "nominal_speed_too_low"
        print(f"Stopping at t = {t_s:.2f} s because nominal speed is too low for this simplified regime")
        break

    if alt_now_m <= 0.0:
        termination_reason = "nominal_ground_reached"
        print(f"Stopping at t = {t_s:.2f} s because the nominal trajectory reached the ground")
        break

    guidance_law.set_prediction_context(
        x_interval_old=x_interval if interval_active else None,
        heat_shield=interval_heat_shield if interval_active else None,
        supervisor_cfg=supervisor_cfg,
        trajectory_id=trajectory_id,
        guidance_cycle_index=k,
        interval_active=bool(interval_active),
    )

    step_result = step_closed_loop_milestone1(
        t_s=float(t_s),
        x_trans=list(x_nom),
        att=att_nom,
        params=params,
        control_step_fn=control_step_fn_main,
        rcs=rcs_system,
        dt_s=float(dt_s),
        Izz_kgm2=float(Izz_kgm2),
    )

    ctrl_out = step_result.control_out
    roll_step = step_result.roll_step
    x_nom_next = list(step_result.x_new)
    dx_nom = list(step_result.dx_trans)

    sigma_actual_next_rad = float(step_result.sigma_actual_after_rad)
    roll_rate_next_rad_s = float(step_result.att_new.omega_b_rad_s[2])
    roll_accel_rad_s2 = (roll_rate_next_rad_s - roll_rate_rad_s) / float(dt_s)

    guidance_debug = dict(ctrl_out.guidance_debug)
    candidate_dicts = guidance_debug.get("candidate_dicts", [])
    selected_candidate_index = int(guidance_debug.get("selected_candidate_index", -1))

    if 0 <= selected_candidate_index < len(candidate_dicts):
        selected_candidate_debug = dict(candidate_dicts[selected_candidate_index])
    else:
        selected_candidate_debug = {}

    selected_interval_failure_kind = str(guidance_debug.get("selected_interval_failure_kind", ""))

    nominal_step_aero = nominal_aero_forces_from_state(
        x=x_nom,
        sigma_rad=float(sigma_actual_next_rad),
        params=params,
    )

    gravity_mps2 = float(constants.gravity(float(x_nom[0])))

    nominal_heat_info = nominal_heating_envelope_from_state(
        x_nominal=list(x_nom),
        dt_s=float(dt_s),
        heat_shield=nominal_heat_shield,
    )
    nominal_heat_shield = nominal_heat_info["heat_shield"]
    nominal_heat_qdot_max_interval = nominal_heat_info["qdot_max"]
    nominal_heat_Q_max_interval = nominal_heat_info["Q_max"]

    nominal_heat_violation = False
    nominal_heat_violation_reason = ""

    if nominal_heat_qdot_max_interval.hi > float(supervisor_cfg.heat_rate_limit):
        nominal_heat_violation = True
        nominal_heat_violation_reason = "nominal_heat_rate_limit"

    if nominal_heat_Q_max_interval.hi > float(supervisor_cfg.heat_load_limit):
        nominal_heat_violation = True
        if len(nominal_heat_violation_reason) > 0:
            nominal_heat_violation_reason += ","
        nominal_heat_violation_reason += "nominal_heat_load_limit"

    annotation = None
    interval_step_failure_kind = ""
    interval_step_failure_reason = ""
    interval_step_heat_violation = False
    interval_step_numerical_failure = False
    interval_step_recentered = False
    interval_step_recenter_reason = ""
    interval_step_split = False
    interval_step_split_reason = ""
    interval_step_split_depth = 0

    if interval_active:
        annotation = annotate_nominal_state_with_interval_supervisor(
            x_nominal_old=x_nom,
            params=params,
            sigma_actual_after_rad=float(sigma_actual_next_rad),
            x_interval_old=x_interval,
            supervisor_cfg=supervisor_cfg,
            dt_s=dt_s,
            heat_shield=interval_heat_shield,
            t_s=t_s,
        )

        interval_step_failure_kind = str(annotation.interval_failure_kind)
        interval_step_failure_reason = str(annotation.interval_failure_reason)
        interval_step_heat_violation = bool(annotation.interval_heat_violation)
        interval_step_numerical_failure = bool(annotation.interval_numerical_failure)
        interval_step_recentered = bool(annotation.recentered_this_step)
        interval_step_recenter_reason = str(getattr(annotation, "recenter_reason", ""))
        interval_step_split = bool(annotation.split_this_step)
        interval_step_split_reason = str(getattr(annotation, "split_reason", ""))
        interval_step_split_depth = int(annotation.split_depth_used)

        if annotation.last_valid_interval_state is not None:
            last_valid_interval_state = [Interval(iv.lo, iv.hi) for iv in annotation.last_valid_interval_state]
            last_valid_interval_time_s = float(t_s)

        if annotation.interval_valid and not annotation.interval_heat_violation:
            x_interval = [Interval(iv.lo, iv.hi) for iv in annotation.x_interval_new]
            interval_heat_shield = getattr(annotation, "heat_shield", interval_heat_shield)
        else:
            interval_failed = True
            interval_active = False
            interval_failure_time_s = float(t_s + dt_s)
            interval_failure_reason = str(annotation.interval_failure_reason)
            interval_failure_kind = str(annotation.interval_failure_kind)
            interval_heat_shield = None
            x_interval = None

    atm_nom = atmosphere_from_state_vector(x_nom)

    active_names = roll_step.fire_cmd.active_names()
    obs = ctrl_out.obs

    # --- Derived physics telemetry (g-load, energy, Mach, FMV) ---
    derived = telemetry.compute_derived_metrics(
        r_m=float(x_nom[0]),
        V_mps=float(x_nom[3]),
        altitude_m=float(atm_nom["alt_m"]),
        T_K=atm_nom.get("T_K", None),
        rho_kgm3=float(atm_nom["rho_kgm3"]),
        drag_mag_N=float(nominal_step_aero["drag_mag_N"]),
        lift_mag_N=float(nominal_step_aero["lift_mag_N"]),
        mass_kg=float(params["mass_kg"]),
        fmv=float(nominal_step_aero.get("FMV", float("nan"))),
    )

    # --- RCS propellant model ---
    # thruster-seconds this step = (active thrusters) * (fired on-time)
    fired_on_time_s = float(roll_step.num_fired_internal_steps) * float(constants.RCS_INTERNAL_DT_S)
    fired_thruster_seconds = float(len(active_names)) * fired_on_time_s
    rcs_fuel_rate_kg_s = telemetry.rcs_propellant_kg(fired_thruster_seconds) / float(dt_s) if dt_s > 0 else 0.0
    rcs_fuel_step_kg = telemetry.rcs_propellant_kg(fired_thruster_seconds)
    rcs_fuel_used_kg_total += rcs_fuel_step_kg

    # --- RL reward ingredients + composite default ---
    reward = telemetry.composite_reward(
        range_to_go_m=float(obs.get("range_to_go_m", np.nan)),
        qdot_w_m2=float(nominal_heat_qdot_max_interval.hi),
        load_factor_g=float(derived["load_factor_g"]),
        rcs_fuel_rate_kg_s=float(rcs_fuel_rate_kg_s),
    )

    candidate_sigma_deg_json = json.dumps([
        float(v.get("sigma_cmd_deg", np.nan))
        for v in candidate_dicts
    ])

    candidate_cost_json = json.dumps([
        float(v.get("total_cost", np.nan))
        for v in candidate_dicts
    ])

    candidate_heat_flag_json = json.dumps([
        int(bool(v.get("fully_feasible", False)))
        for v in candidate_dicts
    ])

    candidate_failure_reason_json = json.dumps([
        str(v.get("failure_reason", ""))
        for v in candidate_dicts
    ])

    interval_sigma_lo = np.nan
    interval_sigma_hi = np.nan
    interval_rho_lo = np.nan
    interval_rho_hi = np.nan
    interval_q_lo = np.nan
    interval_q_hi = np.nan
    interval_alt_lo = np.nan
    interval_alt_hi = np.nan
    width_r_m = np.nan
    width_phi_rad = np.nan
    width_lam_rad = np.nan
    width_V_mps = np.nan
    width_gamma_rad = np.nan
    width_chi_rad = np.nan
    dx_width_r = np.nan
    dx_width_phi = np.nan
    dx_width_lam = np.nan
    dx_width_V = np.nan
    dx_width_gamma = np.nan
    dx_width_chi = np.nan
    heating_qdot_max_lo = np.nan
    heating_qdot_max_hi = np.nan
    heating_qdot_mean_lo = np.nan
    heating_qdot_mean_hi = np.nan
    heating_Q_max_lo = np.nan
    heating_Q_max_hi = np.nan
    safety_status = "interval_inactive"
    safety_checks = "{}"
    layer_indices = ""
    interval_state_lo_hi = {}

    if annotation is not None:
        interval_sigma_lo = float(annotation.sigma_interval_used.lo)
        interval_sigma_hi = float(annotation.sigma_interval_used.hi)
        interval_rho_lo = float(annotation.rho_interval.lo)
        interval_rho_hi = float(annotation.rho_interval.hi)
        interval_q_lo = float(annotation.q_interval.lo)
        interval_q_hi = float(annotation.q_interval.hi)
        interval_alt_lo = float(annotation.altitude_interval.lo)
        interval_alt_hi = float(annotation.altitude_interval.hi)
        width_r_m = float(annotation.state_widths_new["r"])
        width_phi_rad = float(annotation.state_widths_new["phi"])
        width_lam_rad = float(annotation.state_widths_new["lam"])
        width_V_mps = float(annotation.state_widths_new["V"])
        width_gamma_rad = float(annotation.state_widths_new["gamma"])
        width_chi_rad = float(annotation.state_widths_new["chi"])
        dx_width_r = float(annotation.dx_widths["r"])
        dx_width_phi = float(annotation.dx_widths["phi"])
        dx_width_lam = float(annotation.dx_widths["lam"])
        dx_width_V = float(annotation.dx_widths["V"])
        dx_width_gamma = float(annotation.dx_widths["gamma"])
        dx_width_chi = float(annotation.dx_widths["chi"])
        safety_status = str(annotation.safety_status)
        safety_checks = str(annotation.safety_checks)
        layer_indices = ",".join(str(v) for v in annotation.layer_indices)

        if annotation.heating_qdot_max_interval is not None:
            heating_qdot_max_lo = float(annotation.heating_qdot_max_interval.lo)
            heating_qdot_max_hi = float(annotation.heating_qdot_max_interval.hi)
        if annotation.heating_qdot_mean_interval is not None:
            heating_qdot_mean_lo = float(annotation.heating_qdot_mean_interval.lo)
            heating_qdot_mean_hi = float(annotation.heating_qdot_mean_interval.hi)
        if annotation.heating_Q_max_interval is not None:
            heating_Q_max_lo = float(annotation.heating_Q_max_interval.lo)
            heating_Q_max_hi = float(annotation.heating_Q_max_interval.hi)

        for i, name in enumerate(STATE_NAMES):
            interval_state_lo_hi[f"{name}_lo"] = float(annotation.x_interval_new[i].lo)
            interval_state_lo_hi[f"{name}_hi"] = float(annotation.x_interval_new[i].hi)
            interval_state_lo_hi[f"{name}_width"] = float(annotation.x_interval_new[i].width())
    else:
        for name in STATE_NAMES:
            interval_state_lo_hi[f"{name}_lo"] = np.nan
            interval_state_lo_hi[f"{name}_hi"] = np.nan
            interval_state_lo_hi[f"{name}_width"] = np.nan

    nominal_continued_after_interval_failure = int(
        bool(interval_failed)
        and np.isfinite(interval_failure_time_s)
        and float(t_s) >= float(interval_failure_time_s)
        and not bool(nominal_heat_violation)
    )

    row = {
        "step": k,
        "t_s": t_s,
        "guidance_updated": bool(ctrl_out.guidance_updated),
        "trajectory_id": trajectory_id,
        "guidance_cycle_index": int(guidance_debug.get("guidance_cycle_index", k)),

        "sigma_cmd_rad": float(ctrl_out.sigma_cmd_rad),
        "sigma_target_rad": float(ctrl_out.sigma_target_rad),
        "sigma_actual_rad": float(sigma_actual_next_rad),
        "roll_rate_rad_s": float(roll_rate_next_rad_s),
        "roll_accel_rad_s2": float(roll_accel_rad_s2),

        "tau_roll_cmd_Nm": float(ctrl_out.tau_roll_cmd_Nm),
        "tau_roll_capacity_Nm": float(roll_step.tau_roll_capacity_Nm),
        "requested_duty": float(roll_step.requested_duty),
        "fired_this_step": int(bool(roll_step.fired_this_step)),
        "num_internal_steps": int(roll_step.num_internal_steps),
        "num_fired_internal_steps": int(roll_step.num_fired_internal_steps),
        "roll_pos_backlog_s": float(roll_step.roll_pos_backlog_s),
        "roll_neg_backlog_s": float(roll_step.roll_neg_backlog_s),
        "roll_pos_is_on": int(bool(roll_step.roll_pos_is_on)),
        "roll_neg_is_on": int(bool(roll_step.roll_neg_is_on)),
        "active_thrusters": ",".join(active_names),

        "torque_z_from_rcs": float(roll_step.wrench.torque_b_Nm[2]),
        "force_x_from_rcs": float(roll_step.wrench.force_b_N[0]),
        "force_y_from_rcs": float(roll_step.wrench.force_b_N[1]),
        "force_z_from_rcs": float(roll_step.wrench.force_b_N[2]),

        # CPAS (parachute) state and contribution
        "cpas_phase": str(cpas_out.phase),
        "cpas_drag_cdA_m2": float(cpas_out.drag_area_cdA_m2),
        "cpas_lift_scale": float(cpas_out.lift_scale),
        "cpas_open_fraction": float(cpas_out.open_fraction),
        "cpas_force_vertical": int(bool(cpas_out.force_vertical)),
        "cpas_events": ",".join(cpas_out.events) if cpas_out.events else "",
        "cpas_fbc_jettisoned": int(bool(cpas_out.fbc_jettisoned)),
        "cpas_mass_shed_kg": float(cpas_out.mass_shed_kg),
        "cpas_reefing_stage": int(cpas_out.reefing_stage),
        "cpas_num_mains_active": int(cpas_out.num_mains_active),
        "cpas_num_mains_squidding": int(cpas_out.num_mains_squidding),
        "cpas_pendulum_angle_deg": float(math.degrees(cpas_out.pendulum_angle_rad)),
        "cpas_pendulum_rate_deg_s": float(math.degrees(cpas_out.pendulum_rate_rad_s)),
        "cpas_pendulum_lateral_v_mps": float(cpas_out.pendulum_lateral_v_mps),
        "cpas_wind_east_mps": float(cpas_out.wind_east_mps),
        "cpas_wind_north_mps": float(cpas_out.wind_north_mps),
        "cpas_mass_kg_now": float(params["mass_kg"]),

        "r_m": float(x_nom[0]),
        "phi_rad": float(x_nom[1]),
        "lam_rad": float(x_nom[2]),
        "V_mps": float(x_nom[3]),
        "gamma_rad": float(x_nom[4]),
        "chi_rad": float(x_nom[5]),

        "alt_m": float(atm_nom["alt_m"]),
        "rho_kgm3": float(atm_nom["rho_kgm3"]),
        "q_pa": float(atm_nom["q_pa"]),
        "T_K": float(atm_nom["T_K"]) if not pd.isna(atm_nom["T_K"]) else np.nan,
        "p_Pa": float(atm_nom["p_Pa"]) if not pd.isna(atm_nom["p_Pa"]) else np.nan,

        "drag_mag_N": float(nominal_step_aero["drag_mag_N"]),
        "lift_mag_N": float(nominal_step_aero["lift_mag_N"]),
        "lift_vertical_plane_N": float(nominal_step_aero["lift_vertical_plane_N"]),
        "lift_lateral_left_N": float(nominal_step_aero["lift_lateral_left_N"]),
        # Aero coefficients (Bibb 2010 orion_cm_trim, blended-FMV lookup).
        # Logged here so the dashboard / aero_coefficients plot can render
        # the same CD/CL/LD trace the notebook driver already provides.
        "CD": float(nominal_step_aero["CD"]),
        "CL": float(nominal_step_aero["CL"]),
        "LD": float(nominal_step_aero["LD"]),
        "gravity_mps2": float(gravity_mps2),
        # --- Derived physics telemetry (shared telemetry.py) ---
        "mach": float(derived["mach"]),
        "speed_of_sound_mps": float(derived["speed_of_sound_mps"]),
        "fmv": float(derived["fmv"]),
        "specific_kinetic_e_j_kg": float(derived["specific_kinetic_e_j_kg"]),
        "specific_potential_e_j_kg": float(derived["specific_potential_e_j_kg"]),
        "specific_energy_j_kg": float(derived["specific_energy_j_kg"]),
        "energy_height_m": float(derived["energy_height_m"]),
        "aero_force_N": float(derived["aero_force_N"]),
        "aero_decel_mps2": float(derived["aero_decel_mps2"]),
        "drag_decel_mps2": float(derived["drag_decel_mps2"]),
        "load_factor_g": float(derived["load_factor_g"]),
        # --- RCS propellant ---
        "rcs_fuel_step_kg": float(rcs_fuel_step_kg),
        "rcs_fuel_rate_kg_s": float(rcs_fuel_rate_kg_s),
        "rcs_fuel_used_kg": float(rcs_fuel_used_kg_total),
        # --- RL reward ingredients + composite default ---
        "rew_range_term": float(reward["rew_range_term"]),
        "rew_heat_term": float(reward["rew_heat_term"]),
        "rew_gload_term": float(reward["rew_gload_term"]),
        "rew_fuel_term": float(reward["rew_fuel_term"]),
        "reward_default": float(reward["reward_default"]),

        "range_to_go_m": float(obs.get("range_to_go_m", np.nan)),
        "great_circle_range_m": float(obs.get("great_circle_range_m", np.nan)),
        "great_circle_angle_rad": float(obs.get("great_circle_angle_rad", np.nan)),
        "target_course_chi_rad": float(obs.get("target_course_chi_rad", np.nan)),
        "target_azimuth_north_rad": float(obs.get("target_azimuth_north_rad", np.nan)),
        "heading_error_rad": float(obs.get("heading_error_rad", np.nan)),
        "heading_offset_psi_rad": float(obs.get("heading_offset_psi_rad", np.nan)),
        "north_error_m": float(obs.get("north_error_m", np.nan)),
        "east_error_m": float(obs.get("east_error_m", np.nan)),
        "along_track_error_m": float(obs.get("along_track_error_m", np.nan)),
        "cross_track_error_m": float(obs.get("cross_track_error_m", np.nan)),

        "interval_active": int(bool(interval_active)),
        "interval_failed": int(bool(interval_failed)),
        "interval_failure_time_s": float(interval_failure_time_s) if np.isfinite(interval_failure_time_s) else np.nan,
        "interval_failure_reason": str(interval_failure_reason),
        "interval_failure_kind": str(interval_failure_kind),
        "nominal_continued_after_interval_failure": int(nominal_continued_after_interval_failure),
        "last_valid_interval_available": int(last_valid_interval_state is not None),
        "last_valid_interval_time_s": float(last_valid_interval_time_s) if np.isfinite(last_valid_interval_time_s) else np.nan,
        "interval_step_failure_kind": str(interval_step_failure_kind),
        "interval_step_failure_reason": str(interval_step_failure_reason),
        "interval_step_heat_violation": int(bool(interval_step_heat_violation)),
        "interval_step_numerical_failure": int(bool(interval_step_numerical_failure)),
        "interval_step_recentered": int(bool(interval_step_recentered)),
        "interval_step_recenter_reason": str(interval_step_recenter_reason),
        "interval_step_split": int(bool(interval_step_split)),
        "interval_step_split_reason": str(interval_step_split_reason),
        "interval_step_split_depth": int(interval_step_split_depth),

        "interval_sigma_lo": interval_sigma_lo,
        "interval_sigma_hi": interval_sigma_hi,
        "interval_rho_lo": interval_rho_lo,
        "interval_rho_hi": interval_rho_hi,
        "interval_q_lo": interval_q_lo,
        "interval_q_hi": interval_q_hi,
        "interval_alt_lo": interval_alt_lo,
        "interval_alt_hi": interval_alt_hi,

        "width_r_m": width_r_m,
        "width_phi_rad": width_phi_rad,
        "width_lam_rad": width_lam_rad,
        "width_V_mps": width_V_mps,
        "width_gamma_rad": width_gamma_rad,
        "width_chi_rad": width_chi_rad,

        "dx_width_r": dx_width_r,
        "dx_width_phi": dx_width_phi,
        "dx_width_lam": dx_width_lam,
        "dx_width_V": dx_width_V,
        "dx_width_gamma": dx_width_gamma,
        "dx_width_chi": dx_width_chi,

        "heating_qdot_max_lo": heating_qdot_max_lo,
        "heating_qdot_max_hi": heating_qdot_max_hi,
        "heating_qdot_mean_lo": heating_qdot_mean_lo,
        "heating_qdot_mean_hi": heating_qdot_mean_hi,
        "heating_Q_max_lo": heating_Q_max_lo,
        "heating_Q_max_hi": heating_Q_max_hi,

        "nominal_heat_qdot_max_lo": float(nominal_heat_qdot_max_interval.lo),
        "nominal_heat_qdot_max_hi": float(nominal_heat_qdot_max_interval.hi),
        "nominal_heat_Q_max_lo": float(nominal_heat_Q_max_interval.lo),
        "nominal_heat_Q_max_hi": float(nominal_heat_Q_max_interval.hi),
        "nominal_heat_violation": int(bool(nominal_heat_violation)),
        "nominal_heat_violation_reason": str(nominal_heat_violation_reason),
        # Stagnation heating component breakdown (convective vs radiative)
        # + radiative-equilibrium wall temperature.
        "qdot_conv_stag_hi": float(nominal_heat_info["qdot_conv_stag"].hi),
        "qdot_rad_stag_hi": float(nominal_heat_info["qdot_rad_stag"].hi),
        "qdot_total_stag_hi": float(nominal_heat_info["qdot_total_stag"].hi),
        "T_wall_stag_K_hi": float(nominal_heat_info["T_wall_stag_K"].hi),

        "guidance_selected_candidate_index": int(selected_candidate_index),
        "guidance_any_feasible_candidate": int(bool(guidance_debug.get("any_feasible_candidate", True))),
        "guidance_selected_candidate_heat_feasible": int(bool(guidance_debug.get("selected_candidate_heat_feasible", True))),
        "guidance_selected_interval_screen_used": int(bool(guidance_debug.get("selected_interval_screen_used", True))),
        "guidance_selected_interval_failure_kind": str(selected_interval_failure_kind),
        "guidance_chosen_sigma_cmd_deg": float(guidance_debug.get("chosen_sigma_cmd_deg", np.nan)),
        "guidance_chosen_sigma_mag_deg": float(guidance_debug.get("chosen_sigma_mag_deg", np.nan)),
        "guidance_heading_deadband_deg": float(guidance_debug.get("heading_deadband_deg", np.nan)),
        "guidance_current_altitude_rate_mps": float(guidance_debug.get("current_altitude_rate_mps", np.nan)),
        "guidance_selected_geometry_cost": float(guidance_debug.get("selected_geometry_cost", np.nan)),
        "guidance_selected_heat_penalty": float(guidance_debug.get("selected_heat_penalty", np.nan)),
        "guidance_selected_total_cost": float(guidance_debug.get("selected_total_cost", np.nan)),
        "guidance_selected_failure_reason": str(guidance_debug.get("selected_failure_reason", "")),
        "guidance_selected_violation_amount": float(guidance_debug.get("selected_violation_amount", np.nan)),
        "guidance_selected_first_violation_step": int(guidance_debug.get("selected_first_violation_step", -1)),
        "guidance_selected_first_violation_time_s": float(guidance_debug.get("selected_first_violation_time_s", np.nan)),
        "guidance_selected_max_heating_rate_lo": float(guidance_debug.get("selected_max_heating_rate_lo", np.nan)),
        "guidance_selected_max_heating_rate_hi": float(guidance_debug.get("selected_max_heating_rate_hi", np.nan)),
        "guidance_selected_max_heat_load_lo": float(guidance_debug.get("selected_max_heat_load_lo", np.nan)),
        "guidance_selected_max_heat_load_hi": float(guidance_debug.get("selected_max_heat_load_hi", np.nan)),
        "guidance_candidate_count": int(len(candidate_dicts)),
        "guidance_candidate_sigma_deg_json": candidate_sigma_deg_json,
        "guidance_candidate_cost_json": candidate_cost_json,
        "guidance_candidate_heat_flag_json": candidate_heat_flag_json,
        "guidance_candidate_failure_reason_json": candidate_failure_reason_json,

        "safety_status": safety_status,
        "safety_checks": safety_checks,
        "layer_indices": layer_indices,
    }

    row.update(interval_state_lo_hi)
    rows.append(row)

    if bool(ctrl_out.guidance_updated):
        selected_heat_feasible = bool(guidance_debug.get("selected_candidate_heat_feasible", True))
        selected_is_true_heat_failure = (not selected_heat_feasible) and str(selected_interval_failure_kind) != "numerical"

        if selected_is_true_heat_failure:
            interval_propagation_invalid = int(str(selected_interval_failure_kind) == "numerical")
            failed_heat_rows.append(
                {
                    "trajectory_id": trajectory_id,
                    "step": int(k),
                    "guidance_cycle_index": int(guidance_debug.get("guidance_cycle_index", k)),
                    "time_s": float(t_s),
                    "r_m": float(x_nom[0]),
                    "phi_rad": float(x_nom[1]),
                    "lam_rad": float(x_nom[2]),
                    "V_mps": float(x_nom[3]),
                    "gamma_rad": float(x_nom[4]),
                    "chi_rad": float(x_nom[5]),
                    "target_phi_rad": float(target_phi_rad),
                    "target_lam_rad": float(target_lam_rad),
                    "sigma_cmd_candidate_deg_json": candidate_sigma_deg_json,
                    "sigma_cmd_candidate_cost_json": candidate_cost_json,
                    "sigma_cmd_candidate_heat_flag_json": candidate_heat_flag_json,
                    "sigma_cmd_candidate_failure_reason_json": candidate_failure_reason_json,
                    "selected_sigma_cmd_rad": float(ctrl_out.sigma_cmd_rad),
                    "selected_sigma_cmd_deg": float(guidance_debug.get("chosen_sigma_cmd_deg", np.nan)),
                    "sigma_actual_rad": float(sigma_actual_rad),
                    "range_to_go_m": float(obs.get("range_to_go_m", np.nan)),
                    "heading_error_rad": float(obs.get("heading_error_rad", np.nan)),
                    "cross_track_error_m": float(obs.get("cross_track_error_m", np.nan)),
                    "dynamic_pressure_pa": float(obs.get("dynamic_pressure_pa", np.nan)),
                    "density_kgm3": float(obs.get("density_kgm3", np.nan)),
                    "predicted_max_heating_rate_lo": float(guidance_debug.get("selected_max_heating_rate_lo", np.nan)),
                    "predicted_max_heating_rate_hi": float(guidance_debug.get("selected_max_heating_rate_hi", np.nan)),
                    "predicted_max_heat_load_lo": float(guidance_debug.get("selected_max_heat_load_lo", np.nan)),
                    "predicted_max_heat_load_hi": float(guidance_debug.get("selected_max_heat_load_hi", np.nan)),
                    "heat_feasible_flag": int(False),
                    "violation_amount": float(guidance_debug.get("selected_violation_amount", np.nan)),
                    "first_violation_step": int(guidance_debug.get("selected_first_violation_step", -1)),
                    "first_violation_time_s": float(guidance_debug.get("selected_first_violation_time_s", np.nan)),
                    "failure_reason": str(guidance_debug.get("selected_failure_reason", "")),
                    "interval_propagation_invalid": int(interval_propagation_invalid),
                    "geometry_cost": float(guidance_debug.get("selected_geometry_cost", np.nan)),
                    "heat_penalty": float(guidance_debug.get("selected_heat_penalty", np.nan)),
                    "total_cost": float(guidance_debug.get("selected_total_cost", np.nan)),
                }
            )

    # Apply wind drift and pendulum lateral velocity under chutes
    # (shared helper -- the notebook driver uses the same call).
    cpas.apply_horizontal_drift_to_state(x_nom_next, cpas_out, dt_s, _R_EARTH)

    x_nom = list(x_nom_next)
    att_nom = step_result.att_new

    if nominal_heat_violation:
        sim_terminated_by_nominal_heat = True
        termination_reason = nominal_heat_violation_reason
        print(f"Stopping at t = {t_s + dt_s:.2f} s because nominal heating violated the hard limit")
        break

# build the main logged trajectory dataframe
df = pd.DataFrame(rows)

# build the failed heat dataframe with fixed columns
failed_heat_steps_df = pd.DataFrame(failed_heat_rows, columns=failed_heat_columns)

# --- Early persistence ------------------------------------------------------
# Save the most important outputs immediately after the dataframe is built,
# so an exception in the (large) downstream plotting/diagnostic section does
# not lose the actual simulation data. The full save block lower in the file
# still runs on success and will simply overwrite these files.
df.to_csv(success_trajectory_csv_path, index=False)
print("Saved (early)", success_trajectory_csv_path)
failed_heat_steps_df.to_csv(failed_heat_step_csv_path, index=False)
print("Saved (early)", failed_heat_step_csv_path)

_cpas_events_df_early = (
    pd.DataFrame(cpas_events_log)
    if cpas_events_log
    else pd.DataFrame(columns=["t_s", "step", "alt_m", "V_mps", "mach", "phase", "event"])
)
_cpas_events_df_early.to_csv(cpas_events_csv_path, index=False)
print("Saved (early)", cpas_events_csv_path)

_cpas_summary_early = cpas_inst.summary()
_run_summary_early = {
    "trajectory_id": trajectory_id,
    "num_logged_steps": int(len(df)),
    "termination_reason": termination_reason,
    "final_t_s": float(df["t_s"].iloc[-1]) if len(df) else 0.0,
    "final_alt_m": float(df["alt_m"].iloc[-1]) if len(df) else 0.0,
    "final_V_mps": float(df["V_mps"].iloc[-1]) if len(df) else 0.0,
    "final_gamma_deg": float(math.degrees(df["gamma_rad"].iloc[-1])) if len(df) else 0.0,
    "aero_model": params["aero_model"],
    "cpas_summary": _cpas_summary_early,
    "cpas_drogue_deploy_t_s": _cpas_summary_early["drogue_deploy_t_s"],
    "cpas_pilot_deploy_t_s": _cpas_summary_early["pilot_deploy_t_s"],
    "cpas_main_deploy_t_s": _cpas_summary_early["main_deploy_t_s"],
    "cpas_landed_t_s": _cpas_summary_early["landed_t_s"],
    "output_files": {
        "trajectory_csv": success_trajectory_csv_path,
        "failed_heat_step_csv": failed_heat_step_csv_path,
        "cpas_events_csv": cpas_events_csv_path,
    },
}
with open(run_summary_json_path, "w", encoding="utf-8") as _f:
    json.dump(_run_summary_early, _f, indent=2, default=str)
print("Saved (early)", run_summary_json_path)

_landing_summary_early = {
    "termination_reason": termination_reason,
    "final_t_s": float(df["t_s"].iloc[-1]) if len(df) else 0.0,
    "final_alt_km": float(df["alt_m"].iloc[-1]) / 1000.0 if len(df) else 0.0,
    "final_V_mps": float(df["V_mps"].iloc[-1]) if len(df) else 0.0,
    "final_gamma_deg": float(math.degrees(df["gamma_rad"].iloc[-1])) if len(df) else 0.0,
    "final_chi_deg": float(math.degrees(df["chi_rad"].iloc[-1])) if len(df) else 0.0,
    "final_phi_deg": float(math.degrees(df["phi_rad"].iloc[-1])) if len(df) else 0.0,
    "final_lam_deg": float(math.degrees(df["lam_rad"].iloc[-1])) if len(df) else 0.0,
    "cpas_drogue_deploy_t_s": _cpas_summary_early["drogue_deploy_t_s"],
    "cpas_pilot_deploy_t_s": _cpas_summary_early["pilot_deploy_t_s"],
    "cpas_main_deploy_t_s": _cpas_summary_early["main_deploy_t_s"],
    "cpas_landed_t_s": _cpas_summary_early["landed_t_s"],
}
with open(landing_summary_json_path, "w", encoding="utf-8") as _f:
    json.dump(_landing_summary_early, _f, indent=2, default=str)
print("Saved (early)", landing_summary_json_path)
# ---------------------------------------------------------------------------

# --- Unified plotting via plotting.py --------------------------------------
# Both the notebook and run_sim.py call into this module so the two parallel
# drivers produce the same set of named figures.
PY_FIG_DIR = PY_OUTPUT_DIR / "figures"
PY_FIG_DIR.mkdir(parents=True, exist_ok=True)

# Build a simple thruster_fires dataframe from the active_thrusters column
# so plot_thruster_raster gets the same input shape as in the notebook.
_fire_records = []
if "active_thrusters" in df.columns:
    for _row in df.itertuples(index=False):
        _names_str = getattr(_row, "active_thrusters", "")
        if _names_str:
            for _n in str(_names_str).split(","):
                _n = _n.strip()
                if _n:
                    _fire_records.append({
                        "step": int(getattr(_row, "step", -1)),
                        "t_s": float(getattr(_row, "t_s", 0.0)),
                        "thruster": _n,
                    })
thruster_fires_df = pd.DataFrame(_fire_records)

def _run_sim_save_fig(name):
    path = PY_FIG_DIR / f"{name}.png"
    plt.savefig(path, bbox_inches="tight")
    print(f"  saved {path}")

import os as _os_for_plot_categories
_pc_env = _os_for_plot_categories.environ.get("PLOT_CATEGORIES", "").strip()
PLOT_CATEGORIES = set(c.strip() for c in _pc_env.split(",") if c.strip()) if _pc_env else None
if PLOT_CATEGORIES:
    print(f"Rendering figures into {PY_FIG_DIR} -- filtered to categories: {sorted(PLOT_CATEGORIES)}")
else:
    print(f"Rendering figures into {PY_FIG_DIR} (all categories)")
plotting.render_all_figures(
    df,
    _run_sim_save_fig,
    only=PLOT_CATEGORIES,
    target_phi_rad=float(target_phi_rad),
    target_lam_rad=float(target_lam_rad),
    thruster_fires_df=thruster_fires_df,
    rcs_system=rcs_system,
    nominal_heat_shield=nominal_heat_shield,
    dt_s=float(dt_s),
)
print("Plotting complete.")

# Stop the script here; everything below was the legacy plotting/analysis
# block that lived in this file when it was first converted from a notebook.
# The shared plotting.py is now the source of truth.
import sys as _sys
_sys.exit(0)
# ---------------------------------------------------------------------------

interval_failure_detected = bool(len(df) > 0 and df["interval_failed"].max() > 0)
first_interval_failure_time_s = np.nan
interval_failure_kind_summary = ""
interval_failure_reason_summary = ""
nominal_continuation_duration_s = 0.0
nominal_continuation_steps = 0

if interval_failure_detected:
    failure_rows = df[df["interval_failed"] == 1]
    first_failure_row = failure_rows.iloc[0]
    first_interval_failure_time_s = float(first_failure_row["interval_failure_time_s"])
    interval_failure_kind_summary = str(first_failure_row["interval_failure_kind"])
    interval_failure_reason_summary = str(first_failure_row["interval_failure_reason"])
    nominal_continuation_duration_s = max(0.0, float(df["t_s"].iloc[-1]) - float(first_interval_failure_time_s))
    nominal_continuation_steps = int((df["nominal_continued_after_interval_failure"] == 1).sum())

interval_recenter_count = int((df["interval_step_recentered"] == 1).sum()) if len(df) > 0 else 0
interval_recenter_width_count = int(((df["interval_step_recentered"] == 1) & (df["interval_step_recenter_reason"] == "width_threshold")).sum()) if len(df) > 0 else 0
interval_recenter_cadence_count = int(((df["interval_step_recentered"] == 1) & (df["interval_step_recenter_reason"] == "cadence")).sum()) if len(df) > 0 else 0

interval_split_count = int((df["interval_step_split"] == 1).sum()) if len(df) > 0 else 0
interval_split_width_count = int(((df["interval_step_split"] == 1) & (df["interval_step_split_reason"] == "width_threshold")).sum()) if len(df) > 0 else 0
interval_split_denominator_v_count = int(((df["interval_step_split"] == 1) & (df["interval_step_split_reason"] == "denominator_V")).sum()) if len(df) > 0 else 0
interval_split_denominator_cos_gamma_count = int(((df["interval_step_split"] == 1) & (df["interval_step_split_reason"] == "denominator_cos_gamma")).sum()) if len(df) > 0 else 0
interval_split_denominator_cos_phi_count = int(((df["interval_step_split"] == 1) & (df["interval_step_split_reason"] == "denominator_cos_phi")).sum()) if len(df) > 0 else 0

failed_episode_summary_df = pd.DataFrame(
    [
        {
            "trajectory_id": trajectory_id,
            "num_logged_steps": int(len(df)),
            "num_heat_infeasible_guidance_cycles": int(len(failed_heat_steps_df)),
            "any_heat_infeasible_guidance_cycle": int(len(failed_heat_steps_df) > 0),
            "interval_failure_detected": int(interval_failure_detected),
            "interval_failure_time_s": float(first_interval_failure_time_s) if np.isfinite(first_interval_failure_time_s) else np.nan,
            "interval_failure_kind": interval_failure_kind_summary,
            "interval_failure_reason": interval_failure_reason_summary,
            "nominal_continuation_duration_s": float(nominal_continuation_duration_s),
            "nominal_continuation_steps": int(nominal_continuation_steps),
            "interval_recenter_count": int(interval_recenter_count),
            "interval_recenter_width_count": int(interval_recenter_width_count),
            "interval_recenter_cadence_count": int(interval_recenter_cadence_count),
            "interval_split_count": int(interval_split_count),
            "interval_split_width_count": int(interval_split_width_count),
            "interval_split_denominator_v_count": int(interval_split_denominator_v_count),
            "interval_split_denominator_cos_gamma_count": int(interval_split_denominator_cos_gamma_count),
            "interval_split_denominator_cos_phi_count": int(interval_split_denominator_cos_phi_count),
            "sim_terminated_by_nominal_heat": int(bool(sim_terminated_by_nominal_heat)),
            "termination_reason": str(termination_reason),
            "success_csv_path": success_trajectory_csv_path,
            "failed_heat_step_csv_path": failed_heat_step_csv_path,
        }
    ]
)

print("Number of logged steps", len(df))
print("Number of logged heat infeasible guidance cycles", len(failed_heat_steps_df))
print("Interval failure detected", interval_failure_detected)
print("Simulation terminated by nominal heat", bool(sim_terminated_by_nominal_heat))
print("Termination reason", termination_reason)
print("Interval recenter count", interval_recenter_count)
print("Interval recenter width count", interval_recenter_width_count)
print("Interval recenter cadence count", interval_recenter_cadence_count)
print("Interval split count", interval_split_count)
print("Interval split width count", interval_split_width_count)
print("Interval split denominator v count", interval_split_denominator_v_count)
print("Interval split denominator cos gamma count", interval_split_denominator_cos_gamma_count)
print("Interval split denominator cos phi count", interval_split_denominator_cos_phi_count)

if interval_failure_detected:
    print("First interval failure time s", first_interval_failure_time_s)
    print("Interval failure kind", interval_failure_kind_summary)
    print("Interval failure reason", interval_failure_reason_summary)
    print("Nominal continuation duration after interval failure s", nominal_continuation_duration_s)
    print("Nominal continuation steps after interval failure", nominal_continuation_steps)

if len(df) > 0:
    last = df.iloc[-1]
    print("Last logged time s", last["t_s"])
    print("Last logged altitude m", last["alt_m"])
    print("Last logged speed m per s", last["V_mps"])
    print("Last logged gamma deg", deg(last["gamma_rad"]))
    print("Last logged sigma actual deg", deg(last["sigma_actual_rad"]))
    print("Last logged roll rate deg per s", deg(last["roll_rate_rad_s"]))
    print("Last logged heading error deg", deg(last["heading_error_rad"]))
    print("Last logged cross track m", last["cross_track_error_m"])
    print("Last logged range to go m", last["range_to_go_m"])
    print("Last selected predictor cost", last["guidance_selected_total_cost"])
    print("Last selected heat feasible flag", last["guidance_selected_candidate_heat_feasible"])
    print("Last nominal heat qdot hi", last["nominal_heat_qdot_max_hi"])
    print("Last nominal heat Q hi", last["nominal_heat_Q_max_hi"])
    print("Last interval q high Pa", last["interval_q_hi"])
    print("Last interval density high", last["interval_rho_hi"])
    print("Last interval active flag", last["interval_active"])
    print("Last interval recenter reason", last["interval_step_recenter_reason"])
    print("Last interval split reason", last["interval_step_split_reason"])

df.head()

# Guidance candidate choice diagnostics
# Put this cell right after the DataFrame build cell

guidance_choice_df = df[df["guidance_updated"]].copy()

if len(guidance_choice_df) == 0:
    print("No guidance update rows were logged")
else:
    print("Number of guidance updates", len(guidance_choice_df))

    display(
        guidance_choice_df[
            [
                "t_s",
                "guidance_selected_candidate_index",
                "guidance_chosen_sigma_mag_deg",
                "guidance_chosen_sigma_cmd_deg",
                "guidance_selected_total_cost",
                "guidance_selected_candidate_heat_feasible",
            ]
        ].head(20)
    )

    chosen_counts = (
        guidance_choice_df["guidance_chosen_sigma_mag_deg"]
        .value_counts(dropna=False)
        .sort_index()
    )

    print("Chosen candidate counts")
    display(chosen_counts.rename("count").to_frame())

    plt.figure(figsize=(11, 4))
    plt.plot(
        guidance_choice_df["t_s"],
        guidance_choice_df["guidance_chosen_sigma_mag_deg"],
        marker="o",
        linewidth=1.5,
    )
    plt.xlabel("Time s")
    plt.ylabel("Chosen bank magnitude deg")
    plt.title("Chosen predictor candidate over guidance updates")
    plt.grid(True)
    plt.show()

    plt.figure(figsize=(8, 4))
    chosen_counts.plot(kind="bar")
    plt.xlabel("Chosen bank magnitude deg")
    plt.ylabel("Count")
    plt.title("Distribution of chosen predictor candidates")
    plt.grid(True, axis="y")
    plt.show()

# Quick summary table of the conditional uncertainty outputs

cols = [
    "t_s",
    "sigma_cmd_rad",
    "sigma_target_rad",
    "sigma_actual_rad",
    "roll_rate_rad_s",
    "tau_roll_cmd_Nm",
    "requested_duty",
    "fired_this_step",
    "alt_m",
    "interval_alt_lo",
    "interval_alt_hi",
    "V_mps",
    "V_mps_lo",
    "V_mps_hi",
    "interval_rho_lo",
    "interval_rho_hi",
    "interval_q_lo",
    "interval_q_hi",
    "T_K",
    "p_Pa",
    "safety_status",
]

df[cols].head(12)

# Ground track, target comparison, and landing error summary

if len(df) == 0:
    raise ValueError("The dataframe is empty. Run the main loop cell first.")

# Convert trajectory states to degrees for plotting
track_df = df.copy()
track_df["phi_deg"] = np.degrees(track_df["phi_rad"])
track_df["lam_deg"] = np.degrees(track_df["lam_rad"])

# Start and end points from the nominal path
start_phi_deg = float(track_df["phi_deg"].iloc[0])
start_lam_deg = float(track_df["lam_deg"].iloc[0])

end_phi_deg = float(track_df["phi_deg"].iloc[-1])
end_lam_deg = float(track_df["lam_deg"].iloc[-1])

# Target point in degrees
target_phi_deg = float(np.degrees(target_phi_rad))
target_lam_deg = float(np.degrees(target_lam_rad))

# Final nominal state error relative to target
R_earth_m = float(constants.RADIUS_EARTH)

dphi_rad = float(track_df["phi_rad"].iloc[-1] - target_phi_rad)
dlam_rad = float(track_df["lam_rad"].iloc[-1] - target_lam_rad)

north_error_m = R_earth_m * dphi_rad
east_error_m = R_earth_m * dlam_rad * math.cos(target_phi_rad)
ground_error_m = math.sqrt(north_error_m**2 + east_error_m**2)

# Great circle style angular distance and surface distance
sin_dphi_2 = math.sin(dphi_rad / 2.0)
sin_dlam_2 = math.sin(dlam_rad / 2.0)
a = (
    sin_dphi_2**2
    + math.cos(float(track_df["phi_rad"].iloc[-1])) * math.cos(target_phi_rad) * sin_dlam_2**2
)
a = min(1.0, max(0.0, a))
central_angle_rad = 2.0 * math.asin(math.sqrt(a))
great_circle_error_m = R_earth_m * central_angle_rad

# Plot the nominal ground track and mark start, end, and target
plt.figure(figsize=(10, 7))
plt.plot(track_df["lam_deg"], track_df["phi_deg"], label="Nominal ground track")
plt.scatter(start_lam_deg, start_phi_deg, s=100, marker="o", label="Start")
plt.scatter(end_lam_deg, end_phi_deg, s=120, marker="s", label="End")
plt.scatter(target_lam_deg, target_phi_deg, s=180, marker="x", linewidths=3, label="Target")

plt.xlabel("Longitude deg")
plt.ylabel("Latitude deg")
plt.title("Ground track with start, end, and target")
plt.legend()
plt.grid(True)
plt.axis("equal")
plt.show()

# Optional zoomed view around the target and end point
lon_vals = [start_lam_deg, end_lam_deg, target_lam_deg]
lat_vals = [start_phi_deg, end_phi_deg, target_phi_deg]

lon_pad = max(0.5, 0.15 * (max(lon_vals) - min(lon_vals) + 1.0))
lat_pad = max(0.5, 0.15 * (max(lat_vals) - min(lat_vals) + 1.0))

# Print a compact landing style summary
print("Start latitude deg", start_phi_deg)
print("Start longitude deg", start_lam_deg)
print("End latitude deg", end_phi_deg)
print("End longitude deg", end_lam_deg)
print("Target latitude deg", target_phi_deg)
print("Target longitude deg", target_lam_deg)
print("North error m", north_error_m)
print("East error m", east_error_m)
print("Planar ground error m", ground_error_m)
print("Great circle surface error m", great_circle_error_m)

# Tabular summary
landing_summary_df = pd.DataFrame(
    [
        {
            "start_lat_deg": start_phi_deg,
            "start_lon_deg": start_lam_deg,
            "end_lat_deg": end_phi_deg,
            "end_lon_deg": end_lam_deg,
            "target_lat_deg": target_phi_deg,
            "target_lon_deg": target_lam_deg,
            "north_error_m": north_error_m,
            "east_error_m": east_error_m,
            "planar_ground_error_m": ground_error_m,
            "great_circle_error_m": great_circle_error_m,
        }
    ]
)

landing_summary_df

# Approximate total nominal path length traveled

if len(df) < 2:
    raise ValueError("Need at least two trajectory rows in df")

R_earth_m = float(constants.RADIUS_EARTH)

segment_lengths_m = []

for i in range(1, len(df)):
    r1 = float(df.iloc[i - 1]["r_m"])
    r2 = float(df.iloc[i]["r_m"])

    phi1 = float(df.iloc[i - 1]["phi_rad"])
    phi2 = float(df.iloc[i]["phi_rad"])

    lam1 = float(df.iloc[i - 1]["lam_rad"])
    lam2 = float(df.iloc[i]["lam_rad"])

    x1 = r1 * math.cos(phi1) * math.cos(lam1)
    y1 = r1 * math.cos(phi1) * math.sin(lam1)
    z1 = r1 * math.sin(phi1)

    x2 = r2 * math.cos(phi2) * math.cos(lam2)
    y2 = r2 * math.cos(phi2) * math.sin(lam2)
    z2 = r2 * math.sin(phi2)

    ds = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    segment_lengths_m.append(ds)

total_distance_m = float(np.sum(segment_lengths_m))
total_distance_km = total_distance_m / 1000.0

print("Approximate total nominal path length m", total_distance_m)
print("Approximate total nominal path length km", total_distance_km)

# 3D trajectory from four viewing angles with conservative interval corner curves

from mpl_toolkits.mplot3d import Axes3D

def spherical_to_cartesian(r, phi, lam):
    x = r * np.cos(phi) * np.cos(lam)
    y = r * np.cos(phi) * np.sin(lam)
    z = r * np.sin(phi)
    return x, y, z

# Nominal trajectory
r_nom = df["r_m"].to_numpy(dtype=float)
phi_nom = df["phi_rad"].to_numpy(dtype=float)
lam_nom = df["lam_rad"].to_numpy(dtype=float)

x_nom, y_nom, z_nom = spherical_to_cartesian(r_nom, phi_nom, lam_nom)

# Conservative interval corner trajectories
r_lo = df["r_m_lo"].to_numpy(dtype=float)
r_hi = df["r_m_hi"].to_numpy(dtype=float)

phi_lo = df["phi_rad_lo"].to_numpy(dtype=float)
phi_hi = df["phi_rad_hi"].to_numpy(dtype=float)

lam_lo = df["lam_rad_lo"].to_numpy(dtype=float)
lam_hi = df["lam_rad_hi"].to_numpy(dtype=float)

x_lo, y_lo, z_lo = spherical_to_cartesian(r_lo, phi_lo, lam_lo)
x_hi, y_hi, z_hi = spherical_to_cartesian(r_hi, phi_hi, lam_hi)

# Earth sphere for reference
u = np.linspace(0.0, 2.0 * np.pi, 60)
v = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 30)
uu, vv = np.meshgrid(u, v)

x_e = R_earth_m * np.cos(vv) * np.cos(uu)
y_e = R_earth_m * np.cos(vv) * np.sin(uu)
z_e = R_earth_m * np.sin(vv)

fig = plt.figure(figsize=(14, 12))
views = [(20, 35), (20, 125), (50, 215), (80, 300)]

for idx, (elev, azim) in enumerate(views, start=1):
    ax = fig.add_subplot(2, 2, idx, projection="3d")

    ax.plot_surface(x_e, y_e, z_e, alpha=0.15, linewidth=0)

    ax.plot(x_nom, y_nom, z_nom, label="Nominal trajectory")
    ax.plot(x_lo, y_lo, z_lo, linestyle="--", label="Interval low corner")
    ax.plot(x_hi, y_hi, z_hi, linestyle="--", label="Interval high corner")

    ax.scatter(x_nom[0], y_nom[0], z_nom[0], marker="o", s=40, label="Start" if idx == 1 else None)
    ax.scatter(x_nom[-1], y_nom[-1], z_nom[-1], marker="s", s=40, label="End" if idx == 1 else None)

    ax.set_xlabel("X m")
    ax.set_ylabel("Y m")
    ax.set_zlabel("Z m")
    ax.set_title(f"3D trajectory view {idx}")
    ax.view_init(elev=elev, azim=azim)
    ax.set_box_aspect((1, 1, 1))

    max_extent = np.max(np.abs(np.concatenate([x_nom, y_nom, z_nom, x_hi, y_hi, z_hi])))
    ax.set_xlim(-max_extent, max_extent)
    ax.set_ylim(-max_extent, max_extent)
    ax.set_zlim(-max_extent, max_extent)

    if idx == 1:
        ax.legend()

plt.tight_layout()
plt.show()

# 3D descent tube view centered on the trajectory instead of the whole Earth

# A 3D path needs at least two trajectory samples so there is a start and an end to draw.
if len(df) < 2:
    raise ValueError("Need at least two rows in df to plot the 3D trajectory")

# Read the Earth radius once so the latitude and longitude offsets can be converted into meters.
R_earth_m = float(constants.RADIUS_EARTH)

# Use the first trajectory point as the local reference origin for the tangent plane view.
phi0 = float(df["phi_rad"].iloc[0])
lam0 = float(df["lam_rad"].iloc[0])

# This helper maps latitude, longitude, and altitude into a local east north up frame.
# The approximation is local, so it is most useful when the path is viewed relative to the start point.
def local_track_from_lat_lon_alt(phi_rad, lam_rad, alt_m, phi_ref_rad, lam_ref_rad, R_m):
    # East distance comes from longitude change scaled by the cosine of the reference latitude.
    east_m = R_m * (lam_rad - lam_ref_rad) * math.cos(phi_ref_rad)

    # North distance comes directly from latitude change.
    north_m = R_m * (phi_rad - phi_ref_rad)

    # Up is just the altitude above the surface.
    up_m = alt_m

    # Return the local coordinates in meters.
    return east_m, north_m, up_m

# Store the nominal trajectory in local east north up coordinates.
east_nom = []
north_nom = []
up_nom = []

# Store the lower interval corner path so the uncertainty envelope can be drawn.
east_lo = []
north_lo = []
up_lo = []

# Store the upper interval corner path for the other side of the envelope.
east_hi = []
north_hi = []
up_hi = []

# Walk through every row so the nominal path and both interval corners are converted into the same local frame.
for _, row in df.iterrows():
    # Convert the nominal latitude, longitude, and altitude at this time step.
    e, n, u = local_track_from_lat_lon_alt(
        float(row["phi_rad"]),
        float(row["lam_rad"]),
        float(row["alt_m"]),
        phi0,
        lam0,
        R_earth_m,
    )

    # Save the nominal local coordinates.
    east_nom.append(e)
    north_nom.append(n)
    up_nom.append(u)

    # Convert the lower bound trajectory values into the same local frame.
    e_lo, n_lo, u_lo = local_track_from_lat_lon_alt(
        float(row["phi_rad_lo"]),
        float(row["lam_rad_lo"]),
        float(row["interval_alt_lo"]),
        phi0,
        lam0,
        R_earth_m,
    )

    # Save the lower interval corner coordinates.
    east_lo.append(e_lo)
    north_lo.append(n_lo)
    up_lo.append(u_lo)

    # Convert the upper bound trajectory values into the same local frame.
    e_hi, n_hi, u_hi = local_track_from_lat_lon_alt(
        float(row["phi_rad_hi"]),
        float(row["lam_rad_hi"]),
        float(row["interval_alt_hi"]),
        phi0,
        lam0,
        R_earth_m,
    )

    # Save the upper interval corner coordinates.
    east_hi.append(e_hi)
    north_hi.append(n_hi)
    up_hi.append(u_hi)

# Convert the nominal coordinate lists into numpy arrays for plotting.
east_nom = np.array(east_nom)
north_nom = np.array(north_nom)
up_nom = np.array(up_nom)

# Convert the lower interval corner lists into numpy arrays.
east_lo = np.array(east_lo)
north_lo = np.array(north_lo)
up_lo = np.array(up_lo)

# Convert the upper interval corner lists into numpy arrays.
east_hi = np.array(east_hi)
north_hi = np.array(north_hi)
up_hi = np.array(up_hi)

# Import the 3D plotting toolkit used by matplotlib.
from mpl_toolkits.mplot3d import Axes3D

# Create a larger figure so the 3D views are easier to read.
fig = plt.figure(figsize=(14, 12))

# Define several camera angles, even though only the first two will be shown here.
views = [(20, 35), (20, 125), (35, 215), (70, 300)]

# Only show the first two viewing angles so the figure stays compact.
for idx, (elev, azim) in enumerate(views[:2], start=1):
    # Create a 3D subplot for this view.
    ax = fig.add_subplot(2, 2, idx, projection="3d")

    # Draw the nominal trajectory as the center line of the descent path.
    ax.plot(east_nom, north_nom, up_nom, linewidth=2, label="Nominal trajectory")

    # Draw the lower interval corner as one side of the uncertainty envelope.
    ax.plot(east_lo, north_lo, up_lo, linestyle="--", label="Interval low corner")

    # Draw the upper interval corner as the other side of the uncertainty envelope.
    ax.plot(east_hi, north_hi, up_hi, linestyle="--", label="Interval high corner")

    # Choose a stride so only some cross connections are drawn.
    # This keeps the envelope readable instead of cluttered.
    stride = max(1, len(df) // 30)

    # Connect matching low and high corner points so the interval looks more like a tube.
    for i in range(0, len(df), stride):
        ax.plot(
            [east_lo[i], east_hi[i]],
            [north_lo[i], north_hi[i]],
            [up_lo[i], up_hi[i]],
            alpha=0.5,
        )

    # Mark the first nominal point so the start of the descent is obvious.
    ax.scatter(east_nom[0], north_nom[0], up_nom[0], s=60, marker="o", label="Start" if idx == 1 else None)

    # Mark the last nominal point so the end of the descent is obvious.
    ax.scatter(east_nom[-1], north_nom[-1], up_nom[-1], s=60, marker="s", label="End" if idx == 1 else None)

    # Label the local coordinate axes in meters.
    ax.set_xlabel("East m")
    ax.set_ylabel("North m")
    ax.set_zlabel("Altitude m")

    # Give each panel a simple title based on the view number.
    ax.set_title(f"3D descent tube view {idx}")

    # Set the camera angle for this subplot.
    ax.view_init(elev=elev, azim=azim)

    # Show the legend only once so the figure does not repeat it in every panel.
    if idx == 1:
        ax.legend()

# Tighten layout spacing so the subplots fit cleanly in the figure.
plt.tight_layout()

# Render the completed 3D figure.
plt.show()

# Diagnose whether guidance is permanently disabled by the velocity gate

print("Minimum logged speed m/s:", df["V_mps"].min())
print("Maximum logged speed m/s:", df["V_mps"].max())
print("Any rows below 6000 m/s:", bool((df["V_mps"] <= 6000.0).any()))
print("Unique fired_this_step values:", sorted(df["fired_this_step"].unique()))
print("Max abs tau_roll_cmd_Nm:", float(np.abs(df["tau_roll_cmd_Nm"]).max()))
print("Max abs sigma_cmd_rad:", float(np.abs(df["sigma_cmd_rad"]).max()))
print("Max abs sigma_target_rad:", float(np.abs(df["sigma_target_rad"]).max()))

# Check whether the guidance command is identically zero

df[
    [
        "t_s",
        "V_mps",
        "sigma_cmd_rad",
        "sigma_target_rad",
        "sigma_actual_rad",
        "tau_roll_cmd_Nm",
        "requested_duty",
        "fired_this_step",
    ]
].head(20)

# Plot flight path angle and heading interval envelopes

plt.figure(figsize=(10, 5))
plt.plot(df["t_s"], np.degrees(df["gamma_rad"]), label="Nominal gamma")
plt.plot(df["t_s"], np.degrees(df["gamma_rad_lo"]), label="Interval gamma low")
plt.plot(df["t_s"], np.degrees(df["gamma_rad_hi"]), label="Interval gamma high")
plt.xlabel("Time s")
plt.ylabel("Gamma deg")
plt.title("Conditional uncertainty propagation for flight path angle")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(df["t_s"], np.degrees(df["chi_rad"]), label="Nominal chi")
plt.plot(df["t_s"], np.degrees(df["chi_rad_lo"]), label="Interval chi low")
plt.plot(df["t_s"], np.degrees(df["chi_rad_hi"]), label="Interval chi high")
plt.xlabel("Time s")
plt.ylabel("Chi deg")
plt.title("Conditional uncertainty propagation for heading angle")
plt.legend()
plt.grid(True)
plt.show()

# Plot density and dynamic pressure interval envelopes

plt.figure(figsize=(10, 5))
plt.plot(df["t_s"], df["rho_kgm3"], label="Nominal density")
plt.plot(df["t_s"], df["interval_rho_lo"], label="Interval density low")
plt.plot(df["t_s"], df["interval_rho_hi"], label="Interval density high")
plt.xlabel("Time s")
plt.ylabel("Density kg per m^3")
plt.title("Conditional uncertainty propagation for density")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(df["t_s"], df["q_pa"], label="Nominal dynamic pressure")
plt.plot(df["t_s"], df["interval_q_lo"], label="Interval q low")
plt.plot(df["t_s"], df["interval_q_hi"], label="Interval q high")
plt.xlabel("Time s")
plt.ylabel("Dynamic pressure Pa")
plt.title("Conditional uncertainty propagation for dynamic pressure")
plt.legend()
plt.grid(True)
plt.show()

# Plot interval and nominal heating with the interval failure marker

plt.figure(figsize=(10, 5))
plt.plot(df["t_s"], df["nominal_heat_qdot_max_hi"], label="Nominal heating qdot high")
plt.plot(df["t_s"], df["heating_qdot_max_hi"], label="Interval heating qdot high")

if len(df) > 0 and df["interval_failed"].max() > 0:
    t_fail = float(df[df["interval_failed"] == 1]["interval_failure_time_s"].iloc[0])
    plt.axvline(t_fail, linestyle="--", label="Interval failure time")

plt.xlabel("Time s")
plt.ylabel("Heating proxy")
plt.title("Nominal and interval heating rate")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(df["t_s"], df["nominal_heat_Q_max_hi"], label="Nominal integrated heating high")
plt.plot(df["t_s"], df["heating_Q_max_hi"], label="Interval integrated heating high")

if len(df) > 0 and df["interval_failed"].max() > 0:
    t_fail = float(df[df["interval_failed"] == 1]["interval_failure_time_s"].iloc[0])
    plt.axvline(t_fail, linestyle="--", label="Interval failure time")

plt.xlabel("Time s")
plt.ylabel("Integrated heating proxy")
plt.title("Nominal and interval integrated heating")
plt.legend()
plt.grid(True)
plt.show()

# Heat shield face maps for total integrated heat load
#
# This cell replaces the old ring only plotting logic with a ring and sector plot.
# It assumes the HeatShield class has already been upgraded so that the shield is
# split into both radial rings and angular sectors.
#
# What this cell shows
# The first panel shows the nominal final integrated heat load by shield cell.
# The second panel shows the interval lower bound final integrated heat load by shield cell.
# The third panel shows the interval upper bound final integrated heat load by shield cell.
#
# Why this cell is different from the earlier version
# The earlier version drew each ring as a full circular annulus.
# That was correct only for the old radial only shield model.
# This new version draws each shield cell as an annular sector.
# That allows the face map to show both radial and angular heating structure.
#
# Important modeling note
# This plotting cell can only display directional heating that already exists in the
# underlying heat shield data. If the simulation updated the heat shield with a fixed
# hot direction for every step, then the sectors will still look symmetric over time.
# For a truly moving hot region, the main heating update path must pass a changing
# hot_theta_rad during the simulation.
#
# Design choices in this cell
# One shared color scale is used across all three panels.
# That makes the colors directly comparable between nominal lower bound and upper bound.
#
# The colorbar is placed in its own dedicated axis on the far right.
# That prevents any overlap with the third panel.
#
# Labels are only drawn when the number of shield cells is modest.
# For a coarse grid such as 2 rings by 4 sectors, labels remain readable.
# For finer grids, the plot stays clean and the numeric values are still available
# in the summary DataFrame shown below the figure.

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from interval_math import promote

# Validate that the required simulation objects exist

# The DataFrame is needed so the nominal shield can be reconstructed by replaying
# the nominal trajectory heating one time step at a time.
if "df" not in globals() or len(df) == 0:
    raise ValueError("df is missing or empty. Run the main simulation loop first.")

# The interval heat shield must already exist because it carries the interval based
# accumulated heat loads that are visualized in the second and third panels.
if "interval_heat_shield" not in globals() or interval_heat_shield is None:
    raise ValueError("interval_heat_shield is missing. Run the main simulation loop first with heating enabled.")

# This revised plotting cell expects the upgraded HeatShield class that stores
# ring sector geometry. If these attributes are missing, the notebook is still
# using the old ring only HeatShield implementation.
required_sector_attrs = [
    "num_sectors",
    "theta_edges",
    "theta_centers",
    "cell_ring_index",
    "cell_sector_index",
    "cell_r_inner",
    "cell_r_outer",
    "cell_r_center",
    "cell_theta_lo",
    "cell_theta_hi",
    "cell_theta_center",
    "cell_area",
]
missing_sector_attrs = [name for name in required_sector_attrs if not hasattr(interval_heat_shield, name)]
if missing_sector_attrs:
    raise ValueError(
        "interval_heat_shield is missing sector geometry attributes. "
        "Make sure the upgraded HeatShield class has been loaded before running this cell. "
        f"Missing attributes: {missing_sector_attrs}"
    )

# Small helpers for direction replay and label formatting

def get_hot_theta_from_row(row):
    """
    Return the heating direction angle to use for nominal replay.

    Preferred behavior
    If the main loop logged a dedicated heating direction column, use it.

    Fallback behavior
    If no such column exists, fall back to 0.0 so the nominal replay still runs.
    That fallback keeps the code stable, but it does not create a moving hot region.
    """
    candidate_columns = [
        "hot_theta_rad",
        "heating_hot_theta_rad",
        "shield_hot_theta_rad",
    ]

    for col in candidate_columns:
        if col in row.index and pd.notna(row[col]):
            return float(row[col]), col

    return 0.0, "fallback_constant_zero"

def choose_text_color(value, norm):
    """
    Choose a text color that stays readable against the cell color.
    """
    return "white" if norm(value) < 0.58 else "black"

def should_draw_labels(n_cells_total):
    """
    Draw labels only when the shield resolution is still visually readable.
    """
    return n_cells_total <= 16

# Build a fresh nominal heat shield and replay the nominal heating history

# Pull the discretization from constants when available.
# If the constants file is not yet updated, fall back to the geometry already
# present on the interval shield so that nominal and interval plots still align.
num_rings_to_use = int(getattr(constants, "HEAT_SHIELD_NUM_RINGS", interval_heat_shield.num_rings))
num_sectors_to_use = int(getattr(constants, "HEAT_SHIELD_NUM_SECTORS", interval_heat_shield.num_sectors))
radial_exp_to_use = float(getattr(constants, "HEAT_SHIELD_RADIAL_EXP", 1.0))
azimuthal_gain_to_use = float(getattr(constants, "HEAT_SHIELD_AZIMUTHAL_GAIN", 0.0))

# Rebuild a nominal shield using the same geometry settings that the upgraded class expects.
# This nominal shield is separate from the interval shield because the nominal replay uses
# point values rather than interval bounds.
nominal_heat_shield = constants.HeatShield(
    radius_m=float(constants.HEAT_SHIELD_RADIUS_M),
    nose_radius_m=float(constants.HEAT_SHIELD_NOSE_RADIUS_M),
    num_rings=num_rings_to_use,
    num_sectors=num_sectors_to_use,
    radial_exp=radial_exp_to_use,
    azimuthal_gain=azimuthal_gain_to_use,
)

# Confirm that the nominal replay shield and the interval shield describe the same
# number of cells. The plotting logic assumes a one to one correspondence.
if len(nominal_heat_shield.Q) != len(interval_heat_shield.Q):
    raise ValueError(
        "Nominal shield cell count does not match interval shield cell count. "
        "Make sure the notebook is using the same ring and sector settings for both."
    )

# Replay the nominal heat history step by step.
# Each time step uses the logged nominal density and speed.
# If a heating direction column exists, it is also replayed here.
direction_source_name = None

for _, row in df.iterrows():
    rho_nom = promote(float(row["rho_kgm3"]))
    V_nom = promote(float(row["V_mps"]))
    hot_theta_rad, direction_source_name = get_hot_theta_from_row(row)

    nominal_heat_shield.update(
        rho=rho_nom,
        V=V_nom,
        dt=float(dt_s),
        hot_theta_rad=float(hot_theta_rad),
    )

# Extract final cellwise integrated heat values for all three cases

# The nominal shield stores punctual intervals, so mid returns the scalar value cleanly.
# The interval shield stores lower and upper bounds for each cell.
nominal_Q = np.array([float(q.mid()) for q in nominal_heat_shield.Q], dtype=float)
lower_Q = np.array([float(q.lo) for q in interval_heat_shield.Q], dtype=float)
upper_Q = np.array([float(q.hi) for q in interval_heat_shield.Q], dtype=float)

# Extract shield cell geometry from the interval shield

# These arrays describe the exact annular sector geometry for every shield cell.
# They are used to draw the face maps.
cell_ring_index = np.array(interval_heat_shield.cell_ring_index, dtype=int)
cell_sector_index = np.array(interval_heat_shield.cell_sector_index, dtype=int)
cell_r_inner = np.array(interval_heat_shield.cell_r_inner, dtype=float)
cell_r_outer = np.array(interval_heat_shield.cell_r_outer, dtype=float)
cell_r_center = np.array(interval_heat_shield.cell_r_center, dtype=float)
cell_theta_lo = np.array(interval_heat_shield.cell_theta_lo, dtype=float)
cell_theta_hi = np.array(interval_heat_shield.cell_theta_hi, dtype=float)
cell_theta_center = np.array(interval_heat_shield.cell_theta_center, dtype=float)
cell_area = np.array(interval_heat_shield.cell_area, dtype=float)

shield_radius = float(interval_heat_shield.radius)
n_cells_total = len(cell_area)

# Build one shared color scale across all three panels

# A single normalization makes color comparisons fair across all three scenarios.
all_vals = np.concatenate([nominal_Q, lower_Q, upper_Q])
vmin = float(np.min(all_vals))
vmax = float(np.max(all_vals))

# Avoid a zero width normalization range in the degenerate case where all values match.
if np.isclose(vmin, vmax):
    vmax = vmin + 1.0

norm = Normalize(vmin=vmin, vmax=vmax)

# Colormap design
# Lower values move through blue and purple
# Higher values move toward bright yellow
cmap = LinearSegmentedColormap.from_list(
    "blue_purple_yellow",
    [
        "#163b8c",
        "#3b4cc0",
        "#5e2b97",
        "#8e44ad",
        "#c77dff",
        "#f4d35e",
        "#ffe45e",
    ],
    N=256,
)

# Helper function to draw one full shield face map using annular sectors

def draw_sector_heatmap(ax, values, title):
    """
    Draw one shield face map.

    Each shield cell is drawn as a wedge with:
    inner radius from the ring geometry
    outer radius from the ring geometry
    starting angle from the sector geometry
    ending angle from the sector geometry

    This is the key difference from the old ring only plot.
    Instead of drawing full circles, this draws individual ring sector cells.
    """
    ax.set_aspect("equal")
    ax.set_xlim(-shield_radius * 1.10, shield_radius * 1.10)
    ax.set_ylim(-shield_radius * 1.10, shield_radius * 1.10)

    # Draw every ring sector cell.
    for cell_idx, value in enumerate(values):
        r_inner = cell_r_inner[cell_idx]
        r_outer = cell_r_outer[cell_idx]
        theta_lo_deg = math.degrees(cell_theta_lo[cell_idx])
        theta_hi_deg = math.degrees(cell_theta_hi[cell_idx])

        patch = Wedge(
            center=(0.0, 0.0),
            r=r_outer,
            theta1=theta_lo_deg,
            theta2=theta_hi_deg,
            width=(r_outer - r_inner),
            facecolor=cmap(norm(value)),
            edgecolor="black",
            linewidth=1.1,
        )
        ax.add_patch(patch)

    # Add the outer shield boundary so the full face remains visually crisp
    outer_boundary = plt.Circle(
        (0.0, 0.0),
        shield_radius,
        fill=False,
        color="black",
        linewidth=1.5,
    )
    ax.add_artist(outer_boundary)
    
    # Draw labels only for modest cell counts.
    # Fine grids become unreadable if every cell is labeled on the plot itself
    if should_draw_labels(n_cells_total):
        for cell_idx, value in enumerate(values):
            theta = cell_theta_center[cell_idx]
            r_label = 0.72 * cell_r_center[cell_idx]
            x_label = r_label * math.cos(theta)
            y_label = r_label * math.sin(theta)

            ring_id = cell_ring_index[cell_idx] + 1
            sector_id = cell_sector_index[cell_idx] + 1

            ax.text(
                x_label,
                y_label,
                f"R{ring_id} S{sector_id}\n{value:.2e}",
                ha="center",
                va="center",
                fontsize=9,
                color=choose_text_color(value, norm),
                fontweight="bold",
            )

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=16)

# Create the figure and draw the three panels

fig, axes = plt.subplots(1, 3, figsize=(20, 7))

# Leave space on the right for a dedicated colorbar axis
fig.subplots_adjust(left=0.04, right=0.90, wspace=0.12)

draw_sector_heatmap(axes[0], nominal_Q, "Nominal total heat load")
draw_sector_heatmap(axes[1], lower_Q, "Interval lower bound total heat load")
draw_sector_heatmap(axes[2], upper_Q, "Interval upper bound total heat load")

# Build a separate colorbar axis so it never overlaps the third panel
cbar_ax = fig.add_axes([0.92, 0.18, 0.015, 0.64])

sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label("Integrated heat load Q", fontsize=12)

fig.suptitle("Heat shield face maps by ring and sector", fontsize=18, y=0.96)

plt.show()

# Build a summary DataFrame with one row per shield cell

# This table is useful when the grid becomes fine enough that text labels are
# intentionally suppressed on the figure
heat_cells_df = pd.DataFrame({
    "ring_index": cell_ring_index,
    "sector_index": cell_sector_index,
    "ring_label": [f"R{i + 1}" for i in cell_ring_index],
    "sector_label": [f"S{i + 1}" for i in cell_sector_index],
    "r_inner_m": cell_r_inner,
    "r_outer_m": cell_r_outer,
    "r_center_m": cell_r_center,
    "theta_lo_deg": np.degrees(cell_theta_lo),
    "theta_hi_deg": np.degrees(cell_theta_hi),
    "theta_center_deg": np.degrees(cell_theta_center),
    "cell_area_m2": cell_area,
    "Q_nominal": nominal_Q,
    "Q_interval_lo": lower_Q,
    "Q_interval_hi": upper_Q,
})

# Print the direction source used for the nominal replay.
# This makes it easy to see whether the notebook is replaying a logged heating
# direction or just using the fixed fallback angle.
print(f"Nominal replay hot direction source: {direction_source_name}")

heat_cells_df

# plot state widths and mark health intervention steps

active_df = df[df["interval_active"] == 1].copy()

if len(active_df) == 0:
    print("No active interval rows were logged")
else:
    plt.figure(figsize=(10, 5))
    plt.plot(active_df["t_s"], active_df["width_r_m"], label="r width")
    plt.plot(active_df["t_s"], active_df["width_V_mps"], label="V width")
    plt.plot(active_df["t_s"], active_df["width_gamma_rad"], label="gamma width")
    plt.plot(active_df["t_s"], active_df["width_chi_rad"], label="chi width")

    recenter_width_df = active_df[
        (active_df["interval_step_recentered"] == 1)
        & (active_df["interval_step_recenter_reason"] == "width_threshold")
    ]
    if len(recenter_width_df) > 0:
        plt.scatter(
            recenter_width_df["t_s"],
            recenter_width_df["width_r_m"],
            marker="o",
            s=30,
            label="recenter width",
        )

    split_width_df = active_df[
        (active_df["interval_step_split"] == 1)
        & (active_df["interval_step_split_reason"] == "width_threshold")
    ]
    if len(split_width_df) > 0:
        plt.scatter(
            split_width_df["t_s"],
            split_width_df["width_r_m"],
            marker="x",
            s=35,
            label="split width",
        )

    split_denom_df = active_df[
        (active_df["interval_step_split"] == 1)
        & (active_df["interval_step_split_reason"] != "width_threshold")
        & (active_df["interval_step_split_reason"] != "")
    ]
    if len(split_denom_df) > 0:
        plt.scatter(
            split_denom_df["t_s"],
            split_denom_df["width_r_m"],
            marker="^",
            s=35,
            label="split denominator",
        )

    if df["interval_failed"].max() > 0:
        t_fail = float(df[df["interval_failed"] == 1]["interval_failure_time_s"].iloc[0])
        plt.axvline(t_fail, linestyle="--", label="Interval failure time")

    plt.xlabel("Time s")
    plt.ylabel("Interval width")
    plt.title("Growth of state uncertainty with health events")
    plt.legend()
    plt.grid(True)
    plt.show()

# plot derivative widths and mark health intervention steps

active_df = df[df["interval_active"] == 1].copy()

if len(active_df) == 0:
    print("No active interval rows were logged")
else:
    plt.figure(figsize=(10, 5))
    plt.plot(active_df["t_s"], active_df["dx_width_r"], label="r dot width")
    plt.plot(active_df["t_s"], active_df["dx_width_V"], label="V dot width")
    plt.plot(active_df["t_s"], active_df["dx_width_gamma"], label="gamma dot width")
    plt.plot(active_df["t_s"], active_df["dx_width_chi"], label="chi dot width")

    recenter_width_df = active_df[
        (active_df["interval_step_recentered"] == 1)
        & (active_df["interval_step_recenter_reason"] == "width_threshold")
    ]
    if len(recenter_width_df) > 0:
        plt.scatter(
            recenter_width_df["t_s"],
            recenter_width_df["dx_width_r"],
            marker="o",
            s=30,
            label="recenter width",
        )

    split_width_df = active_df[
        (active_df["interval_step_split"] == 1)
        & (active_df["interval_step_split_reason"] == "width_threshold")
    ]
    if len(split_width_df) > 0:
        plt.scatter(
            split_width_df["t_s"],
            split_width_df["dx_width_r"],
            marker="x",
            s=35,
            label="split width",
        )

    split_denom_df = active_df[
        (active_df["interval_step_split"] == 1)
        & (active_df["interval_step_split_reason"] != "width_threshold")
        & (active_df["interval_step_split_reason"] != "")
    ]
    if len(split_denom_df) > 0:
        plt.scatter(
            split_denom_df["t_s"],
            split_denom_df["dx_width_r"],
            marker="^",
            s=35,
            label="split denominator",
        )

    if df["interval_failed"].max() > 0:
        t_fail = float(df[df["interval_failed"] == 1]["interval_failure_time_s"].iloc[0])
        plt.axvline(t_fail, linestyle="--", label="Interval failure time")

    plt.xlabel("Time s")
    plt.ylabel("Derivative interval width")
    plt.title("Growth of derivative uncertainty with health events")
    plt.legend()
    plt.grid(True)
    plt.show()

# Quick interval status summary

print(df["safety_status"].value_counts(dropna=False))
print()
print(df["interval_failure_kind"].value_counts(dropna=False))

# Show rows where the safety status is not fully inside

df[df["safety_status"] != "inside"][
    [
        "t_s",
        "alt_m",
        "V_mps",
        "q_pa",
        "interval_alt_lo",
        "interval_alt_hi",
        "interval_q_lo",
        "interval_q_hi",
        "safety_status",
        "safety_checks",
    ]
].head(20)

# Save outputs for later review

# Save the main success trajectory log
df.to_csv(success_trajectory_csv_path, index=False)
print("Saved", success_trajectory_csv_path)

# Save the failed heat step log for RL teaching
failed_heat_steps_df.to_csv(failed_heat_step_csv_path, index=False)
print("Saved", failed_heat_step_csv_path)

# Save the compact failed episode summary
failed_episode_summary_df.to_csv(failed_heat_episode_csv_path, index=False)
print("Saved", failed_heat_episode_csv_path)

# Also save a legacy named copy of the main log for compatibility
output_csv = str(PY_OUTPUT_DIR / "conditional_uncertainty_propagation_log.csv")
df.to_csv(output_csv, index=False)
print("Saved", output_csv)

# --- CPAS outputs ---
if cpas_events_log:
    cpas_events_df = pd.DataFrame(cpas_events_log)
else:
    cpas_events_df = pd.DataFrame(columns=["t_s", "step", "alt_m", "V_mps", "mach", "phase", "event"])
cpas_events_df.to_csv(cpas_events_csv_path, index=False)
print("Saved", cpas_events_csv_path)

cpas_summary = cpas_inst.summary()

# Run summary JSON (mirrors the notebook's revision_v1/run_summary.json)
run_summary = {
    "trajectory_id": trajectory_id,
    "num_logged_steps": int(len(df)),
    "termination_reason": termination_reason,
    "final_t_s": float(df["t_s"].iloc[-1]) if len(df) else 0.0,
    "final_alt_m": float(df["alt_m"].iloc[-1]) if len(df) else 0.0,
    "final_V_mps": float(df["V_mps"].iloc[-1]) if len(df) else 0.0,
    "final_gamma_deg": float(math.degrees(df["gamma_rad"].iloc[-1])) if len(df) else 0.0,
    "aero_model": params["aero_model"],
    "cpas_summary": cpas_summary,
    "cpas_drogue_deploy_t_s": cpas_summary["drogue_deploy_t_s"],
    "cpas_pilot_deploy_t_s": cpas_summary["pilot_deploy_t_s"],
    "cpas_main_deploy_t_s": cpas_summary["main_deploy_t_s"],
    "cpas_landed_t_s": cpas_summary["landed_t_s"],
    "output_files": {
        "trajectory_csv": success_trajectory_csv_path,
        "failed_heat_step_csv": failed_heat_step_csv_path,
        "failed_heat_episode_csv": failed_heat_episode_csv_path,
        "cpas_events_csv": cpas_events_csv_path,
    },
}
with open(run_summary_json_path, "w", encoding="utf-8") as f:
    json.dump(run_summary, f, indent=2, default=str)
print("Saved", run_summary_json_path)

# Landing summary JSON (final touchdown metrics)
landing_summary = {
    "termination_reason": termination_reason,
    "final_t_s": float(df["t_s"].iloc[-1]) if len(df) else 0.0,
    "final_alt_km": float(df["alt_m"].iloc[-1]) / 1000.0 if len(df) else 0.0,
    "final_V_mps": float(df["V_mps"].iloc[-1]) if len(df) else 0.0,
    "final_gamma_deg": float(math.degrees(df["gamma_rad"].iloc[-1])) if len(df) else 0.0,
    "final_chi_deg": float(math.degrees(df["chi_rad"].iloc[-1])) if len(df) else 0.0,
    "final_phi_deg": float(math.degrees(df["phi_rad"].iloc[-1])) if len(df) else 0.0,
    "final_lam_deg": float(math.degrees(df["lam_rad"].iloc[-1])) if len(df) else 0.0,
    "cpas_drogue_deploy_t_s": cpas_summary["drogue_deploy_t_s"],
    "cpas_pilot_deploy_t_s": cpas_summary["pilot_deploy_t_s"],
    "cpas_main_deploy_t_s": cpas_summary["main_deploy_t_s"],
    "cpas_landed_t_s": cpas_summary["landed_t_s"],
}
with open(landing_summary_json_path, "w", encoding="utf-8") as f:
    json.dump(landing_summary, f, indent=2, default=str)
print("Saved", landing_summary_json_path)

# Check for heading singularity using the current logged gamma column

plt.figure(figsize=(10, 5))
plt.plot(df["t_s"], np.cos(df["gamma_rad"]), label="cos gamma nominal")
plt.axhline(0.0)
plt.xlabel("Time s")
plt.ylabel("cos gamma")
plt.title("Check for heading singularity")
plt.legend()
plt.grid(True)
plt.show()

# plot interval activity and health event flags on one timeline

plt.figure(figsize=(10, 5))
plt.plot(df["t_s"], df["alt_m"], label="Nominal altitude")

active_df = df[df["interval_active"] == 1].copy()
if len(active_df) > 0:
    plt.fill_between(
        active_df["t_s"],
        active_df["interval_alt_lo"],
        active_df["interval_alt_hi"],
        alpha=0.25,
        label="Interval altitude band while active",
    )

if df["interval_failed"].max() > 0:
    t_fail = float(df[df["interval_failed"] == 1]["interval_failure_time_s"].iloc[0])
    plt.axvline(t_fail, linestyle="--", label="Interval failure time")

plt.xlabel("Time s")
plt.ylabel("Altitude m")
plt.title("Nominal continuation after interval failure")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 4))
plt.step(df["t_s"], df["interval_active"], where="post", label="Interval active")
plt.step(df["t_s"], df["nominal_continued_after_interval_failure"], where="post", label="Nominal only continuation")
plt.step(df["t_s"], df["interval_step_recentered"], where="post", label="Recenter step")
plt.step(df["t_s"], df["interval_step_split"], where="post", label="Split step")
plt.xlabel("Time s")
plt.ylabel("Flag")
plt.title("Interval status and health event flags")
plt.yticks([0, 1])
plt.legend()
plt.grid(True)
plt.show()

event_cols = [
    "t_s",
    "interval_active",
    "interval_failed",
    "interval_step_recentered",
    "interval_step_recenter_reason",
    "interval_step_split",
    "interval_step_split_reason",
    "interval_step_split_depth",
    "nominal_continued_after_interval_failure",
]

df[
    (df["interval_step_recentered"] == 1)
    | (df["interval_step_split"] == 1)
    | (df["interval_failed"] == 1)
][event_cols].head(25)

# show interval widths while the interval layer is active

active_df = df[df["interval_active"] == 1].copy()

if len(active_df) == 0:
    print("No active interval rows were logged")
else:
    plt.figure(figsize=(10, 5))
    plt.plot(active_df["t_s"], active_df["width_r_m"], label="r width")
    plt.plot(active_df["t_s"], active_df["width_V_mps"], label="V width")
    plt.plot(active_df["t_s"], active_df["width_gamma_rad"], label="gamma width")
    plt.plot(active_df["t_s"], active_df["width_chi_rad"], label="chi width")
    plt.xlabel("Time s")
    plt.ylabel("Interval width")
    plt.title("Interval widths while interval is active")
    plt.legend()
    plt.grid(True)
    plt.show()

    print(
        active_df[
            [
                "t_s",
                "width_r_m",
                "width_V_mps",
                "width_gamma_rad",
                "width_chi_rad",
                "interval_step_recentered",
                "interval_step_recenter_reason",
                "interval_step_split",
                "interval_step_split_reason",
            ]
        ].tail(10)
    )

# Add degree based convenience columns for control diagnostics

df["sigma_cmd_deg"] = np.degrees(df["sigma_cmd_rad"])
df["sigma_target_deg"] = np.degrees(df["sigma_target_rad"])
df["sigma_actual_deg"] = np.degrees(df["sigma_actual_rad"])
df["roll_rate_deg_s"] = np.degrees(df["roll_rate_rad_s"])
df["roll_accel_deg_s2"] = np.degrees(df["roll_accel_rad_s2"])

df[
    [
        "t_s",
        "sigma_cmd_deg",
        "sigma_target_deg",
        "sigma_actual_deg",
        "roll_rate_deg_s",
        "tau_roll_cmd_Nm",
        "requested_duty",
        "fired_this_step",
        "active_thrusters",
        "torque_z_from_rcs",
    ]
].head(12)

# Initialize a new figure with specific dimensions for the plot
plt.figure(figsize=(11, 5))

# Plot the commanded bank angle with a dashed line and increased width to prevent it from being hidden
plt.plot(df["t_s"], df["sigma_cmd_deg"], label="Sigma cmd", linestyle="dashed", linewidth=3.0, alpha=0.7)

# Plot the smoothed target bank angle with a dotted line
plt.plot(df["t_s"], df["sigma_target_deg"], label="Sigma target", linewidth=2.5, alpha=0.7)

# Plot the actual flown bank angle with a standard solid line
plt.plot(df["t_s"], df["sigma_actual_deg"], label="Sigma actual", linewidth=1.5)

# Apply text labels to both the horizontal and vertical axes
plt.xlabel("Time s")
plt.ylabel("Bank angle deg")

# Set the main title for the graph
plt.title("Closed loop bank angle behavior")

# Insert a legend to clearly identify each variable
plt.legend()

# Enable the background grid for easier visual tracking of values
plt.grid(True)

# Render and display the completed plot
plt.show()

# Plot actual bank tracking error

tracking_error_deg = df["sigma_target_deg"] - df["sigma_actual_deg"]

plt.figure(figsize=(11, 5))
plt.plot(df["t_s"], tracking_error_deg, label="Sigma target minus sigma actual")
plt.xlabel("Time s")
plt.ylabel("Tracking error deg")
plt.title("Bank angle tracking error")
plt.legend()
plt.grid(True)
plt.show()

# Plot roll rate and roll acceleration

plt.figure(figsize=(11, 5))
plt.plot(df["t_s"], df["roll_rate_deg_s"], label="Roll rate")
plt.xlabel("Time s")
plt.ylabel("Roll rate deg per s")
plt.title("Roll rate over time")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(11, 5))
plt.plot(df["t_s"], df["roll_accel_deg_s2"], label="Roll acceleration")
plt.xlabel("Time s")
plt.ylabel("Roll acceleration deg per s squared")
plt.title("Roll acceleration from RCS torque")
plt.legend()
plt.grid(True)
plt.show()

# Plot commanded roll torque and realized RCS roll torque

plt.figure(figsize=(11, 5))
plt.plot(df["t_s"], df["tau_roll_cmd_Nm"], label="Commanded roll torque")
plt.plot(df["t_s"], df["torque_z_from_rcs"], label="Realized RCS z torque")
plt.xlabel("Time s")
plt.ylabel("Torque N m")
plt.title("Roll torque command versus realized RCS torque")
plt.legend()
plt.grid(True)
plt.show()

# Plot requested duty and whether a firing actually happened this step

plt.figure(figsize=(11, 5))
plt.plot(df["t_s"], df["requested_duty"], label="Requested duty")
plt.xlabel("Time s")
plt.ylabel("Duty")
plt.title("Requested roll channel duty")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(11, 3))
plt.step(df["t_s"], df["fired_this_step"], where="post", label="Fired this step")
plt.xlabel("Time s")
plt.ylabel("Fired")
plt.title("Binary RCS firing events")
plt.yticks([0, 1])
plt.legend()
plt.grid(True)
plt.show()

# Expand active thrusters into one row per thruster firing event

# Build a flat event table so each individual thruster firing can be analyzed on its own row
thruster_event_rows = []

# Walk through the trajectory log row by row to inspect which thrusters were active at each time step
for _, row in df.iterrows():
    # Read the raw active thruster field for this time step
    active_str = row["active_thrusters"]

    # Skip rows where no thruster activity was logged
    if pd.isna(active_str):
        continue

    # Clean up whitespace so empty strings do not get treated as valid activity
    active_str = str(active_str).strip()

    # Skip rows where the field exists but contains no actual thruster names
    if active_str == "":
        continue

    # Split the comma separated thruster list into individual thruster names
    active_list = [name.strip() for name in active_str.split(",") if name.strip()]

    # Create one output row per active thruster so each firing becomes its own event record
    for thr_name in active_list:
        thruster_event_rows.append(
            {
                # Save the simulation time of this firing event
                "t_s": float(row["t_s"]),

                # Save the specific thruster name that was active
                "thruster": thr_name,

                # Keep the requested duty value to compare command strength against actual firing behavior
                "requested_duty": float(row["requested_duty"]),

                # Keep the commanded roll torque so each firing can be tied back to the controller demand
                "tau_roll_cmd_Nm": float(row["tau_roll_cmd_Nm"]),

                # Keep the actual roll axis torque generated by the RCS system at this step
                "torque_z_from_rcs": float(row["torque_z_from_rcs"]),
            }
        )

# Convert the collected event records into a dataframe for inspection and plotting
thruster_events_df = pd.DataFrame(thruster_event_rows)

# Print the number of individual firing events extracted from the original log
print("Number of thruster firing records", len(thruster_events_df))

# Show the first several event rows to verify the expansion worked as expected
thruster_events_df.head(20)

# Scatter plot of which thrusters fired and when

if len(thruster_events_df) == 0:
    print("No thruster firing events were logged")
else:
    thruster_names = sorted(thruster_events_df["thruster"].unique())
    thruster_to_y = {name: idx for idx, name in enumerate(thruster_names)}

    y_vals = thruster_events_df["thruster"].map(thruster_to_y)

    plt.figure(figsize=(12, 6))
    plt.scatter(thruster_events_df["t_s"], y_vals)
    plt.xlabel("Time s")
    plt.ylabel("Thruster")
    plt.title("RCS firing events by thruster")
    plt.yticks(list(thruster_to_y.values()), list(thruster_to_y.keys()))
    plt.grid(True)
    plt.show()

# Estimated on time per active thruster during each step

# This is estimated on time per active thruster during the step.
# In this milestone 1 model, a firing step means the channel fired for the full dt.
# The requested duty still shows the fractional demand before pulse accumulation.

if len(thruster_events_df) == 0:
    print("No thruster firing events were logged")
else:
    thruster_events_df["estimated_on_time_s"] = dt_s

    plt.figure(figsize=(12, 6))
    plt.scatter(
        thruster_events_df["t_s"],
        thruster_events_df["estimated_on_time_s"],
    )
    plt.xlabel("Time s")
    plt.ylabel("Estimated on time s")
    plt.title("Estimated thruster on time per firing event")
    plt.grid(True)
    plt.show()

# Requested duty versus actual firing events

if len(thruster_events_df) == 0:
    print("No thruster firing events were logged")
else:
    plt.figure(figsize=(11, 5))
    plt.scatter(df["t_s"], df["requested_duty"], label="Requested duty")
    plt.scatter(
        thruster_events_df["t_s"],
        thruster_events_df["requested_duty"],
        label="Duty on firing steps",
    )
    plt.xlabel("Time s")
    plt.ylabel("Duty")
    plt.title("Requested duty and actual firing steps")
    plt.legend()
    plt.grid(True)
    plt.show()

# Rolling summary of firing frequency

window_steps = 20

df["firing_rate_window"] = df["fired_this_step"].rolling(window_steps, min_periods=1).mean()

plt.figure(figsize=(11, 5))
plt.plot(df["t_s"], df["firing_rate_window"], label="Rolling firing fraction")
plt.xlabel("Time s")
plt.ylabel("Fraction of steps firing")
plt.title("RCS firing activity over time")
plt.legend()
plt.grid(True)
plt.show()

# final compact view with interval health semantics

summary_cols = [
    "t_s",
    "alt_m",
    "V_mps",
    "gamma_rad",
    "sigma_actual_deg",
    "interval_active",
    "interval_failed",
    "interval_failure_kind",
    "interval_step_recentered",
    "interval_step_recenter_reason",
    "interval_step_split",
    "interval_step_split_reason",
    "interval_step_split_depth",
    "nominal_continued_after_interval_failure",
    "nominal_heat_qdot_max_hi",
    "nominal_heat_Q_max_hi",
    "width_r_m",
    "width_V_mps",
    "width_gamma_rad",
    "width_chi_rad",
    "interval_rho_lo",
    "interval_rho_hi",
    "interval_q_lo",
    "interval_q_hi",
    "heating_qdot_max_hi",
    "safety_status",
]

summary_df = df[summary_cols].copy()
summary_df["gamma_deg"] = np.degrees(summary_df["gamma_rad"])

summary_df[
    [
        "t_s",
        "alt_m",
        "V_mps",
        "gamma_deg",
        "sigma_actual_deg",
        "interval_active",
        "interval_failed",
        "interval_failure_kind",
        "interval_step_recentered",
        "interval_step_recenter_reason",
        "interval_step_split",
        "interval_step_split_reason",
        "interval_step_split_depth",
        "nominal_continued_after_interval_failure",
        "nominal_heat_qdot_max_hi",
        "nominal_heat_Q_max_hi",
        "width_r_m",
        "width_V_mps",
        "width_gamma_rad",
        "width_chi_rad",
        "interval_rho_lo",
        "interval_rho_hi",
        "interval_q_lo",
        "interval_q_hi",
        "heating_qdot_max_hi",
        "safety_status",
    ]
].tail(15)

# Plot RCS timing behavior from the logged dataframe

rcs_plot_df = df.copy()

# Realized fired fraction over the outer step
# This shows how many internal substeps actually fired
if "num_fired_internal_steps" in rcs_plot_df.columns and "num_internal_steps" in rcs_plot_df.columns:
    rcs_plot_df["realized_fire_fraction"] = (
        rcs_plot_df["num_fired_internal_steps"] / rcs_plot_df["num_internal_steps"].replace(0, np.nan)
    )
else:
    rcs_plot_df["realized_fire_fraction"] = np.nan

# Convert backlog from seconds to milliseconds for easier reading
if "roll_pos_backlog_s" in rcs_plot_df.columns:
    rcs_plot_df["roll_pos_backlog_ms"] = 1000.0 * rcs_plot_df["roll_pos_backlog_s"]

if "roll_neg_backlog_s" in rcs_plot_df.columns:
    rcs_plot_df["roll_neg_backlog_ms"] = 1000.0 * rcs_plot_df["roll_neg_backlog_s"]

# Requested duty versus realized fired fraction
plt.figure(figsize=(11, 4))
plt.plot(rcs_plot_df["t_s"], rcs_plot_df["requested_duty"], label="Requested duty")
plt.plot(rcs_plot_df["t_s"], rcs_plot_df["realized_fire_fraction"], label="Realized fire fraction")
plt.xlabel("Time s")
plt.ylabel("Fraction")
plt.title("Requested duty versus realized internal firing fraction")
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()

# Positive and negative channel backlog
plt.figure(figsize=(11, 4))
if "roll_pos_backlog_ms" in rcs_plot_df.columns:
    plt.plot(rcs_plot_df["t_s"], rcs_plot_df["roll_pos_backlog_ms"], label="Positive channel backlog ms")
if "roll_neg_backlog_ms" in rcs_plot_df.columns:
    plt.plot(rcs_plot_df["t_s"], rcs_plot_df["roll_neg_backlog_ms"], label="Negative channel backlog ms")
plt.xlabel("Time s")
plt.ylabel("Backlog ms")
plt.title("RCS accumulated backlog by roll channel")
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()

# Positive and negative channel on off state
plt.figure(figsize=(11, 4))
if "roll_pos_is_on" in rcs_plot_df.columns:
    plt.step(rcs_plot_df["t_s"], rcs_plot_df["roll_pos_is_on"], where="post", label="Positive channel on")
if "roll_neg_is_on" in rcs_plot_df.columns:
    plt.step(rcs_plot_df["t_s"], rcs_plot_df["roll_neg_is_on"], where="post", label="Negative channel on")
plt.xlabel("Time s")
plt.ylabel("Channel state")
plt.title("RCS channel on off state")
plt.yticks([0, 1])
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()

# Number of fired internal substeps inside each outer step
plt.figure(figsize=(11, 4))
if "num_fired_internal_steps" in rcs_plot_df.columns:
    plt.step(
        rcs_plot_df["t_s"],
        rcs_plot_df["num_fired_internal_steps"],
        where="post",
        label="Fired internal substeps",
    )
plt.xlabel("Time s")
plt.ylabel("Count")
plt.title("Number of internal RCS firings per outer step")
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()

# Torque command and realized torque with timing context
plt.figure(figsize=(11, 4))
plt.plot(rcs_plot_df["t_s"], rcs_plot_df["tau_roll_cmd_Nm"], label="Commanded roll torque")
plt.plot(rcs_plot_df["t_s"], rcs_plot_df["torque_z_from_rcs"], label="Realized RCS z torque")
plt.xlabel("Time s")
plt.ylabel("Torque N m")
plt.title("Commanded torque versus realized torque after timing logic")
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()

