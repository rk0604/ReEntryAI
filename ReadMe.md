ReEntryAI Interval ReEntry Notebook Usage Guide

Project overview

This codebase models a milestone 1 atmospheric reentry simulation with three tightly connected layers.

1. A nominal translational path that advances the main vehicle state.
2. A physical bank control chain that turns guidance decisions into bank motion through actuator shaping, roll torque control, and RCS firing.
3. An interval annotation layer that propagates uncertainty, tracks conservative heating bounds, and helps the predictor corrector guidance reject bad candidates.

The main notebook is interval_ReEntry3.ipynb. It ties the modules together, runs the simulation loop, logs every step into tables, saves csv outputs, and creates the diagnostic plots.

Codebase map

1. AtmosphereModel.py
   1. Implements the 1976 US Standard Atmosphere model for the float path up to 86 km.
   2. Provides the interval atmosphere functions used by the uncertainty layer when altitude bounds cross one or more atmospheric layers.
   3. Converts geometric altitude into geopotential altitude so the atmosphere calculations stay consistent with the layer definitions.
   4. Returns temperature, pressure, density, and layer index information that is used by both nominal and interval dynamics.
   5. This file is the main place to edit atmosphere behavior, density bounds, or layer handling.

2. constants.py
   1. Stores the shared simulation constants so the notebook and modules all use the same Earth model, timing values, vehicle properties, control limits, and heat limits.
   2. Defines the Earth radius, gravity constants, atmosphere layer tables, capsule mass, aerodynamic reference area, RCS geometry, roll controller gains, and predictor corrector tuning.
   3. Contains the HeatShield class, which tracks interval heat rate and accumulated heat load over shield cells.
   4. Defines the default output csv names used by the notebook save cells.
   5. This file is the main place to change global tuning values without rewriting module logic.

3. interval_math.py
   1. Provides the Interval class and the basic interval arithmetic used across the uncertainty propagation path.
   2. Includes interval versions of log, sqrt, exp, sin, cos, integer powers, and safe division logic.
   3. Supplies box helpers for interval state vectors, including box addition, scalar multiplication, midpoints, widths, and interval Euler stepping.
   4. Makes it possible for the rest of the code to reuse the same formulas with conservative bounds rather than only point values.
   5. This file is the first place to inspect when an interval propagation failure comes from division, trigonometric overestimation, or a state width blow up.

4. math_3d.py
   1. Implements the nominal and interval translational reentry dynamics for the state vector [r, phi, lam, V, gamma, chi].
   2. Contains the corrected spherical Earth equations of motion used by the notebook.
   3. Builds interval annotations for each nominal step, including state width summaries, density bounds, dynamic pressure bounds, altitude bounds, and heating summaries.
   4. Provides the short horizon rollout utility used by the predictor corrector guidance to test candidate bank commands.
   5. This file is the main physics bridge between the notebook, the guidance layer, and the heat supervisor layer.

5. control.py
   1. Holds the high level control stack for the capsule.
   2. Defines the control facing state object, the observation provider, the guidance scheduler, the heat aware predictor corrector guidance law, the bank actuator, the roll controller, and the full capsule control stack.
   3. Converts raw vehicle state into guidance features such as range to go, heading error, target azimuth, and cross track error.
   4. Evaluates multiple bank magnitude candidates, rolls each one forward through the corrected dynamics and interval heating model, scores them, and selects sigma_cmd.
   5. This file is the main place to modify target logic, candidate search behavior, bank schedules, or control shaping.

6. ReactionControl.py
   1. Implements the body frame RCS model for milestone 1.
   2. Defines fixed thrusters on the capsule, computes body force and body torque for each thruster, and allocates roll channel firings from commanded roll torque.
   3. Uses a fixed step pulse model so each selected jet is either on for the whole simulation step or off for the whole step.
   4. Supports duty accumulation across steps so small continuous roll demands can still produce realistic discrete firings.
   5. This file is the main place to adjust thruster layout, pulse logic, thrust capacity, or future pitch and yaw allocation.

7. Pasted code.py
   1. Contains the milestone1_nominal.py style module that owns the nominal milestone 1 simulation path in a standalone module form.
   2. Includes float translational dynamics, quaternion and frame helpers, attitude state objects, aerodynamic force evaluation, and the closed loop nominal step.
   3. Organizes the nominal path cleanly so math_3d.py can stay focused on interval propagation and predictor support.
   4. Serves as a reusable module level version of logic that the notebook currently wires together cell by cell.
   5. A cleaner future project layout would rename this file to milestone1_nominal.py and import it directly.

8. interval_ReEntry3.ipynb
   1. This is the main driver notebook for the full simulation and analysis workflow.
   2. Imports the shared modules, sets up the vehicle and controller, initializes the nominal and interval states, runs the closed loop physical simulation, logs every step, exports csv files, and creates plots.
   3. Also includes RCS diagnostics, interval width plots, heating plots, heat shield maps, landing error summaries, and final compact tables.
   4. This is the first file to run when checking whether the full system still works after code changes.
   5. Most debugging and visual analysis currently happens here rather than in a separate application entry point.

How the main notebook is organized

1. Import and helper setup
   The early cells import numpy, pandas, matplotlib, the project modules, and helper plotting functions.

2. Vehicle and controller setup
   The setup cells build the params dictionary, the guidance stack, the bank actuator, the roll controller, and the physical RCS system.

3. Initial conditions
   The notebook defines the initial translational state, the initial bank and roll rate, the target latitude and longitude, and the interval supervisor widths.

4. Closed loop simulation loop
   The main loop synchronizes the control state, asks guidance for sigma_cmd, runs the RCS allocator, updates actual bank motion, advances the nominal state, annotates the same step with interval propagation, and writes one log row into the dataframe.

5. Analysis and export
   The later cells inspect the dataframe, build trajectory plots, show uncertainty growth, save csv files, and analyze the RCS firing history.

How to run interval_ReEntry3.ipynb

1. Place the following files in the same working folder.
   interval_ReEntry3.ipynb
   AtmosphereModel.py
   constants.py
   control.py
   interval_math.py
   math_3d.py
   ReactionControl.py
   Pasted code.py if the milestone1_nominal module version is being kept alongside the notebook

2. Open the folder in Jupyter Notebook, JupyterLab, or Visual Studio Code with a Python notebook kernel.

3. Make sure the Python environment has numpy, pandas, matplotlib, and jupyter available.

4. Start the notebook kernel from the same folder that contains the project files. The notebook imports the modules by filename, so the working directory matters.

5. Run the notebook from top to bottom in order.

6. Watch the early setup cells for any import failure. If one module fails to import, the later control and simulation cells will not run correctly.

7. After the setup cells, run the main simulation loop cell. This is the core cell that creates df, failed_heat_steps_df, and failed_episode_summary_df.

8. Run the later plotting and export cells only after df exists and contains rows.

Important notebook run order

1. Run the import cell first.
2. Run the helper cell that defines plotting and conversion utilities.
3. Run the setup cell that builds params, guidance, actuator, controller, and RCS objects.
4. Run the initial condition cell that builds x_nominal and state_ctrl.
5. Run the supervisor config cell that defines the uncertainty widths.
6. Run the wrapper cell that defines nominal_step_closed_loop.
7. Run the heat shield activation and short diagnostic cells if needed.
8. Run the main conditional uncertainty propagation loop.
9. Run the dataframe summary and plotting cells.
10. Run the csv save cell near the end if exported logs are needed.

What the main simulation loop does

1. Copies the latest nominal state into the control facing ReentryState object.
2. Applies pre step guards to avoid singular regions such as cos gamma near zero, speed near zero, or invalid interval denominators.
3. Computes nominal aerodynamic telemetry from the current flown bank.
4. Passes live interval and heat shield context into the guidance layer.
5. Lets guidance evaluate short horizon bank candidates.
6. Converts the chosen sigma_cmd into sigma_target through the bank actuator.
7. Converts sigma tracking error into roll torque through the roll controller.
8. Allocates real RCS firings for that torque demand.
9. Updates actual bank angle and roll rate from the realized RCS torque.
10. Advances the nominal translational state using the actual flown bank angle.
11. Builds interval annotations for the same step, including density bounds, dynamic pressure bounds, altitude bounds, state widths, and heating bounds.
12. Logs all nominal, control, RCS, interval, and guidance debug values into one dataframe row.

Expected notebook outputs

1. df
   This is the main step by step trajectory log.

2. failed_heat_steps_df
   This stores guidance cycles where the selected candidate was heat infeasible or interval invalid. It is useful for later RL style labeling.

3. failed_episode_summary_df
   This gives a compact episode level summary for the current run.

4. Saved csv files
   trajectory_success.csv
   trajectory_failed_heat_steps.csv
   trajectory_failed_heat_episodes.csv
   conditional_uncertainty_propagation_log.csv

5. Diagnostic plots
   Ground track plots.
   Three dimensional trajectory views.
   Density and dynamic pressure envelopes.
   Heating envelopes.
   Heat shield face maps.
   State width growth plots.
   RCS requested duty and firing plots.
   Thruster event plots.

Useful columns in df

1. Control path columns
   sigma_cmd_rad
   sigma_target_rad
   sigma_actual_rad
   roll_rate_rad_s
   tau_roll_cmd_Nm
   requested_duty
   fired_this_step
   active_thrusters
   torque_z_from_rcs

2. Nominal flight columns
   r_m
   phi_rad
   lam_rad
   V_mps
   gamma_rad
   chi_rad
   alt_m
   rho_kgm3
   q_pa

3. Interval supervisor columns
   interval_rho_lo
   interval_rho_hi
   interval_q_lo
   interval_q_hi
   interval_alt_lo
   interval_alt_hi
   width_r_m
   width_V_mps
   width_gamma_rad
   width_chi_rad
   safety_status

4. Heating columns
   heating_qdot_max_lo
   heating_qdot_max_hi
   heating_qdot_mean_lo
   heating_qdot_mean_hi
   heating_Q_max_lo
   heating_Q_max_hi

5. Guidance predictor columns
   guidance_selected_candidate_index
   guidance_any_feasible_candidate
   guidance_selected_candidate_heat_feasible
   guidance_chosen_sigma_cmd_deg
   guidance_selected_total_cost
   guidance_selected_failure_reason
   guidance_selected_violation_amount
   guidance_candidate_sigma_deg_json
   guidance_candidate_cost_json
   guidance_candidate_heat_flag_json
   guidance_candidate_failure_reason_json

Common stopping conditions in the notebook

1. The notebook may stop the main loop when nominal cos gamma becomes too small.
   This is a guard against heading propagation singular behavior.

2. The notebook may stop when speed becomes too low for the simplified regime.
   This is an intentional modeling boundary rather than a crash.

3. The notebook may stop when altitude reaches the low altitude regime.
   This is a practical stop point for the current milestone 1 scope.

4. The notebook may stop when the interval state becomes clearly invalid.
   Examples include radius bounds becoming nonphysical or interval denominators crossing zero.

5. These stop conditions usually indicate the current model boundary or a numeric guard, not a notebook bug by themselves.

Practical debugging order

1. If the notebook fails before the main loop, check imports and module names first.

2. If the notebook fails inside interval propagation, inspect interval_math.py and math_3d.py together. Most interval failures come from denominator crossings, interval overestimation, or invalid atmosphere inputs.

3. If guidance seems inactive, inspect control.py and then the notebook guidance diagnostic cells.

4. If sigma_cmd changes but sigma_actual does not, inspect ReactionControl.py, requested_duty, fired_this_step, active_thrusters, and torque_z_from_rcs.

5. If the trajectory shape looks wrong even when control is active, inspect nominal_eom_step in math_3d.py and the observation geometry in control.py.

6. If heating looks symmetric when directional structure is expected, inspect the HeatShield update path and the hot_theta_rad handling in constants.py.

Good workflow for code changes

1. Change one module at a time.
2. Restart the notebook kernel after editing shared modules.
3. Re run the notebook from the import cell downward.
4. Confirm that df is populated before trusting any later plots.
5. Check the final compact summary table and the RCS diagnostic plots after each major change.
6. Save csv outputs after a stable run so later comparisons stay easy.

Suggested future cleanup

1. Rename Pasted code.py to milestone1_nominal.py.
2. Move notebook helper functions into a dedicated utilities module.
3. Add one small runner script that mirrors the notebook simulation loop for batch runs.
4. Keep the notebook for diagnostics and plotting while the runner script handles repeatable experiments.
