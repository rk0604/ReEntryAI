import sys
import os
sys.path.append(r"C:\Users\Risha\desktop\ReEntryAI\StochasticEntrySim")
import math_3d

params_trim = {
    "aero_model": "orion_cm_trim"
}

params_scheduled = {
    "aero_model": "scheduled_orion_like"
}

params_constant = {
    "aero_model": "constant",
    "CD_const": 1.15,
    "CL_const": 0.28
}

test_speeds = [100.0, 300.0, 600.0, 1500.0, 2440.0, 3050.0, 5000.0, 7500.0, 10500.0, 12000.0]

print("--- Orion CM Trim Database Model (Bibb et al. hypersonic) ---")
print(f"{'V_mps':>7} | {'FMV':>6} | {'CD':>5} | {'L/D':>5} | {'CL':>5} | region")
for v in test_speeds:
    res = math_3d.aero_coefficients_from_speed_altitude(v, 50000.0, params_trim)
    print(f"{v:>7.0f} | {res['FMV']:>6.2f} | {res['CD']:>5.3f} | {res['LD']:>5.3f} | {res['CL']:>5.3f} | {res['schedule_region']}")

print("\n--- Old Scheduled Orion-like Model ---")
for v in test_speeds:
    res = math_3d.aero_coefficients_from_speed_altitude(v, 50000.0, params_scheduled)
    print(f"V: {v:>7} m/s | CD: {res['CD']:.3f} | L/D: {res['LD']:.3f} | CL: {res['CL']:.3f} | Region: {res['schedule_region']}")

print("\n--- Constant Model Fallback ---")
res = math_3d.aero_coefficients_from_speed_altitude(3500.0, 50000.0, params_constant)
print(f"V: 3500.0 m/s | CD: {res['CD']:.3f} | L/D: {res['LD']:.3f} | CL: {res['CL']:.3f} | Region: {res['schedule_region']}")
