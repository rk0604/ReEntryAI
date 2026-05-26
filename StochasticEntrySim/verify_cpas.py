"""
Standalone verification of the CPAS module.

Runs a synthetic vertical descent from 8 km / 200 m/s with no aero
lift and no guidance, just gravity + capsule aero drag + chute drag.
Prints phase transitions and final terminal velocity.
"""

import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import math
import cpas
import constants

# --- Capsule parameters (match the main sim) ---
mass_kg = float(constants.CAPSULE_MASS_KG)
ref_area_m2 = float(constants.CAPSULE_REFERENCE_AREA_M2)
cd_capsule = 1.20  # representative subsonic CD for the bare capsule

# --- Simple exponential atmosphere for the demo (sea-level rho, scale height) ---
RHO0 = 1.225
SCALE_HEIGHT_M = 8500.0

def rho_at(alt_m: float) -> float:
    return RHO0 * math.exp(-max(0.0, alt_m) / SCALE_HEIGHT_M)

# --- Initial conditions ---
alt_m = 8000.0
V_mps = 200.0   # falling vertically downward
t_s = 0.0
dt = 0.1
g = 9.81

c = cpas.CPAS()

print(f"{'t_s':>6} {'alt_m':>8} {'V_mps':>7} {'phase':>8} {'open':>5} "
      f"{'cdA':>8} {'F_drag_N':>10} {'events'}")
print("-" * 90)

last_print_t = -1.0
events_seen = []
while alt_m > 0.0 and t_s < 600.0:
    out = c.step(t_s=t_s, alt_m=alt_m, V_mps=V_mps, mach=None)

    rho = rho_at(alt_m)
    q_bar = 0.5 * rho * V_mps * V_mps

    # Bare capsule drag + chute drag, both opposite to velocity (downward fall).
    drag_capsule_N = q_bar * cd_capsule * ref_area_m2
    drag_chute_N = q_bar * out.drag_area_cdA_m2
    drag_total_N = drag_capsule_N + drag_chute_N

    # Vertical EOM: a = g - D/m  (V_mps is the downward fall speed)
    accel = g - drag_total_N / mass_kg
    V_mps = V_mps + accel * dt
    if V_mps < 0.0:
        V_mps = 0.0
    alt_m = alt_m - V_mps * dt
    t_s += dt

    if out.events:
        events_seen.extend([(t_s, alt_m, V_mps, e) for e in out.events])
        print(f"{t_s:6.1f} {alt_m:8.1f} {V_mps:7.2f} {out.phase:>8} "
              f"{out.open_fraction:5.2f} {out.drag_area_cdA_m2:8.2f} "
              f"{drag_total_N:10.1f}  {' '.join(out.events)}")
        last_print_t = t_s
    elif t_s - last_print_t >= 10.0:
        print(f"{t_s:6.1f} {alt_m:8.1f} {V_mps:7.2f} {out.phase:>8} "
              f"{out.open_fraction:5.2f} {out.drag_area_cdA_m2:8.2f} "
              f"{drag_total_N:10.1f}")
        last_print_t = t_s

print("-" * 90)
print(f"Landing summary: t = {t_s:.2f} s, V_impact = {V_mps:.2f} m/s")

# Analytic terminal velocity under mains, sea-level rho:
v_t = cpas.terminal_velocity_mps(c.config, mass_kg=mass_kg, rho_kgm3=RHO0)
print(f"Analytic terminal V under mains at sea level: {v_t:.2f} m/s "
      f"(target ~7.6 m/s)")

print()
print("CPAS summary:", c.summary())
