"""
Standalone verification of the high-fidelity CPAS module.

Runs a synthetic vertical descent from 10 km and demonstrates every
modeled feature:
  - FBC jettison + mass shed
  - 2-stage drogue reefing
  - 3-stage main reefing
  - Mortar / line-stretch initial delay
  - Opening-shock transients
  - Pendulum oscillation
  - Squidding events (random)
  - Wind drift (lateral position offset)
  - 2-of-3 main failure mode (configurable)

Prints a timeline of every CPAS event, periodic state samples, and a
landing summary that should match real Orion within a few percent.
"""

import math
import os
import random
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import constants
import cpas

mass_kg_initial = float(constants.CAPSULE_MASS_KG)
ref_area_m2 = float(constants.CAPSULE_REFERENCE_AREA_M2)
cd_capsule = 1.20

RHO0 = 1.225
SCALE_HEIGHT_M = 8500.0
def rho_at(alt_m):
    return RHO0 * math.exp(-max(0.0, alt_m) / SCALE_HEIGHT_M)


def run_scenario(label, num_mains_operational=3, wind_east=5.0, wind_north=2.0,
                 squid_prob=0.005, rng_seed=0, force_deterministic_failure=False):
    print()
    print("=" * 90)
    print(f"SCENARIO: {label}")
    print("=" * 90)

    cfg = cpas.CPASConfig(
        # Disable stochastic roll so num_mains_operational is honored
        # deterministically. This makes the test repeatable for whatever
        # failure mode the scenario wants to exercise.
        enable_stochastic_failure=not force_deterministic_failure
            and num_mains_operational == 3,
        num_mains_operational=num_mains_operational,
        wind_east_mps=wind_east,
        wind_north_mps=wind_north,
        squid_probability_per_second=squid_prob,
    )
    c = cpas.CPAS(config=cfg, rng=random.Random(rng_seed))

    # Initial state: 10 km, falling at 200 m/s, at (0, 0) lat/lon.
    alt_m = 10000.0
    V_mps = 200.0
    t_s = 0.0
    dt = 0.1
    g = 9.81
    east_m = 0.0      # accumulated horizontal drift (wind + pendulum)
    north_m = 0.0
    mass_kg = mass_kg_initial

    print(f"{'t':>6}  {'alt':>7}  {'V':>6}  {'phase':>8}  {'stage':>5}  "
          f"{'cdA':>8}  {'pend_deg':>8}  {'mass':>7}  events")
    print("-" * 90)

    last_print_t = -10.0
    max_pendulum_amplitude = 0.0
    while alt_m > 0.0 and t_s < 1200.0:
        out = c.step(t_s=t_s, alt_m=alt_m, V_mps=V_mps, dt_s=dt)
        max_pendulum_amplitude = max(max_pendulum_amplitude,
                                     abs(out.pendulum_angle_rad))

        # Vertical descent dynamics
        rho = rho_at(alt_m)
        q_bar = 0.5 * rho * V_mps * V_mps
        drag_total_N = q_bar * (cd_capsule * ref_area_m2 + out.drag_area_cdA_m2)
        accel = g - drag_total_N / mass_kg
        V_mps = max(0.0, V_mps + accel * dt)
        alt_m -= V_mps * dt

        # Horizontal motion (just for diagnostics here)
        east_m += dt * (out.wind_east_mps
                        + out.pendulum_lateral_v_mps * math.sin(out.pendulum_azimuth_rad))
        north_m += dt * (out.wind_north_mps
                         + out.pendulum_lateral_v_mps * math.cos(out.pendulum_azimuth_rad))

        # Apply mass shed from FBC
        if out.mass_shed_delta_kg > 0.0:
            mass_kg = max(1.0, mass_kg - out.mass_shed_delta_kg)

        t_s += dt
        if out.events or t_s - last_print_t >= 15.0:
            print(f"{t_s:6.1f}  {alt_m:7.0f}  {V_mps:6.2f}  {out.phase:>8}  "
                  f"{out.reefing_stage:>5}  {out.drag_area_cdA_m2:8.1f}  "
                  f"{math.degrees(out.pendulum_angle_rad):>8.2f}  "
                  f"{mass_kg:7.0f}  {' '.join(out.events)}")
            last_print_t = t_s

    print("-" * 90)
    print(f"Landing: t={t_s:.1f}s  V_impact={V_mps:.2f} m/s")
    print(f"Drift:   east={east_m:7.1f} m   north={north_m:7.1f} m   "
          f"total={(east_m**2 + north_m**2) ** 0.5:.1f} m")
    print(f"Max pendulum amplitude observed: {math.degrees(max_pendulum_amplitude):.2f} deg")
    print(f"Mass shed total: {mass_kg_initial - mass_kg:.1f} kg")
    print(f"Summary: {c.summary()}")

    # Analytic terminal velocity for this configuration
    n_ops = cfg.num_mains_operational
    v_t = cpas.terminal_velocity_mps(cfg, mass_kg=mass_kg, rho_kgm3=RHO0,
                                     n_mains_active=n_ops)
    print(f"Analytic terminal V (sea-level, {n_ops} mains): {v_t:.2f} m/s")


if __name__ == "__main__":
    # Stochastic mode: per-mission Fuqua 2010 prior. Nominal in 99%+ of seeds.
    run_scenario("Stochastic prior: 3 mains expected (Fuqua 2010 prior)",
                 num_mains_operational=3,
                 wind_east=5.0, wind_north=2.0, squid_prob=0.005, rng_seed=0)
    # Deterministically forced 2-of-3 failure -- exercises the pendulum mode.
    run_scenario("Forced 2-of-3 main failure (pendulum failure mode)",
                 num_mains_operational=2,
                 wind_east=8.0, wind_north=-3.0, squid_prob=0.005, rng_seed=1,
                 force_deterministic_failure=True)
    # Calm air, 3 mains -- shows the baseline descent profile.
    run_scenario("Calm air, 3 mains (baseline)",
                 num_mains_operational=3, wind_east=0.0, wind_north=0.0,
                 squid_prob=0.0, rng_seed=2)
