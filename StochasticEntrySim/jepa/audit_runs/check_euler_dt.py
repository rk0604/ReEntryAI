"""Integration error of the simulator's forward-Euler step (0.25 s) against 0.05 s and 0.01 s, constant bank.
Run from StochasticEntrySim/: python jepa/audit_runs/check_euler_dt.py (jepa/simulator_fidelity_audit.md)."""
import math
import sys

sys.path.insert(0, ".")
from check_fnpeg_predictor import sim  # noqa: E402  (same folder as this script)
from jepa import physics_v4 as PH
from jepa.data import R_EARTH


def main():
    p = PH.nominal_params()
    print(f"{'gam0':>5s} {'bank':>5s} {'R(0.25)':>9s} {'R(0.05)':>9s} {'R(0.01)':>9s} {'err 0.25 km':>11s} {'err 0.05 km':>11s}")
    for gam0 in (-1.2, -2.0, -3.5):
        for bank in (30.0, 60.0, 80.0):
            x = [R_EARTH + 121_920.0, 0.0, 0.0, 7_700.0, math.radians(gam0), math.radians(90.0)]
            r = {dt: sim(x, math.radians(bank), p, dt=dt)[0] for dt in (0.25, 0.05, 0.01)}
            print(f"{gam0:5.1f} {bank:5.0f} {r[0.25]:9.2f} {r[0.05]:9.2f} {r[0.01]:9.2f} "
                  f"{r[0.25] - r[0.01]:11.2f} {r[0.05] - r[0.01]:11.2f}")


if __name__ == "__main__":
    main()
