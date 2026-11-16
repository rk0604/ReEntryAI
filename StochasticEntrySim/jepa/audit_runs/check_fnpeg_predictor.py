"""Independent check: FNPEG's predictor (RK4, 2 s, filters' K) against the simulator's own equations
(math_3d.nominal_eom_step, Euler at 0.05 s) for constant-bank flights, with the true parameter errors."""
import math
import sys

import numpy as np

sys.path.insert(0, ".")
import constants
from math_3d import nominal_eom_step
from jepa import physics_v4 as PH
from jepa.data import apply_domain, R_EARTH
from jepa.guidance_fnpeg import FNPEGConfig, predict, energy
from jepa.testset_v3 import ATMOSPHERES


def sim(x, sig, p, dt=0.05, hf=7600.0):
    x = list(x)
    s = 0.0
    t = 0.0
    while x[0] - R_EARTH > hf and t < 2500:
        r, V, g = x[0], x[3], x[4]
        out = nominal_eom_step(x=x, sigma_rad=sig, params=p, dt_s=dt)
        xn = out["x_next"]
        s += dt * 0.5 * (V * math.cos(g) / r + xn[3] * math.cos(xn[4]) / xn[0])
        if xn[0] - R_EARTH <= hf:
            f = (x[0] - R_EARTH - hf) / (x[0] - xn[0])
            return s * R_EARTH / 1000 - (1 - f) * dt * xn[3] * math.cos(xn[4]) / xn[0] * R_EARTH / 1000, t + f * dt
        if abs(math.cos(xn[4])) < PH.COS_GAMMA_STOP:
            return s * R_EARTH / 1000, t
        x, t = list(xn), t + dt
    return s * R_EARTH / 1000, t


def main():
    p0 = PH.nominal_params()
    cfg = FNPEGConfig(load_relief=False, gamma_safe_deg=None, dt_s=2.0)
    print(f"{'atm':14s} {'gam0':>5s} {'bank':>5s} {'sim range km':>12s} {'pred km':>9s} {'diff km':>8s} {'t sim':>6s} {'t pred':>6s}")
    worst = 0.0
    for atm, d in ATMOSPHERES.items():
        p = apply_domain(dict(p0), d)
        KD, KL = d["density_scale"] / d["mass_scale"], d["density_scale"] * d["lift_scale"] / d["mass_scale"]
        for gam0 in (-1.2, -2.0, -3.5):
            for bank in (30.0, 60.0, 80.0):
                h0, V0 = 121_920.0, 7_700.0
                x = [R_EARTH + h0, 0.0, 0.0, V0, math.radians(gam0), math.radians(90.0)]
                rs, ts = sim(x, math.radians(bank), p)
                sg = np.array([math.radians(bank)])
                ef = np.array([energy(R_EARTH + cfg.hf_m, 1.0)])     # stop on altitude, not energy
                o = predict(np.array([[R_EARTH + h0, V0, math.radians(gam0), 0.0]]), sg, sg, sg, ef,
                            np.array([KD]), np.array([KL]), cfg)
                rp = -o["z_km"][0]
                worst = max(worst, abs(rp - rs))
                print(f"{atm:14s} {gam0:5.1f} {bank:5.0f} {rs:12.2f} {rp:9.2f} {rp - rs:8.2f} {ts:6.0f} {o['t_end_s'][0]:6.0f}")
    print(f"worst |diff| = {worst:.2f} km")


if __name__ == "__main__":
    main()
