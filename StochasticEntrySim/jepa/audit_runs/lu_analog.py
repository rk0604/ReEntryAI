"""Lu (2014) load-relief result (Figs 14 and 17) reproduced on a Lu-analog mission in our simulator.
Lu's own CEV case (rotating Earth, 5.5 m / 8,383 kg capsule, 3,769 km to Edwards) is out of reach of our vehicle on
a non-rotating Earth (lu_cev_search.json: full lift-up falls 16-880 km short for gamma0 -1.0 to -2.0 deg). Here:
Lu's entry state (121.9 km, 7,623.5 m/s, Vf 150 m/s) with ORB-1's gamma0 -1.40 deg, and target ranges short enough
that the unrelieved flight peaks near Lu's 5.4 g.
Run from StochasticEntrySim/ (results are written next to this script; see jepa/guidance_audit.md section 3):
       python jepa/audit_runs/lu_analog.py stage1        (1,700-2,450 km x {none, 4 g, 3 g}, 180 deg cap)
       python jepa/audit_runs/lu_analog.py stage2 2450   (sigma_f 30/60 deg at 180 deg; constant and 30 deg at 90 deg)"""
import json
import os
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
GAM0 = -1.40


def target_cfg_path(rng_km):
    from jepa import benchmark_v5 as B
    ic = dict(B.ORB1Spec().nominal_ic(), gamma_deg=GAM0)
    lat, lon = B.great_circle_point(ic, float(rng_km))
    cfg = B.target_config({"target_index": 9, "frac": 0.0, "downrange_km": float(rng_km), "phi_deg": lat,
                           "lambda_deg": lon})
    p = Path(tempfile.gettempdir()) / f"lu_analog_target_{int(rng_km)}_{os.getpid()}.json"
    p.write_text(json.dumps(cfg, indent=1))
    return str(p), ic


def fly(args):
    rng_km, amax, cap = args[:3]
    sigf = args[3] if len(args) > 3 else None
    relief = amax is not None
    from jepa import fly_v5 as F5
    from jepa.data import NOMINAL_DOMAIN, physical_state
    from jepa.guidance_fnpeg import FNPEGConfig
    cpath, ic = target_cfg_path(rng_km)
    job = {"ic": ic, "domain": dict(NOMINAL_DOMAIN), "config_path": cpath, "sigma_max_deg": float(cap)}
    law = F5.FNPEGLaw(FNPEGConfig(sigma_max_deg=float(cap), load_relief=relief, sigma_f_deg=sigf,
                                  amax_g=float(amax if relief else 4.5), vf_mps=150.0))
    f = F5.GuidanceFlight(job)
    f.state = law.new_flight(job, f.env)
    while not f.done:
        c = law.decide([f.state], [physical_state(f.env)], [float(f.env.t_s)], [f.env], [f.k])[0]
        f.step(c)
    cmds = np.degrees(np.abs(np.asarray(f.cmds)))
    return {"range_km": rng_km, "amax": amax, "cap": cap, "sigma_f": sigf, "outcome": f.info.get("outcome"),
            "miss_km": float(f.info.get("miss_km", float("nan"))), "peak_g": float(f.env.peak_g),
            "t_s": float(f.env.t_s), "max_bank": float(cmds.max()), "n_over_90": int((cmds > 90.5).sum())}


def main():
    sys.path.insert(0, os.getcwd())
    phase = sys.argv[1]
    if phase == "stage1":
        runs = [(r, a, 180.0) for r in (1700.0, 1950.0, 2200.0, 2450.0) for a in (None, 4.0, 3.0)]
    else:
        r = float(sys.argv[2])
        runs = [(r, a, 180.0, 30.0) for a in (None, 4.0, 3.0, 2.7)] + [(r, a, 180.0, 60.0) for a in (4.0, 3.0)]
        runs += [(r, a, 90.0, None) for a in (None, 4.0, 3.0)] + [(r, a, 90.0, 30.0) for a in (4.0, 3.0)]
    with Pool(min(len(runs), max(1, (os.cpu_count() or 2) - 1))) as pool:
        res = pool.map(fly, runs)
    for x in res:
        print(json.dumps(x), flush=True)
    (OUT / f"lu_analog_{phase}.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
