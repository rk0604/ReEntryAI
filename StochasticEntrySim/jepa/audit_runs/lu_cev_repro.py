"""Lu (2014) CEV mission (Table 1: h0 121.9 km, V0 7623.5 m/s, s0 2035 nmi, hf 7.6 km, Vf 150 m/s) flown by our
FNPEG in the full simulator. Lu gives no gamma0, so phase 'search' finds the gamma0 whose unconstrained flight peaks
near Lu's 5.4 g; phase 'sweep' then flies Lu's load limits (none, 4, 3, 2.7, 2.5, 2.3 g) with 90 and 180 deg caps.
Only 'search' was run (2026-11-15): no gamma0 from -1.0 to -2.0 deg reaches the target, so 'sweep' has nothing to
reproduce (jepa/guidance_audit.md section 2).
Run from StochasticEntrySim/: python jepa/audit_runs/lu_cev_repro.py search | sweep <gamma0_deg>"""
import json
import os
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
S0_KM = 2035 * 1.852


def target_cfg_path(gam0):
    from jepa import benchmark_v5 as B
    ic = dict(B.ORB1Spec().nominal_ic(), gamma_deg=float(gam0))
    lat, lon = B.great_circle_point(ic, S0_KM)
    cfg = B.target_config({"target_index": 9, "frac": 0.0, "downrange_km": S0_KM, "phi_deg": lat, "lambda_deg": lon})
    p = Path(tempfile.gettempdir()) / f"lu_cev_target_{os.getpid()}.json"
    p.write_text(json.dumps(cfg, indent=1))
    return str(p), ic


def fly(args):
    gam0, amax, cap, relief = args
    from jepa import fly_v5 as F5
    from jepa.data import NOMINAL_DOMAIN
    from jepa.guidance_fnpeg import FNPEGConfig
    cpath, ic = target_cfg_path(gam0)
    job = {"ic": ic, "domain": dict(NOMINAL_DOMAIN), "config_path": cpath, "sigma_max_deg": float(cap)}
    law = F5.FNPEGLaw(FNPEGConfig(sigma_max_deg=float(cap), load_relief=bool(relief),
                                  amax_g=float(amax if relief else 4.5), vf_mps=150.0))
    f = F5.GuidanceFlight(job)
    f.state = law.new_flight(job, f.env)
    from jepa.data import physical_state
    while not f.done:
        c = law.decide([f.state], [physical_state(f.env)], [float(f.env.t_s)], [f.env], [f.k])[0]
        f.step(c)
    cmds = np.degrees(np.abs(np.asarray(f.cmds)))
    return {"gamma0": gam0, "amax": amax if relief else None, "cap": cap, "outcome": f.info.get("outcome"),
            "miss_km": float(f.info.get("miss_km", float("nan"))), "peak_g": float(f.env.peak_g),
            "t_s": float(f.env.t_s), "max_bank": float(cmds.max()), "s_over_90": int((cmds > 90.5).sum())}


def main():
    import os
    sys.path.insert(0, os.getcwd())
    phase = sys.argv[1]
    if phase == "search":
        runs = [(g, None, 180.0, False) for g in (-1.0, -1.1, -1.2, -1.3, -1.4, -1.5, -1.7, -2.0)]
    else:
        g = float(sys.argv[2])
        runs = [(g, None, cap, False) for cap in (90.0, 180.0)]
        runs += [(g, a, cap, True) for a in (4.0, 3.0, 2.7, 2.5, 2.3) for cap in (90.0, 180.0)]
    with Pool(min(len(runs), max(1, (os.cpu_count() or 2) - 1))) as pool:
        res = pool.map(fly, runs)
    for r in res:
        print(json.dumps(r), flush=True)
    (OUT / f"lu_cev_{phase}.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
