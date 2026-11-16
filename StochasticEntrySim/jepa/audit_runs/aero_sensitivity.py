"""Range effect of our orion_cm_trim schedule vs Bibb Part I's nominal trim L/D (Fig. 27) and trim CD (Fig. 26):
ORB-1 nominal entry (121.92 km, 7,623.5 m/s, -1.40 deg) at constant bank, full simulator."""
import json
import math
import os
import sys

sys.path.insert(0, os.getcwd())
import numpy as np  # noqa: E402

import math_3d  # noqa: E402,F401
from jepa import benchmark_v5 as B  # noqa: E402
from jepa import testset_v3 as T3  # noqa: E402
from jepa.data import SIGMA_MAX_RAD, NOMINAL_DOMAIN, apply_domain  # noqa: E402

OURS = [(0.0, 0.83, 0.10), (0.5, 0.85, 0.15), (0.7, 0.88, 0.20), (0.9, 1.02, 0.23), (1.1, 1.28, 0.25),
        (1.5, 1.38, 0.28), (2.0, 1.43, 0.30), (3.0, 1.47, 0.33), (5.0, 1.50, 0.36), (8.0, 1.52, 0.40),
        (10.0, 1.50, 0.39), (15.0, 1.52, 0.33), (20.0, 1.54, 0.30), (25.0, 1.55, 0.27), (30.0, 1.55, 0.27),
        (35.0, 1.55, 0.28), (40.0, 1.55, 0.30)]
BIBB_LD = {8.0: 0.280, 10.0: 0.277, 15.0: 0.263, 20.0: 0.257, 25.0: 0.250, 30.0: 0.245, 35.0: 0.244, 40.0: 0.245}


def schedule(ld: bool, cd: bool):
    out = []
    for f, c, l in OURS:
        if f >= 8.0:
            out.append((f, 1.38 if cd else c, BIBB_LD[f] if ld else l))
        else:
            out.append((f, c, l))
    return out


def fly(sigma, sch):
    import rl_env
    ic = B.ORB1Spec().nominal_ic()
    env = rl_env.OrionEntryEnv(config_path="configs/default.json", controller_mode="policy", dispersion=False,
                               max_episode_steps=1400, fast_atmosphere=True)
    apply_domain(env.params, dict(NOMINAL_DOMAIN))
    if sch is not None:
        env.params["orion_trim_schedule"] = sch
    env.reset(seed=0, options={"ic_override": dict(ic)})
    a = np.array([math.radians(sigma) / SIGMA_MAX_RAD], dtype=np.float32)
    done, info = False, {}
    while not done:
        _, _, term, trunc, info = env.step(a)
        done = bool(term or trunc)
    Q = T3.entry_frame(math.radians(ic["phi_deg"]), math.radians(ic["lambda_deg"]), math.radians(ic["chi_deg"]))
    dr, _ = T3.to_local(Q, float(env.x[1]), float(env.x[2]))
    return {"downrange_km": round(dr, 1), "peak_g": round(float(info.get("peak_g", float("nan"))), 2),
            "t_s": float(env.t_s), "outcome": str(info.get("outcome"))}


if __name__ == "__main__":
    res = {}
    for name, sch in (("ours (default)", None), ("ours, re-entered", OURS), ("Bibb L/D", schedule(True, False)),
                      ("Bibb L/D and CD 1.38", schedule(True, True))):
        for sigma in (0.0, 60.0):
            r = fly(sigma, sch)
            res[f"{name} | bank {sigma:g}"] = r
            print(name, sigma, r, flush=True)
    open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "aero_sensitivity.json"), "w").write(
        json.dumps(res, indent=1))
