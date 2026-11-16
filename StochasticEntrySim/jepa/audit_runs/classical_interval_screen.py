"""Does the classical baseline's interval heat screen change its flight? rl_env's baseline mode never calls
SimpleBankGuidance.set_prediction_context, so the guidance context keeps its default interval_active=True and
control._compute_sigma_cmd_continuous runs rollout_predictor_corrector_candidate with the interval screen on,
after the bank command is already chosen. This flies the same entries with the screen on (as every study flew
them) and off, and compares every state and the wall time.
Run from StochasticEntrySim/: python jepa/audit_runs/classical_interval_screen.py (jepa/simulator_fidelity_audit.md)"""
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, os.getcwd())
import numpy as np  # noqa: E402

OUT = Path(__file__).resolve().parent


def fly(ic, screen: bool):
    import rl_env
    from jepa.data import NOMINAL_DOMAIN, apply_domain
    env = rl_env.OrionEntryEnv(config_path="configs/default.json", controller_mode="baseline", dispersion=False,
                               max_episode_steps=1400, fast_atmosphere=True)
    apply_domain(env.params, dict(NOMINAL_DOMAIN))
    env.reset(seed=0, options={"ic_override": dict(ic)})
    env.guidance._prediction_context.interval_active = bool(screen)
    used, xs, done, info, t0 = [], [], False, {}, time.time()
    a = np.zeros(1, dtype=np.float32)
    while not done:
        _, _, term, trunc, info = env.step(a)
        used.append(bool(env.guidance.get_last_debug().get("selected_interval_screen_used", False)))
        xs.append(list(env.x) + [float(env.att.sigma_rel_rad)])
        done = bool(term or trunc)
    return {"screen": screen, "screen_used_frac": float(np.mean(used)), "steps": len(xs),
            "wall_s": time.time() - t0, "miss_km": float(info["miss_km"]), "peak_g": float(info["peak_g"]),
            "outcome": info["outcome"]}, np.asarray(xs)


if __name__ == "__main__":
    from jepa import data as D
    ics = D.load_test_ics("configs/heldout_testset.csv", n=3)
    res = []
    for k, ic in enumerate(ics):
        on, xo = fly(ic, True)
        off, xf = fly(ic, False)
        same = bool(xo.shape == xf.shape and np.array_equal(xo, xf))
        r = {"case": k, "on": on, "off": off, "states_bit_identical": same,
             "max_abs_state_diff": float(np.abs(xo - xf).max()) if xo.shape == xf.shape else None}
        print(json.dumps(r), flush=True)
        res.append(r)
    (OUT / "classical_interval_screen.json").write_text(json.dumps(res, indent=1))
