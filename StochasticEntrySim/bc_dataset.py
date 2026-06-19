"""
bc_dataset.py

Collect a behavior-cloning dataset by running OrionEntryEnv in *baseline* mode
(the classical predictor-corrector) and recording, at every guidance step:

    observation (exactly what the RL policy sees)  ->  classical bank command

The bank command is the guidance's chosen bank for that step
(`SimpleBankGuidance.get_last_debug()["chosen_sigma_cmd_deg"]`), mapped into the
policy's action space the same way the env maps it back:

    env.step:   sigma_cmd_rad = action * (pi/2)        # action in [-1, 1]
    so:         action        = sigma_cmd_rad / (pi/2) = sigma_cmd_deg / 90

These (obs, action) pairs let us supervised-pretrain (behavior-clone) the policy
so it *starts* knowing how to steer, instead of having to discover steering by
exploration. Used by the Colab imitation-warmstart cells.

Standalone and importable; no dependency on the RL training loop.
"""

from __future__ import annotations

import math
from typing import List, Optional, Tuple

import numpy as np

# A dummy action is required by the Gym step() signature but is IGNORED in
# baseline mode (the predictor computes the bank).
_DUMMY_ACTION = np.zeros((1,), dtype=np.float32)
_SIGMA_MAX_DEG = 90.0   # action +-1 maps to bank +-90 deg


def _run_one_episode(args: Tuple[str, int, bool]):
    """Run one baseline episode; return (obs[N,13], act[N,1], info)."""
    import rl_env   # imported inside the worker (picklable across processes)
    config_path, seed, dispersion = args
    env = rl_env.OrionEntryEnv(config_path=config_path, controller_mode="baseline",
                               dispersion=dispersion, fast_atmosphere=True)
    obs, info = env.reset(seed=seed)
    obs_rows: List[np.ndarray] = []
    act_rows: List[List[float]] = []
    done = False
    while not done:
        obs_at_decision = np.asarray(obs, dtype=np.float32)
        obs, _r, term, trunc, info = env.step(_DUMMY_ACTION)
        dbg = env.guidance.get_last_debug() or {}
        sigma_cmd_deg = float(dbg.get("chosen_sigma_cmd_deg", 0.0))
        action = max(-1.0, min(1.0, sigma_cmd_deg / _SIGMA_MAX_DEG))
        obs_rows.append(obs_at_decision)
        act_rows.append([action])
        done = bool(term or trunc)
    return (np.asarray(obs_rows, dtype=np.float32),
            np.asarray(act_rows, dtype=np.float32), info)


def collect_bc_dataset(config_path: str = "configs/default.json",
                       n_episodes: int = 60, dispersion: bool = True,
                       workers: int = 1, seed0: int = 12345,
                       out_path: Optional[str] = None):
    """Collect (obs, action) pairs from `n_episodes` baseline episodes.

    Args:
        config_path: mission config the baseline flies.
        n_episodes:  number of dispersed baseline trajectories to imitate.
        dispersion:  draw a fresh dispersed IC each episode (recommended True so
                     the clone generalizes across entry conditions).
        workers:     parallel processes (each runs whole episodes); 1 = serial.
        out_path:    if given, save to this .npz (keys: obs, act).

    Returns (obs[Total,13], act[Total,1]).
    """
    args = [(config_path, seed0 + i, dispersion) for i in range(int(n_episodes))]

    if workers and workers > 1:
        import multiprocessing as mp
        with mp.Pool(int(workers)) as pool:
            results = pool.map(_run_one_episode, args)
    else:
        results = [_run_one_episode(a) for a in args]

    obs = np.concatenate([r[0] for r in results], axis=0)
    act = np.concatenate([r[1] for r in results], axis=0)

    # Quick report: did the baseline episodes actually succeed (good labels)?
    succ = float(np.mean([1.0 if r[2].get("success") else 0.0 for r in results]))
    miss = float(np.median([r[2].get("miss_km", float("nan")) for r in results]))
    print(f"[bc_dataset] {len(results)} episodes | baseline success {100*succ:.0f}%"
          f" | median miss {miss:.2f} km | {len(obs)} (obs, action) pairs")

    if out_path:
        np.savez_compressed(out_path, obs=obs, act=act)
        print(f"[bc_dataset] saved -> {out_path}")
    return obs, act


if __name__ == "__main__":
    import sys
    cfg = sys.argv[1] if len(sys.argv) > 1 else "configs/default.json"
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    O, A = collect_bc_dataset(cfg, n_episodes=n, dispersion=True, workers=1)
    print("obs", O.shape, "act", A.shape,
          "| action range", float(A.min()), float(A.max()))
