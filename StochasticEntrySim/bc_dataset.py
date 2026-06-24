"""
bc_dataset.py

Collect a behavior cloning dataset by running OrionEntryEnv in baseline mode
(the classical predictor corrector) and recording, at every guidance step:

    observation (exactly what the RL policy sees) and the classical bank command

The bank command is the guidance chosen bank for that step, read from
SimpleBankGuidance.get_last_debug()["chosen_sigma_cmd_deg"]. It is mapped into
the policy action space the same way the environment maps it back:

    env.step:   sigma_cmd_rad = action * (pi/2)        with action in [-1, 1]
    therefore:  action        = sigma_cmd_rad / (pi/2) = sigma_cmd_deg / 90

These (observation, action) pairs allow a supervised pretrain step (behavior
cloning) of the policy, so it starts already knowing how to steer instead of
discovering steering by exploration. The Colab imitation warmstart cells use
this dataset.

This file is standalone and importable. It does not depend on the RL training
loop.
"""

from __future__ import annotations

import math
from typing import List, Optional, Tuple

import numpy as np

# A dummy action is required by the Gym step() signature but is ignored in
# baseline mode, where the predictor computes the bank.
_DUMMY_ACTION = np.zeros((1,), dtype=np.float32)
_SIGMA_MAX_DEG = 90.0   # action of plus or minus 1 maps to bank of plus or minus 90 deg


def _run_one_episode(args: Tuple[str, int, bool]):
    """
    Run one baseline episode and return its recorded pairs.

    Parameters:
        args: a tuple (config_path, seed, dispersion).
            config_path: mission config the baseline flies.
            seed:        reset seed for this episode.
            dispersion:  draw a dispersed initial condition when True.

    Returns:
        A tuple (obs, act, info):
            obs:  float32 array of shape [N, 13], one observation per guidance
                  step, captured before the step is taken.
            act:  float32 array of shape [N, 1], the classical bank for each
                  step mapped into the [-1, 1] action space.
            info: the final info dictionary from the environment.
    """
    import rl_env   # imported inside the worker so it is picklable across processes
    config_path, seed, dispersion = args
    env = rl_env.OrionEntryEnv(config_path=config_path, controller_mode="baseline",
                               dispersion=dispersion, fast_atmosphere=True)
    obs, info = env.reset(seed=seed)
    obs_rows: List[np.ndarray] = []
    act_rows: List[List[float]] = []
    done = False
    while not done:
        # Capture the observation the policy would act on, before stepping.
        obs_at_decision = np.asarray(obs, dtype=np.float32)
        obs, _r, term, trunc, info = env.step(_DUMMY_ACTION)
        # Read the bank the classical guidance actually chose this step.
        dbg = env.guidance.get_last_debug() or {}
        sigma_cmd_deg = float(dbg.get("chosen_sigma_cmd_deg", 0.0))
        # Map degrees into the [-1, 1] action space and clamp.
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
    """
    Collect observation and action pairs from several baseline episodes.

    Parameters:
        config_path: mission config the baseline flies.
        n_episodes:  number of dispersed baseline trajectories to imitate.
        dispersion:  draw a fresh dispersed initial condition each episode.
                     True is recommended so the clone generalizes across entry
                     conditions.
        workers:     number of parallel worker processes. Each worker runs whole
                     episodes. A value of 1 runs everything in one process.
        seed0:       base seed. Episode i uses seed0 + i.
        out_path:    when given, the dataset is saved to this .npz file with
                     keys obs and act.

    Returns:
        A tuple (obs, act):
            obs: float32 array of shape [Total, 13].
            act: float32 array of shape [Total, 1].
    """
    # One work item per episode, with a distinct seed.
    args = [(config_path, seed0 + i, dispersion) for i in range(int(n_episodes))]

    # Run in a process pool when more than one worker is requested.
    if workers and workers > 1:
        import multiprocessing as mp
        with mp.Pool(int(workers)) as pool:
            results = pool.map(_run_one_episode, args)
    else:
        results = [_run_one_episode(a) for a in args]

    # Stack every episode into one dataset.
    obs = np.concatenate([r[0] for r in results], axis=0)
    act = np.concatenate([r[1] for r in results], axis=0)

    # Quick report so it is clear whether the baseline produced good labels.
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
