"""
What v5's RL work costs on the current runtime (the Colab-vs-Kaggle decision, jepa/kaggle_rl.md).

Times short pieces of the real training code with the study's frozen settings, then projects every Tier 1 RL
tuning configuration and training run:

    PPO          decisions per second (ppo_hist; bcppo_fnpeg's PPO with its KL reference). A PPO run stops at
                 budget_decisions or budget_hours of training, whichever first, so its hours are
                 min(budget_hours, budget_decisions / rate).
    FNPEG expert one chunk of expert-labelled fast-env episodes (DAgger's iteration, BC's dataset chunk) at the
                 training size (dagger_rollouts envs) and the tuning size (10 % of it). Not wall-clock capped.
    fits         the supervised fit on one chunk's states at dagger_epochs (DAgger refits on all data so far
                 before every iteration; BC fits once on its whole dataset at bc_epochs).
    flights      DEV48 policy flights in the full simulator (48 per tuning configuration), on every CPU core.
    fast eval    each training run's fast-environment evaluation (fast_eval.json).

Writes only into its own scratch folder (never into a study folder) and returns a dict; the notebook
ReEntryAI_JEPA_v5_rl_cost.ipynb saves it to Drive. FNPEG flies its default settings (tuning may change its RK4
step, so the expert's cost is approximate), and every tuning configuration is assumed to run at the default
configuration's speed.
"""

from __future__ import annotations

import math
import time
from dataclasses import replace
from pathlib import Path
from typing import Dict

import numpy as np
import torch

from jepa import rl_v5 as RL
from jepa import study_v5 as S5


def units_per_hour(gpu_name: str) -> float:
    for k, v in S5.UNITS_PER_HOUR.items():
        if k != "cpu" and k.lower() in gpu_name.lower():
            return float(v)
    return float(S5.UNITS_PER_HOUR["cpu"])


def _cfg(cfg: S5.StudyConfigV5, arm: str, **kw) -> RL.RLCfg:
    d = {**cfg.rl, **cfg.rl_arms.get(arm, {}), **kw}
    if "hidden" in d:
        d["hidden"] = tuple(d["hidden"])
    return RL.RLCfg(**d)


def measure(cfg: S5.StudyConfigV5, L: S5.LayoutV5, device: str, workers: int, minutes: float = 5.0,
            chunk: int = 0, tune_chunk: int = 0, n_flights: int = 0, fast_eval_episodes: int = 0,
            ppo_kw: Dict = None, log=print) -> dict:
    """Time each piece on this runtime. chunk/tune_chunk/n_flights/fast_eval_episodes 0 = the study's sizes."""
    dev = torch.device(device)
    ppo_kw = dict(ppo_kw or {})
    root = Path(L.root) / "rl_cost"
    fn = dict(cfg.fnpeg)
    base = _cfg(cfg, "ppo_hist")
    chunk = int(chunk or base.dagger_rollouts)
    tune_chunk = int(tune_chunk or max(2, int(base.dagger_rollouts * cfg.rl_tune_frac)))
    out = {"hardware": RL.hardware_info(device), "minutes_per_ppo_measurement": minutes, "chunk": chunk,
           "tune_chunk": tune_chunk}

    def ppo_rate(arm: str, actor=None, ref=None) -> dict:
        rc = _cfg(cfg, arm, budget_hours=minutes / 60.0, curve_points=0, save_every_s=1e9, seed=0, **ppo_kw)
        d = root / arm
        t0 = time.time()
        info = RL.ppo_train(L, d, rc, device, f"measure {arm}", actor=actor, ref=ref, log=log)
        r = {"decisions": info["decisions"], "train_s": info["train_s"], "wall_s": round(time.time() - t0, 1),
             "decisions_per_s": info["decisions"] / max(info["train_s"], 1e-9), "n_envs": rc.n_envs}
        log(f"[{arm}] {r['decisions']:,} decisions in {r['train_s']:.0f} s = {r['decisions_per_s']:,.0f} decisions/s")
        return r

    out["ppo_hist"] = ppo_rate("ppo_hist")
    bc = _cfg(cfg, "bcppo_fnpeg")
    clone = RL.Actor(hidden=bc.hidden, kind="magsign").to(dev)
    ref = RL.Actor(hidden=bc.hidden, kind="magsign").to(dev)
    ref.load_state_dict(clone.state_dict())
    for p in ref.parameters():
        p.requires_grad_(False)
    out["bcppo_fnpeg_ppo"] = ppo_rate("bcppo_fnpeg", actor=clone, ref=ref)

    def expert_chunk(k: int) -> dict:
        env = RL.make_env(base, 0, dev, n_envs=k)
        expert = RL.FastExpert(env, base, fn)
        gen = torch.Generator(device=dev).manual_seed(9)
        t0 = time.time()
        O, A, s = RL._rollouts(env, expert, None, 1.0, k, gen, dev)
        r = {"episodes": int(s["episodes"]), "states": int(len(O)), "s": round(time.time() - t0, 1)}
        log(f"[FNPEG expert] {r['episodes']} episodes on {k} envs: {r['states']:,} labelled states in {r['s']:.0f} s")
        return r, O, A

    out["expert_chunk"], O, A = expert_chunk(chunk)
    out["expert_tune_chunk"], _, _ = expert_chunk(tune_chunk)

    g = torch.Generator(device=dev).manual_seed(1)
    actor = RL.Actor(hidden=base.hidden, kind="magsign").to(dev)
    t0 = time.time()
    RL.fit_supervised(actor, O, A, base.dagger_epochs, base.sup_batch, base.lr, g, dev, log, "measure fit")
    out["fit_chunk"] = {"states": int(len(O)), "epochs": base.dagger_epochs, "s": round(time.time() - t0, 1)}

    pol_path = root / "ppo_hist" / "policy.pt"
    conds = S5.conditions(cfg, "DEV48")
    n = int(n_flights or min(len(conds), 4 * max(1, workers)))
    js = [S5.policy_job(c, "rlcost_ppo", "rl_cost", str(pol_path), 0) for c in conds[:n]]
    t0 = time.time()
    S5.run_cpu(cfg, js, L, workers, log, "measure policy flights")
    out["flights"] = {"n": n, "workers": workers, "s": round(time.time() - t0, 1),
                      "s_per_flight": round((time.time() - t0) / max(n, 1), 2)}
    log(f"[policy flights] {n} DEV48 flights on {workers} workers in {out['flights']['s']:.0f} s")

    pol = RL.load_policy(pol_path, device)
    t0 = time.time()
    RL.fast_success(pol, int(fast_eval_episodes or cfg.fast_eval_episodes), seed=12345, device=device)
    RL.fast_success_on(pol, S5.conditions(cfg, "A256"), device)
    out["fast_eval"] = {"s": round(time.time() - t0, 1)}
    return out


def project(cfg: S5.StudyConfigV5, m: dict, arms=("dagger_fnpeg", "bcppo_fnpeg", "ppo_hist")) -> dict:
    """Hours for every Tier 1 RL tuning configuration and training run from the measurements m."""
    H = 3600.0
    f = cfg.rl_tune_frac
    n_dev = len(S5.conditions(cfg, "DEV48"))
    flight_h = m["flights"]["s_per_flight"] * n_dev / H            # one configuration's DEV48 flights (all cores)
    fe_h = m["fast_eval"]["s"] / H
    chunk_h = m["expert_chunk"]["s"] / H
    tchunk_h = m["expert_tune_chunk"]["s"] / H
    fit_h = m["fit_chunk"]["s"] / H                                 # one chunk's states at dagger_epochs

    def ppo_h(rate, dec, hours):
        return min(hours, dec / max(rate, 1e-9) / H)

    rows = {}
    for arm in arms:
        c = _cfg(cfg, arm)
        rate = m["bcppo_fnpeg_ppo" if arm == "bcppo_fnpeg" else "ppo_hist"]["decisions_per_s"]
        n_iter, epochs_ratio = c.dagger_iters, c.dagger_epochs / max(m["fit_chunk"]["epochs"], 1)
        if arm == "ppo_hist":
            train = ppo_h(rate, c.budget_decisions, c.budget_hours)
            tune = ppo_h(rate, c.budget_decisions * f, c.budget_hours * f)
        elif arm == "bcppo_fnpeg":
            n_chunks = math.ceil(c.bc_episodes / c.dagger_rollouts)
            bc_fit = fit_h * n_chunks * c.bc_epochs / max(m["fit_chunk"]["epochs"], 1)
            train = n_chunks * chunk_h + bc_fit + ppo_h(rate, c.budget_decisions, c.budget_hours)
            t_rollouts = max(2, int(c.dagger_rollouts * f))
            t_chunks = math.ceil(max(2, int(c.bc_episodes * f)) / t_rollouts)
            tune = (t_chunks * tchunk_h + bc_fit * f + ppo_h(rate, c.budget_decisions * f, c.budget_hours * f))
        elif arm == "dagger_fnpeg":
            fits = sum(range(1, n_iter)) + n_iter                   # before iterations 1..n-1, then the final fit
            train = n_iter * chunk_h + fits * fit_h * epochs_ratio
            tune = n_iter * tchunk_h + fits * fit_h * epochs_ratio * f
        else:
            raise ValueError(arm)
        n_tune, n_seed = int(cfg.rl_tune_n.get(arm, 0)), len(cfg.rl_seeds)
        rows[arm] = {"train_h_per_run": round(train + fe_h, 3), "tune_h_per_config": round(tune + flight_h, 3),
                     "runs": n_seed, "configs": n_tune,
                     "total_h": round(n_seed * (train + fe_h) + n_tune * (tune + flight_h), 2),
                     "ppo_hits_time_cap": (c.budget_decisions / max(rate, 1e-9) / H) >= c.budget_hours
                     if arm != "dagger_fnpeg" else None}
    total_h = sum(r["total_h"] for r in rows.values())
    u = units_per_hour(m["hardware"].get("gpu") or "cpu")
    return {"arms": rows, "total_h": round(total_h, 1), "units_per_hour": u, "units": round(total_h * u, 1),
            "kaggle_weeks_at_30h": round(total_h / 30.0, 1)}


def summary_text(m: dict, p: dict) -> str:
    lines = [f"GPU {m['hardware'].get('gpu')}, {m['hardware'].get('cpu_count')} CPU cores",
             f"PPO: ppo_hist {m['ppo_hist']['decisions_per_s']:,.0f} decisions/s, bcppo_fnpeg "
             f"{m['bcppo_fnpeg_ppo']['decisions_per_s']:,.0f} decisions/s",
             f"FNPEG expert: {m['expert_chunk']['s']:.0f} s per {m['chunk']}-episode chunk, "
             f"{m['expert_tune_chunk']['s']:.0f} s per {m['tune_chunk']}-episode chunk; fit {m['fit_chunk']['s']:.0f} s",
             f"policy flights: {m['flights']['s_per_flight']:.1f} s each (wall, {m['flights']['workers']} workers); "
             f"fast eval {m['fast_eval']['s']:.0f} s per run", ""]
    for arm, r in p["arms"].items():
        cap = "" if r["ppo_hits_time_cap"] is None else (", PPO stops at the 3 h cap" if r["ppo_hits_time_cap"]
                                                         else ", PPO reaches its decision budget before 3 h")
        lines.append(f"{arm}: {r['runs']} runs x {r['train_h_per_run']:.2f} h + {r['configs']} tuning configs x "
                     f"{r['tune_h_per_config']:.2f} h = {r['total_h']:.1f} h{cap}")
    lines += ["", f"ALL TIER 1 RL ON THIS RUNTIME: {p['total_h']:.0f} hours = about {p['units']:.0f} Colab units "
                  f"(at {p['units_per_hour']}/h); on Kaggle's free 30 h/week: about {p['kaggle_weeks_at_30h']} weeks"]
    return "\n".join(lines)
