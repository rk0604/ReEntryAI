"""
jepa/mppi.py

Receding-horizon MPPI on a learned model (SkyJEPA Sec. IV-D), and the
closed-loop episode runner for the head-to-head against the classical
predictor-corrector.

Planner, once per macro step (5 s):
  1. sample S bank sequences = nominal + smooth noise. The noise lives on a few
     knots spread over the predicted time-to-go and is linearly interpolated,
     so a ~100-step horizon costs only n_knots random dimensions.
  2. roll each one through the world model + prober to the drogue altitude
     (all S samples in one batch; stops once every sample has crossed).
  3. cost = predicted great-circle miss at the drogue crossing (km)
          + smoothness (km per rad of command change, a proxy for RCS fuel)
          + g-load over the limit (read off the prober's specific forces)
          + flight-path angle steeper than gamma_limit_deg before drogue (the
            env terminates at |gamma| > 87 deg; the classical controller avoids
            it with a hard wings-level rule, here it is a soft state constraint
            on the prober's interpretable output)
          + a large penalty if the sample never reaches drogue altitude.
  4. softmax-weight the samples (temperature set from the cost spread, so the
     same settings work at 500 km and at 2 km), update the nominal, execute
     its first command, and shift.
No sign logic and no deadband: bank reversals emerge from the cost.

The classical baseline runs through the same runner (OrionEntryEnv baseline
mode). Its predictor keeps the nominal aero model unless classical_knows_truth
is set, which is what makes the domain-mismatch condition a fair question:
does history conditioning let the learned model adapt where a fixed nominal
model cannot?
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional

import numpy as np
import torch

from jepa.data import (G0, NOMINAL_DOMAIN, SIGMA_MAX_RAD, MACRO_DT_S, apply_domain,
                       physical_state)
from jepa.evaluate import crossing_point, haversine_km
from jepa.models import jepa_predict, load_checkpoint, predictive_predict


@dataclass
class MPPIConfig:
    """Planner settings. Costs are in km-equivalents so they share the miss's units."""
    n_samples: int = 256
    horizon: int = 110               # macro steps (550 s), an upper bound on time-to-go
    n_knots: int = 6
    noise_std_deg: float = 30.0
    n_iters: int = 1
    n_iters_first: int = 4           # cold start at entry interface
    temp_quantile: float = 0.1       # temperature = kappa * (J at this quantile - J_min)
    temp_kappa: float = 0.5
    w_smooth_km_per_rad: float = 2.0
    w_g_km: float = 50.0
    g_limit: float = 8.0
    gamma_limit_deg: float = 82.0    # classical reaches drogue at ~79 deg; the gate is 87 deg
    w_gamma_km_per_deg2: float = 2.0
    nocross_km: float = 3000.0
    n_sub: int = 2                   # Euler substeps per macro step in the prober
    seed: int = 0


class MPPI:
    """Model predictive path integral control over a bank-angle sequence."""

    def __init__(self, rollout_fn: Callable[[torch.Tensor], Dict[str, torch.Tensor]],
                 cfg: MPPIConfig, target_phi: float, target_lam: float, h_drogue: float):
        """rollout_fn maps commands [S, K] (rad) to {"y6": [S, K'+1, 6], "gload": [S, K']}."""
        self.rollout_fn, self.cfg = rollout_fn, cfg
        self.target = (float(target_phi), float(target_lam))
        self.h_drogue = float(h_drogue)
        self.gen = torch.Generator().manual_seed(cfg.seed)
        self.reset()

    def reset(self) -> None:
        """Forget the warm start."""
        self.u_nom: Optional[torch.Tensor] = None
        self.k_rem = self.cfg.horizon
        self.last: Dict[str, float] = {}

    def _noise_basis(self, K: int) -> torch.Tensor:
        """Hat-function interpolation matrix [n_knots, K] over the predicted time-to-go."""
        n = self.cfg.n_knots
        k_rem = float(min(max(self.k_rem, n), K))
        knots = torch.linspace(0.0, k_rem - 1.0, n)
        k = torch.arange(K, dtype=torch.float32).clamp_max(k_rem - 1.0)
        spacing = (k_rem - 1.0) / (n - 1)
        return (1.0 - (k[None, :] - knots[:, None]).abs() / spacing).clamp_min(0.0)

    def _cost(self, out: Dict[str, torch.Tensor], u: torch.Tensor, u_prev: torch.Tensor):
        """Per-sample cost, crossing step, and predicted miss."""
        c = self.cfg
        y6 = out["y6"]
        Kp = y6.shape[1] - 1
        phi, lam, crossed, kc = crossing_point(y6, self.h_drogue)
        miss = haversine_km(phi, lam, torch.full_like(phi, self.target[0]), torch.full_like(lam, self.target[1]))
        valid = (torch.arange(Kp)[None, :] < kc[:, None]).float()
        du = torch.diff(torch.cat([u_prev.expand(u.shape[0], 1), u[:, :Kp]], 1), dim=1).abs()
        J = miss + (~crossed).float() * c.nocross_km + c.w_smooth_km_per_rad * (du * valid).sum(1)
        if out.get("gload") is not None:
            over = torch.relu(out["gload"][:, :Kp] - c.g_limit)
            J = J + c.w_g_km * (over ** 2 * valid).sum(1)
        # States 1..Kp that are still above drogue altitude.
        above = (torch.arange(1, Kp + 1)[None, :] < kc[:, None]).float()
        steep = torch.relu(torch.rad2deg(y6[:, 1:, 4].abs()) - c.gamma_limit_deg)
        J = J + c.w_gamma_km_per_deg2 * (steep ** 2 * above).sum(1)
        return J, kc, miss

    @torch.no_grad()
    def plan(self, u_prev: torch.Tensor) -> float:
        """One receding-horizon decision; returns the bank command (rad) to hold."""
        c = self.cfg
        K, S = c.horizon, c.n_samples
        first = self.u_nom is None
        if first:
            self.u_nom = torch.zeros(K)
        M = self._noise_basis(K)
        sd = math.radians(c.noise_std_deg)
        for _ in range(c.n_iters_first if first else c.n_iters):
            eps = (torch.randn(S, c.n_knots, generator=self.gen) * sd) @ M
            eps[0] = 0.0                                    # keep the warm start as a candidate
            u = (self.u_nom + eps).clamp(-SIGMA_MAX_RAD, SIGMA_MAX_RAD)
            J, kc, miss = self._cost(self.rollout_fn(u), u, u_prev)
            J_min = J.min()
            temp = max(c.temp_kappa * float(torch.quantile(J, c.temp_quantile) - J_min), 1e-6)
            w = torch.softmax(-(J - J_min) / temp, 0)
            self.u_nom = (w[:, None] * u).sum(0)
        best = int(torch.argmax(w))
        cmd = float(self.u_nom[0])
        self.last = {"J_min": float(J_min), "miss_best_km": float(miss[best]),
                     "ess": float(1.0 / (w ** 2).sum()), "k_rem": int(kc[best])}
        self.u_nom = torch.cat([self.u_nom[1:], self.u_nom[-1:]])
        self.k_rem = max(2, int(kc[best]) - 1)
        return cmd


class LearnedModelController:
    """
    MPPI on the JEPA world model + PI prober ("jepa") or on the predictive
    baseline ("pred"). Call act() once per macro step with the current
    physical state; it keeps the state/command history the encoders need.
    """

    def __init__(self, kind: str, models: dict, cfg: MPPIConfig, target_phi: float,
                 target_lam: float, h_drogue: float, macro_dt: float = MACRO_DT_S):
        """models is the dict returned by load_checkpoint."""
        self.kind, self.m, self.cfg, self.dt = kind, models, cfg, float(macro_dt)
        self.stats = models["stats"]
        self.rf = models["wm"].rf if kind == "jepa" else models["pred"].rf
        self.h_drogue = float(h_drogue)
        self.mppi = MPPI(self._rollout, cfg, target_phi, target_lam, h_drogue)

    def reset(self, y0: np.ndarray) -> None:
        """Front-pad the history with the entry state and wings-level commands."""
        self.y_hist = torch.as_tensor(y0, dtype=torch.float32)[None].repeat(self.rf, 1)
        self.u_hist = torch.zeros(self.rf - 1)
        self.mppi.reset()

    def _rollout(self, u: torch.Tensor) -> Dict[str, torch.Tensor]:
        """Batched model rollout of candidate command sequences."""
        S = u.shape[0]
        if self.kind == "jepa":
            out = jepa_predict(self.m["wm"], self.m["pi"], self.stats.norm_feats(self.y_hist)[None],
                               self.y_hist[-1, :6].expand(S, 6), self.u_hist.expand(S, -1), u,
                               self.dt, self.cfg.n_sub, h_stop=self.h_drogue)
            out["gload"] = out["acc"].norm(dim=-1) / G0
        else:
            out = predictive_predict(self.m["pred"], self.stats, self.y_hist[None],
                                     self.u_hist.expand(S, -1), u, h_stop=self.h_drogue)
            y = out["y"][:, 1:]
            out["gload"] = torch.sqrt(y[..., 8].clamp_min(0) ** 2 + y[..., 9].clamp_min(0) ** 2) / G0
        return out

    def act(self, y_now: np.ndarray) -> float:
        """Append the current state, plan, append the chosen command, return it."""
        self.y_hist = torch.cat([self.y_hist[1:], torch.as_tensor(y_now, dtype=torch.float32)[None]])
        cmd = self.mppi.plan(self.u_hist[-1:])
        self.u_hist = torch.cat([self.u_hist[1:], torch.tensor([cmd])])
        return cmd


# ============================================================================
# Closed-loop episode runner (worker-safe)
# ============================================================================
_MODEL_CACHE: Dict[str, dict] = {}


def _models(path: str) -> dict:
    """Load a checkpoint once per process."""
    if path not in _MODEL_CACHE:
        _MODEL_CACHE[path] = load_checkpoint(path)
    return _MODEL_CACHE[path]


def run_episode(job: dict) -> dict:
    """
    Fly one entry from the given IC under the given truth domain.

    job keys:
        controller   "jepa", "pred", or "classical"
        ic           ic_override dict (frozen test-set row), or None for nominal
        domain       truth perturbation (see data.apply_domain)
        ckpt         checkpoint path (learned controllers)
        mppi         dict of MPPIConfig fields
        config_path  mission config (default configs/default.json)
        classical_knows_truth  give the classical predictor the true aero too
        record       keep a 5 s trajectory for plotting
        tag, case    free-form labels copied into the result
        threads      torch threads (default 1, right for pool workers)

    Returns a metric dict (miss, success, peak g, fuel, planning time, ...).
    """
    torch.set_num_threads(int(job.get("threads", 1)))
    import rl_env
    kind = job["controller"]
    classical = kind == "classical"
    env = rl_env.OrionEntryEnv(config_path=job.get("config_path", "configs/default.json"),
                               controller_mode="baseline" if classical else "policy",
                               dispersion=False, max_episode_steps=1400, fast_atmosphere=True)
    domain = dict(job.get("domain") or NOMINAL_DOMAIN)
    apply_domain(env.params, domain)
    if classical and job.get("classical_knows_truth"):
        apply_domain(env.guidance.params, domain)
    env.reset(seed=int(job.get("seed", 0)),
              options={"ic_override": job["ic"]} if job.get("ic") else None)
    n_sub = int(round(MACRO_DT_S / env.guidance_period_s))

    ctrl = None
    if not classical:
        cfg = MPPIConfig(**job.get("mppi", {}))
        ctrl = LearnedModelController(kind, _models(job["ckpt"]), cfg, env.target_phi_rad,
                                      env.target_lam_rad, env.drogue_alt_m)
        ctrl.reset(physical_state(env))

    t0 = time.time()
    plan_ms: List[float] = []
    cmds: List[float] = []
    traj: List[np.ndarray] = []
    sigma, k, done, info = 0.0, 0, False, {}
    while not done:
        if ctrl is not None and k % n_sub == 0:
            tp = time.perf_counter()
            sigma = ctrl.act(physical_state(env))
            plan_ms.append(1000.0 * (time.perf_counter() - tp))
        if job.get("record") and k % n_sub == 0:
            traj.append(np.append(physical_state(env), sigma))
        _, _, term, trunc, info = env.step(np.array([sigma / SIGMA_MAX_RAD], dtype=np.float32))
        if classical:
            sigma = math.radians(float((env.guidance.get_last_debug() or {}).get("chosen_sigma_cmd_deg", 0.0)))
        cmds.append(sigma)
        k += 1
        done = bool(term or trunc)
    if job.get("record"):
        traj.append(np.append(physical_state(env), sigma))

    c = np.asarray(cmds)
    s = np.sign(c[np.abs(c) > math.radians(5.0)])
    return {
        "tag": job.get("tag", ""), "case": job.get("case", -1), "controller": kind,
        "outcome": info.get("outcome"), "miss_km": float(info.get("miss_km", float("nan"))),
        "success": bool(info.get("success", False)), "peak_g": float(info.get("peak_g", float("nan"))),
        "peak_qdot_MWm2": float(info.get("peak_qdot_MWm2", float("nan"))),
        "rcs_fuel_kg": float(info.get("rcs_fuel_kg", float("nan"))),
        "t_s": float(info.get("t_s", float("nan"))), "wall_s": round(time.time() - t0, 2),
        "plan_ms_mean": float(np.mean(plan_ms)) if plan_ms else float("nan"),
        "plan_ms_max": float(np.max(plan_ms)) if plan_ms else float("nan"),
        "n_reversals": int((s[1:] != s[:-1]).sum()) if len(s) > 1 else 0,
        **{f"dom_{k_}": v for k_, v in domain.items()},
        "traj": np.asarray(traj, dtype=np.float32) if traj else None,
        "target": (float(env.target_phi_rad), float(env.target_lam_rad)),
    }


def summarize(results: List[dict]):
    """Per (tag, controller) table: success, miss percentiles, peak g, fuel, planning time."""
    import pandas as pd
    df = pd.DataFrame([{k: v for k, v in r.items() if k not in ("traj", "target")} for r in results])
    g = df.groupby(["tag", "controller"])
    return pd.DataFrame({
        "n": g.size(),
        "success_%": 100.0 * g["success"].mean(),
        "miss_med_km": g["miss_km"].median(),
        "miss_p90_km": g["miss_km"].quantile(0.9),
        "peak_g_mean": g["peak_g"].mean(),
        "rcs_fuel_kg": g["rcs_fuel_kg"].mean(),
        "reversals": g["n_reversals"].mean(),
        "plan_ms": g["plan_ms_mean"].mean(),
        "wall_s": g["wall_s"].mean(),
    }).round(2)
