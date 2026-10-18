"""
jepa/mppi_v2.py

MPPI on the v2 ensemble, and a lockstep closed-loop runner for the GPU.

The planner is jepa.mppi.MPPI with two changes, both in EnsembleMPPI:
  * the cost of a sample is the mean of v1's cost over the ensemble members
    (each member flies its own trajectory through its own prober), plus
    w_unc_km x the members' spread of predicted landing points (km). With
    one member and w_unc_km = 0 it is exactly v1's planner.
  * plan() is split into propose() / update(), so many flights can share one
    batched model rollout per decision.

run_lockstep() flies many episodes in lockstep: every 5 s, every deciding
episode proposes its samples, all of them go through the ensemble in one
batched pass on the GPU, and each episode updates its own MPPI state. Env
construction, truth perturbation, the 1 Hz loop, and the result dict mirror
jepa.mppi.run_episode, so jepa.mppi.summarize() tables v1, v2, and classical
results together.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Sequence

import numpy as np
import torch

from jepa.data import G0, MACRO_DT_S, NOMINAL_DOMAIN, SIGMA_MAX_RAD, apply_domain, physical_state
from jepa.evaluate import crossing_point, haversine_km
from jepa.models_v2 import amp_dtype, ensemble_rollout, load_checkpoint_v2
from jepa.mppi import MPPI, MPPIConfig


@dataclass
class MPPIConfigV2(MPPIConfig):
    """v1's MPPIConfig (same defaults) plus the ensemble-uncertainty weight."""
    w_unc_km: float = 1.0              # km of cost per km of landing-point spread across members
    trim_horizon: bool = True          # after the first decision, roll out ~1.5x the predicted time-to-go


class EnsembleMPPI(MPPI):
    """MPPI whose rollouts carry a member dimension: y6 [S, M, K'+1, 6]."""

    def __init__(self, cfg: MPPIConfigV2, target_phi: float, target_lam: float, h_drogue: float):
        super().__init__(rollout_fn=None, cfg=cfg, target_phi=target_phi, target_lam=target_lam,
                         h_drogue=h_drogue)

    def _cost(self, out: Dict[str, torch.Tensor], u: torch.Tensor, u_prev: torch.Tensor):
        """Mean member cost + w_unc * landing spread. Returns (J [S], k_cross [S], miss [S], spread [S])."""
        y6 = out["y6"]
        S, Mm = y6.shape[:2]
        flat = {"y6": y6.flatten(0, 1),
                "gload": out["gload"].flatten(0, 1) if out.get("gload") is not None else None}
        J, kc, miss = super()._cost(flat, u.repeat_interleave(Mm, 0), u_prev)
        phi, lam, _, _ = crossing_point(flat["y6"], self.h_drogue)
        phi, lam = phi.view(S, Mm), lam.view(S, Mm)
        d = haversine_km(phi, lam, phi.mean(1, keepdim=True), lam.mean(1, keepdim=True))
        spread = torch.sqrt((d ** 2).mean(1))
        J = J.view(S, Mm).mean(1) + self.cfg.w_unc_km * spread
        kc = kc.view(S, Mm).float().median(1).values.long()
        return J, kc, miss.view(S, Mm).mean(1), spread

    # ---- plan() split in two (same arithmetic and random stream as MPPI.plan) ----
    def begin(self) -> int:
        """Start a decision; returns how many MPPI iterations it takes."""
        first = self.u_nom is None
        if first:
            self.u_nom = torch.zeros(self.cfg.horizon)
        self._basis = self._noise_basis(self.cfg.horizon)
        return self.cfg.n_iters_first if first else self.cfg.n_iters

    def propose(self) -> torch.Tensor:
        """Sample S command sequences [S, K] around the warm start."""
        c = self.cfg
        eps = (torch.randn(c.n_samples, c.n_knots, generator=self.gen) * math.radians(c.noise_std_deg)) @ self._basis
        eps[0] = 0.0
        return (self.u_nom + eps).clamp(-SIGMA_MAX_RAD, SIGMA_MAX_RAD)

    def update(self, out: Dict[str, torch.Tensor], u: torch.Tensor, u_prev: torch.Tensor,
               final: bool) -> Optional[float]:
        """Weight the samples; on the decision's last iteration return the command and shift."""
        c = self.cfg
        J, kc, miss, spread = self._cost(out, u, u_prev)
        J_min = J.min()
        temp = max(c.temp_kappa * float(torch.quantile(J, c.temp_quantile) - J_min), 1e-6)
        w = torch.softmax(-(J - J_min) / temp, 0)
        self.u_nom = (w[:, None] * u).sum(0)
        if not final:
            return None
        best = int(torch.argmax(w))
        cmd = float(self.u_nom[0])
        self.last = {"J_min": float(J_min), "miss_best_km": float(miss[best]), "spread_km": float(spread[best]),
                     "ess": float(1.0 / (w ** 2).sum()), "k_rem": int(kc[best])}
        self.u_nom = torch.cat([self.u_nom[1:], self.u_nom[-1:]])
        self.k_rem = max(2, int(kc[best]) - 1)
        return cmd


class _Flight:
    """One closed-loop flight in the lockstep runner (env, histories, planner, logs)."""

    def __init__(self, job: dict, rf: int, cfg: MPPIConfigV2):
        import rl_env
        self.job = job
        self.env = rl_env.OrionEntryEnv(config_path=job.get("config_path", "configs/default.json"),
                                        controller_mode="policy", dispersion=False,
                                        max_episode_steps=1400, fast_atmosphere=True)
        self.domain = dict(job.get("domain") or NOMINAL_DOMAIN)
        apply_domain(self.env.params, self.domain)
        self.env.reset(seed=int(job.get("seed", 0)),
                       options={"ic_override": job["ic"]} if job.get("ic") else None)
        self.n_sub = int(round(MACRO_DT_S / self.env.guidance_period_s))
        self.mppi = EnsembleMPPI(cfg, self.env.target_phi_rad, self.env.target_lam_rad, self.env.drogue_alt_m)
        y0 = torch.as_tensor(physical_state(self.env), dtype=torch.float32)
        self.y_hist = y0[None].repeat(rf, 1)
        self.u_hist = torch.zeros(rf - 1)
        self.sigma, self.k, self.done, self.info = 0.0, 0, False, {}
        self.cmds: List[float] = []
        self.plan_ms: List[float] = []
        self.spread: List[float] = []
        self.traj: List[np.ndarray] = []
        self.t0 = time.time()

    def observe(self) -> None:
        """Append the current physical state to the encoder history."""
        self._y_now = physical_state(self.env)
        self.y_hist = torch.cat([self.y_hist[1:], torch.as_tensor(self._y_now, dtype=torch.float32)[None]])

    def commit(self, cmd: float) -> None:
        """Hold cmd for the next macro step (and log the state it was chosen at, like run_episode)."""
        self.sigma = cmd
        self.u_hist = torch.cat([self.u_hist[1:], torch.tensor([cmd])])
        self.spread.append(self.mppi.last.get("spread_km", float("nan")))
        if self.job.get("record"):
            self.traj.append(np.append(self._y_now, cmd))

    def env_step(self) -> None:
        """One 1 Hz guidance step of the real closed loop."""
        _, _, term, trunc, self.info = self.env.step(np.array([self.sigma / SIGMA_MAX_RAD], dtype=np.float32))
        self.cmds.append(self.sigma)
        self.k += 1
        self.done = bool(term or trunc)

    def result(self) -> dict:
        """Metric dict in jepa.mppi.run_episode's format (plus the mean planner spread)."""
        if self.job.get("record"):
            self.traj.append(np.append(physical_state(self.env), self.sigma))
        info, c = self.info, np.asarray(self.cmds)
        s = np.sign(c[np.abs(c) > math.radians(5.0)])
        return {
            "tag": self.job.get("tag", ""), "case": self.job.get("case", -1),
            "controller": self.job.get("label", "v2"),
            "outcome": info.get("outcome"), "miss_km": float(info.get("miss_km", float("nan"))),
            "success": bool(info.get("success", False)), "peak_g": float(info.get("peak_g", float("nan"))),
            "peak_qdot_MWm2": float(info.get("peak_qdot_MWm2", float("nan"))),
            "rcs_fuel_kg": float(info.get("rcs_fuel_kg", float("nan"))),
            "t_s": float(info.get("t_s", float("nan"))), "wall_s": round(time.time() - self.t0, 2),
            "plan_ms_mean": float(np.mean(self.plan_ms)) if self.plan_ms else float("nan"),
            "plan_ms_max": float(np.max(self.plan_ms)) if self.plan_ms else float("nan"),
            "n_reversals": int((s[1:] != s[:-1]).sum()) if len(s) > 1 else 0,
            "spread_km_mean": float(np.nanmean(self.spread)) if self.spread else float("nan"),
            **{f"dom_{k_}": v for k_, v in self.domain.items()},
            "traj": np.asarray(self.traj, dtype=np.float32) if self.traj else None,
            "target": (float(self.env.target_phi_rad), float(self.env.target_lam_rad)),
        }


@torch.no_grad()
def run_lockstep(jobs: Sequence[dict], models: dict, mppi: Optional[dict] = None, n_sub: int = 2,
                 chunk: int = 4096, verbose: bool = True) -> List[dict]:
    """
    Fly every job (same keys as jepa.mppi.run_episode: ic, domain, tag, case,
    seed, record, config_path, plus optional label and per-job mppi overrides)
    with MPPI on the v2 ensemble, all episodes in lockstep.

    models: load_checkpoint_v2(...) dict on any device. mppi: MPPIConfigV2
    fields. plan_ms per episode is the batched round time divided by the
    number of episodes deciding in it (amortized).
    """
    wm, pi, stats = models["wm"], models["pi"], models["stats"]
    dev = next(wm.parameters()).device
    dt_amp = amp_dtype(dev)
    flights = [_Flight(j, wm.rf, MPPIConfigV2(**{**(mppi or {}), **j.get("mppi", {})})) for j in jobs]
    t0, n_round = time.time(), 0
    while True:
        active = [f for f in flights if not f.done]
        if not active:
            break
        deciding = [f for f in active if f.k % f.n_sub == 0]
        if deciding:
            tr = time.perf_counter()
            for f in deciding:
                f.observe()
            feats = torch.stack([stats.norm_feats(f.y_hist) for f in deciding]).to(dev)
            with torch.autocast(dev.type, dtype=dt_amp, enabled=dt_amp is not None):
                s_ctx = wm.encode_states(feats)[:, :, -1].float()                 # [M, E, D]
            iters = [f.mppi.begin() for f in deciding]
            for it in range(max(iters)):
                group = [i for i, f in enumerate(deciding) if it < iters[i]]
                us = [deciding[i].mppi.propose() for i in group]
                S = us[0].shape[0]
                K = flights[0].mppi.cfg.horizon
                if all(deciding[i].mppi.cfg.trim_horizon and deciding[i].mppi.u_nom is not None
                       and deciding[i].mppi.k_rem < K for i in group):
                    K = min(K, int(1.5 * max(deciding[i].mppi.k_rem for i in group)) + 12)
                u = torch.cat(us)[:, :K]
                s0 = s_ctx[:, group].repeat_interleave(S, 1)
                y0 = torch.stack([deciding[i].y_hist[-1, :6] for i in group]).repeat_interleave(S, 0)
                uh = torch.stack([deciding[i].u_hist for i in group]).repeat_interleave(S, 0)
                with torch.autocast(dev.type, dtype=dt_amp, enabled=dt_amp is not None):
                    out = ensemble_rollout(wm, pi, s0, y0.to(dev), uh.to(dev), u.to(dev), MACRO_DT_S, n_sub,
                                           h_stop=float(deciding[0].env.drogue_alt_m), chunk=chunk)
                y6 = out["y6"].cpu()
                gl = (out["acc"].norm(dim=-1) / G0).cpu()
                for j, i in enumerate(group):
                    f = deciding[i]
                    sl = slice(j * S, (j + 1) * S)
                    cmd = f.mppi.update({"y6": y6[sl], "gload": gl[sl]}, us[j], f.u_hist[-1:],
                                        final=it == iters[i] - 1)
                    if cmd is not None:
                        f.commit(cmd)
            ms = 1000.0 * (time.perf_counter() - tr) / len(deciding)
            for f in deciding:
                f.plan_ms.append(ms)
            n_round += 1
            if verbose and n_round % 20 == 0:
                print(f"  [lockstep] round {n_round}: {len(active)}/{len(flights)} flying, "
                      f"{ms:.0f} ms/decision, {time.time() - t0:.0f}s", flush=True)
        for f in active:
            f.env_step()
    return [f.result() for f in flights]


def run_lockstep_job(job: dict) -> List[dict]:
    """
    Pool worker: {"jobs": [...], "ckpt": path, "device": "cuda"|"cpu", "mppi": {...},
    "threads": 1}. Loads the v2 checkpoint once and flies its share of episodes.
    """
    torch.set_num_threads(int(job.get("threads", 1)))
    models = load_checkpoint_v2(job["ckpt"], job.get("device", "cpu"))
    return run_lockstep(job["jobs"], models, job.get("mppi"), n_sub=int(job.get("n_sub", 2)),
                        chunk=int(job.get("chunk", 4096)), verbose=bool(job.get("verbose", False)))


def shard(jobs: Sequence[dict], n: int) -> List[List[dict]]:
    """Split jobs into n round-robin shards (similar flight lengths per shard)."""
    return [list(jobs[i::n]) for i in range(n) if jobs[i::n]]
