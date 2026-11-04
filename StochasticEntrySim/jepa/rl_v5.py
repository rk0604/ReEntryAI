"""
jepa/rl_v5.py

The v5 RL baselines (jepa/v5_plan.md section 4.4), in-house torch, resumable (section 9.1).

    Actor          "gauss": Gaussian over the action a in [-1, 1] (sigma = 90 deg x a), state-independent log std
                   (PPO); "magsign": |a| ~ Gaussian on a sigmoid mean and sign ~ Bernoulli (DAgger, BC -> PPO:
                   reversals are a class, never averaged away); "tanh": tanh-squashed Gaussian (SAC)
    ObsNorm        running mean and variance of the observation (saved with every checkpoint)
    ppo_train      CleanRL-style PPO (GAE 0.95, gamma 0.99, clip 0.2, 4 epochs, value and entropy terms),
                   optional KL penalty to a frozen reference actor (the clone of bcppo_fnpeg) and a critic-only
                   warm-up on frozen-actor rollouts; budget in decisions and in wall hours
    sac_train      SAC with twin critics and automatic temperature; the replay buffer is written as append-only
                   shards with each checkpoint and reloaded up to its capacity
    ppo_full_train ppo_matched: PPO in the full simulator (rl_env, 5 s decisions, the same observation and reward),
                   6,000 episodes, each update a set of complete episodes flown in a process pool
    dagger_train   DAgger with fnpeg as the expert, queried on the learner's own states (the expert's filters
                   run on every flight whoever acts), beta_i = 0.5^i, rollouts per iteration, the aggregated
                   data in data/iter_<i>.npz, supervised fit of the magsign actor (MSE on |a|, BCE on the sign)
    bc_dataset     fnpeg episodes in the fast env (states, expert actions) for the clone
    train_state    every save_every_s seconds and at the end: actors, critics, optimizers, normalizer, counters,
                   the fast env's state and RNG, every RNG state, the learning curve; .prev kept; a run's
                   policy.pt (actor + normalizer) marks it done
    legacy         list MyDrive/ReEntryAI/{models,checkpoints,data}; SB3 zips evaluated with their own
                   settings (1 Hz, rl_env's 13 features, deterministic)
"""

from __future__ import annotations

import io
import json
import math
import os
import time
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as Fnn

from jepa import rlenv_fast_v5 as RF
from jepa.models_v5 import fail_point, load_state, rng_state, save_state, set_rng_state

RL_ARMS = ("dagger_fnpeg", "bcppo_fnpeg", "ppo_hist", "sac_hist", "ppo_matched")
LOG_STD_MIN, LOG_STD_MAX = -5.0, 1.0


# ============================================================================
# Networks
# ============================================================================
def mlp(d_in: int, hidden: Sequence[int], d_out: int, act=nn.Tanh, out_std: float = 0.01) -> nn.Sequential:
    layers: List[nn.Module] = []
    d = d_in
    for h in hidden:
        lin = nn.Linear(d, h)
        nn.init.orthogonal_(lin.weight, math.sqrt(2))
        nn.init.zeros_(lin.bias)
        layers += [lin, act()]
        d = h
    out = nn.Linear(d, d_out)
    nn.init.orthogonal_(out.weight, out_std)
    nn.init.zeros_(out.bias)
    layers.append(out)
    return nn.Sequential(*layers)


class ObsNorm(nn.Module):
    def __init__(self, dim: int, clip: float = 10.0):
        super().__init__()
        self.register_buffer("mean", torch.zeros(dim, dtype=torch.float64))
        self.register_buffer("var", torch.ones(dim, dtype=torch.float64))
        self.register_buffer("count", torch.tensor(1e-4, dtype=torch.float64))
        self.clip = clip

    @torch.no_grad()
    def update(self, x: torch.Tensor) -> None:
        x = x.double().reshape(-1, x.shape[-1])
        b_mean, b_var, n = x.mean(0), x.var(0, unbiased=False), x.shape[0]
        delta = b_mean - self.mean
        tot = self.count + n
        self.mean += delta * n / tot
        self.var = (self.var * self.count + b_var * n + delta ** 2 * self.count * n / tot) / tot
        self.count = tot

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return ((x.double() - self.mean) / torch.sqrt(self.var + 1e-8)).clamp(-self.clip, self.clip).float()


class Actor(nn.Module):
    """kind: gauss | magsign | tanh. act(obs) -> (env action in [-1, 1], raw sample, log prob)."""

    def __init__(self, obs_dim: int = RF.OBS_DIM, hidden: Sequence[int] = (256, 256), kind: str = "gauss"):
        super().__init__()
        self.kind, self.hparams = kind, {"obs_dim": obs_dim, "hidden": list(hidden), "kind": kind}
        self.norm = ObsNorm(obs_dim)
        n_out = {"gauss": 1, "magsign": 2, "tanh": 2}[kind]
        self.net = mlp(obs_dim, hidden, n_out, act=nn.Tanh if kind != "tanh" else nn.ReLU)
        self.log_std = nn.Parameter(torch.full((1,), -0.5)) if kind != "tanh" else None

    def _heads(self, obs):
        o = self.net(self.norm(obs))
        return o

    def dist_params(self, obs):
        o = self._heads(obs)
        if self.kind == "gauss":
            return o[:, 0], self.log_std.clamp(LOG_STD_MIN, LOG_STD_MAX).expand_as(o[:, 0]), None
        if self.kind == "magsign":
            return torch.sigmoid(o[:, 0]), self.log_std.clamp(LOG_STD_MIN, LOG_STD_MAX).expand_as(o[:, 0]), o[:, 1]
        return o[:, 0], o[:, 1].clamp(LOG_STD_MIN, LOG_STD_MAX), None

    def deterministic(self, obs) -> torch.Tensor:
        mu, _, logit = self.dist_params(obs)
        if self.kind == "gauss":
            return mu.clamp(-1, 1)
        if self.kind == "magsign":
            return torch.where(logit >= 0, 1.0, -1.0) * mu.clamp(0, 1)
        return torch.tanh(mu)

    def sample(self, obs, gen: Optional[torch.Generator] = None):
        """(env action, raw [B, 2] = (value, sign or 0), log prob)."""
        mu, ls, logit = self.dist_params(obs)
        eps = torch.randn(mu.shape, generator=gen, device=mu.device)
        x = mu + ls.exp() * eps
        if self.kind == "gauss":
            raw = torch.stack([x, torch.zeros_like(x)], -1)
            return x.clamp(-1, 1), raw, self.log_prob(obs, raw)
        if self.kind == "magsign":
            u = torch.rand(mu.shape, generator=gen, device=mu.device)
            s = (u < torch.sigmoid(logit)).float() * 2 - 1
            raw = torch.stack([x, s], -1)
            return s * x.clamp(0, 1), raw, self.log_prob(obs, raw)
        a = torch.tanh(x)
        lp = (-0.5 * eps ** 2 - ls - 0.5 * math.log(2 * math.pi)) - torch.log(1 - a * a + 1e-6)
        return a, torch.stack([x, torch.zeros_like(x)], -1), lp

    def log_prob(self, obs, raw):
        mu, ls, logit = self.dist_params(obs)
        x = raw[:, 0]
        lp = -0.5 * ((x - mu) / ls.exp()) ** 2 - ls - 0.5 * math.log(2 * math.pi)
        if self.kind == "magsign":
            lp = lp - Fnn.binary_cross_entropy_with_logits(logit, (raw[:, 1] > 0).float(), reduction="none")
        return lp

    def entropy(self, obs):
        mu, ls, logit = self.dist_params(obs)
        e = ls + 0.5 * math.log(2 * math.pi * math.e)
        if self.kind == "magsign":
            p = torch.sigmoid(logit)
            e = e - (p * torch.log(p + 1e-8) + (1 - p) * torch.log(1 - p + 1e-8))
        return e

    def kl_to(self, ref: "Actor", obs) -> torch.Tensor:
        """KL(self || ref), per row (same kind)."""
        mu, ls, lg = self.dist_params(obs)
        with torch.no_grad():
            mu0, ls0, lg0 = ref.dist_params(obs)
        kl = ls0 - ls + (torch.exp(2 * ls) + (mu - mu0) ** 2) / (2 * torch.exp(2 * ls0)) - 0.5
        if self.kind == "magsign":
            p, p0 = torch.sigmoid(lg), torch.sigmoid(lg0)
            kl = kl + p * torch.log((p + 1e-8) / (p0 + 1e-8)) + (1 - p) * torch.log((1 - p + 1e-8) / (1 - p0 + 1e-8))
        return kl


def critic(obs_dim: int = RF.OBS_DIM, hidden=(256, 256), act_dim: int = 0) -> nn.Module:
    return mlp(obs_dim + act_dim, hidden, 1, out_std=1.0)


def save_policy(L, path, actor: Actor, extra: dict) -> None:
    buf = io.BytesIO()
    torch.save({"hparams": actor.hparams, "state": {k: v.cpu() for k, v in actor.state_dict().items()},
                "extra": extra}, buf)
    L.atomic_bytes(Path(path), buf.getvalue())


def load_policy(path, device="cpu") -> Actor:
    ck = torch.load(path, map_location="cpu", weights_only=False)
    a = Actor(**ck["hparams"])
    a.load_state_dict(ck["state"])
    return a.to(device).eval()


# ============================================================================
# Configuration
# ============================================================================
@dataclass
class RLCfg:
    n_envs: int = 4096
    hidden: Tuple[int, ...] = (256, 256)
    lr: float = 3e-4
    gamma: float = 0.99
    gae_lambda: float = 0.95
    n_steps: int = 64                       # decisions per env per PPO rollout
    epochs: int = 4
    minibatch: int = 8192
    clip: float = 0.2
    vf_coef: float = 0.5
    ent_coef: float = 0.0
    max_grad_norm: float = 0.5
    kl_coef: float = 0.0                    # bcppo: KL penalty to the clone
    critic_warmup_updates: int = 0          # bcppo: critic-only updates on frozen-actor rollouts
    budget_decisions: int = 200_000_000
    budget_hours: float = 3.0
    # SAC
    buffer: int = 2_000_000
    batch: int = 1024
    tau: float = 0.005
    learning_starts: int = 20_000
    updates_per_step: float = 0.25          # gradient steps per environment decision (per env)
    # DAgger / BC
    dagger_iters: int = 10
    dagger_rollouts: int = 200
    dagger_epochs: int = 30
    bc_episodes: int = 2000
    bc_epochs: int = 40
    sup_batch: int = 4096
    expert_beta_1hz: float = 0.9            # fnpeg's filter at 1 Hz; ^5 at the fast env's 5 s decisions
    # bookkeeping
    save_every_s: float = 300.0
    save_every_updates: int = 0             # > 0: also every n updates (the smoke test's deterministic resume)
    curve_points: int = 5                   # policy snapshots for the full-simulator learning curve
    env_pool: int = 2048                    # fast env: episodes drawn (footprints flown) per refill
    # ppo_matched (PPO in the full simulator, the world models' data budget)
    full_episodes: int = 6000
    full_per_update: int = 50               # complete episodes per PPO update (no episode spans a checkpoint)
    full_workers: int = 0                   # simulator processes (0: every core)
    seed: int = 0


def make_env(cfg: RLCfg, seed: int, device, n_envs: Optional[int] = None, fast: Optional[dict] = None):
    fc = RF.FastEnvCfg(n_envs=int(n_envs or cfg.n_envs), gamma_rl=cfg.gamma, seed=int(seed),
                       **{"pool": int(cfg.env_pool), **(fast or {})})
    return RF.FastEntryEnv(fc, device)


class Curve:
    """Learning curve rows (written to learning_curve.csv with every checkpoint)."""

    def __init__(self, rows=None):
        if isinstance(rows, dict):              # a saved state: rows and the running sums of the open window
            self.rows = list(rows["rows"])
            self.ret, self.succ, self.miss, self.n = rows["acc"]
        else:
            self.rows = list(rows or [])
            self.ret, self.succ, self.miss, self.n = 0.0, 0, 0.0, 0

    def state(self) -> dict:
        return {"rows": list(self.rows), "acc": (self.ret, self.succ, self.miss, self.n)}

    def add(self, info: dict) -> None:
        d = info["done"]
        if bool(d.any()):
            self.ret += float(info["ep_return"][d].sum())
            self.succ += int(info["success_5km"][d].sum())
            self.miss += float(info["miss_km"][d].clamp_max(5000).sum())
            self.n += int(d.sum())

    def flush(self, step: int, update: int, sec: float, **extra) -> None:
        n = max(self.n, 1)
        self.rows.append({"decisions": int(step), "update": int(update), "sec": round(sec, 1), "episodes": self.n,
                          "ep_return": self.ret / n, "success_5km_fast": self.succ / n, "miss_km_fast": self.miss / n,
                          **extra})
        self.ret, self.succ, self.miss, self.n = 0.0, 0, 0.0, 0


def _write_curve(L, d: Path, curve: Curve, title: str) -> None:
    import pandas as pd
    cols = ["decisions", "update", "sec", "episodes", "ep_return", "success_5km_fast", "miss_km_fast"]
    df = pd.DataFrame(curve.rows) if curve.rows else pd.DataFrame(columns=cols)
    L.write_csv(df, d / "learning_curve.csv")
    if not curve.rows:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(8, 3))
    ax[0].plot(df.decisions, df.ep_return)
    ax[0].set_xlabel("decisions")
    ax[0].set_ylabel("episode return (fast env)")
    ax[1].plot(df.decisions, df.success_5km_fast)
    ax[1].set_xlabel("decisions")
    ax[1].set_ylabel("success@5km (fast env)")
    fig.suptitle(title)
    for ext in (".png", ".pdf"):
        L.save_fig(fig, d / f"F_curves{ext}", dpi=200)
    plt.close(fig)


# ============================================================================
# PPO (ppo_hist, bcppo_fnpeg's RL stage)
# ============================================================================
def ppo_train(L, d: Path, cfg: RLCfg, device, label: str, actor: Optional[Actor] = None,
              ref: Optional[Actor] = None, log=print, stage: str = "rl_train", fast: Optional[dict] = None) -> dict:
    """PPO from `actor` (None: a fresh gauss actor). Resumable from d/train_state.pt; writes d/policy.pt."""
    dev = torch.device(device)
    torch.manual_seed(cfg.seed)
    actor = (actor or Actor(hidden=cfg.hidden, kind="gauss")).to(dev)
    vf = critic(hidden=cfg.hidden).to(dev)
    opt = torch.optim.Adam(list(actor.parameters()) + list(vf.parameters()), lr=cfg.lr, eps=1e-5)
    env = make_env(cfg, cfg.seed, dev, fast=fast)
    gen = torch.Generator(device=dev).manual_seed(cfg.seed * 7 + 3)
    curve, step, update, elapsed = Curve(), 0, 0, 0.0
    st = load_state(d / "train_state.pt", log)
    if st is not None:
        actor.load_state_dict(st["actor"])
        vf.load_state_dict(st["critic"])
        opt.load_state_dict(st["opt"])
        env.set_state(st["env"])
        step, update, elapsed, curve = st["step"], st["update"], st["elapsed"], Curve(st["curve"])
        set_rng_state(st["rng"], {"gen": gen})
        log(f"  [{label}] resumed at update {update} ({step:,} decisions)")
    obs = env.obs()
    N, T = env.cfg.n_envs, cfg.n_steps
    t_save, t0 = time.time(), time.time()
    snaps = _snapshot_marks(cfg)

    def save():
        save_state(L, {"actor": actor.state_dict(), "critic": vf.state_dict(), "opt": opt.state_dict(),
                       "env": env.get_state(), "step": step, "update": update, "elapsed": elapsed,
                       "curve": curve.state(), "rng": rng_state({"gen": gen}), "cfg": asdict(cfg)}, d / "train_state.pt")
        _write_curve(L, d, curve, label)

    while step < cfg.budget_decisions and elapsed + (time.time() - t0) < cfg.budget_hours * 3600:
        warm = update < cfg.critic_warmup_updates
        buf_o = torch.zeros(T, N, RF.OBS_DIM, device=dev)
        buf_raw = torch.zeros(T, N, 2, device=dev)
        buf_lp, buf_r, buf_d, buf_v = (torch.zeros(T, N, device=dev) for _ in range(4))
        with torch.no_grad():
            for k in range(T):
                actor.norm.update(obs) if not warm and ref is None else None
                a, raw, lp = actor.sample(obs, gen)
                buf_o[k], buf_raw[k], buf_lp[k] = obs, raw, lp
                buf_v[k] = vf(actor.norm(obs)).squeeze(-1)
                obs, r, done, info = env.step(a)
                obs = obs.to(dev)
                buf_r[k], buf_d[k] = r.to(dev), done.to(dev).float()
                curve.add(info)
            nv = vf(actor.norm(obs)).squeeze(-1)
            adv = torch.zeros_like(buf_r)
            last = torch.zeros(N, device=dev)
            for k in reversed(range(T)):
                nxt = nv if k == T - 1 else buf_v[k + 1]
                nonterm = 1.0 - buf_d[k]
                delta = buf_r[k] + cfg.gamma * nxt * nonterm - buf_v[k]
                last = delta + cfg.gamma * cfg.gae_lambda * nonterm * last
                adv[k] = last
            ret = adv + buf_v
        _ppo_update(actor, vf, opt, cfg, gen, ref, warm, buf_o.flatten(0, 1), buf_raw.flatten(0, 1),
                    buf_lp.flatten(0, 1), adv.flatten(), ret.flatten())
        step += N * T
        update += 1
        curve.flush(step, update, elapsed + time.time() - t0, warmup=int(warm))
        if update % 10 == 0 or update == 1:
            r = curve.rows[-1]
            log(f"  [{label}] update {update} {step:,} decisions: return {r['ep_return']:.3f} "
                f"success@5km {r['success_5km_fast']:.3f} ({(elapsed + time.time() - t0) / 60:.1f} min)")
        _maybe_snapshot(L, d, actor, step, cfg, snaps)
        due = time.time() - t_save > cfg.save_every_s or (cfg.save_every_updates and update % cfg.save_every_updates == 0)
        if due:
            elapsed += time.time() - t0
            t0 = time.time()
            save()
            t_save = time.time()
            fail_point(stage, update)
    elapsed += time.time() - t0
    t0 = time.time()
    save()
    save_policy(L, d / "policy.pt", actor, {"label": label, "decisions": step, "updates": update,
                                            "train_s": round(elapsed, 1), "cfg": asdict(cfg)})
    return {"decisions": step, "updates": update, "train_s": round(elapsed, 1)}


def _ppo_update(actor: Actor, vf: nn.Module, opt, cfg: RLCfg, gen, ref: Optional[Actor], warm: bool,
                b_o, b_raw, b_lp, b_adv, b_ret) -> None:
    """cfg.epochs passes of clipped-ratio PPO (critic only during a warm-up) over one batch of decisions."""
    n = len(b_o)
    mb = min(cfg.minibatch, n)
    with torch.no_grad():                   # the normalizer moved during the rollout: the ratio starts at 1
        b_lp = actor.log_prob(b_o, b_raw)
    for _ep in range(cfg.epochs):
        perm = torch.randperm(n, generator=gen, device=b_o.device)
        for i in range(0, n, mb):
            idx = perm[i:i + mb]
            v = vf(actor.norm(b_o[idx])).squeeze(-1)
            v_loss = 0.5 * ((v - b_ret[idx]) ** 2).mean()
            if warm:
                loss = cfg.vf_coef * v_loss
            else:
                lp = actor.log_prob(b_o[idx], b_raw[idx])
                ratio = torch.exp(lp - b_lp[idx])
                a_n = (b_adv[idx] - b_adv[idx].mean()) / (b_adv[idx].std() + 1e-8)
                pg = torch.max(-a_n * ratio, -a_n * ratio.clamp(1 - cfg.clip, 1 + cfg.clip)).mean()
                loss = pg + cfg.vf_coef * v_loss - cfg.ent_coef * actor.entropy(b_o[idx]).mean()
                if ref is not None and cfg.kl_coef > 0:
                    loss = loss + cfg.kl_coef * actor.kl_to(ref, b_o[idx]).mean()
            opt.zero_grad(set_to_none=True)
            loss.backward()
            nn.utils.clip_grad_norm_(list(actor.parameters()) + list(vf.parameters()), cfg.max_grad_norm)
            opt.step()


def _snapshot_marks(cfg: RLCfg) -> List[int]:
    if cfg.curve_points <= 0:
        return []
    return [int(cfg.budget_decisions * (i + 1) / cfg.curve_points) for i in range(cfg.curve_points)]


def _maybe_snapshot(L, d: Path, actor: Actor, step: int, cfg: RLCfg, marks: List[int]) -> None:
    for i, m in enumerate(marks):
        p = d / "snapshots" / f"policy_{i + 1}of{len(marks)}.pt"
        if step >= m and not p.exists():
            save_policy(L, p, actor, {"decisions": step, "point": i + 1})


# ============================================================================
# SAC (sac_hist)
# ============================================================================
def sac_train(L, d: Path, cfg: RLCfg, device, label: str, log=print, stage: str = "rl_train") -> dict:
    dev = torch.device(device)
    torch.manual_seed(cfg.seed)
    actor = Actor(hidden=cfg.hidden, kind="tanh").to(dev)
    q1, q2 = critic(hidden=cfg.hidden, act_dim=1).to(dev), critic(hidden=cfg.hidden, act_dim=1).to(dev)
    q1t, q2t = critic(hidden=cfg.hidden, act_dim=1).to(dev), critic(hidden=cfg.hidden, act_dim=1).to(dev)
    q1t.load_state_dict(q1.state_dict())
    q2t.load_state_dict(q2.state_dict())
    log_alpha = torch.zeros(1, device=dev, requires_grad=True)
    oa = torch.optim.Adam(actor.parameters(), lr=cfg.lr)
    oq = torch.optim.Adam(list(q1.parameters()) + list(q2.parameters()), lr=cfg.lr)
    oal = torch.optim.Adam([log_alpha], lr=cfg.lr)
    env = make_env(cfg, cfg.seed, dev)
    gen = torch.Generator(device=dev).manual_seed(cfg.seed * 7 + 5)
    cap, D = int(cfg.buffer), RF.OBS_DIM
    B = {"o": torch.zeros(cap, D), "a": torch.zeros(cap, 1), "r": torch.zeros(cap), "o2": torch.zeros(cap, D),
         "d": torch.zeros(cap)}
    ptr, full, n_shard, shard_from = 0, False, 0, 0
    curve, step, update, elapsed = Curve(), 0, 0, 0.0
    st = load_state(d / "train_state.pt", log)
    if st is not None:
        for m, k in ((actor, "actor"), (q1, "q1"), (q2, "q2"), (q1t, "q1t"), (q2t, "q2t")):
            m.load_state_dict(st[k])
        oa.load_state_dict(st["oa"])
        oq.load_state_dict(st["oq"])
        oal.load_state_dict(st["oal"])
        with torch.no_grad():
            log_alpha.copy_(st["log_alpha"])
        env.set_state(st["env"])
        step, update, elapsed, curve = st["step"], st["update"], st["elapsed"], Curve(st["curve"])
        n_shard = st["n_shard"]
        set_rng_state(st["rng"], {"gen": gen})
        ptr, full = _load_shards(d / "buffer", B, cap, n_shard)
        shard_from = ptr
        log(f"  [{label}] resumed at {step:,} decisions with {cap if full else ptr:,} transitions")
    obs = env.obs()
    N = env.cfg.n_envs
    t_save, t0 = time.time(), time.time()
    tgt_ent = -1.0
    snaps = _snapshot_marks(cfg)

    def save():
        nonlocal n_shard, shard_from
        n_shard = _write_shard(L, d / "buffer", B, shard_from, ptr, full, n_shard, cap)
        shard_from = ptr
        save_state(L, {"actor": actor.state_dict(), "q1": q1.state_dict(), "q2": q2.state_dict(),
                       "q1t": q1t.state_dict(), "q2t": q2t.state_dict(), "oa": oa.state_dict(), "oq": oq.state_dict(),
                       "oal": oal.state_dict(), "log_alpha": log_alpha.detach().clone(), "env": env.get_state(),
                       "step": step, "update": update, "elapsed": elapsed, "curve": curve.state(), "n_shard": n_shard,
                       "rng": rng_state({"gen": gen}), "cfg": asdict(cfg)}, d / "train_state.pt")
        _write_curve(L, d, curve, label)

    while step < cfg.budget_decisions and elapsed + (time.time() - t0) < cfg.budget_hours * 3600:
        with torch.no_grad():
            if step < cfg.learning_starts:
                a = torch.rand(N, generator=gen, device=dev) * 2 - 1
            else:
                actor.norm.update(obs)
                a, _, _ = actor.sample(obs, gen)
            o2, r, done, info = env.step(a)
            curve.add(info)
        idx = (torch.arange(N) + ptr) % cap
        B["o"][idx], B["a"][idx, 0], B["r"][idx] = obs.cpu(), a.cpu(), r.cpu()
        B["o2"][idx], B["d"][idx] = o2.cpu(), done.float().cpu()
        if ptr + N >= cap:
            full = True
        ptr = (ptr + N) % cap
        obs = o2.to(dev)
        step += N
        n_size = cap if full else ptr
        if step >= cfg.learning_starts:
            for _ in range(max(1, int(cfg.updates_per_step * N / 64))):
                j = torch.randint(n_size, (cfg.batch,), generator=torch.Generator().manual_seed(update * 1009 + cfg.seed))
                o, aa, rr, oo, dd = (B[k][j].to(dev) for k in ("o", "a", "r", "o2", "d"))
                alpha = log_alpha.exp().detach()
                with torch.no_grad():
                    a2, raw2, lp2 = actor.sample(oo, gen)
                    xo = torch.cat([actor.norm(oo), a2[:, None]], 1)
                    qt = torch.min(q1t(xo), q2t(xo)).squeeze(-1) - alpha * lp2
                    y = rr + cfg.gamma * (1 - dd) * qt
                xq = torch.cat([actor.norm(o), aa], 1)
                lq = ((q1(xq).squeeze(-1) - y) ** 2).mean() + ((q2(xq).squeeze(-1) - y) ** 2).mean()
                oq.zero_grad(set_to_none=True)
                lq.backward()
                oq.step()
                an, _, lpn = actor.sample(o, gen)
                xa = torch.cat([actor.norm(o), an[:, None]], 1)
                la = (alpha * lpn - torch.min(q1(xa), q2(xa)).squeeze(-1)).mean()
                oa.zero_grad(set_to_none=True)
                la.backward()
                oa.step()
                lal = -(log_alpha * (lpn.detach() + tgt_ent)).mean()
                oal.zero_grad(set_to_none=True)
                lal.backward()
                oal.step()
                with torch.no_grad():
                    for p, pt in zip(list(q1.parameters()) + list(q2.parameters()),
                                     list(q1t.parameters()) + list(q2t.parameters())):
                        pt.mul_(1 - cfg.tau).add_(cfg.tau * p)
                update += 1
        if curve.n >= N:
            curve.flush(step, update, elapsed + time.time() - t0, alpha=float(log_alpha.exp()))
        _maybe_snapshot(L, d, actor, step, cfg, snaps)
        due = time.time() - t_save > cfg.save_every_s or (cfg.save_every_updates and update > 0
                                                          and update % cfg.save_every_updates == 0)
        if due:
            elapsed += time.time() - t0
            t0 = time.time()
            save()
            t_save = time.time()
            fail_point(stage, update)
    elapsed += time.time() - t0
    t0 = time.time()
    save()
    save_policy(L, d / "policy.pt", actor, {"label": label, "decisions": step, "updates": update,
                                            "train_s": round(elapsed, 1), "cfg": asdict(cfg)})
    return {"decisions": step, "updates": update, "train_s": round(elapsed, 1)}


def _write_shard(L, folder: Path, B: dict, frm: int, to: int, full: bool, n_shard: int, cap: int) -> int:
    """Append the transitions added since the last checkpoint as shard_<n>.npz (ring order)."""
    if to == frm and not full:
        return n_shard
    idx = np.arange(frm, to) if to >= frm else np.concatenate([np.arange(frm, cap), np.arange(0, to)])
    if not len(idx):
        return n_shard
    buf = io.BytesIO()
    np.savez_compressed(buf, idx=idx, **{k: v[idx].numpy() for k, v in B.items()})
    path = folder / f"shard_{n_shard:05d}.npz"
    if hasattr(L, "save_large"):
        L.save_large(buf.getvalue(), path)
    else:
        L.atomic_bytes(path, buf.getvalue())
    return n_shard + 1


def _load_shards(folder: Path, B: dict, cap: int, n_shard: int) -> Tuple[int, bool]:
    """Replay the shards in order into the ring buffer (later shards overwrite earlier slots)."""
    ptr, full, total = 0, False, 0
    for i in range(n_shard):
        p = folder / f"shard_{i:05d}.npz"
        if not p.exists():
            continue
        z = np.load(p)
        idx = torch.from_numpy(z["idx"]).long()
        for k in B:
            B[k][idx] = torch.from_numpy(z[k])
        total += len(idx)
        ptr = int((int(idx[-1]) + 1) % cap)
    full = total >= cap
    return ptr, full


# ============================================================================
# The expert in the fast env, DAgger, BC
# ============================================================================
class FastExpert:
    """fnpeg on fast-env flights at the env's 5 s decisions (filters at beta^5); one FNPEG flight per env row."""

    def __init__(self, env: RF.FastEntryEnv, cfg: RLCfg, fnpeg_cfg: Optional[dict] = None):
        from jepa.guidance_fnpeg import FNPEG, FNPEGConfig
        fc = dict(fnpeg_cfg or {})
        fc["filter_beta"] = float(cfg.expert_beta_1hz) ** float(env.cfg.dt_s)
        self.law = FNPEG(FNPEGConfig(**fc))
        self.env = env
        self.flights = [self.law.new_flight(float(p), float(l)) for p, l in zip(env.phi_t.cpu(), env.lam_t.cpu())]

    def reset_rows(self, rows) -> None:
        for i in rows:
            i = int(i)
            self.flights[i] = self.law.new_flight(float(self.env.phi_t[i]), float(self.env.lam_t[i]))

    def label(self) -> torch.Tensor:
        """The expert's action in [-1, 1] for every env row (also advances its filters and warm start)."""
        ys = self.env.expert_rows()
        cmds = self.law.decide(self.flights, list(ys), list(self.env.t.cpu().numpy()))
        return torch.as_tensor(np.asarray(cmds) / RF.SIGMA_MAX_RAD, dtype=torch.float32).clamp(-1, 1)


def _rollouts(env: RF.FastEntryEnv, expert: FastExpert, actor: Optional[Actor], beta: float, n_episodes: int,
              gen: torch.Generator, dev) -> Tuple[np.ndarray, np.ndarray, dict]:
    """Episodes until n_episodes finished; the expert labels every visited state. Returns (obs, labels, stats)."""
    O, A = [], []
    done_n, succ, miss = 0, 0, []
    obs = env.obs()
    while done_n < n_episodes:
        lab = expert.label().to(dev)
        O.append(obs.cpu().numpy())
        A.append(lab.cpu().numpy())
        if actor is None:
            act = lab
        else:
            with torch.no_grad():
                pa = actor.deterministic(obs.to(dev))
            mix = torch.rand(len(lab), generator=gen, device=dev) < beta
            act = torch.where(mix, lab, pa)
        obs, r, done, info = env.step(act)
        idx = torch.nonzero(done).flatten().cpu()
        if len(idx):
            expert.reset_rows(idx)
            done_n += len(idx)
            succ += int(info["success_5km"][done].sum())
            miss += info["miss_km"][done].cpu().tolist()
    return (np.concatenate(O).astype(np.float32), np.concatenate(A).astype(np.float32),
            {"episodes": done_n, "success_5km_fast": succ / max(done_n, 1),
             "miss_km_fast": float(np.median(miss)) if miss else float("nan")})


def fit_supervised(actor: Actor, O: np.ndarray, A: np.ndarray, epochs: int, batch: int, lr: float, gen, dev,
                   log=print, label: str = "") -> float:
    """MSE on |a| and BCE on sign(a) (magsign), or MSE on a (gauss); refits the normalizer on the data."""
    o = torch.from_numpy(O).to(dev)
    a = torch.from_numpy(A).to(dev)
    with torch.no_grad():
        actor.norm.mean.zero_()
        actor.norm.var.fill_(1.0)
        actor.norm.count.fill_(1e-4)
        actor.norm.update(o)
    opt = torch.optim.Adam([p for n, p in actor.named_parameters() if n != "log_std"], lr=lr)
    n = len(o)
    tot = float("nan")
    for ep in range(epochs):
        perm = torch.randperm(n, generator=gen, device=dev)
        tot = 0.0
        for i in range(0, n, batch):
            j = perm[i:i + batch]
            mu, _, logit = actor.dist_params(o[j])
            if actor.kind == "magsign":
                loss = Fnn.mse_loss(mu, a[j].abs()) + Fnn.binary_cross_entropy_with_logits(logit, (a[j] >= 0).float())
            else:
                loss = Fnn.mse_loss(mu, a[j])
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()
            tot += float(loss) * len(j)
        tot /= n
    log(f"  [{label}] supervised fit on {n:,} states: loss {tot:.4f}")
    return tot


def dagger_train(L, d: Path, cfg: RLCfg, device, label: str, fnpeg_cfg: Optional[dict] = None, log=print,
                 stage: str = "rl_train") -> dict:
    """DAgger iterations; data/iter_<i>.npz marks an iteration done; policy.pt marks the run done."""
    dev = torch.device(device)
    torch.manual_seed(cfg.seed)
    gen = torch.Generator(device=dev).manual_seed(cfg.seed * 7 + 9)
    env = make_env(cfg, cfg.seed, dev, n_envs=min(cfg.n_envs, cfg.dagger_rollouts))
    actor = Actor(hidden=cfg.hidden, kind="magsign").to(dev)
    rows = []
    t0 = time.time()
    for i in range(cfg.dagger_iters):
        p = d / "data" / f"iter_{i:02d}.npz"
        if p.exists():
            continue
        # rebuild the learner from the data so far (deterministic: same seed, same data)
        actor, _ = _fit_from_iters(d, i, cfg, dev, label, log)
        if i > 0:
            save_policy(L, d / "snapshots" / f"policy_iter{i:02d}.pt", actor, {"iteration": i, "point": i})
        beta = 0.5 ** i
        st_path = d / "data" / f"env_{i:02d}.pt"
        if i > 0 and (d / "data" / f"env_{i - 1:02d}.pt").exists():
            env.set_state(torch.load(d / "data" / f"env_{i - 1:02d}.pt", map_location="cpu", weights_only=False))
        expert = FastExpert(env, cfg, fnpeg_cfg)
        gen = torch.Generator(device=dev).manual_seed(cfg.seed * 7919 + 9 + 1000 * i)
        O, A, s = _rollouts(env, expert, actor if i > 0 else None, beta, cfg.dagger_rollouts, gen, dev)
        buf = io.BytesIO()
        torch.save(env.get_state(), buf)
        L.atomic_bytes(st_path, buf.getvalue())
        buf = io.BytesIO()
        np.savez_compressed(buf, obs=O, act=A, iteration=i, beta=beta, **{k: np.asarray(v) for k, v in s.items()})
        L.atomic_bytes(p, buf.getvalue())
        log(f"  [{label}] iteration {i}: beta {beta:.3f}, {len(O):,} labeled states, learner+expert success@5km "
            f"{s['success_5km_fast']:.3f} ({(time.time() - t0) / 60:.1f} min)")
        fail_point(stage, i)
    actor, loss = _fit_from_iters(d, cfg.dagger_iters, cfg, dev, label, log)
    for i in range(cfg.dagger_iters):
        z = np.load(d / "data" / f"iter_{i:02d}.npz")
        rows.append({"iteration": i, "beta": float(z["beta"]), "states": int(len(z["obs"])),
                     "success_5km_fast": float(z["success_5km_fast"]), "miss_km_fast": float(z["miss_km_fast"])})
    import pandas as pd
    L.write_csv(pd.DataFrame(rows), d / "learning_curve.csv")
    save_policy(L, d / "policy.pt", actor, {"label": label, "iterations": cfg.dagger_iters, "fit_loss": loss,
                                            "cfg": asdict(cfg)})
    return {"iterations": cfg.dagger_iters, "states": int(sum(r["states"] for r in rows))}


def _fit_from_iters(d: Path, n_iter: int, cfg: RLCfg, dev, label: str, log) -> Tuple[Actor, float]:
    torch.manual_seed(cfg.seed + 100 * n_iter)
    actor = Actor(hidden=cfg.hidden, kind="magsign").to(dev)
    if n_iter == 0:
        return actor, float("nan")
    zs = [np.load(d / "data" / f"iter_{i:02d}.npz") for i in range(n_iter)]
    O = np.concatenate([z["obs"] for z in zs])
    A = np.concatenate([z["act"] for z in zs])
    g = torch.Generator(device=dev).manual_seed(cfg.seed * 31 + n_iter)
    loss = fit_supervised(actor, O, A, cfg.dagger_epochs, cfg.sup_batch, cfg.lr, g, dev, log, f"{label} fit {n_iter}")
    return actor, loss


def bc_dataset(L, d: Path, cfg: RLCfg, device, label: str, fnpeg_cfg: Optional[dict] = None, log=print) -> Path:
    """cfg.bc_episodes fnpeg episodes in the fast env, in chunks (data/bc_<j>.npz; resumable per chunk)."""
    dev = torch.device(device)
    n_chunk = max(1, int(math.ceil(cfg.bc_episodes / cfg.dagger_rollouts)))
    env = make_env(cfg, cfg.seed + 1000, dev, n_envs=min(cfg.n_envs, cfg.dagger_rollouts))
    gen = torch.Generator(device=dev).manual_seed(cfg.seed)
    for j in range(n_chunk):
        p = d / "data" / f"bc_{j:03d}.npz"
        if p.exists():
            continue
        if j > 0 and (d / "data" / f"bc_env_{j - 1:03d}.pt").exists():
            env.set_state(torch.load(d / "data" / f"bc_env_{j - 1:03d}.pt", map_location="cpu", weights_only=False))
        expert = FastExpert(env, cfg, fnpeg_cfg)
        O, A, s = _rollouts(env, expert, None, 1.0, cfg.dagger_rollouts, gen, dev)
        buf = io.BytesIO()
        torch.save(env.get_state(), buf)
        L.atomic_bytes(d / "data" / f"bc_env_{j:03d}.pt", buf.getvalue())
        buf = io.BytesIO()
        np.savez_compressed(buf, obs=O, act=A, **{k: np.asarray(v) for k, v in s.items()})
        L.atomic_bytes(p, buf.getvalue())
        log(f"  [{label}] expert chunk {j + 1}/{n_chunk}: {len(O):,} states, expert success@5km "
            f"{s['success_5km_fast']:.3f}")
    return d / "data"


def bcppo_train(L, d: Path, cfg: RLCfg, device, label: str, fnpeg_cfg: Optional[dict] = None, log=print,
                stage: str = "rl_train") -> dict:
    """BC from fnpeg episodes, then critic warm-up on frozen-actor rollouts, then PPO with a KL penalty."""
    dev = torch.device(device)
    clone_p = d / "clone.pt"
    if not clone_p.exists():
        bc_dataset(L, d, cfg, device, label, fnpeg_cfg, log)
        files = sorted((d / "data").glob("bc_*.npz"))
        O = np.concatenate([np.load(f)["obs"] for f in files])
        A = np.concatenate([np.load(f)["act"] for f in files])
        torch.manual_seed(cfg.seed)
        actor = Actor(hidden=cfg.hidden, kind="magsign").to(dev)
        g = torch.Generator(device=dev).manual_seed(cfg.seed * 13 + 1)
        loss = fit_supervised(actor, O, A, cfg.bc_epochs, cfg.sup_batch, cfg.lr, g, dev, log, f"{label} clone")
        save_policy(L, clone_p, actor, {"label": f"{label} clone", "states": int(len(O)), "fit_loss": loss})
    clone = load_policy(clone_p, dev)
    ref = load_policy(clone_p, dev)
    for p_ in ref.parameters():
        p_.requires_grad_(False)
    return ppo_train(L, d, cfg, device, label, actor=clone, ref=ref, log=log, stage=stage)


# ============================================================================
# One RL arm, one seed
# ============================================================================
def train_rl_dir(L, d: Path, arm: str, cfg: RLCfg, device, fnpeg_cfg: Optional[dict] = None, log=print,
                 label: Optional[str] = None) -> dict:
    """Train (or resume, or skip) one RL run into folder d; policy.pt marks it done, run.json describes it."""
    d = Path(d)
    if (d / "policy.pt").exists():
        p = d / "run.json"
        return json.loads(p.read_text()) if p.exists() else {"arm": arm, "seed": cfg.seed}
    L.guard(d).mkdir(parents=True, exist_ok=True)
    label = label or f"{arm} seed {cfg.seed}"
    t0 = time.time()
    if arm == "dagger_fnpeg":
        out = dagger_train(L, d, cfg, device, label, fnpeg_cfg, log)
    elif arm == "bcppo_fnpeg":
        out = bcppo_train(L, d, cfg, device, label, fnpeg_cfg, log)
    elif arm == "ppo_hist":
        out = ppo_train(L, d, cfg, device, label, log=log)
    elif arm == "sac_hist":
        out = sac_train(L, d, cfg, device, label, log=log)
    elif arm == "ppo_matched":
        out = ppo_full_train(L, d, cfg, device, label, log=log)
    else:
        raise ValueError(f"unknown RL arm {arm!r}")
    info = {"arm": arm, "seed": int(cfg.seed), **out, "wall_s": round(time.time() - t0, 1), "cfg": asdict(cfg)}
    L.write_json(d / "run.json", info)
    return info


def train_rl_arm(L, arm: str, seed: int, cfg: RLCfg, device, fnpeg_cfg: Optional[dict] = None, log=print) -> dict:
    return train_rl_dir(L, L.policy_dir(arm, seed), arm, replace(cfg, seed=int(seed)), device, fnpeg_cfg, log)


# ============================================================================
# PPO in the full simulator (ppo_matched): the world models' data budget, 6,000 episodes
# ============================================================================
def full_episode_specs(n: int, seed: int) -> List[dict]:
    """n training episodes drawn as the fast env draws them (entry, truth parameters, target on the footprint)."""
    from jepa.data import R_EARTH
    env = RF.FastEntryEnv(RF.FastEnvCfg(n_envs=1, pool=max(64, min(2048, n)), seed=int(seed)))
    out: List[dict] = []
    while len(out) < n:
        P = env._pool
        for i in range(len(P["x0"])):
            x = [float(v) for v in P["x0"][i]]
            out.append({"ic": {"altitude_m": x[0] - R_EARTH, "phi_deg": math.degrees(x[1]),
                               "lambda_deg": math.degrees(x[2]), "velocity_mps": x[3],
                               "gamma_deg": math.degrees(x[4]), "chi_deg": math.degrees(x[5])},
                        "domain": dict(zip(("density_scale", "lift_scale", "mass_scale"),
                                           [float(v) for v in P["dom"][i]])),
                        "target_phi": float(P["phi_t"][i]), "target_lam": float(P["lam_t"][i])})
        env._refill()
    return out[:n]


def target_config_file(phi_rad: float, lam_rad: float) -> str:
    """configs/default.json with the mission target moved, as a local cache file (never under configs/)."""
    import tempfile
    import zlib
    key = f"{phi_rad:.12f}_{lam_rad:.12f}"
    d = Path(tempfile.gettempdir()) / "jepa_v5_targets"
    path = d / f"t_{zlib.crc32(key.encode()):08x}.json"
    if not path.exists():
        from jepa.benchmark_v5 import SIM_DIR
        cfg = json.loads((SIM_DIR / "configs" / "default.json").read_text())
        cfg["mission"]["target_phi_deg"] = math.degrees(phi_rad)
        cfg["mission"]["target_lambda_deg"] = math.degrees(lam_rad)
        d.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
        tmp.write_text(json.dumps(cfg))
        os.replace(tmp, path)
    return str(path)


def _full_episode(arg) -> dict:
    """One training episode in the full simulator (rl_env), 5 s decisions, the fast env's reward."""
    hp, state, spec, seed, fast_cfg = arg
    torch.set_num_threads(1)
    import rl_env
    from jepa.data import G0, SIGMA_MAX_RAD, apply_domain, physical_state
    from jepa.fly_v3 import great_circle_km
    actor = Actor(**hp)
    actor.load_state_dict(state)
    actor.eval()
    env = rl_env.OrionEntryEnv(config_path=target_config_file(spec["target_phi"], spec["target_lam"]),
                               controller_mode="policy", dispersion=False, max_episode_steps=1400,
                               fast_atmosphere=True)
    apply_domain(env.params, spec["domain"])
    env.reset(seed=0, options={"ic_override": spec["ic"]})
    c = RF.FastEnvCfg(**fast_cfg)
    fo = RF.FullObs()
    gen = torch.Generator().manual_seed(int(seed))
    rtg = lambda: great_circle_km(float(env.x[1]), float(env.x[2]), spec["target_phi"], spec["target_lam"])
    O, RAW, LP, R = [], [], [], []
    u_prev, done, info = 0.0, False, {}
    while not done:
        o = fo(env, physical_state(env))
        with torch.no_grad():
            a, raw, lp = actor.sample(o[None], gen)
        a0 = float(a[0])
        cmd = a0 * SIGMA_MAX_RAD
        phi0 = -rtg()
        gmax = 0.0
        for _ in range(int(c.dt_s)):
            _, _, term, trunc, info = env.step(np.array([a0], dtype=np.float32))
            y = physical_state(env)
            gmax = max(gmax, math.hypot(float(y[8]), float(y[9])) / G0)
            if term or trunc:
                done = True
                break
        crossed = done and str(info.get("outcome")) == "drogue"
        r = -c.w_smooth_km_per_rad * abs(cmd - u_prev) - c.w_g_km * max(0.0, gmax - c.g_limit) ** 2
        if not crossed:
            r -= c.w_gamma_km_per_deg2 * max(0.0, abs(math.degrees(float(env.x[4]))) - c.gamma_limit_deg) ** 2
        if done:
            r += -rtg() - (0.0 if crossed else c.nocross_km)
        if c.shaping:
            r += c.gamma_rl * (0.0 if done else -rtg()) - phi0
        O.append(o)
        RAW.append(raw[0])
        LP.append(lp[0])
        R.append(r * c.reward_scale)
        u_prev = cmd
    miss = rtg()
    return {"obs": torch.stack(O), "raw": torch.stack(RAW), "lp": torch.stack(LP),
            "r": torch.tensor(R, dtype=torch.float32), "miss_km": miss, "crossed": bool(crossed),
            "success_5km": bool(crossed and miss < 5.0), "ep_return": float(sum(R)), "outcome": str(info.get("outcome"))}


def ppo_full_train(L, d: Path, cfg: RLCfg, device, label: str, log=print, stage: str = "rl_train") -> dict:
    """
    PPO in the full simulator for cfg.full_episodes episodes. Every update collects cfg.full_per_update complete
    episodes (in a process pool, each with its own noise stream), so a checkpoint never cuts an episode and a
    resumed run continues exactly. Resumable from d/train_state.pt; writes d/policy.pt.
    """
    import multiprocessing as mp
    import zlib
    dev = torch.device(device)
    torch.manual_seed(cfg.seed)
    actor = Actor(hidden=cfg.hidden, kind="gauss").to(dev)
    vf = critic(hidden=cfg.hidden).to(dev)
    opt = torch.optim.Adam(list(actor.parameters()) + list(vf.parameters()), lr=cfg.lr, eps=1e-5)
    gen = torch.Generator(device=dev).manual_seed(cfg.seed * 7 + 13)
    specs = full_episode_specs(cfg.full_episodes, cfg.seed + 5000)
    curve, step, update, n_ep, elapsed = Curve(), 0, 0, 0, 0.0
    st = load_state(d / "train_state.pt", log)
    if st is not None:
        actor.load_state_dict(st["actor"])
        vf.load_state_dict(st["critic"])
        opt.load_state_dict(st["opt"])
        step, update, n_ep, elapsed, curve = st["step"], st["update"], st["episodes"], st["elapsed"], Curve(st["curve"])
        set_rng_state(st["rng"], {"gen": gen})
        log(f"  [{label}] resumed at update {update} ({n_ep} episodes)")
    fast_cfg = asdict(RF.FastEnvCfg(n_envs=1, gamma_rl=cfg.gamma))
    workers = int(cfg.full_workers or os.cpu_count() or 1)
    pool = mp.get_context("spawn").Pool(workers) if workers > 1 else None
    t_save, t0 = time.time(), time.time()

    def save():
        save_state(L, {"actor": actor.state_dict(), "critic": vf.state_dict(), "opt": opt.state_dict(),
                       "step": step, "update": update, "episodes": n_ep, "elapsed": elapsed, "curve": curve.state(),
                       "rng": rng_state({"gen": gen}), "cfg": asdict(cfg)}, d / "train_state.pt")
        _write_curve(L, d, curve, label)

    try:
        while n_ep < cfg.full_episodes and elapsed + (time.time() - t0) < cfg.budget_hours * 3600:
            batch = specs[n_ep:n_ep + cfg.full_per_update]
            sd = {k: v.detach().cpu() for k, v in actor.state_dict().items()}
            args = [(actor.hparams, sd, sp, zlib.crc32(f"ppo_matched|{cfg.seed}|{n_ep + i}".encode()), fast_cfg)
                    for i, sp in enumerate(batch)]
            res = pool.map(_full_episode, args) if pool is not None else [_full_episode(a) for a in args]
            O = torch.cat([r["obs"] for r in res]).to(dev)
            actor.norm.update(O)
            ADV, RET = [], []
            with torch.no_grad():
                for r in res:
                    v = vf(actor.norm(r["obs"].to(dev))).squeeze(-1)
                    rr = r["r"].to(dev)
                    adv = torch.zeros_like(rr)
                    last = torch.zeros((), device=dev)
                    for k in reversed(range(len(rr))):
                        nxt = v[k + 1] if k + 1 < len(rr) else torch.zeros((), device=dev)
                        delta = rr[k] + cfg.gamma * nxt - v[k]
                        last = delta + cfg.gamma * cfg.gae_lambda * last
                        adv[k] = last
                    ADV.append(adv)
                    RET.append(adv + v)
            _ppo_update(actor, vf, opt, cfg, gen, None, False, O, torch.cat([r["raw"] for r in res]).to(dev),
                        torch.cat([r["lp"] for r in res]).to(dev), torch.cat(ADV), torch.cat(RET))
            n_ep += len(batch)
            step += int(len(O))
            update += 1
            for r in res:
                curve.ret += r["ep_return"]
                curve.succ += int(r["success_5km"])
                curve.miss += min(r["miss_km"], 5000.0)
                curve.n += 1
            curve.flush(step, update, elapsed + time.time() - t0, episodes=n_ep)
            rw = curve.rows[-1]
            log(f"  [{label}] update {update}: {n_ep}/{cfg.full_episodes} episodes, return {rw['ep_return']:.3f} "
                f"success@5km {rw['success_5km_fast']:.3f} ({(elapsed + time.time() - t0) / 60:.1f} min)")
            due = time.time() - t_save > cfg.save_every_s or (cfg.save_every_updates and update % cfg.save_every_updates == 0)
            if due:
                elapsed += time.time() - t0
                t0 = time.time()
                save()
                t_save = time.time()
                fail_point(stage, update)
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    elapsed += time.time() - t0
    t0 = time.time()
    save()
    save_policy(L, d / "policy.pt", actor, {"label": label, "decisions": step, "updates": update, "episodes": n_ep,
                                            "train_s": round(elapsed, 1), "cfg": asdict(cfg)})
    return {"decisions": step, "updates": update, "episodes": n_ep, "train_s": round(elapsed, 1)}


@torch.no_grad()
def fast_success(policy: Actor, n_episodes: int = 512, seed: int = 12345, device="cpu") -> dict:
    """success@5km of a policy in the fast env (the fast side of the fast-vs-full gap)."""
    env = RF.FastEntryEnv(RF.FastEnvCfg(n_envs=min(256, n_episodes), seed=seed), device)
    obs = env.obs()
    n, succ, miss = 0, 0, []
    while n < n_episodes:
        obs, r, done, info = env.step(policy.deterministic(obs.to(device)))
        if bool(done.any()):
            n += int(done.sum())
            succ += int(info["success_5km"][done].sum())
            miss += info["miss_km"][done].cpu().tolist()
    return {"episodes": n, "success_5km_fast": succ / n, "miss_km_fast_median": float(np.median(miss))}


def fast_success_on(policy: Actor, conds: Sequence[dict], device="cpu") -> dict:
    """
    success@5km of a policy in the fast env on given conditions (entry state, truth parameters and the target of
    each condition's config): the fast side of the fast-vs-full gap on the same cases as the full-simulator flights.
    """
    sim = Path(__file__).resolve().parents[1]
    n = len(conds)
    env = RF.FastEntryEnv(RF.FastEnvCfg(n_envs=n, seed=0, pool=n), device)
    f64 = dict(dtype=torch.float64, device=env.dev)
    x0 = env._x_from_ic([c["ic"] for c in conds])
    dom = torch.as_tensor([[c["domain"]["density_scale"], c["domain"]["lift_scale"], c["domain"]["mass_scale"]]
                           for c in conds], **f64)
    tg = []
    for c in conds:
        m = json.loads((sim / (c.get("config_path") or "configs/default.json")).read_text())["mission"]
        tg.append([math.radians(m["target_phi_deg"]), math.radians(m["target_lambda_deg"])])
    tg = torch.as_tensor(tg, **f64)
    env._pool = {"x0": x0, "S_ref": env.area0 * dom[:, 0], "mass": env.mass0 * dom[:, 2], "ls": dom[:, 1],
                 "phi_t": tg[:, 0], "lam_t": tg[:, 1], "dom": dom}
    env._pool_i = 0
    env.reset_rows(torch.arange(n, device=env.dev))
    obs = env.obs()
    first = np.zeros(n, dtype=bool)
    succ, miss = np.zeros(n, dtype=bool), np.full(n, np.nan)
    with torch.no_grad():
        while not first.all():
            obs, _, done, info = env.step(policy.deterministic(obs.to(device)))
            d = done.cpu().numpy() & ~first
            if d.any():
                succ[d] = info["success_5km"].cpu().numpy()[d]
                miss[d] = info["miss_km"].cpu().numpy()[d]
                first |= d
    return {"n": n, "success_5km_fast": float(succ.mean()), "miss_km_fast_median": float(np.nanmedian(miss)),
            "success_by_case": succ.astype(int).tolist()}


# ============================================================================
# Legacy RL (section 4.6)
# ============================================================================
def legacy_inventory(root: Path) -> List[dict]:
    """Every file under ReEntryAI/{models,checkpoints,data} (os.scandir per folder)."""
    import re
    rows = []
    for sub in ("models", "checkpoints", "data"):
        base = Path(root) / sub
        if not base.exists():
            rows.append({"path": str(base), "run_name": None, "steps": None, "size_bytes": None, "date": None,
                         "kind": "not found"})
            continue
        stack = [base]
        while stack:
            cur = stack.pop()
            with os.scandir(cur) as it:
                for e in it:
                    if e.is_dir():
                        stack.append(Path(e.path))
                        continue
                    name = e.name
                    m = re.search(r"_(\d+)_steps", name)
                    st = e.stat()
                    rows.append({"path": e.path, "run_name": Path(cur).name if sub == "checkpoints" else Path(name).stem,
                                 "steps": int(m.group(1)) if m else None, "size_bytes": st.st_size,
                                 "date": time.strftime("%Y-%m-%d %H:%M", time.localtime(st.st_mtime)),
                                 "kind": "final model" if (sub == "models" and name.endswith(".zip")) else sub})
    return rows


class LegacyPolicy:
    """An SB3 policy (zip) with its own settings: 1 Hz, rl_env's observation, deterministic."""
    period = 1
    legacy = True

    def __init__(self, path: str):
        from stable_baselines3 import PPO
        self.model = PPO.load(path, device="cpu")
        self.obs_dim = int(self.model.observation_space.shape[0])

    def act(self, obs: np.ndarray) -> float:
        a, _ = self.model.predict(obs, deterministic=True)
        return float(np.asarray(a).reshape(-1)[0])
