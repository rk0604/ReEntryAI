"""
jepa/rl2_v5.py

Notebook 2's RL arms (jepa/plan_rl_v5.md): residual RL on FNPEG, Liu-style reference tracking, and the cheap pilot
gates that every arm passes before any L4 time. ppo_hist itself is rl_v5's, unchanged.

    ResidExpert    FNPEG (the tuned configuration) on every fast-env row at the 5 s decisions, returning its command
                   in float64 with its filter state (K_D, K_L) and sign; in-process or sharded over worker processes
                   (one FNPEG flight per row, so both give the same numbers); its flights' state round-trips for an
                   exact resume. FNPEG's memory is never synced: the flown sign is always its own.
    ResidEnv       the fast env driven by sigma = s_F clip(|sigma_F| + b a, 0, 90 deg) (plan section 1.1); observation
                   = the 181 common features + |sigma_F|/90 deg, s_F, K_D, K_L; reward = the fast env's
    LiuEnv         the fast env with Liu et al.'s reference velocity vector (Eqs. 16-21), heading deadband (Eq. 19,
                   FNPEG's tuned deadband), reward (Eqs. 31-37) and training terminations (load > 4.5 g, heading away);
                   action = bank magnitude only (plan sections 1.2-1.3)
    ppo_train_env  rl_v5.ppo_train for any of these environments (its update, curve, checks and snapshots), with the
                   normalizer also updated during the critic burn-in
    gates          RI (residual at init equals FNPEG bit for bit) and RP (pilot competence on 512 fixed episodes
                   against FNPEG alone on the same episodes), cost per million decisions, routing
    ResidLaw,      the trained arms in the full simulator behind fly_v5's 1 Hz runner: FNPEG every second as the fnpeg
    LiuLaw         arm, the policy every 5 s (held)
"""

from __future__ import annotations

import copy
import io
import json
import math
import os
import time
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch

from jepa import rl_v5 as RL
from jepa import rlenv_fast_v5 as RF
from jepa.data import G0, R_EARTH, SIGMA_MAX_RAD
from jepa.models_v5 import fail_point, load_state, rng_state, save_state, set_rng_state

ARMS = ("resid_fnpeg", "liu_track", "ppo_hist")
RESID_EXTRA = 4
RESID_OBS_DIM = RF.OBS_DIM + RESID_EXTRA
LIU_EXTRA = 11
LIU_OBS_DIM = RF.OBS_DIM + LIU_EXTRA

# arm settings that are not RLCfg fields (plan section 4.2); a policy file stores the ones it was trained with
ARM_DEFAULTS = {"resid_fnpeg": {"bound_deg": 20.0, "init_std": 0.1},
                "liu_track": {"ds_frac": 0.3},
                "ppo_hist": {}}


# ============================================================================
# The residual: FNPEG's command, corrected in magnitude
# ============================================================================
def compose(cmd: float, sign: float, delta: float, smax: float = SIGMA_MAX_RAD) -> float:
    """s_F clip(|sigma_F| + delta, 0, smax). With delta == 0 it returns cmd exactly (copysign(|c|, sign(c)) = c)."""
    return math.copysign(min(smax, max(0.0, abs(cmd) + delta)), sign)


def compose_t(cmd: torch.Tensor, sign: torch.Tensor, delta: torch.Tensor, smax: float = SIGMA_MAX_RAD) -> torch.Tensor:
    return torch.copysign((cmd.abs() + delta).clamp(0.0, smax), sign)


def sign_eff(cmd: float, mem: float) -> float:
    """The sign the residual keeps: FNPEG's command's, else its sign memory, else +."""
    if cmd != 0.0:
        return 1.0 if cmd > 0.0 else -1.0
    return 1.0 if mem >= 0.0 else -1.0


def _query_rows(law, flights: Sequence, ys: np.ndarray, ts: Sequence[float]) -> np.ndarray:
    cmds = law.decide(list(flights), list(ys), list(ts))
    return np.array([[c, f.K_D, f.K_L, sign_eff(float(c), f.sign)] for c, f in zip(cmds, flights)],
                    dtype=np.float64).reshape(-1, 4)


def _resid_worker(conn, fc: dict) -> None:
    """A ResidExpert process: the FNPEG flights of its rows."""
    torch.set_num_threads(1)
    from jepa.guidance_fnpeg import FNPEG, FNPEGConfig
    law = FNPEG(FNPEGConfig(**fc))
    flights: Dict[int, object] = {}
    while True:
        msg = conn.recv()
        if msg is None:
            break
        op = msg[0]
        if op == "q":
            _, resets, rows, ys, ts = msg
            for r, ph, la in resets:
                flights[int(r)] = law.new_flight(float(ph), float(la))
            conn.send(_query_rows(law, [flights[int(r)] for r in rows], ys, ts) if len(rows) else
                      np.zeros((0, 4)))
        elif op == "get":
            conn.send(copy.deepcopy(flights))
        elif op == "set":
            flights = dict(msg[1])
            conn.send(True)
    conn.close()


class ResidExpert:
    """FNPEG on every env row; query(active) -> [n, 4] = (command rad, K_D, K_L, sign) (zeros on inactive rows)."""

    def __init__(self, env: RF.FastEntryEnv, cfg: RL.RLCfg, fnpeg_cfg: Optional[dict], workers: int = 1):
        from jepa.guidance_fnpeg import FNPEG, FNPEGConfig
        self.env = env
        self.fc = RL._expert_fc(env, cfg, fnpeg_cfg)
        n = int(env.cfg.n_envs)
        w = max(1, min(int(workers or (os.cpu_count() or 2)), n))
        self.parallel = w > 1
        self.pending_resets: List[List[tuple]] = []
        if not self.parallel:
            self.law = FNPEG(FNPEGConfig(**self.fc))
            self.flights: List[object] = [None] * n
        else:
            import multiprocessing as mp
            self.shards = [a for a in np.array_split(np.arange(n), w) if len(a)]
            self.owner = np.zeros(n, dtype=np.int64)
            for k, a in enumerate(self.shards):
                self.owner[a] = k
            ctx = mp.get_context("spawn")
            self.conns, self.procs = [], []
            for _ in self.shards:
                a, b = ctx.Pipe()
                pr = ctx.Process(target=_resid_worker, args=(b, self.fc), daemon=True)
                pr.start()
                self.conns.append(a)
                self.procs.append(pr)
            self.pending_resets = [[] for _ in self.shards]
        self.reset_rows(range(n))

    def reset_rows(self, rows) -> None:
        for i in rows:
            i = int(i)
            ph, la = float(self.env.phi_t[i]), float(self.env.lam_t[i])
            if self.parallel:
                self.pending_resets[int(self.owner[i])].append((i, ph, la))
            else:
                self.flights[i] = self.law.new_flight(ph, la)

    def query(self, active: Optional[np.ndarray] = None) -> np.ndarray:
        n = int(self.env.cfg.n_envs)
        act = np.ones(n, dtype=bool) if active is None else np.asarray(active, dtype=bool)
        ys = self.env.expert_rows()
        ts = self.env.t.cpu().numpy()
        out = np.zeros((n, 4), dtype=np.float64)
        out[:, 1:3] = 1.0
        out[:, 3] = 1.0
        if not self.parallel:
            rows = np.nonzero(act)[0]
            if len(rows):
                out[rows] = _query_rows(self.law, [self.flights[i] for i in rows], ys[rows], ts[rows])
        else:
            for k, a in enumerate(self.shards):
                rows = a[act[a]]
                self.conns[k].send(("q", self.pending_resets[k], rows, ys[rows], ts[rows]))
                self.pending_resets[k] = []
            for k, a in enumerate(self.shards):
                rows = a[act[a]]
                res = self.conns[k].recv()
                if len(rows):
                    out[rows] = res
        try:
            from jepa import runguard_v5 as RG
            RG.touch("fnpeg labels")
        except Exception:
            pass
        return out

    def get_state(self) -> dict:
        if not self.parallel:
            return {"flights": copy.deepcopy(self.flights)}
        fl: Dict[int, object] = {}
        for k, c in enumerate(self.conns):
            c.send(("q", self.pending_resets[k], np.zeros(0, dtype=np.int64), np.zeros((0, 10)), np.zeros(0)))
            self.pending_resets[k] = []
            c.recv()
            c.send(("get",))
            fl.update(c.recv())
        return {"flights": [fl.get(i) for i in range(int(self.env.cfg.n_envs))]}

    def set_state(self, st: dict) -> None:
        fl = st["flights"]
        if not self.parallel:
            self.flights = copy.deepcopy(list(fl))
            return
        for k, a in enumerate(self.shards):
            self.conns[k].send(("set", {int(i): fl[int(i)] for i in a}))
            self.pending_resets[k] = []
            self.conns[k].recv()

    def close(self) -> None:
        if not self.parallel:
            return
        for c in self.conns:
            try:
                c.send(None)
            except (OSError, BrokenPipeError):
                pass
        for pr in self.procs:
            pr.join(timeout=10)
            if pr.is_alive():
                pr.terminate()
        self.conns, self.procs = [], []

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


def resid_features(q: np.ndarray) -> np.ndarray:
    """|sigma_F| / 90 deg, s_F, K_D, K_L."""
    return np.stack([np.abs(q[:, 0]) / SIGMA_MAX_RAD, q[:, 3], q[:, 1], q[:, 2]], -1)


def set_pool(env: RF.FastEntryEnv, conds: Sequence[dict]) -> None:
    """The env's rows reset to given conditions (entry state, truth parameters, the target of each condition's
    config), as rl_v5.fast_success_on: the fast side of the fast-vs-full gap on the full simulator's cases."""
    sim = Path(__file__).resolve().parents[1]
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
    env.reset_rows(torch.arange(len(conds), device=env.dev))


class ResidEnv:
    """The fast env flown by FNPEG plus a bounded magnitude correction. step(a in [-1, 1]) -> (obs, r, done, info)."""

    def __init__(self, fast: RF.FastEnvCfg, device, cfg: RL.RLCfg, fnpeg_cfg: Optional[dict], bound_deg: float,
                 workers: int = 1, conds: Optional[Sequence[dict]] = None):
        self.env = RF.FastEntryEnv(fast, device)
        if conds:
            set_pool(self.env, conds)
        self.cfg = self.env.cfg
        self.dev = self.env.dev
        self.obs_dim = RESID_OBS_DIM
        self.bound = math.radians(float(bound_deg))
        self.expert = ResidExpert(self.env, cfg, fnpeg_cfg, workers)
        self.q = self.expert.query()
        self.last_cmd = torch.zeros(int(fast.n_envs), dtype=torch.float64, device=self.dev)

    def obs(self) -> torch.Tensor:
        f = torch.as_tensor(resid_features(self.q), dtype=torch.float32, device=self.dev)
        return torch.cat([self.env.obs(), f], 1)

    def command(self, action: torch.Tensor) -> torch.Tensor:
        q = torch.as_tensor(self.q, dtype=torch.float64, device=self.dev)
        a = action.to(self.dev, torch.float64).clamp(-1.0, 1.0)
        return compose_t(q[:, 0], q[:, 3], self.bound * a)

    @torch.no_grad()
    def step(self, action: torch.Tensor, active: Optional[np.ndarray] = None):
        cmd = self.command(action)
        if active is not None:
            cmd = torch.where(torch.as_tensor(active, device=self.dev), cmd, torch.zeros_like(cmd))
        self.last_cmd = cmd
        _, r, done, info = self.env.step(cmd / SIGMA_MAX_RAD)
        idx = torch.nonzero(done).flatten().cpu()
        if len(idx):
            self.expert.reset_rows(idx)
        self.q = self.expert.query(active)
        return self.obs(), r, done, info

    def get_state(self) -> dict:
        return {"env": self.env.get_state(), "expert": self.expert.get_state(), "q": self.q.copy()}

    def set_state(self, st: dict) -> None:
        self.env.set_state(st["env"])
        self.expert.set_state(st["expert"])
        self.q = np.array(st["q"])

    def close(self) -> None:
        self.expert.close()


# ============================================================================
# Liu-style reference tracking
# ============================================================================
@dataclass
class LiuCfg:
    ds_frac: float = 0.3                     # delta_s = ds_frac * s0 (Eq. 16), tuned in {0.2, 0.3, 0.4}
    vf_mps: float = 147.7                    # FNPEG's tuned Vf at the drogue altitude
    hf_m: float = 7_600.0
    db_c0_deg: float = 1.0                   # the heading deadband (Eq. 19's delta_psi) = FNPEG's tuned one
    db_c1_deg_per_mps: float = 5.2e-4
    n_max_g: float = 4.5
    b_shaping: float = 0.25                  # Liu et al. Table 2
    b_alt: float = 0.05
    b_ctrl: float = 0.02
    b_path: float = 0.05
    b_rev: float = 0.1
    bonus: float = 10.0
    kappa_v: float = 0.1                     # kappa_v = 0.1 v_ref
    kappa_alt: float = 0.1                   # kappa_alt = 0.1 max(|sin gamma|, sin floor)
    kappa_alt_floor_deg: float = 0.5
    success_km: float = 5.0
    train_terms: bool = True                 # load > n_max and heading away end a training episode (never an eval)


def liu_cfg_from(fnpeg_cfg: Optional[dict], **kw) -> LiuCfg:
    from jepa.guidance_fnpeg import FNPEGConfig
    f = FNPEGConfig(**(fnpeg_cfg or {}))
    return LiuCfg(vf_mps=float(f.vf_mps), hf_m=float(f.hf_m), db_c0_deg=float(f.deadband_c0_deg),
                  db_c1_deg_per_mps=float(f.deadband_c1_deg_per_mps), **kw)


def liu_terms(x: torch.Tensor, phi_t, lam_t, s0, v0, lc: LiuCfg) -> dict:
    """Liu's reference and errors at states x [B, 6] = (r, phi, lam, V, gamma, chi), float64."""
    r, phi, lam, V, gam, chi = x.unbind(-1)
    h = r - R_EARTH
    ang, course = RF.range_course(phi, lam, phi_t, lam_t)
    s = ang * R_EARTH
    ds = (lc.ds_frac * s0).clamp_min(1.0)
    v_ref = v0 - (v0 - lc.vf_mps) * torch.exp(-s / ds)                                   # Eq. 16
    g_ref = torch.asin(((lc.hf_m - h) / s.clamp_min(1.0)).clamp(-1.0, 1.0))               # Eq. 18
    dpsi = RF.wrap(chi - course)                                                           # FNPEG's dpsi = chi - Psi
    db = torch.deg2rad(lc.db_c0_deg + lc.db_c1_deg_per_mps * V)
    psi_ref = torch.where(dpsi.abs() <= db, chi, course)                                   # Eq. 19
    vc = torch.stack([V * torch.cos(gam) * torch.cos(chi), V * torch.sin(gam), V * torch.cos(gam) * torch.sin(chi)], -1)
    vr = torch.stack([v_ref * torch.cos(g_ref) * torch.cos(psi_ref), v_ref * torch.sin(g_ref),
                      v_ref * torch.cos(g_ref) * torch.sin(psi_ref)], -1)                  # Eq. 21
    return {"s": s, "v_ref": v_ref, "g_ref": g_ref, "dpsi": dpsi, "db": db, "verr": vc - vr, "h": h, "V": V, "gam": gam}


def liu_shaping(T: dict, n_g: torch.Tensor, lc: LiuCfg) -> torch.Tensor:
    """r_shaping + r_alt + r_path (Eqs. 31, 32, 34) at a state."""
    kv = (lc.kappa_v * T["v_ref"]).clamp_min(1.0)
    r = lc.b_shaping * torch.exp(-(T["verr"] ** 2).sum(-1) / kv ** 2)
    ka = lc.kappa_alt * torch.maximum(torch.sin(T["gam"]).abs(),
                                      torch.full_like(T["gam"], math.sin(math.radians(lc.kappa_alt_floor_deg))))
    r = r + lc.b_alt * torch.exp(-((torch.sin(T["gam"]) - torch.sin(T["g_ref"])) / ka) ** 2)
    r = r - lc.b_path * torch.exp((n_g - lc.n_max_g) / lc.n_max_g)
    return r


def liu_features(T: dict, n_g, s0, sign, nrev, lc: LiuCfg) -> torch.Tensor:
    vr = T["v_ref"].clamp_min(1.0)
    return torch.cat([T["verr"] / vr[:, None], torch.stack([
        T["s"] / s0.clamp_min(1.0), T["V"] / vr - 1.0, T["dpsi"] / math.pi, (T["h"] - lc.hf_m) / 1e5,
        n_g / lc.n_max_g, T["g_ref"] / (math.pi / 2), sign, nrev / 10.0], -1)], -1)


def liu_sign(sign: torch.Tensor, dpsi: torch.Tensor, db: torch.Tensor) -> torch.Tensor:
    """Eq. 19 applied to the bank sign (FNPEG's convention: a positive bank turns right): set at the first decision,
    flipped toward the line of sight once |dpsi| exceeds the deadband, kept inside it."""
    s = torch.where(sign == 0, torch.where(dpsi >= 0, 1.0, -1.0).to(sign), sign)
    s = torch.where(dpsi > db, torch.ones_like(s), s)
    return torch.where(dpsi < -db, -torch.ones_like(s), s)


class LiuEnv:
    """The fast env with Liu-style reference tracking: action a -> |sigma| = 45 deg (a + 1), sign by Eq. 19."""

    def __init__(self, fast: RF.FastEnvCfg, device, lc: LiuCfg, conds: Optional[Sequence[dict]] = None):
        self.env = RF.FastEntryEnv(fast, device)
        if conds:
            set_pool(self.env, conds)
        self.cfg = self.env.cfg
        self.dev = self.env.dev
        self.lc = lc
        self.obs_dim = LIU_OBS_DIM
        n = int(fast.n_envs)
        f64 = dict(dtype=torch.float64, device=self.dev)
        self.s0, self.v0, self.sign, self.nrev, self.mag_prev, self.ep_ret = (torch.zeros(n, **f64) for _ in range(6))
        self._start(torch.arange(n, device=self.dev))

    def _start(self, rows) -> None:
        x = self.env.x[rows]
        ang, _ = RF.range_course(x[:, 1], x[:, 2], self.env.phi_t[rows], self.env.lam_t[rows])
        self.s0[rows] = ang * R_EARTH
        self.v0[rows] = x[:, 3]
        for t in (self.sign, self.nrev, self.mag_prev, self.ep_ret):
            t[rows] = 0.0

    def _terms(self) -> dict:
        return liu_terms(self.env.x, self.env.phi_t, self.env.lam_t, self.s0, self.v0, self.lc)

    def obs(self) -> torch.Tensor:
        T = self._terms()
        f = liu_features(T, self.env.gload, self.s0, self.sign, self.nrev, self.lc).float()
        return torch.cat([self.env.obs(), f], 1)

    @torch.no_grad()
    def step(self, action: torch.Tensor):
        lc = self.lc
        a = action.to(self.dev, torch.float64).clamp(-1.0, 1.0)
        mag = 0.5 * (a + 1.0) * SIGMA_MAX_RAD
        T = self._terms()
        new = liu_sign(self.sign, T["dpsi"], T["db"])
        rev = ((new != self.sign) & (self.sign != 0)).double()
        self.sign, self.nrev = new, self.nrev + rev
        r_fix = -lc.b_ctrl * ((mag - self.mag_prev) / SIGMA_MAX_RAD) ** 2 - lc.b_rev * rev
        self.mag_prev = mag
        _, _, done, info = self.env.step(self.sign * mag / SIGMA_MAX_RAD)
        succ = info["success_5km"] if lc.success_km == 5.0 else done & info["crossed"] & (info["miss_km"] < lc.success_km)
        miss = info["miss_km"].clone()
        T2 = self._terms()
        live = ~done
        shaped = liu_shaping(T2, self.env.gload, lc)
        forced = torch.zeros_like(done)
        if lc.train_terms:
            forced = live & ((self.env.gload > lc.n_max_g) | (T2["dpsi"].abs() > math.pi / 2))
        r = r_fix + torch.where(live, shaped, torch.zeros_like(shaped)) + lc.bonus * succ.double()
        self.ep_ret += r
        done_all = done | forced
        miss = torch.where(forced, T2["s"] / 1000.0, miss)
        info = {"done": done_all.clone(), "success_5km": succ & done, "miss_km": miss,
                "ep_return": self.ep_ret.clone(), "crossed": info["crossed"], "forced": forced}
        if bool(forced.any()):
            self.env.reset_rows(torch.nonzero(forced).flatten())
        idx = torch.nonzero(done_all).flatten()
        if len(idx):
            self._start(idx)
        return self.obs(), r.float(), done_all, info

    def get_state(self) -> dict:
        return {"env": self.env.get_state(), **{k: getattr(self, k).clone().cpu() for k in
                ("s0", "v0", "sign", "nrev", "mag_prev", "ep_ret")}}

    def set_state(self, st: dict) -> None:
        self.env.set_state(st["env"])
        for k in ("s0", "v0", "sign", "nrev", "mag_prev", "ep_ret"):
            setattr(self, k, st[k].to(self.dev))

    def close(self) -> None:
        pass


# ============================================================================
# Actors and environments per arm
# ============================================================================
def make_actor(arm: str, cfg: RL.RLCfg, arm_kw: dict) -> RL.Actor:
    if arm == "resid_fnpeg":
        a = RL.Actor(obs_dim=RESID_OBS_DIM, hidden=cfg.hidden, kind="gauss")
        with torch.no_grad():
            a.net[-1].weight.zero_()            # mean exactly 0 at init: the arm is FNPEG (gate RI)
            a.net[-1].bias.zero_()
            a.log_std.fill_(math.log(float(arm_kw.get("init_std", 0.1))))
        return a
    if arm == "liu_track":
        return RL.Actor(obs_dim=LIU_OBS_DIM, hidden=cfg.hidden, kind="gauss")
    return RL.Actor(hidden=cfg.hidden, kind="gauss")


def make_env(arm: str, cfg: RL.RLCfg, seed: int, device, fnpeg_cfg: Optional[dict], arm_kw: dict,
             n_envs: Optional[int] = None, pool: Optional[int] = None, train: bool = True,
             conds: Optional[Sequence[dict]] = None):
    fc = RF.FastEnvCfg(n_envs=int(n_envs or cfg.n_envs), gamma_rl=cfg.gamma, seed=int(seed),
                       pool=int(pool or cfg.env_pool))
    if arm == "resid_fnpeg":
        return ResidEnv(fc, device, cfg, fnpeg_cfg, float(arm_kw.get("bound_deg", 20.0)), int(cfg.label_workers),
                        conds)
    if arm == "liu_track":
        return LiuEnv(fc, device, liu_cfg_from(fnpeg_cfg, ds_frac=float(arm_kw.get("ds_frac", 0.3)),
                                               train_terms=bool(train)), conds)
    env = RF.FastEntryEnv(fc, device)
    if conds:
        set_pool(env, conds)
    env.obs_dim = RF.OBS_DIM
    return env


def _close(env) -> None:
    if hasattr(env, "close"):
        env.close()


# ============================================================================
# PPO on any of the environments (rl_v5.ppo_train's loop)
# ============================================================================
def ppo_train_env(L, d: Path, cfg: RL.RLCfg, device, label: str, arm: str, fnpeg_cfg: Optional[dict], arm_kw: dict,
                  log=print, stage: str = "rl2_train") -> dict:
    """PPO for a notebook-2 arm. Resumable from d/train_state.pt (the environment's state, FNPEG's flights included);
    writes d/policy.pt with the arm and its settings in the policy's extra."""
    dev = torch.device(device)
    torch.manual_seed(cfg.seed)
    actor = make_actor(arm, cfg, arm_kw).to(dev)
    env = make_env(arm, cfg, cfg.seed, dev, fnpeg_cfg, arm_kw)
    try:
        vf = RL.critic(obs_dim=env.obs_dim, hidden=cfg.hidden).to(dev)
        opt = torch.optim.Adam(list(actor.parameters()) + list(vf.parameters()), lr=cfg.lr, eps=1e-5)
        gen = torch.Generator(device=dev).manual_seed(cfg.seed * 7 + 3)
        curve, step, update, elapsed = RL.Curve(), 0, 0, 0.0
        st = load_state(Path(d) / "train_state.pt", log)
        if st is not None:
            actor.load_state_dict(st["actor"])
            vf.load_state_dict(st["critic"])
            opt.load_state_dict(st["opt"])
            env.set_state(st["env"])
            step, update, elapsed, curve = st["step"], st["update"], st["elapsed"], RL.Curve(st["curve"])
            set_rng_state(st["rng"], {"gen": gen})
            log(f"  [{label}] resumed at update {update} ({step:,} decisions)")
        obs = env.obs()
        N, T = int(env.cfg.n_envs), int(cfg.n_steps)
        t_save, t0 = time.time(), time.time()
        snaps = RL._snapshot_marks(cfg)
        checks = [f for f in RL.CHECK_AT if RL._used(step, elapsed, cfg) < f]
        meta = (arm, dict(arm_kw), fnpeg_cfg)

        def save():
            save_state(L, {"actor": actor.state_dict(), "critic": vf.state_dict(), "opt": opt.state_dict(),
                           "env": env.get_state(), "step": step, "update": update, "elapsed": elapsed,
                           "curve": curve.state(), "rng": rng_state({"gen": gen}), "cfg": asdict(cfg),
                           "arm": arm, "arm_kw": dict(arm_kw)}, Path(d) / "train_state.pt")
            RL._write_curve(L, Path(d), curve, label)

        while step < cfg.budget_decisions and elapsed + (time.time() - t0) < cfg.budget_hours * 3600:
            warm = update < cfg.critic_warmup_updates
            buf_o = torch.zeros(T, N, env.obs_dim, device=dev)
            buf_raw = torch.zeros(T, N, 2, device=dev)
            buf_lp, buf_r, buf_d, buf_v = (torch.zeros(T, N, device=dev) for _ in range(4))
            with torch.no_grad():
                for k in range(T):
                    actor.norm.update(obs)              # also during the critic burn-in (the critic learns on it)
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
            RL._ppo_update(actor, vf, opt, cfg, gen, None, warm, buf_o.flatten(0, 1), buf_raw.flatten(0, 1),
                           buf_lp.flatten(0, 1), adv.flatten(), ret.flatten())
            step += N * T
            update += 1
            curve.flush(step, update, elapsed + time.time() - t0, warmup=int(warm))
            if update % 10 == 0 or update == 1:
                rw = curve.rows[-1]
                log(f"  [{label}] update {update} {step:,} decisions: return {rw['ep_return']:.3f} "
                    f"success@5km {rw['success_5km_fast']:.3f} ({(elapsed + time.time() - t0) / 60:.1f} min)")
            _maybe_snapshot(L, Path(d), actor, step, cfg, snaps, elapsed + time.time() - t0, meta)
            while checks and RL._used(step, elapsed + time.time() - t0, cfg) >= checks[0]:
                why = RL.training_check(actor, curve.rows, checks.pop(0))
                if why:
                    elapsed += time.time() - t0
                    t0 = time.time()
                    save()
                    raise RL.RLStopped(f"[{label}] {why}", step)
            due = time.time() - t_save > cfg.save_every_s or (cfg.save_every_updates and
                                                              update % cfg.save_every_updates == 0)
            if due:
                elapsed += time.time() - t0
                t0 = time.time()
                save()
                t_save = time.time()
                fail_point(stage, update)
        elapsed += time.time() - t0
        t0 = time.time()
        _maybe_snapshot(L, Path(d), actor, step, cfg, snaps, elapsed, meta, final=True)
        save()
        save_policy(L, Path(d) / "policy.pt", actor, arm, arm_kw, fnpeg_cfg,
                    {"label": label, "decisions": step, "updates": update, "train_s": round(elapsed, 1),
                     "cfg": asdict(cfg)})
        return {"decisions": step, "updates": update, "train_s": round(elapsed, 1)}
    finally:
        _close(env)


def _maybe_snapshot(L, d: Path, actor: RL.Actor, step: int, cfg: RL.RLCfg, marks, elapsed_s: float, meta: tuple,
                    final: bool = False) -> None:
    """rl_v5._maybe_snapshot, with the arm and its settings in each snapshot (a snapshot flies as its arm)."""
    for i, m in enumerate(marks):
        p = d / "snapshots" / f"policy_{i + 1}of{len(marks)}.pt"
        if p.exists():
            continue
        by_dec = step >= int(cfg.budget_decisions * m)
        by_time = elapsed_s >= cfg.budget_hours * 3600 * m
        why = "decisions" if by_dec else ("wall_clock" if by_time else ("final" if final else None))
        if why:
            save_policy(L, p, actor, *meta, {"decisions": step, "point": i + 1, "trigger": why,
                                             "elapsed_s": round(float(elapsed_s), 1)})


def save_policy(L, path, actor: RL.Actor, arm: str, arm_kw: dict, fnpeg_cfg: Optional[dict], extra: dict) -> None:
    RL.save_policy(L, path, actor, {**extra, "arm": arm, "arm_kw": dict(arm_kw), "fnpeg_cfg": dict(fnpeg_cfg or {})})


def policy_meta(path) -> dict:
    ck = torch.load(path, map_location="cpu", weights_only=False)
    return dict(ck.get("extra") or {})


def train_dir(L, d: Path, arm: str, cfg: RL.RLCfg, device, fnpeg_cfg: Optional[dict], arm_kw: dict, log=print,
              label: Optional[str] = None) -> dict:
    """Train (or resume, or skip) one run of a notebook-2 arm into d; policy.pt marks it done, run.json describes it."""
    d = Path(d)
    label = label or f"{arm} seed {cfg.seed}"
    if (d / "policy.pt").exists():
        RL.ensure_curves(L, d, label)
        p = d / "run.json"
        return json.loads(p.read_text()) if p.exists() else {"arm": arm, "seed": cfg.seed}
    L.guard(d).mkdir(parents=True, exist_ok=True)
    hw = RL.check_hardware(L, d, device, log)
    t0 = time.time()
    if arm == "ppo_hist":
        out = RL.ppo_train(L, d, cfg, device, label, log=log, stage="rl2_train")     # preregistered, unchanged
    else:
        out = ppo_train_env(L, d, cfg, device, label, arm, fnpeg_cfg, arm_kw, log)
    info = {"arm": arm, "seed": int(cfg.seed), **out, "wall_s": round(time.time() - t0, 1), "cfg": asdict(cfg),
            "arm_kw": dict(arm_kw), "hardware": hw}
    L.write_json(d / "run.json", info)
    return info


# ============================================================================
# Evaluation on fixed episodes (the gates' "same episodes")
# ============================================================================
def load_actor(path, device="cpu") -> RL.Actor:
    return RL.load_policy(path, device)


@torch.no_grad()
def eval_fixed(arm: str, actor: Optional[RL.Actor], cfg: RL.RLCfg, fnpeg_cfg: Optional[dict], arm_kw: dict,
               n: int = 512, seed: int = 12345, device="cpu", conds: Optional[Sequence[dict]] = None) -> dict:
    """success@5km on n fixed fast-env episodes: n environment rows, each counted at its first episode's end (the
    first n valid draws of the env's seed, whatever flies them). arm "fnpeg": FNPEG alone (the residual env with a
    zero correction). Training terminations are off (an evaluation ends where physics ends it). conds: the
    episodes are these conditions instead (n = len(conds); the fast side of R-gap)."""
    if conds:
        n = len(conds)
    kind = "resid_fnpeg" if arm == "fnpeg" else arm
    kw = dict(arm_kw, bound_deg=0.0) if arm == "fnpeg" else arm_kw
    env = make_env(kind, cfg, seed, device, fnpeg_cfg, kw, n_envs=n, pool=n, train=False, conds=conds)
    first = np.zeros(n, dtype=bool)
    succ, miss = np.zeros(n, dtype=bool), np.full(n, np.nan)
    steps = 0
    try:
        obs = env.obs()
        while not first.all():
            a = torch.zeros(n, device=env.dev) if actor is None else actor.deterministic(obs.to(device))
            if kind == "resid_fnpeg":
                obs, _, done, info = env.step(a, active=~first)
            else:
                obs, _, done, info = env.step(a)
            steps += 1
            dn = done.cpu().numpy() & ~first
            if dn.any():
                succ[dn] = info["success_5km"].cpu().numpy()[dn]
                miss[dn] = info["miss_km"].cpu().numpy()[dn]
                first |= dn
            if steps > 10_000:
                raise RuntimeError("eval_fixed: episodes did not end")
    finally:
        _close(env)
    return {"arm": arm, "n": int(n), "seed": int(seed), "success_5km_fast": float(succ.mean()),
            "miss_km_fast_median": float(np.nanmedian(miss)), "success_by_episode": succ.astype(int).tolist()}


# ============================================================================
# Gates (plan section 3)
# ============================================================================
RP_RESID_MARGIN = 0.02            # resid_fnpeg: success >= FNPEG's - 2 pp
RP_FRAC = 0.5                     # liu_track, ppo_hist: success >= 50 % of FNPEG's


def _rises(rows: Sequence[dict]) -> Tuple[bool, Optional[float], Optional[float], bool]:
    ret = [float(r["ep_return"]) for r in rows if r.get("episodes")]
    finite = all(math.isfinite(x) for x in ret)
    if len(ret) < 3:
        return False, (float(np.mean(ret)) if ret else None), None, finite
    k = max(1, len(ret) // 3)
    a, b = float(np.mean(ret[:k])), float(np.mean(ret[-k:]))
    return b > a, a, b, finite


def rp_verdict(arm: str, rows: Sequence[dict], actor: torch.nn.Module, succ: float, ref: float) -> dict:
    """Gate RP from a pilot's learning curve and its success on the fixed episodes against FNPEG's."""
    rise, r0, r1, finite = _rises(rows)
    finite = finite and all(bool(torch.isfinite(p).all()) for p in actor.parameters())
    best = max([float(r.get("success_5km_fast", 0.0)) for r in rows] or [0.0])
    if arm == "resid_fnpeg":
        bar = ref - RP_RESID_MARGIN
        ok = finite and rise and succ >= bar
        need = f"need >= FNPEG's - {100 * RP_RESID_MARGIN:.0f} pp = {bar:.3f}, a rising return, no NaN"
    else:
        bar = RP_FRAC * ref
        ok = finite and succ >= bar and ref > 0
        need = f"need >= {RP_FRAC:.0%} of FNPEG's = {bar:.3f}, no NaN"
        if arm == "ppo_hist":                        # the old RV3 as well
            ok = ok and rise and best > 0
            need += ", a rising return and success > 0 during training (old RV3)"
    ret = (f"return {r0:.3f} -> {r1:.3f} ({'rises' if rise else 'does not rise'})" if r1 is not None
           else "too few learning-curve rows for a trend")
    return {"passed": bool(ok), "success_5km_fixed": float(succ), "fnpeg_success_5km_fixed": float(ref),
            "bar": float(bar), "finite": bool(finite), "return_first": r0, "return_last": r1,
            "best_success_5km_training": best, "curve_rows": len(rows),
            "reason": f"pilot {succ:.3f} success@5km on the fixed episodes, FNPEG {ref:.3f} ({need}); {ret}"}


@torch.no_grad()
def ri_fast(cfg: RL.RLCfg, fnpeg_cfg: Optional[dict], arm_kw: dict, n: int = 64, steps: int = 400, seed: int = 4242,
            device="cpu") -> dict:
    """Gate RI part 1: the initialized residual actor's mean is exactly 0 on the fast env's observations along whole
    episodes, and the composed command equals FNPEG's bit for bit."""
    actor = make_actor("resid_fnpeg", cfg, arm_kw).to(device)
    env = make_env("resid_fnpeg", cfg, seed, device, fnpeg_cfg, arm_kw, n_envs=n, pool=max(64, n))
    bad_a, bad_c, k, eps = 0, 0, 0, 0
    try:
        obs = env.obs()
        for k in range(steps):
            actor.norm.update(obs)
            a = actor.deterministic(obs.to(device))
            q = torch.as_tensor(env.q[:, 0], dtype=torch.float64)
            bad_a += int((a != 0).sum())
            cmd = env.command(a).cpu()
            bad_c += int((cmd != q).sum())
            obs, _, done, _ = env.step(a)
            eps += int(done.sum())
    finally:
        _close(env)
    return {"passed": bool(bad_a == 0 and bad_c == 0), "decisions": int(n * (k + 1)), "episodes_finished": eps,
            "nonzero_actions": bad_a, "commands_differing": bad_c}


def flights_identical(a: dict, b: dict) -> Tuple[bool, List[str]]:
    """Two flight records (result + telemetry) equal: every telemetry array of fnpeg's exactly, and miss, outcome,
    peak g. Returns (same, what differs)."""
    diff = []
    for k in ("miss_km", "outcome", "peak_g", "flight_time_s"):
        if k in a["result"] or k in b["result"]:
            va, vb = a["result"].get(k), b["result"].get(k)
            if not (va == vb or (isinstance(va, float) and isinstance(vb, float) and math.isnan(va) and math.isnan(vb))):
                diff.append(k)
    ta, tb = a["telemetry"], b["telemetry"]
    for k in ta:
        if k not in tb:
            diff.append(f"tel:{k} missing")
        elif not np.array_equal(np.asarray(ta[k]), np.asarray(tb[k]), equal_nan=True):
            diff.append(f"tel:{k}")
    return (not diff), diff


# ============================================================================
# Cost per million decisions and routing (plan section 3)
# ============================================================================
def units_per_million(decisions_per_s: float, units_per_hour: float) -> float:
    return 1e6 / max(float(decisions_per_s), 1e-9) / 3600.0 * float(units_per_hour)


def route(rates: Dict[str, dict]) -> Optional[dict]:
    """The cheaper runtime from {runtime: {"decisions_per_s", "units_per_hour"}}; None until "cpu" and "L4" exist."""
    if not all(k in rates for k in ("cpu", "L4")):
        return None
    cost = {k: units_per_million(v["decisions_per_s"], v["units_per_hour"]) for k, v in rates.items()}
    best = min(("cpu", "L4"), key=lambda k: cost[k])
    return {"runtime": best, "units_per_million": cost, "rates": rates}


# ============================================================================
# The arms in the full simulator (fly_v5's 1 Hz runner)
# ============================================================================
class ResidLaw:
    """resid_fnpeg in the full simulator: FNPEG every second exactly as the fnpeg arm, the policy's correction b a
    every 5 s (held), sigma = s_F clip(|sigma_F| + b a, 0, 90 deg)."""
    period = 5

    def __init__(self, path: Optional[str], fnpeg_cfg: dict, bound_deg: Optional[float] = None, device="cpu",
                 actor: Optional[RL.Actor] = None):
        from jepa.fly_v5 import FNPEGLaw
        from jepa.guidance_fnpeg import FNPEGConfig
        meta = policy_meta(path) if path else {}
        self.fn = FNPEGLaw(FNPEGConfig(**fnpeg_cfg), device)
        self.actor = actor if actor is not None else load_actor(path, device)
        kw = meta.get("arm_kw") or {}
        self.bound = math.radians(float(bound_deg if bound_deg is not None else kw.get("bound_deg", 20.0)))
        self.dev = torch.device(device)

    def new_flight(self, job, env):
        from types import SimpleNamespace
        return SimpleNamespace(fn=self.fn.new_flight(job, env), obs=RF.FullObs(), delta=0.0, decided=False, n_dec=0,
                               deltas=[])

    def decide(self, states, ys, t_s, envs, ks):
        cmds = self.fn.decide([s.fn for s in states], ys, t_s, envs, ks)
        idx = [i for i, k in enumerate(ks) if int(k) % self.period == 0]
        if idx:
            O = []
            for i in idx:
                f = states[i].fn
                q = np.array([[cmds[i], f.K_D, f.K_L, sign_eff(float(cmds[i]), f.sign)]])
                O.append(torch.cat([states[i].obs(envs[i], ys[i]), torch.as_tensor(resid_features(q)[0]).float()]))
            a = self.actor.deterministic(torch.stack(O).to(self.dev)).double().cpu().numpy()
            for j, i in enumerate(idx):
                states[i].delta = float(self.bound * float(np.clip(a[j], -1.0, 1.0)))
                states[i].deltas.append(states[i].delta)
                states[i].n_dec += 1
        out = []
        for st, c in zip(states, cmds):
            st.decided = bool(st.fn.decided)
            out.append(compose(float(c), sign_eff(float(c), st.fn.sign), st.delta))
        return out

    def flight_summary(self, st) -> dict:
        d = np.degrees(np.asarray(st.deltas, dtype=float))
        return {**self.fn.flight_summary(st.fn), "policy_decisions": int(st.n_dec), "policy_period_s": self.period,
                "resid_mean_abs_deg": float(np.abs(d).mean()) if len(d) else 0.0,
                "resid_max_abs_deg": float(np.abs(d).max()) if len(d) else 0.0}

    def decisions(self, st) -> Dict[str, np.ndarray]:
        return self.fn.decisions(st.fn)


class LiuLaw:
    """liu_track in the full simulator: every 5 s the policy's |sigma| and Eq. 19's sign (held between decisions)."""
    period = 5

    def __init__(self, path: str, fnpeg_cfg: dict, device="cpu"):
        meta = policy_meta(path)
        self.actor = load_actor(path, device)
        self.lc = liu_cfg_from(fnpeg_cfg, ds_frac=float((meta.get("arm_kw") or {}).get("ds_frac", 0.3)),
                               train_terms=False)
        self.dev = torch.device(device)

    def new_flight(self, job, env):
        from types import SimpleNamespace
        return SimpleNamespace(obs=RF.FullObs(), cmd=0.0, decided=False, n_dec=0, s0=None, v0=None, sign=0.0, nrev=0.0,
                               tphi=float(env.target_phi_rad), tlam=float(env.target_lam_rad))

    def decide(self, states, ys, t_s, envs, ks):
        idx = [i for i, k in enumerate(ks) if int(k) % self.period == 0]
        for st in states:
            st.decided = False
        for i in idx:
            st, y = states[i], np.asarray(ys[i], dtype=np.float64)
            x = torch.tensor([[y[0] + R_EARTH, y[1], y[2], y[3], y[4], y[5]]], dtype=torch.float64)
            pt, lt = torch.tensor([st.tphi], dtype=torch.float64), torch.tensor([st.tlam], dtype=torch.float64)
            if st.s0 is None:
                ang, _ = RF.range_course(x[:, 1], x[:, 2], pt, lt)
                st.s0, st.v0 = float(ang[0]) * R_EARTH, float(y[3])
            s0, v0 = torch.tensor([st.s0], dtype=torch.float64), torch.tensor([st.v0], dtype=torch.float64)
            T = liu_terms(x, pt, lt, s0, v0, self.lc)
            n_g = torch.tensor([math.hypot(float(y[8]), float(y[9])) / G0], dtype=torch.float64)
            sg = torch.tensor([st.sign], dtype=torch.float64)
            f = liu_features(T, n_g, s0, sg, torch.tensor([st.nrev], dtype=torch.float64), self.lc).float()[0]
            o = torch.cat([st.obs(envs[i], ys[i]), f])
            a = float(self.actor.deterministic(o[None].to(self.dev))[0])
            new = float(liu_sign(sg, T["dpsi"], T["db"])[0])
            if st.sign != 0.0 and new != st.sign:
                st.nrev += 1
            st.sign = new
            st.cmd = st.sign * 0.5 * (min(1.0, max(-1.0, a)) + 1.0) * SIGMA_MAX_RAD
            st.decided = True
            st.n_dec += 1
        return [st.cmd for st in states]

    def flight_summary(self, st) -> dict:
        return {"policy_decisions": int(st.n_dec), "policy_period_s": self.period, "liu_reversals": int(st.nrev)}

    def decisions(self, st) -> Dict[str, np.ndarray]:
        return {}


def make_law(model: dict, device="cpu"):
    """The law of a notebook-2 flight job (model["rl2"]), or rl_v5's PolicyLaw for ppo_hist."""
    from jepa import fly_v5 as F5
    r2 = model.get("rl2") or {}
    arm = r2.get("arm")
    if arm == "resid_fnpeg":
        return ResidLaw(model.get("path"), r2["fnpeg_cfg"], r2.get("bound_deg"), device)
    if arm == "liu_track":
        return LiuLaw(model["path"], r2["fnpeg_cfg"], device)
    if arm == "fnpeg":
        return F5.make_law("fnpeg", r2["fnpeg_cfg"], device)
    return F5.make_law("policy", {"path": model["path"], "legacy": False}, device)


class ZeroResidLaw(ResidLaw):
    """resid_fnpeg at init (gate RI part 2): a freshly initialized residual actor, no policy file."""

    def __init__(self, fnpeg_cfg: dict, cfg: RL.RLCfg, arm_kw: dict, device="cpu"):
        super().__init__(None, fnpeg_cfg, float(arm_kw.get("bound_deg", 20.0)), device,
                         actor=make_actor("resid_fnpeg", cfg, arm_kw).eval())
