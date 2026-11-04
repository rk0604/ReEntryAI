"""
jepa/rlenv_fast_v5.py

The batched fast training environment of the v5 RL baselines (jepa/v5_plan.md section 4.4).

    dynamics     physics_v4's model with each environment's true density, L/D and mass (as mppi_truth),
                 float64: 5 s decisions in (2, 2, 1) s explicit Euler substeps, the bank moving toward the
                 command at 20 deg/s (rollout_physics' inner loop, bit for bit: smoke check 6)
    episodes     entry state from the world models' training envelope (data.EnvelopeSpec), truth parameters
                 from data.WIDE_DOMAIN, the target on the case's own footprint (downrange at a uniform fraction
                 of the way from the 90 deg to the 0 deg constant-bank range of that entry and parameters,
                 crossrange uniform in +-60 km); terminates at the drogue altitude (crossing interpolated),
                 |cos gamma| < 0.05, skip-out (above 160 km) or 1,400 s; finished environments reset in place
    reward       minus MPPI's own cost, realized (mppi.MPPI._cost): at the end -miss_km, -3,000 km if the
                 flight never crossed the drogue altitude; per decision -2 km/rad |d sigma|, -50 km/g^2
                 max(0, g - 4.5)^2 (the step's peak load), -2 km/deg^2 max(0, |gamma| - 82 deg)^2; plus the
                 potential-based shaping gamma_rl Phi(s') - Phi(s), Phi = -range to go (km), Phi = 0 at the
                 end; everything x reward_scale
    observation  rl_env's 13 features (computed here in torch; the heat-rate margin uses Sutton-Graves
                 convective heating only, the full simulator's interval upper bound adds a radiative term) +
                 the last 24 decisions (120 s) of [h, V, gamma, a_D, cD_meas, cL_meas, sigma], normalized;
                 obs_full() builds the same vector in the full simulator from rl_env and data.physical_state
    expert_rows  the physical_state rows (h, phi, lam, V, gamma, chi, sigma, roll rate, a_D, a_L) FNPEG reads
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import torch

import constants
from jepa import physics_v4 as PH
from jepa.data import G0, R_EARTH, SIGMA_MAX_RAD, WIDE_DOMAIN, EnvelopeSpec, sample_domain, sample_entry_ic

HIST_LEN = 24
HIST_NAMES = ("h", "V", "gamma", "a_D", "cD_meas", "cL_meas", "sigma")
HIST_SCALE = (120_000.0, 11_000.0, math.pi / 2, 5.0 * G0, 1.0, 1.0, math.pi / 2)
N_BASE = 13
OBS_DIM = N_BASE + HIST_LEN * len(HIST_NAMES)
NORM = {"alt": 120_000.0, "vel": 11_000.0, "gamma": math.pi / 2.0, "range": 2_500_000.0, "heading": math.pi,
        "cross": 500_000.0, "energy": 6.0e7, "sigma": math.pi / 2.0, "hdot": 300.0}    # rl_env._NORM
SKIP_ALT_M = 160_000.0


@dataclass
class FastEnvCfg:
    n_envs: int = 4096
    dt_s: float = 5.0
    h_drogue_m: float = 7600.0
    t_max_s: float = 1400.0
    g_limit: float = 4.5
    w_smooth_km_per_rad: float = 2.0
    w_g_km: float = 50.0
    gamma_limit_deg: float = 82.0
    w_gamma_km_per_deg2: float = 2.0
    nocross_km: float = 3000.0
    reward_scale: float = 0.01
    shaping: bool = True
    gamma_rl: float = 0.99
    target_frac: Tuple[float, float] = (0.15, 0.85)
    crossrange_km: float = 60.0
    gload_limit_obs_g: float = 8.0          # rl_env's default gload_limit_g (the observation's g margin)
    pool: int = 2048                        # episodes drawn (and their footprints flown) per refill
    seed: int = 0


def hist_row(h, V, gam, aD, aL, sig) -> torch.Tensor:
    """One normalized history row [.., 7] (models_v2.raw_features_v2's cD_meas, cL_meas)."""
    from jepa.models_v2 import BETA_NOM, H_SCALE_M, LD_NOM, RHO0_KGM3
    aD_ = aD.clamp_min(1e-6)
    aL_ = aL.clamp_min(1e-6)
    aD_nom = 0.5 * RHO0_KGM3 * torch.exp(-h.clamp_min(0.0) / H_SCALE_M) * V * V / BETA_NOM
    cD = torch.log(aD_ / aD_nom.clamp_min(1e-12)).clamp(-6.0, 6.0)
    cL = torch.log(aL_ / aD_ / LD_NOM).clamp(-6.0, 6.0)
    sc = torch.as_tensor(HIST_SCALE, dtype=h.dtype, device=h.device)
    return torch.stack([h, V, gam, aD, cD, cL, sig], -1) / sc


def range_course(phi, lam, phi_t, lam_t):
    """Central angle to the target and its course in the chi convention (from east toward north)."""
    dlam = lam_t - lam
    a = torch.sin((phi_t - phi) / 2) ** 2 + torch.cos(phi) * torch.cos(phi_t) * torch.sin(dlam / 2) ** 2
    ang = 2 * torch.asin(torch.sqrt(a.clamp(0.0, 1.0)))
    y = torch.sin(dlam) * torch.cos(phi_t)
    x = torch.cos(phi) * torch.sin(phi_t) - torch.sin(phi) * torch.cos(phi_t) * torch.cos(dlam)
    bearing = torch.atan2(y, x)
    return ang, math.pi / 2 - bearing


def wrap(x):
    return torch.atan2(torch.sin(x), torch.cos(x))


def base_obs(x, sig, gload, heat_limit: float, g_lim: float, phi_t, lam_t, atm) -> torch.Tensor:
    """rl_env._build_obs' 13 features for a batch (x = r, phi, lam, V, gamma, chi)."""
    r, phi, lam, V, gam, chi = x.unbind(-1)
    alt = (r - R_EARTH).clamp_min(0.0)
    rho, _ = atm.state(alt)
    qdot = constants.SUTTON_GRAVES_K_EARTH * torch.sqrt(rho / constants.HEAT_SHIELD_NOSE_RADIUS_M) * V ** 3
    ang, course = range_course(phi, lam, phi_t, lam_t)
    rng = ang * R_EARTH
    herr = torch.where(rng > 1e-9, wrap(course - chi), torch.zeros_like(chi))
    spec_e = 0.5 * V * V + G0 * (r - R_EARTH)
    return torch.stack([(r - R_EARTH) / NORM["alt"], V / NORM["vel"], gam / NORM["gamma"], torch.sin(chi),
                        torch.cos(chi), rng / NORM["range"], herr / NORM["heading"],
                        rng * torch.sin(herr) / NORM["cross"], spec_e / NORM["energy"], sig / NORM["sigma"],
                        (heat_limit - qdot) / heat_limit, (g_lim - gload) / g_lim, V * torch.sin(gam) / NORM["hdot"]],
                       -1)


def aero_accels(x, sig, S_ref, mass, ls, atm):
    """(a_D, a_L, load in g) of the true model at x."""
    alt = (x[:, 0] - R_EARTH).clamp_min(0.0)
    rho, T = atm.state(alt)
    q = 0.5 * rho * x[:, 3] ** 2
    CD, CL = PH.aero_coefficients(x[:, 3], T, ls)
    aD, aL = q * S_ref * CD / mass, q * S_ref * CL / mass
    return aD, aL, torch.sqrt(aD * aD + aL * aL) / G0


def target_from_local(phi0, lam0, chi0, dr_km, cr_km):
    """Entry-frame (downrange, crossrange) -> target latitude/longitude (testset_v3.entry_frame inverse)."""
    R = R_EARTH / 1000.0
    E = torch.stack([torch.cos(phi0) * torch.cos(lam0), torch.cos(phi0) * torch.sin(lam0), torch.sin(phi0)], -1)
    east = torch.stack([-torch.sin(lam0), torch.cos(lam0), torch.zeros_like(lam0)], -1)
    north = torch.cross(E, east, dim=-1)
    h = torch.cos(chi0)[:, None] * east + torch.sin(chi0)[:, None] * north
    n = torch.cross(E, h, dim=-1)
    a, c = dr_km / R, cr_km / R
    p = (torch.cos(c) * torch.cos(a))[:, None] * E + (torch.cos(c) * torch.sin(a))[:, None] * h + torch.sin(c)[:, None] * n
    return torch.asin(p[:, 2].clamp(-1, 1)), torch.atan2(p[:, 1], p[:, 0])


class FastEntryEnv:
    """N environments in lockstep on one device. step(actions in [-1, 1]) -> (obs, reward, done, info)."""

    def __init__(self, cfg: Optional[FastEnvCfg] = None, device="cpu"):
        self.cfg = cfg or FastEnvCfg()
        self.dev = torch.device(device)
        self.atm = PH.atm_table(self.dev)
        p = PH.nominal_params()
        self.area0, self.mass0 = float(p["ref_area_m2"]), float(p["mass_kg"])
        self.heat_limit = float(constants.HEAT_RATE_LIMIT_DEFAULT)
        self.rng = np.random.default_rng(self.cfg.seed)
        N = self.cfg.n_envs
        f64 = dict(dtype=torch.float64, device=self.dev)
        self.x = torch.zeros(N, 6, **f64)
        self.sig = torch.zeros(N, **f64)
        self.u_prev = torch.zeros(N, **f64)
        self.S_ref, self.mass, self.ls = (torch.ones(N, **f64) for _ in range(3))
        self.phi_t, self.lam_t = torch.zeros(N, **f64), torch.zeros(N, **f64)
        self.t = torch.zeros(N, **f64)
        self.hist = torch.zeros(N, HIST_LEN, len(HIST_NAMES), **f64)
        self.gload = torch.zeros(N, **f64)
        self.domain = np.zeros((N, 3))
        self.ep_return = torch.zeros(N, **f64)
        self.ep_len = torch.zeros(N, dtype=torch.long, device=self.dev)
        self._pool: Dict[str, torch.Tensor] = {}
        self._pool_i = 0
        self.reset_rows(torch.arange(N, device=self.dev))

    # ---- RNG (checkpoint) ----
    def get_state(self) -> dict:
        return {"rng": self.rng.bit_generator.state, "pool": {k: v.clone().cpu() for k, v in self._pool.items()},
                "pool_i": self._pool_i, **{k: getattr(self, k).clone().cpu() for k in
                ("x", "sig", "u_prev", "S_ref", "mass", "ls", "phi_t", "lam_t", "t", "hist", "gload",
                 "ep_return", "ep_len")}, "domain": self.domain.copy()}

    def set_state(self, st: dict) -> None:
        self.rng.bit_generator.state = st["rng"]
        self._pool = {k: v.to(self.dev) for k, v in st["pool"].items()}
        self._pool_i = int(st["pool_i"])
        for k in ("x", "sig", "u_prev", "S_ref", "mass", "ls", "phi_t", "lam_t", "t", "hist", "gload", "ep_return",
                  "ep_len"):
            setattr(self, k, st[k].to(self.dev))
        self.domain = st["domain"].copy()

    # ---- episodes ----
    def _draw(self, n: int):
        ics = [sample_entry_ic(self.rng, EnvelopeSpec()) for _ in range(n)]
        doms = [sample_domain(self.rng, WIDE_DOMAIN) for _ in range(n)]
        fr = self.rng.uniform(*self.cfg.target_frac, size=n)
        cr = self.rng.uniform(-self.cfg.crossrange_km, self.cfg.crossrange_km, size=n)
        return ics, doms, fr, cr

    def _x_from_ic(self, ics) -> torch.Tensor:
        a = np.array([[ic["altitude_m"] + R_EARTH, math.radians(ic["phi_deg"]), math.radians(ic["lambda_deg"]),
                       ic["velocity_mps"], math.radians(ic["gamma_deg"]), math.radians(ic["chi_deg"])] for ic in ics])
        return torch.as_tensor(a, dtype=torch.float64, device=self.dev)

    def footprint_ranges(self, x0, S_ref, mass, ls) -> Tuple[torch.Tensor, torch.Tensor]:
        """Downrange (km, along the entry heading) of the 90 deg and 0 deg constant-bank flights."""
        n = len(x0)
        y0 = torch.cat([x0, x0]).clone()
        y0[:, 0] -= R_EARTH
        u = torch.cat([torch.full((n,), SIGMA_MAX_RAD, dtype=torch.float64, device=self.dev),
                       torch.zeros(n, dtype=torch.float64, device=self.dev)])
        K = int(self.cfg.t_max_s / self.cfg.dt_s)
        y6, _ = PH.rollout_physics(y0, u, u[:, None].expand(2 * n, K), torch.cat([S_ref, S_ref]),
                                   torch.cat([mass, mass]), torch.cat([ls, ls]), self.cfg.h_drogue_m, self.atm)
        from jepa.evaluate import crossing_point
        phi, lam, crossed, _ = crossing_point(y6, self.cfg.h_drogue_m)
        ang, _ = range_course(torch.cat([x0[:, 1], x0[:, 1]]), torch.cat([x0[:, 2], x0[:, 2]]), phi, lam)
        d = ang * R_EARTH / 1000.0
        # a flight that turned near vertical (|cos gamma| < 0.05) low in the atmosphere lands where it stopped
        # (as benchmark_v5's S0 rule); one that ends high (skip-out, timeout) has no landing point
        low = y6[:, -1, 0] < 30_000.0
        d = torch.where(crossed | low, d, torch.full_like(d, float("nan")))
        return d[:n], d[n:]

    def _refill(self) -> None:
        """Draw cfg.pool episodes (entry, truth, target on the footprint); keep the ones with a footprint."""
        f64 = dict(dtype=torch.float64, device=self.dev)
        n = int(self.cfg.pool)
        ics, doms, fr, cr = self._draw(n)
        x0 = self._x_from_ic(ics)
        dom = torch.as_tensor([[d["density_scale"], d["lift_scale"], d["mass_scale"]] for d in doms], **f64)
        S_ref, mass, ls = self.area0 * dom[:, 0], self.mass0 * dom[:, 2], dom[:, 1]
        d90, d0 = self.footprint_ranges(x0, S_ref, mass, ls)
        ok = torch.isfinite(d90) & torch.isfinite(d0) & (d0 > d90 + 20.0)
        dr = d90 + torch.as_tensor(fr, **f64) * (d0 - d90)
        phi_t, lam_t = target_from_local(x0[:, 1], x0[:, 2], x0[:, 5], torch.nan_to_num(dr), torch.as_tensor(cr, **f64))
        self._pool = {"x0": x0[ok], "S_ref": S_ref[ok], "mass": mass[ok], "ls": ls[ok], "phi_t": phi_t[ok],
                      "lam_t": lam_t[ok], "dom": dom[ok]}
        self._pool_i = 0
        if not bool(ok.any()):
            raise RuntimeError("fast env: no entry of the pool has a footprint")

    def reset_rows(self, rows: torch.Tensor) -> None:
        rows = rows.to(self.dev)
        done = 0
        while done < len(rows):
            if not self._pool or self._pool_i >= len(self._pool["x0"]):
                self._refill()
            take = min(len(rows) - done, len(self._pool["x0"]) - self._pool_i)
            r = rows[done:done + take]
            sl = slice(self._pool_i, self._pool_i + take)
            P = self._pool
            self.x[r], self.S_ref[r], self.mass[r], self.ls[r] = P["x0"][sl], P["S_ref"][sl], P["mass"][sl], P["ls"][sl]
            self.phi_t[r], self.lam_t[r] = P["phi_t"][sl], P["lam_t"][sl]
            self.domain[r.cpu().numpy()] = P["dom"][sl].cpu().numpy()
            self._pool_i += take
            done += take
        self.sig[rows] = 0.0
        self.u_prev[rows] = 0.0
        self.t[rows] = 0.0
        self.ep_return[rows] = 0.0
        self.ep_len[rows] = 0
        aD, aL, g = aero_accels(self.x[rows], self.sig[rows], self.S_ref[rows], self.mass[rows], self.ls[rows], self.atm)
        self.gload[rows] = g
        h = self.x[rows, 0] - R_EARTH
        row = hist_row(h, self.x[rows, 3], self.x[rows, 4], aD, aL, self.sig[rows])
        self.hist[rows] = row[:, None].expand(-1, HIST_LEN, -1)

    # ---- observation ----
    def obs(self) -> torch.Tensor:
        b = base_obs(self.x, self.sig, self.gload, self.heat_limit, self.cfg.gload_limit_obs_g, self.phi_t,
                     self.lam_t, self.atm)
        return torch.cat([b, self.hist.flatten(1)], 1).float()

    def potential(self, x=None) -> torch.Tensor:
        x = self.x if x is None else x
        ang, _ = range_course(x[:, 1], x[:, 2], self.phi_t, self.lam_t)
        return -ang * R_EARTH / 1000.0

    def expert_rows(self) -> np.ndarray:
        """physical_state-format rows (roll rate 0: the fast model has no roll dynamics)."""
        aD, aL, _ = aero_accels(self.x, self.sig, self.S_ref, self.mass, self.ls, self.atm)
        x = self.x
        return torch.stack([x[:, 0] - R_EARTH, x[:, 1], x[:, 2], x[:, 3], x[:, 4], x[:, 5], self.sig,
                            torch.zeros_like(self.sig), aD, aL], -1).cpu().numpy()

    # ---- step ----
    @torch.no_grad()
    def step(self, action: torch.Tensor):
        """One 5 s decision for every environment. Returns (obs, reward, done, info) as float32 tensors."""
        c = self.cfg
        cmd = action.to(self.dev, torch.float64).clamp(-1.0, 1.0) * SIGMA_MAX_RAD
        phi0 = self.potential()
        x, sig = self.x, self.sig
        gmax = torch.zeros_like(sig)
        Re = R_EARTH
        crossed = torch.zeros_like(sig, dtype=torch.bool)
        vert = torch.zeros_like(crossed)
        cross_phi, cross_lam = x[:, 1].clone(), x[:, 2].clone()
        frozen = torch.zeros_like(crossed)
        for ds in PH.SUBSTEPS_S:
            xn, gl = PH.eom_step(x, sig, self.S_ref, self.mass, self.ls, float(ds), self.atm)
            live = ~frozen
            hit = live & (xn[:, 0] - Re <= c.h_drogue_m)
            f = ((x[:, 0] - Re - c.h_drogue_m) / (x[:, 0] - xn[:, 0]).clamp_min(1e-9)).clamp(0, 1)
            cross_phi = torch.where(hit, x[:, 1] + f * (xn[:, 1] - x[:, 1]), cross_phi)
            cross_lam = torch.where(hit, x[:, 2] + f * (xn[:, 2] - x[:, 2]), cross_lam)
            x = torch.where(live[:, None], xn, x)
            gmax = torch.where(live, torch.maximum(gmax, gl), gmax)
            sig = sig + (cmd - sig).clamp(-PH.RATE_MAX_RPS * ds, PH.RATE_MAX_RPS * ds)
            crossed |= hit
            v = live & ~hit & (torch.cos(x[:, 4]).abs() < PH.COS_GAMMA_STOP)
            vert |= v
            frozen = frozen | hit | v
        self.x, self.sig = x, sig
        self.t = self.t + c.dt_s
        aD, aL, g_now = aero_accels(x, sig, self.S_ref, self.mass, self.ls, self.atm)
        self.gload = g_now
        skip = (x[:, 0] - Re > SKIP_ALT_M) & (x[:, 4] > 0)
        timeout = self.t >= c.t_max_s
        done = crossed | vert | skip | timeout | ~torch.isfinite(x).all(1)
        # MPPI's cost, realized
        r = -c.w_smooth_km_per_rad * (cmd - self.u_prev).abs()
        r = r - c.w_g_km * torch.relu(gmax - c.g_limit) ** 2
        steep = torch.relu(torch.rad2deg(x[:, 4].abs()) - c.gamma_limit_deg)
        r = r - c.w_gamma_km_per_deg2 * torch.where(crossed, torch.zeros_like(steep), steep ** 2)
        end_phi = torch.where(crossed, cross_phi, x[:, 1])
        end_lam = torch.where(crossed, cross_lam, x[:, 2])
        ang, _ = range_course(end_phi, end_lam, self.phi_t, self.lam_t)
        miss = ang * R_EARTH / 1000.0
        term = -miss - torch.where(crossed, torch.zeros_like(miss), torch.full_like(miss, c.nocross_km))
        r = r + torch.where(done, term, torch.zeros_like(term))
        if c.shaping:
            phi1 = torch.where(done, torch.zeros_like(phi0), self.potential())
            r = r + c.gamma_rl * phi1 - phi0
        r = r * c.reward_scale
        self.u_prev = cmd
        h = x[:, 0] - Re
        self.hist = torch.cat([self.hist[:, 1:], hist_row(h, x[:, 3], x[:, 4], aD, aL, sig)[:, None]], 1)
        self.ep_return += r
        self.ep_len += 1
        info = {"miss_km": miss, "success_5km": done & crossed & (miss < 5.0), "crossed": crossed,
                "outcome_drogue": crossed, "ep_return": self.ep_return.clone(), "ep_len": self.ep_len.clone(),
                "done": done.clone()}
        idx = torch.nonzero(done).flatten()
        if len(idx):
            self.reset_rows(idx)
        return self.obs(), r.float(), done, info


# ============================================================================
# The same observation in the full simulator
# ============================================================================
class FullObs:
    """rl_env._build_obs() + the 24-decision history, for a policy flying the full simulator at 5 s."""

    def __init__(self):
        self.hist: Optional[torch.Tensor] = None

    def __call__(self, env, y: np.ndarray) -> torch.Tensor:
        yt = torch.as_tensor(np.asarray(y, dtype=np.float64))
        row = hist_row(yt[0:1], yt[3:4], yt[4:5], yt[8:9], yt[9:10], yt[6:7])[0]
        if self.hist is None:
            self.hist = row[None].repeat(HIST_LEN, 1)
        else:
            self.hist = torch.cat([self.hist[1:], row[None]])
        base = torch.as_tensor(env._build_obs()[:N_BASE], dtype=torch.float64)
        return torch.cat([base, self.hist.flatten()]).float()


def check_matches_rollout_physics(n: int = 16, K: int = 40, seed: int = 0) -> float:
    """Smoke check 6: the env's dynamics equal physics_v4.rollout_physics for the same commands (max |diff|)."""
    env = FastEntryEnv(FastEnvCfg(n_envs=n, seed=seed, shaping=False))
    g = torch.Generator().manual_seed(seed)
    u = (torch.rand(n, K, generator=g, dtype=torch.float64) * 2 - 1) * SIGMA_MAX_RAD
    y0 = env.x.clone()
    y0[:, 0] -= R_EARTH
    ref, _ = PH.rollout_physics(y0, env.sig.clone(), u, env.S_ref, env.mass, env.ls, -1e9, env.atm)
    err = 0.0
    for k in range(K):
        x_before = env.x.clone()
        env.cfg.h_drogue_m = -1e9                                # never terminate on altitude here
        env.step(u[:, k] / SIGMA_MAX_RAD)
        live = torch.isfinite(ref[:, k + 1]).all(1) & (torch.cos(x_before[:, 4]).abs() >= PH.COS_GAMMA_STOP)
        got = env.x.clone()
        got[:, 0] -= R_EARTH
        d = (got - ref[:, k + 1]).abs()[live & (env.t > 0)]
        if len(d):
            err = max(err, float(d.max()))
        if (env.ep_len == 0).any():
            break
    return err
