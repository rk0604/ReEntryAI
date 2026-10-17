"""
jepa/data.py

Dataset synthesis for the JEPA world-model experiment (see jepa/PLAN.md).

SkyJEPA's data recipe (Sec. IV-E), mapped onto entry guidance:
  Diverse excitation    random smooth bank profiles, Gaussian-process and
                        reversal-structured, in place of GP position references.
  Dynamically feasible  every profile is flown through the real closed loop
                        (bank actuator -> roll PD -> 12 jet RCS -> 3DOF EOM) by
                        OrionEntryEnv in policy mode, so actuator lag and
                        rate limits are in the data.
  Robust                per episode domain randomization of density (q*S
                        scale), L/D, and mass, on top of the config's entry
                        interface dispersion.

Everything is logged at the 1 Hz guidance rate and subsampled to the world
model macro step afterwards, so the macro step can change without recollecting.

Physical state layout y (Y_NAMES), one row per logged step:
    h, phi, lam, V, gamma, chi    the translational state (altitude, not radius)
    sigma_act, roll_rate          realized bank and roll rate from the attitude
    a_drag, a_lift                aero specific force magnitudes (what an IMU
                                  senses), in m/s^2

Domain randomization uses knobs the simulator already reads from params, so no
core physics module changes:
    ref_area_m2 * k_rho     q*S scales exactly like a density scale for forces
    cpas_lift_scale = k_LD  scales CL (inactive above drogue, so free to reuse)
    mass_kg * k_m
The classical guidance holds its own params copy made at env construction, so
mutating env.params changes the truth but not the predictor-corrector's model.
"""

from __future__ import annotations

import math
import multiprocessing as mp
import pickle
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import torch

import constants

# World model / planner step. One bank command is held for this long.
MACRO_DT_S = 5.0
SIGMA_MAX_RAD = math.pi / 2.0     # same handle as the RL env: action 1 -> 90 deg
SKIP_OUT_ALT_M = 160_000.0       # top of the atmosphere table; a flight above it has skipped out
G0 = float(constants.STANDARD_GRAVITY_MPS2)
R_EARTH = float(constants.RADIUS_EARTH)

Y_NAMES = ("h_m", "phi_rad", "lam_rad", "V_mps", "gamma_rad", "chi_rad",
           "sigma_act_rad", "roll_rate_rps", "a_drag_mps2", "a_lift_mps2")
# Encoder features derived from y. Longitude is dropped: on a non-rotating
# Earth with the altitude-only atmosphere table the dynamics do not depend on it.
FEAT_NAMES = ("h_km", "V_kms", "gamma", "sin_chi", "cos_chi", "phi",
              "sigma_act", "roll_rate", "a_drag_g", "a_lift_g")
N_Y = len(Y_NAMES)
N_FEAT = len(FEAT_NAMES)


# ============================================================================
# Domain randomization
# ============================================================================
@dataclass
class DomainSpec:
    """Uniform ranges for the per-episode truth perturbations."""
    density_scale: tuple = (0.85, 1.15)   # applied through ref_area_m2
    lift_scale: tuple = (0.85, 1.15)      # applied through cpas_lift_scale
    mass_scale: tuple = (0.95, 1.05)


NOMINAL_DOMAIN = {"density_scale": 1.0, "lift_scale": 1.0, "mass_scale": 1.0}

# Wider truth randomization for the wide-envelope dataset.
WIDE_DOMAIN = DomainSpec(density_scale=(0.70, 1.30), lift_scale=(0.70, 1.30), mass_scale=(0.90, 1.10))


@dataclass
class EnvelopeSpec:
    """
    Uniform ranges for the entry-interface state of wide-envelope data.

    The default config only disperses the entry by about 1 km, 100 m/s, 0.4 deg
    of flight-path angle and 3 deg of heading, so every flight starts in the same
    corridor. These ranges cover entry altitude, sub-circular to super-circular
    speed, shallow to steep flight-path angle, any heading, and a band of
    latitudes. Longitude barely matters (the dynamics do not depend on it).

    Fast, shallow entries skip back out of the atmosphere, so the shallowest
    allowed flight-path angle steepens with speed: from gamma_deg's upper bound
    at skip_guard_speed_mps[0] down to skip_guard_gamma_deg at [1] (checked with
    bank 0, thin air, high L/D, the worst case for skipping).
    """
    altitude_m: tuple = (100_000.0, 140_000.0)
    velocity_mps: tuple = (6_800.0, 8_600.0)
    gamma_deg: tuple = (-6.0, -1.0)
    chi_deg: tuple = (0.0, 360.0)
    phi_deg: tuple = (-30.0, 30.0)
    lambda_deg: tuple = (-20.0, 20.0)
    skip_guard_speed_mps: tuple = (7_600.0, 8_600.0)
    skip_guard_gamma_deg: float = -4.0


def sample_entry_ic(rng: np.random.Generator, spec: EnvelopeSpec) -> Dict[str, float]:
    """Draw one entry-interface condition (an ic_override dict) from the envelope."""
    V = float(rng.uniform(*spec.velocity_mps))
    v0, v1 = spec.skip_guard_speed_mps
    frac = min(max((V - v0) / (v1 - v0), 0.0), 1.0)
    shallowest = spec.gamma_deg[1] + frac * (spec.skip_guard_gamma_deg - spec.gamma_deg[1])
    return {"altitude_m": float(rng.uniform(*spec.altitude_m)), "velocity_mps": V,
            "gamma_deg": float(rng.uniform(spec.gamma_deg[0], shallowest)),
            "chi_deg": float(rng.uniform(*spec.chi_deg)), "phi_deg": float(rng.uniform(*spec.phi_deg)),
            "lambda_deg": float(rng.uniform(*spec.lambda_deg))}


def sample_domain(rng: np.random.Generator, spec: DomainSpec) -> Dict[str, float]:
    """Draw one set of truth perturbation factors from the spec ranges."""
    return {k: float(rng.uniform(*v)) for k, v in asdict(spec).items()}


def apply_domain(params: dict, domain: Dict[str, float]) -> dict:
    """
    Set the perturbed aero and mass values on a params dict, in place.

    Values are set from the nominal constants, not multiplied onto the current
    ones, so applying a domain is idempotent and safe on a reused env.
    """
    params["ref_area_m2"] = float(constants.CAPSULE_REFERENCE_AREA_M2) * float(domain["density_scale"])
    params["cpas_lift_scale"] = float(domain["lift_scale"])
    params["mass_kg"] = float(constants.CAPSULE_MASS_KG) * float(domain["mass_scale"])
    return params


# ============================================================================
# Bank profile generators (the excitation signal)
# ============================================================================
def _gp_sample(rng: np.random.Generator, n: int, length_scale: float) -> np.ndarray:
    """One zero-mean unit-variance draw from an RBF Gaussian process on a unit grid."""
    t = np.arange(n, dtype=np.float64)
    K = np.exp(-0.5 * ((t[:, None] - t[None, :]) / length_scale) ** 2)
    L = np.linalg.cholesky(K + 1e-6 * np.eye(n))
    return L @ rng.standard_normal(n)


def sample_bank_profile(rng: np.random.Generator, n_macro: int, kind: str) -> np.ndarray:
    """
    Draw one bank command profile, one value per macro step, in radians.

    Parameters:
        rng:     random generator.
        n_macro: number of macro steps to cover.
        kind:    "gp" for a squashed multi-length-scale GP (broad coverage of
                 magnitude and sign), or "reversal" for a slowly varying
                 magnitude with random sign flips (the structure real entry
                 guidance produces).

    Returns:
        Float array of shape [n_macro] inside [-90 deg, +90 deg].
    """
    if kind == "gp":
        # Slow drift + medium wander + fast wiggle, in macro steps (5 s each).
        u = sum(rng.uniform(0.3, 1.5) * _gp_sample(rng, n_macro, ls) for ls in (3.0, 10.0, 30.0))
        return SIGMA_MAX_RAD * np.tanh(u + rng.uniform(-1.0, 1.0))
    if kind == "reversal":
        mag0 = math.radians(rng.uniform(15.0, 85.0))
        mag = np.clip(mag0 + math.radians(15.0) * _gp_sample(rng, n_macro, 30.0),
                      math.radians(5.0), SIGMA_MAX_RAD)
        # Reversals arrive as a Poisson process, one every ~60 to ~200 s.
        flips = rng.random(n_macro) < rng.uniform(1.0 / 40.0, 1.0 / 12.0)
        sign = rng.choice([-1.0, 1.0]) * np.where(np.cumsum(flips) % 2 == 0, 1.0, -1.0)
        return sign * mag
    raise ValueError(f"unknown profile kind {kind!r}")


# ============================================================================
# Collection (one worker-safe episode, and a parallel driver)
# ============================================================================
def physical_state(env) -> np.ndarray:
    """
    Read the physical state y (Y_NAMES order) from an OrionEntryEnv.

    The aero specific forces are recomputed from the env's own (possibly
    perturbed) params, at the realized bank, so they are what an onboard
    accelerometer would sense.
    """
    from math_3d import nominal_aero_forces_from_state
    r, phi, lam, V, gamma, chi = (float(v) for v in env.x)
    sigma = float(env.att.sigma_rel_rad)
    aero = nominal_aero_forces_from_state(x=list(env.x), sigma_rad=sigma, params=env.params)
    m = float(env.params["mass_kg"])
    return np.array([r - R_EARTH, phi, lam, V, gamma, chi, sigma,
                     float(env.att.omega_b_rad_s[2]),
                     float(aero["drag_mag_N"]) / m, float(aero["lift_mag_N"]) / m],
                    dtype=np.float64)


_WORKER_ENVS: Dict[tuple, object] = {}


def _policy_env(config_path: str, max_steps: int):
    """Build (once per process) a policy mode env with dispersion on."""
    key = (config_path, int(max_steps))
    if key not in _WORKER_ENVS:
        import rl_env
        _WORKER_ENVS[key] = rl_env.OrionEntryEnv(
            config_path=config_path, controller_mode="policy", dispersion=True,
            max_episode_steps=int(max_steps), fast_atmosphere=True)
    return _WORKER_ENVS[key]


def collect_episode(job: dict) -> dict:
    """
    Fly one random bank profile through the real closed loop and log it.

    Parameters:
        job: dict with config_path, seed, kind ("gp" or "reversal"),
             randomize (bool), spec (DomainSpec fields as a dict), max_steps,
             and optionally envelope (EnvelopeSpec fields as a dict) to start
             from a wide-envelope entry condition instead of the config's.

    Returns:
        Episode dict: y [N+1, N_Y] float32 at 1 Hz, cmd [N] float32 (rad, the
        command flown during each 1 s step), domain, kind, seed, outcome,
        miss_km, t_final_s, ic.
    """
    rng = np.random.default_rng(int(job["seed"]))
    env = _policy_env(job["config_path"], job["max_steps"])
    domain = (sample_domain(rng, DomainSpec(**job["spec"])) if job["randomize"]
              else dict(NOMINAL_DOMAIN))
    apply_domain(env.params, domain)
    # Drawn after the domain so jobs without an envelope keep their old seeds' results.
    options = ({"ic_override": sample_entry_ic(rng, EnvelopeSpec(**job["envelope"]))}
               if job.get("envelope") else None)
    _, info = env.reset(seed=int(job["seed"]), options=options)

    n_sub = int(round(MACRO_DT_S / env.guidance_period_s))
    profile = sample_bank_profile(rng, job["max_steps"] // n_sub + 2, job["kind"])
    ys = [physical_state(env)]
    cmds = []
    done, k, outcome = False, 0, None
    while not done:
        sigma = float(profile[k // n_sub])
        _, _, term, trunc, info = env.step(np.array([sigma / SIGMA_MAX_RAD], dtype=np.float32))
        cmds.append(sigma)
        ys.append(physical_state(env))
        k += 1
        done = bool(term or trunc)
        # Above the atmosphere table there is nothing left to learn, only coasting.
        if ys[-1][0] > SKIP_OUT_ALT_M:
            done, outcome = True, "skip_out"
    return {
        "y": np.asarray(ys, dtype=np.float32), "cmd": np.asarray(cmds, dtype=np.float32),
        "domain": domain, "kind": job["kind"], "seed": int(job["seed"]),
        "outcome": outcome or info["outcome"], "miss_km": float(info["miss_km"]),
        "t_final_s": float(info["t_s"]), "ic": dict(info["ic"]),
    }


def run_parallel(fn, jobs: Sequence[dict], workers: int, label: str = "jobs") -> List[dict]:
    """
    Map fn over jobs in a spawn process pool (Windows safe), keeping order.

    fn must be importable from a module (not defined in a notebook), which is
    why the collection and episode workers live in this package.
    """
    t0 = time.time()
    out: List[Optional[dict]] = [None] * len(jobs)
    if workers <= 1:
        for i, j in enumerate(jobs):
            out[i] = fn(j)
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(int(workers)) as pool:
            for n_done, (i, r) in enumerate(pool.imap_unordered(_indexed_call, [(fn, i, j) for i, j in enumerate(jobs)]), 1):
                out[i] = r
                if n_done % max(1, len(jobs) // 10) == 0 or n_done == len(jobs):
                    print(f"  [{label}] {n_done}/{len(jobs)}  ({time.time() - t0:.0f}s)")
    return out  # type: ignore[return-value]


def _indexed_call(arg):
    """Pool helper: call fn(job) and return it tagged with its index."""
    fn, i, job = arg
    return i, fn(job)


def collect_dataset(n_episodes: int, config_path: str = "configs/default.json",
                    seed0: int = 1_000, workers: int = 1, max_steps: int = 900,
                    kinds: Sequence[str] = ("gp", "reversal"), randomize: bool = True,
                    spec: Optional[DomainSpec] = None,
                    cache_path: Optional[str] = None) -> List[dict]:
    """
    Collect (or load from cache) a set of randomized closed-loop episodes.

    Seeds start at seed0 so they never coincide with the frozen test set
    (whose ICs are replayed by value, not by seed, anyway). A smaller cache is
    extended rather than recollected.
    """
    eps: List[dict] = []
    if cache_path and Path(cache_path).exists():
        with open(cache_path, "rb") as f:
            eps = pickle.load(f)
        print(f"[data] loaded {len(eps)} cached episodes <- {cache_path}")
        if len(eps) >= n_episodes:
            return eps[:n_episodes]
    # Only collect what the cache is missing; episode i always uses seed0 + i.
    spec = spec or DomainSpec()
    jobs = [{"config_path": config_path, "seed": seed0 + i, "kind": kinds[i % len(kinds)],
             "randomize": randomize, "spec": asdict(spec), "max_steps": max_steps}
            for i in range(len(eps), n_episodes)]
    eps = eps + run_parallel(collect_episode, jobs, workers, label="collect")
    if cache_path:
        Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "wb") as f:
            pickle.dump(eps, f)
        print(f"[data] saved {len(eps)} episodes -> {cache_path}")
    return eps


def collect_shards(n_total: int, out_dir: str, shard_size: int = 250, seed0: int = 100_000,
                   workers: int = 1, max_steps: int = 1200, kinds: Sequence[str] = ("gp", "reversal"),
                   spec: Optional[DomainSpec] = None, envelope: Optional[EnvelopeSpec] = None,
                   config_path: str = "configs/default.json") -> None:
    """
    Collect a large dataset as numbered shard files, resumable across sessions.

    Shard k holds episodes k*shard_size .. (k+1)*shard_size-1, and episode i
    always uses seed seed0 + i, so a rerun skips finished shards and produces
    the same data. Each shard is written to a temp file and renamed, so a
    disconnect mid-write never leaves a corrupt shard behind.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    spec = spec or DomainSpec()
    n_shards = (n_total + shard_size - 1) // shard_size
    t0, done_now = time.time(), 0
    for k in range(n_shards):
        path = out / f"shard_{k:03d}.pkl"
        if path.exists():
            continue
        ids = range(k * shard_size, min(n_total, (k + 1) * shard_size))
        jobs = [{"config_path": config_path, "seed": seed0 + i, "kind": kinds[i % len(kinds)],
                 "randomize": True, "spec": asdict(spec), "max_steps": max_steps,
                 "envelope": asdict(envelope) if envelope else None} for i in ids]
        eps = run_parallel(collect_episode, jobs, workers, label=f"shard {k + 1}/{n_shards}")
        tmp = path.with_suffix(".tmp")
        with open(tmp, "wb") as f:
            pickle.dump(eps, f)
        tmp.replace(path)
        done_now += 1
        left = sum(not (out / f"shard_{j:03d}.pkl").exists() for j in range(n_shards))
        per = (time.time() - t0) / done_now
        print(f"[data] shard {k + 1}/{n_shards} saved ({len(eps)} episodes) | "
              f"{left} left, ~{per * left / 60:.0f} min remaining")
    print(f"[data] all {n_shards} shards present in {out}")


def load_shards(out_dir: str) -> List[dict]:
    """Concatenate every shard in a directory, in shard order."""
    eps: List[dict] = []
    for path in sorted(Path(out_dir).glob("shard_*.pkl")):
        with open(path, "rb") as f:
            eps += pickle.load(f)
    return eps


# ============================================================================
# Macro-step view of an episode
# ============================================================================
def to_macro(ep: dict, macro_dt_s: float = MACRO_DT_S, guidance_dt_s: float = 1.0) -> dict:
    """
    Subsample a 1 Hz episode onto the macro grid.

    Returns:
        dict with y [L+1, N_Y] (state at each macro boundary), u [L] (command
        held over each macro step), acc [L, 3] (mean drag, vertical lift, and
        lateral lift specific force over each macro step, the prober's
        auxiliary target), y_final [N_Y] (the drogue-crossing state), plus the
        episode's domain and outcome. The trailing partial macro step is dropped.
    """
    n = int(round(macro_dt_s / guidance_dt_s))
    y, cmd = ep["y"].astype(np.float64), ep["cmd"]
    L = len(cmd) // n
    acc = np.zeros((L, 3))
    for k in range(L):
        seg = y[k * n:(k + 1) * n + 1]
        aL, sig = seg[:, 9], seg[:, 6]
        acc[k] = (seg[:, 8].mean(), (aL * np.cos(sig)).mean(), (aL * np.sin(sig)).mean())
    return {"y": y[0:L * n + 1:n].astype(np.float32), "u": cmd[0:L * n:n].astype(np.float32),
            "acc": acc.astype(np.float32), "y_final": y[-1].astype(np.float32),
            "domain": ep["domain"], "outcome": ep["outcome"], "kind": ep.get("kind", "")}


def raw_features(y: torch.Tensor) -> torch.Tensor:
    """Map physical states [..., N_Y] to encoder features [..., N_FEAT] (unnormalized)."""
    chi = y[..., 5]
    return torch.stack([y[..., 0] / 1000.0, y[..., 3] / 1000.0, y[..., 4],
                        torch.sin(chi), torch.cos(chi), y[..., 1],
                        y[..., 6], y[..., 7], y[..., 8] / G0, y[..., 9] / G0], dim=-1)


@dataclass
class Stats:
    """Normalization statistics, computed on the training split only."""
    feat_mean: torch.Tensor
    feat_std: torch.Tensor
    dy_std: torch.Tensor      # [N_Y] std of one macro-step state change

    def norm_feats(self, y: torch.Tensor) -> torch.Tensor:
        """Physical states -> normalized encoder features."""
        return (raw_features(y) - self.feat_mean) / self.feat_std

    def to_dict(self) -> dict:
        """Plain-tensor dict for checkpointing."""
        return {"feat_mean": self.feat_mean, "feat_std": self.feat_std, "dy_std": self.dy_std}


def compute_stats(macro_eps: Sequence[dict]) -> Stats:
    """Feature mean/std and per-step state-change std over a list of macro episodes."""
    Y = torch.from_numpy(np.concatenate([e["y"] for e in macro_eps])).double()
    dY = torch.from_numpy(np.concatenate([np.diff(e["y"], axis=0) for e in macro_eps])).double()
    f = raw_features(Y)
    return Stats(feat_mean=f.mean(0).float(), feat_std=f.std(0).clamp_min(1e-6).float(),
                 dy_std=dY.std(0).clamp_min(1e-9).float())


class MacroDataset:
    """
    Fixed-length training windows over macro episodes.

    Each episode is front-padded with P = rf - 1 copies of its first state and
    P zero commands (the capsule starts wings level), so windows can start at
    entry interface just like the planner does. A window starting at real
    macro index t yields:
        feats [P+T+1, N_FEAT]  normalized features of y[t-P .. t+T]
        u     [P+T]            commands u[t-P .. t+T-1]
        y     [T+1, N_Y]       physical states y[t .. t+T]
        acc   [T, 3]           macro-mean specific forces over steps t .. t+T-1
    """

    def __init__(self, macro_eps: Sequence[dict], rf: int, horizon: int, stats: Stats):
        self.P, self.T, self.stats = int(rf) - 1, int(horizon), stats
        P, T = self.P, self.T
        Ys, Us, As, y_starts, u_starts = [], [], [], [], []
        y_off = u_off = 0
        for ep in macro_eps:
            y, u, acc = ep["y"], ep["u"], ep["acc"]
            L = len(u)
            if L < T:
                continue
            Ys.append(np.concatenate([np.repeat(y[:1], P, 0), y]))
            Us.append(np.concatenate([np.zeros(P, np.float32), u]))
            As.append(np.concatenate([np.zeros((P, 3), np.float32), acc]))
            t = np.arange(L - T + 1)
            y_starts.append(y_off + t)
            u_starts.append(u_off + t)
            y_off += L + 1 + P
            u_off += L + P
        self.Y = torch.from_numpy(np.concatenate(Ys)).float()
        self.F = stats.norm_feats(self.Y)
        self.U = torch.from_numpy(np.concatenate(Us)).float()
        self.A = torch.from_numpy(np.concatenate(As)).float()
        self.y_start = torch.from_numpy(np.concatenate(y_starts)).long()
        self.u_start = torch.from_numpy(np.concatenate(u_starts)).long()

    def __len__(self) -> int:
        return len(self.y_start)

    def batch(self, idx: torch.Tensor) -> Dict[str, torch.Tensor]:
        """Gather the windows with the given indices."""
        P, T = self.P, self.T
        iy = self.y_start[idx, None] + torch.arange(P + T + 1)
        iu = self.u_start[idx, None] + torch.arange(P + T)
        return {"feats": self.F[iy], "u": self.U[iu],
                "y": self.Y[iy[:, P:]], "acc": self.A[iu[:, P:]]}

    def sample(self, batch_size: int, gen: Optional[torch.Generator] = None) -> Dict[str, torch.Tensor]:
        """A random batch of windows."""
        return self.batch(torch.randint(len(self), (batch_size,), generator=gen))


def split_episodes(macro_eps: Sequence[dict], fracs=(0.8, 0.1, 0.1), seed: int = 0):
    """Shuffle and split episodes (not windows) into train / val / test."""
    idx = np.random.default_rng(seed).permutation(len(macro_eps))
    n_tr = int(round(fracs[0] * len(idx)))
    n_va = int(round(fracs[1] * len(idx)))
    pick = lambda ii: [macro_eps[i] for i in ii]
    return pick(idx[:n_tr]), pick(idx[n_tr:n_tr + n_va]), pick(idx[n_tr + n_va:])


def load_test_ics(path: str = "configs/heldout_testset.csv", n: Optional[int] = None) -> List[dict]:
    """The frozen held-out entry-interface conditions, as ic_override dicts."""
    import pandas as pd
    df = pd.read_csv(path, comment="#")
    cols = ["altitude_m", "velocity_mps", "gamma_deg", "chi_deg", "phi_deg", "lambda_deg"]
    rows = df[cols].to_dict("records")
    return [{k: float(v) for k, v in r.items()} for r in (rows[:n] if n else rows)]
