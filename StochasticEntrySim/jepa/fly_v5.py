"""
jepa/fly_v5.py

Closed-loop flights for the v5 study (jepa/v5_plan.md sections 4.1-4.4).

    flight_id_v5            <suite>/<arm>/<cond> for deterministic arms (guidance laws, classical, RL
                            policies), <suite>/<arm>/<load>/seed<k>/mppi<m>/<cond> for MPPI arms;
                            cond = <atm>_case<NN> (A256, DEV48) or case<NNNN> (ORB-1)
    run_guidance_lockstep   the 1 Hz lockstep runner for guidance laws and policies: every flight of a
                            batch steps 1 s at a time, and the law decides for all of them at once (FNPEG's
                            predictions run as one batch). Bank limit 90 deg through the simulator's own
                            action clip; 180 deg arms replace the clipped command on the env's own
                            PolicyGuidance instance (the actuator wraps to +-pi), no simulator change
    PolicyLaw               RL policies in the same runner: v5 policies every 5 s on rlenv_fast_v5.FullObs,
                            legacy SB3 policies at 1 Hz on rl_env's observation (job["env_kwargs"] for their env)
    run_mppi_v5             fly_v4.run_lockstep_v4 unchanged (so E0 holds by construction), each record
                            relabelled with the v5 flight id and metrics
    FilteredPhysicsAdapter  mppi_filtered: v4's PhysicsAdapter("nominal") whose density and lift scales
                            come from FNPEG's fading-memory filters every decision (density_scale = K_D,
                            lift_scale = K_L / K_D, mass 1); frozen=True keeps them at 1 (= mppi_nominal)
    load_adapter_v5         v4's adapters, physics_filtered, and the E3 models (models_v5)
    v5_metrics              success at 2/10 km, NASA-STD-3000 g-duration compliance, authority
                            saturation, planning time median/p99
    fnpeg_peak_g            per-condition matched g limits from fnpeg's flights
"""

from __future__ import annotations

import math
import time
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
import torch

from jepa import fly_v3 as F
from jepa import fly_v4 as F4
from jepa.data import G0, NOMINAL_DOMAIN, SIGMA_MAX_RAD, apply_domain, physical_state

SAT_TOL_DEG = 1.0          # "at the limit": within 1 deg of 0 or of the bank limit


# ============================================================================
# Keys
# ============================================================================
def cond_key(job: dict) -> str:
    if job.get("suite") == "ORB1":
        return f"case{int(job['case']):04d}"
    return f"{job['atmosphere']}_case{int(job['case']):02d}"


def deterministic(job: dict) -> bool:
    return (job.get("model") or {}).get("kind") in ("classical", "fnpeg", "predguid", "policy", "legacy_rl")


def flight_id_v5(job: dict) -> str:
    """Relative path (no extension) of one v5 flight's result and telemetry."""
    if deterministic(job):
        ts = job.get("train_seed") if (job.get("model") or {}).get("kind") == "policy" else None
        if ts is not None:                       # an RL policy: one deterministic flight per training seed
            return f"{job['suite']}/{job['controller']}/seed{int(ts)}/{cond_key(job)}"
        return f"{job['suite']}/{job['controller']}/{cond_key(job)}"
    return (f"{job['suite']}/{job['controller']}/{job['load']}/seed{int(job['seed'])}/"
            f"mppi{int(job.get('mppi_rep', 0))}/{cond_key(job)}")


def result_extras_v5(job: dict) -> dict:
    return {"flight_id": flight_id_v5(job), "suite": job.get("suite"), "arm": job.get("controller"),
            "stage": job.get("stage"), "block": job.get("block"), "target_index": job.get("target_index"),
            "train_seed": job.get("train_seed"), "feasible": job.get("feasible"),
            "in_training_domain": job.get("in_training_domain"), "model_kind": (job.get("model") or {}).get("kind"),
            "sigma_max_deg": float(job.get("sigma_max_deg", 90.0)), "source": "v5"}


# ============================================================================
# Per-flight metrics (section 5)
# ============================================================================
def v5_metrics(result: dict, tel: Dict[str, np.ndarray], sigma_max_deg: float = 90.0,
               plan_ms: Sequence[float] = ()) -> dict:
    from jepa.benchmark_v5 import g_duration
    t = tel["t_s"].astype(np.float64)
    g = tel["g"].astype(np.float64)
    dt = float(np.median(np.diff(t))) if len(t) > 1 else 1.0
    gd = g_duration(g[1:], dt)
    cmd = np.abs(tel["sigma_cmd_deg"][1:].astype(np.float64))
    guided = g[1:] >= 1.52 / G0                                  # "guided time": aero acceleration >= 5 ft/s^2
    sat = ((cmd <= SAT_TOL_DEG) | (cmd >= sigma_max_deg - SAT_TOL_DEG)) & guided
    landed = result.get("outcome") == "drogue"
    miss = float(result.get("miss_km", np.nan))
    pm = np.asarray(plan_ms, dtype=np.float64)
    return {"success_2km": bool(landed and miss < 2.0), "success_10km": bool(landed and miss <= 10.0),
            "g_duration_compliant": gd["compliant"], "g_duration_margin_g": gd["margin_g"],
            "sat_frac": float(sat.sum() / max(int(guided.sum()), 1)),
            "plan_ms_median": float(np.median(pm)) if len(pm) else float("nan"),
            "plan_ms_p99": float(np.percentile(pm, 99)) if len(pm) else float("nan")}


# ============================================================================
# The 1 Hz guidance runner
# ============================================================================
class GuidanceFlight:
    """One closed-loop flight flown by a guidance law (or policy) through the full simulator."""

    def __init__(self, job: dict):
        import rl_env
        self.job = job
        kw = {"max_episode_steps": 1400, "fast_atmosphere": True, **(job.get("env_kwargs") or {})}
        self.env = rl_env.OrionEntryEnv(config_path=job.get("config_path", "configs/default.json"),
                                        controller_mode="policy", dispersion=False, **kw)
        self.domain = dict(job.get("domain") or NOMINAL_DOMAIN)
        apply_domain(self.env.params, self.domain)
        self.env.reset(seed=int(job.get("env_seed", 0)), options={"ic_override": job["ic"]} if job.get("ic") else None)
        self.sigma_max = math.radians(float(job.get("sigma_max_deg", 90.0)))
        self._cmd = 0.0
        if self.sigma_max > SIGMA_MAX_RAD + 1e-9:
            orig = self.env.guidance.set_action
            # rl_env clips the action to +-90 deg and then calls set_action(clipped); this flight's own
            # PolicyGuidance instance receives the unclipped command instead.
            self.env.guidance.set_action = lambda _clipped, _o=orig, _f=self: _o(_f._cmd)
        self.rec = F.Recorder(self.env)
        self.cmds: List[float] = []
        self.plan_ms: List[float] = []
        self.done, self.info, self.k = False, {}, 0
        self.t0 = time.time()
        self.state = None                       # the law's per-flight state

    def step(self, cmd: float) -> None:
        cmd = float(max(-self.sigma_max, min(self.sigma_max, cmd)))
        self._cmd = cmd
        _, _, term, trunc, self.info = self.env.step(np.array([cmd / SIGMA_MAX_RAD], dtype=np.float32))
        self.cmds.append(cmd)
        self.rec.add(self.env, cmd)
        self.k += 1
        self.done = bool(term or trunc)

    def record(self, law=None) -> dict:
        tel = self.rec.arrays()
        # fly_v3.finish builds a v3 id from load and seed, which guidance arms do not have; the v5 id replaces it
        shim = {**self.job, "load": self.job.get("load") or "-", "seed": self.job.get("seed") or 0}
        res = F.finish(shim, self.env, self.info, self.cmds, tel, time.time() - self.t0, self.plan_ms)
        res.update({"load": self.job.get("load"), "seed": self.job.get("seed")})
        res.update(result_extras_v5(self.job))
        res.update(v5_metrics(res, tel, math.degrees(self.sigma_max), self.plan_ms))
        if law is not None and hasattr(law, "flight_summary"):
            res.update(law.flight_summary(self.state))
            tel.update(law.decisions(self.state))
        return {"result": res, "telemetry": tel}


@torch.no_grad()
def run_guidance_lockstep(jobs: Sequence[dict], law, on_done: Optional[Callable[[dict], None]] = None,
                          progress: Optional[Callable[[dict], None]] = None, verbose: bool = True) -> List[dict]:
    """
    Fly every job with `law` (FNPEGLaw, PredGuidLaw, PolicyLaw), all flights in lockstep at 1 Hz.
    law.new_flight(job, env) -> state; law.decide(states, ys, t_s, envs, ks) -> bank commands (rad).
    plan_ms per flight = the batched decision's wall time / number of deciding flights.
    """
    flights = [GuidanceFlight(j) for j in jobs]
    for f in flights:
        f.state = law.new_flight(f.job, f.env)
    out: List[Optional[dict]] = [None] * len(flights)
    t0, n_round = time.time(), 0
    while True:
        active = [i for i, f in enumerate(flights) if not f.done]
        if not active:
            break
        fs = [flights[i] for i in active]
        tr = time.perf_counter()
        ys = [physical_state(f.env) for f in fs]
        cmds = law.decide([f.state for f in fs], ys, [float(f.env.t_s) for f in fs], [f.env for f in fs],
                          [f.k for f in fs])
        ms = 1000.0 * (time.perf_counter() - tr) / len(fs)
        for f, c in zip(fs, cmds):
            if getattr(f.state, "decided", True):
                f.plan_ms.append(ms)
            f.step(c)
        n_round += 1
        if progress is not None and n_round % 10 == 0:
            progress({"round": n_round, "flying": len(active), "total": len(flights),
                      "sec": round(time.time() - t0, 1), "ms_per_decision": round(ms, 2)})
        if verbose and n_round % 100 == 0:
            print(f"  [guidance] t={n_round}s: {len(active)}/{len(flights)} flying, {ms:.1f} ms/flight-decision, "
                  f"{time.time() - t0:.0f}s", flush=True)
        for i in active:
            f = flights[i]
            if f.done:
                out[i] = f.record(law)
                if on_done is not None:
                    on_done(out[i])
    return out  # type: ignore[return-value]


class FNPEGLaw:
    """guidance_fnpeg.FNPEG behind the runner's interface."""

    def __init__(self, cfg=None, device="cpu"):
        from jepa.guidance_fnpeg import FNPEG
        self.g = FNPEG(cfg, device)

    def new_flight(self, job, env):
        return self.g.new_flight(env.target_phi_rad, env.target_lam_rad, job.get("domain"))

    def decide(self, states, ys, t_s, envs, ks):
        return self.g.decide(states, ys, t_s)

    def flight_summary(self, st) -> dict:
        it = np.asarray(st.iters, dtype=float)
        return {"fnpeg_cycles": int(len(it)), "fnpeg_iter_mean": float(it.mean()) if len(it) else float("nan"),
                "fnpeg_iter_max": int(it.max()) if len(it) else -1,
                "fnpeg_iter_le6_frac": float((it <= 6).mean()) if len(it) else float("nan"),
                "fnpeg_K_D_final": float(st.K_D), "fnpeg_K_L_final": float(st.K_L)}

    def decisions(self, st) -> Dict[str, np.ndarray]:
        from jepa.guidance_fnpeg import decisions_array
        return decisions_array(st)


def fly_guidance(jobs: Sequence[dict], law_name: str, law_cfg: dict, device="cpu", **kw) -> List[dict]:
    """Build the law from its name and settings, then fly the jobs."""
    return run_guidance_lockstep(jobs, make_law(law_name, law_cfg, device), **kw)


class PolicyLaw:
    """
    An RL policy behind the runner's interface. New v5 policies (rl_v5.Actor, observation normalizer inside)
    decide every 5 s on rlenv_fast_v5.FullObs (rl_env's 13 features + the last 24 decisions), deterministic;
    the command is held between decisions. Legacy SB3 policies decide at 1 Hz on rl_env's own observation.
    """

    def __init__(self, path: str, legacy: bool = False, device="cpu"):
        from jepa import rl_v5 as RL
        self.legacy = bool(legacy)
        if self.legacy:
            self.pol = RL.LegacyPolicy(path)
            self.period = 1
        else:
            self.actor = RL.load_policy(path, device)
            self.period = 5
        self.dev = torch.device(device)

    def new_flight(self, job, env):
        from types import SimpleNamespace
        from jepa.rlenv_fast_v5 import FullObs
        return SimpleNamespace(obs=FullObs(), cmd=0.0, decided=False, n_dec=0)

    def decide(self, states, ys, t_s, envs, ks):
        idx = [i for i, k in enumerate(ks) if int(k) % self.period == 0]
        for st in states:
            st.decided = False
        if idx:
            if self.legacy:
                a = np.array([self.pol.act(envs[i]._build_obs()) for i in idx])
            else:
                O = torch.stack([states[i].obs(envs[i], ys[i]) for i in idx]).to(self.dev)
                a = self.actor.deterministic(O).double().cpu().numpy()
            for j, i in enumerate(idx):
                states[i].cmd = float(np.clip(a[j], -1.0, 1.0)) * SIGMA_MAX_RAD
                states[i].decided = True
                states[i].n_dec += 1
        return [st.cmd for st in states]

    def flight_summary(self, st) -> dict:
        return {"policy_decisions": int(st.n_dec), "policy_period_s": self.period}

    def decisions(self, st) -> Dict[str, np.ndarray]:
        return {}


def make_law(name: str, cfg: dict, device="cpu"):
    if name == "policy":
        return PolicyLaw(cfg["path"], legacy=bool(cfg.get("legacy")), device=device)
    if name == "fnpeg":
        from jepa.guidance_fnpeg import FNPEGConfig
        return FNPEGLaw(FNPEGConfig(**cfg), device)
    if name == "predguid":
        from jepa.guidance_predguid import PredGuidConfig, PredGuidLaw
        return PredGuidLaw(PredGuidConfig(**cfg), device)
    raise ValueError(f"unknown guidance law {name!r}")


# ============================================================================
# MPPI arms (v4's runner unchanged)
# ============================================================================
def run_mppi_v5(jobs: Sequence[dict], adapter, mppi: Optional[dict] = None, n_sub: int = 2, max_rows: int = 2048,
                on_done: Optional[Callable[[dict], None]] = None, verbose: bool = True,
                progress: Optional[Callable[[dict], None]] = None) -> List[dict]:
    """fly_v4.run_lockstep_v4 on the jobs; every record gets the v5 flight id and section-5 metrics."""
    by_v4 = {F.flight_id(j): j for j in jobs}
    if len(by_v4) != len(jobs):
        raise ValueError("run_mppi_v5: two jobs share a v4 flight id (fly one suite per batch)")

    def relabel(rec: dict) -> dict:
        res, tel = rec["result"], rec["telemetry"]
        job = by_v4[res["flight_id"]]
        res["v4_flight_id"] = res["flight_id"]
        res.update(result_extras_v5(job))
        res["source"] = "v5"
        res.update(v5_metrics(res, tel, 90.0, ()))
        # run_lockstep_v4 keeps the mean and max planning time per flight; RQ5's median and p99 per decision
        # come from the dedicated timing stage (study_v5.compute_table), not from batched flights
        res["plan_ms_median"], res["plan_ms_p99"] = res.get("plan_ms_mean"), res.get("plan_ms_max")
        return rec

    def done(rec):
        relabel(rec)
        if on_done is not None:
            on_done(rec)

    out = F4.run_lockstep_v4(jobs, adapter, mppi, n_sub=n_sub, max_rows=max_rows, on_done=done, verbose=verbose,
                             progress=progress)
    return out


class FilteredPhysicsAdapter:
    """
    mppi_filtered: the classical controller's physics model (physics_v4.PhysicsAdapter, nominal parameters)
    with the drag and lift scales of FNPEG's fading-memory filters, updated at every decision from the
    flight's latest measured accelerations (y[8], y[9]). The filter runs at the decision cadence with
    beta^(dt) so its time constant equals FNPEG's 1 Hz filter. frozen=True keeps K_D = K_L = 1, which makes
    the context identical to PhysicsAdapter("nominal")'s (smoke check 4).
    """
    M = 1
    rf = 2

    def __init__(self, device="cpu", beta: float = 0.9, on_g: float = 0.05, dt_s: float = 5.0, frozen: bool = False):
        from jepa import physics_v4 as PH
        from jepa.guidance_fnpeg import long_model
        self.base = PH.PhysicsAdapter("nominal", device)
        self.model = long_model()
        self.beta = float(beta) ** float(dt_s)
        self.on = float(on_g) * G0
        self.frozen = bool(frozen)
        self.K: Dict[str, List[float]] = {}
        self.dev = self.base.dev

    def context(self, y_hists, jobs):
        rows = []
        for y, j in zip(y_hists, jobs):
            key = F.flight_id(j)
            K = self.K.setdefault(key, [1.0, 1.0])
            yl = y[-1].double().cpu().numpy()
            aD, aL = float(yl[8]), float(yl[9])
            if not self.frozen and math.hypot(aD, aL) >= self.on:
                D0, L0, _ = self.model.aero(yl[0:1], yl[3:4])
                D0, L0 = float(D0[0]), float(L0[0])
                if D0 > 0 and L0 > 0:
                    K[0] = self.beta * K[0] + (1.0 - self.beta) * aD / D0
                    K[1] = self.beta * K[1] + (1.0 - self.beta) * aL / L0
            rows.append([float(y[-1, 6]), float(y[-1, 7]), K[0], K[1] / K[0], 1.0])
        return torch.tensor(rows, dtype=torch.float64)[None]

    def rollout(self, s0, y0, uh, u, n_sub, h_stop):
        return self.base.rollout(s0, y0, uh, u, n_sub, h_stop)


_ADAPTERS: Dict[tuple, object] = {}


def load_adapter_v5(kind: str, ckpt: Optional[str] = None, device: str = "cpu", members=None, chunk: int = 512,
                    keep_only: bool = False, **kw):
    """fly_v4.load_adapter_v4's kinds, plus physics_filtered(_frozen) and the E3 model kinds."""
    if kind in ("v1", "v2", "physics_nominal", "physics_truth", "prior"):
        return F4.load_adapter_v4(kind, ckpt, device, members, chunk, keep_only)
    key = (kind, str(ckpt), str(device))
    if keep_only:
        for k in [k for k in _ADAPTERS if k != key]:
            del _ADAPTERS[k]
    if key not in _ADAPTERS:
        if kind in ("physics_filtered", "physics_filtered_frozen"):
            _ADAPTERS[key] = FilteredPhysicsAdapter(device, frozen=kind.endswith("frozen"), **kw)
        else:
            from jepa import models_v5 as M5
            _ADAPTERS[key] = M5.load_adapter(kind, ckpt, device)
    if kind.startswith("physics_filtered"):
        _ADAPTERS[key].K = {}                   # filters start at 1 for every new batch of flights
    return _ADAPTERS[key]


# ============================================================================
# Matched limits from fnpeg flights
# ============================================================================
def fnpeg_peak_g(layout, suite: str, conds: Sequence[dict], arm: str = "fnpeg") -> Dict[tuple, float]:
    """(case, atmosphere) -> realized peak g of the fnpeg flight on that condition (must exist)."""
    import json
    out, missing = {}, []
    for c in conds:
        fid = f"{suite}/{arm}/{cond_key({**c, 'suite': suite})}"
        p = layout.result_path(fid)
        if p.exists():
            out[(int(c["case"]), c["atmosphere"])] = float(json.loads(p.read_text())["peak_g"])
        else:
            missing.append(fid)
    if missing:
        raise FileNotFoundError(f"matched limits need {len(missing)} fnpeg flights first, e.g. {missing[0]}")
    return out
