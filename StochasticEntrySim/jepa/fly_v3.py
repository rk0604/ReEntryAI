"""
jepa/fly_v3.py

Closed-loop flights for the v3 study (jepa/v3_plan.md sections 4-6).

    Recorder            1 Hz telemetry read from the env after every step:
                        state, commanded and actual bank, roll rate, RCS fuel
                        and jet-on seconds, g-load, heat rate. No simulator change.
    flight_summary      per-flight metrics from the telemetry (heat load, fuel
                        by phase, bank tracking, roll-rate saturation, ...)
    crn_seed            common random numbers: the MPPI noise stream of a
                        flight depends only on (seed, case, atmosphere, rep),
                        so every learned controller sees the same noise
    V1Adapter/V2Adapter one interface over the v1 (TCN + GRU) and v2
                        (ensemble transformer) world models, so all learned
                        controllers fly through the same lockstep MPPI runner
    run_lockstep_v3     jepa.mppi_v2.run_lockstep with a model adapter,
                        telemetry, and an on_done callback per flight
    fly_classical       jepa.mppi.run_episode's classical branch + telemetry
    flight_id           the job key -> results/<controller>/<load>/seed<k>/mppi<m>/<atm>_case<NN>

Every flight's result dict is flat and JSON-serializable; its telemetry is a
dict of float32 arrays saved as .npz next to it (see study_v3.write_flight).
"""

from __future__ import annotations

import math
import time
import zlib
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
import torch

from jepa.data import G0, MACRO_DT_S, NOMINAL_DOMAIN, SIGMA_MAX_RAD, apply_domain, physical_state
from jepa.mppi_v2 import EnsembleMPPI, MPPIConfigV2, _Flight

SIGMA_RATE_MAX_DPS = 20.0          # rl_env: sigma_rate_max_rps = radians(20)
ROLL_SAT_FRAC = 0.9                # "saturated" = |roll rate| >= 90% of the limit
PHASE_EDGES_M = (60_000.0, 25_000.0)
TELEM_FIELDS = ("t_s", "h_m", "V_mps", "gamma_deg", "phi_deg", "lambda_deg", "chi_deg", "sigma_cmd_deg",
                "sigma_deg", "roll_rate_dps", "rcs_fuel_kg", "g", "qdot_MWm2")


# ============================================================================
# Keys and seeds
# ============================================================================
def flight_id(job: dict) -> str:
    """Relative path (no extension) of one flight's result and telemetry files."""
    case = f"{job['atmosphere']}_case{int(job['case']):02d}"
    if job["controller"].startswith("classical"):
        return f"{job['controller']}/{case}"
    return f"{job['controller']}/{job['load']}/seed{int(job['seed'])}/mppi{int(job.get('mppi_rep', 0))}/{case}"


def crn_seed(seed: int, case: int, atmosphere: str, rep: int = 0) -> int:
    """MPPI noise seed shared by every learned controller on the same (seed, case, atmosphere, rep)."""
    return zlib.crc32(f"v3|{int(seed)}|{int(case)}|{atmosphere}|{int(rep)}".encode()) & 0x7FFFFFFF


def great_circle_km(phi1: float, lam1: float, phi2: float, lam2: float) -> float:
    """Haversine distance on the simulator's sphere (km)."""
    import constants
    R = float(constants.RADIUS_EARTH) / 1000.0
    a = math.sin((phi2 - phi1) / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin((lam2 - lam1) / 2) ** 2
    return 2.0 * R * math.asin(min(1.0, math.sqrt(a)))


# ============================================================================
# Telemetry and per-flight metrics
# ============================================================================
class Recorder:
    """1 Hz telemetry, one row at entry interface and one after every env.step."""

    def __init__(self, env):
        import telemetry
        self.kg_per_jet_s = float(telemetry.rcs_propellant_kg(1.0))     # linear in on-time (exact)
        self.rows: List[tuple] = []
        self.add(env, 0.0)

    def add(self, env, sigma_cmd_rad: float) -> None:
        """Append the env's current state (call right after env.step)."""
        import constants
        r, phi, lam, V, gam, chi = (float(v) for v in env.x)
        d = math.degrees
        self.rows.append((float(env.t_s), r - float(constants.RADIUS_EARTH), V, d(gam), d(phi), d(lam), d(chi),
                          d(sigma_cmd_rad), d(float(env.att.sigma_rel_rad)), d(float(env.att.omega_b_rad_s[2])),
                          float(env.rcs_fuel_kg), float(env._last_g), float(env._last_qdot) / 1e6))

    def arrays(self) -> Dict[str, np.ndarray]:
        """Telemetry columns as float32 arrays (+ cumulative jet-on seconds)."""
        a = np.asarray(self.rows, dtype=np.float64)
        out = {k: a[:, i].astype(np.float32) for i, k in enumerate(TELEM_FIELDS)}
        out["jet_s"] = (a[:, TELEM_FIELDS.index("rcs_fuel_kg")] / self.kg_per_jet_s).astype(np.float32)
        return out


def _wrap_deg(x: np.ndarray) -> np.ndarray:
    return (x + 180.0) % 360.0 - 180.0


def flight_summary(tel: Dict[str, np.ndarray]) -> Dict[str, float]:
    """v3_plan.md section 5/6 metrics that come from the 1 Hz telemetry."""
    t = tel["t_s"].astype(np.float64)
    q = tel["qdot_MWm2"].astype(np.float64)
    fuel = tel["rcs_fuel_kg"].astype(np.float64)
    h0 = tel["h_m"][:-1]
    dfuel = np.diff(fuel)
    hi, lo = PHASE_EDGES_M
    rr = np.abs(tel["roll_rate_dps"][1:].astype(np.float64))
    err = _wrap_deg(tel["sigma_cmd_deg"][1:].astype(np.float64) - tel["sigma_deg"][1:].astype(np.float64))
    n = max(len(rr), 1)
    return {
        "heat_load_MJm2": float(((q[1:] + q[:-1]) * 0.5 * np.diff(t)).sum()),
        "jet_s": float(tel["jet_s"][-1]),
        "fuel_kg_above60km": float(dfuel[h0 > hi].sum()),
        "fuel_kg_60_25km": float(dfuel[(h0 <= hi) & (h0 > lo)].sum()),
        "fuel_kg_below25km": float(dfuel[h0 <= lo].sum()),
        "roll_rate_max_dps": float(rr.max()) if len(rr) else 0.0,
        "roll_sat_frac": float((rr >= ROLL_SAT_FRAC * SIGMA_RATE_MAX_DPS).sum() / n),
        "bank_track_rms_deg": float(np.sqrt((err ** 2).mean())) if len(err) else 0.0,
        "telemetry_rows": int(len(t)),
    }


def n_reversals(cmds_rad: Sequence[float]) -> int:
    """Sign changes of the bank command, counting only |sigma| > 5 deg (as v1/v2)."""
    c = np.asarray(cmds_rad, dtype=np.float64)
    s = np.sign(c[np.abs(c) > math.radians(5.0)])
    return int((s[1:] != s[:-1]).sum()) if len(s) > 1 else 0


def finish(job: dict, env, info: dict, cmds: Sequence[float], tel: Dict[str, np.ndarray],
           wall_s: float, plan_ms: Sequence[float], spread_km: Sequence[float] = ()) -> dict:
    """The flat result record of one flight (JSON-serializable)."""
    miss = great_circle_km(float(env.x[1]), float(env.x[2]), float(env.target_phi_rad), float(env.target_lam_rad))
    outcome = str(info.get("outcome"))
    keep = ("controller", "load", "g_limit", "seed", "mppi_rep", "mppi_seed", "case", "set", "atmosphere",
            "model_id", "horizon")
    return {
        "flight_id": flight_id(job), **{k: job.get(k) for k in keep},
        "outcome": outcome, "miss_km": miss,
        "success_5km": bool(outcome == "drogue" and miss < 5.0),
        "success_25km": bool(outcome == "drogue" and miss < 25.0),
        "peak_g": float(info.get("peak_g", float("nan"))),
        "peak_qdot_MWm2": float(info.get("peak_qdot_MWm2", float("nan"))),
        "rcs_fuel_kg": float(info.get("rcs_fuel_kg", float("nan"))),
        "t_s": float(info.get("t_s", float("nan"))),
        **flight_summary(tel),
        "n_reversals": n_reversals(cmds),
        "plan_ms_mean": float(np.mean(plan_ms)) if len(plan_ms) else float("nan"),
        "plan_ms_max": float(np.max(plan_ms)) if len(plan_ms) else float("nan"),
        "spread_km_mean": float(np.nanmean(spread_km)) if len(spread_km) else float("nan"),
        "wall_s": round(float(wall_s), 2),
        **{f"dom_{k}": float(v) for k, v in (job.get("domain") or NOMINAL_DOMAIN).items()},
    }


# ============================================================================
# Classical predictor-corrector (unchanged guidance, plus telemetry)
# ============================================================================
def fly_classical(job: dict) -> dict:
    """
    One classical flight, exactly as jepa.mppi.run_episode flies it (baseline
    mode, dispersion off, truth applied to the env, the guidance keeps the
    nominal aero unless job["classical_knows_truth"]). Returns
    {"result": record, "telemetry": arrays}. Pool-worker safe.
    """
    torch.set_num_threads(int(job.get("threads", 1)))
    import rl_env
    env = rl_env.OrionEntryEnv(config_path=job.get("config_path", "configs/default.json"),
                               controller_mode="baseline", dispersion=False, max_episode_steps=1400,
                               fast_atmosphere=True)
    domain = dict(job.get("domain") or NOMINAL_DOMAIN)
    apply_domain(env.params, domain)
    if job.get("classical_knows_truth"):
        apply_domain(env.guidance.params, domain)
    env.reset(seed=int(job.get("env_seed", 0)), options={"ic_override": job["ic"]} if job.get("ic") else None)
    rec = Recorder(env)
    t0 = time.time()
    cmds: List[float] = []
    step_ms: List[float] = []
    done, info = False, {}
    while not done:
        ts = time.perf_counter()
        _, _, term, trunc, info = env.step(np.array([0.0], dtype=np.float32))     # action ignored (baseline)
        step_ms.append(1000.0 * (time.perf_counter() - ts))
        sigma = math.radians(float((env.guidance.get_last_debug() or {}).get("chosen_sigma_cmd_deg", 0.0)))
        cmds.append(sigma)
        rec.add(env, sigma)
        done = bool(term or trunc)
    tel = rec.arrays()
    # plan_ms for the classical controller is the whole 1 s guidance step (predictor + physics).
    return {"result": finish(job, env, info, cmds, tel, time.time() - t0, step_ms[::5]), "telemetry": tel}


# ============================================================================
# Model adapters: one interface for v1 and v2 inside the lockstep runner
# ============================================================================
class V2Adapter:
    """The v2 ensemble (or a member subset) on its device, with AMP."""

    def __init__(self, models: dict, chunk: int = 512):
        from jepa.models_v2 import amp_dtype
        self.wm, self.pi, self.stats = models["wm"], models["pi"], models["stats"]
        self.rf, self.M = self.wm.rf, self.wm.M
        self.dev = next(self.wm.parameters()).device
        self.amp = amp_dtype(self.dev)
        self.chunk = int(chunk)

    def context(self, y_hists: Sequence[torch.Tensor]) -> torch.Tensor:
        """Encoder histories [rf, NY] per flight -> context latents [M, E, D] (fp32)."""
        feats = torch.stack([self.stats.norm_feats(y) for y in y_hists]).to(self.dev)
        with torch.autocast(self.dev.type, dtype=self.amp, enabled=self.amp is not None):
            return self.wm.encode_states(feats)[:, :, -1].float()

    def rollout(self, s0: torch.Tensor, y0: torch.Tensor, uh: torch.Tensor, u: torch.Tensor,
                n_sub: int, h_stop: float):
        """Samples -> (y6 [B, M, K'+1, 6], gload [B, M, K']) on the CPU."""
        from jepa.models_v2 import ensemble_rollout
        with torch.autocast(self.dev.type, dtype=self.amp, enabled=self.amp is not None):
            out = ensemble_rollout(self.wm, self.pi, s0, y0.to(self.dev), uh.to(self.dev), u.to(self.dev),
                                   MACRO_DT_S, n_sub, h_stop=h_stop, chunk=self.chunk)
        return out["y6"].cpu(), (out["acc"].norm(dim=-1) / G0).cpu()


class V1Adapter:
    """
    The v1 model (models.JEPAWorldModel + PIProber) as a one-member ensemble.
    rollout() is models.jepa_predict's loop with the context latent passed in
    per row (so many flights share one batch), in fp32 on any device.
    """
    M = 1

    def __init__(self, models: dict, device="cpu"):
        self.dev = torch.device(device)
        self.wm, self.pi, self.stats = models["wm"].to(self.dev).eval(), models["pi"].to(self.dev).eval(), models["stats"]
        self.rf = self.wm.rf

    def context(self, y_hists: Sequence[torch.Tensor]) -> torch.Tensor:
        feats = torch.stack([self.stats.norm_feats(y) for y in y_hists]).to(self.dev)
        return self.wm.encode_states(feats)[:, -1][None].float()           # [1, E, D]

    def rollout(self, s0: torch.Tensor, y0: torch.Tensor, uh: torch.Tensor, u: torch.Tensor,
                n_sub: int, h_stop: float):
        u = u.to(self.dev)
        y = y0.to(self.dev).float()
        z = self.wm.encode_actions(torch.cat([uh.to(self.dev), u], 1))
        hid = self.wm.init_hidden(s0[0])
        ys, accs = [y], []
        for k in range(u.shape[1]):
            s, hid = self.wm.step(z[:, k], hid)
            y, a = self.pi.step(s, y, u[:, k], MACRO_DT_S, n_sub)
            ys.append(y)
            accs.append(a)
            if k % 5 == 4 and bool((y[:, 0] < h_stop).all()):
                break
        y6 = torch.stack(ys, 1)[:, None]
        gl = (torch.stack(accs, 1).norm(dim=-1) / G0)[:, None]
        return y6.cpu(), gl.cpu()


_ADAPTERS: Dict[tuple, object] = {}


def load_adapter(kind: str, ckpt: str, device: str = "cpu", members: Optional[Sequence[int]] = None,
                 chunk: int = 512):
    """Build (once per process) the adapter of a checkpoint: kind 'v1' or 'v2'."""
    key = (kind, str(ckpt), str(device), tuple(members) if members is not None else None, int(chunk))
    if key not in _ADAPTERS:
        if kind == "v1":
            from jepa.models import load_checkpoint
            _ADAPTERS[key] = V1Adapter(load_checkpoint(ckpt), device)
        else:
            from jepa.models_v2 import load_checkpoint_v2, select_members
            m = load_checkpoint_v2(ckpt, device)
            _ADAPTERS[key] = V2Adapter(select_members(m, members) if members is not None else m, chunk)
    return _ADAPTERS[key]


# ============================================================================
# Lockstep runner (all learned controllers)
# ============================================================================
class FlightV3(_Flight):
    """mppi_v2._Flight + 1 Hz telemetry."""

    def __init__(self, job: dict, rf: int, cfg: MPPIConfigV2):
        super().__init__(job, rf, cfg)
        self.rec = Recorder(self.env)

    def env_step(self) -> None:
        super().env_step()
        self.rec.add(self.env, self.sigma)

    def record(self) -> dict:
        tel = self.rec.arrays()
        return {"result": finish(self.job, self.env, self.info, self.cmds, tel, time.time() - self.t0,
                                 self.plan_ms, self.spread), "telemetry": tel}


@torch.no_grad()
def run_lockstep_v3(jobs: Sequence[dict], adapter, mppi: Optional[dict] = None, n_sub: int = 2,
                    max_rows: int = 2048, on_done: Optional[Callable[[dict], None]] = None,
                    verbose: bool = True, progress: Optional[Callable[[dict], None]] = None) -> List[dict]:
    """
    Fly every job with MPPI on the adapter's model, all flights in lockstep
    (same loop and arithmetic as mppi_v2.run_lockstep). Job keys: ic, domain,
    controller/load/seed/mppi_rep/case/atmosphere labels, and job["mppi"]
    overrides (its "seed" is the flight's CRN noise seed). on_done(record)
    is called as soon as each flight lands, so results persist mid-batch;
    progress(dict) is called every 5 decision rounds.
    Returns the records ({"result", "telemetry"}) in job order.
    """
    flights = [FlightV3(j, adapter.rf, MPPIConfigV2(**{**(mppi or {}), **j.get("mppi", {})})) for j in jobs]
    out: List[Optional[dict]] = [None] * len(flights)
    t0, n_round = time.time(), 0
    while True:
        active = [i for i, f in enumerate(flights) if not f.done]
        if not active:
            break
        deciding = [flights[i] for i in active if flights[i].k % flights[i].n_sub == 0]
        if deciding:
            tr = time.perf_counter()
            for f in deciding:
                f.observe()
            s_ctx = adapter.context([f.y_hist for f in deciding])                  # [M, E, D]
            iters = [f.mppi.begin() for f in deciding]
            for it in range(max(iters)):
                group = [i for i, f in enumerate(deciding) if it < iters[i]]
                us = {i: deciding[i].mppi.propose() for i in group}
                S = us[group[0]].shape[0]
                per = max(1, max_rows // S)
                for g0 in range(0, len(group), per):
                    sub = group[g0:g0 + per]
                    K = deciding[sub[0]].mppi.cfg.horizon
                    u = torch.cat([us[i] for i in sub])[:, :K]
                    s0 = s_ctx[:, sub].repeat_interleave(S, 1)
                    y0 = torch.stack([deciding[i].y_hist[-1, :6] for i in sub]).repeat_interleave(S, 0)
                    uh = torch.stack([deciding[i].u_hist for i in sub]).repeat_interleave(S, 0)
                    y6, gl = adapter.rollout(s0, y0, uh, u, n_sub, float(deciding[sub[0]].env.drogue_alt_m))
                    for j, i in enumerate(sub):
                        f = deciding[i]
                        sl = slice(j * S, (j + 1) * S)
                        cmd = f.mppi.update({"y6": y6[sl], "gload": gl[sl]}, us[i], f.u_hist[-1:],
                                            final=it == iters[i] - 1)
                        if cmd is not None:
                            f.commit(cmd)
            ms = 1000.0 * (time.perf_counter() - tr) / len(deciding)
            for f in deciding:
                f.plan_ms.append(ms)
            n_round += 1
            if progress is not None and n_round % 5 == 0:
                progress({"round": n_round, "flying": len(active), "total": len(flights),
                          "sec": round(time.time() - t0, 1), "ms_per_decision": round(ms, 1)})
            if verbose and n_round % 25 == 0:
                print(f"  [lockstep] round {n_round}: {len(active)}/{len(flights)} flying, "
                      f"{ms:.0f} ms/decision, {time.time() - t0:.0f}s", flush=True)
        for i in active:
            f = flights[i]
            f.env_step()
            if f.done:
                out[i] = f.record()
                if on_done is not None:
                    on_done(out[i])
    return out  # type: ignore[return-value]
