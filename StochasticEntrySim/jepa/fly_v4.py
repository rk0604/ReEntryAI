"""
jepa/fly_v4.py

Closed-loop flights for the v4 study (jepa/v4_plan.md sections 4.0-4.2, 4.4).

    run_lockstep_v4     fly_v3.run_lockstep_v3 with two changes only: the
                        adapter's context() also gets the deciding flights' jobs
                        (a physics model needs each flight's realized bank and,
                        for the oracle, its true parameters), and every decision
                        is logged (mppi.last, the command, the predicted peak g
                        and landing point of the best sample) into the flight's
                        telemetry .npz as dec_* arrays (section 4.4)
    V1AdapterV4, V2AdapterV4
                        fly_v3's adapters, accepting and ignoring the jobs
    load_adapter_v4     one adapter per (kind, checkpoint, members): "v1" (v1 and
                        the B variants), "v2" (v2, v2single, R0 and the A
                        variants), "physics_nominal", "physics_truth", "prior"
    matched_limits      the per-condition g limit of the matched arms: the
                        classical controller's peak g on that (case, atmosphere)
                        in v3's classical flights
    e0a_table, e0b_table, e0c_check
                        the reproducibility checks of section 4.0
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
import torch

from jepa import fly_v3 as F
from jepa.data import G0, MACRO_DT_S, NOMINAL_DOMAIN
from jepa.evaluate import crossing_point, haversine_km
from jepa.mppi_v2 import MPPIConfigV2

DEC_FIELDS = ("t_s", "h_m", "V_mps", "cmd_deg", "J_min", "miss_best_km", "spread_km", "ess", "k_rem",
              "peak_g_pred", "phi_pred_deg", "lambda_pred_deg")


# ============================================================================
# Adapters
# ============================================================================
class V2AdapterV4(F.V2Adapter):
    """fly_v3.V2Adapter; context() accepts and ignores the jobs."""

    def context(self, y_hists, jobs=None):
        return super().context(y_hists)


class V1AdapterV4(F.V1Adapter):
    """fly_v3.V1Adapter; context() accepts and ignores the jobs."""

    def context(self, y_hists, jobs=None):
        return super().context(y_hists)


def load_v1_family(ckpt: str) -> dict:
    """
    models.load_checkpoint, plus the v2 input features of B1: its checkpoint
    records feat_set "v2", and its stats are then StatsV2 (15 features).
    """
    from jepa.models import load_checkpoint
    from jepa.models_v2 import StatsV2
    m = load_checkpoint(ckpt)
    if (m.get("extra") or {}).get("feat_set") == "v2":
        ck = torch.load(ckpt, map_location="cpu", weights_only=False)
        m["stats"] = StatsV2(**ck["stats"])
    return m


_ADAPTERS: Dict[tuple, object] = {}


def load_adapter_v4(kind: str, ckpt: Optional[str] = None, device: str = "cpu",
                    members: Optional[Sequence[int]] = None, chunk: int = 512, keep_only: bool = False):
    """
    Build (once per process) the adapter of a model. keep_only=True first drops every other cached
    adapter (a long-lived pool worker flies dozens of models one after another).
    """
    key = (kind, str(ckpt), str(device), tuple(members) if members is not None else None, int(chunk))
    if keep_only and any(k != key for k in _ADAPTERS):
        for k in [k for k in _ADAPTERS if k != key]:
            del _ADAPTERS[k]
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
    if key not in _ADAPTERS:
        from jepa import physics_v4 as PH
        if kind == "v1":
            _ADAPTERS[key] = V1AdapterV4(load_v1_family(ckpt), device)
        elif kind == "v2":
            from jepa.models_v2 import load_checkpoint_v2, select_members
            m = load_checkpoint_v2(ckpt, device)
            _ADAPTERS[key] = V2AdapterV4(select_members(m, members) if members is not None else m, chunk)
        elif kind == "physics_nominal":
            _ADAPTERS[key] = PH.PhysicsAdapter("nominal", device)
        elif kind == "physics_truth":
            _ADAPTERS[key] = PH.PhysicsAdapter("truth", device)
        elif kind == "prior":
            _ADAPTERS[key] = PH.PriorAdapter(device)
        else:
            raise ValueError(f"unknown model kind {kind!r}")
    return _ADAPTERS[key]


# ============================================================================
# Lockstep runner with the decision log
# ============================================================================
class FlightV4(F.FlightV3):
    """fly_v3.FlightV3 + the per-decision log."""

    def __init__(self, job: dict, rf: int, cfg: MPPIConfigV2):
        super().__init__(job, rf, cfg)
        self.dec: List[tuple] = []

    def record(self) -> dict:
        rec = super().record()
        d = np.asarray(self.dec, dtype=np.float64).reshape(-1, len(DEC_FIELDS))
        rec["telemetry"].update({f"dec_{k}": d[:, i].astype(np.float32) for i, k in enumerate(DEC_FIELDS)})
        rec["result"].update(result_extras(self.job))
        return rec


def result_extras(job: dict) -> dict:
    """The v4 labels every result record carries next to fly_v3.finish's fields."""
    return {"arm": job.get("controller"), "stage": job.get("stage"), "variant": job.get("variant"),
            "model_kind": (job.get("model") or {}).get("kind"), "source": "v4"}


def _best_sample(f: FlightV4, o: Dict[str, torch.Tensor], u: torch.Tensor, u_prev: torch.Tensor) -> tuple:
    """
    The sample EnsembleMPPI.update will pick (same arithmetic: argmax of the
    softmax weights): its predicted peak g (each member's largest predicted
    load before its own drogue crossing, mean over members) and landing point.
    Pure function of the rollout: no random numbers are drawn.
    """
    c = f.mppi.cfg
    J, _, _, _ = f.mppi._cost(o, u, u_prev)
    J_min = J.min()
    temp = max(c.temp_kappa * float(torch.quantile(J, c.temp_quantile) - J_min), 1e-6)
    best = int(torch.argmax(torch.softmax(-(J - J_min) / temp, 0)))
    y6 = o["y6"][best]                                                 # [M, K'+1, 6]
    phi, lam, _, kc = crossing_point(y6, f.mppi.h_drogue)
    gl = o["gload"][best]                                              # [M, K']
    peaks = [float(gl[m, :max(1, int(kc[m]))].max()) for m in range(gl.shape[0])]
    return float(np.mean(peaks)), math.degrees(float(phi.mean())), math.degrees(float(lam.mean()))


@torch.no_grad()
def run_lockstep_v4(jobs: Sequence[dict], adapter, mppi: Optional[dict] = None, n_sub: int = 2,
                    max_rows: int = 2048, on_done: Optional[Callable[[dict], None]] = None,
                    verbose: bool = True, progress: Optional[Callable[[dict], None]] = None) -> List[dict]:
    """
    fly_v3.run_lockstep_v3 (same loop and arithmetic) with adapter.context(y_hists, jobs)
    and the per-decision log. Returns the records ({"result", "telemetry"}) in job order.
    """
    flights = [FlightV4(j, adapter.rf, MPPIConfigV2(**{**(mppi or {}), **j.get("mppi", {})})) for j in jobs]
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
            s_ctx = adapter.context([f.y_hist for f in deciding], [f.job for f in deciding])   # [M, E, D]
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
                        o = {"y6": y6[sl], "gload": gl[sl]}
                        final = it == iters[i] - 1
                        best = _best_sample(f, o, us[i], f.u_hist[-1:]) if final else None
                        cmd = f.mppi.update(o, us[i], f.u_hist[-1:], final=final)
                        if cmd is not None:
                            f.commit(cmd)
                            L = f.mppi.last
                            f.dec.append((float(f.env.t_s), float(f._y_now[0]), float(f._y_now[3]), math.degrees(cmd),
                                          L.get("J_min", np.nan), L.get("miss_best_km", np.nan),
                                          L.get("spread_km", np.nan), L.get("ess", np.nan), L.get("k_rem", np.nan),
                                          *best))
            ms = 1000.0 * (time.perf_counter() - tr) / len(deciding)
            for f in deciding:
                f.plan_ms.append(ms)
            n_round += 1
            if progress is not None and n_round % 5 == 0:
                progress({"round": n_round, "flying": len(active), "total": len(flights),
                          "sec": round(time.time() - t0, 1), "ms_per_decision": round(ms, 1)})
            if verbose and n_round % 25 == 0:
                print(f"  [lockstep v4] round {n_round}: {len(active)}/{len(flights)} flying, "
                      f"{ms:.0f} ms/decision, {time.time() - t0:.0f}s", flush=True)
        for i in active:
            f = flights[i]
            f.env_step()
            if f.done:
                out[i] = f.record()
                if on_done is not None:
                    on_done(out[i])
    return out  # type: ignore[return-value]


def fly_classical_v4(job: dict) -> dict:
    """fly_v3.fly_classical (unchanged guidance) with the v4 labels on the result."""
    rec = F.fly_classical(job)
    rec["result"].update(result_extras(job))
    return rec


# ============================================================================
# Matched load limits
# ============================================================================
def classical_peak_g(v3_root, cases: Sequence[int], atmospheres: Sequence[str]) -> Dict[tuple, float]:
    """(case, atmosphere) -> peak g of v3's classical flight (results/classical/<atm>_case<NN>.json)."""
    out, missing = {}, []
    for c in cases:
        for a in atmospheres:
            p = Path(v3_root) / "04_flights" / "results" / "classical" / f"{a}_case{int(c):02d}.json"
            if p.exists():
                out[(int(c), a)] = float(json.loads(p.read_text())["peak_g"])
            else:
                missing.append(p.name)
    if missing:
        raise FileNotFoundError(f"v3 classical flights missing for the matched load: {missing[:6]}"
                                f"{' ...' if len(missing) > 6 else ''} ({len(missing)} in all)")
    return out


# ============================================================================
# E0: reproducibility checks
# ============================================================================
def _read(path: Path) -> Optional[dict]:
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def e0a_table(v3_root, v4_root, cases: Sequence[int], atmospheres: Sequence[str]):
    """
    E0a: v3's classical flights re-flown (arm classical_E0a) against v3's own,
    plus classical_oracle against classical in the nominal atmosphere.
    """
    import pandas as pd
    rows = []
    for c in cases:
        for a in atmospheres:
            name = f"{a}_case{int(c):02d}.json"
            r3 = _read(Path(v3_root) / "04_flights" / "results" / "classical" / name)
            r4 = _read(Path(v4_root) / "04_flights" / "results" / "classical_E0a" / name)
            if r3 is None or r4 is None:
                continue
            rows.append({"check": "E0a", "case": int(c), "atmosphere": a, "outcome_v3": r3["outcome"],
                         "outcome_v4": r4["outcome"], "miss_v3_km": r3["miss_km"], "miss_v4_km": r4["miss_km"],
                         "d_miss_km": abs(r4["miss_km"] - r3["miss_km"]), "fuel_v3_kg": r3["rcs_fuel_kg"],
                         "fuel_v4_kg": r4["rcs_fuel_kg"], "d_fuel_kg": abs(r4["rcs_fuel_kg"] - r3["rcs_fuel_kg"]),
                         "identical": bool(r3["outcome"] == r4["outcome"] and abs(r4["miss_km"] - r3["miss_km"]) <= 1e-6
                                           and r4["rcs_fuel_kg"] == r3["rcs_fuel_kg"])})
    for c in cases:
        r3 = _read(Path(v3_root) / "04_flights" / "results" / "classical" / f"nominal_case{int(c):02d}.json")
        ro = _read(Path(v4_root) / "04_flights" / "results" / "classical_oracle" / f"nominal_case{int(c):02d}.json")
        if r3 is None or ro is None:
            continue
        rows.append({"check": "oracle_vs_classical_nominal", "case": int(c), "atmosphere": "nominal",
                     "outcome_v3": r3["outcome"], "outcome_v4": ro["outcome"], "miss_v3_km": r3["miss_km"],
                     "miss_v4_km": ro["miss_km"], "d_miss_km": abs(ro["miss_km"] - r3["miss_km"]),
                     "fuel_v3_kg": r3["rcs_fuel_kg"], "fuel_v4_kg": ro["rcs_fuel_kg"],
                     "d_fuel_kg": abs(ro["rcs_fuel_kg"] - r3["rcs_fuel_kg"]),
                     "identical": bool(r3["outcome"] == ro["outcome"] and abs(ro["miss_km"] - r3["miss_km"]) <= 1e-6
                                       and ro["rcs_fuel_kg"] == r3["rcs_fuel_kg"])})
    return pd.DataFrame(rows)


def e0b_table(v3_root, v4_root, pairs: Sequence[tuple], seed: int, cases: Sequence[int], atmospheres: Sequence[str],
              load: str = "4.5g"):
    """E0b: per condition, v3's learned flight against the v4 runner's re-fly. pairs: (v3 arm, v4 arm)."""
    import pandas as pd
    rows = []
    for a3, a4 in pairs:
        for c in cases:
            for a in atmospheres:
                rel = f"{load}/seed{int(seed)}/mppi0/{a}_case{int(c):02d}.json"
                r3 = _read(Path(v3_root) / "04_flights" / "results" / a3 / rel)
                r4 = _read(Path(v4_root) / "04_flights" / "results" / a4 / rel)
                if r3 is None or r4 is None:
                    continue
                l3, l4 = math.log(max(r3["miss_km"], 0.01)), math.log(max(r4["miss_km"], 0.01))
                rows.append({"v3_arm": a3, "v4_arm": a4, "seed": int(seed), "case": int(c), "atmosphere": a,
                             "outcome_v3": r3["outcome"], "outcome_v4": r4["outcome"],
                             "same_outcome": r3["outcome"] == r4["outcome"],
                             "same_success5": bool(r3["success_5km"]) == bool(r4["success_5km"]),
                             "miss_v3_km": r3["miss_km"], "miss_v4_km": r4["miss_km"], "d_log_miss": l4 - l3})
    return pd.DataFrame(rows)


def e0c_check(jobs: Sequence[dict], v1_ckpt: str, v2_ckpt: str, mppi: dict) -> List[dict]:
    """
    E0c: run_lockstep_v4 against run_lockstep_v3 on the CPU, for the v1 adapter
    and for v2's member 0, on the same jobs: every result field (except wall
    time and planning time) and every telemetry array v3 records must be equal.
    """
    from jepa.models import load_checkpoint
    from jepa.models_v2 import load_checkpoint_v2, select_members
    out = []
    n_threads = torch.get_num_threads()
    torch.set_num_threads(1)
    try:
        pairs = [("v1", F.V1Adapter(load_checkpoint(v1_ckpt), "cpu"), V1AdapterV4(load_checkpoint(v1_ckpt), "cpu")),
                 ("v2 member 0", F.V2Adapter(select_members(load_checkpoint_v2(v2_ckpt, "cpu"), [0])),
                  V2AdapterV4(select_members(load_checkpoint_v2(v2_ckpt, "cpu"), [0])))]
        runs = [(kind, F.run_lockstep_v3(jobs, a3, mppi, verbose=False), run_lockstep_v4(jobs, a4, mppi, verbose=False))
                for kind, a3, a4 in pairs]
    finally:
        torch.set_num_threads(n_threads)
    for kind, r3, r4 in runs:
        for j, x, y in zip(jobs, r3, r4):
            skip = {"wall_s", "plan_ms_mean", "plan_ms_max"}
            bad = [k for k, v in x["result"].items() if k not in skip and
                   not (v == y["result"].get(k) or (isinstance(v, float) and math.isnan(v) and
                                                    math.isnan(y["result"].get(k, 0.0))))]
            bad += [k for k, v in x["telemetry"].items() if not np.array_equal(v, y["telemetry"][k])]
            out.append({"adapter": kind, "case": j["case"], "atmosphere": j["atmosphere"],
                        "miss_v3_km": x["result"]["miss_km"], "miss_v4_km": y["result"]["miss_km"],
                        "decisions": int(len(y["telemetry"]["dec_t_s"])), "identical": not bad,
                        "differs_in": ",".join(bad)})
    return out
