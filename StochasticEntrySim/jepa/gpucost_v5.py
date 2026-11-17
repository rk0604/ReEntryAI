"""
What v5's GPU work (everything except RL) costs on the current runtime: E3 training, MPPI flights, open loop.

Times short pieces of the real study code with the study's frozen settings, then projects the work the study
folder still has to do (read only: its status, nothing is written there):

    E3 training  each E3 arm trains `train_steps` steps of its real training loop (the curriculum runs by fraction
                 of the budget, so the short run sees the same mix of 24/48/96-step windows as a full run); its
                 seconds per step x the full budget = one model. An arm still needing models also needs its
                 learning-rate pilots (len(lr_mults) x lr_pilot_frac of a full run).
    MPPI flights the study's run_mppi on the study's number of GPU workers and lockstep sizes, one model run per
                 model kind: mppi_filtered (physics), v2single (v3's v2, one member), the three E3 arms and
                 r0_nojepa (the short-trained models), on A256 conditions; mppi_filtered and v2single also on
                 ORB-1 conditions. Each model run's seconds come from its flights' result-file times (the pool's
                 start-up is paid by a small warm-up run first), divided by its flights.
    open loop    one footprint evaluation per E3 arm (open_loop_seed) x the seeds left.

Writes only into its own scratch folder; ReEntryAI_JEPA_v5_gpu_cost.ipynb saves the result to Drive. The short
models fly worse than trained ones, which changes flight lengths somewhat, so treat MPPI numbers as +-20 %.
"""

from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import torch

from jepa import fly_v5 as F5
from jepa import models_v5 as M5
from jepa import study_v5 as S5

GPU_TYPES = ("mppi_physics", "mppi_v2", "mppi_e3")


def units_per_hour(gpu_name: str) -> float:
    for k, v in S5.UNITS_PER_HOUR.items():
        if k != "cpu" and k.lower() in gpu_name.lower():
            return float(v)
    return float(S5.UNITS_PER_HOUR["cpu"])


def gpu_workers(device: str) -> int:
    """The study notebook's GPU worker count (build_notebook_v5.py SETUP)."""
    if device != "cuda":
        return 1
    mem = torch.cuda.get_device_properties(0).total_memory / 1e9
    return max(1, min(4, os.cpu_count() or 2, int((mem - 2) // 6)))


def _full_steps(cfg: S5.StudyConfigV5) -> int:
    return int(M5.E3TrainCfg(**cfg.e3_train).steps)


def measure_training(cfg: S5.StudyConfigV5, L: S5.LayoutV5, data: dict, device: str, steps: int,
                     log=print) -> Dict[str, dict]:
    """Seed 0 of every E3 arm for `steps` steps into the scratch folder (also the models the MPPI timing flies)."""
    out = {}
    for arm in cfg.e3_arms:
        t0 = time.time()
        info = M5.train_arm_seed(cfg, data, L, arm, 0, device, 1.0, log, steps=steps)
        wall = time.time() - t0
        s = float(info.get("train_s", wall))
        out[arm] = {"steps": int(steps), "train_s": round(s, 1), "wall_s": round(wall, 1),
                    "s_per_step": s / max(steps, 1), "params": info.get("params")}
        log(f"[train {arm}] {steps} steps in {s:.0f} s = {steps / max(s, 1e-9):.1f} steps/s; full run "
            f"({_full_steps(cfg):,} steps) ~ {s / max(steps, 1) * _full_steps(cfg) / 60:.0f} min")
    return out


def _mtimes(L: S5.LayoutV5, jobs: Sequence[dict]) -> List[float]:
    return sorted(os.stat(L.result_path(F5.flight_id_v5(j))).st_mtime for j in jobs
                  if L.result_path(F5.flight_id_v5(j)).exists())


def measure_mppi(cfg: S5.StudyConfigV5, L: S5.LayoutV5, L3, device: str, workers: int, n_model: int = 0,
                 n_physics: int = 0, log=print) -> Dict[str, dict]:
    """
    One model run per (kind, suite) through the study's run_mppi. n_model / n_physics flights per run (0 = one
    lockstep chunk per worker: max_rows / n_samples flights for a model, max_rows_physics / n_samples for physics).
    """
    S = int(cfg.mppi()["n_samples"])
    n_model = int(n_model or workers * max(1, cfg.max_rows // S))
    n_physics = int(n_physics or workers * max(1, cfg.max_rows_physics // S))
    a256, orb1 = S5.conditions(cfg, "A256"), S5.conditions(cfg, "ORB1")
    s0 = cfg.seeds[0]

    def pick(conds, n):
        step = max(1, len(conds) // max(n, 1))          # spread over the suite (cases and atmospheres)
        return (conds[::step] + conds)[:n]

    ds1 = S5.conditions(cfg, "DS1")
    runs = [("warmup", [S5.arm_job(cfg, L, L3, c, "mppi_filtered", "DS1", s0) for c in ds1[:2]])]
    runs.append(("mppi_filtered@A256", [S5.arm_job(cfg, L, L3, c, "mppi_filtered", "A256", s0)
                                        for c in pick(a256, n_physics)]))
    runs.append(("mppi_filtered@ORB1", [S5.arm_job(cfg, L, L3, c, "mppi_filtered", "ORB1", s0)
                                        for c in pick(orb1, n_physics)]))
    runs.append(("v2single@A256", [S5.arm_job(cfg, L, L3, c, "v2single", "A256", s0) for c in pick(a256, n_model)]))
    runs.append(("v2single@ORB1", [S5.arm_job(cfg, L, L3, c, "v2single", "ORB1", s0) for c in pick(orb1, n_model)]))
    for arm in cfg.e3_arms:
        runs.append((f"{arm}@A256", [S5.arm_job(cfg, L, L3, c, arm, "E3", 0) for c in pick(a256, n_model)]))
    jobs = [j for _, g in runs for j in g]
    t0 = time.time()
    S5.run_mppi(cfg, jobs, L, device, workers, log, "measure mppi", report_every_s=120.0)
    total = time.time() - t0
    out, prev_end = {}, None
    for name, g in runs:
        t = _mtimes(L, g)
        if len(t) < len(g):
            raise RuntimeError(f"measure mppi: {name}: {len(t)}/{len(g)} flights flown")
        if prev_end is not None:
            s = t[-1] - prev_end
            out[name] = {"flights": len(g), "s": round(s, 1), "s_per_flight": s / len(g),
                         "kind": g[0]["model"]["kind"]}
            log(f"[{name}] {len(g)} flights in {s:.0f} s = {s / len(g):.1f} s per flight ({workers} GPU workers)")
        prev_end = t[-1]
    out["_wall_s"] = round(total, 1)
    return out


def measure_open_loop(cfg: S5.StudyConfigV5, L: S5.LayoutV5, data: dict, device: str, log=print) -> Dict[str, dict]:
    out = {}
    for arm in cfg.e3_arms:
        t0 = time.time()
        M5.open_loop_seed(cfg, L, data, arm, 0, device)
        out[arm] = {"s": round(time.time() - t0, 1)}
        log(f"[open loop {arm}] {out[arm]['s']:.0f} s")
    return out


def study_left(cfg: S5.StudyConfigV5, study_root: Optional[str], L3, scratch: S5.LayoutV5, log=print) -> dict:
    """
    What the study folder still has to do, read only: its status (flights per stage and arm), models trained,
    learning rates chosen, open-loop files. No study folder: everything is left (a fresh study).
    """
    Ls, src = scratch, "fresh study (no study folder found)"
    if study_root and Path(study_root, "00_study", "study_config.json").exists():
        Ls, src = S5.LayoutV5(study_root), study_root

        def refuse(path):                 # read only: any write into the study folder raises
            raise S5.ReadOnlyViolation(f"gpu cost measurement refuses to write into the study folder: {path}")
        Ls.guard = refuse
    J = S5.all_stage_jobs(cfg, Ls, L3)
    st = S5.status(cfg, Ls, J, "T4", write=False)
    matched_missing = not S5.jobs_matched(cfg, Ls, L3)
    rows = [r for r in st["arms"] if r["type"] in GPU_TYPES]
    if matched_missing:                   # v2single@fnpeg_matched appears once fnpeg has flown A256
        n = len(S5.conditions(cfg, "A256")) * len(cfg.seeds)
        rows.append({"stage": "A256", "suite": "A256", "arm": "v2single", "load": "fnpeg_matched", "type": "mppi_v2",
                     "done": 0, "total": n, "note": "added: waits for fnpeg's A256 flights"})
    models = {a: [s for s in cfg.seeds if not Ls.ckpt(a, s).exists()] for a in cfg.e3_arms}
    pilots = {a: bool(models[a]) and not (Ls.models / a / "lr_choice.json").exists() for a in cfg.e3_arms}
    open_loop = {a: [s for s in cfg.seeds if not (Ls.open_loop / f"footprint_{a}_seed{s}.csv").exists()]
                 for a in cfg.e3_arms}
    cpu_h = sum(r["hours_left"] or 0.0 for r in st["arms"] if r["type"] not in GPU_TYPES)
    cpu_unknown = sum(r["total"] - r["done"] for r in st["arms"]
                      if r["type"] not in GPU_TYPES and r["hours_left"] is None and r["total"] > r["done"])
    log(f"[study] {src}: {st['flights_done']}/{st['flights_total']} flights done, "
        f"{st['trainings_done']}/{st['trainings_total']} trainings done")
    return {"source": src, "gpu_rows": rows, "models_left": models, "pilots_left": pilots,
            "open_loop_left": open_loop, "cpu_hours_left_status": round(cpu_h, 2),
            "cpu_flights_without_rate": int(cpu_unknown), "flights_done": st["flights_done"],
            "flights_total": st["flights_total"]}


def _rate_key(r: dict) -> str:
    arm, suite = r["arm"], r["suite"]
    if r["type"] == "mppi_physics":
        return "mppi_filtered@ORB1" if suite == "ORB1" else "mppi_filtered@A256"
    if arm == "v2single":
        return "v2single@ORB1" if suite == "ORB1" else "v2single@A256"
    return f"{arm}@A256"


def project(cfg: S5.StudyConfigV5, m: dict, left: dict) -> dict:
    """Hours of GPU work left, from the measurements m and the study's remaining work."""
    H = 3600.0
    full = _full_steps(cfg)
    pilot_runs = len(cfg.lr_mults) * float(cfg.lr_pilot_frac)
    train = {}
    for arm, r in m["training"].items():
        run_h = r["s_per_step"] * full / H
        n = len(left["models_left"][arm])
        p = pilot_runs if left["pilots_left"][arm] else 0.0
        train[arm] = {"h_per_model": round(run_h, 3), "models_left": n, "pilot_runs": p,
                      "h": round(run_h * (n + p), 2)}
    flights: Dict[str, dict] = {}
    for r in left["gpu_rows"]:
        k = _rate_key(r)
        rate = m["mppi"].get(k) or m["mppi"].get(k.split("@")[0] + "@A256")
        n = int(r["total"] - r["done"])
        d = flights.setdefault(k, {"flights_left": 0, "s_per_flight": round(rate["s_per_flight"], 2), "h": 0.0})
        d["flights_left"] += n
        d["h"] = round(d["flights_left"] * d["s_per_flight"] / H, 2)
    ol = {a: round(m["open_loop"][a]["s"] * len(left["open_loop_left"][a]) / H, 3) for a in m["open_loop"]}
    data_h = m.get("data_load_s", 0.0) / H
    h = {"e3_training": round(sum(v["h"] for v in train.values()), 2),
         "mppi_flights": round(sum(v["h"] for v in flights.values()), 2),
         "open_loop": round(sum(ol.values()), 2)}
    total = sum(h.values())
    u = units_per_hour(m["hardware"])
    sessions = max(1, int(-(-total // 6)))          # GPU sessions of up to ~6 h (each loads the data once)
    total += sessions * data_h
    cpu_u = left["cpu_hours_left_status"] * S5.UNITS_PER_HOUR["cpu"]
    return {"training": train, "flights": flights, "open_loop_h": ol, "hours": h, "data_load_h_per_session":
            round(data_h, 3), "sessions_est": sessions, "total_h": round(total, 1), "units_per_hour": u,
            "units": round(total * u, 1), "cpu_hours_left_status": left["cpu_hours_left_status"],
            "cpu_units_left_status": round(cpu_u, 1), "cpu_flights_without_rate": left["cpu_flights_without_rate"]}


def summary_text(m: dict, left: dict, p: dict) -> str:
    lines = [f"GPU {m['hardware']}, {m['workers']} GPU worker(s), {m['cpus']} CPU cores; study: {left['source']}",
             f"study flights done: {left['flights_done']:,}/{left['flights_total']:,}", "",
             "E3 training (per full model; models left + learning-rate pilots):"]
    for arm, r in p["training"].items():
        lines.append(f"  {arm}: {r['h_per_model'] * 60:.0f} min per model x ({r['models_left']} models + "
                     f"{r['pilot_runs']:.1f} pilot runs) = {r['h']:.1f} h")
    lines.append("MPPI flights left (s per flight on this GPU, all workers):")
    for k, r in p["flights"].items():
        lines.append(f"  {k}: {r['flights_left']:,} flights x {r['s_per_flight']:.1f} s = {r['h']:.1f} h")
    lines += [f"open loop: {p['hours']['open_loop']:.2f} h; data loading {p['data_load_h_per_session'] * 60:.0f} min "
              f"per GPU session x ~{p['sessions_est']} sessions", "",
              f"ALL REMAINING GPU WORK EXCEPT RL: {p['total_h']:.0f} hours = about {p['units']:.0f} Colab units "
              f"(at {p['units_per_hour']}/h)",
              f"CPU-pool flights left (the study's own measured rates): {p['cpu_hours_left_status']:.0f} h = about "
              f"{p['cpu_units_left_status']:.0f} units on a CPU runtime"
              + (f" (+ {p['cpu_flights_without_rate']:,} flights of arms with no rate yet)"
                 if p["cpu_flights_without_rate"] else "")]
    return "\n".join(lines)
