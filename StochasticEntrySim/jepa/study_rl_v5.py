"""
jepa/study_rl_v5.py

Notebook 2's study (jepa/plan_rl_v5.md): RQ3 with RL baselines that never decide the bank-reversal sign. The arms,
environments, gates and flight laws are jepa/rl2_v5.py's; this module is the study around them:

    config, layout   RLStudyConfig (frozen to 00_study/study_config.json), the output folder MyDrive/ReEntryAI_rl_v5
                     with notebook 1's folder (and v3/v4) read-only and guarded; FNPEG's tuned configuration copied
                     with its SHA-256 into 00_study/inputs.json (a later change stops Setup)
    stop rules       STOPPED.json / ACK_STOP / 00_study/stops/ (plan section 6): stop on a technical failure, never
                     on a bad result
    gates            RI (resid_fnpeg at init is FNPEG, fast env and DEV48), RP pilots on 512 fixed episodes, the
                     per-runtime rates and the routing to the cheaper runtime
    tuning, training DEV48 random search (<= 12 configurations), seeds 0-4, fast-env evaluation (R-gap)
    flights          A256 (every seed), DS1, AS1, ORB-1 block 1 (seed = case mod 5), CPU pool, flight-error stop
    status, analysis what is left and where; RQ3 with Holm over 3, estimates vs fnpeg, DS1/AS1/ORB-1, R-gap, gates,
                     dropped arms, tuning, cost (tables in 05_analysis/tables/, summary 05_analysis/RQ3_summary.md)

Uses study_v5 and analyze_v5 helpers as they are; its one hook there is analyze_v5.load_flights_v5(extra_v5=...).
"""

from __future__ import annotations

import hashlib
import json
import math
import multiprocessing as mp
import os
import time
import traceback
import zlib
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch

from jepa import fly_v5 as F5
from jepa import rl2_v5 as R2
from jepa import rl_v5 as RL
from jepa import runguard_v5 as RG
from jepa import study_v5 as S5
from jepa.models_v5 import fail_point

NOTEBOOK_REV_RL = 2           # ReEntryAI_RL_v5.ipynb's revision; an older copy stops at Setup
ARMS = R2.ARMS
SUITES = ("A256", "DS1", "AS1", "ORB1")

# plan section 2: the arms notebook 2 replaces, with their gate numbers (laptop pilots, commit 5471ceb)
DROPPED = {
    "dagger_fnpeg": {"RV1": "fail: held-out sign accuracy 0.919 vs best trivial 0.909, |bank| error 9.7 deg",
                     "RV2": "fail: success@5km 0.014 vs FNPEG 0.805 on the same fast-env episodes",
                     "pilot": "pilot_dagger2 (Part B session, commit 5471ceb)"},
    "bcppo_fnpeg": {"RV1": "pass: 0.957 vs 0.818, 5.8 deg",
                    "RV2": "fail: success@5km 0.031 vs FNPEG 0.805 on the same fast-env episodes",
                    "pilot": "pilot_bc2 (Part B session, commit 5471ceb)"},
}

# plan section 4.2: tuning spaces (keys that are not RLCfg fields are the arm's settings, rl2_v5.ARM_DEFAULTS)
SPACES = {
    "resid_fnpeg": {"bound_deg": ("choice", [10.0, 20.0, 30.0]), "lr": ("loguniform", 3e-5, 3e-4),
                    "init_std": ("choice", [0.1, 0.3]), "critic_warmup_updates": ("choice", [5, 20])},
    "liu_track": {"lr": ("loguniform", 1e-4, 1e-3), "ds_frac": ("choice", [0.2, 0.3, 0.4]),
                  "ent_coef": ("choice", [0.0, 1e-3, 1e-2]), "clip": ("choice", [0.1, 0.2, 0.3])},
    "ppo_hist": dict(S5.RL_SPACES["ppo_hist"]),
}
KW_KEYS = {"resid_fnpeg": ("bound_deg", "init_std"), "liu_track": ("ds_frac",), "ppo_hist": ()}
ZERO_STOP = S5.TUNE_ZERO_STOP

STOP_TRIGGERS = {
    "rl_nonfinite": "an RL run's weights or return went non-finite, or its training check stopped it (no flight "
                    "within 5 km and no rise in the return by 50 % of its budget)",
    "ri_failed": "gate RI failed: resid_fnpeg at init is not FNPEG bit for bit",
    "flight_errors": "5 or more of an arm's first 20 flights raised or ended non-finite",
    "failures": "3 or more failure records in one session",
    "stage_exception": "a stage raised",
}
GATES_PASS_ENV = "JEPA_RL2_GATES_PASS"   # smoke only ("1"): a failed RP still lets its arm tune, train and fly
RUNTIME_KEYS = ("guidance_batch",)
# Settings added after a folder may have been frozen: a frozen study_config.json without the key takes the new value
# and records the amendment (00_study/study_config.json "amendments"); a different frozen value still refuses.
AMENDED_KEYS = {"rpx_arms": "2026-11-21: extended pilot RPX (plan_rl_v5.md section 3.1, v5_plan.md section 15)",
                "rpx_snapshots": "2026-11-21: extended pilot RPX (plan_rl_v5.md section 3.1)"}


# ============================================================================
# Configuration
# ============================================================================
def _rl_arms_default() -> dict:
    return {"resid_fnpeg": {"n_envs": 2048, "n_steps": 16, "budget_decisions": 4_000_000, "minibatch": 4096,
                            "lr": 1e-4, "critic_warmup_updates": 5},
            "liu_track": {},
            "ppo_hist": {}}


@dataclass
class RLStudyConfig:
    """Every number of plan_rl_v5.md. Frozen to 00_study/study_config.json on first use."""
    mode: str = "full"                                  # "full" | "smoke" (every stage, tiny, numbers meaningless)
    rl_seeds: Tuple[int, ...] = (0, 1, 2, 3, 4)
    arms: Tuple[str, ...] = ARMS
    rl: dict = field(default_factory=dict)              # RLCfg overrides for every arm
    rl_arms: dict = field(default_factory=_rl_arms_default)
    arm_kw: dict = field(default_factory=lambda: json.loads(json.dumps(R2.ARM_DEFAULTS)))
    # gates (plan section 3)
    pilot_decisions: int = 20_000_000
    pilot_hours: float = 1.0
    fixed_n: int = 512
    fixed_seed: int = 12345
    ri_n: int = 64
    ri_steps: int = 400
    timing_min: float = 5.0
    rpx_arms: Tuple[str, ...] = ("liu_track", "ppo_hist")   # extended pilot at the full budget after a failed RP (3.1)
    rpx_snapshots: int = 10                              # RPX policies scored on the fixed episodes along the run
    # tuning (section 4.2)
    tune_n: int = 12
    tune_frac: float = 0.1
    tune_seed: int = 20261119
    # budgets (section 4.1)
    wall_min_h: float = 3.0
    wall_mult: float = 1.5
    fast_eval_episodes: int = 2048
    # stop rules (section 6)
    flight_err_first: int = 20
    flight_err_max: int = 5
    failures_stop: int = 3
    # analysis
    n_boot: int = 10_000
    holm_m: int = 3
    # runtime
    guidance_batch: int = 16
    units_per_hour: dict = field(default_factory=lambda: {"L4": 1.98, "cpu": RG.UNITS_PER_HOUR_CPU_HIGH})

    @staticmethod
    def smoke() -> "RLStudyConfig":
        rl = dict(n_envs=16, hidden=(32, 32), n_steps=8, epochs=2, minibatch=64, budget_decisions=256,
                  budget_hours=0.5, save_every_s=30.0, curve_points=2, env_pool=64, critic_warmup_updates=1)
        return RLStudyConfig(
            mode="smoke", rl_seeds=(0,), rl=rl,
            rl_arms={"resid_fnpeg": {"n_envs": 8, "n_steps": 4, "budget_decisions": 96, "critic_warmup_updates": 1},
                     "liu_track": {}, "ppo_hist": {}},
            pilot_decisions=128, pilot_hours=0.2, rpx_snapshots=2, fixed_n=8, ri_n=4, ri_steps=30, timing_min=0.1, tune_n=2,
            tune_frac=0.5, wall_min_h=0.5, fast_eval_episodes=8, n_boot=200, guidance_batch=4)

    def v5(self) -> S5.StudyConfigV5:
        """The v5 configuration whose suites (A256, DEV48, DS1, AS1, ORB-1) and analysis settings notebook 2 uses."""
        base = S5.StudyConfigV5.smoke() if self.mode == "smoke" else S5.StudyConfigV5()
        return replace(base, rl_seeds=tuple(self.rl_seeds), n_boot=int(self.n_boot),
                       guidance_batch=int(self.guidance_batch))


# ============================================================================
# Layout, inputs, frozen configuration
# ============================================================================
def layout(root, jepa_root, v3_root=None, v4_root=None) -> Tuple[S5.LayoutV5, Optional[S5.LayoutV5]]:
    """(notebook 2's folder, notebook 1's folder read as a layout); every write is checked, notebook 1's folder and
    the v3/v4 folders are read-only."""
    ro = [p for p in (jepa_root, v3_root, v4_root) if p]
    L = S5.LayoutV5(root, read_only=ro).make()
    J = S5.LayoutV5(jepa_root) if jepa_root and os.path.isdir(jepa_root) else None
    return L, J


def fnpeg_input(L: S5.LayoutV5, J: Optional[S5.LayoutV5], log=print) -> dict:
    """FNPEG's tuned configuration from notebook 1 (03_tuning/fnpeg/chosen.json), copied with its SHA-256 at first
    use; a later change of notebook 1's file stops here (flights under two FNPEG configurations must not mix)."""
    ip = L.study / "inputs.json"
    have = json.loads(ip.read_text()) if ip.exists() else None
    src = (J.tuning / "fnpeg" / "chosen.json") if J is not None else None
    if src is None or not src.exists():
        if have is not None:
            log(f"[inputs] notebook 1's FNPEG chosen.json not found; using the copy from {have['copied']}")
            return dict(have["fnpeg_cfg"])
        raise RuntimeError(f"FNPEG's tuned configuration is missing ({src}): notebook 1 (ReEntryAI_JEPA_v5) tunes it; "
                           "notebook 2 starts after that.")
    raw = src.read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    fn = dict(json.loads(raw)["config"])
    if have is not None:
        if have["fnpeg_sha256"] != sha:
            raise RuntimeError(f"notebook 1's FNPEG configuration changed since notebook 2 copied it ({src}: sha256 "
                               f"{sha[:12]}, copied {have['fnpeg_sha256'][:12]} on {have['copied']}). Notebook 2's "
                               "policies and flights use the copied one; start a new notebook-2 folder to use the "
                               "new one.")
        return dict(have["fnpeg_cfg"])
    L.atomic_bytes(L.study / "fnpeg_chosen.json", raw)
    L.write_json(ip, {"fnpeg_chosen": str(src), "fnpeg_sha256": sha, "fnpeg_cfg": fn,
                      "copied": time.strftime("%Y-%m-%d %H:%M:%S"), "jepa_root": str(J.root),
                      "git_commit": S5._GIT.get()})
    log(f"[inputs] FNPEG's tuned configuration copied (sha256 {sha[:12]}): {fn}")
    return fn


def freeze(cfg: RLStudyConfig, L: S5.LayoutV5, inputs: Optional[dict] = None) -> dict:
    path = L.study / "study_config.json"
    want = json.loads(json.dumps(asdict(cfg), default=list))
    if path.exists():
        have = json.loads(path.read_text())
        new = {k: v for k, v in want.items() if k in AMENDED_KEYS and k not in have["config"]}
        diff = {k: (have["config"].get(k), v) for k, v in want.items()
                if k not in RUNTIME_KEYS and k not in new and have["config"].get(k) != v}
        if not diff and new:
            have["config"].update(new)
            have.setdefault("amendments", []).append({"keys": new, "why": {k: AMENDED_KEYS[k] for k in new},
                                                      "time": time.strftime("%Y-%m-%d %H:%M:%S"),
                                                      "git_commit": S5._GIT.get()})
            L.write_json(path, have)
        if diff:
            raise ValueError(f"settings differ from the frozen {path}: {diff}. This folder was started with git commit "
                             f"{str(have.get('git_commit', '?'))[:10]}; continue it with that commit or start a new "
                             "folder.")
        return have
    doc = {"config": want, "git_commit": S5._GIT.get(), "created": time.strftime("%Y-%m-%d %H:%M:%S"),
           "inputs": inputs or {}, "notebook_rev": NOTEBOOK_REV_RL}
    L.write_json(path, doc)
    return doc


def notebook_check(rev) -> None:
    if rev != NOTEBOOK_REV_RL:
        raise RuntimeError(f"this notebook is revision {rev}, the code is revision {NOTEBOOK_REV_RL}: it is an old "
                           "copy. Open ReEntryAI_RL_v5.ipynb from GitHub (its Open in Colab badge), then Run all.")


# ============================================================================
# Stop rules (plan section 6)
# ============================================================================
class Stopped(RuntimeError):
    """A stop rule fired: 00_study/STOPPED.json is written and the runtime is released (the notebook's hook)."""


def stop(L: S5.LayoutV5, trigger: str, evidence: dict, log=print) -> None:
    doc = {"trigger": trigger, "description": STOP_TRIGGERS.get(trigger, trigger),
           "evidence": RG.json_safe(evidence), "time": time.strftime("%Y-%m-%d %H:%M:%S"),
           "git_commit": S5._GIT.get()}
    try:
        L.write_json(L.study / "STOPPED.json", doc)
    except Exception as e:                                   # the stop still happens
        log(f"[stop] STOPPED.json not written: {e}")
    RG.SUMMARY.failure("STOPPED", trigger, f"{doc['description']}: {json.dumps(doc['evidence'], default=str)[:400]}",
                       "00_study/STOPPED.json")
    RG.SUMMARY.write()
    msg = (f"STOPPED ({trigger}): {doc['description']}. Evidence in 00_study/STOPPED.json. Read it, fix or accept the "
           f"cause, then set ACK_STOP = \"{trigger}\" in the Settings cell and Run all.")
    log(f"[stop] {msg}")
    raise Stopped(msg)


def check_stopped(L: S5.LayoutV5, ack: Optional[str], log=print) -> Optional[dict]:
    """At Setup: a STOPPED.json refuses the session unless ACK_STOP names its trigger; acknowledged, the record moves
    to 00_study/stops/ and the notebook goes on."""
    p = L.study / "STOPPED.json"
    if not p.exists():
        return None
    doc = json.loads(p.read_text())
    if ack and ack == doc["trigger"]:
        dest = L.study / "stops" / f"{time.strftime('%Y%m%d_%H%M%S')}_{doc['trigger']}.json"
        L.write_json(dest, {**doc, "acknowledged": time.strftime("%Y-%m-%d %H:%M:%S")})
        os.remove(L.guard(p))
        log(f"[stop] {doc['trigger']} acknowledged (ACK_STOP); record moved to {dest.relative_to(L.root)}")
        return doc
    raise Stopped(f"this study STOPPED on {doc['time']} ({doc['trigger']}: {doc['description']}). Evidence: "
                  f"{json.dumps(doc['evidence'], default=str)[:800]}. Set ACK_STOP = \"{doc['trigger']}\" in the "
                  "Settings cell to continue.")


def check_failures(L: S5.LayoutV5, cfg: RLStudyConfig, log=print) -> None:
    fails = [f for f in RG.SUMMARY.failures if f.get("kind") != "STOPPED"]
    if len(fails) >= cfg.failures_stop:
        stop(L, "failures", {"n": len(fails), "failures": fails[-10:]}, log)


@contextmanager
def stage(L: S5.LayoutV5, cfg: RLStudyConfig, name: str, log=print):
    """A notebook stage: any exception is the stop rule stage_exception (a stop rule's own Stopped passes through)."""
    RG.STATE["stop_on_failure"] = True                        # notebook 2 stops on a technical failure
    try:
        with RG.stage(name, log) as rec:
            yield rec
    except (Stopped, KeyboardInterrupt):
        raise
    except Exception as e:
        stop(L, "stage_exception", {"stage": name, "error": f"{type(e).__name__}: {str(e)[:600]}",
                                    "traceback": traceback.format_exc()[-3000:]}, log)
    check_failures(L, cfg, log)


# ============================================================================
# Arm settings
# ============================================================================
def tmethod(arm: str) -> str:
    return f"rl2_{arm}"


def chosen(L: S5.LayoutV5, arm: str) -> Optional[dict]:
    p = L.tuning / tmethod(arm) / "chosen.json"
    return json.loads(p.read_text()) if p.exists() else None


def tuning_stopped(L: S5.LayoutV5, arm: str) -> Optional[dict]:
    p = L.tuning / tmethod(arm) / "stopped.json"
    return json.loads(p.read_text()) if p.exists() else None


def arm_settings(cfg: RLStudyConfig, L: Optional[S5.LayoutV5], arm: str, extra: Optional[dict] = None,
                 tuned: bool = True) -> Tuple[RL.RLCfg, dict]:
    """(RLCfg, the arm's own settings): common overrides, the arm's, the tuned configuration (tuned=True and
    chosen), then extra. A key of KW_KEYS[arm] goes to the arm's settings, every other to RLCfg."""
    rc = {**cfg.rl, **cfg.rl_arms.get(arm, {})}
    kw = dict(cfg.arm_kw.get(arm, {}))
    t = (chosen(L, arm) or {}).get("config", {}) if (tuned and L is not None) else {}
    for src in (t, extra or {}):
        for k, v in src.items():
            (kw if k in KW_KEYS.get(arm, ()) else rc)[k] = v
    if "hidden" in rc:
        rc["hidden"] = tuple(rc["hidden"])
    return RL.RLCfg(**rc), kw


def _rate_path(L, arm: str, runtime: str) -> Path:
    return L.checks / "rl2_rates" / f"{arm}_{runtime}.json"


def rates(L: S5.LayoutV5, arm: str) -> Dict[str, dict]:
    out = {}
    for rt in ("cpu", "L4"):
        p = _rate_path(L, arm, rt)
        if p.exists():
            out[rt] = json.loads(p.read_text())
    return out


def route(L: S5.LayoutV5, arm: str) -> Optional[dict]:
    """01_checks/route_<arm>.json (frozen once both the CPU and the L4 rate exist), else None (waits)."""
    p = L.checks / f"route_{arm}.json"
    if p.exists():
        return json.loads(p.read_text())
    r = R2.route(rates(L, arm))
    if r is None:
        return None
    r = {"arm": arm, **r, "time": time.strftime("%Y-%m-%d %H:%M:%S")}
    L.write_json(p, r)
    return r


def wall_cap_h(cfg: RLStudyConfig, L: S5.LayoutV5, arm: str, decisions: int) -> float:
    """The budget rule (plan section 4.1, change): decisions bind; the wall-clock cap is max(3 h, 1.5 x the time
    projected at the routed runtime's measured rate)."""
    r = route(L, arm)
    rt = (r or {}).get("runtime")
    rate = rates(L, arm).get(rt, {}).get("decisions_per_s") if rt else None
    if not rate:
        rr = [v.get("decisions_per_s") for v in rates(L, arm).values() if v.get("decisions_per_s")]
        rate = min(rr) if rr else None
    proj = float(decisions) / float(rate) / 3600.0 if rate else 0.0
    return float(max(cfg.wall_min_h, cfg.wall_mult * proj))


# ============================================================================
# Gates (plan section 3)
# ============================================================================
def _pilots(L) -> Path:
    return L.checks / "rl2_pilots"


def fnpeg_fixed(cfg: RLStudyConfig, L: S5.LayoutV5, fn: dict, workers: int, log=print) -> dict:
    """FNPEG alone on the 512 fixed fast-env episodes: every RP gate's reference (computed once)."""
    p = _pilots(L) / "fnpeg_fixed.json"
    if p.exists():
        return json.loads(p.read_text())
    rc, kw = arm_settings(cfg, None, "resid_fnpeg", {"label_workers": int(workers)}, tuned=False)
    t0 = time.time()
    with RG.allow("FNPEG on the fixed episodes", 90):
        r = R2.eval_fixed("fnpeg", None, rc, fn, kw, cfg.fixed_n, cfg.fixed_seed)
    r["wall_s"] = round(time.time() - t0, 1)
    L.write_json(p, r)
    log(f"[gates] FNPEG on the {cfg.fixed_n} fixed episodes: success@5km {r['success_5km_fast']:.3f} "
        f"({r['wall_s']:.0f} s)")
    return r


def _ri_worker(arg) -> List[dict]:
    conds, fn, rc, kw = arg
    torch.set_num_threads(1)
    rc = RL.RLCfg(**{**rc, "hidden": tuple(rc["hidden"])})
    out = []
    for c in conds:
        j = S5.guidance_job(c, "gate_ri", "gates", "fnpeg", fn)
        ra, rb = [], []
        F5.run_guidance_lockstep([j], F5.make_law("fnpeg", fn), on_done=ra.append, verbose=False)
        F5.run_guidance_lockstep([j], R2.ZeroResidLaw(fn, rc, kw), on_done=rb.append, verbose=False)
        same, diff = R2.flights_identical(ra[0], rb[0])
        out.append({"case": int(c["case"]), "atmosphere": c["atmosphere"], "identical": bool(same), "differs": diff[:10],
                    "miss_km_fnpeg": float(ra[0]["result"]["miss_km"]),
                    "miss_km_resid0": float(rb[0]["result"]["miss_km"])})
    return out


def gate_ri(cfg: RLStudyConfig, L: S5.LayoutV5, fn: dict, workers: int, log=print) -> dict:
    """Gate RI: (1) the initialized residual actor's mean is exactly 0 and its command FNPEG's, bit for bit, along
    whole fast-env episodes; (2) DEV48 flown by fnpeg and by resid_fnpeg@init in this session: identical flights.
    A failure is the stop rule ri_failed (01_checks/RI_resid_fnpeg.json; deleting it reruns the gate)."""
    p = L.checks / "RI_resid_fnpeg.json"
    if p.exists():
        return json.loads(p.read_text())
    rc, kw = arm_settings(cfg, None, "resid_fnpeg", {"label_workers": int(workers)}, tuned=False)
    t0 = time.time()
    fast = R2.ri_fast(rc, fn, kw, cfg.ri_n, cfg.ri_steps)
    log(f"[gates] RI part 1 (fast env): {fast}")
    conds = S5.conditions(cfg.v5(), "DEV48")
    w = max(1, min(int(workers), len(conds)))
    shards = [list(a) for a in np.array_split(np.arange(len(conds)), w) if len(a)]
    args = [([conds[i] for i in s], fn, asdict(rc), kw) for s in shards]
    with RG.allow("gate RI's DEV48 flights", 120):
        if w > 1:
            with mp.get_context("spawn").Pool(w) as pool:
                parts = pool.map(_ri_worker, args)
        else:
            parts = [_ri_worker(a) for a in args]
    rows = [r for part in parts for r in part]
    ok = bool(fast["passed"] and rows and all(r["identical"] for r in rows))
    doc = {"gate": "RI", "arm": "resid_fnpeg", "passed": ok, "fast": fast, "dev48": rows,
           "n_dev48_conditions": len(rows), "n_flights": 2 * len(rows),
           "differing": [r for r in rows if not r["identical"]], "arm_kw": kw, "wall_s": round(time.time() - t0, 1),
           "time": time.strftime("%Y-%m-%d %H:%M:%S"), "git_commit": S5._GIT.get()}
    L.write_json(p, doc)
    log(f"[gates] RI: {'PASS' if ok else 'FAIL'} ({len(rows)} DEV48 conditions, "
        f"{sum(not r['identical'] for r in rows)} differ)")
    if not ok:
        stop(L, "ri_failed", {"fast": fast, "differing": doc["differing"][:10]}, log)
    return doc


def _curve_rows(d: Path) -> List[dict]:
    st = torch.load(Path(d) / "train_state.pt", map_location="cpu", weights_only=False)
    return list(RL.Curve(st["curve"]).rows), st


def _pilot_run(cfg: RLStudyConfig, L: S5.LayoutV5, arm: str, fn: dict, device, workers: int, d: Path, extra: dict,
               label: str, log=print):
    """One untuned seed-0 run of an arm into d (resumable): (curve rows, train state, final actor, decisions, seconds,
    the training check's reason or None, RLCfg, arm settings)."""
    rc, kw = arm_settings(cfg, None, arm, {"seed": 0, "label_workers": int(workers), **extra}, tuned=False)
    why = None
    try:
        info = R2.train_dir(L, d, arm, rc, device, fn, kw, log, label=label)
    except RL.RLStopped as e:
        why, info = str(e), {}
    rows, st = _curve_rows(d)
    dec, secs = int(info.get("decisions", st["step"])), float(info.get("train_s", st["elapsed"]))
    if (d / "policy.pt").exists():
        actor = R2.load_actor(d / "policy.pt", device)
    else:
        actor = R2.make_actor(arm, rc, kw).to(device)
        actor.load_state_dict(st["actor"])
        actor.eval()
    return rows, st, actor, dec, secs, why, rc, kw


def _finite(actor) -> bool:
    return all(bool(torch.isfinite(q).all()) for q in actor.parameters())


def gate_rp(cfg: RLStudyConfig, L: S5.LayoutV5, arm: str, fn: dict, device, workers: int, runtime: str,
            log=print) -> Optional[dict]:
    """Gate RP: the arm's untuned pilot (seed 0, 2e7 decisions or 1 h) against FNPEG on the fixed episodes. A failed
    RP is a result (the arm gets no tuning, training or flights, unless its extended pilot RPX passes); non-finite
    weights or returns are rl_nonfinite. Also the pilot's decisions per second on this runtime
    (01_checks/rl2_rates/<arm>_<runtime>.json)."""
    p = L.checks / f"RP_{arm}.json"
    if p.exists():
        return json.loads(p.read_text())
    if arm == "resid_fnpeg":
        ri = L.checks / "RI_resid_fnpeg.json"
        if not ri.exists() or not json.loads(ri.read_text())["passed"]:
            log("[gates] resid_fnpeg's pilot waits for gate RI")
            return None
    ref = fnpeg_fixed(cfg, L, fn, workers, log)
    rows, st, actor, dec, secs, why, rc, kw = _pilot_run(
        cfg, L, arm, fn, device, workers, _pilots(L) / arm,
        {"budget_decisions": int(cfg.pilot_decisions), "budget_hours": float(cfg.pilot_hours), "curve_points": 0},
        f"pilot {arm}", log)
    finite = _finite(actor)
    rate = {"arm": arm, "runtime": runtime, "decisions_per_s": dec / max(secs, 1e-9), "decisions": dec,
            "train_s": secs, "units_per_hour": float(cfg.units_per_hour[runtime]), "source": "pilot",
            "n_envs": int(rc.n_envs), "label_workers": int(workers), "device": str(device),
            "time": time.strftime("%Y-%m-%d %H:%M:%S")}
    L.write_json(_rate_path(L, arm, runtime), rate)
    if not finite or (why and "non-finite" in why):
        doc = {"gate": "RP", "arm": arm, "passed": False, "finite": False, "reason": why or "non-finite weights",
               "rate": rate, "time": time.strftime("%Y-%m-%d %H:%M:%S"), "git_commit": S5._GIT.get()}
        L.write_json(p, doc)
        stop(L, "rl_nonfinite", {"arm": arm, "stage": "pilot", "reason": doc["reason"]}, log)
    ev = R2.eval_fixed(arm, actor, rc, fn, kw, cfg.fixed_n, cfg.fixed_seed, device)
    v = R2.rp_verdict(arm, rows, actor, ev["success_5km_fast"], float(ref["success_5km_fast"]))
    if why:
        v["passed"] = False
        v["reason"] += f"; the pilot's training check stopped it: {why}"
    doc = {"gate": "RP", "arm": arm, **v, "pilot": {"decisions": dec, "train_s": round(secs, 1), "cfg": asdict(rc),
                                                     "arm_kw": kw, "stopped": why},
           "fixed_eval": {k: ev[k] for k in ("n", "seed", "success_5km_fast", "miss_km_fast_median")},
           "rate": rate, "time": time.strftime("%Y-%m-%d %H:%M:%S"), "git_commit": S5._GIT.get()}
    L.write_json(p, doc)
    log(f"[gates] RP {arm}: {'PASS' if v['passed'] else 'FAIL (a result: no tuning, training or flights)'}: "
        f"{v['reason']}; {rate['decisions_per_s']:.0f} decisions/s on {runtime}")
    return doc


def _rp_doc(L, arm) -> Optional[dict]:
    p = L.checks / f"RP_{arm}.json"
    return json.loads(p.read_text()) if p.exists() else None


def _rpx_doc(L, arm) -> Optional[dict]:
    p = L.checks / f"RPX_{arm}.json"
    return json.loads(p.read_text()) if p.exists() else None


def _smoke_pass(cfg: RLStudyConfig) -> bool:
    return cfg.mode == "smoke" and os.environ.get(GATES_PASS_ENV) == "1"


def rpx_due(cfg: RLStudyConfig, L: S5.LayoutV5, arm: str) -> bool:
    """Does the arm get the extended pilot (plan section 3.1)? Its RP failed with finite weights and it is in
    cfg.rpx_arms."""
    g = _rp_doc(L, arm)
    return bool(arm in cfg.rpx_arms and g is not None and not g["passed"] and g.get("finite", True))


def rpx_pending(cfg: RLStudyConfig, L: S5.LayoutV5, arm: str) -> bool:
    return rpx_due(cfg, L, arm) and _rpx_doc(L, arm) is None


def exploratory(cfg: RLStudyConfig, L: S5.LayoutV5, arm: str) -> bool:
    """An arm that failed RP and passed RPX: it trains and flies, and its comparisons are exploratory only."""
    x = _rpx_doc(L, arm)
    return bool(rpx_due(cfg, L, arm) and x is not None and x["passed"])


def gate_rpx(cfg: RLStudyConfig, L: S5.LayoutV5, arm: str, fn: dict, device, workers: int, runtime: str,
             log=print) -> Optional[dict]:
    """Extended pilot RPX (plan section 3.1): an arm that failed RP runs once more untuned, seed 0, at its full budget
    (the budget rule's wall-clock cap, the same training check as a full run), on its routed runtime; the final
    policy is judged by the arm's RP rule on the same fixed episodes against the same FNPEG reference, and
    cfg.rpx_snapshots snapshots are scored there too (the learning curve). Pass: the arm goes on as exploratory."""
    p = L.checks / f"RPX_{arm}.json"
    if p.exists():
        return json.loads(p.read_text())
    if not rpx_due(cfg, L, arm):
        return None
    ref = fnpeg_fixed(cfg, L, fn, workers, log)
    full, _ = arm_settings(cfg, L, arm, tuned=False)
    dec_budget = int(full.budget_decisions)
    hours = wall_cap_h(cfg, L, arm, dec_budget)
    d = _pilots(L) / f"{arm}_full"
    lim = S5.SESSION["limit"]
    if lim is not None and not (d / "policy.pt").exists() and not lim.may_start(
            hours, log, f"extended pilot {arm}", extra_s=600.0):
        return None
    log(f"[gates] RPX {arm}: the full budget, {dec_budget:,} decisions (wall-clock cap {hours:.2f} h) on {runtime}")
    rows, st, actor, dec, secs, why, rc, kw = _pilot_run(
        cfg, L, arm, fn, device, workers, d,
        {"budget_decisions": dec_budget, "budget_hours": hours, "curve_points": int(cfg.rpx_snapshots)},
        f"extended pilot {arm}", log)
    if lim is not None:
        lim.done()
    if not _finite(actor) or (why and "non-finite" in why):
        doc = {"gate": "RPX", "arm": arm, "passed": False, "finite": False, "reason": why or "non-finite weights",
               "time": time.strftime("%Y-%m-%d %H:%M:%S"), "git_commit": S5._GIT.get()}
        L.write_json(p, doc)
        stop(L, "rl_nonfinite", {"arm": arm, "stage": "extended pilot", "reason": doc["reason"]}, log)
    curve = []
    snaps = sorted((d / "snapshots").glob("policy_*of*.pt"), key=lambda q: int(q.stem.split("_")[1].split("of")[0]))
    for sp in snaps:
        meta = R2.policy_meta(sp)
        e = R2.eval_fixed(arm, R2.load_actor(sp, device), rc, fn, kw, cfg.fixed_n, cfg.fixed_seed, device)
        curve.append({"point": sp.stem, "decisions": int(meta.get("decisions", 0)), "trigger": meta.get("trigger"),
                      "success_5km_fixed": float(e["success_5km_fast"]),
                      "miss_km_fixed_median": float(e["miss_km_fast_median"])})
    ev = R2.eval_fixed(arm, actor, rc, fn, kw, cfg.fixed_n, cfg.fixed_seed, device)
    v = R2.rp_verdict(arm, rows, actor, ev["success_5km_fast"], float(ref["success_5km_fast"]))
    v["reason"] = "full-budget run" + v["reason"][len("pilot"):] if v["reason"].startswith("pilot") else v["reason"]
    if why:
        v["passed"] = False
        v["reason"] += f"; the run's training check stopped it: {why}"
    reached = "the full budget" if dec >= dec_budget else f"{dec:,} of {dec_budget:,} decisions"
    doc = {"gate": "RPX", "arm": arm, **v, "reached": reached,
           "run": {"decisions": dec, "budget_decisions": dec_budget, "wall_cap_h": hours, "train_s": round(secs, 1),
                   "runtime": runtime, "cfg": asdict(rc), "arm_kw": kw, "stopped": why},
           "fixed_eval": {k: ev[k] for k in ("n", "seed", "success_5km_fast", "miss_km_fast_median")},
           "fixed_curve": curve, "rp": {k: (_rp_doc(L, arm) or {}).get(k) for k in ("success_5km_fixed", "reason")},
           "time": time.strftime("%Y-%m-%d %H:%M:%S"), "git_commit": S5._GIT.get()}
    L.write_json(p, doc)
    log(f"[gates] RPX {arm}: {'PASS (the arm goes on as exploratory)' if v['passed'] else 'FAIL (a result)'} "
        f"after {reached}: {v['reason']}; fixed-episode curve: "
        + ", ".join(f"{c['decisions']:,}: {c['success_5km_fixed']:.3f}" for c in curve))
    return doc


def gate_passed(cfg: RLStudyConfig, L: S5.LayoutV5, arm: str) -> Optional[bool]:
    """Does the arm tune, train and fly? RP's verdict; for an arm that failed RP and gets the extended pilot, RPX's
    (None: not decided yet). Smoke mode with JEPA_RL2_GATES_PASS=1 lets a failed arm go on (after its RPX, if due)."""
    g = _rp_doc(L, arm)
    if g is None:
        return None
    if g["passed"]:
        return True
    if rpx_due(cfg, L, arm):
        x = _rpx_doc(L, arm)
        if x is None:
            return None
        ok = bool(x["passed"])
    else:
        ok = False
    return True if (not ok and _smoke_pass(cfg)) else ok


def timing_run(cfg: RLStudyConfig, L: S5.LayoutV5, arm: str, fn: dict, device, workers: int, runtime: str = "L4",
               log=print) -> Optional[dict]:
    """A cfg.timing_min-minute run of a gated arm on this runtime: its decisions per second for the routing."""
    p = _rate_path(L, arm, runtime)
    if p.exists():
        return json.loads(p.read_text())
    if not (gate_passed(cfg, L, arm) or rpx_pending(cfg, L, arm)):
        return None
    rc, kw = arm_settings(cfg, None, arm, {"seed": 0, "label_workers": int(workers), "budget_decisions": 10 ** 12,
                                           "budget_hours": cfg.timing_min / 60.0, "curve_points": 0,
                                           "save_every_s": 1e9}, tuned=False)
    d = L.checks / "rl2_rates" / f"{arm}_{runtime}_run"
    try:
        info = R2.train_dir(L, d, arm, rc, device, fn, kw, log, label=f"timing {arm} on {runtime}")
    except RL.RLStopped:
        info = {}
    _, st = _curve_rows(d)
    dec, secs = int(info.get("decisions", st["step"])), float(info.get("train_s", st["elapsed"]))
    rate = {"arm": arm, "runtime": runtime, "decisions_per_s": dec / max(secs, 1e-9), "decisions": dec,
            "train_s": secs, "units_per_hour": float(cfg.units_per_hour[runtime]), "source": "timing",
            "n_envs": int(rc.n_envs), "label_workers": int(workers), "device": str(device),
            "time": time.strftime("%Y-%m-%d %H:%M:%S")}
    L.write_json(p, rate)
    log(f"[routing] {arm} on {runtime}: {rate['decisions_per_s']:.0f} decisions/s "
        f"({R2.units_per_million(rate['decisions_per_s'], rate['units_per_hour']):.3f} units per million)")
    return rate


# ============================================================================
# Flights
# ============================================================================
def seed_of(cfg: RLStudyConfig, cond: dict) -> int:
    """DS1, AS1, ORB-1: case i is flown by training seed i mod 5 (v5's ORB-1 rule)."""
    s = tuple(cfg.rl_seeds)
    return int(s[int(cond.get("train_seed", cond["case"])) % len(s)])


def rl_job(cond: dict, arm: str, stage: str, path, seed: int, kw: dict, fn: dict, kind: Optional[str] = None) -> dict:
    """A flight of controller `arm` (its name in the results, e.g. tune_resid_fnpeg_c00) flown as notebook-2 arm
    `kind` (default: arm), which selects the law (rl2_v5.make_law)."""
    j = S5.policy_job(cond, arm, stage, str(path), seed)
    j["model"]["rl2"] = {"arm": kind or arm, "fnpeg_cfg": dict(fn), **{k: v for k, v in kw.items()}}
    return j


def arm_ready(cfg: RLStudyConfig, L: S5.LayoutV5, arm: str) -> bool:
    return bool(gate_passed(cfg, L, arm)) and chosen(L, arm) is not None and tuning_stopped(L, arm) is None


def flight_jobs(cfg: RLStudyConfig, L: S5.LayoutV5, arm: str, suite: str, fn: dict) -> List[dict]:
    """The test flights of an arm (empty for an arm that failed its gate or stopped its tuning)."""
    if arm not in cfg.arms or not arm_ready(cfg, L, arm):
        return []
    _, kw = arm_settings(cfg, L, arm)
    conds = S5.conditions(cfg.v5(), suite)
    if suite == "ORB1":
        conds = [c for c in conds if int(c["block"]) == 0]            # ORB-1 block 1
    if suite == "A256":
        return [rl_job(c, arm, suite, L.policy_dir(arm, s) / "policy.pt", s, kw, fn) for s in cfg.rl_seeds
                for c in conds]
    return [rl_job(c, arm, suite, L.policy_dir(arm, seed_of(cfg, c)) / "policy.pt", seed_of(cfg, c), kw, fn)
            for c in conds]


def _fly(shard: Sequence[dict], L: S5.LayoutV5) -> None:
    S5._fault(shard)
    law = R2.make_law(shard[0]["model"])
    F5.run_guidance_lockstep(shard, law, on_done=S5._checked_writer(L), verbose=False)


def _flight_worker(arg) -> dict:
    """Pool task: one lockstep batch; a batch that raises is retried flight by flight (the laws are
    batch-independent); a flight that raises again is returned as an error, not written."""
    shard, root, ro = arg
    torch.set_num_threads(1)
    L = S5.LayoutV5(root, ro)
    try:
        _fly(shard, L)
        return {"n": len(shard), "errors": [], "arm": shard[0]["controller"]}
    except Exception:
        errors = []
        for j in S5._not_done(L, shard):
            try:
                _fly([j], L)
            except Exception as e:
                errors.append(S5._err(j, e))
        return {"n": len(shard) - len(errors), "errors": errors, "arm": shard[0]["controller"]}


_COUNT = {"flights": 0}


def run_flights(cfg: RLStudyConfig, L: S5.LayoutV5, jobs: Sequence[dict], workers: int, log=print,
                label: str = "flights") -> int:
    """Fly the jobs not done whose policy exists (CPU pool). >= flight_err_max of an arm's first flight_err_first
    flights raising (or ending non-finite) is the stop rule flight_errors."""
    done = S5.done_ids(L, jobs)
    todo = [j for j in jobs if F5.flight_id_v5(j) not in done and os.path.exists(j["model"]["path"])]
    wait = sum(1 for j in jobs if F5.flight_id_v5(j) not in done and not os.path.exists(j["model"]["path"]))
    log(f"[{label}] {len(done)}/{len(jobs)} done, {len(todo)} to fly on {workers} worker(s)"
        + (f", {wait} wait for training" if wait else ""))
    if not todo:
        return 0
    groups: Dict[Tuple[str, str], List[dict]] = {}
    for j in todo:
        groups.setdefault((j["controller"], j["model"]["path"]), []).append(j)
    n_b = max(1, int(cfg.guidance_batch))
    tasks = [g[i:i + n_b] for g in groups.values() for i in range(0, len(g), n_b)]
    lim = S5.SESSION["limit"]
    if lim is not None and not lim.may_start(0.0, log, f"{label} flights"):
        return 0
    ro = [str(p) for p in L.read_only]
    frac = (cfg.flight_err_max - 0.5) / cfg.flight_err_first          # > frac * n_first  <=>  >= flight_err_max
    errs: Dict[str, S5._FlightErrors] = {}
    pool = mp.get_context("spawn").Pool(int(workers)) if workers > 1 else None
    it = pool.imap(_flight_worker, [(t, str(L.root), ro) for t in tasks]) if pool else \
        map(_flight_worker, [(t, str(L.root), ro) for t in tasks])
    n, t0, bad = 0, time.time(), None
    step = max(1, len(tasks) // 20)
    try:
        for i, k in enumerate(it, 1):
            n += k["n"]
            e = errs.setdefault(k["arm"], S5._FlightErrors(L, f"{label} {k['arm']}", log, cfg.flight_err_first, frac))
            e.add(k)
            _COUNT["flights"] += 1
            fail_point("rl2_flights", _COUNT["flights"])
            if e.too_many():
                bad = k["arm"]
                break
            if i % step == 0 or i == len(tasks):
                el = time.time() - t0
                log(f"  [{label}] {n}/{len(todo)} flights in {el / 60:.1f} min "
                    f"(~{el / max(n, 1) * (len(todo) - n) / 60:.0f} min left)")
    finally:
        if pool is not None:
            if bad:
                pool.terminate()
            else:
                pool.close()
            pool.join()
    for arm, e in errs.items():
        try:
            e.close(finished=True)
        except RG.StageFailed as ex:
            bad = bad or arm
            stop(L, "flight_errors", {"arm": arm, "label": label, "reason": str(ex), "errors": e.errors[:20]}, log)
    if bad:
        stop(L, "flight_errors", {"arm": bad, "label": label, "errors": errs[bad].errors[:20]}, log)
    return n


# ============================================================================
# Tuning and training
# ============================================================================
def _zero_stop(L: S5.LayoutV5, arm: str, log=print) -> bool:
    """True (stopped.json written, a result, not a failure) when the first ZERO_STOP configurations all scored 0 %."""
    m = L.tuning / tmethod(arm)
    if (m / "stopped.json").exists():
        return True
    z = [m / f"config_{k:02d}.json" for k in range(ZERO_STOP)]
    if not all(p.exists() for p in z):
        return False
    if not all(float(json.loads(p.read_text())["score"]["success_5km"]) == 0.0 for p in z):
        return False
    L.write_json(m / "stopped.json", {"method": tmethod(arm), "status": "stopped",
                                      "reason": f"its first {ZERO_STOP} configurations all scored 0 % success@5km on "
                                                "DEV48", "time": time.strftime("%Y-%m-%d %H:%M:%S"),
                                      "git_commit": S5._GIT.get()})
    log(f"[tuning] {arm}: stopped (first {ZERO_STOP} configurations at 0 %); reported as a result")
    return True


def tune(cfg: RLStudyConfig, L: S5.LayoutV5, arm: str, fn: dict, device, workers: int, log=print) -> Optional[dict]:
    """DEV48 random search: config 0 = the defaults, cfg.tune_n configurations, seed 0, tune_frac of the budget, each
    scored by its 48 full-simulator flights; then S5's choice rule (success@5km, mean log miss, peak g <= 4.6 g)."""
    m = tmethod(arm)
    c0 = chosen(L, arm)
    if c0 is not None or tuning_stopped(L, arm) or not gate_passed(cfg, L, arm):
        return c0
    full, _ = arm_settings(cfg, L, arm, tuned=False)
    dec = max(1, int(full.budget_decisions * cfg.tune_frac))
    dev = S5.conditions(cfg.v5(), "DEV48")
    lim = S5.SESSION["limit"]
    configs = S5.sample_configs(SPACES[arm], cfg.tune_n, cfg.tune_seed + zlib.crc32(arm.encode()) % 1000)
    for j, c in enumerate(configs):
        if j >= ZERO_STOP and _zero_stop(L, arm, log):
            return None
        cp = L.tuning / m / f"config_{j:02d}.json"
        if cp.exists():
            continue
        d = L.tuning / m / f"config_{j:02d}"
        rc, kw = arm_settings(cfg, L, arm, {**c, "seed": 0, "label_workers": int(workers), "budget_decisions": dec,
                                            "budget_hours": wall_cap_h(cfg, L, arm, dec), "curve_points": 0},
                              tuned=False)
        if lim is not None and not (d / "policy.pt").exists() and not lim.may_start(
                rc.budget_hours, log, f"tune {arm} config {j}", extra_s=len(dev) * 90.0 / max(1, int(workers))):
            return None
        note = None
        try:
            R2.train_dir(L, d, arm, rc, device, fn, kw, log, label=f"tune {arm} config {j}")
        except RL.RLStopped as e:
            if "non-finite" in str(e):
                RG.record_failure(L, d / "failed.json", {"arm": arm, "config": j, "reason": str(e)}, "RL tuning",
                                  f"{arm} config {j}", log)
                stop(L, "rl_nonfinite", {"arm": arm, "stage": f"tuning config {j}", "reason": str(e)}, log)
            note = str(e)
        if note:                         # a configuration that does not learn scores 0 (a result, not a failure)
            score = {"n": 0, "success_5km": 0.0, "mean_log_miss": 99.0, "mean_peak_g": 99.0,
                     "median_miss_km": float("nan"), "stopped": note}
        else:
            js = [rl_job(cd, f"tune_{arm}_c{j:02d}", "tuning", d / "policy.pt", 0, kw, fn, kind=arm) for cd in dev]
            run_flights(cfg, L, js, workers, log, f"tuning {arm} config {j}")
            rs = S5._results(L, js)
            if len(rs) < len(js):
                return None
            score = S5._score_dev(rs)
        L.write_json(cp, {"method": m, "config_index": j, "config": dict(c), "sampled": dict(c),
                          "budget": {"budget_decisions": dec, "budget_hours": rc.budget_hours}, "score": score})
        log(f"[tuning] {arm} config {j}: {score}")
        if lim is not None:
            lim.done()
    if len(configs) >= ZERO_STOP and _zero_stop(L, arm, log):
        return None
    return S5._choose_method(cfg.v5(), L, m, len(configs), log)


def fast_eval(cfg: RLStudyConfig, L: S5.LayoutV5, arm: str, seed: int, fn: dict, device) -> Optional[dict]:
    """A trained policy in the fast env: the fixed-episode distribution and A256's cases (the fast side of R-gap)."""
    d = L.policy_dir(arm, seed)
    p = d / "fast_eval.json"
    if p.exists() or not (d / "policy.pt").exists():
        return json.loads(p.read_text()) if p.exists() else None
    rc, kw = arm_settings(cfg, L, arm, {"seed": int(seed)})
    actor = R2.load_actor(d / "policy.pt", device)
    conds = S5.conditions(cfg.v5(), "A256")
    tr = R2.eval_fixed(arm, actor, rc, fn, kw, cfg.fast_eval_episodes, cfg.fixed_seed, device)
    a2 = R2.eval_fixed(arm, actor, rc, fn, kw, device=device, conds=conds)
    doc = {"train_distribution": {k: tr[k] for k in ("n", "seed", "success_5km_fast", "miss_km_fast_median")},
           "A256": {"n": a2["n"], "success_5km_fast": a2["success_5km_fast"],
                    "miss_km_fast_median": a2["miss_km_fast_median"], "success_by_case": a2["success_by_episode"]},
           "A256_cases": [[int(c["case"]), c["atmosphere"]] for c in conds]}
    L.write_json(p, doc)
    return doc


def train(cfg: RLStudyConfig, L: S5.LayoutV5, arm: str, fn: dict, device, workers: int, log=print) -> List[dict]:
    """Every seed of a tuned arm at the full budget (resumable); then its fast evaluation. A run stopped by its
    training check gets failed.json and is the stop rule rl_nonfinite."""
    out = []
    if chosen(L, arm) is None:
        log(f"[train] {arm} waits for its tuning")
        return out
    lim = S5.SESSION["limit"]
    for s in cfg.rl_seeds:
        d = L.policy_dir(arm, s)
        if (d / "failed.json").exists():
            log(f"[train] {arm} seed {s}: stopped earlier ({json.loads((d / 'failed.json').read_text())['reason']}); "
                "delete failed.json to train it again")
            continue
        rc, kw = arm_settings(cfg, L, arm, {"seed": int(s), "label_workers": int(workers)})
        rc = replace(rc, budget_hours=wall_cap_h(cfg, L, arm, rc.budget_decisions))
        fresh = not (d / "fast_eval.json").exists()
        if lim is not None and fresh and not lim.may_start(0.0 if (d / "policy.pt").exists() else rc.budget_hours,
                                                           log, f"{arm} seed {s}", extra_s=600.0):
            return out
        try:
            info = R2.train_dir(L, d, arm, rc, device, fn, kw, log, label=f"{arm} seed {s}")
        except (RL.RLStopped, FloatingPointError) as e:
            RG.record_failure(L, d / "failed.json", {"arm": arm, "seed": int(s), "status": "stopped", "reason": str(e),
                                                     "git_commit": S5._GIT.get()}, "RL run", f"{arm} seed {s}", log)
            stop(L, "rl_nonfinite", {"arm": arm, "seed": int(s), "reason": str(e)}, log)
        fast_eval(cfg, L, arm, s, fn, device)
        out.append(info)
        if lim is not None and fresh:
            lim.done()
    return out


# ============================================================================
# Status and the next pass
# ============================================================================
def runtime_kind(device: str, gpu_name: str = "") -> str:
    if str(device).startswith("cuda"):
        return "L4" if "L4" in str(gpu_name).split() or str(gpu_name).endswith("L4") else str(gpu_name)
    return "cpu"


def status(cfg: RLStudyConfig, L: S5.LayoutV5, fn: Optional[dict]) -> dict:
    rows = []
    for arm in cfg.arms:
        g = _rp_doc(L, arm)
        x = _rpx_doc(L, arm)
        r = route(L, arm)
        n_cfg = sum(1 for k in range(cfg.tune_n) if (L.tuning / tmethod(arm) / f"config_{k:02d}.json").exists())
        trained = sum(1 for s in cfg.rl_seeds if (L.policy_dir(arm, s) / "policy.pt").exists())
        row = {"arm": arm, "RP": "not run" if g is None else ("pass" if g["passed"] else "FAIL (result)"),
               "RPX": ("pending" if x is None else ("pass (exploratory)" if x["passed"] else "FAIL (result)"))
               if rpx_due(cfg, L, arm) else "-",
               "rates": ",".join(sorted(rates(L, arm))), "route": (r or {}).get("runtime", "waits"),
               "tuning": "stopped" if tuning_stopped(L, arm) else ("chosen" if chosen(L, arm) else f"{n_cfg}/{cfg.tune_n}"),
               "seeds_trained": f"{trained}/{len(cfg.rl_seeds)}"}
        for su in SUITES:
            js = flight_jobs(cfg, L, arm, su, fn) if fn else []
            row[su] = f"{len(S5.done_ids(L, js))}/{len(js)}" if js else "-"
        rows.append(row)
    ri = L.checks / "RI_resid_fnpeg.json"
    doc = {"arms": rows, "RI": json.loads(ri.read_text())["passed"] if ri.exists() else None,
           "stopped": (L.study / "STOPPED.json").exists(), "next_pass": next_pass(cfg, L, fn),
           "time": time.strftime("%Y-%m-%d %H:%M:%S")}
    L.write_json(L.study / "status.json", doc)
    return doc


def next_pass(cfg: RLStudyConfig, L: S5.LayoutV5, fn: Optional[dict]) -> str:
    """What the next Run all has to do, and on which runtime."""
    cpu, l4 = [], []
    for arm in cfg.arms:
        if _rp_doc(L, arm) is None:
            cpu.append(f"gate RP {arm}")
            continue
        if rpx_pending(cfg, L, arm):
            if "L4" not in rates(L, arm):
                l4.append(f"timing {arm}")
            else:
                r = route(L, arm)
                (l4 if r and r["runtime"] == "L4" else cpu).append(f"extended pilot {arm} (full budget)")
            continue
        g = gate_passed(cfg, L, arm)
        if not g:
            continue
        rt = rates(L, arm)
        if "L4" not in rt:
            l4.append(f"timing {arm}")
            continue
        r = route(L, arm)
        where = l4 if r and r["runtime"] == "L4" else cpu
        if tuning_stopped(L, arm):
            continue
        if chosen(L, arm) is None:
            where.append(f"tuning {arm}")
        elif any(not (L.policy_dir(arm, s) / "policy.pt").exists() and not (L.policy_dir(arm, s) / "failed.json").exists()
                 for s in cfg.rl_seeds):
            where.append(f"training {arm}")
        elif fn:
            left = sum(len(js) - len(S5.done_ids(L, js)) for js in (flight_jobs(cfg, L, arm, su, fn) for su in SUITES))
            if left:
                cpu.append(f"{left} flights of {arm}")
    if l4:
        return "L4 pass (GPU_TYPE L4): " + "; ".join(l4) + ("; then a CPU pass: " + "; ".join(cpu) if cpu else "")
    if cpu:
        return "CPU pass (High-RAM CPU runtime): " + "; ".join(cpu)
    return "nothing left to train or fly: the analysis is complete once notebook 1's paired flights exist"


# ============================================================================
# Analysis (plan section 4.5)
# ============================================================================
def load_flights(cfg: RLStudyConfig, L: S5.LayoutV5, J: Optional[S5.LayoutV5], L3=None, L4=None):
    """Notebook 2's flights joined with notebook 1's (read only, from its flights.csv) and v3/v4's. Notebook 1's own
    copies of notebook 2's arm names (if any) are dropped: RQ3's arms are notebook 2's."""
    from jepa import analyze_v5 as A5
    df = A5.load_flights_v5(cfg.v5(), L, L3, L4, extra_v5=[J] if J is not None else [])
    if J is not None and "v5_folder" in df:
        df = df[~((df.v5_folder == str(J.root)) & df.arm.isin(list(ARMS)))].reset_index(drop=True)
    return df


def criteria_rq3(cfg: RLStudyConfig, L: S5.LayoutV5, t: "pd.DataFrame") -> "pd.DataFrame":
    """R-each for the three arms, Holm over m = 3 (an arm without a comparison enters with p = 1)."""
    import pandas as pd
    from jepa import analyze_v3 as A3
    out = []
    for arm in ARMS:
        r = None
        if len(t):
            q = t[(t.id == f"R-{arm}") & (t.metric == "success5") & (t.atmosphere == "all")]
            r = q.iloc[0].to_dict() if len(q) else None
        g = _rp_doc(L, arm)
        failed_rp = g is not None and not g["passed"] and not _smoke_pass(cfg)
        row = {"id": f"R-{arm}", "question": f"Does v2 beat {arm}? (v2single - {arm}, success@5km, A256)",
               "gate_RP": None if g is None else bool(g["passed"]), "diff_pp": np.nan, "ci_lo": np.nan,
               "ci_hi": np.nan, "p": 1.0, "seeds_holding": None, "n_pairs": 0, "n_expected": 0, "met": None}
        if failed_rp:                         # confirmatory: p = 1 whatever the extended pilot shows (section 3.1)
            row["status"] = ("did not reach the gate (RPX passed: exploratory comparison in RQ3_exploratory.csv)"
                             if exploratory(cfg, L, arm) else "did not reach the gate")
        elif r is None:
            row["status"] = "tuning stopped" if tuning_stopped(L, arm) else "not run"
        else:
            n = int(r["n_seeds"])
            need = 4 if n >= 5 else n
            done = bool(r["complete"])
            row.update({"diff_pp": r["diff"], "ci_lo": r["ci_lo"], "ci_hi": r["ci_hi"], "p": float(r["p_boot"]),
                        "seeds_holding": f"{int(r['seeds_holding'])}/{n}", "n_pairs": int(r["n_pairs"]),
                        "n_expected": int(r["n_expected"]), "t_lo": r.get("t_lo"), "t_hi": r.get("t_hi")})
            row["met"] = bool(r["diff"] > 0 and r["ci_lo"] > 0 and int(r["seeds_holding"]) >= need) if done else None
            row["status"] = "pending" if not done else None
        out.append(row)
    ps = [o["p"] if np.isfinite(o["p"]) else 1.0 for o in out] + [1.0] * max(0, cfg.holm_m - len(out))
    adj = A3.holm(ps)[:len(out)]
    for o, pa in zip(out, adj):
        o["p_holm"] = float(pa)
        if o["met"] is not None:
            o["met"] = bool(o["met"] and pa < 0.05)
            o["status"] = "met" if o["met"] else "not met"
    return pd.DataFrame(out)


def exploratory_rq3(cfg: RLStudyConfig, L: S5.LayoutV5, t: "pd.DataFrame") -> "pd.DataFrame":
    """v2single - arm for the arms that failed RP and passed RPX: exploratory, unadjusted p, not in Holm's family."""
    import pandas as pd
    out = []
    for arm in ARMS:
        if not exploratory(cfg, L, arm):
            continue
        q = t[(t.id == f"R-{arm}") & (t.metric == "success5") & (t.atmosphere == "all")] if len(t) else t
        row = {"id": f"X-{arm}", "question": f"v2single - {arm}, success@5km, A256 (exploratory: {arm} failed RP "
                                             "and passed the extended pilot RPX; not in RQ3's Holm family)",
               "diff_pp": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "p_unadjusted": np.nan, "seeds_holding": None,
               "n_pairs": 0, "n_expected": 0, "status": "pending"}
        if len(q):
            r = q.iloc[0].to_dict()
            row.update({"diff_pp": r["diff"], "ci_lo": r["ci_lo"], "ci_hi": r["ci_hi"], "p_unadjusted": float(r["p_boot"]),
                        "seeds_holding": f"{int(r['seeds_holding'])}/{int(r['n_seeds'])}", "n_pairs": int(r["n_pairs"]),
                        "n_expected": int(r["n_expected"]),
                        "status": "exploratory" if bool(r["complete"]) else "pending"})
        out.append(row)
    return pd.DataFrame(out)


def analyze(cfg: RLStudyConfig, L: S5.LayoutV5, J: Optional[S5.LayoutV5], L3=None, L4=None, log=print) -> dict:
    """RQ3's tables (05_analysis/tables/RQ3_*.csv) and summary (05_analysis/RQ3_summary.md) from whatever exists;
    a criterion is decided only on complete data (pending otherwise)."""
    import pandas as pd
    from jepa import analyze_v5 as A5
    from jepa import benchmark_v5 as B
    v5 = cfg.v5()
    df = load_flights(cfg, L, J, L3, L4)
    main = A5.main_flights(df) if len(df) else df
    cases = sorted({int(c["case"]) for c in S5.conditions(v5, "A256")})
    seeds = list(cfg.rl_seeds)
    out: Dict[str, object] = {"n_flights": int(len(df))}
    # ---- A256: R-each and estimates vs fnpeg ----
    rows = []
    if len(main):
        comp = A5.Comparer(main, A5.Boot5(L, cfg.n_boot), cases, A5.ClusterBoot(L, cfg.n_boot))
        for arm in ARMS:
            rows += comp(f"R-{arm}", A5.V2S, (arm, A5.OWN), "success5", seeds)
            for m in ("success5", "log_miss"):
                rows += comp(f"E-{arm}-fnpeg", (arm, A5.OWN), A5.FN, m, seeds)
    t = pd.DataFrame(rows)
    crit = criteria_rq3(cfg, L, t)
    L.write_csv(t, L.tables / "RQ3_paired_A256.csv")
    L.write_csv(crit, L.tables / "RQ3_criteria.csv")
    out["paired"], out["criteria"] = t, crit
    out["exploratory"] = exploratory_rq3(cfg, L, t)
    L.write_csv(out["exploratory"], L.tables / "RQ3_exploratory.csv")
    # ---- DS1 / AS1 ----
    for su, base in (("DS1", None), ("AS1", ("nominal",))):
        comps = [(f"{su[0]}-{a}", "v2single", a) for a in ARMS] + [(f"{su[0]}-{a}-fnpeg", a, "fnpeg") for a in ARMS]
        sh = A5.shift_analysis(v5, L, main, su, comps, base_atmospheres=base) if len(main) else {"comp": [], "drop": []}
        c = pd.DataFrame(sh["comp"])
        L.write_csv(c, L.tables / f"RQ3_{su}.csv")
        out[su] = c
    # ---- ORB-1 block 1: arm - fnpeg paired by case, success@10km with its Clopper-Pearson interval ----
    orows = []
    o = main[(main.suite == "ORB1") & (main.block == 0)] if len(main) else main
    for arm in ARMS:
        if not len(o):
            break
        piv = o.pivot_table(index="case", columns="key", values="success_10km", aggfunc="first")
        if arm not in piv:
            continue
        x = piv[arm].dropna()
        k, n = int(x.sum()), int(len(x))
        lo, hi = B.clopper_pearson(k, n)
        row = {"arm": arm, "n": n, "success_10km_%": 100 * k / max(n, 1), "cp_lo_%": 100 * lo, "cp_hi_%": 100 * hi}
        if "fnpeg" in piv:
            p = piv[[arm, "fnpeg"]].astype(float).dropna()
            r = A5.orb_paired(L, f"rl2|ORB1|{arm}-fnpeg|success_10km", p[arm].to_numpy(), p["fnpeg"].to_numpy(),
                              p.index.to_numpy() % len(seeds), cfg.n_boot)
            row.update({"n_pairs": int(len(p)), "diff_vs_fnpeg_pp": 100 * r["diff"], "ci_lo_pp": 100 * r["ci_lo"],
                        "ci_hi_pp": 100 * r["ci_hi"], "p_boot": r["p_boot"]})
        orows.append(row)
    out["ORB1"] = pd.DataFrame(orows)
    L.write_csv(out["ORB1"], L.tables / "RQ3_ORB1_block1.csv")
    # ---- R-gap: fast env vs full simulator on the same A256 cases, per seed ----
    grows = []
    a = main[main.suite == "A256"] if len(main) else main
    for arm in ARMS:
        for s in seeds:
            fe = L.policy_dir(arm, s) / "fast_eval.json"
            if not fe.exists() or not len(a):
                continue
            doc = json.loads(fe.read_text())
            full = a[(a.arm == arm) & (a.seed == s)].set_index(["case", "atmosphere"]).success_5km
            pairs = [(f, full.get((c, at))) for f, (c, at) in zip(doc["A256"]["success_by_case"], doc["A256_cases"])]
            pairs = [(f, v) for f, v in pairs if v is not None and not pd.isna(v)]
            if pairs:
                grows.append({"arm": arm, "seed": s, "n": len(pairs),
                              "fast_success_5km_%": 100 * float(np.mean([f for f, _ in pairs])),
                              "full_success_5km_%": 100 * float(np.mean([float(v) for _, v in pairs])),
                              "gap_pp": 100 * float(np.mean([f - float(v) for f, v in pairs]))})
    out["gap"] = pd.DataFrame(grows)
    L.write_csv(out["gap"], L.tables / "RQ3_gap.csv")
    # ---- gates, dropped arms, tuning, cost ----
    gates = []
    ri = L.checks / "RI_resid_fnpeg.json"
    if ri.exists():
        d = json.loads(ri.read_text())
        gates.append({"gate": "RI", "arm": "resid_fnpeg", "passed": d["passed"],
                      "reason": f"fast env: {d['fast']}; DEV48: {d['n_dev48_conditions']} conditions, "
                                f"{len(d['differing'])} differ"})
    for arm in ARMS:
        g = _rp_doc(L, arm)
        if g is not None:
            gates.append({"gate": "RP", "arm": arm, "passed": g["passed"], "reason": g.get("reason"),
                          "success_5km_fixed": g.get("success_5km_fixed"),
                          "fnpeg_success_5km_fixed": g.get("fnpeg_success_5km_fixed"), "bar": g.get("bar")})
        x = _rpx_doc(L, arm)
        if x is not None:
            gates.append({"gate": "RPX (full budget)", "arm": arm, "passed": x["passed"],
                          "reason": f"{x.get('reason')}; reached {x.get('reached')}; fixed-episode curve "
                                    + ", ".join(f"{c['decisions']:,}: {c['success_5km_fixed']:.3f}"
                                                for c in x.get("fixed_curve", [])),
                          "success_5km_fixed": x.get("success_5km_fixed"),
                          "fnpeg_success_5km_fixed": x.get("fnpeg_success_5km_fixed"), "bar": x.get("bar")})
    for arm, v in DROPPED.items():
        gates.append({"gate": "RV1/RV2 (notebook 1's design, dropped)", "arm": arm, "passed": False,
                      "reason": f"RV1 {v['RV1']}; RV2 {v['RV2']} ({v['pilot']})"})
    out["gates"] = pd.DataFrame(gates)
    L.write_csv(out["gates"], L.tables / "RQ3_gates.csv")
    trows = []
    for arm in ARMS:
        for k in range(cfg.tune_n):
            p = L.tuning / tmethod(arm) / f"config_{k:02d}.json"
            if p.exists():
                d = json.loads(p.read_text())
                trows.append({"arm": arm, "config_index": k, **d["score"],
                              **{f"cfg_{a}": b for a, b in d["config"].items()},
                              "chosen": (chosen(L, arm) or {}).get("config_index") == k})
    out["tuning"] = pd.DataFrame(trows)
    L.write_csv(out["tuning"], L.tables / "RQ3_tuning.csv")
    crows = []
    for arm in ARMS:
        rt, r = rates(L, arm), route(L, arm)
        row = {"arm": arm, "route": (r or {}).get("runtime")}
        for k, v in rt.items():
            row[f"{k}_decisions_per_s"] = v["decisions_per_s"]
            row[f"{k}_units_per_million"] = R2.units_per_million(v["decisions_per_s"], v["units_per_hour"])
        secs = [json.loads((L.policy_dir(arm, s) / "run.json").read_text()).get("train_s", 0.0) for s in seeds
                if (L.policy_dir(arm, s) / "run.json").exists()]
        row["train_hours"] = sum(secs) / 3600.0
        if r:
            row["train_units"] = row["train_hours"] * float(rt[r["runtime"]]["units_per_hour"])
        crows.append(row)
    out["cost"] = pd.DataFrame(crows)
    L.write_csv(out["cost"], L.tables / "RQ3_cost.csv")
    md = summary_md(cfg, out)
    L.atomic_bytes(L.analysis / "RQ3_summary.md", md.encode("utf-8"))
    out["summary"] = md
    log(f"[analysis] {len(df)} flights; RQ3: " + "; ".join(f"{r.id} {r.status}" for r in crit.itertuples()))
    return out


def summary_md(cfg: RLStudyConfig, out: dict) -> str:
    def tab(df):
        if df is None or not len(df):
            return "(none yet)\n"
        cols = list(df.columns)
        lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
        for r in df.itertuples(index=False):
            lines.append("| " + " | ".join(f"{v:.3g}" if isinstance(v, float) else str(v) for v in r) + " |")
        return "\n".join(lines) + "\n"
    c = out["criteria"][["id", "status", "diff_pp", "ci_lo", "ci_hi", "p", "p_holm", "seeds_holding", "gate_RP"]]
    return (f"# RQ3: RL baselines (notebook 2, jepa/plan_rl_v5.md)\n\n{out['n_flights']} flights read. "
            "Holm over m = 3 (an arm without a comparison enters with p = 1). `liu_track` is Liu-style: its "
            "deviations from Liu et al. (2026) are listed in plan_rl_v5.md section 1.3.\n\n"
            "## Criteria\n\n" + tab(c)
            + "\n## Exploratory (arms that failed RP and passed the extended pilot RPX; unadjusted, not in the "
              "Holm family)\n\n" + tab(out.get("exploratory"))
            + "\n## Gates (failed gates are results)\n\n" + tab(out["gates"])
            + "\n## ORB-1 block 1\n\n" + tab(out["ORB1"]) + "\n## Sim-to-sim gap\n\n" + tab(out["gap"])
            + "\n## Cost\n\n" + tab(out["cost"]))
