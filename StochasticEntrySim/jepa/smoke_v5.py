"""
jepa/smoke_v5.py

CPU smoke test of the v5 study code (jepa/v5_plan.md section 12). The numbers it produces are meaningless (tiny
models, 2 cases, 2 seeds); the checks are what matter. It runs on the v4 smoke fixtures (smoke_v4.py --work,
default <tempdir>/jepa_v4_smoke/: wide_smoke/, v3_smoke/, v4_smoke/; v4_smoke may be partial).

  1  suites: ORB-1, DEV48 and DS1 regenerated equal the frozen files byte for byte; a changed file fails the
     checksum; ORB-1's simulator knobs reproduce the dispersed drag and lift exactly (100 draws, the simulator's
     aero); the density hook scales drag by DS1's factor, and a flight with a profile of ones is bit-identical to
     one without; GPU passes route no CPU-pool flight (GPU_PASS_SKIPS_CPU)
  2  FNPEG: predictor = truth lands < 1 km; batched == unbatched; load relief keeps a 3 g case under 3.15 g; the
     filters converge to the true drag scale
  3  PredGuid: phase sequences on a short and a long case; the final-phase design is reproducible
  4  mppi_filtered with frozen filters == mppi_nominal, and the v5 runner == run_lockstep_v4 (v2, physics) (E0)
  5  every E3 model builds at full size within ±10% of R0, trains a few steps, saves, resumes
  6  fast env == physics_v4.rollout_physics; 10,000 batched random-policy steps; PPO, SAC, DAgger, BC -> PPO and
     full-simulator PPO run, checkpoint and resume
  7  legacy RL: a dummy SB3 zip is found by the inventory and flown
  8  job queue: finished flights skipped, a deleted result re-flown, ORB-1 stopping rule applied from files
  9  analysis: every table and figure, a complete criteria table, the bootstrap cache reused; the report builds
  10 the notebook runs end to end in smoke mode; a second run finds everything done in under a minute; v5 refuses
     to write under the v3, v4, data and legacy RL folders
  11 v3 and v4 notebooks and code unchanged (math_3d.py and data.py: only the default-off density hook added)
  12 resume equivalence (bit-exact on CPU): E3, PPO, SAC (with buffer shards), DAgger, BC -> PPO
  13 interrupted writes: stale .tmp, truncated train state, corrupt flight JSON, artifact missing behind .done
  14 simulated disconnects (JEPA_V5_FAIL_AT) in the notebook's long stages, then an unchanged rerun: completes,
     finished flights keep their modification times, outputs equal the uninterrupted run's
  15 outputs: every table and figure (PNG and PDF), curves, captions, logs, status files and the report exist
  16 RL on Kaggle: inputs exported, two notebook sessions (restore, resume), one platform enforced, import into
     Drive (finished runs only, checked, nothing overwritten)

    python jepa/smoke_v5.py [--work <dir>] [--workers N] [--only 1,2,...]
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))
os.chdir(SIM_DIR)
os.environ.setdefault("TMPDIR", tempfile.gettempdir())

import numpy as np
import pandas as pd
import torch

RESULTS = []
NOTES = {}


def check(name: str, ok: bool, detail: str = "") -> bool:
    RESULTS.append((name, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {detail}", flush=True)
    return bool(ok)


def section(t: str) -> None:
    print(f"\n{t}", flush=True)


# ============================================================================
# Children (separate processes, so a simulated disconnect kills a real process)
# ============================================================================
class _L:
    """A minimal layout for the training children (unguarded; temp folders only)."""

    def __init__(self, root):
        self.root = Path(root)
        self.models = self.root / "02_models"
        self.open_loop = self.root / "03_open_loop"

    def guard(self, p):
        return Path(p)

    def atomic_bytes(self, p, b):
        p = Path(p)
        p.parent.mkdir(parents=True, exist_ok=True)
        tmp = p.with_name(p.name + ".tmp")
        tmp.write_bytes(b)
        os.replace(tmp, p)

    def save_large(self, b, p):
        self.atomic_bytes(p, b)

    def write_json(self, p, o):
        self.atomic_bytes(p, json.dumps(o, default=str).encode())

    def write_csv(self, df, p):
        self.atomic_bytes(p, df.to_csv(index=False).encode())

    def save_fig(self, fig, p, dpi=150):
        import io
        b = io.BytesIO()
        fig.savefig(b, format=Path(p).suffix[1:], dpi=dpi)
        self.atomic_bytes(p, b.getvalue())

    def model_dir(self, a, s):
        return self.models / a / f"seed_{s}"

    def ckpt(self, a, s):
        return self.model_dir(a, s) / "ckpt.pt"

    def policy_dir(self, a, s):
        return self.root / "02_policies" / a / f"seed_{s}"


TINY_ARCH = dict(d_latent=16, enc_dim=32, enc_layers=2, enc_heads=2, enc_band=5, pred_dim=48, pred_layers=2,
                 pred_heads=2, act_chunk=4, mlp_ratio=2)
E3_TRAIN = dict(steps=8, batch=8, warmup=2, curriculum=((0.0, 8), (0.5, 16)), log_every=4, save_every_s=1e9,
                save_every_steps=3, amp=False, val_rows=8)
E3_SIZES = {"mlp_direct": {"width": 32, "depth": 2}, "gru_direct": {"hidden": 16},
            "gru_pi": {"hidden": 16, "prober_hidden": 16}}
RL_TINY = dict(n_envs=16, hidden=(32, 32), n_steps=8, epochs=2, minibatch=64, budget_decisions=512, budget_hours=1.0,
               buffer=4096, batch=64, learning_starts=64, updates_per_step=4.0, dagger_iters=2, dagger_rollouts=4,
               dagger_epochs=2, bc_episodes=4, bc_epochs=2, sup_batch=256, save_every_s=1e9, save_every_updates=2,
               curve_points=2, env_pool=64, full_episodes=4, full_per_update=2, full_workers=1,
               critic_warmup_updates=1, kl_coef=0.1)
RL_EXTRA = {"sac_hist": {"budget_decisions": 2400, "save_every_updates": 40}, "ppo_matched": {"save_every_updates": 1}}
RL_FAIL = {"ppo_hist": "rl_train:2", "bcppo_fnpeg": "rl_train:2", "sac_hist": "rl_train:80",
           "dagger_fnpeg": "rl_train:1", "ppo_matched": "rl_train:1"}


def _e3_cfg():
    class C5:
        pass
    c = C5()
    c.e3_train, c.e3_sizes, c.r0_arch = dict(E3_TRAIN), dict(E3_SIZES), dict(TINY_ARCH)
    c.lr_mults, c.lr_pilot_frac, c.ttg_s = (1.0, 3.0), 0.5, (100, 50)
    return c


def child(kind: str, root: str, arm: str, wide: str) -> None:
    if kind == "e3":
        from jepa import models_v5 as M5
        from jepa import train_v2 as T2
        eps = T2.load_episodes(wide)[:60]
        data = T2.prepare(eps, 16, "cpu", 25, split_seed=0)
        M5.train_arm_seed(_e3_cfg(), data, _L(root), arm, 0, "cpu", 1.0, log=print)
    else:
        from jepa import rl_v5 as RL
        cfg = dict(RL_TINY, **RL_EXTRA.get(arm, {}))
        info = RL.train_rl_arm(_L(root), arm, 0, RL.RLCfg(**cfg), "cpu", log=print)
        print("DONE", json.dumps({k: v for k, v in info.items() if k != "cfg"}, default=str))
    sys.exit(0)


def _run_child(kind, root, arm, wide, fail=None):
    env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
    env.pop("JEPA_V5_FAIL_AT", None)
    if fail:
        env["JEPA_V5_FAIL_AT"] = fail
    return subprocess.run([sys.executable, __file__, "--child", kind, str(root), arm, str(wide)], env=env,
                          capture_output=True, text=True, encoding="utf-8", errors="replace")


# ============================================================================
# 1. Suites
# ============================================================================
def t1_suites(tmp: Path, workers: int) -> None:
    section("1. suites")
    from jepa import benchmark_v5 as B
    import rl_env
    from jepa.data import apply_domain
    import math_3d
    sums = B.verify()
    check("frozen suites verify against SHA256SUMS", len(sums) > 5, f"({len(sums)} files)")
    spec = json.loads((B.BENCH_DIR / "ORB1_spec.json").read_text())
    fp = spec["footprint"]
    t0 = time.time()
    orb, _ = B.build_orb1(fp["sigma90"]["downrange_km"], fp["sigma0"]["downrange_km"], B.ORB1Spec(), 12, workers,
                          log=lambda *a: None)
    got = B._csv_bytes(orb).decode().splitlines()
    have = (B.BENCH_DIR / "ORB1.csv").read_text().splitlines()[:13]
    check("ORB-1 regenerated (first 12 cases, footprints flown again) == the frozen ORB1.csv bytes", got == have,
          f"({time.time() - t0:.0f}s)")
    t0 = time.time()
    dev, _ = B.build_dev48(workers, log=lambda *a: None)
    check("DEV48 regenerated == the frozen DEV48.csv bytes",
          B._csv_bytes(dev) == (B.BENCH_DIR / "DEV48.csv").read_bytes(), f"({time.time() - t0:.0f}s)")
    bad = tmp / "bench_bad"
    shutil.copytree(B.BENCH_DIR, bad)
    b = bytearray((bad / "DEV48.csv").read_bytes())
    b[200] = ord("7") if b[200] != ord("7") else ord("8")
    (bad / "DEV48.csv").write_bytes(bytes(b))
    try:
        B.verify(bad)
        ok = False
    except B.ChecksumError:
        ok = True
    check("a changed benchmark file fails the checksum", ok)
    # ORB-1's knobs: the simulator's drag and lift with the dispersed values == nominal x rho*CD/m and rho*CL/m
    env = rl_env.OrionEntryEnv(config_path="configs/default.json", controller_mode="policy", dispersion=False,
                               max_episode_steps=10, fast_atmosphere=True)
    nom = dict(env.params)
    x = [6_378_137.0 + 60_000.0, 0.1, 0.2, 6500.0, math.radians(-2.0), 0.7]
    worst = 0.0
    for r in B.orb1_dispersions(B.ORB1Spec(), 100):
        p0 = apply_domain(dict(nom), {"density_scale": 1.0, "lift_scale": 1.0, "mass_scale": 1.0})
        p1 = apply_domain(dict(nom), B.domain_of(r))
        a0 = math_3d.nominal_aero_forces_from_state(x, 0.3, p0)
        a1 = math_3d.nominal_aero_forces_from_state(x, 0.3, p1)
        d0, l0 = a0["drag_mag_N"] / p0["mass_kg"], a0["lift_mag_N"] / p0["mass_kg"]
        d1, l1 = a1["drag_mag_N"] / p1["mass_kg"], a1["lift_mag_N"] / p1["mass_kg"]
        want_d = r["rho_mult"] * r["cd_mult"] / r["mass_mult"]
        want_l = r["rho_mult"] * r["cl_mult"] / r["mass_mult"]
        worst = max(worst, abs(d1 / d0 / want_d - 1.0), abs(l1 / l0 / want_l - 1.0))
    check("ORB-1 knobs reproduce the dispersed drag and lift (100 draws, simulator aero)", worst < 1e-12,
          f"(max relative error {worst:.1e})")
    t1_ds1(nom)
    t1_as1(nom)
    t1_gpu_pass()
    t1_gpu_lock()


def t1_gpu_pass() -> None:
    """GPU_PASS_SKIPS_CPU: run_jobs on a GPU pass flies no CPU-pool flight and counts it; a CPU pass no MPPI flight."""
    from jepa import study_v5 as S5
    calls = []
    saved = (S5.run_cpu, S5.run_mppi, S5.done_ids)
    S5.run_cpu = lambda cfg, jobs, L, *a, **k: calls.append(("cpu", len(jobs))) or 0
    S5.run_mppi = lambda cfg, jobs, L, *a, **k: calls.append(("mppi", len(jobs))) or 0
    S5.done_ids = lambda L, jobs: set()
    try:
        jobs = [{"model": {"kind": "fnpeg"}}] * 3 + [{"model": {"kind": "policy"}}] * 2 + [{"model": {"kind": "v2"}}] * 4
        before = dict(S5.WAITING)
        msgs = []
        S5.run_jobs(None, jobs, None, "cuda", 1, 1, True, msgs.append, "t", can_cpu=False)
        gpu_pass, w1 = list(calls), dict(S5.WAITING)
        calls.clear()
        S5.run_jobs(None, jobs, None, "cpu", 1, 1, False, msgs.append, "t")
        cpu_pass, w2 = list(calls), dict(S5.WAITING)
        calls.clear()
        S5.run_jobs(None, jobs, None, "cuda", 1, 1, True, msgs.append, "t")
        both = list(calls)
    finally:
        S5.run_cpu, S5.run_mppi, S5.done_ids = saved
        S5.WAITING.update(before)
    nb = (SIM_DIR / "ReEntryAI_JEPA_v5.ipynb").read_text(encoding="utf-8")
    check("GPU pass skips the CPU-pool flights (5 counted as waiting), a CPU pass the MPPI flights (4), the default "
          "flies both; every notebook flight stage passes can_cpu",
          gpu_pass == [("mppi", 4)] and w1["cpu"] - before["cpu"] == 5 and cpu_pass == [("cpu", 9)]
          and w2["gpu"] - w1["gpu"] == 4 and both == [("cpu", 9), ("mppi", 4)]
          and nb.count("S5.run_jobs(") == nb.count("can_cpu=CAN_CPU") == 10 and "GPU_PASS_SKIPS_CPU = True" in nb,
          f"(gpu pass {gpu_pass}, cpu pass {cpu_pass}, default {both})")


def t1_gpu_lock() -> None:
    """GPU_TYPE: the first GPU session fixes 00_study/gpu_platform.json; another GPU type, or one that is not the
    required type, is refused without touching the file."""
    from jepa import study_v5 as S5
    with tempfile.TemporaryDirectory() as td:
        L = S5.LayoutV5(td).make()
        f = L.root / S5.GPU_PLATFORM
        refused = []

        def refuses(*a):
            before = f.read_bytes() if f.exists() else None
            try:
                S5.check_gpu_platform(L, *a, log=lambda m: None)
                return False
            except RuntimeError:
                return (f.read_bytes() if f.exists() else None) == before

        refused.append(refuses("NVIDIA L4", "T4"))                       # not a T4, nothing recorded yet
        none_yet = not f.exists()
        d1 = S5.check_gpu_platform(L, "Tesla T4", "T4", {"torch": "2.4"}, log=lambda m: None)
        d2 = S5.check_gpu_platform(L, "Tesla T4", "T4", {"torch": "2.5"}, log=lambda m: None)   # torch update: allowed
        refused.append(refuses("NVIDIA L4", None))                       # no type required, but the study is on T4
        refused.append(refuses("NVIDIA A100-SXM4-40GB", "T4"))
        doc = json.loads(f.read_text())
    nb = (SIM_DIR / "ReEntryAI_JEPA_v5.ipynb").read_text(encoding="utf-8")
    check("GPU lock: the first T4 session records it, a later T4 session (even on newer torch) passes, an L4 or A100 "
          "is refused with the file unchanged; the notebook requires a T4 and opens on one",
          all(refused) and none_yet and d1["gpu"] == d2["gpu"] == doc["gpu"] == "Tesla T4"
          and len(doc["sessions"]) == 2 and 'GPU_TYPE = \\"T4\\"' in nb and "S5.check_gpu_platform(" in nb
          and '"gpuType": "T4"' in nb, f"(refused {refused}, sessions {len(doc['sessions'])})")


def t1_as1(nom: dict) -> None:
    """AS1 and the aero-schedule hook (v5_plan.md section 3.4b)."""
    from jepa import benchmark_v5 as B
    from jepa import physics_v4 as PH
    from jepa.data import apply_domain
    import math_3d
    spec = B.AS1Spec()
    doc = B.build_as1_doc(spec, workers=max(1, min(8, (os.cpu_count() or 2) - 1)))
    check("AS1 regenerated == the frozen AS1.json bytes (schedules and reachability flags)",
          (json.dumps(doc, indent=1) + "\n").encode() == (B.BENCH_DIR / "AS1.json").read_bytes(),
          f"(reachable: {[doc['vehicles'][v]['n_feasible'] for v in B.AS1_VEHICLES]} of {len(doc['cases'])})")
    conds = B.load_suite("AS1")
    check("AS1 loads 64 cases x 2 vehicles in nominal air; every domain value a float; its copy of the simulator's "
          "schedule equals the port's", len(conds) == 128 and {c["atmosphere"] for c in conds} == set(B.AS1_VEHICLES)
          and all(isinstance(float(v), float) for c in conds for v in c["domain"].values())
          and all(isinstance(c["feasible_as1"], bool) for c in conds)
          and [tuple(r) for r in spec.ours] == [tuple(float(v) for v in r) for r in PH.TRIM_SCHEDULE])
    # the hook: at a speed where the schedules differ, lift/drag of the true vehicle == the AS1 schedule's L/D; a
    # later domain without the id puts the parameters back exactly
    base = {"density_scale": 1.0, "lift_scale": 1.0, "mass_scale": 1.0}
    x = [6_378_100.0 + 60_000.0, 0.1, 0.2, 4500.0, math.radians(-2.0), 0.7]
    f0 = math_3d.nominal_aero_forces_from_state(x, 0.3, apply_domain(dict(nom), base))
    ok, shown = True, {}
    for v in B.AS1_VEHICLES:
        f1 = math_3d.nominal_aero_forces_from_state(x, 0.3, apply_domain(dict(nom), B.as1_domain(v)))
        ld0, ld1 = f0["lift_mag_N"] / f0["drag_mag_N"], f1["lift_mag_N"] / f1["drag_mag_N"]
        sch = dict((r[0], r) for r in spec.schedule(v))
        lo, hi = min(r[2] for r in sch.values()), max(r[2] for r in sch.values())
        shown[v] = (round(ld0, 4), round(ld1, 4), round(f1["drag_mag_N"] / f0["drag_mag_N"], 4))
        ok &= ld1 < ld0 - 0.02 and lo - 1e-9 <= ld1 <= hi + 1e-9
        ok &= (abs(f1["drag_mag_N"] / f0["drag_mag_N"] - 1.0) < 1e-12) == (v == "bibb_ld")
    fresh = apply_domain(dict(nom), base)
    reused = apply_domain(apply_domain(dict(nom), B.as1_domain("bibb_ld_cd")), base)
    check("aero hook: the true vehicle flies the AS1 schedule (lower L/D; drag changes only for bibb_ld_cd); a later "
          "domain without it restores the parameters exactly", ok and reused == fresh
          and "aero_schedule_before" not in reused, f"(L/D ours, AS1, drag ratio: {shown})")
    ic = dict(B.ORB1Spec().nominal_ic())
    a = B.fly_constant_bank_sim(ic, 45.0, dict(base))
    b = B.fly_constant_bank_sim(ic, 45.0, B.as1_domain("bibb_ld"))
    check("a full flight on the AS1 vehicle lands shorter than the simulator's own vehicle at the same bank",
          b["downrange_km"] < a["downrange_km"] - 10.0, f"({a['downrange_km']:.0f} km vs {b['downrange_km']:.0f} km)")


def t1_ds1(nom: dict) -> None:
    """DS1 and the density hook (v5_plan.md section 3.4)."""
    from jepa import benchmark_v5 as B
    from jepa.data import apply_domain
    import math_3d
    doc = B.build_ds1_doc()
    check("DS1 regenerated == the frozen DS1.json bytes",
          (json.dumps(doc, indent=1) + "\n").encode() == (B.BENCH_DIR / "DS1.json").read_bytes())
    conds = B.load_suite("DS1")
    ids_ok = all(B.ds1_profile(c["domain"]["density_profile_id"])["factor"]
                 == doc["atmospheres"][c["atmosphere"]]["profiles"][str(int(c["case"]))] for c in conds)
    check("DS1 loads 64 cases x 4 atmospheres; each domain names its profile by a number (every domain value a "
          "float, as flight records need) that resolves to that case's profile", len(conds) == 256 and ids_ok
          and all(isinstance(float(v), float) for c in conds for v in c["domain"].values())
          and {c["atmosphere"] for c in conds} == set(B.DS1_ATMOSPHERES))
    # the hook: drag with the profile / drag without == the profile's factor at that altitude; absent = untouched
    worst, h_list = 0.0, (5_000.0, 41_000.0, 55_000.0, 69_000.0, 85_000.0, 130_000.0)
    for c in conds[::37]:
        for h in h_list:
            x = [6_378_137.0 + h, 0.1, 0.2, 6500.0, math.radians(-2.0), 0.7]
            p0 = apply_domain(dict(nom), {k: v for k, v in c["domain"].items() if k != "density_profile_id"})
            p1 = apply_domain(dict(nom), c["domain"])
            want = math_3d.density_profile_factor(p1["density_profile"], max(0.0, math_3d.constants.altitude(x[0])))
            r = (math_3d.nominal_aero_forces_from_state(x, 0.3, p1)["drag_mag_N"]
                 / math_3d.nominal_aero_forces_from_state(x, 0.3, p0)["drag_mag_N"])
            worst = max(worst, abs(r / want - 1.0))
    reused = apply_domain(apply_domain(dict(nom), conds[0]["domain"]), {"density_scale": 1.0, "lift_scale": 1.0,
                                                                         "mass_scale": 1.0})
    check("density hook: drag scales by the profile's factor at each altitude; a later domain without a profile "
          "clears it", worst < 1e-12 and "density_profile" not in reused, f"(max relative error {worst:.1e})")
    # bit-exact: a whole full-simulator flight with a profile of ones equals the flight without one
    ic = dict(B.ORB1Spec().nominal_ic())
    flat = {"alt_m": list(doc["alt_m"]), "factor": [1.0] * len(doc["alt_m"])}
    base_dom = {"density_scale": 1.0, "lift_scale": 1.0, "mass_scale": 1.0}
    t0 = time.time()
    a = B.fly_constant_bank_sim(ic, 45.0, dict(base_dom))
    b = B.fly_constant_bank_sim(ic, 45.0, dict(base_dom, density_profile=flat))
    s = B.fly_constant_bank_sim(ic, 45.0, next(c["domain"] for c in conds if c["atmosphere"] == "shear_dense_low"))
    keys = sorted(k for k in a if isinstance(a[k], (int, float)))
    check("a full flight with a profile of ones is bit-identical to one without; a shear profile changes it",
          all(a[k] == b[k] for k in keys) and any(a[k] != s[k] for k in keys),
          f"(compared {keys}; {time.time() - t0:.0f}s)")


# ============================================================================
# 2-3. Guidance laws
# ============================================================================
def _a256(case: int, atm: str = "nominal") -> dict:
    from jepa import study_v5 as S5
    cfg = S5.StudyConfigV5()
    return next(c for c in S5.conditions(cfg, "A256") if int(c["case"]) == case and c["atmosphere"] == atm)


def t2_fnpeg(tmp: Path) -> None:
    section("2. FNPEG")
    from jepa import fly_v5 as F5
    from jepa import study_v5 as S5
    t0 = time.time()
    orc = {"oracle": True, "filters": False, "load_relief": False}
    jobs = [S5.guidance_job(_a256(c), "smoke_fnpeg", "smoke", "fnpeg", orc) for c in (0, 16)]
    both = F5.run_guidance_lockstep(jobs, F5.make_law("fnpeg", orc), verbose=False)
    one = [F5.run_guidance_lockstep([j], F5.make_law("fnpeg", orc), verbose=False)[0] for j in jobs]
    miss = [r["result"]["miss_km"] for r in both]
    check("predictor = truth: lands < 1 km", both[0]["result"]["outcome"] == "drogue" and miss[0] < 1.0,
          f"(misses {[round(m, 3) for m in miss]} km, {time.time() - t0:.0f}s for 4 flights)")
    same = all(a["result"]["miss_km"] == b["result"]["miss_km"] and
               all(np.array_equal(a["telemetry"][k], b["telemetry"][k], equal_nan=True) for k in a["telemetry"])
               for a, b in zip(both, one))
    check("batched == unbatched (results and telemetry bit-identical)", same)
    NOTES["fnpeg_ms"] = float(np.median([r["result"]["plan_ms_median"] for r in one]))
    j3 = S5.guidance_job(_a256(16), "smoke_fnpeg3", "smoke", "fnpeg", {"oracle": True, "filters": False, "amax_g": 3.0})
    r3 = F5.run_guidance_lockstep([j3], F5.make_law("fnpeg", j3["model"]["cfg"]), verbose=False)[0]["result"]
    check("load relief at amax 3 g keeps the peak under 3.15 g", r3["peak_g"] <= 3.15,
          f"(peak {r3['peak_g']:.2f} g, miss {r3['miss_km']:.2f} km)")
    jf = S5.guidance_job(_a256(0, "mismatch_high"), "smoke_fnpeg_f", "smoke", "fnpeg", {})
    rf = F5.run_guidance_lockstep([jf], F5.make_law("fnpeg", {}), verbose=False)[0]
    true = rf["result"].get("dom_density_scale", 1.0) / rf["result"].get("dom_mass_scale", 1.0)
    t = rf["telemetry"]
    on = np.nonzero(t["dec_a_g"] >= 0.05)[0]
    k = np.nonzero(t["dec_t_s"] >= t["dec_t_s"][on[0]] + 150.0)[0] if len(on) else []
    kd = float(t["dec_K_D"][k[0]]) if len(k) else float("nan")
    check("filters converge to the true drag scale (within 3%, 150 s after 0.05 g; gate V3's rule)",
          abs(kd / true - 1.0) < 0.03, f"(K_D {kd:.4f} vs true {true:.4f})")


def t3_predguid(tmp: Path) -> None:
    section("3. PredGuid")
    from jepa import fly_v5 as F5
    from jepa import guidance_predguid as PG
    from jepa import study_v5 as S5
    conds = {c["case"]: c for c in S5.conditions(S5.StudyConfigV5(), "A256") if c["atmosphere"] == "nominal"}
    dr = sorted(conds, key=lambda k: float(conds[k].get("target_downrange_km", conds[k].get("downrange_km", k))))
    short, long_ = dr[0], dr[-1]
    t0 = time.time()
    jobs = [S5.guidance_job(conds[k], "smoke_predguid", "smoke", "predguid", {}) for k in (short, long_)]
    rs = F5.run_guidance_lockstep(jobs, F5.make_law("predguid", {}), verbose=False)
    ph = [r["result"]["predguid_phases"] for r in rs]
    check("phase sequences logged and starting with the initial roll",
          all(p.startswith("initial_roll") and len(p.split("-")) >= 2 for p in ph),
          f"(case {short}: {ph[0]}; case {long_}: {ph[1]}; misses {[round(r['result']['miss_km'], 2) for r in rs]} km, "
          f"{time.time() - t0:.0f}s)")
    NOTES["predguid_ms"] = float(np.median([r["result"]["plan_ms_median"] for r in rs]))
    import rl_env
    env = rl_env.OrionEntryEnv(config_path="configs/default.json", controller_mode="policy", dispersion=False,
                               max_episode_steps=10, fast_atmosphere=True)
    tp, tl = float(env.target_phi_rad), float(env.target_lam_rad)
    cfg = PG.PredGuidConfig()
    a = PG.design_final_phase(conds[short]["ic"], tp, tl, cfg)
    b = PG.design_final_phase(conds[short]["ic"], tp, tl, cfg)
    same = a.summary() == b.summary() and all(np.array_equal(a.ftable[k], b.ftable[k]) for k in a.ftable)
    check("final-phase reference and gain design reproducible", same, f"({a.summary()['range_ref_km']:.1f} km ref)")


# ============================================================================
# 4. E0 (runner equivalence, frozen filters)
# ============================================================================
def t4_e0(tmp: Path, L3) -> None:
    section("4. mppi_filtered frozen == mppi_nominal; v5 runner == run_lockstep_v4 (E0)")
    from jepa import study_v5 as S5
    cfg = S5.StudyConfigV5.smoke()
    L = S5.LayoutV5(tmp / "e0", read_only=[L3.root]).make()
    t0 = time.time()
    r = S5.gate_e0(cfg, L, L3, log=lambda *a: None)
    df = pd.read_csv(L.checks / "E0.csv")
    for chk, g in df.groupby("check"):
        check(chk, bool(g.identical.all()), f"({int(g.identical.sum())}/{len(g)} flights identical)")
    check("E0 gate passes", r["pass"] is True, f"({time.time() - t0:.0f}s)")


# ============================================================================
# 5 + 12 + 13b. E3 models
# ============================================================================
def t5_e3(tmp: Path, wide: Path) -> None:
    section("5. E3 models: full size, a few steps, save, resume (12: bit-exact; 13: truncated train state)")
    from jepa import models_v5 as M5
    r0 = M5.r0_param_count()
    sizes = {a: M5.n_params(M5.build_e3(a)) for a in M5.ARMS_E3}
    check("every E3 model at full size within ±10% of R0's parameters",
          all(abs(n / r0 - 1) <= 0.10 for n in sizes.values()),
          f"(R0 {r0:,}; " + ", ".join(f"{a} {n:,} ({100 * (n / r0 - 1):+.2f}%)" for a, n in sizes.items()) + ")")
    NOTES["e3_sizes"] = {**sizes, "R0": r0}
    for arm in M5.ARMS_E3:
        a, b = tmp / "e3" / arm / "a", tmp / "e3" / arm / "b"
        t0 = time.time()
        r = _run_child("e3", a, arm, wide)
        ok = r.returncode == 0
        ta = time.time() - t0
        r1 = _run_child("e3", b, arm, wide, fail="e3_train:3")
        r2 = _run_child("e3", b, arm, wide)
        ok = ok and r1.returncode == 17 and r2.returncode == 0 and "resumed at step 3" in r2.stdout
        if not ok:
            print(r.stdout[-1500:], r.stderr[-2500:], r1.stdout[-800:], r2.stdout[-1500:], r2.stderr[-2500:])
        same, hist = _same_e3(_L(a), _L(b), arm) if ok else (False, False)
        check(f"12. {arm}: killed after the step-3 checkpoint, resumed == uninterrupted",
              ok and same and hist, f"(weights {same}, history {hist}; {ta:.0f}s per run)")
    # 13b: the newest train state truncated -> falls back to .prev, same result
    arm = "gru_direct"
    c = tmp / "e3" / arm / "c"
    r1 = _run_child("e3", c, arm, wide, fail="e3_train:6")
    st = _L(c).model_dir(arm, 0) / "train_state.pt"
    prev = st.with_name("train_state.prev.pt")
    ok = r1.returncode == 17 and st.exists() and prev.exists()
    if ok:
        data = st.read_bytes()
        st.write_bytes(data[: len(data) // 3])
    r2 = _run_child("e3", c, arm, wide)
    same, hist = _same_e3(_L(tmp / "e3" / arm / "a"), _L(c), arm) if r2.returncode == 0 else (False, False)
    check("13. truncated newest train state: resumes from .prev (logged) and ends identical",
          ok and r2.returncode == 0 and "does not load" in r2.stdout and same,
          f"(weights {same}, history {hist})")


def _same_e3(La, Lb, arm):
    a = torch.load(La.ckpt(arm, 0), weights_only=False)
    b = torch.load(Lb.ckpt(arm, 0), weights_only=False)
    if arm == "r0_nojepa":
        same = all(torch.equal(a["wm"][1][k], b["wm"][1][k]) for k in a["wm"][1]) and \
            all(torch.equal(a["pi"][1][k], b["pi"][1][k]) for k in a["pi"][1])
    else:
        same = all(torch.equal(a["state"][k], b["state"][k]) for k in a["state"])
    ha = pd.read_csv(La.model_dir(arm, 0) / "train_log.csv")
    hb = pd.read_csv(Lb.model_dir(arm, 0) / "train_log.csv")
    cols = [c for c in ("step", "state", "acc", "val") if c in ha]
    return same, bool(np.array_equal(ha[cols].values, hb[cols].values))


# ============================================================================
# 6 + 12. Fast env and RL
# ============================================================================
def t6_rl(tmp: Path, wide: Path) -> None:
    section("6. fast env and RL methods (12: checkpoint/resume bit-exact)")
    from jepa import rlenv_fast_v5 as RF
    err = RF.check_matches_rollout_physics()
    check("fast env == physics_v4.rollout_physics for the same commands", err == 0.0, f"(max |diff| {err:.1e})")
    env = RF.FastEntryEnv(RF.FastEnvCfg(n_envs=1000, seed=1, pool=1024))
    g = torch.Generator().manual_seed(0)
    t0 = time.time()
    n_done = 0
    for _ in range(10):
        _, r, d, _ = env.step(torch.rand(1000, generator=g) * 2 - 1)
        n_done += int(d.sum())
    el = time.time() - t0
    check("a random policy runs 10,000 batched steps", bool(torch.isfinite(r).all()),
          f"({el:.1f}s, {10_000 / el:.0f} steps/s on the CPU, {n_done} episodes ended)")
    NOTES["fast_env_steps_per_s_cpu"] = round(10_000 / el)
    for arm in ("ppo_hist", "sac_hist", "dagger_fnpeg", "bcppo_fnpeg", "ppo_matched"):
        a, b = tmp / "rl" / arm / "a", tmp / "rl" / arm / "b"
        t0 = time.time()
        r = _run_child("rl", a, arm, wide)
        ta = time.time() - t0
        r1 = _run_child("rl", b, arm, wide, fail=RL_FAIL[arm])
        r2 = _run_child("rl", b, arm, wide)
        ok = r.returncode == 0 and r1.returncode == 17 and r2.returncode == 0
        if not ok:
            print(r.stdout[-1500:], r.stderr[-2500:], r1.stdout[-800:], r1.stderr[-1500:], r2.stdout[-1500:],
                  r2.stderr[-2500:])
        w = o = c = False
        if ok:
            La, Lb = _L(a), _L(b)
            pa = torch.load(La.policy_dir(arm, 0) / "policy.pt", weights_only=False)
            pb = torch.load(Lb.policy_dir(arm, 0) / "policy.pt", weights_only=False)
            w = all(torch.equal(pa["state"][k], pb["state"][k]) for k in pa["state"])
            o = True
            sa, sb = La.policy_dir(arm, 0) / "train_state.pt", Lb.policy_dir(arm, 0) / "train_state.pt"
            if sa.exists() and sb.exists():
                xa, xb = torch.load(sa, weights_only=False), torch.load(sb, weights_only=False)
                for k in ("opt", "oq", "oa", "oal"):
                    if k in xa:
                        for u, v in zip(xa[k]["state"].values(), xb[k]["state"].values()):
                            o &= all(torch.equal(u[q], v[q]) for q in u)
            ca = pd.read_csv(La.policy_dir(arm, 0) / "learning_curve.csv")
            cb = pd.read_csv(Lb.policy_dir(arm, 0) / "learning_curve.csv")
            cols = [x for x in ca.columns if x != "sec"]
            c = ca[cols].equals(cb[cols])
        tag = "12. " if arm in ("ppo_hist", "sac_hist", "dagger_fnpeg", "bcppo_fnpeg") else ""
        check(f"{tag}{arm}: runs, checkpoints, killed and resumed == uninterrupted", ok and w and o and c,
              f"(weights {w}, optimizer {o}, curve {c}; {ta:.0f}s per run)")
        if arm == "sac_hist" and ok:
            shards = sorted((_L(b).policy_dir(arm, 0) / "buffer").glob("shard_*.npz"))
            check("12. SAC replay buffer saved as append-only shards", len(shards) >= 2, f"({len(shards)} shards)")


# ============================================================================
# 7. Legacy RL
# ============================================================================
def make_legacy(root: Path) -> Path:
    """MyDrive/ReEntryAI/models/<name>.zip with a dummy SB3 PPO on rl_env's 13-feature observation."""
    p = root / "models" / "dummy_ppo_policy.zip"
    ck = root / "checkpoints" / "run_x" / "ppo_1000_steps.zip"
    if not ck.exists():                     # a temp-folder cleanup can remove one file and leave the other
        ck.parent.mkdir(parents=True, exist_ok=True)
        ck.write_bytes(b"checkpoint")
    if p.exists():
        return p
    import gymnasium as gym
    from stable_baselines3 import PPO

    class Env(gym.Env):
        observation_space = gym.spaces.Box(-np.inf, np.inf, (13,), np.float32)
        action_space = gym.spaces.Box(-1.0, 1.0, (1,), np.float32)

        def reset(self, seed=None, options=None):
            return np.zeros(13, np.float32), {}

        def step(self, a):
            return np.zeros(13, np.float32), 0.0, True, False, {}

    p.parent.mkdir(parents=True, exist_ok=True)
    PPO("MlpPolicy", Env(), seed=0, n_steps=8, batch_size=8).save(str(p))
    return p


def t7_legacy(tmp: Path, legacy: Path) -> None:
    section("7. legacy RL")
    from jepa import study_v5 as S5
    cfg = S5.StudyConfigV5.smoke()
    L = S5.LayoutV5(tmp / "legacy_study", read_only=[legacy]).make()
    fins = S5.legacy_inventory(cfg, L, str(legacy), log=lambda *a: None)
    inv = pd.read_csv(L.study / "legacy_rl_inventory.csv")
    check("the dummy SB3 zip is found by the inventory (and the checkpoint listed)",
          len(fins) == 1 and (inv.kind == "checkpoints").any() and (inv.kind == "not found").any(),
          f"({len(inv)} rows: {inv.kind.value_counts().to_dict()})")
    jobs = [j for j in S5.legacy_jobs(cfg, L, fins)][:2]
    n = S5.run_cpu(cfg, jobs, L, 1, lambda *a: None, "legacy")
    r = [json.loads(L.result_path(S5.F5.flight_id_v5(j)).read_text()) for j in jobs]
    check("the legacy policy is flown with its own settings (1 Hz, rl_env's observation)",
          n == 2 and all(x["policy_period_s"] == 1 for x in r),
          f"({[x['suite'] for x in r]}, misses {[round(x['miss_km'], 1) for x in r]} km)")


# ============================================================================
# 8. Job queue
# ============================================================================
def t8_queue(tmp: Path, L3) -> None:
    section("8. job queue")
    from jepa import fly_v5 as F5
    from jepa import study_v5 as S5
    cfg = S5.StudyConfigV5.smoke()
    L = S5.LayoutV5(tmp / "queue", read_only=[L3.root]).make()
    two = [S5.guidance_job(c, "smoke_q", "smoke", "predguid", {}) for c in S5.conditions(cfg, "A256")[:2]]
    S5.run_cpu(cfg, two, L, 1, lambda *a: None, "q")
    p0, p1 = (L.result_path(F5.flight_id_v5(j)) for j in two)
    m0 = p0.stat().st_mtime_ns
    p1.unlink()
    n = S5.run_cpu(cfg, two, L, 1, lambda *a: None, "q")
    check("finished flight skipped, deleted flight re-flown", n == 1 and p0.stat().st_mtime_ns == m0 and p1.exists())
    # ORB-1 stopping rule from files: fake block-1 results for an arm
    c3 = S5.StudyConfigV5.smoke()
    c3.orb1_per_block = 3
    by_block = {}
    for c in S5.conditions(c3, "ORB1"):
        by_block.setdefault(int(c["block"]), []).append(c)
    out = []
    for n_bad in (1, 2):
        Lq = S5.LayoutV5(tmp / f"orb_{n_bad}").make()
        for i, c in enumerate(by_block[0]):
            j = S5.arm_job(c3, Lq, L3, c, "predguid", "ORB1")
            rec = {"flight_id": F5.flight_id_v5(j), "outcome": "drogue", "miss_km": 25.0 if i < n_bad else 1.0}
            Lq.write_json(Lq.result_path(rec["flight_id"]), rec)
        ok, why = S5.orb1_eligible_files(c3, Lq, L3, "predguid", 1, by_block)
        js = S5.orb1_jobs(c3, Lq, L3, ["predguid"], [1], "ORB1", log=lambda *a: None)
        out.append((ok, len(js), why))
    check("ORB-1 stopping rule from files: 1 miss > 10 km keeps the arm, 2 stop it",
          out[0][0] and out[0][1] == len(by_block[1]) and not out[1][0] and out[1][1] == 0,
          f"({out[0][2]}; {out[1][2]})")
    Lw = S5.LayoutV5(tmp / "orb_wait").make()
    ok, why = S5.orb1_eligible_files(c3, Lw, L3, "predguid", 1, by_block)
    check("ORB-1 block 2 waits until block 1 is complete", not ok and "waits" in why, f"({why})")
    # tuning round 2: drawn only once round 1 is scored, around its best, inside the ranges, the same on a rerun
    full = S5.StudyConfigV5()
    Lt = S5.LayoutV5(tmp / "tune_plan").make()
    n1, n2 = full.tune_n, full.tune_refine_n
    before = S5.guidance_tune_plan(full, Lt, "fnpeg", S5.FNPEG_SPACE)
    for j, _rnd, c in before:
        Lt.write_json(Lt.tuning / "fnpeg" / f"config_{j:02d}.json",
                      {"method": "fnpeg", "config_index": j, "config": dict(c), "sampled": c,
                       "score": {"success_5km": 0.9 if j == 7 else 0.5, "mean_log_miss": 0.0, "mean_peak_g": 4.0}})
    after = S5.guidance_tune_plan(full, Lt, "fnpeg", S5.FNPEG_SPACE)
    again = S5.guidance_tune_plan(full, Lt, "fnpeg", S5.FNPEG_SPACE)
    doc = json.loads((Lt.tuning / "fnpeg" / "refine.json").read_text())
    best = before[7][2]
    near = all(abs(c["deadband_c0_deg"] - best["deadband_c0_deg"]) <= 0.25 * 2.5 + 1e-9 and 0.5 <= c["deadband_c0_deg"] <= 3.0
               and c["dt_s"] in (1.0, 2.0) for _, r, c in after if r == 2)
    check("tuning round 2: drawn once round 1 is scored, around its best, recorded, identical on a rerun",
          len(before) == n1 and len(after) == n1 + n2 and after == again and doc["center_config_index"] == 7 and near,
          f"({len(before)} then {len(after)} configurations; center = configuration {doc['center_config_index']})")
    st = full.stages()
    gj = S5.gate_jobs(full, Lt)
    check("Tier 1 flies the 180° and oracle arms; gates V1-V4 fly DEV48 only; RL tuning uses 12 configurations",
          st["explore"] and not st["ORB1_more"] and {j["suite"] for v in gj.values() for j in v} == {"DEV48"}
          and len(gj["V1"]) == 12 and all(v == 12 for v in full.rl_tune_n.values()),
          f"({len(gj['V1'])} V1 flights on {sorted({j['suite'] for j in gj['V1']})})")


# ============================================================================
# 10 (guard). v5 refuses to write under its inputs
# ============================================================================
def t10_guard(tmp: Path, ro: dict) -> None:
    section("10. v5 refuses to write under the v3, v4, data and legacy RL folders")
    from jepa import study_v5 as S5
    L = S5.LayoutV5(tmp / "guard", read_only=list(ro.values())).make()
    for name, root in ro.items():
        p = Path(root) / "x_v5_probe.json"
        try:
            L.write_json(p, {"x": 1})
            ok = False
        except S5.S4.ReadOnlyViolation:
            ok = True
        check(f"write into the {name} folder refused", ok and not p.exists())
    for what, p in {"outside every folder": tmp / "elsewhere.json",
                    "a path escaping the v5 folder": L.root / ".." / "escape.json",
                    "a Windows extended-length path inside v3": "\\\\?\\" + str(Path(ro["v3"]).resolve() / "y.json")}.items():
        try:
            L.write_json(p, {"x": 1})
            ok = False
        except S5.S4.ReadOnlyViolation:
            ok = True
        check(f"write to {what} refused", ok)
    rec = {"result": {"flight_id": "../../../escape"}, "telemetry": {"t_s": np.zeros(1, np.float32)}}
    try:
        L.write_flight(rec)
        ok = False
    except S5.S4.ReadOnlyViolation:
        ok = True
    check("a flight id escaping the results folder refused", ok)
    p = L.root / "inside.json"
    L.write_json(p, {"x": 1})
    check("an extended-length path inside the v5 folder is accepted (Windows guard fix)",
          os.name != "nt" or L.guard("\\\\?\\" + str(p.resolve())) is not None)


# ============================================================================
# 9, 10, 13, 14, 15: the notebook
# ============================================================================
def nb_env(root: Path, fx: dict, fail: str = None) -> dict:
    env = {**os.environ, "JEPA_V5_MODE": "smoke", "JEPA_V5_ROOT": str(root), "JEPA_V3_ROOT": str(fx["v3"]),
           "JEPA_V4_ROOT": str(fx["v4"]), "JEPA_WIDE_DATA": str(fx["data"]), "JEPA_V5_LEGACY": str(fx["legacy"]),
           "PYTHONIOENCODING": "utf-8"}
    env.pop("JEPA_V5_FAIL_AT", None)
    if fail:
        env["JEPA_V5_FAIL_AT"] = fail
    return env


def _newest(root: Path) -> float:
    m = 0.0
    for dp, _, fs in os.walk(root):
        for f in fs:
            try:
                m = max(m, os.stat(os.path.join(dp, f)).st_mtime)
            except OSError:
                pass
    return m


def run_nb(root: Path, fx: dict, tag: str, fail: str = None, stall_min: float = 25.0) -> tuple:
    """Execute the notebook with nbconvert (output and errors to nb_runs/<tag>.out/.err). A run that writes nothing
    under its study folder for stall_min minutes is killed and reported as stalled (a pause of the whole machine,
    such as a laptop's sleep, restarts that clock). Run the smoke test from a console that stays open: when the
    console of this process is closed, the notebook's worker processes can no longer start and the run hangs."""
    out = root.parent / "nb_runs"
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    fo, fe = out / f"{tag}.out", out / f"{tag}.err"
    with open(fo, "w", encoding="utf-8") as so, open(fe, "w", encoding="utf-8") as se:
        p = subprocess.Popen([sys.executable, "-m", "nbconvert", "--to", "notebook", "--execute",
                              "--ExecutePreprocessor.timeout=-1", "--ExecutePreprocessor.startup_timeout=180",
                              "--output-dir", str(out), "--output", f"{tag}.ipynb", "ReEntryAI_JEPA_v5.ipynb"],
                             env=nb_env(root, fx, fail), stdout=so, stderr=se)
        stalled = False
        last_poll, t_ref = time.time(), t0
        while p.poll() is None:
            try:
                p.wait(timeout=20)
            except subprocess.TimeoutExpired:
                now = time.time()
                if now - last_poll > 300:           # the machine slept: that pause is not a stall
                    t_ref = now
                last_poll = now
                if now - max(t_ref, _newest(root)) > stall_min * 60:
                    stalled = True
                    subprocess.run(["taskkill", "/F", "/T", "/PID", str(p.pid)] if os.name == "nt"
                                   else ["kill", "-9", str(p.pid)], capture_output=True)
                    p.wait()
    rd = lambda f: f.read_text(encoding="utf-8", errors="replace")
    r = subprocess.CompletedProcess(p.args, p.returncode, rd(fo), rd(fe) + (" STALLED" if stalled else ""))
    if stalled:
        NOTES.setdefault("stalled_runs", []).append(tag)
    return r, time.time() - t0, out / f"{tag}.ipynb"


def _nb_text(path: Path) -> str:
    try:
        nb = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return ""
    return "".join("".join(o.get("text", "")) for c in nb["cells"] if c["cell_type"] == "code"
                   for o in c.get("outputs", []))


def t10_notebook(work: Path, fx: dict, keep: bool = False) -> Path:
    """keep: the smoke study folder exists from an earlier run 1 (same code); only the reruns are checked."""
    section("10. the notebook in smoke mode (first run = the whole smoke study)")
    subprocess.run([sys.executable, "jepa/build_notebook_v5.py"], check=True, capture_output=True)
    root = work / "v5_smoke"
    m_in = _mtimes([fx["v3"], fx["v4"], fx["legacy"]])
    if keep and (root / "00_study" / "study_config.json").exists():
        r, t1, _ = run_nb(root, fx, "run1b")
        if r.returncode != 0:
            print(r.stdout[-3000:], r.stderr[-8000:])
        check("notebook rerun on the kept smoke study ends without error", r.returncode == 0, f"({t1:.0f}s)")
    else:
        if root.exists():
            shutil.rmtree(root)
        r, t1, out = run_nb(root, fx, "run1")
        if r.returncode != 0:
            print(r.stdout[-3000:], r.stderr[-8000:])
        check("notebook run 1 ends without error", r.returncode == 0, f"({t1 / 60:.1f} min)")
        NOTES["nb_run1_min"] = round(t1 / 60, 1)
    r, t2, out2 = run_nb(root, fx, "run2")
    check("notebook run 2 ends without error", r.returncode == 0, f"({t2:.0f}s)")
    check("run 2 finds everything done in under a minute", t2 < 60.0, f"({t2:.0f}s)")
    txt = _nb_text(out2)
    flown = re.findall(r"\[[^\]]+\] \d+/\d+ done, (\d+) to fly", txt)
    check("run 2: no flight flown, no model trained", all(int(x) == 0 for x in flown) and "ms/step" not in txt,
          f"({len(flown)} queues checked)")
    check("v3, v4 and legacy folders untouched by the whole smoke study", m_in == _mtimes([fx["v3"], fx["v4"],
                                                                                         fx["legacy"]]))
    return root


def _mtimes(roots) -> dict:
    out = {}
    for r in roots:
        for dp, _, fs in os.walk(r):
            for f in fs:
                p = os.path.join(dp, f)
                out[p] = os.stat(p).st_mtime_ns
    return out


def t9_analysis(root: Path, fx: dict) -> None:
    section("9. analysis and report on the smoke data")
    from jepa import analyze_v5 as A5
    from jepa import report_v5 as R5
    from jepa import study_v4 as S4
    from jepa import study_v5 as S5
    cfg = S5.StudyConfigV5.smoke()
    L = S5.LayoutV5(root, read_only=[fx["v3"], fx["v4"], fx["legacy"], fx["data"]])
    need_t = [f"T{i:02d}" for i in range(1, 15)]
    need_f = [f"F{i:02d}" for i in range(1, 11)]
    have_t = {p.stem.split("_")[0] for p in L.tables.glob("T*.csv")}
    have_f = {p.stem.split("_")[0] for p in L.figures.glob("F*.png")}
    have_pdf = {p.stem.split("_")[0] for p in L.figures.glob("F*.pdf")}
    check("every table T01-T14", set(need_t) <= have_t, f"(missing {sorted(set(need_t) - have_t)})")
    check("every figure F01-F10 as PNG and PDF", set(need_f) <= have_f and set(need_f) <= have_pdf,
          f"(missing {sorted((set(need_f) - have_f) | (set(need_f) - have_pdf))})")
    caps = json.loads((L.analysis / "captions.json").read_text())
    check("every table and figure has a caption and a data source",
          all(any(k.startswith(x) for k in caps) for x in need_t + need_f) and all(v.get("data") for v in caps.values()))
    t05 = pd.read_csv(next(L.tables.glob("T05_*.csv")))
    ids = ["G1", "G2", "G3", "G4", "G5", "G6", "G7a", "G7b", "G7c", "G7d", "G7e", "A1", "A2", "A-shift", "J-mlp_direct", "J-gru_direct", "J-gru_pi",
           "J-r0_nojepa", "J-equiv", "J-mech", "R-dagger_fnpeg", "R-bcppo_fnpeg", "R-ppo_hist", "R-imit", "R-gap",
           "Q3a", "Q3b", "Q3a-sens", "Q3b-sens", "D1", "D2", "D-shift"]
    check("criteria table complete (G1-G7, J-each/equiv/mech, R-each/imit/gap, Q1-Q3, D1-D2)",
          set(ids) <= set(t05.id) and t05.id.str.startswith("Q1-").any() and t05.id.str.startswith("Q2-").any(),
          f"({len(t05)} rows; missing {sorted(set(ids) - set(t05.id))}; not run: "
          f"{list(t05[t05.measured == 'not run'].id)})")
    NOTES["criteria_not_run"] = list(t05[t05.measured == "not run"].id)
    pend = list(t05[t05.status == "pending"].id) if "status" in t05 else ["no status column"]
    check("every criterion has a status; on the complete smoke study only Q1 can be pending (too few cases)",
          "status" in t05 and set(t05.status) <= {"met", "not met", "pending", "estimate", "not run"}
          and all(str(x).startswith("Q1-") for x in pend), f"(pending: {pend})")
    t04 = pd.read_csv(next(L.tables.glob("T04_*.csv")))
    pooled = t04[t04.atmosphere == "all"]
    check("paired-flight counts and clustered-case intervals in T04; Q3 sensitivity table T08b",
          {"n_pairs", "n_expected", "complete", "ci_lo_clustered", "ci_hi_clustered"} <= set(t04.columns)
          and bool(pooled.complete.all()) and pooled.ci_lo_clustered.notna().any() and any(L.tables.glob("T08b_*.csv")),
          f"({len(pooled)} pooled comparisons, {int(pooled.complete.sum())} complete)")
    t13 = pd.read_csv(next(L.tables.glob("T13_*.csv")))
    t13b = pd.read_csv(next(L.tables.glob("T13b_*.csv")))
    d1 = t05.set_index("id").loc[["D1", "D2"]]
    check("DS1: T13 has D1 and D2 complete, T13b every DS1 arm against A256, D1/D2 decided in their own family",
          {"D1", "D2"} <= set(t13.id) and bool(t13.set_index("id").loc[["D1", "D2"], "complete"].all())
          and set(cfg.ds1_arms) <= set(t13b.arm) and set(d1.family) == {"DS1"}
          and set(d1.status) <= {"met", "not met"},
          f"(T13 {list(t13.id)}; T13b {list(t13b.arm)}; D1/D2 {list(d1.status)})")
    t14 = pd.read_csv(next(L.tables.glob("T14_*.csv")))
    t14b = pd.read_csv(next(L.tables.glob("T14b_*.csv")))
    a1 = t05.set_index("id").loc[["A1", "A2"]]
    check("AS1: T14 has A1 and A2 complete, T14b every AS1 arm against A256's nominal air, A1/A2 decided in their own "
          "family", {"A1", "A2"} <= set(t14.id) and bool(t14.set_index("id").loc[["A1", "A2"], "complete"].all())
          and set(cfg.as1_arms) <= set(t14b.arm) and set(a1.family) == {"AS1"} and set(a1.status) <= {"met", "not met"},
          f"(T14 {list(t14.id)}; T14b {list(t14b.arm)}; A1/A2 {list(a1.status)})")
    fl = pd.read_csv(L.flights / "flights.csv", low_memory=False)
    v5f = fl[fl.source == "v5"] if "source" in fl else fl
    iv_cols = {"use_intervals", "shield_interventions", "shield_steps", "interval_violations"} <= set(fl.columns)
    legacy = v5f.arm.astype(str).str.startswith("rl_legacy_")
    check("every v5 flight record carries the interval fields: shield never stepped, no interventions, intervals off "
          "except legacy policies with the interval observation",
          iv_cols and len(v5f) > 0 and bool(v5f.shield_interventions.notna().all())
          and int(v5f.shield_interventions.sum()) == 0 and int(v5f.shield_steps.sum()) == 0
          and not bool(v5f.use_intervals[~legacy].astype(str).str.lower().isin(["true", "1", "1.0"]).any()),
          f"({len(v5f)} v5 flights, columns present: {iv_cols})")
    t09 = pd.read_csv(next(L.tables.glob("T09_*.csv")))
    check("compute table times the classical controller as flown and with its interval heat screen off",
          {"classical", "classical_noscreen"} <= set(t09.arm), f"({sorted(set(t09.arm))})")
    L3 = S4.v3_layout(fx["v3"])
    L4 = S5.v4_layout(fx["v4"])
    try:
        R = A5.analyze(cfg, L, L3, L4, log=lambda *a: None, force=True)
        b, c = R["boot"], R["boot_clustered"]
        check("bootstrap cache reused on a forced rerun (no recomputation)", b.computed == 0 and b.reused > 0
              and c.computed == 0 and c.reused > 0, f"({b.reused} + {c.reused} clustered reused)")
        R = A5.analyze(cfg, L, L3, L4, log=lambda *a: None)
        check("an unforced rerun with unchanged inputs skips the analysis", bool(R.get("skipped")))
    except Exception as e:
        check("analysis reruns", False, f"({type(e).__name__}: {e})")
    pdf = L.report / "SpaceJEPA_v5_report.pdf"
    check("report written by the notebook", pdf.exists() and pdf.stat().st_size > 50_000,
          f"({pdf.stat().st_size // 1000 if pdf.exists() else 0} kB)")
    from jepa import make_report_v2 as R2
    old = R2.FONT_DIR
    R2.FONT_DIR = str(root / "no_system_fonts_here")
    try:
        p = R5.build_report(L, out_pdf=str(L.report / "smoke_fallback_fonts.pdf"))
        ok = Path(p).stat().st_size > 50_000
        Path(p).unlink()
    except Exception as e:
        ok = False
        print(e)
    finally:
        R2.FONT_DIR = old
    check("report builds with matplotlib's DejaVu (no system fonts)", ok)
    # an incomplete arm gets no verdict: remove one FNPEG flight from a copy and analyze it
    c = root.parent / "v5_smoke_pending"
    if c.exists():
        shutil.rmtree(c)
    shutil.copytree(root, c)
    try:
        Lc = S5.LayoutV5(c, read_only=[fx["v3"], fx["v4"], fx["legacy"], fx["data"]])
        sorted((Lc.flights / "results" / "A256" / "fnpeg").rglob("*.json"))[0].unlink()
        t = A5.analyze(cfg, Lc, L3, L4, log=lambda *a: None, force=True)["T05"].set_index("id")
        pdf2 = Path(R5.build_report(Lc, out_pdf=str(Lc.report / "pending.pdf")))
        check("a criterion with a missing flight is pending (no verdict, numbers marked provisional); others decided",
              t.status["G1"] == "pending" and pd.isna(t.met["G1"]) and "PROVISIONAL" in t.measured["G1"]
              and t.n_pairs["G1"] < t.n_expected["G1"] and t.status["G2"] != "pending" and pdf2.stat().st_size > 50_000,
              f"(G1: {int(t.n_pairs['G1'])} of {int(t.n_expected['G1'])} paired flights, status {t.status['G1']}; "
              f"G2: {t.status['G2']})")
    except Exception as e:
        check("a criterion with a missing flight is pending", False, f"({type(e).__name__}: {e})")
    finally:
        shutil.rmtree(c, ignore_errors=True)


def _flights(root: Path) -> pd.DataFrame:
    rows = []
    for p in sorted((root / "04_flights" / "results").rglob("*.json")):
        d = json.loads(p.read_text())
        d["_path"] = str(p.relative_to(root))
        d["_mtime"] = p.stat().st_mtime_ns
        rows.append(d)
    return pd.DataFrame(rows)


VOLATILE = ("wall_s", "plan_ms_mean", "plan_ms_max", "plan_ms_median", "plan_ms_p99", "git_commit", "_mtime")


def _compare_flights(a: pd.DataFrame, b: pd.DataFrame) -> tuple:
    """(number compared, identical, list of differing flight ids) on every non-volatile result field."""
    a, b = a.set_index("_path"), b.set_index("_path")
    common = sorted(set(a.index) & set(b.index))
    cols = [c for c in a.columns if c in b.columns and c not in VOLATILE]
    diff = []
    for p in common:
        x, y = a.loc[p, cols], b.loc[p, cols]
        same = all((u == v) or (isinstance(u, float) and isinstance(v, float) and math.isnan(u) and math.isnan(v))
                   or (pd.isna(u) and pd.isna(v)) for u, v in zip(x, y))
        if not same:
            diff.append(p)
    return len(common), len(common) - len(diff), diff, sorted(set(a.index) ^ set(b.index))


def _weights_equal(ra: Path, rb: Path) -> tuple:
    n, bad = 0, []
    for p in sorted(list((ra / "02_models").rglob("ckpt.pt")) + list((ra / "02_policies").rglob("policy.pt"))):
        q = rb / p.relative_to(ra)
        if "03_tuning" in str(p):
            continue
        n += 1
        if not q.exists():
            bad.append(str(p.relative_to(ra)) + " (missing)")
            continue
        x, y = torch.load(p, weights_only=False), torch.load(q, weights_only=False)

        def flat(d, pre=""):
            out = {}
            if isinstance(d, dict):
                for k, v in d.items():
                    out.update(flat(v, f"{pre}/{k}"))
            elif isinstance(d, (list, tuple)):
                for i, v in enumerate(d):
                    out.update(flat(v, f"{pre}/{i}"))
            elif torch.is_tensor(d):
                out[pre] = d
            return out
        fx, fy = flat(x), flat(y)
        if set(fx) != set(fy) or not all(torch.equal(fx[k], fy[k]) for k in fx):
            bad.append(str(p.relative_to(ra)))
    return n, bad


def t13_repair(root: Path, work: Path, fx: dict) -> None:
    section("13. interrupted writes: Setup repairs each case and the rerun ends with the same outputs")
    from jepa import study_v5 as S5
    c = work / "v5_smoke_repair"
    if c.exists():
        shutil.rmtree(c)
    shutil.copytree(root, c)
    L = S5.LayoutV5(c)
    f0 = _flights(c)
    (L.flights / "results" / ".stale.json.123.tmp").write_bytes(b"half")
    victim = sorted((L.flights / "results" / "A256" / "predguid").glob("*.json"))[0]
    victim.write_text('{"flight_id": "A256/predguid/nominal_case00", "miss_km": ')
    os.utime(victim, None)
    e0csv = L.checks / "E0.csv"
    e0csv.unlink()
    rep = S5.repair(L, log=lambda *a: None)
    check("repair: stale .tmp deleted, corrupt flight JSON deleted, .done with a missing artifact removed",
          rep["tmp_deleted"] >= 1 and rep["json_deleted"] == 1 and rep["done_removed"] >= 1
          and not victim.exists() and not (L.checks / "E0.done").exists(), f"({rep})")
    r, t, _ = run_nb(c, fx, "repair_run")
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-5000:])
    f1 = _flights(c)
    n, same, diff, extra = _compare_flights(f0, f1)
    check("the rerun re-flies the corrupt flight, recomputes E0, and every flight equals the untouched run's",
          r.returncode == 0 and victim.exists() and e0csv.exists() and (L.checks / "E0.done").exists()
          and same == n and not extra, f"({same}/{n} identical, {t:.0f}s; differing {diff[:3]})")
    shutil.rmtree(c, ignore_errors=True)


# in the notebook's stage order, so each lands in a stage the earlier runs have not finished: legacy flights,
# RL tuning (DAgger), the A256 MPPI arms, E3 training, E3 flights, E4 training, E4 flights
FAILS = ["cpu_flights:3", "rl_train:1", "mppi_flights:2", "e3_train:4", "mppi_flights:1", "rl_train:1",
         "cpu_flights:2"]


def t14_disconnect(root: Path, work: Path, fx: dict) -> None:
    section("14. simulated disconnects in the notebook's long stages, then an unchanged rerun")
    b = work / "v5_smoke_b"
    if b.exists():
        shutil.rmtree(b)
    died = []
    for k, f in enumerate(FAILS):
        r, t, _ = run_nb(b, fx, f"fail{k}", fail=f)
        hit = "simulated disconnect" in (r.stdout + r.stderr) or "DeadKernel" in r.stderr or "died" in r.stderr
        died.append((f, r.returncode != 0 and hit, round(t)))
        if r.returncode != 0 and not hit:
            print(r.stdout[-2000:], r.stderr[-5000:])
    check("each forced failure stopped the run (simulated disconnect)", all(d for _, d, _ in died), f"({died})")
    before = _flights(b)
    r, t, _ = run_nb(b, fx, "after_fails")
    _wait_reference(root)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-5000:])
    check("the unchanged rerun completes", r.returncode == 0, f"({t / 60:.1f} min)")
    after = _flights(b)
    m0 = before.set_index("_path")._mtime
    m1 = after.set_index("_path")._mtime
    kept = all(m1.get(p) == m for p, m in m0.items())
    check("finished flights keep their modification times", kept, f"({len(m0)} flights before the rerun)")
    ref = _flights(root)
    n, same, diff, extra = _compare_flights(ref, after)
    check("every flight equals the uninterrupted run's (bit-exact on CPU)", same == n and not extra,
          f"({same}/{n} identical; differing {diff[:4]}; only in one run {extra[:4]})")
    NOTES["t14_diff"] = diff
    nw, bad = _weights_equal(root, b)
    check("every trained model and policy equals the uninterrupted run's (bit-exact on CPU)", not bad,
          f"({nw - len(bad)}/{nw} identical; differing {bad[:4]})")


T_START = time.time()


def _wait_reference(root: Path, hours: float = 8.0) -> None:
    """When check 14 runs in parallel with check 10 (another process), wait for that process's run 2."""
    done = root.parent / "nb_runs" / "run2.ipynb"
    t0 = time.time()
    while time.time() - t0 < hours * 3600:
        if done.exists() and done.stat().st_mtime > T_START - 4 * 3600:
            return
        time.sleep(30)


def t15_outputs(root: Path) -> None:
    section("15. outputs in the smoke study folder")
    from jepa import study_v5 as S5
    L = S5.LayoutV5(root)
    cfg = S5.StudyConfigV5.smoke()
    need = [L.study / "status.json", L.study / "status_history.csv", L.study / "README.md",
            L.study / "study_config.json", L.analysis / "captions.json", L.report / "SpaceJEPA_v5_report.pdf",
            L.study / "legacy_rl_inventory.csv", L.checks / "compute_per_decision.csv"]
    missing = [str(p.relative_to(root)) for p in need if not p.exists()]
    check("status.json, status_history.csv, README, frozen config, captions, report, inventory, timing", not missing,
          f"(missing {missing})")
    hist = pd.read_csv(L.study / "status_history.csv") if (L.study / "status_history.csv").exists() else []
    check("status history appended by every run", len(hist) >= 4, f"({len(hist)} rows)")
    logs = {p.stem for p in L.logs.glob("*.log")}
    want = {"setup", "benchmark", "gates", "tuning", "A256", "ORB1", "E3_train", "E3_fly", "E4_train", "E4_fly",
            "explore", "DS1", "AS1", "open_loop", "analysis"}
    check("a log per stage in 00_study/logs/", want <= logs, f"(missing {sorted(want - logs)})")
    curves = [L.model_dir(a, s) / "F_curves.png" for a in cfg.e3_arms for s in cfg.seeds] + \
        [L.policy_dir(a, s) / "F_curves.png" for a in S5.RL_TIER1 for s in cfg.rl_seeds] + \
        [L.policy_dir(a, s) / "learning_curve.csv" for a in S5.RL_TIER1 for s in cfg.rl_seeds]
    miss = [str(p.relative_to(root)) for p in curves if not p.exists()]
    check("curves of every trained model and policy", not miss, f"(missing {miss[:5]})")
    res = L.flights / "results"
    g_dev = sorted(p.name for p in (res / "DEV48").glob("gate_*")) if (res / "DEV48").exists() else []
    g_test = sorted(p.name for p in (res / "A256").glob("gate_*")) if (res / "A256").exists() else []
    check("gate V1-V4 flights are on the dev suite, none on a test suite", bool(g_dev) and not g_test,
          f"(DEV48: {g_dev}; A256: {g_test})")
    tn = {}
    for m in ("fnpeg", "predguid"):
        d = L.tuning / m
        ch = json.loads((d / "chosen.json").read_text()) if (d / "chosen.json").exists() else {}
        cs = pd.read_csv(d / "configs.csv") if (d / "configs.csv").exists() else pd.DataFrame(columns=["round"])
        tn[m] = ((d / "refine.json").exists(), ch.get("n_configs"), sorted(cs["round"].unique().tolist()))
    want_n = cfg.tune_n + cfg.tune_refine_n
    check("tuning: round 2 recorded (refine.json) and both rounds scored before the choice",
          all(v == (True, want_n, [1, 2]) for v in tn.values()), f"({tn})")
    ex = sorted(p.name for p in (res / "A256").glob("*") if p.name in S5.EXPLORE_ARMS) if (res / "A256").exists() else []
    check("the 180°, oracle and no-relief arms flew in Tier 1", ex == sorted(S5.EXPLORE_ARMS), f"({ex})")
    tabs = len(list(L.tables.glob("T*.csv")))
    figs = len(list(L.figures.glob("F*.png"))), len(list(L.figures.glob("F*.pdf")))
    check("tables and figures (PNG and PDF) on disk", tabs >= 12 and min(figs) >= 10, f"({tabs} tables, {figs} figures)")


# ============================================================================
# 11. v3 / v4 unchanged
# ============================================================================
V34 = ["jepa/study_v3.py", "jepa/fly_v3.py", "jepa/analyze_v3.py", "jepa/testset_v3.py", "jepa/build_notebooks_v3.py",
       "jepa/study_v4.py", "jepa/fly_v4.py", "jepa/physics_v4.py", "jepa/analyze_v4.py", "jepa/report_v4.py",
       "jepa/ablate_v4.py", "jepa/smoke_v4.py", "jepa/build_notebook_v4.py", "jepa/models_v2.py", "jepa/models.py",
       "jepa/train_v2.py", "jepa/mppi_v2.py", "jepa/data.py", "rl_env.py", "math_3d.py", "ReEntryAI_JEPA_v4.ipynb",
       "jepa/v3_plan.md", "jepa/v4_plan.md", "ReEntryAI_JEPA_v3_0_testset.ipynb", "ReEntryAI_JEPA_v3_1_train.ipynb",
       "ReEntryAI_JEPA_v3_2_fly.ipynb", "ReEntryAI_JEPA_v3_3_analyze.ipynb", "ReEntryAI_JEPA_v3_resultspdf.ipynb"]


# The one change v5 makes to the simulator (approved by the user 2026-11-15, v5_plan.md section 3.4): an optional
# altitude-dependent density multiplier, off unless the parameters carry a density_profile.
HOOK_FILES = {"math_3d.py": 2, "jepa/data.py": 1}     # file -> number of added hunks


def _hook_only(base: str) -> Tuple[bool, str]:
    """The diff of the hook files against base adds lines only, in exactly the expected hunks, each about the
    density profile."""
    r = subprocess.run(["git", "diff", "-U0", base, "--", *HOOK_FILES], capture_output=True, text=True)
    if r.returncode:
        return False, r.stderr.strip()
    hunks: Dict[str, List[List[str]]] = {f: [] for f in HOOK_FILES}
    removed, cur = [], None
    for ln in r.stdout.splitlines():
        if ln.startswith("+++ "):
            cur = ln[6:] if ln.startswith("+++ b/") else None
            cur = next((f for f in HOOK_FILES if cur and cur.endswith(f)), None)
        elif ln.startswith("@@") and cur:
            hunks[cur].append([])
        elif ln.startswith("-") and not ln.startswith("--- "):
            removed.append(ln)
        elif ln.startswith("+") and cur and hunks[cur]:
            hunks[cur][-1].append(ln[1:])
    n = {f: len(h) for f, h in hunks.items()}
    about = all("density_profile" in "\n".join(h) or "aero_schedule" in "\n".join(h)
                for hs in hunks.values() for h in hs)
    return (not removed and n == HOOK_FILES and about), f"(removed lines {len(removed)}, added hunks {n})"


def t11_unchanged(tmp: Path) -> None:
    section("11. v3 and v4 code, plans and notebooks unchanged")
    base = "2565808"                                   # the v4 study commit (v5 work starts after it)
    others = [p for p in V34 if p not in HOOK_FILES]
    r = subprocess.run(["git", "diff", "--name-only", base, "--", *others], capture_output=True, text=True)
    check("v3/v4 modules, plans, simulator and the v4 notebook identical to the v4 commit (except the density hook)",
          r.returncode == 0 and not r.stdout.strip(), f"(changed: {r.stdout.split()})")
    ok, why = _hook_only(base)
    check("math_3d.py and data.py differ from the v4 commit only by the default-off density hook (lines added, none "
          "removed)", ok, why)
    before = (SIM_DIR / "ReEntryAI_JEPA_v4.ipynb").read_bytes()
    subprocess.run([sys.executable, "jepa/build_notebook_v4.py"], check=True, capture_output=True)
    after = (SIM_DIR / "ReEntryAI_JEPA_v4.ipynb").read_bytes()
    check("the v4 notebook builder still writes the committed notebook byte for byte", before == after)
    if before != after:
        (SIM_DIR / "ReEntryAI_JEPA_v4.ipynb").write_bytes(before)
    d = tmp / "nb3" / "jepa"
    d.mkdir(parents=True, exist_ok=True)
    shutil.copy(SIM_DIR / "jepa" / "build_notebooks_v3.py", d / "build_notebooks_v3.py")
    subprocess.run([sys.executable, str(d / "build_notebooks_v3.py")], cwd=d.parent, check=True, capture_output=True)
    # cell types and sources (the committed v3 analysis notebook also carries the outputs of an executed run)
    cells = lambda p: [(c["cell_type"], "".join(c["source"]))
                       for c in json.loads(p.read_text(encoding="utf-8"))["cells"]]
    names = sorted(p.name for p in d.parent.glob("ReEntryAI_JEPA_v3_*.ipynb"))
    same = [cells(d.parent / n) == cells(SIM_DIR / n) for n in names]
    # a notebook the user last saved from Colab ("Created using Colab") carries their edits; git diff covers it
    colab = [n for n, s in zip(names, same) if not s and "Created using Colab" in subprocess.run(
        ["git", "log", "-1", "--format=%s", "--", n], capture_output=True, text=True).stdout]
    bad = [n for n, s in zip(names, same) if not s and n not in colab]
    check(f"the v3 builder reproduces the cells of the {len(names)} committed v3 notebooks", len(names) == 5 and not bad,
          f"(differ: {bad}; saved from Colab by the user, checked by git diff only: {colab})")


# ============================================================================
def throughput(root: Path) -> pd.DataFrame:
    """Measured smoke throughput per arm type (for v5_plan.md section 14)."""
    f = _flights(root)
    if not len(f):
        return pd.DataFrame()
    f["type"] = f.model_kind.fillna("classical")
    f = f[~f.stage.isin(["gates"])]
    g = f.groupby("type").agg(flights=("wall_s", "size"), wall_s_median=("wall_s", "median"),
                              plan_ms_median=("plan_ms_median", "median")).reset_index()
    return g.round(2)


# ============================================================================
# 16. RL on Kaggle (jepa/kaggle_rl.md)
# ============================================================================
def run_kaggle_nb(base: Path, live: Path, out: Path, inputs: Sequence[Path], tag: str,
                  max_runs: Optional[int] = None) -> tuple:
    nbdir = base / "nb"
    nbdir.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "JEPA_V5_MODE": "smoke", "JEPA_V5_KAGGLE_LIVE": str(live), "JEPA_V5_KAGGLE_OUT": str(out),
           "JEPA_V5_KAGGLE_INPUTS": os.pathsep.join(map(str, inputs)), "PYTHONIOENCODING": "utf-8"}
    env.pop("JEPA_V5_KAGGLE_MAX_RUNS", None)
    if max_runs is not None:
        env["JEPA_V5_KAGGLE_MAX_RUNS"] = str(max_runs)
    t0 = time.time()
    r = subprocess.run([sys.executable, "-m", "nbconvert", "--to", "notebook", "--execute",
                        "--ExecutePreprocessor.timeout=-1", "--ExecutePreprocessor.startup_timeout=180",
                        "--output-dir", str(nbdir), "--output", f"{tag}.ipynb", "ReEntryAI_JEPA_v5_kaggle_rl.ipynb"],
                       env=env, capture_output=True, text=True, encoding="utf-8", errors="replace")
    return r, time.time() - t0, _nb_text(nbdir / f"{tag}.ipynb")


def t16_kaggle(work: Path) -> None:
    section("16. RL on Kaggle: inputs out, two sessions (restore, resume), one platform, import into Drive")
    from jepa import kaggle_v5 as K
    from jepa import rl_v5 as RL
    from jepa import study_v5 as S5
    subprocess.run([sys.executable, "jepa/build_notebook_v5_kaggle.py"], check=True, capture_output=True)
    nb = (SIM_DIR / "ReEntryAI_JEPA_v5.ipynb").read_text(encoding="utf-8")
    check("Colab notebook: RL_ON_KAGGLE in full mode; RL tuning and training then only export and import",
          'RL_ON_KAGGLE = MODE == \\"full\\"' in nb and "K.export_inputs(L, log)" in nb and "K.import_rl(L, K_FROM, log)"
          in nb and "if RL_ON_KAGGLE:" in nb)
    src = work / "v5_smoke"
    if not (src / "03_tuning" / "fnpeg" / "chosen.json").exists():
        check("Kaggle check needs check 10's smoke study (FNPEG tuned)", False, f"({src})")
        return
    base = work / "kaggle_smoke"
    shutil.rmtree(base, ignore_errors=True)
    drive = base / "drive"                                    # a Drive study folder holding only what Colab has
    Ld = S5.LayoutV5(drive).make()
    for rel in K.INPUT_FILES:
        if (src / rel).exists():
            Ld.atomic_bytes(drive / rel, (src / rel).read_bytes())
    inp = K.export_inputs(Ld)
    ds = base / "input1" / "reentryai-v5-rl-inputs"
    ds.mkdir(parents=True)
    shutil.copy2(inp, ds / K.INPUTS)
    cfg = S5.StudyConfigV5.smoke()
    n_tune = sum(int(cfg.rl_tune_n.get(a, 1)) for a in S5.RL_TIER1)
    n_runs = n_tune + len(S5.RL_TIER1) * len(cfg.rl_seeds)

    out1 = base / "out1"
    r1, t1, txt1 = run_kaggle_nb(base, base / "live1", out1, [base / "input1"], "k1", max_runs=2)
    m1 = json.loads((out1 / K.MANIFEST).read_text()) if (out1 / K.MANIFEST).exists() else {}
    if r1.returncode != 0:
        print(r1.stdout[-2000:], r1.stderr[-6000:])
    check("Kaggle session 1 ends without error, stops after its 2 runs and saves its state",
          r1.returncode == 0 and m1.get("seq") == 1 and len(m1.get("runs_done", [])) == 2 and "not started" in txt1,
          f"({t1:.0f}s; runs done {m1.get('runs_done')})")

    in2 = base / "input2"                                     # session 2: a new machine, session 1's output attached
    shutil.copytree(base / "input1", in2)
    shutil.copytree(out1, in2 / "kaggle-rl-session-1")
    out2 = base / "out2"
    r2, t2, txt2 = run_kaggle_nb(base, base / "live2", out2, [in2], "k2")
    m2 = json.loads((out2 / K.MANIFEST).read_text()) if (out2 / K.MANIFEST).exists() else {}
    if r2.returncode != 0:
        print(r2.stdout[-2000:], r2.stderr[-6000:])
    check("Kaggle session 2 restores session 1, finishes every RL tuning and training run, prints RL COMPLETE",
          r2.returncode == 0 and m2.get("seq") == 2 and "restored session 1" in txt2 and "RL COMPLETE" in txt2
          and len(m2.get("runs_done", [])) == n_runs, f"({t2:.0f}s; {len(m2.get('runs_done', []))}/{n_runs} runs)")

    L2 = S5.LayoutV5(base / "live2")
    hw = json.loads((base / "live2" / K.PLATFORM).read_text())
    refused = []
    try:
        K.check_platform(L2, {**hw, "gpu": "Tesla P100-PCIE-16GB"})
    except RuntimeError:
        refused.append("session")
    run_d = base / "live1" / "03_tuning" / "rl_dagger_fnpeg" / "config_00"
    if (run_d / "hardware.json").exists():
        h = json.loads((run_d / "hardware.json").read_text())
        (run_d / "hardware.json").write_text(json.dumps({**h, "gpu": "Tesla P100-PCIE-16GB"}))
        (run_d / "policy.pt").unlink(missing_ok=True)
        try:
            RL.check_hardware(S5.LayoutV5(base / "live1"), run_d, "cpu")
        except RuntimeError:
            refused.append("run")
    check("a session or a run on another GPU type is refused", refused == ["session", "run"], f"({refused})")

    fk = drive / "00_study" / "kaggle" / "from_kaggle"
    fk.mkdir(parents=True)
    for f in (K.STATE_TAR, K.MANIFEST):
        shutil.copy2(out2 / f, fk / f)
    ok_imp, detail = False, ""
    try:
        res = K.import_rl(Ld, fk, log=lambda *a: None)
        pols = all((Ld.policy_dir(a, s) / "policy.pt").exists() for a in S5.RL_TIER1 for s in cfg.rl_seeds)
        chosen = all((Ld.tuning / f"rl_{a}" / "chosen.json").exists() for a in S5.RL_TIER1)
        no_state = not any(p.name.startswith("train_state") for p in drive.rglob("*"))
        tune_fl = any((drive / "04_flights" / "results" / "DEV48").glob("tune_*"))
        res2 = K.import_rl(Ld, fk, log=lambda *a: None)
        ok_imp = pols and chosen and no_state and tune_fl and (drive / K.PLATFORM).exists() and res2["files_new"] == 0
        detail = f"({res}; again: {res2['files_new']} new)"
    except Exception as e:                                    # noqa: BLE001
        detail = f"({type(e).__name__}: {e})"
    check("import: every finished policy and RL tuning result on Drive, no train state, the platform recorded; "
          "importing again writes nothing", ok_imp, detail)

    p = Ld.policy_dir(S5.RL_TIER1[0], cfg.rl_seeds[0]) / "run.json"
    keep = p.read_bytes()
    p.write_bytes(keep + b" ")
    n_before = sum(1 for _ in drive.rglob("*"))
    shutil.rmtree(Ld.policy_dir(S5.RL_TIER1[1], cfg.rl_seeds[0]))
    conflict = False
    try:
        K.import_rl(Ld, fk, log=lambda *a: None)
    except ValueError:
        conflict = True
    untouched = not Ld.policy_dir(S5.RL_TIER1[1], cfg.rl_seeds[0]).exists()
    p.write_bytes(keep)
    bad_tar = fk / K.STATE_TAR
    b = bytearray(bad_tar.read_bytes())
    b[len(b) // 2] ^= 0xFF
    bad_tar.write_bytes(bytes(b))
    corrupt = False
    try:
        K.import_rl(Ld, fk, log=lambda *a: None)
    except ValueError:
        corrupt = True
    check("import refuses a file that differs from Drive's (writing nothing) and a corrupted download",
          conflict and untouched and corrupt, f"(conflict {conflict}, nothing written {untouched}, corrupt {corrupt}; "
                                              f"{n_before} files)")


def main(work: Path, workers: int, only=None, keep: bool = False) -> None:
    t_all = time.time()
    torch.set_num_threads(max(1, (os.cpu_count() or 2) // 2))
    fx = {"data": work / "wide_smoke", "v3": work / "v3_smoke", "v4": work / "v4_smoke", "legacy": work / "legacy_v5"}
    for k in ("data", "v3", "v4"):
        if not fx[k].exists():
            raise SystemExit(f"fixture {fx[k]} missing: run python jepa/smoke_v4.py first")
    make_legacy(fx["legacy"])
    from jepa import study_v4 as S4
    L3 = S4.v3_layout(fx["v3"])
    tmp = Path(tempfile.mkdtemp(prefix="smoke_v5_", dir=work))
    run = lambda k: only is None or k in only
    if run(11):
        t11_unchanged(tmp)
    if run(10):
        t10_guard(tmp, fx)
    if run(1):
        t1_suites(tmp, workers)
    if run(2):
        t2_fnpeg(tmp)
    if run(3):
        t3_predguid(tmp)
    if run(4):
        t4_e0(tmp, L3)
    if run(8):
        t8_queue(tmp, L3)
    if run(7):
        t7_legacy(tmp, fx["legacy"])
    if run(5):
        t5_e3(tmp, fx["data"])
    if run(6):
        t6_rl(tmp, fx["data"])
    root = work / "v5_smoke"
    if run(10) or run(9) or run(13) or run(14) or run(15):
        if run(10) or not root.exists():
            root = t10_notebook(work, fx, keep)
        if run(9):
            t9_analysis(root, fx)
        if run(15):
            t15_outputs(root)
        if run(13):
            t13_repair(root, work, fx)
        if run(14):
            t14_disconnect(root, work, fx)
        tp = throughput(root)
        if len(tp):
            print("\nmeasured smoke throughput per arm type (CPU):\n" + tp.to_string(index=False))
    if run(16):
        t16_kaggle(work)
    shutil.rmtree(tmp, ignore_errors=True)
    n_ok = sum(ok for _, ok, _ in RESULTS)
    print(f"\nnotes: {json.dumps(NOTES, default=str)}")
    print(f"\nv5 smoke test: {n_ok}/{len(RESULTS)} checks passed in {(time.time() - t_all) / 60:.1f} min. "
          f"Study folders in {work}")
    for name, ok, detail in RESULTS:
        if not ok:
            print(f"  FAILED: {name} {detail}")
    if n_ok != len(RESULTS):
        raise SystemExit(1)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
    except (AttributeError, ValueError):
        pass
    if len(sys.argv) > 2 and sys.argv[1] == "--child":
        child(*sys.argv[2:6])
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=os.environ.get("JEPA_SMOKE_WORK", str(Path(tempfile.gettempdir()) / "jepa_v4_smoke")))
    ap.add_argument("--workers", type=int, default=max(1, min(8, (os.cpu_count() or 2) - 1)))
    ap.add_argument("--only", default="", help="comma-separated check numbers (default: all)")
    ap.add_argument("--keep", action="store_true", help="check 10 on the existing smoke study (no fresh run 1)")
    a = ap.parse_args()
    main(Path(a.work), a.workers, {int(x) for x in a.only.split(",") if x} or None, a.keep)
