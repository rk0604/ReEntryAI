"""
jepa/smoke_v3.py

CPU smoke test of the v3 study code (jepa/v3_plan.md section 11.2). The
numbers it prints are meaningless (tiny models, a few flights); the checks
are what matter:

  1. test set: the footprint sweep matches real env flights, the entry
     placement round-trips, and the frozen CSV is reproduced from its Sobol
     indices
  2. the lockstep runner flies the same decisions as the unchanged
     mppi.run_episode for v1 (shipped checkpoint, 20 decisions, 3 cases);
     for single-member v2 it is identical to mppi_v2.run_lockstep over whole
     flights, its parallel rollout matches the step-by-step one, and its first
     decisions match mppi.run_episode
  3. common random numbers: every learned controller gets the same MPPI
     noise for the same (seed, case, atmosphere, rep), and different
     conditions get different noise
  4. telemetry: fields, row count, fuel -> jet-seconds identity, fuel by
     phase sums to the total, final fuel equals the simulator's
  5. the job queue skips finished flights and re-flies a deleted one
  6. resumable training: interrupt + resume, the convergence-extension
     rule, per-seed v1/v2 training on a tiny config
  7. statistics: Holm, McNemar, sign test, bootstrap on known inputs

    python jepa/smoke_v3.py [--wide <wide-envelope shard dir or pickle>]
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))
os.chdir(SIM_DIR)

import numpy as np
import pandas as pd
import torch

from jepa import analyze_v3 as A
from jepa import data as D
from jepa import fly_v3 as F
from jepa import mppi as P
from jepa import study_v3 as S
from jepa import testset_v3 as TS

V1_CKPT = "jepa_runs/full/jepa_ckpt.pt"


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {detail}", flush=True)
    if not ok:
        raise SystemExit(1)


def t1_testset() -> None:
    print("1. test set")
    import rl_env
    cases = pd.read_csv("configs/v3_testset.csv")
    spec = TS.TestsetSpec()
    params = {t: TS.truth_params(d) for t, d in TS.ATMOSPHERES.items()}
    sob = __import__("torch").quasirandom.SobolEngine(dimension=6, scramble=True, seed=spec.sobol_seed)
    pts = sob.draw(int(cases.sobol_index.max()) + 1).double().numpy()
    new = cases[cases.set == "new"].head(2)
    for r in new.itertuples():
        cand, why = TS._candidate(pts[int(r.sobol_index)], spec, params)
        same = cand is not None and all(abs(cand["row"][k] - getattr(r, k)) < 1e-6 for k in TS.IC_FIELDS)
        check(f"case {r.case} rebuilt from Sobol index {r.sobol_index}", same)
    r = new.iloc[0]
    Q = TS.entry_frame(math.radians(r.phi_deg), math.radians(r.lambda_deg), math.radians(r.chi_deg))
    back = TS.to_local(Q, math.radians(spec.target_phi_deg), math.radians(spec.target_lambda_deg))
    check("entry placement round trip", abs(back[0] - r.downrange_km) < 0.1 and abs(back[1] - r.crossrange_km) < 0.1)
    fp = TS.footprint(r.altitude_m, r.velocity_mps, math.radians(r.gamma_deg), params["mismatch_high"])
    for sdeg in (0.0, 30.0):
        env = rl_env.OrionEntryEnv(controller_mode="policy", dispersion=False, fast_atmosphere=True)
        D.apply_domain(env.params, TS.ATMOSPHERES["mismatch_high"])
        env.reset(seed=0, options={"ic_override": TS.ic_of(r)})
        done = False
        while not done:
            _, _, te, tr, info = env.step(np.array([sdeg / 90.0], dtype=np.float32))
            done = te or tr
        dr, cr = TS.to_local(Q, env.x[1], env.x[2])
        p = next(q for q in fp if q["sigma_deg"] == sdeg)
        err = math.hypot(dr - p["downrange_km"], cr - p["crossrange_km"])
        check(f"footprint vs env at bank {sdeg:g} deg", err < 10.0,
              f"(env {dr:.0f}/{cr:.0f} km, sweep {p['downrange_km']:.0f}/{p['crossrange_km']:.0f} km, {err:.1f} km apart)")
    check("bridge cases 5 and 14 flagged infeasible at x0.78", not cases.set_index("case").loc[[5, 14],
                                                                                               "feasible_out_of_range"].any())


def t2_equivalence(v2_ckpt: str) -> None:
    print("2. lockstep runner == the unchanged runners (3 cases)")
    from jepa import models as M1
    from jepa import models_v2 as V
    from jepa import mppi_v2 as P2
    cases = pd.read_csv("configs/v3_testset.csv").set_index("case")
    cfg_m = {"n_samples": 48, "n_sub": 2, "g_limit": 4.5, "horizon": 120}

    def job(c, seed, ctrl):
        return {"ic": TS.ic_of(cases.loc[c]), "domain": dict(D.NOMINAL_DOMAIN), "atmosphere": "nominal", "case": c,
                "set": "x", "controller": ctrl, "load": "4.5g", "seed": 0, "mppi_rep": 0, "record": True,
                "mppi": {**cfg_m, "seed": seed, "w_unc_km": 0.0}}

    # v1: the v3 runner against the unchanged mppi.run_episode, closed loop (same arithmetic -> exact).
    ad1 = F.load_adapter("v1", V1_CKPT, "cpu")
    worst = 0.0
    for c in (0, 16, 40):
        seed = F.crn_seed(0, c, "nominal", 0)
        ref = P.run_episode({"controller": "jepa", "ic": TS.ic_of(cases.loc[c]), "domain": dict(D.NOMINAL_DOMAIN),
                             "ckpt": V1_CKPT, "mppi": {**cfg_m, "seed": seed}, "record": True})
        got = F.run_lockstep_v3([job(c, seed, "v1")], ad1, {}, verbose=False)[0]
        n = min(20, len(ref["traj"]) - 1)
        worst = max(worst, float(np.abs(np.degrees(ref["traj"][:n, 10])
                                        - got["telemetry"]["sigma_cmd_deg"][1::5][:n]).max()))
    check("v1: same decisions as mppi.run_episode over the first 20 decisions", worst < 0.5,
          f"(max |d cmd| {worst:.2e} deg)")
    # v2 single member: (a) the v3 runner against mppi_v2.run_lockstep (the v2 study's runner), whole flight.
    one = V.select_members(V.load_checkpoint_v2(v2_ckpt), [0])
    ad2 = F.V2Adapter(one)
    worst, n_cmp = 0.0, 0
    for c in (0, 16, 40):
        seed = F.crn_seed(0, c, "nominal", 0)
        old = P2.run_lockstep([{"ic": TS.ic_of(cases.loc[c]), "domain": dict(D.NOMINAL_DOMAIN), "record": True}], one,
                              {**cfg_m, "seed": seed, "w_unc_km": 0.0}, verbose=False)[0]
        got = F.run_lockstep_v3([job(c, seed, "v2single")], ad2, {}, verbose=False)[0]
        a, b = np.degrees(old["traj"][:-1, 10]), got["telemetry"]["sigma_cmd_deg"][1::5]
        n = min(len(a), len(b))
        n_cmp += n
        worst = max(worst, float(np.abs(a[:n] - b[:n]).max()), abs(old["miss_km"] - got["result"]["miss_km"]))
    check("v2 single: identical to mppi_v2.run_lockstep for the whole flight", worst < 1e-4,   # float32 telemetry
          f"({n_cmp} decisions, max diff {worst:.1e})")
    # (b) the parallel predictor (lockstep) against the KV-cached step-by-step path run_episode uses.
    import rl_env
    env = rl_env.OrionEntryEnv(controller_mode="policy", dispersion=False, fast_atmosphere=True)
    env.reset(seed=0, options={"ic_override": TS.ic_of(cases.loc[16])})
    for _ in range(60):
        env.step(np.array([0.3], dtype=np.float32))
    y = torch.as_tensor(D.physical_state(env), dtype=torch.float32)
    yh, uh = y[None].repeat(one["wm"].rf, 1), torch.full((one["wm"].rf - 1,), 0.4)
    u = (torch.randn(64, 110, generator=torch.Generator().manual_seed(0)) * 0.6).clamp(-1.5, 1.5)
    from jepa.evaluate import crossing_point
    view = V.as_v1(one)
    with torch.no_grad():
        ref = M1.jepa_predict(view["wm"], view["pi"], one["stats"].norm_feats(yh)[None], yh[-1, :6].expand(64, 6),
                              uh.expand(64, -1), u, D.MACRO_DT_S, 2, h_stop=7600.0)["y6"]
        got, _ = ad2.rollout(ad2.context([yh]).repeat_interleave(64, 1), yh[-1, :6].expand(64, 6), uh.expand(64, -1),
                             u, 2, 7600.0)
    p1, p2 = crossing_point(ref.double(), 7600.0), crossing_point(got[:, 0].double(), 7600.0)
    d_km = float(torch.maximum((p1[0] - p2[0]).abs(), (p1[1] - p2[1]).abs()).max()) * 6371.0
    check("v2: parallel rollout == step-by-step rollout (landing points)", d_km < 0.01, f"(max {1000 * d_km:.2f} m)")
    # (c) closed loop against mppi.run_episode: the first decisions agree; float32 round-off between the two
    # predictor paths (above: sub-meter) is then amplified by the closed loop, so no long-horizon claim.
    P._MODEL_CACHE["smoke-v3-v2"] = view
    seed = F.crn_seed(0, 0, "nominal", 0)
    ref = P.run_episode({"controller": "jepa", "ic": TS.ic_of(cases.loc[0]), "domain": dict(D.NOMINAL_DOMAIN),
                         "ckpt": "smoke-v3-v2", "mppi": {**cfg_m, "seed": seed}, "record": True})
    got = F.run_lockstep_v3([job(0, seed, "v2single")], ad2, {}, verbose=False)[0]
    d = np.abs(np.degrees(ref["traj"][:3, 10]) - got["telemetry"]["sigma_cmd_deg"][1::5][:3])
    check("v2 single: first 3 decisions match mppi.run_episode", float(d.max()) < 0.5, f"(max {d.max():.1e} deg)")


def t3_crn() -> None:
    print("3. common random numbers")
    cfg = S.StudyConfig.smoke()
    jobs = S.build_jobs(cfg, S.Layout(tempfile.mkdtemp()))
    by = {}
    for j in jobs:
        if j["controller"] in S.LEARNED:
            by.setdefault((j["seed"], j["case"], j["atmosphere"], j["mppi_rep"]), set()).add(j["mppi"]["seed"])
    check("one noise seed per (seed, case, atmosphere, rep) across controllers", all(len(v) == 1 for v in by.values()))
    seeds = [next(iter(v)) for v in by.values()]
    check("different conditions get different noise", len(set(seeds)) == len(seeds), f"({len(seeds)} streams)")
    from jepa.mppi_v2 import EnsembleMPPI, MPPIConfigV2
    a = EnsembleMPPI(MPPIConfigV2(seed=seeds[0], n_samples=16), 0.26, 0.26, 7600.0)
    b = EnsembleMPPI(MPPIConfigV2(seed=seeds[0], n_samples=16, w_unc_km=0.0), 0.26, 0.26, 7600.0)
    a.begin(), b.begin()
    check("identical first proposals", torch.equal(a.propose(), b.propose()))


def t4_telemetry() -> dict:
    print("4. telemetry and RCS fuel")
    cases = pd.read_csv("configs/v3_testset.csv")
    row = cases.set_index("case").loc[16]
    job = {"ic": TS.ic_of(row), "domain": dict(TS.ATMOSPHERES["mismatch"]), "atmosphere": "mismatch", "case": 16,
           "set": "new", "controller": "classical", "load": "own"}
    t0 = time.time()
    rec = F.fly_classical(job)
    tel, r = rec["telemetry"], rec["result"]
    print(f"  classical case 16 mismatch: {r['outcome']} miss {r['miss_km']:.2f} km, fuel {r['rcs_fuel_kg']:.1f} kg, "
          f"{r['telemetry_rows']} rows ({time.time() - t0:.0f}s)")
    check("all telemetry fields", set(F.TELEM_FIELDS) | {"jet_s"} == set(tel))
    check("one row per second + entry", r["telemetry_rows"] == int(round(r["t_s"])) + 1 or
          abs(r["telemetry_rows"] - (r["t_s"] + 1)) <= 1)
    import telemetry
    check("jet-seconds x kg/jet-s == fuel", abs(tel["jet_s"][-1] * telemetry.rcs_propellant_kg(1.0)
                                                 - tel["rcs_fuel_kg"][-1]) < 1e-3)
    check("final fuel == simulator's", abs(tel["rcs_fuel_kg"][-1] - r["rcs_fuel_kg"]) < 1e-3)
    phases = r["fuel_kg_above60km"] + r["fuel_kg_60_25km"] + r["fuel_kg_below25km"]
    check("fuel by phase sums to the total", abs(phases - r["rcs_fuel_kg"]) < 1e-3)
    check("result JSON-serializable", bool(__import__("json").dumps(r, default=S._json_default)))
    return r


def t5_queue(cfg: S.StudyConfig, layout: S.Layout) -> None:
    print("5. job queue resume")
    jobs = [j for j in S.build_jobs(cfg, layout) if j["controller"] == "v1" and j["load"] == "4.5g"
            and j["seed"] == 0][:2]
    S.run_learned(cfg, jobs, layout, device="cpu", workers=1)
    check("both flights written", all(S.flight_done(layout, j) for j in jobs))
    first = layout.result_path(F.flight_id(jobs[0]))
    m0 = first.stat().st_mtime
    gone = layout.result_path(F.flight_id(jobs[1]))
    gone.unlink()
    S.run_learned(cfg, jobs, layout, device="cpu", workers=1)
    check("finished flight skipped", first.stat().st_mtime == m0)
    check("deleted flight re-flown", gone.exists())
    df = S.collect_flights(layout)
    check("flights.csv rebuilt", (layout.flights / "flights.csv").exists() and len(df) == 2)


def t6_training(cfg: S.StudyConfig, layout: S.Layout, data: dict) -> str:
    print("6. resumable training")
    m = torch.nn.Linear(4, 1)
    x = torch.randn(256, 4)
    y = x @ torch.tensor([[1.0], [-2.0], [0.5], [3.0]])
    calls = {"n": 0, "stop": None}

    def loss_fn():
        calls["n"] += 1
        if calls["stop"] is not None and calls["n"] > calls["stop"]:
            raise KeyboardInterrupt
        loss = ((m(x) - y) ** 2).mean()
        return loss, {"mse": float(loss.detach())}

    tmp = Path(tempfile.mkdtemp()) / "state.pt"
    c = S.StudyConfig.smoke()
    c.save_every_s, c.v1_log_every, c.v1_extend_tol = 0.0, 10, -1.0      # tol < 0: always extend
    calls["stop"] = 37                     # call 38 raises: 37 steps done and checkpointed
    try:
        S.fit_resumable(m, loss_fn, 60, tmp, lambda: float(((m(x) - y) ** 2).mean()), c, "unit")
    except KeyboardInterrupt:
        pass
    calls["stop"], calls["n"] = None, 0
    h = S.fit_resumable(m, loss_fn, 60, tmp, lambda: float(((m(x) - y) ** 2).mean()), c, "unit")
    check("resumed where it stopped", calls["n"] == 90 - 37, f"({calls['n']} steps after resume)")
    check("extension rule ran once", h["extended"]["extend"] and h["step"][-1] == 90)
    c.v1_extend_tol = 10.0
    tmp2 = Path(tempfile.mkdtemp()) / "state.pt"
    h = S.fit_resumable(torch.nn.Linear(4, 1), lambda: (torch.zeros((), requires_grad=True), {}), 30, tmp2,
                        lambda: 1.0, c, "flat")
    check("no extension when converged", not h["extended"]["extend"] and h["step"][-1] == 30)
    r1 = S.train_v1_seed(cfg, data, layout, 0, "cpu")
    r2 = S.train_v2_seed(cfg, data, layout, 0, "cpu")
    check("v1 and v2 seed checkpoints", layout.ckpt("v1", 0).exists() and layout.ckpt("v2", 0).exists())
    check("train state files cleaned up", not (layout.model_dir("v2", 0) / "train_state_wm.pt").exists())
    check("cached on rerun", S.train_v2_seed(cfg, data, layout, 0, "cpu")["cached"])
    torch.manual_seed(0)
    from jepa import models_v2 as V
    w0 = V.WorldModelV2(n_members=cfg.v2_members, **cfg.v2_arch).state_dict()
    torch.manual_seed(1)
    w1 = V.WorldModelV2(n_members=cfg.v2_members, **cfg.v2_arch).state_dict()
    check("seeds give different initializations", not all(torch.equal(w0[k], w1[k]) for k in w0))
    return str(layout.ckpt("v2", 0))


def t7_stats() -> None:
    print("7. statistics")
    check("Holm", np.allclose(A.holm([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06]))
    a = np.array([1.0] * 10 + [0.0] * 5)
    b = np.array([0.0] * 10 + [0.0] * 5)
    n10, n01, p = A.mcnemar_p(a, b)
    check("McNemar exact", (n10, n01) == (10, 0) and abs(p - 2 / 1024) < 1e-12)
    w, n, p = A.sign_test(np.arange(10.0), np.arange(10.0) + 1)
    check("sign test", (w, n) == (10, 10) and abs(p - 2 / 1024) < 1e-12)
    rng = np.random.default_rng(0)
    X = rng.random((3, 4, 20))
    r = A.hboot(X, X, "median", 500)
    check("bootstrap of identical arms is 0", r["diff"] == 0 and r["ci_lo"] == 0 and r["ci_hi"] == 0)
    r = A.hboot(X + 1.0, X, "mean", 500)
    check("bootstrap recovers a shift", abs(r["diff"] - 1.0) < 1e-9 and r["ci_lo"] > 0.99)
    r = A.hboot(X[:1] + 0.5, X, "mean", 500)
    check("deterministic arm broadcast over seeds", r["ci_lo"] < 0.5 + 0.2 and r["ci_hi"] > 0.5 - 0.2)


def main(wide: str) -> None:
    t_all = time.time()
    torch.set_num_threads(max(1, (os.cpu_count() or 2) // 2))
    root = Path(tempfile.mkdtemp(prefix="jepa_v3_smoke_"))
    cfg = S.StudyConfig.smoke()
    layout = S.Layout(root).make()
    S.freeze_config(cfg, layout)
    t1_testset()
    t3_crn()
    t7_stats()
    t4_telemetry()
    data = S.load_data(cfg, wide, "cpu")
    v2 = t6_training(cfg, layout, data)
    t2_equivalence(v2)
    t5_queue(cfg, layout)
    S.update_manifest(layout)
    check("manifest and README", (layout.study / "manifest.csv").exists() and (layout.study / "README.md").exists())
    print(f"v3 smoke test passed in {time.time() - t_all:.0f}s (scratch: {root})")
    shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--wide", default=os.environ.get("JEPA_WIDE_DATA", "jepa_runs/wide_smoke"),
                    help="wide-envelope episodes (pickle or shard directory); a small subset is enough")
    main(ap.parse_args().wide)
