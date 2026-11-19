"""
jepa/smoke_rl_v5.py: notebook 2's smoke checks (jepa/plan_rl_v5.md section 7). CPU, numbers meaningless.

    python jepa/smoke_rl_v5.py            (from StochasticEntrySim/; ~1 h: FNPEG's labels are slow at tiny batches)
    python jepa/smoke_rl_v5.py --units    (checks 1-11 without executing the notebook)

Works in a temporary folder (JEPA_RL2_SMOKE_DIR, default <tmp>/jepa_rl2_smoke): a fixture of notebook 1's folder
(FNPEG's tuned chosen.json and a synthetic flights.csv of fnpeg and v2single on the smoke suites, read only) and
notebook 2's smoke folder. Never touches Drive.
"""
from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import replace
from pathlib import Path

SIM = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SIM))

import numpy as np
import pandas as pd
import torch

from jepa import fly_v5 as F5
from jepa import rl2_v5 as R2
from jepa import rl_v5 as RL
from jepa import rlenv_fast_v5 as RF
from jepa import runguard_v5 as RG
from jepa import study_rl_v5 as SR
from jepa import study_v5 as S5

FNPEG_TUNED = {"vf_mps": 147.7, "sigma_f_deg": 34.01292100529351, "deadband_c0_deg": 1.3211845169029366,
               "deadband_c1_deg_per_mps": 0.0003967847332051459, "delta_nd": 0.02349383591168709,
               "k0_nd": 50.748696110460386, "filter_beta": 0.8820230264084713, "dt_s": 1.0}   # notebook 1's config 28
RESULTS = []


def check(n, name):
    def deco(fn):
        def run(*a, **k):
            t0 = time.time()
            try:
                fn(*a, **k)
                RESULTS.append((n, name, "PASS", round(time.time() - t0, 1), ""))
                print(f"[{n}] PASS {name} ({time.time() - t0:.0f} s)", flush=True)
            except Exception as e:
                RESULTS.append((n, name, "FAIL", round(time.time() - t0, 1), f"{type(e).__name__}: {e}"))
                print(f"[{n}] FAIL {name}: {type(e).__name__}: {e}", flush=True)
                import traceback
                traceback.print_exc()
        return run
    return deco


def tiny(**kw) -> RL.RLCfg:
    base = dict(n_envs=8, hidden=(32, 32), n_steps=4, epochs=2, minibatch=32, budget_decisions=64, budget_hours=0.5,
                env_pool=64, curve_points=0, save_every_s=1e9, critic_warmup_updates=1)
    base.update(kw)
    return RL.RLCfg(**base)


# ============================================================================
# Fixture: notebook 1's folder (read only for notebook 2)
# ============================================================================
def make_fixture(root: Path) -> Path:
    jroot = root / "jepa_v5_fixture"
    if jroot.exists():
        shutil.rmtree(jroot)
    J = S5.LayoutV5(jroot).make()
    J.write_json(J.tuning / "fnpeg" / "chosen.json", {"method": "fnpeg", "config": FNPEG_TUNED, "config_index": 28})
    cfg = SR.RLStudyConfig.smoke()
    v5 = cfg.v5()
    rng = np.random.default_rng(0)
    rows = []
    for su in SR.SUITES:
        conds = S5.conditions(v5, su)
        if su == "ORB1":
            conds = [c for c in conds if int(c["block"]) == 0]
        for c in conds:
            for arm, load, seeds in (("fnpeg", "own", [None]), ("v2single", "4.5g", list(cfg.rl_seeds))):
                for s in seeds:
                    miss = float(rng.lognormal(1.0, 1.0))
                    rows.append({"flight_id": f"{su}/{arm}/{s}/{c['case']}/{c['atmosphere']}", "arm": arm,
                                 "controller": arm, "load": load, "suite": su, "seed": -1 if s is None else s,
                                 "case": int(c["case"]), "atmosphere": c["atmosphere"], "miss_km": miss,
                                 "outcome": "drogue", "peak_g": 4.0, "stage": su, "mppi_rep": -1,
                                 "feasible": c.get("feasible"), "block": c.get("block"),
                                 "target_index": c.get("target_index")})
    J.write_csv(pd.DataFrame(rows), J.flights / "flights.csv")
    return jroot


# ============================================================================
# Unit checks (plan section 7)
# ============================================================================
@check(1, "resid_fnpeg at init: mean exactly 0, command = FNPEG's bit for bit (fast env and a full-sim flight)")
def c1(fn):
    kw = {"bound_deg": 20.0, "init_std": 0.1}
    r = R2.ri_fast(tiny(), fn, kw, n=4, steps=30)
    assert r["passed"], r
    for c in (0.0, 0.3, -1.2, R2.SIGMA_MAX_RAD, -R2.SIGMA_MAX_RAD):
        assert R2.compose(c, R2.sign_eff(c, -1.0), 0.0) == c
    cond = S5.conditions(SR.RLStudyConfig.smoke().v5(), "DEV48")[0]
    j = S5.guidance_job(cond, "gate_ri", "gates", "fnpeg", fn)
    ra, rb = [], []
    F5.run_guidance_lockstep([j], F5.make_law("fnpeg", fn), on_done=ra.append, verbose=False)
    F5.run_guidance_lockstep([j], R2.ZeroResidLaw(fn, tiny(), kw), on_done=rb.append, verbose=False)
    same, diff = R2.flights_identical(ra[0], rb[0])
    assert same, diff
    # and a nonzero correction does change the flight (the comparison can fail)
    actor = R2.make_actor("resid_fnpeg", tiny(), kw)
    with torch.no_grad():
        actor.net[-1].bias.fill_(0.5)
    rc = []
    F5.run_guidance_lockstep([j], R2.ResidLaw(None, fn, 20.0, actor=actor.eval()), on_done=rc.append, verbose=False)
    assert not R2.flights_identical(ra[0], rc[0])[0]


@check(2, "residual env: FNPEG labels = rl_v5.FastExpert's, in-process = workers, state round-trips")
def c2(fn):
    kw = {"bound_deg": 20.0}

    def run(workers, steps=14):
        e = R2.make_env("resid_fnpeg", tiny(label_workers=workers), 3, "cpu", fn, kw)
        ref = RL.FastExpert(e.env, tiny(), fn) if workers == 1 else None
        if ref is not None:
            lab = ref.label()
            assert np.allclose(lab.numpy(), np.clip(e.q[:, 0] / R2.SIGMA_MAX_RAD, -1, 1), atol=1e-6)
        qs = []
        for k in range(steps):
            _, _, d, _ = e.step(torch.zeros(8))
            qs.append(e.q.copy())
            if ref is not None:
                ref.reset_rows(torch.nonzero(d).flatten())
                lab = ref.label()
                assert np.allclose(lab.numpy(), np.clip(e.q[:, 0] / R2.SIGMA_MAX_RAD, -1, 1), atol=1e-6), k
        st = e.get_state()
        a = [e.step(torch.full((8,), 0.4))[1].clone() for _ in range(6)]
        e.set_state(st)
        b = [e.step(torch.full((8,), 0.4))[1].clone() for _ in range(6)]
        assert all(torch.equal(x, y) for x, y in zip(a, b)), f"state round trip, workers={workers}"
        e.close()
        return qs
    q1, q2 = run(1), run(2)
    assert all(np.array_equal(x, y) for x, y in zip(q1, q2)), "workers differ from in-process"


@check(3, "liu_track: reference at the start, sign only outside the deadband, finite reward, terminations")
def c3(fn):
    lc = R2.liu_cfg_from(fn, ds_frac=0.3)
    e = R2.LiuEnv(RF.FastEnvCfg(n_envs=8, seed=5, pool=64), "cpu", lc)
    T = e._terms()
    bound = (e.v0 - lc.vf_mps) * math.exp(-1 / lc.ds_frac)
    assert bool(((T["v_ref"] - e.v0).abs() <= bound + 1e-6).all()), "v_ref at s0"
    g = torch.Generator().manual_seed(0)
    for _ in range(150):
        T = e._terms()
        old = e.sign.clone()
        _, r, d, info = e.step(torch.rand(8, generator=g) * 2 - 1)
        assert torch.isfinite(r).all()
        ch = (old != 0) & (e.sign != old) & ~d
        assert bool((T["dpsi"][ch].abs() > T["db"][ch]).all()), "a reversal inside the deadband"
    ev = R2.LiuEnv(e.env.cfg, "cpu", replace(lc, train_terms=False))
    tr = R2.LiuEnv(e.env.cfg, "cpu", replace(lc, n_max_g=0.05))
    f_ev = f_tr = 0
    for _ in range(40):
        f_ev += int(ev.step(torch.zeros(8))[3]["forced"].sum())
        f_tr += int(tr.step(torch.zeros(8))[3]["forced"].sum())
    assert f_ev == 0 and f_tr > 0, (f_ev, f_tr)


@check(4, "PPO on the new environments trains and resumes to the same weights")
def c4(fn, root):
    L = S5.LayoutV5(root / "c4").make()
    for arm, kw in (("liu_track", {"ds_frac": 0.3}), ("resid_fnpeg", {"bound_deg": 20.0, "init_std": 0.1})):
        full, half = tiny(budget_decisions=4 * 32), tiny(budget_decisions=2 * 32)
        da, db = L.policies / arm / "a", L.policies / arm / "b"
        R2.ppo_train_env(L, da, full, "cpu", arm, arm, fn, kw)
        R2.ppo_train_env(L, db, half, "cpu", arm, arm, fn, kw)
        os.remove(db / "policy.pt")
        R2.ppo_train_env(L, db, full, "cpu", arm, arm, fn, kw)
        A, B = (torch.load(d / "policy.pt", weights_only=False)["state"] for d in (da, db))
        assert all(torch.equal(A[k], B[k]) for k in A), f"{arm}: resumed weights differ"
        assert R2.policy_meta(da / "policy.pt")["arm"] == arm


@check(5, "the fixed episodes are the same entries for every arm")
def c5(fn):
    xs = []
    for arm, kw in (("resid_fnpeg", {"bound_deg": 20.0}), ("liu_track", {}), ("ppo_hist", {})):
        e = R2.make_env(arm, tiny(), 12345, "cpu", fn, kw, n_envs=8, pool=8, train=False)
        inner = getattr(e, "env", e)
        xs.append((inner.x.clone(), inner.phi_t.clone(), inner.S_ref.clone()))
        R2._close(e)
    assert all(all(torch.equal(a, b) for a, b in zip(xs[0], x)) for x in xs[1:])


@check(6, "gate verdicts on forced passes and fails; a failed RP is a result")
def c6():
    actor = R2.make_actor("ppo_hist", tiny(), {})
    up = [{"episodes": 5, "ep_return": float(i), "success_5km_fast": 0.1} for i in range(9)]
    flat = [{"episodes": 5, "ep_return": 1.0, "success_5km_fast": 0.0} for _ in range(9)]
    assert R2.rp_verdict("resid_fnpeg", up, actor, 0.79, 0.80)["passed"]
    assert not R2.rp_verdict("resid_fnpeg", up, actor, 0.77, 0.80)["passed"]
    assert not R2.rp_verdict("resid_fnpeg", flat, actor, 0.85, 0.80)["passed"]
    assert R2.rp_verdict("liu_track", flat, actor, 0.41, 0.80)["passed"]
    assert not R2.rp_verdict("liu_track", up, actor, 0.39, 0.80)["passed"]
    assert R2.rp_verdict("ppo_hist", up, actor, 0.41, 0.80)["passed"]
    assert not R2.rp_verdict("ppo_hist", flat, actor, 0.41, 0.80)["passed"]     # old RV3: no rise, no success
    with torch.no_grad():
        actor.net[-1].bias.fill_(float("nan"))
    assert not R2.rp_verdict("liu_track", flat, actor, 0.9, 0.8)["passed"]


@check(7, "routing chooses the cheaper runtime and waits while a rate is missing")
def c7(root):
    assert R2.route({"cpu": {"decisions_per_s": 300, "units_per_hour": 0.26}}) is None
    r = R2.route({"cpu": {"decisions_per_s": 300, "units_per_hour": 0.26},
                  "L4": {"decisions_per_s": 400, "units_per_hour": 1.98}})
    assert r["runtime"] == "cpu", r
    r = R2.route({"cpu": {"decisions_per_s": 3700, "units_per_hour": 0.26},
                  "L4": {"decisions_per_s": 48500, "units_per_hour": 1.98}})
    assert r["runtime"] == "L4", r
    L = S5.LayoutV5(root / "c7").make()
    assert SR.route(L, "ppo_hist") is None
    L.write_json(L.checks / "rl2_rates" / "ppo_hist_cpu.json", {"decisions_per_s": 3700, "units_per_hour": 0.26})
    assert SR.route(L, "ppo_hist") is None
    L.write_json(L.checks / "rl2_rates" / "ppo_hist_L4.json", {"decisions_per_s": 48500, "units_per_hour": 1.98})
    assert SR.route(L, "ppo_hist")["runtime"] == "L4" and (L.checks / "route_ppo_hist.json").exists()
    cfg = SR.RLStudyConfig()
    assert SR.wall_cap_h(cfg, L, "ppo_hist", 200_000_000) == 3.0          # 1.1 h projected: the 3 h floor
    assert abs(SR.wall_cap_h(cfg, L, "ppo_hist", 2_000_000_000) - 1.5 * 2e9 / 48500 / 3600) < 1e-9


@check(8, "jobs: flight ids, seed rule, pairing keys; a dropped or gated-out arm has no jobs")
def c8(fn, root):
    L = S5.LayoutV5(root / "c8").make()
    cfg = SR.RLStudyConfig()
    for arm in ("resid_fnpeg", "liu_track"):
        L.write_json(L.checks / f"RP_{arm}.json", {"passed": arm == "resid_fnpeg"})
        L.write_json(L.tuning / SR.tmethod(arm) / "chosen.json", {"config": {"bound_deg": 10.0} if arm ==
                                                                  "resid_fnpeg" else {}})
    v5 = cfg.v5()
    a = SR.flight_jobs(cfg, L, "resid_fnpeg", "A256", fn)
    assert len(a) == 5 * len(S5.conditions(v5, "A256"))
    j = a[0]
    assert F5.flight_id_v5(j).startswith("A256/resid_fnpeg/seed0/"), F5.flight_id_v5(j)
    assert j["model"]["rl2"]["bound_deg"] == 10.0 and j["model"]["rl2"]["fnpeg_cfg"] == fn
    fj = S5.guidance_job(S5.conditions(v5, "A256")[0], "fnpeg", "A256", "fnpeg", fn)
    assert F5.cond_key(fj) == F5.cond_key(j)
    o = SR.flight_jobs(cfg, L, "resid_fnpeg", "ORB1", fn)
    assert len(o) == 1000 and all(int(x["seed"]) == int(x["case"]) % 5 for x in o)
    assert all(int(x["block"]) == 0 for x in o)
    d = SR.flight_jobs(cfg, L, "resid_fnpeg", "DS1", fn)
    assert all(int(x["seed"]) == int(x["case"]) % 5 for x in d)
    assert SR.flight_jobs(cfg, L, "liu_track", "A256", fn) == []          # failed its gate
    assert SR.flight_jobs(cfg, L, "dagger_fnpeg", "A256", fn) == []       # dropped
    assert SR.flight_jobs(cfg, L, "ppo_hist", "A256", fn) == []           # not gated yet


@check(9, "stop rules: each trigger writes STOPPED.json; Setup refuses without ACK_STOP, continues with it")
def c9(root):
    L = S5.LayoutV5(root / "c9").make()
    cfg = SR.RLStudyConfig.smoke()
    for trig in SR.STOP_TRIGGERS:
        try:
            SR.stop(L, trig, {"forced": True}, log=lambda *a: None)
            raise AssertionError("stop did not raise")
        except SR.Stopped:
            pass
        assert json.loads((L.study / "STOPPED.json").read_text())["trigger"] == trig
        for ack in (None, "wrong"):
            try:
                SR.check_stopped(L, ack, log=lambda *a: None)
                raise AssertionError("Setup went on without ACK_STOP")
            except SR.Stopped:
                pass
        assert SR.check_stopped(L, trig, log=lambda *a: None)["trigger"] == trig
        assert not (L.study / "STOPPED.json").exists()
    assert len(list((L.study / "stops").glob("*.json"))) == len(SR.STOP_TRIGGERS)
    try:                                                     # stage_exception from a raising stage
        with SR.stage(L, cfg, "boom", log=lambda *a: None):
            raise ValueError("boom")
        raise AssertionError("no stop")
    except SR.Stopped:
        assert json.loads((L.study / "STOPPED.json").read_text())["trigger"] == "stage_exception"
    SR.check_stopped(L, "stage_exception", log=lambda *a: None)
    keep = list(RG.SUMMARY.failures)
    RG.SUMMARY.failures = [{"kind": "x", "what": str(i), "reason": "r"} for i in range(3)]
    try:
        SR.check_failures(L, cfg, log=lambda *a: None)
        raise AssertionError("no stop on 3 failures")
    except SR.Stopped:
        assert json.loads((L.study / "STOPPED.json").read_text())["trigger"] == "failures"
    finally:
        RG.SUMMARY.failures = keep
    SR.check_stopped(L, "failures", log=lambda *a: None)
    # flight_errors: 5 of an arm's first 20 flights raising
    L.write_json(L.checks / "RP_ppo_hist.json", {"passed": True})
    L.write_json(L.tuning / SR.tmethod("ppo_hist") / "chosen.json", {"config": {}})
    pol = L.policy_dir("ppo_hist", 0)
    RL.save_policy(L, pol / "policy.pt", R2.make_actor("ppo_hist", tiny(), {}), {})
    js = SR.flight_jobs(replace(cfg, rl_seeds=(0,)), L, "ppo_hist", "A256", FNPEG_TUNED)[:2]
    os.environ[S5.FAULT_ENV] = "flight:ppo_hist"
    try:
        try:
            SR.run_flights(replace(cfg, flight_err_first=2, flight_err_max=2), L, js, 1, lambda *a: None, "c9")
            raise AssertionError("no flight_errors stop")
        except SR.Stopped:
            assert json.loads((L.study / "STOPPED.json").read_text())["trigger"] == "flight_errors"
    finally:
        os.environ.pop(S5.FAULT_ENV, None)


@check(10, "read-only: a write into notebook 1's folder raises")
def c10(root, jroot):
    L, J = SR.layout(root / "c10", str(jroot))
    for f in (lambda: L.write_json(Path(jroot) / "x.json", {}), lambda: L.atomic_bytes(J.tuning / "y.bin", b"1")):
        try:
            f()
            raise AssertionError("write into notebook 1's folder did not raise")
        except S5.ReadOnlyViolation:
            pass


@check(11, "the notebook: CLAUDE.md release code, opens on a CPU runtime, an old revision stops")
def c11():
    import nbformat
    from jepa import build_notebook_v5_rl as BN
    nb = nbformat.read(str(BN.build()), as_version=4)
    src = "\n".join(c.source for c in nb.cells)
    for s in ("flush_and_unmount", "runtime.unassign", "post_run_cell", "KeyboardInterrupt", "_ip.events.unregister",
              "AUTO_DISCONNECT = False", "disconnects itself", "ACK_STOP"):
        assert s in src, s
    assert nb.cells[-1].source.strip() == 'release_runtime("all cells finished")'
    assert "accelerator" not in nb.metadata and "gpuType" not in nb.metadata.get("colab", {})
    assert BN.REV == SR.NOTEBOOK_REV_RL
    try:
        SR.notebook_check(SR.NOTEBOOK_REV_RL - 1)
        raise AssertionError("old revision accepted")
    except RuntimeError:
        pass


@check(12, "the analysis builds the RQ3 table with Holm over 3 (a missing arm counts as p = 1)")
def c12(root, run_dir=None):
    L = S5.LayoutV5(root / "c12").make()
    cfg = SR.RLStudyConfig.smoke()
    c = SR.criteria_rq3(cfg, L, pd.DataFrame())
    assert list(c.id) == [f"R-{a}" for a in SR.ARMS] and (c.p == 1.0).all() and (c.p_holm == 1.0).all()
    t = pd.DataFrame([{"id": "R-ppo_hist", "metric": "success5", "atmosphere": "all", "diff": 20.0, "ci_lo": 5.0,
                       "ci_hi": 30.0, "p_boot": 0.01, "n_seeds": 5, "seeds_holding": 5, "complete": True,
                       "n_pairs": 1280, "n_expected": 1280}])
    c = SR.criteria_rq3(cfg, L, t)
    r = c[c.id == "R-ppo_hist"].iloc[0]
    assert abs(r.p_holm - 0.03) < 1e-12 and r.status == "met", r.to_dict()  # 3 x 0.01 with two p = 1 arms
    if run_dir is not None:
        crit = pd.read_csv(Path(run_dir) / "05_analysis" / "tables" / "RQ3_criteria.csv")
        assert len(crit) == 3 and crit.p_holm.notna().all(), crit
        assert (Path(run_dir) / "05_analysis" / "RQ3_summary.md").exists()


# ============================================================================
# The notebook, executed (smoke mode)
# ============================================================================
def run_notebook(root: Path, jroot: Path, tag: str) -> Path:
    rl_root = root / "rl_v5_smoke"
    out = root / "nb_runs"
    out.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "JEPA_RL2_MODE": "smoke", "JEPA_RL2_ROOT": str(rl_root), "JEPA_V5_ROOT": str(jroot),
           SR.GATES_PASS_ENV: "1", "PYTHONIOENCODING": "utf-8"}
    env.pop("JEPA_V3_ROOT", None)
    env.pop("JEPA_V4_ROOT", None)
    t0 = time.time()
    with open(out / f"{tag}.out", "w", encoding="utf-8") as so, open(out / f"{tag}.err", "w", encoding="utf-8") as se:
        r = subprocess.run([sys.executable, "-m", "nbconvert", "--to", "notebook", "--execute",
                            "--ExecutePreprocessor.timeout=-1", "--ExecutePreprocessor.startup_timeout=180",
                            "--output-dir", str(out), "--output", f"{tag}.ipynb", "ReEntryAI_RL_v5.ipynb"],
                           cwd=str(SIM), env=env, stdout=so, stderr=se)
    print(f"[notebook] {tag}: exit {r.returncode} in {(time.time() - t0) / 60:.1f} min ({out / tag}.ipynb)",
          flush=True)
    if r.returncode != 0:
        raise RuntimeError(f"notebook run {tag} failed: see {out / tag}.err")
    return rl_root


@check(13, "smoke run of the notebook: every stage, then a second Run all that only skips and reports")
def c13(root, jroot, box):
    rl_root = run_notebook(root, jroot, "run1")
    L = S5.LayoutV5(rl_root)
    assert not (L.study / "STOPPED.json").exists(), (L.study / "STOPPED.json").read_text()
    ri = json.loads((L.checks / "RI_resid_fnpeg.json").read_text())
    assert ri["passed"] and ri["n_dev48_conditions"] >= 1, ri
    for arm in SR.ARMS:
        assert (L.checks / f"RP_{arm}.json").exists(), arm
        assert (L.checks / f"route_{arm}.json").exists(), arm
        assert (L.tuning / SR.tmethod(arm) / "chosen.json").exists() or SR.tuning_stopped(L, arm), arm
    st = json.loads((L.study / "status.json").read_text())
    print(pd.DataFrame(st["arms"]).to_string(), flush=True)
    t0 = time.time()
    run_notebook(root, jroot, "run2")
    assert time.time() - t0 < 1800
    box["rl_root"] = rl_root


def main() -> int:
    units = "--units" in sys.argv
    root = Path(os.environ.get("JEPA_RL2_SMOKE_DIR") or Path(tempfile.gettempdir()) / "jepa_rl2_smoke")
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    os.environ.setdefault("JEPA_V5_GIT", "smoke")
    jroot = make_fixture(root)
    fn = dict(FNPEG_TUNED)
    c11()
    c6()
    c7(root)
    c9(root)
    c10(root, jroot)
    c12(root)
    c5(fn)
    c8(fn, root)
    c3(fn)
    c1(fn)
    c4(fn, root)
    c2(fn)
    box = {}
    if not units:
        c13(root, jroot, box)
        if box.get("rl_root"):
            c12(root, box["rl_root"])
    print("\n".join(f"[{n}] {s} {name} ({t} s) {e}" for n, name, s, t, e in sorted(RESULTS, key=lambda r: r[0])))
    bad = [r for r in RESULTS if r[2] != "PASS"]
    print(f"{len(RESULTS) - len(bad)}/{len(RESULTS)} checks passed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
