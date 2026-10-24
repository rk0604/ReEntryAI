"""
jepa/smoke_v4.py

CPU smoke test of the v4 study code (jepa/v4_plan.md section 12). The numbers
it produces are meaningless (tiny models, 2 cases, 2 seeds); the checks are what
matter. It builds everything it needs:

    data      ~200 wide-envelope flights (data.collect_shards, WIDE_DOMAIN), cached
    v3        a smoke v3 study folder made by the v3 code (study_v3.StudyConfig.smoke():
              tiny v1 and v2, 2 seeds, 2 cases, its classical and learned flights), cached
    v4        the v4 notebook executed in smoke mode on top of it (check 9)

Checks (numbered as in the plan):
  1  default-argument equivalence: every changed v3 function (models_v2, models)
     returns bit-identical outputs to the pristine jepa-v3 code (git show)
  2  E0c: run_lockstep_v4 == run_lockstep_v3 on the CPU, v1 and v2 adapters
  3  physics port vs math_3d.nominal_eom_step; a whole rollout step (substeps,
     bank-rate limit, loads, units) against math_3d; mppi_truth with nominal
     parameters == mppi_nominal; mppi_prior == an untrained PIProberV2
  4  CRN: physics arm index s and v2 seed s get the same noise stream
  5  every variant builds, trains a few steps, saves, resumes after an
     interruption (and after a torn save), with its claimed parameter count
     (A7a, A7b within 10% of target at full size); a resumed run equals an
     uninterrupted one; sliced windows equal rebuilt ones
  6  job queue: finished flights skipped, a deleted result re-flown, the matched
     g limits read from v3's classical flights
  7  analysis: every table and figure, a complete criteria table, the bootstrap
     cache reused on a rerun, recomputed (only where needed) after a result changes
  8  the report PDF builds without system DejaVu fonts
  9  the notebook runs end to end in smoke mode, then a second run finds
     everything done in under a minute
  10 the v4 code refuses to write under the v3 folder
  +  the v3 notebooks still build identically (except nbformat's random cell ids)

    python jepa/smoke_v4.py [--work <dir>] [--wide <wide-envelope shards>] [--workers N]
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

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))
os.chdir(SIM_DIR)
os.environ.setdefault("TMPDIR", tempfile.gettempdir())       # study_v3's progress files

import numpy as np
import pandas as pd
import torch

from jepa import ablate_v4 as AB
from jepa import data as D
from jepa import fly_v3 as F
from jepa import fly_v4 as F4
from jepa import models as M
from jepa import models_v2 as V
from jepa import physics_v4 as PH
from jepa import study_v3 as S
from jepa import study_v4 as S4

RESULTS = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {detail}", flush=True)
    if not ok:
        raise SystemExit(1)


# ============================================================================
# 1. Default-argument equivalence with the pristine jepa-v3 code
# ============================================================================
def _pristine(tmp: Path):
    """Import jepa-v3's models.py and models_v2.py side by side, as package jepa_v3ref."""
    pkg = tmp / "jepa_v3ref"
    pkg.mkdir(parents=True, exist_ok=True)
    (pkg / "__init__.py").write_text("")
    for name in ("models", "models_v2"):
        src = None
        for ref in ("jepa-v3", "origin/jepa-v3"):
            r = subprocess.run(["git", "-C", str(SIM_DIR), "show", f"{ref}:StochasticEntrySim/jepa/{name}.py"],
                               capture_output=True, text=True, encoding="utf-8")
            if r.returncode == 0:
                src = r.stdout
                break
        if src is None:
            return None
        src = src.replace("from jepa.models import", "from jepa_v3ref.models import")
        (pkg / f"{name}.py").write_text(src, encoding="utf-8")
    sys.path.insert(0, str(tmp))
    import importlib
    return importlib.import_module("jepa_v3ref.models"), importlib.import_module("jepa_v3ref.models_v2")


def t1_defaults(tmp: Path) -> None:
    print("1. default-argument equivalence with jepa-v3 (bit-identical)")
    ref = _pristine(tmp)
    if ref is None:
        check("pristine jepa-v3 sources available (git show jepa-v3:...)", False)
    M0, V0 = ref
    g = torch.Generator().manual_seed(0)
    u = (torch.rand(3, 30, generator=g) - 0.5) * 3.0
    check("models_v2.action_features", torch.equal(V.action_features(u), V0.action_features(u)))
    for kw in ({"n_members": 2, "d_latent": 16, "enc_dim": 32, "enc_layers": 2, "enc_heads": 2, "enc_band": 5,
                "pred_dim": 48, "pred_layers": 2, "pred_heads": 2, "act_chunk": 4, "mlp_ratio": 2}, {"n_members": 1}):
        torch.manual_seed(3)
        a = V0.WorldModelV2(**kw)
        torch.manual_seed(3)
        b = V.WorldModelV2(**kw)
        same_init = all(torch.equal(x, y) for x, y in zip(a.state_dict().values(), b.state_dict().values()))
        check(f"WorldModelV2({len(kw)} args): same init from the same seed, same hparams",
              same_init and a.hparams == b.hparams and list(a.state_dict()) == list(b.state_dict()))
        rf, M_ = b.rf, b.M
        feats = torch.randn(4, rf + 12, V.N_FEAT_V2, generator=g)
        uu = (torch.rand(4, rf - 1 + 12, generator=g) - 0.5) * 3
        with torch.no_grad():
            Sa, Sb = a.encode_states(feats), b.encode_states(feats)
            Ta, Tb = a.action_tokens(uu), b.action_tokens(uu)
            Pa, Pb = a.predict(Sa[:, :, 0], Ta), b.predict(Sb[:, :, 0], Tb)
        check("  encode_states, action_tokens, predict", torch.equal(Sa, Sb) and torch.equal(Ta, Tb) and torch.equal(Pa, Pb))
        pa, pb = V0.PIProberV2(M_, a.d_latent), V.PIProberV2(M_, b.d_latent)
        pb.load_state_dict(pa.state_dict())
        for p in (pa.net[-1].weight, pb.net[-1].weight):
            torch.nn.init.normal_(p, std=0.1, generator=torch.Generator().manual_seed(5))
        y6 = torch.tensor([60e3, 0.2, 0.3, 4000.0, -0.05, 0.7]).expand(M_, 4, 6).contiguous()
        with torch.no_grad():
            ya, aa = pa.step(Sa[:, :, -1], y6, uu[:, -1].expand(M_, 4), 5.0, 2)
            yb, ab = pb.step(Sb[:, :, -1], y6, uu[:, -1].expand(M_, 4), 5.0, 2)
            ea = V0.ensemble_rollout(a, pa, Sa[:, :, -1], y6[0], uu[:, :rf - 1], uu[:, rf - 1:], 5.0, 2, h_stop=7600.0)
            eb = V.ensemble_rollout(b, pb, Sb[:, :, -1], y6[0], uu[:, :rf - 1], uu[:, rf - 1:], 5.0, 2, h_stop=7600.0)
        check("  PIProberV2.step, ensemble_rollout", torch.equal(ya, yb) and torch.equal(aa, ab) and
              all(torch.equal(ea[k], eb[k]) for k in ea))
        ds = {"feats": torch.randn(M_, 6, rf + 9, V.N_FEAT_V2, generator=g),
              "u": (torch.rand(M_, 6, rf - 1 + 9, generator=g) - 0.5) * 3,
              "mask": (torch.rand(M_, 6, 9, generator=g) > 0.2).float()}
        sig0, sig1 = V0.SIGRegV2(), V.SIGRegV2()
        torch.manual_seed(11)
        la, loga = V0.wm_loss_v2(a, sig0, ds, 0.08, sig_min_rows=2)
        rng_a = torch.get_rng_state()
        torch.manual_seed(11)
        lb, logb = V.wm_loss_v2(b, sig1, ds, 0.08, sig_min_rows=2)
        check("  wm_loss_v2 (per-step SIGReg): same loss, same random numbers drawn",
              torch.equal(la, lb) and loga["pred"] == logb["pred"] and loga["sigreg"] == logb["sigreg"]
              and torch.equal(rng_a, torch.get_rng_state()))
        if kw.get("n_members") == 1:
            c = V0.WorldModelV2(**b.hparams)
            c.load_state_dict(b.state_dict())
            check("  a default model's hparams and weights load in the pristine code", True)
    torch.manual_seed(7)
    p0 = M0.PIProber(32)
    torch.manual_seed(7)
    p1 = M.PIProber(32)
    torch.nn.init.normal_(p0.net[-1].weight, std=0.1, generator=torch.Generator().manual_seed(1))
    p1.load_state_dict(p0.state_dict())
    s = torch.randn(5, 32, generator=g)
    y6 = torch.tensor([60e3, 0.2, 0.3, 4000.0, -0.05, 0.7]).expand(5, 6)
    with torch.no_grad():
        same = torch.equal(p0.corrections(s, u[0, :5]), p1.corrections(s, u[0, :5]))
        a0, a1 = p0.step(s, y6, u[0, :5], 5.0, 2), p1.step(s, y6, u[0, :5], 5.0, 2)
    check("models.PIProber (default): init, corrections, step", same and p0.hparams == p1.hparams and
          torch.equal(a0[0], a1[0]) and torch.equal(a0[1], a1[1]))
    w0 = M0.JEPAWorldModel()
    with torch.no_grad():
        o0 = M0.jepa_predict(w0, p0, torch.randn(1, w0.rf, 10, generator=torch.Generator().manual_seed(2)), y6,
                             torch.zeros(5, w0.rf - 1), u[0:1, :20].expand(5, 20).contiguous(), 5.0, 2)
        o1 = M.jepa_predict(w0, p1, torch.randn(1, w0.rf, 10, generator=torch.Generator().manual_seed(2)), y6,
                            torch.zeros(5, w0.rf - 1), u[0:1, :20].expand(5, 20).contiguous(), 5.0, 2)
    check("models.jepa_predict through the default PIProber", all(torch.equal(o0[k], o1[k]) for k in ("y6", "acc")))


# ============================================================================
# 3. Physics
# ============================================================================
def t3_physics() -> None:
    print("3. physics model")
    r = PH.port_check(1000)
    check("torch port vs math_3d.nominal_eom_step (1,000 states, float64)", r["max_rel_err"] < 1e-9,
          f"(max relative error per step {r['max_rel_err']:.1e}, increment {r['max_rel_err_increment']:.1e}, "
          f"load {r['max_rel_err_gload']:.1e})")
    from math_3d import nominal_aero_forces_from_state, nominal_eom_step
    p = D.apply_domain(dict(PH.nominal_params()), {"density_scale": 1.1, "lift_scale": 0.9, "mass_scale": 1.03})
    Re = float(__import__("constants").RADIUS_EARTH)
    y0 = torch.tensor([[70_000.0, 0.2, 0.25, 6000.0, math.radians(-2.0), 0.8]], dtype=torch.float64)
    s0, cmd = math.radians(-30.0), math.radians(70.0)
    trace = []
    atm = PH.atm_table("cpu")
    f64 = lambda v: torch.tensor(v, dtype=torch.float64)
    y6, gl = PH.rollout_physics(y0, f64([s0]), f64([[cmd, cmd, -cmd]]), f64([p["ref_area_m2"]]), f64([p["mass_kg"]]),
                                f64([p["cpas_lift_scale"]]), 7600.0, atm, bank_trace=trace)
    x, sig, err, gerr = [y0[0, 0].item() + Re] + y0[0, 1:].tolist(), s0, 0.0, 0.0
    for k, c in enumerate((cmd, cmd, -cmd)):
        gm = 0.0
        for ds in PH.SUBSTEPS_S:
            a = nominal_aero_forces_from_state(x=x, sigma_rad=sig, params=p)
            gm = max(gm, math.hypot(a["drag_mag_N"], a["lift_mag_N"]) / p["mass_kg"] / D.G0)
            x = nominal_eom_step(x=x, sigma_rad=sig, params=p, dt_s=ds)["x_next"]
            sig = sig + max(-PH.RATE_MAX_RPS * ds, min(PH.RATE_MAX_RPS * ds, c - sig))
        got = y6[0, k + 1].tolist()
        err = max(err, max(abs(got[i] - (x[0] - Re if i == 0 else x[i])) / max(abs(x[i] - (Re if i == 0 else 0.0)), 1.0)
                           for i in range(6)))
        gerr = max(gerr, abs(float(gl[0, k]) - gm) / gm)
    check("one rollout vs math_3d stepped by hand (2+2+1 s Euler, lagged bank, max load per 5 s)", err < 1e-9 and gerr < 1e-9,
          f"(state {err:.1e}, load {gerr:.1e}; gload shape {tuple(gl.shape)}, in g: {float(gl[0, 0]):.3f})")
    steps = [(ds, float(b[0])) for ds, b in trace]
    rate_ok = all(abs(b1 - b0) <= PH.RATE_MAX_RPS * ds + 1e-12 for (ds, b0), (_, b1) in zip(steps, steps[1:]))
    check("bank moves at most 20 deg/s toward each command",
          rate_ok and abs(steps[0][1] - s0) < 1e-15 and abs(steps[5][1] - cmd) < 1e-12,
          f"(bank per substep, deg: {[round(math.degrees(b), 1) for _, b in steps]})")
    nom, tru = PH.PhysicsAdapter("nominal"), PH.PhysicsAdapter("truth")
    yh = torch.zeros(3, 10)
    yh[:, 6] = 0.3
    u = (torch.rand(32, 60, generator=torch.Generator().manual_seed(0)) - 0.5) * 3
    y0b = torch.tensor([95e3, 0.1, 0.2, 7500.0, math.radians(-1.8), 0.4]).expand(32, 6)
    job = {"domain": dict(D.NOMINAL_DOMAIN)}
    a = nom.rollout(nom.context([yh], [job]).repeat_interleave(32, 1), y0b, None, u, 2, 7600.0)
    b = tru.rollout(tru.context([yh], [job]).repeat_interleave(32, 1), y0b, None, u, 2, 7600.0)
    check("mppi_truth with nominal parameters == mppi_nominal (exactly)", torch.equal(a[0], b[0]) and torch.equal(a[1], b[1]))
    c = tru.rollout(tru.context([yh], [{"domain": {"density_scale": 1.12, "lift_scale": 1.12, "mass_scale": 1.0}}])
                    .repeat_interleave(32, 1), y0b, None, u, 2, 7600.0)
    check("  ... and the truth parameters change the prediction", not torch.equal(a[0], c[0]))
    pr = PH.PriorAdapter()
    yp, gp = pr.rollout(pr.context([yh], [job]), y0b, None, u, 2, 7600.0)
    zero = V.PIProberV2(1, 8)
    with torch.no_grad():
        ref = V.prober_rollout_v2(zero, torch.randn(1, 32, 60, 8), y0b.float()[None], u.float()[None], 5.0, 2, 7600.0)
    check("mppi_prior == an untrained (zero-output) PIProberV2", torch.equal(yp[:, 0], ref[0][0]) and
          torch.equal(gp[:, 0], ref[1][0].norm(dim=-1) / D.G0))


# ============================================================================
# 4. CRN, 10. write guard
# ============================================================================
def t4_crn(L3: S.Layout, tmp: Path) -> None:
    print("4. common random numbers")
    cfg = S4.StudyConfigV4.smoke()
    L = S4.LayoutV4(tmp / "crn", read_only=[L3.root]).make()
    jobs = S4.build_jobs(cfg, L, L3)
    by = {}
    for j in jobs:
        if j["model"]["kind"] != "classical":
            by.setdefault((j["seed"], j["case"], j["atmosphere"]), set()).add(j["mppi"]["seed"])
    check("one noise stream per (seed, case, atmosphere) across every v4 arm and load", all(len(v) == 1 for v in by.values()),
          f"({len(by)} conditions)")
    v3 = {(j["seed"], j["case"], j["atmosphere"]): j["mppi"]["seed"] for j in S.build_jobs(S.StudyConfig.smoke(), L3)
          if j["controller"] == "v2" and j["mppi_rep"] == 0}
    phys = {(j["seed"], j["case"], j["atmosphere"]): j["mppi"]["seed"] for j in jobs if j["model"]["kind"] == "physics_truth"}
    check("physics noise index s == v3's v2 seed s, condition by condition", phys and all(phys[k] == v3[k] for k in phys))
    check("different conditions get different streams", len(set(phys.values())) == len(phys))


def t10_guard(L3: S.Layout, tmp: Path) -> None:
    print("10. the v4 code refuses to write under the v3 folder")
    L = S4.LayoutV4(tmp / "guard", read_only=[L3.root]).make()
    tries = {"v3 folder": L3.study / "x.json", "outside both": tmp / "elsewhere.json",
             "path escaping the v4 folder": L.root / ".." / "escape.json"}
    for what, p in tries.items():
        try:
            L.write_json(p, {"x": 1})
            ok = False
        except S4.ReadOnlyViolation:
            ok = True
        check(f"write into the {what} refused", ok and not Path(p).exists())
    try:
        S4.LayoutV4(L3.root / "inside_v3", read_only=[L3.root])
        ok = False
    except S4.ReadOnlyViolation:
        ok = True
    check("a v4 folder inside the v3 folder refused", ok)
    rec = {"result": {"flight_id": "../../../escape"}, "telemetry": {"t_s": np.zeros(1, np.float32)}}
    try:
        L.write_flight(rec)
        ok = False
    except S4.ReadOnlyViolation:
        ok = True
    check("a flight id escaping the results folder refused", ok)
    m3 = {p: p.stat().st_mtime for p in L3.root.rglob("*") if p.is_file()}
    check("v3 folder untouched so far (no file added or modified)", m3 == {p: p.stat().st_mtime for p in
                                                                         L3.root.rglob("*") if p.is_file()})


# ============================================================================
# Fixtures: data and the smoke v3 folder
# ============================================================================
def ensure_data(path: Path, workers: int) -> Path:
    if path.exists() and len(list(path.glob("shard_*.pkl"))) >= 4:
        return path
    print(f"collecting ~200 wide-envelope flights into {path} ...")
    D.collect_shards(200, str(path), shard_size=50, seed0=100_000, workers=workers, max_steps=1200,
                     spec=D.WIDE_DOMAIN, envelope=D.EnvelopeSpec())
    return path


def ensure_v3(root: Path, wide: Path, workers: int) -> S.Layout:
    """A smoke v3 study made by the v3 code (resumable: finished work is skipped)."""
    cfg = S.StudyConfig.smoke()
    L3 = S.Layout(root).make()
    S.freeze_config(cfg, L3)
    t0 = time.time()
    data = None
    for seed in cfg.seeds:
        if not (L3.ckpt("v1", seed).exists() and L3.ckpt("v2", seed).exists()):
            data = data or S.load_data(cfg, str(wide), "cpu")
            S.train_v2_seed(cfg, data, L3, seed, "cpu")
            S.train_v1_seed(cfg, data, L3, seed, "cpu")
    jobs = S.build_jobs(cfg, L3)
    S.run_classical(jobs, L3, workers=min(workers, 8))
    S.run_learned(cfg, jobs, L3, device="cpu", workers=min(4, workers))
    df = S.collect_flights(L3)
    print(f"smoke v3 folder ready: {len(df)} flights ({time.time() - t0:.0f}s)")
    return L3


# ============================================================================
# 2. E0c
# ============================================================================
def t2_e0c(L3: S.Layout, tmp: Path) -> None:
    print("2. run_lockstep_v4 == run_lockstep_v3 on the CPU (E0c)")
    cfg = S4.StudyConfigV4.smoke()
    L = S4.LayoutV4(tmp / "e0c", read_only=[L3.root]).make()
    S4.freeze_config_v4(cfg, L, L3)
    t0 = time.time()
    ok = S4.e0c_run(cfg, L, L3)
    df = pd.read_csv(L.checks / "E0c.csv")
    check("v1 and v2 adapters: every result field and telemetry array identical", ok,
          f"({len(df)} flights, {int(df.decisions.sum())} logged decisions, {time.time() - t0:.0f}s)")


# ============================================================================
# 5. Variants and resumable training
# ============================================================================
class _Stop(Exception):
    pass


def t5_variants(wide: Path, L3: S.Layout, tmp: Path) -> None:
    print("5. variants: parameters, build, train, save, interrupt, resume")
    full = S4.StudyConfigV4()
    for name, target in (("A7a_size_250k", 250_000), ("A7b_size_33k", 32_675)):
        n = AB.param_count(full, name)["per_member"]
        check(f"{name} at full size: {n:,} parameters", abs(n - target) / target <= 0.10,
              f"({100 * (n - target) / target:+.1f}% of {target:,})")
    check("B3 receptive field 25 (120 s)", AB.param_count(full, "B3_v1_plus_memory_120s")["rf"] == 25)
    check("A2 receptive field 7 (35 s)", AB.param_count(full, "A2_memory_35s")["rf"] == 7)
    cfg = S4.StudyConfigV4.smoke()
    cfg.save_every_s = -1.0                   # save after every step, so an interruption loses nothing
    data = AB.load_data_v4(cfg, str(wide), "cpu")
    ds7 = V.WindowsV2(data["train"], 5, data["ds_tr"].max_T, data["stats"])
    sl = AB.SlicedWindows(data["ds_tr"], 5)
    idx = torch.randint(len(sl), (64,), generator=torch.Generator().manual_seed(0))
    b1, b2 = sl.batch(idx, 8), ds7.batch(idx, 8)
    check("sliced windows == windows rebuilt with the smaller receptive field",
          len(sl) == len(ds7) and all(torch.equal(b1[k], b2[k]) for k in b1))
    L = S4.LayoutV4(tmp / "train", read_only=[L3.root]).make()
    real = {"wm": V.wm_loss_v2, "v1": M.jepa_loss}
    for name in cfg.variants:
        stop_after = {"n": 7}

        def flaky(fn):
            def f(*a, **k):
                stop_after["n"] -= 1
                if stop_after["n"] < 0:
                    raise _Stop()
                return fn(*a, **k)
            return f

        fam = AB.family(name)
        if fam == "v2":
            V.wm_loss_v2 = flaky(real["wm"])
        else:
            M.jepa_loss = flaky(real["v1"])
        try:
            AB.train_variant_seed(cfg, data, L, name, 0, "cpu")
            interrupted = False
        except _Stop:
            interrupted = True
        finally:
            V.wm_loss_v2, M.jepa_loss = real["wm"], real["v1"]
        st = torch.load(L.model_dir(name, 0) / "train_state_wm.pt", map_location="cpu", weights_only=False)
        step_saved = st["step"] if fam == "v2" else st["st"]["step"]
        # a torn save: the temp file of an interrupted write is left behind; the last complete state is untouched
        (L.model_dir(name, 0) / ".train_state_wm.pt.99999.tmp").write_bytes(b"torn")
        r = AB.train_variant_seed(cfg, data, L, name, 0, "cpu")
        info = json.loads(L.variant_json(name, 0).read_text())
        n_total = cfg.r0_wm["steps"] if fam == "v2" else cfg.v3_training().v1_steps
        ck = L.ckpt(name, 0)
        loads = F4.load_adapter_v4("v1" if fam == "v1" else "v2", str(ck), "cpu") is not None
        cnt = AB.param_count(cfg, name)["per_member"]
        check(f"{name}: interrupted at step {step_saved}, resumed, finished, saved, reloads",
              interrupted and step_saved == 7 and ck.exists() and loads and info["status"] == "done" and
              not (L.model_dir(name, 0) / "train_state_wm.pt").exists() and info["params"]["per_member"] == cnt,
              f"({cnt:,} parameters, {n_total} stage-1 steps, {r['train_s']}s)")
    # A resumed run continues the same random stream: identical to an uninterrupted run (CPU).
    for name in ("R0", "B1_v1_plus_inputs"):
        L2 = S4.LayoutV4(tmp / f"train_ref_{name}", read_only=[L3.root]).make()
        AB.train_variant_seed(cfg, data, L2, name, 0, "cpu")
        a, b = torch.load(L.ckpt(name, 0), weights_only=False), torch.load(L2.ckpt(name, 0), weights_only=False)
        same = all(torch.equal(a["wm"][1][k], b["wm"][1][k]) for k in a["wm"][1]) and \
            all(torch.equal(a["pi"][1][k], b["pi"][1][k]) for k in a["pi"][1])
        check(f"{name}: interrupted + resumed == uninterrupted, bit for bit", same)


# ============================================================================
# 6. Job queue
# ============================================================================
def t6_queue(L3: S.Layout, tmp: Path) -> None:
    print("6. job queue")
    cfg = S4.StudyConfigV4.smoke()
    L = S4.LayoutV4(tmp / "queue", read_only=[L3.root]).make()
    S4.freeze_config_v4(cfg, L, L3)
    jobs = S4.build_jobs(cfg, L, L3)
    m = [j for j in jobs if j["controller"] == "mppi_nominal" and j["load"] == "matched"]
    ok = all(abs(j["g_limit"] - json.loads((L3.flights / "results" / "classical" /
                                            f"{j['atmosphere']}_case{j['case']:02d}.json").read_text())["peak_g"]) == 0
             for j in m)
    check("matched limits == v3 classical peak g per condition", ok and len(m) > 0,
          f"({len(m)} jobs, limits {sorted(set(round(j['g_limit'], 2) for j in m))[:6]} ...)")
    two = [j for j in jobs if j["controller"] == "mppi_truth" and j["seed"] == 0][:2]
    S4.run_flights_v4(cfg, two, L, device="cpu", workers=1)
    check("both flights written", all(S4.flight_done(L, j) for j in two))
    first, second = L.result_path(F.flight_id(two[0])), L.result_path(F.flight_id(two[1]))
    t1 = first.stat().st_mtime_ns
    second.unlink()
    n = S4.run_flights_v4(cfg, two, L, device="cpu", workers=1)
    check("finished flight skipped, deleted flight re-flown", first.stat().st_mtime_ns == t1 and second.exists() and n == 1)
    tel = S.load_telemetry(L, F.flight_id(two[0]))
    check("decision log in the telemetry", all(f"dec_{k}" in tel for k in F4.DEC_FIELDS) and len(tel["dec_t_s"]) > 10,
          f"({len(tel['dec_t_s'])} decisions)")


# ============================================================================
# 9. The notebook; 7. analysis; 8. report
# ============================================================================
def t9_notebook(v4_root: Path, v3_root: Path, wide: Path) -> float:
    print("9. the notebook, executed in smoke mode (first run = the whole smoke study)")
    subprocess.run([sys.executable, "jepa/build_notebook_v4.py"], check=True)
    env = {**os.environ, "JEPA_V4_MODE": "smoke", "JEPA_V4_ROOT": str(v4_root), "JEPA_V3_ROOT": str(v3_root),
           "JEPA_WIDE_DATA": str(wide)}
    out = v4_root.parent / "nb_runs"
    out.mkdir(parents=True, exist_ok=True)
    times = []
    for k in (1, 2):
        t0 = time.time()
        r = subprocess.run([sys.executable, "-m", "nbconvert", "--to", "notebook", "--execute",
                            "--ExecutePreprocessor.timeout=-1", "--output-dir", str(out), "--output", f"run{k}.ipynb",
                            "ReEntryAI_JEPA_v4.ipynb"], env=env, capture_output=True, text=True)
        times.append(time.time() - t0)
        if r.returncode != 0:
            print(r.stdout[-3000:], r.stderr[-6000:])
        check(f"notebook run {k} ends without error", r.returncode == 0, f"({times[-1]:.0f}s)")
    check("second run finds everything done in under a minute", times[1] < 60.0, f"({times[1]:.0f}s)")
    nb = json.loads((out / "run2.ipynb").read_text(encoding="utf-8"))
    text = "".join("".join(o.get("text", "")) for c in nb["cells"] if c["cell_type"] == "code" for o in c.get("outputs", []))
    check("second run: no flight flown, no model trained", "flights in " not in text and "ms/step" not in text)
    return times[0]


def t7_analysis(v4_root: Path, v3_root: Path) -> None:
    print("7. analysis on the smoke data")
    from jepa import analyze_v3 as A3
    from jepa import analyze_v4 as A4
    rng = np.random.default_rng(0)
    XA, XB, XC = rng.random((5, 4, 64)), rng.random((5, 4, 64)), rng.random((1, 4, 64))
    XB[0, 1, 3] = np.nan
    bt = A4.Boot.__new__(A4.Boot)
    bt.n_boot = 2000
    same = all(abs(A3.hboot(X, Y, "mean", 2000, seed=17, atm=a)[k] - bt.run([X, Y], "diff_mean", a, seed=17)[k]) < 1e-12
               for X, Y in ((XA, XB), (XA, XC)) for a in (None, 2) for k in ("diff", "ci_lo", "ci_hi", "p_boot"))
    check("v4 bootstrap == analyze_v3.hboot (mean difference; with a deterministic arm; per atmosphere)", same)
    cfg = S4.StudyConfigV4.smoke()
    L3 = S4.v3_layout(v3_root)
    L = S4.LayoutV4(v4_root, read_only=[v3_root])
    need_t = [f"T{i:02d}" for i in range(1, 14)]
    need_f = [f"F{i}" for i in range(10, 19)]
    have_t = {p.stem.split("_")[0] for p in L.tables.glob("T*.csv")}
    have_f = {p.stem.split("_")[0] for p in L.figures.glob("F*.png")}
    check("every table T01-T13", set(need_t) <= have_t, f"(missing {sorted(set(need_t) - have_t)})")
    check("every figure F10-F18", set(need_f) <= have_f, f"(missing {sorted(set(need_f) - have_f)})")
    caps = json.loads((L.analysis / "captions.json").read_text())
    check("every table and figure has a caption and a data source",
          all(any(k.startswith(x) for k in caps) for x in need_t + need_f) and all(v.get("data") for v in caps.values()))
    t05 = pd.read_csv(next(L.tables.glob("T05_*.csv")))
    ids = ["P1", "P2", "P3", "P4", "P5", "P6", "L1", "L2", "L3", "R0-check", "M"] + [v for v in cfg.variants if v != "R0"]
    check("criteria table complete (P1-P6, L1-L3, 13 variants, R0-check, M)", set(ids) <= set(t05.id),
          f"({len(t05)} rows; not run: {list(t05[t05.measured == 'not run'].id)})")
    R = A4.analyze(cfg, L, L3, log=lambda *a: None)
    b = R["boot"]
    check("bootstrap cache reused on a rerun (no recomputation)", b.computed == 0 and b.reused > 0,
          f"({b.reused} reused)")
    fid = next(p for p in (L.flights / "results" / "mppi_truth").rglob("*.json"))
    orig = fid.read_bytes()
    rec = json.loads(orig)
    rec["miss_km"] = rec["miss_km"] + 3.0
    L.write_json(fid, rec)
    try:
        R = A4.analyze(cfg, L, L3, log=lambda *a: None)
        b2 = R["boot"]
        check("one changed result recomputes only the comparisons that use it", 0 < b2.computed < b.reused,
              f"({b2.computed} recomputed, {b2.reused} reused)")
    finally:
        L.atomic_bytes(fid, orig)
        A4.analyze(cfg, L, L3, log=lambda *a: None)


def t8_report(v4_root: Path, v3_root: Path) -> None:
    print("8. report")
    from jepa import make_report_v2 as R2
    from jepa import report_v4 as R4
    L = S4.LayoutV4(v4_root, read_only=[v3_root])
    pdf = L.report / "SpaceJEPA_v4_report.pdf"
    check("report written by the notebook", pdf.exists() and pdf.stat().st_size > 50_000, f"({pdf.stat().st_size // 1000} kB)")
    old = R2.FONT_DIR
    R2.FONT_DIR = str(v4_root / "no_system_fonts_here")
    try:
        p = R4.build_report(L, out_pdf=str(L.report / "smoke_fallback_fonts.pdf"))
    finally:
        R2.FONT_DIR = old
    check("report builds with matplotlib's DejaVu (no system fonts)", Path(p).stat().st_size > 50_000)
    Path(p).unlink()


def t_v3_notebooks(tmp: Path) -> None:
    print("+. v3 notebooks still build identically")
    for side, src in (("old", None), ("new", SIM_DIR / "jepa" / "build_notebooks_v3.py")):
        d = tmp / "nb" / side / "jepa"
        d.mkdir(parents=True, exist_ok=True)
        if src is None:
            txt = subprocess.run(["git", "-C", str(SIM_DIR), "show", "jepa-v3:StochasticEntrySim/jepa/build_notebooks_v3.py"],
                                 capture_output=True, text=True, encoding="utf-8", check=True).stdout
            (d / "build_notebooks_v3.py").write_text(txt, encoding="utf-8")
        else:
            shutil.copy(src, d / "build_notebooks_v3.py")
        subprocess.run([sys.executable, str(d / "build_notebooks_v3.py")], cwd=d.parent, check=True, capture_output=True)
    strip = lambda s: re.sub(r'"id": "[0-9a-f]+"', '"id": ""', s)
    names = sorted(p.name for p in (tmp / "nb" / "old").glob("*.ipynb"))
    same = [strip((tmp / "nb" / "old" / n).read_text(encoding="utf-8")) == strip((tmp / "nb" / "new" / n).read_text(encoding="utf-8"))
            for n in names]
    check(f"{len(names)} v3 notebooks byte-identical except nbformat's random cell ids", len(names) == 5 and all(same))
    import importlib
    sys.path.insert(0, str(SIM_DIR / "jepa"))
    before = {p: p.stat().st_mtime for p in SIM_DIR.glob("ReEntryAI_JEPA_v3_*.ipynb")}
    importlib.import_module("build_notebooks_v3")
    check("importing build_notebooks_v3 writes nothing", before == {p: p.stat().st_mtime for p in
                                                                    SIM_DIR.glob("ReEntryAI_JEPA_v3_*.ipynb")})


# ============================================================================
def main(work: Path, wide: Path, workers: int) -> None:
    t_all = time.time()
    torch.set_num_threads(max(1, (os.cpu_count() or 2) // 2))
    work.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix="smoke_v4_", dir=work))
    t1_defaults(tmp)
    t3_physics()
    wide = ensure_data(wide, workers)
    L3 = ensure_v3(work / "v3_smoke", wide, workers)
    t4_crn(L3, tmp)
    t10_guard(L3, tmp)
    t2_e0c(L3, tmp)
    t5_variants(wide, L3, tmp)
    t6_queue(L3, tmp)
    v4_root = work / "v4_smoke"
    if v4_root.exists():
        shutil.rmtree(v4_root)
    m3 = {p: p.stat().st_mtime for p in L3.root.rglob("*") if p.is_file()}
    t_nb = t9_notebook(v4_root, L3.root, wide)
    check("10. v3 folder unchanged by the whole v4 smoke study", m3 == {p: p.stat().st_mtime for p in
                                                                         L3.root.rglob("*") if p.is_file()})
    t7_analysis(v4_root, L3.root)
    t8_report(v4_root, L3.root)
    t_v3_notebooks(tmp)
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"\nv4 smoke test passed: {sum(ok for _, ok, _ in RESULTS)}/{len(RESULTS)} checks in "
          f"{(time.time() - t_all) / 60:.1f} min (notebook first run {t_nb / 60:.1f} min). Study folders in {work}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=os.environ.get("JEPA_SMOKE_WORK", str(Path(tempfile.gettempdir()) / "jepa_v4_smoke")))
    ap.add_argument("--wide", default=os.environ.get("JEPA_WIDE_DATA", ""),
                    help="wide-envelope shard directory (default: <work>/wide_smoke, collected if missing)")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()
    w = Path(a.work)
    main(w, Path(a.wide) if a.wide else w / "wide_smoke", a.workers)
