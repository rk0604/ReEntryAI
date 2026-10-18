"""
jepa/smoke_v2.py

CPU smoke test of the whole v2 pipeline on a small data subset. It proves the
code runs end to end and checks the invariants the design relies on; the
numbers it prints are meaningless (a few hundred steps of a shrunken model).

    python jepa/smoke_v2.py [--episodes PATH_OR_SHARD_DIR] [--n 120] [--steps 300]

Run from StochasticEntrySim/. Checks:
  1. encoder: one pass over a long window == encoding each rf window alone
  2. predictor: KV-cached step-by-step == parallel pass
  3. per-step SIGReg == v1's SIGReg on the same rows and projections
  4. training: stage 1 (curriculum) and stage 2 run and reduce their losses;
     checkpoint save -> load round trip; resume from a stage checkpoint
  5. evaluate.footprint_errors (unchanged) runs on v1, v2 members, v2 ensemble
  6. the lockstep MPPI runner with one member and no uncertainty term flies
     the same commands as the unchanged mppi.run_episode on the same model
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
import time
from pathlib import Path

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))
os.chdir(SIM_DIR)

import numpy as np
import torch

from jepa import data as D, evaluate as E, models as M1, mppi as P
from jepa import models_v2 as V, mppi_v2 as P2, train_v2 as T2, evaluate_v2 as E2

TINY = dict(d_latent=16, enc_dim=32, enc_layers=2, enc_heads=2, enc_band=5, pred_dim=48, pred_layers=2,
            pred_heads=2, act_chunk=4, mlp_ratio=2)


def check(name: str, ok: bool, detail: str = "") -> None:
    """Print a PASS/FAIL line and stop on failure."""
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        raise SystemExit(1)


def load(src: str, n: int):
    """A shard directory (collect_shards output) or one pickle of episodes."""
    eps = D.load_shards(src) if os.path.isdir(src) else T2.load_episodes(src)
    return eps[:n]


def main(src: str, n: int, steps: int) -> None:
    torch.set_num_threads(max(1, (os.cpu_count() or 2) // 2))
    torch.manual_seed(0)
    dev = torch.device("cpu")
    tmp = tempfile.mkdtemp(prefix="jepa_v2_smoke_")
    t_all = time.time()

    print("1-3. invariants")
    wm = V.WorldModelV2(n_members=2, **TINY).double()
    torch.nn.init.normal_(wm.readout.weight, std=0.1)
    f = torch.randn(2, 3, 30, V.N_FEAT_V2, dtype=torch.double)
    S = wm.encode_states(f)
    S2 = torch.stack([wm.encode_states(f[:, :, j:j + wm.rf])[:, :, 0] for j in range(30 - wm.rf + 1)], 2)
    check("sliding encoder", (S - S2).abs().max() < 1e-8, f"(max diff {(S - S2).abs().max():.1e})")
    a = wm.action_tokens(torch.randn(2, 3, wm.rf - 1 + 10, dtype=torch.double))
    par = wm.predict(S[:, :, 0], a)
    cache, inc = wm.start(S[:, :, 0]), []
    for k in range(a.shape[2]):
        s, cache = wm.step_token(a[:, :, k], cache)
        inc.append(s)
    check("KV cache == parallel", (par - torch.stack(inc, 2)).abs().max() < 1e-10)
    x = torch.randn(512, 16)
    torch.manual_seed(5)
    v1s = M1.SIGReg()(x)
    torch.manual_seed(5)
    v2s = V.SIGRegV2()(x[None])[0][0]
    check("SIGReg v2 == v1", abs(float(v1s - v2s)) < 1e-4 * float(v1s), f"({float(v1s):.3f} vs {float(v2s):.3f})")

    print("4. training")
    eps = load(src, n)
    wm = V.WorldModelV2(n_members=3, **TINY)
    data = T2.prepare(eps, max_T=24, device=dev, rf=wm.rf)
    print(f"  {len(data['train'])}/{len(data['val'])}/{len(data['test'])} episodes, "
          f"{len(data['ds_tr'])} train windows, rf={wm.rf}")
    cfg = V.TrainCfgV2(steps=steps, batch=32, warmup=20, curriculum=((0.0, 8), (0.3, 16), (0.6, 24)),
                       log_every=max(1, steps // 4), amp=False)
    ck = os.path.join(tmp, "wm_stage.pt")
    vb = data["ds_va"].sample((wm.M, 128), 24, torch.Generator().manual_seed(9))
    sig = V.SIGRegV2()
    with torch.no_grad():
        l0 = V.wm_loss_v2(wm, sig, vb, 0.0)[1]["pred"]
    T2.train_world_model(wm, data, cfg, [0.02, 0.05, 0.1], dev, ck)
    with torch.no_grad():
        l1 = V.wm_loss_v2(wm, sig, vb, 0.0)[1]["pred"]
    check("stage 1 val latent MSE at T=24 fell", l1 < l0, f"({l0:.3f} -> {l1:.3f})")
    wm_b = V.WorldModelV2(n_members=3, **TINY)
    T2.train_world_model(wm_b, data, cfg, [0.02, 0.05, 0.1], dev, ck)   # resumes at the end: no steps run
    check("resume restores weights", all(torch.equal(p, q) for p, q in zip(wm.state_dict().values(),
                                                                         wm_b.state_dict().values())))
    pi = V.PIProberV2(wm.M, wm.d_latent, hidden=64, depth=2)
    pcfg = V.TrainCfgV2(steps=max(40, steps // 3), batch=32, warmup=10, lr=2e-3,
                        curriculum=((0.0, 8), (0.5, 16)), log_every=max(1, steps // 6), amp=False)
    vbp = data["ds_va"].sample((wm.M, 128), 16, torch.Generator().manual_seed(10))
    s0 = V.prober_loss_v2(wm, pi, vbp, D.MACRO_DT_S, 2, 0.0)[1]["state"]       # = physics prior only
    T2.train_prober(wm, pi, data, pcfg, dev)
    s1 = V.prober_loss_v2(wm, pi, vbp, D.MACRO_DT_S, 2, 0.0)[1]["state"]
    check("stage 2 val state loss at T=16 below physics prior only", s1 < s0, f"({s0:.2f} -> {s1:.2f})")
    print("  latent health:", [round(r["eff_rank"], 1) for r in T2.latent_health_v2(wm, data["ds_va"], n=512)])
    path = os.path.join(tmp, "v2.pt")
    V.save_checkpoint_v2(path, wm, pi, data["stats"], extra={"smoke": True})
    mv2 = V.load_checkpoint_v2(path)
    b = data["ds_va"].sample((4,), 8)
    with torch.no_grad():
        same = torch.allclose(wm.encode_states(b["feats"]), mv2["wm"].encode_states(b["feats"]))
    check("checkpoint round trip", same)
    print("  params:", V.param_report(mv2["wm"], mv2["pi"]))

    print("5. open loop through the unchanged evaluate.py")
    v1 = M1.load_checkpoint("jepa_runs/full/jepa_ckpt.pt")
    test = data["test"] + data["val"]
    rows = E.footprint_errors({"v1": {"kind": "jepa", "wm": v1["wm"], "prober": v1["pi"]}}, test, v1["stats"],
                              v1["wm"].rf, ttg_s=(200, 100), h_cross=7600.0)
    rows += E.footprint_errors({"v2 ens": V.as_v1(mv2), "v2 m0": V.as_v1(mv2, [0])}, test, mv2["stats"],
                               mv2["wm"].rf, ttg_s=(200, 100), h_cross=7600.0)
    import pandas as pd
    fp = pd.DataFrame(rows)
    print(fp.groupby(["model", "ttg_s"]).err_km.median().unstack(0).round(1).to_string())
    check("footprint rows for every model", fp.model.nunique() == 3)
    pts = pd.DataFrame(E2.footprint_points({"v2 ens": {**V.as_v1(mv2), "rf": mv2["wm"].rf}}, test, (200,), 7600.0))
    ref = fp[(fp.model == "v2 ens") & (fp.ttg_s == 200)].err_km.values
    check("evaluate_v2 matches evaluate.footprint_errors", np.allclose(np.sort(pts.err_km.values), np.sort(ref), atol=1e-3))
    sp = pd.DataFrame(E2.member_spread(mv2, test, (200,), 7600.0))
    print(f"  member spread median {sp.spread_km.median():.1f} km, ens-mean error median {sp.err_mean_km.median():.1f} km")
    sid = pd.DataFrame(E2.system_id_probe({"v2 m0": {"wm": V.as_v1(mv2, [0])["wm"], "stats": mv2["stats"], "rf": mv2["wm"].rf},
                                           "features": {"features": mv2["stats"]}},
                                          data["train"], test, max_eps=40))
    print(sid.round(2).to_string(index=False))
    rr = E.rollout_errors({"v2 ens": V.as_v1(mv2)}, test, mv2["stats"], mv2["wm"].rf, K=30, n_windows=64)
    print(f"  rollout pos RMSE at 150 s: {rr['v2 ens']['pos_km'][-1]:.1f} km")

    print("6. closed loop: lockstep runner vs the unchanged mppi.run_episode (1 member, no uncertainty)")
    one = V.select_members(mv2, [0])
    ic = D.load_test_ics("configs/heldout_testset.csv", n=1)[0]
    cfg_m = {"n_samples": 64, "n_sub": 2}
    P._MODEL_CACHE["smoke-v2"] = V.as_v1(one)
    t0 = time.time()
    ref = P.run_episode({"controller": "jepa", "ic": ic, "domain": dict(D.NOMINAL_DOMAIN), "ckpt": "smoke-v2",
                         "mppi": cfg_m, "tag": "nominal", "case": 0, "record": True})
    t_ref = time.time() - t0
    t0 = time.time()
    got = P2.run_lockstep([{"ic": ic, "domain": dict(D.NOMINAL_DOMAIN), "tag": "nominal", "case": 0, "record": True}],
                          one, {**cfg_m, "w_unc_km": 0.0, "trim_horizon": False}, verbose=False)[0]
    t_got = time.time() - t0
    n_cmp = min(len(ref["traj"]), len(got["traj"]), 20)
    dcmd = np.abs(ref["traj"][:n_cmp, 10] - got["traj"][:n_cmp, 10]).max()
    print(f"  run_episode: {ref['outcome']} miss {ref['miss_km']:.1f} km ({t_ref:.0f}s) | "
          f"lockstep: {got['outcome']} miss {got['miss_km']:.1f} km ({t_got:.0f}s) | "
          f"max |d cmd| over first {n_cmp} decisions {np.degrees(dcmd):.3f} deg")
    check("same first decisions", np.degrees(dcmd) < 0.5)
    # Full ensemble, uncertainty on, two truths in one lockstep batch, plus v1 via run_episode.
    jobs = [{"ic": ic, "domain": dom, "tag": tag, "case": 0} for tag, dom in
            (("nominal", dict(D.NOMINAL_DOMAIN)), ("mismatch_high", {"density_scale": 1.12, "lift_scale": 1.12,
                                                                    "mass_scale": 1.0}))]
    res = P2.run_lockstep(jobs, mv2, {"n_samples": 64, "n_sub": 2, "w_unc_km": 1.0}, verbose=False)
    res.append(P.run_episode({"controller": "jepa", "ic": ic, "domain": dict(D.NOMINAL_DOMAIN),
                              "ckpt": "jepa_runs/full/jepa_ckpt.pt", "mppi": {"n_samples": 64, "n_sub": 2},
                              "tag": "nominal", "case": 0}))
    print(P.summarize(res)[["n", "success_%", "miss_med_km", "peak_g_mean", "plan_ms"]].to_string())
    print(f"smoke test passed in {time.time() - t_all:.0f}s (scratch: {tmp})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--episodes", default="jepa_runs/episodes_900.pkl",
                    help="pickle of episodes or a collect_shards directory")
    ap.add_argument("--n", type=int, default=120)
    ap.add_argument("--steps", type=int, default=300)
    a = ap.parse_args()
    main(a.episodes, a.n, a.steps)
