"""
jepa/size_sweep.py

Three JEPA world models of different size, each with one distinguishing idea,
trained on the same data with the same budget and compared open loop and in
closed loop against the classical predictor-corrector.

    A  paper-scale   SkyJEPA's own sizes (latent 16, GRU 24), ~4k params
    B  base          the notebook's model (latent 32, GRU 64), 26k params
    C  long-memory   75 s of history instead of 35 s, latent 64, GRU 128

Run from StochasticEntrySim/ (reuses jepa_runs/episodes_900.pkl and the
classical-result cache):

    python jepa/size_sweep.py [--cases 16]

Writes checkpoints, CSVs, and plots to jepa_runs/sweep/.
"""

from __future__ import annotations

import argparse
import pickle
import sys
import time
from pathlib import Path

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

import numpy as np
import torch

from jepa import data as D, evaluate as E, models as M, mppi as P

OUT = SIM_DIR / "jepa_runs" / "sweep"
EPISODES = SIM_DIR / "jepa_runs" / "episodes_900.pkl"
CLASSICAL_CACHE = SIM_DIR / "jepa_runs" / "classical_cache.pkl"

VARIANTS = {
    "A_paper_scale": dict(d_latent=16, z_dim=8, hidden=24, state_ch=(8, 16), act_ch=(4, 8)),
    "B_base": dict(d_latent=32, z_dim=16, hidden=64, state_ch=(32, 32), act_ch=(16, 16)),
    "C_long_memory": dict(d_latent=64, z_dim=32, hidden=128, state_ch=(64, 64, 64), act_ch=(32, 32, 32),
                          dilations=(1, 2, 4)),
}
# Tags match the classical cache, so cached classical runs are reused.
CONDITIONS = {
    "nominal": dict(D.NOMINAL_DOMAIN),
    "mismatch": {"density_scale": 0.88, "lift_scale": 0.88, "mass_scale": 1.0},
    "mismatch_high": {"density_scale": 1.12, "lift_scale": 1.12, "mass_scale": 1.0},
    "out_of_range": {"density_scale": 0.78, "lift_scale": 0.80, "mass_scale": 1.0},
}
STEPS_JEPA, STEPS_PI, BATCH, HORIZON, PROBER_HORIZON = 3000, 1500, 256, 24, 16
MACRO_DT, N_SUB, LAMBDA_SIG = D.MACRO_DT_S, 2, 0.02


def _splits():
    """Same macro episodes, split, and normalization stats for every variant."""
    with open(EPISODES, "rb") as f:
        eps = pickle.load(f)
    mac = [D.to_macro(e, MACRO_DT) for e in eps]
    tr, va, te = D.split_episodes(mac, (0.8, 0.1, 0.1), seed=0)
    return tr, va, te, D.compute_stats(tr)


def train_variant(job: dict) -> dict:
    """Train one variant's world model and PI prober (pool worker)."""
    torch.set_num_threads(1)
    torch.manual_seed(0)
    name, arch = job["name"], job["arch"]
    tr, va, _, stats = _splits()
    wm = M.JEPAWorldModel(**arch)
    sigreg = M.SIGReg()
    ds, dsv = D.MacroDataset(tr, wm.rf, HORIZON, stats), D.MacroDataset(va, wm.rf, HORIZON, stats)
    dsp, dspv = D.MacroDataset(tr, wm.rf, PROBER_HORIZON, stats), D.MacroDataset(va, wm.rf, PROBER_HORIZON, stats)
    vb = dsv.sample(512, torch.Generator().manual_seed(1))
    vbp = dspv.sample(512, torch.Generator().manual_seed(2))
    t0 = time.time()
    M.fit(wm.parameters(), lambda: M.jepa_loss(wm, sigreg, ds.sample(BATCH), LAMBDA_SIG), STEPS_JEPA,
          log_every=1000, label=f"{name} jepa")
    wm.eval()
    pi = M.PIProber(wm.d_latent)
    M.fit(pi.parameters(), lambda: M.prober_loss(wm, pi, dsp.sample(BATCH), MACRO_DT, N_SUB, w_acc=1.0),
          STEPS_PI, log_every=500, label=f"{name} prober")
    pi.eval()
    with torch.no_grad():
        val_pred = M.jepa_loss(wm, sigreg, vb, 0.0)[1]["pred"]
        val_state = M.prober_loss(wm, pi, vbp, MACRO_DT, N_SUB)[1]["state"]
    health = E.latent_health(wm, dsv)
    path = OUT / f"{name}.pt"
    M.save_checkpoint(str(path), wm, pi, None, None, stats, extra={"variant": name, "arch": arch})
    return {"variant": name, "params_wm": M.n_params(wm), "params_prober": M.n_params(pi),
            "history_s": wm.rf * MACRO_DT, "train_s": round(time.time() - t0),
            "val_latent_mse": round(val_pred, 4), "val_state_loss": round(val_state, 3),
            "eff_rank": round(health["eff_rank"], 1), "d_latent": wm.d_latent, "ckpt": str(path)}


def main(n_cases: int, workers: int) -> None:
    """Train, evaluate open loop, fly closed loop, and write the comparison."""
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    OUT.mkdir(parents=True, exist_ok=True)

    # 1. Train the three variants in parallel.
    t0 = time.time()
    train = D.run_parallel(train_variant, [{"name": n, "arch": a} for n, a in VARIANTS.items()],
                           min(workers, len(VARIANTS)), label="train")
    train_df = pd.DataFrame(train).set_index("variant")
    print(train_df.drop(columns="ckpt"), f"\ntraining: {time.time() - t0:.0f}s")
    train_df.to_csv(OUT / "training.csv")

    # 2. Open loop: footprint error vs time-to-go on held-out episodes.
    _, va, te, stats = _splits()
    nominal = E.nominal_params()
    rows = []
    for name in VARIANTS:
        m = M.load_checkpoint(train_df.loc[name, "ckpt"])
        rows += E.footprint_errors({name: {"kind": "jepa", "wm": m["wm"], "prober": m["pi"]}}, te + va, stats,
                                   m["wm"].rf, ttg_s=(400, 300, 200, 100, 50), h_cross=7600.0, dt=MACRO_DT,
                                   n_sub=N_SUB, nominal=nominal if name == "B_base" else None)
    fp = pd.DataFrame(rows)
    fp.to_csv(OUT / "footprint.csv", index=False)
    print("\nmedian footprint error [km]\n", fp.groupby(["model", "ttg_s"]).err_km.median().unstack(0).round(2))

    # 3. Closed loop: every variant on the same ICs and truths; classical reused from cache.
    ics = D.load_test_ics(str(SIM_DIR / "configs" / "heldout_testset.csv"), n=n_cases)
    cache = pickle.load(open(CLASSICAL_CACHE, "rb")) if CLASSICAL_CACHE.exists() else {}
    jobs = []
    for tag, dom in CONDITIONS.items():
        for i, ic in enumerate(ics):
            if (tag, i) not in cache:
                jobs.append({"controller": "classical", "ic": ic, "domain": dom, "tag": tag, "case": i,
                             "record": True, "config_path": str(SIM_DIR / "configs" / "default.json")})
            for name in VARIANTS:
                jobs.append({"controller": "jepa", "ic": ic, "domain": dom, "ckpt": train_df.loc[name, "ckpt"],
                             "mppi": {"n_sub": N_SUB}, "tag": tag, "case": i, "variant": name,
                             "config_path": str(SIM_DIR / "configs" / "default.json")})
    jobs.sort(key=lambda j: j["controller"] != "classical")
    t0 = time.time()
    new = D.run_parallel(_run_tagged, jobs, workers, label="closed loop")
    print(f"closed loop: {len(jobs)} episodes in {time.time() - t0:.0f}s")
    for r in new:
        if r["controller"] == "classical":
            cache[(r["tag"], r["case"])] = r
    pickle.dump(cache, open(CLASSICAL_CACHE, "wb"))
    results = [r for r in new if r["controller"] != "classical"]
    for tag in CONDITIONS:
        for i in range(n_cases):
            results.append({**cache[(tag, i)], "controller": "classical"})
    df = pd.DataFrame([{k: v for k, v in r.items() if k not in ("traj", "target")} for r in results])
    df.to_csv(OUT / "closed_loop.csv", index=False)
    g = df.groupby(["tag", "controller"])
    summary = pd.DataFrame({"success_%": 100 * g.success.mean(), "miss_med_km": g.miss_km.median(),
                            "miss_p90_km": g.miss_km.quantile(0.9), "peak_g": g.peak_g.mean(),
                            "fuel_kg": g.rcs_fuel_kg.mean(), "plan_ms": g.plan_ms_mean.mean()}).round(2)
    summary.to_csv(OUT / "summary.csv")
    print(summary)

    # 4. Plots: footprint vs time-to-go, and miss CDF per truth.
    col = {"A_paper_scale": "#1baf7a", "B_base": "#2a78d6", "C_long_memory": "#eb6834",
           "physics-nominal": "#4a3aa7", "classical": "#4a3aa7"}
    med = fp.groupby(["model", "ttg_s"]).err_km.median().unstack(0)
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    for n in med.columns:
        ax.plot(med.index, med[n], marker="o", ms=5, color=col[n], label=n)
    ax.set(xlabel="time-to-go when predicting [s]", ylabel="median footprint error [km]", yscale="log",
           title="Landing-point prediction by model size")
    ax.invert_xaxis()
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "footprint_by_size.png", dpi=200)
    fig, axes = plt.subplots(1, len(CONDITIONS), figsize=(16, 3.4), sharey=True)
    for ax, tag in zip(axes, CONDITIONS):
        for c in ["classical", *VARIANTS]:
            m = np.sort(df[(df.tag == tag) & (df.controller == c)].miss_km.values)
            ax.step(m, np.arange(1, len(m) + 1) / len(m), where="post", color=col[c], label=c)
        ax.axvline(25, color="#a3a29c", lw=1, ls=":")
        ax.set(xscale="log", xlabel="miss at drogue [km]", title=tag)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("fraction of cases")
    axes[0].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "miss_cdf_by_size.png", dpi=200)
    print(f"wrote {OUT}")


def _run_tagged(job: dict) -> dict:
    """run_episode, with the variant name used as the controller label."""
    r = P.run_episode(job)
    if job.get("variant"):
        r["controller"] = job["variant"]
    return r


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", type=int, default=16)
    ap.add_argument("--workers", type=int, default=11)
    a = ap.parse_args()
    main(a.cases, a.workers)
