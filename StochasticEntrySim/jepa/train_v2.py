"""
jepa/train_v2.py

Training orchestration for SpaceJEPA v2 (jepa/models_v2.py), kept out of the
notebook so it stays a thin configure-call-plot layer like v1's.

    find_wide_dataset     locate episodes_wide_6000.pkl (Drive, local, or the
                          GitHub release data-wide-6000)
    prepare               macro steps, episode split, v2 stats, GPU windows
    train_world_model     stage 1: ensemble JEPA, per-step SIGReg, curriculum
    train_prober          stage 2: state-aware PI probers on frozen latents
    lambda_sweep          one member per SIGReg weight, scored on val footprint
    train_v1_baseline     v1's architecture and recipe, retrained on the wide
                          data (separates "more data" from "new model")
"""

from __future__ import annotations

import os
import pickle
import subprocess
import time
from dataclasses import asdict, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import torch

from jepa import data as D
from jepa import models_v2 as V

RELEASE_TAG = "data-wide-6000"
RELEASE_FILE = "episodes_wide_6000.pkl"
REPO = "rk0604/ReEntryAI"


# ============================================================================
# Data
# ============================================================================
def _has_shards(path: str) -> bool:
    """True if path is a directory holding collect_shards output."""
    return bool(path) and os.path.isdir(path) and any(Path(path).glob("shard_*.pkl"))


def find_wide_dataset(candidates: Sequence[str], download_dir: str) -> str:
    """
    Return a path to the wide-envelope episodes: a pickle, or a shard directory.

    Tries each candidate in order (a merged episodes_wide_6000.pkl, or the
    shard directory ReEntryAI_JEPA_data.ipynb writes to Drive), then the
    GitHub release (gh CLI if installed and authenticated, else plain HTTPS,
    which works for a public repo), saving into download_dir.
    load_episodes() reads either kind.
    """
    for c in candidates:
        if c and (os.path.isfile(c) or _has_shards(c)):
            return c
    os.makedirs(download_dir, exist_ok=True)
    dst = os.path.join(download_dir, RELEASE_FILE)
    if os.path.exists(dst):
        return dst
    url = f"https://github.com/{REPO}/releases/download/{RELEASE_TAG}/{RELEASE_FILE}"
    cmds = [["gh", "release", "download", RELEASE_TAG, "-R", REPO, "-p", RELEASE_FILE, "-D", download_dir],
            ["wget", "-q", "-O", dst + ".part", url], ["curl", "-sSLf", "-o", dst + ".part", url]]
    for cmd in cmds:
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
        if os.path.exists(dst + ".part") and r.returncode == 0:
            os.replace(dst + ".part", dst)
        if os.path.exists(dst) and os.path.getsize(dst) > 1_000_000:
            return dst
        print(f"[data] {' '.join(cmd[:3])}: exit {r.returncode} {r.stderr.strip()[:200]}")
    raise FileNotFoundError(
        f"No wide-envelope data in {[c for c in candidates if c]} and the {RELEASE_TAG} release could not be "
        f"downloaded. Collect it with ReEntryAI_JEPA_data.ipynb (it writes shards and the merged file to Drive).")


def local_copy(path: str, local_dir: str) -> str:
    """Copy a (Drive-mounted) pickle to local disk once; faster to read, and pool workers reread it."""
    if os.path.isdir(path):
        return path
    os.makedirs(local_dir, exist_ok=True)
    dst = os.path.join(local_dir, os.path.basename(path))
    if not os.path.exists(dst) or os.path.getsize(dst) != os.path.getsize(path):
        import shutil
        shutil.copyfile(path, dst + ".part")
        os.replace(dst + ".part", dst)
    return dst


def load_episodes(path: str) -> List[dict]:
    """A pickle of collect_episode dicts, or a collect_shards directory."""
    if os.path.isdir(path):
        return D.load_shards(path)
    with open(path, "rb") as f:
        return pickle.load(f)


def prepare(episodes: Sequence[dict], max_T: int, device, rf: int, split_seed: int = 0,
            fracs=(0.8, 0.1, 0.1)) -> dict:
    """
    Macro-step the episodes, split by episode, fit v2 stats on train, and
    build the train/val window sets on the device.
    """
    macro = [D.to_macro(e, D.MACRO_DT_S) for e in episodes]
    tr, va, te = D.split_episodes(macro, fracs, seed=split_seed)
    stats = V.compute_stats_v2(tr)
    return {"train": tr, "val": va, "test": te, "stats": stats, "rf": rf,
            "ds_tr": V.WindowsV2(tr, rf, max_T, stats, device),
            "ds_va": V.WindowsV2(va, rf, max_T, stats, device)}


# ============================================================================
# Stage 1 and stage 2
# ============================================================================
def train_world_model(wm: V.WorldModelV2, data: dict, cfg: V.TrainCfgV2, lam, device,
                      ckpt_path: Optional[str] = None, sig_proj: int = 32, sig_steps: int = 16,
                      sig_min_rows: int = 32, label: str = "wm") -> Dict[str, list]:
    """
    Stage 1. lam: float, or a per-member list/tensor (lambda sweep).
    Validation: masked latent MSE per member on a fixed val batch at the
    current curriculum horizon (lambda = 0 so it is comparable across lambdas).
    """
    sig = V.SIGRegV2(n_proj=sig_proj).to(device)
    lam_t = torch.as_tensor(lam, dtype=torch.float32, device=device)
    g = torch.Generator(device=data["ds_va"].device).manual_seed(1)
    vbs = {T: data["ds_va"].sample((wm.M, 256), T, g) for _, T in cfg.curriculum}
    state = {"T": cfg.curriculum[0][1]}

    def loss_fn(step, T):
        state["T"] = T
        return V.wm_loss_v2(wm, sig, data["ds_tr"].sample((wm.M, cfg.batch), T), lam_t,
                            sig_steps=sig_steps, sig_min_rows=sig_min_rows)

    def val_fn():
        _, logs = V.wm_loss_v2(wm, sig, vbs[state["T"]], 0.0, sig_steps=sig_steps, sig_min_rows=sig_min_rows)
        return {"pred": logs["pred"]}

    torch.manual_seed(cfg.seed)
    return V.train_stage(wm, loss_fn, cfg, torch.device(device), ckpt_path, val_fn, label)


def train_prober(wm: V.WorldModelV2, pi: V.PIProberV2, data: dict, cfg: V.TrainCfgV2, device,
                 ckpt_path: Optional[str] = None, dt: float = D.MACRO_DT_S, n_sub: int = 2,
                 w_acc: float = 1.0, label: str = "prober") -> Dict[str, list]:
    """Stage 2 on frozen latents (the world model is put in eval mode and not optimized)."""
    wm.eval()
    for p in wm.parameters():
        p.requires_grad_(False)
    g = torch.Generator(device=data["ds_va"].device).manual_seed(2)
    vbs = {T: data["ds_va"].sample((pi.M, 256), T, g) for _, T in cfg.curriculum}
    state = {"T": cfg.curriculum[0][1]}

    def loss_fn(step, T):
        state["T"] = T
        return V.prober_loss_v2(wm, pi, data["ds_tr"].sample((pi.M, cfg.batch), T), dt, n_sub, w_acc)

    def val_fn():
        return {"state": V.prober_loss_v2(wm, pi, vbs[state["T"]], dt, n_sub, 0.0)[1]["state"]}

    torch.manual_seed(cfg.seed + 1)
    try:
        return V.train_stage(pi, loss_fn, cfg, torch.device(device), ckpt_path, val_fn, label)
    finally:
        for p in wm.parameters():
            p.requires_grad_(True)


@torch.no_grad()
def latent_health_v2(wm: V.WorldModelV2, ds: V.WindowsV2, n: int = 4096, seed: int = 0) -> List[dict]:
    """Per-member collapse check (as evaluate.latent_health): per-dim std, effective rank, |corr|."""
    g = torch.Generator(device=ds.device).manual_seed(seed)
    b = ds.sample((n,), 1, g)
    s = wm.encode_states(b["feats"][:, :wm.rf]).float()[:, :, 0]          # [M, n, D]
    out = []
    for m in range(wm.M):
        x = s[m] - s[m].mean(0)
        sv = torch.linalg.svdvals(x)
        p = sv / sv.sum()
        corr = torch.corrcoef(x.T)
        off = corr[~torch.eye(corr.shape[0], dtype=torch.bool, device=corr.device)]
        out.append({"member": m, "eff_rank": float(torch.exp(-(p * torch.log(p + 1e-12)).sum())),
                    "std_min": float(x.std(0).min()), "std_med": float(x.std(0).median()),
                    "std_max": float(x.std(0).max()), "mean_abs_corr": float(off.abs().mean())})
    return out


# ============================================================================
# SIGReg lambda sweep (one ensemble member per lambda, trained in one pass)
# ============================================================================
def lambda_sweep(data: dict, lams: Sequence[float], hparams: dict, wm_cfg: V.TrainCfgV2,
                 pi_cfg: V.TrainCfgV2, device, run_dir: str, ttg_s=(300, 200, 100),
                 h_cross: float = 7600.0, n_eval: int = 150) -> dict:
    """
    Train len(lams) members at once, member i with SIGReg weight lams[i], then
    their probers, and score each on validation footprint error (median over
    the given times-to-go, via evaluate.footprint_errors). The latent MSE is
    not comparable across lambdas (a more isotropic latent is harder to
    predict), so the downstream footprint is the selection metric.
    """
    import pandas as pd
    from jepa import evaluate as E
    os.makedirs(run_dir, exist_ok=True)
    hp = dict(hparams, n_members=len(lams))
    wm = V.WorldModelV2(**hp).to(device)
    pi = V.PIProberV2(len(lams), wm.d_latent).to(device)
    h1 = train_world_model(wm, data, wm_cfg, list(lams), device, os.path.join(run_dir, "sweep_wm.pt"),
                           label="sweep wm")
    h2 = train_prober(wm, pi, data, pi_cfg, device, os.path.join(run_dir, "sweep_pi.pt"), label="sweep prober")
    models = {"wm": wm, "pi": pi, "stats": data["stats"]}
    health = latent_health_v2(wm, data["ds_va"])
    val_eps = [e for e in data["val"] if e["outcome"] == "drogue"][:n_eval]
    rows = []
    for i, lam in enumerate(lams):
        fp = pd.DataFrame(E.footprint_errors({f"lam={lam:g}": V.as_v1(models, [i])}, val_eps, data["stats"],
                                             wm.rf, ttg_s=ttg_s, h_cross=h_cross))
        med = fp.groupby("ttg_s").err_km.median()
        rows.append({"lambda": lam, **{f"fp_med_km@{int(t)}s": med.get(float(t), np.nan) for t in ttg_s},
                     "score_km": float(np.exp(np.log(med.values).mean())),
                     "val_latent_mse": float(h1["val"][-1]["pred"]) if h1["val"] else np.nan,
                     "val_state": h2["val"][-1]["state"] if h2["val"] else np.nan,
                     "eff_rank": health[i]["eff_rank"], "std_med": health[i]["std_med"]})
    table = pd.DataFrame(rows).set_index("lambda")
    best = float(table["score_km"].idxmin())
    return {"table": table, "best": best, "hist_wm": h1, "hist_pi": h2, "models": models}


# ============================================================================
# v1 architecture retrained on the wide data (ablation: data vs model)
# ============================================================================
def train_v1_baseline(job: dict) -> dict:
    """
    v1's notebook recipe (models.JEPAWorldModel 26k params, T = 24, PI prober
    T = 16, lambda = 0.02, models.fit on CPU) on the v2 training split.
    Pool-worker safe: job = {"episodes_path", "ckpt", "steps_jepa", "steps_pi",
    "split_seed", "threads"}. Returns the training record.
    """
    from jepa import models as M
    torch.set_num_threads(int(job.get("threads", 1)))
    torch.manual_seed(0)
    if os.path.exists(job["ckpt"]):
        return {"ckpt": job["ckpt"], "cached": True}
    t0 = time.time()
    eps = load_episodes(job["episodes_path"])
    if job.get("n_episodes"):
        eps = eps[:int(job["n_episodes"])]
    macro = [D.to_macro(e, D.MACRO_DT_S) for e in eps]
    tr, va, _ = D.split_episodes(macro, (0.8, 0.1, 0.1), seed=int(job.get("split_seed", 0)))
    stats = D.compute_stats(tr)
    wm = M.JEPAWorldModel()
    sig = M.SIGReg()
    ds = D.MacroDataset(tr, wm.rf, 24, stats)
    dsp = D.MacroDataset(tr, wm.rf, 16, stats)
    M.fit(wm.parameters(), lambda: M.jepa_loss(wm, sig, ds.sample(256), 0.02), int(job["steps_jepa"]),
          log_every=max(1, int(job["steps_jepa"]) // 5), label="v1-wide jepa")
    wm.eval()
    pi = M.PIProber(wm.d_latent)
    M.fit(pi.parameters(), lambda: M.prober_loss(wm, pi, dsp.sample(256), D.MACRO_DT_S, 2, w_acc=1.0),
          int(job["steps_pi"]), log_every=max(1, int(job["steps_pi"]) // 5), label="v1-wide prober")
    pi.eval()
    tmp = job["ckpt"] + ".tmp"
    M.save_checkpoint(tmp, wm, pi, None, None, stats, extra={"variant": "v1-wide"})
    os.replace(tmp, job["ckpt"])
    return {"ckpt": job["ckpt"], "cached": False, "train_s": round(time.time() - t0)}
