"""
jepa/evaluate_v2.py

Diagnostics that evaluate.py does not have, aimed at v2's questions. The
headline open-loop metric stays evaluate.footprint_errors (unchanged); these
explain it:

    footprint_points   predicted vs true drogue point per (episode, time-to-go),
                       split into along-track / cross-track components. A
                       systematic along-track bias that tracks the true density
                       is the v1 dense-air failure.
    member_spread      per-member landing points of a v2 ensemble: spread
                       (the MPPI uncertainty signal) and its calibration
                       against the actual error.
    system_id_probe    how well a linear readout of the latent recovers the
                       truth perturbation (density/mass, L/D), by time since
                       entry: implicit system ID, the claim behind robustness.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np
import torch

from jepa.data import MACRO_DT_S, R_EARTH
from jepa.evaluate import crossing_point, haversine_km
from jepa.models import jepa_predict
from jepa import models_v2 as V


def _windows(eps: Sequence[dict], rf: int, stats, ttg: float, dt: float, margin: int):
    """Inputs at time-to-go ttg, built exactly as evaluate.footprint_errors builds them."""
    P = rf - 1
    k_back = int(round(ttg / dt))
    sel = [e for e in eps if e["outcome"] == "drogue" and len(e["u"]) - k_back >= 0]
    if not sel:
        return None
    K = k_back + margin
    feats, y0, u_hist, u_plan, truth = [], [], [], [], []
    for e in sel:
        t = len(e["u"]) - k_back
        y = torch.from_numpy(e["y"]).float()
        u = torch.from_numpy(e["u"]).float()
        yp = torch.cat([y[:1].repeat(P, 1), y])
        up = torch.cat([torch.zeros(P), u])
        feats.append(stats.norm_feats(yp[t:t + P + 1]))
        y0.append(y[t, :6])
        u_hist.append(up[t:t + P])
        fut = u[t:]
        fut = torch.cat([fut, (fut[-1:] if len(fut) else u[-1:]).repeat(K - len(fut))])
        u_plan.append(fut[:K])
        truth.append(e["y_final"])
    return sel, torch.stack(feats), torch.stack(y0), torch.stack(u_hist), torch.stack(u_plan), \
        np.stack(truth).astype(np.float64)


def _decompose(phi, lam, truth):
    """Error of (phi, lam) vs truth rows, in km: total, along-track, cross-track (true heading)."""
    phi_t, lam_t, chi = truth[:, 1], truth[:, 2], truth[:, 5]
    dn = (phi - phi_t) * R_EARTH / 1000.0
    de = (lam - lam_t) * R_EARTH * np.cos(phi_t) / 1000.0
    # chi is measured from east toward north (math_3d convention).
    along = de * np.cos(chi) + dn * np.sin(chi)
    cross = -de * np.sin(chi) + dn * np.cos(chi)
    return haversine_km(phi, lam, phi_t, lam_t), along, cross


@torch.no_grad()
def footprint_points(models: Dict[str, dict], eps: Sequence[dict], ttg_s: Sequence[float], h_cross: float,
                     dt: float = MACRO_DT_S, n_sub: int = 2, margin_steps: int = 30) -> List[dict]:
    """
    Signed footprint errors. models: label -> {"wm", "prober", "stats", "rf"}
    in the v1 interface (a v1 checkpoint, or models_v2.as_v1(...) plus its
    stats and rf). Rows: model, ttg_s, err_km, along_km (+ = predicted beyond
    the true point), cross_km, and the episode's truth factors.
    """
    rows = []
    for name, m in models.items():
        for ttg in ttg_s:
            w = _windows(eps, m["rf"], m["stats"], ttg, dt, margin_steps)
            if w is None:
                continue
            sel, feats, y0, uh, up, truth = w
            y6 = jepa_predict(m["wm"], m["prober"], feats, y0, uh, up, dt, n_sub)["y6"]
            phi, lam, _, _ = crossing_point(y6.double().cpu(), h_cross)
            err, along, cross = _decompose(phi.numpy(), lam.numpy(), truth)
            for i, e in enumerate(sel):
                rows.append({"model": name, "ttg_s": float(ttg), "err_km": float(err[i]),
                             "along_km": float(along[i]), "cross_km": float(cross[i]),
                             **{f"dom_{k}": v for k, v in e["domain"].items()}})
    return rows


@torch.no_grad()
def member_spread(models_v2: dict, eps: Sequence[dict], ttg_s: Sequence[float], h_cross: float,
                  dt: float = MACRO_DT_S, n_sub: int = 2, margin_steps: int = 30) -> List[dict]:
    """
    Every member's landing point through models_v2.ensemble_rollout. Rows: ttg_s,
    spread_km (rms distance of members from their mean point), err_mean_km
    (error of the mean point), and err_m{i}_km per member.
    """
    wm, pi, stats = models_v2["wm"], models_v2["pi"], models_v2["stats"]
    dev = next(wm.parameters()).device
    rows = []
    for ttg in ttg_s:
        w = _windows(eps, wm.rf, stats, ttg, dt, margin_steps)
        if w is None:
            continue
        sel, feats, y0, uh, up, truth = w
        s0 = wm.encode_states(feats.to(dev))[:, :, -1].float()
        out = V.ensemble_rollout(wm, pi, s0, y0.to(dev), uh.to(dev), up.to(dev), dt, n_sub)
        y6 = out["y6"].double().cpu()                                         # [B, M, K+1, 6]
        B, Mm = y6.shape[:2]
        phi, lam, _, _ = crossing_point(y6.flatten(0, 1), h_cross)
        phi, lam = phi.view(B, Mm).numpy(), lam.view(B, Mm).numpy()
        pm, lm = phi.mean(1), lam.mean(1)
        spread = np.sqrt((haversine_km(phi, lam, pm[:, None], lm[:, None]) ** 2).mean(1))
        err_mean, along, _ = _decompose(pm, lm, truth)
        err_m = haversine_km(phi, lam, truth[:, 1:2], truth[:, 2:3])
        for i, e in enumerate(sel):
            rows.append({"ttg_s": float(ttg), "spread_km": float(spread[i]), "err_mean_km": float(err_mean[i]),
                         "along_mean_km": float(along[i]),
                         **{f"err_m{j}_km": float(err_m[i, j]) for j in range(Mm)},
                         **{f"dom_{k}": v for k, v in e["domain"].items()}})
    return rows


# ============================================================================
# Implicit system ID
# ============================================================================
def _truth_targets(ep: dict) -> np.ndarray:
    """What is identifiable from the forces: log(density/mass) (ballistic coefficient) and log(L/D)."""
    d = ep["domain"]
    return np.array([np.log(d["density_scale"] / d["mass_scale"]), np.log(d["lift_scale"])])


@torch.no_grad()
def _latents_along(wm, stats, rf: int, ep: dict) -> np.ndarray:
    """Encoded latent at every macro step of an episode [L+1, D] (v1 interface, one pass)."""
    y = torch.from_numpy(ep["y"]).float()
    yp = torch.cat([y[:1].repeat(rf - 1, 1), y])
    return wm.encode_states(stats.norm_feats(yp)[None])[0].float().cpu().numpy()


def system_id_probe(models: Dict[str, dict], train_eps: Sequence[dict], test_eps: Sequence[dict],
                    bins_s=(0, 60, 180, 1e9), ridge: float = 1e-2, max_eps: int = 600) -> List[dict]:
    """
    Ridge regression from the encoded latent at time t to the episode's truth
    factors, fit on train episodes, R^2 on test episodes, per time-since-entry
    bin. models: label -> {"wm", "stats", "rf"} (v1 interface), or
    {"features": stats} for the raw normalized features at t (a no-memory reference).
    """
    def collect(eps, m):
        X, Y, T = [], [], []
        for ep in eps[:max_eps]:
            if "features" in m:
                lat = m["features"].norm_feats(torch.from_numpy(ep["y"]).float()).numpy()
            else:
                lat = _latents_along(m["wm"], m["stats"], m["rf"], ep)
            X.append(lat)
            Y.append(np.repeat(_truth_targets(ep)[None], len(lat), 0))
            T.append(np.arange(len(lat)) * MACRO_DT_S)
        return np.concatenate(X), np.concatenate(Y), np.concatenate(T)

    rows = []
    for name, m in models.items():
        Xtr, Ytr, _ = collect(train_eps, m)
        Xte, Yte, Tte = collect(test_eps, m)
        mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-6
        A = np.c_[(Xtr - mu) / sd, np.ones(len(Xtr))]
        W = np.linalg.solve(A.T @ A + ridge * len(A) * np.eye(A.shape[1]), A.T @ Ytr)
        pred = np.c_[(Xte - mu) / sd, np.ones(len(Xte))] @ W
        for lo, hi in zip(bins_s[:-1], bins_s[1:]):
            sel = (Tte >= lo) & (Tte < hi)
            if sel.sum() < 10:
                continue
            ss_res = ((pred[sel] - Yte[sel]) ** 2).sum(0)
            ss_tot = ((Yte[sel] - Yte[sel].mean(0)) ** 2).sum(0)
            r2 = 1.0 - ss_res / np.maximum(ss_tot, 1e-12)
            rows.append({"model": name, "t_since_entry": f"{lo:.0f}-{hi:.0f} s" if hi < 1e8 else f">{lo:.0f} s",
                         "R2_density/mass": float(r2[0]), "R2_L/D": float(r2[1]), "n": int(sel.sum())})
    return rows
