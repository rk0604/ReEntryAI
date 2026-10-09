"""
jepa/evaluate.py

Open-loop evaluation of the learned models (SkyJEPA Sec. VI-A/B, adapted).

    rollout_errors     RMSE of altitude, speed, gamma, heading, and position vs.
                       unroll length, for JEPA+PI, JEPA+MLP, and Predictive
                       (paper Fig. 8 / Table III), plus per-step error growth
                       (paper Eq. 30).
    footprint_errors   the guidance-relevant metric: given the state at some
                       time-to-go and the commands actually flown afterwards,
                       how far is the predicted drogue-deploy point from the
                       true one? Includes the nominal-physics predictor (what
                       the classical predictor-corrector's model would say:
                       nominal aero, instantaneous bank) as the reference.
    latent_health      collapse diagnostics (per-dim std, effective rank).
    straightness       temporal straightening of a latent rollout (Sec. VI-A).
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence

import numpy as np
import torch

import constants
from jepa.data import MacroDataset, Stats, R_EARTH
from jepa.models import jepa_predict, predictive_predict


def haversine_km(phi1, lam1, phi2, lam2):
    """Great-circle distance in km; works on floats, numpy arrays, or tensors."""
    lib = torch if isinstance(phi1, torch.Tensor) else np
    a = lib.sin((phi2 - phi1) / 2) ** 2 + lib.cos(phi1) * lib.cos(phi2) * lib.sin((lam2 - lam1) / 2) ** 2
    if lib is torch:
        a = a.clamp(0.0, 1.0)
    else:
        a = np.clip(a, 0.0, 1.0)
    return 2.0 * R_EARTH / 1000.0 * lib.arcsin(lib.sqrt(a))


def crossing_point(y6: torch.Tensor, h_cross: float):
    """
    First downward crossing of h_cross along each rollout, linearly interpolated.

    y6 [B, K+1, 6]. Returns (phi [B], lam [B], crossed [B] bool, k_cross [B]
    long, the index of the first step at or below h_cross; K if never). Samples
    that never cross report their final position.
    """
    h = y6[..., 0]
    below = h <= h_cross
    crossed = below.any(1)
    K = y6.shape[1] - 1
    k = torch.where(crossed, below.float().argmax(1), torch.full_like(crossed, K, dtype=torch.long))
    k_prev = (k - 1).clamp_min(0)
    ar = torch.arange(y6.shape[0])
    y1, y0 = y6[ar, k], y6[ar, k_prev]
    dh = (y0[:, 0] - y1[:, 0]).clamp_min(1e-6)
    f = torch.where(crossed & (k > 0), ((y0[:, 0] - h_cross) / dh).clamp(0, 1), torch.ones_like(dh))
    phi = y0[:, 1] + f * (y1[:, 1] - y0[:, 1])
    lam = y0[:, 2] + f * (y1[:, 2] - y0[:, 2])
    return phi, lam, crossed, k


# ============================================================================
# Rollout error vs. unroll length
# ============================================================================
@torch.no_grad()
def rollout_errors(models: Dict[str, dict], test_eps: Sequence[dict], stats: Stats, rf: int,
                   K: int, n_windows: int = 256, dt: float = 5.0, n_sub: int = 2,
                   seed: int = 0) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Open-loop RMSE per unroll step k = 1..K on held-out episodes.

    models maps a label to {"kind": "jepa", "wm", "prober"} or
    {"kind": "pred", "pm"}. Returns label -> {quantity: [K] RMSE}, with
    quantities h_km, V_mps, gamma_deg, chi_deg, pos_km, plus err_growth_pos_km
    (mean increase in position error per step).
    """
    ds = MacroDataset(test_eps, rf=rf, horizon=K, stats=stats)
    g = torch.Generator().manual_seed(seed)
    idx = torch.randperm(len(ds), generator=g)[:n_windows]
    b = ds.batch(idx)
    P = rf - 1
    y = b["y"]
    y_ctx = ds.Y[ds.y_start[idx, None] + torch.arange(P + 1)]
    out = {}
    for name, m in models.items():
        if m["kind"] == "jepa":
            yhat = jepa_predict(m["wm"], m["prober"], b["feats"][:, :P + 1], y[:, 0, :6],
                                b["u"][:, :P], b["u"][:, P:], dt, n_sub)["y6"]
        else:
            yhat = predictive_predict(m["pm"], stats, y_ctx, b["u"][:, :P], b["u"][:, P:])["y6"]
        out[name] = _error_curves(yhat[:, 1:], y[:, 1:, :6])
    return out


def _error_curves(yhat6: torch.Tensor, y6: torch.Tensor) -> Dict[str, np.ndarray]:
    """RMSE over the batch per step, in physical units."""
    dn = (yhat6[..., 1] - y6[..., 1]) * R_EARTH / 1000.0
    de = (yhat6[..., 2] - y6[..., 2]) * R_EARTH * torch.cos(y6[..., 1]) / 1000.0
    pos = torch.sqrt(dn ** 2 + de ** 2)
    dchi = torch.atan2(torch.sin(yhat6[..., 5] - y6[..., 5]), torch.cos(yhat6[..., 5] - y6[..., 5]))
    rms = lambda e: torch.sqrt((e ** 2).mean(0)).numpy()
    return {"h_km": rms((yhat6[..., 0] - y6[..., 0]) / 1000.0),
            "V_mps": rms(yhat6[..., 3] - y6[..., 3]),
            "gamma_deg": rms(torch.rad2deg(yhat6[..., 4] - y6[..., 4])),
            "chi_deg": rms(torch.rad2deg(dchi)),
            "pos_km": rms(pos),
            "err_growth_pos_km": torch.diff(pos, dim=1, prepend=torch.zeros_like(pos[:, :1])).mean(0).numpy()}


# ============================================================================
# Footprint (drogue-deploy point) prediction
# ============================================================================
def nominal_params(config_path: str = "configs/default.json") -> dict:
    """The nominal aero/mass params the env and classical guidance start from."""
    import AtmosphereModel
    import mission_config
    if not AtmosphereModel.atmosphere_table_active():
        AtmosphereModel.build_atmosphere_table()
    p = {"mass_kg": float(constants.CAPSULE_MASS_KG),
         "ref_area_m2": float(constants.CAPSULE_REFERENCE_AREA_M2),
         "aero_model": "orion_cm_trim", "CD_const": 1.15, "CL_const": 0.28}
    mission_config.apply_aero_to_params(mission_config.load_config(config_path), p)
    return p


def physics_nominal_crossing(y0: np.ndarray, u_macro: np.ndarray, params: dict, h_cross: float,
                             macro_dt: float = 5.0, dt: float = 1.0, max_s: float = 1500.0):
    """
    Nominal-model footprint prediction: the classical predictor's view.

    Integrates math_3d.nominal_eom_step with nominal params and the commanded
    bank applied instantly (no actuator lag), holding the last command past
    the end of u_macro. Returns (phi, lam) at the h_cross crossing.
    """
    from math_3d import nominal_eom_step
    x = [float(y0[0]) + R_EARTH, float(y0[1]), float(y0[2]), float(y0[3]), float(y0[4]), float(y0[5])]
    n_per = int(round(macro_dt / dt))
    for i in range(int(max_s / dt)):
        sigma = float(u_macro[min(i // n_per, len(u_macro) - 1)])
        x_new = nominal_eom_step(x=x, sigma_rad=sigma, params=params, dt_s=dt)["x_next"]
        h0, h1 = x[0] - R_EARTH, x_new[0] - R_EARTH
        if h1 <= h_cross:
            f = (h0 - h_cross) / max(h0 - h1, 1e-6)
            return x[1] + f * (x_new[1] - x[1]), x[2] + f * (x_new[2] - x[2])
        x = x_new
    return x[1], x[2]


@torch.no_grad()
def footprint_errors(models: Dict[str, dict], test_eps: Sequence[dict], stats: Stats, rf: int,
                     ttg_s: Sequence[float], h_cross: float, dt: float = 5.0, n_sub: int = 2,
                     nominal: Optional[dict] = None, margin_steps: int = 30) -> List[dict]:
    """
    Predicted vs. true drogue-deploy point, from several times-to-go.

    For each held-out episode that reached drogue deploy and each time-to-go,
    every model is given the true history up to that point and the commands
    actually flown afterwards (last command held past the end). nominal, when
    given, adds the nominal-physics predictor under the label "physics-nominal".
    Returns rows {model, ttg_s, err_km, domain...}.
    """
    P = rf - 1
    eps = [e for e in test_eps if e["outcome"] == "drogue"]
    rows: List[dict] = []
    for ttg in ttg_s:
        k_back = int(round(ttg / dt))
        sel = [e for e in eps if len(e["u"]) - k_back >= 0]
        if not sel:
            continue
        K = k_back + margin_steps
        feats, y_ctx, y0, u_hist, u_plan, truth = [], [], [], [], [], []
        for e in sel:
            t = len(e["u"]) - k_back
            y = torch.from_numpy(e["y"]).float()
            u = torch.from_numpy(e["u"]).float()
            yp = torch.cat([y[:1].repeat(P, 1), y])           # front pad like MacroDataset
            up = torch.cat([torch.zeros(P), u])
            hist = yp[t:t + P + 1]
            y_ctx.append(hist)
            feats.append(stats.norm_feats(hist))
            y0.append(y[t, :6])
            u_hist.append(up[t:t + P])
            fut = u[t:]
            fut = torch.cat([fut, (fut[-1:] if len(fut) else u[-1:]).repeat(K - len(fut))])
            u_plan.append(fut[:K])
            truth.append(e["y_final"])
        feats, y_ctx, y0 = torch.stack(feats), torch.stack(y_ctx), torch.stack(y0)
        u_hist, u_plan = torch.stack(u_hist), torch.stack(u_plan)
        truth = np.stack(truth).astype(np.float64)
        preds = {}
        for name, m in models.items():
            if m["kind"] == "jepa":
                y6 = jepa_predict(m["wm"], m["prober"], feats, y0, u_hist, u_plan, dt, n_sub)["y6"]
            else:
                y6 = predictive_predict(m["pm"], stats, y_ctx, u_hist, u_plan)["y6"]
            phi, lam, _, _ = crossing_point(y6.double(), h_cross)
            preds[name] = (phi.numpy(), lam.numpy())
        if nominal is not None:
            pl = [physics_nominal_crossing(y0[i].double().numpy(), u_plan[i].numpy(), nominal, h_cross, dt)
                  for i in range(len(sel))]
            preds["physics-nominal"] = (np.array([p[0] for p in pl]), np.array([p[1] for p in pl]))
        for name, (phi, lam) in preds.items():
            err = haversine_km(phi, lam, truth[:, 1], truth[:, 2])
            for i, e in enumerate(sel):
                rows.append({"model": name, "ttg_s": float(ttg), "err_km": float(err[i]),
                             **{f"dom_{k}": v for k, v in e["domain"].items()}})
    return rows


# ============================================================================
# Latent diagnostics
# ============================================================================
@torch.no_grad()
def latent_health(wm, ds: MacroDataset, n: int = 4096, seed: int = 0) -> Dict[str, object]:
    """
    Collapse check on encoded latents.

    Returns per-dim std (SIGReg targets 1.0), the effective rank (exp of the
    entropy of the normalized singular values; D means fully spread, 1 means
    collapsed to a line), and the mean absolute off-diagonal correlation.
    """
    g = torch.Generator().manual_seed(seed)
    b = ds.batch(torch.randint(len(ds), (n,), generator=g))
    s = wm.encode_states(b["feats"])[:, 0]
    sv = torch.linalg.svdvals(s - s.mean(0))
    p = sv / sv.sum()
    corr = torch.corrcoef(s.T)
    off = corr[~torch.eye(corr.shape[0], dtype=torch.bool)]
    return {"std": s.std(0).numpy(), "eff_rank": float(torch.exp(-(p * torch.log(p + 1e-12)).sum())),
            "mean_abs_corr": float(off.abs().mean()), "latents": s.numpy()}


def straightness(seq: torch.Tensor) -> np.ndarray:
    """Mean cosine similarity of consecutive steps of seq [B, K, D] (1 = straight line)."""
    v = torch.diff(seq, dim=1)
    cos = torch.nn.functional.cosine_similarity(v[:, 1:], v[:, :-1], dim=-1)
    return cos.mean(1).numpy()
