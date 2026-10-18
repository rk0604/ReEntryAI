"""
jepa/models_v2.py

SpaceJEPA v2: an ensemble of attention-based JEPA world models for the wide
entry envelope. Design and rationale: jepa/PLAN.md section 11. v1
(jepa/models.py) is untouched; v2 reuses its integrator (integrate_eom), its
dimensionless state errors, and its nominal aero constants.

What changes from v1, and why (one line each, the PLAN has the argument):

    inputs      + physics-residual features: the drag and L/D corrections the
                  prober needs, *measured* from the IMU at the current state
                  (log a_D / a_D,nominal(h, V), log (L/D) / 0.30). System ID
                  becomes reading a feature instead of inferring a product of
                  normalized inputs. Targets the dense-air bias.
    encoder     banded causal transformer, rf = 25 macro steps (120 s, v1: 35 s).
                  Exact receptive field and relative (RoPE) positions, so one
                  pass over a window still yields every sliding encoding, like
                  v1's valid-padding TCN.
    actions     one token per macro step built from the last c = 4 commands
                  (multi-step action encoding: the actuator lag lives here).
    predictor   causal transformer over [context latent, action tokens]
                  instead of a GRU. The latent at step k only sees s_t and
                  commands up to k, exactly the GRU's causal structure, but
                  training is one parallel pass, so a 96-step horizon is cheap.
                  Planning runs it in parallel too (one pass for all samples);
                  a KV cache gives the step-by-step v1 API for evaluate.py.
    loss        masked multi-step latent MSE + per-step SIGReg (the paper's
                  form) on encoded and predicted latents, horizon curriculum
                  24 -> 48 -> 96. Windows run past the end of short flights with
                  a mask, so short flights are not dropped from long horizons.
    prober      the v1 physics-inspired prober (corrections to a crude
                  nominal aero model through the exact 3DOF EOM), plus the
                  integrated altitude and speed as inputs: the density-profile
                  and Mach parts of the correction are functions of where the
                  *rollout* is, which the latent cannot know. Pseudo-Huber
                  state loss so long prober horizons stay stable.
    ensemble    M members trained in one batched pass. Every weight has a
                  leading member dimension (ELinear), so the whole ensemble is
                  one set of batched matmuls in training and in MPPI. Member
                  disagreement on the landing point is the planner's
                  epistemic-uncertainty signal (jepa/mppi_v2.py).

Shapes: M members, B batch, L sequence, T horizon, D latent, F features.
"""

from __future__ import annotations

import math
import os
import time
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.checkpoint import checkpoint as grad_checkpoint

from jepa.data import FEAT_NAMES, N_Y, SIGMA_MAX_RAD, raw_features
from jepa.models import (BETA_NOM, H_SCALE_M, LD_NOM, RHO0_KGM3, A_SCALE, integrate_eom,
                         state_errors)

# ============================================================================
# Features (v1 features + physics residuals)
# ============================================================================
FEAT_NAMES_V2 = tuple(FEAT_NAMES) + ("sin_sigma", "cos_sigma", "log_a_drag", "cD_meas", "cL_meas")
N_FEAT_V2 = len(FEAT_NAMES_V2)
A_FLOOR = 1e-6                   # m/s^2, floor before taking logs of sensed forces


def raw_features_v2(y: torch.Tensor) -> torch.Tensor:
    """
    Physical states [..., N_Y] -> v2 encoder features [..., N_FEAT_V2] (unnormalized).

    The last two are the corrections the PI prober would need under its crude
    nominal model, measured at the current state from the sensed forces:
        cD_meas = log(a_D / (0.5 rho_exp(h) V^2 / beta_nom))
                = log(density x / mass x) + log(rho_table / rho_exp)(h) + log(CD(M) / 1.5)
        cL_meas = log((a_L / a_D) / 0.30)
    """
    h, V, sig = y[..., 0], y[..., 3], y[..., 6]
    aD = y[..., 8].clamp_min(A_FLOOR)
    aL = y[..., 9].clamp_min(A_FLOOR)
    aD_nom = 0.5 * RHO0_KGM3 * torch.exp(-h.clamp_min(0.0) / H_SCALE_M) * V * V / BETA_NOM
    cD = torch.log(aD / aD_nom.clamp_min(1e-12)).clamp(-6.0, 6.0)
    cL = torch.log(aL / aD / LD_NOM).clamp(-6.0, 6.0)
    extra = torch.stack([torch.sin(sig), torch.cos(sig), torch.log(aD), cD, cL], -1)
    return torch.cat([raw_features(y), extra], -1)


@dataclass
class StatsV2:
    """Normalization statistics (training split only). Duck-types data.Stats."""
    feat_mean: torch.Tensor
    feat_std: torch.Tensor
    dy_std: torch.Tensor

    def norm_feats(self, y: torch.Tensor) -> torch.Tensor:
        """Physical states -> normalized v2 features (on y's device)."""
        return (raw_features_v2(y) - self.feat_mean.to(y.device)) / self.feat_std.to(y.device)

    def to_dict(self) -> dict:
        """Plain-tensor dict for checkpointing."""
        return {"feat_mean": self.feat_mean.cpu(), "feat_std": self.feat_std.cpu(), "dy_std": self.dy_std.cpu()}


def compute_stats_v2(macro_eps: Sequence[dict]) -> StatsV2:
    """Feature mean/std and per-step state-change std over macro episodes."""
    Y = torch.from_numpy(np.concatenate([e["y"] for e in macro_eps])).double()
    dY = torch.from_numpy(np.concatenate([np.diff(e["y"], axis=0) for e in macro_eps])).double()
    f = raw_features_v2(Y)
    return StatsV2(feat_mean=f.mean(0).float(), feat_std=f.std(0).clamp_min(1e-6).float(),
                   dy_std=dY.std(0).clamp_min(1e-9).float())


# ============================================================================
# Batched-ensemble building blocks (every weight has a leading member dim)
# ============================================================================
class ELinear(nn.Module):
    """M independent linear layers applied to x [M, ..., d_in] in one bmm."""

    def __init__(self, M: int, d_in: int, d_out: int, std: Optional[float] = None, bias: bool = True):
        """std=None -> 0.02 (transformer default); std=0 -> zero init."""
        super().__init__()
        self.weight = nn.Parameter(torch.empty(M, d_in, d_out))
        std = 0.02 if std is None else float(std)
        if std > 0:
            nn.init.trunc_normal_(self.weight, std=std, a=-2 * std, b=2 * std)
        else:
            nn.init.zeros_(self.weight)
        self.bias = nn.Parameter(torch.zeros(M, d_out)) if bias else None

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """x [M, ..., d_in] -> [M, ..., d_out]."""
        M, shp = x.shape[0], x.shape
        x2 = x.reshape(M, -1, shp[-1])
        y = torch.bmm(x2, self.weight.to(x2.dtype))
        if self.bias is not None:
            y = y + self.bias.to(y.dtype)[:, None]
        return y.reshape(*shp[:-1], y.shape[-1])


class ELayerNorm(nn.Module):
    """LayerNorm over the last dim with per-member affine parameters."""

    def __init__(self, M: int, d: int):
        super().__init__()
        self.weight = nn.Parameter(torch.ones(M, d))
        self.bias = nn.Parameter(torch.zeros(M, d))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """x [M, ..., d]."""
        y = F.layer_norm(x, (x.shape[-1],))
        view = (x.shape[0],) + (1,) * (x.dim() - 2) + (x.shape[-1],)
        return y * self.weight.view(view).to(y.dtype) + self.bias.view(view).to(y.dtype)


def _rope(L: int, dh: int, pos0: int, device, dtype) -> Tuple[torch.Tensor, torch.Tensor]:
    """Rotary-embedding cos/sin tables [L, dh/2] for positions pos0 .. pos0+L-1."""
    inv = 1.0 / (10000.0 ** (torch.arange(0, dh, 2, device=device, dtype=torch.float32) / dh))
    f = torch.outer(torch.arange(pos0, pos0 + L, device=device, dtype=torch.float32), inv)
    return f.cos().to(dtype), f.sin().to(dtype)


def _apply_rope(x: torch.Tensor, cos: torch.Tensor, sin: torch.Tensor) -> torch.Tensor:
    """Rotate consecutive channel pairs of x [..., L, dh]."""
    x1, x2 = x[..., 0::2], x[..., 1::2]
    return torch.stack([x1 * cos - x2 * sin, x1 * sin + x2 * cos], -1).flatten(-2)


def band_mask(L: int, band: int, device) -> torch.Tensor:
    """Boolean [L, L], True where query i may attend key j: i - band < j <= i."""
    i = torch.arange(L, device=device)
    d = i[:, None] - i[None, :]
    return (d >= 0) & (d < band)


class EAttention(nn.Module):
    """Multi-head self-attention with RoPE, per-member weights, optional KV cache."""

    def __init__(self, M: int, d: int, n_heads: int, out_std: float):
        super().__init__()
        assert d % n_heads == 0 and (d // n_heads) % 2 == 0
        self.h, self.dh = n_heads, d // n_heads
        self.qkv = ELinear(M, d, 3 * d)
        self.proj = ELinear(M, d, d, std=out_std)

    def forward(self, x: torch.Tensor, mask: Optional[torch.Tensor] = None, causal: bool = False,
                pos0: int = 0, cache: Optional[dict] = None) -> torch.Tensor:
        """
        x [M, B, L, d]. mask: boolean [L, L] (True = attend). causal: lower
        triangular. cache: dict holding this layer's rotated keys and values;
        new tokens attend to everything cached (used one token at a time).
        """
        M, B, L, d = x.shape
        qkv = self.qkv(x).view(M * B, L, 3, self.h, self.dh).permute(2, 0, 3, 1, 4)
        q, k, v = qkv[0], qkv[1], qkv[2]                       # [MB, H, L, dh]
        cos, sin = _rope(L, self.dh, pos0, x.device, q.dtype)
        q, k = _apply_rope(q, cos, sin), _apply_rope(k, cos, sin)
        if cache is not None:
            if "k" in cache:
                assert L == 1, "incremental decoding feeds one token at a time"
                k = torch.cat([cache["k"], k], 2)
                v = torch.cat([cache["v"], v], 2)
            else:
                causal = L > 1
            cache["k"], cache["v"] = k, v
            o = F.scaled_dot_product_attention(q, k, v, is_causal=causal and L > 1)
        else:
            o = F.scaled_dot_product_attention(q, k, v, attn_mask=mask, is_causal=causal and mask is None)
        return self.proj(o.transpose(1, 2).reshape(M, B, L, d))


class EBlock(nn.Module):
    """Pre-LN transformer block with per-member weights."""

    def __init__(self, M: int, d: int, n_heads: int, mlp_ratio: int, n_layers: int, dropout: float = 0.0):
        super().__init__()
        out_std = 0.02 / math.sqrt(2 * n_layers)
        self.ln1, self.ln2 = ELayerNorm(M, d), ELayerNorm(M, d)
        self.attn = EAttention(M, d, n_heads, out_std)
        self.fc1 = ELinear(M, d, mlp_ratio * d)
        self.fc2 = ELinear(M, mlp_ratio * d, d, std=out_std)
        self.dropout = float(dropout)

    def forward(self, x: torch.Tensor, **kw) -> torch.Tensor:
        """x [M, B, L, d]; keyword args go to the attention."""
        a = self.attn(self.ln1(x), **kw)
        x = x + F.dropout(a, self.dropout, self.training)
        m = self.fc2(F.gelu(self.fc1(self.ln2(x))))
        return x + F.dropout(m, self.dropout, self.training)


def _run_blocks(blocks: nn.ModuleList, x: torch.Tensor, use_ckpt: bool, **kw) -> torch.Tensor:
    """Apply blocks, optionally with activation checkpointing (training only)."""
    for blk in blocks:
        if use_ckpt and torch.is_grad_enabled():
            x = grad_checkpoint(lambda t, b=blk: b(t, **kw), x, use_reentrant=False)
        else:
            x = blk(x, **kw)
    return x


def _full(x: torch.Tensor) -> torch.Tensor:
    """Upcast half-precision (autocast) outputs to fp32; leave fp32/fp64 alone."""
    return x.float() if x.dtype in (torch.float16, torch.bfloat16) else x


def action_features(u: torch.Tensor) -> torch.Tensor:
    """Bank commands (rad) [..., L] -> [..., L, 3]: scaled command and its sine/cosine."""
    return torch.stack([u / SIGMA_MAX_RAD, torch.sin(u), torch.cos(u)], -1)


# ============================================================================
# World model
# ============================================================================
class WorldModelV2(nn.Module):
    """Banded-attention state encoder + chunked action tokens + causal transformer predictor."""

    def __init__(self, n_members: int = 5, n_feat: int = N_FEAT_V2, d_latent: int = 64,
                 enc_dim: int = 128, enc_layers: int = 3, enc_heads: int = 4, enc_band: int = 9,
                 pred_dim: int = 192, pred_layers: int = 4, pred_heads: int = 6, act_chunk: int = 4,
                 mlp_ratio: int = 4, dropout: float = 0.0):
        """Defaults: rf = 1 + 3 * 8 = 25 macro steps (120 s), ~2.4M parameters per member."""
        super().__init__()
        self.hparams = dict(n_members=n_members, n_feat=n_feat, d_latent=d_latent, enc_dim=enc_dim,
                            enc_layers=enc_layers, enc_heads=enc_heads, enc_band=enc_band,
                            pred_dim=pred_dim, pred_layers=pred_layers, pred_heads=pred_heads,
                            act_chunk=act_chunk, mlp_ratio=mlp_ratio, dropout=dropout)
        M = n_members
        self.M, self.d_latent, self.enc_band, self.c = M, d_latent, enc_band, act_chunk
        self.rf = 1 + enc_layers * (enc_band - 1)
        assert act_chunk <= self.rf
        self.grad_ckpt = False
        self.enc_in = ELinear(M, n_feat, enc_dim, std=n_feat ** -0.5)
        self.enc_blocks = nn.ModuleList(EBlock(M, enc_dim, enc_heads, mlp_ratio, enc_layers, dropout)
                                        for _ in range(enc_layers))
        self.enc_norm = ELayerNorm(M, enc_dim)
        self.enc_out = ELinear(M, enc_dim, d_latent, std=enc_dim ** -0.5)
        self.ctx_in = ELinear(M, d_latent, pred_dim, std=d_latent ** -0.5)
        self.act_in = ELinear(M, 3 * act_chunk, pred_dim, std=(3 * act_chunk) ** -0.5)
        self.pred_blocks = nn.ModuleList(EBlock(M, pred_dim, pred_heads, mlp_ratio, pred_layers, dropout)
                                         for _ in range(pred_layers))
        self.pred_norm = ELayerNorm(M, pred_dim)
        # Zero init: before training the prediction is "the latent does not change".
        self.readout = ELinear(M, pred_dim, d_latent, std=0.0)

    def _members(self, x: torch.Tensor, member_dims: int) -> torch.Tensor:
        """Broadcast an input without a member dim (x.dim() == member_dims - 1) to [M, ...]."""
        if x.dim() == member_dims - 1:
            return x.unsqueeze(0).expand(self.M, *x.shape)
        return x

    def encode_states(self, feats: torch.Tensor) -> torch.Tensor:
        """Normalized features [M or -, B, L, F] -> latents [M, B, L - rf + 1, D]."""
        x = self.enc_in(self._members(feats, 4))
        x = _run_blocks(self.enc_blocks, x, self.grad_ckpt,
                        mask=band_mask(x.shape[2], self.enc_band, x.device))
        return _full(self.enc_out(self.enc_norm(x[:, :, self.rf - 1:])))   # latents never in half precision

    def action_tokens(self, u: torch.Tensor) -> torch.Tensor:
        """
        Commands (rad) [M or -, B, L] -> tokens [M, B, L - rf + 1, d_pred].

        Token j covers the command applied over the (rf-1+j)-th macro step and
        the c-1 before it, so with L = rf-1 history + K planned commands there
        is one token per planned step (the same alignment as v1's encode_actions).
        """
        u = self._members(u, 3)
        a = action_features(u).unfold(2, self.c, 1)           # [M, B, L-c+1, 3, c]
        a = a.transpose(-1, -2).flatten(-2)                   # [M, B, L-c+1, 3c]
        return self.act_in(a[:, :, self.rf - self.c:])

    def predict(self, s0: torch.Tensor, a: torch.Tensor) -> torch.Tensor:
        """Parallel rollout: s0 [M, B, D], tokens a [M, B, K, d] -> latents s~_t+1..t+K [M, B, K, D]."""
        ctx = self.ctx_in(s0)[:, :, None]
        x = _run_blocks(self.pred_blocks, torch.cat([ctx, a.to(ctx.dtype)], 2), self.grad_ckpt, causal=True)
        return _full(s0)[:, :, None] + _full(self.readout(self.pred_norm(x[:, :, 1:])))

    # ---- incremental (KV-cached) predictor, the v1 step-by-step API -----------
    def start(self, s0: torch.Tensor) -> dict:
        """Feed the context token; returns the KV cache."""
        cache = {"s0": s0, "pos": 1, "layers": [dict() for _ in self.pred_blocks]}
        x = self.ctx_in(s0)[:, :, None]
        for blk, c in zip(self.pred_blocks, cache["layers"]):
            x = blk(x, pos0=0, cache=c)
        return cache

    def step_token(self, a_k: torch.Tensor, cache: dict) -> Tuple[torch.Tensor, dict]:
        """One action token a_k [M, B, d] -> (next latent [M, B, D], cache)."""
        x = a_k[:, :, None]
        for blk, c in zip(self.pred_blocks, cache["layers"]):
            x = blk(x, pos0=cache["pos"], cache=c)
        cache["pos"] += 1
        return _full(cache["s0"]) + _full(self.readout(self.pred_norm(x[:, :, 0]))), cache


# ============================================================================
# Per-step SIGReg
# ============================================================================
class SIGRegV2(nn.Module):
    """
    v1's SIGReg (Epps-Pulley on random 1-D projections, same knots and
    weighting, statistic x N), evaluated separately per group of rows with an
    optional row mask. Grouping by rollout step gives the paper's per-step
    SIGReg; each member gets its own statistic.
    """

    def __init__(self, n_proj: int = 32, n_knots: int = 17, t_max: float = 4.0):
        super().__init__()
        t = torch.linspace(0.0, t_max, n_knots)
        w = torch.full((n_knots,), t_max / (n_knots - 1))
        w[0] = w[-1] = 0.5 * w[0]
        phi0 = torch.exp(-0.5 * t ** 2)
        self.n_proj = n_proj
        self.register_buffer("t", t)
        self.register_buffer("phi0", phi0)
        self.register_buffer("w", 2.0 * w * phi0)

    def forward(self, x: torch.Tensor, mask: Optional[torch.Tensor] = None,
                A: Optional[torch.Tensor] = None) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        x [G, N, D] (G groups of N rows), mask [G, N] (1 = row counts).
        Returns (statistic per group [G], valid row count per group [G]).
        """
        x = x.float()
        if A is None:
            A = torch.randn(x.shape[-1], self.n_proj, device=x.device)
            A = A / A.norm(dim=0, keepdim=True)
        xt = (x @ A).unsqueeze(-1) * self.t                   # [G, N, P, K]
        if mask is None:
            n = torch.full(x.shape[:1], float(x.shape[1]), device=x.device)
            c, s = torch.cos(xt).mean(1), torch.sin(xt).mean(1)
        else:
            m = mask.float()
            n = m.sum(1)
            wm = (m / n.clamp_min(1.0)[:, None])[:, :, None, None]
            c, s = (torch.cos(xt) * wm).sum(1), (torch.sin(xt) * wm).sum(1)
        err = (c - self.phi0) ** 2 + s ** 2                   # [G, P, K]
        return (err @ self.w).mean(-1) * n, n


def per_step_sigreg(sig: SIGRegV2, x: torch.Tensor, mask: torch.Tensor, steps: torch.Tensor,
                    min_rows: int) -> torch.Tensor:
    """
    Mean over the chosen rollout steps of SIGReg on x[:, :, k] (one statistic
    per member per step). x [M, B, T, D], mask [M, B, T]. Steps with fewer than
    min_rows valid rows (late steps of short flights) are skipped. Returns [M].
    """
    M, B, _, D = x.shape
    xs = x[:, :, steps].permute(0, 2, 1, 3).reshape(M * len(steps), B, D)
    ms = mask[:, :, steps].permute(0, 2, 1).reshape(M * len(steps), B)
    stat, n = sig(xs, ms)
    ok = (n >= min_rows).float()
    stat, ok = stat.view(M, -1), ok.view(M, -1)
    return (stat * ok).sum(1) / ok.sum(1).clamp_min(1.0)


# ============================================================================
# Physics-inspired prober (v1 structure, state-aware, ensemble-batched)
# ============================================================================
H_NORM = (60_000.0, 30_000.0)    # (center, scale) of the altitude input
V_NORM = (4_000.0, 2_500.0)      # (center, scale) of the speed input


class PIProberV2(nn.Module):
    """
    Latent (+ integrated h, V) -> (c_D, c_L, bank offset) -> v1's crude nominal
    aero model -> exact 3DOF EOM (models.integrate_eom). Zero-initialized, so
    before training it is the "physics prior only" model, as in v1.
    """

    def __init__(self, n_members: int, d_latent: int, hidden: int = 128, depth: int = 3,
                 state_aware: bool = True):
        super().__init__()
        self.hparams = dict(n_members=n_members, d_latent=d_latent, hidden=hidden, depth=depth,
                            state_aware=state_aware)
        self.M, self.state_aware = n_members, state_aware
        d_in = d_latent + (2 if state_aware else 0)
        layers: List[nn.Module] = []
        for i in range(depth):
            layers += [ELinear(n_members, d_in if i == 0 else hidden, hidden, std=(d_in if i == 0 else hidden) ** -0.5),
                       nn.GELU()]
        layers.append(ELinear(n_members, hidden, 3, std=0.0))
        self.net = nn.Sequential(*layers)

    def corrections(self, s: torch.Tensor, y6: torch.Tensor, u: torch.Tensor) -> torch.Tensor:
        """s [M, B, D], y6 [M, B, 6], u [M, B] -> (c_D, c_L, sigma_eff) [M, B, 3]."""
        x = s.float()
        if self.state_aware:
            x = torch.cat([x, ((y6[..., 0] - H_NORM[0]) / H_NORM[1])[..., None],
                           ((y6[..., 3] - V_NORM[0]) / V_NORM[1])[..., None]], -1)
        o = self.net(x)
        return torch.stack([o[..., 0].clamp(-5.0, 5.0), o[..., 1].clamp(-3.0, 3.0),
                            u + SIGMA_MAX_RAD * torch.tanh(o[..., 2])], -1)

    def step(self, s: torch.Tensor, y6: torch.Tensor, u: torch.Tensor, dt: float, n_sub: int):
        """One macro step for every member: returns (y6_next [M, B, 6], mean forces [M, B, 3])."""
        return integrate_eom(y6, self.corrections(s, y6, u), dt, n_sub)


def prober_rollout_v2(prober: PIProberV2, lat: torch.Tensor, y0: torch.Tensor, u: torch.Tensor,
                      dt: float, n_sub: int, h_stop: Optional[float] = None):
    """
    Integrate latents lat [M, B, T, D] from y0 [M, B, 6] under commands u [M, B, T].
    The correction for step k comes from s~_t+k+1, as in v1. Returns
    (y6 [M, B, T'+1, 6], acc [M, B, T', 3]); T' < T only if h_stop ended early.
    """
    ys, accs, y = [y0], [], y0
    for k in range(lat.shape[2]):
        y, a = prober.step(lat[:, :, k], y, u[:, :, k], dt, n_sub)
        ys.append(y)
        accs.append(a)
        if h_stop is not None and k % 5 == 4 and bool((y[..., 0] < h_stop).all()):
            break
    return torch.stack(ys, 2), torch.stack(accs, 2)


# ============================================================================
# Training windows on the GPU
# ============================================================================
class WindowsV2:
    """
    All training windows of a set of macro episodes, resident on one device.

    Like data.MacroDataset (front pad with rf-1 copies of the entry state and
    wings-level commands), but every start t in [0, L-1] is a window for every
    horizon: episodes are also tail-padded and a mask marks the real future
    steps, so short flights still train the long horizons (wide-envelope
    flights are 46-150 macro steps; dropping the short ones from T = 96 windows
    would leave mostly long glides there).
    """

    def __init__(self, macro_eps: Sequence[dict], rf: int, max_T: int, stats: StatsV2,
                 device: Union[str, torch.device] = "cpu"):
        self.P, self.max_T, self.device = int(rf) - 1, int(max_T), torch.device(device)
        P, Tm = self.P, self.max_T
        Ys, Us, As, ys0, us0, rem, ep_id = [], [], [], [], [], [], []
        y_off = u_off = 0
        for i, ep in enumerate(macro_eps):
            y, u, acc = ep["y"], ep["u"], ep["acc"]
            L = len(u)
            if L < 1:
                continue
            Ys.append(np.concatenate([np.repeat(y[:1], P, 0), y, np.repeat(y[-1:], Tm, 0)]))
            Us.append(np.concatenate([np.zeros(P, np.float32), u, np.repeat(u[-1:], Tm)]))
            As.append(np.concatenate([np.zeros((P, 3), np.float32), acc, np.zeros((Tm, 3), np.float32)]))
            t = np.arange(L)
            ys0.append(y_off + t)
            us0.append(u_off + t)
            rem.append(L - t)
            ep_id.append(np.full(L, i))
            y_off += P + L + 1 + Tm
            u_off += P + L + Tm
        dev = self.device
        self.Y = torch.from_numpy(np.concatenate(Ys)).float().to(dev)
        self.F = stats.norm_feats(self.Y)
        self.U = torch.from_numpy(np.concatenate(Us)).float().to(dev)
        self.A = torch.from_numpy(np.concatenate(As)).float().to(dev)
        self.y_start = torch.from_numpy(np.concatenate(ys0)).long().to(dev)
        self.u_start = torch.from_numpy(np.concatenate(us0)).long().to(dev)
        self.rem = torch.from_numpy(np.concatenate(rem)).long().to(dev)
        self.ep_id = torch.from_numpy(np.concatenate(ep_id)).long().to(dev)

    def __len__(self) -> int:
        return len(self.y_start)

    def batch(self, idx: torch.Tensor, T: int) -> Dict[str, torch.Tensor]:
        """
        Windows idx (any shape S, e.g. [M, B]) at horizon T:
            feats [S, P+T+1, F], u [S, P+T], y [S, T+1, N_Y], acc [S, T, 3], mask [S, T]
        """
        assert T <= self.max_T
        P = self.P
        iy = self.y_start[idx][..., None] + torch.arange(P + T + 1, device=self.device)
        iu = self.u_start[idx][..., None] + torch.arange(P + T, device=self.device)
        k = torch.arange(1, T + 1, device=self.device)
        return {"feats": self.F[iy], "u": self.U[iu], "y": self.Y[iy[..., P:]],
                "acc": self.A[iu[..., P:]], "mask": (k <= self.rem[idx][..., None]).float()}

    def sample(self, shape: Tuple[int, ...], T: int, gen: Optional[torch.Generator] = None):
        """Uniformly random windows (one independent draw per member and row)."""
        idx = torch.randint(len(self), shape, generator=gen, device=self.device)
        return self.batch(idx, T)


# ============================================================================
# Losses
# ============================================================================
def wm_loss_v2(wm: WorldModelV2, sig: SIGRegV2, b: Dict[str, torch.Tensor],
               lam: Union[float, torch.Tensor], sig_steps: int = 16, sig_min_rows: int = 32):
    """
    Masked multi-step latent prediction + per-step SIGReg, for all members.

    b holds [M, B, ...] windows (WindowsV2.sample). lam is a float or a
    per-member tensor [M] (a lambda sweep trains one member per value).
    SIGReg is evaluated at sig_steps random rollout steps per iteration
    (an unbiased estimate of the mean over all steps; all of them is ~T/16 x
    the memory), on the encoded targets and on the predictions.
    """
    S = wm.encode_states(b["feats"])                           # [M, B, T+1, D]
    pred = wm.predict(S[:, :, 0], wm.action_tokens(b["u"]))    # [M, B, T, D]
    m = b["mask"]
    T = m.shape[-1]
    err = (pred.float() - S[:, :, 1:].float()).pow(2).mean(-1)
    l_pred = (err * m).sum((1, 2)) / m.sum((1, 2)).clamp_min(1.0)        # [M]
    dev = m.device
    k = torch.randperm(T, device=dev)[:max(1, sig_steps - 1)] + 1             # rollout steps 1..T
    k_tgt = torch.cat([torch.zeros(1, dtype=torch.long, device=dev), k])     # + the context latent
    k_pred = k - 1                                                           # pred[:, :, j] is step j+1
    m_tgt = torch.cat([torch.ones_like(m[:, :, :1]), m], 2)
    l_sig = (per_step_sigreg(sig, S, m_tgt, k_tgt, sig_min_rows)
             + per_step_sigreg(sig, pred, m, k_pred, sig_min_rows))         # [M]
    lam_t = torch.as_tensor(lam, dtype=torch.float32, device=dev)
    loss = (l_pred + lam_t * l_sig).sum()
    return loss, {"pred": float(l_pred.mean().detach()), "sigreg": float(l_sig.mean().detach()),
                  "pred_m": l_pred.detach().cpu(), "sigreg_m": l_sig.detach().cpu()}


def pseudo_huber(e: torch.Tensor) -> torch.Tensor:
    """sqrt(1 + e^2) - 1: quadratic below 1, linear above (robust long-horizon errors)."""
    return torch.sqrt(1.0 + e * e) - 1.0


def prober_loss_v2(wm: WorldModelV2, prober: PIProberV2, b: Dict[str, torch.Tensor], dt: float,
                   n_sub: int, w_acc: float = 1.0):
    """
    Supervised rollout loss on frozen latents, every member on its own windows.
    Pseudo-Huber on v1's dimensionless state errors (1 km, 5 km, 50 m/s,
    0.5 deg, 0.5 deg), masked past the end of the flight, plus v1's auxiliary
    match to the sensed macro-mean specific forces.
    """
    with torch.no_grad():
        S0 = wm.encode_states(b["feats"][:, :, :wm.rf])[:, :, 0]
        lat = wm.predict(S0, wm.action_tokens(b["u"])).float()
    m = b["mask"]
    T = m.shape[-1]
    y6 = b["y"][..., :6]
    with torch.autocast(lat.device.type, enabled=False):            # physics in fp32, always
        yhat, acc = prober_rollout_v2(prober, lat, y6[:, :, 0], b["u"][:, :, -T:], dt, n_sub)
    e = pseudo_huber(state_errors(yhat[:, :, 1:], y6[:, :, 1:])).mean(-1)
    denom = m.sum((1, 2)).clamp_min(1.0)
    l_state = (e * m).sum((1, 2)) / denom
    l_acc = ((((acc - b["acc"]) / A_SCALE) ** 2).mean(-1) * m).sum((1, 2)) / denom
    loss = (l_state + w_acc * l_acc).sum()
    return loss, {"state": float(l_state.mean().detach()), "acc": float(l_acc.mean().detach()),
                  "state_m": l_state.detach().cpu()}


# ============================================================================
# Training loop (resumable, AMP, horizon curriculum)
# ============================================================================
@dataclass
class TrainCfgV2:
    """One training stage. curriculum: (fraction of steps, horizon) pairs."""
    steps: int = 20_000
    batch: int = 256                  # windows per member per step
    lr: float = 6e-4
    lr_min: float = 2e-5
    warmup: int = 500
    weight_decay: float = 0.05
    clip: float = 1.0
    curriculum: Tuple[Tuple[float, int], ...] = ((0.0, 24), (0.2, 48), (0.5, 96))
    log_every: int = 250
    save_every_s: float = 600.0
    amp: bool = True                  # bf16 if supported, else fp16 + GradScaler (CUDA only)
    seed: int = 0


def horizon_at(step: int, cfg: TrainCfgV2) -> int:
    """Current curriculum horizon."""
    T = cfg.curriculum[0][1]
    for frac, h in cfg.curriculum:
        if step >= frac * cfg.steps:
            T = h
    return T


def _lr(step: int, cfg: TrainCfgV2) -> float:
    """Linear warmup, cosine decay."""
    if step < cfg.warmup:
        return cfg.lr * (step + 1) / cfg.warmup
    p = (step - cfg.warmup) / max(1, cfg.steps - cfg.warmup)
    return cfg.lr_min + 0.5 * (cfg.lr - cfg.lr_min) * (1.0 + math.cos(math.pi * min(p, 1.0)))


def _param_groups(module: nn.Module, wd: float):
    """Weight decay on matrices only (not biases, norms)."""
    decay, no = [], []
    for n, p in module.named_parameters():
        if not p.requires_grad:
            continue
        (decay if p.dim() >= 3 and "ln" not in n and "norm" not in n else no).append(p)
    return [{"params": decay, "weight_decay": wd}, {"params": no, "weight_decay": 0.0}]


def amp_dtype(device: torch.device) -> Optional[torch.dtype]:
    """bf16 on Ampere+ GPUs, fp16 on older ones (T4), None on CPU."""
    if device.type != "cuda":
        return None
    return torch.bfloat16 if torch.cuda.is_bf16_supported() else torch.float16


def train_stage(module: nn.Module, loss_fn: Callable[[int, int], tuple], cfg: TrainCfgV2,
                device: torch.device, ckpt_path: Optional[str] = None,
                val_fn: Optional[Callable[[], Dict[str, float]]] = None, label: str = "") -> Dict[str, list]:
    """
    Generic resumable AdamW loop.

    loss_fn(step, T) draws its own batch at horizon T and returns (loss, logs).
    Every save_every_s seconds (and at the end) the module, optimizer, scaler,
    step, and history are written to ckpt_path atomically; a rerun resumes
    from there, so a Colab disconnect costs at most save_every_s of training.
    """
    opt = torch.optim.AdamW(_param_groups(module, cfg.weight_decay), lr=cfg.lr, betas=(0.9, 0.98))
    dt_amp = amp_dtype(device) if cfg.amp else None
    scaler = torch.amp.GradScaler(device.type, enabled=dt_amp == torch.float16)
    hist: Dict[str, list] = {"step": [], "T": [], "val": [], "sec": []}
    start, elapsed = 0, 0.0
    if ckpt_path and os.path.exists(ckpt_path):
        st = torch.load(ckpt_path, map_location=device, weights_only=False)
        module.load_state_dict(st["module"])
        opt.load_state_dict(st["opt"])
        scaler.load_state_dict(st["scaler"])
        start, hist, elapsed = st["step"], st["hist"], st.get("elapsed", 0.0)
        print(f"  [{label}] resumed at step {start}/{cfg.steps} from {ckpt_path}")

    def save(step):
        if not ckpt_path:
            return
        tmp = ckpt_path + ".tmp"
        torch.save({"module": module.state_dict(), "opt": opt.state_dict(), "scaler": scaler.state_dict(),
                    "step": step, "hist": hist, "elapsed": elapsed}, tmp)
        os.replace(tmp, ckpt_path)

    params = [p for g in opt.param_groups for p in g["params"]]
    acc: Dict[str, float] = {}
    n_acc, t_log, t_save = 0, time.time(), time.time()
    module.train()
    for step in range(start, cfg.steps):
        for g in opt.param_groups:
            g["lr"] = _lr(step, cfg)
        T = horizon_at(step, cfg)
        with torch.autocast(device.type, dtype=dt_amp, enabled=dt_amp is not None):
            loss, logs = loss_fn(step, T)
        opt.zero_grad(set_to_none=True)
        scaler.scale(loss).backward()
        scaler.unscale_(opt)
        torch.nn.utils.clip_grad_norm_(params, cfg.clip)
        scaler.step(opt)
        scaler.update()
        for k, v in logs.items():
            if isinstance(v, float):
                acc[k] = acc.get(k, 0.0) + v
        n_acc += 1
        last = step + 1 == cfg.steps
        if (step + 1) % cfg.log_every == 0 or last:
            if device.type == "cuda":
                torch.cuda.synchronize()
            now = time.time()
            elapsed += now - t_log
            per = (now - t_log) / n_acc
            t_log = now
            hist["step"].append(step + 1)
            hist["T"].append(T)
            hist["sec"].append(round(elapsed, 1))
            for k, v in acc.items():
                hist.setdefault(k, []).append(v / n_acc)
            val = {}
            if val_fn is not None:
                module.eval()
                with torch.no_grad(), torch.autocast(device.type, dtype=dt_amp, enabled=dt_amp is not None):
                    val = val_fn()
                module.train()
            hist["val"].append(val)
            eta = per * (cfg.steps - step - 1) / 60
            msg = " ".join(f"{k}={v / n_acc:.4g}" for k, v in acc.items())
            vmsg = " ".join(f"val_{k}={v:.4g}" for k, v in val.items())
            print(f"  [{label}] step {step + 1:>6}/{cfg.steps} T={T:<3} {msg} {vmsg} "
                  f"({per * 1000:.0f} ms/step, {elapsed / 60:.1f} min, ~{eta:.0f} min left)", flush=True)
            acc, n_acc = {}, 0
        if last or (ckpt_path and time.time() - t_save > cfg.save_every_s):
            save(step + 1)
            t_save = time.time()
    module.eval()
    return hist


# ============================================================================
# v1-compatible views (so jepa.models.jepa_predict, evaluate.py and mppi.py run v2)
# ============================================================================
class V1WorldModelView:
    """
    Presents a WorldModelV2 through v1's JEPAWorldModel interface
    (rf, d_latent, encode_states, encode_actions, init_hidden, step), with the
    members flattened into the latent dimension. Inputs may live on the CPU;
    the work runs on the model's device.

    members: which members to expose (default all). Paired with V1ProberView,
    the exposed members' prober corrections are averaged, which turns the
    ensemble into one deterministic model ("ensemble mean").
    """

    def __init__(self, wm: WorldModelV2, members: Optional[Sequence[int]] = None):
        self.wm = wm.eval()
        self.members = list(range(wm.M)) if members is None else list(members)
        self.rf, self.D = wm.rf, wm.d_latent
        self.d_latent = self.D * len(self.members)
        self.device = next(wm.parameters()).device

    def _flat(self, x: torch.Tensor) -> torch.Tensor:
        """[M, B, ..., X] -> [B, ..., m*X] for the exposed members."""
        x = x[self.members]
        return x.movedim(0, -2).flatten(-2)

    def _unflat(self, x: torch.Tensor) -> torch.Tensor:
        """[B, m*X] -> [m, B, X]."""
        return x.view(x.shape[0], len(self.members), -1).movedim(1, 0)

    @torch.no_grad()
    def encode_states(self, feats: torch.Tensor) -> torch.Tensor:
        return self._flat(self.wm.encode_states(feats.to(self.device)))

    @torch.no_grad()
    def encode_actions(self, u: torch.Tensor) -> torch.Tensor:
        return self._flat(self.wm.action_tokens(u.to(self.device)))

    def _full(self, x: torch.Tensor) -> torch.Tensor:
        """Exposed-member tensor [m, B, X] -> all-member tensor [M, B, X] (others zero)."""
        if len(self.members) == self.wm.M:
            return x
        out = x.new_zeros((self.wm.M,) + tuple(x.shape[1:]))
        out[self.members] = x
        return out

    @torch.no_grad()
    def init_hidden(self, s0: torch.Tensor) -> dict:
        return self.wm.start(self._full(self._unflat(s0)))

    @torch.no_grad()
    def step(self, z_k: torch.Tensor, hid: dict) -> Tuple[torch.Tensor, dict]:
        s, hid = self.wm.step_token(self._full(self._unflat(z_k)), hid)
        return self._flat(s), hid


class V1ProberView:
    """v1 prober interface over a PIProberV2: averages the exposed members' corrections."""

    def __init__(self, prober: PIProberV2, members: Sequence[int]):
        self.p = prober.eval()
        self.members = list(members)
        self.device = next(prober.parameters()).device

    @torch.no_grad()
    def step(self, s: torch.Tensor, y6: torch.Tensor, u: torch.Tensor, dt: float, n_sub: int):
        dev_in = y6.device
        m = len(self.members)
        s = s.to(self.device).view(s.shape[0], m, -1).movedim(1, 0)          # [m, B, D]
        full = s.new_zeros((self.p.M,) + tuple(s.shape[1:]))
        full[self.members] = s
        y = y6.to(self.device).float()
        corr = self.p.corrections(full, y.expand(self.p.M, *y.shape), u.to(self.device).float()
                                  .expand(self.p.M, *u.shape))[self.members]
        # Average the bank through its unit vector (a plain mean is fine away from +-180 deg).
        c = corr.mean(0)
        sig = torch.atan2(torch.sin(corr[..., 2]).mean(0), torch.cos(corr[..., 2]).mean(0))
        c = torch.stack([c[..., 0], c[..., 1], sig], -1)
        y_next, a = integrate_eom(y, c, dt, n_sub)
        return y_next.to(dev_in), a.to(dev_in)


def as_v1(models: dict, members: Optional[Sequence[int]] = None) -> dict:
    """
    {"kind": "jepa", "wm", "prober"} entry for evaluate.footprint_errors /
    rollout_errors, or the models dict for mppi.LearnedModelController.
    members=None -> the ensemble mean; [i] -> member i alone.
    """
    wm = V1WorldModelView(models["wm"], members)
    return {"kind": "jepa", "wm": wm, "prober": V1ProberView(models["pi"], wm.members),
            "pi": V1ProberView(models["pi"], wm.members), "stats": models["stats"]}


# ============================================================================
# Batched ensemble rollout (MPPI v2 and evaluation)
# ============================================================================
@torch.no_grad()
def ensemble_rollout(wm: WorldModelV2, prober: PIProberV2, s0: torch.Tensor, y0: torch.Tensor,
                     u_hist: torch.Tensor, u_plan: torch.Tensor, dt: float, n_sub: int,
                     h_stop: Optional[float] = None, chunk: int = 4096) -> Dict[str, torch.Tensor]:
    """
    Roll every member forward under planned commands, all on the model's device.

        s0     [M, B, D]  encoded context latent (WorldModelV2.encode_states(...)[:, :, -1])
        y0     [B, 6]     current translational state
        u_hist [B, rf-1]  previous commands;  u_plan [B, K] planned commands (rad)

    The predictor runs as one parallel pass (in chunks of `chunk` rows, under
    the ambient autocast), the prober integrates in fp32 with early stop at
    h_stop. Returns y6 [B, M, K'+1, 6], acc [B, M, K', 3], latents [M, B, K', D].
    """
    M, B = s0.shape[:2]
    u_all = torch.cat([u_hist, u_plan], 1)
    lats = []
    for i in range(0, B, chunk):
        lats.append(wm.predict(s0[:, i:i + chunk], wm.action_tokens(u_all[i:i + chunk])))
    lat = torch.cat(lats, 1).float()
    with torch.autocast(lat.device.type, enabled=False):
        y6, acc = prober_rollout_v2(prober, lat, y0.float().expand(M, *y0.shape), u_plan.float()
                                    .expand(M, *u_plan.shape), dt, n_sub, h_stop)
    K = acc.shape[2]
    return {"y6": y6.transpose(0, 1), "acc": acc.transpose(0, 1), "latents": lat[:, :, :K]}


# ============================================================================
# Checkpoints and parameter counts
# ============================================================================
def save_checkpoint_v2(path: str, wm: WorldModelV2, prober: PIProberV2, stats: StatsV2,
                       extra: Optional[dict] = None) -> None:
    """World model + prober + stats, with constructor args (atomic write)."""
    tmp = path + ".tmp"
    torch.save({"version": 2, "wm": (wm.hparams, {k: v.cpu() for k, v in wm.state_dict().items()}),
                "pi": (prober.hparams, {k: v.cpu() for k, v in prober.state_dict().items()}),
                "stats": stats.to_dict(), "extra": extra or {}}, tmp)
    os.replace(tmp, path)


def load_checkpoint_v2(path: str, device: Union[str, torch.device] = "cpu") -> dict:
    """Rebuild the models saved by save_checkpoint_v2 (eval mode, on device)."""
    ck = torch.load(path, map_location="cpu", weights_only=False)
    wm = WorldModelV2(**ck["wm"][0])
    wm.load_state_dict(ck["wm"][1])
    pi = PIProberV2(**ck["pi"][0])
    pi.load_state_dict(ck["pi"][1])
    return {"wm": wm.to(device).eval(), "pi": pi.to(device).eval(), "stats": StatsV2(**ck["stats"]),
            "extra": ck["extra"]}


def select_members(models: dict, members: Sequence[int]) -> dict:
    """A new models dict holding only the given ensemble members (e.g. one lambda of a sweep)."""
    members = list(members)
    hp = dict(models["wm"].hparams, n_members=len(members))
    wm = WorldModelV2(**hp)
    wm.load_state_dict({k: v[members] for k, v in models["wm"].state_dict().items()})
    php = dict(models["pi"].hparams, n_members=len(members))
    pi = PIProberV2(**php)
    pi.load_state_dict({k: v[members] for k, v in models["pi"].state_dict().items()})
    dev = next(models["wm"].parameters()).device
    return {"wm": wm.to(dev).eval(), "pi": pi.to(dev).eval(), "stats": models["stats"],
            "extra": models.get("extra", {})}


def param_report(wm: WorldModelV2, prober: Optional[PIProberV2] = None) -> Dict[str, int]:
    """Parameter counts per component, per member and for the whole ensemble."""
    M = wm.M
    groups = {"state encoder": [wm.enc_in, wm.enc_blocks, wm.enc_norm, wm.enc_out],
              "action tokens": [wm.act_in], "predictor": [wm.ctx_in, wm.pred_blocks, wm.pred_norm, wm.readout]}
    out = {k: sum(p.numel() for m in v for p in m.parameters()) // M for k, v in groups.items()}
    if prober is not None:
        out["PI prober"] = sum(p.numel() for p in prober.parameters()) // prober.M
    out["per member"] = sum(v for k, v in out.items())
    out[f"ensemble (x{M})"] = out["per member"] * M
    return out


def train_flops_per_step(wm: WorldModelV2, B: int, T: int) -> float:
    """Approximate training FLOPs (forward + backward = 3x forward matmuls) for one step."""
    hp = wm.hparams
    P = wm.rf - 1
    enc_tok, pred_tok = B * (P + T + 1), B * (T + 1)
    e, p = hp["enc_dim"], hp["pred_dim"]
    r = hp["mlp_ratio"]
    enc = enc_tok * hp["enc_layers"] * (2 * (4 + 2 * r) * e * e + 4 * hp["enc_band"] * e)
    pred = pred_tok * hp["pred_layers"] * (2 * (4 + 2 * r) * p * p + 4 * (T + 1) / 2 * p)
    return 3.0 * wm.M * (enc + pred)
