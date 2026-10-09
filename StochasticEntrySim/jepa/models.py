"""
jepa/models.py

The networks, losses, training loops, and batched rollouts.

JEPA world model (SkyJEPA Sec. IV-A/B):
    s_t   = Enc_theta(X_t)            causal TCN over the last rf macro states
    z_t   = Enc_phi(A_t)              causal TCN over the last rf bank commands
    s_t+1 = Pred(s_t, z_t)            GRU, unrolled T steps on encoded actions
    loss  = mean ||s~_t+k - s_t+k||^2 + lambda * SIGReg     (no decoder, no EMA)

Probers (Sec. IV-C), trained second on frozen latents:
    PIProber   latent -> three corrections (log drag factor, log L/D factor,
               bank offset) to a crude nominal aero model (exponential
               atmosphere, one ballistic coefficient, lift along the commanded
               bank), integrated through the exact 3DOF kinematics + gravity of
               math_3d.nominal_eom_step. The latent only carries what is hard
               to model (real density profile, Mach-scheduled CD and L/D,
               actuator lag, the randomized truth); the geometry is free.
    MLPProber  ablation: latent -> free-form per-step state increment, summed.

PredictiveModel is the paper's "Predictive" baseline: same encoders, but it
regresses the next physical state directly and is rolled out on its own outputs.

Everything is small on purpose (~25k parameters for the world model) so it
trains on a laptop CPU in minutes and MPPI can afford 256 x 100-step rollouts
per decision.
"""

from __future__ import annotations

import math
import time
from typing import Callable, Dict, Optional, Sequence, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F

import constants
from jepa.data import N_FEAT, N_Y, R_EARTH, SIGMA_MAX_RAD, Stats

MU = float(constants.MU)
A_SCALE = 10.0                      # m/s^2, scale of the auxiliary force loss
# Error scales that make the prober's state loss dimensionless. Position is
# looser than altitude because footprint error is what guidance cares about
# and a 5 km miss is already a good landing.
ERR_SCALE = {"h_m": 1000.0, "pos_m": 5000.0, "V_mps": 50.0,
             "gamma_rad": math.radians(0.5), "chi_rad": math.radians(0.5)}


# ============================================================================
# Building blocks
# ============================================================================
class CausalTCN(nn.Module):
    """
    Dilated causal convolutions with no padding.

    [B, C_in, L] -> [B, C_out, L - rf + 1]. Output j only sees inputs j .. j+rf-1,
    so running it over a whole sequence yields every sliding-window encoding in
    one pass (used to encode all T+1 targets, and all planned actions, at once).
    """

    def __init__(self, c_in: int, channels: Sequence[int], kernel: int = 3,
                 dilations: Sequence[int] = (1, 2)):
        """Stack one conv + GELU per (channel, dilation) pair."""
        super().__init__()
        layers, c = [], c_in
        for ch, d in zip(channels, dilations):
            layers += [nn.Conv1d(c, ch, kernel, dilation=d), nn.GELU()]
            c = ch
        self.net = nn.Sequential(*layers)
        self.rf = 1 + sum((kernel - 1) * d for d in dilations)
        self.c_out = c

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Apply the conv stack."""
        return self.net(x)


def _mlp(d_in: int, hidden: int, d_out: int) -> nn.Sequential:
    """Two hidden layer GELU MLP."""
    return nn.Sequential(nn.Linear(d_in, hidden), nn.GELU(),
                         nn.Linear(hidden, hidden), nn.GELU(), nn.Linear(hidden, d_out))


# ============================================================================
# JEPA world model
# ============================================================================
class JEPAWorldModel(nn.Module):
    """State encoder, action encoder, and GRU latent predictor."""

    def __init__(self, n_feat: int = N_FEAT, d_latent: int = 32, z_dim: int = 16,
                 hidden: int = 64, state_ch: Sequence[int] = (32, 32),
                 act_ch: Sequence[int] = (16, 16), kernel: int = 3,
                 dilations: Sequence[int] = (1, 2)):
        """Build the encoders and predictor. The two TCNs share a receptive field."""
        super().__init__()
        self.hparams = dict(n_feat=n_feat, d_latent=d_latent, z_dim=z_dim, hidden=hidden,
                            state_ch=tuple(state_ch), act_ch=tuple(act_ch), kernel=kernel,
                            dilations=tuple(dilations))
        self.state_tcn = CausalTCN(n_feat, state_ch, kernel, dilations)
        self.state_head = nn.Conv1d(self.state_tcn.c_out, d_latent, 1)
        self.act_tcn = CausalTCN(1, act_ch, kernel, dilations)
        self.act_head = nn.Conv1d(self.act_tcn.c_out, z_dim, 1)
        self.rf = self.state_tcn.rf
        self.d_latent = d_latent
        # The GRU carries a bounded hidden state; the latent is a linear readout
        # of it, so predictions can reach the unbounded N(0, I) targets SIGReg
        # shapes the encoder toward.
        self.h_init = nn.Linear(d_latent, hidden)
        self.gru = nn.GRUCell(z_dim, hidden)
        self.readout = nn.Linear(hidden, d_latent)

    def encode_states(self, feats: torch.Tensor) -> torch.Tensor:
        """Normalized features [B, L, F] -> latents [B, L - rf + 1, D]."""
        return self.state_head(self.state_tcn(feats.transpose(1, 2))).transpose(1, 2)

    def encode_actions(self, u: torch.Tensor) -> torch.Tensor:
        """Bank commands in rad [B, L] -> action codes [B, L - rf + 1, Z]."""
        return self.act_head(self.act_tcn((u / SIGMA_MAX_RAD).unsqueeze(1))).transpose(1, 2)

    def init_hidden(self, s0: torch.Tensor) -> torch.Tensor:
        """Predictor hidden state from the current encoded latent."""
        return torch.tanh(self.h_init(s0))

    def step(self, z_k: torch.Tensor, hid: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """One latent step: returns (next latent, next hidden)."""
        hid = self.gru(z_k, hid)
        return self.readout(hid), hid

    def rollout(self, s0: torch.Tensor, z: torch.Tensor) -> torch.Tensor:
        """Unroll on encoded actions z [B, T, Z] from s0 [B, D] -> [B, T, D]."""
        hid = self.init_hidden(s0)
        out = []
        for k in range(z.shape[1]):
            s, hid = self.step(z[:, k], hid)
            out.append(s)
        return torch.stack(out, 1)


class SIGReg(nn.Module):
    """
    Sketched Isotropic Gaussian Regularization (Balestriero & LeCun, LeJEPA).

    Projects the latents onto M random unit directions and, for each, measures
    the weighted L2 distance between the empirical characteristic function and
    that of N(0, 1) (Epps-Pulley), by trapezoid quadrature on n_knots points.
    By Cramer-Wold, matching every 1-D marginal matches the joint, so this is
    the only anti-collapse term needed: no stop-gradient, EMA target, or decoder.
    """

    def __init__(self, n_proj: int = 32, n_knots: int = 17, t_max: float = 4.0):
        """Precompute the quadrature grid and the Gaussian weighting."""
        super().__init__()
        t = torch.linspace(0.0, t_max, n_knots)
        w = torch.full((n_knots,), t_max / (n_knots - 1))
        w[0] = w[-1] = 0.5 * w[0]
        phi0 = torch.exp(-0.5 * t ** 2)
        self.n_proj = n_proj
        self.register_buffer("t", t)
        self.register_buffer("phi0", phi0)
        # The integrand is even in t, so integrate [0, t_max] and double it.
        self.register_buffer("w", 2.0 * w * phi0)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Mean Epps-Pulley statistic of x [N, D] over the random projections."""
        A = torch.randn(x.shape[1], self.n_proj, device=x.device)
        A = A / A.norm(dim=0, keepdim=True)
        xt = (x @ A).unsqueeze(-1) * self.t                     # [N, M, K]
        err = (torch.cos(xt).mean(0) - self.phi0) ** 2 + torch.sin(xt).mean(0) ** 2
        return (err @ self.w).mean() * x.shape[0]


# ============================================================================
# Physics-inspired prober
# ============================================================================
# Nominal structure the prober corrects. Deliberately crude and parameter-light:
# a single exponential atmosphere and one nominal ballistic coefficient. The
# truth has the NRLMSISE table, a Mach-scheduled CD and L/D, actuator lag, and
# randomized density / L/D / mass; the latent has to supply all of that.
RHO0_KGM3 = 1.225
H_SCALE_M = 7200.0
BETA_NOM = float(constants.CAPSULE_MASS_KG) / (1.5 * float(constants.CAPSULE_REFERENCE_AREA_M2))
LD_NOM = 0.30


def aero_from_corrections(h: torch.Tensor, V: torch.Tensor, corr: torch.Tensor):
    """
    Specific forces from the nominal structure and the latent's corrections.

        a_D   = 0.5 rho_exp(h) V^2 / beta_nom * exp(c_D)
        a_L   = a_D * L/D_nom * exp(c_L)
        a_Lv  = a_L cos(sigma_eff),   a_Lh = a_L sin(sigma_eff)

    corr [..., 3] = (c_D, c_L, sigma_eff). Evaluated at the integrated h and V,
    so drag responds to the rollout's own altitude (the analog of the paper's
    nominal thrust term using the integrated attitude R_t).
    """
    rho = RHO0_KGM3 * torch.exp(-h.clamp_min(0.0) / H_SCALE_M)
    aD = 0.5 * rho * V * V / BETA_NOM * torch.exp(corr[..., 0])
    aL = aD * LD_NOM * torch.exp(corr[..., 1])
    return aD, aL * torch.cos(corr[..., 2]), aL * torch.sin(corr[..., 2])


def integrate_eom(y6: torch.Tensor, corr: torch.Tensor, dt: float, n_sub: int):
    """
    Advance (h, phi, lam, V, gamma, chi) one macro step.

    Same equations as math_3d.nominal_eom_step (chi measured from east toward
    north, non-rotating spherical Earth, inverse-square gravity), with the aero
    terms from aero_from_corrections:
        V_dot     = -a_D - g sin(gamma)
        gamma_dot =  a_Lv / V + (V/r - g/V) cos(gamma)
        chi_dot   = -a_Lh / (V cos(gamma)) - (V/r) cos(gamma) cos(chi) tan(phi)
    The corrections are held over the macro step, the forces are re-evaluated
    every one of the n_sub explicit Euler substeps. The prober is trained
    through this exact integrator, so its corrections also absorb the
    discretization error (the same trick as the paper's prober).

    Returns (y6_next, mean specific forces [..., 3] over the substeps).
    """
    h, phi, lam, V, gam, chi = y6.unbind(-1)
    ds = float(dt) / int(n_sub)
    f_sum = 0.0
    for _ in range(int(n_sub)):
        aD, aLv, aLh = aero_from_corrections(h, V, corr)
        f_sum = f_sum + torch.stack([aD, aLv, aLh], -1)
        r = R_EARTH + h
        g = MU / (r * r)
        sg, cg = torch.sin(gam), torch.cos(gam)
        sc, cc = torch.sin(chi), torch.cos(chi)
        Vs = V.clamp_min(50.0)
        h_n = h + ds * V * sg
        phi_n = phi + ds * V * cg * sc / r
        lam_n = lam + ds * V * cg * cc / (r * torch.cos(phi))
        V_n = V + ds * (-aD - g * sg)
        gam_n = gam + ds * (aLv / Vs + (Vs / r - g / Vs) * cg)
        chi_n = chi + ds * (-aLh / (Vs * cg.clamp_min(0.05)) - (Vs / r) * cg * cc * torch.tan(phi))
        h, phi, lam, V, gam, chi = h_n, phi_n, lam_n, V_n, gam_n, chi_n
    return torch.stack([h, phi, lam, V, gam, chi], -1), f_sum / int(n_sub)


class PIProber(nn.Module):
    """
    Frozen latent -> three corrections to a nominal aero model -> exact kinematics.

    The MLP sees only the latent (paper Eq. 14). Its outputs are a log drag
    factor, a log L/D factor, and a bank offset added to the commanded bank
    (the paper's residual torque is likewise linear in the action, K_t a_t).
    Zero-initialized, so before training the prober is the crude nominal
    physics model with instantaneous bank: the "physics prior only" ablation.
    """

    def __init__(self, d_latent: int, hidden: int = 64):
        """Zero-init the output layer."""
        super().__init__()
        self.hparams = dict(d_latent=d_latent, hidden=hidden)
        self.net = _mlp(d_latent, hidden, 3)
        nn.init.zeros_(self.net[-1].weight)
        nn.init.zeros_(self.net[-1].bias)

    def corrections(self, s: torch.Tensor, u: torch.Tensor) -> torch.Tensor:
        """Latent [..., D] and command [...] -> (c_D, c_L, sigma_eff)."""
        o = self.net(s)
        return torch.stack([o[..., 0].clamp(-5.0, 5.0), o[..., 1].clamp(-3.0, 3.0),
                            u + SIGMA_MAX_RAD * torch.tanh(o[..., 2])], -1)

    def step(self, s: torch.Tensor, y6: torch.Tensor, u: torch.Tensor, dt: float, n_sub: int):
        """Advance the physical state one macro step. Returns (y6_next, mean forces)."""
        return integrate_eom(y6, self.corrections(s, u), dt, n_sub)


class MLPProber(nn.Module):
    """Ablation: latent -> free-form state increment per macro step (no physics)."""

    def __init__(self, d_latent: int, dy_scale6: torch.Tensor, hidden: int = 64):
        """dy_scale6 is the per-dimension std of one macro-step state change."""
        super().__init__()
        self.hparams = dict(d_latent=d_latent, hidden=hidden)
        self.net = _mlp(d_latent, hidden, 6)
        self.register_buffer("scale", dy_scale6.clone().float())

    def step(self, s: torch.Tensor, y6: torch.Tensor, u: torch.Tensor = None,
             dt: float = 0.0, n_sub: int = 1):
        """Advance the physical state one macro step. Returns (y6_next, None)."""
        return y6 + self.net(s) * self.scale, None


def prober_rollout(prober: nn.Module, s_seq: torch.Tensor, y0: torch.Tensor, u_seq: torch.Tensor,
                   dt: float, n_sub: int) -> Tuple[torch.Tensor, Optional[torch.Tensor]]:
    """
    Map predicted latents s~_t+1..t+T [B, T, D] to states y6 [B, T+1, 6].

    The correction for macro step k (t+k -> t+k+1, command u_seq[:, k]) comes
    from s~_t+k+1, the first latent that has seen that command.
    """
    ys, accs, y = [y0], [], y0
    for k in range(s_seq.shape[1]):
        y, a = prober.step(s_seq[:, k], y, u_seq[:, k], dt, n_sub)
        ys.append(y)
        if a is not None:
            accs.append(a)
    return torch.stack(ys, 1), (torch.stack(accs, 1) if accs else None)


def state_errors(yhat6: torch.Tensor, y6: torch.Tensor) -> torch.Tensor:
    """Dimensionless errors [..., 6]: altitude, north, east, speed, gamma, heading."""
    dn = (yhat6[..., 1] - y6[..., 1]) * R_EARTH
    de = (yhat6[..., 2] - y6[..., 2]) * R_EARTH * torch.cos(y6[..., 1])
    dchi = torch.atan2(torch.sin(yhat6[..., 5] - y6[..., 5]), torch.cos(yhat6[..., 5] - y6[..., 5]))
    return torch.stack([(yhat6[..., 0] - y6[..., 0]) / ERR_SCALE["h_m"],
                        dn / ERR_SCALE["pos_m"], de / ERR_SCALE["pos_m"],
                        (yhat6[..., 3] - y6[..., 3]) / ERR_SCALE["V_mps"],
                        (yhat6[..., 4] - y6[..., 4]) / ERR_SCALE["gamma_rad"],
                        dchi / ERR_SCALE["chi_rad"]], -1)


# ============================================================================
# Predictive (autoregressive state-space) baseline
# ============================================================================
class PredictiveModel(nn.Module):
    """History of states and commands -> next macro-step state change (normalized)."""

    def __init__(self, dy_std: torch.Tensor, n_feat: int = N_FEAT, hidden: int = 64,
                 state_ch: Sequence[int] = (32, 32), act_ch: Sequence[int] = (16, 16),
                 kernel: int = 3, dilations: Sequence[int] = (1, 2)):
        """Same encoder sizes as the world model, with a pointwise MLP head."""
        super().__init__()
        self.hparams = dict(n_feat=n_feat, hidden=hidden, state_ch=tuple(state_ch),
                            act_ch=tuple(act_ch), kernel=kernel, dilations=tuple(dilations))
        self.state_tcn = CausalTCN(n_feat, state_ch, kernel, dilations)
        self.act_tcn = CausalTCN(1, act_ch, kernel, dilations)
        self.rf = self.state_tcn.rf
        c = self.state_tcn.c_out + self.act_tcn.c_out
        self.head = nn.Sequential(nn.Conv1d(c, hidden, 1), nn.GELU(),
                                  nn.Conv1d(hidden, hidden, 1), nn.GELU(), nn.Conv1d(hidden, N_Y, 1))
        self.register_buffer("dy_std", dy_std.clone().float())

    def forward(self, feats: torch.Tensor, u: torch.Tensor) -> torch.Tensor:
        """feats [B, L, F], u [B, L] -> normalized state change [B, L - rf + 1, N_Y]."""
        hs = self.state_tcn(feats.transpose(1, 2))
        ha = self.act_tcn((u / SIGMA_MAX_RAD).unsqueeze(1))
        return self.head(torch.cat([hs, ha], 1)).transpose(1, 2)


# ============================================================================
# Batched rollouts (shared by evaluation and MPPI)
# ============================================================================
@torch.no_grad()
def jepa_predict(wm: JEPAWorldModel, prober: nn.Module, feats_hist: torch.Tensor,
                 y0_6: torch.Tensor, u_hist: torch.Tensor, u_plan: torch.Tensor,
                 dt: float, n_sub: int, h_stop: Optional[float] = None) -> Dict[str, torch.Tensor]:
    """
    Roll the JEPA model + prober forward under planned commands.

    Parameters:
        feats_hist: [B or 1, rf, F] normalized feature history (last = now).
        y0_6:       [B, 6] current translational state.
        u_hist:     [B, rf-1] commands of the previous rf-1 macro steps.
        u_plan:     [B, K] planned commands, rad.
        h_stop:     stop early once every sample is below this altitude.

    Returns:
        dict with y6 [B, K'+1, 6], acc [B, K', 3] (PI prober only), and
        latents [B, K', D], where K' <= K if h_stop ended the rollout early.
    """
    B, K = u_plan.shape
    s0 = wm.encode_states(feats_hist)[:, -1].expand(B, -1)
    z = wm.encode_actions(torch.cat([u_hist, u_plan], 1))
    hid = wm.init_hidden(s0)
    y = y0_6
    ys, accs, lats = [y], [], []
    for k in range(K):
        s, hid = wm.step(z[:, k], hid)
        y, a = prober.step(s, y, u_plan[:, k], dt, n_sub)
        ys.append(y)
        lats.append(s)
        if a is not None:
            accs.append(a)
        if h_stop is not None and k % 5 == 4 and bool((y[:, 0] < h_stop).all()):
            break
    return {"y6": torch.stack(ys, 1), "acc": torch.stack(accs, 1) if accs else None,
            "latents": torch.stack(lats, 1)}


@torch.no_grad()
def predictive_predict(pm: PredictiveModel, stats: Stats, y_hist: torch.Tensor,
                       u_hist: torch.Tensor, u_plan: torch.Tensor,
                       h_stop: Optional[float] = None) -> Dict[str, torch.Tensor]:
    """
    Autoregressive rollout of the predictive baseline on its own outputs.

    y_hist [B or 1, rf, N_Y] physical history, u_hist [B, rf-1], u_plan [B, K].
    Returns dict with y [B, K'+1, N_Y] and y6 (its first six columns).
    """
    B, K = u_plan.shape
    yw = y_hist.expand(B, -1, -1)
    uw = u_hist
    ys = [yw[:, -1]]
    for k in range(K):
        uk = torch.cat([uw, u_plan[:, k:k + 1]], 1)
        dy = pm(stats.norm_feats(yw), uk)[:, -1] * pm.dy_std
        y_new = yw[:, -1] + dy
        ys.append(y_new)
        yw = torch.cat([yw[:, 1:], y_new[:, None]], 1)
        uw = uk[:, 1:]
        if h_stop is not None and k % 5 == 4 and bool((y_new[:, 0] < h_stop).all()):
            break
    y = torch.stack(ys, 1)
    return {"y": y, "y6": y[..., :6]}


# ============================================================================
# Losses
# ============================================================================
def jepa_loss(wm: JEPAWorldModel, sigreg: SIGReg, b: Dict[str, torch.Tensor],
              lam_sig: float, sig_rows: int = 1024, target_stopgrad: bool = False):
    """
    Multi-step latent prediction + SIGReg (paper Eq. 9, 12, 13).

    SIGReg is applied to a random subset of the pooled predicted and encoded
    latents (the paper applies it per rollout step; pooling is ~T x cheaper).
    Measured on a laptop CPU, 1024 rows x 32 projections costs ~5 ms per step
    against ~40 ms for 2048 x 64; LeJEPA reports the projection count matters little.
    """
    S = wm.encode_states(b["feats"])                 # [B, T+1, D]
    Z = wm.encode_actions(b["u"])                    # [B, T, Z]
    tgt = S[:, 1:].detach() if target_stopgrad else S[:, 1:]
    pred = wm.rollout(S[:, 0], Z)
    l_pred = (pred - tgt).pow(2).mean()
    lat = torch.cat([pred.reshape(-1, wm.d_latent), S.reshape(-1, wm.d_latent)])
    lat = lat[torch.randperm(lat.shape[0])[:sig_rows]]
    l_sig = sigreg(lat)
    return l_pred + lam_sig * l_sig, {"pred": float(l_pred.detach()), "sigreg": float(l_sig.detach())}


def prober_loss(wm: JEPAWorldModel, prober: nn.Module, b: Dict[str, torch.Tensor],
                dt: float, n_sub: int, w_acc: float = 0.0):
    """
    Supervised rollout loss on frozen latents (paper Eq. 17).

    w_acc > 0 adds an auxiliary match to the measured macro-mean specific
    forces (PI prober only). Not in the paper; it conditions the force outputs
    so the g-load the planner reads off the prober is physical.
    """
    with torch.no_grad():
        S = wm.encode_states(b["feats"])
        pred = wm.rollout(S[:, 0], wm.encode_actions(b["u"]))
    y6 = b["y"][..., :6]
    T = pred.shape[1]
    yhat, acc = prober_rollout(prober, pred, y6[:, 0], b["u"][:, -T:], dt, n_sub)
    l_state = state_errors(yhat[:, 1:], y6[:, 1:]).pow(2).mean()
    logs = {"state": float(l_state.detach())}
    loss = l_state
    if acc is not None and w_acc > 0.0:
        l_acc = ((acc - b["acc"]) / A_SCALE).pow(2).mean()
        loss = loss + w_acc * l_acc
        logs["acc"] = float(l_acc.detach())
    return loss, logs


def predictive_loss(pm: PredictiveModel, b: Dict[str, torch.Tensor]):
    """One-step teacher-forced loss at every position of the window (paper Eq. 4)."""
    out = pm(b["feats"][:, :-1], b["u"])
    dy = (b["y"][:, 1:] - b["y"][:, :-1]) / pm.dy_std
    l = (out - dy).pow(2).mean()
    return l, {"onestep": float(l.detach())}


# ============================================================================
# Training loop
# ============================================================================
def lr_at(step: int, n_steps: int, lr_max: float, lr_min: float, warmup: int) -> float:
    """Linear warmup then cosine decay (the paper's schedule, scaled down)."""
    if step < warmup:
        return lr_max * (step + 1) / warmup
    p = (step - warmup) / max(1, n_steps - warmup)
    return lr_min + 0.5 * (lr_max - lr_min) * (1.0 + math.cos(math.pi * p))


def fit(params, loss_fn: Callable[[], tuple], n_steps: int, val_fn: Optional[Callable[[], float]] = None,
        lr_max: float = 3e-3, lr_min: float = 1e-4, warmup_frac: float = 0.05,
        clip: float = 0.5, weight_decay: float = 1e-5, log_every: int = 200, label: str = "") -> Dict[str, list]:
    """
    Generic AdamW loop with the paper's clip (0.5) and weight decay (1e-5).

    loss_fn() draws its own batch and returns (loss, logs dict). val_fn()
    returns a scalar validation loss. Returns a history of averaged logs.
    """
    params = [p for p in params if p.requires_grad]
    opt = torch.optim.AdamW(params, lr=lr_max, weight_decay=weight_decay)
    warmup = max(1, int(warmup_frac * n_steps))
    hist: Dict[str, list] = {"step": [], "val": []}
    acc: Dict[str, float] = {}
    n_acc, t0 = 0, time.time()
    for step in range(n_steps):
        for g in opt.param_groups:
            g["lr"] = lr_at(step, n_steps, lr_max, lr_min, warmup)
        loss, logs = loss_fn()
        opt.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(params, clip)
        opt.step()
        for k, v in logs.items():
            acc[k] = acc.get(k, 0.0) + v
        n_acc += 1
        if (step + 1) % log_every == 0 or step + 1 == n_steps:
            hist["step"].append(step + 1)
            for k, v in acc.items():
                hist.setdefault(k, []).append(v / n_acc)
            val = float(val_fn()) if val_fn else float("nan")
            hist["val"].append(val)
            msg = " ".join(f"{k}={v / n_acc:.4g}" for k, v in acc.items())
            print(f"  [{label}] step {step + 1:>5}/{n_steps}  {msg}  val={val:.4g}  ({time.time() - t0:.0f}s)")
            acc, n_acc = {}, 0
    return hist


def n_params(m: nn.Module) -> int:
    """Trainable parameter count."""
    return sum(p.numel() for p in m.parameters() if p.requires_grad)


# ============================================================================
# Checkpoints
# ============================================================================
def save_checkpoint(path: str, wm: JEPAWorldModel, pi_prober: PIProber,
                    mlp_prober: Optional[MLPProber], pm: Optional[PredictiveModel],
                    stats: Stats, extra: Optional[dict] = None) -> None:
    """Save every model, its constructor args, and the normalization stats."""
    torch.save({
        "wm": (wm.hparams, wm.state_dict()),
        "pi": (pi_prober.hparams, pi_prober.state_dict()),
        "mlp": (mlp_prober.hparams, mlp_prober.state_dict()) if mlp_prober else None,
        "pred": (pm.hparams, pm.state_dict()) if pm else None,
        "stats": stats.to_dict(), "extra": extra or {},
    }, path)


def load_checkpoint(path: str) -> dict:
    """Rebuild the models saved by save_checkpoint (CPU, eval mode)."""
    ck = torch.load(path, map_location="cpu", weights_only=False)
    stats = Stats(**ck["stats"])
    wm = JEPAWorldModel(**ck["wm"][0])
    wm.load_state_dict(ck["wm"][1])
    pi = PIProber(**ck["pi"][0])
    pi.load_state_dict(ck["pi"][1])
    out = {"wm": wm.eval(), "pi": pi.eval(), "stats": stats, "extra": ck["extra"],
           "mlp": None, "pred": None}
    if ck["mlp"]:
        mlp = MLPProber(dy_scale6=stats.dy_std[:6], **ck["mlp"][0])
        mlp.load_state_dict(ck["mlp"][1])
        out["mlp"] = mlp.eval()
    if ck["pred"]:
        pm = PredictiveModel(dy_std=stats.dy_std, **ck["pred"][0])
        pm.load_state_dict(ck["pred"][1])
        out["pred"] = pm.eval()
    return out
