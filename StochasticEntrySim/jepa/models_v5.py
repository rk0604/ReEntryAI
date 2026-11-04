"""
jepa/models_v5.py

E3 of the v5 study (jepa/v5_plan.md section 4.3): four dynamics models without the JEPA objective, trained on
R0's data, split, normalization and inputs, and flown under v2's MPPI through fly_v4's runner, unchanged.

    mlp_direct   MLP over the flattened 120 s window (25 rows of v2's 15 features, the 24 commands of the
                 history) and the next command's v2 action features -> normalized 5 s change of the six
                 translational states and the step's mean aerodynamic g-load. Rolled out on its own outputs:
                 the window slides over predicted rows (predicted state, the commanded bank as the realized
                 bank, zero roll rate, drag and lift from the predicted g-load split by the last measured L/D)
    gru_direct   GRU encoder (2 layers) over the 25 feature rows -> hidden state; in the rollout a GRU cell
                 driven by the action features; linear head -> normalized state change and g-load
    gru_pi       gru_direct's recurrent core feeding v2's PIProberV2 (the GRU state is the prober's latent)
    r0_nojepa    R0's architecture exactly (WorldModelV2 + PIProberV2, one member), trained end to end on the
                 prober's state loss only: no latent prediction loss, no SIGReg. Saved as a v2 checkpoint, so it
                 flies with fly_v4's V2 adapter like R0

    training     free-running rollouts through the model's own predictions (never teacher forcing), loss =
                 the prober's dimensionless pseudo-Huber state loss (1 km, 5 km, 50 m/s, 0.5 deg, 0.5 deg),
                 masked past the end of the flight, + v2's auxiliary force term at weight 1 (PI models: the
                 macro-mean force vector; direct models: its magnitude = the g-load the planner uses); R0's
                 AdamW (betas 0.9/0.98, weight decay 0.05 on matrices), warmup + cosine, clip 1, AMP; 36,000
                 steps x 256 with R0's stage-1 horizon curriculum 24 -> 48 -> 96 at 20 % and 50 %
    resume       train_state.pt (+ .prev) every cfg.save_every_s seconds, at every curriculum switch and at
                 the end: model, optimizer, scaler, step, history, elapsed time, every RNG state (Python,
                 NumPy, torch CPU and CUDA) and the window sampler's generator; a resumed run continues the
                 same random stream (bit-exact on CPU, smoke check 12). ckpt.pt is written last and marks done
    lr choice    {1, 0.3, 3} x R0's rate, each trained on seed 0 for lr_pilot_frac of the budget; the lowest
                 validation open-loop state loss (horizon 96) wins and is used for every seed
    E3Adapter    the fly_v4 adapter interface (rf, M, context, rollout) for the three non-v2 models
    open loop    footprint error at 400-50 s to go on the test split (as v4 T10), any adapter
"""

from __future__ import annotations

import io
import json
import math
import os
import random
import time
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as Fnn

from jepa import models_v2 as V
from jepa.data import G0, MACRO_DT_S, N_Y, SIGMA_MAX_RAD
from jepa.models import A_SCALE, state_errors

ARMS_E3 = ("mlp_direct", "gru_direct", "gru_pi", "r0_nojepa")
DESCRIPTIONS = {
    "mlp_direct": "MLP over the flattened 120 s window; black-box feedforward dynamics",
    "gru_direct": "GRU encoder + action-driven GRU cell; black-box recurrent dynamics",
    "gru_pi": "GRU core + v2's physics-informed prober; physics decoder without the JEPA objective",
    "r0_nojepa": "R0's architecture trained on the prober's state loss only (no latent loss, no SIGReg)",
}
RF = 25                                  # 120 s at 5 s: R0's receptive field (1 + 3 x 8)
P_HIST = RF - 1
N_OUT_DIRECT = 7                         # six state changes + g-load
FAIL_ENV = "JEPA_V5_FAIL_AT"


# ============================================================================
# Shared persistence helpers (also used by rl_v5 and study_v5)
# ============================================================================
def rng_state(gens: Optional[Dict[str, torch.Generator]] = None) -> dict:
    st = {"python": random.getstate(), "numpy": np.random.get_state(), "torch": torch.get_rng_state()}
    if torch.cuda.is_available():
        st["cuda"] = torch.cuda.get_rng_state_all()
    if gens:
        st["gens"] = {k: g.get_state() for k, g in gens.items()}
    return st


def set_rng_state(st: dict, gens: Optional[Dict[str, torch.Generator]] = None) -> None:
    random.setstate(st["python"])
    np.random.set_state(st["numpy"])
    torch.set_rng_state(st["torch"])
    if "cuda" in st and torch.cuda.is_available():
        torch.cuda.set_rng_state_all(st["cuda"])
    for k, g in (gens or {}).items():
        if k in st.get("gens", {}):
            g.set_state(st["gens"][k])


def fail_point(stage: str, n: int) -> None:
    """The smoke test's simulated disconnect: JEPA_V5_FAIL_AT=<stage>:<n> kills the process at that point."""
    want = os.environ.get(FAIL_ENV, "")
    if want and want == f"{stage}:{int(n)}":
        print(f"[fail_point] {stage}:{n}: simulated disconnect", flush=True)
        os._exit(17)


def save_state(L, obj: dict, path) -> None:
    """Train state to Drive: the previous one becomes .prev, then the new one is written (large-file path)."""
    path = Path(path)
    buf = io.BytesIO()
    torch.save(obj, buf)
    if path.exists():
        prev = path.with_name(path.stem + ".prev" + path.suffix)
        L.guard(prev)
        os.replace(path, prev)
    if hasattr(L, "save_large"):
        L.save_large(buf.getvalue(), path)
    else:
        L.atomic_bytes(path, buf.getvalue())


def load_state(path, log=print) -> Optional[dict]:
    """The newest train state that loads cleanly (else .prev, logged); None if neither exists or loads."""
    path = Path(path)
    for p in (path, path.with_name(path.stem + ".prev" + path.suffix)):
        if p.exists():
            try:
                return torch.load(p, map_location="cpu", weights_only=False)
            except Exception as e:                       # truncated or corrupt: fall back
                log(f"  [resume] {p.name} does not load ({type(e).__name__}); trying the previous state")
    return None


# ============================================================================
# Models
# ============================================================================
class _Norm(nn.Module):
    """v2's normalization statistics as buffers (the checkpoint is self-contained)."""

    def __init__(self, stats: Optional[V.StatsV2] = None, n_feat: int = V.N_FEAT_V2):
        super().__init__()
        self.register_buffer("feat_mean", torch.zeros(n_feat) if stats is None else stats.feat_mean.clone())
        self.register_buffer("feat_std", torch.ones(n_feat) if stats is None else stats.feat_std.clone())
        self.register_buffer("dy_std", torch.ones(N_Y) if stats is None else stats.dy_std.clone())

    def feats(self, y: torch.Tensor) -> torch.Tensor:
        return (V.raw_features_v2(y) - self.feat_mean) / self.feat_std


def _mlp(d_in: int, width: int, depth: int, d_out: int) -> nn.Sequential:
    layers: List[nn.Module] = []
    for i in range(depth):
        layers += [nn.Linear(d_in if i == 0 else width, width), nn.GELU()]
    head = nn.Linear(width, d_out)
    nn.init.zeros_(head.weight)
    nn.init.zeros_(head.bias)
    layers.append(head)
    return nn.Sequential(*layers)


class MLPDirect(nn.Module):
    """Flattened-window MLP; the window slides over its own predicted rows."""

    def __init__(self, stats: Optional[V.StatsV2] = None, width: int = 832, depth: int = 4, rf: int = RF):
        super().__init__()
        self.hparams = {"arm": "mlp_direct", "width": width, "depth": depth, "rf": rf}
        self.rf, self.P = rf, rf - 1
        self.norm = _Norm(stats)
        d_in = rf * V.N_FEAT_V2 + self.P * 3 + 3
        self.net = _mlp(d_in, width, depth, N_OUT_DIRECT)

    def context(self, feats_hist: torch.Tensor, y_now: torch.Tensor) -> torch.Tensor:
        """feats [B, rf, F] (normalized), y_now [B, N_Y] -> [B, rf F + 2] (window, realized bank, L/D)."""
        ld = y_now[:, 9].clamp_min(V.A_FLOOR) / y_now[:, 8].clamp_min(V.A_FLOOR)
        return torch.cat([feats_hist.flatten(1), y_now[:, 6:7], ld[:, None].clamp(0.0, 2.0)], 1)

    def rollout(self, ctx, y0, uh, u, n_sub=2, h_stop=None):
        B, K = u.shape
        F_ = V.N_FEAT_V2
        win = ctx[:, :self.rf * F_].view(B, self.rf, F_)
        ld = ctx[:, -1]
        uw = uh
        y = y0.float()
        ys, gs = [y], []
        for k in range(K):
            x = torch.cat([win.flatten(1), V.action_features(uw).flatten(1), V.action_features(u[:, k])], 1)
            o = V._full(self.net(x))
            y = y + o[:, :6].float() * self.norm.dy_std[:6]
            g = Fnn.softplus(o[:, 6].float())
            a = g * G0
            aD = a / torch.sqrt(1.0 + ld * ld)
            row = torch.cat([y, u[:, k:k + 1], torch.zeros_like(y[:, :1]), aD[:, None], (aD * ld)[:, None]], 1)
            win = torch.cat([win[:, 1:], self.norm.feats(row)[:, None].to(win.dtype)], 1)
            uw = torch.cat([uw[:, 1:], u[:, k:k + 1]], 1)
            ys.append(y)
            gs.append(g)
            if h_stop is not None and k % 5 == 4 and bool((y[:, 0] < h_stop).all()):
                break
        return torch.stack(ys, 1), torch.stack(gs, 1), None


class GRUCore(nn.Module):
    def __init__(self, hidden: int, layers: int = 2):
        super().__init__()
        self.enc = nn.GRU(V.N_FEAT_V2, hidden, num_layers=layers, batch_first=True)
        self.cell = nn.GRUCell(3, hidden)

    def encode(self, feats_hist: torch.Tensor) -> torch.Tensor:
        _, h = self.enc(feats_hist)
        return h[-1]

    def step(self, u_k: torch.Tensor, s: torch.Tensor) -> torch.Tensor:
        return self.cell(V.action_features(u_k).to(s.dtype), s)


class GRUDirect(nn.Module):
    def __init__(self, stats: Optional[V.StatsV2] = None, hidden: int = 452, layers: int = 2, rf: int = RF):
        super().__init__()
        self.hparams = {"arm": "gru_direct", "hidden": hidden, "layers": layers, "rf": rf}
        self.rf, self.P = rf, rf - 1
        self.norm = _Norm(stats)
        self.core = GRUCore(hidden, layers)
        self.head = nn.Linear(hidden, N_OUT_DIRECT)
        nn.init.zeros_(self.head.weight)
        nn.init.zeros_(self.head.bias)

    def context(self, feats_hist, y_now):
        return self.core.encode(feats_hist)

    def rollout(self, ctx, y0, uh, u, n_sub=2, h_stop=None):
        s, y = ctx, y0.float()
        ys, gs = [y], []
        for k in range(u.shape[1]):
            s = self.core.step(u[:, k], s)
            o = V._full(self.head(s)).float()
            y = y + o[:, :6] * self.norm.dy_std[:6]
            ys.append(y)
            gs.append(Fnn.softplus(o[:, 6]))
            if h_stop is not None and k % 5 == 4 and bool((y[:, 0] < h_stop).all()):
                break
        return torch.stack(ys, 1), torch.stack(gs, 1), None


class GRUPI(nn.Module):
    def __init__(self, stats: Optional[V.StatsV2] = None, hidden: int = 447, layers: int = 2, prober_hidden: int = 128,
                 rf: int = RF):
        super().__init__()
        self.hparams = {"arm": "gru_pi", "hidden": hidden, "layers": layers, "prober_hidden": prober_hidden, "rf": rf}
        self.rf, self.P = rf, rf - 1
        self.norm = _Norm(stats)
        self.core = GRUCore(hidden, layers)
        self.pi = V.PIProberV2(1, hidden, hidden=prober_hidden, state_aware=True)

    def context(self, feats_hist, y_now):
        return self.core.encode(feats_hist)

    def rollout(self, ctx, y0, uh, u, n_sub=2, h_stop=None):
        s, y = ctx, y0.float()[None]
        ys, accs = [y[0]], []
        for k in range(u.shape[1]):
            s = self.core.step(u[:, k], s)
            with torch.autocast(s.device.type, enabled=False):
                y, a = self.pi.step(V._full(s).float()[None], y, u[None, :, k].float(), MACRO_DT_S, n_sub)
            ys.append(y[0])
            accs.append(a[0])
            if h_stop is not None and k % 5 == 4 and bool((y[0, :, 0] < h_stop).all()):
                break
        acc = torch.stack(accs, 1)
        return torch.stack(ys, 1), acc.norm(dim=-1) / G0, acc


class R0NoJEPA(nn.Module):
    """R0 (WorldModelV2 + PIProberV2, one member) as one module; trained end to end on the state loss."""

    def __init__(self, stats: Optional[V.StatsV2] = None, arch: Optional[dict] = None, prober_hidden: int = 128):
        super().__init__()
        arch = dict(arch or {})
        self.hparams = {"arm": "r0_nojepa", "arch": arch, "prober_hidden": prober_hidden}
        self.norm = _Norm(stats)
        self.wm = V.WorldModelV2(n_members=1, **arch)
        self.pi = V.PIProberV2(1, self.wm.d_latent, hidden=prober_hidden, state_aware=True)
        self.rf, self.P = self.wm.rf, self.wm.rf - 1

    def context(self, feats_hist, y_now):
        return self.wm.encode_states(feats_hist[None])[0, :, -1].float()

    def rollout(self, ctx, y0, uh, u, n_sub=2, h_stop=None):
        lat = self.wm.predict(ctx[None], self.wm.action_tokens(torch.cat([uh, u], 1)[None])).float()
        with torch.autocast(lat.device.type, enabled=False):
            y6, acc = V.prober_rollout_v2(self.pi, lat, y0.float()[None], u.float()[None], MACRO_DT_S, n_sub, h_stop)
        return y6[0], acc[0].norm(dim=-1) / G0, acc[0]


def n_params(m: nn.Module) -> int:
    return int(sum(p.numel() for p in m.parameters()))


def r0_param_count(arch: Optional[dict] = None) -> int:
    """R0's parameters per member (world model + prober): 2,454,531 at v2's defaults."""
    wm = V.WorldModelV2(n_members=1, **(arch or {}))
    return n_params(wm) + n_params(V.PIProberV2(1, wm.d_latent, hidden=128, state_aware=True))


def _size(arm: str, target: int) -> dict:
    """The width (MLP) or hidden size (GRU) whose parameter count is closest to target."""
    lo, hi = 8, 4096
    build = {"mlp_direct": lambda w: MLPDirect(width=w), "gru_direct": lambda w: GRUDirect(hidden=w),
             "gru_pi": lambda w: GRUPI(hidden=w)}[arm]
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if n_params(build(mid)) < target:
            lo = mid
        else:
            hi = mid
    w = min((lo, hi), key=lambda w: abs(n_params(build(w)) - target))
    return {"width": w} if arm == "mlp_direct" else {"hidden": w}


def build_e3(arm: str, stats: Optional[V.StatsV2] = None, sizes: Optional[dict] = None,
             r0_arch: Optional[dict] = None) -> nn.Module:
    """One E3 model; sizes = {arm: {width|hidden: n}} (None: sized to R0's parameter count within 10 %)."""
    if arm == "r0_nojepa":
        return R0NoJEPA(stats, r0_arch)
    kw = dict((sizes or {}).get(arm) or _size(arm, r0_param_count(r0_arch)))
    return {"mlp_direct": MLPDirect, "gru_direct": GRUDirect, "gru_pi": GRUPI}[arm](stats, **kw)


def rebuild(hp: dict) -> nn.Module:
    hp = dict(hp)
    arm = hp.pop("arm")
    if arm == "r0_nojepa":
        return R0NoJEPA(None, hp["arch"], hp["prober_hidden"])
    return {"mlp_direct": MLPDirect, "gru_direct": GRUDirect, "gru_pi": GRUPI}[arm](None, **hp)


# ============================================================================
# Loss, training
# ============================================================================
@dataclass
class E3TrainCfg:
    steps: int = 36_000
    batch: int = 256
    lr: float = 6e-4                         # R0's stage-1 rate; x {1, 0.3, 3} chosen on validation
    lr_min: float = 2e-5
    warmup: int = 1_000
    weight_decay: float = 0.05
    clip: float = 1.0
    curriculum: Tuple[Tuple[float, int], ...] = ((0.0, 24), (0.2, 48), (0.5, 96))
    log_every: int = 500
    save_every_s: float = 300.0
    save_every_steps: int = 0                # > 0: also save every n steps (the smoke test's deterministic resume)
    amp: bool = True
    w_acc: float = 1.0
    n_sub: int = 2
    val_rows: int = 256


def _cfg_v2(cfg: E3TrainCfg) -> V.TrainCfgV2:
    return V.TrainCfgV2(steps=cfg.steps, batch=cfg.batch, lr=cfg.lr, lr_min=cfg.lr_min, warmup=cfg.warmup,
                        weight_decay=cfg.weight_decay, clip=cfg.clip, curriculum=tuple(cfg.curriculum))


def e3_loss(model: nn.Module, b: Dict[str, torch.Tensor], w_acc: float, n_sub: int):
    """Free-running rollout from the window's current state; masked pseudo-Huber state loss + force term."""
    P = model.P
    feats, u, y, acc, m = (b[k][0] for k in ("feats", "u", "y", "acc", "mask"))   # drop the member dim
    ctx = model.context(feats[:, :P + 1], y[:, 0])
    y6, g, a3 = model.rollout(ctx, y[:, 0, :6], u[:, :P], u[:, P:], n_sub, None)
    e = V.pseudo_huber(state_errors(y6[:, 1:], y[:, 1:, :6])).mean(-1)
    den = m.sum().clamp_min(1.0)
    l_state = (e * m).sum() / den
    if a3 is not None:
        l_acc = ((((a3 - acc) / A_SCALE) ** 2).mean(-1) * m).sum() / den
    else:
        l_acc = ((((g * G0 - acc.norm(dim=-1)) / A_SCALE) ** 2) * m).sum() / den
    return l_state + w_acc * l_acc, {"state": float(l_state.detach()), "acc": float(l_acc.detach())}


def param_groups(module: nn.Module, wd: float):
    """Weight decay on matrices (Linear/GRU weights, ensemble weights), not on biases or norms."""
    decay, no = [], []
    for n, p in module.named_parameters():
        if not p.requires_grad:
            continue
        mat = (p.dim() >= 3) or (p.dim() == 2 and not n.startswith(("wm.", "pi.")) and "bias" not in n)
        (decay if mat and "norm" not in n and "ln" not in n else no).append(p)
    return [{"params": decay, "weight_decay": wd}, {"params": no, "weight_decay": 0.0}]


def train_e3(model: nn.Module, ds_tr, ds_va, cfg: E3TrainCfg, device, L, state_path, seed: int,
             label: str = "", log=print, stage: str = "e3_train") -> Dict[str, list]:
    """
    Resumable AdamW loop (models_v2.train_stage's arithmetic) with the full train state of v5_plan.md 9.1.
    Returns the history (step, T, sec, state, acc, val).
    """
    device = torch.device(device)
    model.to(device)
    c2 = _cfg_v2(cfg)
    opt = torch.optim.AdamW(param_groups(model, cfg.weight_decay), lr=cfg.lr, betas=(0.9, 0.98))
    dt_amp = V.amp_dtype(device) if cfg.amp else None
    scaler = torch.amp.GradScaler(device.type, enabled=dt_amp == torch.float16)
    gen = torch.Generator(device=ds_tr.device).manual_seed(int(seed) * 7919 + 11)
    gv = torch.Generator(device=ds_va.device).manual_seed(1)
    vbs = {T: ds_va.sample((1, cfg.val_rows), T, gv) for _, T in cfg.curriculum}
    hist: Dict[str, list] = {"step": [], "T": [], "sec": [], "lr": [], "state": [], "acc": [], "val": []}
    start, elapsed = 0, 0.0
    acc0, n_acc0 = {}, 0
    torch.manual_seed(int(seed))
    st = load_state(state_path, log) if state_path else None
    if st is not None:
        model.load_state_dict(st["module"])
        opt.load_state_dict(st["opt"])
        scaler.load_state_dict(st["scaler"])
        start, hist, elapsed = int(st["step"]), st["hist"], float(st.get("elapsed", 0.0))
        acc0, n_acc0 = dict(st.get("acc", {})), int(st.get("n_acc", 0))
        set_rng_state(st["rng"], {"gen": gen})
        log(f"  [{label}] resumed at step {start}/{cfg.steps}")

    def save(step: int):
        if state_path:
            save_state(L, {"module": model.state_dict(), "opt": opt.state_dict(), "scaler": scaler.state_dict(),
                           "step": step, "hist": hist, "elapsed": elapsed, "acc": acc, "n_acc": n_acc,
                           "rng": rng_state({"gen": gen}), "cfg": asdict(cfg)}, state_path)

    def val_at(T: int) -> float:
        model.eval()
        with torch.no_grad(), torch.autocast(device.type, dtype=dt_amp, enabled=dt_amp is not None):
            v = e3_loss(model, vbs[T], 0.0, cfg.n_sub)[1]["state"]
        model.train()
        return v

    params = [p for g in opt.param_groups for p in g["params"]]
    acc: Dict[str, float] = acc0               # the log window's running sums (saved, so a resume logs the same)
    n_acc, t_log, t_save = n_acc0, time.time(), time.time()
    model.train()
    for step in range(start, cfg.steps):
        lr = V._lr(step, c2)
        for g in opt.param_groups:
            g["lr"] = lr
        T = V.horizon_at(step, c2)
        b = ds_tr.sample((1, cfg.batch), T, gen)
        with torch.autocast(device.type, dtype=dt_amp, enabled=dt_amp is not None):
            loss, logs = e3_loss(model, b, cfg.w_acc, cfg.n_sub)
        opt.zero_grad(set_to_none=True)
        scaler.scale(loss).backward()
        scaler.unscale_(opt)
        torch.nn.utils.clip_grad_norm_(params, cfg.clip)
        scaler.step(opt)
        scaler.update()
        for k, v in logs.items():
            acc[k] = acc.get(k, 0.0) + v
        n_acc += 1
        last = step + 1 == cfg.steps
        if (step + 1) % cfg.log_every == 0 or last:
            now = time.time()
            elapsed += now - t_log
            per = (now - t_log) / max(n_acc, 1)
            t_log = now
            hist["step"].append(step + 1)
            hist["T"].append(T)
            hist["sec"].append(round(elapsed, 1))
            hist["lr"].append(lr)
            for k in ("state", "acc"):
                hist[k].append(acc.get(k, 0.0) / max(n_acc, 1))
            hist["val"].append(val_at(T))
            log(f"  [{label}] step {step + 1:>6}/{cfg.steps} T={T:<3} state={hist['state'][-1]:.4g} "
                f"acc={hist['acc'][-1]:.4g} val={hist['val'][-1]:.4g} ({per * 1000:.0f} ms/step, "
                f"~{per * (cfg.steps - step - 1) / 60:.0f} min left)")
            acc, n_acc = {}, 0
        switch = not last and V.horizon_at(step + 1, c2) != T
        due = state_path and (time.time() - t_save > cfg.save_every_s
                              or (cfg.save_every_steps and (step + 1) % cfg.save_every_steps == 0))
        if last or switch or due:
            save(step + 1)
            t_save = time.time()
            fail_point(stage, step + 1)
    model.eval()
    return hist


def val_open_loop(model: nn.Module, ds_va, T: int, rows: int = 512, n_sub: int = 2) -> float:
    """Validation open-loop state loss at horizon T (the learning-rate choice)."""
    g = torch.Generator(device=ds_va.device).manual_seed(3)
    b = ds_va.sample((1, rows), T, g)
    model.eval()
    with torch.no_grad():
        return float(e3_loss(model, b, 0.0, n_sub)[1]["state"])


# ============================================================================
# Checkpoints and adapters
# ============================================================================
def save_e3(L, path, model: nn.Module, stats: V.StatsV2, extra: dict) -> None:
    """r0_nojepa as a models_v2 checkpoint (flies with the V2 adapter); the others as an e3 checkpoint."""
    path = Path(path)
    if isinstance(model, R0NoJEPA):
        buf = io.BytesIO()
        torch.save({"version": 2, "wm": (model.wm.hparams, {k: v.cpu() for k, v in model.wm.state_dict().items()}),
                    "pi": (model.pi.hparams, {k: v.cpu() for k, v in model.pi.state_dict().items()}),
                    "stats": stats.to_dict(), "extra": {**extra, "e3_hparams": model.hparams}}, buf)
    else:
        buf = io.BytesIO()
        torch.save({"version": "e3", "hparams": model.hparams,
                    "state": {k: v.cpu() for k, v in model.state_dict().items()},
                    "stats": stats.to_dict(), "extra": extra}, buf)
    L.atomic_bytes(path, buf.getvalue())


def load_e3(path, device="cpu") -> Tuple[nn.Module, V.StatsV2, dict]:
    ck = torch.load(path, map_location="cpu", weights_only=False)
    stats = V.StatsV2(**ck["stats"])
    if ck["version"] == 2:
        hp = ck["extra"].get("e3_hparams") or {"arm": "r0_nojepa", "arch": ck["wm"][0], "prober_hidden": 128}
        arch = {k: v for k, v in ck["wm"][0].items() if k != "n_members"}
        m = R0NoJEPA(stats, arch, ck["pi"][0]["hidden"])
        m.wm.load_state_dict(ck["wm"][1])
        m.pi.load_state_dict(ck["pi"][1])
    else:
        m = rebuild(ck["hparams"])
        m.load_state_dict(ck["state"])
    return m.to(device).eval(), stats, ck["extra"]


class E3Adapter:
    """fly_v4's adapter interface for an E3 model (one member)."""
    M = 1

    def __init__(self, model: nn.Module, device="cpu"):
        self.dev = torch.device(device)
        self.m = model.to(self.dev).eval()
        self.rf = model.rf
        self.amp = V.amp_dtype(self.dev)

    @torch.no_grad()
    def context(self, y_hists, jobs=None) -> torch.Tensor:
        Y = torch.stack([torch.as_tensor(y) for y in y_hists]).float().to(self.dev)
        with torch.autocast(self.dev.type, dtype=self.amp, enabled=self.amp is not None):
            return self.m.context(self.m.norm.feats(Y), Y[:, -1])[None].float()

    @torch.no_grad()
    def rollout(self, s0, y0, uh, u, n_sub, h_stop):
        with torch.autocast(self.dev.type, dtype=self.amp, enabled=self.amp is not None):
            y6, g, _ = self.m.rollout(s0[0].to(self.dev), y0.to(self.dev), uh.to(self.dev), u.to(self.dev), n_sub,
                                      h_stop)
        return y6.float()[:, None].cpu(), g.float()[:, None].cpu()


def load_adapter(kind: str, ckpt: str, device="cpu"):
    """E3 model kinds for fly_v5.load_adapter_v5 ("e3"; r0_nojepa flies as kind "v2")."""
    if kind != "e3":
        raise ValueError(f"unknown model kind {kind!r}")
    m, _, _ = load_e3(ckpt, device)
    return E3Adapter(m, device)


# ============================================================================
# One arm, one seed (resumable; ckpt.pt marks done)
# ============================================================================
def train_arm_seed(cfg5, data: dict, L, arm: str, seed: int, device, lr_mult: float, log=print,
                   train_dir: Optional[Path] = None, steps: Optional[int] = None) -> dict:
    """
    Train (or resume, or skip) one E3 model. cfg5: study_v5.StudyConfigV5 (e3_train, e3_sizes, r0_arch).
    train_dir/steps: the learning-rate pilot's folder and budget. Returns the run's summary.
    """
    d = Path(train_dir) if train_dir is not None else L.model_dir(arm, seed)
    ckpt = d / "ckpt.pt"
    if ckpt.exists():
        info = d / "model.json"
        return json.loads(info.read_text()) if info.exists() else {"arm": arm, "seed": seed, "ckpt": str(ckpt)}
    L.guard(d).mkdir(parents=True, exist_ok=True)
    tc = E3TrainCfg(**cfg5.e3_train)
    tc = replace(tc, lr=tc.lr * float(lr_mult), steps=int(steps or tc.steps))
    torch.manual_seed(int(seed))
    model = build_e3(arm, data["stats"], cfg5.e3_sizes, cfg5.r0_arch)
    ds_tr, ds_va = data["ds_tr"], data["ds_va"]
    if model.rf != ds_tr.P + 1:
        from jepa.ablate_v4 import SlicedWindows
        ds_tr, ds_va = SlicedWindows(ds_tr, model.rf), SlicedWindows(ds_va, model.rf)
    t0 = time.time()
    hist = train_e3(model, ds_tr, ds_va, tc, device, L, d / "train_state.pt", seed,
                    label=f"{arm} seed {seed}" + (f" lr x{lr_mult:g}" if train_dir is not None else ""), log=log)
    T_max = max(T for _, T in tc.curriculum)
    info = {"arm": arm, "seed": int(seed), "lr_mult": float(lr_mult), "lr": tc.lr, "steps": tc.steps,
            "params": n_params(model), "r0_params": r0_param_count(cfg5.r0_arch), "hparams": model.hparams,
            "val_open_loop_T%d" % T_max: val_open_loop(model, ds_va, T_max, n_sub=tc.n_sub),
            "train_s": round(time.time() - t0, 1), "train_cfg": asdict(tc), "ckpt": str(ckpt)}
    import pandas as pd
    L.write_csv(pd.DataFrame({k: v for k, v in hist.items()}), d / "train_log.csv")
    _curves(L, hist, f"{arm} seed {seed}", d / "F_curves.png")
    save_e3(L, ckpt, model, data["stats"], info)
    L.write_json(d / "model.json", info)
    for f in ("train_state.pt", "train_state.prev.pt"):
        if (d / f).exists():
            L.guard(d / f).unlink()
    return info


def _curves(L, hist: dict, title: str, path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not hist.get("step"):
        return
    fig, ax = plt.subplots(1, 1, figsize=(6, 3.4))
    ax.semilogy(hist["step"], hist["state"], label="train state loss")
    ax.semilogy(hist["step"], hist["val"], label="val state loss")
    ax.semilogy(hist["step"], hist["acc"], label="train force term", alpha=0.6)
    ax.set_xlabel("step")
    ax.set_title(title)
    ax.legend(fontsize=7)
    for ext in (".png", ".pdf"):
        L.save_fig(fig, Path(path).with_suffix(ext), dpi=200)
    plt.close(fig)


def choose_lr(cfg5, data: dict, L, arm: str, device, log=print) -> float:
    """The learning-rate multiplier of an arm: pilots on seed 0, chosen on validation open-loop loss."""
    out = L.models / arm / "lr_choice.json"
    if out.exists():
        return float(json.loads(out.read_text())["lr_mult"])
    steps = max(1, int(round(cfg5.lr_pilot_frac * E3TrainCfg(**cfg5.e3_train).steps)))
    T_max = max(T for _, T in E3TrainCfg(**cfg5.e3_train).curriculum)
    rows = []
    for mult in cfg5.lr_mults:
        info = train_arm_seed(cfg5, data, L, arm, 0, device, mult, log, L.models / arm / "lr_pilot" / f"x{mult:g}",
                              steps)
        rows.append({"lr_mult": mult, "val": info["val_open_loop_T%d" % T_max]})
    best = min(rows, key=lambda r: (r["val"] if math.isfinite(r["val"]) else float("inf")))
    L.write_json(out, {"arm": arm, "lr_mult": best["lr_mult"], "pilots": rows, "pilot_steps": steps,
                       "criterion": f"validation open-loop state loss at T={T_max}"})
    log(f"  [{arm}] learning rate x{best['lr_mult']:g} (pilots: {rows})")
    return float(best["lr_mult"])


# ============================================================================
# Open loop (v4 T10's footprint errors, any adapter)
# ============================================================================
@torch.no_grad()
def footprint_rows(adapter, name: str, eps: Sequence[dict], ttg_s: Sequence[float], h_cross: float = 7600.0,
                   dt: float = MACRO_DT_S, n_sub: int = 2, margin_steps: int = 30) -> List[dict]:
    """Predicted vs true drogue point from several times to go (evaluate_v2.footprint_points' windows)."""
    from jepa.evaluate import crossing_point
    from jepa.evaluate_v2 import _decompose
    P = adapter.rf - 1
    eps = [e for e in eps if e["outcome"] == "drogue"]
    rows: List[dict] = []
    for ttg in ttg_s:
        k_back = int(round(ttg / dt))
        sel = [e for e in eps if len(e["u"]) - k_back >= 0]
        if not sel:
            continue
        K = k_back + margin_steps
        yh, y0, uh, up, truth = [], [], [], [], []
        for e in sel:
            t = len(e["u"]) - k_back
            y = torch.from_numpy(e["y"]).float()
            u = torch.from_numpy(e["u"]).float()
            yp = torch.cat([y[:1].repeat(P, 1), y])
            upad = torch.cat([torch.zeros(P), u])
            yh.append(yp[t:t + P + 1])
            y0.append(y[t, :6])
            uh.append(upad[t:t + P])
            fut = u[t:]
            fut = torch.cat([fut, (fut[-1:] if len(fut) else u[-1:]).repeat(K - len(fut))])
            up.append(fut[:K])
            truth.append(e["y_final"])
        ctx = adapter.context(yh)
        y6, _ = adapter.rollout(ctx, torch.stack(y0), torch.stack(uh), torch.stack(up), n_sub, None)
        y6 = y6[:, 0].double()
        phi, lam, _, _ = crossing_point(y6, h_cross)
        err, along, cross = _decompose(phi.numpy(), lam.numpy(), np.stack(truth).astype(np.float64))
        for i, e in enumerate(sel):
            rows.append({"model": name, "ttg_s": float(ttg), "err_km": float(err[i]), "along_km": float(along[i]),
                         "cross_km": float(cross[i]), **{f"dom_{k}": v for k, v in e["domain"].items()}})
    return rows


def open_loop_seed(cfg5, L, data: dict, arm: str, seed: int, device="cpu"):
    """03_open_loop/footprint_<arm>_seed<k>.csv (test and validation splits); cached."""
    import pandas as pd
    path = L.open_loop / f"footprint_{arm}_seed{seed}.csv"
    if path.exists():
        return pd.read_csv(path)
    ck = L.ckpt(arm, seed)
    if not ck.exists():
        return None
    m, _, _ = load_e3(str(ck), device)
    ad = E3Adapter(m, device)
    rows = []
    for split in ("test",):
        r = pd.DataFrame(footprint_rows(ad, arm, data[split], cfg5.ttg_s))
        r["split"] = split
        rows.append(r)
    df = pd.concat(rows, ignore_index=True)
    df["seed"] = seed
    L.write_csv(df, path)
    return df
