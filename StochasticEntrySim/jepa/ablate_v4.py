"""
jepa/ablate_v4.py

E3 of the v4 study (jepa/v4_plan.md section 4.3): the reference R0 (v2's recipe
with one member) and the 13 single-change variants, as model and training
configurations, and their resumable training.

    variant_spec       name -> the full configuration (architecture, prober,
                       loss, curricula, budget); what differs from the reference
                       is exactly the one change in the plan's tables
    param_count        parameters per member (world model + prober), without training
    SlicedWindows      the shared GPU training windows seen with a shorter
                       encoder history (A2's 35 s receptive field)
    train_stage_v4     models_v2.train_stage (same arithmetic) whose train state
                       also holds the RNG state and is written atomically through
                       the guarded layout, so a resumed run continues the same
                       random stream
    train_variant_seed one variant and seed: train (or resume), checkpoint, log,
                       curves, latent health, variant.json; ckpt.pt marks it done

v2 family (R0, A1-A8): models_v2.WorldModelV2 / PIProberV2 / wm_loss_v2 with the
v4 options. v1 family (B1-B3): models.JEPAWorldModel / PIProber trained with
study_v3.fit_resumable on v3's fair budget (v1's recipe, convergence rule).
"""

from __future__ import annotations

import io
import json
import math
import os
import time
from dataclasses import asdict
from typing import Callable, Dict, Optional, Sequence, Tuple

import numpy as np
import torch

from jepa import data as D
from jepa import models as M
from jepa import models_v2 as V
from jepa import study_v3 as S

V1_FEATS = tuple(range(len(D.FEAT_NAMES)))                                        # A1: v1's 10 inputs
NO_MEASURED_AERO = tuple(i for i, n in enumerate(V.FEAT_NAMES_V2) if n not in ("cD_meas", "cL_meas"))   # A1b

DESCRIPTIONS = {
    "R0": "v2's recipe with one member (reference)",
    "A1_inputs_v1": "encoder inputs = v1's 10 features (no sin/cos bank, log a_D, measured aero)",
    "A1b_no_measured_aero": "no measured-aero inputs (cD_meas, cL_meas)",
    "A2_memory_35s": "encoder receptive field 120 s -> 35 s (band 3)",
    "A3_blind_prober": "prober sees the latent only, not the rollout's h and V",
    "A4_train_T24": "v1's training horizons (stage 1 T = 24, stage 2 T = 16)",
    "A5_sigreg_pooled": "v1's SIGReg: pooled over time, lambda 0.02",
    "A6_latent_32": "latent 64 -> 32",
    "A7a_size_250k": "about 250k parameters (widths scaled ~10x down)",
    "A7b_size_33k": "about 33k parameters (v1's size)",
    "A8_actions_v1": "v1's action input: 7-command window, u / sigma_max only",
    "B1_v1_plus_inputs": "v1 + the five v2 input features",
    "B2_v1_plus_state_prober": "v1's prober also sees the rollout's h and V",
    "B3_v1_plus_memory_120s": "v1's TCNs with a 120 s receptive field",
}


def family(name: str) -> str:
    return "v1" if name.startswith("B") else "v2"


# ============================================================================
# Variant definitions
# ============================================================================
def variant_spec(cfg, name: str) -> dict:
    """The full model and training configuration of one variant (cfg: study_v4.StudyConfigV4)."""
    if name not in DESCRIPTIONS:
        raise KeyError(name)
    if family(name) == "v1":
        v3 = cfg.v3_training()
        wm = {}
        feat_set, state_aware = "v1", False
        if name == "B1_v1_plus_inputs":
            wm["n_feat"], feat_set = V.N_FEAT_V2, "v2"
        elif name == "B2_v1_plus_state_prober":
            state_aware = True
        elif name == "B3_v1_plus_memory_120s":
            wm.update(state_ch=(32, 32, 32, 32), act_ch=(16, 16, 16, 16), dilations=(1, 2, 3, 6))
        return {"name": name, "family": "v1", "change": DESCRIPTIONS[name], "wm": wm,
                "prober": {"state_aware": state_aware}, "feat_set": feat_set, "lambda": v3.v1_lambda,
                "T": v3.v1_T, "T_pi": v3.v1_T_pi,
                "budget": {k: getattr(v3, k) for k in ("v1_steps", "v1_batch", "v1_pi_steps", "v1_pi_batch", "v1_lr",
                                                       "v1_lr_min", "v1_warmup_frac", "v1_clip", "v1_weight_decay",
                                                       "v1_extend_tol", "v1_extend_frac", "v1_extend_lr_frac")}}
    arch = {"n_members": 1, **cfg.r0_arch}
    base = V.WorldModelV2(n_members=1, **cfg.r0_arch).hparams
    prober = {"hidden": 128, "state_aware": True}
    wm_train, pi_train = dict(cfg.r0_wm), dict(cfg.r0_pi)
    lam, mode = cfg.r0_lambda, "per_step"
    if name == "A1_inputs_v1":
        arch["feat_mask"] = V1_FEATS
    elif name == "A1b_no_measured_aero":
        arch["feat_mask"] = NO_MEASURED_AERO
    elif name == "A2_memory_35s":
        arch["enc_band"] = 3
    elif name == "A3_blind_prober":
        prober["state_aware"] = False
    elif name == "A4_train_T24":
        wm_train["curriculum"], pi_train["curriculum"] = cfg.a4_curricula
    elif name == "A5_sigreg_pooled":
        lam, mode = cfg.a5_lambda, "pooled"
    elif name == "A6_latent_32":
        arch["d_latent"] = base["d_latent"] // 2
    elif name in ("A7a_size_250k", "A7b_size_33k"):
        w = dict(cfg.a7a if name == "A7a_size_250k" else cfg.a7b)
        prober["hidden"] = w.pop("hidden")
        arch.update(w)
    elif name == "A8_actions_v1":
        arch.update(act_chunk=7, action_mode="v1")
    return {"name": name, "family": "v2", "change": DESCRIPTIONS[name], "arch": arch, "prober": prober,
            "lambda": lam, "sigreg_mode": mode, "wm_train": wm_train, "pi_train": pi_train}


def build_v2(spec: dict, device="cpu") -> Tuple[V.WorldModelV2, V.PIProberV2]:
    wm = V.WorldModelV2(**spec["arch"]).to(device)
    pi = V.PIProberV2(wm.M, wm.d_latent, hidden=spec["prober"]["hidden"],
                      state_aware=spec["prober"]["state_aware"]).to(device)
    return wm, pi


def build_v1(spec: dict, device="cpu") -> Tuple[M.JEPAWorldModel, M.PIProber]:
    wm = M.JEPAWorldModel(**spec["wm"]).to(device)
    pi = M.PIProber(wm.d_latent, state_aware=spec["prober"]["state_aware"]).to(device)
    return wm, pi


_PARAMS: Dict[str, dict] = {}


def param_count(cfg, name: str) -> dict:
    """Parameters per member (world model + prober) and the receptive field, without training."""
    spec = variant_spec(cfg, name)
    key = json.dumps(spec, sort_keys=True, default=str)
    if key not in _PARAMS:
        if spec["family"] == "v1":
            wm, pi = build_v1(spec)
            _PARAMS[key] = {"variant": name, "world_model": M.n_params(wm), "prober": M.n_params(pi),
                            "per_member": M.n_params(wm) + M.n_params(pi), "rf": wm.rf}
        else:
            wm, pi = build_v2(spec)
            r = V.param_report(wm, pi)
            _PARAMS[key] = {"variant": name, "world_model": r["per member"] - r["PI prober"], "prober": r["PI prober"],
                            "per_member": r["per member"], "rf": wm.rf}
    return dict(_PARAMS[key])


# ============================================================================
# Training windows with a shorter encoder history
# ============================================================================
class SlicedWindows:
    """
    models_v2.WindowsV2 seen by a model with a smaller receptive field rf: the
    last rf-1 of its P history steps. Front padding repeats the entry state and
    zero commands either way, so every window equals WindowsV2(..., rf, ...)'s
    (smoke_v4 checks this), with the same indexing and random draws.
    """

    def __init__(self, ds: V.WindowsV2, rf: int):
        self.ds, self.P, self.max_T, self.device = ds, int(rf) - 1, ds.max_T, ds.device
        self.cut = ds.P - self.P
        assert self.cut >= 0, "the shared windows must have at least this model's history"

    def __len__(self) -> int:
        return len(self.ds)

    def batch(self, idx: torch.Tensor, T: int) -> Dict[str, torch.Tensor]:
        b = self.ds.batch(idx, T)
        b["feats"], b["u"] = b["feats"][..., self.cut:, :], b["u"][..., self.cut:]
        return b

    def sample(self, shape, T: int, gen: Optional[torch.Generator] = None):
        idx = torch.randint(len(self), shape, generator=gen, device=self.device)
        return self.batch(idx, T)


def windows_for(data: dict, rf: int) -> Tuple[object, object]:
    """(train, val) windows for a v2-family model with receptive field rf."""
    if rf == data["ds_tr"].P + 1:
        return data["ds_tr"], data["ds_va"]
    return SlicedWindows(data["ds_tr"], rf), SlicedWindows(data["ds_va"], rf)


# ============================================================================
# Resumable training loop (models_v2.train_stage + RNG state)
# ============================================================================
def _rng_state() -> dict:
    st = {"cpu": torch.get_rng_state()}
    if torch.cuda.is_available():
        st["cuda"] = torch.cuda.get_rng_state_all()
    return st


def _set_rng_state(st: dict) -> None:
    torch.set_rng_state(st["cpu"])
    if "cuda" in st and torch.cuda.is_available():
        torch.cuda.set_rng_state_all(st["cuda"])


def train_stage_v4(module: torch.nn.Module, loss_fn: Callable[[int, int], tuple], cfg: V.TrainCfgV2,
                   device: torch.device, L, state_path=None, val_fn: Optional[Callable[[], Dict[str, float]]] = None,
                   label: str = "") -> Dict[str, list]:
    """
    models_v2.train_stage (AdamW, warmup + cosine, AMP, clip, curriculum; same
    arithmetic), whose train state (module, optimizer, scaler, step, history,
    elapsed time and RNG state) is written every cfg.save_every_s seconds and
    at the end, atomically through L's guarded writer. A rerun resumes there.
    """
    opt = torch.optim.AdamW(V._param_groups(module, cfg.weight_decay), lr=cfg.lr, betas=(0.9, 0.98))
    dt_amp = V.amp_dtype(device) if cfg.amp else None
    scaler = torch.amp.GradScaler(device.type, enabled=dt_amp == torch.float16)
    hist: Dict[str, list] = {"step": [], "T": [], "val": [], "sec": []}
    start, elapsed = 0, 0.0
    if state_path and os.path.exists(state_path):
        st = torch.load(state_path, map_location="cpu", weights_only=False)
        module.load_state_dict(st["module"])
        opt.load_state_dict(st["opt"])
        scaler.load_state_dict(st["scaler"])
        start, hist, elapsed = st["step"], st["hist"], st.get("elapsed", 0.0)
        _set_rng_state(st["rng"])
        print(f"  [{label}] resumed at step {start}/{cfg.steps} from {state_path}", flush=True)

    def save(step):
        if state_path:
            L.torch_save({"module": module.state_dict(), "opt": opt.state_dict(), "scaler": scaler.state_dict(),
                          "step": step, "hist": hist, "elapsed": elapsed, "rng": _rng_state()}, state_path)

    params = [p for g in opt.param_groups for p in g["params"]]
    acc: Dict[str, float] = {}
    n_acc, t_log, t_save = 0, time.time(), time.time()
    module.train()
    for step in range(start, cfg.steps):
        for g in opt.param_groups:
            g["lr"] = V._lr(step, cfg)
        T = V.horizon_at(step, cfg)
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
        if last or (state_path and time.time() - t_save > cfg.save_every_s):
            save(step + 1)
            t_save = time.time()
    module.eval()
    return hist


def train_world_model_v4(wm: V.WorldModelV2, ds_tr, ds_va, cfg: V.TrainCfgV2, lam: float, sigreg_mode: str,
                         device, L, state_path=None, label: str = "wm") -> Dict[str, list]:
    """train_v2.train_world_model (same validation batches and seeding) with the SIGReg mode, resumable via L."""
    sig = V.SIGRegV2(n_proj=32).to(device)
    lam_t = torch.as_tensor(lam, dtype=torch.float32, device=device)
    g = torch.Generator(device=ds_va.device).manual_seed(1)
    vbs = {T: ds_va.sample((wm.M, 256), T, g) for _, T in cfg.curriculum}
    state = {"T": cfg.curriculum[0][1]}

    def loss_fn(step, T):
        state["T"] = T
        return V.wm_loss_v2(wm, sig, ds_tr.sample((wm.M, cfg.batch), T), lam_t, sigreg_mode=sigreg_mode)

    def val_fn():
        _, logs = V.wm_loss_v2(wm, sig, vbs[state["T"]], 0.0, sigreg_mode=sigreg_mode)
        return {"pred": logs["pred"]}

    torch.manual_seed(cfg.seed)
    return train_stage_v4(wm, loss_fn, cfg, torch.device(device), L, state_path, val_fn, label)


def train_prober_v4(wm: V.WorldModelV2, pi: V.PIProberV2, ds_tr, ds_va, cfg: V.TrainCfgV2, device, L,
                    state_path=None, label: str = "prober") -> Dict[str, list]:
    """train_v2.train_prober (frozen world model, same seeding), resumable via L."""
    wm.eval()
    for p in wm.parameters():
        p.requires_grad_(False)
    g = torch.Generator(device=ds_va.device).manual_seed(2)
    vbs = {T: ds_va.sample((pi.M, 256), T, g) for _, T in cfg.curriculum}
    state = {"T": cfg.curriculum[0][1]}

    def loss_fn(step, T):
        state["T"] = T
        return V.prober_loss_v2(wm, pi, ds_tr.sample((pi.M, cfg.batch), T), D.MACRO_DT_S, 2, 1.0)

    def val_fn():
        return {"state": V.prober_loss_v2(wm, pi, vbs[state["T"]], D.MACRO_DT_S, 2, 0.0)[1]["state"]}

    torch.manual_seed(cfg.seed + 1)
    try:
        return train_stage_v4(pi, loss_fn, cfg, torch.device(device), L, state_path, val_fn, label)
    finally:
        for p in wm.parameters():
            p.requires_grad_(True)


# ============================================================================
# One variant, one seed
# ============================================================================
def _write_logs(L, d, name: str, seed: int, h1: dict, h2: dict, health) -> None:
    import pandas as pd
    log = pd.concat([S._hist_frame("world model", h1), S._hist_frame("prober", h2)], ignore_index=True)
    L.write_csv(log, d / "train_log.csv")
    L.write_csv(pd.DataFrame(health), d / "latent_health.csv")
    L.guard(d / "F_curves.png")
    S._curves_figure(log, f"{name} seed {seed}", d / "F_curves.png")


def train_variant_seed(cfg, data: dict, L, name: str, seed: int, device, grad_ckpt: bool = False,
                       keep_train_state: bool = False) -> dict:
    """Train (or resume, or skip if done) one variant for one seed; seed sets init, data order, projections."""
    d = L.model_dir(name, seed)
    ckpt = L.ckpt(name, seed)
    if ckpt.exists():
        return {"variant": name, "seed": seed, "ckpt": str(ckpt), "cached": True}
    L.guard(d).mkdir(parents=True, exist_ok=True)
    spec = variant_spec(cfg, name)
    info = {**spec, "seed": int(seed), "params": param_count(cfg, name), "status": "training"}
    L.write_json(L.variant_json(name, seed), info)
    t0 = time.time()
    device = torch.device(device)
    if spec["family"] == "v2":
        extra = _train_v2_family(cfg, spec, data, L, d, seed, device, grad_ckpt, ckpt)
    else:
        extra = _train_v1_family(cfg, spec, data, L, d, seed, device, ckpt)
    if not keep_train_state:
        for f in ("train_state_wm.pt", "train_state_pi.pt"):
            if (d / f).exists():
                L.guard(d / f).unlink()
    info.update({"status": "done", "train_s": round(time.time() - t0), **extra})
    L.write_json(L.variant_json(name, seed), info)
    return {"variant": name, "seed": seed, "ckpt": str(ckpt), "cached": False, "train_s": round(time.time() - t0)}


def _train_v2_family(cfg, spec, data, L, d, seed, device, grad_ckpt, ckpt) -> dict:
    from jepa import train_v2 as T2
    torch.manual_seed(seed)
    wm, pi = build_v2(spec, device)
    wm.grad_ckpt = bool(grad_ckpt)
    ds_tr, ds_va = windows_for(data, wm.rf)
    wm_cfg = V.TrainCfgV2(**{**spec["wm_train"], "seed": seed, "save_every_s": cfg.save_every_s})
    pi_cfg = V.TrainCfgV2(**{**spec["pi_train"], "seed": seed, "save_every_s": cfg.save_every_s})
    t0 = time.time()
    h1 = train_world_model_v4(wm, ds_tr, ds_va, wm_cfg, spec["lambda"], spec["sigreg_mode"], device, L,
                              str(L.guard(d / "train_state_wm.pt")), label=f"{spec['name']} seed {seed} wm")
    h2 = train_prober_v4(wm, pi, ds_tr, ds_va, pi_cfg, device, L, str(L.guard(d / "train_state_pi.pt")),
                         label=f"{spec['name']} seed {seed} prober")
    health = T2.latent_health_v2(wm, ds_va)
    _write_logs(L, d, spec["name"], seed, h1, h2, health)
    extra = {"variant": spec["name"], "seed": seed, "spec": json.loads(json.dumps(spec, default=str)),
             "lambda": spec["lambda"], "wm_cfg": asdict(wm_cfg), "pi_cfg": asdict(pi_cfg),
             "train_s": round(time.time() - t0), "n_train_episodes": len(data["train"]), "health": health}
    V.save_checkpoint_v2(str(L.guard(ckpt)), wm, pi, data["stats"], extra=extra)
    return {"min_eff_rank": min(h["eff_rank"] for h in health),
            "final_val_pred": (h1["val"][-1] or {}).get("pred") if h1["val"] else None,
            "final_val_state": (h2["val"][-1] or {}).get("state") if h2["val"] else None}


@torch.no_grad()
def _latent_health_v1(wm: M.JEPAWorldModel, ds, n: int = 4096, seed: int = 0) -> list:
    """train_v2.latent_health_v2's collapse check for a v1-family model (one member)."""
    wm.eval()
    b = ds.sample(n, torch.Generator(device=ds.device).manual_seed(seed))
    s = wm.encode_states(b["feats"])[:, 0].float()
    x = s - s.mean(0)
    sv = torch.linalg.svdvals(x)
    p = sv / sv.sum()
    corr = torch.corrcoef(x.T)
    off = corr[~torch.eye(corr.shape[0], dtype=torch.bool, device=corr.device)]
    return [{"member": 0, "eff_rank": float(torch.exp(-(p * torch.log(p + 1e-12)).sum())),
             "std_min": float(x.std(0).min()), "std_med": float(x.std(0).median()), "std_max": float(x.std(0).max()),
             "mean_abs_corr": float(off.abs().mean())}]


def _train_v1_family(cfg, spec, data, L, d, seed, device, ckpt) -> dict:
    c3 = cfg.v3_training()
    torch.manual_seed(seed)
    wm, pi = build_v1(spec, device)
    sig = M.SIGReg().to(device)
    stats = V.compute_stats_v2(data["train"]) if spec["feat_set"] == "v2" else D.compute_stats(data["train"])
    ds = S.V1Windows(data["train"], wm.rf, spec["T"], stats, device)
    dsp = S.V1Windows(data["train"], wm.rf, spec["T_pi"], stats, device)
    g = torch.Generator(device=device).manual_seed(1)
    vb = S.V1Windows(data["val"], wm.rf, spec["T"], stats, device).sample(1024, g)
    vbp = S.V1Windows(data["val"], wm.rf, spec["T_pi"], stats, device).sample(1024, g)
    print(f"  [{spec['name']} seed {seed}] {len(ds):,} windows at T={spec['T']}, rf {wm.rf}, {M.n_params(wm):,} + "
          f"{M.n_params(pi):,} parameters, device {device}", flush=True)
    t0 = time.time()
    h1 = S.fit_resumable(wm, lambda: M.jepa_loss(wm, sig, ds.sample(c3.v1_batch), spec["lambda"]), c3.v1_steps,
                         L.guard(d / "train_state_wm.pt"), lambda: M.jepa_loss(wm, sig, vb, 0.0)[1]["pred"], c3,
                         label=f"{spec['name']} seed {seed} wm")
    wm.eval()
    for p in wm.parameters():
        p.requires_grad_(False)
    h2 = S.fit_resumable(pi, lambda: M.prober_loss(wm, pi, dsp.sample(c3.v1_pi_batch), D.MACRO_DT_S, 2, w_acc=1.0),
                         c3.v1_pi_steps, L.guard(d / "train_state_pi.pt"),
                         lambda: M.prober_loss(wm, pi, vbp, D.MACRO_DT_S, 2, w_acc=0.0)[1]["state"], c3,
                         label=f"{spec['name']} seed {seed} prober")
    for p in wm.parameters():
        p.requires_grad_(True)
    health = _latent_health_v1(wm, S.V1Windows(data["val"], wm.rf, spec["T"], stats, device))
    _write_logs(L, d, spec["name"], seed, h1, h2, health)
    extra = {"variant": spec["name"], "seed": seed, "feat_set": spec["feat_set"],
             "spec": json.loads(json.dumps(spec, default=str)), "lambda": spec["lambda"],
             "steps": [c3.v1_steps, c3.v1_pi_steps], "batch": [c3.v1_batch, c3.v1_pi_batch],
             "extended": {"wm": h1["extended"], "prober": h2["extended"]}, "train_s": round(time.time() - t0),
             "n_train_episodes": len(data["train"]), "health": health}
    tmp = L.guard(str(ckpt) + ".tmp")
    M.save_checkpoint(str(tmp), wm.cpu(), pi.cpu(), None, None, stats, extra=extra)
    os.replace(tmp, ckpt)
    return {"min_eff_rank": health[0]["eff_rank"], "extended": extra["extended"],
            "final_val_pred": h1["val"][-1] if h1.get("val") else None,
            "final_val_state": h2["val"][-1] if h2.get("val") else None}


# ============================================================================
# Data
# ============================================================================
def load_data_v4(cfg, data_path: str, device) -> dict:
    """The 6,000 flights, the v3 split (seed 0, 80/10/10), v2 stats and windows long enough for every variant."""
    from jepa import train_v2 as T2
    eps = T2.load_episodes(data_path)
    if cfg.n_episodes:
        eps = eps[:cfg.n_episodes]
    rf = max(V.WorldModelV2(**variant_spec(cfg, n)["arch"]).rf for n in cfg.variants if family(n) == "v2")
    max_T = max(T for n in cfg.variants if family(n) == "v2" for k in ("wm_train", "pi_train")
                for _, T in variant_spec(cfg, n)[k]["curriculum"])
    data = T2.prepare(eps, max_T, device, rf, split_seed=cfg.split_seed)
    data["n_episodes"] = len(eps)
    return data
