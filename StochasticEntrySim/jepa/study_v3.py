"""
jepa/study_v3.py

The v3 study's plumbing (jepa/v3_plan.md sections 2-4, 9):

    StudyConfig        every setting of the study, frozen to 00_study/study_config.json
                       together with the git commit and package versions
    Layout             the Google Drive folder layout (00_study ... 06_report)
    atomic writes      temp file + rename, so a disconnect never leaves half a file
    manifest/README    an index of every artifact (path, stage, size, sha1)
    train_v2_seed      one v2 ensemble per seed (lambda fixed at 0.08), resumable
    train_v1_seed      v1's architecture and recipe with v2's per-member budget,
                       on the GPU, resumable, with the plan's convergence-extension rule
    build_jobs         every flight of the study, keyed (controller, load, seed,
                       mppi_rep, case, atmosphere), in seed-major order
    run_classical / run_learned
                       the resumable job queue: a flight whose result file
                       exists is skipped; results and 1 Hz telemetry are written
                       to Drive as soon as each flight lands
    collect_flights    rebuild 04_flights/flights.csv from the result files
"""

from __future__ import annotations

import hashlib
import json
import math
import multiprocessing as mp
import os
import platform
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

from jepa import data as D
from jepa import fly_v3 as F
from jepa import testset_v3 as TS

LEARNED = ("v1", "v2", "v2single")


# ============================================================================
# Configuration
# ============================================================================
@dataclass
class StudyConfig:
    """Every number of the study (v3_plan.md). Frozen to study_config.json on first use."""
    mode: str = "full"                                  # "full" | "smoke" (every stage, tiny, numbers meaningless)
    seeds: Tuple[int, ...] = (0, 1, 2, 3, 4)
    split_seed: int = 0
    n_episodes: Optional[int] = None                    # None = all 6,000 flights
    testset_csv: str = "configs/v3_testset.csv"
    cases: Optional[Tuple[int, ...]] = None             # None = all 64
    atmospheres: Tuple[str, ...] = tuple(TS.ATMOSPHERES)
    # ---- v2 (as in the first full run, lambda fixed) ----
    v2_members: int = 5
    v2_arch: dict = field(default_factory=dict)         # models_v2 defaults (2.45M parameters per member)
    v2_lambda: float = 0.08
    v2_wm: dict = field(default_factory=lambda: dict(
        steps=24_000, batch=256, warmup=1_000, lr=6e-4, curriculum=((0.0, 24), (0.2, 48), (0.5, 96)), log_every=500))
    v2_pi: dict = field(default_factory=lambda: dict(
        steps=6_000, batch=512, warmup=300, lr=1e-3, weight_decay=1e-4,
        curriculum=((0.0, 16), (0.3, 32), (0.6, 48)), log_every=250))
    # ---- v1 with a fair budget (v1's recipe, one v2 member's steps and batch) ----
    v1_lambda: float = 0.02
    v1_T: int = 24
    v1_T_pi: int = 16
    v1_steps: int = 24_000
    v1_batch: int = 256
    v1_pi_steps: int = 6_000
    v1_pi_batch: int = 512
    v1_lr: float = 3e-3
    v1_lr_min: float = 1e-4
    v1_warmup_frac: float = 0.05
    v1_clip: float = 0.5
    v1_weight_decay: float = 1e-5
    v1_log_every: int = 500
    v1_extend_tol: float = 0.02                         # extend once if val improved > 2% over the last 20%
    v1_extend_frac: float = 0.5                         # ... by 50% more steps
    v1_extend_lr_frac: float = 0.3                      # warm restart at 30% of the peak lr, cosine to lr_min
    save_every_s: float = 600.0
    # ---- flights ----
    g_primary: float = 4.5
    g_legacy: float = 8.0
    legacy_seeds: Tuple[int, ...] = (0,)
    noise_seed: int = 0
    noise_reps: Tuple[int, ...] = (1, 2)
    run_oracle: bool = False
    mppi: dict = field(default_factory=lambda: dict(n_samples=256, n_knots=6, noise_std_deg=30.0, n_iters=1,
                                                    n_iters_first=4, n_sub=2))
    horizon: Optional[int] = None                       # None = from the test set (plan 4.2)
    w_unc_km: float = 1.0
    max_rows: int = 4096                                # samples per batched rollout (latents: M x rows x K x D floats)
    v2_chunk: int = 2048                                # rows per predictor pass (activations ~ M x chunk x K x 768)
    # ---- analysis ----
    success_km: float = 5.0
    legacy_success_km: float = 25.0
    n_boot: int = 10_000
    ttg_s: Tuple[int, ...] = (400, 300, 200, 100, 50)

    @staticmethod
    def smoke() -> "StudyConfig":
        """Every stage on a CPU in minutes: tiny models, 2 seeds, 2 cases, 32 samples."""
        tiny = dict(d_latent=16, enc_dim=32, enc_layers=2, enc_heads=2, enc_band=5, pred_dim=48, pred_layers=2,
                    pred_heads=2, act_chunk=4, mlp_ratio=2)
        return StudyConfig(
            mode="smoke", seeds=(0, 1), n_episodes=160, cases=(0, 16), v2_members=2, v2_arch=tiny,
            v2_wm=dict(steps=120, batch=32, warmup=10, curriculum=((0.0, 8), (0.5, 16)), log_every=60, amp=True),
            v2_pi=dict(steps=60, batch=32, warmup=5, lr=2e-3, curriculum=((0.0, 8), (0.5, 16)), log_every=30),
            v1_steps=120, v1_batch=64, v1_pi_steps=60, v1_pi_batch=64, v1_log_every=30, save_every_s=30.0,
            mppi=dict(n_samples=32, n_knots=6, noise_std_deg=30.0, n_iters=1, n_iters_first=2, n_sub=2),
            noise_reps=(1,), n_boot=200, ttg_s=(200, 100))


def git_commit() -> str:
    """The checked-out commit (or 'unknown')."""
    try:
        return subprocess.run(["git", "-C", str(SIM_DIR), "rev-parse", "HEAD"], capture_output=True, text=True,
                              timeout=10).stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def environment_info() -> dict:
    """Package versions and hardware, recorded next to the config."""
    import pandas as pd
    info = {"python": platform.python_version(), "torch": torch.__version__, "numpy": np.__version__,
            "pandas": pd.__version__, "platform": platform.platform(), "cpus": os.cpu_count(),
            "cuda": torch.version.cuda, "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None}
    try:
        import pymsis
        info["pymsis"] = getattr(pymsis, "__version__", "?")
    except Exception:
        pass
    return info


# ============================================================================
# Drive layout and atomic writes
# ============================================================================
STAGES = {
    "00_study": "frozen settings (study_config.json), this index (README.md), and manifest.csv",
    "01_testset": "the 64 entry cases, their reachable footprints, coverage figures",
    "02_models": "v1 and v2 checkpoints per seed, training logs and curves",
    "03_open_loop": "open-loop footprint errors per model and seed (CSV)",
    "04_flights": "one JSON result + one 1 Hz telemetry .npz per closed-loop flight; flights.csv merges them",
    "05_analysis": "numbered tables (T##) and figures (F##) with captions.json",
    "06_report": "SpaceJEPA_v3_report.pdf, generated from 05_analysis",
}


class Layout:
    """Paths of the study folder (v3_plan.md section 9)."""

    def __init__(self, root):
        self.root = Path(root)
        self.study = self.root / "00_study"
        self.testset = self.root / "01_testset"
        self.models = self.root / "02_models"
        self.open_loop = self.root / "03_open_loop"
        self.flights = self.root / "04_flights"
        self.analysis = self.root / "05_analysis"
        self.tables = self.analysis / "tables"
        self.figures = self.analysis / "figures"
        self.report = self.root / "06_report"

    def make(self) -> "Layout":
        for p in (self.study, self.testset, self.models, self.open_loop, self.flights / "results",
                  self.flights / "telemetry", self.tables, self.figures, self.report):
            p.mkdir(parents=True, exist_ok=True)
        return self

    def model_dir(self, kind: str, seed: int) -> Path:
        return self.models / kind / f"seed_{int(seed)}"

    def ckpt(self, kind: str, seed: int) -> Path:
        return self.model_dir(kind, seed) / "ckpt.pt"

    def result_path(self, fid: str) -> Path:
        return self.flights / "results" / f"{fid}.json"

    def telemetry_path(self, fid: str) -> Path:
        return self.flights / "telemetry" / f"{fid}.npz"


def _tmp(path: Path) -> Path:
    return path.with_name(f".{path.name}.{os.getpid()}.tmp")


def atomic_bytes(path, data: bytes) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = _tmp(path)
    with open(tmp, "wb") as f:
        f.write(data)
    os.replace(tmp, path)


def write_json(path, obj) -> None:
    atomic_bytes(path, json.dumps(obj, indent=1, default=_json_default).encode())


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (np.ndarray,)):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def write_csv(df, path, **kw) -> None:
    atomic_bytes(path, df.to_csv(index=kw.pop("index", False), **kw).encode())


def save_fig(fig, path, dpi: int = 150) -> None:
    import io
    buf = io.BytesIO()
    fig.savefig(buf, format=Path(path).suffix[1:] or "png", dpi=dpi, bbox_inches="tight")
    atomic_bytes(path, buf.getvalue())


def write_flight(layout: Layout, rec: dict) -> None:
    """Telemetry first, then the result JSON (its presence marks the flight done)."""
    import io
    fid = rec["result"]["flight_id"]
    buf = io.BytesIO()
    np.savez_compressed(buf, **rec["telemetry"])
    atomic_bytes(layout.telemetry_path(fid), buf.getvalue())
    write_json(layout.result_path(fid), rec["result"])


def flight_done(layout: Layout, job: dict) -> bool:
    return layout.result_path(F.flight_id(job)).exists()


def sha1(path, block: int = 1 << 20) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(block), b""):
            h.update(b)
    return h.hexdigest()


def freeze_config(cfg: StudyConfig, layout: Layout, extra: Optional[dict] = None) -> dict:
    """
    Write 00_study/study_config.json on first use. Later sessions must use the
    same settings: a difference raises (start a new study folder instead).
    """
    path = layout.study / "study_config.json"
    want = json.loads(json.dumps(asdict(cfg), default=_json_default))
    if path.exists():
        have = json.loads(path.read_text())
        diff = {k: (have["config"].get(k), v) for k, v in want.items() if have["config"].get(k) != v}
        if diff:
            raise ValueError(f"settings differ from the frozen {path}: {diff}. Use a new study folder.")
        return have
    doc = {"config": want, "git_commit": git_commit(), "created": time.strftime("%Y-%m-%d %H:%M:%S"),
           "environment": environment_info(), **(extra or {})}
    write_json(path, doc)
    return doc


def update_manifest(layout: Layout) -> "object":
    """
    Rebuild 00_study/manifest.csv (path, stage, created, size, sha1) and
    README.md. sha1s of unchanged files (same size and mtime) are reused.
    """
    import pandas as pd
    mpath = layout.study / "manifest.csv"
    old = {}
    if mpath.exists():
        o = pd.read_csv(mpath)
        old = {r.path: r for r in o.itertuples()}
    rows = []
    for p in sorted(layout.root.rglob("*")):
        if not p.is_file() or p.name.startswith(".") or p.name == "manifest.csv":
            continue
        rel = str(p.relative_to(layout.root))
        st = p.stat()
        created = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(st.st_mtime))
        o = old.get(rel)
        digest = o.sha1 if o is not None and o.size == st.st_size and o.created == created else sha1(p)
        rows.append({"path": rel, "stage": rel.split("/")[0], "created": created, "size": st.st_size, "sha1": digest})
    df = pd.DataFrame(rows, columns=["path", "stage", "created", "size", "sha1"])
    write_csv(df, mpath)
    _write_readme(layout, df)
    return df


def _write_readme(layout: Layout, man) -> None:
    lines = ["# SpaceJEPA v3 study folder", "",
             "Auto-written by `jepa/study_v3.py` (do not edit). The plan is `jepa/v3_plan.md` in the repo.", ""]
    cfg = layout.study / "study_config.json"
    if cfg.exists():
        c = json.loads(cfg.read_text())
        lines += [f"* git commit: `{c.get('git_commit')}`", f"* created: {c.get('created')}",
                  f"* mode: {c['config'].get('mode')}", ""]
    lines += ["| folder | what | files | MB |", "|---|---|---|---|"]
    for stage, what in STAGES.items():
        sub = man[man.stage == stage]
        lines.append(f"| `{stage}/` | {what} | {len(sub)} | {sub['size'].sum() / 1e6:.1f} |")
    caps = layout.analysis / "captions.json"
    if caps.exists():
        lines += ["", "## Figures and tables", ""]
        for k, v in sorted(json.loads(caps.read_text()).items()):
            lines.append(f"* `{k}`: {v.get('caption', '')} (data: `{v.get('data', '')}`)")
    atomic_bytes(layout.study / "README.md", ("\n".join(lines) + "\n").encode())


# ============================================================================
# Data
# ============================================================================
def load_data(cfg: StudyConfig, data_path: str, device) -> dict:
    """The 6,000 flights, macro-stepped, one fixed split, v2 windows on the device."""
    from jepa import models_v2 as V
    from jepa import train_v2 as T2
    eps = T2.load_episodes(data_path)
    if cfg.n_episodes:
        eps = eps[:cfg.n_episodes]
    rf = V.WorldModelV2(n_members=1, **cfg.v2_arch).rf
    max_T = max(T for _, T in tuple(cfg.v2_wm["curriculum"]) + tuple(cfg.v2_pi["curriculum"]))
    data = T2.prepare(eps, max_T, device, rf, split_seed=cfg.split_seed)
    data["n_episodes"] = len(eps)
    return data


class V1Windows(D.MacroDataset):
    """data.MacroDataset (v1's windows, unchanged recipe) resident on one device."""

    def __init__(self, macro_eps, rf: int, horizon: int, stats, device):
        super().__init__(macro_eps, rf, horizon, stats)
        self.device = torch.device(device)
        for k in ("Y", "F", "U", "A", "y_start", "u_start"):
            setattr(self, k, getattr(self, k).to(self.device))

    def batch(self, idx: torch.Tensor) -> Dict[str, torch.Tensor]:
        P, T, dev = self.P, self.T, self.device
        idx = idx.to(dev)
        iy = self.y_start[idx, None] + torch.arange(P + T + 1, device=dev)
        iu = self.u_start[idx, None] + torch.arange(P + T, device=dev)
        return {"feats": self.F[iy], "u": self.U[iu], "y": self.Y[iy[:, P:]], "acc": self.A[iu[:, P:]]}

    def sample(self, batch_size: int, gen: Optional[torch.Generator] = None) -> Dict[str, torch.Tensor]:
        dev = gen.device if gen is not None else self.device
        return self.batch(torch.randint(len(self), (batch_size,), generator=gen, device=dev))


# ============================================================================
# Training
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


def fit_resumable(module: torch.nn.Module, loss_fn, n_steps: int, state_path, val_fn, cfg: StudyConfig,
                  label: str = "") -> dict:
    """
    v1's models.fit (AdamW, warmup + cosine, clip 0.5, wd 1e-5) made resumable,
    plus the plan's convergence rule: if the validation loss improved by more
    than v1_extend_tol over the last 20% of the steps, train once more for
    v1_extend_frac x the steps (warm restart at v1_extend_lr_frac x peak lr,
    cosine to lr_min). Returns the history {step, lr, val, sec, logs..., extended}.
    """
    from jepa.models import lr_at
    params = [p for p in module.parameters() if p.requires_grad]
    opt = torch.optim.AdamW(params, lr=cfg.v1_lr, weight_decay=cfg.v1_weight_decay)
    st = {"step": 0, "n_ext": 0, "extended": None, "elapsed": 0.0,
          "hist": {"step": [], "lr": [], "val": [], "sec": []}}
    state_path = Path(state_path)
    if state_path.exists():
        ck = torch.load(state_path, map_location="cpu", weights_only=False)
        module.load_state_dict(ck["module"])
        opt.load_state_dict(ck["opt"])
        st = ck["st"]
        _set_rng_state(ck["rng"])
        print(f"  [{label}] resumed at step {st['step']} from {state_path}")

    def save():
        import io
        buf = io.BytesIO()
        torch.save({"module": module.state_dict(), "opt": opt.state_dict(), "st": st, "rng": _rng_state()}, buf)
        atomic_bytes(state_path, buf.getvalue())

    def lr(step):
        if step < n_steps:
            return lr_at(step, n_steps, cfg.v1_lr, cfg.v1_lr_min, max(1, int(cfg.v1_warmup_frac * n_steps)))
        k, n = step - n_steps, st["n_ext"]
        return lr_at(k, n, cfg.v1_extend_lr_frac * cfg.v1_lr, cfg.v1_lr_min, max(1, int(cfg.v1_warmup_frac * n)))

    def decide_extension():
        h = st["hist"]
        steps, val = np.asarray(h["step"]), np.asarray(h["val"], dtype=float)
        i80 = max(0, int(np.searchsorted(steps, 0.8 * n_steps, side="right")) - 1)   # last log at <= 80%
        v80, v_end = float(val[i80]), float(val[-1])
        gain = (v80 - v_end) / max(abs(v80), 1e-12)
        st["extended"] = {"val_at_80pct": v80, "val_at_end": v_end, "gain": gain,
                          "extend": bool(gain > cfg.v1_extend_tol)}
        st["n_ext"] = int(round(cfg.v1_extend_frac * n_steps)) if gain > cfg.v1_extend_tol else 0
        print(f"  [{label}] val gain over the last 20%: {100 * gain:.2f}% -> "
              f"{'extend by ' + str(st['n_ext']) + ' steps' if st['n_ext'] else 'converged, no extension'}")

    dev = next(module.parameters()).device
    acc: Dict[str, float] = {}
    n_acc, t_log, t_save = 0, time.time(), time.time()
    module.train()
    while True:
        total = n_steps + st["n_ext"]
        if st["step"] >= n_steps and st["extended"] is None:
            decide_extension()
            save()
            total = n_steps + st["n_ext"]
        if st["step"] >= total:
            break
        step = st["step"]
        for g in opt.param_groups:
            g["lr"] = lr(step)
        loss, logs = loss_fn()
        opt.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(params, cfg.v1_clip)
        opt.step()
        for k, v in logs.items():
            acc[k] = acc.get(k, 0.0) + float(v)
        n_acc += 1
        st["step"] = step + 1
        end_of_phase = st["step"] in (n_steps, total)
        if st["step"] % cfg.v1_log_every == 0 or end_of_phase:
            if dev.type == "cuda":
                torch.cuda.synchronize()
            now = time.time()
            st["elapsed"] += now - t_log
            per = (now - t_log) / max(n_acc, 1)
            t_log = now
            module.eval()
            with torch.no_grad():
                val = float(val_fn())
            module.train()
            h = st["hist"]
            h["step"].append(st["step"])
            h["lr"].append(lr(step))
            h["val"].append(val)
            h["sec"].append(round(st["elapsed"], 1))
            for k, v in acc.items():
                h.setdefault(k, []).append(v / n_acc)
            msg = " ".join(f"{k}={v / n_acc:.4g}" for k, v in acc.items())
            eta = per * (total - st["step"]) / 60
            print(f"  [{label}] step {st['step']:>6}/{total} {msg} val={val:.4g} "
                  f"({per * 1000:.0f} ms/step, {st['elapsed'] / 60:.1f} min, ~{eta:.0f} min left)", flush=True)
            acc, n_acc = {}, 0
        if end_of_phase or time.time() - t_save > cfg.save_every_s:
            save()
            t_save = time.time()
    module.eval()
    return {**st["hist"], "extended": st["extended"], "elapsed_s": st["elapsed"]}


def _hist_frame(stage: str, hist: dict):
    """A training history (either loop's) as a long-format DataFrame."""
    import pandas as pd
    rows = []
    for i, s in enumerate(hist.get("step", [])):
        r = {"stage": stage, "step": s}
        for k, v in hist.items():
            if k in ("step", "extended", "elapsed_s") or not isinstance(v, list) or i >= len(v):
                continue
            if isinstance(v[i], dict):
                r.update({f"val_{kk}": vv for kk, vv in v[i].items()})
            elif k == "val":
                r["val"] = v[i]
            else:
                r[k] = v[i]
        rows.append(r)
    return pd.DataFrame(rows)


def _curves_figure(log, title: str, path) -> None:
    import matplotlib.pyplot as plt
    stages = list(dict.fromkeys(log.stage))
    fig, axes = plt.subplots(1, len(stages), figsize=(5.2 * len(stages), 3.2), squeeze=False)
    for ax, stg in zip(axes[0], stages):
        d = log[log.stage == stg]
        for col in [c for c in d.columns if c not in ("stage", "step", "T", "sec", "lr") and d[c].notna().any()]:
            if d[col].dtype.kind in "fi" and (d[col] > 0).all():
                ax.plot(d.step, d[col], label=col)
        ax.set(yscale="log", xlabel="step", title=f"{title}: {stg}")
        ax.legend(fontsize=7)
    save_fig(fig, path)
    plt.close(fig)


def train_v2_seed(cfg: StudyConfig, data: dict, layout: Layout, seed: int, device, grad_ckpt: bool = False,
                  keep_train_state: bool = False) -> dict:
    """One v2 ensemble for one seed (init, data order, SIGReg projections all from `seed`)."""
    from jepa import models_v2 as V
    from jepa import train_v2 as T2
    d = layout.model_dir("v2", seed)
    d.mkdir(parents=True, exist_ok=True)
    ckpt = layout.ckpt("v2", seed)
    if ckpt.exists():
        return {"seed": seed, "ckpt": str(ckpt), "cached": True}
    t0 = time.time()
    torch.manual_seed(seed)
    wm = V.WorldModelV2(n_members=cfg.v2_members, **cfg.v2_arch).to(device)
    wm.grad_ckpt = bool(grad_ckpt)
    pi = V.PIProberV2(wm.M, wm.d_latent).to(device)
    wm_cfg = V.TrainCfgV2(**{**cfg.v2_wm, "seed": seed, "save_every_s": cfg.save_every_s})
    pi_cfg = V.TrainCfgV2(**{**cfg.v2_pi, "seed": seed, "save_every_s": cfg.save_every_s})
    h1 = T2.train_world_model(wm, data, wm_cfg, cfg.v2_lambda, device, str(d / "train_state_wm.pt"),
                              label=f"v2 seed {seed} wm")
    h2 = T2.train_prober(wm, pi, data, pi_cfg, device, str(d / "train_state_pi.pt"), label=f"v2 seed {seed} prober")
    health = T2.latent_health_v2(wm, data["ds_va"])
    import pandas as pd
    log = pd.concat([_hist_frame("world model", h1), _hist_frame("prober", h2)], ignore_index=True)
    write_csv(log, d / "train_log.csv")
    write_csv(pd.DataFrame(health), d / "latent_health.csv")
    _curves_figure(log, f"v2 seed {seed}", d / "F_curves.png")
    V.save_checkpoint_v2(str(ckpt), wm, pi, data["stats"],
                         extra={"variant": "v2", "seed": seed, "lambda": cfg.v2_lambda, "wm_cfg": asdict(wm_cfg),
                                "pi_cfg": asdict(pi_cfg), "train_s": round(time.time() - t0),
                                "train_min_wm": h1["sec"][-1] / 60 if h1.get("sec") else None,
                                "train_min_pi": h2["sec"][-1] / 60 if h2.get("sec") else None,
                                "n_train_episodes": len(data["train"]), "health": health})
    if not keep_train_state:
        for f in ("train_state_wm.pt", "train_state_pi.pt"):
            (d / f).unlink(missing_ok=True)
    return {"seed": seed, "ckpt": str(ckpt), "cached": False, "train_s": round(time.time() - t0)}


def train_v1_seed(cfg: StudyConfig, data: dict, layout: Layout, seed: int, device,
                  keep_train_state: bool = False) -> dict:
    """
    v1's architecture and recipe (26k-parameter TCN + GRU, pooled SIGReg 0.02,
    T = 24, PI prober at T = 16) with one v2 member's budget, on the same
    split, on the GPU. Validation = latent MSE (stage 1) and state loss (stage 2)
    on fixed validation batches.
    """
    from jepa import models as M
    d = layout.model_dir("v1", seed)
    d.mkdir(parents=True, exist_ok=True)
    ckpt = layout.ckpt("v1", seed)
    if ckpt.exists():
        return {"seed": seed, "ckpt": str(ckpt), "cached": True}
    t0 = time.time()
    device = torch.device(device)
    torch.manual_seed(seed)
    wm = M.JEPAWorldModel().to(device)
    pi = M.PIProber(wm.d_latent).to(device)
    sig = M.SIGReg().to(device)
    stats = D.compute_stats(data["train"])
    ds = V1Windows(data["train"], wm.rf, cfg.v1_T, stats, device)
    dsp = V1Windows(data["train"], wm.rf, cfg.v1_T_pi, stats, device)
    g = torch.Generator(device=device).manual_seed(1)
    vb = V1Windows(data["val"], wm.rf, cfg.v1_T, stats, device).sample(1024, g)
    vbp = V1Windows(data["val"], wm.rf, cfg.v1_T_pi, stats, device).sample(1024, g)
    print(f"  [v1 seed {seed}] {len(ds):,} training windows at T={cfg.v1_T}, {M.n_params(wm):,} + "
          f"{M.n_params(pi):,} parameters, device {device}")
    h1 = fit_resumable(wm, lambda: M.jepa_loss(wm, sig, ds.sample(cfg.v1_batch), cfg.v1_lambda), cfg.v1_steps,
                       d / "train_state_wm.pt", lambda: M.jepa_loss(wm, sig, vb, 0.0)[1]["pred"], cfg,
                       label=f"v1 seed {seed} wm")
    wm.eval()
    for p in wm.parameters():
        p.requires_grad_(False)
    h2 = fit_resumable(pi, lambda: M.prober_loss(wm, pi, dsp.sample(cfg.v1_pi_batch), D.MACRO_DT_S, 2, w_acc=1.0),
                       cfg.v1_pi_steps, d / "train_state_pi.pt",
                       lambda: M.prober_loss(wm, pi, vbp, D.MACRO_DT_S, 2, w_acc=0.0)[1]["state"], cfg,
                       label=f"v1 seed {seed} prober")
    for p in wm.parameters():
        p.requires_grad_(True)
    import pandas as pd
    log = pd.concat([_hist_frame("world model", h1), _hist_frame("prober", h2)], ignore_index=True)
    write_csv(log, d / "train_log.csv")
    _curves_figure(log, f"v1 seed {seed}", d / "F_curves.png")
    tmp = str(ckpt) + ".tmp"
    M.save_checkpoint(tmp, wm.cpu(), pi.cpu(), None, None, stats,
                      extra={"variant": "v1-fair", "seed": seed, "lambda": cfg.v1_lambda,
                             "steps": [cfg.v1_steps, cfg.v1_pi_steps], "batch": [cfg.v1_batch, cfg.v1_pi_batch],
                             "extended": {"wm": h1["extended"], "prober": h2["extended"]},
                             "train_s": round(time.time() - t0), "n_train_episodes": len(data["train"])})
    os.replace(tmp, ckpt)
    if not keep_train_state:
        for f in ("train_state_wm.pt", "train_state_pi.pt"):
            (d / f).unlink(missing_ok=True)
    return {"seed": seed, "ckpt": str(ckpt), "cached": False, "train_s": round(time.time() - t0),
            "extended": {"wm": h1["extended"], "prober": h2["extended"]}}


# ============================================================================
# Flight jobs
# ============================================================================
def load_testset(cfg: StudyConfig):
    """The frozen cases (optionally a subset) as a DataFrame."""
    import pandas as pd
    cases = pd.read_csv(SIM_DIR / cfg.testset_csv if not os.path.isabs(cfg.testset_csv) else cfg.testset_csv)
    if cfg.cases is not None:
        cases = cases[cases.case.isin(cfg.cases)].reset_index(drop=True)
    return cases


def horizon_steps(cfg: StudyConfig) -> int:
    """The MPPI horizon of every learned controller (plan 4.2), from the full test set."""
    if cfg.horizon:
        return int(cfg.horizon)
    import pandas as pd
    return TS.planning_horizon_steps(pd.read_csv(SIM_DIR / cfg.testset_csv), D.MACRO_DT_S)


def model_ckpt(layout: Layout, controller: str, seed: int) -> Tuple[str, str, Optional[Tuple[int, ...]]]:
    """(adapter kind, checkpoint path, members) of a learned controller."""
    if controller == "v1":
        return "v1", str(layout.ckpt("v1", seed)), None
    if controller == "v2single":
        return "v2", str(layout.ckpt("v2", seed)), (0,)
    return "v2", str(layout.ckpt("v2", seed)), None


def build_jobs(cfg: StudyConfig, layout: Layout) -> List[dict]:
    """Every flight of the study, classical first, then seed-major (v3_plan.md section 3, 4.3)."""
    cases = load_testset(cfg)
    H = horizon_steps(cfg)
    conds = [(c, a) for c in cases.to_dict("records") for a in cfg.atmospheres]

    def base(c, atm):
        return {"case": int(c["case"]), "set": c["set"], "atmosphere": atm, "ic": TS.ic_of(c),
                "domain": dict(TS.ATMOSPHERES[atm])}

    jobs = [{**base(c, a), "controller": "classical", "load": "own", "g_limit": None, "seed": None,
             "mppi_rep": None} for c, a in conds]
    if cfg.run_oracle:
        jobs += [{**base(c, a), "controller": "classical_oracle", "load": "own", "g_limit": None, "seed": None,
                  "mppi_rep": None, "classical_knows_truth": True} for c, a in conds if a != "nominal"]

    def learned(ctrl, load, g, seed, rep):
        out = []
        for c, a in conds:
            s = F.crn_seed(seed, int(c["case"]), a, rep)
            out.append({**base(c, a), "controller": ctrl, "load": load, "g_limit": g, "seed": seed,
                        "mppi_rep": rep, "mppi_seed": s, "horizon": H,
                        "mppi": {"seed": s, "g_limit": g, "horizon": H,
                                 "w_unc_km": cfg.w_unc_km if ctrl == "v2" else 0.0}})
        return out

    for seed in cfg.seeds:
        for ctrl in LEARNED:
            jobs += learned(ctrl, f"{cfg.g_primary:g}g", cfg.g_primary, seed, 0)
        if seed in cfg.legacy_seeds:
            for ctrl in ("v1", "v2"):
                jobs += learned(ctrl, f"{cfg.g_legacy:g}g", cfg.g_legacy, seed, 0)
        if seed == cfg.noise_seed:
            for rep in cfg.noise_reps:
                jobs += learned("v2", f"{cfg.g_primary:g}g", cfg.g_primary, seed, rep)
    return jobs


def groups(jobs: Sequence[dict]) -> List[Tuple[tuple, List[dict]]]:
    """Learned jobs grouped by model run (controller, load, seed, rep), in queue order."""
    out: Dict[tuple, List[dict]] = {}
    for j in jobs:
        if j["controller"] in LEARNED:
            out.setdefault((j["controller"], j["load"], j["seed"], j["mppi_rep"]), []).append(j)
    return list(out.items())


# ============================================================================
# The job queue
# ============================================================================
def _classical_worker(arg) -> dict:
    job, root = arg
    rec = F.fly_classical(job)
    write_flight(Layout(root), rec)
    return rec["result"]


def run_classical(jobs: Sequence[dict], layout: Layout, workers: int = 1, log=print) -> List[dict]:
    """Fly every pending classical job in a CPU process pool; each result is written as it lands."""
    todo = [j for j in jobs if j["controller"].startswith("classical") and not flight_done(layout, j)]
    n_all = sum(j["controller"].startswith("classical") for j in jobs)
    log(f"[classical] {n_all - len(todo)}/{n_all} done, {len(todo)} to fly on {workers} worker(s)")
    if not todo:
        return []
    # Longest flights first (mismatch_high glides farthest), so the pool drains evenly.
    order = {"mismatch_high": 0, "nominal": 1, "mismatch": 2, "out_of_range": 3}
    todo.sort(key=lambda j: order.get(j["atmosphere"], 4))
    out, t0 = [], time.time()
    args = [(j, str(layout.root)) for j in todo]
    if workers <= 1:
        it = map(_classical_worker, args)
        pool = None
    else:
        pool = mp.get_context("spawn").Pool(int(workers))
        it = pool.imap_unordered(_classical_worker, args)
    try:
        for i, r in enumerate(it, 1):
            out.append(r)
            if i % max(1, len(todo) // 20) == 0 or i == len(todo):
                el = time.time() - t0
                log(f"  [classical] {i}/{len(todo)}  last: {r['flight_id']} {r['outcome']} {r['miss_km']:.2f} km  "
                    f"({el / 60:.1f} min, ~{el / i * (len(todo) - i) / 60:.0f} min left)")
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    return out


def _learned_worker(arg) -> int:
    """Pool worker: fly one shard of a model run in lockstep; results are written by on_done."""
    shard, kind, ckpt, members, device, mppi, max_rows, chunk, root, threads, verbose = arg
    torch.set_num_threads(int(threads))
    adapter = F.load_adapter(kind, ckpt, device, members, chunk)
    layout = Layout(root)
    F.run_lockstep_v3(shard, adapter, mppi, n_sub=int(mppi.get("n_sub", 2)), max_rows=max_rows,
                      on_done=lambda rec: write_flight(layout, rec), verbose=verbose)
    return len(shard)


def run_learned(cfg: StudyConfig, jobs: Sequence[dict], layout: Layout, device: str = "cuda", workers: int = 1,
                only: Optional[Sequence[tuple]] = None, log=print) -> None:
    """
    Fly every pending learned job, one model run (controller, load, seed, rep)
    at a time in queue order, its flights in lockstep, split over `workers`
    processes that share the GPU (each steps its own envs on one CPU core).
    """
    gs = [(k, g) for k, g in groups(jobs) if only is None or k in only]
    ready = [(k, g) for k, g in gs if os.path.exists(model_ckpt(layout, k[0], k[2])[1])]
    n_all = sum(len(g) for _, g in ready)
    n_done0 = sum(flight_done(layout, j) for _, g in ready for j in g)
    log(f"[learned] {n_done0}/{n_all} flights done in {len(ready)} model runs with a checkpoint "
        f"({len(gs) - len(ready)} runs wait for training); {workers} worker(s) on {device}")
    pool = mp.get_context("spawn").Pool(int(workers)) if workers > 1 else None
    t_all, n_flown = time.time(), 0
    try:
        for key, g in gs:
            todo = [j for j in g if not flight_done(layout, j)]
            if not todo:
                continue
            ctrl, load, seed, rep = key
            kind, ckpt, members = model_ckpt(layout, ctrl, seed)
            if not os.path.exists(ckpt):
                log(f"  [learned] skip {key}: no checkpoint {ckpt} (train seed {seed} first)")
                continue
            t0 = time.time()
            shards = [todo[i::max(1, workers)] for i in range(max(1, workers))]
            args = [(s, kind, ckpt, members, device, dict(cfg.mppi), cfg.max_rows, cfg.v2_chunk, str(layout.root), 1,
                     pool is None) for s in shards if s]
            if pool is None:
                for a in args:
                    _learned_worker(a)
            else:
                pool.map(_learned_worker, args)
            n_flown += len(todo)
            el = time.time() - t_all
            left = n_all - n_done0 - n_flown
            log(f"  [learned] {ctrl} {load} seed {seed} mppi {rep}: {len(todo)} flights in "
                f"{(time.time() - t0) / 60:.1f} min; {left} flights left (~{el / n_flown * left / 60:.0f} min)")
    finally:
        if pool is not None:
            pool.close()
            pool.join()


def collect_flights(layout: Layout, write: bool = True):
    """Rebuild 04_flights/flights.csv from every result JSON."""
    import pandas as pd
    rows = []
    for p in sorted((layout.flights / "results").rglob("*.json")):
        try:
            rows.append(json.loads(p.read_text()))
        except (json.JSONDecodeError, OSError):
            print(f"  [collect] unreadable {p}, skipped")
    df = pd.DataFrame(rows)
    if write and len(df):
        write_csv(df, layout.flights / "flights.csv")
    return df


def load_telemetry(layout: Layout, flight_id: str) -> Dict[str, np.ndarray]:
    with np.load(layout.telemetry_path(flight_id)) as z:
        return {k: z[k] for k in z.files}


def eta_report(cfg: StudyConfig, jobs: Sequence[dict], layout: Layout, s_per_classical: float,
               s_per_learned_batch_flight: float, workers_cpu: int, workers_gpu: int) -> dict:
    """Remaining-work estimate for the notebook (inputs are measured rates)."""
    n_c = sum(j["controller"].startswith("classical") and not flight_done(layout, j) for j in jobs)
    n_l = sum(j["controller"] in LEARNED and not flight_done(layout, j) for j in jobs)
    return {"classical_left": n_c, "learned_left": n_l,
            "classical_h": n_c * s_per_classical / max(1, workers_cpu) / 3600,
            "learned_h": n_l * s_per_learned_batch_flight / max(1, workers_gpu) / 3600}
