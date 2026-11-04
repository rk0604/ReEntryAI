"""
jepa/study_v5.py

The v5 study's plumbing (jepa/v5_plan.md sections 3, 4, 9).

    StudyConfigV5      every number of the plan, the compute tier (D1), smoke() for the CPU smoke run
    LayoutV5           the v5 Drive folder (00_study ... 06_report, 01_benchmark, 01_checks, 02_policies,
                       03_tuning). Every write goes through its guard (refuses anything outside the v5 folder or
                       under a read-only input: the v3, v4, data and legacy RL folders), is fsync'd and size-checked
                       before the rename; large files are staged on the local disk first (save_large)
    repair             Setup's self-repair (section 9.1 rule 8): stale *.tmp deleted, unreadable flight JSONs
                       deleted (re-flown), .done markers whose artifacts are missing removed
    freeze_config_v5   00_study/study_config.json, frozen at first use (a later session with other settings raises)
    conditions, *_jobs every flight of every stage, with the v4 job fields plus suite, block, target range, training
                       seed; ORB-1's sequential stopping rule is recomputed from the files on disk
    run_jobs           the resumable queues: guidance laws, RL policies and classical on a CPU process pool (one
                       lockstep batch of flights per task), MPPI arms one model run at a time, its flights in
                       lockstep shards; every flight is written the moment it lands, finished flights are skipped
    gates              S0, V1-V4, E0 (01_checks/<gate>.json and <gate>.done)
    tuning             DEV48 only (section 4.5): Vf from classical, random search for FNPEG and PredGuid
                       (03_tuning/<method>/config_<j>.json), RL configurations at a reduced budget
    train_e3_stage, train_rl_stage, open_loop_stage
                       wrappers around models_v5 and rl_v5 (resumable, done markers ckpt.pt / policy.pt)
    legacy_inventory   MyDrive/ReEntryAI/{models,checkpoints,data} -> 00_study/legacy_rl_inventory.csv
    status             done / total per stage and arm, rates measured on the first 20 flights, hours and Colab
                       units left (00_study/status.json, status_history.csv)
    timing_stage       compute per decision (RQ5): median and p99 ms, CPU single thread and GPU
"""

from __future__ import annotations

import io
import json
import math
import multiprocessing as mp
import os
import pickle
import shutil
import sys
import tempfile
import time
import zlib
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

from jepa import benchmark_v5 as B
from jepa import fly_v3 as F
from jepa import fly_v4 as F4
from jepa import fly_v5 as F5
from jepa import models_v5 as M5
from jepa import study_v3 as S
from jepa import study_v4 as S4
from jepa import testset_v3 as TS
from jepa.data import NOMINAL_DOMAIN
from jepa.models_v5 import fail_point

ATM = tuple(TS.ATMOSPHERES)
G = "4.5g"
GUIDANCE_KINDS = ("classical", "fnpeg", "predguid", "policy", "legacy_rl")
RL_TIER1 = ("dagger_fnpeg", "bcppo_fnpeg", "ppo_hist")
RL_TIER2 = ("sac_hist", "ppo_matched")
ORB1_ARMS = ("classical", "fnpeg", "predguid", "mppi_filtered", "v2single")
EXPLORE_ARMS = ("fnpeg_180", "predguid_180", "fnpeg_oracle")
# Batching and worker knobs: change speed and memory, never results; free to differ between sessions.
RUNTIME_KEYS = ("max_rows", "max_rows_physics", "v2_chunk", "guidance_batch", "units_per_hour", "tier")
# Approximate Colab compute units per hour by runtime (the status cell's projection only).
UNITS_PER_HOUR = {"A100": 11.8, "L4": 4.8, "T4": 1.8, "V100": 4.9, "cpu": 0.08}


# ============================================================================
# Configuration
# ============================================================================
def _e3_default() -> dict:
    return asdict(M5.E3TrainCfg())


@dataclass
class StudyConfigV5:
    """Every number of v5_plan.md. Frozen to 00_study/study_config.json on first use."""
    mode: str = "full"                                  # "full" | "smoke" (every stage, tiny, numbers meaningless)
    tier: int = 1                                       # D1: Tier 1 first; Tier 2 after the user read Tier 1
    seeds: Tuple[int, ...] = (0, 1, 2, 3, 4)            # learned models (v3's v2 seeds, E3 seeds, MPPI noise indices)
    rl_seeds: Tuple[int, ...] = (0, 1, 2)
    a256_cases: Optional[Tuple[int, ...]] = None        # None: all 64
    dev48_cases: Optional[Tuple[int, ...]] = None       # None: all 12
    orb1_per_block: Optional[int] = None                # None: all 1,000 per block (smoke: the first k)
    atmospheres: Tuple[str, ...] = ATM
    g_primary: float = 4.5
    split_seed: int = 0
    n_episodes: Optional[int] = None
    # ---- guidance laws: the settings before tuning (FNPEGConfig / PredGuidConfig overrides) ----
    fnpeg: dict = field(default_factory=dict)
    predguid: dict = field(default_factory=dict)
    tune_n: int = 24                                    # configurations per guidance law (config 0 = defaults)
    tune_seed: int = 20261104
    tune_peak_g_max: float = 4.6
    rl_tune_n: dict = field(default_factory=lambda: {"dagger_fnpeg": 4, "bcppo_fnpeg": 4, "ppo_hist": 6,
                                                     "sac_hist": 4})
    rl_tune_frac: float = 0.1                           # budget fraction of an RL tuning run
    # ---- gates (section 4.0) ----
    gate_cases: Tuple[int, ...] = tuple(range(0, 64, 4))
    gates_stop: bool = True                             # a failed gate stops the notebook (smoke: reported only)
    v1_miss_km: float = 2.0
    v1_min_frac: float = 15 / 16
    v1_iter_max: int = 6
    v1_iter_frac: float = 0.99
    v2_amax_g: float = 3.0
    v2_trigger_g: float = 3.5
    v2_tol_g: float = 0.15
    v2_sweep_g: Tuple[float, ...] = (2.7, 2.5, 2.3)
    v2_sweep_n: int = 2
    v3_n_cases: int = 8
    v3_tol: float = 0.03
    v3_after_s: float = 150.0
    v4_miss_km: float = 3.5
    v4_min_frac: float = 15 / 16
    e0_conditions: Tuple[Tuple[int, str], ...] = ((0, "nominal"), (16, "nominal"), (40, "mismatch_high"))
    e0_samples: int = 32
    s0_feasible_min: float = 0.97
    s0_horizon_s: float = 940.0
    # ---- E3 ----
    e3_arms: Tuple[str, ...] = M5.ARMS_E3
    e3_train: dict = field(default_factory=_e3_default)
    e3_sizes: dict = field(default_factory=dict)        # {} = each arm sized to R0's parameter count
    r0_arch: dict = field(default_factory=dict)         # models_v2 defaults (R0)
    lr_mults: Tuple[float, ...] = (1.0, 0.3, 3.0)
    lr_pilot_frac: float = 0.1
    ttg_s: Tuple[int, ...] = (400, 300, 200, 100, 50)
    # ---- E4 ----
    rl: dict = field(default_factory=dict)              # RLCfg overrides for every RL arm
    rl_arms: dict = field(default_factory=lambda: {
        "bcppo_fnpeg": {"lr": 3e-5, "kl_coef": 0.1, "critic_warmup_updates": 20},
        "sac_hist": {"budget_decisions": 20_000_000, "n_envs": 256},
        "ppo_matched": {"budget_hours": 24.0}})
    rl_curve_cases: int = 12                            # DEV48 nominal cases for the full-simulator learning curve
    legacy_heldout_n: int = 256
    fast_eval_episodes: int = 2048
    # ---- analysis ----
    n_boot: int = 10_000
    success_km: float = 5.0
    log_miss_floor_km: float = 0.01
    equiv_log: float = 0.15
    thresholds_km: Tuple[float, ...] = tuple(np.round(np.arange(0.5, 25.01, 0.5), 2).tolist())
    # ---- runtime ----
    max_rows: int = 2048
    max_rows_physics: int = 16384
    v2_chunk: int = 512
    guidance_batch: int = 16
    units_per_hour: dict = field(default_factory=lambda: dict(UNITS_PER_HOUR))

    @staticmethod
    def smoke() -> "StudyConfigV5":
        """Every stage on a CPU in minutes, on top of the v3 and v4 smoke folders (smoke_v4.py)."""
        v3 = S.StudyConfig.smoke()
        tiny = dict(d_latent=16, enc_dim=32, enc_layers=2, enc_heads=2, enc_band=5, pred_dim=48, pred_layers=2,
                    pred_heads=2, act_chunk=4, mlp_ratio=2)
        rl = dict(n_envs=16, hidden=(32, 32), n_steps=8, epochs=2, minibatch=64, budget_decisions=512,
                  budget_hours=1.0, buffer=4096, batch=64, learning_starts=64, updates_per_step=4.0, dagger_iters=2,
                  dagger_rollouts=4, dagger_epochs=2, bc_episodes=4, bc_epochs=2, sup_batch=256, save_every_s=30.0,
                  curve_points=2, env_pool=64, full_episodes=2, full_per_update=2, full_workers=2,
                  critic_warmup_updates=1)
        return StudyConfigV5(
            mode="smoke", seeds=v3.seeds, rl_seeds=(0,), a256_cases=v3.cases, dev48_cases=(0,), orb1_per_block=1,
            n_episodes=v3.n_episodes, tune_n=2, rl_tune_n={"dagger_fnpeg": 1, "bcppo_fnpeg": 1, "ppo_hist": 2,
                                                           "sac_hist": 1},
            rl_tune_frac=0.5, gate_cases=v3.cases, gates_stop=False, v2_sweep_g=(2.7,), v2_sweep_n=1, v3_n_cases=1,
            e0_conditions=((0, "nominal"), (16, "mismatch_high")), e0_samples=16,
            e3_train=dict(_e3_default(), steps=8, batch=8, warmup=2, curriculum=((0.0, 8), (0.5, 16)), log_every=4,
                          save_every_s=30.0, amp=False, val_rows=8),
            e3_sizes={"mlp_direct": {"width": 32, "depth": 2}, "gru_direct": {"hidden": 16},
                      "gru_pi": {"hidden": 16, "prober_hidden": 16}},
            r0_arch=tiny, lr_mults=(1.0, 3.0), lr_pilot_frac=0.5, ttg_s=(200, 100), rl=rl,
            rl_arms={"bcppo_fnpeg": {"lr": 3e-5, "kl_coef": 0.1}, "sac_hist": {"budget_decisions": 1024}},
            rl_curve_cases=1, legacy_heldout_n=2, fast_eval_episodes=16, n_boot=200, guidance_batch=4)

    def v3_expected(self) -> S.StudyConfig:
        return S.StudyConfig.smoke() if self.mode == "smoke" else S.StudyConfig()

    def mppi(self) -> dict:
        return dict(self.v3_expected().mppi)

    def horizon(self) -> int:
        return S.horizon_steps(self.v3_expected())

    def stages(self) -> Dict[str, bool]:
        """RUN per notebook stage, derived from the tier (the notebook may override any entry)."""
        t2 = self.tier >= 2
        return {"benchmark": True, "legacy": True, "gates": True, "tuning": True, "A256": True, "ORB1_block1": True,
                "E3_train": True, "E3_fly": True, "E4_train": True, "E4_fly": True, "ORB1_more": t2,
                "tier2_rl": t2, "explore": t2, "open_loop": True, "timing": True, "analysis": True, "report": True}


def git_commit() -> str:
    return S.git_commit()


# ============================================================================
# Layout: guarded, fsync'd, size-checked writes
# ============================================================================
class ReadOnlyViolation(S4.ReadOnlyViolation):
    """A v5 write aimed outside the v5 study folder (or into a read-only input)."""


def _norm(p) -> str:
    """A path as a normalized absolute string (no filesystem access; Windows' \\\\?\\ prefix removed)."""
    s = str(p)
    if s.startswith("\\\\?\\"):
        s = s[4:]
    return os.path.normcase(os.path.abspath(s))


def _inside_s(p: str, root: str) -> bool:
    return p == root or p.startswith(root.rstrip("\\/") + os.sep)


STAGES = {
    "00_study": "frozen settings (study_config.json), status.json, status_history.csv, logs/, legacy RL inventory",
    "01_benchmark": "a copy of configs/benchmark_v5 (A256, ORB-1, DEV48) with SHA256SUMS",
    "01_checks": "gates S0, V1-V4, E0 (<gate>.json + .done) and compute per decision (RQ5)",
    "02_models": "E3 models per arm and seed: ckpt.pt, model.json, train_log.csv, F_curves.png/pdf",
    "02_policies": "RL policies per arm and seed: policy.pt, run.json, learning_curve.csv, F_curves.png/pdf",
    "03_tuning": "DEV48 tuning: config_<j>.json per configuration, chosen.json per method",
    "03_open_loop": "open-loop footprint errors of the E3 models (CSV)",
    "04_flights": "one JSON result + one 1 Hz telemetry .npz per flight; flights.csv",
    "05_analysis": "tables T01-T11, figures F01-F10 (PNG and PDF) with captions.json; cache/",
    "06_report": "SpaceJEPA_v5_report.pdf",
}


class LayoutV5(S4.LayoutV4):
    """Paths of the v5 study folder (v5_plan.md section 9) and its guarded write helpers."""

    def __init__(self, root, read_only: Sequence = ()):
        super().__init__(root, read_only)
        self.bench = self.root / "01_benchmark"
        self.policies = self.root / "02_policies"
        self.tuning = self.root / "03_tuning"
        self.logs = self.study / "logs"
        self._root_s = _norm(self.root)
        self._ro_s = [_norm(p) for p in self.read_only]

    def make(self) -> "LayoutV5":
        super().make()
        for p in (self.bench, self.policies, self.tuning, self.logs):
            self.guard(p)
            p.mkdir(parents=True, exist_ok=True)
        return self

    def guard(self, path) -> Path:
        p = _norm(path)
        if not _inside_s(p, self._root_s):
            raise ReadOnlyViolation(f"v5 refuses to write outside its study folder: {p} (folder {self.root})")
        for ro in self._ro_s:
            if _inside_s(p, ro):
                raise ReadOnlyViolation(f"v5 refuses to write under the read-only input {ro}: {p}")
        return Path(path)

    def atomic_bytes(self, path, data: bytes) -> None:
        """<name>.tmp next to the target, flush + fsync, size check, then os.replace (section 9.1 rule 2)."""
        path = self.guard(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
        with open(tmp, "wb") as f:
            f.write(data)
            f.flush()
            os.fsync(f.fileno())
        if os.path.getsize(tmp) != len(data):
            raise OSError(f"short write: {tmp} has {os.path.getsize(tmp)} of {len(data)} bytes")
        os.replace(tmp, path)

    def save_large(self, data: bytes, path) -> None:
        """Large files: local disk first, then a copy to the Drive .tmp name, size check, rename."""
        path = self.guard(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        stage = Path(tempfile.gettempdir()) / "jepa_v5_stage"
        stage.mkdir(parents=True, exist_ok=True)
        local = stage / f"{os.getpid()}_{path.name}"
        with open(local, "wb") as f:
            f.write(data)
            f.flush()
            os.fsync(f.fileno())
        tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
        shutil.copyfile(local, tmp)
        with open(tmp, "rb+") as f:
            os.fsync(f.fileno())
        if os.path.getsize(tmp) != len(data):
            raise OSError(f"short copy: {tmp} has {os.path.getsize(tmp)} of {len(data)} bytes")
        os.replace(tmp, path)
        local.unlink(missing_ok=True)

    def torch_save(self, obj, path) -> None:
        buf = io.BytesIO()
        torch.save(obj, buf)
        self.save_large(buf.getvalue(), path)

    def write_flight(self, rec: dict) -> None:
        rec["result"].setdefault("git_commit", _GIT.get())
        super().write_flight(rec)

    def policy_dir(self, arm: str, seed: int) -> Path:
        return self.policies / arm / f"seed_{int(seed)}"

    def logger(self, stage: str) -> Callable:
        """print + append to 00_study/logs/<stage>.log (Colab's cell output is lost on a disconnect, the log is not)."""
        path = self.guard(self.logs / f"{stage}.log")

        def log(*args):
            msg = " ".join(str(a) for a in args)
            print(msg, flush=True)
            try:
                path.parent.mkdir(parents=True, exist_ok=True)
                with open(path, "a", encoding="utf-8") as f:
                    f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}\n")
            except OSError:
                pass
        return log

    def mark_done(self, path, artifacts: Sequence = (), info: Optional[dict] = None) -> None:
        """A .done marker listing its artifacts (relative paths); repair removes it when one goes missing."""
        rel = [os.path.relpath(_norm(a), self._root_s) for a in artifacts]
        self.write_json(path, {"artifacts": rel, "time": time.strftime("%Y-%m-%d %H:%M:%S"),
                               "git_commit": _GIT.get(), **({"info": info} if info else {})})


class _Git:
    def __init__(self):
        self.v = None

    def get(self) -> str:
        if self.v is None:
            self.v = os.environ.get("JEPA_V5_GIT") or git_commit()
            os.environ["JEPA_V5_GIT"] = self.v                  # inherited by the pool workers
        return self.v


_GIT = _Git()


def repair(L: LayoutV5, log=print) -> dict:
    """
    Section 9.1 rule 8, at Setup: delete stale *.tmp files; delete flight JSONs written since the last repair that
    do not parse (they are re-flown); remove .done markers whose artifacts are missing. Checkpoints that do not load
    fall back to .prev when they are loaded (models_v5.load_state).
    """
    t0 = time.time()
    sp = L.study / "repair_state.json"
    last = 0.0
    if sp.exists():
        try:
            last = float(json.loads(sp.read_text())["time"])
        except (OSError, ValueError, KeyError):
            last = 0.0
    out = {"tmp_deleted": 0, "json_deleted": 0, "done_removed": 0, "files_seen": 0}
    res_root = _norm(L.flights / "results")
    for dirpath, _dirs, files in os.walk(L.root):
        for name in files:
            out["files_seen"] += 1
            p = os.path.join(dirpath, name)
            if name.endswith(".tmp"):
                L.guard(p)
                os.remove(p)
                out["tmp_deleted"] += 1
            elif name.endswith(".done"):
                try:
                    arts = json.loads(Path(p).read_text()).get("artifacts", [])
                    ok = all((L.root / a).exists() for a in arts)
                except (OSError, ValueError):
                    ok = False
                if not ok:
                    os.remove(L.guard(p))
                    out["done_removed"] += 1
            elif name.endswith(".json") and _inside_s(_norm(p), res_root) and os.path.getmtime(p) >= last - 600:
                try:
                    json.loads(Path(p).read_text())
                except (OSError, ValueError):
                    os.remove(L.guard(p))
                    out["json_deleted"] += 1
    L.write_json(sp, {"time": t0, "result": out})
    if out["tmp_deleted"] or out["json_deleted"] or out["done_removed"]:
        log(f"[repair] {out}")
    out["sec"] = round(time.time() - t0, 1)
    return out


def freeze_config_v5(cfg: StudyConfigV5, L: LayoutV5, inputs: Optional[dict] = None) -> dict:
    """00_study/study_config.json on first use; later sessions must use the same settings (except RUNTIME_KEYS)."""
    path = L.study / "study_config.json"
    want = json.loads(json.dumps(asdict(cfg), default=S._json_default))
    if path.exists():
        have = json.loads(path.read_text())
        diff = {k: (have["config"].get(k), v) for k, v in want.items()
                if k not in RUNTIME_KEYS and have["config"].get(k) != v}
        if diff:
            raise ValueError(f"settings differ from the frozen {path}: {diff}. Use a new study folder.")
        return have
    doc = {"config": want, "git_commit": _GIT.get(), "created": time.strftime("%Y-%m-%d %H:%M:%S"),
           "environment": S.environment_info(), "inputs": inputs or {}, "benchmark_sha256": B.verify()}
    L.write_json(path, doc)
    return doc


def copy_benchmark(L: LayoutV5) -> None:
    """The frozen suites to 01_benchmark/ (checksums verified first)."""
    sums = B.verify()
    for name in list(sums) + ["SHA256SUMS"]:
        dst = L.bench / name
        if not dst.exists():
            L.atomic_bytes(dst, (B.BENCH_DIR / name).read_bytes())


def v4_layout(root) -> Optional[S4.LayoutV4]:
    """The v4 folder (R0's flights and open loop for RQ2, v2's load curve), or None if it is not there. Read only:
    v5 never writes through it (its own guard refuses any path under it)."""
    if root and (Path(root) / "00_study" / "study_config.json").exists():
        return S4.LayoutV4(root)
    return None


# ============================================================================
# Conditions and jobs
# ============================================================================
_SUITES: Dict[str, List[dict]] = {}


def conditions(cfg: StudyConfigV5, suite: str) -> List[dict]:
    """The conditions of a suite this config flies (A256, DEV48, ORB1, HELD256)."""
    if suite not in _SUITES:
        if suite == "HELD256":
            import pandas as pd
            df = pd.read_csv(SIM_DIR / "configs" / "heldout_testset.csv", comment="#")
            _SUITES[suite] = [{"case": int(r["idx"]), "suite": "HELD256", "set": "heldout", "atmosphere": "nominal",
                               "ic": TS.ic_of(r), "domain": dict(NOMINAL_DOMAIN), "config_path": "configs/default.json"}
                              for r in df.to_dict("records")]
        else:
            _SUITES[suite] = B.load_suite(suite)
    out = []
    for c in _SUITES[suite]:
        c = dict(c)
        k = int(c["case"])
        if suite == "A256":
            if (cfg.a256_cases is not None and k not in cfg.a256_cases) or c["atmosphere"] not in cfg.atmospheres:
                continue
            c["feasible"] = c.get(f"feasible_{c['atmosphere']}")
        elif suite == "DEV48":
            if cfg.dev48_cases is not None and k not in cfg.dev48_cases:
                continue
            c["feasible"] = c.get(f"feasible_{c['atmosphere']}")
        elif suite == "ORB1":
            c["block"] = B.block_of(k)
            if cfg.orb1_per_block is not None and k % (B.ORB1Spec().n_cases // B.ORB1Spec().n_blocks) >= cfg.orb1_per_block:
                continue
        elif suite == "HELD256":
            if k >= cfg.legacy_heldout_n:
                continue
        out.append(c)
    return out


def _base(cond: dict, arm: str, stage: str) -> dict:
    f = cond.get("feasible")
    return {"case": int(cond["case"]), "set": cond.get("set", cond["suite"]), "atmosphere": cond["atmosphere"],
            "ic": dict(cond["ic"]), "domain": dict(cond["domain"]),
            "config_path": cond.get("config_path") or "configs/default.json", "suite": cond["suite"],
            "controller": arm, "stage": stage, "block": cond.get("block"),
            "target_index": None if cond.get("target_index") is None else int(cond["target_index"]),
            "feasible": None if f is None or (isinstance(f, float) and math.isnan(f)) else bool(f),
            "in_training_domain": cond.get("in_training_domain"), "sigma_max_deg": 90.0}


def guidance_job(cond: dict, arm: str, stage: str, law: str, law_cfg: dict, sigma_max: float = 90.0) -> dict:
    j = _base(cond, arm, stage)
    j.update({"load": "own", "g_limit": None, "seed": None, "mppi_rep": None, "sigma_max_deg": float(sigma_max),
              "model": {"kind": law, "cfg": dict(law_cfg)}})
    if law == "predguid" and cond["suite"] == "ORB1":
        j["design_ic"] = B.ORB1Spec().nominal_ic()      # preflight design: the mission's nominal entry, nominal model
    return j


def classical_job(cond: dict, stage: str, arm: str = "classical") -> dict:
    j = _base(cond, arm, stage)
    j.update({"load": "own", "g_limit": None, "seed": None, "mppi_rep": None, "classical_knows_truth": False,
              "model": {"kind": "classical"}})
    return j


def policy_job(cond: dict, arm: str, stage: str, path: str, seed: Optional[int], legacy: bool = False,
               env_kwargs: Optional[dict] = None) -> dict:
    j = _base(cond, arm, stage)
    j.update({"load": "own", "g_limit": None, "seed": None if seed is None else int(seed), "mppi_rep": None,
              "train_seed": None if seed is None else int(seed),
              "model": {"kind": "legacy_rl" if legacy else "policy", "path": str(path)}})
    if env_kwargs:
        j["env_kwargs"] = dict(env_kwargs)
    return j


def mppi_job(cfg: StudyConfigV5, cond: dict, arm: str, load: str, g: float, seed: int, model: dict, stage: str) -> dict:
    H = cfg.horizon()
    s = F.crn_seed(seed, int(cond["case"]), cond["atmosphere"], 0)
    j = _base(cond, arm, stage)
    j.update({"load": load, "g_limit": float(g), "seed": int(seed), "mppi_rep": 0, "mppi_seed": s, "horizon": H,
              "train_seed": int(seed), "mppi": {"seed": s, "g_limit": float(g), "horizon": H, "w_unc_km": 0.0},
              "model": dict(model)})
    return j


def v2single_model(L3: S.Layout, seed: int) -> dict:
    return {"kind": "v2", "ckpt": str(L3.ckpt("v2", seed)), "members": [0]}


def e3_model(L: LayoutV5, arm: str, seed: int) -> dict:
    ck = str(L.ckpt(arm, seed))
    return {"kind": "v2", "ckpt": ck, "members": None} if arm == "r0_nojepa" else {"kind": "e3", "ckpt": ck,
                                                                                    "members": None}


def orb1_seed(cfg: StudyConfigV5, cond: dict) -> int:
    """ORB-1: case i is flown by training seed i mod 5 (the CSV's train_seed), noise index likewise."""
    return int(cfg.seeds[int(cond.get("train_seed", cond["case"])) % len(cfg.seeds)])


def tuned(L: LayoutV5, method: str, cfg: StudyConfigV5) -> dict:
    """The chosen configuration of a tuned method (defaults until tuning has chosen; is_tuned tells which)."""
    p = L.tuning / method / "chosen.json"
    if p.exists():
        return dict(json.loads(p.read_text())["config"])
    if method == "fnpeg":
        return dict(cfg.fnpeg, **_vf(L))
    if method == "predguid":
        return dict(cfg.predguid)
    return {}


def is_tuned(L: LayoutV5, method: str) -> bool:
    return (L.tuning / method / "chosen.json").exists()


def _vf(L: LayoutV5) -> dict:
    p = L.tuning / "fnpeg" / "vf.json"
    return {"vf_mps": float(json.loads(p.read_text())["vf_mps"])} if p.exists() else {}


def arm_job(cfg: StudyConfigV5, L: LayoutV5, L3: S.Layout, cond: dict, arm: str, stage: str,
            seed: Optional[int] = None) -> dict:
    """One flight of a named arm on a condition (guidance laws and MPPI arms of E1/E2)."""
    if arm == "classical":
        return classical_job(cond, stage)
    if arm in ("fnpeg", "fnpeg_180", "fnpeg_oracle"):
        c = tuned(L, "fnpeg", cfg)
        if arm == "fnpeg_180":
            c["sigma_max_deg"] = 180.0
        if arm == "fnpeg_oracle":
            c.update({"oracle": True, "filters": False})
        j = guidance_job(cond, arm, stage, "fnpeg", c, 180.0 if arm.endswith("180") else 90.0)
        j["untuned"] = not is_tuned(L, "fnpeg")          # test flights wait until tuning has chosen (section 4.5)
        return j
    if arm in ("predguid", "predguid_180"):
        c = tuned(L, "predguid", cfg)
        if arm == "predguid_180":
            c["sigma_max_deg"] = 180.0
        j = guidance_job(cond, arm, stage, "predguid", c, 180.0 if arm.endswith("180") else 90.0)
        j["untuned"] = not is_tuned(L, "predguid")
        return j
    if seed is None and arm in RL_TIER1 + RL_TIER2:      # ORB-1: case i flies RL training seed i mod 3
        seed = int(cfg.rl_seeds[int(cond.get("train_seed", cond["case"])) % len(cfg.rl_seeds)])
    if seed is None:
        seed = orb1_seed(cfg, cond)
    if arm == "mppi_filtered":
        return mppi_job(cfg, cond, arm, G, cfg.g_primary, seed, {"kind": "physics_filtered", "ckpt": None,
                                                                 "members": None}, stage)
    if arm == "mppi_truth":
        return mppi_job(cfg, cond, arm, G, cfg.g_primary, seed, {"kind": "physics_truth", "ckpt": None,
                                                                 "members": None}, stage)
    if arm == "v2single":
        return mppi_job(cfg, cond, arm, G, cfg.g_primary, seed, v2single_model(L3, seed), stage)
    if arm in M5.ARMS_E3:
        return mppi_job(cfg, cond, arm, G, cfg.g_primary, seed, e3_model(L, arm, seed), stage)
    if arm in RL_TIER1 + RL_TIER2:
        return policy_job(cond, arm, stage, str(L.policy_dir(arm, seed) / "policy.pt"), seed)
    raise ValueError(f"unknown arm {arm!r}")


def jobs_a256(cfg: StudyConfigV5, L: LayoutV5, L3: S.Layout, matched: bool = True) -> List[dict]:
    """Tier 1 A256: fnpeg, predguid, mppi_filtered (5 noise indices), v2single@fnpeg_matched (once fnpeg has flown)."""
    conds = conditions(cfg, "A256")
    jobs = [arm_job(cfg, L, L3, c, a, "A256") for a in ("fnpeg", "predguid") for c in conds]
    for s in cfg.seeds:
        jobs += [arm_job(cfg, L, L3, c, "mppi_filtered", "A256", s) for c in conds]
    if matched:
        jobs += jobs_matched(cfg, L, L3)
    return jobs


def jobs_matched(cfg: StudyConfigV5, L: LayoutV5, L3: S.Layout) -> List[dict]:
    """v2single with each condition's soft g limit = fnpeg's realized peak g there ([] until fnpeg has flown A256)."""
    conds = conditions(cfg, "A256")
    try:
        pk = F5.fnpeg_peak_g(L, "A256", conds)
    except FileNotFoundError:
        return []
    return [mppi_job(cfg, c, "v2single", "fnpeg_matched", pk[(int(c["case"]), c["atmosphere"])], s,
                     v2single_model(L3, s), "A256") for s in cfg.seeds for c in conds]


def orb1_jobs(cfg: StudyConfigV5, L: LayoutV5, L3: S.Layout, arms: Sequence[str], blocks: Sequence[int],
              stage: str, log=print) -> List[dict]:
    """ORB-1 flights of `arms` in `blocks`; block b > 0 only once every earlier block is flown and the arm is eligible."""
    conds = conditions(cfg, "ORB1")
    by_block: Dict[int, List[dict]] = {}
    for c in conds:
        by_block.setdefault(int(c["block"]), []).append(c)
    jobs = []
    for arm in arms:
        for b in sorted(blocks):
            js = [arm_job(cfg, L, L3, c, arm, stage) for c in by_block.get(b, [])]
            if b > 0:
                ok, why = orb1_eligible_files(cfg, L, L3, arm, b, by_block)
                if not ok:
                    log(f"  [ORB-1] {arm} block {b + 1}: {why}")
                    continue
            jobs += js
    return jobs


def orb1_eligible_files(cfg, L, L3, arm, b, by_block) -> Tuple[bool, str]:
    """Sequential stopping (section 3.2) from the result files: earlier blocks complete and < 2 non-successes."""
    bad, missing = 0, 0
    for bb in range(b):
        for c in by_block.get(bb, []):
            p = L.result_path(F5.flight_id_v5(arm_job(cfg, L, L3, c, arm, "ORB1")))
            r = _read_json(p)
            if r is None:
                missing += 1
            elif not (r.get("outcome") == "drogue" and float(r.get("miss_km", 1e9)) <= 10.0):
                bad += 1
    if missing:
        return False, f"waits for {missing} flights of earlier blocks"
    if bad >= 2:
        return False, f"stopped: {bad} flights of earlier blocks missed 10 km (failed Q1)"
    return True, "eligible"


_JSON_CACHE: Dict[str, Tuple[float, dict]] = {}


def _read_json(p: Path) -> Optional[dict]:
    try:
        m = os.path.getmtime(p)
    except OSError:
        return None
    k = str(p)
    if k in _JSON_CACHE and _JSON_CACHE[k][0] == m:
        return _JSON_CACHE[k][1]
    try:
        d = json.loads(Path(p).read_text())
    except (OSError, ValueError):
        return None
    _JSON_CACHE[k] = (m, d)
    return d


def jobs_e3(cfg: StudyConfigV5, L: LayoutV5, L3: S.Layout) -> List[dict]:
    conds = conditions(cfg, "A256")
    return [arm_job(cfg, L, L3, c, arm, "E3", s) for arm in cfg.e3_arms for s in cfg.seeds for c in conds]


def jobs_e4(cfg: StudyConfigV5, L: LayoutV5, L3: S.Layout, arms: Sequence[str]) -> List[dict]:
    """RL policies on A256 (seed-major per arm) and their snapshots on the DEV48 learning-curve cases."""
    conds = conditions(cfg, "A256")
    jobs = [arm_job(cfg, L, L3, c, arm, "E4", s) for arm in arms for s in cfg.rl_seeds for c in conds]
    dev = [c for c in conditions(cfg, "DEV48") if c["atmosphere"] == "nominal"][:cfg.rl_curve_cases]
    for arm in arms:
        for s in cfg.rl_seeds:
            snap = L.policy_dir(arm, s) / "snapshots"
            if not snap.exists():
                continue
            for p in sorted(snap.glob("policy_*.pt")):
                jobs += [policy_job(c, f"{arm}@{p.stem.replace('policy_', '')}", "E4_curve", str(p), s) for c in dev]
    return jobs


def jobs_explore(cfg: StudyConfigV5, L: LayoutV5, L3: S.Layout) -> List[dict]:
    conds = conditions(cfg, "A256")
    return [arm_job(cfg, L, L3, c, a, "explore") for a in EXPLORE_ARMS for c in conds]


def done_ids(L: LayoutV5, jobs: Sequence[dict]) -> set:
    """Flight ids among jobs whose result exists (one os.scandir per result folder)."""
    by_dir: Dict[Path, Dict[str, str]] = {}
    for j in jobs:
        fid = F5.flight_id_v5(j)
        p = L.result_path(fid)
        by_dir.setdefault(p.parent, {})[p.name] = fid
    out = set()
    for d, names in by_dir.items():
        try:
            with os.scandir(d) as it:
                for e in it:
                    if e.name in names:
                        out.add(names[e.name])
        except FileNotFoundError:
            pass
    return out


def flight_done(L: LayoutV5, job: dict) -> bool:
    return L.result_path(F5.flight_id_v5(job)).exists()


# ============================================================================
# Queues
# ============================================================================
def _relabel(job: dict, rec: dict) -> dict:
    """A classical flight's record (fly_v4) with the v5 id, labels and section-5 metrics."""
    res, tel = rec["result"], rec["telemetry"]
    res["v4_flight_id"] = res["flight_id"]
    res.update(F5.result_extras_v5(job))
    res.update(F5.v5_metrics(res, tel, 90.0, ()))
    res["plan_ms_median"], res["plan_ms_p99"] = res.get("plan_ms_mean"), res.get("plan_ms_max")
    return rec


def _cpu_worker(arg) -> int:
    """Pool task: classical flights one by one, or one lockstep batch of guidance-law / policy flights."""
    shard, root, ro = arg
    torch.set_num_threads(1)
    L = LayoutV5(root, ro)
    m = shard[0]["model"]
    if m["kind"] == "classical":
        for j in shard:
            L.write_flight(_relabel(j, F4.fly_classical_v4(j)))
        return len(shard)
    if m["kind"] in ("policy", "legacy_rl"):
        law = F5.make_law("policy", {"path": m["path"], "legacy": m["kind"] == "legacy_rl"})
    else:
        law = F5.make_law(m["kind"], m["cfg"])
    F5.run_guidance_lockstep(shard, law, on_done=L.write_flight, verbose=False)
    return len(shard)


def _mppi_worker(arg) -> int:
    shard, model, device, mppi, max_rows, chunk, root, ro, threads, verbose, prog = arg
    torch.set_num_threads(int(threads))
    adapter = F5.load_adapter_v5(model["kind"], model.get("ckpt"), device, model.get("members"), chunk, keep_only=True)
    L = LayoutV5(root, ro)

    def report(p):
        if prog:
            try:
                S.atomic_bytes(prog, json.dumps(p).encode())
            except OSError:
                pass

    F5.run_mppi_v5(shard, adapter, mppi, n_sub=int(mppi.get("n_sub", 2)), max_rows=max_rows, on_done=L.write_flight,
                   verbose=verbose, progress=report)
    return len(shard)


def _cpu_key(j: dict) -> str:
    return json.dumps([j["controller"], j["model"], j.get("env_kwargs")], sort_keys=True, default=str)


_COUNT = {"cpu": 0, "mppi": 0}      # tasks finished in this process (the smoke test's fail points)


def run_cpu(cfg: StudyConfigV5, jobs: Sequence[dict], L: LayoutV5, workers: int = 1, log=print,
            label: str = "cpu") -> int:
    """Guidance laws, policies and classical flights (CPU pool); returns flights flown."""
    jobs = [j for j in jobs if j["model"]["kind"] in GUIDANCE_KINDS]
    done = done_ids(L, jobs)
    todo = [j for j in jobs if F5.flight_id_v5(j) not in done]
    wait = [j for j in todo if j.get("untuned") or (j["model"]["kind"] in ("policy", "legacy_rl")
                                                    and not os.path.exists(j["model"]["path"]))]
    todo = [j for j in todo if not (j.get("untuned") or (j["model"]["kind"] in ("policy", "legacy_rl")
                                                         and not os.path.exists(j["model"]["path"])))]
    log(f"[{label}] {len(jobs) - len(todo) - len(wait)}/{len(jobs)} done, {len(todo)} to fly on {workers} CPU "
        f"worker(s)" + (f", {len(wait)} wait for tuning or training" if wait else ""))
    if not todo:
        return 0
    groups: Dict[str, List[dict]] = {}
    for j in todo:
        groups.setdefault(_cpu_key(j), []).append(j)
    tasks = []
    for g in groups.values():
        n = 1 if g[0]["model"]["kind"] == "classical" else max(1, int(cfg.guidance_batch))
        tasks += [g[i:i + n] for i in range(0, len(g), n)]
    ro = [str(p) for p in L.read_only]
    args = [(t, str(L.root), ro) for t in tasks]
    t0, n = time.time(), 0
    pool = mp.get_context("spawn").Pool(int(workers)) if workers > 1 else None
    it = pool.imap(_cpu_worker, args) if pool else map(_cpu_worker, args)
    step = max(1, len(tasks) // 20)
    try:
        for i, k in enumerate(it, 1):
            n += k
            _COUNT["cpu"] += 1
            fail_point("cpu_flights", _COUNT["cpu"])
            if i % step == 0 or i == len(tasks):
                el = time.time() - t0
                log(f"  [{label}] {n}/{len(todo)} flights in {el / 60:.1f} min (~{el / n * (len(todo) - n) / 60:.0f} min left)")
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    return n


def run_mppi(cfg: StudyConfigV5, jobs: Sequence[dict], L: LayoutV5, device: str = "cuda", workers: int = 1,
             log=print, label: str = "mppi", report_every_s: float = 300.0) -> int:
    """MPPI arms, one model run (arm, load, seed, model) at a time, its flights in lockstep shards."""
    jobs = [j for j in jobs if j["model"]["kind"] not in GUIDANCE_KINDS]
    done = done_ids(L, jobs)
    gs: Dict[tuple, List[dict]] = {}
    for j in jobs:
        m = j["model"]
        gs.setdefault((j["suite"], j["controller"], j["load"], j["seed"], m["kind"], m.get("ckpt"),
                       tuple(m["members"]) if m.get("members") is not None else None), []).append(j)
    ready = [(k, g) for k, g in gs.items() if g[0]["model"].get("ckpt") is None or os.path.exists(g[0]["model"]["ckpt"])]
    n_all = sum(len(g) for _, g in ready)
    n_done0 = sum(F5.flight_id_v5(j) in done for _, g in ready for j in g)
    log(f"[{label}] {n_done0}/{n_all} done in {len(ready)} model runs ({len(gs) - len(ready)} wait for training); "
        f"{workers} worker(s) on {device}")
    mppi = cfg.mppi()
    pool, t_all, n_flown = None, time.time(), 0
    ro = [str(p) for p in L.read_only]
    try:
        for key, g in ready:
            todo = [j for j in g if F5.flight_id_v5(j) not in done]
            if not todo:
                continue
            if pool is None and workers > 1:
                pool = mp.get_context("spawn").Pool(int(workers))
            suite, arm, load, seed, kind = key[:5]
            t0 = time.time()
            shards = [sh for sh in (todo[i::max(1, workers)] for i in range(max(1, workers))) if sh]
            S4.PROGRESS_DIR.mkdir(parents=True, exist_ok=True)
            progs = [S4.PROGRESS_DIR / f"v5_worker{i}.json" for i in range(len(shards))]
            for f in progs:
                f.unlink(missing_ok=True)
            rows = cfg.max_rows_physics if kind.startswith("physics") else cfg.max_rows
            threads = 1 if workers > 1 or device == "cuda" else max(1, (os.cpu_count() or 2) // 2)
            args = [(sh, g[0]["model"], device, mppi, rows, cfg.v2_chunk, str(L.root), ro, threads, pool is None,
                     str(pf) if pool is not None else None) for sh, pf in zip(shards, progs)]
            if pool is None:
                for a in args:
                    _mppi_worker(a)
            else:
                res = pool.map_async(_mppi_worker, args)
                while True:
                    try:
                        res.get(timeout=report_every_s)
                        break
                    except mp.TimeoutError:
                        log(f"    ... {suite} {arm} {load} seed {seed}: {(time.time() - t0) / 60:.0f} min, "
                            f"{S._progress_line(progs, len(todo))}")
            n_flown += len(todo)
            _COUNT["mppi"] += 1
            fail_point("mppi_flights", _COUNT["mppi"])
            el = time.time() - t_all
            left = n_all - n_done0 - n_flown
            log(f"  [{label}] {suite} {arm} {load} seed {seed}: {len(todo)} flights in {(time.time() - t0) / 60:.1f} "
                f"min; {left} left in this call (~{el / n_flown * left / 60:.0f} min)")
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    return n_flown


def run_jobs(cfg: StudyConfigV5, jobs: Sequence[dict], L: LayoutV5, device: str = "cpu", cpu_workers: int = 1,
             gpu_workers: int = 1, can_gpu: bool = True, log=print, label: str = "flights") -> int:
    n = run_cpu(cfg, jobs, L, cpu_workers, log, label)
    mj = [j for j in jobs if j["model"]["kind"] not in GUIDANCE_KINDS]
    if mj and can_gpu:
        n += run_mppi(cfg, mj, L, device, gpu_workers, log, label)
    elif mj:
        log(f"[{label}] {len(mj) - len(done_ids(L, mj))} MPPI flights wait for a GPU runtime")
    return n


# ============================================================================
# Flights table (incremental: only new or changed result files are parsed)
# ============================================================================
def collect_flights(L: LayoutV5, write: bool = True):
    import pandas as pd
    idx_p = L.cache / "flights_index.pkl"
    index: Dict[str, Tuple[float, dict]] = {}
    if idx_p.exists():
        try:
            index = pickle.loads(idx_p.read_bytes())
        except Exception:
            index = {}
    seen, changed = set(), 0
    root = L.flights / "results"
    stack = [root] if root.exists() else []
    while stack:
        d = stack.pop()
        with os.scandir(d) as it:
            for e in it:
                if e.is_dir():
                    stack.append(Path(e.path))
                elif e.name.endswith(".json"):
                    rel = os.path.relpath(e.path, root)
                    seen.add(rel)
                    m = e.stat().st_mtime
                    if rel not in index or index[rel][0] != m:
                        try:
                            index[rel] = (m, json.loads(Path(e.path).read_text()))
                            changed += 1
                        except (OSError, ValueError):
                            index.pop(rel, None)
    gone = set(index) - seen
    for k in gone:
        del index[k]
    if changed or gone or not idx_p.exists():
        L.atomic_bytes(idx_p, pickle.dumps(index))
    df = pd.DataFrame([v[1] for v in index.values()])
    if write and len(df) and (changed or gone or not (L.flights / "flights.csv").exists()):
        L.write_csv(df, L.flights / "flights.csv")
    return df


# ============================================================================
# Gates (section 4.0)
# ============================================================================
def _gate_done(L: LayoutV5, gate: str) -> Optional[dict]:
    p = L.checks / f"{gate}.done"
    j = L.checks / f"{gate}.json"
    if p.exists() and j.exists():
        return json.loads(j.read_text())
    return None


def _gate_write(L: LayoutV5, gate: str, doc: dict, extra: Sequence[Path] = ()) -> dict:
    doc = {"gate": gate, **doc, "git_commit": _GIT.get(), "time": time.strftime("%Y-%m-%d %H:%M:%S")}
    L.write_json(L.checks / f"{gate}.json", doc)
    L.mark_done(L.checks / f"{gate}.done", [L.checks / f"{gate}.json", *extra])
    return doc


def _results(L: LayoutV5, jobs: Sequence[dict]) -> List[dict]:
    return [r for r in (_read_json(L.result_path(F5.flight_id_v5(j))) for j in jobs) if r is not None]


def gate_s0(cfg: StudyConfigV5, L: LayoutV5, data: Optional[dict], log=print) -> dict:
    """Benchmark: checksums, targets inside the footprint, ORB-1 feasibility, flight time < horizon, feature ranges."""
    have = _gate_done(L, "S0")
    if have is not None:
        return have
    import pandas as pd
    out = {}
    try:
        B.verify()
        out["checksums"] = True
    except B.ChecksumError as e:
        out["checksums"] = False
        out["checksum_error"] = str(e)
    spec = json.loads((B.BENCH_DIR / "ORB1_spec.json").read_text())
    fp = spec["footprint"]
    d90, d0 = fp["sigma90"]["downrange_km"], fp["sigma0"]["downrange_km"]
    out["footprint_km"] = [d90, d0]
    out["targets_km"] = [t["downrange_km"] for t in spec["targets"]]
    out["targets_inside"] = all(d90 < t < d0 for t in out["targets_km"])
    orb = pd.read_csv(B.BENCH_DIR / "ORB1.csv")
    out["feasible_frac"] = float(orb["feasible"].mean())
    out["feasible_4p5g_frac"] = float(orb["feasible_4p5g"].mean())
    out["longest_nominal_s"] = float(max(fp["sigma90"]["t_s"], fp["sigma0"]["t_s"]))
    out["horizon_s"] = cfg.s0_horizon_s
    ok = out["checksums"] and out["targets_inside"] and out["feasible_frac"] >= cfg.s0_feasible_min and \
        out["longest_nominal_s"] < cfg.s0_horizon_s
    if data is not None:
        ic = B.ORB1Spec().nominal_ic()
        cols = {"h_m": 0, "V_mps": 3, "gamma": 4, "a_D": 8, "a_L": 9}
        Y = np.concatenate([e["y"] for e in data["train"]])
        lo, hi = Y.min(0), Y.max(0)
        rows = []
        for sig in (90.0, 0.0):
            y = B.fly_constant_bank_sim(ic, sig, record=True)["y"]
            for name, k in cols.items():
                outside = float(((y[:, k] < lo[k]) | (y[:, k] > hi[k])).mean())
                rows.append({"sigma_deg": sig, "feature": name, "min": float(y[:, k].min()), "max": float(y[:, k].max()),
                             "train_min": float(lo[k]), "train_max": float(hi[k]), "frac_outside": outside})
        L.write_csv(pd.DataFrame(rows), L.checks / "S0_feature_ranges.csv")
        out["feature_ranges_inside"] = all(r["frac_outside"] == 0.0 for r in rows)
        ok = ok and out["feature_ranges_inside"]
    else:
        out["feature_ranges_inside"] = None
    out["pass"] = bool(ok)
    copy_benchmark(L)
    if data is None:                                       # not final until the feature ranges are checked
        L.write_json(L.checks / "S0.json", {"gate": "S0", **out})
        return out
    return _gate_write(L, "S0", out, [L.checks / "S0_feature_ranges.csv"])


def gate_jobs(cfg: StudyConfigV5, L: LayoutV5) -> Dict[str, List[dict]]:
    """The flights of V1-V4 (V2's depend on V1's realized peaks)."""
    conds = {(int(c["case"]), c["atmosphere"]): c for c in conditions(cfg, "A256")}
    cases = [k for k in cfg.gate_cases if (k, "nominal") in conds]
    fn = dict(cfg.fnpeg)                                   # the defaults: gates run before tuning
    out = {"V1": [guidance_job(conds[(k, "nominal")], "gate_V1_fnpeg", "gates", "fnpeg",
                               dict(fn, oracle=True, filters=False, load_relief=False)) for k in cases]}
    v1 = {r["case"]: r for r in _results(L, out["V1"])}
    if len(v1) == len(cases):
        heavy = [k for k in cases if v1[k]["peak_g"] > cfg.v2_trigger_g]
        out["V2"] = [guidance_job(conds[(k, "nominal")], "gate_V2_fnpeg_3g", "gates", "fnpeg",
                                  dict(fn, oracle=True, filters=False, amax_g=cfg.v2_amax_g)) for k in heavy]
        for g in cfg.v2_sweep_g:
            out["V2"] += [guidance_job(conds[(k, "nominal")], f"gate_V2_fnpeg_{g:g}g", "gates", "fnpeg",
                                       dict(fn, oracle=True, filters=False, amax_g=g))
                          for k in heavy[:cfg.v2_sweep_n]]
    out["V3"] = [guidance_job(conds[(k, a)], "gate_V3_fnpeg", "gates", "fnpeg", fn)
                 for a in ("mismatch", "mismatch_high") for k in cases[:cfg.v3_n_cases] if (k, a) in conds]
    out["V4"] = [guidance_job(conds[(k, "nominal")], "gate_V4_predguid", "gates", "predguid", dict(cfg.predguid))
                 for k in cases]
    return out


def gates(cfg: StudyConfigV5, L: LayoutV5, L3: S.Layout, workers: int = 1, log=print) -> dict:
    """Fly and evaluate V1-V4 (CPU), then E0. Each passed or failed gate is cached (01_checks/<gate>.done)."""
    out = {}
    for _ in range(2):                                     # V2's jobs need V1's flights: two passes
        js = gate_jobs(cfg, L)
        run_cpu(cfg, [j for v in js.values() for j in v], L, workers, log, "gates")
    js = gate_jobs(cfg, L)
    out["V1"] = _gate_v1(cfg, L, js["V1"])
    out["V2"] = _gate_v2(cfg, L, js.get("V2", []))
    out["V3"] = _gate_v3(cfg, L, js["V3"])
    out["V4"] = _gate_v4(cfg, L, js["V4"])
    out["E0"] = gate_e0(cfg, L, L3, log)
    return out


def _gate_v1(cfg, L, jobs) -> dict:
    have = _gate_done(L, "V1")
    if have is not None:
        return have
    rs = _results(L, jobs)
    if len(rs) < len(jobs):
        return {"pass": None}
    miss = np.array([r["miss_km"] for r in rs])
    land = np.array([r["outcome"] == "drogue" for r in rs])
    n_ok = int((land & (miss < cfg.v1_miss_km)).sum())
    it_frac = [r.get("fnpeg_iter_le6_frac", np.nan) for r in rs]
    tel = [S.load_telemetry(L, r["flight_id"]) for r in rs]
    its = np.concatenate([t["dec_n_iter"][t["dec_n_iter"] >= 0] for t in tel]) if tel else np.array([])
    frac_le = float((its <= cfg.v1_iter_max).mean()) if len(its) else float("nan")
    ok = n_ok >= math.ceil(cfg.v1_min_frac * len(rs) - 1e-9) and frac_le >= cfg.v1_iter_frac
    return _gate_write(L, "V1", {"pass": bool(ok), "n": len(rs), "n_miss_lt_2km": n_ok, "miss_km": miss.tolist(),
                                 "cycles": int(len(its)), "frac_cycles_iter_le_6": frac_le,
                                 "per_flight_iter_le6": it_frac})


def _gate_v2(cfg, L, jobs) -> dict:
    have = _gate_done(L, "V2")
    if have is not None:
        return have
    if not (L.checks / "V1.done").exists():
        return {"pass": None}
    rs = _results(L, jobs)
    if len(rs) < len(jobs):
        return {"pass": None}
    main = [r for r in rs if r["arm"] == "gate_V2_fnpeg_3g"]
    peak_ok = all(r["peak_g"] <= cfg.v2_amax_g + cfg.v2_tol_g for r in main)
    miss_ok = [r["outcome"] == "drogue" and r["miss_km"] < cfg.v1_miss_km for r in main]
    sweep = {}
    for r in rs:
        sweep.setdefault(int(r["case"]), {})[r["arm"].replace("gate_V2_fnpeg_", "")] = (r["peak_g"], r["miss_km"])
    # "where feasible": the 3 g miss criterion is required where the 3 g flight lands (an infeasible 3 g case
    # cannot reach the target at that load; it is reported, as Lu's Fig. 17 sweep is)
    ok = peak_ok and all(m for r, m in zip(main, miss_ok) if r["outcome"] == "drogue")
    return _gate_write(L, "V2", {"pass": bool(ok), "n_cases": len(main), "peak_ok": bool(peak_ok),
                                 "miss_lt_2km": miss_ok, "flights": {r["flight_id"]: {"peak_g": r["peak_g"],
                                                                                    "miss_km": r["miss_km"]} for r in rs},
                                 "sweep": sweep})


def _gate_v3(cfg, L, jobs) -> dict:
    have = _gate_done(L, "V3")
    if have is not None:
        return have
    rs = _results(L, jobs)
    if len(rs) < len(jobs):
        return {"pass": None}
    rows = []
    for r in rs:
        t = S.load_telemetry(L, r["flight_id"])
        ds, ms = r.get("dom_density_scale", 1.0), r.get("dom_mass_scale", 1.0)
        true = ds / ms
        tt, a, kd = t["dec_t_s"], t["dec_a_g"], t["dec_K_D"]
        on = np.nonzero(a >= 0.05)[0]
        if not len(on):
            rows.append({"flight": r["flight_id"], "ok": False, "err": float("nan")})
            continue
        t_chk = tt[on[0]] + cfg.v3_after_s
        k = np.nonzero(tt >= t_chk)[0]
        err = float(abs(kd[k[0]] / true - 1.0)) if len(k) else float("nan")
        rows.append({"flight": r["flight_id"], "true_K_D": true, "K_D": float(kd[k[0]]) if len(k) else float("nan"),
                     "rel_err": err, "ok": bool(np.isfinite(err) and err <= cfg.v3_tol)})
    ok = all(x["ok"] for x in rows)
    return _gate_write(L, "V3", {"pass": bool(ok), "n": len(rows), "flights": rows})


def _gate_v4(cfg, L, jobs) -> dict:
    have = _gate_done(L, "V4")
    if have is not None:
        return have
    rs = _results(L, jobs)
    if len(rs) < len(jobs):
        return {"pass": None}
    n_ok = int(sum(r["outcome"] == "drogue" and r["miss_km"] < cfg.v4_miss_km for r in rs))
    seqs = {r["flight_id"]: r.get("predguid_phases") for r in rs}
    ok = n_ok >= math.ceil(cfg.v4_min_frac * len(rs) - 1e-9)
    return _gate_write(L, "V4", {"pass": bool(ok), "n": len(rs), "n_miss_lt_3p5km": n_ok,
                                 "miss_km": [r["miss_km"] for r in rs], "phase_sequences": seqs})


def gate_e0(cfg: StudyConfigV5, L: LayoutV5, L3: S.Layout, log=print) -> dict:
    """The v5 runner equals run_lockstep_v4 bit for bit on the CPU (v2 member 0, physics), and mppi_filtered with
    frozen filters equals mppi_nominal; 3 flights each, never written to the study's flights."""
    have = _gate_done(L, "E0")
    if have is not None:
        return have
    conds = {(int(c["case"]), c["atmosphere"]): c for c in conditions(cfg, "A256")}
    mppi = dict(cfg.mppi(), n_samples=cfg.e0_samples)
    rows = []
    for kind, model in (("v2", v2single_model(L3, cfg.seeds[0])),
                        ("physics_nominal", {"kind": "physics_nominal", "ckpt": None, "members": None})):
        jobs = [mppi_job(cfg, conds[k], f"E0_{kind}", G, cfg.g_primary, cfg.seeds[0], model, "gates")
                for k in cfg.e0_conditions if k in conds]
        for j in jobs:
            j["mppi"]["n_samples"] = cfg.e0_samples
        ad = F5.load_adapter_v5(model["kind"], model.get("ckpt"), "cpu", model.get("members"), cfg.v2_chunk)
        a = F4.run_lockstep_v4(jobs, ad, mppi, verbose=False)
        ad = F5.load_adapter_v5(model["kind"], model.get("ckpt"), "cpu", model.get("members"), cfg.v2_chunk)
        b = F5.run_mppi_v5(jobs, ad, mppi, verbose=False)
        rows += [_same(x, y, f"v5 runner == run_lockstep_v4 ({kind})") for x, y in zip(a, b)]
        if kind == "physics_nominal":
            ad = F5.load_adapter_v5("physics_filtered_frozen", None, "cpu")
            c = F5.run_mppi_v5(jobs, ad, mppi, verbose=False)
            rows += [_same(x, y, "mppi_filtered (filters frozen at 1) == mppi_nominal") for x, y in zip(b, c)]
    import pandas as pd
    L.write_csv(pd.DataFrame(rows), L.checks / "E0.csv")
    ok = all(r["identical"] for r in rows)
    log(f"[E0] {sum(r['identical'] for r in rows)}/{len(rows)} identical")
    return _gate_write(L, "E0", {"pass": bool(ok), "n": len(rows)}, [L.checks / "E0.csv"])


_SKIP = ("flight_id", "v4_flight_id", "wall_s", "arm", "source", "stage", "plan_ms_mean", "plan_ms_max",
         "plan_ms_median", "plan_ms_p99", "suite", "block", "target_index", "train_seed", "feasible",
         "in_training_domain", "model_kind", "sigma_max_deg", "controller", "git_commit", "variant")


def _same(a: dict, b: dict, check: str) -> dict:
    ra, rb = a["result"], b["result"]
    keys = [k for k in ra if k not in _SKIP and k in rb]
    diff = [k for k in keys if not (ra[k] == rb[k] or (isinstance(ra[k], float) and isinstance(rb[k], float)
                                                     and math.isnan(ra[k]) and math.isnan(rb[k])))]
    ta, tb = a["telemetry"], b["telemetry"]
    tdiff = [k for k in ta if k in tb and not np.array_equal(ta[k], tb[k], equal_nan=True)]
    return {"check": check, "case": ra.get("case"), "atmosphere": ra.get("atmosphere"), "miss_km": ra.get("miss_km"),
            "identical": not diff and not tdiff, "diff_fields": ";".join(diff + tdiff)}


# ============================================================================
# Tuning (section 4.5, DEV48 only)
# ============================================================================
FNPEG_SPACE = {"sigma_f_deg": ("none_or_uniform", 20.0, 80.0, 0.25), "deadband_c0_deg": ("uniform", 0.5, 3.0),
               "deadband_c1_deg_per_mps": ("uniform", 2e-4, 1e-3), "delta_nd": ("uniform", 0.01, 0.04),
               "k0_nd": ("loguniform", 50.0, 200.0), "filter_beta": ("uniform", 0.8, 0.97),
               "dt_s": ("choice", [1.0, 2.0])}
PREDGUID_SPACE = {"LATBIAS": ("uniform", 0.2, 0.8), "KLAT": ("loguniform", 0.006, 0.025), "KA": ("uniform", 45.0, 85.0),
                  "Q7F": ("uniform", 3.0, 9.0), "Q19": ("uniform", 0.1, 0.4), "overcontrol": ("uniform", 2.0, 6.0),
                  "K_RHO_FILTER_GAIN": ("uniform", 0.02, 0.10), "gmax_g": ("uniform", 4.2, 4.5)}
RL_SPACES = {
    "ppo_hist": {"lr": ("loguniform", 1e-4, 1e-3), "n_steps": ("choice", [32, 64, 128]),
                 "ent_coef": ("choice", [0.0, 1e-3, 1e-2]), "gae_lambda": ("uniform", 0.9, 0.98),
                 "clip": ("choice", [0.1, 0.2, 0.3])},
    "bcppo_fnpeg": {"lr": ("loguniform", 1e-5, 1e-4), "kl_coef": ("loguniform", 0.01, 1.0),
                    "critic_warmup_updates": ("choice", [10, 20, 50])},
    "dagger_fnpeg": {"lr": ("loguniform", 1e-4, 3e-3), "dagger_epochs": ("choice", [20, 30, 50]),
                     "hidden": ("choice", [[256, 256], [512, 512]])},
    "sac_hist": {"lr": ("loguniform", 1e-4, 1e-3), "tau": ("choice", [0.005, 0.01]),
                 "batch": ("choice", [512, 1024, 2048])},
}


def sample_configs(space: dict, n: int, seed: int) -> List[dict]:
    """Config 0 = {} (the defaults); 1..n-1 drawn from the documented ranges with a fixed seed."""
    rng = np.random.default_rng(seed)
    out = [{}]
    for _ in range(max(0, n - 1)):
        c = {}
        for k, spec in space.items():
            kind = spec[0]
            if kind == "uniform":
                c[k] = float(rng.uniform(spec[1], spec[2]))
            elif kind == "loguniform":
                c[k] = float(math.exp(rng.uniform(math.log(spec[1]), math.log(spec[2]))))
            elif kind == "choice":
                c[k] = spec[1][int(rng.integers(len(spec[1])))]
            elif kind == "none_or_uniform":
                c[k] = None if rng.uniform() < spec[3] else float(rng.uniform(spec[1], spec[2]))
        out.append(c)
    return out


def _score_dev(rs: List[dict]) -> dict:
    land = np.array([r["outcome"] == "drogue" for r in rs])
    miss = np.array([float(r["miss_km"]) for r in rs])
    return {"n": len(rs), "success_5km": float((land & (miss < 5.0)).mean()) if len(rs) else float("nan"),
            "mean_log_miss": float(np.mean(np.log(np.maximum(miss, 0.01)))) if len(rs) else float("nan"),
            "mean_peak_g": float(np.mean([r["peak_g"] for r in rs])) if len(rs) else float("nan"),
            "median_miss_km": float(np.median(miss)) if len(rs) else float("nan")}


def _choose(rows: List[dict], peak_max: float) -> dict:
    """success@5km, then mean log miss, subject to mean peak g <= peak_max (else the lowest peak g)."""
    ok = [r for r in rows if r["score"]["mean_peak_g"] <= peak_max]
    if not ok:
        return min(rows, key=lambda r: r["score"]["mean_peak_g"])
    return max(ok, key=lambda r: (r["score"]["success_5km"], -r["score"]["mean_log_miss"]))


def tune_vf(cfg: StudyConfigV5, L: LayoutV5, workers: int, log=print) -> Optional[float]:
    """Vf = median speed at 7.6 km of the DEV48 classical flights (Lu's CEV: 150 m/s)."""
    p = L.tuning / "fnpeg" / "vf.json"
    if p.exists():
        return float(json.loads(p.read_text())["vf_mps"])
    jobs = [classical_job(c, "tuning") for c in conditions(cfg, "DEV48")]
    run_cpu(cfg, jobs, L, workers, log, "tuning classical")
    rs = _results(L, jobs)
    if len(rs) < len(jobs):
        return None
    v = [float(S.load_telemetry(L, r["flight_id"])["V_mps"][-1]) for r in rs if r["outcome"] == "drogue"]
    vf = float(round(np.median(v), 1)) if v else 150.0
    L.write_json(p, {"vf_mps": vf, "n_flights": len(v), "rule": "median speed at the drogue altitude over the "
                     "DEV48 classical flights that landed"})
    log(f"[tuning] Vf = {vf} m/s ({len(v)} classical DEV48 flights)")
    return vf


def tune_guidance(cfg: StudyConfigV5, L: LayoutV5, workers: int, log=print) -> Dict[str, dict]:
    """Random search for FNPEG and PredGuid on DEV48 (cached per configuration), then the choice."""
    vf = tune_vf(cfg, L, workers, log)
    if vf is None:
        return {}
    out = {}
    spaces = {"fnpeg": (FNPEG_SPACE, dict(cfg.fnpeg, vf_mps=vf)), "predguid": (PREDGUID_SPACE, dict(cfg.predguid))}
    dev = conditions(cfg, "DEV48")
    jobs, plan = [], {}
    for m, (space, base) in spaces.items():
        if (L.tuning / m / "chosen.json").exists():
            continue
        for j, c in enumerate(sample_configs(space, cfg.tune_n, cfg.tune_seed + (0 if m == "fnpeg" else 1))):
            if (L.tuning / m / f"config_{j:02d}.json").exists():
                continue
            full = dict(base, **c)
            js = [guidance_job(cd, f"tune_{m}_c{j:02d}", "tuning", m, full) for cd in dev]
            plan[(m, j)] = (full, c, js)
            jobs += js
    run_cpu(cfg, jobs, L, workers, log, "tuning")
    for (m, j), (full, c, js) in plan.items():
        rs = _results(L, js)
        if len(rs) == len(js):
            L.write_json(L.tuning / m / f"config_{j:02d}.json", {"method": m, "config_index": j, "config": full,
                                                                "sampled": c, "score": _score_dev(rs)})
    for m in spaces:
        out[m] = _choose_method(cfg, L, m, cfg.tune_n, log)
    return out


def _choose_method(cfg, L, m: str, n: int, log=print) -> Optional[dict]:
    p = L.tuning / m / "chosen.json"
    if p.exists():
        return json.loads(p.read_text())
    rows = [json.loads((L.tuning / m / f"config_{j:02d}.json").read_text()) for j in range(n)
            if (L.tuning / m / f"config_{j:02d}.json").exists()]
    if len(rows) < n:
        return None
    best = _choose(rows, cfg.tune_peak_g_max)
    doc = {"method": m, "config": best["config"], "config_index": best["config_index"], "score": best["score"],
           "n_configs": len(rows), "objective": "success@5km, then mean log miss, subject to mean peak g <= "
                                                f"{cfg.tune_peak_g_max} g (DEV48)"}
    L.write_json(p, doc)
    import pandas as pd
    L.write_csv(pd.DataFrame([{"config_index": r["config_index"], **r["score"],
                               **{f"cfg_{k}": v for k, v in (r.get("sampled") or {}).items()}} for r in rows]),
                L.tuning / m / "configs.csv")
    log(f"[tuning] {m}: config {best['config_index']} chosen ({best['score']})")
    return doc


def rl_cfg(cfg: StudyConfigV5, L: LayoutV5, arm: str, extra: Optional[dict] = None):
    from jepa import rl_v5 as RL
    base_arm = "ppo_hist" if arm == "ppo_matched" else arm        # ppo_matched reuses ppo_hist's tuned settings
    t = tuned(L, f"rl_{base_arm}", cfg) if base_arm in RL_SPACES else {}
    kw = {**cfg.rl, **cfg.rl_arms.get(arm, {}), **t, **(extra or {})}
    if "hidden" in kw:
        kw["hidden"] = tuple(kw["hidden"])
    return RL.RLCfg(**kw)


def tune_rl(cfg: StudyConfigV5, L: LayoutV5, arms: Sequence[str], device: str, workers: int, log=print) -> dict:
    """RL configurations (seed 0, cfg.rl_tune_frac of the budget) scored on DEV48 in the full simulator."""
    from jepa import rl_v5 as RL
    if not is_tuned(L, "fnpeg"):
        log("[tuning] RL tuning waits for FNPEG's tuning (its expert and clone source)")
        return {}
    fn = tuned(L, "fnpeg", cfg)
    out = {}
    for arm in arms:
        if arm not in RL_SPACES:
            continue
        m = f"rl_{arm}"
        n = int(cfg.rl_tune_n.get(arm, 1))
        if (L.tuning / m / "chosen.json").exists():
            out[arm] = json.loads((L.tuning / m / "chosen.json").read_text())
            continue
        base = asdict(RL.RLCfg(**{**cfg.rl, **cfg.rl_arms.get(arm, {})}))
        for j, c in enumerate(sample_configs(RL_SPACES[arm], n, cfg.tune_seed + zlib.crc32(arm.encode()) % 1000)):
            if (L.tuning / m / f"config_{j:02d}.json").exists():
                continue
            f = cfg.rl_tune_frac
            red = {"budget_decisions": max(1, int(base["budget_decisions"] * f)), "budget_hours": base["budget_hours"] * f,
                   "dagger_rollouts": max(2, int(base["dagger_rollouts"] * f)),
                   "bc_episodes": max(2, int(base["bc_episodes"] * f)), "curve_points": 0}
            kw = {**base, **c, **red, "seed": 0}
            kw["hidden"] = tuple(kw["hidden"])
            d = L.tuning / m / f"config_{j:02d}"
            RL.train_rl_dir(L, d, arm, RL.RLCfg(**kw), device, fn, log, label=f"tune {arm} config {j}")
            js = [policy_job(cd, f"tune_{arm}_c{j:02d}", "tuning", str(d / "policy.pt"), 0)
                  for cd in conditions(cfg, "DEV48")]
            run_cpu(cfg, js, L, workers, log, f"tuning {arm}")
            rs = _results(L, js)
            if len(rs) == len(js):
                L.write_json(L.tuning / m / f"config_{j:02d}.json", {"method": m, "config_index": j,
                                                                    "config": {**c}, "sampled": c,
                                                                    "budget": red, "score": _score_dev(rs)})
        out[arm] = _choose_method(cfg, L, m, n, log)
    return out


# ============================================================================
# Training stages
# ============================================================================
def load_data_v5(cfg: StudyConfigV5, data_path: str, device) -> dict:
    """The 6,000 flights, the v3 split, v2 stats, windows for R0's receptive field and the E3 curriculum."""
    from jepa import models_v2 as V
    from jepa import train_v2 as T2
    eps = T2.load_episodes(data_path)
    if cfg.n_episodes:
        eps = eps[:cfg.n_episodes]
    rf = max(M5.RF, V.WorldModelV2(n_members=1, **cfg.r0_arch).rf)
    max_T = max(T for _, T in M5.E3TrainCfg(**cfg.e3_train).curriculum)
    data = T2.prepare(eps, max_T, device, rf, split_seed=cfg.split_seed)
    data["n_episodes"] = len(eps)
    return data


def train_e3_stage(cfg: StudyConfigV5, L: LayoutV5, data_fn: Callable[[], dict], device, log=print) -> List[dict]:
    """Every E3 arm: learning-rate pilots on seed 0, then seeds (seed-major within the arm)."""
    out = []
    todo = [(a, s) for a in cfg.e3_arms for s in cfg.seeds if not L.ckpt(a, s).exists()]
    log(f"[E3] {len(cfg.e3_arms) * len(cfg.seeds) - len(todo)}/{len(cfg.e3_arms) * len(cfg.seeds)} models trained")
    for arm in cfg.e3_arms:
        if all(L.ckpt(arm, s).exists() for s in cfg.seeds):
            continue
        data = data_fn()
        lr = M5.choose_lr(cfg, data, L, arm, device, log)
        for s in cfg.seeds:
            t0 = time.time()
            info = M5.train_arm_seed(cfg, data, L, arm, s, device, lr, log)
            log(f"== {arm} seed {s}: {info.get('params')} parameters, {(time.time() - t0) / 60:.1f} min")
            out.append(info)
    return out


def train_rl_stage(cfg: StudyConfigV5, L: LayoutV5, arms: Sequence[str], device, log=print) -> List[dict]:
    """Every RL arm and seed (seed-major within the arm), then its fast-env evaluation (fast-vs-full gap)."""
    from jepa import rl_v5 as RL
    fn = tuned(L, "fnpeg", cfg)
    out = []
    for arm in arms:
        base = "ppo_hist" if arm == "ppo_matched" else arm
        if not is_tuned(L, "fnpeg") or (base in RL_SPACES and not is_tuned(L, f"rl_{base}")):
            log(f"[E4] {arm} waits for tuning (03_tuning/)")
            continue
        for s in cfg.rl_seeds:
            info = RL.train_rl_arm(L, arm, s, rl_cfg(cfg, L, arm), device, fn, log)
            d = L.policy_dir(arm, s)
            fe = d / "fast_eval.json"
            if (d / "policy.pt").exists() and not fe.exists():
                pol = RL.load_policy(d / "policy.pt", device)
                L.write_json(fe, {"train_distribution": RL.fast_success(pol, cfg.fast_eval_episodes, seed=12345,
                                                                         device=device),
                                  "A256": RL.fast_success_on(pol, conditions(cfg, "A256"), device),
                                  "A256_cases": [[int(c["case"]), c["atmosphere"]] for c in conditions(cfg, "A256")]})
            out.append(info)
    return out


def open_loop_stage(cfg: StudyConfigV5, L: LayoutV5, data_fn, device, log=print) -> int:
    todo = [(a, s) for a in cfg.e3_arms for s in cfg.seeds
            if L.ckpt(a, s).exists() and not (L.open_loop / f"footprint_{a}_seed{s}.csv").exists()]
    log(f"[open loop] {len(todo)} model-seeds to evaluate")
    for a, s in todo:
        t0 = time.time()
        M5.open_loop_seed(cfg, L, data_fn(), a, s, device)
        log(f"  {a} seed {s}: {time.time() - t0:.0f}s")
    return len(todo)


# ============================================================================
# Legacy RL (section 4.6)
# ============================================================================
def legacy_inventory(cfg: StudyConfigV5, L: LayoutV5, legacy_root: Optional[str], log=print) -> List[dict]:
    """00_study/legacy_rl_inventory.csv and the final models to evaluate ([] if the folder is not found)."""
    import pandas as pd
    from jepa import rl_v5 as RL
    rows = RL.legacy_inventory(Path(legacy_root)) if legacy_root else [
        {"path": None, "run_name": None, "steps": None, "size_bytes": None, "date": None, "kind": "not found"}]
    L.write_csv(pd.DataFrame(rows), L.study / "legacy_rl_inventory.csv")
    finals = [r for r in rows if r["kind"] == "final model"]
    log(f"[legacy RL] {len(rows)} files listed, {len(finals)} final models")
    return finals


def legacy_jobs(cfg: StudyConfigV5, L: LayoutV5, finals: Sequence[dict]) -> List[dict]:
    import re
    jobs = []
    for r in finals:
        name = "rl_legacy_" + re.sub(r"[^A-Za-z0-9_]+", "_", Path(r["path"]).stem)
        kw = _legacy_env_kwargs(r["path"])
        for suite in ("HELD256", "A256"):
            jobs += [policy_job(c, name, "legacy", r["path"], None, legacy=True, env_kwargs=kw)
                     for c in conditions(cfg, suite)]
    return jobs


def _legacy_env_kwargs(path: str) -> dict:
    """rl_env's interval observation if the policy's input is wider than the 13 basic features."""
    try:
        from jepa.rl_v5 import LegacyPolicy
        return {"interval_obs": True} if LegacyPolicy(path).obs_dim > 13 else {}
    except Exception:
        return {}


def legacy_finals(L: LayoutV5) -> List[dict]:
    """The final legacy models listed by the inventory ([] before the legacy stage ran)."""
    p = L.study / "legacy_rl_inventory.csv"
    if not p.exists():
        return []
    import pandas as pd
    d = pd.read_csv(p)
    return d[d.kind == "final model"].to_dict("records") if "kind" in d else []


def best_arms(cfg: StudyConfigV5, L: LayoutV5, L3: S.Layout) -> List[str]:
    """Tier 2: the non-JEPA (E3) arm and the RL arm with the highest A256 success@5km (once all have flown)."""
    out = []
    for group, kind in ((list(cfg.e3_arms), "E3"), (list(RL_TIER1), "E4")):
        best, best_s = None, -1.0
        for arm in group:
            js = [j for j in (jobs_e3(cfg, L, L3) if kind == "E3" else jobs_e4(cfg, L, L3, [arm]))
                  if j["controller"] == arm and j["suite"] == "A256"]
            rs = _results(L, js)
            if not js or len(rs) < len(js):
                return []                                   # not every arm has flown A256 yet
            s = float(np.mean([r["outcome"] == "drogue" and r["miss_km"] < 5.0 for r in rs]))
            if s > best_s:
                best, best_s = arm, s
        if best:
            out.append(best)
    return out


def all_stage_jobs(cfg: StudyConfigV5, L: LayoutV5, L3: S.Layout, log=lambda *a: None) -> Dict[str, List[dict]]:
    """Every flight stage of this tier (status and the notebook's flight cells)."""
    t2 = cfg.tier >= 2
    rl = list(RL_TIER1) + (list(RL_TIER2) if t2 else [])
    J = {"A256": jobs_a256(cfg, L, L3), "ORB1_block1": orb1_jobs(cfg, L, L3, ORB1_ARMS, [0], "ORB1", log),
         "E3": jobs_e3(cfg, L, L3), "E4": jobs_e4(cfg, L, L3, rl)}
    fin = legacy_finals(L)
    if fin:
        J["legacy"] = legacy_jobs(cfg, L, fin)
    if t2:
        J["ORB1_more"] = orb1_jobs(cfg, L, L3, ORB1_ARMS, [1, 2], "ORB1", log)
        J["ORB1_extra"] = orb1_jobs(cfg, L, L3, ["mppi_truth"] + best_arms(cfg, L, L3), [0], "ORB1", log)
        J["explore"] = jobs_explore(cfg, L, L3)
    return J


# ============================================================================
# Status
# ============================================================================
def _times(L: LayoutV5, jobs: Sequence[dict]) -> List[float]:
    by_dir: Dict[Path, set] = {}
    for j in jobs:
        p = L.result_path(F5.flight_id_v5(j))
        by_dir.setdefault(p.parent, set()).add(p.name)
    out = []
    for d, names in by_dir.items():
        try:
            with os.scandir(d) as it:
                out += [e.stat().st_mtime for e in it if e.name in names]
        except FileNotFoundError:
            pass
    return sorted(out)


def _arm_type(j: dict) -> str:
    k = j["model"]["kind"]
    return {"physics_filtered": "mppi_physics", "physics_truth": "mppi_physics", "e3": "mppi_e3",
            "v2": "mppi_v2"}.get(k, k)


def status(cfg: StudyConfigV5, L: LayoutV5, stage_jobs: Dict[str, List[dict]], runtime: str = "cpu",
           write: bool = True) -> dict:
    """Done / total per stage and arm; flights per hour per arm type from the first 20 flights; hours and units left."""
    rows, rate_by_type = [], {}
    for stage, js in stage_jobs.items():
        arms: Dict[tuple, List[dict]] = {}
        for j in js:
            arms.setdefault((j["suite"], j["controller"], j.get("load")), []).append(j)
        for (suite, arm, load), g in arms.items():
            t = _times(L, g)
            first = t[:20]
            rate = (len(first) - 1) / max(first[-1] - first[0], 1.0) * 3600 if len(first) >= 2 else None
            typ = _arm_type(g[0])
            if rate and (typ not in rate_by_type or len(first) > rate_by_type[typ][1]):
                rate_by_type[typ] = (rate, len(first))
            rows.append({"stage": stage, "suite": suite, "arm": arm, "load": load, "type": typ, "done": len(t),
                         "total": len(g), "flights_per_h": round(rate, 1) if rate else None})
    for r in rows:
        rt = r["flights_per_h"] or (rate_by_type.get(r["type"], (None,))[0])
        r["hours_left"] = round((r["total"] - r["done"]) / rt, 2) if rt else None
    stages = []
    for stage in stage_jobs:
        rr = [r for r in rows if r["stage"] == stage]
        hl = [r["hours_left"] for r in rr if r["hours_left"] is not None]
        stages.append({"stage": stage, "done": sum(r["done"] for r in rr), "total": sum(r["total"] for r in rr),
                       "hours_left": round(sum(hl), 2) if hl else None})
    train = []
    for a in cfg.e3_arms:
        for s in cfg.seeds:
            d = L.model_dir(a, s)
            st = "done" if (d / "ckpt.pt").exists() else ("partial" if (d / "train_state.pt").exists() else "todo")
            train.append({"kind": "E3", "arm": a, "seed": s, "state": st})
    for a in RL_TIER1 + (RL_TIER2 if cfg.tier >= 2 else ()):
        for s in cfg.rl_seeds:
            d = L.policy_dir(a, s)
            st = "done" if (d / "policy.pt").exists() else ("partial" if d.exists() else "todo")
            train.append({"kind": "E4", "arm": a, "seed": s, "state": st})
    hours = sum(s["hours_left"] or 0.0 for s in stages)
    key = next((k for k in cfg.units_per_hour if k.lower() in runtime.lower()), "cpu")
    doc = {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "tier": cfg.tier, "runtime": runtime, "stages": stages,
           "arms": rows, "training": train, "rates_by_type": {k: round(v[0], 1) for k, v in rate_by_type.items()},
           "flights_done": int(sum(r["done"] for r in rows)), "flights_total": int(sum(r["total"] for r in rows)),
           "trainings_done": int(sum(t["state"] == "done" for t in train)), "trainings_total": len(train),
           "flight_hours_left": round(hours, 2), "units_per_hour": cfg.units_per_hour[key],
           "flight_units_left": round(hours * cfg.units_per_hour[key], 1)}
    if write:
        import pandas as pd
        L.write_json(L.study / "status.json", doc)
        hp = L.study / "status_history.csv"
        new = pd.DataFrame([{"time": doc["time"], "runtime": runtime, "flights_done": doc["flights_done"],
                             "flights_total": doc["flights_total"], "trainings_done": doc["trainings_done"],
                             "trainings_total": doc["trainings_total"], "flight_hours_left": doc["flight_hours_left"],
                             **{f"{s['stage']}_done": s["done"] for s in stages},
                             **{f"{s['stage']}_total": s["total"] for s in stages},
                             **{f"{s['stage']}_hours_left": s["hours_left"] for s in stages}}])
        old = pd.read_csv(hp) if hp.exists() else None
        L.write_csv(pd.concat([old, new], ignore_index=True) if old is not None else new, hp)
    return doc


def update_manifest(L: LayoutV5) -> None:
    """00_study/README.md: the folder index (counts with os.scandir, no tree walk)."""
    lines = ["# SpaceJEPA v5 study folder", "", "Auto-written by `jepa/study_v5.py` (do not edit). The plan is "
             "`jepa/v5_plan.md` in the repo.", ""]
    cfg = L.study / "study_config.json"
    if cfg.exists():
        c = json.loads(cfg.read_text())
        lines += [f"* git commit: `{c.get('git_commit')}`", f"* created: {c.get('created')}",
                  f"* mode: {c['config'].get('mode')}, tier {c['config'].get('tier')}", ""]
    st = L.study / "status.json"
    if st.exists():
        s = json.loads(st.read_text())
        lines += [f"Status ({s.get('time')}): {s.get('flights_done')} / {s.get('flights_total')} flights, "
                  f"{s.get('trainings_done')} / {s.get('trainings_total')} trainings.", ""]
    lines += ["| folder | what |", "|---|---|"] + [f"| `{k}/` | {v} |" for k, v in STAGES.items()]
    caps = L.analysis / "captions.json"
    if caps.exists():
        lines += ["", "## Figures and tables", ""]
        for k, v in sorted(json.loads(caps.read_text()).items()):
            lines.append(f"* `{k}`: {v.get('caption', '')}")
    L.atomic_bytes(L.study / "README.md", ("\n".join(lines) + "\n").encode())


# ============================================================================
# Compute per decision (RQ5)
# ============================================================================
def timing_stage(cfg: StudyConfigV5, L: LayoutV5, L3: S.Layout, device: str, log=print) -> Optional[object]:
    """
    01_checks/compute_per_decision.csv: median and p99 ms per guidance decision for each arm type, one flight
    (A256 case 0 or the first case, nominal air), single CPU thread, and on the GPU for the MPPI arms.
    Parameter counts of every learned model. Cached per (arm, device).
    """
    import pandas as pd
    path = L.checks / "compute_per_decision.csv"
    have = pd.read_csv(path) if path.exists() else pd.DataFrame(columns=["arm", "device"])
    done = set(zip(have.arm, have.device))
    cond = [c for c in conditions(cfg, "A256") if c["atmosphere"] == "nominal"][0]
    rows = []
    old_threads = torch.get_num_threads()
    torch.set_num_threads(1)
    try:
        for arm in ("fnpeg", "predguid", "classical"):
            if (arm, "cpu") in done:
                continue
            if arm == "classical":
                ms = _classical_step_ms(classical_job(cond, "timing"))
            else:
                j = arm_job(cfg, L, L3, cond, arm, "timing")
                law = F5.make_law(j["model"]["kind"], j["model"]["cfg"])
                ms = []
                rec = F5.run_guidance_lockstep([j], law, verbose=False)[0]
                ms = [rec["result"]["plan_ms_median"], rec["result"]["plan_ms_p99"]]
                rows.append({"arm": arm, "device": "cpu", "ms_median": ms[0], "ms_p99": ms[1], "params": 0})
                continue
            rows.append({"arm": arm, "device": "cpu", "ms_median": float(np.median(ms)),
                         "ms_p99": float(np.percentile(ms, 99)), "params": 0})
        mppi_arms = [("mppi_filtered", None), ("v2single", None)] + [(a, None) for a in cfg.e3_arms
                                                                     if L.ckpt(a, cfg.seeds[0]).exists()]
        devs = ["cpu"] + (["cuda"] if device == "cuda" else [])
        for arm, _ in mppi_arms:
            for dv in devs:
                if (arm, dv) in done:
                    continue
                j = arm_job(cfg, L, L3, cond, arm, "timing", cfg.seeds[0])
                ms = _mppi_plan_ms(cfg, j, dv)
                rows.append({"arm": arm, "device": dv, "ms_median": float(np.median(ms)),
                             "ms_p99": float(np.percentile(ms, 99)), "params": _params(j["model"])})
        for arm in RL_TIER1 + RL_TIER2:
            p = L.policy_dir(arm, cfg.rl_seeds[0]) / "policy.pt"
            if (arm, "cpu") in done or not p.exists():
                continue
            from jepa import rl_v5 as RL
            pol = RL.load_policy(p)
            o = torch.zeros(1, pol.hparams["obs_dim"])
            ts = []
            for _ in range(200):
                t0 = time.perf_counter()
                with torch.no_grad():
                    pol.deterministic(o)
                ts.append(1000 * (time.perf_counter() - t0))
            rows.append({"arm": arm, "device": "cpu", "ms_median": float(np.median(ts)),
                         "ms_p99": float(np.percentile(ts, 99)), "params": int(sum(x.numel() for x in pol.parameters()))})
    finally:
        torch.set_num_threads(old_threads)
    if rows:
        df = pd.concat([have, pd.DataFrame(rows)], ignore_index=True)
        L.write_csv(df, path)
        log(f"[timing] {len(rows)} arm-device timings added")
        return df
    return have


def _params(model: dict) -> int:
    if model["kind"] == "e3":
        m, _, _ = M5.load_e3(model["ckpt"])
        return M5.n_params(m)
    if model["kind"] == "v2":
        from jepa import models_v2 as V
        ck = V.load_checkpoint_v2(model["ckpt"], "cpu")
        n_wm = M5.n_params(ck["wm"])
        k = len(model["members"]) if model.get("members") else ck["wm"].M
        return int(n_wm * k / ck["wm"].M + M5.n_params(ck["pi"]))
    return 0


def _classical_step_ms(job: dict) -> List[float]:
    """Wall time of each 1 s step of the classical controller (predictor-corrector + simulator)."""
    import rl_env
    from jepa.data import apply_domain
    env = rl_env.OrionEntryEnv(config_path=job["config_path"], controller_mode="baseline", dispersion=False,
                               max_episode_steps=1400, fast_atmosphere=True)
    apply_domain(env.params, job["domain"])
    env.reset(seed=0, options={"ic_override": job["ic"]})
    ms, done = [], False
    while not done:
        t0 = time.perf_counter()
        _, _, term, trunc, _ = env.step(np.array([0.0], dtype=np.float32))
        ms.append(1000 * (time.perf_counter() - t0))
        done = bool(term or trunc)
    return ms[::5]


def _mppi_plan_ms(cfg: StudyConfigV5, job: dict, device: str) -> List[float]:
    """Per-decision planning ms of one MPPI flight: fly_v3.finish is wrapped only here, to keep the list."""
    got = {}
    orig = F.finish

    def keep(job_, env, info, cmds, tel, wall_s, plan_ms, spread_km=()):
        got["ms"] = list(plan_ms)
        return orig(job_, env, info, cmds, tel, wall_s, plan_ms, spread_km)

    F.finish = keep
    try:
        m = job["model"]
        ad = F5.load_adapter_v5(m["kind"], m.get("ckpt"), device, m.get("members"), cfg.v2_chunk)
        F4.run_lockstep_v4([job], ad, cfg.mppi(), verbose=False)
    finally:
        F.finish = orig
    return got.get("ms") or [float("nan")]
