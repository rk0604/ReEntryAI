"""
jepa/study_v4.py

The v4 study's plumbing (jepa/v4_plan.md sections 2-4, 9).

    StudyConfigV4      every number of the plan (smoke() for the CPU smoke run)
    LayoutV4           the v4 Drive folder (00_study ... 06_report, 01_checks).
                       Every v4 write goes through its guarded helpers, which
                       refuse a path outside the v4 folder or inside a read-only
                       input (the v3 folder, the 6,000 flights)
    freeze_config_v4   check the v3 folder's frozen settings against what v4
                       expects, then freeze 00_study/study_config.json
    ARMS, build_jobs   the arm registry and every flight of the study, in run
                       order (E0 -> E1 -> E2 -> E4 -> E3), seed-major per stage
    run_classical_v4, run_flights_v4
                       the resumable queues (a flight whose result exists is
                       skipped; each one is written the moment it lands)
    status             done / total per stage and arm, measured rates, hours
                       left (00_study/status.json); counts files with os.scandir
    e0_checks          the section 4.0 checks into 01_checks/
"""

from __future__ import annotations

import io
import json
import math
import multiprocessing as mp
import os
import sys
import time
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch

SIM_DIR = Path(__file__).resolve().parents[1]
if str(SIM_DIR) not in sys.path:
    sys.path.insert(0, str(SIM_DIR))

from jepa import fly_v3 as F
from jepa import fly_v4 as F4
from jepa import study_v3 as S
from jepa import testset_v3 as TS

ATM = tuple(TS.ATMOSPHERES)
VARIANTS = ("R0", "A1_inputs_v1", "A1b_no_measured_aero", "A2_memory_35s", "A3_blind_prober", "A4_train_T24",
            "A5_sigreg_pooled", "A6_latent_32", "A7a_size_250k", "A7b_size_33k", "A8_actions_v1",
            "B1_v1_plus_inputs", "B2_v1_plus_state_prober", "B3_v1_plus_memory_120s")
STAGE_ORDER = ("E0", "E1", "E2", "E4", "E3")
# v3 settings v4 relies on (v4_plan.md section 2, "Frozen from v3"); the v3 folder must match them.
V3_FROZEN = ("seeds", "split_seed", "n_episodes", "testset_csv", "cases", "atmospheres", "v2_members", "v2_arch",
             "v2_lambda", "v2_wm", "v2_pi", "v1_lambda", "v1_T", "v1_T_pi", "v1_steps", "v1_batch", "v1_pi_steps",
             "v1_pi_batch", "v1_lr", "v1_lr_min", "v1_warmup_frac", "v1_clip", "v1_weight_decay", "v1_extend_tol",
             "v1_extend_frac", "v1_extend_lr_frac", "g_primary", "g_legacy", "legacy_seeds", "mppi", "horizon",
             "w_unc_km", "success_km")
# Batching and worker knobs: change speed and memory, never results; free to differ between sessions.
RUNTIME_KEYS = ("max_rows", "max_rows_physics", "v2_chunk")


# ============================================================================
# Configuration
# ============================================================================
@dataclass
class StudyConfigV4:
    """Every number of v4_plan.md. Frozen to 00_study/study_config.json on first use."""
    mode: str = "full"                                  # "full" | "smoke" (every stage, tiny, numbers meaningless)
    seeds: Tuple[int, ...] = (0, 1, 2, 3, 4)
    cases: Optional[Tuple[int, ...]] = None             # None = all 64
    atmospheres: Tuple[str, ...] = ATM
    testset_csv: str = "configs/v3_testset.csv"
    split_seed: int = 0
    n_episodes: Optional[int] = None
    # ---- loads ----
    g_primary: float = 4.5
    g_extra: Tuple[float, ...] = (4.0, 3.5)             # E2
    g_legacy: float = 8.0                               # E4 (and v3's legacy arm)
    # ---- E0 ----
    e0a_cases: Tuple[int, ...] = (0, 8, 16, 24, 32, 40, 48, 56)
    e0b_seed: int = 0
    e0b_pairs: Tuple[Tuple[str, str], ...] = (("v2single", "v2single_v4run"), ("v1", "v1_v4run"))
    e0b_agree_min: float = 0.98
    e0c_conditions: Tuple[Tuple[int, str], ...] = ((0, "nominal"), (16, "nominal"), (40, "nominal"))
    e0c_samples: int = 32
    fidelity_cases: Tuple[int, ...] = tuple(range(0, 64, 2))
    fidelity_atmospheres: Tuple[str, ...] = ("nominal", "mismatch_high")
    fidelity_ttg_s: Tuple[int, ...] = (300, 200, 100)
    port_check_states: int = 1000
    # ---- E4 ----
    diag_case: int = 10
    diag_atmospheres: Tuple[str, ...] = ("nominal", "mismatch_high")
    # ---- E3: the reference R0 (v2's recipe, one member) and v1's budget for the B variants ----
    variants: Tuple[str, ...] = VARIANTS
    r0_arch: dict = field(default_factory=dict)         # models_v2 defaults (2.45M parameters)
    r0_lambda: float = 0.08
    r0_wm: dict = field(default_factory=lambda: dict(
        steps=24_000, batch=256, warmup=1_000, lr=6e-4, curriculum=((0.0, 24), (0.2, 48), (0.5, 96)), log_every=500))
    r0_pi: dict = field(default_factory=lambda: dict(
        steps=6_000, batch=512, warmup=300, lr=1e-3, weight_decay=1e-4,
        curriculum=((0.0, 16), (0.3, 32), (0.6, 48)), log_every=250))
    a4_curricula: Tuple[tuple, tuple] = (((0.0, 24),), ((0.0, 16),))
    a5_lambda: float = 0.02
    a7a: dict = field(default_factory=lambda: dict(enc_dim=40, pred_dim=60, hidden=40))
    a7b: dict = field(default_factory=lambda: dict(enc_dim=16, pred_dim=18, pred_heads=3, hidden=16))
    v1_budget: dict = field(default_factory=dict)       # {} = v3's v1 settings (the fair budget); smoke overrides
    save_every_s: float = 300.0
    # ---- flights ----
    max_rows: int = 2048
    max_rows_physics: int = 16384                       # physics models are elementwise: bigger batches, fewer launches
    v2_chunk: int = 512
    # ---- analysis ----
    success_km: float = 5.0
    success2_km: float = 2.0
    log_miss_floor_km: float = 0.01
    thresholds_km: Tuple[float, ...] = tuple(np.round(np.arange(0.5, 25.01, 0.5), 2).tolist())
    n_boot: int = 10_000
    equiv_pp: float = 5.0
    r0_check_pp: float = 3.0
    l1_tol_g: float = 0.1
    gap_closed_min: float = 0.5
    ttg_s: Tuple[int, ...] = (400, 300, 200, 100, 50)
    mechanism_ttg_s: int = 300
    separation_km: float = 5.0                          # case 10: "separates" = ground tracks 5 km apart

    @staticmethod
    def smoke() -> "StudyConfigV4":
        """Every stage on a CPU in minutes, on top of a v3 folder built with study_v3.StudyConfig.smoke()."""
        tiny = dict(d_latent=16, enc_dim=32, enc_layers=2, enc_heads=2, enc_band=5, pred_dim=48, pred_layers=2,
                    pred_heads=2, act_chunk=4, mlp_ratio=2)
        v3 = S.StudyConfig.smoke()
        return StudyConfigV4(
            mode="smoke", seeds=v3.seeds, cases=v3.cases, n_episodes=v3.n_episodes, e0a_cases=v3.cases,
            e0c_conditions=((0, "nominal"), (16, "nominal"), (16, "mismatch")), e0c_samples=16,
            fidelity_cases=v3.cases, port_check_states=200, diag_case=0, r0_arch=tiny,
            r0_wm=dict(steps=40, batch=16, warmup=5, curriculum=((0.0, 8), (0.5, 16)), log_every=20),
            r0_pi=dict(steps=20, batch=16, warmup=3, lr=2e-3, curriculum=((0.0, 8), (0.5, 16)), log_every=10),
            a4_curricula=(((0.0, 8),), ((0.0, 8),)),
            a7a=dict(enc_dim=16, pred_dim=24, hidden=32), a7b=dict(enc_dim=8, pred_dim=12, hidden=16),
            v1_budget=dict(v1_steps=40, v1_batch=32, v1_pi_steps=20, v1_pi_batch=32, v1_log_every=10),
            save_every_s=30.0, n_boot=200, ttg_s=(200, 100), mechanism_ttg_s=200, fidelity_ttg_s=(200, 100))

    def v3_expected(self) -> S.StudyConfig:
        """The v3 settings this v4 config expects in the v3 folder."""
        return S.StudyConfig.smoke() if self.mode == "smoke" else S.StudyConfig()

    def v3_training(self) -> S.StudyConfig:
        """The v3 config fit_resumable trains the v1-family variants with (v1's recipe and budget)."""
        return replace(self.v3_expected(), save_every_s=self.save_every_s, **self.v1_budget)


def git_commit() -> str:
    return S.git_commit()


# ============================================================================
# Layout and guarded writes
# ============================================================================
class ReadOnlyViolation(AssertionError):
    """A v4 write aimed outside the v4 study folder (or into a read-only input)."""


STAGES = {
    "00_study": "frozen settings (study_config.json), status.json, this index (README.md), manifest.csv",
    "01_checks": "E0 reproducibility checks and physics-model fidelity (CSV/JSON)",
    "02_models": "R0 and the 13 ablation variants per seed: checkpoint, variant.json, training log and curves",
    "03_open_loop": "open-loop footprint errors and along-track bias per variant and seed (CSV)",
    "04_flights": "one JSON result + one 1 Hz telemetry .npz (with the decision log) per flight; flights.csv",
    "05_analysis": "numbered tables (T##) and figures (F##) with captions.json; cache/ holds bootstrap results",
    "06_report": "SpaceJEPA_v4_report.pdf, generated from 05_analysis",
}


def _inside(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


class LayoutV4(S.Layout):
    """Paths of the v4 study folder (v4_plan.md section 3) and its guarded write helpers."""

    def __init__(self, root, read_only: Sequence = ()):
        super().__init__(root)
        self.checks = self.root / "01_checks"
        self.cache = self.analysis / "cache"
        self.testset = self.checks                  # (v3's 01_testset slot; v4 has no test-set stage)
        self.read_only = [Path(p).resolve() for p in read_only if p]
        r = self.root.resolve()
        for ro in self.read_only:
            if _inside(r, ro) or _inside(ro, r):
                raise ReadOnlyViolation(f"the v4 folder {r} and the read-only input {ro} overlap")

    def make(self) -> "LayoutV4":
        for p in (self.study, self.checks, self.models, self.open_loop, self.flights / "results",
                  self.flights / "telemetry", self.tables, self.figures, self.cache, self.report):
            self.guard(p)
            p.mkdir(parents=True, exist_ok=True)
        return self

    # ---- the guard ----
    def guard(self, path) -> Path:
        """Raise unless path is inside the v4 folder and outside every read-only input."""
        p = Path(path).resolve()
        if not _inside(p, self.root.resolve()):
            raise ReadOnlyViolation(f"v4 refuses to write outside its study folder: {p} (folder {self.root})")
        for ro in self.read_only:
            if _inside(p, ro):
                raise ReadOnlyViolation(f"v4 refuses to write under the read-only input {ro}: {p}")
        return Path(path)

    def atomic_bytes(self, path, data: bytes) -> None:
        S.atomic_bytes(self.guard(path), data)

    def write_json(self, path, obj) -> None:
        self.atomic_bytes(path, json.dumps(obj, indent=1, default=S._json_default).encode())

    def write_csv(self, df, path, **kw) -> None:
        self.atomic_bytes(path, df.to_csv(index=kw.pop("index", False), **kw).encode())

    def save_fig(self, fig, path, dpi: int = 150) -> None:
        buf = io.BytesIO()
        fig.savefig(buf, format=Path(path).suffix[1:] or "png", dpi=dpi, bbox_inches="tight")
        self.atomic_bytes(path, buf.getvalue())

    def torch_save(self, obj, path) -> None:
        buf = io.BytesIO()
        torch.save(obj, buf)
        self.atomic_bytes(path, buf.getvalue())

    def write_flight(self, rec: dict) -> None:
        """Telemetry first, then the result JSON (its presence marks the flight done)."""
        fid = rec["result"]["flight_id"]
        buf = io.BytesIO()
        np.savez_compressed(buf, **rec["telemetry"])
        self.atomic_bytes(self.telemetry_path(fid), buf.getvalue())
        self.write_json(self.result_path(fid), rec["result"])

    def variant_json(self, name: str, seed: int) -> Path:
        return self.model_dir(name, seed) / "variant.json"


def v3_layout(root) -> S.Layout:
    """The read-only v3 folder."""
    L3 = S.Layout(root)
    if not (L3.study / "study_config.json").exists():
        raise FileNotFoundError(f"no v3 study at {root} (00_study/study_config.json missing)")
    return L3


def update_manifest(L: LayoutV4):
    """00_study/manifest.csv (study_v3.update_manifest) and a v4 README.md."""
    L.guard(L.study / "manifest.csv")
    man = S.update_manifest(L)
    lines = ["# SpaceJEPA v4 study folder", "",
             "Auto-written by `jepa/study_v4.py` (do not edit). The plan is `jepa/v4_plan.md` in the repo.", ""]
    cfg = L.study / "study_config.json"
    if cfg.exists():
        c = json.loads(cfg.read_text())
        lines += [f"* git commit: `{c.get('git_commit')}`", f"* created: {c.get('created')}",
                  f"* mode: {c['config'].get('mode')}", f"* v3 folder checked: `{c.get('v3', {}).get('root')}` "
                  f"(git `{str(c.get('v3', {}).get('git_commit'))[:10]}`)", ""]
    st = L.study / "status.json"
    if st.exists():
        s = json.loads(st.read_text())
        lines += [f"Status ({s.get('time')}): {s.get('flights_done')} / {s.get('flights_total')} flights, "
                  f"{s.get('models_done')} / {s.get('models_total')} trained models.", ""]
    lines += ["| folder | what | files | MB |", "|---|---|---|---|"]
    for stage, what in STAGES.items():
        sub = man[man.stage == stage]
        lines.append(f"| `{stage}/` | {what} | {len(sub)} | {sub['size'].sum() / 1e6:.1f} |")
    caps = L.analysis / "captions.json"
    if caps.exists():
        lines += ["", "## Figures and tables", ""]
        for k, v in sorted(json.loads(caps.read_text()).items()):
            lines.append(f"* `{k}`: {v.get('caption', '')} (data: `{v.get('data', '')}`)")
    L.atomic_bytes(L.study / "README.md", ("\n".join(lines) + "\n").encode())
    return man


def freeze_config_v4(cfg: StudyConfigV4, L: LayoutV4, L3: S.Layout) -> dict:
    """
    Check the v3 folder's study_config.json against the v3 settings v4 expects
    (V3_FROZEN), then write 00_study/study_config.json on first use. Later
    sessions must use the same settings (except RUNTIME_KEYS): a difference
    raises (start a new study folder instead).
    """
    v3doc = json.loads((L3.study / "study_config.json").read_text())
    want3 = json.loads(json.dumps(asdict(cfg.v3_expected()), default=S._json_default))
    diff3 = {k: (v3doc["config"].get(k), want3.get(k)) for k in V3_FROZEN if v3doc["config"].get(k) != want3.get(k)}
    if diff3:
        raise ValueError(f"the v3 folder {L3.root} was run with other settings than v4 expects: {diff3}")
    if cfg.mode == "full":
        same = {"r0_wm": (cfg.r0_wm, want3["v2_wm"]), "r0_pi": (cfg.r0_pi, want3["v2_pi"]),
                "r0_lambda": (cfg.r0_lambda, want3["v2_lambda"]), "r0_arch": (cfg.r0_arch, want3["v2_arch"])}
        bad = {k: v for k, v in same.items() if json.loads(json.dumps(v[0], default=S._json_default)) != v[1]}
        if bad:
            raise ValueError(f"R0 must be v2's recipe (v4_plan.md section 4.3): {bad}")
    path = L.study / "study_config.json"
    want = json.loads(json.dumps(asdict(cfg), default=S._json_default))
    if path.exists():
        have = json.loads(path.read_text())
        diff = {k: (have["config"].get(k), v) for k, v in want.items()
                if k not in RUNTIME_KEYS and have["config"].get(k) != v}
        if diff:
            raise ValueError(f"settings differ from the frozen {path}: {diff}. Use a new study folder.")
        return have
    doc = {"config": want, "git_commit": git_commit(), "created": time.strftime("%Y-%m-%d %H:%M:%S"),
           "environment": S.environment_info(),
           "v3": {"root": str(L3.root), "git_commit": v3doc.get("git_commit"), "created": v3doc.get("created"),
                  "config": {k: v3doc["config"].get(k) for k in V3_FROZEN}}}
    L.write_json(path, doc)
    return doc


# ============================================================================
# Arms and jobs
# ============================================================================
def load_testset(cfg: StudyConfigV4):
    import pandas as pd
    cases = pd.read_csv(SIM_DIR / cfg.testset_csv)
    if cfg.cases is not None:
        cases = cases[cases.case.isin(cfg.cases)].reset_index(drop=True)
    return cases


def gname(g: float) -> str:
    return f"{g:g}g"


def _base(c: dict, atm: str) -> dict:
    return {"case": int(c["case"]), "set": c["set"], "atmosphere": atm, "ic": TS.ic_of(c),
            "domain": dict(TS.ATMOSPHERES[atm])}


def _learned(cfg, H, c, atm, arm, load, g, seed, model, stage, w_unc=0.0, variant=None) -> dict:
    s = F.crn_seed(seed, int(c["case"]), atm, 0)
    return {**_base(c, atm), "controller": arm, "load": load, "g_limit": g, "seed": int(seed), "mppi_rep": 0,
            "mppi_seed": s, "horizon": H, "mppi": {"seed": s, "g_limit": g, "horizon": H, "w_unc_km": w_unc},
            "model": model, "stage": stage, "variant": variant}


def _classical(c, atm, arm, stage, knows_truth: bool) -> dict:
    return {**_base(c, atm), "controller": arm, "load": "own", "g_limit": None, "seed": None, "mppi_rep": None,
            "classical_knows_truth": bool(knows_truth), "model": {"kind": "classical"}, "stage": stage}


def v3_model(L3: S.Layout, arm: str, seed: int) -> dict:
    """A v3 checkpoint as a v4 model spec: v1, v2 (ensemble) or v2single (member 0)."""
    if arm == "v1":
        return {"kind": "v1", "ckpt": str(L3.ckpt("v1", seed)), "members": None}
    return {"kind": "v2", "ckpt": str(L3.ckpt("v2", seed)), "members": [0] if arm == "v2single" else None}


def variant_model(L: LayoutV4, name: str, seed: int) -> dict:
    from jepa import ablate_v4 as AB
    return {"kind": "v1" if AB.family(name) == "v1" else "v2", "ckpt": str(L.ckpt(name, seed)), "members": None}


def e0b_failed(L: LayoutV4) -> bool:
    """True once E0b has been evaluated and failed (then the v3 learned arms are re-flown, section 4.0)."""
    p = L.checks / "E0b_summary.json"
    return p.exists() and not json.loads(p.read_text()).get("pass", True)


def build_jobs(cfg: StudyConfigV4, L: LayoutV4, L3: S.Layout, stages: Sequence[str] = STAGE_ORDER) -> List[dict]:
    """Every flight of the study in run order (E0 -> E1 -> E2 -> E4 -> E3), seed-major within each stage."""
    cases = load_testset(cfg)
    H = S.horizon_steps(cfg.v3_expected())
    rows = cases.to_dict("records")
    conds = [(c, a) for c in rows for a in cfg.atmospheres]
    by_case = {int(c["case"]): c for c in rows}
    gp = cfg.g_primary
    matched = None
    jobs: List[dict] = []

    def need_matched():
        nonlocal matched
        if matched is None:
            matched = F4.classical_peak_g(L3.root, [int(c["case"]) for c in rows], cfg.atmospheres)
        return matched

    for stage in stages:
        if stage == "E0":
            jobs += [_classical(by_case[k], a, "classical_E0a", "E0", False)
                     for k in cfg.e0a_cases if k in by_case for a in cfg.atmospheres]
            refly = e0b_failed(L)
            for seed in cfg.seeds:
                arms = [(a3, a4) for a3, a4 in cfg.e0b_pairs if seed == cfg.e0b_seed]
                if refly:
                    arms = [(a3, f"{a3}_v4run") for a3 in ("v1", "v2", "v2single")]
                for a3, a4 in arms:
                    jobs += [_learned(cfg, H, c, a, a4, gname(gp), gp, seed, v3_model(L3, a3, seed), "E0",
                                      cfg.v3_expected().w_unc_km if a3 == "v2" else 0.0) for c, a in conds]
        elif stage == "E1":
            jobs += [_classical(c, a, "classical_oracle", "E1", True) for c, a in conds]
            for seed in cfg.seeds:
                for arm, kind, load in (("mppi_nominal", "physics_nominal", gname(gp)),
                                        ("mppi_nominal", "physics_nominal", "matched"),
                                        ("mppi_truth", "physics_truth", gname(gp)),
                                        ("mppi_prior", "prior", gname(gp))):
                    m = {"kind": kind, "ckpt": None, "members": None}
                    for c, a in conds:
                        g = need_matched()[(int(c["case"]), a)] if load == "matched" else gp
                        jobs.append(_learned(cfg, H, c, a, arm, load, g, seed, m, "E1"))
        elif stage == "E2":
            w = cfg.v3_expected().w_unc_km
            for seed in cfg.seeds:
                for load in [gname(g) for g in cfg.g_extra] + ["matched"]:
                    for c, a in conds:
                        g = need_matched()[(int(c["case"]), a)] if load == "matched" else float(load[:-1])
                        jobs.append(_learned(cfg, H, c, a, "v2", load, g, seed, v3_model(L3, "v2", seed), "E2", w))
        elif stage == "E4":
            c = by_case.get(cfg.diag_case)
            if c is None:
                continue
            w = cfg.v3_expected().w_unc_km
            for seed in cfg.seeds:
                for arm in ("v2", "v2single"):
                    for g in (gp, cfg.g_legacy):
                        for a in cfg.diag_atmospheres:
                            jobs.append(_learned(cfg, H, c, a, "diag_case10", f"{arm}_{gname(g)}", g, seed,
                                                 v3_model(L3, arm, seed), "E4", w if arm == "v2" else 0.0))
        elif stage == "E3":
            for seed in cfg.seeds:
                for name in cfg.variants:
                    m = variant_model(L, name, seed)
                    jobs += [_learned(cfg, H, c, a, name, gname(gp), gp, seed, m, "E3", variant=name)
                             for c, a in conds]
    return jobs


def flight_done(L: S.Layout, job: dict) -> bool:
    return L.result_path(F.flight_id(job)).exists()


def group_key(job: dict) -> tuple:
    m = job["model"]
    return (job["controller"], job["load"], job["seed"], m["kind"], m.get("ckpt"),
            tuple(m["members"]) if m.get("members") is not None else None)


def groups(jobs: Sequence[dict]) -> List[Tuple[tuple, List[dict]]]:
    """Model-run groups (arm, load, seed, model) of the non-classical jobs, in queue order."""
    out: Dict[tuple, List[dict]] = {}
    for j in jobs:
        if j["model"]["kind"] != "classical":
            out.setdefault(group_key(j), []).append(j)
    return list(out.items())


# ============================================================================
# Queues
# ============================================================================
def _classical_worker(arg) -> dict:
    job, root, ro = arg
    rec = F4.fly_classical_v4(job)
    LayoutV4(root, ro).write_flight(rec)
    return rec["result"]


def run_classical_v4(jobs: Sequence[dict], L: LayoutV4, workers: int = 1, log=print) -> List[dict]:
    """Fly every pending classical job (E0a re-fly, classical_oracle) in a CPU process pool."""
    todo = [j for j in jobs if j["model"]["kind"] == "classical" and not flight_done(L, j)]
    n_all = sum(j["model"]["kind"] == "classical" for j in jobs)
    log(f"[classical] {n_all - len(todo)}/{n_all} done, {len(todo)} to fly on {workers} worker(s)")
    if not todo:
        return []
    order = {"mismatch_high": 0, "nominal": 1, "mismatch": 2, "out_of_range": 3}
    todo.sort(key=lambda j: order.get(j["atmosphere"], 4))
    args = [(j, str(L.root), [str(p) for p in L.read_only]) for j in todo]
    out, t0 = [], time.time()
    pool = mp.get_context("spawn").Pool(int(workers)) if workers > 1 else None
    it = pool.imap_unordered(_classical_worker, args) if pool else map(_classical_worker, args)
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


PROGRESS_DIR = Path(os.environ.get("TMPDIR", os.environ.get("TEMP", "/tmp"))) / "jepa_v4_progress"


def _flight_worker(arg) -> int:
    """Pool worker: one shard of a model run in lockstep, each flight written as it lands."""
    shard, model, device, mppi, max_rows, chunk, root, ro, threads, verbose, prog = arg
    torch.set_num_threads(int(threads))
    adapter = F4.load_adapter_v4(model["kind"], model.get("ckpt"), device, model.get("members"), chunk, keep_only=True)
    L = LayoutV4(root, ro)

    def report(p):
        if prog:
            try:
                S.atomic_bytes(prog, json.dumps(p).encode())
            except OSError:
                pass

    F4.run_lockstep_v4(shard, adapter, mppi, n_sub=int(mppi.get("n_sub", 2)), max_rows=max_rows,
                       on_done=L.write_flight, verbose=verbose, progress=report)
    return len(shard)


def run_flights_v4(cfg: StudyConfigV4, jobs: Sequence[dict], L: LayoutV4, device: str = "cuda", workers: int = 1,
                   log=print, report_every_s: float = 300.0) -> int:
    """
    Fly every pending learned or physics job, one model run (arm, load, seed,
    model) at a time in queue order, its flights in lockstep, split over
    `workers` processes sharing the device. A run whose checkpoint does not
    exist yet (E3 before training) is skipped with a note. Returns flights flown.
    """
    mppi = dict(cfg.v3_expected().mppi)
    gs = groups(jobs)
    ready = [(k, g) for k, g in gs if g[0]["model"].get("ckpt") is None or os.path.exists(g[0]["model"]["ckpt"])]
    n_all = sum(len(g) for _, g in ready)
    n_done0 = sum(flight_done(L, j) for _, g in ready for j in g)
    log(f"[flights] {n_done0}/{n_all} done in {len(ready)} model runs ({len(gs) - len(ready)} runs wait for "
        f"training); {workers} worker(s) on {device}")
    pool = None
    t_all, n_flown = time.time(), 0
    ro = [str(p) for p in L.read_only]
    try:
        for key, g in ready:
            todo = [j for j in g if not flight_done(L, j)]
            if not todo:
                continue
            if pool is None and workers > 1:
                pool = mp.get_context("spawn").Pool(int(workers))
            arm, load, seed, kind = key[:4]
            model = g[0]["model"]
            t0 = time.time()
            shards = [sh for sh in (todo[i::max(1, workers)] for i in range(max(1, workers))) if sh]
            PROGRESS_DIR.mkdir(parents=True, exist_ok=True)
            progs = [PROGRESS_DIR / f"worker{i}.json" for i in range(len(shards))]
            for f in progs:
                f.unlink(missing_ok=True)
            rows = cfg.max_rows_physics if kind in ("physics_nominal", "physics_truth", "prior") else cfg.max_rows
            threads = 1 if workers > 1 or device == "cuda" else max(1, (os.cpu_count() or 2) // 2)
            args = [(sh, model, device, mppi, rows, cfg.v2_chunk, str(L.root), ro, threads, pool is None,
                     str(pf) if pool is not None else None) for sh, pf in zip(shards, progs)]
            if pool is None:
                for a in args:
                    _flight_worker(a)
            else:
                res = pool.map_async(_flight_worker, args)
                while True:
                    try:
                        res.get(timeout=report_every_s)
                        break
                    except mp.TimeoutError:
                        log(f"    ... {arm} {load} seed {seed}: {(time.time() - t0) / 60:.0f} min, "
                            f"{S._progress_line(progs, len(todo))}")
            n_flown += len(todo)
            el = time.time() - t_all
            left = n_all - n_done0 - n_flown
            log(f"  [flights] {arm} {load} seed {seed}: {len(todo)} flights in {(time.time() - t0) / 60:.1f} min; "
                f"{left} left in this call (~{el / n_flown * left / 60:.0f} min)")
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    return n_flown


def collect_flights(L: LayoutV4, write: bool = True):
    """Rebuild 04_flights/flights.csv from every v4 result JSON."""
    df = S.collect_flights(L, write=False)
    if write and len(df):
        L.write_csv(df, L.flights / "flights.csv")
    return df


# ============================================================================
# Status
# ============================================================================
def _count_done(L: S.Layout, jobs: Sequence[dict]) -> Tuple[int, List[float]]:
    """Finished flights among jobs (one os.scandir per result folder) and their file times."""
    by_dir: Dict[Path, set] = {}
    for j in jobs:
        p = L.result_path(F.flight_id(j))
        by_dir.setdefault(p.parent, set()).add(p.name)
    n, times = 0, []
    for d, names in by_dir.items():
        try:
            with os.scandir(d) as it:
                for e in it:
                    if e.name in names:
                        n += 1
                        times.append(e.stat().st_mtime)
        except FileNotFoundError:
            pass
    return n, times


def status(cfg: StudyConfigV4, L: LayoutV4, jobs: Sequence[dict], write: bool = True) -> dict:
    """Done / total per stage and arm, the measured rate, and hours left; also every model's training state."""
    from jepa import ablate_v4 as AB
    arms: Dict[tuple, List[dict]] = {}
    for j in jobs:
        arms.setdefault((j["stage"], j["controller"], j["load"]), []).append(j)
    rows = []
    for (stage, arm, load), js in arms.items():
        n, times = _count_done(L, js)
        rate = (len(times) - 1) / max(max(times) - min(times), 1.0) * 3600 if len(times) >= 2 else float("nan")
        rows.append({"stage": stage, "arm": arm, "load": load, "done": n, "total": len(js),
                     "flights_per_h": round(rate, 1) if np.isfinite(rate) else None})
    train = []
    for name in cfg.variants:
        for seed in cfg.seeds:
            d = L.model_dir(name, seed)
            st = "done" if (d / "ckpt.pt").exists() else (
                "partial" if any((d / f).exists() for f in ("train_state_wm.pt", "train_state_pi.pt")) else "todo")
            sec = None
            if st == "done" and L.variant_json(name, seed).exists():
                sec = json.loads(L.variant_json(name, seed).read_text()).get("train_s")
            train.append({"variant": name, "seed": seed, "family": AB.family(name), "state": st, "train_s": sec})
    stages = []
    for stage in STAGE_ORDER:
        r = [x for x in rows if x["stage"] == stage]
        done, total = sum(x["done"] for x in r), sum(x["total"] for x in r)
        rates = [x["flights_per_h"] for x in r if x["flights_per_h"]]
        rate = float(np.median(rates)) if rates else None
        stages.append({"stage": stage, "flights_done": done, "flights_total": total,
                       "flights_per_h": rate, "hours_left": round((total - done) / rate, 2) if rate else None})
    secs = [t["train_s"] for t in train if t["train_s"]]
    left = sum(t["state"] != "done" for t in train)
    doc = {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "stages": stages, "arms": rows, "training": train,
           "flights_done": int(sum(x["done"] for x in rows)), "flights_total": int(sum(x["total"] for x in rows)),
           "models_done": int(sum(t["state"] == "done" for t in train)), "models_total": len(train),
           "training_hours_left": round(left * float(np.mean(secs)) / 3600, 2) if secs else None,
           "e0b_failed": e0b_failed(L)}
    if write:
        L.write_json(L.study / "status.json", doc)
    return doc


# ============================================================================
# E0 checks (v4_plan.md section 4.0) and the physics-model fidelity (4.1.3)
# ============================================================================
def e0_checks(cfg: StudyConfigV4, L: LayoutV4, L3: S.Layout, log=print) -> dict:
    """
    Evaluate whatever E0 flights exist: 01_checks/E0a.csv, E0b.csv and
    E0b_summary.json (pass/fail once every E0b flight has landed). Returns
    {"E0a": pass|None, "E0b": pass|None, ...}; None = flights still missing.
    """
    from jepa import analyze_v3 as A3
    out = {}
    cases = [int(c) for c in load_testset(cfg).case]
    a = F4.e0a_table(L3.root, L.root, [c for c in cfg.e0a_cases if c in cases], cfg.atmospheres)
    n_a = len([c for c in cfg.e0a_cases if c in cases]) * len(cfg.atmospheres)
    if len(a):
        L.write_csv(a, L.checks / "E0a.csv")
    ea = a[a.check == "E0a"] if len(a) else a
    out["E0a"] = bool(ea.identical.all()) if len(ea) == n_a else None
    orc = a[a.check == "oracle_vs_classical_nominal"] if len(a) else a
    out["oracle_equals_classical_nominal"] = bool(orc.identical.all()) if len(orc) else None
    out["oracle_vs_classical_nominal_max_d_miss_km"] = float(orc.d_miss_km.max()) if len(orc) else None
    b = F4.e0b_table(L3.root, L.root, cfg.e0b_pairs, cfg.e0b_seed, cases, cfg.atmospheres, gname(cfg.g_primary))
    n_b = len(cfg.e0b_pairs) * len(cases) * len(cfg.atmospheres)
    if len(b):
        L.write_csv(b, L.checks / "E0b.csv")
    if len(b) == n_b and n_b:
        summ = {"n": int(len(b)), "outcome_agreement": float(b.same_outcome.mean()),
                "success5_agreement": float(b.same_success5.mean())}
        per = []
        for (a3, a4), g in b.groupby(["v3_arm", "v4_arm"]):
            X = np.full((1, len(cfg.atmospheres), len(cases)), np.nan)
            for r in g.itertuples():
                X[0, cfg.atmospheres.index(r.atmosphere), cases.index(r.case)] = r.d_log_miss
            bt = A3.hboot(X, np.zeros_like(X), "mean", cfg.n_boot, seed=4)
            d = g.d_log_miss.to_numpy()
            per.append({"v3_arm": a3, "v4_arm": a4, "n": int(len(g)), "outcome_agreement": float(g.same_outcome.mean()),
                        "mean_d_log_miss": bt["diff"], "ci_lo": bt["ci_lo"], "ci_hi": bt["ci_hi"],
                        "sd_d_log_miss": float(np.std(d)), "p05": float(np.percentile(d, 5)),
                        "p95": float(np.percentile(d, 95)), "max_abs": float(np.abs(d).max()),
                        "frac_identical_miss": float((g.miss_v3_km == g.miss_v4_km).mean())})
        ok = all(p["outcome_agreement"] >= cfg.e0b_agree_min and p["ci_lo"] <= 0.0 <= p["ci_hi"] for p in per)
        summ.update({"pairs": per, "pass": bool(ok)})
        if not (L.checks / "E0b_summary.json").exists() or \
                json.loads((L.checks / "E0b_summary.json").read_text()).get("pass") != ok:
            log(f"[E0b] {'pass' if ok else 'FAIL: the v3 learned arms will be re-flown with the v4 runner (section 4.0); record it in v4_plan.md section 15'}")
        L.write_json(L.checks / "E0b_summary.json", summ)
        out["E0b"] = bool(ok)
    else:
        out["E0b"] = None
    c = L.checks / "E0c.csv"
    if c.exists():
        import pandas as pd
        out["E0c"] = bool(pd.read_csv(c).identical.all())
    else:
        out["E0c"] = None
    L.write_json(L.checks / "E0_summary.json", out)
    return out


def e0c_run(cfg: StudyConfigV4, L: LayoutV4, L3: S.Layout, log=print) -> bool:
    """E0c on the CPU with seed 0's v3 checkpoints (cached in 01_checks/E0c.csv)."""
    import pandas as pd
    p = L.checks / "E0c.csv"
    if p.exists():
        return bool(pd.read_csv(p).identical.all())
    cases = load_testset(cfg).set_index("case")
    H = S.horizon_steps(cfg.v3_expected())
    jobs = []
    for c, a in cfg.e0c_conditions:
        row = {**cases.loc[c].to_dict(), "case": c}
        s = F.crn_seed(0, c, a, 0)
        jobs.append({**_base(row, a), "controller": "E0c", "load": gname(cfg.g_primary), "seed": 0, "mppi_rep": 0,
                     "mppi_seed": s, "horizon": H, "mppi": {"seed": s, "g_limit": cfg.g_primary, "horizon": H,
                                                            "w_unc_km": 0.0, "n_samples": cfg.e0c_samples}})
    t0 = time.time()
    rows = F4.e0c_check(jobs, str(L3.ckpt("v1", 0)), str(L3.ckpt("v2", 0)), dict(cfg.v3_expected().mppi))
    df = pd.DataFrame(rows)
    L.write_csv(df, p)
    log(f"[E0c] run_lockstep_v4 == run_lockstep_v3 on the CPU: {bool(df.identical.all())} "
        f"({len(df)} flights, {time.time() - t0:.0f}s)")
    return bool(df.identical.all())


def physics_checks(cfg: StudyConfigV4, L: LayoutV4, L3: S.Layout, log=print) -> dict:
    """Section 4.1.3 checks 1 and 2 (cached in 01_checks/)."""
    import pandas as pd
    from jepa import physics_v4 as PH
    out = {}
    p1 = L.checks / "E0_physics_port.json"
    if not p1.exists():
        L.write_json(p1, PH.port_check(cfg.port_check_states))
    out["port"] = json.loads(p1.read_text())
    p2 = L.checks / "E0_physics_fidelity.csv"
    if not p2.exists():
        tels = []
        for c in cfg.fidelity_cases:
            for a in cfg.fidelity_atmospheres:
                fid = f"v2/{gname(cfg.g_primary)}/seed0/mppi0/{a}_case{int(c):02d}"
                if L3.result_path(fid).exists() and L3.telemetry_path(fid).exists():
                    tels.append((json.loads(L3.result_path(fid).read_text()), S.load_telemetry(L3, fid)))
        if tels:
            L.write_csv(pd.DataFrame(PH.open_loop_fidelity(tels, cfg.fidelity_ttg_s)), p2)
            log(f"[physics] open-loop fidelity on {len(tels)} recorded v3 flights")
    if p2.exists():
        fd = pd.read_csv(p2)
        out["fidelity"] = fd.groupby("ttg_s")[["err_truth_km", "err_nominal_km"]].median().round(3).to_dict()
    return out
