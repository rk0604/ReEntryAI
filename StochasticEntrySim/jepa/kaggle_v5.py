"""
jepa/kaggle_v5.py

Every v5 RL run (tuning and training) trains on Kaggle's free GPU, on one platform and one GPU type
(jepa/kaggle_rl.md, the user's decision of 2026-11-16). This module moves the files between the two places.

    export_inputs   Colab: 00_study/kaggle/to_kaggle/v5_rl_inputs.json, the frozen settings and FNPEG's tuning in one
                    file, for the user to upload as a Kaggle dataset
    copy_inputs     Kaggle: that file's contents into the session's study folder (an existing different file raises)
    check_platform  Kaggle: 00_study/rl_platform.json at the first session; a later session on other hardware
                    (GPU type, torch, CUDA, platform) raises, so no RL run of the study trains elsewhere
    restore         Kaggle: the newest earlier session's state (v5_rl_state.tar + v5_rl_manifest.json, found among
                    the attached inputs), checked and unpacked into the session's study folder
    pack            Kaggle: the study folder as one tar in the session's output (Kaggle keeps few files per output),
                    with a manifest of every file's size and sha256; called after each run and at the end
    import_rl       Colab: a session's state into the Drive study folder, only finished RL runs, every file checked,
                    nothing overwritten (a different file already on Drive stops the import before any write)
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tarfile
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

STATE_TAR = "v5_rl_state.tar"
MANIFEST = "v5_rl_manifest.json"
INPUTS = "v5_rl_inputs.json"
INPUT_FILES = ("00_study/study_config.json", "03_tuning/fnpeg/chosen.json", "03_tuning/fnpeg/vf.json")
PLATFORM = "00_study/rl_platform.json"
HW_SAME = ("gpu", "torch", "cuda", "platform")
RL_ARMS = ("dagger_fnpeg", "bcppo_fnpeg", "ppo_hist", "sac_hist", "ppo_matched")


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def _rel(p: Path, root: Path) -> str:
    return p.relative_to(root).as_posix()


def _walk(root: Path, max_depth: int = 5):
    root = Path(root)
    if not root.exists():
        return
    base = len(root.parts)
    for dp, dirs, files in os.walk(root):
        if len(Path(dp).parts) - base >= max_depth:
            dirs[:] = []
        for f in files:
            yield Path(dp) / f


# ============================================================================
# Colab side: inputs out, results in
# ============================================================================
def kaggle_dir(L) -> Path:
    return L.study / "kaggle"


def export_inputs(L, log=print) -> Optional[Path]:
    """The frozen settings and FNPEG's tuning as one JSON (Colab, once FNPEG is tuned); None until then."""
    files = {}
    for rel in INPUT_FILES:
        p = L.root / rel
        if not p.exists():
            if rel != "03_tuning/fnpeg/vf.json":
                log(f"[kaggle] the Kaggle inputs wait for {rel} (FNPEG's tuning on a CPU runtime)")
                return None
            continue
        files[rel] = p.read_text(encoding="utf-8")
    doc = {"what": "v5 RL inputs for Kaggle (jepa/kaggle_rl.md)", "files": files,
           "sha256": {k: hashlib.sha256(v.encode()).hexdigest() for k, v in files.items()},
           "time": time.strftime("%Y-%m-%d %H:%M:%S")}
    out = kaggle_dir(L) / "to_kaggle" / INPUTS
    data = json.dumps(doc, indent=1).encode()
    if not out.exists() or out.read_bytes() != data:
        L.atomic_bytes(out, data)
    return out


def import_rl(L, src: Path, log=print) -> dict:
    """A Kaggle session's state (the folder holding v5_rl_state.tar and v5_rl_manifest.json) into the study folder."""
    src = Path(src)
    man_p, tar_p = src / MANIFEST, src / STATE_TAR
    if not (man_p.exists() and tar_p.exists()):
        raise FileNotFoundError(f"{src} must hold {STATE_TAR} and {MANIFEST} (a Kaggle session's output)")
    man = json.loads(man_p.read_text())
    if sha256_file(tar_p) != man["tar_sha256"]:
        raise ValueError(f"{tar_p} does not match its manifest (sha256): download both files of one session again")
    with tempfile.TemporaryDirectory() as td:
        _extract(tar_p, Path(td))
        root = Path(td)
        files = man["files"]
        bad = [r for r, (size, sha) in files.items() if (root / r).exists() and sha256_file(root / r) != sha]
        if bad:
            raise ValueError(f"files in {tar_p} do not match the manifest: {bad[:5]}")
        done_runs = {str(Path(r).parent.as_posix()) for r in files if r.endswith("/policy.pt")}
        plan: List[Tuple[str, str]] = []                       # (source rel, target rel)
        for r in sorted(files):
            t = _import_target(r, done_runs)
            if t:
                plan.append((r, t))
        # every finished run trained on the platform the study's RL is bound to
        plat = json.loads((root / PLATFORM).read_text()) if (root / PLATFORM).exists() else None
        if plat is None:
            raise ValueError(f"{tar_p} has no {PLATFORM}: not a v5 Kaggle RL session")
        have_p = L.root / PLATFORM
        if have_p.exists():
            have = json.loads(have_p.read_text())
            diff = {k: (have.get(k), plat.get(k)) for k in HW_SAME if have.get(k) != plat.get(k)}
            if diff:
                raise ValueError(f"this Kaggle state trained on other hardware than the RL already imported: {diff}")
        for run in done_runs:
            hp = root / run / "hardware.json"
            hw = json.loads(hp.read_text()) if hp.exists() else {}
            diff = {k: (plat.get(k), hw.get(k)) for k in HW_SAME if plat.get(k) != hw.get(k)}
            if diff:
                raise ValueError(f"RL run {run} did not train on the study's RL platform: {diff}")
        conflicts, new, same = [], [], 0
        for r, t in plan:
            dst = L.root / t
            if dst.exists():
                if sha256_file(dst) == files[r][1]:
                    same += 1
                elif not t.startswith("00_study/kaggle/logs/"):
                    conflicts.append(t)
            else:
                new.append((r, t))
        if conflicts:
            raise ValueError(f"{len(conflicts)} file(s) on Drive differ from this Kaggle state; nothing was imported "
                             f"(one RL run trained twice?): {conflicts[:5]}")
        logs = [(r, t) for r, t in plan if t.startswith("00_study/kaggle/logs/") and (r, t) not in new]
        for r, t in new + logs:                                # logs: the newest session's copy replaces the older
            data = (root / r).read_bytes()
            if len(data) > (4 << 20):
                L.save_large(data, L.root / t)
            else:
                L.atomic_bytes(L.root / t, data)
    out = {"session": man.get("seq"), "runs_finished": len(done_runs), "files_new": len(new), "files_same": same,
           "gpu": plat.get("gpu")}
    log(f"[kaggle] imported Kaggle session {man.get('seq')}: {out}")
    return out


def _import_target(rel: str, done_runs: set) -> Optional[str]:
    """Where a file of the Kaggle state goes in the Drive study folder (None: not imported)."""
    p = Path(rel)
    parts = p.parts
    name = p.name
    if name.startswith(".") or name.endswith(".tmp") or name.startswith("train_state"):
        return None
    if rel == PLATFORM:
        return rel
    if parts[:2] == ("00_study", "logs"):
        return f"00_study/kaggle/logs/{name}"
    if parts[0] == "02_policies" and len(parts) >= 4:
        return rel if "/".join(parts[:3]) in done_runs else None
    if parts[0] == "03_tuning" and len(parts) >= 2 and parts[1].startswith("rl_"):
        if len(parts) >= 4:                                    # a configuration's run folder
            return rel if "/".join(parts[:3]) in done_runs else None
        return rel                                             # config_<j>.json, chosen.json, configs.csv
    if parts[0] == "04_flights" and len(parts) >= 4 and parts[1] in ("results", "telemetry"):
        if any(parts[3].startswith(f"tune_{a}_c") for a in RL_ARMS):
            return rel
    return None


# ============================================================================
# Kaggle side
# ============================================================================
def find_inputs(search_roots: Sequence) -> Optional[Path]:
    hits = [p for r in search_roots for p in _walk(Path(r)) if p.name == INPUTS]
    return sorted(hits)[-1] if hits else None


def copy_inputs(L, search_roots: Sequence, log=print) -> dict:
    """v5_rl_inputs.json (an attached Kaggle dataset) into the study folder; a different existing file raises."""
    p = find_inputs(search_roots)
    if p is None:
        raise FileNotFoundError(f"no {INPUTS} among the inputs {list(map(str, search_roots))}: add the Kaggle dataset "
                                "made from the Colab study's 00_study/kaggle/to_kaggle/ folder")
    doc = json.loads(p.read_text(encoding="utf-8"))
    for rel, text in doc["files"].items():
        if hashlib.sha256(text.encode()).hexdigest() != doc["sha256"][rel]:
            raise ValueError(f"{p}: {rel} does not match its checksum")
        dst = L.root / rel
        if dst.exists():
            if dst.read_text(encoding="utf-8") != text:
                raise ValueError(f"{rel} in this session's state differs from {p}: the inputs changed after RL started")
        else:
            L.atomic_bytes(dst, text.encode("utf-8"))
    log(f"[kaggle] inputs from {p}: {sorted(doc['files'])}")
    return doc


def check_platform(L, hw: dict, log=print) -> dict:
    """The study's RL platform: written by the first session, every later session must match it."""
    p = L.root / PLATFORM
    if p.exists():
        have = json.loads(p.read_text())
        diff = {k: (have.get(k), hw.get(k)) for k in HW_SAME if have.get(k) != hw.get(k)}
        if diff:
            raise RuntimeError(f"this session's hardware differs from the study's RL platform {have}: {diff}. Every RL "
                               "run must train on the same platform (jepa/kaggle_rl.md): choose that accelerator and "
                               "pin the Kaggle environment (Settings -> Environment -> pin to original).")
        return have
    doc = {**hw, "since": time.strftime("%Y-%m-%d %H:%M:%S")}
    L.write_json(p, doc)
    log(f"[kaggle] RL platform fixed for the study: {doc}")
    return doc


def find_states(search_roots: Sequence) -> List[Tuple[int, Path, dict]]:
    out = []
    for r in search_roots:
        for p in _walk(Path(r)):
            if p.name == MANIFEST and (p.parent / STATE_TAR).exists():
                try:
                    man = json.loads(p.read_text())
                    out.append((int(man["seq"]), p.parent, man))
                except (OSError, ValueError, KeyError):
                    pass
    return sorted(out, key=lambda x: x[0])


def _extract(tar_p: Path, dst: Path) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_p) as tf:
        for m in tf.getmembers():
            if not (m.isfile() or m.isdir()) or m.name.startswith("/") or ".." in Path(m.name).parts:
                raise ValueError(f"{tar_p}: unexpected member {m.name!r}")
        tf.extractall(dst, **({"filter": "data"} if hasattr(tarfile, "data_filter") else {}))


def restore(L, search_roots: Sequence, log=print) -> int:
    """The newest earlier session's state into the study folder; returns its session number (0: none)."""
    states = find_states(search_roots)
    if not states:
        log("[kaggle] no earlier session's state among the inputs: this is session 1")
        return 0
    seq, d, man = states[-1]
    if sha256_file(d / STATE_TAR) != man["tar_sha256"]:
        raise ValueError(f"{d / STATE_TAR} does not match its manifest: attach that session's output again")
    with tempfile.TemporaryDirectory() as td:
        _extract(d / STATE_TAR, Path(td))
        n = 0
        for p in _walk(Path(td), 50):
            rel = _rel(p, Path(td))
            dst = L.root / rel
            if not dst.exists():
                L.guard(dst).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(p, dst)
                n += 1
    others = [s for s, _, _ in states[:-1]]
    log(f"[kaggle] restored session {seq} ({n} files; {len(man.get('runs_done', []))} RL runs finished)"
        + (f"; older states also attached: {others}" if others else ""))
    return seq


def pack(L, out_dir: Path, seq: int, log=print) -> dict:
    """The study folder as out_dir/v5_rl_state.tar + v5_rl_manifest.json (session number seq)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    root = Path(L.root)
    files: Dict[str, list] = {}
    tmp_tar = out_dir / f".{STATE_TAR}.tmp"
    with tarfile.open(tmp_tar, "w") as tf:
        for p in sorted(_walk(root, 50)):
            rel = _rel(p, root)
            if p.name.endswith(".tmp") or rel.startswith("05_analysis/cache/"):
                continue
            files[rel] = [p.stat().st_size, sha256_file(p)]
            tf.add(p, arcname=rel, recursive=False)
    os.replace(tmp_tar, out_dir / STATE_TAR)
    runs = sorted(str(Path(r).parent.as_posix()) for r in files if r.endswith("/policy.pt"))
    man = {"seq": int(seq), "time": time.strftime("%Y-%m-%d %H:%M:%S"), "tar_sha256": sha256_file(out_dir / STATE_TAR),
           "tar_bytes": (out_dir / STATE_TAR).stat().st_size, "runs_done": runs, "files": files}
    tmp = out_dir / f".{MANIFEST}.tmp"
    tmp.write_text(json.dumps(man, indent=1))
    os.replace(tmp, out_dir / MANIFEST)
    log(f"[kaggle] saved session {seq}: {len(files)} files, {len(runs)} finished RL runs, "
        f"{man['tar_bytes'] / 1e6:.1f} MB")
    return man
