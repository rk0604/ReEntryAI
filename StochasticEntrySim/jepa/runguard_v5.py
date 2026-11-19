"""
Unattended-run safety for the v5 notebook (v5_plan.md section 15, 2026-11-18): a stall watchdog, the session's unit
budget, failure records and 00_study/RUN_SUMMARY.txt.

Progress. Every study log line, every Drive write of the study code (LayoutV5.atomic_bytes, save_large) and every
worker's progress report calls touch(). A touch in any process (pool workers included) also refreshes a small beat
file in the local temp folder, so the watchdog in the notebook process sees the workers' progress too. A waiting
loop that only says "still waiting" logs with quiet(log), which does not count as progress. A known long silent
step runs inside allow(label, minutes), which raises the stall limit to that many minutes while it runs.

Watchdog. A daemon thread started at Setup. With no progress for STALL_MIN minutes it writes the reason to
00_study/logs/watchdog.log (and a local copy) and RUN_SUMMARY.txt, then calls the notebook's release_runtime
(Drive flush, then runtime.unassign). Each Drive step runs in its own thread with a timeout, so a hung Drive mount
cannot stop the release: after the timeout google.colab.runtime.unassign() is called directly.

Failures. stage(name) wraps a notebook stage: StageFailed (a stage stopped by its own failure rule) and any other
exception are recorded (00_study/logs/failures.log, RUN_SUMMARY.txt) and the notebook goes on with the next stage,
unless STOP_ON_FAILURE is set, which re-raises (the post_run_cell hook then releases the runtime). After every
stage RUN_SUMMARY.txt is rewritten, and when the session's unit budget is used up the runtime is released.

Stop rules (plan_jepa_needed_v5.md 4b, v5_plan.md 15 of 2026-11-19; on when STATE["stop_rules"], which only
ReEntryAI_JEPA_v5_needed.ipynb sets). A run that more compute will not fix stops the notebook: stop_run() writes
00_study/STOPPED.json (trigger id and evidence), records the failure, releases the runtime and raises StopRun. The
triggers: an E3 arm that diverged at every multiplier, an E3 validation-loss stop, an MPPI arm whose first flights
fail, FAILURE_LIMIT failure records in one session, and any stage that raises. The next Run all stops at Setup
(check_stopped) until ACK_STOP = "<trigger id>" is set in Settings. Watchdog, budget and 24 h releases are not stop
rules: those runs resume normally.

Nothing here imports the study modules at import time (study_v5 imports this module).
"""

from __future__ import annotations

import json
import os
import tempfile
import threading
import time
import traceback
from contextlib import contextmanager
from pathlib import Path
from typing import Callable, List, Optional

BEAT_DIR = Path(os.environ.get("JEPA_V5_BEAT_DIR") or Path(tempfile.gettempdir()) / "jepa_v5_beat")   # env: tests
_LAST = {"t": time.time(), "what": "start", "file_t": 0.0}
_ALLOW: List[float] = []                 # minutes allowed by the allow() blocks now running (innermost last)


class StageFailed(RuntimeError):
    """A stage stopped by its own failure rule (e.g. too many flight errors); recorded, then the notebook goes on."""


class StopRun(RuntimeError):
    """A stop rule fired (stop_run): STOPPED.json is written and the runtime released; nothing else runs."""

    def __init__(self, trigger: str, evidence: str):
        super().__init__(f"stop rule {trigger}: {evidence}")
        self.trigger, self.evidence = trigger, evidence


STOP_FILE = "STOPPED.json"           # in 00_study
FAILURE_LIMIT = 3                    # failure records in one session that stop the run (stop rules only)
NOT_COUNTED = ("cell", "stop", "E3 divergence")   # failure kinds not counted: the crash hook's own record, a stop
                                     # rule's (it already stopped), and an E3 divergence the learning-rate fallback
                                     # handles (its end, "did not train", is a stop rule of its own)


# ============================================================================
# Progress
# ============================================================================
def touch(what: str = "") -> None:
    """Progress. Cheap: the beat file is refreshed at most every 10 s per process."""
    now = time.time()
    _LAST["t"] = now
    if what:
        _LAST["what"] = what
    if now - _LAST["file_t"] > 10.0:
        _LAST["file_t"] = now
        try:
            BEAT_DIR.mkdir(parents=True, exist_ok=True)
            p = BEAT_DIR / f"{os.getpid()}"
            with open(p, "w") as f:
                f.write(f"{now:.0f} {what[:200]}")
        except OSError:
            pass


def last_progress() -> float:
    """The newest progress time of this process or of any worker's beat file."""
    t = _LAST["t"]
    try:
        with os.scandir(BEAT_DIR) as it:
            for e in it:
                try:
                    t = max(t, e.stat().st_mtime)
                except OSError:
                    pass
    except OSError:
        pass
    return t


def quiet(log: Callable) -> Callable:
    """The log function without the progress touch (for "still waiting" lines)."""
    return getattr(log, "quiet", log)


@contextmanager
def allow(label: str, minutes: float):
    """A known long silent step: the stall limit is at least `minutes` while it runs."""
    _ALLOW.append(float(minutes))
    touch(f"{label} (up to {minutes:.0f} min)")
    try:
        yield
    finally:
        _ALLOW.remove(float(minutes))
        touch(f"{label} done")


# ============================================================================
# Session: rate, budget, summary
# ============================================================================
UNITS_PER_HOUR_CPU = 0.08            # a standard CPU runtime (2 cores)
UNITS_PER_HOUR_CPU_HIGH = 0.26       # the High-RAM CPU runtime (8 cores; the user's Colab, 2026-11-17)


def units_per_hour(runtime: str, cpus: int, table: dict, override: Optional[float] = None) -> float:
    if override:
        return float(override)
    if runtime and runtime != "cpu":
        for k, v in table.items():
            if k != "cpu" and k.lower() in runtime.lower():
                return float(v)
    return UNITS_PER_HOUR_CPU_HIGH if int(cpus or 2) >= 8 else UNITS_PER_HOUR_CPU


class Summary:
    """00_study/RUN_SUMMARY.txt: this session's stages, failures, units and release reason (read on a phone)."""

    def __init__(self):
        self.L = None
        self.t0 = time.time()
        self.info: dict = {}
        self.stages: List[dict] = []
        self.failures: List[dict] = []
        self.released_why: Optional[str] = None
        self.next_note: Callable[[], str] = lambda: ""
        self.rate = 0.0
        self.budget_units: Optional[float] = None

    def configure(self, L, info: dict, rate: float, budget_units: Optional[float], next_note: Callable[[], str]):
        self.L, self.info, self.rate, self.budget_units, self.next_note = L, dict(info), float(rate), budget_units, \
            next_note
        self.t0, self.stages, self.failures, self.released_why = time.time(), [], [], None

    def units(self) -> float:
        return (time.time() - self.t0) / 3600.0 * self.rate

    def failure(self, kind: str, what: str, reason: str, path: Optional[str] = None) -> None:
        self.failures.append({"time": time.strftime("%Y-%m-%d %H:%M:%S"), "kind": kind, "what": what,
                              "reason": str(reason)[:600], "record": path})
        if self.L is not None:
            try:
                p = self.L.logs / "failures.log"
                self.L.guard(p).parent.mkdir(parents=True, exist_ok=True)
                with open(p, "a", encoding="utf-8") as f:
                    f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} [{kind}] {what}: {str(reason)[:600]}"
                            + (f" ({path})" if path else "") + "\n")
            except OSError:
                pass

    def text(self) -> str:
        h = (time.time() - self.t0) / 3600.0
        lines = ["SpaceJEPA v5: summary of the latest session (rewritten after every stage)", "",
                 f"session started {time.strftime('%Y-%m-%d %H:%M', time.localtime(self.t0))}, "
                 f"updated {time.strftime('%Y-%m-%d %H:%M')}",
                 *(f"{k}: {v}" for k, v in self.info.items()),
                 f"running for {h:.1f} h, about {self.units():.1f} units used at {self.rate:g} units/h"
                 + (f" (budget {self.budget_units:g} units)" if self.budget_units else ""), ""]
        lines.append("RELEASED: " + self.released_why if self.released_why else "status: running")
        lines += ["", "stages this session:"]
        for s in self.stages:
            lines.append(f"  {s['name']}: {s['status']}" + (f" ({s['minutes']:.0f} min)" if s.get("minutes") else "")
                         + (f" - {s['note']}" if s.get("note") else ""))
        if not self.stages:
            lines.append("  (none yet)")
        lines += ["", f"failures this session: {len(self.failures) or 'none'}"]
        for f in self.failures:
            lines.append(f"  {f['time']} [{f['kind']}] {f['what']}: {f['reason']}"
                         + (f" -> {f['record']}" if f.get("record") else ""))
        try:
            note = self.next_note() or ""
        except Exception as e:                  # never let the summary itself fail
            note = f"(next-pass note unavailable: {e})"
        lines += ["", "next pass: " + (note or "nothing reported yet")]
        return "\n".join(lines) + "\n"

    def write(self) -> None:
        if self.L is None:
            return
        try:
            self.L.atomic_bytes(self.L.study / "RUN_SUMMARY.txt", self.text().encode("utf-8"))
        except Exception as e:
            print(f"[summary] RUN_SUMMARY.txt not written: {e}", flush=True)

    def released(self, why: str) -> None:
        self.released_why = f"{time.strftime('%Y-%m-%d %H:%M')}: {why}"


SUMMARY = Summary()
STATE = {"release": None, "limit": None, "stop_on_failure": False, "stall_min": 45.0, "stop_rules": False,
         "stopped": None}


def configure(L, info: dict, rate: float, release: Callable[[str], None], limit=None,
              budget_units: Optional[float] = None, stop_on_failure: bool = False,
              next_note: Callable[[], str] = lambda: "", stop_rules: bool = False) -> None:
    STATE.update({"release": release, "limit": limit, "stop_on_failure": bool(stop_on_failure),
                  "stop_rules": bool(stop_rules), "stopped": None})
    SUMMARY.configure(L, info, rate, budget_units, next_note)
    SUMMARY.write()


def deadline(budget_units: Optional[float], rate: float) -> Optional[float]:
    """MAX_SESSION_UNITS as a wall-clock deadline at the session's rate (None: no budget)."""
    if not budget_units or rate <= 0:
        return None
    return time.time() + float(budget_units) / float(rate) * 3600.0


def release(why: str) -> None:
    SUMMARY.released(why)
    SUMMARY.write()
    fn = STATE.get("release")
    if fn is not None:
        fn(why)


def budget_stopped() -> bool:
    lim = STATE.get("limit")
    return bool(lim is not None and getattr(lim, "stopped", False))


def record_failure(L, path, doc: dict, kind: str, what: str, log=print) -> None:
    """A failure record (JSON at `path` inside the study folder), the failures log and RUN_SUMMARY.txt."""
    doc = {**doc, "time": time.strftime("%Y-%m-%d %H:%M:%S")}
    rel = None
    try:
        L.write_json(path, doc)
        rel = os.path.relpath(str(path), str(L.root))
    except Exception as e:
        log(f"[failure] record {path} not written: {e}")
    SUMMARY.failure(kind, what, doc.get("reason", ""), rel)
    log(f"[failure] {kind} {what}: {doc.get('reason', '')}" + (f" ({rel})" if rel else ""))
    SUMMARY.write()
    check_failure_count(L, log)


def counted_failures() -> List[dict]:
    return [f for f in SUMMARY.failures if f["kind"] not in NOT_COUNTED]


def check_failure_count(L, log=print) -> None:
    """Stop rule: FAILURE_LIMIT failure records of any kind in one session."""
    fs = counted_failures()
    if STATE.get("stop_rules") and not STATE.get("stopped") and len(fs) >= FAILURE_LIMIT:
        ev = "; ".join(f"[{f['kind']}] {f['what']}: {f['reason'][:160]}" for f in fs)
        stop_run(L, "failures_%d" % FAILURE_LIMIT, f"{len(fs)} failure records this session: {ev}", log)


def stop_run(L, trigger: str, evidence: str, log=print, extra: Optional[dict] = None) -> None:
    """A stop rule fired: 00_study/STOPPED.json, the failure record, RUN_SUMMARY.txt, the release, then StopRun."""
    if STATE.get("stopped"):                  # already stopped (and released) in this session
        raise StopRun(STATE["stopped"], evidence)
    STATE["stopped"] = trigger
    doc = {"trigger": trigger, "evidence": str(evidence)[:4000], "time": time.strftime("%Y-%m-%d %H:%M:%S"),
           "session": SUMMARY.info, "failures_this_session": SUMMARY.failures, **(extra or {}),
           "to_continue": f'look at the evidence, then set ACK_STOP = "{trigger}" in the Settings cell and Run all'}
    rel = None
    try:
        L.write_json(L.study / STOP_FILE, json_safe(doc))
        rel = os.path.relpath(str(L.study / STOP_FILE), str(L.root))
    except Exception as e:
        log(f"[stop] {STOP_FILE} not written: {e}")
    SUMMARY.failure("stop", trigger, evidence, rel)
    try:
        log(f"[stop] stop rule {trigger}: {evidence}")
    except Exception:
        pass
    release(f"stop rule {trigger} ({str(evidence)[:300]}); see 00_study/{STOP_FILE}, then set "
            f'ACK_STOP = "{trigger}" in Settings to continue')
    raise StopRun(trigger, evidence)


def check_stopped(L, ack: Optional[str], log=print) -> Optional[dict]:
    """Setup: a STOPPED.json from an earlier session stops this Run all unless ACK_STOP is its trigger id; then it
    moves to 00_study/stopped_history/ and the run goes on. Returns the acknowledged record (or None)."""
    p = L.study / STOP_FILE
    if not p.exists():
        return None
    doc = json.loads(p.read_text(encoding="utf-8"))
    trig = doc.get("trigger", "?")
    show = json.dumps({k: doc.get(k) for k in ("trigger", "time", "evidence")}, indent=1)
    if ack != trig:
        raise RuntimeError(f"an earlier session stopped by its stop rule {trig!r} (00_study/{STOP_FILE}):\n{show}\n"
                           f'Look at the evidence, then set ACK_STOP = "{trig}" in the Settings cell and Run all '
                           f"(ACK_STOP is now {ack!r}).")
    dst = L.study / "stopped_history" / f"STOPPED_{time.strftime('%Y%m%d_%H%M%S')}_{_safe_id(trig)}.json"
    L.write_json(dst, {**doc, "acknowledged": time.strftime("%Y-%m-%d %H:%M:%S")})
    L.guard(p).unlink()
    log(f"[stop] stop rule {trig!r} acknowledged (ACK_STOP); its record moved to {dst.relative_to(L.root)}")
    return doc


def _safe_id(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in str(s))[:80]


@contextmanager
def stage(name: str, log=print):
    """A notebook stage: failures are recorded and the notebook goes on (STOP_ON_FAILURE re-raises)."""
    t0 = time.time()
    rec = {"name": name, "status": "running"}
    SUMMARY.stages.append(rec)
    touch(f"stage {name}")
    try:
        yield rec
        rec["status"] = rec.get("result") or "done"
    except KeyboardInterrupt:
        rec["status"] = "interrupted"
        raise
    except StopRun:
        rec["status"] = "STOPPED"
        raise
    except Exception as e:
        rec["status"] = "FAILED"
        rec["note"] = f"{type(e).__name__}: {str(e)[:300]}"
        SUMMARY.failure("stage", name, rec["note"])
        try:
            log(f"[stage {name}] FAILED: {type(e).__name__}: {e}\n{traceback.format_exc()[-3000:]}")
        except Exception:
            pass
        if STATE.get("stop_rules"):           # stop rule: a stage that raises stops the run
            rec["minutes"] = (time.time() - t0) / 60.0
            rec["status"] = "STOPPED"
            stop_run(SUMMARY.L, f"stage_failed:{_safe_id(name)}", rec["note"], log)
        if STATE["stop_on_failure"]:
            rec["minutes"] = (time.time() - t0) / 60.0
            SUMMARY.write()
            raise
        log(f"[stage {name}] recorded; going on with the next stage (STOP_ON_FAILURE = False)")
    finally:
        rec["minutes"] = (time.time() - t0) / 60.0
        if budget_stopped() and rec["status"] in ("done", "running"):
            rec["note"] = (rec.get("note", "") + " session budget reached").strip()
        SUMMARY.write()
    if STATE.get("stop_rules") and SUMMARY.L is not None:
        check_failure_count(SUMMARY.L, log)
    if budget_stopped():
        release("the session's unit budget (MAX_SESSION_UNITS) is used up; Run all again to continue")


# ============================================================================
# Watchdog
# ============================================================================
def _bounded(fn: Callable, timeout_s: float) -> bool:
    """Run fn in a daemon thread; True if it returned within timeout_s."""
    t = threading.Thread(target=fn, daemon=True)
    t.start()
    t.join(timeout_s)
    return not t.is_alive()


class Watchdog(threading.Thread):
    """Releases the runtime after stall_min minutes without progress (see the module docstring)."""

    def __init__(self, stall_min: float, on_stall: Callable[[str], None], check_s: float = 30.0):
        super().__init__(daemon=True, name="jepa_v5_watchdog")
        self.stall_s = float(stall_min) * 60.0
        self.on_stall = on_stall
        self.check_s = float(check_s)
        self._halt = threading.Event()
        self.fired: Optional[str] = None        # why it fired (set before it acts)
        self.acted = threading.Event()          # set once on_stall has returned

    def limit_s(self) -> float:
        return max([self.stall_s] + [m * 60.0 for m in _ALLOW])

    def run(self) -> None:
        while not self._halt.wait(self.check_s):
            idle = time.time() - last_progress()
            if idle > self.limit_s():
                self.fired = (f"no progress for {idle / 60:.1f} min (limit {self.limit_s() / 60:.1f} min); last "
                              f"progress: {_LAST['what']!r}")
                try:
                    self.on_stall(self.fired)
                finally:
                    self.acted.set()
                    return

    def stop(self) -> None:
        self._halt.set()


_WD: List[Watchdog] = []


def _local_log(msg: str) -> None:
    try:
        BEAT_DIR.parent.mkdir(parents=True, exist_ok=True)
        with open(Path(tempfile.gettempdir()) / "jepa_v5_watchdog.log", "a", encoding="utf-8") as f:
            f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}\n")
    except OSError:
        pass


def stall_action(L, release_fn: Callable[[str], None], drive_timeout_s: float = 120.0,
                 release_timeout_s: float = 180.0) -> Callable[[str], None]:
    """What the watchdog does on a stall: log (local and Drive), RUN_SUMMARY.txt, release (bounded)."""

    def act(reason: str) -> None:
        msg = f"[watchdog] {reason}: releasing the runtime"
        print(msg, flush=True)
        _local_log(msg)

        def drive_logs():
            p = L.logs / "watchdog.log"
            L.guard(p).parent.mkdir(parents=True, exist_ok=True)
            with open(p, "a", encoding="utf-8") as f:
                f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}\n")
            SUMMARY.failure("watchdog", "session", reason)
            SUMMARY.released(f"watchdog: {reason}")
            SUMMARY.write()

        if not _bounded(drive_logs, drive_timeout_s):
            _local_log("[watchdog] Drive did not answer within the timeout; releasing without the Drive logs")
        if not _bounded(lambda: release_fn(f"watchdog: {reason}"), release_timeout_s):
            _local_log("[watchdog] release_runtime did not return; calling runtime.unassign() directly")
            try:
                from google.colab import runtime          # noqa: WPS433 (Colab only)
                runtime.unassign()
            except Exception as e:
                _local_log(f"[watchdog] runtime.unassign() failed: {e}")
    return act


def start_watchdog(stall_min: Optional[float], on_stall: Callable[[str], None], check_s: float = 30.0,
                   log=print) -> Optional[Watchdog]:
    """Start (or restart, when Setup runs again) the watchdog; stall_min None or 0 switches it off."""
    for w in _WD:
        w.stop()
    _WD.clear()
    touch("watchdog started")
    if not stall_min:
        log("watchdog off (STALL_MIN = None)")
        return None
    STATE["stall_min"] = float(stall_min)
    w = Watchdog(stall_min, on_stall, check_s)
    w.start()
    _WD.append(w)
    log(f"watchdog on: the runtime is released after {stall_min:g} min without progress (STALL_MIN)")
    return w


def json_safe(o):
    return json.loads(json.dumps(o, default=str))
