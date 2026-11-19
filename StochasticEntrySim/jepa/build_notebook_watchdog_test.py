"""
Build ReEntryAI_watchdog_test.ipynb: a 5-minute check, on a CPU runtime, that the v5 notebook's stall watchdog
(jepa/runguard_v5.py) really releases a Colab runtime from its background thread. Edit this builder, never the .ipynb:

    python jepa/build_notebook_watchdog_test.py          (from StochasticEntrySim/)

It writes nothing to Drive (a scratch folder on the runtime's local disk) and does not mount it. Outside Colab (the
smoke test, JEPA_WD_TEST_MODE=smoke) the release only prints, and the hanging cell checks that the watchdog fired.
"""
import os
import sys

import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from build_notebooks_v3 import AUTO_DISCONNECT_CODE, DONE, DONE_MD  # noqa: E402

NAME = "ReEntryAI_watchdog_test.ipynb"
BADGE = ("[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)]"
         "(https://colab.research.google.com/github/rk0604/ReEntryAI/blob/space-jepa5/StochasticEntrySim/" + NAME + ")")

INTRO = f"""
# Watchdog test: does a stalled notebook release its runtime?

{BADGE}

Run this **once before leaving a long v5 pass unattended**. It opens on a **CPU runtime** (about 5 minutes,
≈ 0.01 units). **Runtime → Run all**, then watch: cell 3 hangs on purpose, printing nothing. With
`STALL_MIN = 2` the watchdog (the same code the v5 notebook uses, `jepa/runguard_v5.py`) sees no progress for
2 minutes and releases the runtime from its background thread.

**Pass:** about 2–3 minutes into cell 3 the output shows `[watchdog] no progress for 2 min ... releasing the
runtime`, then Colab shows the runtime as disconnected (Runtime → "Connect" appears again). **Fail:** cell 3 is
still running 2 minutes after that line, or it ends with `WATCHDOG TEST FAILED`: then tell Claude, and do not leave
the v5 pass unattended.

Nothing is written to Drive (Drive is not even mounted). **The notebook disconnects itself** in every case: the
watchdog releases it, a failed cell releases it, and the last cell releases it; `AUTO_DISCONNECT = False` in the
setup cell turns that off (the test then only prints what it would do).
"""

SETTINGS = r"""
import os, sys
IN_COLAB = "google.colab" in sys.modules
MODE = os.environ.get("JEPA_WD_TEST_MODE", "colab" if IN_COLAB else "smoke")
STALL_MIN = 2 if MODE == "colab" else 0.05       # the v5 notebook uses 45
CHECK_S = 15 if MODE == "colab" else 1           # how often the watchdog looks
HANG_S = 600 if MODE == "colab" else 30          # cell 3 hangs this long at most
"""

SETUP = r"""
BRANCH = "space-jepa5"
if IN_COLAB:
    REPO = "/content/ReEntryAI_v5"
    if not os.path.exists(REPO):
        !git clone -q -b {BRANCH} https://github.com/rk0604/ReEntryAI.git {REPO}
    else:
        !git -C {REPO} pull -q origin {BRANCH}
    os.chdir(f"{REPO}/StochasticEntrySim")
    from google.colab import drive      # not mounted: release_runtime's Drive flush only prints that it had nothing
    SCRATCH = "/content/watchdog_test"
else:
    import tempfile
    SCRATCH = os.environ.get("JEPA_WD_TEST_ROOT", os.path.join(tempfile.gettempdir(), "jepa_watchdog_test"))
""" + AUTO_DISCONNECT_CODE + r"""
import time
sys.path.insert(0, os.getcwd())
from jepa import runguard_v5 as RG, study_v5 as S5
L = S5.LayoutV5(SCRATCH).make()          # local scratch folder: watchdog.log and RUN_SUMMARY.txt go here
RG.configure(L, {"test": "watchdog"}, 0.08, release_runtime)
WATCHDOG = RG.start_watchdog(STALL_MIN, RG.stall_action(L, release_runtime), check_s=CHECK_S)
print(f"watchdog armed: release after {STALL_MIN:g} min without progress (checked every {CHECK_S} s)")
"""

HANG = r"""
print(f"hanging now, printing nothing for up to {HANG_S} s; the watchdog should release the runtime after "
      f"{STALL_MIN:g} min ...", flush=True)
t0 = time.time()
while time.time() - t0 < HANG_S and WATCHDOG.fired is None:
    time.sleep(1)                                     # no touch(): no progress
fired_at = time.time() - t0
if WATCHDOG.fired is None:
    raise RuntimeError(f"WATCHDOG TEST FAILED: the watchdog did not fire within {HANG_S} s")
print(f"the watchdog fired after {fired_at / 60:.1f} min: {WATCHDOG.fired}", flush=True)
WATCHDOG.acted.wait(330)                              # its Drive logs and the release take a moment
if IN_COLAB and AUTO_DISCONNECT:
    time.sleep(120)                                   # the runtime should be gone long before this ends
    raise RuntimeError("WATCHDOG TEST FAILED: the watchdog fired but the runtime is still running 2 minutes later")
print("watchdog test passed (outside Colab the release only prints)")
log = (L.logs / "watchdog.log").read_text()
assert "no progress" in log and "releasing the runtime" in log, log
assert "watchdog" in (L.study / "RUN_SUMMARY.txt").read_text()
print("watchdog.log and RUN_SUMMARY.txt record the release")
"""

CELLS = [
    ("md", INTRO),
    ("md", "## 1. Settings"),
    ("code", SETTINGS),
    ("md", "## 2. Setup\n\nClone the `space-jepa5` code, the auto-disconnect hook (CLAUDE.md), and the watchdog with "
           "`STALL_MIN = 2`. Seconds."),
    ("code", SETUP),
    ("md", "## 3. A stalled cell\n\nPrints one line, then nothing. The watchdog should release the runtime about 2–3 "
           "minutes in."),
    ("code", HANG),
    ("md", DONE_MD),
    ("code", DONE),
]


def build() -> str:
    nb = nbf.v4.new_notebook()
    nb["cells"] = [nbf.v4.new_markdown_cell(t.strip("\n")) if k == "md" else nbf.v4.new_code_cell(t.strip("\n"))
                   for k, t in CELLS]
    for i, c in enumerate(nb["cells"]):
        c["id"] = f"wd-cell-{i:02d}"
    nb["metadata"] = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                      "language_info": {"name": "python"}, "colab": {"provenance": []}}
    path = os.path.join(SIM, NAME)
    nbf.validate(nb)
    nbf.write(nb, path)
    print(f"wrote {path} ({len(CELLS)} cells)")
    return path


if __name__ == "__main__":
    build()
