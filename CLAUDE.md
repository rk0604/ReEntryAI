# Working rules for this repo

## Google Colab notebooks: always release the runtime automatically

Every Colab notebook, new or updated, must free its runtime on its own, so it never sits idle on
paid compute units (a finished run once idled overnight and wasted credits).

- **Setup cell:** define `AUTO_DISCONNECT = IN_COLAB` and a `release_runtime(reason)` helper that
  first calls `google.colab.drive.flush_and_unmount()`, so nothing written to Drive is lost, and then
  calls `google.colab.runtime.unassign()`.
- **Release on failure too:** register an IPython `post_run_cell` hook that calls `release_runtime`
  when a cell raises. Run all stops at an error and would otherwise leave the runtime idling. Ignore
  `KeyboardInterrupt`, which means the user pressed stop. Unregister any previous copy of the hook, so
  rerunning setup doesn't add a second one.
- **Last cell:** call `release_runtime("all cells finished")`.
- **Before releasing:** everything worth keeping (checkpoints, results, figures) is already on Drive.
- **Notebook text:** say that the notebook disconnects itself, and that `AUTO_DISCONNECT = False` in
  the setup cell turns it off.
- **Runtime type:** request a GPU (and the cheapest GPU that fits) only for notebooks that use one.
  CPU-only notebooks open on a CPU runtime.

Reference implementation: `AUTO_DISCONNECT_CODE`, `DONE_MD`, and `DONE` in
`StochasticEntrySim/jepa/build_notebooks_v3.py`.
