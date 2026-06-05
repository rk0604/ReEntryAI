# Colab — RL training & evaluation

All **heavy compute** (baseline characterization, RL training, evaluation) runs
here, not on the local machine.

- **`ReEntryAI_RL_Colab.ipynb`** — the notebook. Open it in Google Colab via
  *File → Open notebook → GitHub →* `rk0604/ReEntryAI`, or:

  `https://colab.research.google.com/github/rk0604/ReEntryAI/blob/main/StochasticEntrySim/colab/ReEntryAI_RL_Colab.ipynb`

- **`build_colab_notebook.py`** — generator for the notebook. Edit this, then
  `python colab/build_colab_notebook.py` to regenerate, and commit both files.
  (The `.ipynb` is generated; don't hand-edit it.)

### What the notebook does
1. Clone the repo, install deps (`pymsis`, `gymnasium`, `stable-baselines3`,
   `wandb`), optionally mount Drive for checkpoints.
2. **Phase 3** — classical predictor-corrector over the frozen held-out set
   (`configs/heldout_testset.csv`) → `baseline_metrics.csv`.
3. **Phase 4** — train PPO on `OrionEntryEnv` (policy mode, dispersed ICs),
   logging reward curves + per-episode mission metrics to Weights & Biases.
4. **Phase 5** — run the trained policy on the **same** frozen ICs and compare
   the two metric distributions head-to-head.

Set `Runtime → GPU`. Throughput is CPU-bound (the env step), so more vCPUs help
training more than the GPU does.

### Persistence & caching (survives Colab disconnects)
Everything is written to **`MyDrive/ReEntryAI/`** on Google Drive:

```
MyDrive/ReEntryAI/
  data/         baseline_metrics.csv, policy_metrics.csv
  checkpoints/  periodic PPO checkpoints  (training auto-resumes from the latest)
  models/       final saved models
  tb/           tensorboard logs
```

Caching means heavy compute isn't repeated:
- **Baseline** — skipped if `data/baseline_metrics.csv` already exists
  (`FORCE_BASELINE=True` to recompute).
- **Training** — auto-resumes from the newest checkpoint in `checkpoints/`, and
  the W&B run (`resume="allow"`) continues the same curves.
- **In-sim** — each worker builds the atmosphere lookup table once and reuses it
  (the optimization that made runs fast in the first place).

