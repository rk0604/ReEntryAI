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
