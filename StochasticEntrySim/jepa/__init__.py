"""
jepa: a SkyJEPA-style latent world model for Orion entry guidance.

Rao et al., "SkyJEPA: Learning Long-Horizon World Models for Zero-Shot
Sim-to-Real Control of Quadrotors" (arXiv 2606.23444) learns a JEPA latent
dynamics model, maps frozen latent rollouts to physical state with a
physics-inspired prober, and plans through it with MPPI. This package applies
the same recipe to bank-angle entry guidance. See jepa/PLAN.md for the design
and ReEntryAI_JEPA.ipynb for the end-to-end experiment.

Modules:
    data      dataset synthesis (random bank profiles flown through the real
              closed loop, with domain randomization) and macro-step windowing
    models    JEPA world model, SIGReg, physics-inspired and MLP probers, the
              autoregressive "predictive" baseline, training loops, rollouts
    evaluate  open-loop metrics (rollout error, footprint prediction, latent
              health) and the nominal-physics footprint predictor
    mppi      MPPI planner on the learned model and the closed-loop episode
              runner (JEPA-MPPI, predictive-MPPI, classical predictor corrector)

All modules expect StochasticEntrySim/ on sys.path, like the rest of the sim.
"""
