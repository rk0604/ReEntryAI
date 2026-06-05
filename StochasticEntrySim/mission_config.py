"""
Mission configuration loader for the entry simulation.

A single JSON file controls every parameter we routinely vary for RL
training data generation. The same config is consumed by both
run_sim.py and ReEntryAI_run2.ipynb so the two drivers always produce
identical outputs for identical configs.

JSON schema
-----------

    {
      "run_id":     "nominal_001",       // short unique identifier
      "run_name":   "Baseline mission",  // human-readable label
      "output_dir": "runs/nominal_001",  // where to write outputs (relative
                                         //   to StochasticEntrySim/)
      "seed":       42,                  // master RNG seed

      "initial_state": {
        "altitude_m":   115000.0,
        "velocity_mps": 10500.0,
        "gamma_deg":    -5.5,
        "chi_deg":      90.0,
        "phi_deg":      0.0,
        "lambda_deg":   0.0
      },

      "aerodynamics": {
        "model":         "orion_cm_trim",     // "orion_cm_trim" | "scheduled_orion_like" | "constant"
        "cd_const":      1.15,                // used only when model="constant"
        "cl_const":      0.28,
        "trim_schedule": null                 // optional override of the FMV->(CD, L/D) table
      },

      "cpas": {
        // Any field from cpas.CPASConfig is overridable here.
        // Examples:
        "enable_stochastic_failure":            true,
        "main_failure_probability_per_chute":   0.00241,
        "num_mains_operational":                3,
        "wind_east_mps":                        5.0,
        "wind_north_mps":                       2.0,
        "pendulum_target_amplitude_deg":        25.0,
        "main_deploy_alt_m":                    2400.0
      },

      "guidance": {
        "gamma_safety_threshold_deg": 30.0,
        "candidate_bank_magnitudes_deg": [0, 25, 45, 70, 90]
      },

      "interval_supervisor": {
        "heat_rate_limit_w_per_m2":  15.0e6,
        "heat_load_limit_j_per_m2":  2.5e9,
        "cos_gamma_safety":          0.2
      },

      "mission": {
        "target_phi_deg":               15.0,
        "target_lambda_deg":            15.0,
        "cos_gamma_termination_gate":   0.2
      }
    }

Every section is optional. Anything you omit falls back to the in-code
defaults of the relevant subsystem, so the smallest valid config is just
``{"run_id": "x"}``.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import cpas


THIS_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH = THIS_DIR / "configs" / "default.json"


# =============================================================================
# Loader
# =============================================================================

@dataclass
class MissionConfig:
    """In-memory view of a loaded mission config.

    The raw dict is kept around so subsystems that don't have a typed
    accessor can still reach into it.
    """
    raw: Dict[str, Any]
    source_path: Optional[Path] = None

    # ---- top-level metadata ----
    @property
    def run_id(self) -> str:
        return str(self.raw.get("run_id", "unnamed_run"))

    @property
    def run_name(self) -> str:
        return str(self.raw.get("run_name", self.run_id))

    @property
    def output_dir(self) -> Path:
        sub = self.raw.get("output_dir", f"runs/{self.run_id}")
        p = Path(sub)
        if not p.is_absolute():
            p = THIS_DIR / p
        return p

    @property
    def seed(self) -> int:
        return int(self.raw.get("seed", 42))


def load_config(path: Optional[str | Path] = None) -> MissionConfig:
    """Load a mission config from JSON.

    If ``path`` is None, falls back to (in order):
      1. The SIM_CONFIG environment variable
      2. configs/default.json next to this module
    """
    if path is None:
        env_path = os.environ.get("SIM_CONFIG", "").strip()
        if env_path:
            path = env_path
        else:
            path = DEFAULT_CONFIG_PATH

    p = Path(path)
    if not p.is_absolute() and not p.exists():
        # Try resolving relative to StochasticEntrySim/
        p_alt = THIS_DIR / p
        if p_alt.exists():
            p = p_alt
    if not p.exists():
        raise FileNotFoundError(f"Mission config not found: {path}")

    with open(p, "r", encoding="utf-8") as f:
        data = json.load(f)
    return MissionConfig(raw=data, source_path=p)


# =============================================================================
# Subsystem appliers
# =============================================================================

def initial_state_vector(cfg: MissionConfig, radius_earth_m: float = None) -> List[float]:
    """Return the 6-element [r, phi, lam, V, gamma, chi] entry-interface vector.

    radius_earth_m is accepted for backwards compatibility with the existing
    call sites, but defaults to `constants.RADIUS_EARTH` (the only value any
    caller ever passes).
    """
    if radius_earth_m is None:
        import constants  # deferred to avoid an import cycle at module load
        radius_earth_m = float(constants.RADIUS_EARTH)
    ic = cfg.raw.get("initial_state", {})
    return [
        float(radius_earth_m + float(ic.get("altitude_m", 115_000.0))),
        math.radians(float(ic.get("phi_deg", 0.0))),
        math.radians(float(ic.get("lambda_deg", 0.0))),
        float(ic.get("velocity_mps", 10_500.0)),
        math.radians(float(ic.get("gamma_deg", -5.5))),
        math.radians(float(ic.get("chi_deg", 90.0))),
    ]


# =============================================================================
# Initial-condition dispersion (Phase 2)
# =============================================================================
# The six entry-interface variables that may be dispersed for RL training and
# for building the frozen held-out test set. Values are sampled in *config
# units* (m, m/s, deg) around the nominal `initial_state`.
_DISPERSION_VARS = (
    "altitude_m", "velocity_mps", "gamma_deg", "chi_deg", "phi_deg", "lambda_deg",
)

_IC_DEFAULTS = {
    "altitude_m": 115_000.0, "velocity_mps": 10_500.0, "gamma_deg": -5.5,
    "chi_deg": 90.0, "phi_deg": 0.0, "lambda_deg": 0.0,
}


def dispersion_enabled(cfg: MissionConfig) -> bool:
    """Whether the config's dispersion block is turned on."""
    return bool(cfg.raw.get("dispersion", {}).get("enabled", False))


def dispersion_spec(cfg: MissionConfig) -> Dict[str, dict]:
    """Per-variable dispersion ranges from the config's `dispersion` block.

    Each variable maps to one of:
      * {"delta": d}          -> uniform on [nominal - d, nominal + d]
      * {"low": a, "high": b} -> uniform on [a, b] (absolute)
      * {"sigma": s, "trunc": k} -> Gaussian(nominal, s), truncated at +/- k*s
    """
    d = dict(cfg.raw.get("dispersion", {}))
    d.pop("enabled", None)
    return {k: v for k, v in d.items() if k in _DISPERSION_VARS and isinstance(v, dict)}


def _sample_one(rng, nominal: float, spec: dict) -> float:
    if "low" in spec or "high" in spec:
        lo = float(spec.get("low", nominal))
        hi = float(spec.get("high", nominal))
        return float(rng.uniform(lo, hi))
    if "sigma" in spec:
        s = float(spec["sigma"])
        k = float(spec.get("trunc", 3.0))
        z = float(rng.normal(0.0, 1.0))
        z = max(-k, min(k, z))
        return float(nominal + s * z)
    delta = float(spec.get("delta", 0.0))
    return float(rng.uniform(nominal - delta, nominal + delta))


def nominal_initial_conditions(cfg: MissionConfig) -> Dict[str, float]:
    """Return the nominal six IC values in config units (m, m/s, deg)."""
    ic = dict(cfg.raw.get("initial_state", {}))
    return {k: float(ic.get(k, _IC_DEFAULTS[k])) for k in _DISPERSION_VARS}


def initial_state_vector_from_ic(ic: Dict[str, float],
                                 radius_earth_m: float = None) -> List[float]:
    """Build the [r, phi, lam, V, gamma, chi] vector from a dict of IC values
    (config units). Used for both the nominal and any dispersed/override IC."""
    if radius_earth_m is None:
        import constants
        radius_earth_m = float(constants.RADIUS_EARTH)
    g = lambda k: float(ic.get(k, _IC_DEFAULTS[k]))
    return [
        float(radius_earth_m + g("altitude_m")),
        math.radians(g("phi_deg")),
        math.radians(g("lambda_deg")),
        g("velocity_mps"),
        math.radians(g("gamma_deg")),
        math.radians(g("chi_deg")),
    ]


def sample_dispersed_ic(cfg: MissionConfig, rng) -> Dict[str, float]:
    """Draw one dispersed IC (config units) per the config's dispersion spec.

    Variables without a spec entry keep their nominal value. Returns the drawn
    dict so the exact conditions can be logged / frozen into a test set.
    """
    nominal = nominal_initial_conditions(cfg)
    spec = dispersion_spec(cfg)
    return {
        k: (_sample_one(rng, v, spec[k]) if k in spec else v)
        for k, v in nominal.items()
    }


def apply_aero_to_params(cfg: MissionConfig, params: Dict[str, Any]) -> Dict[str, Any]:
    """Mutate and return the params dict with aero overrides applied."""
    aero = cfg.raw.get("aerodynamics", {})
    if "model" in aero:
        params["aero_model"] = str(aero["model"])
    if "cd_const" in aero:
        params["CD_const"] = float(aero["cd_const"])
    if "cl_const" in aero:
        params["CL_const"] = float(aero["cl_const"])
    if aero.get("trim_schedule") is not None:
        params["orion_trim_schedule"] = aero["trim_schedule"]
    return params


def build_cpas_config(cfg: MissionConfig) -> cpas.CPASConfig:
    """Construct a CPASConfig honoring whatever fields the config specifies.

    Two special conveniences vs the raw dataclass:
      * pendulum_target_amplitude_deg overrides pendulum_target_amplitude_rad
      * pendulum_initial_amplitude_deg also accepted (alias, legacy)
    """
    c = dict(cfg.raw.get("cpas", {}))

    # deg->rad conveniences
    if "pendulum_target_amplitude_deg" in c:
        c["pendulum_target_amplitude_rad"] = math.radians(
            float(c.pop("pendulum_target_amplitude_deg"))
        )
    if "pendulum_initial_amplitude_deg" in c:
        # treat legacy name as target amplitude
        c["pendulum_target_amplitude_rad"] = math.radians(
            float(c.pop("pendulum_initial_amplitude_deg"))
        )

    # only pass through known dataclass fields
    known = set(cpas.CPASConfig.__dataclass_fields__.keys())
    kwargs = {k: v for k, v in c.items() if k in known}
    unknown = set(c.keys()) - known
    if unknown:
        print(f"  [config] WARN: unknown cpas fields ignored: {sorted(unknown)}")
    return cpas.CPASConfig(**kwargs)


def guidance_overrides(cfg: MissionConfig) -> Dict[str, Any]:
    """Return guidance overrides (consumed manually by the main loop)."""
    g = dict(cfg.raw.get("guidance", {}))
    out: Dict[str, Any] = {}
    if "gamma_safety_threshold_deg" in g:
        out["gamma_safety_threshold_rad"] = math.radians(float(g["gamma_safety_threshold_deg"]))
    if "candidate_bank_magnitudes_deg" in g:
        out["candidate_bank_magnitudes_deg"] = list(g["candidate_bank_magnitudes_deg"])
    if "velocity_enable_mps" in g:
        out["velocity_enable_mps"] = float(g["velocity_enable_mps"])
    return out


def supervisor_overrides(cfg: MissionConfig) -> Dict[str, Any]:
    """Return interval supervisor overrides."""
    s = dict(cfg.raw.get("interval_supervisor", {}))
    out: Dict[str, Any] = {}
    if "heat_rate_limit_w_per_m2" in s:
        out["heat_rate_limit"] = float(s["heat_rate_limit_w_per_m2"])
    if "heat_load_limit_j_per_m2" in s:
        out["heat_load_limit"] = float(s["heat_load_limit_j_per_m2"])
    if "cos_gamma_safety" in s:
        out["interval_denominator_safety_cos_gamma"] = float(s["cos_gamma_safety"])
    return out


def mission_targets(cfg: MissionConfig) -> Dict[str, float]:
    """Return target latitude/longitude in radians and the termination gate."""
    m = cfg.raw.get("mission", {})
    return {
        "target_phi_rad":  math.radians(float(m.get("target_phi_deg",    15.0))),
        "target_lam_rad":  math.radians(float(m.get("target_lambda_deg", 15.0))),
        "cos_gamma_termination_gate": float(m.get("cos_gamma_termination_gate", 0.2)),
    }


# =============================================================================
# Helpful summary
# =============================================================================

def summarize(cfg: MissionConfig) -> str:
    """One-line summary for the console banner at the top of a run."""
    src = str(cfg.source_path) if cfg.source_path else "<inline>"
    return (
        f"[mission_config] {cfg.run_id} ({cfg.run_name})\n"
        f"                 source = {src}\n"
        f"                 output = {cfg.output_dir}\n"
        f"                 seed   = {cfg.seed}"
    )
