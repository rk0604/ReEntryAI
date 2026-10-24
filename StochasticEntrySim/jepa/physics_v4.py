"""
jepa/physics_v4.py

The non-learned models v4 flies MPPI through (jepa/v4_plan.md section 4.1).

    AtmTable           the simulator's tabulated atmosphere (AtmosphereModel's
                       fast table: log-density and temperature against altitude)
    aero_coefficients  the orion_cm_trim CD and L/D schedule of
                       math_3d.aero_coefficients_from_speed_altitude
    eom_step           batched torch port of math_3d.nominal_eom_step (one
                       explicit Euler step of the translational model the
                       classical predictor-corrector predicts with)
    rollout_physics    that model as an MPPI rollout: the classical predictor's
                       2 s Euler step, the realized bank moving toward each 5 s
                       command at the 20 deg/s roll-rate limit, g-load = the
                       largest aerodynamic load over the substeps of each 5 s step
    PhysicsAdapter     mppi_nominal (nominal parameters, the classical
                       controller's own model) and mppi_truth (the flight's true
                       density, L/D and mass) for fly_v4.run_lockstep_v4
    PriorAdapter       mppi_prior: v2's own physics backbone (models.integrate_eom:
                       exponential atmosphere, one ballistic coefficient, the
                       commanded bank applied at once) with zero learned
                       corrections, i.e. v2 with the learning removed
    port_check / open_loop_fidelity
                       the section 4.1.3 checks (stored in 01_checks/)

Nothing here changes the simulator: math_3d, AtmosphereModel and the guidance
are only read.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch

import constants
from jepa.data import G0, MACRO_DT_S, NOMINAL_DOMAIN, R_EARTH
from jepa.models import integrate_eom

MU = float(constants.MU)
RATE_MAX_RPS = math.radians(20.0)         # rl_env BankActuatorLimits.sigma_rate_max_rps (fly_v3.SIGMA_RATE_MAX_DPS)
PREDICTOR_DT_S = 2.0                      # control.SimpleBankGuidance.long_prediction_dt_s
SUBSTEPS_S = (2.0, 2.0, 1.0)              # Euler steps per 5 s macro step: the predictor's 2 s, cut to land on 5 s
COS_GAMMA_STOP = 0.05                     # math_3d.rollout_nominal_to_terminal stops here ("singular")
EPS = 1.0e-9                              # math_3d.nominal_eom_step's guards
# math_3d.aero_coefficients_from_speed_altitude, default_trim_schedule (FMV, CD_trim, LD_trim); port_check
# compares every step against the original, so a change there shows up as a failed check.
TRIM_SCHEDULE = ((0.0, 0.83, 0.10), (0.5, 0.85, 0.15), (0.7, 0.88, 0.20), (0.9, 1.02, 0.23), (1.1, 1.28, 0.25),
                 (1.5, 1.38, 0.28), (2.0, 1.43, 0.30), (3.0, 1.47, 0.33), (5.0, 1.50, 0.36), (8.0, 1.52, 0.40),
                 (10.0, 1.50, 0.39), (15.0, 1.52, 0.33), (20.0, 1.54, 0.30), (25.0, 1.55, 0.27), (30.0, 1.55, 0.27),
                 (35.0, 1.55, 0.28), (40.0, 1.55, 0.30))


def nominal_params() -> dict:
    """The parameters the classical guidance predicts with (evaluate.nominal_params)."""
    from jepa.evaluate import nominal_params as _np_
    p = _np_()
    if p.get("aero_model") != "orion_cm_trim":
        raise ValueError(f"physics_v4 ports the orion_cm_trim model, the config uses {p.get('aero_model')!r}")
    return p


# ============================================================================
# Atmosphere and aerodynamics
# ============================================================================
class AtmTable:
    """AtmosphereModel's fast table on a torch device (float64), with np.interp's arithmetic."""

    def __init__(self, device="cpu"):
        import AtmosphereModel as AM
        if not AM.atmosphere_table_active():
            AM.build_atmosphere_table()
        t = AM._ATM_TABLE
        kw = dict(dtype=torch.float64, device=torch.device(device))
        self.alt = torch.as_tensor(t["alt"], **kw).contiguous()
        self.log_rho = torch.as_tensor(t["log_rho"], **kw)
        self.T = torch.as_tensor(t["T"], **kw)
        self.alt_max = float(t["alt_max"])

    def _interp(self, x: torch.Tensor, fp: torch.Tensor) -> torch.Tensor:
        shp = x.shape
        x = x.reshape(-1).contiguous()
        j = (torch.searchsorted(self.alt, x, right=True) - 1).clamp(0, len(self.alt) - 2)
        x0, x1, y0, y1 = self.alt[j], self.alt[j + 1], fp[j], fp[j + 1]
        return ((y1 - y0) / (x1 - x0) * (x - x0) + y0).reshape(shp)

    def state(self, h_m: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """AtmosphereModel.get_atmosphere_state with the table active: (rho, T; NaN where undefined)."""
        h = h_m.clamp_min(0.0)
        hc = h.clamp_max(self.alt_max)
        rho = torch.exp(self._interp(hc, self.log_rho))
        T = self._interp(hc, self.T)
        space = h / 1000.0 > 1000.0
        rho = torch.where(space, torch.zeros_like(rho), rho)
        T = torch.where((T > 0.0) & ~space, T, torch.full_like(T, float("nan")))
        return rho, T


_ATM: Dict[str, AtmTable] = {}


def atm_table(device) -> AtmTable:
    key = str(torch.device(device))
    if key not in _ATM:
        _ATM[key] = AtmTable(device)
    return _ATM[key]


def aero_coefficients(V: torch.Tensor, T: torch.Tensor, lift_scale: torch.Tensor,
                      schedule: Sequence[Sequence[float]] = TRIM_SCHEDULE) -> Tuple[torch.Tensor, torch.Tensor]:
    """(CD, CL) of math_3d.aero_coefficients_from_speed_altitude for orion_cm_trim, with cpas_lift_scale."""
    sch = torch.as_tensor(sorted(schedule, key=lambda r: r[0]), dtype=V.dtype, device=V.device)
    f, cd, ld = sch[:, 0].contiguous(), sch[:, 1], sch[:, 2]
    V_kfps = V * 0.00328084
    t_ok = torch.isfinite(T) & (T > 0.0)
    mach = V / torch.sqrt(1.4 * 287.053 * torch.where(t_ok, T, torch.ones_like(T)))
    a = (V_kfps - 8.8) / (9.8 - 8.8)
    fmv = torch.where(V_kfps >= 9.8, V_kfps,
                      torch.where(~t_ok, V_kfps, torch.where(V_kfps <= 8.8, mach, (1.0 - a) * mach + a * V_kfps)))
    shp = fmv.shape
    x = fmv.reshape(-1).contiguous()
    i = (torch.searchsorted(f, x) - 1).clamp(0, len(f) - 2)      # first segment with f_i <= FMV <= f_i+1
    f0, f1 = f[i], f[i + 1]
    al = (x - f0) / (f1 - f0)
    CD = cd[i] + al * (cd[i + 1] - cd[i])
    LD = ld[i] + al * (ld[i + 1] - ld[i])
    lo, hi = x <= f[0], x >= f[-1]
    CD = torch.where(lo, cd[0], torch.where(hi, cd[-1], CD)).reshape(shp)
    LD = torch.where(lo, ld[0], torch.where(hi, ld[-1], LD)).reshape(shp)
    return CD, CD * LD * lift_scale


def eom_step(x: torch.Tensor, sigma: torch.Tensor, S_ref: torch.Tensor, mass: torch.Tensor,
             lift_scale: torch.Tensor, dt: float, atm: AtmTable) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    math_3d.nominal_eom_step for a batch: x [B, 6] = (r, phi, lam, V, gamma, chi),
    bank sigma [B], per-row reference area, mass and lift scale [B].
    Returns (x_next [B, 6], aerodynamic load [B] in g = |drag + lift| / (m g0)).
    """
    r, phi, lam, V, gam, chi = x.unbind(-1)
    alt = (r - float(constants.RADIUS_EARTH)).clamp_min(0.0)
    rho, T = atm.state(alt)
    q = 0.5 * rho * V * V
    CD, CL = aero_coefficients(V, T, lift_scale)
    drag = q * S_ref * CD
    lift = q * S_ref * CL
    g = MU / (r * r)
    cos_phi, sin_phi = torch.cos(phi), torch.sin(phi)
    cos_g, sin_g = torch.cos(gam), torch.sin(gam)
    sin_c, cos_c = torch.sin(chi), torch.cos(chi)
    safe_r = r.clamp_min(float(constants.RADIUS_EARTH) + 1.0)
    safe_V = V.abs().clamp_min(EPS)
    tiny = lambda c: torch.where(c.abs() < EPS, torch.where(c >= 0.0, torch.full_like(c, EPS),
                                                             torch.full_like(c, -EPS)), c)
    cos_phi, cos_g = tiny(cos_phi), tiny(cos_g)
    tan_phi = sin_phi / cos_phi
    L_v = lift * torch.cos(sigma)
    L_lat = -lift * torch.sin(sigma)
    r_dot = V * sin_g
    phi_dot = (V * cos_g * sin_c) / safe_r
    lam_dot = (V * cos_g * cos_c) / (safe_r * cos_phi)
    V_dot = -drag / mass - g * sin_g
    gam_dot = L_v / (mass * safe_V) + (safe_V / safe_r - g / safe_V) * cos_g
    chi_dot = L_lat / (mass * safe_V * cos_g) - (safe_V / safe_r) * cos_g * cos_c * tan_phi
    x_dot = torch.stack([r_dot, phi_dot, lam_dot, V_dot, gam_dot, chi_dot], -1)
    return x + x_dot * dt, torch.sqrt(drag * drag + lift * lift) / (mass * G0)


@torch.no_grad()
def rollout_physics(y0: torch.Tensor, sigma0: torch.Tensor, u: torch.Tensor, S_ref: torch.Tensor,
                    mass: torch.Tensor, lift_scale: torch.Tensor, h_stop: float, atm: AtmTable,
                    substeps: Sequence[float] = SUBSTEPS_S, rate: float = RATE_MAX_RPS,
                    bank_trace: Optional[list] = None):
    """
    Roll the classical predictor's model forward under planned 5 s commands.

        y0 [B, 6]   (h, phi, lam, V, gamma, chi), altitude convention (as fly_v3's flights pass it)
        sigma0 [B]  realized bank now; u [B, K] planned commands (rad)

    Explicit Euler with the predictor's 2 s step, cut to (2, 2, 1) s so every
    command starts on the 5 s grid. The bank used over a substep is the one at
    its start; after the substep it moves toward the command by at most
    rate x substep. A row is frozen once it is below h_stop or its flight path
    is near vertical (|cos gamma| < 0.05, where the predictor stops), so it
    never produces non-finite values. Stops early when every row is below
    h_stop, checked every 5 steps (as fly_v3.V1Adapter). Returns
    (y6 [B, K'+1, 6], gload [B, K']) in float64, gload = the largest
    aerodynamic load over the substeps of each 5 s step. bank_trace, if a
    list, receives the bank used in every substep (for the smoke test).
    """
    Re = float(constants.RADIUS_EARTH)
    x = y0.to(torch.float64).clone()
    x[:, 0] = x[:, 0] + Re
    sig = sigma0.to(torch.float64).clone()
    u = u.to(torch.float64)
    frozen = (x[:, 0] - Re <= h_stop) | (torch.cos(x[:, 4]).abs() < COS_GAMMA_STOP)
    ys, gls = [x.clone()], []
    for k in range(u.shape[1]):
        cmd = u[:, k]
        gmax = torch.zeros_like(sig)
        for ds in substeps:
            if bank_trace is not None:
                bank_trace.append((float(ds), sig.clone()))
            xn, gl = eom_step(x, sig, S_ref, mass, lift_scale, float(ds), atm)
            live = ~frozen
            x = torch.where(live[:, None], xn, x)
            gmax = torch.where(live, torch.maximum(gmax, gl), gmax)
            sig = sig + (cmd - sig).clamp(-rate * ds, rate * ds)
            frozen = frozen | (x[:, 0] - Re <= h_stop) | (torch.cos(x[:, 4]).abs() < COS_GAMMA_STOP)
        ys.append(x.clone())
        gls.append(gmax)
        if k % 5 == 4 and bool((x[:, 0] - Re < h_stop).all()):
            break
    y6 = torch.stack(ys, 1)
    y6[..., 0] = y6[..., 0] - Re
    return y6, torch.stack(gls, 1)


# ============================================================================
# MPPI adapters (fly_v4.run_lockstep_v4 interface)
# ============================================================================
class PhysicsAdapter:
    """
    MPPI through the classical controller's physics model. mode "nominal": the
    nominal parameters (mppi_nominal); "truth": each flight's true density, L/D
    and mass from its job's domain, applied as data.apply_domain does (mppi_truth).

    context() returns [1, E, 5] = (realized bank, roll rate, density scale,
    lift scale, mass scale) per flight; the runner repeats it per sample row and
    rollout() reads each row's parameters from s0[0]. The roll rate is carried
    for the record; the model's bank dynamics are the rate limit alone.
    """
    M = 1
    rf = 2            # one previous command (the planner's smoothness term needs it), no history model

    def __init__(self, mode: str, device="cpu"):
        assert mode in ("nominal", "truth"), mode
        self.mode, self.dev = mode, torch.device(device)
        p = nominal_params()
        self.area0, self.mass0 = float(p["ref_area_m2"]), float(p["mass_kg"])
        self.atm = atm_table(self.dev)

    def context(self, y_hists: Sequence[torch.Tensor], jobs: Sequence[dict]) -> torch.Tensor:
        rows = []
        for y, j in zip(y_hists, jobs):
            d = dict(j.get("domain") or NOMINAL_DOMAIN) if self.mode == "truth" else dict(NOMINAL_DOMAIN)
            rows.append([float(y[-1, 6]), float(y[-1, 7]), float(d["density_scale"]), float(d["lift_scale"]),
                         float(d["mass_scale"])])
        return torch.tensor(rows, dtype=torch.float64)[None]

    def rollout(self, s0: torch.Tensor, y0: torch.Tensor, uh: torch.Tensor, u: torch.Tensor,
                n_sub: int, h_stop: float):
        """Samples -> (y6 [B, 1, K'+1, 6], gload [B, 1, K']) as float32 on the CPU, like the learned adapters."""
        c = s0[0].to(self.dev, torch.float64)
        S_ref = self.area0 * c[:, 2]
        mass = self.mass0 * c[:, 4]
        y6, gl = rollout_physics(y0.to(self.dev), c[:, 0], u.to(self.dev), S_ref, mass, c[:, 3], h_stop, self.atm)
        return y6.float().cpu()[:, None], gl.float().cpu()[:, None]


class PriorAdapter:
    """
    mppi_prior: models.integrate_eom with zero corrections, the "physics prior
    only" model of PLAN.md section 3.4 and exactly what an untrained (zero-init)
    PIProberV2 integrates. Bank as in v2's prober: the commanded bank at once
    (zero learned offset). fp32, like v2's prober.
    """
    M = 1
    rf = 2

    def __init__(self, device="cpu"):
        self.dev = torch.device(device)

    def context(self, y_hists: Sequence[torch.Tensor], jobs: Sequence[dict]) -> torch.Tensor:
        return torch.zeros(1, len(y_hists), 1)

    def rollout(self, s0: torch.Tensor, y0: torch.Tensor, uh: torch.Tensor, u: torch.Tensor,
                n_sub: int, h_stop: float):
        u = u.to(self.dev).float()
        y = y0.to(self.dev).float()
        zero = torch.zeros_like(u[:, 0])
        ys, accs = [y], []
        for k in range(u.shape[1]):
            y, a = integrate_eom(y, torch.stack([zero, zero, u[:, k]], -1), MACRO_DT_S, n_sub)
            ys.append(y)
            accs.append(a)
            if k % 5 == 4 and bool((y[:, 0] < h_stop).all()):
                break
        y6 = torch.stack(ys, 1)[:, None]
        gl = (torch.stack(accs, 1).norm(dim=-1) / G0)[:, None]
        return y6.cpu(), gl.cpu()


# ============================================================================
# Checks (v4_plan.md section 4.1.3)
# ============================================================================
def port_check(n: int = 1000, seed: int = 0, dt: float = PREDICTOR_DT_S) -> Dict[str, float]:
    """
    Check 1: the torch step against math_3d.nominal_eom_step on n random states,
    banks and parameters (float64). Relative error of x_next per component
    (|a - b| / max(|b|, 1): angles in radians are compared absolutely) and of
    the increment x_next - x (relative to |increment|, floor 1e-12 |x|), plus the load.
    """
    from math_3d import nominal_aero_forces_from_state, nominal_eom_step
    from jepa.data import apply_domain
    rng = np.random.default_rng(seed)
    p0 = nominal_params()
    atm = atm_table("cpu")
    X, sig, dom, ref, gref = [], [], [], [], []
    for _ in range(n):
        h = rng.uniform(0.0, 165_000.0)
        x = [h + float(constants.RADIUS_EARTH), math.radians(rng.uniform(-40, 40)), math.radians(rng.uniform(-60, 60)),
             rng.uniform(50.0, 8_600.0), math.radians(rng.uniform(-60, 10)), math.radians(rng.uniform(0, 360))]
        s = math.radians(rng.uniform(-90, 90))
        d = {"density_scale": rng.uniform(0.7, 1.3), "lift_scale": rng.uniform(0.7, 1.3),
             "mass_scale": rng.uniform(0.9, 1.1)}
        p = apply_domain(dict(p0), d)
        out = nominal_eom_step(x=x, sigma_rad=s, params=p, dt_s=dt)
        a = nominal_aero_forces_from_state(x=x, sigma_rad=s, params=p)
        X.append(x), sig.append(s), dom.append(d), ref.append(out["x_next"])
        gref.append(math.hypot(a["drag_mag_N"], a["lift_mag_N"]) / p["mass_kg"] / G0)
    X = torch.tensor(X, dtype=torch.float64)
    S_ref = torch.tensor([float(constants.CAPSULE_REFERENCE_AREA_M2) * d["density_scale"] for d in dom], dtype=torch.float64)
    mass = torch.tensor([float(constants.CAPSULE_MASS_KG) * d["mass_scale"] for d in dom], dtype=torch.float64)
    ls = torch.tensor([d["lift_scale"] for d in dom], dtype=torch.float64)
    xn, gl = eom_step(X, torch.tensor(sig, dtype=torch.float64), S_ref, mass, ls, dt, atm)
    xn, X, ref, gref = xn.numpy(), X.numpy(), np.asarray(ref), np.asarray(gref)
    rel = np.abs(xn - ref) / np.maximum(np.abs(ref), 1.0)
    inc_ref, inc = ref - X, xn - X
    rel_inc = np.abs(inc - inc_ref) / np.maximum(np.abs(inc_ref), 1e-12 * np.maximum(np.abs(X), 1.0))
    rel_g = np.abs(gl.numpy() - gref) / np.maximum(gref, 1e-12)
    names = ("r", "phi", "lam", "V", "gamma", "chi")
    return {"n": int(n), "dt_s": float(dt), "max_rel_err": float(rel.max()), "max_rel_err_increment": float(rel_inc.max()),
            "max_rel_err_gload": float(rel_g.max()),
            **{f"max_rel_err_{k}": float(rel[:, i].max()) for i, k in enumerate(names)}}


def _crossing(t: np.ndarray, h: np.ndarray, phi: np.ndarray, lam: np.ndarray, h_cross: float):
    """First downward crossing of h_cross in 1 Hz telemetry, linearly interpolated (lat, lon in rad)."""
    below = np.nonzero(h <= h_cross)[0]
    k = int(below[0]) if len(below) else len(h) - 1
    if k == 0:
        return phi[0], lam[0]
    f = np.clip((h[k - 1] - h_cross) / max(h[k - 1] - h[k], 1e-6), 0.0, 1.0)
    return phi[k - 1] + f * (phi[k] - phi[k - 1]), lam[k - 1] + f * (lam[k] - lam[k - 1])


def open_loop_fidelity(tels: Sequence[Tuple[dict, Dict[str, np.ndarray]]], ttg_s: Sequence[float] = (300, 200, 100),
                       h_cross: float = 7600.0, device="cpu") -> List[dict]:
    """
    Check 2: from a recorded flight's state at time-to-go ttg (on its 5 s
    decision grid), with the commands it actually flew afterwards (last one held),
    the physics model's predicted drogue point against the real one, with the
    flight's true parameters (the oracle's model error) and with nominal ones
    (the classical controller's). tels: (result row, telemetry arrays) pairs.
    """
    from jepa.evaluate import crossing_point, haversine_km
    adapters = {"truth": PhysicsAdapter("truth", device), "nominal": PhysicsAdapter("nominal", device)}
    rows = []
    for res, tel in tels:
        t = tel["t_s"].astype(np.float64)
        h = tel["h_m"].astype(np.float64)
        phi, lam = np.radians(tel["phi_deg"].astype(np.float64)), np.radians(tel["lambda_deg"].astype(np.float64))
        p_true = _crossing(t, h, phi, lam, h_cross)
        for ttg in ttg_s:
            i0 = int(5 * round((t[-1] - ttg) / 5.0))
            if i0 <= 0 or i0 >= len(t) - 1:
                continue
            y0 = torch.tensor([[h[i0], phi[i0], lam[i0], tel["V_mps"][i0], math.radians(tel["gamma_deg"][i0]),
                                math.radians(tel["chi_deg"][i0])]], dtype=torch.float64)
            cmds = np.radians(tel["sigma_cmd_deg"][i0 + 1::5].astype(np.float64))
            K = int(ttg // 5) + 30
            u = torch.tensor(np.concatenate([cmds, np.repeat(cmds[-1:], max(0, K - len(cmds)))])[:K])[None]
            job = {"domain": {k[4:]: res[k] for k in ("dom_density_scale", "dom_lift_scale", "dom_mass_scale")}}
            y_hist = torch.zeros(1, 10, dtype=torch.float64)
            y_hist[0, 6] = math.radians(float(tel["sigma_deg"][i0]))
            row = {"flight_id": res["flight_id"], "case": int(res["case"]), "atmosphere": res["atmosphere"],
                   "ttg_s": float(ttg)}
            for name, ad in adapters.items():
                y6, _ = ad.rollout(ad.context([y_hist], [job]), y0, None, u, 2, h_cross)
                pp, pl, crossed, _ = crossing_point(y6[:, 0].double(), h_cross)
                row[f"err_{name}_km"] = float(haversine_km(pp, pl, torch.tensor([p_true[0]]), torch.tensor([p_true[1]]))[0])
                row[f"crossed_{name}"] = bool(crossed[0])
            rows.append(row)
    return rows
