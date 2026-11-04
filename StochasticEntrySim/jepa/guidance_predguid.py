"""
jepa/guidance_predguid.py

PredGuid: the enhanced Apollo entry guidance of Bairstow (2006), Appendix B, with the final-phase
reference trajectory and gains designed preflight by his section 5.5 procedure (jepa/v5_plan.md section
4.1), for the `predguid` and `predguid_180` arms.

    executive   B.2: targeting, then the phase named by SELECTOR (1 initial roll / pre-entry attitude hold,
                2 Huntest with the constant-drag controller, 3 down/upcontrol (low-loft version), 4 ballistic,
                5 final), then the lateral logic. Run every 2 s, Apollo's guidance cycle (Bairstow ch. 4);
                the command is held in between
    targeting   B.9-B.17 (non-rotating Earth: WT = 0, RELVELSW irrelevant) and the density-bias estimator
                B.12, which feeds the PredGuid predictor
    PredGuid    B.56-B.59 with pc_sequencer: a numerical predictor (our nominal physics model, density times
                K_rho_est, Table B.10 stop conditions: crash, capture, escape) and the six-method corrector of
                Table B.11, in bank angle
    final       B.54-B.60: PREDANGL from the reference table, L/D = LOD + 4 (THETNM - PREDANGL)/F3, the g
                limiter at GMAX (4.5 g here, Apollo's 10 g), GONEPAST, VQUIT
    lateral     B.61-B.67; the command's sign maps Apollo's K2ROLL onto the simulator's bank convention
    design      section 5.5: the final-phase initial state from a nominal run of the executive (nominal
                parameters, the flight's entry state); the reference L/D (LOD) by root finding so the
                reference trajectory (Bairstow Eq. 5.3, constant vertical L/D, our model's drag and trim L/D)
                reaches the target at 7.6 km; the adjoint (lambda' = -A^T lambda, numerical Jacobian) backward
                from lambda(tf) = [1, 0, 0, -cot(gamma_f) RE/(RE+Hf), 0]; F1 = -HS lambda_H / D,
                F2 = lambda_gamma / (V cos gamma), F3 = lambda_LD, RTOGO; Q2, Q3, Q5, Q6 (B.5-B.8)
    |sigma|     <= 90 deg (D2): an L/D below 0 (lift down) saturates at 90 deg, counted per flight

Apollo's constants are kept in Apollo's units (ft, fps, fpss, nm); the navigation state is converted on
the way in. Planet constants follow the simulator's sphere (ATK = R_E in nm per rad, RE, GS = mu/R_E^2,
VSAT = sqrt(GS RE), WIE = KWE = 0); vehicle constants follow our capsule (LAD, LDCMINR = LAD cos 15 deg,
CD_EST, MASS_EI, S_REF); GMAX = 4.5 g; every other constant is Bairstow's Table B.4 / B.7 value.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import constants
from jepa.data import G0, R_EARTH
from jepa.guidance_fnpeg import MU, long_model, range_and_azimuth

FT = 0.3048
NM = 1852.0
LAD_OURS = 0.28                       # our trim L/D at Mach 20-40 is 0.27-0.30 (physics_v4.TRIM_SCHEDULE); Lu's CEV 0.28
PHASES = {1: "initial_roll", 2: "huntest", 3: "upcontrol", 4: "ballistic", 5: "final"}
DEC_FIELDS_PG = ("t_s", "h_m", "V_mps", "D_fpss", "selector", "LD", "cmd_deg", "THETNM", "PREDANGL", "LATANG",
                 "K2ROLL", "Krho_est", "sat_liftdown", "gonepast", "quit")


@dataclass
class PredGuidConfig:
    sigma_max_deg: float = 90.0
    cycle_s: float = 2.0                     # Apollo's guidance cycle
    # ---- Table B.4 (Bairstow's values unless noted) ----
    C1: float = 1.25
    C16: float = 0.01
    C17: float = 0.001
    C18: float = 500.0
    C20: float = 175.0
    CHOOK: float = 0.25
    CH1: float = 0.75
    gmax_g: float = 4.5                      # GMAX: Apollo 322 fpss (10 g); the v5 load limit here
    KA: float = 64.4
    KDMIN: float = 0.5
    KLAT: float = 0.0125
    K44: float = 44389312.0
    LAD: float = LAD_OURS
    LATBIAS: float = 0.4
    VFINAL: float = 29000.0                  # Bairstow 25,000 fps; raised so LEO entries fly initial roll -> final (Putnam)
    VLMIN: float = 18000.0
    VRCONTRL: float = 700.0
    VCORLIM: float = 1000.0
    TOL: float = 25.0
    VQUIT: float = 1000.0
    HS: float = 28500.0
    LEWD1: float = 0.1
    LEWD2: float = 0.2
    D2: float = 175.0
    D0: float = 130.0
    Q7F: float = 6.0
    Q19: float = 0.2
    ASENSIBLE: float = 1.6
    overcontrol: float = 4.0                 # the "4" of B.57
    upcontrol: str = "low_loft"              # "low_loft" (B.50-B.53) or "high_loft" (B.54-B.55)
    # ---- Table B.7 ----
    APOGEE_EPSILON_TIGHTER: float = 40.0
    APOGEE_EPSILON_TIGHT: float = 25.0
    ALT_MAX_FT: float = 600000.0
    CD_EST_INITIAL: float = 1.55             # our trim CD above Mach 25
    DELTA_T_PRED: float = 10.0
    K_RHO_FILTER_GAIN: float = 0.05
    LIFT_INC_CAPTURE_DEG: float = -3.0
    MAX_NUMBER_RUNS: int = 10
    VI_LOOSE_APOGEE: float = 26600.0
    INFINITY: float = 99999999.0
    # ---- design (section 5.5) ----
    design_dt_s: float = 1.0
    hf_m: float = 7600.0
    lod_fixed_frac: Optional[float] = None   # None: LOD by root finding to the target; else LOD = frac x LAD (Bairstow 0.6)
    gamma_safe_deg: Optional[float] = 30.0   # |gamma| beyond it: wings level (the classical controller's override,
                                             # needed by the simulator's |cos gamma| < 0.05 gate); in flight, in the
                                             # reference trajectory and in the PredGuid predictor

    def to_dict(self) -> dict:
        return asdict(self)


def planet() -> dict:
    """The planet constants on the simulator's sphere (Apollo units)."""
    RE = R_EARTH / FT
    GS = MU / R_EARTH ** 2 / FT
    return {"ATK": R_EARTH / NM, "RE": RE, "GS": GS, "VSAT": math.sqrt(GS * RE)}


# ============================================================================
# Reference trajectory and gains (section 5.5)
# ============================================================================
def _ref_deriv(m, X, V, gam, H, p, lad, gs=None):
    """Bairstow Eq. (5.3) with our model: constant vertical L/D p = cos(sigma) LAD, the exact ground-range rate."""
    if gs is not None:
        p = np.where(gam < -gs, lad, p)
    D0, L0, _ = m.aero(H, V)
    r = R_EARTH + H
    g = MU / (r * r)
    Vs = np.maximum(V, 1e-6)
    cg = np.cos(gam)
    return (V * cg * R_EARTH / r, -D0 - g * np.sin(gam), (p / lad) * L0 / Vs + (Vs / r - g / Vs) * cg,
            V * np.sin(gam))


def fly_reference(x0: Sequence[float], p: np.ndarray, lad: float, hf: float, dt: float = 1.0,
                  t_max: float = 3000.0, record: bool = False, gs: Optional[float] = None):
    """RK4 from x0 = (X, V, gamma, H) at constant vertical L/D p (vector). Returns range at hf (m) and the path."""
    m = long_model()
    p = np.atleast_1d(np.asarray(p, dtype=np.float64))
    B = len(p)
    X, V, gam, H = (np.full(B, float(v)) for v in x0)
    Xf = np.full(B, np.nan)
    live = np.ones(B, dtype=bool)
    path = []
    for k in range(int(t_max / dt)):
        if record:
            path.append(np.stack([np.full(B, k * dt), X, V, gam, H], -1))
        k1 = _ref_deriv(m, X, V, gam, H, p, lad, gs)
        k2 = _ref_deriv(m, X + .5 * dt * k1[0], V + .5 * dt * k1[1], gam + .5 * dt * k1[2], H + .5 * dt * k1[3], p, lad, gs)
        k3 = _ref_deriv(m, X + .5 * dt * k2[0], V + .5 * dt * k2[1], gam + .5 * dt * k2[2], H + .5 * dt * k2[3], p, lad, gs)
        k4 = _ref_deriv(m, X + dt * k3[0], V + dt * k3[1], gam + dt * k3[2], H + dt * k3[3], p, lad, gs)
        Xn = X + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        Vn = V + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        Gn = gam + dt / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        Hn = H + dt / 6 * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3])
        cross = live & (Hn <= hf)
        if cross.any():
            f = (H[cross] - hf) / np.maximum(H[cross] - Hn[cross], 1e-9)
            Xf[cross] = X[cross] + f * (Xn[cross] - X[cross])
            live &= ~cross
        X, V, gam, H = np.where(live, Xn, X), np.where(live, Vn, V), np.where(live, Gn, gam), np.where(live, Hn, H)
        if record and not live.any():
            path.append(np.stack([np.full(B, (k + 1) * dt), Xn, Vn, Gn, Hn], -1))
        if not live.any() or (H > 200_000.0).all():
            break
    return Xf, (np.stack(path, 1) if record else None)


@dataclass
class Design:
    ftable: Dict[str, np.ndarray]            # VREF ascending (fps) and the other columns (Apollo units)
    LOD: float
    c_ref: float
    Q2: float
    Q3: float
    Q5: float
    Q6: float
    ic_ref: dict                             # final-phase initial state of the design run (SI)
    range_ref_km: float
    saturated: bool
    t_final_start_s: float
    phases: List[Tuple[float, int]]

    def summary(self) -> dict:
        return {"LOD": self.LOD, "c_ref": self.c_ref, "Q2": self.Q2, "Q3": self.Q3, "Q5": self.Q5, "Q6": self.Q6,
                "range_ref_km": self.range_ref_km, "design_saturated": self.saturated,
                "t_final_start_s": self.t_final_start_s, "V_final_start_fps": float(self.ftable["VREF"][-1])}


def _jacobian(m, x, p, lad, gs=None):
    """Numerical A = d f / d (X, V, gamma, H, p) of the reference dynamics at one point (central differences)."""
    X, V, gam, H = x
    steps = (1.0, max(1e-3 * abs(V), 1e-3), 1e-6, 1.0, 1e-5)
    base = np.array([X, V, gam, H, p], dtype=np.float64)
    A = np.zeros((5, 5))
    pts = []
    for j in range(5):
        for sgn in (1.0, -1.0):
            b = base.copy()
            b[j] += sgn * steps[j]
            pts.append(b)
    P = np.asarray(pts)
    f = np.stack(_ref_deriv(m, P[:, 0], P[:, 1], P[:, 2], P[:, 3], P[:, 4], lad, gs), -1)     # [10, 4]
    for j in range(5):
        A[:4, j] = (f[2 * j] - f[2 * j + 1]) / (2 * steps[j])
    return A


def design_final_phase(ic: dict, target_phi: float, target_lam: float, cfg: PredGuidConfig,
                       sigma_max_deg: Optional[float] = None) -> Design:
    """
    Section 5.5 for one mission: nominal executive run to the final-phase start, LOD by root finding (or a fixed
    fraction of LAD), the reference trajectory, the adjoint gains, and the FTABLE. Nominal parameters only.
    """
    m = long_model()
    lad = cfg.LAD
    smax = math.radians(sigma_max_deg if sigma_max_deg is not None else cfg.sigma_max_deg)
    gs = math.radians(cfg.gamma_safe_deg) if cfg.gamma_safe_deg is not None else None
    # ---- 1. final-phase initial state from a nominal run of the executive (longitudinal) ----
    st = _dry_run_to_final(ic, target_phi, target_lam, cfg)
    x0 = (0.0, st["V"], st["gamma"], st["h"])
    s_go_m = st["s_go"] * R_EARTH
    # ---- 2. LOD ----
    c_lo = max(math.cos(smax), -1.0)
    if cfg.lod_fixed_frac is not None:
        c, sat = float(cfg.lod_fixed_frac), False
    else:
        grid = np.linspace(c_lo, 1.0, 41)
        Xf, _ = fly_reference(x0, grid * lad, lad, cfg.hf_m, cfg.design_dt_s, gs=gs)
        ok = np.isfinite(Xf)
        err = np.where(ok, Xf - s_go_m, np.nan)
        sat = False
        idx = np.nonzero(ok[:-1] & ok[1:] & (np.sign(err[:-1]) != np.sign(err[1:])))[0]
        if len(idx):
            a, b = grid[idx[0]], grid[idx[0] + 1]
            ea = err[idx[0]]
            for _ in range(40):
                mid = 0.5 * (a + b)
                em = float(fly_reference(x0, [mid * lad], lad, cfg.hf_m, cfg.design_dt_s, gs=gs)[0][0] - s_go_m)
                if np.sign(em) == np.sign(ea):
                    a, ea = mid, em
                else:
                    b = mid
                if b - a < 1e-7:
                    break
            c = 0.5 * (a + b)
        else:
            j = int(np.nanargmin(np.abs(err))) if ok.any() else len(grid) - 1
            c, sat = float(grid[j]), True
    p = c * lad
    Xf, path = fly_reference(x0, [p], lad, cfg.hf_m, cfg.design_dt_s, record=True, gs=gs)
    path = path[0]
    keep = np.isfinite(path[:, 1])
    path = path[keep]
    # truncate after the hf crossing (last recorded point is the first one below hf)
    below = np.nonzero(path[:, 4] <= cfg.hf_m)[0]
    if len(below):
        path = path[:below[0] + 1]
    t, X, V, gam, H = path.T
    Xend = float(Xf[0]) if np.isfinite(Xf[0]) else float(X[-1])
    # ---- 3. adjoint, backward, with the numerical Jacobian at every point ----
    r_f = R_EARTH + cfg.hf_m
    lam = np.array([1.0, 0.0, 0.0, -math.cos(gam[-1]) / math.sin(gam[-1]) * R_EARTH / r_f, 0.0])
    lams = np.zeros((len(t), 5))
    lams[-1] = lam
    As = [_jacobian(m, (X[k], V[k], gam[k], H[k]), p, lad, gs) for k in range(len(t))]
    for k in range(len(t) - 1, 0, -1):
        h = t[k] - t[k - 1]
        A1, A0 = As[k], As[k - 1]
        Am = 0.5 * (A0 + A1)
        f = lambda A, l: A.T @ l                  # d lambda / d(-t) = A^T lambda
        k1 = f(A1, lam)
        k2 = f(Am, lam + 0.5 * h * k1)
        k3 = f(Am, lam + 0.5 * h * k2)
        k4 = f(A0, lam + h * k3)
        lam = lam + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        lams[k - 1] = lam
    D0, L0, _ = m.aero(H, V)
    Dtot = np.hypot(D0, L0)
    conv = FT / NM                                # m per (m/s or m/s^2) -> nm per (fps or fpss)
    F1 = -(cfg.HS * FT) * lams[:, 3] / np.maximum(Dtot, 1e-9) * conv
    F2 = lams[:, 2] / np.maximum(V * np.cos(gam), 1e-6) * conv
    F3 = lams[:, 4] / NM
    RTOGO = (Xend - X) / NM
    tab = {"VREF": V / FT, "RDOTREF": V * np.sin(gam) / FT, "DREFR": Dtot / FT, "F2": F2, "F1": F1,
           "RTOGO": RTOGO, "F3": F3}
    # V decreases along the reference until the terminal descent; keep the strictly decreasing part, ascending
    vmin = np.minimum.accumulate(V)
    mono = np.concatenate([[True], V[1:] < vmin[:-1]])
    order = np.nonzero(mono)[0][::-1]
    ftable = {k: np.asarray(v)[order].astype(np.float64) for k, v in tab.items()}
    Q3 = float(lams[0, 1] * conv)
    Q5 = float(lams[0, 2] / NM)
    Q2 = float(RTOGO[0] - Q3 * V[0] / FT)
    Q6 = float(-gam[0])
    return Design(ftable=ftable, LOD=float(p), c_ref=float(c), Q2=Q2, Q3=Q3, Q5=Q5, Q6=Q6,
                  ic_ref={"h": float(H[0]), "V": float(V[0]), "gamma": float(gam[0]), "s_go_km": s_go_m / 1000.0},
                  range_ref_km=(Xend - X[0]) / 1000.0, saturated=bool(sat), t_final_start_s=float(st["t"]),
                  phases=st["phases"])


def _dry_run_to_final(ic: dict, target_phi: float, target_lam: float, cfg: PredGuidConfig) -> dict:
    """Fly the executive on the nominal longitudinal model from the entry state until SELECTOR = 5."""
    m = long_model()
    s_go, _ = range_and_azimuth(math.radians(ic["phi_deg"]), math.radians(ic["lambda_deg"]), target_phi, target_lam)
    h, V, gam, s = float(ic["altitude_m"]), float(ic["velocity_mps"]), math.radians(ic["gamma_deg"]), s_go
    pg = PredGuidState(cfg)
    sig = 0.0
    t, dt = 0.0, 1.0
    rate = math.radians(20.0)
    while t < 3000.0:
        if int(round(t)) % int(cfg.cycle_s) == 0:
            D0, L0, _ = m.aero(np.array([h]), np.array([V]))
            nav = {"V": V / FT, "RDOT": V * math.sin(gam) / FT, "D": float(np.hypot(D0, L0)[0]) / FT,
                   "THETNM": s * planet()["ATK"], "LATANG": 0.0, "gonepast": s < 0.0, "alt": h / FT,
                   "long": (R_EARTH + h, V, gam, s)}
            ld = executive(pg, nav, cfg, lateral=False)
            if pg.SELECTOR == 5:
                return {"h": h, "V": V, "gamma": gam, "s_go": s, "t": t, "phases": pg.phase_log}
            cmd = math.acos(max(-1.0, min(1.0, ld / cfg.LAD)))
            cmd = min(cmd, math.radians(cfg.sigma_max_deg))
        sig = sig + max(-rate * dt, min(rate * dt, cmd - sig))
        x = np.array([R_EARTH + h, V, gam, s])
        D0, L0, _ = m.aero(np.array([h]), np.array([V]))
        r = R_EARTH + h
        g = MU / r ** 2
        cs = math.cos(sig)
        # one RK4 step of the longitudinal model (as FNPEG's predictor)
        def fd(xx):
            hh = xx[0] - R_EARTH
            d0, l0, _ = m.aero(np.array([hh]), np.array([xx[1]]))
            gg = MU / xx[0] ** 2
            return np.array([xx[1] * math.sin(xx[2]), -d0[0] - gg * math.sin(xx[2]),
                             l0[0] * cs / max(xx[1], 1e-6) + (xx[1] / xx[0] - gg / max(xx[1], 1e-6)) * math.cos(xx[2]),
                             -xx[1] * math.cos(xx[2]) / xx[0]])
        k1 = fd(x); k2 = fd(x + .5 * dt * k1); k3 = fd(x + .5 * dt * k2); k4 = fd(x + dt * k3)
        x = x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        h, V, gam, s = x[0] - R_EARTH, x[1], x[2], x[3]
        t += dt
        if h <= cfg.hf_m:
            break
    raise RuntimeError(f"PredGuid design: the nominal run never reached the final phase (t={t:.0f}s, h={h:.0f} m)")


# ============================================================================
# The executive (B.2-B.10), one flight
# ============================================================================
@dataclass
class PredGuidState:
    cfg: PredGuidConfig
    design: Optional[Design] = None
    SELECTOR: int = 1
    INRLSW: int = 0
    HUNTIND: int = 0
    HIND: int = 0
    EGSW: int = 0
    GONEPAST: int = 0
    PIND: int = 0
    QUITIND: int = 0
    K1ROLL: int = 0
    K2ROLL: int = 0
    LD: float = 0.0
    A0: float = 0.0
    A1: float = 0.0
    V1: float = 0.0
    V1OLD: float = 0.0
    VL: float = 0.0
    Q7: float = 0.0
    DIFF: float = 0.0
    DIFFOLD: float = 0.0
    GAMMAL: float = 0.0
    GAMMAL1: float = 0.0
    ALP: float = 0.0
    FACT1: float = 0.0
    FACT2: float = 0.0
    VCORR: float = 0.0
    PREDANGL: float = float("nan")
    Krho_est: float = 1.0
    pg_last_bank: Optional[float] = None
    phase_log: List[Tuple[float, int]] = field(default_factory=list)
    t: float = 0.0
    # flight bookkeeping (not App. B)
    target_phi: float = 0.0
    target_lam: float = 0.0
    last_cmd: float = 0.0
    decided: bool = False
    n_cycles: int = 0
    n_sat_liftdown: int = 0
    dec: List[tuple] = field(default_factory=list)


def _interp(tab: Dict[str, np.ndarray], key: str, V: float) -> float:
    return float(np.interp(V, tab["VREF"], tab[key]))


def _final_range(pg: PredGuidState, V: float, gammal: float, q7: float, cfg: PredGuidConfig) -> float:
    """ASPF (B.35/B.36 = B.56/B.57): the final-phase range from (VL, GAMMAL, Q7)."""
    d = pg.design
    if d is None:
        raise RuntimeError("PredGuid: Huntest/PredGuid needs the final-phase design (super-orbital entry: not built)")
    if V < float(d.ftable["VREF"][-1]):
        tab = d.ftable
        return (_interp(tab, "RTOGO", V) + _interp(tab, "F2", V) * (-gammal * V - _interp(tab, "RDOTREF", V))
                + _interp(tab, "F1", V) * (q7 - _interp(tab, "DREFR", V)))
    return d.Q2 + d.Q3 * V + d.Q5 * (d.Q6 - gammal)


def executive(pg: PredGuidState, nav: dict, cfg: PredGuidConfig, lateral: bool = True) -> float:
    """
    One guidance cycle. nav (Apollo units): V, RDOT, D (total aerodynamic acceleration), THETNM, LATANG, alt,
    gonepast (target behind the horizontal velocity), long = (r, V, gamma, s) in SI for the PredGuid predictor.
    Returns LD (the commanded vertical L/D); pg.K2ROLL / K1ROLL carry the bank's sign.
    """
    P = planet()
    V, RDOT, D = nav["V"], nav["RDOT"], nav["D"]
    VSQ = V * V / (P["VSAT"] ** 2)
    LEQ = (VSQ - 1.0) * P["GS"]
    THETNM, LATANG = nav["THETNM"], nav["LATANG"]
    # density-bias estimator (B.12)
    if D > 0.0 and nav.get("alt") is not None:
        m = long_model()
        rho_std_D, _, _ = m.aero(np.array([nav["alt"] * FT]), np.array([V * FT]))
        # rho = 2 M D / (V^2 S CD sqrt(1 + LAD^2)) over the standard density, i.e. D / D_std at CD_EST
        cd_tab = rho_std_D[0] / max(0.5 * (V * FT) ** 2 * m.S / m.m, 1e-30)      # rho_std x CD_tab
        d_std = 0.5 * (V * FT) ** 2 * m.S / m.m * cd_tab * math.sqrt(1 + cfg.LAD ** 2) / FT
        if d_std > 0.0 and D >= cfg.ASENSIBLE:
            k_now = D / d_std * (m.S * 0 + 1.0) * (cd_tab / max(cd_tab, 1e-30))
            pg.Krho_est = cfg.K_RHO_FILTER_GAIN * k_now + (1.0 - cfg.K_RHO_FILTER_GAIN) * pg.Krho_est
    gmax = cfg.gmax_g * G0 / FT
    nxt = pg.SELECTOR
    for _guard in range(8):                 # phase transitions inside one cycle ("go to")
        pg.SELECTOR = nxt
        if nxt == 1:
            nxt = _initial_roll(pg, V, RDOT, D, cfg)
        elif nxt == 2:
            nxt = _huntest(pg, V, RDOT, D, LEQ, THETNM, cfg, P)
        elif nxt == 3:
            nxt = _upcontrol(pg, nav, V, RDOT, D, LEQ, THETNM, cfg, P)
        elif nxt == 4:
            nxt = _ballistic(pg, nav, V, D, THETNM, cfg, P)
        elif nxt == 5:
            if pg.design is None:           # the preflight dry run stops where the final phase starts
                pg.SELECTOR = 5
                nxt = 0
                break
            nxt = _final(pg, nav, V, RDOT, D, LEQ, THETNM, gmax, cfg)
        if nxt in (0, -1):
            break
    if not pg.phase_log or pg.phase_log[-1][1] != pg.SELECTOR:
        pg.phase_log.append((pg.t, pg.SELECTOR))
    if nxt == -1:                           # QUIT: B.54 sends ROLLC directly to output
        return pg.LD
    if lateral:
        _lateral(pg, V, VSQ, LATANG, P, cfg)
    else:
        pg.LD = max(-cfg.LAD, min(cfg.LAD, pg.LD))
    return pg.LD


def _initial_roll(pg, V, RDOT, D, cfg) -> int:
    if pg.INRLSW == 0:
        if D < cfg.ASENSIBLE:
            return 0
        pg.INRLSW = 1
        if V < cfg.VFINAL:
            pg.LD = cfg.LAD
            return_sel = 4
            pg.SELECTOR = return_sel
            return 0                        # B.18: "skip to step 5" (lateral logic); ballistic next cycle
        pg.LD = -cfg.LAD if V > cfg.VFINAL - cfg.K44 * (RDOT / V) ** 3 else cfg.LAD
        return 0
    if D > cfg.KA:
        pg.LD = cfg.LAD
    if RDOT < -cfg.VRCONTRL:
        return 0
    return 2


def _huntest(pg, V, RDOT, D, LEQ, THETNM, cfg, P) -> int:
    LEWD = cfg.LEWD1 if D > cfg.D2 else cfg.LEWD2
    if RDOT < 0.0:
        V1 = V + RDOT / LEWD
        A0 = (V1 / V) ** 2 * (D + RDOT ** 2 / (2 * cfg.HS * LEWD))
        A1 = D
    else:
        V1 = V + RDOT / cfg.LAD
        A0 = (V1 / V) ** 2 * (D + RDOT ** 2 / (2 * cfg.HS * cfg.LAD))
        A1 = A0
    pg.V1, pg.A0, pg.A1 = V1, A0, A1
    if pg.HUNTIND == 0:
        pg.HUNTIND, pg.DIFFOLD, pg.V1OLD, pg.Q7 = 1, 0.0, V1 + cfg.C18, cfg.Q7F
    pg.HIND = 0
    for _it in range(25):
        ALP = 2 * cfg.C1 * pg.A0 * cfg.HS / (LEWD * pg.V1 ** 2)
        FACT1 = pg.V1 / (1 - ALP)
        FACT2 = ALP * (ALP - 1) / pg.A0
        VL = FACT1 * (1 - math.sqrt(max(FACT2 * pg.Q7 + ALP, 0.0)))
        pg.ALP, pg.FACT1, pg.FACT2, pg.VL = ALP, FACT1, FACT2, VL
        if VL < cfg.VLMIN:
            pg.EGSW = 1
            return 5
        if VL > P["VSAT"]:
            return _constant_drag(pg, V, RDOT, D, LEQ, cfg)
        VS1 = min(pg.V1, P["VSAT"])
        DVL = VS1 - VL
        DHOOK = ((1 - VS1 / FACT1) ** 2 - ALP) / FACT2
        AHOOK = cfg.CHOOK * (DHOOK / pg.Q7 - 1) / DVL if DVL != 0 else 0.0
        GAMMAL1 = LEWD * (pg.V1 - VL) / VL
        GAMMAL = GAMMAL1 - (cfg.CH1 * P["GS"] * DVL ** 2 * (1 + AHOOK * DVL)) / (DHOOK * VL ** 2)
        if GAMMAL < 0.0:
            VL = VL + GAMMAL * VL / (LEWD - (3 * AHOOK * DVL ** 2 + 2 * DVL) * (cfg.CH1 * P["GS"]) / (DHOOK * VL ** 2))
            pg.Q7 = ((1 - VL / FACT1) ** 2 - ALP) / FACT2
            GAMMAL = 0.0
            pg.VL = VL
        GAMMAL1 = GAMMAL1 * (1 - cfg.Q19) + cfg.Q19 * GAMMAL
        pg.GAMMAL, pg.GAMMAL1 = GAMMAL, GAMMAL1
        VBARS = VL ** 2 / P["VSAT"] ** 2
        COSG = 1 - GAMMAL ** 2 / 2
        E = math.sqrt(max(1 + (VBARS - 2) * COSG ** 2 * VBARS, 1e-12))
        ASKEP = 2 * P["ATK"] * math.asin(max(-1.0, min(1.0, VBARS * COSG * GAMMAL / E)))
        ASPF = _final_range(pg, VL, GAMMAL, pg.Q7, cfg)
        ASPUP = P["ATK"] / P["RE"] * (cfg.HS / max(GAMMAL1, 1e-9)) * math.log(max(pg.A0 * VL ** 2 / (pg.Q7 * pg.V1 ** 2), 1e-12))
        ASPDWN = -RDOT * V * P["ATK"] / (pg.A0 * cfg.LAD * P["RE"])
        ASP = ASKEP + ASPF + ASPUP + ASPDWN
        pg.DIFF = THETNM - ASP
        if abs(pg.DIFF) < cfg.TOL:
            return 3
        if pg.HIND == 0:
            if pg.DIFF < 0.0:
                pg.DIFFOLD, pg.V1OLD = pg.DIFF, pg.V1
                return _constant_drag(pg, V, RDOT, D, LEQ, cfg)
            pg.VCORR = pg.V1 - pg.V1OLD
        den = pg.DIFFOLD - pg.DIFF
        pg.VCORR = pg.VCORR * pg.DIFF / den if den != 0 else cfg.VCORLIM
        pg.VCORR = min(pg.VCORR, cfg.VCORLIM)
        if P["VSAT"] - VL < pg.VCORR:
            pg.VCORR /= 2.0
        pg.V1 += pg.VCORR
        pg.HIND, pg.DIFFOLD = 1, pg.DIFF
    return _constant_drag(pg, V, RDOT, D, LEQ, cfg)


def _constant_drag(pg, V, RDOT, D, LEQ, cfg) -> int:
    pg.LD = -LEQ / cfg.D0 + cfg.C16 * (D - cfg.D0) - cfg.C17 * (RDOT + 2 * cfg.HS * cfg.D0 / V)
    _negtest(pg, D, cfg)
    return 0


def _negtest(pg, D, cfg) -> None:
    if pg.LD < 0.0 and D > cfg.C20:
        pg.LD = 0.0


def _upcontrol(pg, nav, V, RDOT, D, LEQ, THETNM, cfg, P) -> int:
    if cfg.upcontrol == "low_loft":
        if V > pg.V1:
            RDTR = cfg.LAD * (pg.V1 - V)
            DR = (V / pg.V1) ** 2 * pg.A0 - RDTR ** 2 / (2 * cfg.HS * cfg.LAD)
            pg.LD = cfg.LAD + cfg.C16 * (D - DR) - cfg.C17 * (RDOT - RDTR)
            return 0
    if V < pg.VL + cfg.C18 and RDOT < 0.0:
        pg.EGSW = 1
        return 5
    if D < pg.Q7:
        return 4
    if cfg.upcontrol == "low_loft" and D > pg.A0:
        pg.LD = cfg.LAD
        return 0
    _predguid(pg, nav, V, THETNM, cfg, P)
    _negtest(pg, D, cfg)
    return 0


def _ballistic(pg, nav, V, D, THETNM, cfg, P) -> int:
    if D > pg.Q7 + cfg.KDMIN:
        pg.EGSW = 1
        return 5
    _predguid(pg, nav, V, THETNM, cfg, P)
    return 0


def _final(pg, nav, V, RDOT, D, LEQ, THETNM, gmax, cfg) -> int:
    if V < cfg.VQUIT:
        pg.LD = 0.0
        pg.QUITIND = 1
        return -1
    if pg.GONEPAST or nav.get("gonepast"):
        pg.GONEPAST = 1
        pg.LD = -cfg.LAD
    else:
        tab = pg.design.ftable
        pg.PREDANGL = (_interp(tab, "RTOGO", V) + _interp(tab, "F2", V) * (RDOT - _interp(tab, "RDOTREF", V))
                       + _interp(tab, "F1", V) * (D - _interp(tab, "DREFR", V)))
        F3 = _interp(tab, "F3", V)
        pg.LD = pg.design.LOD + cfg.overcontrol * (THETNM - pg.PREDANGL) / (F3 if abs(F3) > 1e-9 else 1e-9)
    if D < gmax / 2.0:
        return 0
    if D > gmax:
        pg.LD = cfg.LAD
        return 0
    X = math.sqrt(max(2 * cfg.HS * (gmax - D) * (LEQ / gmax + cfg.LAD) + (2 * cfg.HS * gmax / V) ** 2, 0.0))
    if RDOT < -X:
        pg.LD = cfg.LAD
    return 0


def _lateral(pg, V, VSQ, LATANG, P, cfg) -> None:
    if not pg.GONEPAST:
        Y = cfg.KLAT * VSQ + cfg.LATBIAS / P["ATK"]
        LDCMINR = cfg.LAD * math.cos(math.radians(15.0))
        skip = False
        if abs(pg.LD) > LDCMINR:
            Y = Y / 2.0
            if pg.K2ROLL * LATANG > 0:
                pg.LD = LDCMINR * math.copysign(1.0, pg.LD)
                skip = True
        if not skip and pg.K2ROLL * LATANG > Y:
            pg.K2ROLL = -pg.K2ROLL
            if pg.LD < 0:
                pg.K1ROLL = pg.K1ROLL - pg.K2ROLL
    if abs(pg.LD) > cfg.LAD:
        pg.LD = math.copysign(cfg.LAD, pg.LD)


# ============================================================================
# PredGuid (B.7): pc_sequencer with the skip predictor and the six-method corrector
# ============================================================================
def _pg_predict(long0, bank: float, krho: float, pg_q7: float, ind_ini: int, cfg: PredGuidConfig):
    """Table B.10: (range nm or +-INF, end speed fps, end rdot fps). Constant bank, density x krho."""
    m = long_model()
    r, V, gam, s0 = long0
    s = s0
    dt = cfg.DELTA_T_PRED
    ind = ind_ini
    cs0 = math.cos(bank)
    gs = math.radians(cfg.gamma_safe_deg) if cfg.gamma_safe_deg is not None else None
    P = planet()
    for _ in range(int(4000 / dt)):
        h = r - R_EARTH
        if h < 0:
            return -cfg.INFINITY, V / FT, V * math.sin(gam) / FT
        if h / FT > cfg.ALT_MAX_FT:
            return cfg.INFINITY, V / FT, V * math.sin(gam) / FT
        D0, L0, _ = m.aero(np.array([h]), np.array([V]))
        amag = krho * float(np.hypot(D0, L0)[0]) / FT
        if ind == 0 and V * math.sin(gam) >= 0:
            ind = 1
        if ind == 1 and amag > pg_q7 and V * math.sin(gam) < 0:
            return (s0 - s) * P["ATK"], V / FT, V * math.sin(gam) / FT

        cs = 1.0 if (gs is not None and gam < -gs) else cs0

        def fd(x):
            d0, l0, _ = m.aero(np.array([x[0] - R_EARTH]), np.array([x[1]]))
            gg = MU / x[0] ** 2
            vv = max(x[1], 1e-6)
            return np.array([x[1] * math.sin(x[2]), -krho * d0[0] - gg * math.sin(x[2]),
                             krho * l0[0] * cs / vv + (vv / x[0] - gg / vv) * math.cos(x[2]),
                             -x[1] * math.cos(x[2]) / x[0]])
        x = np.array([r, V, gam, s])
        k1 = fd(x); k2 = fd(x + .5 * dt * k1); k3 = fd(x + .5 * dt * k2); k4 = fd(x + dt * k3)
        r, V, gam, s = x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return cfg.INFINITY, V / FT, V * math.sin(gam) / FT


def _predguid(pg, nav, V, THETNM, cfg, P) -> None:
    gammal = pg.GAMMAL
    aspf = _final_range(pg, pg.VL if pg.VL > 0 else V, gammal, pg.Q7, cfg)
    target = THETNM - aspf
    if not (target > 100.0 or pg.PIND == 0):
        return
    lift_inc = math.radians(abs(cfg.LIFT_INC_CAPTURE_DEG))
    runs = cfg.MAX_NUMBER_RUNS
    smax = math.radians(cfg.sigma_max_deg)
    if pg.PIND == 0:
        pg.PIND = 1
        bank = math.pi / 2.0                       # CPhi_Desired = 0
        lift_inc, runs = math.radians(10.0), 20
    else:
        bank = math.acos(max(-1.0, min(1.0, pg.LD / cfg.LAD)))
    bank = min(max(bank, 0.0), smax)
    pg_q7 = pg.Q7 + cfg.KDMIN
    ind_ini = 0 if (nav["RDOT"] < 0 and nav["D"] > pg_q7) else 1
    tol = cfg.APOGEE_EPSILON_TIGHTER if V > cfg.VI_LOOSE_APOGEE else cfg.APOGEE_EPSILON_TIGHT
    sols = []                                      # (bank, range, vend, rdend)
    for _ in range(runs):
        rng, vend, rdend = _pg_predict(nav["long"], bank, pg.Krho_est, pg_q7, ind_ini, cfg)
        sols.append((bank, rng, vend, rdend))
        miss = rng - target
        if abs(miss) < tol:
            break
        good = [x for x in sols if abs(x[1]) < cfg.INFINITY]
        hi = [x for x in good if x[1] > target]            # too far: more bank
        lo = [x for x in good if x[1] <= target]           # too short: less bank
        esc = [x for x in sols if x[1] >= cfg.INFINITY]
        crash = [x for x in sols if x[1] <= -cfg.INFINITY]
        if hi and lo:                                        # method 2
            a, b = min(hi, key=lambda x: x[1] - target), max(lo, key=lambda x: x[1] - target)
            bank = a[0] + (target - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
        elif lo and esc:                                     # method 3 (low + escape bracket): bisect
            bank = 0.5 * (max(lo, key=lambda x: x[1])[0] + min(esc, key=lambda x: x[0])[0])
        elif len(good) >= 2:                                 # method 6: extrapolate from two on one side
            a, b = good[-1], good[-2]
            bank = a[0] + (target - a[1]) * (b[0] - a[0]) / (b[1] - a[1]) if b[1] != a[1] else a[0]
        elif esc and not good:                               # method 4: march out of escape (more bank)
            bank = bank + lift_inc
        else:                                                # method 4/5: crash only or one capture: less bank
            bank = bank - lift_inc if (crash or (good and good[-1][1] < target)) else bank + lift_inc
        bank = min(max(bank, 0.0), smax)
    best = min(sols, key=lambda x: abs(x[1] - target))
    pg.LD = math.cos(best[0]) * cfg.LAD
    if abs(best[1] - target) < 1000.0 and best[2] > 0:
        pg.VL = best[2]
        pg.GAMMAL = -best[3] / best[2]


# ============================================================================
# The law behind fly_v5's runner
# ============================================================================
_DESIGNS: Dict[tuple, Design] = {}


def design_for(ic: dict, target_phi: float, target_lam: float, cfg: PredGuidConfig) -> Design:
    key = (tuple(round(float(ic[k]), 6) for k in ("altitude_m", "velocity_mps", "gamma_deg", "chi_deg", "phi_deg",
                                                  "lambda_deg")), round(target_phi, 9), round(target_lam, 9),
           cfg.sigma_max_deg, cfg.VFINAL, cfg.LAD, cfg.lod_fixed_frac, cfg.K44, cfg.KA)
    if key not in _DESIGNS:
        _DESIGNS[key] = design_final_phase(ic, target_phi, target_lam, cfg)
    return _DESIGNS[key]


def _vec(phi, lam):
    return np.array([math.cos(phi) * math.cos(lam), math.cos(phi) * math.sin(lam), math.sin(phi)])


def nav_from_state(y: np.ndarray, target_phi: float, target_lam: float) -> dict:
    """data.physical_state row -> App. B navigation quantities (Apollo units) plus the SI longitudinal state."""
    h, phi, lam, V, gam, chi = (float(v) for v in y[:6])
    u = _vec(phi, lam)
    east = np.array([-math.sin(lam), math.cos(lam), 0.0])
    north = np.cross(u, east)
    vh = math.cos(chi) * east + math.sin(chi) * north
    vel = V * (math.sin(gam) * u + math.cos(gam) * vh)
    URT = _vec(target_phi, target_lam)
    uni = np.cross(vel, u)
    uni = uni / max(np.linalg.norm(uni), 1e-12)
    THETA = math.acos(max(-1.0, min(1.0, float(URT @ u))))
    P = planet()
    D = math.hypot(float(y[8]), float(y[9]))
    return {"V": V / FT, "RDOT": V * math.sin(gam) / FT, "D": D / FT, "THETNM": THETA * P["ATK"],
            "LATANG": float(URT @ uni), "alt": h / FT, "gonepast": float(URT @ vh) < 0.0,
            "long": (R_EARTH + h, V, gam, THETA)}


class PredGuidLaw:
    """App. B behind fly_v5.run_guidance_lockstep; one guidance cycle every cfg.cycle_s seconds."""

    def __init__(self, cfg: Optional[PredGuidConfig] = None, device="cpu"):
        self.cfg = cfg or PredGuidConfig()

    def new_flight(self, job: dict, env) -> PredGuidState:
        c = self.cfg
        st = PredGuidState(c)
        st.target_phi, st.target_lam = float(env.target_phi_rad), float(env.target_lam_rad)
        ic = job.get("design_ic") or job["ic"]
        st.design = design_for(ic, st.target_phi, st.target_lam, c)
        # B.2.4: run the targeting twice, then K2ROLL = -sign(LATANG)
        from jepa.data import physical_state
        nav = nav_from_state(physical_state(env), st.target_phi, st.target_lam)
        st.K2ROLL = -1 if nav["LATANG"] >= 0 else 1
        st.last_cmd = float(env.att.sigma_rel_rad)
        return st

    def decide(self, states, ys, t_s, envs, ks) -> List[float]:
        c = self.cfg
        smax = math.radians(c.sigma_max_deg)
        out = []
        for st, y, t, k in zip(states, ys, t_s, ks):
            st.t = float(t)
            st.decided = False
            if st.QUITIND or (k % int(c.cycle_s)) != 0:
                out.append(st.last_cmd)
                continue
            nav = nav_from_state(np.asarray(y), st.target_phi, st.target_lam)
            if c.gamma_safe_deg is not None and abs(float(y[4])) > math.radians(c.gamma_safe_deg):
                st.last_cmd = 0.0                   # steep descent: all lift vertical (cos gamma gate)
                out.append(0.0)
                continue
            ld = executive(st, nav, c, lateral=True)
            mag = math.acos(max(-1.0, min(1.0, ld / c.LAD)))
            sat = mag > smax + 1e-12
            if sat:
                st.n_sat_liftdown += 1
            mag = min(mag, smax)
            # Apollo's K2ROLL = -sign(LATANG) turns toward the target; in the simulator a positive bank turns
            # right (chi decreases), which is the opposite sign for UNI = unit(V x R)
            cmd = -st.K2ROLL * mag if st.K2ROLL != 0 else mag
            st.last_cmd = cmd
            st.decided = True
            st.n_cycles += 1
            st.dec.append((float(t), float(y[0]), float(y[3]), nav["D"], st.SELECTOR, ld, math.degrees(cmd),
                           nav["THETNM"], st.PREDANGL, nav["LATANG"], st.K2ROLL, st.Krho_est, float(sat),
                           st.GONEPAST, st.QUITIND))
            out.append(cmd)
        return out

    def flight_summary(self, st: PredGuidState) -> dict:
        seq = []
        for _, p in st.phase_log:
            if not seq or seq[-1] != p:
                seq.append(p)
        return {"predguid_phases": "-".join(PHASES[p] for p in seq), "predguid_cycles": st.n_cycles,
                "predguid_sat_liftdown": st.n_sat_liftdown,
                "predguid_sat_frac": st.n_sat_liftdown / max(st.n_cycles, 1),
                "predguid_Krho_final": st.Krho_est, **{f"predguid_{k}": v for k, v in st.design.summary().items()}}

    def decisions(self, st: PredGuidState) -> Dict[str, np.ndarray]:
        d = np.asarray(st.dec, dtype=np.float64).reshape(-1, len(DEC_FIELDS_PG))
        return {f"dec_{k}": d[:, i].astype(np.float32) for i, k in enumerate(DEC_FIELDS_PG)}
