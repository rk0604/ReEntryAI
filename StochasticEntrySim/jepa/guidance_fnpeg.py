"""
jepa/guidance_fnpeg.py

FNPEG: Lu's (2014) fully numerical predictor-corrector entry guidance, batched over flights
(jepa/v5_plan.md section 4.1), for the `fnpeg`, `fnpeg_180`, `fnpeg_oracle` arms and the DAgger expert.

    predictor   physics_v4's port of the classical predictor's equations (non-rotating Earth, the
                simulator's atmosphere table and trim schedule, same arithmetic), longitudinal states
                (r, V, gamma) plus range-to-go s carried by Lu's Eq. (17), s' = -V cos(gamma) / r, signed;
                RK4 in time (step cfg.dt_s), float64, vectorized over every prediction of a corrector round
                (numpy: at these batch sizes it is several times faster than torch on a CPU); the bank
                magnitude follows the profile through the simulator's 20 deg/s rate limit; modeled drag
                and lift multiplied by the filters' K_D and K_L along the whole prediction
    profile     |sigma(e)| = sigma0 + (e - e0)/(ef - e0) (sigma_f - sigma0), e = mu/r - V^2/2 (Eq. 19);
                ef from hf = 7.6 km and Vf; the prediction stops at h <= hf or e >= ef (interpolated
                inside the step), or when the path is near vertical (|cos gamma| < 0.05)
    corrector   Eq. (24) step-size-controlled Newton with a finite-difference derivative, Eq. (25) secant
                after the first iteration, lambda = 1/2^i until f = z^2/2 decreases, stop on Eq. (23)
                (|z dz/dsigma0| <= eps) or |z| <= z_tol; warm start from the previous cycle
    lateral     heading offset dpsi = chi - Psi (Eq. 26) against a deadband linear in V; the sign flips
                as soon as |dpsi| exceeds it on the side the current sign makes worse; the sign flip is
                instantaneous in the prediction (it does not enter the longitudinal channel)
    load relief Eqs. (36)-(41) inside every predicted trajectory, and once more on navigation data for
                the output command (Lu's point 3, Eq. 42); Lu's CEV values delta = 0.02 and k0 = 100
                (nondimensional) converted to seconds
    filters     Brunner and Lu (2008) first-order fading-memory filters K_D, K_L on measured over
                modeled drag and lift acceleration, active from 0.05 g
    cadence     a solution every guidance cycle (1 Hz) while the aerodynamic acceleration is
                >= 1.52 m/s^2 (Lu); below it the command holds the current bank

The flight runner (fly_v5.run_guidance_lockstep) calls FNPEG.decide once per second with every
flight's navigation state; all predictions of a corrector round run as one vectorized batch.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

import constants
from jepa import physics_v4 as PH
from jepa.data import G0, R_EARTH

MU = float(constants.MU)
LU_TIME_S = math.sqrt(R_EARTH / G0)             # Lu's time unit sqrt(R0/g0), ~806 s
DEC_FIELDS_FNPEG = ("t_s", "h_m", "V_mps", "a_g", "K_D", "K_L", "guided", "sigma0_deg", "z_km", "n_iter",
                    "n_pred", "sigma_base_deg", "cmd_deg", "dpsi_deg", "deadband_deg", "relief")


def lu_k0(k0_nd: float) -> float:
    """Lu's nondimensional load-relief gain -> 1/s (L cos sigma in m/s^2, hdot in m/s)."""
    return float(k0_nd) * math.sqrt(G0 / R_EARTH)


def lu_delta(delta_nd: float) -> float:
    """Lu's nondimensional prediction horizon -> s."""
    return float(delta_nd) * LU_TIME_S


@dataclass
class FNPEGConfig:
    sigma_max_deg: float = 90.0             # D2: 90 deg for the primary arm, 180 for fnpeg_180
    sigma_f_deg: Optional[float] = None     # None: sigma_f = sigma0 (constant profile), Lu's option
    dt_s: float = 2.0                       # RK4 step, tuned in {1, 2}
    hf_m: float = 7_600.0                   # drogue altitude = Lu's CEV final altitude
    vf_mps: float = 150.0                   # set from DEV48 classical flights by the study (Lu's CEV: 150)
    deadband_c0_deg: float = 1.0            # deadband(V) = c0 + c1 V: 1 deg at rest, 5 deg at 7,700 m/s (tuned)
    deadband_c1_deg_per_mps: float = 5.2e-4
    gamma_safe_deg: Optional[float] = 30.0  # |gamma| beyond it: wings level (the classical controller's override)
    load_relief: bool = True
    amax_g: float = 4.5
    margin_g: float = 0.1                   # Lu: "a small margin in amax of 0.1 g or so"
    delta_nd: float = 0.02                  # Lu's CEV value (16 s)
    k0_nd: float = 100.0                    # Lu's CEV value
    filters: bool = True
    filter_beta: float = 0.9                # K <- beta K + (1 - beta) measured / modeled, at 1 Hz
    filter_on_g: float = 0.05
    oracle: bool = False                    # predictor gets the true K_D, K_L; filters off
    a_guid_mps2: float = 1.52               # Lu: guidance only while sqrt(L^2 + D^2) >= 5 ft/s^2
    sigma_init_deg: float = 45.0            # first guess of the first cycle
    fd_deg: float = 1.0                     # finite-difference step of the first Newton iteration
    max_iter: int = 8
    max_halvings: int = 5
    z_tol_km: float = 0.05
    eps_km2: float = 5.0                    # Eq. (23): |z dz/dsigma0| <= eps (z in km, sigma0 in rad)
    t_max_s: float = 2_500.0                # longest prediction
    t_hold_s: float = 10.0                  # hold the command once the predicted time to go is shorter
    rate_dps: float = 20.0                  # the simulator's bank rate limit (also inside the predictor)

    def to_dict(self) -> dict:
        return asdict(self)


# ============================================================================
# The longitudinal model (vectorized numpy, float64)
# ============================================================================
class LongModel:
    """physics_v4's equations (longitudinal part) with the atmosphere table and trim schedule as arrays."""

    def __init__(self):
        atm = PH.atm_table("cpu")
        self.alt = atm.alt.numpy().astype(np.float64)
        self.log_rho = atm.log_rho.numpy().astype(np.float64)
        self.T = atm.T.numpy().astype(np.float64)
        self.alt_max = float(atm.alt_max)
        self.n_seg = len(self.alt) - 1
        self.dlog = np.diff(self.log_rho) / np.diff(self.alt)
        sch = np.asarray(sorted(PH.TRIM_SCHEDULE, key=lambda r: r[0]), dtype=np.float64)
        self.f, self.cd, self.ld = sch[:, 0].copy(), sch[:, 1].copy(), sch[:, 2].copy()
        p = PH.nominal_params()
        self.S, self.m = float(p["ref_area_m2"]), float(p["mass_kg"])

    def _seg(self, hh: np.ndarray) -> np.ndarray:
        """physics_v4.AtmTable._interp's segment: searchsorted(right) - 1, clamped."""
        return np.clip(np.searchsorted(self.alt, hh, side="right") - 1, 0, self.n_seg - 1)

    def aero(self, h: np.ndarray, V: np.ndarray):
        """(drag accel, lift accel at lift scale 1, d ln(rho)/dh) with nominal parameters: m/s^2, m/s^2, 1/m."""
        hc = np.maximum(h, 0.0)
        hh = np.minimum(hc, self.alt_max)
        j = self._seg(hh)
        x0 = self.alt[j]
        w = (hh - x0) / (self.alt[j + 1] - x0)
        rho = np.exp(self.log_rho[j] + w * (self.log_rho[j + 1] - self.log_rho[j]))
        T = self.T[j] + w * (self.T[j + 1] - self.T[j])
        space = hc / 1000.0 > 1000.0
        rho = np.where(space, 0.0, rho)
        t_ok = (T > 0.0) & ~space
        # physics_v4.aero_coefficients (same arithmetic), lift scale 1
        V_kfps = V * 0.00328084
        mach = V / np.sqrt(1.4 * 287.053 * np.where(t_ok, T, 1.0))
        a = (V_kfps - 8.8) / (9.8 - 8.8)
        fmv = np.where(V_kfps >= 9.8, V_kfps,
                       np.where(~t_ok, V_kfps, np.where(V_kfps <= 8.8, mach, (1.0 - a) * mach + a * V_kfps)))
        i = np.clip(np.searchsorted(self.f, fmv, side="left") - 1, 0, len(self.f) - 2)
        f0 = self.f[i]
        al = (fmv - f0) / (self.f[i + 1] - f0)
        CD = self.cd[i] + al * (self.cd[i + 1] - self.cd[i])
        LD = self.ld[i] + al * (self.ld[i + 1] - self.ld[i])
        lo, hi = fmv <= self.f[0], fmv >= self.f[-1]
        CD = np.where(lo, self.cd[0], np.where(hi, self.cd[-1], CD))
        LD = np.where(lo, self.ld[0], np.where(hi, self.ld[-1], LD))
        qS_m = 0.5 * rho * V * V * self.S / self.m
        return qS_m * CD, qS_m * CD * LD, self.dlog[j]


_MODEL: Optional[LongModel] = None


def long_model() -> LongModel:
    global _MODEL
    if _MODEL is None:
        _MODEL = LongModel()
    return _MODEL


def energy(r, V):
    """Lu's specific energy, dimensional: e = mu/r - V^2/2 (increases monotonically during entry)."""
    return MU / r - 0.5 * V * V


def relief_cos(cos_base, L, D, a, V, sin_gam, beta_h, amax, delta, k0):
    """
    Eqs. (38)-(41): cos(sigma_cmd) = cos(sigma_base) - k0 (hdot - hdot_ref) / L, hdot_ref = V max(sin gamma, Ua),
    Ua = [amax - a (1 - 2 D delta / V)] / (a V beta delta), beta = d ln(rho)/dh < 0. Returns (cos, active).
    """
    Vs = np.maximum(V, 1.0)
    den = a * Vs * beta_h * delta
    ok = (a > 1e-6) & (np.abs(den) > 1e-12)
    Ua = np.where(ok, (amax - a * (1.0 - 2.0 * D * delta / Vs)) / np.where(ok, den, 1.0), -2.0)
    active = Ua > sin_gam
    corr = k0 * Vs * (sin_gam - np.maximum(sin_gam, Ua)) / np.maximum(L, 1e-9)
    return cos_base - np.where(active, corr, 0.0), active


def _deriv(m: LongModel, r, V, gam, cs, KD, KL):
    D0, L0, _ = m.aero(r - R_EARTH, V)
    return _deriv_from(D0, L0, r, V, gam, cs, KD, KL)


def _deriv_from(D0, L0, r, V, gam, cs, KD, KL):
    g = MU / (r * r)
    Vs = np.maximum(V, 1e-9)
    cg, sg = np.cos(gam), np.sin(gam)
    return V * sg, -KD * D0 - g * sg, KL * L0 * cs / Vs + (Vs / r - g / Vs) * cg, -V * cg / r


def predict(x0: np.ndarray, sig_now: np.ndarray, sigma0: np.ndarray, sigma_f: np.ndarray, ef: np.ndarray,
            KD: np.ndarray, KL: np.ndarray, cfg: FNPEGConfig, record: bool = False) -> dict:
    """
    One prediction per row. x0 [B, 4] = (r, V, gamma, s) (s in rad, > 0 toward the target), sig_now [B] the
    realized |bank|, sigma0/sigma_f [B] the profile (rad), ef [B]. Returns z_km [B] (range to go at the end:
    > 0 undershoot), t_end_s [B], a_peak_g [B] (predicted peak aerodynamic load), and, with record, the path.
    """
    m = long_model()
    dt = float(cfg.dt_s)
    x0 = np.asarray(x0, dtype=np.float64)
    r, V, gam, s = (x0[:, i].copy() for i in range(4))
    e0 = energy(r, V)
    den = ef - e0
    den = np.where(np.abs(den) < 1e-6, 1e-6, den)
    smax = math.radians(cfg.sigma_max_deg)
    rate = math.radians(cfg.rate_dps) * dt
    sig = np.clip(np.asarray(sig_now, dtype=np.float64), 0.0, smax)
    amax = (cfg.amax_g - cfg.margin_g) * G0
    delta, k0 = lu_delta(cfg.delta_nd), lu_k0(cfg.k0_nd)
    hf = R_EARTH + cfg.hf_m
    B = len(r)
    done = np.zeros(B, dtype=bool)
    z = np.zeros(B)
    t_end = np.full(B, float(cfg.t_max_s))
    a_peak = np.zeros(B)
    path = [] if record else None
    idx = np.arange(B)
    # integrate only the live rows (rows drop out as they finish)
    R, Vv, Gm, Sg, Sv = r, V, gam, sig, s
    S0, SF, EF, E0, DEN, kd, kl = sigma0, sigma_f, ef, e0, den, KD, KL
    n_steps = int(math.ceil(cfg.t_max_s / dt))
    for k in range(n_steps):
        e = energy(R, Vv)
        prof = np.clip(S0 + (e - E0) / DEN * (SF - S0), 0.0, smax)
        D0, L0, beta = m.aero(R - R_EARTH, Vv)
        D, L = kd * D0, kl * L0
        a = np.sqrt(D * D + L * L)
        a_peak[idx] = np.maximum(a_peak[idx], a)
        cos_cmd = np.cos(prof)
        if cfg.load_relief:
            cos_cmd, _ = relief_cos(cos_cmd, L, D, a, Vv, np.sin(Gm), beta, amax, delta, k0)
        cmd = np.arccos(np.clip(cos_cmd, math.cos(smax), 1.0))
        if cfg.gamma_safe_deg is not None:
            cmd = np.where(Gm < -math.radians(cfg.gamma_safe_deg), 0.0, cmd)
        cs = np.cos(Sg)
        # RK4 with the bank held over the step (as the simulator's actuator holds it between guidance calls)
        k1 = _deriv_from(D0, L0, R, Vv, Gm, cs, kd, kl)
        k2 = _deriv(m, R + 0.5 * dt * k1[0], Vv + 0.5 * dt * k1[1], Gm + 0.5 * dt * k1[2], cs, kd, kl)
        k3 = _deriv(m, R + 0.5 * dt * k2[0], Vv + 0.5 * dt * k2[1], Gm + 0.5 * dt * k2[2], cs, kd, kl)
        k4 = _deriv(m, R + dt * k3[0], Vv + dt * k3[1], Gm + dt * k3[2], cs, kd, kl)
        Rn = R + dt / 6.0 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        Vn = Vv + dt / 6.0 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        Gn = Gm + dt / 6.0 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        Sn = Sv + dt / 6.0 * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3])
        en = energy(Rn, Vn)
        with np.errstate(divide="ignore", invalid="ignore"):
            f_h = np.where(Rn < hf, (R - hf) / np.maximum(R - Rn, 1e-9), 2.0)
            f_e = np.where(en >= EF, (EF - e) / np.maximum(en - e, 1e-12), 2.0)
        f = np.clip(np.minimum(f_h, f_e), 0.0, 2.0)
        vert = np.abs(np.cos(Gn)) < PH.COS_GAMMA_STOP
        bad = ~np.isfinite(Rn) | ~np.isfinite(Sn) | (Vn < 1.0)
        end = (f <= 1.0) | vert | bad
        if record:
            row = np.full((B, 5), np.nan)
            row[idx] = np.stack([R - R_EARTH, Vv, Gm, Sv, Sg], -1)
            path.append(row)
        if end.any():
            fe = np.where(f <= 1.0, f, 1.0)[end]
            ie = idx[end]
            z[ie] = (Sv[end] + fe * (Sn[end] - Sv[end])) * (R_EARTH / 1000.0)
            t_end[ie] = (k + fe) * dt
            done[ie] = True
            keep = ~end
            idx = idx[keep]
            if not len(idx):
                break
            Rn, Vn, Gn, Sn = Rn[keep], Vn[keep], Gn[keep], Sn[keep]
            Sg, cmd = Sg[keep], cmd[keep]
            S0, SF, EF, E0, DEN, kd, kl = S0[keep], SF[keep], EF[keep], E0[keep], DEN[keep], kd[keep], kl[keep]
        R, Vv, Gm, Sv = Rn, Vn, Gn, Sn
        Sg = Sg + np.clip(cmd - Sg, -rate, rate)
    if len(idx) and not done[idx].all():
        z[idx] = Sv * (R_EARTH / 1000.0)               # ran out of time: range to go at t_max
    out = {"z_km": z, "t_end_s": t_end, "a_peak_g": a_peak / G0}
    if record:
        out["path"] = np.stack(path, 1)
    return out


# ============================================================================
# Navigation helpers
# ============================================================================
def range_and_azimuth(phi: float, lam: float, phi_t: float, lam_t: float):
    """Great-circle central angle to the target (rad) and its azimuth in the simulator's chi convention."""
    dl = lam_t - lam
    a = math.sin((phi_t - phi) / 2) ** 2 + math.cos(phi) * math.cos(phi_t) * math.sin(dl / 2) ** 2
    s = 2.0 * math.asin(min(1.0, math.sqrt(a)))
    brg = math.atan2(math.sin(dl) * math.cos(phi_t),
                     math.cos(phi) * math.sin(phi_t) - math.sin(phi) * math.cos(phi_t) * math.cos(dl))
    return s, math.pi / 2.0 - brg          # bearing from north clockwise -> chi (from east toward north)


def wrap_pi(x: float) -> float:
    return (x + math.pi) % (2.0 * math.pi) - math.pi


@dataclass
class FNPEGFlight:
    """Per-flight guidance state."""
    target_phi: float
    target_lam: float
    K_D: float = 1.0
    K_L: float = 1.0
    sigma0: Optional[float] = None
    sign: float = 0.0
    hold: bool = False
    last_cmd: float = 0.0
    decided: bool = False
    true_K: Optional[tuple] = None          # (K_D, K_L) for the oracle
    dec: List[tuple] = field(default_factory=list)
    iters: List[int] = field(default_factory=list)


class FNPEG:
    """Batched FNPEG. One instance flies any number of flights; decide() is one 1 Hz guidance cycle."""

    def __init__(self, cfg: Optional[FNPEGConfig] = None, device="cpu"):
        self.cfg = cfg or FNPEGConfig()
        self.model = long_model()
        self.ef = float(energy(R_EARTH + self.cfg.hf_m, self.cfg.vf_mps))

    def new_flight(self, target_phi: float, target_lam: float, domain: Optional[dict] = None) -> FNPEGFlight:
        f = FNPEGFlight(float(target_phi), float(target_lam))
        if self.cfg.oracle and domain is not None:
            ds, ls, ms = float(domain["density_scale"]), float(domain["lift_scale"]), float(domain["mass_scale"])
            f.true_K = (ds / ms, ds * ls / ms)
            f.K_D, f.K_L = f.true_K
        return f

    def deadband(self, V: float) -> float:
        return math.radians(self.cfg.deadband_c0_deg + self.cfg.deadband_c1_deg_per_mps * V)

    def update_filters(self, flights: Sequence[FNPEGFlight], Y: np.ndarray) -> np.ndarray:
        """Brunner-Lu filters on measured / modeled acceleration. Returns the measured load (m/s^2)."""
        c = self.cfg
        D0, L0, _ = self.model.aero(Y[:, 0], Y[:, 3])
        a_meas = np.hypot(Y[:, 8], Y[:, 9])
        if c.filters and not c.oracle:
            b = c.filter_beta
            for i, f in enumerate(flights):
                if a_meas[i] >= c.filter_on_g * G0 and D0[i] > 0 and L0[i] > 0:
                    f.K_D = b * f.K_D + (1.0 - b) * Y[i, 8] / D0[i]
                    f.K_L = b * f.K_L + (1.0 - b) * Y[i, 9] / L0[i]
        return a_meas

    def decide(self, flights: Sequence[FNPEGFlight], ys: Sequence[np.ndarray], t_s: Sequence[float]) -> List[float]:
        """
        One guidance cycle for every flight in `flights`. ys: data.physical_state rows (h, phi, lam, V, gamma,
        chi, sigma, roll rate, a_drag, a_lift), the measured accelerations from the true parameters. Returns
        the bank commands (rad, signed, |cmd| <= sigma_max).
        """
        c = self.cfg
        n = len(flights)
        Y = np.asarray(ys, dtype=np.float64).reshape(n, -1)
        a_meas = self.update_filters(flights, Y)
        _, _, beta = self.model.aero(Y[:, 0], Y[:, 3])
        guided = [i for i, f in enumerate(flights) if a_meas[i] >= c.a_guid_mps2 and not f.hold]
        smax = math.radians(c.sigma_max_deg)
        sol = self._solve([flights[i] for i in guided], Y[guided]) if guided else {}
        sol = {guided[j]: v for j, v in sol.items()}
        amax = (c.amax_g - c.margin_g) * G0
        cmds = []
        for i, f in enumerate(flights):
            _, psi = range_and_azimuth(Y[i, 1], Y[i, 2], f.target_phi, f.target_lam)
            dpsi = wrap_pi(Y[i, 5] - psi)
            db = self.deadband(Y[i, 3])
            row = sol.get(i)
            relief = False
            f.decided = row is not None
            if c.gamma_safe_deg is not None and abs(Y[i, 4]) > math.radians(c.gamma_safe_deg):
                cmd, sb = 0.0, 0.0                   # steep descent: all lift vertical (cos gamma gate)
                f.decided = False
            elif f.hold:
                cmd = f.last_cmd
                sb = abs(cmd)
            elif row is not None:
                sb = float(row["sigma_base"])
                cb = math.cos(sb)
                if c.load_relief:
                    cc, act = relief_cos(np.array([cb]), np.array([Y[i, 9]]), np.array([Y[i, 8]]),
                                         np.array([a_meas[i]]), np.array([Y[i, 3]]), np.array([math.sin(Y[i, 4])]),
                                         np.array([beta[i]]), amax, lu_delta(c.delta_nd), lu_k0(c.k0_nd))
                    cb, relief = float(cc[0]), bool(act[0])
                mag = math.acos(min(1.0, max(math.cos(smax), cb)))
                if f.sign == 0.0:
                    f.sign = 1.0 if dpsi >= 0.0 else -1.0
                if dpsi > db:
                    f.sign = 1.0                     # a positive bank turns right (chi decreases)
                elif dpsi < -db:
                    f.sign = -1.0
                cmd = f.sign * mag
                if row["t_end_s"] < c.t_hold_s:
                    f.hold = True
            else:
                cmd = float(Y[i, 6])                 # Lu: below 1.52 m/s^2 the command is the current bank
                sb = abs(cmd)
            f.last_cmd = cmd
            cmds.append(cmd)
            f.dec.append((float(t_s[i]), Y[i, 0], Y[i, 3], a_meas[i] / G0, f.K_D, f.K_L, float(row is not None),
                          math.degrees(f.sigma0) if f.sigma0 is not None else np.nan,
                          row["z_km"] if row is not None else np.nan, row["n_iter"] if row is not None else -1,
                          row["n_pred"] if row is not None else 0, math.degrees(sb), math.degrees(cmd),
                          math.degrees(dpsi), math.degrees(db), float(relief)))
            if row is not None:
                f.iters.append(int(row["n_iter"]))
        return cmds

    def _solve(self, flights: Sequence[FNPEGFlight], Y: np.ndarray) -> Dict[int, dict]:
        """Batched corrector over the guided flights. Returns {row: sigma_base, z_km, t_end_s, n_iter, n_pred}."""
        c = self.cfg
        n = len(flights)
        smax = math.radians(c.sigma_max_deg)
        s_go = np.array([range_and_azimuth(Y[i, 1], Y[i, 2], f.target_phi, f.target_lam)[0]
                         for i, f in enumerate(flights)])
        x0 = np.stack([Y[:, 0] + R_EARTH, Y[:, 3], Y[:, 4], s_go], -1)
        sig_now = np.abs(Y[:, 6])
        KD = np.array([f.K_D for f in flights])
        KL = np.array([f.K_L for f in flights])
        ef = np.full(n, self.ef)
        sf_fixed = None if c.sigma_f_deg is None else math.radians(c.sigma_f_deg)

        def run(rows: np.ndarray, s0: np.ndarray):
            sf = s0 if sf_fixed is None else np.full_like(s0, sf_fixed)
            o = predict(x0[rows], sig_now[rows], s0, sf, ef[rows], KD[rows], KL[rows], c)
            return o["z_km"], o["t_end_s"]

        sigma = np.clip(np.array([f.sigma0 if f.sigma0 is not None else math.radians(c.sigma_init_deg)
                                  for f in flights]), 0.0, smax)
        h_fd = math.radians(c.fd_deg)
        pert = np.where(sigma + h_fd <= smax, sigma + h_fd, sigma - h_fd)
        rows = np.arange(n)
        zz, tt = run(np.concatenate([rows, rows]), np.concatenate([sigma, pert]))
        z, t_end = zz[:n].copy(), tt[:n].copy()
        dz = (zz[n:] - z) / (pert - sigma)
        n_pred = np.full(n, 2)
        n_iter = np.zeros(n, dtype=int)
        conv = (np.abs(z) <= c.z_tol_km) | (np.abs(z * dz) <= c.eps_km2)
        for _ in range(c.max_iter):
            act = np.nonzero(~conv)[0]
            if not len(act):
                break
            lam = np.ones(len(act))
            slope = np.where(np.abs(dz[act]) > 1e-9, dz[act], 1e-9)
            cand = np.clip(sigma[act] - lam * z[act] / slope, 0.0, smax)
            zc, tc = run(act, cand)
            n_pred[act] += 1
            for _h in range(c.max_halvings):
                worse = np.abs(zc) >= np.abs(z[act])
                if not worse.any():
                    break
                wi = np.nonzero(worse)[0]
                lam[wi] *= 0.5
                cand[wi] = np.clip(sigma[act][wi] - lam[wi] * z[act][wi] / slope[wi], 0.0, smax)
                zc[wi], tc[wi] = run(act[wi], cand[wi])
                n_pred[act[wi]] += 1
            improved = np.abs(zc) < np.abs(z[act])
            moved = np.abs(cand - sigma[act]) > 1e-12
            new_dz = np.where(moved, (zc - z[act]) / np.where(moved, cand - sigma[act], 1.0), dz[act])
            n_iter[act] += 1
            up = act[improved]
            sigma[up], z[up], t_end[up], dz[up] = cand[improved], zc[improved], tc[improved], new_dz[improved]
            conv[up] = (np.abs(z[up]) <= c.z_tol_km) | (np.abs(z[up] * dz[up]) <= c.eps_km2)
            conv[act[~improved]] = True               # no decrease of f along the step: a stationary point
        out = {}
        for i, f in enumerate(flights):
            f.sigma0 = float(sigma[i])
            out[i] = {"sigma_base": float(sigma[i]), "z_km": float(z[i]), "t_end_s": float(t_end[i]),
                      "n_iter": int(n_iter[i]), "n_pred": int(n_pred[i])}
        return out


def decisions_array(f: FNPEGFlight) -> Dict[str, np.ndarray]:
    d = np.asarray(f.dec, dtype=np.float64).reshape(-1, len(DEC_FIELDS_FNPEG))
    return {f"dec_{k}": d[:, i].astype(np.float32) for i, k in enumerate(DEC_FIELDS_FNPEG)}
