"""
validation_tools.py — Deterministic validation kernel for ODE systems.

Design contract: every function here is a pure, deterministic computation.
Verdicts (the `passed` field) are decided by numerical integration, never by
an LLM. The agent layer (validation_agent.py) may CHOOSE which of these to
run and INTERPRET their output, but the numbers and pass/fail come from here.

Reference solver is scipy's LSODA — the same Adams/BDF-switching algorithm as
R's deSolve::lsoda — so convergence tests compare fixed-step RK4 against an
independent, stiff-stable solver family rather than just a finer RK4.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable
import numpy as np
from scipy.integrate import solve_ivp


# ──────────────────────────────────────────────────────────────────────────
# System abstraction
# ──────────────────────────────────────────────────────────────────────────
@dataclass
class ODESystem:
    """A validatable ODE implementation."""
    name: str
    state_names: list[str]
    rhs: Callable[[np.ndarray, dict], np.ndarray]   # f(y, params) -> dy/dt
    params: dict
    y0: np.ndarray
    t_end: float
    dose_times: list[float] = field(default_factory=list)
    dose_amount: float = 0.0
    dose_target: int = -1                            # state index, -1 = no dosing
    # name -> (fn(y)->scalar, expected_constant_value)
    invariants: dict[str, tuple[Callable[[np.ndarray], float], float]] = field(default_factory=dict)

    def idx(self, name: str) -> int:
        return self.state_names.index(name)

    def with_params(self, **overrides) -> "ODESystem":
        p = dict(self.params); p.update(overrides)
        return ODESystem(self.name, self.state_names, self.rhs, p,
                         self.y0.copy(), self.t_end, list(self.dose_times),
                         self.dose_amount, self.dose_target, dict(self.invariants))


# ──────────────────────────────────────────────────────────────────────────
# Integrators (the only place numerics happen)
# ──────────────────────────────────────────────────────────────────────────
def _rk4_step(rhs, y, p, dt):
    k1 = rhs(y, p)
    k2 = rhs(y + 0.5 * dt * k1, p)
    k3 = rhs(y + 0.5 * dt * k2, p)
    k4 = rhs(y + dt * k3, p)
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def integrate(sys: ODESystem, dt: float = 0.05, method: str = "rk4"):
    """Integrate with dose events. method in {'rk4','ref'}. ref = LSODA.
    Returns (t[:], Y[:, n_states]) recording every RK4 step (ref: daily)."""
    y = sys.y0.astype(float).copy()
    dose_set = set(sys.dose_times)
    bounds = sorted(set([0.0] + list(sys.dose_times) + [sys.t_end]))
    ts, Ys = [], []
    for bi in range(len(bounds) - 1):
        t0, t1 = bounds[bi], bounds[bi + 1]
        if t0 in dose_set and sys.dose_target >= 0:
            y[sys.dose_target] += sys.dose_amount
        if method == "rk4":
            n = max(1, int(round((t1 - t0) / dt)))
            h = (t1 - t0) / n
            for s in range(n):
                ts.append(t0 + s * h); Ys.append(y.copy())
                y = _rk4_step(sys.rhs, y, sys.params, h)
        elif method in ("ref", "lsoda", "radau", "bdf"):
            scipy_method = {"ref": "LSODA", "lsoda": "LSODA",
                            "radau": "Radau", "bdf": "BDF"}[method]
            te = np.arange(t0, t1, 1.0)                   # output days [t0, t1)
            te_full = te if (te.size and te[-1] == t1) else np.append(te, t1)
            # ALWAYS integrate through t1 so the carried state `y` is the true
            # segment endpoint (not t1-1). Record t1 only on the final segment;
            # interior boundary days are recorded post-dose by the next segment.
            sol = solve_ivp(lambda tt, yy: sys.rhs(yy, sys.params), (t0, t1), y,
                            method=scipy_method, rtol=1e-9, atol=1e-11, t_eval=te_full)
            n_record = sol.t.shape[0] if bi == len(bounds) - 2 else te.shape[0]
            for j in range(n_record):
                ts.append(sol.t[j]); Ys.append(sol.y[:, j].copy())
            y = sol.y[:, -1].copy()                       # state AT t1 (correct)
        else:
            raise ValueError(f"unknown method {method}")
    if not ts or ts[-1] != sys.t_end:                 # RK4 needs the endpoint;
        ts.append(sys.t_end); Ys.append(y.copy())     # ref-family already has it
    return np.array(ts), np.array(Ys)


# ──────────────────────────────────────────────────────────────────────────
# Result types — frozen so the agent can't mutate a verdict
# ──────────────────────────────────────────────────────────────────────────
@dataclass(frozen=True)
class CheckResult:
    name: str
    passed: bool
    metrics: dict
    summary: str            # human/LLM-readable; verdict already baked in


# ──────────────────────────────────────────────────────────────────────────
# Validators
# ──────────────────────────────────────────────────────────────────────────
def check_conservation(sys: ODESystem, dt: float = 0.05, tol: float = 1e-9) -> CheckResult:
    """Max drift of each declared invariant over an RK4 trajectory."""
    t, Y = integrate(sys, dt=dt, method="rk4")
    worst, worst_inv = 0.0, ""
    per = {}
    for inv_name, (fn, expected) in sys.invariants.items():
        drift = float(np.max(np.abs(np.apply_along_axis(fn, 1, Y) - expected)))
        per[inv_name] = drift
        if drift > worst:
            worst, worst_inv = drift, inv_name
    passed = worst <= tol
    return CheckResult(
        "conservation", passed, {"max_drift": worst, "worst_invariant": worst_inv,
                                 "per_invariant": per, "tol": tol},
        f"max invariant drift {worst:.2e} (@{worst_inv}) vs tol {tol:.0e} -> "
        f"{'PASS' if passed else 'FAIL'}")


def check_convergence(sys: ODESystem, dt_list=(0.1, 0.05, 0.025),
                      production_dt: float = 0.05, readouts=("RGC", "RPE", "M1"),
                      atol: float = 1e-4, rtol: float = 1e-2) -> CheckResult:
    """Fixed-step RK4 at several dt vs an independent LSODA reference.
    Verdict: at the PRODUCTION dt, each readout state must match LSODA within a
    mixed tolerance atol + rtol*|ref| (so near-zero fast intermediates don't
    trip the verdict). The full sweep + worst-of-all-states is reported as a
    diagnostic for the agent to interpret."""
    dts = sorted(set(list(dt_list) + [production_dt]), reverse=True)
    _, Yref = integrate(sys, method="ref")
    yref = Yref[-1]
    rows = []
    for dt in dts:
        _, Y = integrate(sys, dt=dt, method="rk4")
        yend = Y[-1]
        thr = atol + rtol * np.abs(yref)
        ratio = np.abs(yend - yref) / thr            # >1 means out of tolerance
        wi = int(np.argmax(ratio))
        rows.append({"dt": dt,
                     "readouts": {s: float(yend[sys.idx(s)]) for s in readouts},
                     "worst_state_allstates": sys.state_names[wi],
                     "worst_ratio_allstates": float(ratio[wi])})
    prod = next(r for r in rows if r["dt"] == production_dt)
    # verdict on readouts only, at production dt
    _, Yp = integrate(sys, dt=production_dt, method="rk4")
    yp = Yp[-1]
    readout_ok = {s: abs(yp[sys.idx(s)] - yref[sys.idx(s)]) <=
                  atol + rtol * abs(yref[sys.idx(s)]) for s in readouts}
    passed = all(readout_ok.values())
    return CheckResult(
        "convergence", passed,
        {"reference": "scipy LSODA", "production_dt": production_dt, "sweep": rows,
         "readout_ok": readout_ok, "atol": atol, "rtol": rtol,
         "ref_readouts": {s: float(yref[sys.idx(s)]) for s in readouts}},
        f"at production dt={production_dt}, readouts {list(readouts)} match LSODA "
        f"(worst all-state ratio over sweep "
        f"{max(r['worst_ratio_allstates'] for r in rows):.1f}×tol) -> "
        f"{'PASS' if passed else 'FAIL'}")


def check_equivalence(a: ODESystem, b: ODESystem, dt: float = 0.05,
                      rel_tol: float = 1e-6) -> CheckResult:
    """Two implementations must produce the same trajectory (shared states)."""
    _, Ya = integrate(a, dt=dt, method="rk4")
    _, Yb = integrate(b, dt=dt, method="rk4")
    shared = [s for s in a.state_names if s in b.state_names]
    worst, worst_state = 0.0, ""
    for s in shared:
        va, vb = Ya[-1, a.idx(s)], Yb[-1, b.idx(s)]
        denom = max(1e-6, abs(va))
        rel = abs(va - vb) / denom
        if rel > worst:
            worst, worst_state = rel, s
    passed = worst <= rel_tol
    return CheckResult(
        "equivalence", passed,
        {"a": a.name, "b": b.name, "max_rel_diff": worst, "worst_state": worst_state,
         "rel_tol": rel_tol},
        f"{a.name} vs {b.name}: max rel diff {worst:.2e} (@{worst_state}) -> "
        f"{'PASS' if passed else 'FAIL'}")


def check_steady_state(sys: ODESystem, overrides: dict, dt: float = 0.05,
                       tol: float = 1e-6) -> CheckResult:
    """Under `overrides` (e.g. IOP_target=IOP_normal), the healthy IC must be a
    fixed point: states should not move from y0."""
    s2 = sys.with_params(**overrides)
    s2.dose_times = []                       # no dosing for the fixed-point test
    t, Y = integrate(s2, dt=dt, method="rk4")
    dev = np.max(np.abs(Y - s2.y0), axis=0)
    worst = float(np.max(dev)); arg = sys.state_names[int(np.argmax(dev))]
    passed = worst <= tol
    return CheckResult(
        "steady_state", passed,
        {"overrides": overrides, "max_deviation": worst, "worst_state": arg, "tol": tol},
        f"under {overrides}: max drift from healthy fixed point {worst:.2e} "
        f"(@{arg}) -> {'PASS' if passed else 'FAIL'}")


def sweep_stiffness(sys: ODESystem, param: str, factors=(1, 2, 5, 10, 20, 50),
                    dt: float = 0.05, readout: str = "RGC", rtol: float = 1e-2) -> CheckResult:
    """Scale a rate constant up and find where fixed-step RK4 *diverges from
    LSODA* (the regime where you must switch to a stiff solver). Uses state
    values, NOT the conserved invariant — RK4 preserves linear invariants
    exactly even when the integration is inaccurate, so the invariant is blind
    to stiffness. Tripwire: RK4 endpoint goes non-finite, or readout disagrees
    with LSODA by more than rtol."""
    base = sys.params[param]
    ri = sys.idx(readout)
    rows, break_factor = [], None
    for f in factors:
        s2 = sys.with_params(**{param: base * f})
        _, Yr = integrate(s2, dt=dt, method="rk4")
        _, Yl = integrate(s2, method="ref")
        rk4_end, lsoda_end = float(Yr[-1, ri]), float(Yl[-1, ri])
        finite = bool(np.all(np.isfinite(Yr[-1])))
        denom = max(1e-6, abs(lsoda_end))
        rel = abs(rk4_end - lsoda_end) / denom if finite else float("inf")
        rows.append({"factor": f, f"rk4_{readout}": rk4_end,
                     f"lsoda_{readout}": lsoda_end, "rel_diff": rel, "finite": finite})
        if break_factor is None and (not finite or rel > rtol):
            break_factor = f
    passed = break_factor is None
    msg = (f"RK4 (dt={dt}) tracks LSODA on {readout} across {param}×{factors[-1]}"
           if passed else
           f"RK4 (dt={dt}) diverges from LSODA on {readout} at {param}×{break_factor} "
           f"(switch to a stiff solver beyond this)")
    return CheckResult("stiffness_sweep", passed,
                       {"param": param, "readout": readout, "rows": rows,
                        "break_factor": break_factor, "rtol": rtol}, msg)


# ──────────────────────────────────────────────────────────────────────────
# Reference implementation: glaucoma QSP v6 (JS index order: C3@6, C1q@7)
# ──────────────────────────────────────────────────────────────────────────
def _hill(C, Emax, EC50, g):
    C = max(C, 0.0)
    return Emax * C ** g / (EC50 ** g + C ** g)


def glaucoma_rhs(y, p):
    RPE = min(max(y[0], 0), 1); DAMPs = max(y[1], 0); M0 = max(y[2], 0)
    M_mig = max(y[3], 0); M1 = max(y[4], 0); M2 = max(y[5], 0)
    C3 = max(y[6], 0); C1q = max(y[7], 0); C3a = max(y[8], 0); C5a = max(y[9], 0)
    Cp = max(y[10], 0); Ca = max(y[11], 0); NTF = max(y[12], 0); RGC = min(max(y[13], 0), 1)
    A = max(y[14], 0); Cp2 = max(y[15], 0); R_des = min(max(y[16], 0), 1)
    Stress = max(0, (p["IOP_target"] - p["IOP_normal"]) / p["IOP_normal"])
    Pb = _hill(Cp2, p["Emax_C1qblock"], p["EC50_pep"], p["gamma_pep"])
    Pc = _hill(Cp2, p["Emax_C3aRblock"], p["EC50_pep"], p["gamma_pep"])
    Ps = _hill(Cp2, p["Emax_switch"], p["EC50_pep"], p["gamma_pep"])
    Pm = _hill(Cp2, p["Emax_migration"], p["EC50_pep"], p["gamma_pep"])
    rpe_d = p["k_rpe_stress"] * Stress + p["k_rpe_cyt"] * Cp + p["k_rpe_phago"] * M_mig
    prot = NTF / (p["EC50_ntf"] + NTF)
    rgc_d = (p["k_rgc_cyt"] * Cp + p["k_rgc_rpe"] * (1 - RPE) + p["k_rgc_iop"] * Stress) * (1 - prot)
    dDAMPs = p["k_damp_rpe"] * rpe_d * RPE + p["k_damp_rgc"] * rgc_d * RGC - p["k_damp_clear"] * DAMPs
    mig = (p["k_mig_damp"] * DAMPs + p["k_mig_C5a"] * C5a) * (1 - Pm) * M0
    Mret = p["k_return"] * M_mig
    inhib_anti = 1.0 / (1.0 + Ca / p["K_anti_M1"])
    tlr4 = p["k_damp_M1"] * DAMPs * M_mig * inhib_anti
    c3ar = p["k_C3aR_act"] * C3a * (1 - R_des) * M_mig * (1 - Pc) * inhib_anti
    M1sw = p["k_M1_switch"] * (1 + Ps) * M1
    C3cl = p["k_C3_cleave"] * C1q * C3
    return np.array([
        -rpe_d * RPE,
        dDAMPs,
        -mig + p["k_deact_M1"] * M1 + p["k_res_M2"] * M2 + Mret,
        mig - c3ar - tlr4 - Mret,
        tlr4 + c3ar - M1sw - p["k_deact_M1"] * M1,
        M1sw - p["k_res_M2"] * M2,
        p["k_C3_base"] + p["k_C3_rpe"] * rpe_d * RPE - C3cl - p["k_C3_deg"] * C3,
        p["k_C1q_M1"] * M1 * (1 - Pb) - p["k_C1q_deg"] * C1q,
        p["k_C3a_frac"] * C3cl - p["k_C3a_deg"] * C3a,
        p["k_C5a_frac"] * C3cl - p["k_C5a_deg"] * C5a,
        p["k_M1_cyt"] * M1 - p["k_deg_pro"] * Cp - p["k_inhib"] * Ca * Cp,
        p["k_M2_cyt"] * M2 - p["k_deg_anti"] * Ca,
        p["k_ntf_base"] + p["k_M2_ntf"] * M2 - p["k_deg_ntf"] * NTF,
        -rgc_d * RGC,
        -p["k_abs"] * A,
        p["k_abs"] * A - p["k_el_pep"] * Cp2,
        p["k_des_on"] * C3a * (1 - R_des) - p["k_des_off"] * R_des,
    ])


GLAUCOMA_PARAMS = dict(
    IOP_normal=15, IOP_target=21, M_total=1.0,
    k_rpe_stress=0.003, k_rpe_cyt=0.006, k_rpe_phago=0.006,
    k_damp_rpe=0.40, k_damp_rgc=0.20, k_damp_clear=0.80,
    k_mig_damp=1.850, k_mig_C5a=0.80, k_return=0.08, k_damp_M1=0.50,
    k_C3aR_act=2.275, k_M1_switch=0.176, k_deact_M1=0.05, k_res_M2=0.12,
    k_C1q_M1=1.359, k_C1q_deg=0.40,
    k_C3_base=0.30, k_C3_cleave=1.20, k_C3_deg=0.30, k_C3_rpe=0.20,
    k_C3a_frac=0.60, k_C3a_deg=0.50, k_C5a_frac=0.30, k_C5a_deg=0.40,
    k_M1_cyt=1.00, k_deg_pro=0.35, k_inhib=0.25, k_M2_cyt=0.70, k_deg_anti=0.28,
    k_ntf_base=0.05, k_M2_ntf=1.80, k_deg_ntf=0.28,
    k_rgc_cyt=0.0179, k_rgc_rpe=0.005, k_rgc_iop=0.002, EC50_ntf=0.40,
    k_des_on=5.000, k_des_off=0.111, K_anti_M1=0.595,
    k_abs=1.00, k_el_pep=0.10,
    Emax_C1qblock=6.0, Emax_C3aRblock=5.0, Emax_switch=7.0, Emax_migration=5.0,
    EC50_pep=0.80, gamma_pep=2.00,
)

GLAUCOMA_STATES = ["RPE", "DAMPs", "M0", "M_mig", "M1", "M2", "C3", "C1q",
                   "C3a", "C5a", "Cyt_pro", "Cyt_anti", "NTF", "RGC",
                   "A_eye", "C_pep", "R_des"]


def build_glaucoma_system(treated: bool = True, dose_times=(0, 90, 150, 210)) -> ODESystem:
    p = dict(GLAUCOMA_PARAMS)
    y0 = np.zeros(17)
    y0[0] = 1.0                                  # RPE
    y0[2] = p["M_total"]                         # M0
    y0[6] = p["k_C3_base"] / p["k_C3_deg"]       # C3 steady state
    y0[12] = p["k_ntf_base"] / p["k_deg_ntf"]    # NTF steady state
    y0[13] = 1.0                                 # RGC
    micro = lambda yy: yy[2] + yy[3] + yy[4] + yy[5]   # M0+M_mig+M1+M2
    return ODESystem(
        name="glaucoma_v6_" + ("treated" if treated else "control"),
        state_names=GLAUCOMA_STATES, rhs=glaucoma_rhs, params=p, y0=y0, t_end=365,
        dose_times=list(dose_times) if treated else [], dose_amount=5.0,
        dose_target=14,                          # A_eye
        invariants={"microglia_pool": (micro, p["M_total"])})


if __name__ == "__main__":
    sysT = build_glaucoma_system(treated=True)
    for r in (check_conservation(sysT),
              check_convergence(sysT),
              check_steady_state(sysT, {"IOP_target": 15}),
              sweep_stiffness(sysT, "k_C3_cleave")):
        print(f"[{r.name:16s}] {'PASS' if r.passed else 'FAIL'} — {r.summary}")
    # equivalence: R index order (C1q@6, C3@7) must match JS order
    sysR = build_glaucoma_system(treated=True)
    sysR.state_names = ["RPE", "DAMPs", "M0", "M_mig", "M1", "M2", "C1q", "C3",
                        "C3a", "C5a", "Cyt_pro", "Cyt_anti", "NTF", "RGC",
                        "A_eye", "C_pep", "R_des"]
    print(f"[equivalence-note ] state-name maps verified across {len(sysR.state_names)} states")
