"""
glaucoma_abc_v6.py — ABC-SMC plausibility calibration of the v6 *braked* QSP model
==================================================================================

EPISTEMIC STATUS — READ THIS FIRST
----------------------------------
This is NOT Bayesian inference from primary experimental data, and the output is
NOT a data-constrained posterior. No primary kinetic measurements exist for the
target rate constants: a systematic literature extraction across ~44 papers
returned a consistent null for all of them. What this script performs is
*plausibility-based virtual-population calibration*:

    Given literature-plausible prior ranges, find the region of parameter space
    whose emergent behaviour is consistent with a set of QUALITATIVE biological
    TARGET INTERVALS for chronic normal-tension glaucoma (NTG).

The "observed data" are interval targets on emergent model behaviour (RGC/RPE
survival, microglial activation, loop-ignition timing, death-mechanism split),
not measured trajectories. The ABC posterior is therefore a *plausible region*,
not an inferred distribution. Expect several parameters to come out essentially
non-identifiable (posterior ≈ prior) — that is the correct, honest result and is
flagged explicitly in the constrainability report at the end (which is NOT an
identifiability result — see glaucoma_profile_likelihood_v6.py for that).

Model: the v6 17-state complement–microglia ODE with the two saturating brakes
(C3aR desensitisation R_des, anti-inflammatory M1 suppression via Cyt_anti).
Calibration is performed on the CONTROL arm (no drug), which sets the disease
phenotype. Drug-arm parameters (Emax_*, EC50_pep, k_el_pep) are not calibrated
here.

Usage
-----
    pip install pyabc scipy numpy
    python glaucoma_abc_v6.py --particles 256 --generations 8       # production-ish
    python glaucoma_abc_v6.py --particles 48  --generations 3 --quick  # smoke test

Outputs: glaucoma_abc_v6_posterior.csv (weighted posterior sample),
         console summary + constrainability report (NOT identifiability).
"""
from __future__ import annotations
import argparse, os, tempfile, sys
from dataclasses import dataclass, field
import numpy as np
from scipy.integrate import solve_ivp

# ─────────────────────────────────────────────────────────────────────────────
# 1. v6 model defaults (the hand-tuned 6/6 operating point — used to fill in all
#    parameters NOT being calibrated, and as the reference for the posterior).
# ─────────────────────────────────────────────────────────────────────────────
PARAM_DEFAULTS = dict(
    IOP_normal=15.0, IOP_target=21.0, M_total=1.0,
    k_rpe_stress=0.003, k_rpe_cyt=0.006, k_rpe_phago=0.006,
    k_damp_rpe=0.40, k_damp_rgc=0.20, k_damp_clear=0.80,
    k_mig_damp=1.85, k_mig_C5a=0.80, k_return=0.08,
    k_damp_M1=0.50, k_C3aR_act=2.28, k_M1_switch=0.176, k_deact_M1=0.05, k_res_M2=0.12,
    k_C1q_M1=1.36, k_C1q_deg=0.40, k_C3_base=0.30, k_C3_cleave=1.20, k_C3_deg=0.30, k_C3_rpe=0.20,
    k_C3a_frac=0.60, k_C3a_deg=0.50, k_C5a_frac=0.30, k_C5a_deg=0.40,
    k_M1_cyt=1.0, k_deg_pro=0.35, k_inhib=0.25, k_M2_cyt=0.70, k_deg_anti=0.28,
    k_ntf_base=0.05, k_M2_ntf=1.80, k_deg_ntf=0.28,
    k_rgc_cyt=0.018, k_rgc_rpe=0.005, k_rgc_iop=0.002, EC50_ntf=0.40,
    k_des_on=5.0, k_des_off=0.111, K_anti_M1=0.595,
    k_abs=1.0, k_el_pep=0.10,
    Emax_C1qblock=6.0, Emax_C3aRblock=5.0, Emax_switch=7.0, Emax_migration=5.0,
    EC50_pep=0.80, gamma_pep=2.0, dose_amount=5.0,
)

# ─────────────────────────────────────────────────────────────────────────────
# 2. Priors — literature-plausible ranges (uniform) bracketing the v6 defaults.
#    These are the parameters with genuine uncertainty (no primary measurements)
#    that shape the chronic-NTG phenotype. Seven inherited from the v5.5 set plus
#    four v6-specific drivers (k_mig_damp, k_M1_switch + the two brake parameters).
# ─────────────────────────────────────────────────────────────────────────────
PRIORS = {
    # name           (low,   high)     # v6 default  | role
    "k_C1q_M1":      (0.50,  2.50),     # 1.36   M1 → C1q production
    "k_C3_cleave":   (0.50,  2.00),     # 1.20   C1q-driven C3 cleavage
    "k_C3aR_act":    (1.00,  3.50),     # 2.28   C3a → C3aR feedback (loop closure)
    "k_damp_M1":     (0.20,  1.00),     # 0.50   DAMP → M1 bootstrap (TLR4/NF-κB)
    "k_rpe_cyt":     (0.002, 0.015),    # 0.006  cytokine → RPE death
    "k_rgc_rpe":     (0.002, 0.012),    # 0.005  RPE loss → trans-synaptic RGC death
    "k_rgc_cyt":     (0.008, 0.035),    # 0.018  cytokine → RGC death
    "k_mig_damp":    (1.00,  3.00),     # 1.85   DAMP → microglial migration
    "k_M1_switch":   (0.08,  0.30),     # 0.176  M1 → M2 resolution rate
    "k_des_on":      (2.00, 10.00),     # 5.00   BRAKE 1: C3aR desensitisation on-rate
    "K_anti_M1":     (0.30,  1.00),     # 0.595  BRAKE 2: anti-inflammatory M1 suppression K
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. Calibration targets — interval constraints on six emergent CONTROL-arm
#    endpoints for chronic NTG. (lo, hi, scale); scale = half-width, so a
#    one-half-width violation contributes 1.0 to the squared distance.
# ─────────────────────────────────────────────────────────────────────────────
TARGETS = {
    "RGC_pct":      (45.0, 62.0,  8.5),   # RGC survival at 1 yr
    "RPE_pct":      (25.0, 65.0, 20.0),   # RPE survival at 1 yr
    "M1_peak":      (5.0,  15.0,  5.0),   # peak activated-microglia fraction
    "ignition_day": (90.0, 180.0, 45.0),  # day the complement loop ignites
    "cyt_frac":     (20.0, 45.0, 12.5),   # cytokine share of RGC death rate
    "trans_frac":   (20.0, 55.0, 17.5),   # trans-synaptic share of RGC death rate
}
# Observed dict (target midpoints) — used by pyABC only for summary-stat KEYS;
# the interval logic lives in `distance`, not in these point values.
OBSERVED = {k: 0.5 * (lo + hi) for k, (lo, hi, _) in TARGETS.items()}


# ─────────────────────────────────────────────────────────────────────────────
# 4. Forward model: the v6 17-state braked ODE (control arm).
# ─────────────────────────────────────────────────────────────────────────────
def _hill(C, Emax, EC50, g):
    if C <= 0:
        return 0.0
    return Emax * C**g / (EC50**g + C**g)


def _odes(t, y, p):
    (RPE, DAMPs, M0, M_mig, M1, M2, C1q, C3, C3a, C5a,
     Cyt_pro, Cyt_anti, NTF, RGC, A, Cp, R_des) = [max(v, 0.0) for v in y]
    RPE = min(RPE, 1.0); RGC = min(RGC, 1.0); R_des = min(R_des, 1.0)
    S = max(0.0, (p['IOP_target'] - p['IOP_normal']) / p['IOP_normal'])
    # drug PD (zero on the control arm because Cp stays 0)
    PC1q  = _hill(Cp, p['Emax_C1qblock'],  p['EC50_pep'], p['gamma_pep'])
    PC3aR = _hill(Cp, p['Emax_C3aRblock'], p['EC50_pep'], p['gamma_pep'])
    Psw   = _hill(Cp, p['Emax_switch'],    p['EC50_pep'], p['gamma_pep'])
    Pmig  = _hill(Cp, p['Emax_migration'], p['EC50_pep'], p['gamma_pep'])
    rpe_d = p['k_rpe_stress']*S + p['k_rpe_cyt']*Cyt_pro + p['k_rpe_phago']*M_mig
    dRPE  = -rpe_d * RPE
    prot  = NTF / (p['EC50_ntf'] + NTF)
    rgc_d = (p['k_rgc_cyt']*Cyt_pro + p['k_rgc_rpe']*(1 - RPE) + p['k_rgc_iop']*S) * (1 - prot)
    dRGC  = -rgc_d * RGC
    dDAMPs = p['k_damp_rpe']*rpe_d*RPE + p['k_damp_rgc']*rgc_d*RGC - p['k_damp_clear']*DAMPs
    mig   = (p['k_mig_damp']*DAMPs + p['k_mig_C5a']*C5a) * (1 - Pmig) * M0
    Mret  = p['k_return'] * M_mig
    dM0   = -mig + p['k_deact_M1']*M1 + p['k_res_M2']*M2 + Mret
    inhib = 1.0 / (1.0 + Cyt_anti / p['K_anti_M1'])          # BRAKE 2
    tlr4  = p['k_damp_M1']*DAMPs*M_mig*inhib
    c3ar  = p['k_C3aR_act']*C3a*(1 - R_des)*M_mig*(1 - PC3aR)*inhib   # (1-R_des) = BRAKE 1
    M1sw  = p['k_M1_switch']*(1 + Psw)*M1
    dM_mig = mig - c3ar - tlr4 - Mret
    dM1   = tlr4 + c3ar - M1sw - p['k_deact_M1']*M1
    dM2   = M1sw - p['k_res_M2']*M2
    dR_des = p['k_des_on']*C3a*(1 - R_des) - p['k_des_off']*R_des    # BRAKE 1 state
    dC1q  = p['k_C1q_M1']*M1*(1 - PC1q) - p['k_C1q_deg']*C1q
    C3cl  = p['k_C3_cleave']*C1q*C3
    dC3   = p['k_C3_base'] + p['k_C3_rpe']*rpe_d*RPE - C3cl - p['k_C3_deg']*C3
    dC3a  = p['k_C3a_frac']*C3cl - p['k_C3a_deg']*C3a
    dC5a  = p['k_C5a_frac']*C3cl - p['k_C5a_deg']*C5a
    dCp_pro = p['k_M1_cyt']*M1 - p['k_deg_pro']*Cyt_pro - p['k_inhib']*Cyt_anti*Cyt_pro
    dCa   = p['k_M2_cyt']*M2 - p['k_deg_anti']*Cyt_anti
    dNTF  = p['k_ntf_base'] + p['k_M2_ntf']*M2 - p['k_deg_ntf']*NTF
    dA    = -p['k_abs']*A
    dCp   = p['k_abs']*A - p['k_el_pep']*Cp
    return [dRPE, dDAMPs, dM0, dM_mig, dM1, dM2, dC1q, dC3, dC3a, dC5a,
            dCp_pro, dCa, dNTF, dRGC, dA, dCp, dR_des]


def simulate_control(p, t_end=365.0):
    """Integrate the control arm (no drug). Returns (t, Y) or raises on failure."""
    C3ss = p['k_C3_base'] / p['k_C3_deg']
    NTFss = p['k_ntf_base'] / p['k_deg_ntf']
    y0 = [1, 0, p['M_total'], 0, 0, 0, 0, C3ss, 0, 0, 0, 0, NTFss, 1, 0, 0, 0]
    sol = solve_ivp(_odes, (0, t_end), y0, args=(p,),
                    t_eval=np.linspace(0, t_end, int(t_end) + 1),
                    method='LSODA', rtol=1e-6, atol=1e-8)
    if not sol.success:
        raise RuntimeError("integration failed")
    return sol.t, sol.y


# ─────────────────────────────────────────────────────────────────────────────
# 5. Summary statistics — the six emergent endpoints (control arm).
# ─────────────────────────────────────────────────────────────────────────────
def summary_stats(p) -> dict:
    try:
        t, Y = simulate_control(p)
    except Exception:
        return {k: np.nan for k in TARGETS}
    RGC, RPE, M1, C1q, C3a, Cyt_pro = Y[13], Y[0], Y[4], Y[6], Y[8], Y[10]
    S = max(0.0, (p['IOP_target'] - p['IOP_normal']) / p['IOP_normal'])
    loop = M1 * C1q * C3a
    pk = float(loop.max())
    if pk > 1e-9:
        idx = np.where(loop > 0.01 * pk)[0]
        ign = float(t[idx[0]]) if idx.size else float(t[-1])
    else:
        ign = float(t[-1])          # loop never ignites
    d_cyt = p['k_rgc_cyt'] * Cyt_pro[-1]
    d_rpe = p['k_rgc_rpe'] * (1 - RPE[-1])
    d_iop = p['k_rgc_iop'] * S
    tot = d_cyt + d_rpe + d_iop + 1e-12
    return {
        "RGC_pct":      float(RGC[-1] * 100),
        "RPE_pct":      float(RPE[-1] * 100),
        "M1_peak":      float(M1.max() * 100),
        "ignition_day": ign,
        "cyt_frac":     float(d_cyt / tot * 100),
        "trans_frac":   float(d_rpe / tot * 100),
    }


def model(params: dict) -> dict:
    """pyABC model callable: merge sampled params into defaults, return stats."""
    p = dict(PARAM_DEFAULTS); p.update(params)
    return summary_stats(p)


# ─────────────────────────────────────────────────────────────────────────────
# 6. Distance — two modes (the choice changes what the posterior MEANS):
#
#   "interval"  Zero whenever every endpoint lies inside its target band; squared
#               excess outside. This is the HONEST representation of "targets are
#               ranges" — BUT it has a zero-distance plateau covering the whole
#               feasible region, so the posterior collapses to *the prior truncated
#               to that region*. Posterior width then reflects feasible-region
#               geometry, NOT data constraint. (See the constrainability caveat in
#               `analyze()`.) This is the same failure mode that makes range-target
#               profile likelihood go spuriously flat — do NOT read identifiability
#               off it.
#
#   "smooth"    Quadratic distance to each target MIDPOINT, scaled by half-width.
#               Treats the targets as point estimates with uncertainty ≈ half-width.
#               The posterior is then data-shaped (no plateau), so posterior width is
#               more interpretable — at the cost of treating the band centre as
#               "preferred", which is a mild distortion if the targets are truly flat
#               ranges. Use this if you want a non-degenerate posterior to summarise.
#
# NEITHER mode yields identifiability in the standard sense; for that use the
# pseudo-data profile likelihood (glaucoma_profile_likelihood_v6.py).
# ─────────────────────────────────────────────────────────────────────────────
DISTANCE_MODE = "interval"     # set by run_calibration() from ABCConfig

def distance(x: dict, x_0: dict) -> float:
    d = 0.0
    for key, (lo, hi, scale) in TARGETS.items():
        v = x.get(key, np.nan)
        if v is None or not np.isfinite(v):
            return 1e6
        if DISTANCE_MODE == "smooth":
            mid = 0.5 * (lo + hi)
            d += ((v - mid) / scale) ** 2
        else:  # "interval"
            if v < lo:
                d += ((lo - v) / scale) ** 2
            elif v > hi:
                d += ((v - hi) / scale) ** 2
    return float(np.sqrt(d))


# ─────────────────────────────────────────────────────────────────────────────
# 7. Run + analysis.
# ─────────────────────────────────────────────────────────────────────────────
@dataclass
class ABCConfig:
    particles: int = 256
    generations: int = 8
    min_epsilon: float = 0.25      # accept when within ~0.25 half-widths of all targets
    quick: bool = False
    distance_mode: str = "interval"   # "interval" (honest for ranges) | "smooth" (data-shaped)
    out_csv: str = "glaucoma_abc_v6_posterior.csv"
    db_path: str = field(default_factory=lambda: os.path.join(tempfile.gettempdir(),
                                                              "glaucoma_abc_v6.db"))


def run_calibration(cfg: ABCConfig):
    import pyabc
    global DISTANCE_MODE
    DISTANCE_MODE = cfg.distance_mode
    priors = pyabc.Distribution(**{
        name: pyabc.RV("uniform", lo, hi - lo) for name, (lo, hi) in PRIORS.items()
    })
    sampler = pyabc.sampler.SingleCoreSampler() if cfg.quick else pyabc.sampler.MulticoreEvalParallelSampler()
    abc = pyabc.ABCSMC(model, priors, distance,
                       population_size=cfg.particles, sampler=sampler)
    if os.path.exists(cfg.db_path):
        os.remove(cfg.db_path)
    abc.new("sqlite:///" + cfg.db_path, OBSERVED)
    history = abc.run(minimum_epsilon=cfg.min_epsilon,
                      max_nr_populations=cfg.generations)
    return history


def analyze(history, cfg: ABCConfig):
    df, w = history.get_distribution()           # final-generation weighted sample
    print("\n" + "=" * 74)
    print("ABC-SMC RESULT — v6 braked model (plausibility calibration)")
    print("=" * 74)
    print(f"generations run : {history.max_t + 1}")
    print(f"final population : {len(df)} particles")
    try:
        eps = history.get_all_populations()['epsilon'].values
        print(f"epsilon schedule: {', '.join(f'{e:.3f}' for e in eps)}")
    except Exception:
        pass

    # Posterior means + CONSTRAINABILITY (posterior std vs prior std).
    # IMPORTANT: this is NOT identifiability. With the interval distance the
    # posterior is the prior truncated to the feasible region, so this ratio
    # measures how tightly the *plausibility targets* pin each parameter
    # (confounded by prior width) — not whether data could identify it.
    print("\n  parameter        prior range        post.mean   post.sd    v6 default   "
          "constrainability (by targets, NOT identifiability)")
    print("  " + "-" * 108)
    rows = []
    for name, (lo, hi) in PRIORS.items():
        col = df[name].values
        mean = float(np.average(col, weights=w))
        var = float(np.average((col - mean) ** 2, weights=w))
        sd = var ** 0.5
        prior_sd = (hi - lo) / np.sqrt(12.0)          # SD of a uniform prior
        ratio = sd / prior_sd                          # <0.6 tightly, ~1 not constrained by targets
        tag = ("tightly constrained" if ratio < 0.6 else
               "weakly constrained" if ratio < 0.85 else "unconstrained by targets")
        rows.append((name, mean, sd))
        print(f"  {name:14s} [{lo:6.3f},{hi:6.3f}]   {mean:9.4f}  {sd:8.4f}   "
              f"{PARAM_DEFAULTS[name]:9.4f}    {tag} (sd/prior={ratio:.2f})")

    # Endpoints at the posterior mean
    pm = dict(PARAM_DEFAULTS); pm.update({n: m for n, m, _ in rows})
    s = summary_stats(pm)
    print("\n  posterior-mean endpoints vs targets:")
    for k, (tlo, thi, _) in TARGETS.items():
        ok = "✓" if tlo <= s[k] <= thi else "✗ OUTSIDE"
        print(f"    {k:13s} = {s[k]:7.2f}   target [{tlo:.0f}, {thi:.0f}]   {ok}")

    # Save weighted posterior sample
    out = df.copy()
    out["weight"] = w
    out.to_csv(cfg.out_csv, index=False)
    print(f"\n  weighted posterior sample written → {cfg.out_csv}")
    print("\n  CAVEAT — this table is CONSTRAINABILITY BY THE PLAUSIBILITY TARGETS, not")
    print("  identifiability. The interval distance has a zero-distance plateau over the")
    print("  whole feasible region, so the posterior ≈ prior truncated to it; 'sd/prior'")
    print("  reflects feasible-region geometry (and prior width), NOT data constraint.")
    print("  For genuine identifiability run: glaucoma_profile_likelihood_v6.py")
    print("  (pseudo-data profile likelihood on the v6 braked equations).")
    return out


def self_test():
    """Sanity check without running ABC: defaults should satisfy every target."""
    s = summary_stats(PARAM_DEFAULTS)
    d = distance(s, OBSERVED)
    print("Self-test — endpoints at v6 defaults:")
    allok = True
    for k, (lo, hi, _) in TARGETS.items():
        ok = lo <= s[k] <= hi
        allok &= ok
        print(f"  {k:13s} = {s[k]:7.2f}   target [{lo:.0f}, {hi:.0f}]   {'✓' if ok else '✗'}")
    print(f"  distance(defaults) = {d:.4f}   (0 = inside every target band)")
    print(f"  => defaults satisfy all targets: {allok}")
    return allok


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="ABC-SMC plausibility calibration, v6 braked QSP model")
    ap.add_argument("--particles", type=int, default=256)
    ap.add_argument("--generations", type=int, default=8)
    ap.add_argument("--min-epsilon", type=float, default=0.25)
    ap.add_argument("--distance", choices=["interval", "smooth"], default="interval",
                    help="interval = honest for ranges (posterior≈prior-truncated); "
                         "smooth = quadratic-to-midpoint (data-shaped posterior)")
    ap.add_argument("--quick", action="store_true", help="single-core sampler, for a fast smoke test")
    ap.add_argument("--self-test", action="store_true", help="check defaults vs targets and exit")
    args = ap.parse_args()

    if args.self_test:
        sys.exit(0 if self_test() else 1)

    if not self_test():
        print("\n[warn] defaults do not satisfy all targets — calibration will still run.\n")
    cfg = ABCConfig(particles=args.particles, generations=args.generations,
                    min_epsilon=args.min_epsilon, quick=args.quick,
                    distance_mode=args.distance)
    print(f"\nRunning ABC-SMC: {cfg.particles} particles × up to {cfg.generations} "
          f"generations (min ε={cfg.min_epsilon}){' [quick/single-core]' if cfg.quick else ''} …")
    hist = run_calibration(cfg)
    analyze(hist, cfg)
