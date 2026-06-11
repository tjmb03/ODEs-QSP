"""
glaucoma_profile_likelihood_v6.py — practical identifiability of the v6 BRAKED model
====================================================================================

WHY THIS EXISTS
---------------
The ABC "constrainability" report (glaucoma_abc_v6.py) is NOT an identifiability
result: with an interval/range distance the posterior is just the prior truncated
to the feasible region. Worse, range-based targets make profile likelihood go
*spuriously flat* — fixing any parameter still lets the others hit the targets, so
every curve sits at zero and every parameter looks (falsely) non-identifiable.

The correct method (Raue et al. 2009; and the fix adopted in the GBM project) is
**profile likelihood against PSEUDO-DATA TRAJECTORY RESIDUALS**:

  1. Run the calibrated v6 model once → reference trajectories of MEASURABLE
     observables, sampled at several timepoints.
  2. Add measurement noise (known σ, 5% CV) → "pseudo-data". This makes
     χ²_min > 0, so the likelihood has a real scale.
  3. For each parameter θ_i: fix it on a grid, RE-OPTIMISE all other parameters
     to minimise χ², and record Δχ²(θ_i) = χ²_profile(θ_i) − χ²_min.
     Fixing θ_i wrong forces trajectory deviations the other parameters cannot
     fully absorb → a genuine U-shaped curve where the data constrain θ_i.
  4. Threshold at χ²_{0.95,1} = 3.84:
        crosses on both sides → IDENTIFIABLE (finite CI)
        crosses on one side   → PARTIALLY identifiable (one-sided bound)
        never crosses         → NON-IDENTIFIABLE (flat)

This is performed on the v6 *braked* equations (R_des desensitisation +
anti-inflammatory M1 suppression), which the earlier SimBiology profile likelihood
never saw — so it tells you whether the brakes changed the identifiability
structure found previously on the pre-brake (v5.5) model.

EPISTEMIC NOTE: the "data" are pseudo-data generated from the model itself, so this
measures *structural/practical identifiability of the model given a measurable
observable set and noise level* — it is NOT a fit to real experimental data (none
exists). It answers: "if you measured these observables at this precision, which
parameters could you pin down?"

Usage
-----
    python glaucoma_profile_likelihood_v6.py                       # all 11 params
    python glaucoma_profile_likelihood_v6.py --quick               # 3 params, coarse
    python glaucoma_profile_likelihood_v6.py --params k_rgc_rpe k_C3aR_act
Outputs: glaucoma_pl_v6_profiles.json, console identifiability table + ABC compare.
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np
from scipy.optimize import minimize
import glaucoma_abc_v6 as G

# ── Observables: measurable species (state index) sampled at these timepoints ──
OBS_TIMES = np.arange(30, 366, 30, dtype=float)        # monthly, days 30..360
OBS_SPECIES = {            # name : (state index, "pct" => ×100)
    "RGC": (13, "pct"),    # RGC survival (OCT / histology)
    "RPE": (0,  "pct"),    # RPE survival
    "M1":  (4,  "abs"),    # activated microglia (IBA1+/CD68+)
    "C1q": (6,  "abs"),    # complement C1q
    "C3a": (8,  "abs"),    # complement C3a
}
NOISE_CV = 0.05            # 5% measurement noise (known σ)
CHI2_THRESH = 3.8415       # χ²_{0.95, 1}
PL_PARAMS = list(G.PRIORS.keys())   # the 11 calibrated parameters


def observables(p) -> dict:
    """Control-arm observable trajectories at OBS_TIMES."""
    t, Y = G.simulate_control(p, t_end=365)
    out = {}
    for name, (idx, kind) in OBS_SPECIES.items():
        series = Y[idx]
        vals = np.interp(OBS_TIMES, t, series)
        out[name] = vals * 100.0 if kind == "pct" else vals
    return out


def make_pseudodata(p_true, seed=0):
    """Reference observables + Gaussian noise (known σ, relative floor)."""
    rng = np.random.default_rng(seed)
    ref = observables(p_true)
    data, sigma = {}, {}
    for name, vals in ref.items():
        floor = 0.05 * (np.max(np.abs(vals)) + 1e-9)
        sd = np.maximum(NOISE_CV * np.abs(vals), floor)
        data[name] = vals + rng.normal(0.0, sd)
        sigma[name] = sd
    n_obs = sum(len(v) for v in data.values())
    return data, sigma, ref, n_obs


def chi2(p, data, sigma) -> float:
    """Weighted SSE (χ²) between model observables and pseudo-data; robust to failures."""
    try:
        sim = observables(p)
    except Exception:
        return 1e12
    s = 0.0
    for name in data:
        r = (sim[name] - data[name]) / sigma[name]
        s += float(np.sum(r * r))
    return s if np.isfinite(s) else 1e12


def fit_free(data, sigma, p_start, fixed=None, maxiter=80):
    """Minimise χ² over all PL_PARAMS except `fixed` (held at p_start[fixed])."""
    free = [k for k in PL_PARAMS if k != fixed]
    x0 = np.array([p_start[k] for k in free])
    bounds = [G.PRIORS[k] for k in free]
    base = dict(p_start)

    def obj(x):
        p = dict(base)
        for k, v in zip(free, x):
            p[k] = v
        return chi2(p, data, sigma)

    res = minimize(obj, x0, method="L-BFGS-B", bounds=bounds,
                   options={"maxiter": maxiter, "maxfun": maxiter * 4, "ftol": 1e-9})
    p_fit = dict(base)
    for k, v in zip(free, res.x):
        p_fit[k] = v
    return res.fun, p_fit


def profile_param(name, data, sigma, p_hat, chi2_min, n_grid=15, maxiter=60):
    """Profile one parameter: fix on grid, re-optimise others (warm-started both ways)."""
    lo, hi = G.PRIORS[name]
    hat = p_hat[name]
    # grid split at θ̂ so we can warm-start outward in each direction
    left = np.linspace(lo, hat, max(2, n_grid // 2))[::-1]      # hat → lo
    right = np.linspace(hat, hi, max(2, n_grid // 2 + 1))        # hat → hi
    def sweep(grid):
        out, warm = [], dict(p_hat)
        for v in grid:
            warm[name] = v
            f, p_fit = fit_free(data, sigma, warm, fixed=name, maxiter=maxiter)
            warm = p_fit                                          # continuation
            out.append((float(v), float(f - chi2_min)))
        return out
    pts = sorted(set(sweep(left) + sweep(right)), key=lambda z: z[0])
    grid = np.array([g for g, _ in pts])
    d2 = np.array([d for _, d in pts])
    d2 = np.maximum(d2, 0.0)                                       # numerical guard
    return grid, d2


def confidence_interval(grid, d2, hat):
    """Where Δχ² crosses CHI2_THRESH on each side of the minimum (linear interp)."""
    imin = int(np.argmin(d2))
    def cross(side_g, side_d):
        for i in range(len(side_g) - 1):
            a, b = side_d[i], side_d[i + 1]
            if (a - CHI2_THRESH) * (b - CHI2_THRESH) < 0:         # sign change
                frac = (CHI2_THRESH - a) / (b - a)
                return side_g[i] + frac * (side_g[i + 1] - side_g[i])
        return None
    lo = cross(grid[imin::-1], d2[imin::-1])                      # toward small θ
    hi = cross(grid[imin:], d2[imin:])                            # toward large θ
    return lo, hi


def classify(lo, hi):
    if lo is not None and hi is not None:
        return "IDENTIFIABLE"
    if lo is not None or hi is not None:
        return "partial"
    return "NON-IDENTIFIABLE"


def main():
    ap = argparse.ArgumentParser(description="Pseudo-data profile likelihood, v6 braked model")
    ap.add_argument("--params", nargs="+", default=None, help="subset of parameters to profile")
    ap.add_argument("--grid", type=int, default=15)
    ap.add_argument("--maxiter", type=int, default=60)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--quick", action="store_true", help="3 params, coarse grid, tight budget")
    args = ap.parse_args()

    params = args.params or PL_PARAMS
    n_grid, maxiter = args.grid, args.maxiter
    if args.quick:
        params = params[:3] if args.params is None else params
        n_grid, maxiter = 7, 25

    # 1) pseudo-data from the calibrated reference point (v6 defaults)
    p_true = dict(G.PARAM_DEFAULTS)
    data, sigma, ref, n_obs = make_pseudodata(p_true, seed=args.seed)
    print(f"Pseudo-data: {len(OBS_SPECIES)} observables × {len(OBS_TIMES)} timepoints "
          f"= {n_obs} obs, {NOISE_CV*100:.0f}% CV noise")

    # 2) global best fit (all params free) → χ²_min, θ̂
    chi2_true = chi2(p_true, data, sigma)
    chi2_min, p_hat = fit_free(data, sigma, p_true, fixed=None, maxiter=maxiter * 2)
    dfe = n_obs - len(PL_PARAMS)
    print(f"χ²(p_true)={chi2_true:.1f}   χ²_min={chi2_min:.1f}   DFE={dfe}   "
          f"(χ²_min>0 ⇒ no spurious-flat plateau)\n")

    # 3) profile each parameter
    print(f"Profiling {len(params)} parameter(s): grid={n_grid}, maxiter={maxiter}"
          f"{'  [QUICK]' if args.quick else ''}")
    print("\n  parameter        θ̂        95% CI (profile)              status")
    print("  " + "-" * 78)
    results = {}
    for name in params:
        grid, d2 = profile_param(name, data, sigma, p_hat, chi2_min,
                                  n_grid=n_grid, maxiter=maxiter)
        lo, hi = confidence_interval(grid, d2, p_hat[name])
        status = classify(lo, hi)
        ci = (f"[{lo:.4g}, {hi:.4g}]" if status == "IDENTIFIABLE" else
              f"[{lo:.4g}, ∞)" if lo is not None else
              f"(-∞, {hi:.4g}]" if hi is not None else "(-∞, ∞)  flat")
        results[name] = {"theta_hat": p_hat[name], "ci_low": lo, "ci_high": hi,
                         "status": status, "max_dchi2": float(d2.max()),
                         "grid": grid.tolist(), "dchi2": d2.tolist()}
        print(f"  {name:14s} {p_hat[name]:8.4f}   {ci:28s}  {status} "
              f"(maxΔχ²={d2.max():.1f})")

    json.dump({"observables": list(OBS_SPECIES), "obs_times": OBS_TIMES.tolist(),
               "n_obs": n_obs, "chi2_min": chi2_min, "results": results},
              open("glaucoma_pl_v6_profiles.json", "w"), indent=2)
    print("\n  profiles written → glaucoma_pl_v6_profiles.json")

    # 4) cross-check against the ABC constrainability report, if present
    if os.path.exists("abc_posterior_summary.json"):
        abc = json.load(open("abc_posterior_summary.json"))["posterior"]
        print("\n  PROFILE LIKELIHOOD  vs  ABC constrainability (do they agree?):")
        print("    parameter        PL status            ABC sd/prior")
        print("    " + "-" * 60)
        for name in params:
            if name in abc:
                v = abc[name]; prior_sd = (v["hi"] - v["lo"]) / np.sqrt(12.0)
                ratio = v["sd"] / prior_sd
                print(f"    {name:14s} {results[name]['status']:18s}   {ratio:.2f} "
                      f"({'tight' if ratio<0.6 else 'weak' if ratio<0.85 else 'loose'})")
        print("\n  Disagreements are expected & informative: ABC constrainability reflects")
        print("  feasible-region width under the prior; PL reflects what the OBSERVABLES")
        print("  actually pin down. A parameter can look 'tight' to ABC yet be flat in PL")
        print("  (the targets constrain a combination the raw observables don't resolve).")

    print("\n  NOTE: pseudo-data identifiability — answers 'if these observables were")
    print("  measured at this precision, which parameters are pinned down', NOT a fit to")
    print("  real data. Re-run with richer/with-treatment observables to test designs.")


if __name__ == "__main__":
    main()
