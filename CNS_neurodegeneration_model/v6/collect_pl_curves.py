import json, sys, os
import numpy as np
import glaucoma_abc_v6 as G
import glaucoma_profile_likelihood_v6 as PL

OUT = "pl_curves_full.json"
p_true = dict(G.PARAM_DEFAULTS)
data, sigma, ref, n = PL.make_pseudodata(p_true, seed=0)
chi2_min, p_hat = PL.fit_free(data, sigma, p_true, fixed=None, maxiter=40)
db = json.load(open(OUT)) if os.path.exists(OUT) else {"chi2_min": chi2_min, "thresh": PL.CHI2_THRESH, "params": {}}
for name in sys.argv[1:]:
    grid, d2 = PL.profile_param(name, data, sigma, p_hat, chi2_min, n_grid=9, maxiter=18)
    lo, hi = PL.confidence_interval(grid, d2, p_hat[name])
    db["params"][name] = {"grid": grid.tolist(), "dchi2": d2.tolist(),
                          "theta_hat": p_hat[name], "ci_low": lo, "ci_high": hi,
                          "bounds": list(G.PRIORS[name]), "status": PL.classify(lo, hi)}
    print(f"  {name}: {PL.classify(lo,hi)}")
json.dump(db, open(OUT, "w"), indent=2)
