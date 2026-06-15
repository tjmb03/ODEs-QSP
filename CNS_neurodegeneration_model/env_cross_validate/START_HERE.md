# START HERE — glaucoma v6 cross-engine validation

A deterministic gate that confirms the glaucoma v6 ODE model is implemented
consistently across four independent stiff solvers (R, Python, MATLAB, Julia),
plus an optional LLM agent for interactive triage. **Build it up locally first;
CI is the last step.**

---

## File inventory

**Model core**
- `glaucoma_model.R` — the 17-state vector field + DEFAULTS + SNAMES, no Shiny
  deps. Your Shiny app should source this too (single source of truth).

**Engines** (each: `<harness> <control|treated> <out.csv>`)
- `harness_python.py` — Python, scipy Radau
- `harness_r.R`       — R, deSolve::lsoda  (sources a model file, arg 3)
- `harness_julia.jl`  — Julia, OrdinaryDiffEq FBDF
- `harness_matlab.m`  — MATLAB, SimBiology ode15s (needs `build_glaucoma_simbiology.m`)

**Comparator + gate** (deterministic, no LLM)
- `validation_tools.py` — single-engine numerics kernel (RK4, conservation, etc.)
- `crossval.py`         — cross-engine comparator (matches CSVs by column name)
- `run_crossval_gate.py`— the gate: conservation + cross-validation, exit codes

**Deployment**
- `requirements.txt`              — numpy, scipy (gate deps only)
- `Dockerfile`                    — reproducible R+Python+Julia image
- `.github/workflows/crossval.yml`— four-engine matrix → gate

**Interactive (separate, local only)**
- `validation_agent.py` — Claude Agent SDK orchestrator for adaptive triage

---

## Step 1 — two-engine local smoke test (≈10 min)

The minimum that exercises the full comparator, using the two engines you
already have. Put the model-core, comparator, gate, and the Python + R files in
one folder.

```bash
pip install -r requirements.txt
Rscript -e 'install.packages("deSolve", repos="https://cloud.r-project.org")'

python  harness_python.py control python_control.csv
python  harness_python.py treated python_treated.csv
Rscript harness_r.R      control r_control.csv glaucoma_model.R
Rscript harness_r.R      treated r_treated.csv glaucoma_model.R

python run_crossval_gate.py --csv-dir .
```

Expect: conservation PASS for both engines, cross-validation PASS on both arms,
exit 0. **If this is green, the architecture works.** If control disagrees, the
vector fields differ; if only treated disagrees, it's the dose handling.

## Step 2 — the trust anchor: does the model core match YOUR model?

Point the R harness at your real Shiny file and confirm it gives the identical
trajectory to `glaucoma_model.R`:

```bash
Rscript harness_r.R control real_control.csv glaucoma_qsp_v6.R
python -c "import crossval as cv, numpy as np; \
a=cv.load_trajectory('r_control.csv'); b=cv.load_trajectory('real_control.csv'); \
print('max diff:', max(float(np.max(np.abs(a[k]-b[k]))) for k in a))"
```

Near machine precision → the extracted core is faithful and you are validating
the right equations. (If `glaucoma_qsp_v6.R` calls `runApp()` at top level,
guard it with `if (interactive())` so sourcing doesn't launch the app.)

## Step 3 — add Julia, then MATLAB

```bash
julia -e 'using Pkg; Pkg.add(["OrdinaryDiffEq","DiffEqCallbacks"])'
julia harness_julia.jl control julia_control.csv
julia harness_julia.jl treated julia_treated.csv

matlab -batch "harness_matlab('control','matlab_control.csv')"   # SimBiology on path
matlab -batch "harness_matlab('treated','matlab_treated.csv')"

python run_crossval_gate.py --csv-dir .          # now 4 engines, pairwise
```

Independent stiff solvers should agree to ~1e-6 or better — once confirmed,
tighten `rtol` in `crossval.cross_validate` from 5e-3 toward 1e-3.

## Step 4 — containerize (reproducibility)

```bash
docker build -t glaucoma-xval .     # R+Python+Julia; first build is slow (Julia precompile)
docker run --rm glaucoma-xval       # runs the 3 in-container engines + the gate
```

## Step 5 — CI (the deployment)

Commit every file above (with `crossval.yml` under `.github/workflows/`) to the
repo. **Private repo:** request a MATLAB Batch Licensing token and add it as a
secret named `MLM_LICENSE_TOKEN`. **Public repo:** nothing to do. Push → the gate
runs on every PR; MATLAB is best-effort, R/Python/Julia are required.

---

## Optional — the interactive agent

Not part of CI (CI must be deterministic). For local "validate this and tell me
what's wrong" sessions: `pip install claude-agent-sdk`, `claude login`, then
`python validation_agent.py`. Uses your subscription; no API key.
