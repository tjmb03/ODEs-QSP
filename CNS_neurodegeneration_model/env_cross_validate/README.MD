# Four-Engine Cross-Validation — Glaucoma v6 QSP ODE

![cross-engine-validation](https://github.com/tjmb03/ODEs-QSP/actions/workflows/crossval.yml/badge.svg)

> All model parameters are pseudodata generated programmatically for methodological
> demonstration. Ranges are literature-informed; no values derive from proprietary,
> internal-assay, or employer-affiliated sources.

A deterministic harness that confirms the **v6 17-state glaucoma QSP ODE** — the
RPE–DAMPs–complement–microglial vicious cycle with intravitreal peptide therapy — is
implemented **consistently across four independent stiff solvers** (R, Python, Julia,
MATLAB), plus an optional LLM agent for interactive triage when a check fails.

## What this establishes — and what it doesn't

This is a **solver-independence / numerical-reproducibility** check. It confirms that four
unrelated implementations integrate the *same* vector field and produce the *same*
trajectories to tight tolerance. It is **not** a claim about biological correctness, and it is
distinct from identifiability analysis (which asks whether parameters can be recovered).
Solver agreement means the math is implemented the same way everywhere — not that the model is
"true." If all four engines shared a bug they would agree; independence across methods and
codebases is what makes that vanishingly unlikely.

## The four engines

| Language | Solver | Method |
|---|---|---|
| R | `deSolve::lsoda` | ODEPACK / Fortran (Adams + BDF) |
| Python | `scipy` Radau | implicit Runge–Kutta |
| Julia | `OrdinaryDiffEq` FBDF | multistep BDF |
| MATLAB | `SimBiology` ode15s | NDF / BDF |

Two deliberate choices. Python uses **Radau, not LSODA** — LSODA shares Fortran lineage with
R's `deSolve`, so using it would make those two engines non-independent. And all four are
**stiff** solvers because the system spans fast complement kinetics and year-scale
neurodegeneration; an explicit method would be unstable without impractically small steps.

## Results

**Right-hand side** — verified **bit-identical** across all four implementations: `0.0`
difference at 2,000 random states, all **50 / 50** parameters matched. The engines integrate
the same function, not four lookalikes.

**Trajectory cross-validation** — max pairwise deviation across engines, gate tolerance `5e-3`:

| Arm | RGC | RPE | M1 | Verdict |
|---|---|---|---|---|
| Control (pure vector field) | `1.6e-8` | `7.8e-8` | `6.9e-8` | **PASS** |
| Treated (with dosing) | `1.2e-8` | `8.7e-9` | `5.4e-15` | **PASS** |

**Conservation** — the microglia pool (`M0 + M_mig + M1 + M2`) is preserved throughout:
`~1e-15` (R, Python) and `~8e-9` (Julia, at its solver tolerance). PASS for all six engine/arm
runs. Every readout sits four to five orders of magnitude inside the gate tolerance.

## How it works

Two arms isolate two concerns. The **control** arm runs with no dosing, so it tests the pure
vector field. The **treated** arm applies the dose schedule `{0, 90, 150, 210}` to the
intravitreal compartment, so it tests dosing semantics — the part easiest to get subtly wrong
across engines.

The **deterministic gate** (`run_crossval_gate.py`) runs each engine, compares trajectories
**by state name** (the engines order their states differently), and exits non-zero on any
disagreement past tolerance or any conservation failure. No LLM is involved; every number
comes from a solver.

The **optional LLM agent** (`validation_agent.py`) is for interactive triage only. It wires
the same deterministic tools as callable functions and orchestrates them — deciding which to
run and explaining the results — but it produces no numbers itself. The pattern is
*LLM → solver → LLM*: the model plans and interprets, the solvers compute. It is not in the CI
path and not required for the validation.

## A bug it caught

Solver-independence is only worth the effort if it surfaces real problems. It did. The gate
failed on the treated arm; the per-state diagnostic localized the disagreement to the
intravitreal peptide `A_eye` at **day 90**, where two engines differed by **exactly 5.0** — the
dose amount. That signature fingerprints a **pre- vs post-dose reporting mismatch**: one engine
reported the compartment value *before* the bolus was applied. Pinning every engine to a
segment-restart dose convention (dose applied at the segment boundary, reported post-dose)
closed it; the arms then agreed to `4e-12`.

## Running it

**Locally** — run each harness for both arms, then the gate:

```bash
python harness_python.py control python_control.csv
python harness_python.py treated python_treated.csv
Rscript harness_r.R control r_control.csv glaucoma_model.R
Rscript harness_r.R treated r_treated.csv glaucoma_model.R
julia  harness_julia.jl control julia_control.csv
julia  harness_julia.jl treated julia_treated.csv
# MATLAB, where SimBiology is licensed:
matlab -batch "harness_matlab('control','matlab_control.csv')"
matlab -batch "harness_matlab('treated','matlab_treated.csv')"
python run_crossval_gate.py --csv-dir .
```

**Docker** — a reproducible R + Python + Julia gate in a clean container (MATLAB excluded,
license-gated):

```bash
docker build -t glaucoma-xval .
docker run --rm glaucoma-xval
```

**CI** — the workflow at the repository root runs all four engines on every push (MATLAB
best-effort, license permitting) and gates on agreement.

## Files

| File | Role |
|---|---|
| `run_crossval_gate.py` | deterministic gate — per-engine conservation + pairwise cross-validation |
| `validation_tools.py` | Python kernel: the validated RHS, parameters, integrators, checks |
| `crossval.py` | trajectory comparator (interpolate to a common grid, compare by name) |
| `harness_python.py` · `harness_r.R` · `harness_julia.jl` · `harness_matlab.m` | one engine each |
| `glaucoma_model.R` | standalone R model core, sourced by the R harness |
| `build_glaucoma_simbiology.m` | programmatic SimBiology build (rate rules) |
| `diagnose_treated.py` | per-state divergence localizer (the triage tool) |
| `validation_agent.py` | optional LLM agent orchestrating the deterministic tools |
| `Dockerfile` | reproducible R + Python + Julia gate |
| `START_HERE.md` | runbook |

---

Part of the [CNS Neurodegeneration Model](../) · [QSPplus on YouTube](https://youtube.com/@QSPplus)
