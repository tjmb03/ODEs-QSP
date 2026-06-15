> **Note:** All model parameters are pseudodata generated programmatically for
> methodological demonstration. Ranges are literature-informed; no values derive from
> proprietary, internal-assay, or employer-affiliated sources.

# CNS Neurodegeneration Model — Complement-Driven NTG (QSP)

A quantitative systems pharmacology model of the RPE–DAMPs–complement–microglial vicious
cycle in **normal-tension glaucoma**, with intravitreal peptide therapy and a full
practical + structural identifiability analysis.

## Versions

| Version | Folder | Highlights |
|---|---|---|
| **v6** (current) | [`v6/`](v6/) | 17-state **braked** ODE; RK4 / RK45 / ROS2 solvers; profile-likelihood + SIAN identifiability (**all parameters structurally identifiable**); 8 result & validation figures |
| v5.5 | [`v5.5/`](v5.5/) | 14-state ODE; in-browser RK4; original three-method identifiability framework |

The v6 folder contains a **"What's New in v6 (vs v5.5)"** changelog and a consolidated
validation summary. Start at [`v6/README.md`](v6/README.md).

## Cross-engine validation

[`env_cross_validate/`](env_cross_validate/) — a solver-independence harness for the v6
17-state ODE. The same system is integrated by **four independent stiff solvers** — R
`deSolve::lsoda`, Python `scipy Radau`, Julia `OrdinaryDiffEq FBDF`, and MATLAB
`SimBiology ode15s` — and trajectories are compared state-by-state (by name, since engines
order states differently). Two checks:

- **Right-hand side** — verified **bit-identical** across all four implementations
  (equations and all 50 parameters), so the engines integrate the same vector field, not
  four lookalikes.
- **Trajectories** — cross-validated on a no-drug (pure vector-field) arm and a dosed arm,
  dose schedule pinned identically across engines. Current agreement is below `1e-7`
  against a `5e-3` gate.

A deterministic CI gate ([`.github/workflows/crossval.yml`](../.github/workflows/crossval.yml))
re-runs the three open-source engines on every push and fails on disagreement; no LLM
produces any number (an optional `validation_agent.py` only orchestrates the same tools).
MATLAB/SimBiology is best-effort, license-gated on CI runners. Start at
[`env_cross_validate/START_HERE.md`](env_cross_validate/START_HERE.md).
