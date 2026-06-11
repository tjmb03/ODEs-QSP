# v6 Braked Model — Identifiability Analysis (results, rev. 2)

Three methods on the **v6 braked** equations (R_des desensitisation +
anti-inflammatory M1 suppression), observable set {RGC, RPE, M1, C1q, C3a} sampled
monthly (12 timepoints), 5% CV. **Rev. 2** incorporates the SIAN algebraic proof,
which overturns the rev. 1 "rank 9/11 ⇒ structurally non-identifiable" reading.

**Epistemic status:** characterises the model + observation map (pseudo-data /
symbolic), not a fit to real data — none exists.

## Headline

| question | method | result |
|---|---|---|
| Are parameters identifiable *in principle*? | **structural** — SIAN `assess_local_identifiability` + generic-point sensitivity rank | **Yes — all 11 params + all 15 states locally structurally identifiable (rank 11/11)** |
| Identifiable *from this observable set + noise*? | **practical** — pseudo-data profile likelihood | 8 identifiable, 3 partial (the death-rate parameters) |
| Why the gap? | operating-point sensitivity rank | natural-history regime is collinear (Cytpro ≈ 1−RPE) → death directions ~1000× weak |

The death-pathway degeneracy is **practical (operating-regime collinearity + weak
excitation), NOT structural.** That distinction is decisive: a practical collinearity
is removable by experimental design; a structural one is not. **Confirmed globally:** on
the reduced subsystem with cytokine as a known input, all three death-rate parameters
are `:globally` identifiable and separate individually (§4).

---

## 1. Structural identifiability — DEFINITIVE (SIAN, algebraic)

`glaucoma_v6_SIAN.jl` → `assess_local_identifiability(ode)` returned **true for all 26
entries** (11 parameters + 15 states). The model is **locally structurally identifiable
and fully observable**: every parameter is uniquely recoverable from the five outputs
up to at most discrete ambiguities. No continuous non-identifiable direction exists.

Numerical confirmation (`glaucoma_structural_id_v6.py`): the normalised sensitivity
matrix has **rank 11/11 at every generic parameter point** (8/8 random draws), with all
singular values real and ≫ the 1e-8 finite-difference floor. The smallest two are only
~2×10⁻⁴ of σ_max — real, but weak.

The rev. 1 "9/11" was a **double artifact**: (a) evaluated at the non-generic
calibrated point, and (b) a 1e-3 cutoff that rounds the two smallest *real* singular
values to zero. Tolerance sweep at a generic point: rank = 9 at tol 1e-3, **11 at tol
≤ 1e-4** — i.e. the directions are present, just below 0.1% relative strength.

**Output-set ablation (vs v5.5).** Re-running SIAN with C3a demoted to a *hidden* state
and only {RGC, RPE, M1, C1q} observed — the v5.5 output set
(`glaucoma_v6_SIAN_4out.jl`) — still returns all-true. So the added C3a output is
sufficient but not necessary; the v6 model is structurally identifiable from the same
four outputs v5.5 reported as having non-identifiable parameters. Since structural
identifiability depends only on the equations, the resolution is attributable to the v6
structural change — the two saturating brakes (R_des + anti-inflammatory suppression),
not any added measurement.

Global (`assess_identifiability`) was killed (OOM) at the IO-equation step — the
documented Gröbner/differential-elimination wall for a 15-state rational system. Local
is the operative result; the global / discrete-alias check is in the reduced subsystem
below.

## 2. Practical identifiability — profile likelihood (pseudo-data)

`glaucoma_profile_likelihood_v6.py` · χ²_min = 45.7 (DFE 49) · threshold Δχ² = 3.84.

| Parameter | θ̂ | 95% CI (profile) | Status |
|---|---|---|---|
| k_C1q_M1    | 1.367  | [1.356, 1.396]  | identifiable |
| k_C3_cleave | 1.204  | [1.190, 1.213]  | identifiable |
| k_C3aR_act  | 2.286  | [2.178, 2.555]  | identifiable |
| k_damp_M1   | 0.502  | [0.380, 0.558]  | identifiable |
| k_mig_damp  | 1.856  | [1.832, 1.865]  | identifiable |
| k_M1_switch | 0.175  | [0.173, 0.177]  | identifiable |
| k_des_on    | 4.984  | [4.941, 5.031]  | identifiable |
| K_anti_M1   | 0.597  | [0.533, 0.784]  | identifiable |
| k_rpe_cyt   | 0.0028 | (−∞, 0.00928]   | partial (upper) |
| k_rgc_rpe   | 0.0040 | (−∞, 0.00573]   | partial (upper) |
| k_rgc_cyt   | 0.0212 | [0.01225, ∞)    | partial (lower) |

8 identifiable · 3 partial. The 3 partials are exactly the weak directions from §1.

## 3. The reconciliation

All three methods agree once read correctly:

- **Structurally** (SIAN, generic rank): all 11 identifiable — `k_rgc_cyt` multiplies
  `Cytpro(t)` and `k_rgc_rpe` multiplies `(1−RPE(t))`, *distinct* signals, algebraically
  separable.
- **Operating-point rank**: the two weakest SVD directions are
  `−0.70·k_rgc_cyt +0.63·k_rgc_rpe` and a `K_anti_M1`/death-rate mix — because in the
  natural-history trajectory `Cytpro` and `(1−RPE)` rise together (collinear), so those
  terms are nearly proportional *there*.
- **Profile likelihood**: the same parameters get one-sided CIs — finite noise cannot
  resolve a near-collinear pair.

→ A **practical** degeneracy of the operating regime, not a property of the model.

## 4. Global / discrete-alias check — reduced subsystem

`glaucoma_v6_SIAN_reduced.jl`: states {RPE, RGC}; cytokine, migratory microglia and
NTF as **known inputs** (the experimental condition that decouples cytokine from RPE
loss); outputs RGC, RPE; parameters `k_rpe_cyt, k_rgc_cyt, k_rgc_rpe`.

**Run (16 s) — confirmed:** `assess_identifiability` returns **all three `:globally`
identifiable** (both states globally observable too) — global, so no discrete aliases
either. And `find_identifiable_functions` returns the **bare parameters**
`k_rpe_cyt, k_rgc_rpe, k_rgc_cyt` — *not* combinations. That is the strongest possible
statement: the three do not confound each other at all once cytokine is a known input.
(Contrast the textbook case where this function returns combinations like `a+b`, `a·b`,
which would signal the individual parameters are *not* separable.)

→ The death-pathway degeneracy is **entirely an artifact of the observation condition**
(cytokine hidden and co-varying with RPE loss in natural history); it is removed by an
experiment that measures or perturbs cytokine independently of RPE loss.

## 5. Why this differs from the earlier (v5.5) profile likelihood

The earlier analysis used a **single endpoint** (t=180) of raw species and found most
parameters flat. Here, **time-series observables (12 timepoints + C3a)** make the
signalling block identifiable and leave only the death-rate collinearity. So the
earlier non-identifiability was mostly **measurement design (single endpoint)**, not
structure — confirmed now from two directions (trajectory PL + SIAN).

## 6. Experimental implication (corrected)

The model is identifiable **in principle** (SIAN). The barrier is practical: the three
death-rate parameters are collinear in the natural-history regime. The fix is **not**
more survival curves — it is an experiment that **breaks the Cytpro ↔ RPE-loss
collinearity**: a cytokine-neutralisation / C3aR-block arm (drives `Cytpro` down while
IOP-driven RPE loss continues), or an RPE-protection arm (isolates `k_rgc_rpe`). The
reduced subsystem shows algebraically that with cytokine measured independently the
three parameters separate. This is the top identifiability-driven experimental
priority, distinct from the sensitivity-driven priorities in the dashboard roadmap.

---

**Artifacts:** `glaucoma_profile_likelihood_v6.py`, `glaucoma_structural_id_v6.py`,
`glaucoma_structural_id_v6.json`, `glaucoma_v6_SIAN.jl` (local proof, run),
`glaucoma_v6_SIAN_reduced.jl` (global check).
