# ODEs-QSP

> **A collection of quantitative systems pharmacology (QSP) models built as fully interactive, in-browser dashboards — no installation required.**

Each model implements a disease-specific ODE system with a paired therapeutic intervention, simulated in real time.

---

## Projects

### 🧠 [Glioblastoma (GBM) — rAAV/shRNA Gene Therapy](./Glioblastoma)

Models the **DNA-PK/MYT1L-CXCR1-ERK1/2 positive feedback loop** in GBM and its disruption by rAAV-delivered shRNA. The loop drives autocrine CXCR1 amplification and ERK-mediated proliferation — a self-sustaining vicious cycle targeted by post-transcriptional gene silencing.

| Feature | Detail |
|---------|--------|
| ODE system | 13 state variables |
| Therapy | rAAV-shCXCR1 / shDNAPK / shMYT1L |
| Key loops | DNA-PK/MYT1L-CXCR1, TAM-TGFβ immune escape, CTL protective |
| Burn-in | 90-day pre-treatment equilibration |
| Outputs | GBM burden, loop suppression, ERK/DNAPK/MYT1L cascade, immune compartment |

---

### 👁 [Glaucoma — Complement-Driven Neurodegeneration](./Glaucoma)

Models the **RPE–DAMPs–complement–microglial vicious cycle** in normal-tension glaucoma (NTG) and its disruption by an intravitreal peptide with three simultaneous PD mechanisms.
Quantifies the drug-accessible cytokine death fraction versus the structural trans-synaptic fraction that defines the therapeutic efficacy ceiling.

| Feature | Detail |
|---------|--------|
| ODE system | 14-state RK4 (RPE, DAMPs, M0/M_mig/M1/M2, C3/C1q/C3a/C5a, Cyt_pro/Cyt_anti, NTF, RGC) + 2 PK compartments |
| Feedback loop | M1 → C1q → C4b2a → C3a (C3aR, closes loop) + C5a (C5aR1, drives migration) → M1 |
| Loop gain | L = ∏k / ∏δ around loop; calibrated L ≈ 5.5× (supercritical, self-sustaining) |
| Therapy | Intravitreal peptide: C1q block (PD1) · C3aR/C5aR antagonism (PD2) · M1→M2 switch (PD3) |
| Key insight | ~42–50% of RGC death is trans-synaptic (structural, drug-independent); ~46–55% is cytokine-mediated (drug-accessible) — sets hard ceiling on achievable protection |
| Outputs | RGC/RPE survival, death pathway decomposition, loop intensity, microglial polarisation, NTF neuroprotection, +11–12 pp protection at 180 days |

---

## Repository Structure

```
ODEs-QSP/
├── README.md               ← this file
├── Glioblastoma/
│   ├── figures/            ← simulation output PNGs
│   └── README.md           ← model documentation & results
└── Glaucoma/
    ├── figures/            ← simulation output PNGs
    └── README.md           ← model documentation & results
```


---

## Contact

**Author:** tjmb03  


> © 2026 tjmb03. Source code for the interactive dashboards is **available on request** for academic and research use.

