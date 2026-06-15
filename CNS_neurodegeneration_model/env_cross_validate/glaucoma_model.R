# glaucoma_model.R — standalone model CORE for the glaucoma QSP v6 system.
#
# Contains ONLY the ODE vector field (glaucoma_odes), default parameters
# (DEFAULTS), and state names (SNAMES) — no Shiny / UI dependencies — so it can
# be sourced by harness_r.R in CI without pulling the full app stack.
#
# This is the same 17-state vector field validated term-by-term against the JS
# and Python implementations. Your Shiny app (glaucoma_qsp_v6.R) should source
# THIS file too, to keep a single source of truth. Diff it against your existing
# glaucoma_odes once to confirm they match.

SNAMES <- c("RPE","DAMPs","M0","M_mig","M1","M2","C1q","C3","C3a","C5a",
            "Cyt_pro","Cyt_anti","NTF","RGC","A_eye","C_pep","R_des")

DEFAULTS <- list(
  IOP_normal = 15, IOP_target = 21, M_total = 1.0,
  k_rpe_stress = 0.003, k_rpe_cyt = 0.006, k_rpe_phago = 0.006,
  k_damp_rpe = 0.40, k_damp_rgc = 0.20, k_damp_clear = 0.80,
  k_mig_damp = 1.850, k_mig_C5a = 0.80, k_return = 0.08, k_damp_M1 = 0.50,
  k_C3aR_act = 2.275, k_M1_switch = 0.176, k_deact_M1 = 0.05, k_res_M2 = 0.12,
  k_C1q_M1 = 1.359, k_C1q_deg = 0.40,
  k_C3_base = 0.30, k_C3_cleave = 1.20, k_C3_deg = 0.30, k_C3_rpe = 0.20,
  k_C3a_frac = 0.60, k_C3a_deg = 0.50, k_C5a_frac = 0.30, k_C5a_deg = 0.40,
  k_M1_cyt = 1.00, k_deg_pro = 0.35, k_inhib = 0.25,
  k_M2_cyt = 0.70, k_deg_anti = 0.28,
  k_ntf_base = 0.05, k_M2_ntf = 1.80, k_deg_ntf = 0.28,
  k_rgc_cyt = 0.0179, k_rgc_rpe = 0.005, k_rgc_iop = 0.002, EC50_ntf = 0.40,
  k_des_on = 5.000, k_des_off = 0.111, K_anti_M1 = 0.595,
  k_abs = 1.00, k_el_pep = 0.10,
  Emax_C1qblock = 6.0, Emax_C3aRblock = 5.0, Emax_switch = 7.0, Emax_migration = 5.0,
  EC50_pep = 0.80, gamma_pep = 2.00,
  dose_amount = 5.0
)

.hill <- function(C, Emax, EC50, g) { C <- max(C, 0); Emax * C^g / (EC50^g + C^g) }

glaucoma_odes <- function(t, state, p) {
  with(as.list(c(state, p)), {
    # defensive read clamps (match Python/Julia; inert under correct dynamics)
    RPE <- min(max(RPE, 0), 1); RGC <- min(max(RGC, 0), 1); R_des <- min(max(R_des, 0), 1)
    DAMPs <- max(DAMPs, 0); M0 <- max(M0, 0); M_mig <- max(M_mig, 0)
    M1 <- max(M1, 0); M2 <- max(M2, 0); C1q <- max(C1q, 0); C3 <- max(C3, 0)
    C3a <- max(C3a, 0); C5a <- max(C5a, 0); Cyt_pro <- max(Cyt_pro, 0)
    Cyt_anti <- max(Cyt_anti, 0); NTF <- max(NTF, 0); A_eye <- max(A_eye, 0); C_pep <- max(C_pep, 0)

    Stress <- max(0, (IOP_target - IOP_normal) / IOP_normal)
    Pb <- .hill(C_pep, Emax_C1qblock, EC50_pep, gamma_pep)
    Pc <- .hill(C_pep, Emax_C3aRblock, EC50_pep, gamma_pep)
    Ps <- .hill(C_pep, Emax_switch, EC50_pep, gamma_pep)
    Pm <- .hill(C_pep, Emax_migration, EC50_pep, gamma_pep)

    rpe_d <- k_rpe_stress * Stress + k_rpe_cyt * Cyt_pro + k_rpe_phago * M_mig
    prot  <- NTF / (EC50_ntf + NTF)
    rgc_d <- (k_rgc_cyt * Cyt_pro + k_rgc_rpe * (1 - RPE) + k_rgc_iop * Stress) * (1 - prot)
    mig   <- (k_mig_damp * DAMPs + k_mig_C5a * C5a) * (1 - Pm) * M0
    Mret  <- k_return * M_mig
    inhib_anti <- 1 / (1 + Cyt_anti / K_anti_M1)
    tlr4  <- k_damp_M1 * DAMPs * M_mig * inhib_anti
    c3ar  <- k_C3aR_act * C3a * (1 - R_des) * M_mig * (1 - Pc) * inhib_anti
    M1sw  <- k_M1_switch * (1 + Ps) * M1
    C3cl  <- k_C3_cleave * C1q * C3

    dRPE      <- -rpe_d * RPE
    dDAMPs    <- k_damp_rpe * rpe_d * RPE + k_damp_rgc * rgc_d * RGC - k_damp_clear * DAMPs
    dM0       <- -mig + k_deact_M1 * M1 + k_res_M2 * M2 + Mret
    dM_mig    <- mig - c3ar - tlr4 - Mret
    dM1       <- tlr4 + c3ar - M1sw - k_deact_M1 * M1
    dM2       <- M1sw - k_res_M2 * M2
    dC1q      <- k_C1q_M1 * M1 * (1 - Pb) - k_C1q_deg * C1q
    dC3       <- k_C3_base + k_C3_rpe * rpe_d * RPE - C3cl - k_C3_deg * C3
    dC3a      <- k_C3a_frac * C3cl - k_C3a_deg * C3a
    dC5a      <- k_C5a_frac * C3cl - k_C5a_deg * C5a
    dCyt_pro  <- k_M1_cyt * M1 - k_deg_pro * Cyt_pro - k_inhib * Cyt_anti * Cyt_pro
    dCyt_anti <- k_M2_cyt * M2 - k_deg_anti * Cyt_anti
    dNTF      <- k_ntf_base + k_M2_ntf * M2 - k_deg_ntf * NTF
    dRGC      <- -rgc_d * RGC
    dA_eye    <- -k_abs * A_eye
    dC_pep    <- k_abs * A_eye - k_el_pep * C_pep
    dR_des    <- k_des_on * C3a * (1 - R_des) - k_des_off * R_des

    list(c(dRPE, dDAMPs, dM0, dM_mig, dM1, dM2, dC1q, dC3, dC3a, dC5a,
           dCyt_pro, dCyt_anti, dNTF, dRGC, dA_eye, dC_pep, dR_des))
  })
}
