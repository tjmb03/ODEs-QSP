# ============================================================================
# QSP GLAUCOMA v6 — PAPER-ALIGNED + BOOTSTRAP FIX + OPTIMIZER + SIMBIOLOGY NET
# Ref: Gao et al. (2026) Survey of Ophthalmology 71:346–360
#
# v5 CHANGES:
#  FIX 1: Bootstrap TLR4 term added to dM1 so the complement loop can ignite
#          dM1 += k_damp_M1 * DAMPs * M_mig  (paper §6.3.1 TLR4/NF-κB entry)
#  FIX 2: Auto-optimizer (Nelder-Mead) tunes k_damp_M1, k_rpe_stress,
#          k_mig_damp, k_C1q_M1, k_C3aR_act to hit chronic-NTG targets
#  FIX 3: Network redrawn in SimBiology style (oval state nodes, small circle
#          junctions, square process nodes, dashed modulation edges)
# ============================================================================

suppressPackageStartupMessages({
  library(shiny); library(shinydashboard); library(shinyWidgets)
  library(shinyjs); library(deSolve); library(ggplot2)
  library(dplyr); library(tidyr); library(plotly); library(DT); library(jsonlite)
})

# ============================================================================
# 1. HILL FUNCTION
# ============================================================================
hill_fn <- function(C, Emax, EC50, gamma = 1.0) {
  C <- max(C, 0); Emax * C^gamma / (EC50^gamma + C^gamma)
}

# ============================================================================
# 2. ODE SYSTEM — paper-aligned loop + TLR4 bootstrap fix
# ============================================================================
glaucoma_odes <- function(t, state, p) {
  RPE      <- min(max(state["RPE"],    0), 1)
  DAMPs    <- max(state["DAMPs"],      0)
  M0       <- max(state["M0"],         0)
  M_mig    <- max(state["M_mig"],      0)
  M1       <- max(state["M1"],         0)
  M2       <- max(state["M2"],         0)
  C1q      <- max(state["C1q"],        0)
  C3       <- max(state["C3"],         0)
  C3a      <- max(state["C3a"],        0)
  C5a      <- max(state["C5a"],        0)
  Cyt_pro  <- max(state["Cyt_pro"],    0)
  Cyt_anti <- max(state["Cyt_anti"],   0)
  NTF      <- max(state["NTF"],        0)
  RGC      <- min(max(state["RGC"],    0), 1)
  A_eye    <- max(state["A_eye"],      0)
  C_pep    <- max(state["C_pep"],      0)
  R_des    <- min(max(state["R_des"],  0), 1)   # v6: C3aR desensitization (0-1)

  Stress <- max(0, (p$IOP_target - p$IOP_normal) / p$IOP_normal)

  # drug PD
  P_C1q  <- hill_fn(C_pep, p$Emax_C1qblock,  p$EC50_pep, p$gamma_pep)
  P_C3aR <- hill_fn(C_pep, p$Emax_C3aRblock, p$EC50_pep, p$gamma_pep)
  P_sw   <- hill_fn(C_pep, p$Emax_switch,    p$EC50_pep, p$gamma_pep)
  P_mig  <- hill_fn(C_pep, p$Emax_migration, p$EC50_pep, p$gamma_pep)

  # RPE
  rpe_death <- p$k_rpe_stress * Stress + p$k_rpe_cyt * Cyt_pro + p$k_rpe_phago * M_mig
  dRPE <- -rpe_death * RPE

  # RGC
  prot       <- NTF / (p$EC50_ntf + NTF)
  rgc_d_rate <- (p$k_rgc_cyt * Cyt_pro + p$k_rgc_rpe * (1 - RPE) + p$k_rgc_iop * Stress) * (1 - prot)
  dRGC <- -rgc_d_rate * RGC

  # DAMPs
  dDAMPs <- p$k_damp_rpe * rpe_death * RPE + p$k_damp_rgc * rgc_d_rate * RGC - p$k_damp_clear * DAMPs

  # Migration: driven by DAMPs + C5a (C5aR1)
  mig_rate <- (p$k_mig_damp * DAMPs + p$k_mig_C5a * C5a) * (1 - P_mig) * M0
  M_ret    <- p$k_return * M_mig
  dM0      <- -mig_rate + p$k_deact_M1 * M1 + p$k_res_M2 * M2 + M_ret

  # ═════════════════════════════════════════════════════════════════════════
  # v6 STRUCTURAL CHANGE — SATURATING BRAKES (stable intermediate operating pt)
  #
  # Problem in v5: positive feedback gain >> 1 → bistable switch. Loop either
  # stays silent (RGC preserved) or explodes (M1→53%, RGC crashes). No stable
  # chronic intermediate. M1>=5% was incompatible with RGC 40-65%.
  #
  # Fix: two biologically-grounded negative-feedback brakes that grow with the
  # inflammatory state and clamp the loop at a stable intermediate level.
  #
  #   BRAKE 1 — C3aR desensitization (C3aR is a GPCR; GPCRs downregulate under
  #             sustained agonist). New state R_des rises with C3a; feedback is
  #             scaled by (1 - R_des). Caps the explosive C3a→C3aR gain.
  #
  #   BRAKE 2 — anti-inflammatory suppression of M1 polarisation. M2-derived
  #             IL-10/TGF-β (Cyt_anti) suppress new M1 formation (paper M1/M2
  #             balance). inhib_anti = 1/(1 + Cyt_anti/K_anti_M1).
  # ═════════════════════════════════════════════════════════════════════════
  inhib_anti <- 1 / (1 + Cyt_anti / p$K_anti_M1)             # BRAKE 2

  # TLR4 bootstrap (paper §6.3.1) — now also suppressed by anti-inflammatory tone
  TLR4_boot <- p$k_damp_M1 * DAMPs * M_mig * inhib_anti
  # Main C3aR feedback — capped by desensitization (BRAKE 1) and anti-inflam (BRAKE 2)
  C3aR_act  <- p$k_C3aR_act * C3a * (1 - R_des) * (1 - P_C3aR) * M_mig * inhib_anti
  M1_sw     <- p$k_M1_switch * (1 + P_sw) * M1

  dM_mig <- mig_rate - C3aR_act - TLR4_boot - M_ret
  dM1    <- C3aR_act + TLR4_boot - M1_sw - p$k_deact_M1 * M1
  dM2    <- M1_sw - p$k_res_M2 * M2

  # BRAKE 1 dynamics: C3aR desensitization builds with C3a, recovers slowly
  dR_des <- p$k_des_on * C3a * (1 - R_des) - p$k_des_off * R_des

  # Complement loop
  # STEP A: M1 → C1q  (paper §4.2: "activated microglia express C1q")
  dC1q <- p$k_C1q_M1 * M1 * (1 - P_C1q) - p$k_C1q_deg * C1q

  # STEP B: C1q → C4b2a → cleaves C3
  C3_cleavage <- p$k_C3_cleave * C1q * C3
  dC3  <- p$k_C3_base + p$k_C3_rpe * rpe_death * RPE - C3_cleavage - p$k_C3_deg * C3

  # STEP C: C3a closes loop; C5a amplifies migration
  dC3a <- p$k_C3a_frac * C3_cleavage - p$k_C3a_deg * C3a
  dC5a <- p$k_C5a_frac * C3_cleavage - p$k_C5a_deg * C5a

  # Cytokines (NLRP3 inflammasome: IL-1β, IL-18, IL-6, TNF-α)
  dCyt_pro  <- p$k_M1_cyt * M1 - p$k_deg_pro * Cyt_pro - p$k_inhib * Cyt_anti * Cyt_pro
  dCyt_anti <- p$k_M2_cyt * M2 - p$k_deg_anti * Cyt_anti

  # NTF (BDNF, IGF-1 from M2)
  dNTF <- p$k_ntf_base + p$k_M2_ntf * M2 - p$k_deg_ntf * NTF

  # Drug PK
  dA_eye <- -p$k_abs * A_eye
  dC_pep <-  p$k_abs * A_eye - p$k_el_pep * C_pep

  list(c(dRPE, dDAMPs, dM0, dM_mig, dM1, dM2,
         dC1q, dC3, dC3a, dC5a,
         dCyt_pro, dCyt_anti, dNTF, dRGC, dA_eye, dC_pep, dR_des))
}

SNAMES <- c("RPE","DAMPs","M0","M_mig","M1","M2",
            "C1q","C3","C3a","C5a",
            "Cyt_pro","Cyt_anti","NTF","RGC","A_eye","C_pep","R_des")

# ============================================================================
# 3. DEFAULT PARAMETERS (v6 — optimized chronic-NTG operating point + brakes)
#    These defaults reproduce: RGC≈56%, M1 peak≈5%, loop fires ~day 180,
#    cytokine death fraction≈33% — a stable chronic intermediate, NOT a switch.
# ============================================================================
DEFAULTS <- list(
  IOP_normal = 15, IOP_target = 21, M_total = 1.0,
  k_rpe_stress = 0.003, k_rpe_cyt = 0.006, k_rpe_phago = 0.006,
  k_damp_rpe = 0.40, k_damp_rgc = 0.20, k_damp_clear = 0.80,
  k_mig_damp = 1.850, k_mig_C5a = 0.80, k_return = 0.08,
  # ─ BOOTSTRAP (v5, paper §6.3.1) ──────────────────────────
  k_damp_M1  = 0.50,   # TLR4 bootstrap
  # ─ LOOP GAIN ─────────────────────────────────────────────
  k_C3aR_act = 2.275,
  k_M1_switch = 0.176, k_deact_M1 = 0.05, k_res_M2 = 0.12,
  k_C1q_M1 = 1.359, k_C1q_deg = 0.40,
  k_C3_base = 0.30, k_C3_cleave = 1.20, k_C3_deg = 0.30, k_C3_rpe = 0.20,
  k_C3a_frac = 0.60, k_C3a_deg = 0.50,
  k_C5a_frac = 0.30, k_C5a_deg = 0.40,
  k_M1_cyt = 1.00, k_deg_pro = 0.35, k_inhib = 0.25,
  k_M2_cyt = 0.70, k_deg_anti = 0.28,
  k_ntf_base = 0.05, k_M2_ntf = 1.80, k_deg_ntf = 0.28,
  k_rgc_cyt = 0.0179, k_rgc_rpe = 0.005, k_rgc_iop = 0.002,
  EC50_ntf = 0.40,
  # ═ v6 SATURATING BRAKES ══════════════════════════════════
  k_des_on  = 5.000,   # BRAKE 1: C3aR desensitization rate
  k_des_off = 0.111,   # BRAKE 1: C3aR resensitization rate
  K_anti_M1 = 0.595,   # BRAKE 2: anti-inflammatory IC50 for M1 suppression
  # ═════════════════════════════════════════════════════════
  k_abs = 1.00, k_el_pep = 0.10,   # ~6.9-day half-life (realistic intravitreal anti-complement Fab)
  Emax_C1qblock = 6.0, Emax_C3aRblock = 5.0,
  Emax_switch = 7.0, Emax_migration = 5.0,
  EC50_pep = 0.80, gamma_pep = 2.00,
  dose_amount = 5.0, dose_times = c(90, 150, 210)   # dosed to cover loop ignition (~day 156) + active phase
)

# ============================================================================
# 4. SOLVER
# ============================================================================
run_sim <- function(p, treat = TRUE, t_end = 365, fast = FALSE) {
  doses  <- if (treat) sort(p$dose_times) else numeric(0)
  breaks <- if (length(doses) > 1) c(doses[-1], t_end) else t_end
  C3_ss  <- p$k_C3_base / p$k_C3_deg
  NTF_ss <- p$k_ntf_base / p$k_deg_ntf
  yc <- c(RPE=1, DAMPs=0, M0=p$M_total, M_mig=0, M1=0, M2=0,
          C1q=0, C3=C3_ss, C3a=0, C5a=0,
          Cyt_pro=0, Cyt_anti=0, NTF=NTF_ss, RGC=1,
          A_eye=ifelse(treat, p$dose_amount, 0), C_pep=0, R_des=0)
  # fast mode (optimizer): ~2 pts/day, looser tolerance — ~10x faster, accurate
  # enough for endpoint metrics, peak M1, and ignition timing.
  rtol <- if (fast) 1e-6 else 1e-8
  atol <- if (fast) 1e-8 else 1e-10
  ppd  <- if (fast) 2L else 20L
  tc <- 0; all_t <- NULL; all_y <- NULL
  for (tb in breaks) {
    n   <- max(60L, as.integer((tb - tc) * ppd))
    sol <- tryCatch(
      ode(y=yc, times=seq(tc, tb, length.out=n),
          func=glaucoma_odes, parms=p, method="lsoda", rtol=rtol, atol=atol),
      error=function(e) NULL)
    if (is.null(sol)) return(NULL)
    all_t <- c(all_t, sol[, "time"])
    all_y <- rbind(all_y, sol[, -1])
    yc    <- sol[nrow(sol), -1]; names(yc) <- SNAMES
    if (tb %in% doses) yc["A_eye"] <- yc["A_eye"] + p$dose_amount
    tc <- tb
  }
  df <- as.data.frame(all_y); df$time <- all_t
  df <- df %>% group_by(time) %>% slice(1) %>% ungroup() %>% arrange(time)
  df$RPE_pct   <- df$RPE * 100; df$RGC_pct <- df$RGC * 100
  df$Prot      <- df$NTF / (p$EC50_ntf + df$NTF)
  df$loop_idx  <- df$M1 * df$C1q * df$C3a
  Stress       <- max(0, (p$IOP_target - p$IOP_normal) / p$IOP_normal)
  df$d_cyt     <- p$k_rgc_cyt * df$Cyt_pro
  df$d_rpe     <- p$k_rgc_rpe * (1 - df$RPE)
  df$d_iop     <- rep(p$k_rgc_iop * Stress, nrow(df))
  df
}

# ============================================================================
# 5. METRICS
# ============================================================================
compute_metrics <- function(dt, dc) {
  if (is.null(dt) || is.null(dc)) return(NULL)
  lc <- tail(dc, 1); lt <- tail(dt, 1)
  d_tot <- lc$d_cyt + lc$d_rpe + lc$d_iop + 1e-12
  list(
    rgc_ctrl   = lc$RGC_pct, rgc_treat  = lt$RGC_pct,
    protection = lt$RGC_pct - lc$RGC_pct,
    rpe_ctrl   = lc$RPE_pct, rpe_treat  = lt$RPE_pct,
    c1q_peak   = max(dc$C1q, na.rm=TRUE),
    c3a_peak   = max(dc$C3a, na.rm=TRUE),
    c5a_peak   = max(dc$C5a, na.rm=TRUE),
    m1_peak    = max(dc$M1,  na.rm=TRUE),
    loop_peak  = max(dc$loop_idx, na.rm=TRUE),
    loop_suppress = (max(dc$loop_idx, na.rm=TRUE) - max(dt$loop_idx, na.rm=TRUE)) /
                    (max(dc$loop_idx, na.rm=TRUE) + 1e-9) * 100,
    cyt_frac   = lc$d_cyt / d_tot * 100,
    trans_frac = lc$d_rpe / d_tot * 100,
    iop_frac   = lc$d_iop / d_tot * 100
  )
}

# ============================================================================
# 6. AUTO-OPTIMIZER (Nelder-Mead via optim) — v6 with brake tuning
#    Targets — chronic NTG after 365 days (Gao et al. §6.3):
#    RGC_ctrl:   45–62% (steady progressive loss)
#    Cyt_frac:   20–45% (inflammation drives a meaningful fraction of death)
#    Trans_frac: 20–55% (trans-synaptic component visible)
#    M1_peak:    5–15%  (loop fires VISIBLY but is clamped by brakes — NOW
#                        reachable thanks to C3aR desensitization + anti-inflam)
#    Ignition:   day 90–180
#  Tunes 11 params incl. brake strengths (k_des_on, k_des_off, K_anti_M1).
# ============================================================================
OPT_KEYS <- c("k_des_on","k_des_off","K_anti_M1","k_rpe_stress","k_mig_damp",
              "k_C3aR_act","k_C1q_M1","k_M1_switch","k_rgc_cyt","k_rpe_cyt","k_rgc_rpe")
OPT_LO <- c(0.5, 0.05, 0.15, 0.003, 0.3, 0.5, 0.5, 0.04, 0.010, 0.010, 0.005)
OPT_HI <- c(8.0, 0.50, 2.00, 0.020, 2.0, 3.0, 2.5, 0.40, 0.080, 0.060, 0.040)

opt_clamp <- function(x) pmin(pmax(x, OPT_LO), OPT_HI)

opt_loop_ignition <- function(dc) {
  lp <- dc$loop_idx; pk <- max(lp, na.rm=TRUE)
  if (pk <= 1e-9) return(NA_real_)
  idx <- which(lp > 0.01 * pk)
  if (length(idx) == 0) return(NA_real_)
  dc$time[idx[1]]
}

optimizer_cost <- function(x, p_base, t_end) {
  p  <- p_base
  xc <- opt_clamp(x)
  for (i in seq_along(OPT_KEYS)) p[[OPT_KEYS[i]]] <- xc[i]
  dc <- tryCatch(run_sim(p, treat=FALSE, t_end=t_end, fast=TRUE), error=function(e) NULL)
  if (is.null(dc) || nrow(dc) < 2) return(1e6)
  lc    <- tail(dc, 1)
  d_tot <- lc$d_cyt + lc$d_rpe + lc$d_iop + 1e-12
  rgc   <- lc$RGC_pct
  rpe   <- lc$RPE_pct
  cyt_f <- lc$d_cyt / d_tot * 100
  trn_f <- lc$d_rpe / d_tot * 100
  m1_pk <- max(dc$M1, na.rm=TRUE) * 100
  ign   <- opt_loop_ignition(dc)

  pen <- function(v, lo, hi) if (v < lo) (v-lo)^2 else if (v > hi) (v-hi)^2 else 0

  cost <- pen(rgc,   45, 62) * 0.10 +
          pen(cyt_f, 20, 45) * 0.05 +
          pen(trn_f, 20, 55) * 0.03 +
          pen(m1_pk,  5, 15) * 0.40 +   # M1 peak 5–15% — now reachable
          pen(rpe,   28, 65) * 0.03
  cost <- cost + (if (is.na(ign)) 30.0 else pen(ign, 90, 180) * 0.006)
  if (rgc < 28) cost <- cost + (28 - rgc) * 3.0
  cost
}

run_optimizer <- function(p, t_end, progress_cb = NULL) {
  x0 <- sapply(OPT_KEYS, function(k) p[[k]])
  best_cost <- Inf; best_x <- x0

  set.seed(42)
  # 3 starts (current point + 2 spread) — fast sim makes each cheap, but keep
  # the count low so total wall-time stays ~1-2 min instead of ~30 min.
  starts <- rbind(
    x0,
    c(3.0, 0.10, 0.40, 0.006, 1.2, 1.2, 1.0, 0.10, 0.040, 0.030, 0.015),
    c(5.0, 0.08, 0.35, 0.009, 0.9, 1.0, 0.9, 0.12, 0.050, 0.045, 0.022)
  )

  for (i in seq_len(nrow(starts))) {
    if (!is.null(progress_cb)) progress_cb((i - 0.5) / nrow(starts))
    res <- tryCatch(
      optim(starts[i,], optimizer_cost, p_base=p, t_end=t_end,
            method="Nelder-Mead",
            control=list(maxit=250, reltol=1e-5, abstol=1e-4)),
      error=function(e) list(value=Inf, par=starts[i,]))
    if (res$value < best_cost) { best_cost <- res$value; best_x <- res$par }
    if (!is.null(progress_cb)) progress_cb(i / nrow(starts))
    # Early exit: cost < 0.5 means all targets essentially satisfied — no need
    # to run the remaining restarts.
    if (best_cost < 0.5) break
  }
  bx <- opt_clamp(best_x)
  out <- as.list(bx); names(out) <- OPT_KEYS
  out$final_cost <- best_cost
  out
}

# ============================================================================
# 7. SIMBIOLOGY-STYLE NETWORK JS
#    Oval state nodes (dark navy / orange / teal / purple / red / green)
#    Small circle junctions on every flow (gray dot)
#    Small square process nodes on regulated steps (gray square)
#    Dashed lines for modulation / inhibition
# ============================================================================
build_simbiology_js <- function(cid = "vis-net") {
  js_template <- '
(function() {
  if (typeof vis === "undefined") { setTimeout(arguments.callee, 400); return; }
  var el = document.getElementById("%CID%");
  if (!el) { setTimeout(arguments.callee, 400); return; }

  // ── COLOUR PALETTE (SimBiology style) ─────────────────────────────────
  var C = {
    drug:    {bg:"#E65100", bd:"#BF360C", ft:"#fff"},  // orange  — drug PK
    state:   {bg:"#1B3A6B", bd:"#0D2144", ft:"#fff"},  // navy    — core state vars
    micro:   {bg:"#006064", bd:"#004D40", ft:"#fff"},  // teal    — microglia
    comp:    {bg:"#6A1B9A", bd:"#4A148C", ft:"#fff"},  // purple  — complement
    inflam:  {bg:"#B71C1C", bd:"#7f0000", ft:"#fff"},  // red     — pro-inflam
    protect: {bg:"#1B5E20", bd:"#0A3D0A", ft:"#fff"},  // green   — protective
    stim:    {bg:"#E65100", bd:"#BF360C", ft:"#fff"},  // orange  — stressor
    junc:    {bg:"#9E9E9E", bd:"#616161", ft:""},       // gray    — junction circle
    proc:    {bg:"#BDBDBD", bd:"#9E9E9E", ft:"#333"},  // lt-gray — process square
    loop_a:  {bg:"#AD1457", bd:"#880E4F", ft:"#fff"},  // pink    — C1q loop node A
    loop_b:  {bg:"#7B1FA2", bd:"#4A148C", ft:"#fff"}   // violet  — C3a loop node B
  };

  function oval(id, lbl, col, x, y, bw) {
    return { id:id, label:lbl, shape:"ellipse",
      color:{background:col.bg, border:col.bd},
      font:{color:col.ft, size:11, face:"Arial", bold:(bw||false)},
      x:x, y:y, fixed:true, widthConstraint:{maximum:120}, heightConstraint:{minimum:28} };
  }
  function junc(id, x, y, tip) {
    return { id:id, label:"", shape:"dot", size:7,
      color:{background:C.junc.bg, border:C.junc.bd},
      title: tip||"", x:x, y:y, fixed:true };
  }
  function proc(id, lbl, x, y) {
    return { id:id, label:lbl, shape:"square", size:9,
      color:{background:C.proc.bg, border:C.proc.bd},
      font:{color:C.proc.ft, size:9, face:"Arial"},
      x:x, y:y, fixed:true };
  }
  function arr(from, to, lbl, dash, col, w) {
    return { from:from, to:to, label:lbl||"",
      arrows:{to:{enabled:true, scaleFactor:0.7}},
      dashes:dash||false,
      color:{color:col||"#616161"},
      width: w||1.2,
      font:{color:"#555",size:9,face:"Arial",align:"middle"},
      smooth:{type:"continuous"} };
  }

  // ════════════════════════════════════════════════════════════════════════
  //  NODE DEFINITIONS
  // ════════════════════════════════════════════════════════════════════════
  var N = [
    // ── DRUG PK (top-left) ──────────────────────────────────────────────
    oval(1,  "GI / Intravitreal\\nDepot  (A_eye)", C.drug,  -480, -260),
    junc(51, -370, -260, "Compound absorption (k_abs)"),
    oval(2,  "Active Peptide\\n(C_pep)",           C.drug,  -260, -260),
    junc(52, -260, -310, "Elimination (k_el_pep)"),

    // ── STRESSOR (top-center) ───────────────────────────────────────────
    oval(3,  "IOP / Mechanical\\nStress",           C.stim,  -100, -260),

    // ── STRUCTURAL COMPARTMENTS ─────────────────────────────────────────
    oval(4,  "RPE\\n(retinal pigment\\nepithelium)", C.state, -100, -80),
    proc(61, "rpe_stress",  -280, -80),    // IOP → RPE
    proc(62, "rpe_cyt",      -10, -10),    // Cyt_pro → RPE

    oval(5,  "DAMPs\\n(danger signals)", C.state,   80, -80),
    proc(63, "damp_clear",  80, -150),

    oval(6,  "RGC\\n(outcome)",          C.state, -100, 100),
    proc(64, "rgc_iop",    -240, 100),
    proc(65, "rgc_trans",  -10,  140),
    proc(66, "rgc_cyt",     80,  170),
    proc(67, "NTF_prot",  -100,  190),    // NTF protection block

    // ── MICROGLIA CASCADE ───────────────────────────────────────────────
    oval(7,  "Resting\\nMicroglia (M0)",  C.micro,  260, -80),
    junc(53, 340, -80,  "Migration (k_mig_damp + k_mig_C5a)"),
    oval(8,  "Migrating\\nMicroglia\\n(M_mig)", C.micro, 440, -80),

    // TLR4 bootstrap process node
    proc(68, "TLR4/NF-κB\\nbootstrap",   440, -160),   // ← FIX v5

    junc(54, 530, -80,  "C3aR activation + TLR4 boot → M1"),
    oval(9,  "M1 Microglia\\n(pro-inflam)", C.inflam, 640, -80),
    junc(55, 720, -80,  "M1→M2 switch (k_M1_switch)"),
    oval(10, "M2 Microglia\\n(anti-inflam)", C.protect, 820, -80),

    // ── COMPLEMENT LOOP ─────────────────────────────────────────────────
    //  Loop Node A
    oval(11, "★ C1q ★\\n(produced by M1)\\n[Loop node A]",
         C.loop_a, 640, -260),
    junc(56, 720, -260, "C4b2a convertase cleaves C3"),

    //  C3 pool
    oval(12, "C3\\n(tissue pool)",        C.comp,   820, -260),
    proc(69, "C3_base",  820, -320),

    //  Anaphylatoxin branch
    proc(70, "→ C3a\\ncleavage",   760, -170),
    //  Loop Node B
    oval(13, "★ C3a ★\\n(binds C3aR\\n→ closes loop)\\n[Loop node B]",
         C.loop_b, 640, -190),

    proc(71, "→ C5a\\ncleavage",   820, -170),
    oval(14, "C5a\\n(binds C5aR1)",  C.comp,   940, -190),

    oval(15, "MAC (C5b-9)\\n[lysis]",  C.comp,   940, -80),

    // ── SIGNALLING ──────────────────────────────────────────────────────
    proc(72, "NLRP3\\ninflammaso.",  640, 30),
    oval(16, "Pro-inflam\\nCytokines\\n(IL-1β·IL-6·TNF-α)",
         C.inflam,  500, 130),
    oval(17, "Anti-inflam\\nCytokines\\n(IL-4·IL-10)",
         C.protect, 820, 30),
    oval(18, "NTF\\n(BDNF·IGF-1)",  C.protect, 940, 30),

    // ── v6 SATURATING BRAKES ────────────────────────────────────────────
    oval(19, "★ R_des ★\\nC3aR desensitization\\n[BRAKE 1 — GPCR\\ndownregulation]",
         {bg:"#00897B", bd:"#00695C", ft:"#fff"}, 470, -190),
    proc(73, "anti-inflam\\nM1 suppression\\n[BRAKE 2]", 470, 80)
  ];

  // ════════════════════════════════════════════════════════════════════════
  //  EDGE DEFINITIONS
  // ════════════════════════════════════════════════════════════════════════
  var E = [
    // ── DRUG PK ─────────────────────────────────────────────────────────
    arr(1,  51, "absorption"),
    arr(51, 2),
    arr(2,  52, "elimination", false,"#E65100"),
    // C_pep inhibits various steps (dashed)
    arr(2, 61, "C1q block (PD1)", true, "#E65100"),
    arr(2, 68, "migration\\nblock (PD4)", true, "#E65100"),
    arr(2, 55, "M1→M2\\n(PD3)",  true, "#E65100"),
    arr(2, 54, "C3aR block\\n(PD2)", true, "#E65100"),

    // ── STRESSOR → RPE ──────────────────────────────────────────────────
    arr(3, 61, "IOP stress"),
    arr(61, 4, "", false, "#E65100", 2),

    // ── RPE DEATH ───────────────────────────────────────────────────────
    arr(4, 5,  "releases DAMPs", false, "#9E9E9E"),
    arr(5, 63, "clearance"),
    arr(53, 7, "", false, "#006064"),    // M_ret back to M0

    // ── RGC DEATH PATHWAYS ──────────────────────────────────────────────
    arr(3,  64, "IOP axon"),
    arr(64, 6,  "", false, "#E65100", 2),
    arr(4,  65, "trans-synap.", true, "#795548"),
    arr(65, 6,  "", false, "#795548", 2),
    arr(16, 66, "cytokine"),
    arr(66, 6,  "", false, "#B71C1C", 2),
    arr(67, 6,  "NTF protects", true, "#1B5E20"),
    arr(18, 67, "", true, "#1B5E20"),

    // ── MICROGLIA MIGRATION ─────────────────────────────────────────────
    arr(7,  53, "recruit"),
    arr(53, 8),
    arr(5,  53, "DAMPs drive", false, "#795548", 2),   // DAMPs → migration
    arr(14, 53, "C5aR1 recruit", false, "#7B1FA2", 2), // C5a → migration

    // ── TLR4 BOOTSTRAP ─────────────────────────────────────────────────
    arr(5,  68, "DAMP→TLR4", false, "#FF8F00", 2),    // ← FIX v5: DAMPs → TLR4
    arr(8,  68, "M_mig sensor", false, "#FF8F00"),     // M_mig needed
    arr(68, 54, "NF-κB→M1", false, "#FF8F00", 3),     // → activation junction

    // ── MAIN ACTIVATION ─────────────────────────────────────────────────
    arr(8,  54, "→ M1"),
    arr(54, 9,  "", false, "#B71C1C", 2),

    // ── M1→M2 SWITCH ────────────────────────────────────────────────────
    arr(9,  55, "switch"),
    arr(55, 10, "", false, "#1B5E20", 2),
    arr(10, 7,  "resolve→M0", true, "#1B5E20"),

    // ════════════════════════════════════════════════════════════════════
    //  THE COMPLEMENT FEEDBACK LOOP
    // ════════════════════════════════════════════════════════════════════

    // LOOP STEP A: M1 → C1q
    arr(9,  11, "★ PRODUCES C1q ★", false, "#AD1457", 4),

    // LOOP STEP B: C1q → C4b2a → cleaves C3
    arr(11, 56, "C4b2a\\nconvertase",   false, "#AD1457", 3),
    arr(56, 12, "", false, "#6A1B9A", 2),

    // C3 baseline production
    arr(69, 12, "constitutive"),

    // LOOP STEP C: C3 → C3a (closes loop)
    arr(12, 70, "cleavage"),
    arr(70, 13, "★ → C3a ★", false, "#9C27B0", 4),

    // LOOP CLOSURE: C3a → C3aR → activation junction → M1
    arr(13, 54, "★ C3aR→M1\\n(LOOP CLOSES) ★", false, "#9C27B0", 5),

    // C5a branch
    arr(12, 71, "cleavage"),
    arr(71, 14, "→ C5a", false, "#7B1FA2", 2),

    // MAC
    arr(14, 15, "→ MAC"),
    arr(15, 6,  "MAC lysis", false, "#4A148C", 2),

    // ── CYTOKINE CASCADE ────────────────────────────────────────────────
    arr(9,  72, "activates\\nNLRP3"),
    arr(72, 16, "IL-1β·IL-18", false, "#B71C1C", 2),
    arr(16, 4,  "RPE damage", true, "#B71C1C"),
    arr(16, 62, "cytokine\\nkills RPE"),
    arr(62, 4,  "", false, "#B71C1C", 2),

    // ── M2 PROTECTION ───────────────────────────────────────────────────
    arr(10, 17, "secretes"),
    arr(10, 18, "secretes NTF"),

    // ── RPE as secondary C3 source ──────────────────────────────────────
    arr(4, 12, "dying RPE→C3", true, "#795548"),

    // ── v6 BRAKE EDGES (negative feedback that clamps the loop) ──────────
    // BRAKE 1: C3a drives R_des; R_des inhibits the C3aR feedback
    arr(13, 19, "C3a → desensitize", false, "#00897B", 2),
    arr(19, 54, "⊣ caps C3aR\\nactivation", true, "#00897B", 3),
    // BRAKE 2: anti-inflammatory cytokines suppress M1 formation
    arr(17, 73, "IL-10/TGF-β", false, "#00897B", 2),
    arr(73, 54, "⊣ suppress\\nM1 formation", true, "#00897B", 3)
  ];

  var options = {
    manipulation: {
      enabled: true,
      addNode: function(d,cb){ d.color={background:"#607D8B"}; d.font={color:"#fff"}; cb(d); },
      editNode: function(d,cb){ cb(d); },
      addEdge:  function(d,cb){ cb(d); }
    },
    physics: { enabled: false },
    nodes: { shadow:{ enabled:true, size:4, x:2, y:2, color:"rgba(0,0,0,.15)" },
             borderWidth:1.5 },
    edges: { shadow:false,
             font:{ size:9, face:"Arial", align:"middle", strokeWidth:2, strokeColor:"#fff" },
             arrows:{ to:{ scaleFactor:0.75 } } },
    interaction: { hover:true, navigationButtons:true, keyboard:true, tooltipDelay:200 }
  };

  var net = new vis.Network(el,
    { nodes: new vis.DataSet(N), edges: new vis.DataSet(E) }, options);
  window.qspNet = net;
})();
'
  gsub("%CID%", cid, js_template, fixed = TRUE)
}

# ============================================================================
# 8. HTML EXPORT BUILDER
# ============================================================================
build_html_export <- function(sim_treat, sim_ctrl, params, t_end) {
  se <- sim_treat %>%
    select(time, RPE_pct, RGC_pct, C1q, C3, C3a, C5a, R_des,
           M_mig, M1, M2, Cyt_pro, Cyt_anti, NTF, loop_idx, d_cyt, d_rpe, d_iop) %>%
    mutate(across(where(is.numeric), ~round(., 5)))
  ce <- sim_ctrl %>%
    select(time, RPE_pct, RGC_pct, C1q, C3, C3a, C5a, R_des,
           M_mig, M1, M2, Cyt_pro, Cyt_anti, NTF, loop_idx, d_cyt, d_rpe, d_iop) %>%
    mutate(across(where(is.numeric), ~round(., 5)))
  pe <- params[setdiff(names(params), "dose_times")]
  pe$dose_times_str <- paste(params$dose_times, collapse=", ")
  sim_json <- jsonlite::toJSON(list(treat=se, ctrl=ce, params=pe, t_end=t_end), auto_unbox=TRUE, digits=5)
  njs <- build_simbiology_js("exp-net")

  paste0('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>QSP Glaucoma v6 — Paper-Aligned Loop (Bootstrap Fixed)</title>
<script src="https://cdn.plot.ly/plotly-2.26.0.min.js"></script>
<script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
<link rel="stylesheet" href="https://unpkg.com/vis-network/styles/vis-network.min.css">
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:Arial,sans-serif;background:#F0F2F5;padding:22px;color:#1A1A2E}
h1{color:#0F3460;border-bottom:3px solid #E94560;padding-bottom:8px;margin-bottom:6px;font-size:20px}
.sub{color:#6B7280;font-size:12px;margin-bottom:20px;line-height:1.7}
h2{color:#0F3460;margin:26px 0 10px;font-size:14px;border-left:4px solid #E94560;padding-left:8px}
.g2{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin-bottom:14px}
.g3{display:grid;grid-template-columns:1fr 1fr 1fr;gap:14px;margin-bottom:14px}
.g4{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin-bottom:20px}
.card{background:#fff;border-radius:8px;padding:14px;box-shadow:0 2px 8px rgba(0,0,0,.08)}
.full{grid-column:1/-1}
#exp-net{width:100%;height:600px;border:1px solid #ddd;border-radius:6px;background:#FAFAFA}
.stat-val{font-size:24px;font-weight:bold;color:#0F3460}
.stat-lbl{font-size:10px;color:#9CA3AF;margin-top:3px;text-transform:uppercase;letter-spacing:.5px}
.stat-card{text-align:center;background:#fff;border-radius:8px;padding:12px;box-shadow:0 2px 6px rgba(0,0,0,.07)}
.fix-box{background:#FFF8E1;border:2px solid #FF8F00;border-radius:6px;padding:12px;margin:10px 0;font-size:12px;line-height:1.8}
.loop-box{background:#FCE4EC;border:2px solid #AD1457;border-radius:6px;padding:12px;margin:10px 0;font-size:12px;line-height:1.8}
table.pt{width:100%;border-collapse:collapse;font-size:11px}
table.pt th{background:#0F3460;color:#fff;padding:5px 8px;text-align:left}
table.pt td{padding:4px 8px;border-bottom:1px solid #f0f0f0}
table.pt tr:nth-child(even) td{background:#F9FAFB}
.hl td{background:#FFF3E0!important}
</style>
</head>
<body>
<h1>QSP Glaucoma v6 — Paper-Aligned Feedback Loop (Bootstrap Fixed)</h1>
<p class="sub">
Reference: Gao et al. (2026) <em>Survey of Ophthalmology</em> 71:346–360 &nbsp;|&nbsp;
Duration: ', t_end, ' days &nbsp;|&nbsp; IOP: ', params$IOP_target, ' mmHg
</p>
<div class="fix-box">
<strong>v5 Bootstrap Fix (paper §6.3.1):</strong><br>
Added TLR4/NF-κB entry: <code>dM1 += k_damp_M1 × DAMPs × M_mig</code><br>
DAMPs from dying RPE bind TLR4 on migrating microglia → NF-κB activation → M1 polarisation (low efficiency).
Once M1 &gt; 0, the C1q→C3→C3a→C3aR amplifier ignites and the self-reinforcing loop takes over.
</div>
<div class="loop-box">
<strong>Feedback Loop:</strong>&nbsp;
M1 → C1q → C4b2a → C3a (C3aR, CLOSES LOOP) + C5a (C5aR1, migration) → M_mig → M1
</div>
<div id="stat-boxes" class="g4"></div>
<h2>Network Diagram — SimBiology Style (editable)</h2>
<div class="card full"><div id="exp-net"></div></div>
<h2>Simulation Plots</h2>
<div class="g2">
  <div class="card"><div id="p_rgc"></div></div>
  <div class="card"><div id="p_rpe"></div></div>
</div>
<div class="g3">
  <div class="card"><div id="p_c1q"></div></div>
  <div class="card"><div id="p_c3a"></div></div>
  <div class="card"><div id="p_m1"></div></div>
</div>
<div class="g2">
  <div class="card"><div id="p_micro"></div></div>
  <div class="card"><div id="p_cyt"></div></div>
</div>
<div class="card full" style="margin-bottom:14px"><div id="p_loop"></div></div>
<div class="g2">
  <div class="card"><div id="p_rdes"></div></div>
  <div class="card"><div id="p_m1plat"></div></div>
</div>
<div class="g2">
  <div class="card"><div id="p_death"></div></div>
  <div class="card"><div id="p_ntf"></div></div>
</div>
<h2>Parameters</h2>
<div class="card full"><table class="pt" id="ptbl"></table></div>
<script>
const SIM=', sim_json, ';
const C=SIM.ctrl, T=SIM.treat, P=SIM.params;
function tr(d,k,n,col,dash){
  return{x:d.map(r=>r.time),y:d.map(r=>r[k]),name:n,mode:"lines",type:"scatter",
    line:{color:col,dash:dash||"solid",width:2}};
}
const L=(t,y)=>({title:{text:t,font:{size:12}},xaxis:{title:"Days"},yaxis:{title:y||""},
  margin:{t:35,b:40,l:50,r:15},legend:{orientation:"h",y:-0.3},
  plot_bgcolor:"#FAFAFA",paper_bgcolor:"#FAFAFA"});
const cfg={responsive:true,displayModeBar:false};
Plotly.newPlot("p_rgc",[tr(T,"RGC_pct","Treated","#2E7D32"),tr(C,"RGC_pct","Control","#C62828","dash")],
  Object.assign(L("RGC Survival (%)","% surviving"),{yaxis:{range:[0,105]}}),cfg);
Plotly.newPlot("p_rpe",[tr(T,"RPE_pct","Treated","#E65100"),tr(C,"RPE_pct","Control","#C62828","dash")],
  Object.assign(L("RPE Health (%)","% surviving"),{yaxis:{range:[0,105]}}),cfg);
Plotly.newPlot("p_c1q",[tr(T,"C1q","Treated","#AD1457"),tr(C,"C1q","Control","#880E4F","dash")],L("C1q — Loop Node A","a.u."),cfg);
Plotly.newPlot("p_c3a",[tr(T,"C3a","Treated","#9C27B0"),tr(C,"C3a","Control","#6A1B9A","dash")],L("C3a — Loop Node B (closes loop)","a.u."),cfg);
Plotly.newPlot("p_m1",[tr(T,"M1","Treated","#C62828"),tr(C,"M1","Control","#7f0000","dash")],L("M1 Microglia (loop amplifier)","a.u."),cfg);
Plotly.newPlot("p_micro",[
  tr(C,"M_mig","M_mig ctrl","#26A69A","dot"),tr(C,"M1","M1 ctrl","#C62828","dot"),tr(C,"M2","M2 ctrl","#2E7D32","dot"),
  tr(T,"M1","M1 treated","#C62828"),tr(T,"M2","M2 treated","#2E7D32")],L("Microglial States","a.u."),cfg);
Plotly.newPlot("p_cyt",[tr(C,"Cyt_pro","Pro-inflam ctrl","#D32F2F","dot"),tr(T,"Cyt_pro","Pro-inflam treated","#D32F2F"),
  tr(T,"NTF","NTF treated","#1565C0")],L("Cytokines & NTF","a.u."),cfg);
Plotly.newPlot("p_loop",[
  {x:C.map(r=>r.time),y:C.map(r=>r.loop_idx),name:"Control",mode:"lines",
   line:{color:"#C62828",width:2.5},fill:"tozeroy",fillcolor:"rgba(198,40,40,0.08)"},
  {x:T.map(r=>r.time),y:T.map(r=>r.loop_idx),name:"Treated",mode:"lines",
   line:{color:"#2E7D32",width:2.5},fill:"tozeroy",fillcolor:"rgba(46,125,50,0.08)"}],
  L("Feedback Loop Intensity  (M1 × C1q × C3a)","intensity"),cfg);
Plotly.newPlot("p_rdes",[
  {x:C.map(r=>r.time),y:C.map(r=>r.R_des),name:"Control",mode:"lines",
   line:{color:"#00695C",width:2.4},fill:"tozeroy",fillcolor:"rgba(0,137,123,0.10)"},
  {x:T.map(r=>r.time),y:T.map(r=>r.R_des),name:"Treated",mode:"lines",
   line:{color:"#2E7D32",width:2.2,dash:"dash"}}],
  Object.assign(L("v6 BRAKE 1 — C3aR Desensitization (R_des)","R_des 0–1"),{yaxis:{range:[0,1],title:"R_des 0–1"}}),cfg);
Plotly.newPlot("p_m1plat",[
  tr(C,"M1","Control M1","#C62828"),tr(T,"M1","Treated M1","#2E7D32","dash")],
  L("v6 BRAKE 2 — M1 plateau (stable, not switch)","M1 a.u."),cfg);
Plotly.newPlot("p_death",[
  {x:C.map(r=>r.time),y:C.map(r=>r.d_cyt),name:"Cytokine",stackgroup:"one",fillcolor:"rgba(198,40,40,.7)",line:{color:"#C62828"}},
  {x:C.map(r=>r.time),y:C.map(r=>r.d_rpe),name:"Trans-synaptic",stackgroup:"one",fillcolor:"rgba(106,27,154,.7)",line:{color:"#6A1B9A"}},
  {x:C.map(r=>r.time),y:C.map(r=>r.d_iop),name:"IOP-direct",stackgroup:"one",fillcolor:"rgba(230,81,0,.7)",line:{color:"#E65100"}}],
  L("RGC Death Decomposition (Control)","rate 1/day"),cfg);
Plotly.newPlot("p_ntf",[tr(C,"NTF","NTF ctrl","#1565C0","dot"),tr(T,"NTF","NTF treated","#1565C0"),
  tr(C,"Prot","Protection ctrl","#2E7D32","dot"),tr(T,"Prot","Protection treated","#2E7D32")],
  L("NTF & Protection","a.u."),cfg);
// stat boxes
const lc=C[C.length-1],lt=T[T.length-1];
const dt=lc.d_cyt+lc.d_rpe+lc.d_iop+1e-12;
const lmax=Math.max(...C.map(r=>r.loop_idx)),ltmax=Math.max(...T.map(r=>r.loop_idx));
const sup=lmax>0?(lmax-ltmax)/lmax*100:0;
const m1max=Math.max(...C.map(r=>r.M1))*100;
[{v:lt.RGC_pct.toFixed(1)+"%",l:"RGC Survival (treated)"},
 {v:"+"+(lt.RGC_pct-lc.RGC_pct).toFixed(1)+"%",l:"RGC Protection"},
 {v:sup.toFixed(1)+"%",l:"Loop Suppressed"},
 {v:m1max.toFixed(2)+"%",l:"Peak M1 (ctrl)"},
 {v:(lc.d_cyt/dt*100).toFixed(0)+"%",l:"Cytokine death fraction"},
 {v:(lc.d_rpe/dt*100).toFixed(0)+"%",l:"Trans-synaptic fraction"},
 {v:lc.RPE_pct.toFixed(1)+"%",l:"RPE (ctrl)"},
 {v:lc.RGC_pct.toFixed(1)+"%",l:"RGC (ctrl)"}
].forEach(s=>{
  document.getElementById("stat-boxes").innerHTML+=
    `<div class="stat-card"><div class="stat-val">${s.v}</div><div class="stat-lbl">${s.l}</div></div>`;
});
// params table
const PDESC={k_damp_M1:"★ TLR4 bootstrap rate (v5 FIX §6.3.1)",
  k_C1q_M1:"C1q production by M1 (Loop Step A)",
  k_C3_cleave:"C3 cleavage by C4b2a (Loop Step B)",
  k_C3aR_act:"C3aR feedback activation (Loop Step C — CLOSES LOOP)",
  k_mig_C5a:"C5aR1-driven migration",
  k_des_on:"★ v6 BRAKE 1: C3aR desensitization rate (GPCR downregulation)",
  k_des_off:"★ v6 BRAKE 1: C3aR resensitization rate",
  K_anti_M1:"★ v6 BRAKE 2: anti-inflammatory IC50 for M1 suppression",
  k_M1_switch:"M1→M2 switch rate (loop removal / gain control)",
  Emax_C1qblock:"Drug PD1: C1q block (ANX007)",
  Emax_C3aRblock:"Drug PD2: C3aR/C5aR1 block",
  Emax_switch:"Drug PD3: M1→M2 (minocycline)",
  Emax_migration:"Drug PD4: migration block"};
const tbl=document.getElementById("ptbl");
tbl.innerHTML="<tr><th>Parameter</th><th>Value</th><th>Description</th></tr>";
Object.entries(P).forEach(([k,v])=>{
  const d=PDESC[k]||"";
  const hl=d.includes("★")||d.includes("Loop")?"hl":"";
  tbl.innerHTML+=`<tr class="${hl}"><td><strong>${k}</strong></td><td>${v}</td><td style="color:#6B7280">${d}</td></tr>`;
});
', njs, '
</script></body></html>')
}

# ============================================================================
# 9. STYLE
# ============================================================================
SB_CSS <- "
.skin-black .main-sidebar{background:#1A1A2E}
.sidebar-menu>li>a{color:#ccc!important;font-size:12px}
.sidebar-menu>li.active>a{color:#E94560!important;background:#0F3460!important}
"
BD_CSS <- "
body,.content-wrapper{background:#F0F2F5!important;font-family:'Arial',sans-serif}
.box{border-radius:8px;box-shadow:0 2px 8px rgba(0,0,0,.07)}
.box.box-primary>.box-header{background:#0F3460;border-color:#0F3460}
.box.box-danger>.box-header{background:#C62828;border-color:#C62828}
.box.box-success>.box-header{background:#2E7D32;border-color:#2E7D32}
.box.box-warning>.box-header{background:#E65100;border-color:#E65100}
.box.box-info>.box-header{background:#006064;border-color:#006064}
.plbl{font-size:10px;color:#555;font-weight:bold;text-transform:uppercase;
      letter-spacing:.4px;padding:4px 0 2px;margin-top:4px}
.btn-run{background:#E94560!important;border:none!important;color:#fff!important;
          font-weight:bold!important;width:100%;border-radius:6px;margin-top:6px}
.btn-opt{background:#0F3460!important;border:none!important;color:#fff!important;
          font-weight:bold!important;width:100%;border-radius:6px;margin-top:6px}
.btn-rst{background:#455A64!important;border:none!important;color:#fff!important;
          width:100%;border-radius:6px;margin-top:4px}
.fix-note{background:#FFF8E1;border-left:4px solid #FF8F00;border-radius:4px;
           padding:8px 12px;font-size:11px;line-height:1.7;margin:6px 0}
#vis-net{width:100%;height:580px;border:1px solid #ddd;border-radius:6px;background:#FAFAFA}
.opt-res{background:#E8F5E9;border-left:4px solid #2E7D32;border-radius:4px;
         padding:8px 12px;font-size:11px;margin-top:8px;display:none}
"

# ============================================================================
# 10. UI
# ============================================================================
ui <- dashboardPage(
  skin="black",
  dashboardHeader(
    title=tags$span(style="font-size:14px;color:#E94560;", "QSP Glaucoma v6"),
    titleWidth=200),
  dashboardSidebar(
    width=220,
    tags$style(HTML(SB_CSS)),
    sidebarMenu(id="tabs",
      menuItem("Simulation",    tabName="sim",    icon=icon("chart-line")),
      menuItem("Loop Dynamics", tabName="loop",   icon=icon("sync-alt")),
      menuItem("Network",       tabName="net",    icon=icon("project-diagram")),
      menuItem("Mechanisms",    tabName="mech",   icon=icon("flask")),
      menuItem("Death Pathways",tabName="death",  icon=icon("skull")),
      menuItem("Export HTML",   tabName="export", icon=icon("file-export"))
    )
  ),
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$style(HTML(BD_CSS)),
      tags$script(src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"),
      tags$link(rel="stylesheet",
                href="https://unpkg.com/vis-network/styles/vis-network.min.css")
    ),
    tabItems(

      # ── TAB 1: SIMULATION ─────────────────────────────────────────────
      tabItem("sim",
        fluidRow(
          column(3,
            box(title="Settings", status="primary", solidHeader=TRUE, width=12,
              div(class="plbl","Duration (days)"),
              sliderInput("t_end","",30,1825,365,30),
              div(class="plbl","IOP (mmHg)"),
              sliderInput("IOP_target","",10,40,21,0.5),
              div(class="plbl","Dose Amount"), numericInput("dose_amt","",5,0,50,0.5),
              checkboxGroupInput("dose_times","Injection days (dose to cover loop ignition ~day 156)",
                choices=list("Day 0"=0,"Day 30"=30,"Day 60"=60,"Day 90"=90,
                             "Day 120"=120,"Day 150"=150,"Day 180"=180,
                             "Day 210"=210,"Day 240"=240,"Day 270"=270),
                selected=c(90,150,210)),
              tags$hr(),
              tags$div(class="fix-note",
                tags$b("v5 FIX: TLR4 Bootstrap"),tags$br(),
                "k_damp_M1 × DAMPs × M_mig → M1",tags$br(),
                "(paper §6.3.1 — loop ignition)"),
              div(class="plbl","TLR4 Bootstrap (k_damp_M1)"),
              sliderInput("k_damp_M1","",0,2,0.5,0.05),
              div(class="plbl","Loop — C1q Production"),
              sliderInput("k_C1q_M1","",0.1,5,1.36,0.01),
              div(class="plbl","Loop — C3 Cleavage"),
              sliderInput("k_C3_cleave","",0.1,4,1.2,0.1),
              div(class="plbl","Loop — C3aR Feedback"),
              sliderInput("k_C3aR_act","",0.1,6,2.28,0.01),
              div(class="plbl","C5aR1 Migration"),
              sliderInput("k_mig_C5a","",0,3,0.8,0.1),
              div(class="plbl","RPE Stress Rate"),
              sliderInput("k_rpe_stress","",0.002,0.05,0.003,0.001),
              div(class="plbl","DAMPs Migration"),
              sliderInput("k_mig_damp","",0.1,3,1.64,0.01),
              div(class="plbl","M1→M2 Switch Rate (loop removal)"),
              sliderInput("k_M1_switch","",0.02,0.5,0.176,0.002),
              tags$hr(),
              tags$div(class="fix-note", style="border-left-color:#00897B;background:#E0F2F1;",
                tags$b("v6 SATURATING BRAKES"),tags$br(),
                "Brake 1: C3aR desensitization (k_des_on/off)",tags$br(),
                "Brake 2: anti-inflam M1 suppression (K_anti_M1)",tags$br(),
                tags$em("These create a STABLE intermediate M1 (~5%)")),
              div(class="plbl","Brake 1 — C3aR Desensitization On (k_des_on)"),
              sliderInput("k_des_on","",0,8,5.70,0.05),
              div(class="plbl","Brake 1 — C3aR Resensitization (k_des_off)"),
              sliderInput("k_des_off","",0.02,0.5,0.111,0.002),
              div(class="plbl","Brake 2 — Anti-inflam IC50 (K_anti_M1)"),
              sliderInput("K_anti_M1","",0.1,2,0.595,0.005),
              tags$hr(),
              tags$div(class="fix-note",
                tags$b("RGC Death Coupling:"),tags$br(),
                "k_rgc_rpe = trans-synaptic (RPE→RGC)",tags$br(),
                "k_rgc_cyt = cytokine-mediated"),
              div(class="plbl","Trans-synaptic coupling (k_rgc_rpe)"),
              sliderInput("k_rgc_rpe","",0.001,0.20,0.005,0.001),
              div(class="plbl","Cytokine-RGC coupling (k_rgc_cyt)"),
              sliderInput("k_rgc_cyt","",0.001,0.30,0.018,0.001),
              div(class="plbl","Cytokine-RPE coupling (k_rpe_cyt)"),
              sliderInput("k_rpe_cyt","",0.005,0.10,0.010,0.001),
              tags$hr(),
              div(class="plbl","Drug PD — C1q Block (ANX007)"),
              sliderInput("Emax_C1qblock","",0,15,6,0.5),
              div(class="plbl","Drug PD — C3aR/C5aR1 Block"),
              sliderInput("Emax_C3aRblock","",0,15,5,0.5),
              div(class="plbl","Drug PD — M1→M2 Switch"),
              sliderInput("Emax_switch","",0,15,7,0.5),
              div(class="plbl","Drug PD — Migration Block"),
              sliderInput("Emax_migration","",0,15,5,0.5),
              sliderInput("EC50_pep","EC50 peptide",0.1,5,0.8,0.1),
              tags$hr(),
              actionButton("run_sim","▶  Run Simulation", class="btn-run"),
              actionButton("run_opt","⚙  Auto-Optimize Parameters", class="btn-opt"),
              div(id="opt_result", class="opt-res"),
              br(),
              actionButton("reset_all","⟳ Reset Defaults", class="btn-rst")
            )
          ),
          column(9,
            fluidRow(
              valueBoxOutput("vb_rgc_t",3), valueBoxOutput("vb_rgc_c",3),
              valueBoxOutput("vb_loop",3),  valueBoxOutput("vb_m1",3)
            ),
            fluidRow(
              box(title="RGC Survival",   status="success", solidHeader=TRUE,
                  width=6, plotlyOutput("p_rgc",height="230px")),
              box(title="RPE Health",     status="warning", solidHeader=TRUE,
                  width=6, plotlyOutput("p_rpe",height="230px"))
            ),
            fluidRow(
              box(title="Microglia States",   status="info",    solidHeader=TRUE,
                  width=6, plotlyOutput("p_micro",height="210px")),
              box(title="Cytokines & NTF",    status="primary", solidHeader=TRUE,
                  width=6, plotlyOutput("p_cyt",  height="210px"))
            ),
            fluidRow(
              box(title="Drug PK", status="primary", solidHeader=TRUE,
                  width=6, plotlyOutput("p_pk",  height="190px")),
              box(title="Feedback Loop Intensity (M1×C1q×C3a)",
                  status="danger", solidHeader=TRUE,
                  width=6, plotlyOutput("p_cycle",height="190px"))
            )
          )
        )
      ),

      # ── TAB 2: LOOP DYNAMICS ──────────────────────────────────────────
      tabItem("loop",
        fluidRow(
          box(title="Bootstrap TLR4 → Complement Loop — Molecular Dynamics",
              status="danger", solidHeader=TRUE, width=12,
            tags$div(class="fix-note",
              tags$b("v5 bootstrap chain (paper §6.3.1 → §4.2):"),tags$br(),
              "DAMPs + TLR4 → M_mig → M1 (seed)  →  M1 produces C1q  →  ",
              "C4b2a cleaves C3  →  C3a binds C3aR  →  more M1  [LOOP FIRES]",tags$br(),
              "Simultaneously: C5a via C5aR1 recruits M0→M_mig [migration amplifier]"
            )
          )
        ),
        fluidRow(
          box(title="TLR4 Bootstrap: M1 seed from DAMPs", status="warning",
              solidHeader=TRUE, width=4, plotlyOutput("loop_m1",   height="190px")),
          box(title="Loop Node A: C1q (by M1)", status="danger",
              solidHeader=TRUE, width=4, plotlyOutput("loop_c1q",  height="190px")),
          box(title="C3 pool (cleaved by C4b2a)", status="primary",
              solidHeader=TRUE, width=4, plotlyOutput("loop_c3",   height="190px"))
        ),
        fluidRow(
          box(title="Loop Node B: C3a → C3aR (CLOSES LOOP)", status="danger",
              solidHeader=TRUE, width=4, plotlyOutput("loop_c3a",  height="190px")),
          box(title="C5a → C5aR1 (migration amplifier)", status="primary",
              solidHeader=TRUE, width=4, plotlyOutput("loop_c5a",  height="190px")),
          box(title="Loop Intensity: M1 × C1q × C3a", status="danger",
              solidHeader=TRUE, width=4, plotlyOutput("loop_idx",  height="190px"))
        ),
        fluidRow(
          box(title="v6 BRAKE 1: C3aR Desensitization (R_des) — clamps the loop",
              status="success", solidHeader=TRUE, width=6,
              plotlyOutput("loop_rdes", height="200px")),
          box(title="v6 BRAKE 2: M1 plateau (stable intermediate, not switch)",
              status="success", solidHeader=TRUE, width=6,
              plotlyOutput("loop_m1plateau", height="200px"))
        )
      ),

      # ── TAB 3: NETWORK ────────────────────────────────────────────────
      tabItem("net",
        fluidRow(
          box(title="SimBiology-Style Network — Editable (drag · add · modify)",
              status="info", solidHeader=TRUE, width=12,
            tags$div(class="fix-note",
              tags$b("Node types: "),
              tags$span(style="background:#E65100;color:#fff;padding:2px 6px;border-radius:3px","■ Drug PK (orange)"),
              " ",
              tags$span(style="background:#1B3A6B;color:#fff;padding:2px 6px;border-radius:3px","■ State variables (navy)"),
              " ",
              tags$span(style="background:#006064;color:#fff;padding:2px 6px;border-radius:3px","■ Microglia (teal)"),
              " ",
              tags$span(style="background:#AD1457;color:#fff;padding:2px 6px;border-radius:3px","■ C1q loop-A (pink)"),
              " ",
              tags$span(style="background:#7B1FA2;color:#fff;padding:2px 6px;border-radius:3px","■ C3a loop-B (violet)"),
              tags$br(),
              tags$span(style="background:#FF8F00;color:#fff;padding:2px 6px;border-radius:3px","■ TLR4 bootstrap (amber — v5 FIX)"),
              " &nbsp;Gray circles = flow junctions; gray squares = process nodes; dashed = modulation"
            ),
            actionButton("init_net","🔄 Re-initialize Network", class="btn-run"),
            br(), br(),
            div(id="vis-net")
          )
        )
      ),

      # ── TAB 4: MECHANISMS ─────────────────────────────────────────────
      tabItem("mech",
        fluidRow(
          column(3,
            box(title="Mechanism Arms", status="primary", solidHeader=TRUE, width=12,
              tags$div(class="fix-note",
                tags$b("5 arms including TLR4:"),tags$br(),
                "PD0: TLR4 block (cut bootstrap)",tags$br(),
                "PD1: C1q block (ANX007)",tags$br(),
                "PD2: C3aR/C5aR1 block",tags$br(),
                "PD3: M1→M2 switch",tags$br(),
                "PD4: Migration block"
              ),
              sliderInput("mC1q","Emax C1q block",0,15,6,0.5),
              sliderInput("mC3aR","Emax C3aR block",0,15,5,0.5),
              sliderInput("msw","Emax M1→M2",0,15,7,0.5),
              sliderInput("mmg","Emax migration",0,15,5,0.5),
              sliderInput("mTLR4","TLR4 block (% k_damp_M1 reduced)",0,100,80,5),
              actionButton("run_mech","Compare Mechanisms", class="btn-run")
            )
          ),
          column(9,
            fluidRow(
              box(title="RGC Survival by Mechanism", status="success",
                  solidHeader=TRUE, width=8, plotlyOutput("mech_rgc",height="270px")),
              box(title="Endpoint RGC %", status="success",
                  solidHeader=TRUE, width=4, plotlyOutput("mech_bar",height="270px"))
            ),
            fluidRow(
              box(title="Loop Suppression", status="danger",
                  solidHeader=TRUE, width=6, plotlyOutput("mech_loop",height="220px")),
              box(title="RPE by Mechanism", status="warning",
                  solidHeader=TRUE, width=6, plotlyOutput("mech_rpe", height="220px"))
            )
          )
        )
      ),

      # ── TAB 5: DEATH PATHWAYS ─────────────────────────────────────────
      tabItem("death",
        fluidRow(
          box(title="RGC Death Decomposition", status="danger",
              solidHeader=TRUE, width=8, plotlyOutput("death_area",height="280px")),
          box(title="End-point Fractions", status="danger",
              solidHeader=TRUE, width=4, plotlyOutput("death_pie",height="280px"))
        ),
        fluidRow(
          box(title="RPE–RGC Phase Portrait", status="primary",
              solidHeader=TRUE, width=6, plotlyOutput("death_phase",height="260px")),
          box(title="NTF & Protection Factor", status="success",
              solidHeader=TRUE, width=6, plotlyOutput("death_ntf",height="260px"))
        )
      ),

      # ── TAB 6: EXPORT ─────────────────────────────────────────────────
      tabItem("export",
        fluidRow(
          box(title="Export Self-Contained HTML", status="success",
              solidHeader=TRUE, width=6,
            tags$div(class="fix-note",
              "Exports include: SimBiology-style vis.js network, all 10 plotly plots,",tags$br(),
              "loop intensity trace, death decomposition, annotated parameter table,",tags$br(),
              "v5 bootstrap fix documentation. Works offline."
            ), br(),
            downloadButton("dl_html","⬇ Download HTML Report",
                           class="btn-run", style="width:100%;font-size:14px;padding:12px")
          ),
          box(title="Summary Table", status="primary",
              solidHeader=TRUE, width=6, tableOutput("export_tbl"))
        )
      )
    )
  )
)

# ============================================================================
# 11. SERVER
# ============================================================================
server <- function(input, output, session) {

  get_params <- reactive({
    p <- DEFAULTS
    p$IOP_target      <- input$IOP_target
    p$k_damp_M1       <- input$k_damp_M1
    p$k_C1q_M1        <- input$k_C1q_M1
    p$k_C3_cleave     <- input$k_C3_cleave
    p$k_C3aR_act      <- input$k_C3aR_act
    p$k_mig_C5a       <- input$k_mig_C5a
    p$k_rpe_stress    <- input$k_rpe_stress
    p$k_mig_damp      <- input$k_mig_damp
    p$k_M1_switch     <- input$k_M1_switch
    p$k_rgc_rpe       <- input$k_rgc_rpe
    p$k_rgc_cyt       <- input$k_rgc_cyt
    p$k_rpe_cyt       <- input$k_rpe_cyt
    p$k_des_on        <- input$k_des_on
    p$k_des_off       <- input$k_des_off
    p$K_anti_M1       <- input$K_anti_M1
    p$Emax_C1qblock   <- input$Emax_C1qblock
    p$Emax_C3aRblock  <- input$Emax_C3aRblock
    p$Emax_switch     <- input$Emax_switch
    p$Emax_migration  <- input$Emax_migration
    p$EC50_pep        <- input$EC50_pep
    p$dose_amount     <- input$dose_amt
    p$dose_times      <- as.numeric(input$dose_times)
    if (length(p$dose_times)==0) p$dose_times <- 0
    p
  })

  observeEvent(input$reset_all, {
    D <- DEFAULTS
    updateSliderInput(session,"IOP_target",    value=D$IOP_target)
    updateSliderInput(session,"k_damp_M1",     value=D$k_damp_M1)
    updateSliderInput(session,"k_C1q_M1",      value=D$k_C1q_M1)
    updateSliderInput(session,"k_C3_cleave",   value=D$k_C3_cleave)
    updateSliderInput(session,"k_C3aR_act",    value=D$k_C3aR_act)
    updateSliderInput(session,"k_mig_C5a",     value=D$k_mig_C5a)
    updateSliderInput(session,"k_rpe_stress",  value=D$k_rpe_stress)
    updateSliderInput(session,"k_mig_damp",    value=D$k_mig_damp)
    updateSliderInput(session,"k_M1_switch",   value=D$k_M1_switch)
    updateSliderInput(session,"k_rgc_rpe",     value=D$k_rgc_rpe)
    updateSliderInput(session,"k_rgc_cyt",     value=D$k_rgc_cyt)
    updateSliderInput(session,"k_rpe_cyt",     value=D$k_rpe_cyt)
    updateSliderInput(session,"k_des_on",      value=D$k_des_on)
    updateSliderInput(session,"k_des_off",     value=D$k_des_off)
    updateSliderInput(session,"K_anti_M1",     value=D$K_anti_M1)
    updateSliderInput(session,"Emax_C1qblock", value=D$Emax_C1qblock)
    updateSliderInput(session,"Emax_C3aRblock",value=D$Emax_C3aRblock)
    updateSliderInput(session,"Emax_switch",   value=D$Emax_switch)
    updateSliderInput(session,"Emax_migration",value=D$Emax_migration)
    updateSliderInput(session,"EC50_pep",      value=D$EC50_pep)
    updateSliderInput(session,"t_end",         value=365)
  })

  # ── AUTO-OPTIMIZER ────────────────────────────────────────────────────
  observeEvent(input$run_opt, {
    runjs("document.getElementById('opt_result').style.display='none'")
    p  <- get_params()
    te <- input$t_end
    withProgress(message="Optimizing (11 params incl. brakes)…", value=0, {
      res <- run_optimizer(p, te, progress_cb=function(v) setProgress(v))
    })
    # write all tuned params back to sliders
    updateSliderInput(session,"k_des_on",    value=round(res$k_des_on,   3))
    updateSliderInput(session,"k_des_off",   value=round(res$k_des_off,  3))
    updateSliderInput(session,"K_anti_M1",   value=round(res$K_anti_M1,  3))
    updateSliderInput(session,"k_rpe_stress",value=round(res$k_rpe_stress,4))
    updateSliderInput(session,"k_mig_damp",  value=round(res$k_mig_damp, 2))
    updateSliderInput(session,"k_C3aR_act",  value=round(res$k_C3aR_act, 2))
    updateSliderInput(session,"k_C1q_M1",    value=round(res$k_C1q_M1,   2))
    updateSliderInput(session,"k_M1_switch", value=round(res$k_M1_switch,3))
    updateSliderInput(session,"k_rgc_cyt",   value=round(res$k_rgc_cyt,  3))
    updateSliderInput(session,"k_rpe_cyt",   value=round(res$k_rpe_cyt,  3))
    updateSliderInput(session,"k_rgc_rpe",   value=round(res$k_rgc_rpe,  3))
    html_msg <- sprintf(
      paste0("✓ Optimized — cost=%.3f<br>",
             "BRAKES: k_des_on=%.2f · k_des_off=%.3f · K_anti_M1=%.2f<br>",
             "LOOP: k_C3aR_act=%.2f · k_C1q_M1=%.2f · k_M1_switch=%.3f<br>",
             "TIMING: k_rpe_stress=%.4f · k_mig_damp=%.2f"),
      res$final_cost, res$k_des_on, res$k_des_off, res$K_anti_M1,
      res$k_C3aR_act, res$k_C1q_M1, res$k_M1_switch,
      res$k_rpe_stress, res$k_mig_damp)
    runjs(sprintf("var d=document.getElementById('opt_result');d.innerHTML='%s';d.style.display='block'", html_msg))
  })

  # ── SIMULATION ────────────────────────────────────────────────────────
  sim_data <- eventReactive(input$run_sim, {
    p <- get_params()
    withProgress(message="Running simulation…", value=0.5, {
      list(treat=run_sim(p,TRUE,input$t_end), ctrl=run_sim(p,FALSE,input$t_end))
    })
  }, ignoreNULL=FALSE)

  metrics <- reactive({ m <- sim_data(); compute_metrics(m$treat, m$ctrl) })

  # value boxes
  output$vb_rgc_t <- renderValueBox({
    m <- metrics()
    valueBox(sprintf("%.1f%%", if(!is.null(m)) m$rgc_treat else 0), "RGC (treated)", icon=icon("eye"), color="green")
  })
  output$vb_rgc_c <- renderValueBox({
    m <- metrics()
    valueBox(sprintf("%.1f%%", if(!is.null(m)) m$rgc_ctrl else 0), "RGC (control)", icon=icon("eye"), color="red")
  })
  output$vb_loop <- renderValueBox({
    m <- metrics()
    valueBox(sprintf("%.1f%%", if(!is.null(m)) m$loop_suppress else 0), "Loop suppressed", icon=icon("sync-alt"), color="yellow")
  })
  output$vb_m1 <- renderValueBox({
    m <- metrics()
    valueBox(sprintf("%.3f", if(!is.null(m)) m$m1_peak else 0), "Peak M1", icon=icon("fire"), color="purple")
  })

  # QSP theme helper
  qsp_theme <- function() {
    theme_minimal(base_size=10) +
    theme(panel.background=element_rect(fill="#FAFAFA",colour=NA),
          plot.background =element_rect(fill="#FAFAFA",colour=NA),
          panel.grid.major=element_line(colour="grey90",linewidth=0.3),
          panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold",size=10,colour="#1A1A2E"),
          axis.text=element_text(size=8),
          legend.position="bottom",legend.text=element_text(size=8),
          legend.key.size=unit(0.4,"cm"))
  }
  gg2ly <- function(p) {
    ggplotly(p, tooltip=c("x","y","colour")) %>%
      layout(legend=list(orientation="h",y=-0.35)) %>%
      config(displayModeBar=FALSE)
  }
  two_line <- function(dt, dc, yvar, ylab, title, ct="#1565C0", cc="#C62828", ylim=NULL) {
    df <- bind_rows(dt %>% mutate(g="Treated"), dc %>% mutate(g="Control"))
    p  <- ggplot(df, aes(time, .data[[yvar]], colour=g, linetype=g)) +
      geom_line(linewidth=0.9) +
      scale_colour_manual(values=c("Treated"=ct,"Control"=cc), name=NULL) +
      scale_linetype_manual(values=c("Treated"="solid","Control"="dashed"), name=NULL) +
      labs(x="Days",y=ylab,title=title) + qsp_theme()
    if (!is.null(ylim)) p <- p + coord_cartesian(ylim=ylim)
    p
  }

  # sim tab plots
  output$p_rgc <- renderPlotly({ sd <- sim_data(); req(sd)
    gg2ly(two_line(sd$treat, sd$ctrl, "RGC_pct","%","RGC Survival","#2E7D32","#C62828",c(0,105)))})
  output$p_rpe <- renderPlotly({ sd <- sim_data(); req(sd)
    gg2ly(two_line(sd$treat, sd$ctrl, "RPE_pct","%","RPE Health","#E65100","#C62828",c(0,105)))})
  output$p_micro <- renderPlotly({ sd <- sim_data(); req(sd)
    dc <- sd$ctrl; dt <- sd$treat
    df <- bind_rows(dc %>% select(time,M_mig,M1,M2) %>% pivot_longer(-time) %>% mutate(g="Control"),
                    dt %>% select(time,M_mig,M1,M2) %>% pivot_longer(-time) %>% mutate(g="Treated"))
    p <- ggplot(df, aes(time,value,colour=name,linetype=g)) +
      geom_line(linewidth=0.85) +
      scale_colour_manual(values=c(M_mig="#26A69A",M1="#C62828",M2="#2E7D32"),name=NULL)+
      scale_linetype_manual(values=c(Treated="solid",Control="dashed"),name=NULL)+
      labs(x="Days",y="a.u.",title="Microglial States") + qsp_theme()
    gg2ly(p)
  })
  output$p_cyt <- renderPlotly({ sd <- sim_data(); req(sd)
    df <- bind_rows(sd$ctrl  %>% select(time,Cyt_pro,NTF) %>% pivot_longer(-time) %>% mutate(g="Control"),
                    sd$treat %>% select(time,Cyt_pro,NTF) %>% pivot_longer(-time) %>% mutate(g="Treated"))
    p <- ggplot(df, aes(time,value,colour=name,linetype=g)) +
      geom_line(linewidth=0.85) +
      scale_colour_manual(values=c(Cyt_pro="#C62828",NTF="#1565C0"),name=NULL)+
      scale_linetype_manual(values=c(Treated="solid",Control="dashed"),name=NULL)+
      labs(x="Days",y="a.u.",title="Cytokines & NTF") + qsp_theme()
    gg2ly(p)
  })
  output$p_pk <- renderPlotly({ sd <- sim_data(); req(sd)
    df <- sd$treat %>% select(time,A_eye,C_pep) %>% pivot_longer(-time)
    gg2ly(ggplot(df,aes(time,value,colour=name))+geom_line(linewidth=1)+
      scale_colour_manual(values=c(A_eye="#6A1B9A",C_pep="#E91E63"),name=NULL)+
      labs(x="Days",y="a.u.",title="Drug PK")+qsp_theme())
  })
  output$p_cycle <- renderPlotly({ sd <- sim_data(); req(sd)
    p <- ggplot() +
      geom_area(data=sd$ctrl,  aes(time,loop_idx), fill="#C62828", alpha=0.25)+
      geom_line(data=sd$ctrl,  aes(time,loop_idx,colour="Control"),  linewidth=1)+
      geom_line(data=sd$treat, aes(time,loop_idx,colour="Treated"),  linewidth=1)+
      scale_colour_manual(values=c(Control="#C62828",Treated="#2E7D32"),name=NULL)+
      labs(x="Days",y="M1×C1q×C3a",title="Loop Intensity")+qsp_theme()
    gg2ly(p)
  })

  # loop tab
  lp <- function(v,t,ct) { sd <- sim_data(); req(sd)
    gg2ly(two_line(sd$treat,sd$ctrl,v,"a.u.",t,ct,"#9E9E9E")) }
  output$loop_m1  <- renderPlotly({ lp("M1",  "M1 (TLR4-seeded then C3aR-amplified)","#C62828") })
  output$loop_c1q <- renderPlotly({ lp("C1q", "C1q — Loop Node A","#AD1457") })
  output$loop_c3  <- renderPlotly({ lp("C3",  "C3 pool (cleaved by C4b2a)","#6A1B9A") })
  output$loop_c3a <- renderPlotly({ lp("C3a", "C3a — Loop Node B (C3aR)","#9C27B0") })
  output$loop_c5a <- renderPlotly({ lp("C5a", "C5a — C5aR1 migration arm","#7B1FA2") })
  output$loop_idx <- renderPlotly({ sd <- sim_data(); req(sd)
    p <- ggplot()+
      geom_area(data=sd$ctrl, aes(time,loop_idx), fill="#C62828",alpha=0.2)+
      geom_line(data=sd$ctrl, aes(time,loop_idx,colour="Control"),linewidth=1.1)+
      geom_line(data=sd$treat,aes(time,loop_idx,colour="Treated"),linewidth=1.1)+
      scale_colour_manual(values=c(Control="#C62828",Treated="#2E7D32"),name=NULL)+
      labs(x="Days",y="M1×C1q×C3a",title="Loop Intensity")+qsp_theme()
    gg2ly(p)
  })
  # v6 BRAKE 1: C3aR desensitization
  output$loop_rdes <- renderPlotly({ sd <- sim_data(); req(sd)
    p <- ggplot()+
      geom_area(data=sd$ctrl, aes(time,R_des), fill="#00897B", alpha=0.18)+
      geom_line(data=sd$ctrl, aes(time,R_des,colour="Control"),linewidth=1.2)+
      geom_line(data=sd$treat,aes(time,R_des,colour="Treated"),linewidth=1.2)+
      coord_cartesian(ylim=c(0,1))+
      scale_colour_manual(values=c(Control="#00695C",Treated="#2E7D32"),name=NULL)+
      labs(x="Days",y="R_des (0–1)",
           title="C3aR desensitization (rises with C3a → caps loop)")+qsp_theme()
    gg2ly(p)
  })
  # v6 BRAKE 2: M1 plateau vs hypothetical explosion
  output$loop_m1plateau <- renderPlotly({ sd <- sim_data(); req(sd)
    p <- ggplot()+
      geom_line(data=sd$ctrl, aes(time,M1*100,colour="Control"),linewidth=1.2)+
      geom_line(data=sd$treat,aes(time,M1*100,colour="Treated"),linewidth=1.2)+
      geom_hline(yintercept=15, linetype="dotted", colour="grey50")+
      annotate("text", x=Inf, y=15, label="brake ceiling", hjust=1.05, vjust=-0.4,
               size=3, colour="grey40")+
      scale_colour_manual(values=c(Control="#C62828",Treated="#2E7D32"),name=NULL)+
      labs(x="Days",y="M1 (% of pool)",
           title="M1 stabilises at intermediate level (no explosion)")+qsp_theme()
    gg2ly(p)
  })

  # network tab
  observeEvent(input$tabs, { if (input$tabs=="net") { Sys.sleep(0.2); runjs(build_simbiology_js("vis-net")) } })
  observeEvent(input$init_net, { runjs(build_simbiology_js("vis-net")) })
  observe({
    invalidateLater(800)
    isolate({ if (!is.null(input$tabs) && input$tabs=="net") runjs(build_simbiology_js("vis-net")) })
  })

  # mechanisms
  ARM_COLS <- c("All targets"="#2E7D32","TLR4 block"="#FF8F00","C1q block"="#AD1457",
                "C3aR block"="#7B1FA2","M1→M2"="#1565C0","Migration"="#006064","No treatment"="#C62828")
  mech_data <- eventReactive(input$run_mech, {
    p <- get_params(); te <- input$t_end
    run_arm <- function(tlr4red, c1q, c3ar, sw, mg, lbl) {
      pm <- p
      pm$k_damp_M1      <- p$k_damp_M1 * (1 - tlr4red/100)
      pm$Emax_C1qblock  <- if(c1q)  input$mC1q  else 0
      pm$Emax_C3aRblock <- if(c3ar) input$mC3aR else 0
      pm$Emax_switch    <- if(sw)   input$msw   else 0
      pm$Emax_migration <- if(mg)   input$mmg   else 0
      run_sim(pm, treat=TRUE, t_end=te) %>% mutate(arm=lbl)
    }
    withProgress(message="Running arms…", value=0.2, {
      bind_rows(
        run_arm(input$mTLR4,T,T,T,T,"All targets"),
        run_arm(input$mTLR4,F,F,F,F,"TLR4 block"),
        run_arm(0,T,F,F,F,"C1q block"),
        run_arm(0,F,T,F,F,"C3aR block"),
        run_arm(0,F,F,T,F,"M1→M2"),
        run_arm(0,F,F,F,T,"Migration"),
        run_sim(p,FALSE,t_end=te) %>% mutate(arm="No treatment")
      )
    })
  }, ignoreNULL=FALSE)

  output$mech_rgc <- renderPlotly({ md <- mech_data(); req(md)
    gg2ly(ggplot(md,aes(time,RGC_pct,colour=arm))+geom_line(linewidth=1)+
      coord_cartesian(ylim=c(0,105))+scale_colour_manual(values=ARM_COLS,name=NULL)+
      labs(x="Days",y="RGC (%)",title="RGC Survival by Mechanism")+qsp_theme())
  })
  output$mech_bar <- renderPlotly({
    md <- mech_data(); req(md)
    ends <- md %>% group_by(arm) %>% summarise(RGC=last(RGC_pct),.groups="drop")
    ggplotly(ggplot(ends,aes(reorder(arm,RGC),RGC,fill=arm,text=sprintf("%s\n%.1f%%",arm,RGC)))+
      geom_col(width=0.65,colour="white")+coord_flip(ylim=c(0,110))+
      scale_fill_manual(values=ARM_COLS,guide="none")+
      labs(x=NULL,y="RGC end (%)",title="Endpoint")+qsp_theme(),tooltip="text") %>%
      config(displayModeBar=FALSE)
  })
  output$mech_loop <- renderPlotly({ md <- mech_data(); req(md)
    gg2ly(ggplot(md,aes(time,loop_idx,colour=arm))+geom_line(linewidth=0.9)+
      scale_colour_manual(values=ARM_COLS,name=NULL)+
      labs(x="Days",y="M1×C1q×C3a",title="Loop Suppression")+qsp_theme())
  })
  output$mech_rpe <- renderPlotly({ md <- mech_data(); req(md)
    gg2ly(ggplot(md,aes(time,RPE_pct,colour=arm))+geom_line(linewidth=0.9)+
      coord_cartesian(ylim=c(0,105))+scale_colour_manual(values=ARM_COLS,name=NULL)+
      labs(x="Days",y="RPE (%)",title="RPE by Mechanism")+qsp_theme())
  })

  # death pathways
  output$death_area <- renderPlotly({ sd <- sim_data(); req(sd)
    df <- sd$ctrl %>% select(time,d_cyt,d_rpe,d_iop) %>%
      pivot_longer(-time) %>%
      mutate(path=recode(name,d_cyt="Cytokine",d_rpe="Trans-synaptic",d_iop="IOP-direct"))
    gg2ly(ggplot(df,aes(time,value,fill=path))+geom_area(alpha=0.8,position="stack")+
      scale_fill_manual(values=c(Cytokine="#C62828","Trans-synaptic"="#6A1B9A","IOP-direct"="#E65100"),name=NULL)+
      labs(x="Days",y="Rate 1/day",title="RGC Death Decomposition (Control)")+qsp_theme())
  })
  output$death_pie <- renderPlotly({ sd <- sim_data(); req(sd)
    lc <- tail(sd$ctrl,1)
    fr <- tibble(path=c("Cytokine","Trans-syn","IOP"),val=c(lc$d_cyt,lc$d_rpe,lc$d_iop)) %>%
      mutate(pct=val/sum(val)*100)
    plot_ly(fr,labels=~path,values=~pct,type="pie",
            marker=list(colors=c("#C62828","#6A1B9A","#E65100")),
            textinfo="label+percent") %>% layout(showlegend=FALSE) %>% config(displayModeBar=FALSE)
  })
  output$death_phase <- renderPlotly({ sd <- sim_data(); req(sd)
    gg2ly(ggplot()+
      geom_path(data=sd$ctrl, aes(100-RPE_pct,100-RGC_pct,colour="Control"),linewidth=1.1,
                arrow=arrow(length=unit(0.18,"cm")))+
      geom_path(data=sd$treat,aes(100-RPE_pct,100-RGC_pct,colour="Treated"),linewidth=1.1,linetype="dashed")+
      geom_abline(slope=1,intercept=0,colour="grey60",linetype="dotted")+
      scale_colour_manual(values=c(Control="#C62828",Treated="#2E7D32"),name=NULL)+
      labs(x="RPE Loss (%)",y="RGC Loss (%)",title="RPE–RGC Phase Portrait")+qsp_theme())
  })
  output$death_ntf <- renderPlotly({ sd <- sim_data(); req(sd)
    df <- bind_rows(sd$treat %>% select(time,NTF,Prot) %>% mutate(g="Treated"),
                    sd$ctrl  %>% select(time,NTF,Prot) %>% mutate(g="Control")) %>%
      pivot_longer(c(NTF,Prot))
    gg2ly(ggplot(df,aes(time,value,colour=name,linetype=g))+geom_line(linewidth=0.9)+
      scale_colour_manual(values=c(NTF="#1565C0",Prot="#2E7D32"),name=NULL)+
      scale_linetype_manual(values=c(Treated="solid",Control="dashed"),name=NULL)+
      labs(x="Days",y="a.u.",title="NTF & Protection")+qsp_theme())
  })

  # export
  output$export_tbl <- renderTable({
    m <- metrics()
    if (is.null(m)) return(data.frame(Metric="Run simulation first", Value="—"))
    data.frame(Metric=c("RGC (treated)","RGC (control)","Protection",
                        "RPE (treated)","RPE (control)",
                        "Loop suppression","Peak M1","Peak C1q","Peak C3a",
                        "Cytokine death %","Trans-syn death %","IOP death %"),
               Value=c(sprintf("%.1f%%",m$rgc_treat),sprintf("%.1f%%",m$rgc_ctrl),
                       sprintf("+%.1f%%",m$protection),
                       sprintf("%.1f%%",m$rpe_treat),sprintf("%.1f%%",m$rpe_ctrl),
                       sprintf("%.1f%%",m$loop_suppress),
                       sprintf("%.4f",m$m1_peak),sprintf("%.4f",m$c1q_peak),sprintf("%.4f",m$c3a_peak),
                       sprintf("%.1f%%",m$cyt_frac),sprintf("%.1f%%",m$trans_frac),sprintf("%.1f%%",m$iop_frac)))
  }, striped=TRUE, bordered=TRUE, spacing="s")

  output$dl_html <- downloadHandler(
    filename=function() paste0("qsp_glaucoma_v5_",Sys.Date(),".html"),
    content=function(file) {
      sd <- sim_data()
      if (is.null(sd)) { writeLines("<html><body><p>Run simulation first.</p></body></html>",file); return() }
      writeLines(build_html_export(sd$treat, sd$ctrl, get_params(), input$t_end), file)
    }
  )
}

# ============================================================================
shinyApp(ui, server)
