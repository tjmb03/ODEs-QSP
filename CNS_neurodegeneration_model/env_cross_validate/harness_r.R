#!/usr/bin/env Rscript
# harness_r.R — R engine for cross-validation, with SEGMENT-RESTART dosing.
#
# Doses are applied by restarting the integration at each dose boundary
# (A_eye += dose at the segment start, so a dose day's output is POST-dose),
# exactly matching the Python engine. This deliberately avoids deSolve's
# event mechanism: events specified at times that coincide with output times
# are returned PRE-event and the dose enters the dynamics offset by a step,
# which shifted the drug trajectory ~1 day relative to Python.
#
# The CONTROL arm is a single segment, identical to a plain ode() call, so its
# (already machine-precision) agreement with Python is unchanged.
#
# Usage:
#   Rscript harness_r.R <control|treated> <out.csv> [path/to/glaucoma_model.R]

suppressMessages(library(deSolve))

args <- commandArgs(trailingOnly = TRUE)
arm  <- if (length(args) >= 1) args[1] else "control"
out  <- if (length(args) >= 2) args[2] else "r_engine.csv"
srcf <- if (length(args) >= 3) args[3] else "glaucoma_model.R"

source(srcf, local = TRUE)            # -> glaucoma_odes, DEFAULTS, SNAMES

p <- DEFAULTS
dose_times  <- c(0, 90, 150, 210)     # CANONICAL — identical across all engines
dose_amount <- p$dose_amount
C3_ss  <- p$k_C3_base  / p$k_C3_deg
NTF_ss <- p$k_ntf_base / p$k_deg_ntf

y <- setNames(numeric(length(SNAMES)), SNAMES)
y["RPE"] <- 1; y["M0"] <- p$M_total; y["C3"] <- C3_ss
y["NTF"] <- NTF_ss; y["RGC"] <- 1
# A_eye starts at 0; the t=0 loading dose is applied in the loop (like Python).

treated  <- identical(arm, "treated")
ds       <- if (treated) sort(unique(dose_times)) else numeric(0)
bounds   <- sort(unique(c(0, ds, 365)))
all_days <- 0:365

rec <- NULL
for (i in seq_len(length(bounds) - 1)) {
  t0 <- bounds[i]; t1 <- bounds[i + 1]
  if (treated && (t0 %in% ds)) y["A_eye"] <- y["A_eye"] + dose_amount  # dose at boundary
  seg_days <- all_days[all_days >= t0 & all_days < t1]                 # [t0, t1)
  if (i == length(bounds) - 1) seg_days <- c(seg_days, t1)            # final: include endpoint
  tt  <- sort(unique(c(t0, seg_days, t1)))                            # integrate t0 -> t1
  sol <- ode(y = y, times = tt, func = glaucoma_odes, parms = p,
             method = "lsoda", rtol = 1e-8, atol = 1e-10)
  rec <- rbind(rec, sol[sol[, "time"] %in% seg_days, , drop = FALSE])  # record requested days
  y   <- sol[sol[, "time"] == t1, SNAMES]                             # carry state to next dose
}

df <- as.data.frame(rec)              # columns: time + SNAMES (named)
write.csv(df, out, row.names = FALSE)
cat(sprintf("R/lsoda (segment-restart dosing) wrote %d rows x %d states (%s arm) to %s\n",
            nrow(df), length(SNAMES), arm, out))
