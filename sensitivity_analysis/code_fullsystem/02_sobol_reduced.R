###############################################################################
# 02_sobol_reduced.R  -  Stage 2: variance-based (Sobol) SA on the reduced set
#
# Uses the Morris screening to fix the 5 negligible parameters at baseline and
# computes first-order (Si) and total-order (Ti) Sobol indices for the 14
# influential parameters, in the focal scenario (forestA_50), for all 12
# system-level outputs. Jansen estimators via sensobol; bootstrap CIs.
#
# Convergence is checked by re-running at two base sample sizes (N1 < N2).
#
# Run from the PROJECT ROOT:
#   source("sensitivity_analysis/code_fullsystem/02_sobol_reduced.R")
###############################################################################

source("sensitivity_analysis/code_fullsystem/00_helpers_fullsystem.R")
suppressPackageStartupMessages(library(sensobol))

set.seed(20260730)

# ---- Config ----------------------------------------------------------------
# Parameters fixed at baseline (Morris: never in any output top-4; max mu* < 0.15)
FIXED_PARAMS <- c("CUP_pind", "g_pbs", "delta_g_sigma", "gamma_ast", "delta_CH")

N_PRIMARY <- 2048   # base sample; total evals = N*(k+2)
N_CONV    <- 1024   # convergence check
BOOT_R    <- 1000

out_dir <- "sensitivity_analysis/outputs_fullsystem"
fig_dir <- "sensitivity_analysis/figures_fullsystem"
tab_dir <- "sensitivity_analysis/tables_fullsystem"
log_dir <- "sensitivity_analysis/logs_fullsystem"
for (d in c(out_dir, fig_dir, tab_dir, log_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(log_dir, "sobol_run.log")
cat(sprintf("Reduced Sobol start: %s\n", Sys.time()), file = log_file)

# ---- Ledger + reduced parameter set ----------------------------------------
led    <- load_ledger()
keep   <- !(led$param %in% FIXED_PARAMS)
pnames <- led$param[keep]
lower  <- setNames(led$lower[keep], pnames)
upper  <- setNames(led$upper[keep], pnames)
k      <- length(pnames)

cat(sprintf(">> Reduced Sobol on %d parameters (fixed at baseline: %s)\n",
            k, paste(FIXED_PARAMS, collapse = ", ")))
cat(sprintf(">> Retained: %s\n", paste(pnames, collapse = ", ")))

# ---- Compile model once -----------------------------------------------------
cat(">> Compiling model under focal scenario (forestA_50)...\n")
sys <- compile_model_focal()

# ---- Sobol engine -----------------------------------------------------------
#' Run one Sobol experiment at base sample size N; returns per-output indices.
run_sobol_at <- function(N, tag) {
  cat(sprintf("\n>> Sobol at N=%d  (%d evaluations)\n", N, N * (k + 2)))
  mat <- sobol_matrices(matrices = c("A", "B", "AB"),
                        N = N, params = pnames, order = "first", type = "QRN")
  # Scale unit hypercube to physical bounds.
  Xphys <- mat
  for (p in pnames) Xphys[, p] <- lower[p] + mat[, p] * (upper[p] - lower[p])

  t0 <- Sys.time()
  Y  <- run_design(sys, Xphys, log_file = log_file, verbose = TRUE)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf(">> N=%d done: %d evals in %.1f s (%.4f s/run)\n",
              N, nrow(Xphys), el, el / nrow(Xphys)))

  saveRDS(list(N = N, mat = mat, Xphys = Xphys, Y = Y, params = pnames),
          file = file.path(out_dir, sprintf("sobol_raw_%s.rds", tag)))

  ind_list <- list()
  for (oname in FS_OUTPUT_NAMES) {
    y <- Y[, oname]
    if (sum(is.finite(y)) < 0.8 * length(y)) {
      cat(sprintf("   [skip] %s: too many non-finite\n", oname)); next
    }
    y[!is.finite(y)] <- median(y, na.rm = TRUE)
    sob <- sobol_indices(matrices = c("A", "B", "AB"), Y = y, N = N,
                         params = pnames, boot = TRUE, R = BOOT_R,
                         order = "first", type = "QRN")
    d <- sob$results
    d$output <- oname
    d$N      <- N
    ind_list[[oname]] <- d
  }
  do.call(rbind, ind_list)
}

ind_primary <- run_sobol_at(N_PRIMARY, "N2048")
ind_conv    <- run_sobol_at(N_CONV,    "N1024")

write.csv(ind_primary, file.path(tab_dir, "sobol_indices_N2048.csv"), row.names = FALSE)
write.csv(ind_conv,    file.path(tab_dir, "sobol_indices_N1024.csv"), row.names = FALSE)

# ---- Convergence: max |Si(N2048) - Si(N1024)| per output -------------------
conv_rows <- list()
for (oname in FS_OUTPUT_NAMES) {
  a <- ind_primary[ind_primary$output == oname & ind_primary$sensitivity == "Si", ]
  b <- ind_conv[ind_conv$output == oname & ind_conv$sensitivity == "Si", ]
  if (nrow(a) == 0 || nrow(b) == 0) next
  m <- merge(a[, c("parameters", "original")], b[, c("parameters", "original")],
             by = "parameters", suffixes = c("_N2048", "_N1024"))
  conv_rows[[oname]] <- data.frame(
    output = oname,
    max_abs_dSi  = max(abs(m$original_N2048 - m$original_N1024), na.rm = TRUE),
    mean_abs_dSi = mean(abs(m$original_N2048 - m$original_N1024), na.rm = TRUE)
  )
}
conv_tab <- do.call(rbind, conv_rows)
write.csv(conv_tab, file.path(tab_dir, "sobol_convergence_N1024_vs_N2048.csv"), row.names = FALSE)

cat("\n=== Sobol convergence (Si, N=1024 vs N=2048) ===\n")
print(conv_tab, digits = 3, row.names = FALSE)

# ---- Console summary: top total-order drivers per output -------------------
cat("\n=== Top total-order (Ti) drivers per output (N=2048) ===\n")
for (oname in FS_OUTPUT_NAMES) {
  s <- ind_primary[ind_primary$output == oname & ind_primary$sensitivity == "Ti", ]
  if (nrow(s) == 0) next
  s <- s[order(-s$original), ]
  cat(sprintf("%-26s: %s\n", oname,
      paste(sprintf("%s(%.2f)", head(s$parameters, 4), head(s$original, 4)), collapse = ", ")))
}

cat(sprintf("\nTables written to %s\n", tab_dir))
