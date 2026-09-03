###############################################################################
# 04_sobol_N512_refined.R  –  Sobol SA at N=512 with the REFINED extractor
#
# Uses the same extract_sa_outputs (with inflation outputs) as the N=2048 run.
# Then builds a proper comparison table: sobol_comparison_N512_vs_N2048.csv
#
# Outputs:
#   tables_final/sobol_indices_N512_refined.csv
#   outputs_refined/sobol_raw_N512.rds
#   outputs_refined/sobol_batch_N512.rds
#   tables_final/sobol_comparison_N512_vs_N2048.csv  (overwritten)
#   logs_refined/sobol_N512.log
###############################################################################

cat("============================================================\n")
cat(" REFINED SOBOL SA  (N = 512)  +  COMPARISON TABLE\n")
cat("============================================================\n")

source("sensitivity_analysis/code_refined/00_helpers_refined.R")

sa_fig  <- "sensitivity_analysis/figures_final"
sa_tab  <- "sensitivity_analysis/tables_final"
sa_out  <- "sensitivity_analysis/outputs_refined"
sa_log  <- "sensitivity_analysis/logs_refined/sobol_N512.log"

cat(sprintf("Sobol SA N=512 (refined) started at %s\n", Sys.time()),
    file = sa_log)

# ---- 1. Compile model -------------------------------------------------------

cat(">> Compiling model...\n")
sys <- compile_model()
cat(">> Model compiled.\n")

# ---- 2. Parameter bounds (±30%) — identical to N=2048 run -------------------

phi1_lo <- BASELINE_PHI1 * 0.70
phi1_hi <- BASELINE_PHI1 * 1.30
phi2_lo <- BASELINE_PHI2 * 0.70
phi2_hi <- BASELINE_PHI2 * 1.30

params_sa <- c("phi1", "phi2")

# ---- 3. Output names — identical to N=2048 run ------------------------------

output_names <- c(
  "MC_2100", "tropical_stock_2100", "sequestration_2100",
  "net_emissions_2100", "temperature_2100", "output_2100",
  "private_debt_ratio_2100", "cumulative_MC_2020_2100", "cumulative_emissions",
  "inflation_2100", "avg_inflation_2050_2100", "min_inflation_post2050"
)

# ---- 4. Sobol run at N = 512 ------------------------------------------------

N_512 <- 512  # 512 * (2+2+1) = 2560 evaluations
cat(sprintf(">> N = %d  =>  %d total evaluations\n", N_512, N_512 * 5))

set.seed(42)
mat_512 <- sobol_matrices(
  matrices = c("A", "B", "AB"),
  N        = N_512,
  params   = params_sa,
  order    = "second",
  type     = "QRN"
)

# Scale to physical bounds
mat_512_phys <- mat_512
mat_512_phys[, "phi1"] <- phi1_lo + mat_512[, "phi1"] * (phi1_hi - phi1_lo)
mat_512_phys[, "phi2"] <- phi2_lo + mat_512[, "phi2"] * (phi2_hi - phi2_lo)

param_df_512 <- data.frame(
  phi1 = mat_512_phys[, "phi1"],
  phi2 = mat_512_phys[, "phi2"]
)

cat(">> Running N=512 Sobol batch...\n")
t0 <- Sys.time()
res_512 <- run_batch(sys, param_df_512, log_file = sa_log, progress = TRUE)
t1 <- Sys.time()
elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
cat(sprintf(">> Completed: %d evaluations in %.1f s (%.4f s/run)\n",
            nrow(param_df_512), elapsed, elapsed / nrow(param_df_512)))
cat(sprintf("N=512 run: %d evals in %.1f s\n", nrow(param_df_512), elapsed),
    file = sa_log, append = TRUE)

n_fail <- sum(is.na(res_512$MC_2100))
cat(sprintf(">> Failures: %d / %d\n", n_fail, nrow(res_512)))

# ---- 5. Compute Sobol indices for N=512 ------------------------------------

cat(sprintf("\n>> Computing Sobol indices for %d outputs...\n",
            length(output_names)))

sobol_results_512 <- list()

for (oname in output_names) {
  if (!oname %in% names(res_512)) {
    cat(sprintf("   [SKIP] %s: column not found in results\n", oname))
    next
  }

  y_vec <- res_512[[oname]]

  if (sum(is.finite(y_vec)) < 0.8 * length(y_vec)) {
    cat(sprintf("   [SKIP] %s: too many NAs (%.0f%%)\n",
                oname, 100 * sum(!is.finite(y_vec)) / length(y_vec)))
    next
  }

  y_clean <- y_vec
  y_clean[!is.finite(y_clean)] <- median(y_vec, na.rm = TRUE)

  sob <- sobol_indices(
    matrices = c("A", "B", "AB"),
    Y        = y_clean,
    N        = N_512,
    params   = params_sa,
    boot     = TRUE,
    R        = 1000,
    order    = "second",
    type     = "QRN"
  )

  sobol_results_512[[oname]] <- sob
  cat(sprintf("   [OK] %s\n", oname))
}

# ---- 6. Save raw results ----------------------------------------------------

saveRDS(sobol_results_512, file.path(sa_out, "sobol_raw_N512.rds"))
saveRDS(res_512,           file.path(sa_out, "sobol_batch_N512.rds"))
cat(">> Raw N=512 objects saved.\n")

# ---- 7. Tidy summary table --------------------------------------------------

dt_512 <- tidy_sobol(sobol_results_512)
write.csv(dt_512,
          file.path(sa_tab, "sobol_indices_N512_refined.csv"),
          row.names = FALSE)
cat(">> sobol_indices_N512_refined.csv saved.\n")

# ---- 8. Build comparison table (N=512 vs N=2048) ----------------------------

cat("\n>> Building comparison table (N=512 refined vs N=2048)...\n")

path_2048 <- file.path(sa_tab, "sobol_indices_N2048.csv")
if (!file.exists(path_2048)) {
  cat("!! sobol_indices_N2048.csv not found — skipping comparison.\n")
} else {
  dt_2048 <- data.table::fread(path_2048)

  dt_compare <- merge(
    dt_512[, .(output, parameters, sensitivity,
               original_512 = original,
               se_512       = std_error,
               lo_512       = low_ci,
               hi_512       = high_ci)],
    dt_2048[, .(output, parameters, sensitivity,
                original_2048 = original,
                se_2048       = std_error,
                lo_2048       = low_ci,
                hi_2048       = high_ci)],
    by = c("output", "parameters", "sensitivity"),
    all = TRUE
  )

  dt_compare$delta <- dt_compare$original_2048 - dt_compare$original_512
  dt_compare$se_reduction_pct <- ifelse(
    !is.na(dt_compare$se_512) & dt_compare$se_512 > 0,
    (1 - dt_compare$se_2048 / dt_compare$se_512) * 100,
    NA_real_
  )

  write.csv(dt_compare,
            file.path(sa_tab, "sobol_comparison_N512_vs_N2048.csv"),
            row.names = FALSE)
  cat(">> sobol_comparison_N512_vs_N2048.csv saved.\n")

  # Print Ti comparison
  cat("\n>> Ti comparison (N=512 refined vs N=2048):\n")
  dt_Ti <- dt_compare[dt_compare$sensitivity == "Ti",
                      c("output", "parameters", "original_512",
                        "original_2048", "delta", "se_reduction_pct")]
  print(dt_Ti, digits = 3)
}

cat(sprintf("\nSobol N=512 (refined) completed at %s\n", Sys.time()),
    file = sa_log, append = TRUE)
cat("\n>> 04_sobol_N512_refined.R finished.\n")
