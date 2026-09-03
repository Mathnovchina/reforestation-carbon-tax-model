###############################################################################
# 01_sobol_refined.R  –  Refined Sobol SA at N=2048 with inflation outputs
#
# Adds:  inflation_2100, avg_inflation_2050_2100, deflation_episodes
# Saves: tables_final/sobol_indices_N2048.csv
#        tables_final/sobol_comparison_N512_vs_N2048.csv
#        outputs_refined/sobol_raw_N2048.rds
#        outputs_refined/sobol_batch_N2048.rds
#        logs_refined/sobol_N2048.log
###############################################################################

cat("============================================================\n")
cat(" REFINED SOBOL SA  (N = 2048)\n")
cat("============================================================\n")

source("sensitivity_analysis/code_refined/00_helpers_refined.R")

sa_fig  <- "sensitivity_analysis/figures_final"
sa_tab  <- "sensitivity_analysis/tables_final"
sa_out  <- "sensitivity_analysis/outputs_refined"
sa_log  <- "sensitivity_analysis/logs_refined/sobol_N2048.log"

cat(sprintf("Sobol SA (refined) started at %s\n", Sys.time()), file = sa_log)

# ---- 1. Compile model -------------------------------------------------------

cat(">> Compiling model...\n")
sys <- compile_model()
cat(">> Model compiled.\n")

# ---- 2. Parameter bounds (±30%) --------------------------------------------

phi1_lo <- BASELINE_PHI1 * 0.70
phi1_hi <- BASELINE_PHI1 * 1.30
phi2_lo <- BASELINE_PHI2 * 0.70
phi2_hi <- BASELINE_PHI2 * 1.30

params_sa <- c("phi1", "phi2")

# ---- 3. Output names (extended with inflation) ------------------------------

output_names <- c(
  "MC_2100", "tropical_stock_2100", "sequestration_2100",
  "net_emissions_2100", "temperature_2100", "output_2100",
  "private_debt_ratio_2100", "cumulative_MC_2020_2100", "cumulative_emissions",
  # New inflation outputs
  "inflation_2100", "avg_inflation_2050_2100", "min_inflation_post2050"
)

# ---- 4. Full Sobol run at N = 2048 -----------------------------------------

N_full <- 2048  # 2048 * (2+2+1) = 10240 evaluations
cat(sprintf(">> N = %d  =>  %d total evaluations\n", N_full, N_full * 5))

set.seed(42)
mat_full <- sobol_matrices(
  matrices = c("A", "B", "AB"),
  N        = N_full,
  params   = params_sa,
  order    = "second",
  type     = "QRN"
)

# Scale to physical bounds
mat_full_phys <- mat_full
mat_full_phys[, "phi1"] <- phi1_lo + mat_full[, "phi1"] * (phi1_hi - phi1_lo)
mat_full_phys[, "phi2"] <- phi2_lo + mat_full[, "phi2"] * (phi2_hi - phi2_lo)

param_df_full <- data.frame(
  phi1 = mat_full_phys[, "phi1"],
  phi2 = mat_full_phys[, "phi2"]
)

cat(">> Running full Sobol batch...\n")
t0 <- Sys.time()
res_full <- run_batch(sys, param_df_full, log_file = sa_log, progress = TRUE)
t1 <- Sys.time()
elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
cat(sprintf(">> Full run completed: %d evaluations in %.1f s (%.4f s/run)\n",
            nrow(param_df_full), elapsed, elapsed / nrow(param_df_full)))
cat(sprintf("Full run: %d evals in %.1f s\n", nrow(param_df_full), elapsed),
    file = sa_log, append = TRUE)

n_fail <- sum(is.na(res_full$MC_2100))
cat(sprintf(">> Failures: %d / %d\n", n_fail, nrow(res_full)))

# ---- 5. Compute Sobol indices -----------------------------------------------

cat("\n>> Computing Sobol indices for %d outputs...\n", length(output_names))

sobol_results <- list()

for (oname in output_names) {
  if (!oname %in% names(res_full)) {
    cat(sprintf("   [SKIP] %s: column not found in results\n", oname))
    next
  }
  
  y_vec <- res_full[[oname]]
  
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
    N        = N_full,
    params   = params_sa,
    boot     = TRUE,
    R        = 1000,          # more bootstrap replicates for tighter CIs
    order    = "second",
    type     = "QRN"
  )
  
  sobol_results[[oname]] <- sob
  cat(sprintf("   [OK] %s\n", oname))
}

# ---- 6. Save raw results ----------------------------------------------------

saveRDS(sobol_results, file.path(sa_out, "sobol_raw_N2048.rds"))
saveRDS(res_full,      file.path(sa_out, "sobol_batch_N2048.rds"))
cat(">> Raw objects saved.\n")

# ---- 7. Tidy summary table --------------------------------------------------

dt_sobol <- tidy_sobol(sobol_results)
write.csv(dt_sobol, file.path(sa_tab, "sobol_indices_N2048.csv"), row.names = FALSE)
cat(">> Sobol indices table saved.\n")

# Print Ti summary for quick inspection
cat("\n>> Total-order indices (Ti):\n")
dt_Ti <- dt_sobol[dt_sobol$sensitivity == "Ti", ]
dt_Ti_wide <- reshape(as.data.frame(dt_Ti[, c("output", "parameters", "original")]),
                      idvar = "output", timevar = "parameters",
                      direction = "wide")
names(dt_Ti_wide) <- gsub("original\\.", "Ti_", names(dt_Ti_wide))
print(dt_Ti_wide, digits = 3)

# ---- 8. Comparison with N=512 results (if available) ------------------------

old_sobol_path <- "sensitivity_analysis/tables/sobol_indices.csv"
if (file.exists(old_sobol_path)) {
  cat("\n>> Building comparison table (N=512 vs N=2048)...\n")
  dt_old <- data.table::fread(old_sobol_path)
  dt_old$N <- 512
  dt_sobol$N <- 2048
  
  # Merge on output + parameters + sensitivity
  dt_compare <- merge(
    dt_old[, .(output, parameters, sensitivity, original_512 = original,
               se_512 = std_error, lo_512 = low_ci, hi_512 = high_ci)],
    dt_sobol[, .(output, parameters, sensitivity, original_2048 = original,
                 se_2048 = std_error, lo_2048 = low_ci, hi_2048 = high_ci)],
    by = c("output", "parameters", "sensitivity"),
    all = TRUE
  )
  dt_compare$delta <- dt_compare$original_2048 - dt_compare$original_512
  dt_compare$se_reduction_pct <- ifelse(
    !is.na(dt_compare$se_512) & dt_compare$se_512 > 0,
    (1 - dt_compare$se_2048 / dt_compare$se_512) * 100,
    NA_real_
  )
  
  write.csv(dt_compare, file.path(sa_tab, "sobol_comparison_N512_vs_N2048.csv"),
            row.names = FALSE)
  cat(">> Comparison table saved.\n")
  
  # Print Ti comparison
  cat("\n>> Ti comparison (change from N=512 to N=2048):\n")
  print(dt_compare[dt_compare$sensitivity == "Ti",
                   c("output", "parameters", "original_512", "original_2048",
                     "delta", "se_reduction_pct")],
        digits = 3)
}

cat(sprintf("\nSobol SA (refined) completed at %s\n", Sys.time()),
    file = sa_log, append = TRUE)
cat("\n>> 01_sobol_refined.R finished.\n")
