###############################################################################
# 02_sobol_global_sa.R  –  Global Sobol sensitivity analysis
#
# Purpose:
#   - Sobol first-order (S1), total-order (ST) and second-order (S2)
#     indices for phi1 and phi2 on all key outputs
#   - Uses the sensobol package (Saltelli 2002 / Jansen estimators)
#   - Bootstrap confidence intervals
#   - Pilot run first (small N), then full run
#
# Outputs:
#   sensitivity_analysis/outputs/sobol_raw.rds
#   sensitivity_analysis/tables/sobol_indices.csv
#   sensitivity_analysis/figures/sobol_barplot_*.png / .pdf
#   sensitivity_analysis/figures/sobol_heatmap_*.png / .pdf
#   sensitivity_analysis/figures/sobol_fan_MC.png / .pdf
#   sensitivity_analysis/logs/sobol_sa.log
###############################################################################

cat("============================================================\n")
cat(" 02  GLOBAL SOBOL SENSITIVITY ANALYSIS\n")
cat("============================================================\n")

# ---- 0. Setup ---------------------------------------------------------------

source("sensitivity_analysis/code/00_helpers.R")

sa_fig  <- "sensitivity_analysis/figures"
sa_tab  <- "sensitivity_analysis/tables"
sa_out  <- "sensitivity_analysis/outputs"
sa_log  <- "sensitivity_analysis/logs/sobol_sa.log"

cat(sprintf("Sobol SA started at %s\n", Sys.time()), file = sa_log)

# ---- 1. Compile model -------------------------------------------------------

cat(">> Compiling model...\n")
sys <- compile_model()
cat(">> Model compiled.\n")

# ---- 2. Define parameter bounds --------------------------------------------

# Uniform distributions: [0.7*baseline, 1.3*baseline]
# Modify these bounds as needed
phi1_lo <- BASELINE_PHI1 * 0.70
phi1_hi <- BASELINE_PHI1 * 1.30
phi2_lo <- BASELINE_PHI2 * 0.70
phi2_hi <- BASELINE_PHI2 * 1.30

cat(sprintf(">> phi1 bounds: [%.5f, %.5f]\n", phi1_lo, phi1_hi))
cat(sprintf(">> phi2 bounds: [%.5f, %.5f]\n", phi2_lo, phi2_hi))

params_sa <- c("phi1", "phi2")

# ---- 3. Pilot run -----------------------------------------------------------

N_pilot <- 100  # gives 100*(2+2) = 400 model evaluations

cat(sprintf("\n>> Pilot run: N = %d (%d total evaluations)...\n",
            N_pilot, N_pilot * (length(params_sa) + 2)))

# Generate Sobol quasi-random matrices using sensobol
set.seed(2024)
mat_pilot <- sobol_matrices(
  matrices = c("A", "B", "AB"),
  N        = N_pilot,
  params   = params_sa,
  order    = "second",
  type     = "QRN"    # quasi-random (Sobol sequences)
)

# Scale from [0,1] to physical bounds
mat_pilot_phys <- mat_pilot
mat_pilot_phys[, "phi1"] <- phi1_lo + mat_pilot[, "phi1"] * (phi1_hi - phi1_lo)
mat_pilot_phys[, "phi2"] <- phi2_lo + mat_pilot[, "phi2"] * (phi2_hi - phi2_lo)

param_df_pilot <- data.frame(
  phi1 = mat_pilot_phys[, "phi1"],
  phi2 = mat_pilot_phys[, "phi2"]
)

# Run batch
t0 <- Sys.time()
res_pilot <- run_batch(sys, param_df_pilot, log_file = sa_log, progress = TRUE)
t1 <- Sys.time()
pilot_time <- as.numeric(difftime(t1, t0, units = "secs"))

cat(sprintf(">> Pilot run completed in %.1f seconds (%.3f s/run).\n",
            pilot_time, pilot_time / nrow(param_df_pilot)))

# Check for failures
n_fail <- sum(is.na(res_pilot$MC_2100))
cat(sprintf(">> Failures: %d / %d (%.1f%%)\n",
            n_fail, nrow(res_pilot), 100 * n_fail / nrow(res_pilot)))

# Quick Sobol estimate on pilot data
output_names <- c(
  "MC_2100", "tropical_stock_2100", "sequestration_2100",
  "net_emissions_2100", "temperature_2100", "output_2100",
  "private_debt_ratio_2100", "cumulative_MC_2020_2100", "cumulative_emissions"
)

# ---- 4. Full Sobol run ------------------------------------------------------

# Estimate N for stable indices: target total evals ~ 2000-5000
# With 2 params: total evals = N * (2 + 2) = 4*N
N_full <- 512   # 512 * 4 = 2048 evaluations; increase to 1024 or 2048 if time permits

# Estimate runtime
est_seconds <- (N_full * 4) * (pilot_time / nrow(param_df_pilot))
cat(sprintf("\n>> Full run: N = %d (%d evaluations), est. %.0f seconds\n",
            N_full, N_full * 4, est_seconds))

set.seed(42)
mat_full <- sobol_matrices(
  matrices = c("A", "B", "AB"),
  N        = N_full,
  params   = params_sa,
  order    = "second",
  type     = "QRN"
)

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
cat(sprintf(">> Full run completed in %.1f seconds.\n",
            as.numeric(difftime(t1, t0, units = "secs"))))

n_fail_full <- sum(is.na(res_full$MC_2100))
cat(sprintf(">> Failures: %d / %d\n", n_fail_full, nrow(res_full)))

# ---- 5. Compute Sobol indices -----------------------------------------------

cat("\n>> Computing Sobol indices...\n")

sobol_results <- list()

for (oname in output_names) {
  
  y_vec <- res_full[[oname]]
  
  # Skip if too many NAs
  if (sum(is.finite(y_vec)) < 0.8 * length(y_vec)) {
    cat(sprintf("   [SKIP] %s: too many NAs (%.0f%%)\n",
                oname, 100 * sum(!is.finite(y_vec)) / length(y_vec)))
    next
  }
  
  # Replace NAs with median (sensobol does not handle NAs)
  y_clean <- y_vec
  y_clean[!is.finite(y_clean)] <- median(y_vec, na.rm = TRUE)
  
  # Compute indices using sensobol
  sob <- sobol_indices(
    matrices = c("A", "B", "AB"),
    Y        = y_clean,
    N        = N_full,
    params   = params_sa,
    boot     = TRUE,
    R        = 500,          # bootstrap replicates
    order    = "second",
    type     = "QRN"
  )
  
  sobol_results[[oname]] <- sob
}

# ---- 6. Save raw results ----------------------------------------------------

saveRDS(sobol_results, file.path(sa_out, "sobol_raw.rds"))
saveRDS(res_full,      file.path(sa_out, "sobol_batch_results.rds"))
cat(">> Raw Sobol objects saved.\n")

# ---- 7. Tidy summary table --------------------------------------------------

tidy_sobol <- function(sobol_results, conf_level = 0.95) {
  z <- qnorm(1 - (1 - conf_level) / 2)
  rows <- list()
  for (oname in names(sobol_results)) {
    sob <- sobol_results[[oname]]
    idx <- sob$results
    
    for (j in seq_len(nrow(idx))) {
      # Use explicit CI columns if available, otherwise compute from std.error
      lo <- idx$`low.ci`[j]
      hi <- idx$`high.ci`[j]
      if (is.na(lo) && !is.na(idx$std.error[j])) {
        lo <- idx$original[j] - z * idx$std.error[j]
        hi <- idx$original[j] + z * idx$std.error[j]
      }
      rows[[length(rows) + 1]] <- data.frame(
        output     = oname,
        parameters = idx$parameters[j],
        sensitivity= idx$sensitivity[j],
        original   = idx$original[j],
        std_error  = idx$std.error[j],
        low_ci     = lo,
        high_ci    = hi,
        stringsAsFactors = FALSE
      )
    }
  }
  data.table::rbindlist(rows)
}

dt_sobol <- tidy_sobol(sobol_results)
write.csv(dt_sobol, file.path(sa_tab, "sobol_indices.csv"), row.names = FALSE)
cat(">> Sobol indices table saved.\n")
print(dt_sobol)

# ---- 8. Sobol bar plots with CI --------------------------------------------

nice_labels <- c(
  MC_2100                 = "MC at 2100",
  tropical_stock_2100     = "Trop. forest 2100",
  sequestration_2100      = "Sequestration 2100",
  net_emissions_2100      = "Net emissions 2100",
  temperature_2100        = "Temperature 2100",
  output_2100             = "Output 2100",
  private_debt_ratio_2100 = "Debt ratio 2100",
  cumulative_MC_2020_2100 = "Cumul. MC 2020-2100",
  cumulative_emissions    = "Cumul. emissions"
)

# First and total order bar plot
dt_ft <- dt_sobol[sensitivity %in% c("Si", "Ti")]
dt_ft[, output_label := nice_labels[output]]
dt_ft[, param_label  := ifelse(parameters == "phi1",
                                expression(phi[1]),
                                expression(phi[2]))]

p_sobol_bar <- ggplot(dt_ft, aes(x = output_label, y = original,
                                  fill = parameters, alpha = sensitivity)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = low_ci, ymax = high_ci),
                position = position_dodge(width = 0.7), width = 0.2) +
  scale_fill_manual(
    values = c(phi1 = "#D55E00", phi2 = "#0072B2"),
    labels = c(phi1 = expression(phi[1]~"(level)"),
               phi2 = expression(phi[2]~"(curvature)"))
  ) +
  scale_alpha_manual(
    values = c(Si = 1.0, Ti = 0.5),
    labels = c(Si = "First-order (Si)", Ti = "Total-order (Ti)")
  ) +
  coord_flip() +
  labs(
    x     = NULL,
    y     = "Sobol index",
    fill  = "Parameter",
    alpha = "Index type",
    title = "Sobol sensitivity indices for MC parameters",
    subtitle = sprintf("N = %d, bootstrap R = 500, uniform bounds ±30%%", N_full)
  ) +
  theme_sa() +
  theme(legend.box = "vertical")

ggsave(file.path(sa_fig, "sobol_barplot_first_total.png"), p_sobol_bar,
       width = 10, height = 7, dpi = 300)
ggsave(file.path(sa_fig, "sobol_barplot_first_total.pdf"), p_sobol_bar,
       width = 10, height = 7)
cat(">> Sobol bar plot saved.\n")

# ---- 9. Second-order interaction plot ---------------------------------------

dt_s2 <- dt_sobol[sensitivity == "Sij"]
if (nrow(dt_s2) > 0) {
  dt_s2[, output_label := nice_labels[output]]
  
  p_s2 <- ggplot(dt_s2, aes(x = output_label, y = original)) +
    geom_col(fill = "#009E73", width = 0.5) +
    geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = 0.2) +
    coord_flip() +
    labs(
      x     = NULL,
      y     = expression("Second-order index  " * S[ij] * "  (" * phi[1] * " × " * phi[2] * ")"),
      title = expression("Interaction between " * phi[1] * " and " * phi[2])
    ) +
    theme_sa()
  
  ggsave(file.path(sa_fig, "sobol_interaction_s2.png"), p_s2,
         width = 8, height = 6, dpi = 300)
  ggsave(file.path(sa_fig, "sobol_interaction_s2.pdf"), p_s2,
         width = 8, height = 6)
  cat(">> Interaction plot saved.\n")
}

# ---- 10. Fan chart of MC(RE) under sampled phi1, phi2 ----------------------

cat(">> Generating MC fan chart...\n")

# Sample a subset of parameter pairs for the fan chart
set.seed(99)
n_fan   <- min(200, nrow(param_df_full))
fan_idx <- sample(nrow(param_df_full), n_fan)

RE_grid <- seq(0, 4, length.out = 200)   # RE in GtCO2-e

fan_data <- data.table::rbindlist(lapply(fan_idx, function(i) {
  p1 <- param_df_full$phi1[i]
  p2 <- param_df_full$phi2[i]
  data.frame(
    RE = RE_grid,
    MC = p1 * exp(p2 * RE_grid),
    sample_id = i
  )
}))

# Percentile ribbons
fan_ribbons <- fan_data[, .(
  med  = median(MC),
  lo90 = quantile(MC, 0.05),
  hi90 = quantile(MC, 0.95),
  lo50 = quantile(MC, 0.25),
  hi50 = quantile(MC, 0.75)
), by = RE]

# Baseline curve
baseline_mc <- data.frame(
  RE = RE_grid,
  MC = BASELINE_PHI1 * exp(BASELINE_PHI2 * RE_grid)
)

p_fan <- ggplot() +
  geom_ribbon(data = fan_ribbons, aes(x = RE, ymin = lo90, ymax = hi90),
              fill = "#0072B2", alpha = 0.2) +
  geom_ribbon(data = fan_ribbons, aes(x = RE, ymin = lo50, ymax = hi50),
              fill = "#0072B2", alpha = 0.3) +
  geom_line(data = fan_ribbons, aes(x = RE, y = med),
            colour = "#0072B2", linewidth = 0.8) +
  geom_line(data = baseline_mc, aes(x = RE, y = MC),
            colour = "black", linewidth = 1, linetype = "dashed") +
  labs(
    x        = "Avoided deforestation RE (GtCO2-e)",
    y        = "Marginal Cost MC (tril. 2015 US$)",
    title    = expression("Fan chart of MC(RE) = " * phi[1] * " exp(" * phi[2] * " RE)"),
    subtitle = "Shaded: 50% and 90% envelopes under uniform ±30% uncertainty; dashed = baseline"
  ) +
  coord_cartesian(ylim = c(0, quantile(fan_data$MC, 0.98))) +
  theme_sa()

ggsave(file.path(sa_fig, "sobol_fan_MC_RE.png"), p_fan,
       width = 8, height = 6, dpi = 300)
ggsave(file.path(sa_fig, "sobol_fan_MC_RE.pdf"), p_fan,
       width = 8, height = 6)
cat(">> Fan chart saved.\n")

# ---- 11. Response-surface heatmaps -----------------------------------------

cat(">> Generating response-surface heatmaps...\n")

heatmap_outputs <- c("output_2100", "private_debt_ratio_2100",
                     "temperature_2100", "tropical_stock_2100",
                     "MC_2100", "net_emissions_2100")

heatmap_labels <- c(
  output_2100             = "Output 2100 (tril. $)",
  private_debt_ratio_2100 = "Private debt ratio 2100",
  temperature_2100        = "Temperature 2100 (°C)",
  tropical_stock_2100     = "Tropical forest 2100 (bn m3)",
  MC_2100                 = "MC 2100 (tril. $)",
  net_emissions_2100      = "Net emissions 2100 (GtCO2-e)"
)

# Use the full batch results — they already have phi1, phi2, and all outputs
res_heat <- as.data.frame(res_full)

heat_plots <- list()
for (oname in heatmap_outputs) {
  p <- ggplot(res_heat, aes(x = phi1, y = phi2, colour = .data[[oname]])) +
    geom_point(size = 1.2, alpha = 0.7) +
    scale_colour_viridis_c(option = "C", name = heatmap_labels[oname]) +
    labs(
      x     = expression(phi[1]~"(level)"),
      y     = expression(phi[2]~"(curvature)"),
      title = heatmap_labels[oname]
    ) +
    theme_sa() +
    theme(legend.position = "right")
  heat_plots[[oname]] <- p
}

p_heat_all <- wrap_plots(heat_plots, ncol = 2) +
  plot_annotation(
    title = "Response surfaces in MC parameter space",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(file.path(sa_fig, "sobol_heatmaps.png"), p_heat_all,
       width = 14, height = 16, dpi = 300)
ggsave(file.path(sa_fig, "sobol_heatmaps.pdf"), p_heat_all,
       width = 14, height = 16)
cat(">> Heatmaps saved.\n")

# ---- 12. Trajectory ribbons from Sobol sample --------------------------------

cat(">> Generating trajectory ribbons...\n")

# Run a subset with time series
set.seed(123)
n_ribbon  <- min(100, nrow(param_df_full))
rib_idx   <- sample(nrow(param_df_full), n_ribbon)
ribbon_ts <- list()

for (i in seq_along(rib_idx)) {
  j <- rib_idx[i]
  out <- safe_run_model(sys,
                        phi1 = param_df_full$phi1[j],
                        phi2 = param_df_full$phi2[j],
                        return_ts = TRUE, log_file = sa_log)
  if (!is.null(out$ts)) {
    out$ts$sample_id <- i
    ribbon_ts[[i]] <- out$ts
  }
}
dt_ribbon <- data.table::rbindlist(ribbon_ts, fill = TRUE)

# Compute percentile envelopes by year
ribbon_vars <- c("MC", "tropical_stock", "net_emissions", "temperature",
                 "output", "debt_ratio")
ribbon_labels <- c(
  MC             = "MC (tril. $)",
  tropical_stock = "Tropical forest stock (bn m3)",
  net_emissions  = "Net emissions (GtCO2-e/yr)",
  temperature    = "Temperature anomaly (°C)",
  output         = "Output (tril. $)",
  debt_ratio     = "Private debt ratio"
)

ribbon_plots <- list()
for (vname in ribbon_vars) {
  env <- dt_ribbon[, .(
    med  = median(get(vname), na.rm = TRUE),
    lo90 = quantile(get(vname), 0.05, na.rm = TRUE),
    hi90 = quantile(get(vname), 0.95, na.rm = TRUE),
    lo50 = quantile(get(vname), 0.25, na.rm = TRUE),
    hi50 = quantile(get(vname), 0.75, na.rm = TRUE)
  ), by = year]
  
  p <- ggplot(env, aes(x = year)) +
    geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "#0072B2", alpha = 0.15) +
    geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = "#0072B2", alpha = 0.30) +
    geom_line(aes(y = med), colour = "#0072B2", linewidth = 0.8) +
    labs(x = "Year", y = ribbon_labels[vname], title = ribbon_labels[vname]) +
    theme_sa()
  
  ribbon_plots[[vname]] <- p
}

p_ribbons <- wrap_plots(ribbon_plots, ncol = 2) +
  plot_annotation(
    title    = expression("Trajectory uncertainty from MC parameters (" * phi[1] * ", " * phi[2] * ")"),
    subtitle = "Shaded: 50% and 90% envelopes; line = median",
    theme    = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(file.path(sa_fig, "sobol_trajectory_ribbons.png"), p_ribbons,
       width = 12, height = 10, dpi = 300)
ggsave(file.path(sa_fig, "sobol_trajectory_ribbons.pdf"), p_ribbons,
       width = 12, height = 10)
cat(">> Trajectory ribbons saved.\n")

cat(sprintf("\nSobol SA completed at %s\n", Sys.time()), file = sa_log, append = TRUE)
cat("\n>> 02_sobol_global_sa.R finished.\n")
