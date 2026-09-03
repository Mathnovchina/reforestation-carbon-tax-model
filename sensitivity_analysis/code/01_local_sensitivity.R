###############################################################################
# 01_local_sensitivity.R  –  Local (one-at-a-time) sensitivity checks
#
# Purpose:
#   - Vary phi1 and phi2 by +/-10%, +/-20%, +/-30%
#   - Verify numerical stability
#   - Produce line plots and summary tables
#   - Help determine plausible bounds for global SA
#
# Outputs:
#   sensitivity_analysis/tables/local_sensitivity_summary.csv
#   sensitivity_analysis/figures/local_OAT_*.png / .pdf
#   sensitivity_analysis/logs/local_sensitivity.log
###############################################################################

cat("============================================================\n")
cat(" 01  LOCAL SENSITIVITY ANALYSIS\n")
cat("============================================================\n")

# ---- 0. Setup ---------------------------------------------------------------

# Source helpers (assumes working directory is the project root)
source("sensitivity_analysis/code/00_helpers.R")

sa_fig  <- "sensitivity_analysis/figures"
sa_tab  <- "sensitivity_analysis/tables"
sa_log  <- "sensitivity_analysis/logs/local_sensitivity.log"

cat(sprintf("Local SA started at %s\n", Sys.time()), file = sa_log)

# ---- 1. Compile model -------------------------------------------------------

cat(">> Compiling model...\n")
sys <- compile_model()
cat(">> Model compiled.\n")

# ---- 2. Define perturbation grid --------------------------------------------

perturbations <- c(-0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30)
pct_labels    <- sprintf("%+.0f%%", perturbations * 100)

# phi1 varied, phi2 at baseline
grid_phi1 <- data.frame(
  phi1     = BASELINE_PHI1 * (1 + perturbations),
  phi2     = BASELINE_PHI2,
  varied   = "phi1",
  pct      = perturbations,
  pct_lab  = pct_labels,
  stringsAsFactors = FALSE
)

# phi2 varied, phi1 at baseline
grid_phi2 <- data.frame(
  phi1     = BASELINE_PHI1,
  phi2     = BASELINE_PHI2 * (1 + perturbations),
  varied   = "phi2",
  pct      = perturbations,
  pct_lab  = pct_labels,
  stringsAsFactors = FALSE
)

grid_local <- rbind(grid_phi1, grid_phi2)

# ---- 3. Run all OAT evaluations --------------------------------------------

cat(sprintf(">> Running %d local evaluations...\n", nrow(grid_local)))

local_results <- list()
local_ts      <- list()

for (i in seq_len(nrow(grid_local))) {
  g <- grid_local[i, ]
  out <- safe_run_model(sys, phi1 = g$phi1, phi2 = g$phi2,
                        return_ts = TRUE, log_file = sa_log)
  local_results[[i]] <- c(as.list(g), as.list(out$summary))
  if (!is.null(out$ts)) {
    out$ts$varied   <- g$varied
    out$ts$pct      <- g$pct
    out$ts$pct_lab  <- g$pct_lab
    out$ts$phi1_val <- g$phi1
    out$ts$phi2_val <- g$phi2
    local_ts[[i]]   <- out$ts
  }
}

dt_local  <- data.table::rbindlist(local_results, fill = TRUE)
dt_ts     <- data.table::rbindlist(local_ts, fill = TRUE)

cat(">> All local runs completed.\n")

# ---- 4. Summary table -------------------------------------------------------

write.csv(dt_local, file.path(sa_tab, "local_sensitivity_summary.csv"), row.names = FALSE)
cat(">> Summary table saved.\n")

# Print a compact view
print_cols <- c("varied", "pct_lab", "MC_2100", "tropical_stock_2100",
                "net_emissions_2100", "temperature_2100", "output_2100",
                "private_debt_ratio_2100")
print(dt_local[, ..print_cols])

# ---- 5. Check numerical stability ------------------------------------------

failed <- dt_local[!is.finite(MC_2100) | !is.finite(temperature_2100)]
if (nrow(failed) > 0) {
  cat("\n[WARNING] Some runs produced non-finite outputs:\n")
  print(failed[, .(varied, pct_lab, MC_2100, temperature_2100)])
} else {
  cat("\n>> All local runs numerically stable.\n")
}

# ---- 6. OAT line plots (each output vs perturbation) -----------------------

output_labels <- c(
  MC_2100                 = "C_F at 2100 (tril. $)",
  tropical_stock_2100     = "Tropical stock 2100 (bn m3)",
  sequestration_2100      = "Sequestration 2100 (GtCO2-e/yr)",
  net_emissions_2100      = "Net emissions 2100 (GtCO2-e/yr)",
  temperature_2100        = "Temperature 2100 (\u00B0C)",
  output_2100             = "Output 2100 (tril. $)",
  private_debt_ratio_2100 = "Private debt ratio 2100",
  cumulative_MC_2020_2100 = "Cumulative C_F 2020-2100 (tril. $\u00B7yr)"
)

dt_long <- dt_local %>%
  tidyr::pivot_longer(
    cols      = all_of(names(output_labels)),
    names_to  = "output",
    values_to = "value"
  ) %>%
  mutate(output_label = output_labels[output])

p_oat <- ggplot(dt_long, aes(x = pct * 100, y = value, colour = varied, shape = varied)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ output_label, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c(phi1 = "#D55E00", phi2 = "#0072B2"),
                      labels = c(phi1 = expression(varphi[1]~"(level)"),
                                 phi2 = expression(varphi[2]~"(curvature)"))) +
  scale_shape_manual(values = c(phi1 = 16, phi2 = 17),
                     labels = c(phi1 = expression(varphi[1]~"(level)"),
                                phi2 = expression(varphi[2]~"(curvature)"))) +
  labs(
    x      = "Perturbation from baseline (%)",
    y      = NULL,
    colour = "Parameter", shape = "Parameter"
  ) +
  theme_sa() +
  theme(strip.text = element_text(size = 9))

ggsave(file.path(sa_fig, "local_OAT_summary.png"), p_oat,
       width = 10, height = 12, dpi = 300)
ggsave(file.path(sa_fig, "local_OAT_summary.pdf"), p_oat,
       width = 10, height = 12)

cat(">> OAT summary plot saved.\n")

# ---- 7. Time-series ribbon plots -------------------------------------------

ts_vars <- c("MC", "tropical_stock", "net_emissions", "temperature",
             "output", "debt_ratio")
ts_labels <- c(
  MC              = "MC (tril. $)",
  tropical_stock  = "Tropical forest stock (bn m3)",
  net_emissions   = "Net emissions (GtCO2-e/yr)",
  temperature     = "Temperature anomaly (°C)",
  output          = "Output (tril. $)",
  debt_ratio      = "Private debt ratio"
)

for (vname in ts_vars) {
  dt_sub <- dt_ts[varied == "phi1", c("year", "pct_lab", vname), with = FALSE]
  setnames(dt_sub, vname, "value")
  
  p1 <- ggplot(dt_sub, aes(x = year, y = value, colour = pct_lab)) +
    geom_line(linewidth = 0.6) +
    scale_colour_viridis_d(option = "C", name = expression(phi[1]~"perturbation")) +
    labs(x = "Year", y = ts_labels[vname],
         title = bquote(.(ts_labels[vname]) ~ " — varying " ~ phi[1])) +
    theme_sa()
  
  dt_sub2 <- dt_ts[varied == "phi2", c("year", "pct_lab", vname), with = FALSE]
  setnames(dt_sub2, vname, "value")
  
  p2 <- ggplot(dt_sub2, aes(x = year, y = value, colour = pct_lab)) +
    geom_line(linewidth = 0.6) +
    scale_colour_viridis_d(option = "D", name = expression(phi[2]~"perturbation")) +
    labs(x = "Year", y = ts_labels[vname],
         title = bquote(.(ts_labels[vname]) ~ " — varying " ~ phi[2])) +
    theme_sa()
  
  p_combined <- p1 / p2 + plot_annotation(
    title = paste("Local sensitivity:", ts_labels[vname]),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )
  
  ggsave(file.path(sa_fig, paste0("local_ts_", vname, ".png")), p_combined,
         width = 10, height = 8, dpi = 300)
  ggsave(file.path(sa_fig, paste0("local_ts_", vname, ".pdf")), p_combined,
         width = 10, height = 8)
}

cat(">> Time-series plots saved.\n")

# ---- 8. Elasticity table ---------------------------------------------------

# Compute approximate elasticities at baseline (central difference, ±10%)
baseline_vals <- as.numeric(dt_local[varied == "phi1" & pct == 0,
                                      names(output_labels), with = FALSE])
names(baseline_vals) <- names(output_labels)

compute_elasticity <- function(param_name, pct_lo = -0.10, pct_hi = 0.10) {
  lo_row <- dt_local[varied == param_name & abs(pct - pct_lo) < 0.001]
  hi_row <- dt_local[varied == param_name & abs(pct - pct_hi) < 0.001]
  
  lo_vals <- as.numeric(lo_row[, names(output_labels), with = FALSE])
  hi_vals <- as.numeric(hi_row[, names(output_labels), with = FALSE])
  
  delta_pct <- pct_hi - pct_lo  # 0.20
  elasticity <- ((hi_vals - lo_vals) / baseline_vals) / delta_pct
  names(elasticity) <- names(output_labels)
  elasticity
}

elast_phi1 <- compute_elasticity("phi1")
elast_phi2 <- compute_elasticity("phi2")

dt_elast <- data.frame(
  output      = names(output_labels),
  description = as.character(output_labels),
  elast_phi1  = round(elast_phi1, 4),
  elast_phi2  = round(elast_phi2, 4),
  dominant    = ifelse(abs(elast_phi1) > abs(elast_phi2), "phi1", "phi2"),
  stringsAsFactors = FALSE
)

write.csv(dt_elast, file.path(sa_tab, "local_elasticities.csv"), row.names = FALSE)
cat("\n>> Elasticity table:\n")
print(dt_elast)

# ---- 9. Suggest bounds for global SA ---------------------------------------

cat("\n>> Suggested bounds for global Sobol SA (±30% of baseline):\n")
cat(sprintf("   phi1 in [%.5f, %.5f]\n",
            BASELINE_PHI1 * 0.70, BASELINE_PHI1 * 1.30))
cat(sprintf("   phi2 in [%.5f, %.5f]\n",
            BASELINE_PHI2 * 0.70, BASELINE_PHI2 * 1.30))

# Check if ±30% runs are stable
stable_30 <- all(is.finite(dt_local[abs(pct) == 0.30]$MC_2100))
if (stable_30) {
  cat("   All ±30% runs are numerically stable. Bounds are suitable.\n")
} else {
  cat("   [WARNING] Some ±30% runs are unstable. Consider narrower bounds.\n")
}

cat(sprintf("\nLocal SA completed at %s\n", Sys.time()), file = sa_log, append = TRUE)
cat("\n>> 01_local_sensitivity.R finished.\n")
