###############################################################################
# 02_convergence_diagnostics.R  –  Sobol index convergence vs sample size
#
# Runs Sobol SA at N = 128, 256, 512, 1024, 2048 and tracks how Si and Ti
# stabilise. Saves convergence plots and table.
#
# Outputs:
#   figures_final/convergence_Ti.{png,pdf}
#   figures_final/convergence_Si.{png,pdf}
#   tables_final/convergence_table.csv
#   logs_refined/convergence.log
###############################################################################

cat("============================================================\n")
cat(" CONVERGENCE DIAGNOSTICS\n")
cat("============================================================\n")

source("sensitivity_analysis/code_refined/00_helpers_refined.R")

sa_fig <- "sensitivity_analysis/figures_final"
sa_tab <- "sensitivity_analysis/tables_final"
sa_log <- "sensitivity_analysis/logs_refined/convergence.log"

cat(sprintf("Convergence diagnostics started at %s\n", Sys.time()), file = sa_log)

# ---- 1. Compile model -------------------------------------------------------

sys <- compile_model()

# ---- 2. Setup ----------------------------------------------------------------

phi1_lo <- BASELINE_PHI1 * 0.70;  phi1_hi <- BASELINE_PHI1 * 1.30
phi2_lo <- BASELINE_PHI2 * 0.70;  phi2_hi <- BASELINE_PHI2 * 1.30
params_sa <- c("phi1", "phi2")

# Key outputs to track convergence for
key_outputs <- c("MC_2100", "temperature_2100", "output_2100",
                 "private_debt_ratio_2100", "inflation_2100")

N_values <- c(128, 256, 512, 1024, 2048)

# ---- 3. Run at each N -------------------------------------------------------

convergence_rows <- list()

for (N_val in N_values) {
  cat(sprintf("\n>> N = %d  (%d evaluations)...\n", N_val, N_val * 5))
  
  set.seed(42)  # same seed for nested QRN sequences
  mat <- sobol_matrices(
    matrices = c("A", "B", "AB"),
    N        = N_val,
    params   = params_sa,
    order    = "second",
    type     = "QRN"
  )
  
  mat_phys <- mat
  mat_phys[, "phi1"] <- phi1_lo + mat[, "phi1"] * (phi1_hi - phi1_lo)
  mat_phys[, "phi2"] <- phi2_lo + mat[, "phi2"] * (phi2_hi - phi2_lo)
  
  param_df <- data.frame(phi1 = mat_phys[, "phi1"], phi2 = mat_phys[, "phi2"])
  
  t0 <- Sys.time()
  res <- run_batch(sys, param_df, log_file = sa_log, progress = FALSE)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("   Completed in %.1f s\n", elapsed))
  
  for (oname in key_outputs) {
    if (!oname %in% names(res)) next
    y_vec   <- res[[oname]]
    y_clean <- y_vec
    y_clean[!is.finite(y_clean)] <- median(y_vec, na.rm = TRUE)
    
    if (sum(is.finite(y_vec)) < 0.8 * length(y_vec)) next
    
    sob <- sobol_indices(
      matrices = c("A", "B", "AB"),
      Y        = y_clean,
      N        = N_val,
      params   = params_sa,
      boot     = TRUE,
      R        = 500,
      order    = "second",
      type     = "QRN"
    )
    
    idx <- sob$results
    z <- qnorm(0.975)
    for (j in seq_len(nrow(idx))) {
      lo <- idx$`low.ci`[j]
      hi <- idx$`high.ci`[j]
      if (is.na(lo) && !is.na(idx$std.error[j])) {
        lo <- idx$original[j] - z * idx$std.error[j]
        hi <- idx$original[j] + z * idx$std.error[j]
      }
      convergence_rows[[length(convergence_rows) + 1]] <- data.frame(
        N          = N_val,
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
}

dt_conv <- data.table::rbindlist(convergence_rows)
write.csv(dt_conv, file.path(sa_tab, "convergence_table.csv"), row.names = FALSE)
cat("\n>> Convergence table saved.\n")

# ---- 4. Convergence plots ---------------------------------------------------

nice_labels <- c(
  MC_2100                 = "MC at 2100",
  temperature_2100        = "Temperature 2100",
  output_2100             = "Output 2100",
  private_debt_ratio_2100 = "Debt ratio 2100",
  inflation_2100          = "Inflation 2100"
)

# --- Ti convergence ---
dt_Ti <- dt_conv[dt_conv$sensitivity == "Ti", ]
dt_Ti$output_label <- nice_labels[dt_Ti$output]
dt_Ti$param_label  <- ifelse(dt_Ti$parameters == "phi1",
                              "phi[1]~(level)", "phi[2]~(curvature)")

p_conv_Ti <- ggplot(dt_Ti, aes(x = N, y = original, colour = parameters)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = low_ci, ymax = high_ci, fill = parameters),
              alpha = 0.15, colour = NA) +
  facet_wrap(~ output_label, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = N_values, trans = "log2",
                     labels = N_values) +
  scale_colour_manual(
    values = c(phi1 = "#D55E00", phi2 = "#0072B2"),
    labels = c(phi1 = expression(varphi[1]), phi2 = expression(varphi[2]))
  ) +
  scale_fill_manual(
    values = c(phi1 = "#D55E00", phi2 = "#0072B2"),
    labels = c(phi1 = expression(varphi[1]), phi2 = expression(varphi[2]))
  ) +
  labs(
    x      = "Sobol sample size N",
    y      = "Total-order index (Ti)",
    colour = "Parameter", fill = "Parameter"
  ) +
  theme_sa(base_size = 10) +
  theme(panel.grid.minor.x = element_blank())

ggsave(file.path(sa_fig, "convergence_Ti.png"), p_conv_Ti,
       width = 12, height = 7, dpi = 300)
ggsave(file.path(sa_fig, "convergence_Ti.pdf"), p_conv_Ti,
       width = 12, height = 7)
cat(">> Ti convergence plot saved.\n")

# --- Si convergence ---
dt_Si <- dt_conv[dt_conv$sensitivity == "Si", ]
dt_Si$output_label <- nice_labels[dt_Si$output]

p_conv_Si <- ggplot(dt_Si, aes(x = N, y = original, colour = parameters)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = low_ci, ymax = high_ci, fill = parameters),
              alpha = 0.15, colour = NA) +
  facet_wrap(~ output_label, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = N_values, trans = "log2",
                     labels = N_values) +
  scale_colour_manual(
    values = c(phi1 = "#D55E00", phi2 = "#0072B2"),
    labels = c(phi1 = expression(varphi[1]), phi2 = expression(varphi[2]))
  ) +
  scale_fill_manual(
    values = c(phi1 = "#D55E00", phi2 = "#0072B2"),
    labels = c(phi1 = expression(varphi[1]), phi2 = expression(varphi[2]))
  ) +
  labs(
    x      = "Sobol sample size N",
    y      = "First-order index (Si)",
    colour = "Parameter", fill = "Parameter"
  ) +
  theme_sa(base_size = 10) +
  theme(panel.grid.minor.x = element_blank())

ggsave(file.path(sa_fig, "convergence_Si.png"), p_conv_Si,
       width = 12, height = 7, dpi = 300)
ggsave(file.path(sa_fig, "convergence_Si.pdf"), p_conv_Si,
       width = 12, height = 7)
cat(">> Si convergence plot saved.\n")

cat(sprintf("Convergence diagnostics completed at %s\n", Sys.time()),
    file = sa_log, append = TRUE)
cat("\n>> 02_convergence_diagnostics.R finished.\n")
