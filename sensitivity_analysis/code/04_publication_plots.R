###############################################################################
# 04_publication_plots.R  –  Assemble all publication-quality figures
#
# This script loads saved outputs and produces or refines:
#   A. Literature MC overlay
#   B. Fan chart of MC(RE)
#   C. Sobol bar plots with CI
#   D. Response-surface heatmaps (composite)
#   E. Scenario trajectory ribbons
#   F. Combined summary figure for the paper
#
# All figures are saved in both PNG 300dpi and PDF.
###############################################################################

cat("============================================================\n")
cat(" 04  PUBLICATION PLOTS\n")
cat("============================================================\n")

source("sensitivity_analysis/code/00_helpers.R")

sa_fig <- "sensitivity_analysis/figures"
sa_out <- "sensitivity_analysis/outputs"
sa_tab <- "sensitivity_analysis/tables"

# ---- 1. Load saved results ---------------------------------------------------

sobol_results <- NULL
res_full      <- NULL

if (file.exists(file.path(sa_out, "sobol_raw.rds"))) {
  sobol_results <- readRDS(file.path(sa_out, "sobol_raw.rds"))
  cat(">> Loaded Sobol results.\n")
}
if (file.exists(file.path(sa_out, "sobol_batch_results.rds"))) {
  res_full <- readRDS(file.path(sa_out, "sobol_batch_results.rds"))
  cat(">> Loaded batch results.\n")
}

dt_sobol <- NULL
if (file.exists(file.path(sa_tab, "sobol_indices.csv"))) {
  dt_sobol <- data.table::fread(file.path(sa_tab, "sobol_indices.csv"))
}

dt_local <- NULL
if (file.exists(file.path(sa_tab, "local_sensitivity_summary.csv"))) {
  dt_local <- data.table::fread(file.path(sa_tab, "local_sensitivity_summary.csv"))
}

# ---- 2. Refined Sobol bar plot (Figure C) ------------------------------------

if (!is.null(dt_sobol)) {
  nice_labels <- c(
    MC_2100                 = "MC at 2100",
    tropical_stock_2100     = "Tropical forest 2100",
    sequestration_2100      = "Sequestration 2100",
    net_emissions_2100      = "Net emissions 2100",
    temperature_2100        = "Temperature 2100",
    output_2100             = "Output 2100",
    private_debt_ratio_2100 = "Debt ratio 2100",
    cumulative_MC_2020_2100 = "Cumul. MC 2020-2100",
    cumulative_emissions    = "Cumul. emissions"
  )
  
  # Separate first and total order into side-by-side panels
  dt_ft <- dt_sobol[sensitivity %in% c("Si", "Ti")]
  dt_ft[, output_label := factor(nice_labels[output],
                                  levels = rev(nice_labels))]
  dt_ft[, index_label := ifelse(sensitivity == "Si",
                                 "First-order (Si)", "Total-order (Ti)")]
  
  p_sobol_pub <- ggplot(dt_ft, aes(x = output_label, y = original, fill = parameters)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.55) +
    geom_errorbar(aes(ymin = pmax(low_ci, 0), ymax = high_ci),
                  position = position_dodge(width = 0.7), width = 0.18,
                  linewidth = 0.4) +
    facet_wrap(~ index_label) +
    coord_flip() +
    scale_fill_manual(
      values = c(phi1 = "#D55E00", phi2 = "#0072B2"),
      labels = c(phi1 = expression(phi[1]~"(level)"),
                 phi2 = expression(phi[2]~"(curvature)"))
    ) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    labs(
      x    = NULL,
      y    = "Sobol index",
      fill = "Parameter",
      title = expression("Sobol indices: sensitivity of model outputs to MC parameters " *
                            phi[1] * " and " * phi[2])
    ) +
    theme_sa(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.text = element_text(face = "bold", size = 11)
    )
  
  ggsave(file.path(sa_fig, "pub_sobol_indices.png"), p_sobol_pub,
         width = 11, height = 6, dpi = 300)
  ggsave(file.path(sa_fig, "pub_sobol_indices.pdf"), p_sobol_pub,
         width = 11, height = 6)
  cat(">> Publication Sobol bar plot saved.\n")
}

# ---- 3. Summary 2x2 composite figure ----------------------------------------

cat(">> Compiling model for composite figure...\n")
sys <- compile_model()

# Baseline run
out_base <- run_model_sa(sys, phi1 = BASELINE_PHI1, phi2 = BASELINE_PHI2, return_ts = TRUE)

# Panel A: MC(RE) static curve with uncertainty band
RE_grid <- seq(0, 4, length.out = 200)

# Simulate envelope from ±30% bounds
set.seed(55)
n_env <- 300
phi1_samp <- runif(n_env, BASELINE_PHI1 * 0.7, BASELINE_PHI1 * 1.3)
phi2_samp <- runif(n_env, BASELINE_PHI2 * 0.7, BASELINE_PHI2 * 1.3)

env_data <- data.table::rbindlist(lapply(seq_len(n_env), function(i) {
  data.frame(RE = RE_grid, MC = phi1_samp[i] * exp(phi2_samp[i] * RE_grid), id = i)
}))

env_stats <- env_data[, .(
  med  = median(MC),
  lo90 = quantile(MC, 0.05),
  hi90 = quantile(MC, 0.95)
), by = RE]

pA <- ggplot() +
  geom_ribbon(data = env_stats, aes(x = RE, ymin = lo90, ymax = hi90),
              fill = "#0072B2", alpha = 0.25) +
  geom_line(data = env_stats, aes(x = RE, y = med),
            colour = "#0072B2", linewidth = 0.7) +
  geom_line(data = data.frame(RE = RE_grid,
                               MC = BASELINE_PHI1 * exp(BASELINE_PHI2 * RE_grid)),
            aes(x = RE, y = MC), colour = "black", linewidth = 1, linetype = "dashed") +
  coord_cartesian(ylim = c(0, quantile(env_data$MC, 0.97))) +
  labs(x = "RE (GtCO2-e)", y = "MC (tril. $)", title = "(A) MC(RE) uncertainty fan") +
  theme_sa(base_size = 9)

# Panel B: Temperature trajectories
# Run a small sample for trajectories
set.seed(77)
n_traj <- 50
traj_ts <- list()
for (i in seq_len(n_traj)) {
  p1 <- runif(1, BASELINE_PHI1 * 0.7, BASELINE_PHI1 * 1.3)
  p2 <- runif(1, BASELINE_PHI2 * 0.7, BASELINE_PHI2 * 1.3)
  out <- safe_run_model(sys, phi1 = p1, phi2 = p2, return_ts = TRUE)
  if (!is.null(out$ts)) {
    out$ts$id <- i
    traj_ts[[i]] <- out$ts
  }
}
dt_traj <- data.table::rbindlist(traj_ts, fill = TRUE)

env_temp <- dt_traj[, .(
  med  = median(temperature, na.rm = TRUE),
  lo   = quantile(temperature, 0.05, na.rm = TRUE),
  hi   = quantile(temperature, 0.95, na.rm = TRUE)
), by = year]

pB <- ggplot(env_temp, aes(x = year)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#CC79A7", alpha = 0.25) +
  geom_line(aes(y = med), colour = "#CC79A7", linewidth = 0.8) +
  geom_line(data = out_base$ts, aes(x = year, y = temperature),
            colour = "black", linewidth = 0.8, linetype = "dashed") +
  labs(x = "Year", y = "°C anomaly", title = "(B) Temperature trajectory uncertainty") +
  theme_sa(base_size = 9)

# Panel C: Tropical forest stock trajectories
env_ftro <- dt_traj[, .(
  med  = median(tropical_stock, na.rm = TRUE),
  lo   = quantile(tropical_stock, 0.05, na.rm = TRUE),
  hi   = quantile(tropical_stock, 0.95, na.rm = TRUE)
), by = year]

pC <- ggplot(env_ftro, aes(x = year)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#009E73", alpha = 0.25) +
  geom_line(aes(y = med), colour = "#009E73", linewidth = 0.8) +
  geom_line(data = out_base$ts, aes(x = year, y = tropical_stock),
            colour = "black", linewidth = 0.8, linetype = "dashed") +
  labs(x = "Year", y = "Billion m3", title = "(C) Tropical forest stock uncertainty") +
  theme_sa(base_size = 9)

# Panel D: Output (GDP) trajectories
env_gdp <- dt_traj[, .(
  med  = median(output, na.rm = TRUE),
  lo   = quantile(output, 0.05, na.rm = TRUE),
  hi   = quantile(output, 0.95, na.rm = TRUE)
), by = year]

pD <- ggplot(env_gdp, aes(x = year)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#E69F00", alpha = 0.25) +
  geom_line(aes(y = med), colour = "#E69F00", linewidth = 0.8) +
  geom_line(data = out_base$ts, aes(x = year, y = output),
            colour = "black", linewidth = 0.8, linetype = "dashed") +
  labs(x = "Year", y = "Tril. 2015 US$", title = "(D) Output trajectory uncertainty") +
  theme_sa(base_size = 9)

# Compose
p_composite <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title    = expression("Sensitivity of model outcomes to avoided-deforestation MC parameters " *
                            phi[1] * " and " * phi[2]),
    subtitle = "Dashed black = baseline; shaded = 90% envelope under ±30% uniform uncertainty",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, colour = "grey30")
    )
  )

ggsave(file.path(sa_fig, "pub_composite_4panel.png"), p_composite,
       width = 12, height = 9, dpi = 300)
ggsave(file.path(sa_fig, "pub_composite_4panel.pdf"), p_composite,
       width = 12, height = 9)
cat(">> Composite 4-panel figure saved.\n")

# ---- 4. Refined heatmap (single page) ---------------------------------------

if (!is.null(res_full)) {
  res_h <- as.data.frame(res_full)
  
  make_heat <- function(oname, lab, pal = "C") {
    ggplot(res_h, aes(x = phi1, y = phi2, colour = .data[[oname]])) +
      geom_point(size = 0.8, alpha = 0.6) +
      scale_colour_viridis_c(option = pal, name = NULL) +
      labs(x = expression(phi[1]), y = expression(phi[2]), title = lab) +
      theme_sa(base_size = 8) +
      theme(legend.position = "right", legend.key.height = unit(0.4, "cm"))
  }
  
  p_h1 <- make_heat("temperature_2100", "Temperature 2100 (°C)", "B")
  p_h2 <- make_heat("tropical_stock_2100", "Tropical forest 2100 (bn m3)", "D")
  p_h3 <- make_heat("output_2100", "Output 2100 (tril. $)", "C")
  p_h4 <- make_heat("private_debt_ratio_2100", "Debt ratio 2100", "A")
  
  p_heat <- (p_h1 | p_h2) / (p_h3 | p_h4) +
    plot_annotation(
      title = "Response surfaces in MC parameter space",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
  
  ggsave(file.path(sa_fig, "pub_heatmaps_2x2.png"), p_heat,
         width = 11, height = 9, dpi = 300)
  ggsave(file.path(sa_fig, "pub_heatmaps_2x2.pdf"), p_heat,
         width = 11, height = 9)
  cat(">> Publication heatmaps saved.\n")
}

cat("\n>> 04_publication_plots.R finished.\n")
