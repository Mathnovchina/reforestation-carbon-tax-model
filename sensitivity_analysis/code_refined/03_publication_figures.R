###############################################################################
# 03_publication_figures.R  –  Final publication-quality figures
#
# Reads refined Sobol outputs and generates:
#   Figure A (main text): Literature MC overlay  [reuses existing]
#   Figure B (main text): Sobol bar chart (Ti prominent, Si lighter)
#   Figure C (main text): 4-panel trajectory bands
#                         (MC, tropical forest, inflation, debt ratio)
#   Figure D (appendix):  Convergence diagnostics  [already produced by 02]
#   Figure E (appendix):  Response-surface heatmaps 2x2
#   Figure F (appendix):  6-panel trajectory ribbons (all variables)
#   Figure G (appendix):  OAT spider plot  [reuses existing]
###############################################################################

cat("============================================================\n")
cat(" PUBLICATION FIGURES (refined)\n")
cat("============================================================\n")

source("sensitivity_analysis/code_refined/00_helpers_refined.R")

sa_fig <- "sensitivity_analysis/figures_final"
sa_tab <- "sensitivity_analysis/tables_final"
sa_out <- "sensitivity_analysis/outputs_refined"

# ---- 0. Load data ------------------------------------------------------------

dt_sobol <- NULL
if (file.exists(file.path(sa_tab, "sobol_indices_N2048.csv"))) {
  dt_sobol <- data.table::fread(file.path(sa_tab, "sobol_indices_N2048.csv"))
  cat(">> Loaded N=2048 Sobol indices.\n")
}

res_full <- NULL
if (file.exists(file.path(sa_out, "sobol_batch_N2048.rds"))) {
  res_full <- readRDS(file.path(sa_out, "sobol_batch_N2048.rds"))
  cat(">> Loaded batch results.\n")
}

# ---- 1. FIGURE A  –  Literature MC overlay (copy existing) ------------------

# The existing literature overlay in sensitivity_analysis/figures/ is already
# publication quality. Copy it to figures_final/ for completeness.
for (ext in c("png", "pdf")) {
  src <- file.path("sensitivity_analysis/figures", paste0("literature_cost_curves_overlay.", ext))
  dst <- file.path(sa_fig, paste0("FigA_literature_overlay.", ext))
  if (file.exists(src)) file.copy(src, dst, overwrite = TRUE)
}
cat(">> Figure A: Literature overlay copied.\n")

# ---- 2. FIGURE B  –  Sobol bar chart ----------------------------------------

if (!is.null(dt_sobol)) {
  
  nice_labels <- c(
    MC_2100                 = "C_F at 2100",
    tropical_stock_2100     = "Tropical forest 2100",
    sequestration_2100      = "Sequestration 2100",
    net_emissions_2100      = "Net emissions 2100",
    temperature_2100        = "Temperature 2100",
    output_2100             = "Output (GDP) 2100",
    private_debt_ratio_2100 = "Debt ratio 2100",
    cumulative_MC_2020_2100 = "Cumul. C_F 2020-2100",
    cumulative_emissions    = "Cumul. emissions",
    inflation_2100          = "Inflation 2100",
    avg_inflation_2050_2100 = "Avg. inflation 2050-2100",
    min_inflation_post2050  = "Min. inflation post-2050"
  )
  
  # Filter to main outputs for the paper figure
  main_outputs <- c("MC_2100", "cumulative_MC_2020_2100", "tropical_stock_2100",
                    "net_emissions_2100", "temperature_2100", "output_2100",
                    "private_debt_ratio_2100", "inflation_2100")
  
  dt_ft <- dt_sobol[sensitivity %in% c("Si", "Ti") & output %in% main_outputs]
  dt_ft[, output_label := factor(nice_labels[output],
                                  levels = rev(nice_labels[main_outputs]))]
  dt_ft[, index_label := ifelse(sensitivity == "Si",
                                 "First-order (Si)", "Total-order (Ti)")]
  
  p_FigB <- ggplot(dt_ft, aes(x = output_label, y = original, fill = parameters)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.55) +
    geom_errorbar(aes(ymin = pmax(low_ci, -0.05), ymax = pmin(high_ci, 1.05)),
                  position = position_dodge(width = 0.7), width = 0.18,
                  linewidth = 0.4) +
    facet_wrap(~ index_label) +
    coord_flip() +
    scale_fill_manual(
      values = c(phi1 = "#D55E00", phi2 = "#0072B2"),
      labels = c(phi1 = expression(varphi[1]~"(level)"),
                 phi2 = expression(varphi[2]~"(curvature)"))
    ) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(-0.05, 1.05)) +
    labs(
      x    = NULL,
      y    = "Sobol index",
      fill = "Parameter"
    ) +
    theme_sa(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.text = element_text(face = "bold", size = 11)
    )
  
  ggsave(file.path(sa_fig, "FigB_sobol_indices.png"), p_FigB,
         width = 11, height = 6.5, dpi = 300)
  ggsave(file.path(sa_fig, "FigB_sobol_indices.pdf"), p_FigB,
         width = 11, height = 6.5)
  cat(">> Figure B: Sobol bar chart saved.\n")
}

# ---- 3. FIGURE C  –  4-panel trajectory bands -------------------------------

cat(">> Compiling model for trajectory bands...\n")
sys <- compile_model()

# Baseline
out_base <- run_model_sa(sys, phi1 = BASELINE_PHI1, phi2 = BASELINE_PHI2,
                         return_ts = TRUE)

# Draw trajectories from the Sobol parameter sample (or fresh sample)
set.seed(77)
n_traj <- 200

if (!is.null(res_full)) {
  # Use a random subset from the Sobol sample
  traj_idx <- sample(nrow(res_full), min(n_traj, nrow(res_full)))
  phi1_samp <- res_full$phi1[traj_idx]
  phi2_samp <- res_full$phi2[traj_idx]
} else {
  phi1_samp <- runif(n_traj, BASELINE_PHI1 * 0.7, BASELINE_PHI1 * 1.3)
  phi2_samp <- runif(n_traj, BASELINE_PHI2 * 0.7, BASELINE_PHI2 * 1.3)
}

cat(sprintf(">> Generating %d trajectory runs...\n", length(phi1_samp)))
traj_ts <- list()
for (i in seq_along(phi1_samp)) {
  out <- safe_run_model(sys, phi1 = phi1_samp[i], phi2 = phi2_samp[i],
                        return_ts = TRUE)
  if (!is.null(out$ts)) {
    out$ts$id <- i
    traj_ts[[i]] <- out$ts
  }
}
dt_traj <- data.table::rbindlist(traj_ts, fill = TRUE)

# Envelope helper
make_envelope <- function(dt, varname) {
  dt[, .(
    med  = median(get(varname), na.rm = TRUE),
    lo90 = quantile(get(varname), 0.05, na.rm = TRUE),
    hi90 = quantile(get(varname), 0.95, na.rm = TRUE),
    lo50 = quantile(get(varname), 0.25, na.rm = TRUE),
    hi50 = quantile(get(varname), 0.75, na.rm = TRUE)
  ), by = year]
}

# Panel A: MC
env_MC <- make_envelope(dt_traj, "MC")
pA <- ggplot(env_MC, aes(x = year)) +
  geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "#0072B2", alpha = 0.15) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = "#0072B2", alpha = 0.30) +
  geom_line(aes(y = med), colour = "#0072B2", linewidth = 0.8) +
  geom_line(data = out_base$ts, aes(x = year, y = MC),
            colour = "black", linewidth = 0.8, linetype = "dashed") +
  labs(x = NULL, y = "Tril. 2015 US$",
       title = "(A) Forest-policy cost") +
  theme_sa(base_size = 9)

# Panel B: Tropical forest stock
env_F <- make_envelope(dt_traj, "tropical_stock")
pB <- ggplot(env_F, aes(x = year)) +
  geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "#009E73", alpha = 0.15) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = "#009E73", alpha = 0.30) +
  geom_line(aes(y = med), colour = "#009E73", linewidth = 0.8) +
  geom_line(data = out_base$ts, aes(x = year, y = tropical_stock),
            colour = "black", linewidth = 0.8, linetype = "dashed") +
  labs(x = NULL, y = "Billion m\u00B3",
       title = "(B) Tropical forest stock") +
  theme_sa(base_size = 9)

# Panel C: Inflation
env_I <- make_envelope(dt_traj, "inflation_pct")
pC <- ggplot(env_I, aes(x = year)) +
  geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "#CC79A7", alpha = 0.15) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = "#CC79A7", alpha = 0.30) +
  geom_line(aes(y = med), colour = "#CC79A7", linewidth = 0.8) +
  geom_line(data = out_base$ts, aes(x = year, y = inflation_pct),
            colour = "black", linewidth = 0.8, linetype = "dashed") +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dotted") +
  labs(x = "Year", y = "Inflation rate (%)",
       title = "(C) Inflation rate") +
  theme_sa(base_size = 9)

# Panel D: Private debt ratio
env_D <- make_envelope(dt_traj, "debt_ratio")
pD <- ggplot(env_D, aes(x = year)) +
  geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "#E69F00", alpha = 0.15) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = "#E69F00", alpha = 0.30) +
  geom_line(aes(y = med), colour = "#E69F00", linewidth = 0.8) +
  geom_line(data = out_base$ts, aes(x = year, y = debt_ratio),
            colour = "black", linewidth = 0.8, linetype = "dashed") +
  labs(x = "Year", y = "Debt / nominal GDP",
       title = "(D) Private debt ratio") +
  theme_sa(base_size = 9)

p_FigC <- (pA | pB) / (pC | pD)

ggsave(file.path(sa_fig, "FigC_4panel_trajectories.png"), p_FigC,
       width = 12, height = 9, dpi = 300)
ggsave(file.path(sa_fig, "FigC_4panel_trajectories.pdf"), p_FigC,
       width = 12, height = 9)
cat(">> Figure C: 4-panel trajectories saved.\n")

# ---- 4. FIGURE E (appendix) –  Response-surface heatmaps --------------------

if (!is.null(res_full)) {
  res_h <- as.data.frame(res_full)
  
  make_heat <- function(oname, lab, pal = "C") {
    ggplot(res_h, aes(x = phi1, y = phi2, colour = .data[[oname]])) +
      geom_point(size = 0.6, alpha = 0.5) +
      scale_colour_viridis_c(option = pal, name = NULL) +
      labs(x = expression(varphi[1]), y = expression(varphi[2]), title = lab) +
      theme_sa(base_size = 8) +
      theme(legend.position = "right", legend.key.height = unit(0.4, "cm"))
  }
  
  p_h1 <- make_heat("temperature_2100", "Temperature 2100 (\u00B0C)", "B")
  p_h2 <- make_heat("tropical_stock_2100", "Tropical forest (bn m\u00B3)", "D")
  p_h3 <- make_heat("output_2100", "Output 2100 (tril. $)", "C")
  p_h4 <- make_heat("inflation_2100", "Inflation 2100 (%)", "A")
  
  p_FigE <- (p_h1 | p_h2) / (p_h3 | p_h4)

  ggsave(file.path(sa_fig, "FigE_heatmaps_2x2.png"), p_FigE,
         width = 11, height = 9, dpi = 300)
  ggsave(file.path(sa_fig, "FigE_heatmaps_2x2.pdf"), p_FigE,
         width = 11, height = 9)
  cat(">> Figure E: Heatmaps saved.\n")
}

# ---- 5. FIGURE F (appendix)  –  6-panel trajectory ribbons -------------------

ribbon_vars   <- c("MC", "tropical_stock", "net_emissions",
                   "temperature", "output", "inflation_pct")
ribbon_labels <- c(
  MC             = "C_F (tril. $)",
  tropical_stock = "Tropical forest (bn m\u00B3)",
  net_emissions  = "Net emissions (GtCO2-e/yr)",
  temperature    = "Temperature (\u00B0C)",
  output         = "Output (tril. $)",
  inflation_pct  = "Inflation rate (%)"
)

ribbon_colours <- c(
  MC = "#0072B2", tropical_stock = "#009E73", net_emissions = "#56B4E9",
  temperature = "#CC79A7", output = "#E69F00", inflation_pct = "#D55E00"
)

ribbon_plots <- list()
for (vname in ribbon_vars) {
  env <- make_envelope(dt_traj, vname)
  base_var <- if (vname %in% names(out_base$ts)) out_base$ts[[vname]] else NULL
  
  p <- ggplot(env, aes(x = year)) +
    geom_ribbon(aes(ymin = lo90, ymax = hi90),
                fill = ribbon_colours[vname], alpha = 0.15) +
    geom_ribbon(aes(ymin = lo50, ymax = hi50),
                fill = ribbon_colours[vname], alpha = 0.30) +
    geom_line(aes(y = med), colour = ribbon_colours[vname], linewidth = 0.7)
  
  if (!is.null(base_var)) {
    p <- p + geom_line(data = out_base$ts, aes(x = year, y = .data[[vname]]),
                       colour = "black", linewidth = 0.7, linetype = "dashed")
  }
  
  p <- p +
    labs(x = "Year", y = ribbon_labels[vname], title = ribbon_labels[vname]) +
    theme_sa(base_size = 8)
  
  ribbon_plots[[vname]] <- p
}

p_FigF <- wrap_plots(ribbon_plots, ncol = 3)

ggsave(file.path(sa_fig, "FigF_6panel_ribbons.png"), p_FigF,
       width = 14, height = 9, dpi = 300)
ggsave(file.path(sa_fig, "FigF_6panel_ribbons.pdf"), p_FigF,
       width = 14, height = 9)
cat(">> Figure F: 6-panel ribbons (2x3) saved.\n")

# ---- 6. Copy existing OAT spider plot as Figure G ---------------------------

for (ext in c("png", "pdf")) {
  src <- file.path("sensitivity_analysis/figures", paste0("local_OAT_summary.", ext))
  dst <- file.path(sa_fig, paste0("FigG_OAT_spider.", ext))
  if (file.exists(src)) file.copy(src, dst, overwrite = TRUE)
}
cat(">> Figure G: OAT spider plot copied.\n")

# ---- 7. MC(RE) fan chart with literature overlay (Figure A alt) --------------

RE_grid <- seq(0, 4, length.out = 200)

# Literature curves
lit_curves <- data.frame(
  source = rep(c("Kindermann (central)", "Busch & Engelmann (2018)",
                  "Havlik (2014)", "Lubowski & Rose (2013)"), each = length(RE_grid)),
  RE = rep(RE_grid, 4),
  MC = c(
    1.527 * exp(1.183 * RE_grid),
    1.1   * exp(0.9   * RE_grid),
    1.8   * exp(1.4   * RE_grid),
    1.2   * exp(0.85  * RE_grid)
  )
)

# Envelope from ±30%
set.seed(55)
n_env <- 500
phi1_s <- runif(n_env, BASELINE_PHI1 * 0.7, BASELINE_PHI1 * 1.3)
phi2_s <- runif(n_env, BASELINE_PHI2 * 0.7, BASELINE_PHI2 * 1.3)

env_RE <- data.table::rbindlist(lapply(seq_len(n_env), function(i) {
  data.frame(RE = RE_grid, MC = phi1_s[i] * exp(phi2_s[i] * RE_grid))
}))
env_RE_stats <- env_RE[, .(
  med  = median(MC),
  lo90 = quantile(MC, 0.05),
  hi90 = quantile(MC, 0.95)
), by = RE]

# Clip for readability
y_clip <- quantile(env_RE$MC, 0.97)

p_FigA_alt <- ggplot() +
  geom_ribbon(data = env_RE_stats, aes(x = RE, ymin = lo90, ymax = pmin(hi90, y_clip)),
              fill = "grey80", alpha = 0.5) +
  geom_line(data = env_RE_stats, aes(x = RE, y = pmin(med, y_clip)),
            colour = "grey40", linewidth = 0.6) +
  geom_line(data = lit_curves, aes(x = RE, y = pmin(MC, y_clip),
                                    colour = source, linetype = source),
            linewidth = 0.7) +
  geom_line(data = data.frame(RE = RE_grid,
                               MC = BASELINE_PHI1 * exp(BASELINE_PHI2 * RE_grid)),
            aes(x = RE, y = pmin(MC, y_clip)),
            colour = "black", linewidth = 1.1, linetype = "solid") +
  coord_cartesian(ylim = c(0, y_clip)) +
  scale_colour_brewer(palette = "Set1", name = "Literature") +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "twodash"),
                        name = "Literature") +
  labs(
    x        = expression("Avoided deforestation RE (GtCO"[2]*"-e)"),
    y        = expression("MC (tril. 2015 US$)"),
    title    = expression("MC(RE) = " * phi[1] * " exp(" * phi[2] * " RE): baseline vs literature"),
    subtitle = "Grey band = 90% envelope under \u00B130% uniform uncertainty"
  ) +
  theme_sa(base_size = 10) +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.8)))

ggsave(file.path(sa_fig, "FigA_MC_literature_fan.png"), p_FigA_alt,
       width = 9, height = 6.5, dpi = 300)
ggsave(file.path(sa_fig, "FigA_MC_literature_fan.pdf"), p_FigA_alt,
       width = 9, height = 6.5)
cat(">> Figure A (alt): MC + literature fan chart saved.\n")

cat("\n>> 03_publication_figures.R finished.\n")
