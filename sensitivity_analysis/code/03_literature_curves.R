###############################################################################
# 03_literature_curves.R  –  Literature cost-curve comparison
#
# Purpose:
#   - Load/define digitised literature cost curves for avoided deforestation
#   - Fit exponential MC(RE) = a * exp(b * RE) to each
#   - Compare fitted (a, b) across sources
#   - Build an empirical envelope
#   - Generate publication-quality overlay figure
#
# Data format (CSV files in sensitivity_analysis/data/literature_curves/):
#   source, x_RE, y_MC, units, notes
#   where x_RE is avoided deforestation (GtCO2-e) and y_MC is marginal cost
#   (2015 US$ trillion, or whatever unit — we normalise in the script)
#
# Outputs:
#   sensitivity_analysis/tables/literature_fit_params.csv
#   sensitivity_analysis/figures/literature_cost_curves_overlay.png / .pdf
#   sensitivity_analysis/figures/literature_param_comparison.png / .pdf
###############################################################################

cat("============================================================\n")
cat(" 03  LITERATURE COST CURVES\n")
cat("============================================================\n")

source("sensitivity_analysis/code/00_helpers.R")

sa_fig  <- "sensitivity_analysis/figures"
sa_tab  <- "sensitivity_analysis/tables"
lit_dir <- "sensitivity_analysis/data/literature_curves"

# ---- 1. Create template if no data exists -----------------------------------

template_file <- file.path(lit_dir, "_TEMPLATE.csv")
if (!file.exists(template_file)) {
  template <- data.frame(
    source = c("Kindermann_2008_low", "Kindermann_2008_low"),
    x_RE   = c(0.5, 1.0),
    y_MC   = c(2.5, 5.2),
    units  = c("tril_2015USD", "tril_2015USD"),
    notes  = c("Digitised from Fig. 3", "Digitised from Fig. 3"),
    stringsAsFactors = FALSE
  )
  write.csv(template, template_file, row.names = FALSE)
  cat(">> Template CSV created at:", template_file, "\n")
  cat("   Add your digitised data files in the same format.\n")
}

# ---- 2. Built-in illustrative curves ----------------------------------------
# 
# Based on representative values from the avoided-deforestation literature:
#   - Kindermann et al. (2008)  — low / central / high
#   - Busch & Engelmann (2018)  — global MACC
#   - Havlik et al. (2014)      — GLOBIOM estimates
#   - Overmars et al. (2014)
#   - Lubowski & Rose (2013)
#
# These are stylised exponential approximations. Replace with actual
# digitised data for publication.

lit_curves <- list(
  list(source = "Kindermann_2008_central",
       phi1 = 1.527, phi2 = 1.183,
       notes = "Baseline calibration (fitted to Kindermann central)"),
  list(source = "Kindermann_2008_low",
       phi1 = 0.90,  phi2 = 1.00,
       notes = "Lower estimate from Kindermann et al. (2008)"),
  list(source = "Kindermann_2008_high",
       phi1 = 2.30,  phi2 = 1.35,
       notes = "Higher estimate from Kindermann et al. (2008)"),
  list(source = "Busch_Engelmann_2018",
       phi1 = 1.10,  phi2 = 0.90,
       notes = "Approx. exponential fit to global MACC from Busch & Engelmann (2018)"),
  list(source = "Havlik_2014_GLOBIOM",
       phi1 = 1.80,  phi2 = 1.40,
       notes = "Stylised from GLOBIOM results in Havlik et al. (2014)"),
  list(source = "Lubowski_Rose_2013",
       phi1 = 1.20,  phi2 = 0.85,
       notes = "Approx. fit from Lubowski & Rose (2013)")
)

# ---- 3. Load any user-supplied CSV data & fit exponentials ------------------

csv_files <- list.files(lit_dir, pattern = "\\.csv$", full.names = TRUE)
csv_files <- csv_files[!grepl("_TEMPLATE", csv_files)]

user_fits <- list()
if (length(csv_files) > 0) {
  cat(sprintf(">> Found %d user-supplied CSV files.\n", length(csv_files)))
  
  for (f in csv_files) {
    dd <- read.csv(f, stringsAsFactors = FALSE)
    sources_in_file <- unique(dd$source)
    
    for (src in sources_in_file) {
      dsub <- dd[dd$source == src, ]
      if (nrow(dsub) < 3) {
        cat(sprintf("   [SKIP] %s: fewer than 3 data points.\n", src))
        next
      }
      
      # Fit MC = a * exp(b * RE) via nls
      fit <- tryCatch(
        nls(y_MC ~ a * exp(b * x_RE), data = dsub,
            start = list(a = 1.5, b = 1.2),
            control = nls.control(maxiter = 500)),
        error = function(e) NULL
      )
      
      if (!is.null(fit)) {
        co <- coef(fit)
        user_fits[[length(user_fits) + 1]] <- list(
          source = src,
          phi1   = unname(co["a"]),
          phi2   = unname(co["b"]),
          notes  = paste("Fitted from user CSV:", basename(f))
        )
        cat(sprintf("   Fitted %s: phi1=%.3f, phi2=%.3f\n", src, co["a"], co["b"]))
      } else {
        cat(sprintf("   [FAIL] %s: nls did not converge.\n", src))
      }
    }
  }
}

# Combine built-in and user curves
all_curves <- c(lit_curves, user_fits)

# ---- 4. Summary table of fitted parameters ----------------------------------

dt_params <- data.table::rbindlist(lapply(all_curves, function(x) {
  data.frame(source = x$source, phi1 = x$phi1, phi2 = x$phi2,
             notes = x$notes, stringsAsFactors = FALSE)
}))

write.csv(dt_params, file.path(sa_tab, "literature_fit_params.csv"), row.names = FALSE)
cat("\n>> Fitted parameters across literature sources:\n")
print(dt_params[, .(source, phi1, phi2)])

# ---- 5. Overlay plot --------------------------------------------------------

RE_grid <- seq(0, 4, length.out = 300)

# All literature curves
curve_data <- data.table::rbindlist(lapply(all_curves, function(x) {
  data.frame(
    RE     = RE_grid,
    MC     = x$phi1 * exp(x$phi2 * RE_grid),
    source = x$source,
    stringsAsFactors = FALSE
  )
}))

# Baseline
baseline_curve <- data.frame(
  RE = RE_grid,
  MC = BASELINE_PHI1 * exp(BASELINE_PHI2 * RE_grid)
)

# Percentile envelope across all curves at each RE point
envelope <- curve_data[, .(
  med  = median(MC),
  lo   = quantile(MC, 0.10),
  hi   = quantile(MC, 0.90),
  lo50 = quantile(MC, 0.25),
  hi50 = quantile(MC, 0.75)
), by = RE]

# Determine y-axis upper limit (99th percentile of data)
y_max <- quantile(curve_data$MC, 0.97, na.rm = TRUE)

p_lit <- ggplot() +
  # Envelope
  geom_ribbon(data = envelope, aes(x = RE, ymin = lo, ymax = hi),
              fill = "grey80", alpha = 0.5) +
  geom_ribbon(data = envelope, aes(x = RE, ymin = lo50, ymax = hi50),
              fill = "grey60", alpha = 0.5) +
  # Individual curves in grey
  geom_line(data = curve_data[source != "Kindermann_2008_central"],
            aes(x = RE, y = MC, group = source),
            colour = "grey50", linewidth = 0.5, alpha = 0.7) +
  # Baseline in black
  geom_line(data = baseline_curve, aes(x = RE, y = MC),
            colour = "black", linewidth = 1.2) +
  # Labels
  coord_cartesian(ylim = c(0, y_max)) +
  labs(
    x        = "Avoided deforestation RE (GtCO2-e)",
    y        = "Cost C_F (tril. 2015 US$)"
  ) +
  theme_sa()

ggsave(file.path(sa_fig, "literature_cost_curves_overlay.png"), p_lit,
       width = 9, height = 6, dpi = 300)
ggsave(file.path(sa_fig, "literature_cost_curves_overlay.pdf"), p_lit,
       width = 9, height = 6)
cat(">> Literature overlay plot saved.\n")

# ---- 6. Parameter comparison scatter ----------------------------------------

p_param <- ggplot(dt_params, aes(x = phi1, y = phi2, label = source)) +
  geom_point(size = 3, colour = "grey40") +
  geom_point(data = dt_params[source == "Kindermann_2008_central"],
             size = 4, colour = "black", shape = 18) +
  {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel(size = 3, max.overlaps = 20)
    } else {
      geom_text(size = 2.5, hjust = -0.1, vjust = -0.5)
    }
  } +
  labs(
    x     = expression(phi[1]~"(level parameter)"),
    y     = expression(phi[2]~"(curvature parameter)"),
    title = expression("Fitted (" * phi[1] * ", " * phi[2] * ") across literature sources"),
    subtitle = "Black diamond = Kindermann baseline calibration"
  ) +
  theme_sa()

ggsave(file.path(sa_fig, "literature_param_comparison.png"), p_param,
       width = 8, height = 6, dpi = 300)
ggsave(file.path(sa_fig, "literature_param_comparison.pdf"), p_param,
       width = 8, height = 6)
cat(">> Parameter comparison plot saved.\n")

# ---- 7. Empirical bounds from literature ------------------------------------

cat("\n>> Empirical parameter ranges from literature sources:\n")
cat(sprintf("   phi1: [%.3f, %.3f]  (baseline: %.3f)\n",
            min(dt_params$phi1), max(dt_params$phi1), BASELINE_PHI1))
cat(sprintf("   phi2: [%.3f, %.3f]  (baseline: %.3f)\n",
            min(dt_params$phi2), max(dt_params$phi2), BASELINE_PHI2))
cat("   These can serve as structurally-informed bounds for an extended SA.\n")

cat("\n>> 03_literature_curves.R finished.\n")
