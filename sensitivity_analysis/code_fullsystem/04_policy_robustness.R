###############################################################################
# 04_policy_robustness.R  -  Stage 3: paired policy sign-robustness
#
# Question answered: does the SIGN of every forest-policy conclusion survive the
# 14-parameter uncertainty? (Sobol/viability say WHAT drives outcomes; this says
# whether the policy CONCLUSION is robust.)
#
# Design (lean; low-damage only, matching the reported Forest A/B scenarios):
#   - 6 systems compiled once each: {Baseline, Forest A, Forest B} x {tax~50, ~300}
#   - one space-filling ensemble of N draws over the 14 retained parameters
#   - the SAME draw is pushed through all 6 systems (paired)
#   - report P( policy - matched-baseline keeps its expected sign ) over the
#     jointly-viable draws, for the policy-relevant outputs
#   - plus a carbon-tax stress test: P( tax~300 worsens tropical stock vs tax~50 )
#
# The damage axis is intentionally NOT swept: Forest A/B are only ever reported
# in the low-damage case, so a damage cross-product would add cost without
# informing any claim in the paper.
#
# Run from the PROJECT ROOT:
#   source("sensitivity_analysis/code_fullsystem/04_policy_robustness.R")
###############################################################################

source("sensitivity_analysis/code_fullsystem/00_helpers_fullsystem.R")
suppressPackageStartupMessages(library(sensobol))

set.seed(20260731)

# ---- Config ----------------------------------------------------------------
FIXED_PARAMS <- c("CUP_pind", "g_pbs", "delta_g_sigma", "gamma_ast", "delta_CH")
N_DRAWS      <- 1024
VIAB_CAP     <- 1e5   # same viability rule as Stage-2 refined

out_dir <- "sensitivity_analysis/outputs_fullsystem"
fig_dir <- "sensitivity_analysis/figures_fullsystem"
tab_dir <- "sensitivity_analysis/tables_fullsystem"
log_dir <- "sensitivity_analysis/logs_fullsystem"
for (d in c(out_dir, fig_dir, tab_dir, log_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(log_dir, "robustness_run.log")
cat(sprintf("Policy robustness start: %s\n", Sys.time()), file = log_file)

# ---- Scenario specs (low damage; mirror run_model.R SCENARIOS) -----------
# g_mid = 0.04 -> ~58 $/tonC ("50");  g_high = 0.06 -> ~308 $/tonC ("300").
SCN <- list(
  baseline_50  = list(g_p_car = 0.04, RD_slope = 0,    MC_on = 0, MC_mode = 1, dam_on = 1, f_K = 0),
  baseline_300 = list(g_p_car = 0.06, RD_slope = 0,    MC_on = 0, MC_mode = 1, dam_on = 1, f_K = 0),
  forestA_50   = list(g_p_car = 0.04, RD_slope = 1/42, MC_on = 1, MC_mode = 1, dam_on = 1, f_K = 0),
  forestA_300  = list(g_p_car = 0.06, RD_slope = 1/42, MC_on = 1, MC_mode = 1, dam_on = 1, f_K = 0),
  forestB_50   = list(g_p_car = 0.04, RD_slope = 1/22, MC_on = 1, MC_mode = 1, dam_on = 1, f_K = 0),
  forestB_300  = list(g_p_car = 0.06, RD_slope = 1/22, MC_on = 1, MC_mode = 1, dam_on = 1, f_K = 0)
)

# Policy-vs-baseline pairings and the tax stress-test pairings.
POLICY_PAIRS <- list(
  c(policy = "forestA_50",  base = "baseline_50"),
  c(policy = "forestA_300", base = "baseline_300"),
  c(policy = "forestB_50",  base = "baseline_50"),
  c(policy = "forestB_300", base = "baseline_300")
)
TAX_PAIRS <- list(
  c(high = "forestA_300", low = "forestA_50"),
  c(high = "forestB_300", low = "forestB_50")
)

# Outputs whose policy sign is a stated conclusion of the paper.
KEY_OUTPUTS <- c("tropical_stock_2100", "total_forest_stock_2100",
                 "cumulative_sequestration", "temperature_2100",
                 "cumulative_emissions", "private_debt_ratio_2100",
                 "avg_inflation_2050_2100", "CF_2100")

# ---- Ledger + 14-parameter ensemble ----------------------------------------
led    <- load_ledger()
keep   <- !(led$param %in% FIXED_PARAMS)
pnames <- led$param[keep]
lower  <- setNames(led$lower[keep], pnames)
upper  <- setNames(led$upper[keep], pnames)
base14 <- setNames(led$baseline[keep], pnames)
k      <- length(pnames)

cat(sprintf(">> Ensemble: %d draws over %d parameters\n", N_DRAWS, k))
cat(sprintf(">> Retained: %s\n", paste(pnames, collapse = ", ")))

# Space-filling A-block (Sobol' QRN), scaled to physical bounds.
A <- sobol_matrices(matrices = "A", N = N_DRAWS, params = pnames, type = "QRN")
X <- A
for (p in pnames) X[, p] <- lower[p] + A[, p] * (upper[p] - lower[p])

# ---- Run the same ensemble through all 6 systems ---------------------------
Ylist <- list()
ref   <- matrix(NA_real_, nrow = length(SCN), ncol = length(FS_OUTPUT_NAMES),
                dimnames = list(names(SCN), FS_OUTPUT_NAMES))
for (nm in names(SCN)) {
  cat(sprintf("\n>> Compiling + running scenario '%s' ...\n", nm))
  sys <- compile_model_focal(SCN[[nm]])
  ref[nm, ] <- safe_run_named(sys, as.list(base14))   # deterministic reference draw
  t0  <- Sys.time()
  Ylist[[nm]] <- run_design(sys, X, log_file = log_file, verbose = TRUE)
  cat(sprintf(">> '%s': %d runs in %.1f s\n", nm,
              nrow(X), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# ---- Viability mask per scenario -------------------------------------------
LEVEL_OUTPUTS <- c("tropical_stock_2100", "total_forest_stock_2100",
                   "net_emissions_2100", "temperature_2100",
                   "private_debt_ratio_2100", "CF_2100")
viable_mask <- function(Y) {
  fin <- apply(Y, 1, function(r) all(is.finite(r)))
  cap <- apply(Y[, LEVEL_OUTPUTS, drop = FALSE], 1,
               function(r) all(abs(r) < VIAB_CAP))
  fin & cap
}
viab <- sapply(Ylist, viable_mask)   # N x 6 logical
cat("\n>> Viable draws per scenario:\n")
print(colSums(viab))

# ---- Policy-vs-baseline sign robustness ------------------------------------
sign_rows <- list()
for (pr in POLICY_PAIRS) {
  pol <- pr[["policy"]]; bas <- pr[["base"]]
  both <- viab[, pol] & viab[, bas]
  nb   <- sum(both)
  for (o in KEY_OUTPUTS) {
    d_ref  <- ref[pol, o] - ref[bas, o]
    esign  <- sign(d_ref)
    d_ens  <- Ylist[[pol]][both, o] - Ylist[[bas]][both, o]
    p_sign <- if (esign == 0) NA_real_ else mean(sign(d_ens) == esign)
    sign_rows[[length(sign_rows) + 1]] <- data.frame(
      comparison = sprintf("%s - %s", pol, bas),
      output = o, label = unname(fs_label(o)),
      expected_sign = ifelse(esign > 0, "+", ifelse(esign < 0, "-", "0")),
      ref_delta = d_ref,
      median_delta = median(d_ens),
      p_sign_robust = p_sign,
      n_joint_viable = nb,
      stringsAsFactors = FALSE
    )
  }
}
sign_tab <- do.call(rbind, sign_rows)

# ---- Tax stress test: does higher tax worsen tropical stock? ---------------
tax_rows <- list()
for (tp in TAX_PAIRS) {
  hi <- tp[["high"]]; lo <- tp[["low"]]
  both <- viab[, hi] & viab[, lo]
  for (o in c("tropical_stock_2100", "total_forest_stock_2100")) {
    d_ref <- ref[hi, o] - ref[lo, o]           # expected negative (tax worsens)
    esign <- sign(d_ref)
    d_ens <- Ylist[[hi]][both, o] - Ylist[[lo]][both, o]
    tax_rows[[length(tax_rows) + 1]] <- data.frame(
      comparison = sprintf("%s - %s", hi, lo),
      output = o, label = unname(fs_label(o)),
      expected_sign = ifelse(esign > 0, "+", ifelse(esign < 0, "-", "0")),
      ref_delta = d_ref,
      median_delta = median(d_ens),
      p_sign_robust = if (esign == 0) NA_real_ else mean(sign(d_ens) == esign),
      n_joint_viable = sum(both),
      stringsAsFactors = FALSE
    )
  }
}
tax_tab <- do.call(rbind, tax_rows)

# ---- Write tables + raw ----------------------------------------------------
write.csv(sign_tab, file.path(tab_dir, "policy_sign_robustness.csv"), row.names = FALSE)
write.csv(tax_tab,  file.path(tab_dir, "tax_effect_robustness.csv"),  row.names = FALSE)
saveRDS(list(X = X, Ylist = Ylist, viab = viab, ref = ref,
             sign_tab = sign_tab, tax_tab = tax_tab, params = pnames),
        file.path(out_dir, "policy_robustness.rds"))

# ---- Figure: sign-robustness heatmap ---------------------------------------
pairs_lab <- unique(sign_tab$comparison)
M <- matrix(NA_real_, nrow = length(KEY_OUTPUTS), ncol = length(pairs_lab),
            dimnames = list(KEY_OUTPUTS, pairs_lab))
for (i in seq_len(nrow(sign_tab)))
  M[sign_tab$output[i], sign_tab$comparison[i]] <- sign_tab$p_sign_robust[i]

save_sa_fig(file.path(fig_dir, "policy_sign_robustness"), w = 11, h = 6, drawfun = function() {
  op <- par(mar = c(9, 17, 3, 2))
  image(x = seq_along(pairs_lab), y = seq_len(nrow(M)), z = t(M),
        col = colorRampPalette(c("white", "steelblue"))(100), zlim = c(0.5, 1),
        axes = FALSE, xlab = "", ylab = "",
        main = "P(policy effect keeps expected sign) across 14-param uncertainty")
  axis(1, at = seq_along(pairs_lab), labels = pairs_lab, las = 2, cex.axis = 0.8)
  axis(2, at = seq_len(nrow(M)), labels = fs_label(rownames(M)), las = 1, cex.axis = 0.75)
  for (ii in seq_along(pairs_lab)) for (jj in seq_len(nrow(M)))
    if (is.finite(M[jj, ii])) text(ii, jj, sprintf("%.2f", M[jj, ii]), cex = 0.7)
  box(); par(op)
})

# ---- Console summary -------------------------------------------------------
cat("\n================ POLICY SIGN ROBUSTNESS ================\n")
print(sign_tab[, c("comparison", "label", "expected_sign",
                   "median_delta", "p_sign_robust", "n_joint_viable")],
      row.names = FALSE, digits = 3)
cat("\n================ TAX STRESS TEST (higher tax vs lower) ================\n")
print(tax_tab[, c("comparison", "label", "expected_sign",
                  "median_delta", "p_sign_robust", "n_joint_viable")],
      row.names = FALSE, digits = 3)
cat("\nTables -> ", tab_dir, "\nFigure -> ", fig_dir, "\n", sep = "")
