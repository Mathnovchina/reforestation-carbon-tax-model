###############################################################################
# 03_viability_conditional.R  -  Stage 2 (refined): viability + conditional SA
#
# Separates two questions:
#   (1) MODEL VIABILITY: which parameters push the coupled economy-climate-forest
#       system into the non-viable (collapse / divergence) basin?  ~11.6% of the
#       plausible 14-D space is non-viable (NaN or non-physical explosion); this
#       is genuine model behaviour (finer RK4 steps do NOT rescue it), echoing the
#       debt-deflation bifurcation of Bovari (2018).
#       -> variance-based Sobol on the binary viability indicator (full cube)
#       -> standardized logistic regression (interpretable complement)
#
#   (2) CONDITIONAL SENSITIVITY: among viable trajectories, which parameters drive
#       each system output?  -> standardized regression coefficients (SRC, with the
#       linear-model R^2 as an adequacy check) and Spearman rank correlation, both
#       computed on the viable subset (imputation-free).
#
# Reuses the design already evaluated in 02_sobol_reduced.R (sobol_raw_N2048.rds).
#
# Run from the PROJECT ROOT:
#   source("sensitivity_analysis/code_fullsystem/03_viability_conditional.R")
###############################################################################

source("sensitivity_analysis/code_fullsystem/00_helpers_fullsystem.R")
suppressPackageStartupMessages(library(sensobol))

set.seed(20260730)

out_dir <- "sensitivity_analysis/outputs_fullsystem"
fig_dir <- "sensitivity_analysis/figures_fullsystem"
tab_dir <- "sensitivity_analysis/tables_fullsystem"
for (d in c(out_dir, fig_dir, tab_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- Load the already-evaluated Sobol design -------------------------------
raw <- readRDS(file.path(out_dir, "sobol_raw_N2048.rds"))
Y      <- raw$Y
X      <- as.data.frame(raw$Xphys)
pnames <- raw$params
N      <- raw$N
k      <- length(pnames)
stopifnot(nrow(Y) == N * (k + 2))

# ---- Viability indicator ---------------------------------------------------
# Non-viable = any non-finite output OR any level output beyond a generous
# physical bound (collapse/divergence). Level outputs can explode to ~1e260.
LEVEL <- c("tropical_stock_2100", "total_forest_stock_2100", "sequestration_2100",
           "cumulative_sequestration", "net_emissions_2100", "cumulative_emissions")
finite_ok <- apply(Y, 1, function(r) all(is.finite(r)))
mag_ok    <- apply(Y[, LEVEL, drop = FALSE], 1, function(r) all(abs(r) < 1e5))
viable    <- finite_ok & mag_ok
frac_nv   <- mean(!viable)

cat(sprintf(">> Viable: %d / %d (%.2f%%)  |  non-viable: %d (%.2f%%)\n",
            sum(viable), length(viable), 100 * mean(viable),
            sum(!viable), 100 * frac_nv))

# =============================================================================
# (1) VIABILITY: variance-based Sobol on the binary indicator (full cube)
# =============================================================================
y_viab <- as.numeric(viable)   # aligned to the A/B/AB design row order
sob_v  <- sobol_indices(matrices = c("A", "B", "AB"), Y = y_viab, N = N,
                        params = pnames, boot = TRUE, R = 1000,
                        order = "first", type = "QRN")
viab_sobol <- sob_v$results
write.csv(viab_sobol, file.path(tab_dir, "viability_sobol_indices.csv"), row.names = FALSE)

cat("\n=== Viability Sobol total-order (Ti), descending ===\n")
vt <- viab_sobol[viab_sobol$sensitivity == "Ti", ]
vt <- vt[order(-vt$original), ]
print(vt[, c("parameters", "original", "low.ci", "high.ci")], digits = 3, row.names = FALSE)

# Standardized logistic regression (interpretable complement)
Xs  <- as.data.frame(scale(X)); names(Xs) <- pnames
nv  <- as.integer(!viable)
fit <- suppressWarnings(glm(nv ~ ., data = cbind(nv = nv, Xs), family = binomial()))
co  <- as.data.frame(summary(fit)$coefficients)
co$parameter <- rownames(co); co <- co[co$parameter != "(Intercept)", ]
co <- co[order(-abs(co$`z value`)), ]
names(co)[names(co) == "z value"] <- "z_value"
write.csv(co[, c("parameter", "Estimate", "Std. Error", "z_value", "Pr(>|z|)")],
          file.path(tab_dir, "viability_logistic.csv"), row.names = FALSE)
cat("\n=== Viability logistic (standardized; +coef => param INCREASES collapse) ===\n")
print(co[, c("parameter", "Estimate", "z_value")], digits = 3, row.names = FALSE)

# =============================================================================
# (2) CONDITIONAL SA on viable runs: SRC (with R^2) + Spearman rho
# =============================================================================
Xv  <- X[viable, , drop = FALSE]
Xvs <- as.data.frame(scale(Xv)); names(Xvs) <- pnames

src_rows  <- list()
spear_rows <- list()
r2_rows   <- list()
for (oname in FS_OUTPUT_NAMES) {
  y <- Y[viable, oname]
  ok <- is.finite(y)
  if (sum(ok) < 0.5 * length(y) || sd(y[ok]) == 0) next
  ys  <- scale(y[ok])
  dat <- cbind(as.data.frame(ys), Xvs[ok, , drop = FALSE]); names(dat)[1] <- "y"
  lmf <- lm(y ~ ., data = dat)
  b   <- coef(lmf)[pnames]                 # standardized betas = SRC
  r2  <- summary(lmf)$r.squared
  src_rows[[oname]] <- data.frame(output = oname, parameter = pnames,
                                  src = as.numeric(b), row.names = NULL)
  r2_rows[[oname]]  <- data.frame(output = oname, R2 = r2)
  rho <- sapply(pnames, function(p) suppressWarnings(
    cor(Xv[ok, p], y[ok], method = "spearman")))
  spear_rows[[oname]] <- data.frame(output = oname, parameter = pnames,
                                    spearman = as.numeric(rho), row.names = NULL)
}
src_tab   <- do.call(rbind, src_rows)
spear_tab <- do.call(rbind, spear_rows)
r2_tab    <- do.call(rbind, r2_rows)
write.csv(src_tab,   file.path(tab_dir, "conditional_src.csv"),      row.names = FALSE)
write.csv(spear_tab, file.path(tab_dir, "conditional_spearman.csv"), row.names = FALSE)
write.csv(r2_tab,    file.path(tab_dir, "conditional_src_R2.csv"),   row.names = FALSE)

cat("\n=== Conditional SA: top-3 |SRC| drivers per output (viable runs; R2) ===\n")
for (oname in FS_OUTPUT_NAMES) {
  s <- src_tab[src_tab$output == oname, ]
  if (nrow(s) == 0) next
  s <- s[order(-abs(s$src)), ]
  r2 <- r2_tab$R2[r2_tab$output == oname]
  cat(sprintf("%-26s (R2=%.2f): %s\n", oname, r2,
      paste(sprintf("%s(%+.2f)", head(s$parameter, 3), head(s$src, 3)), collapse = ", ")))
}

# ---- Figures ---------------------------------------------------------------
# Fig 1: viability Sobol Ti (with 95% CI when available)
vt2 <- vt
has_ci <- all(is.finite(vt2$low.ci)) && all(is.finite(vt2$high.ci))
ytop   <- if (has_ci) max(vt2$high.ci, na.rm = TRUE) else max(vt2$original, na.rm = TRUE)
save_sa_fig(file.path(fig_dir, "viability_sobol"), w = 7, h = 5.5, drawfun = function() {
  op <- par(mar = c(6, 4, 3, 1))
  bp <- barplot(vt2$original, names.arg = param_expr(vt2$parameters), las = 2,
                ylim = c(0, ytop * 1.15),
                col = "firebrick", ylab = expression("Total-order Sobol index" ~ (T[i])),
                cex.names = 0.85)
  if (has_ci) arrows(bp, vt2$low.ci, bp, vt2$high.ci, angle = 90, code = 3, length = 0.04)
  par(op)
})

# Fig 2: SRC heatmap (outputs x parameters)
M <- matrix(NA_real_, nrow = length(FS_OUTPUT_NAMES), ncol = k,
            dimnames = list(FS_OUTPUT_NAMES, pnames))
for (oname in FS_OUTPUT_NAMES) {
  s <- src_tab[src_tab$output == oname, ]
  if (nrow(s)) M[oname, s$parameter] <- s$src
}
save_sa_fig(file.path(fig_dir, "conditional_src_heatmap"), w = 11, h = 6, drawfun = function() {
  lim  <- max(abs(M), na.rm = TRUE)
  cols <- colorRampPalette(c("navy", "white", "firebrick"))(101)
  op   <- par(mar = c(7, 17, 3, 2))
  image(x = seq_len(k), y = seq_len(nrow(M)), z = t(M), col = cols, zlim = c(-lim, lim),
        axes = FALSE, xlab = "", ylab = "")
  axis(1, at = seq_len(k), labels = param_expr(pnames), las = 2, cex.axis = 0.85)
  axis(2, at = seq_len(nrow(M)), labels = fs_label(rownames(M)), las = 1, cex.axis = 0.75)
  box(); par(op)
})

saveRDS(list(viable = viable, frac_nv = frac_nv, viab_sobol = viab_sobol,
             logistic = co, src = src_tab, spearman = spear_tab, r2 = r2_tab),
        file = file.path(out_dir, "viability_conditional.rds"))

cat(sprintf("\nTables -> %s\nFigures -> %s\n", tab_dir, fig_dir))
