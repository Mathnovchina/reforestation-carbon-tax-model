###############################################################################
# 01_morris_screening.R  -  Stage 1: whole-system Morris elementary-effects
#
# Screens all 19 ledger parameters in the focal scenario (forestA_50) to rank
# influence and detect non-linearity/interaction (sigma) before the reduced
# Sobol stage. Outputs mu* / sigma per parameter for each system output.
#
# Run from the PROJECT ROOT:
#   source("sensitivity_analysis/code_fullsystem/01_morris_screening.R")
###############################################################################

source("sensitivity_analysis/code_fullsystem/00_helpers_fullsystem.R")

set.seed(20260730)

# ---- Config ----------------------------------------------------------------
MORRIS_R      <- 20      # trajectories -> N = r*(k+1) = 20*20 = 400 runs
MORRIS_LEVELS <- 6
GRID_JUMP     <- 3

out_dir <- "sensitivity_analysis/outputs_fullsystem"
fig_dir <- "sensitivity_analysis/figures_fullsystem"
log_dir <- "sensitivity_analysis/logs_fullsystem"
for (d in c(out_dir, fig_dir, log_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(log_dir, "morris_run.log")
cat(sprintf("Morris screening start: %s\n", Sys.time()), file = log_file)

# ---- Ledger + design -------------------------------------------------------
led    <- load_ledger()
pnames <- led$param

cat(sprintf("Morris design: %d parameters, r=%d -> %d model runs.\n",
            length(pnames), MORRIS_R, MORRIS_R * (length(pnames) + 1)))

design <- morris_design(
  pnames    = pnames,
  lower     = led$lower,
  upper     = led$upper,
  r         = MORRIS_R,
  levels    = MORRIS_LEVELS,
  grid_jump = GRID_JUMP
)
X <- design$X_real

# ---- Compile + run ---------------------------------------------------------
cat(">> Compiling model under focal scenario (forestA_50)...\n")
sys <- compile_model_focal()

cat(sprintf(">> Running %d Morris design points...\n", nrow(X)))
t0 <- Sys.time()
Y  <- run_design(sys, X, log_file = log_file)
cat(sprintf(">> Done in %.1f s.\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

saveRDS(list(design = design, X = X, Y = Y, ledger = led),
        file = file.path(out_dir, "morris_raw.rds"))

# ---- Elementary effects per output -----------------------------------------
ee_rows <- list()
for (oname in FS_OUTPUT_NAMES) {
  y <- Y[, oname]
  if (all(is.na(y))) { cat(sprintf("  [skip] %s: all NA\n", oname)); next }
  ee <- compute_morris_ee(design, y)
  ee$output <- oname
  ee_rows[[oname]] <- ee[, c("output", "parameter", "mu", "mu_star", "sigma")]
}
ee_tab <- do.call(rbind, ee_rows)
write.csv(ee_tab, file.path(out_dir, "morris_elementary_effects.csv"), row.names = FALSE)

# Rank by mean normalized mu* across outputs (robust cross-output importance)
norm_mu <- do.call(cbind, lapply(FS_OUTPUT_NAMES, function(o) {
  v <- ee_tab$mu_star[ee_tab$output == o]
  if (all(!is.finite(v)) || max(v, na.rm = TRUE) == 0) return(rep(NA_real_, length(pnames)))
  v / max(v, na.rm = TRUE)
}))
rownames(norm_mu) <- ee_tab$parameter[ee_tab$output == FS_OUTPUT_NAMES[1]]
rank_tab <- data.frame(
  parameter    = rownames(norm_mu),
  mean_norm_mu = rowMeans(norm_mu, na.rm = TRUE),
  max_norm_mu  = apply(norm_mu, 1, max, na.rm = TRUE),
  row.names = NULL
)
rank_tab <- rank_tab[order(-rank_tab$mean_norm_mu), ]
write.csv(rank_tab, file.path(out_dir, "morris_importance_ranking.csv"), row.names = FALSE)

cat("\n=== Morris cross-output importance ranking (normalized mu*) ===\n")
print(rank_tab, digits = 3)

# ---- Plots: mu* vs sigma per output ----------------------------------------
# One colour per parameter; shared legend below the 3×4 panel grid.
pnames_all <- unique(ee_tab$parameter)
pal_params <- setNames(
  colorRampPalette(c(
    "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628",
    "#F781BF","#666666","#66C2A5","#FC8D62","#8DA0CB","#E78AC3",
    "#A6D854","#FFD92F","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E"
  ))(length(pnames_all)),
  pnames_all
)

save_sa_fig(file.path(fig_dir, "morris_mustar_sigma"), w = 11, h = 9.5, drawfun = function() {
  layout(matrix(c(1:12, rep(13, 4)), nrow = 4, byrow = TRUE),
         heights = c(1, 1, 1, 0.46))
  op <- par(mar = c(2.2, 3.5, 3.4, 1), mgp = c(2, 0.6, 0), tcl = -0.3)
  for (oname in FS_OUTPUT_NAMES) {
    sub <- ee_tab[ee_tab$output == oname, ]
    if (nrow(sub) == 0 || all(!is.finite(sub$mu_star))) {
      plot.new(); title(main = fs_title(oname), cex.main = 0.8); next
    }
    plot(sub$mu_star, sub$sigma,
         pch = 19, col = pal_params[sub$parameter], cex = 1.1,
         xlab = expression(mu^"*"), ylab = expression(sigma),
         main = fs_title(oname), cex.main = 0.8, cex.axis = 0.75, cex.lab = 0.8)
    # Smart label pos: flip below near top, flip left near right edge
    usr <- par("usr"); xr <- usr[2] - usr[1]; yr <- usr[4] - usr[3]
    pos_vec <- ifelse(sub$sigma   > usr[3] + 0.78 * yr, 1L,
               ifelse(sub$mu_star > usr[1] + 0.78 * xr, 2L, 3L))
    text(sub$mu_star, sub$sigma,
         labels = param_expr(sub$parameter),
         pos = pos_vec, cex = 0.55, xpd = NA)
  }
  # Legend panel — paper symbols via plotmath, close to grid (small top mar)
  par(mar = c(0.2, 0, 0, 0))
  plot.new()
  legend("center",
         legend = param_expr(pnames_all),
         col    = pal_params[pnames_all],
         pch    = 19, pt.cex = 1.3,
         ncol   = ceiling(length(pnames_all) / 2),
         bty    = "n", cex = 1.0, x.intersp = 0.7, y.intersp = 1.15)
  par(op)
  layout(1)
})

cat(sprintf("\nOutputs written to:\n  %s\n  %s\n  %s\n",
            file.path(out_dir, "morris_elementary_effects.csv"),
            file.path(out_dir, "morris_importance_ranking.csv"),
            file.path(fig_dir, "morris_mustar_sigma.pdf")))
