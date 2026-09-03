###############################################################################
# replot_morris_units.R  -  Regenerate ONLY the Morris mu*/sigma figure with
# per-panel physical units, reusing the already-computed elementary effects
# (no model recompile / re-run needed).
#
#   source("sensitivity_analysis/code_fullsystem/replot_morris_units.R")
###############################################################################

source("sensitivity_analysis/code_fullsystem/00_helpers_fullsystem.R")

out_dir <- "sensitivity_analysis/outputs_fullsystem"
fig_dir <- "sensitivity_analysis/figures_fullsystem"

ee_tab <- read.csv(file.path(out_dir, "morris_elementary_effects.csv"),
                   stringsAsFactors = FALSE)

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
    usr <- par("usr"); xr <- usr[2] - usr[1]; yr <- usr[4] - usr[3]
    pos_vec <- ifelse(sub$sigma   > usr[3] + 0.78 * yr, 1L,
               ifelse(sub$mu_star > usr[1] + 0.78 * xr, 2L, 3L))
    text(sub$mu_star, sub$sigma,
         labels = param_expr(sub$parameter),
         pos = pos_vec, cex = 0.55, xpd = NA)
  }
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

cat("Re-plotted:", file.path(fig_dir, "morris_mustar_sigma.pdf"), "+ .eps\n")
