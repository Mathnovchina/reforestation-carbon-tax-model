# =====================================================================================
#  RUN THE MODEL  --  build the RK4 solver from model.R and reproduce the paper figures
# =====================================================================================
#  Set the working directory to the repository root (the folder containing model.R and
#  SourceCode.R), then run:  source("run_model.R")
#  Figures are written to figures/ as both .png (preview) and .eps (publication).

source("SourceCode.R")


# =====================================================================================
#  POLICY DASHBOARD  --  choose the scenario(s) to simulate
# =====================================================================================
#  Every scenario is fully described by six switches that are injected into model.R:
#    g_p_car   : growth rate of the carbon price  -> sets the 2100 carbon-tax level
#    RD_slope  : slope of the deforestation-control rate  RD = RD_slope * t
#                  0    = Baseline (no forest policy)
#                  1/84 = RD reaches 1 (deforestation halted) in 2100
#                  1/42 = RD reaches 2 in 2100            (Forest Policy A)
#                  1/22 = RD reaches ~3.8 in 2100         (Forest Policy B)
#    MC_on     : 1 = forest-policy cost C_F active ; 0 = no forest policy  -> C_F = 0
#    MC_mode   : 0 = level C_F=f1*exp(f2*RE) ; 1 = incremental C_F=f1*(exp(f2*RE)-1)
#    dam_on    : 1 = climate damages active ; 0 = no-damage scenario
#    f_K       : 0 = low damage (output only) ; 1/3 = high damage (output + capital)
#
#  Convert a desired 2100 carbon-tax target (US$/tonC) into the growth rate g_p_car.
#  p_car(t) = p_car0 * exp(g_p_car * t), with p_car0 = 2.0024 and horizon 2016->2100.
g_from_target <- function(target_2100, p_car0 = 2.0024, horizon = 84) {
	log(target_2100 / p_car0) / horizon
}

g_bau <- g_from_target(10)   # ~0.0191 -> 2100 carbon tax ~ 10 $/tonC (baseline path)
g_mid <- 0.04                # -> 2100 carbon tax ~ 58 $/tonC (labelled "50" in the paper)
g_high <- 0.06               # -> 2100 carbon tax ~ 308 $/tonC (labelled "300" in the paper)

# MC comparability rule:
# - Baseline runs always have MC_on=0 (no forest-policy cost).
# - Forest-policy runs use MC_mode=1 (incremental form), so C_F(0)=0 and initial
#   ratio variables (e.g., wage share, inflation) are comparable at t0.
# - If strict continuity with old level-form runs is needed, set policy_mc_mode <- 0.
policy_mc_mode <- 1

# Scenario library (name = list of the six switches + a human-readable label)
SCENARIOS <- list(
	# --- baselines on the BAU carbon path (~10 $/tonC), one per damage case ---
	baseline_nodmg = list(g_p_car = g_bau,  RD_slope = 0,    MC_on = 0, MC_mode = 0, dam_on = 0, f_K = 0,   label = "Baseline (no damage)"),
	baseline_low   = list(g_p_car = g_bau,  RD_slope = 0,    MC_on = 0, MC_mode = 0, dam_on = 1, f_K = 0,   label = "Baseline (low damage)"),
	baseline_high  = list(g_p_car = g_bau,  RD_slope = 0,    MC_on = 0, MC_mode = 0, dam_on = 1, f_K = 1/3, label = "Baseline (high damage)"),

	# --- baselines on the two forest-policy carbon paths ("50" and "300"), one per damage case (Figure 2) ---
	baseline_nodmg_50  = list(g_p_car = g_mid,  RD_slope = 0, MC_on = 0, MC_mode = 0, dam_on = 0, f_K = 0,   label = "Baseline no-damage, tax~50"),
	baseline_nodmg_300 = list(g_p_car = g_high, RD_slope = 0, MC_on = 0, MC_mode = 0, dam_on = 0, f_K = 0,   label = "Baseline no-damage, tax~300"),
	baseline_low_50    = list(g_p_car = g_mid,  RD_slope = 0, MC_on = 0, MC_mode = 0, dam_on = 1, f_K = 0,   label = "Baseline low-damage, tax~50"),
	baseline_low_300   = list(g_p_car = g_high, RD_slope = 0, MC_on = 0, MC_mode = 0, dam_on = 1, f_K = 0,   label = "Baseline low-damage, tax~300"),
	baseline_high_50   = list(g_p_car = g_mid,  RD_slope = 0, MC_on = 0, MC_mode = 0, dam_on = 1, f_K = 1/3, label = "Baseline high-damage, tax~50"),
	baseline_high_300  = list(g_p_car = g_high, RD_slope = 0, MC_on = 0, MC_mode = 0, dam_on = 1, f_K = 1/3, label = "Baseline high-damage, tax~300"),

	# --- legacy low-damage baselines kept for back-compat ---
	baseline_50    = list(g_p_car = g_mid,  RD_slope = 0,    MC_on = 0, MC_mode = 0, dam_on = 1, f_K = 0,   label = "Baseline, tax~50 (low damage)"),
	baseline_300   = list(g_p_car = g_high, RD_slope = 0,    MC_on = 0, MC_mode = 0, dam_on = 1, f_K = 0,   label = "Baseline, tax~300 (low damage)"),

	# --- forest policies A (RD_slope=1/42) and B (RD_slope=1/22) at the two carbon paths ---
	forestA_50     = list(g_p_car = g_mid,  RD_slope = 1/42, MC_on = 1, MC_mode = policy_mc_mode, dam_on = 1, f_K = 0,   label = "Forest Policy A, tax~50 (low damage)"),
	forestA_300    = list(g_p_car = g_high, RD_slope = 1/42, MC_on = 1, MC_mode = policy_mc_mode, dam_on = 1, f_K = 0,   label = "Forest Policy A, tax~300 (low damage)"),
	forestB_50     = list(g_p_car = g_mid,  RD_slope = 1/22, MC_on = 1, MC_mode = policy_mc_mode, dam_on = 1, f_K = 0,   label = "Forest Policy B, tax~50 (low damage)"),
	forestB_300    = list(g_p_car = g_high, RD_slope = 1/22, MC_on = 1, MC_mode = policy_mc_mode, dam_on = 1, f_K = 0,   label = "Forest Policy B, tax~300 (low damage)")
)

# Build + run one scenario: sets the global switches, regenerates the C++ from model.R, integrates.
run_scenario <- function(name) {
	scn <- SCENARIOS[[name]]
	if (is.null(scn)) stop("Unknown scenario '", name, "'. Available: ", paste(names(SCENARIOS), collapse = ", "))
	scenario_g_p_car  <<- scn$g_p_car
	scenario_RD_slope <<- scn$RD_slope
	scenario_MC_on    <<- scn$MC_on
	scenario_MC_mode  <<- scn$MC_mode
	scenario_dam_on   <<- scn$dam_on
	scenario_f_K      <<- scn$f_K
	message("Running scenario: ", scn$label)
	sys <- cppMakeSys("model.R", reportVars = 3)
	cppRK4(sys)
}


# =====================================================================================
#  PAPER FIGURES  --  original-article colour palette
# =====================================================================================
#  Line-type convention (used throughout):
#    solid  (lty = 1)  = carbon tax ~ 50  ($/tonC by 2100)  -- less convex price path
#    dashed (lty = 2)  = carbon tax ~ 300 ($/tonC by 2100)  -- more convex price path
#
#  Colour conventions:
#    Figure 2 (baseline macro): damage case -> blue = no damage, gold = low, red = high.
#    Biomes (Figures 3-5): Tropical = goldenrod, Boreal = darkolivegreen4,
#                          Temperate = cornflowerblue.
LTY50  <- 1   # solid  = carbon tax ~ 50
LTY300 <- 2   # dashed = carbon tax ~ 300

PAL <- list(
	# Figure 2 -- damage cases
	dmg_none = "blue",
	dmg_low  = "goldenrod3",
	dmg_high = "red3",
	# biomes (original article palette)
	biome_tro = "goldenrod",        # Tropical
	biome_bor = "darkolivegreen4",  # Boreal
	biome_tem = "cornflowerblue",   # Temperate
	# forest-grid variable palette (original article palette)
	macro    = "black",
	emp      = "blue",
	wage     = "olivedrab",
	debt     = "blue",
	infl     = "olivedrab",
	fossil   = "orange4",
	bioen    = "olivedrab",
	em_gross = "orange4",
	em_net   = "olivedrab",
	# baseline reference overlaid inside the forest-policy grids
	base_ref = "grey50"
)

# -------------------------------------------------------------------------------------
#  Small plotting helpers
# -------------------------------------------------------------------------------------
tvec <- function(r) r$time + 2016

# Draw one panel from a list of series; each series = list(run=, var=, col=, lty=, scale=).
# ylim is auto-computed across all series when NULL, so nothing gets clipped.
panel <- function(series, main, ylim = NULL, ylab = "") {
	getval <- function(s) {
		sc <- if (is.null(s$scale)) 1 else s$scale
		R[[s$run]][[s$var]] * sc
	}
	if (is.null(ylim)) ylim <- range(unlist(lapply(series, getval)), na.rm = TRUE)
	first <- TRUE
	for (s in series) {
		y <- getval(s)
		if (first) {
			plot(tvec(R[[s$run]]), y, type = "l", col = s$col, lty = s$lty, lwd = 2,
				 xlab = "", ylab = ylab, main = main, ylim = ylim)
			first <- FALSE
		} else {
			lines(tvec(R[[s$run]]), y, col = s$col, lty = s$lty, lwd = 2)
		}
	}
}

# Dual-axis panel (two variables, own scales) with the two tax schedules (solid 50 / dashed 300).
dual_panel <- function(run50, run300, var1, col1, var2, col2, main, lab2, lab4) {
	t50 <- tvec(R[[run50]]); t300 <- tvec(R[[run300]])
	y1 <- range(R[[run50]][[var1]], R[[run300]][[var1]], na.rm = TRUE)
	y2 <- range(R[[run50]][[var2]], R[[run300]][[var2]], na.rm = TRUE)
	plot(t50, R[[run50]][[var1]], type = "l", col = col1, lty = LTY50, lwd = 2,
		 xlab = "", ylab = "", ylim = y1, main = main)
	lines(t300, R[[run300]][[var1]], col = col1, lty = LTY300, lwd = 2)
	par(new = TRUE)
	plot(t50, R[[run50]][[var2]], type = "l", col = col2, lty = LTY50, lwd = 2,
		 axes = FALSE, xlab = "", ylab = "", ylim = y2)
	lines(t300, R[[run300]][[var2]], col = col2, lty = LTY300, lwd = 2)
	axis(side = 4)
	mtext(lab2, side = 2, line = 3,   cex = 0.7)
	mtext(lab4, side = 4, line = 2.4, cex = 0.7)
}

# Write a figure to both PNG (screen/preview) and EPS (publication).
save_fig <- function(path, drawfun, w = 1500, h = 1150, res = 130, ew = 11, eh = 8.2) {
	dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
	png(paste0(path, ".png"), width = w, height = h, res = res); drawfun(); dev.off()
	postscript(paste0(path, ".eps"), onefile = FALSE, horizontal = FALSE,
			 paper = "special", width = ew, height = eh); drawfun(); dev.off()
	message("written: ", path, ".png  +  .eps")
}


# =====================================================================================
#  RUN every scenario needed for the three paper figures (each rebuilds + integrates once)
# =====================================================================================
figure_scenarios <- c(
	"baseline_nodmg_50", "baseline_nodmg_300",
	"baseline_low_50",   "baseline_low_300",
	"baseline_high_50",  "baseline_high_300",
	"forestA_50",        "forestA_300",
	"forestB_50",        "forestB_300"
)
R <- setNames(lapply(figure_scenarios, run_scenario), figure_scenarios)
# derived series used in the figures
for (nm in names(R)) R[[nm]]$div_share <- R[[nm]]$Div / (R[[nm]]$p * R[[nm]]$Y)
# derived output growth rate (%/yr) on the fine integration grid
for (nm in names(R)) {
	dt <- diff(R[[nm]]$time)[1]
	R[[nm]]$g_Y <- 100 * c(NA, diff(log(R[[nm]]$Y))) / dt
	R[[nm]]$g_Y[1:2] <- NA   # drop the first RK4 sub-step (start-up artefact)
}


# =====================================================================================
#  FIGURE 2 -- Baseline main (macro) variables: 3 damage cases x 2 carbon-tax schedules
#  Colour = damage case (blue/gold/red) ; line type = tax schedule (solid 50 / dashed 300)
#  Rationale: only in the baseline do we examine the three damage cases, to (a) reproduce
#  the Bovari et al. (2018) outcomes and (b) show how the damage function shapes stability.
# =====================================================================================
fig2_series <- function(var, scale = 1) list(
	list(run = "baseline_nodmg_50",  var = var, col = PAL$dmg_none, lty = LTY50,  scale = scale),
	list(run = "baseline_nodmg_300", var = var, col = PAL$dmg_none, lty = LTY300, scale = scale),
	list(run = "baseline_low_50",    var = var, col = PAL$dmg_low,  lty = LTY50,  scale = scale),
	list(run = "baseline_low_300",   var = var, col = PAL$dmg_low,  lty = LTY300, scale = scale),
	list(run = "baseline_high_50",   var = var, col = PAL$dmg_high, lty = LTY50,  scale = scale),
	list(run = "baseline_high_300",  var = var, col = PAL$dmg_high, lty = LTY300, scale = scale)
)

draw_fig2 <- function(growth = FALSE) {
	par(mfrow = c(3, 3), mar = c(2.6, 4.4, 3.2, 1.2), xpd = FALSE)

	if (growth) {
		panel(fig2_series("g_Y"), "Output growth rate (%/yr)")
		lpos <- "topright"
	} else {
		panel(fig2_series("Y"),   "Output (2015 US$ tril.)")
		lpos <- "topleft"
	}
	legend(lpos, bty = "n", cex = 0.78,
		 legend = c("No damage", "Low damage", "High damage", "tax ~50 (solid)", "tax ~300 (dashed)"),
		 col = c(PAL$dmg_none, PAL$dmg_low, PAL$dmg_high, "black", "black"),
		 lty = c(1, 1, 1, LTY50, LTY300), lwd = 2)

	panel(fig2_series("lambda"),     "Employment rate")
	panel(fig2_series("omega"),      "Wage share")
	panel(fig2_series("div_share"),  "Dividend share")
	panel(fig2_series("pi"),         "Profit share")
	panel(fig2_series("i", 100),     "Inflation (%)")
	panel(fig2_series("d"),          "Private debt ratio (%)")
	panel(fig2_series("E"),          "Emission (GtCO2-e)")
	panel(fig2_series("Temp"),       "Temperature anomaly (\u00b0C)")
}

save_fig("figures/fig2_baseline_main_3damage_2tax", draw_fig2, w = 1550, h = 1200)
# growth variant: the GDP-level panel is replaced by a GDP-growth panel
save_fig("figures/fig2_baseline_main_3damage_2tax_growth",
		 function() draw_fig2(growth = TRUE), w = 1550, h = 1200)


# =====================================================================================
#  FIGURE 3 -- Baseline forestry-module variables: LOW damage only, 2 carbon-tax schedules
#  Colour = biome ; line type = tax schedule (solid 50 / dashed 300).
#  Rationale: forestry / forest-policy scenarios use the low-damage case only, to avoid
#  interference with the forest policies and allow a clean reading of the dynamics.
# =====================================================================================
biome_series <- function(vtro, vtem, vbor) list(
	list(run = "baseline_low_50",  var = vtro, col = PAL$biome_tro, lty = LTY50),
	list(run = "baseline_low_300", var = vtro, col = PAL$biome_tro, lty = LTY300),
	list(run = "baseline_low_50",  var = vtem, col = PAL$biome_tem, lty = LTY50),
	list(run = "baseline_low_300", var = vtem, col = PAL$biome_tem, lty = LTY300),
	list(run = "baseline_low_50",  var = vbor, col = PAL$biome_bor, lty = LTY50),
	list(run = "baseline_low_300", var = vbor, col = PAL$biome_bor, lty = LTY300)
)

draw_fig3 <- function() {
	par(mfrow = c(2, 3), mar = c(2.6, 4.6, 3.2, 4.1), xpd = FALSE)

	panel(biome_series("F_tro", "F_tem", "F_bor"), "Biomass stocks (billion m3)")
	legend("topright", bty = "n", cex = 0.72,
		 legend = c("Tropical", "Temperate", "Boreal", "tax ~50 (solid)", "tax ~300 (dashed)"),
		 col = c(PAL$biome_tro, PAL$biome_tem, PAL$biome_bor, "black", "black"),
		 lty = c(1, 1, 1, LTY50, LTY300), lwd = 2)

	panel(biome_series("H_tro", "H_tem", "H_bor"), "Total harvest (billion m3)")

	dual_panel("baseline_low_50", "baseline_low_300",
			 "FO_c", PAL$fossil, "BI_m3", PAL$bioen,
			 "Energetic demand", "Fossil (GtC, brown)", "Bioenergy (bn m3, olive)")

	panel(biome_series("HS_tro", "HS_tem", "HS_bor"), "Roundwood harvest (billion m3)")
	panel(biome_series("HB_tro", "HB_tem", "HB_bor"), "Bioenergy harvest (billion m3)")

	emissions_series <- list(
		list(run = "baseline_low_50",  var = "E_ant", col = PAL$em_gross, lty = LTY50),
		list(run = "baseline_low_300", var = "E_ant", col = PAL$em_gross, lty = LTY300),
		list(run = "baseline_low_50",  var = "E",     col = PAL$em_net,   lty = LTY50),
		list(run = "baseline_low_300", var = "E",     col = PAL$em_net,   lty = LTY300)
	)
	emissions_max <- max(unlist(lapply(emissions_series,
		function(s) R[[s$run]][[s$var]])), na.rm = TRUE)
	panel(emissions_series, "Emissions (GtCO2-e)", ylim = c(0, emissions_max))
	axis(side = 4)
	mtext("Gross (brown)", side = 2, line = 3, cex = 0.8)
	mtext("- sequestration (olive)", side = 4, line = 3, cex = 0.75)
}

save_fig("figures/fig3_baseline_forestry_lowdamage_2tax", draw_fig3,
		 w = 1600, h = 1050, ew = 11, eh = 7.2)


# =====================================================================================
#  FIGURES 4 & 5 -- Forest Policy A / B, 3x3 grid (original palette) WITH baselines
#  Forest policy: original per-variable colours ; Baseline: grey50 overlay.
#  Line type = tax schedule (solid 50 / dashed 300) for BOTH policy and baseline.
#  To keep the multi-biome panels legible, only the baseline TROPICAL series is overlaid.
# =====================================================================================
draw_forest_grid <- function(p50, p300, base50, base300, growth = FALSE) {
	tt <- function(nm) R[[nm]]$time + 2016
	g  <- PAL$base_ref
	emissions_max <- max(R[[p50]]$E_ant, R[[p300]]$E_ant,
					 R[[p50]]$E, R[[p300]]$E,
					 R[[base50]]$E_ant, R[[base300]]$E_ant,
					 R[[base50]]$E, R[[base300]]$E, na.rm = TRUE)
	par(mfrow = c(3, 3), mar = c(2.1, 5.1, 4.1, 4.1), xpd = TRUE)

	## 1 -- Output (level) OR Output growth rate
	if (growth) {
		yl <- range(R[[p50]]$g_Y, R[[p300]]$g_Y,
				  R[[base50]]$g_Y, R[[base300]]$g_Y, na.rm = TRUE)
		plot(tt(p50), R[[p50]]$g_Y, type = "l", col = PAL$macro, lty = LTY50, lwd = 2,
			 xlab = "", ylab = "", ylim = yl, main = "Output growth rate (%/yr)")
		lines(tt(p300),    R[[p300]]$g_Y,    col = PAL$macro, lty = LTY300, lwd = 2)
		lines(tt(base50),  R[[base50]]$g_Y,  col = g,         lty = LTY50,  lwd = 2)
		lines(tt(base300), R[[base300]]$g_Y, col = g,         lty = LTY300, lwd = 2)
		lpos <- "topright"
	} else {
		plot(tt(p50), R[[p50]]$Y, type = "l", col = PAL$macro, lty = LTY50, lwd = 2,
			 xlab = "", ylab = "", ylim = c(50, 430), main = "Output (2015 US$ tril.)")
		lines(tt(p300),    R[[p300]]$Y,    col = PAL$macro, lty = LTY300, lwd = 2)
		lines(tt(base50),  R[[base50]]$Y,  col = g,         lty = LTY50,  lwd = 2)
		lines(tt(base300), R[[base300]]$Y, col = g,         lty = LTY300, lwd = 2)
		lpos <- "topleft"
	}
	legend(lpos, bty = "n", cex = 0.68,
		 legend = c("Forest policy", "Baseline", "tax ~50 (solid)", "tax ~300 (dashed)"),
		 col = c(PAL$macro, g, "black", "black"), lty = c(1, 1, LTY50, LTY300), lwd = 2)
	mtext("MACROECONOMIC", side = 2, line = 3, font = 2, cex = 1)

	## 2 -- Labour (employment blue, wage share olive; baseline grey)
	plot(tt(p50), R[[p50]]$lambda, type = "l", col = PAL$emp, lty = LTY50, lwd = 2,
		 xlab = "", ylab = "", ylim = c(0.4, 0.85), main = "Labour")
	lines(tt(p300),    R[[p300]]$lambda, col = PAL$emp,  lty = LTY300, lwd = 2)
	lines(tt(p50),     R[[p50]]$omega,   col = PAL$wage, lty = LTY50,  lwd = 2)
	lines(tt(p300),    R[[p300]]$omega,  col = PAL$wage, lty = LTY300, lwd = 2)
	lines(tt(base50),  R[[base50]]$lambda,  col = g, lty = LTY50,  lwd = 2)
	lines(tt(base300), R[[base300]]$lambda, col = g, lty = LTY300, lwd = 2)
	lines(tt(base50),  R[[base50]]$omega,   col = g, lty = LTY50,  lwd = 2)
	lines(tt(base300), R[[base300]]$omega,  col = g, lty = LTY300, lwd = 2)
	mtext("Employment (blue)",  side = 2, line = 3, cex = 0.7)
	mtext("Wage share (olive)", side = 4, line = 1, cex = 0.7)

	## 3 -- Financial stability (debt blue primary, inflation olive secondary; baseline grey)
	plot(tt(p50), R[[p50]]$d, type = "l", col = PAL$debt, lty = LTY50, lwd = 2,
		 xlab = "", ylab = "", ylim = c(0, 6), main = "Financial stability")
	lines(tt(p300),    R[[p300]]$d,    col = PAL$debt, lty = LTY300, lwd = 2)
	lines(tt(base50),  R[[base50]]$d,  col = g,        lty = LTY50,  lwd = 2)
	lines(tt(base300), R[[base300]]$d, col = g,        lty = LTY300, lwd = 2)
	par(new = TRUE)
	plot(tt(p50), R[[p50]]$i * 100, type = "l", col = PAL$infl, lty = LTY50, lwd = 2,
		 axes = FALSE, xlab = "", ylab = "", ylim = c(-2, 6))
	lines(tt(p300), R[[p300]]$i * 100, col = PAL$infl, lty = LTY300, lwd = 2)
	lines(tt(base50),  R[[base50]]$i * 100, col = g, lty = LTY50,  lwd = 2)
	lines(tt(base300), R[[base300]]$i * 100, col = g, lty = LTY300, lwd = 2)
	axis(side = 4)
	mtext("Debt ratio (blue)",   side = 2, line = 3,   cex = 0.7)
	mtext("Inflation % (olive)", side = 4, line = 2.4, cex = 0.7)

	## 4 -- Biomass stocks (biome palette; baseline tropical grey)
	plot(tt(p50), R[[p50]]$F_tro, type = "l", col = PAL$biome_tro, lty = LTY50, lwd = 2,
		 xlab = "", ylab = "", ylim = c(50, 700), main = "Biomass stocks (billion m3)")
	lines(tt(p300), R[[p300]]$F_tro, col = PAL$biome_tro, lty = LTY300, lwd = 2)
	lines(tt(p50),  R[[p50]]$F_bor,  col = PAL$biome_bor, lty = LTY50,  lwd = 2)
	lines(tt(p300), R[[p300]]$F_bor, col = PAL$biome_bor, lty = LTY300, lwd = 2)
	lines(tt(p50),  R[[p50]]$F_tem,  col = PAL$biome_tem, lty = LTY50,  lwd = 2)
	lines(tt(p300), R[[p300]]$F_tem, col = PAL$biome_tem, lty = LTY300, lwd = 2)
	lines(tt(base50),  R[[base50]]$F_tro,  col = g, lty = LTY50,  lwd = 2)
	lines(tt(base300), R[[base300]]$F_tro, col = g, lty = LTY300, lwd = 2)
	mtext("BIOPHYSICAL", side = 2, line = 3, font = 2, cex = 1)
	legend("topleft", bty = "n", cex = 0.66,
		 legend = c("Tropical", "Temperate", "Boreal", "Baseline (tropical)"),
		 col = c(PAL$biome_tro, PAL$biome_tem, PAL$biome_bor, g), lty = 1, lwd = 2)

	## 5 -- Total harvest (biome palette; baseline tropical grey)
	plot(tt(p50), R[[p50]]$H_tro, type = "l", col = PAL$biome_tro, lty = LTY50, lwd = 2,
		 xlab = "", ylab = "", ylim = c(0, 10), main = "Total harvest (billion m3)")
	lines(tt(p300), R[[p300]]$H_tro, col = PAL$biome_tro, lty = LTY300, lwd = 2)
	lines(tt(p50),  R[[p50]]$H_bor,  col = PAL$biome_bor, lty = LTY50,  lwd = 2)
	lines(tt(p300), R[[p300]]$H_bor, col = PAL$biome_bor, lty = LTY300, lwd = 2)
	lines(tt(p50),  R[[p50]]$H_tem,  col = PAL$biome_tem, lty = LTY50,  lwd = 2)
	lines(tt(p300), R[[p300]]$H_tem, col = PAL$biome_tem, lty = LTY300, lwd = 2)
	lines(tt(base50),  R[[base50]]$H_tro,  col = g, lty = LTY50,  lwd = 2)
	lines(tt(base300), R[[base300]]$H_tro, col = g, lty = LTY300, lwd = 2)

	## 6 -- Energetic demand (fossil orange4 primary, bioenergy olive secondary; baseline grey)
	plot(tt(p50), R[[p50]]$FO_c, type = "l", col = PAL$fossil, lty = LTY50, lwd = 2,
		 xlab = "", ylab = "", ylim = c(0, 75), main = "Energetic demand")
	lines(tt(p300), R[[p300]]$FO_c, col = PAL$fossil, lty = LTY300, lwd = 2)
	lines(tt(p50),  R[[p50]]$BI_m3,  col = PAL$bioen, lty = LTY50,  lwd = 2)
	lines(tt(p300), R[[p300]]$BI_m3, col = PAL$bioen, lty = LTY300, lwd = 2)
	lines(tt(base50),  R[[base50]]$FO_c,  col = g, lty = LTY50,  lwd = 2)
	lines(tt(base300), R[[base300]]$FO_c, col = g, lty = LTY300, lwd = 2)
	lines(tt(base50),  R[[base50]]$BI_m3, col = g, lty = LTY50,  lwd = 2)
	lines(tt(base300), R[[base300]]$BI_m3, col = g, lty = LTY300, lwd = 2)
	mtext("Fossil (GtC, brown)",       side = 2, line = 3, cex = 0.7)
	mtext("Bioenergy (bn m3, olive)",  side = 4, line = 1, cex = 0.7)

	## 7 -- Carbon policy (tax)
	plot(tt(p50), R[[p50]]$tax, type = "l", col = PAL$macro, lty = LTY50, lwd = 2,
		 xlab = "", ylab = "", ylim = c(0, 320), main = "Carbon policy ($/tonC)")
	lines(tt(p300), R[[p300]]$tax, col = PAL$macro, lty = LTY300, lwd = 2)
	mtext("CLIMATE", side = 2, line = 3, font = 2, cex = 1)

	## 8 -- Temperature anomaly
	plot(tt(p50), R[[p50]]$Temp, type = "l", col = PAL$macro, lty = LTY50, lwd = 2,
		 xlab = "", ylab = "", ylim = c(0.5, 4.5), main = "Temperature anomaly (\u00b0C)")
	lines(tt(p300),    R[[p300]]$Temp,    col = PAL$macro, lty = LTY300, lwd = 2)
	lines(tt(base50),  R[[base50]]$Temp,  col = g,         lty = LTY50,  lwd = 2)
	lines(tt(base300), R[[base300]]$Temp, col = g,         lty = LTY300, lwd = 2)

	## 9 -- Emissions (gross orange4, net olive; baseline grey both)
	plot(tt(p50), R[[p50]]$E_ant, type = "l", col = PAL$em_gross, lty = LTY50, lwd = 2,
		 xlab = "", ylab = "", ylim = c(0, emissions_max), main = "Emissions (GtCO2-e)")
	lines(tt(p300), R[[p300]]$E_ant, col = PAL$em_gross, lty = LTY300, lwd = 2)
	lines(tt(p50),  R[[p50]]$E,  col = PAL$em_net, lty = LTY50,  lwd = 2)
	lines(tt(p300), R[[p300]]$E, col = PAL$em_net, lty = LTY300, lwd = 2)
	lines(tt(base50),  R[[base50]]$E_ant, col = g, lty = LTY50,  lwd = 2)
	lines(tt(base300), R[[base300]]$E_ant, col = g, lty = LTY300, lwd = 2)
	lines(tt(base50),  R[[base50]]$E,     col = g, lty = LTY50,  lwd = 2)
	lines(tt(base300), R[[base300]]$E,    col = g, lty = LTY300, lwd = 2)
	axis(side = 4)
	mtext("Gross (brown)", side = 2, line = 3, cex = 0.8)
	mtext("- sequestration (olive)", side = 4, line = 3, cex = 0.75)
}

save_fig("figures/fig4_forestA_with_baseline",
		 function() draw_forest_grid("forestA_50", "forestA_300", "baseline_low_50", "baseline_low_300"),
		 w = 1550, h = 1200)

save_fig("figures/fig5_forestB_with_baseline",
		 function() draw_forest_grid("forestB_50", "forestB_300", "baseline_low_50", "baseline_low_300"),
		 w = 1550, h = 1200)

# growth variants: the GDP-level panel is replaced by a GDP-growth panel
save_fig("figures/fig4_forestA_growth_with_baseline",
		 function() draw_forest_grid("forestA_50", "forestA_300", "baseline_low_50", "baseline_low_300",
								 growth = TRUE),
		 w = 1550, h = 1200)

save_fig("figures/fig5_forestB_growth_with_baseline",
		 function() draw_forest_grid("forestB_50", "forestB_300", "baseline_low_50", "baseline_low_300",
								 growth = TRUE),
		 w = 1550, h = 1200)


# =====================================================================================
#  FIGURE -- Output GROWTH RATE instead of output level, original-article style
#  Colour = scenario (grey = Baseline, blue = Forest A, green = Forest B) ;
#  line type = tax schedule (solid 50 / dashed 300).  Low-damage case, matched to Figs 3-5.
# =====================================================================================
growth_series <- list(
	list(run = "baseline_low_50",  var = "g_Y", col = PAL$base_ref,  lty = LTY50),
	list(run = "baseline_low_300", var = "g_Y", col = PAL$base_ref,  lty = LTY300),
	list(run = "forestA_50",       var = "g_Y", col = PAL$emp,       lty = LTY50),
	list(run = "forestA_300",      var = "g_Y", col = PAL$emp,       lty = LTY300),
	list(run = "forestB_50",       var = "g_Y", col = PAL$biome_bor, lty = LTY50),
	list(run = "forestB_300",      var = "g_Y", col = PAL$biome_bor, lty = LTY300)
)

draw_growth <- function() {
	par(mar = c(4.0, 4.4, 3.2, 1.2), xpd = FALSE)
	panel(growth_series, "Output growth rate (%/yr)", ylab = "%/yr")
	title(xlab = "Year", line = 2.4)
	legend("topright", bty = "n", cex = 0.82,
		 legend = c("Baseline", "Forest policy A", "Forest policy B",
				  "tax ~50 (solid)", "tax ~300 (dashed)"),
		 col = c(PAL$base_ref, PAL$emp, PAL$biome_bor, "black", "black"),
		 lty = c(1, 1, 1, LTY50, LTY300), lwd = 2)
}

save_fig("figures/fig_output_growth", draw_growth, w = 1100, h = 800, ew = 8, eh = 6)


# =====================================================================================
#  SUMMARY -- key 2100 outcomes for every figure scenario
# =====================================================================================
pct_of_first <- function(x) 100 * tail(x, 1) / x[1]

summary_2100 <- data.frame(
	scenario = vapply(figure_scenarios, function(nm) SCENARIOS[[nm]]$label, character(1)),
	Temp_2100     = round(vapply(figure_scenarios, function(nm) tail(R[[nm]]$Temp, 1), numeric(1)), 3),
	Debt_2100     = round(vapply(figure_scenarios, function(nm) tail(R[[nm]]$d, 1), numeric(1)), 3),
	Y_2100        = round(vapply(figure_scenarios, function(nm) tail(R[[nm]]$Y, 1), numeric(1)), 2),
	E_2100        = round(vapply(figure_scenarios, function(nm) tail(R[[nm]]$E, 1), numeric(1)), 2),
	Ftro_pct_2100 = round(vapply(figure_scenarios, function(nm) pct_of_first(R[[nm]]$F_tro), numeric(1)), 1),
	CF_2100       = round(vapply(figure_scenarios, function(nm) tail(R[[nm]]$MC, 1), numeric(1)), 3),
	stringsAsFactors = FALSE, row.names = NULL
)
cat("\n=== 2100 outcomes (figure scenarios) ===\n")
print(summary_2100, row.names = FALSE)
