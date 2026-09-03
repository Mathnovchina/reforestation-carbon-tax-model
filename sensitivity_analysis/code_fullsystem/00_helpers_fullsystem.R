###############################################################################
# 00_helpers_fullsystem.R  -  Full-system sensitivity-analysis helpers
#
# Extends the forest-policy-cost SA to a whole-system exercise:
#   - compiles the model ONCE under a fixed focal scenario (default forestA_50)
#   - runs the model for an arbitrary NAMED parameter vector (ledger names)
#   - handles the two forest multipliers and the HL->delta_CH override
#   - extracts a clean set of system-level scalar outputs (incl. growth)
#
# Ledger (single source of truth): sensitivity_analysis/parameter_ledger.csv
###############################################################################

# ---- 1. Base helpers (packages, compile_model, find_project_root) -----------

source("sensitivity_analysis/code/00_helpers.R")

# NOTE: Stage-1 screening uses a self-contained Morris (1991) elementary-effects
# implementation (see morris_design()/compute_morris_ee() below) rather than the
# `sensitivity` package, to keep the screening fully transparent and dependency-light.
# Stage-2 variance decomposition uses the `sensobol` package.

# ---- 2. Focal scenario ------------------------------------------------------
# Mirrors SCENARIOS$forestA_50 in run_model.R (Forest Policy A, low damage, tax~50).
# g_mid = 0.04 is the mid carbon-tax growth used in the paper figures.

FOCAL_SCENARIO <- list(
  g_p_car  = 0.04,   # ~58 $/tonC by 2100 ("50")
  RD_slope = 1/42,   # Forest Policy A (DC reaches 2 by 2100)
  MC_on    = 1,      # forest-policy cost active
  MC_mode  = 1,      # incremental C_F (zero at t0)
  dam_on   = 1,      # damages active
  f_K      = 0       # low-damage (output only)
)

#' Assign the scenario switches as globals so loadModel() can evaluate model.R's
#' `g_p_car = scenario_g_p_car` etc. at compile time. Must be called BEFORE
#' compile_model().
set_scenario_globals <- function(scn = FOCAL_SCENARIO) {
  scenario_g_p_car  <<- scn$g_p_car
  scenario_RD_slope <<- scn$RD_slope
  scenario_MC_on    <<- scn$MC_on
  scenario_MC_mode  <<- scn$MC_mode
  scenario_dam_on   <<- scn$dam_on
  scenario_f_K      <<- scn$f_K
  invisible(scn)
}

#' Compile the model under a fixed focal scenario (call once per session).
compile_model_focal <- function(scn = FOCAL_SCENARIO, project_root = NULL) {
  if (is.null(project_root)) project_root <- find_project_root()
  set_scenario_globals(scn)
  old_wd <- setwd(project_root)
  on.exit(setwd(old_wd))
  source("SourceCode.R", local = FALSE)
  cppMakeSys("model.R", reportVars = 3)
}

# ---- 3. Ledger --------------------------------------------------------------

load_ledger <- function(path = "sensitivity_analysis/parameter_ledger.csv") {
  led <- read.csv(path, stringsAsFactors = FALSE)
  stopifnot(all(c("param", "baseline", "lower", "upper") %in% names(led)))
  led
}

# Underlying model targets for the two forest multipliers.
FOREST_PHI  <- c("phi_bor", "phi_tem", "phi_tro")   # parameters
FOREST_FMAX <- c("Fmax_bor", "Fmax_tem")            # parameters
FMAX_STATE  <- "Fmax_tro"                           # initial state (y0)

# ---- 4. Generalized runner --------------------------------------------------

#' Run the model for a NAMED parameter vector using ledger names.
#'
#' Recognised special names:
#'   m_phiF  -> multiplies phi_bor/phi_tem/phi_tro
#'   m_Fmax  -> multiplies Fmax_bor/Fmax_tem (parms) AND Fmax_tro (y0)
#' Names present in sys$parms are written to parms; names present only in sys$y0
#' (e.g. g_sigma, an integrated state) are written to the initial state.
#' NB: delta_CH is the live wood-product decay control (the HL intermediate is inert).
run_model_named <- function(sys, pvec, return_ts = FALSE) {
  parms <- sys$parms
  y0    <- sys$y0

  for (nm in names(pvec)) {
    v <- pvec[[nm]]
    if (nm == "m_phiF") {
      parms[FOREST_PHI] <- sys$parms[FOREST_PHI] * v
    } else if (nm == "m_Fmax") {
      parms[FOREST_FMAX] <- sys$parms[FOREST_FMAX] * v
      y0[FMAX_STATE]     <- sys$y0[FMAX_STATE] * v
    } else if (nm %in% names(parms)) {
      parms[nm] <- v
    } else if (nm %in% names(y0)) {
      # State initial value (e.g. g_sigma is integrated, its baseline lives in y0).
      y0[nm] <- v
    } else {
      stop("Unknown parameter '", nm, "' not in sys$parms or sys$y0")
    }
  }

  sys_run       <- sys
  sys_run$parms <- parms
  sys_run$y0    <- y0

  res <- cppRK4(sys_run)
  extract_fs_outputs(res, return_ts = return_ts)
}

# ---- 5. System-level output extraction --------------------------------------

extract_fs_outputs <- function(res, return_ts = FALSE) {
  year <- res$time + BASE_YEAR

  val_at <- function(var, ty) var[which.min(abs(year - ty))]
  cum_trap <- function(var, t0, t1) {
    m <- year >= t0 & year <= t1
    tt <- year[m]; vv <- var[m]
    if (length(tt) < 2) return(NA_real_)
    sum(diff(tt) * (head(vv, -1) + tail(vv, -1)) / 2)
  }

  sequestration <- -(res$E_ant - res$E)          # positive = net forest removal
  total_stock   <- res$F_bor + res$F_tem + res$F_tro
  infl_pct      <- res$i * 100

  # Output growth 2050-2100 as compound annual growth rate (%)
  Y50 <- val_at(res$Y, 2050); Y100 <- val_at(res$Y, 2100)
  growth_50_100 <- if (is.finite(Y50) && is.finite(Y100) && Y50 > 0 && Y100 > 0) {
    ((Y100 / Y50)^(1 / 50) - 1) * 100
  } else NA_real_

  mask_50_100 <- year >= 2050 & year <= 2100
  avg_infl_50_100 <- if (any(mask_50_100)) mean(infl_pct[mask_50_100], na.rm = TRUE) else NA_real_

  summary_vec <- c(
    tropical_stock_2100      = val_at(res$F_tro, 2100),
    total_forest_stock_2100  = val_at(total_stock, 2100),
    sequestration_2100       = val_at(sequestration, 2100),
    cumulative_sequestration = cum_trap(sequestration, 2020, 2100),
    net_emissions_2100       = val_at(res$E, 2100),
    cumulative_emissions     = cum_trap(res$E, 2020, 2100),
    temperature_2100         = val_at(res$Temp, 2100),
    output_growth_2050_2100  = growth_50_100,
    private_debt_ratio_2100  = val_at(res$d, 2100),
    avg_inflation_2050_2100  = avg_infl_50_100,
    CF_2100                  = val_at(res$MC, 2100),
    cumulative_CF            = cum_trap(res$MC, 2020, 2100)
  )

  out <- list(summary = summary_vec)
  if (return_ts) {
    idx <- which(abs(res$time - round(res$time)) < TIME_STEP / 2)
    out$ts <- data.frame(
      year = year[idx], tropical_stock = res$F_tro[idx],
      total_stock = total_stock[idx], sequestration = sequestration[idx],
      net_emissions = res$E[idx], temperature = res$Temp[idx],
      output = res$Y[idx], debt_ratio = res$d[idx],
      inflation_pct = infl_pct[idx], CF = res$MC[idx]
    )
  }
  out
}

FS_OUTPUT_NAMES <- c(
  "tropical_stock_2100", "total_forest_stock_2100", "sequestration_2100",
  "cumulative_sequestration", "net_emissions_2100", "cumulative_emissions",
  "temperature_2100", "output_growth_2050_2100", "private_debt_ratio_2100",
  "avg_inflation_2050_2100", "CF_2100", "cumulative_CF"
)

# Display labels aligned with the paper figure/panel titles (run_model.R).
# The time-window suffix marks the aggregation: "(2100)" = end-of-century level,
# "(avg. 2050-2100)" = second-half mean of a volatile rate, "(cum. 2020-2100)" =
# time-integrated flow.
FS_OUTPUT_LABELS <- c(
  tropical_stock_2100      = "Tropical biomass stock (2100)",
  total_forest_stock_2100  = "Total biomass stock (2100)",
  sequestration_2100       = "Sequestration flow (2100)",
  cumulative_sequestration = "Cumulative sequestration (2020-2100)",
  net_emissions_2100       = "Emissions (2100)",
  cumulative_emissions     = "Cumulative emissions (2020-2100)",
  temperature_2100         = "Temperature anomaly (2100)",
  output_growth_2050_2100  = "Output growth (avg. 2050-2100)",
  private_debt_ratio_2100  = "Private debt ratio (2100)",
  avg_inflation_2050_2100  = "Inflation (avg. 2050-2100)",
  CF_2100                  = "Forest-policy cost C_F (2100)",
  cumulative_CF            = "Cumulative forest-policy cost (2020-2100)"
)

#' Map internal output name(s) to paper-aligned display label(s).
fs_label <- function(x) ifelse(x %in% names(FS_OUTPUT_LABELS), FS_OUTPUT_LABELS[x], x)

# Physical units of each raw system output (plotmath strings). Morris mu* and
# sigma inherit the units of the underlying output, so these are shown beneath
# every panel title. Units match the magnitudes reported in the paper:
# stocks in billion m^3, sequestration in GtC/yr, net emissions in GtCO2e/yr,
# temperature in degrees C, growth/inflation in % per year, the private-debt
# ratio is dimensionless, and the forest-policy cost C_F is in 10^12 USD (the
# model's output unit; the paper additionally reports it as a share of GDP).
FS_OUTPUT_UNITS <- c(
  tropical_stock_2100      = "billion~m^3",
  total_forest_stock_2100  = "billion~m^3",
  sequestration_2100       = "GtC~yr^-1",
  cumulative_sequestration = "GtC",
  net_emissions_2100       = "GtCO[2]*e~yr^-1",
  cumulative_emissions     = "GtCO[2]*e",
  temperature_2100         = "degree*C",
  output_growth_2050_2100  = "'%'~yr^-1",
  private_debt_ratio_2100  = "ratio",
  avg_inflation_2050_2100  = "'%'~yr^-1",
  CF_2100                  = "10^12~USD~yr^-1",
  cumulative_CF            = "10^12~USD"
)

#' Plotmath unit string for an output name ("" if unknown).
fs_unit <- function(x) if (x %in% names(FS_OUTPUT_UNITS)) FS_OUTPUT_UNITS[[x]] else ""

#' Two-line panel title: paper label on top, bracketed physical unit below.
fs_title <- function(oname) {
  lab <- unname(fs_label(oname))
  u   <- fs_unit(oname)
  if (is.null(u) || u == "") return(lab)
  bquote(atop(.(lab), .(str2lang(paste0('"["*', u, '*"]"')))))
}

# Parameter display labels — plotmath strings matching the paper's symbol table.
# phi_1 = Phillips slope (phi_lambda); f1/f2 use varphi to distinguish from phi_lambda.
PARAM_LABELS <- c(
  phi_1        = "phi[lambda]",
  alpha        = "alpha",
  k_pi         = "kappa[pi]",
  div_pi       = "Delta[pi]",
  eta          = "eta",
  g_sigma      = "g[sigma]",
  delta_g_sigma = "delta[g[sigma]]",
  theta        = "theta",
  g_pbs        = "delta[pBS]",
  beta_car     = "beta",
  m_phiF       = "m[psi]",
  m_Fmax       = "m[F^max]",
  delta_E_tro  = "delta[E^Land]",
  delta_CH     = "delta[CH]",
  S            = "S",
  CUP_pind     = "C[pind]^UP",
  gamma_ast    = "gamma^'*'",
  f1           = "varphi[1]",
  f2           = "varphi[2]"
)

#' Write a figure to both PDF (preview) and EPS (publication) — mirrors save_fig() in run_model.R.
save_sa_fig <- function(path, drawfun, w = 9, h = 6) {
  pdf(paste0(path, ".pdf"), width = w, height = h); drawfun(); dev.off()
  postscript(paste0(path, ".eps"), onefile = FALSE, horizontal = FALSE,
             paper = "special", width = w, height = h); drawfun(); dev.off()
  message("written: ", path, ".pdf  +  .eps")
}

#' Return an expression vector of paper symbols for a vector of code-name parameters.
param_expr <- function(pnames) {
  lbls <- ifelse(pnames %in% names(PARAM_LABELS), PARAM_LABELS[pnames], pnames)
  do.call("expression", lapply(lbls, function(s) parse(text = s)[[1]]))
}

# ---- 6. Safe wrapper --------------------------------------------------------

#' Run a single named draw, returning an all-NA output row on failure or on any
#' non-finite output (flagged in log_file).
safe_run_named <- function(sys, pvec, log_file = NULL) {
  tryCatch({
    out <- run_model_named(sys, pvec, return_ts = FALSE)
    s <- out$summary[FS_OUTPUT_NAMES]
    if (any(!is.finite(s))) {
      msg <- sprintf("[WARN] non-finite outputs: %s",
                     paste(sprintf("%s=%.4g", names(pvec), unlist(pvec)), collapse = " "))
      if (!is.null(log_file)) cat(msg, "\n", file = log_file, append = TRUE)
    }
    s
  }, error = function(e) {
    msg <- sprintf("[ERROR] %s | %s",
                   paste(sprintf("%s=%.4g", names(pvec), unlist(pvec)), collapse = " "),
                   e$message)
    if (!is.null(log_file)) cat(msg, "\n", file = log_file, append = TRUE)
    setNames(rep(NA_real_, length(FS_OUTPUT_NAMES)), FS_OUTPUT_NAMES)
  })
}

#' Run a whole design matrix X (rows = draws, cols named by ledger params).
#' Returns a matrix (rows x outputs).
run_design <- function(sys, X, log_file = NULL, verbose = TRUE) {
  n <- nrow(X)
  Y <- matrix(NA_real_, nrow = n, ncol = length(FS_OUTPUT_NAMES),
              dimnames = list(NULL, FS_OUTPUT_NAMES))
  for (i in seq_len(n)) {
    pvec <- as.list(X[i, ]); names(pvec) <- colnames(X)
    Y[i, ] <- safe_run_named(sys, pvec, log_file = log_file)
    if (verbose && (i %% 50 == 0 || i == n)) cat(sprintf("  run %d / %d\n", i, n))
  }
  Y
}

# ---- 7. Morris (1991) elementary-effects design -----------------------------

#' Build r Morris OAT trajectories over k factors on a p-level grid.
#'
#' Each trajectory is a "winding stairs" walk: from a random feasible base point,
#' the factors are perturbed one at a time (random order, random sign) by a fixed
#' step delta = grid_jump/(p-1) on the unit hypercube. Consecutive design points
#' therefore differ in exactly one coordinate. Points are returned BOTH in the
#' normalized [0,1] space and mapped to real units via [lower, upper].
#'
#' @return list(X_real, X_norm, traj_rows, traj_order, traj_sign, delta, pnames)
morris_design <- function(pnames, lower, upper, r = 20, levels = 6, grid_jump = 3) {
  k     <- length(pnames)
  p     <- levels
  delta <- grid_jump / (p - 1)
  grid  <- seq(0, 1, length.out = p)

  X_norm     <- matrix(NA_real_, nrow = r * (k + 1), ncol = k, dimnames = list(NULL, pnames))
  traj_rows  <- vector("list", r)   # row indices in X for each trajectory (length k+1)
  traj_order <- vector("list", r)   # order in which factors are perturbed (length k)
  traj_sign  <- vector("list", r)   # +1/-1 step direction per perturbed factor (length k)

  row <- 0L
  for (m in seq_len(r)) {
    ord  <- sample.int(k)                        # perturbation order
    sgn  <- sample(c(-1, 1), k, replace = TRUE)  # per-factor step direction
    # Feasible base: if we step +delta the base must allow room, and vice-versa.
    base <- numeric(k)
    for (i in seq_len(k)) {
      feas <- if (sgn[i] > 0) grid[grid <= 1 - delta + 1e-9] else grid[grid >= delta - 1e-9]
      base[i] <- sample(feas, 1)
    }

    pts  <- matrix(NA_real_, nrow = k + 1, ncol = k)
    cur  <- base
    pts[1, ] <- cur
    for (j in seq_len(k)) {
      fi        <- ord[j]
      cur[fi]   <- cur[fi] + sgn[fi] * delta
      pts[j + 1, ] <- cur
    }

    idx <- (row + 1L):(row + k + 1L)
    X_norm[idx, ] <- pts
    traj_rows[[m]]  <- idx
    traj_order[[m]] <- ord
    traj_sign[[m]]  <- sgn
    row <- row + k + 1L
  }

  # Map to real units.
  span   <- upper - lower
  X_real <- sweep(sweep(X_norm, 2, span, `*`), 2, lower, `+`)
  colnames(X_real) <- pnames

  list(X_real = X_real, X_norm = X_norm, traj_rows = traj_rows,
       traj_order = traj_order, traj_sign = traj_sign,
       delta = delta, pnames = pnames, r = r)
}

#' Compute elementary effects and Morris measures (mu, mu*, sigma) for one output
#' vector y (aligned to the design's X rows).
compute_morris_ee <- function(design, y) {
  pnames <- design$pnames
  k      <- length(pnames)
  r      <- design$r
  delta  <- design$delta

  EE <- matrix(NA_real_, nrow = r, ncol = k, dimnames = list(NULL, pnames))
  for (m in seq_len(r)) {
    rows <- design$traj_rows[[m]]
    ord  <- design$traj_order[[m]]
    sgn  <- design$traj_sign[[m]]
    yv   <- y[rows]
    for (j in seq_len(k)) {
      fi <- ord[j]
      EE[m, fi] <- (yv[j + 1] - yv[j]) / (sgn[fi] * delta)
    }
  }
  data.frame(
    parameter = pnames,
    mu        = apply(EE, 2, function(e) mean(e,      na.rm = TRUE)),
    mu_star   = apply(EE, 2, function(e) mean(abs(e), na.rm = TRUE)),
    sigma     = apply(EE, 2, function(e) sd(e,        na.rm = TRUE)),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

cat(">> Full-system SA helpers loaded (focal scenario:",
    "forestA_50; 19-param ledger; 12 outputs; manual Morris).\n")
