###############################################################################
# 00_helpers.R  –  Helper functions for MC sensitivity analysis
# 
# This script provides:
#   - Package loading
#   - The model wrapper  run_model_sa()
#   - Output extraction  extract_sa_outputs()
#   - Safe wrapper       safe_run_model()
#   - Utility functions for time-series indexing
#
# Units:
#   MC            : 2015 US$ trillion
#   Forest stocks : billion m3
#   Emissions     : GtCO2-e / yr
#   Temperature   : °C anomaly above pre-industrial
#   Output (Y)    : 2015 US$ trillion
#   Debt ratio (d): dimensionless (D / pY)
#   Inflation (i) : fractional rate per year
###############################################################################

# ---- 1. Packages -----------------------------------------------------------

required_packages <- c(
  "data.table", "ggplot2", "dplyr", "tidyr", "purrr", "patchwork",
  "readr", "scales", "viridis", "lhs", "sensobol", "Rcpp", "parallel"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ---- 2. Project root -------------------------------------------------------

# Detect project root (the folder that contains model.R and SourceCode.R)
find_project_root <- function() {
  # Try working directory first
  candidates <- c(
    getwd(),
    file.path(getwd(), ".."),
    normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."), mustWork = FALSE)
  )
  for (d in candidates) {
    if (file.exists(file.path(d, "model.R")) && file.exists(file.path(d, "SourceCode.R"))) {
      return(normalizePath(d))
    }
  }
  stop("Cannot find project root (model.R + SourceCode.R). Set working directory to the project folder.")
}

# ---- 3. Baseline parameters ------------------------------------------------

BASELINE_PHI1 <- 1.52638276079
BASELINE_PHI2 <- 1.18258619118984

# Time grid (model runs 0->84 with step 0.05, i.e. 2016->2100)
BASE_YEAR   <- 2016
END_YEAR    <- 2100
TIME_STEP   <- 0.05

# ---- 4. Compile model once -------------------------------------------------

#' Compile the C++ RK4 solver (call once per session)
#' @param project_root Path to the folder containing model.R and SourceCode.R
#' @return The sys object from cppMakeSys
compile_model <- function(project_root = NULL) {
  if (is.null(project_root)) project_root <- find_project_root()
  old_wd <- setwd(project_root)
  on.exit(setwd(old_wd))
  # The refactored model.R reads scenario switches from the global env at compile
  # time. The cost-channel SA is nested within the reported forest-policy run
  # (Forest Policy A, low damage, carbon tax ~50, level-form C_F), so we fix those
  # switches here; f1/f2 are perturbed at run time via sys$parms (no recompile).
  scenario_g_p_car  <<- 0.04    # carbon-tax growth (~50 $/tonC by 2100)
  scenario_RD_slope <<- 1/42    # Forest Policy A (DC reaches 2 by 2100)
  scenario_MC_on    <<- 1       # forest-policy cost active
  scenario_MC_mode  <<- 0       # level-form C_F = f1*exp(f2*RE) (as in the original figures)
  scenario_dam_on   <<- 1       # climate damages active
  scenario_f_K      <<- 0       # low-damage case (output only)
  source("SourceCode.R", local = FALSE)
  sys <- cppMakeSys("model.R", reportVars = 3)
  return(sys)
}

# ---- 5. Run model with custom phi1, phi2 -----------------------------------

#' Run the model with given phi1 (f1) and phi2 (f2) values
#'
#' @param sys       Compiled system object from compile_model()
#' @param phi1      Marginal-cost level parameter (default: baseline)
#' @param phi2      Marginal-cost curvature parameter (default: baseline)
#' @param return_ts If TRUE, return full time series; if FALSE, only summary
#' @return A list with elements:
#'   $summary  - named vector of scalar outputs
#'   $ts       - data.frame of annual time series (if return_ts = TRUE)
run_model_sa <- function(sys,
                         phi1 = BASELINE_PHI1,
                         phi2 = BASELINE_PHI2,
                         return_ts = TRUE) {
  
  # Modify parameters
  new_parms <- sys$parms
  new_parms["f1"] <- phi1
  new_parms["f2"] <- phi2
  
  # Run: create a modified sys copy
  sys_run       <- sys
  sys_run$parms <- new_parms
  
  res <- cppRK4(sys_run)
  
  # Extract outputs
  extract_sa_outputs(res, return_ts = return_ts)
}

# ---- 6. Output extraction ---------------------------------------------------

#' Extract summary indicators and (optionally) annual time series
#'
#' @param res       data.frame from cppRK4
#' @param return_ts logical
#' @return list with $summary (named numeric) and optionally $ts (data.frame)
extract_sa_outputs <- function(res, return_ts = TRUE) {
  
  year <- res$time + BASE_YEAR
  
  # Helper: value at nearest year
  val_at <- function(var, target_year) {
    idx <- which.min(abs(year - target_year))
    var[idx]
  }
  
  # Helper: cumulative integral via trapezoidal rule
  cum_trap <- function(var, t_start, t_end) {
    mask <- year >= t_start & year <= t_end
    tt   <- year[mask]
    vv   <- var[mask]
    if (length(tt) < 2) return(NA_real_)
    sum(diff(tt) * (head(vv, -1) + tail(vv, -1)) / 2)
  }
  
  # Sequestration = -(E_ant - E) ; note EF is negative when forests sequester
  sequestration <- -(res$E_ant - res$E)   # positive = net removal by forests
  
  summary_vec <- c(
    MC_2050                  = val_at(res$MC, 2050),
    MC_2100                  = val_at(res$MC, 2100),
    cumulative_MC_2020_2100  = cum_trap(res$MC, 2020, 2100),
    tropical_stock_2050      = val_at(res$F_tro, 2050),
    tropical_stock_2100      = val_at(res$F_tro, 2100),
    sequestration_2100       = val_at(sequestration, 2100),
    net_emissions_2100       = val_at(res$E, 2100),
    cumulative_emissions     = cum_trap(res$E, 2020, 2100),
    output_2100              = val_at(res$Y, 2100),
    private_debt_ratio_2100  = val_at(res$d, 2100),
    min_inflation_post2050   = min(res$i[year >= 2050] * 100, na.rm = TRUE),
    temperature_2100         = val_at(res$Temp, 2100)
  )
  
  out <- list(summary = summary_vec)
  
  if (return_ts) {
    # Subsample to (approximate) annual frequency
    annual_idx <- which(abs(res$time - round(res$time)) < TIME_STEP / 2)
    out$ts <- data.frame(
      year             = year[annual_idx],
      MC               = res$MC[annual_idx],
      tropical_stock   = res$F_tro[annual_idx],
      sequestration    = sequestration[annual_idx],
      net_emissions    = res$E[annual_idx],
      gross_emissions  = res$E_ant[annual_idx],
      output           = res$Y[annual_idx],
      debt_ratio       = res$d[annual_idx],
      inflation_pct    = res$i[annual_idx] * 100,
      temperature      = res$Temp[annual_idx],
      employment       = res$lambda[annual_idx],
      wage_share       = res$omega[annual_idx],
      fossil_demand    = res$FO_c[annual_idx],
      bioenergy_m3     = res$BI_m3[annual_idx]
    )
  }
  
  out
}

# ---- 7. Safe wrapper -------------------------------------------------------

#' Run the model safely, returning NA on failure
#'
#' @param sys   Compiled system
#' @param phi1  f1 parameter
#' @param phi2  f2 parameter
#' @param return_ts logical
#' @param log_file  path to log file (NULL = no logging)
#' @return Same structure as run_model_sa, or list with NAs
safe_run_model <- function(sys, phi1, phi2, return_ts = FALSE, log_file = NULL) {
  
  result <- tryCatch(
    {
      out <- run_model_sa(sys, phi1 = phi1, phi2 = phi2, return_ts = return_ts)
      # Check for non-finite values in summary
      if (any(!is.finite(out$summary))) {
        msg <- sprintf("[WARN] phi1=%.4f phi2=%.4f : non-finite outputs detected", phi1, phi2)
        if (!is.null(log_file)) cat(msg, "\n", file = log_file, append = TRUE)
        # Keep what we have but flag problematic entries
      }
      out
    },
    error = function(e) {
      msg <- sprintf("[ERROR] phi1=%.4f phi2=%.4f : %s", phi1, phi2, e$message)
      if (!is.null(log_file)) cat(msg, "\n", file = log_file, append = TRUE)
      message(msg)
      # Return NA summary
      na_vec <- rep(NA_real_, 12)
      names(na_vec) <- c(
        "MC_2050", "MC_2100", "cumulative_MC_2020_2100",
        "tropical_stock_2050", "tropical_stock_2100",
        "sequestration_2100", "net_emissions_2100", "cumulative_emissions",
        "output_2100", "private_debt_ratio_2100",
        "min_inflation_post2050", "temperature_2100"
      )
      list(summary = na_vec, ts = NULL)
    }
  )
  result
}

# ---- 8. Batch runner -------------------------------------------------------

#' Run a batch of model evaluations (sequential)
#'
#' @param sys       Compiled system object
#' @param param_df  data.frame with columns phi1, phi2
#' @param log_file  path to log file
#' @param progress  logical, print progress
#' @return data.table with param_df columns + output columns
run_batch <- function(sys, param_df, log_file = NULL, progress = TRUE) {
  
  n <- nrow(param_df)
  results <- vector("list", n)
  
  if (!is.null(log_file)) {
    cat(sprintf("=== Batch run: %d evaluations started at %s ===\n",
                n, Sys.time()), file = log_file, append = TRUE)
  }
  
  for (i in seq_len(n)) {
    if (progress && i %% max(1, n %/% 20) == 0) {
      message(sprintf("  run %d / %d (%.0f%%)", i, n, 100 * i / n))
    }
    out <- safe_run_model(sys,
                          phi1 = param_df$phi1[i],
                          phi2 = param_df$phi2[i],
                          return_ts = FALSE,
                          log_file = log_file)
    results[[i]] <- out$summary
  }
  
  res_dt <- data.table::rbindlist(lapply(results, function(x) as.list(x)))
  cbind(data.table::as.data.table(param_df), res_dt)
}

# ---- 9. Utility: nice theme for ggplot2 ------------------------------------

theme_sa <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 2),
      strip.text    = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

cat(">> 00_helpers.R loaded successfully.\n")
