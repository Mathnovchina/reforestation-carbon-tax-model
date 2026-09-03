###############################################################################
# 00_helpers_refined.R  –  Extended helpers for publication-quality SA
#
# EXTENDS the original 00_helpers.R with:
#   - Additional inflation outputs (inflation_2100, avg_inflation_2050_2100,
#     deflation_episodes, first_deflation_year)
#   - Identical compile_model / run infrastructure (sources the original)
#   - Refined output extraction
###############################################################################

# ---- 1. Source the original helpers (packages, compile_model, etc.) ---------

source("sensitivity_analysis/code/00_helpers.R")

# ---- 2. Override extract_sa_outputs with expanded version -------------------

extract_sa_outputs <- function(res, return_ts = TRUE) {
  
  year <- res$time + BASE_YEAR
  
  val_at <- function(var, target_year) {
    idx <- which.min(abs(year - target_year))
    var[idx]
  }
  
  cum_trap <- function(var, t_start, t_end) {
    mask <- year >= t_start & year <= t_end
    tt   <- year[mask]
    vv   <- var[mask]
    if (length(tt) < 2) return(NA_real_)
    sum(diff(tt) * (head(vv, -1) + tail(vv, -1)) / 2)
  }
  
  sequestration <- -(res$E_ant - res$E)
  
  # Inflation in percent
  infl_pct <- res$i * 100
  
  # Average inflation 2050-2100
  mask_50_100 <- year >= 2050 & year <= 2100
  avg_infl_50_100 <- if (sum(mask_50_100) > 0) {
    mean(infl_pct[mask_50_100], na.rm = TRUE)
  } else NA_real_
  
  # Deflation episodes: contiguous runs where inflation < 0, post-2030
  mask_post2030 <- year >= 2030
  infl_post2030 <- infl_pct[mask_post2030]
  year_post2030 <- year[mask_post2030]
  is_deflation  <- infl_post2030 < 0
  
  # Count deflation episodes (transitions from non-deflation to deflation)
  n_defl_episodes <- 0L
  first_defl_year <- NA_real_
  if (any(is_deflation)) {
    runs <- rle(is_deflation)
    n_defl_episodes <- sum(runs$values)
    # First deflation year
    first_defl_idx <- which(is_deflation)[1]
    first_defl_year <- year_post2030[first_defl_idx]
  }
  
  summary_vec <- c(
    # --- Original outputs ---
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
    temperature_2100         = val_at(res$Temp, 2100),
    # --- Extended inflation outputs ---
    min_inflation_post2050   = min(infl_pct[year >= 2050], na.rm = TRUE),
    inflation_2100           = val_at(infl_pct, 2100),
    avg_inflation_2050_2100  = avg_infl_50_100,
    deflation_episodes       = as.numeric(n_defl_episodes),
    first_deflation_year     = first_defl_year
  )
  
  out <- list(summary = summary_vec)
  
  if (return_ts) {
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
      inflation_pct    = infl_pct[annual_idx],
      temperature      = res$Temp[annual_idx],
      employment       = res$lambda[annual_idx],
      wage_share       = res$omega[annual_idx],
      fossil_demand    = res$FO_c[annual_idx],
      bioenergy_m3     = res$BI_m3[annual_idx]
    )
  }
  
  out
}

# ---- 3. Override safe_run_model to use new output count ---------------------

safe_run_model <- function(sys, phi1, phi2, return_ts = FALSE, log_file = NULL) {
  result <- tryCatch(
    {
      out <- run_model_sa(sys, phi1 = phi1, phi2 = phi2, return_ts = return_ts)
      if (any(!is.finite(out$summary[!names(out$summary) %in%
                                      c("first_deflation_year")]))) {
        msg <- sprintf("[WARN] phi1=%.4f phi2=%.4f : non-finite outputs detected",
                       phi1, phi2)
        if (!is.null(log_file)) cat(msg, "\n", file = log_file, append = TRUE)
      }
      out
    },
    error = function(e) {
      msg <- sprintf("[ERROR] phi1=%.4f phi2=%.4f : %s", phi1, phi2, e$message)
      if (!is.null(log_file)) cat(msg, "\n", file = log_file, append = TRUE)
      message(msg)
      na_vec <- rep(NA_real_, 16)
      names(na_vec) <- c(
        "MC_2050", "MC_2100", "cumulative_MC_2020_2100",
        "tropical_stock_2050", "tropical_stock_2100",
        "sequestration_2100", "net_emissions_2100", "cumulative_emissions",
        "output_2100", "private_debt_ratio_2100", "temperature_2100",
        "min_inflation_post2050", "inflation_2100", "avg_inflation_2050_2100",
        "deflation_episodes", "first_deflation_year"
      )
      list(summary = na_vec, ts = NULL)
    }
  )
  result
}

# ---- 4. Tidy Sobol helper (handles sensobol 1.1.6 CI issue) ----------------

tidy_sobol <- function(sobol_results, conf_level = 0.95) {
  z <- qnorm(1 - (1 - conf_level) / 2)
  rows <- list()
  for (oname in names(sobol_results)) {
    sob <- sobol_results[[oname]]
    idx <- sob$results
    for (j in seq_len(nrow(idx))) {
      lo <- idx$`low.ci`[j]
      hi <- idx$`high.ci`[j]
      if (is.na(lo) && !is.na(idx$std.error[j])) {
        lo <- idx$original[j] - z * idx$std.error[j]
        hi <- idx$original[j] + z * idx$std.error[j]
      }
      rows[[length(rows) + 1]] <- data.frame(
        output     = oname,
        parameters = idx$parameters[j],
        sensitivity= idx$sensitivity[j],
        original   = idx$original[j],
        std_error  = idx$std.error[j],
        low_ci     = lo,
        high_ci    = hi,
        stringsAsFactors = FALSE
      )
    }
  }
  data.table::rbindlist(rows)
}

cat(">> Refined helpers loaded (inflation outputs added).\n")
