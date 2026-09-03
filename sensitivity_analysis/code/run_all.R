###############################################################################
# run_all.R  –  Master script: reproduce the full sensitivity analysis
#
# Usage:
#   1. Open R in the project root (the folder with model.R, SourceCode.R)
#   2. source("sensitivity_analysis/code/run_all.R")
#
# This will execute, in order:
#   00_helpers.R           – load packages and define functions
#   01_local_sensitivity.R – OAT perturbation analysis
#   02_sobol_global_sa.R   – global Sobol SA (pilot + full)
#   03_literature_curves.R – literature MC curve overlay
#   04_publication_plots.R – refined publication figures
#
# All outputs are placed under  sensitivity_analysis/
###############################################################################

cat("\n")
cat("##############################################################\n")
cat("##  FULL SENSITIVITY ANALYSIS – MC(RE) = phi1 * exp(phi2*RE) ##\n")
cat("##############################################################\n\n")

t_start <- Sys.time()

# ---- Verify working directory -----------------------------------------------
if (!file.exists("model.R") || !file.exists("SourceCode.R")) {
  stop("Working directory must be the project root containing model.R and SourceCode.R.\n",
       "  Current wd: ", getwd())
}

# ---- Ensure output directories exist ----------------------------------------
dirs <- c(
  "sensitivity_analysis/code",
  "sensitivity_analysis/data",
  "sensitivity_analysis/data/literature_curves",
  "sensitivity_analysis/outputs",
  "sensitivity_analysis/figures",
  "sensitivity_analysis/tables",
  "sensitivity_analysis/logs"
)
for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- Run each step -----------------------------------------------------------

cat("\n===== STEP 1: Local sensitivity analysis =====\n")
source("sensitivity_analysis/code/01_local_sensitivity.R", local = FALSE)

cat("\n===== STEP 2: Global Sobol sensitivity analysis =====\n")
source("sensitivity_analysis/code/02_sobol_global_sa.R", local = FALSE)

cat("\n===== STEP 3: Literature cost curves =====\n")
source("sensitivity_analysis/code/03_literature_curves.R", local = FALSE)

cat("\n===== STEP 4: Publication plots =====\n")
source("sensitivity_analysis/code/04_publication_plots.R", local = FALSE)

# ---- Done --------------------------------------------------------------------

t_end <- Sys.time()
cat(sprintf("\n>> Total runtime: %.1f minutes\n",
            as.numeric(difftime(t_end, t_start, units = "mins"))))

cat("\nAll outputs saved under sensitivity_analysis/\n")
cat("  figures/  – PNG + PDF plots\n")
cat("  tables/   – CSV summary tables\n")
cat("  outputs/  – RDS objects (Sobol results, batch data)\n")
cat("  logs/     – run logs\n")
cat("\nDone.\n")
