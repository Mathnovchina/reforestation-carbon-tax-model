###############################################################################
# run_refined.R  –  Master script: run all refined SA scripts in sequence
###############################################################################

cat("================================================================\n")
cat("  REFINED SENSITIVITY ANALYSIS  –  Master runner\n")
cat("================================================================\n\n")

t_start <- Sys.time()

cat("STEP 1/3: Sobol SA at N = 2048...\n")
source("sensitivity_analysis/code_refined/01_sobol_refined.R")

cat("\n\nSTEP 2/3: Convergence diagnostics...\n")
source("sensitivity_analysis/code_refined/02_convergence_diagnostics.R")

cat("\n\nSTEP 3/3: Publication figures...\n")
source("sensitivity_analysis/code_refined/03_publication_figures.R")

cat(sprintf("\n\n>> ALL DONE in %.1f minutes.\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
