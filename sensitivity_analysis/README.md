# Sensitivity analysis

Three complementary sensitivity-analysis suites for the coupled
economy–climate–forest model. All scripts assume the **working directory is the
repository root** (the folder containing `model.R` and `SourceCode.R`) and use
paths relative to it.

---

## Suites

### 1. Cost-channel (`code/`)

Sensitivity of the model to the two parameters of the forest-policy cost
schedule. Includes one-at-a-time perturbations, a global Sobol decomposition,
and an overlay of the fitted cost curve against literature estimates.

```r
source("sensitivity_analysis/code/run_all.R")
```

Outputs are written under `sensitivity_analysis/`:
`figures/`, `tables/`, `outputs/`, `logs/` (created automatically).

### 2. Refined cost-channel (`code_refined/`)

The same cost-channel exercise at larger sample sizes, with Sobol-index
convergence diagnostics and publication-quality figures.

```r
source("sensitivity_analysis/code_refined/run_refined.R")
```

Outputs: `figures_final/`, `tables_final/`, `outputs_refined/`
(shipped as empty folders; populated when the scripts run).

### 3. Full-system global (`code_fullsystem/`)

A whole-system global sensitivity analysis over 19 structural and behavioural
parameters, drawn from `parameter_ledger.csv`:

1. `01_morris_screening.R`   – Morris elementary-effects screening
2. `02_sobol_reduced.R`      – variance-based (Sobol) decomposition on the retained parameters
3. `03_viability_conditional.R` – model viability and, conditional on viability, output drivers
4. `04_policy_robustness.R`  – paired sign-robustness of the forest-policy conclusions

```r
source("sensitivity_analysis/code_fullsystem/01_morris_screening.R")
source("sensitivity_analysis/code_fullsystem/02_sobol_reduced.R")
source("sensitivity_analysis/code_fullsystem/03_viability_conditional.R")
source("sensitivity_analysis/code_fullsystem/04_policy_robustness.R")
```

Outputs: `figures_fullsystem/`, `tables_fullsystem/`, `outputs_fullsystem/`,
`logs_fullsystem/` (created automatically). Stages 3–4 reuse the design
evaluated in stage 2, so run the scripts in order.

---

## Recommended run order

The refined publication-figure step reuses a figure produced by the full-system
suite. For a complete reproduction from scratch, run the suites in this order:

```r
source("sensitivity_analysis/code/run_all.R")
source("sensitivity_analysis/code_fullsystem/01_morris_screening.R")
source("sensitivity_analysis/code_fullsystem/02_sobol_reduced.R")
source("sensitivity_analysis/code_fullsystem/03_viability_conditional.R")
source("sensitivity_analysis/code_fullsystem/04_policy_robustness.R")
source("sensitivity_analysis/code_refined/run_refined.R")
```

---

## Notes

- `parameter_ledger.csv` is the single source of truth for the full-system SA:
  parameter names, baselines, ranges, distributions, and module assignments.
- `data/literature_curves/` holds the input data for the literature cost-curve
  overlay in `code/03_literature_curves.R`.
- Global sensitivity analysis re-runs the model thousands of times and can take
  a substantial amount of wall-clock time; sample sizes are set at the top of
  each script.
