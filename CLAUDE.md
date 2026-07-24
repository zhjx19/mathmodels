# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment

- R 4.6.1: `D:\R-4.6.1\bin\Rscript.exe`
- Pandoc: `D:\Program Files\RStudio\resources\app\bin\quarto\bin\tools\pandoc.exe`
- Quarto: `D:\Program Files\RStudio\resources\app\bin\quarto\bin\quarto.exe`

## What this project is

`mathmodels` is an R package of mathematical modeling algorithms. The package is organized by modeling family:

- Evaluation pipeline: preprocessing → weights → TOPSIS / RSR / grey relational analysis / fuzzy evaluation / DEA / inequality / regional economics / system evaluation.
- Prediction pipeline: grey prediction and Markov / Grey-Markov models in `R/pred_grey.R` and `R/pred_markov.R`, plus the time-series toolkit in `R/pred_ts.R`.
- Dynamic models: generic ODE solving, population and epidemic models, epidemic plots, and epidemic metrics in `R/diffEq_models.R`.
- Data helpers: `R/read_nbs.R` for National Bureau of Statistics XLS imports and `data/water_quality.rda` as the bundled example dataset.

The package’s user-facing reference is split across `README.md`, `vignettes/welcome.Rmd`, `man/`, and the pkgdown site in `docs/`. The online manual linked from the README is the broadest feature reference.

Generated / derived outputs live in `NAMESPACE`, `man/`, and `docs/`. Change the roxygen source in `R/`, then regenerate those artifacts.

## Common commands

```bash
# Load the package for interactive work
Rscript -e "devtools::load_all()"

# Run the full test suite
Rscript -e "devtools::test()"

# Run tests for a subset of files by name pattern
Rscript -e "devtools::test(filter = '^ts')"

# Run tests for the active source file
Rscript -e "devtools::test_active_file('R/pred_ts.R')"

# Run one named test in the active source file
Rscript -e "devtools::test_active_file('R/pred_ts.R', desc = 'specific test name')"

# Regenerate roxygen documentation, NAMESPACE, and man/
Rscript -e "devtools::document()"

# Build a source package tarball
Rscript -e "devtools::build()"

# Run the package check
Rscript -e "devtools::check()"

# Check pkgdown reference coverage
Rscript -e "pkgdown::check_pkgdown()"

# Rebuild the committed docs/ site
Rscript -e "pkgdown::build_site()"

# Format R code
air format .
```

## Working conventions

- Use the base pipe `|>` rather than `%>%`.
- Use `\() ...` for single-line anonymous functions; use `function() { ... }` for multi-line ones.
- Keep exported functions roxygen-documented, then rerun `devtools::document()`.
- Add user-facing changes to `NEWS.md`.
- Keep tests close to the feature they cover, following `tests/testthat/test-<topic>.R`.
- `tests/testthat.R` only calls `test_check("mathmodels")`.
- The GitHub Pages workflow publishes the committed `docs/` directory directly, so rebuilt docs matter when public reference pages change.
- `README.md` and `vignettes/welcome.Rmd` are public-facing; keep them aligned with the exported surface and docs.

## Architecture notes

- `R/eval_preprocess.R`, `R/eval_weights.R`, `R/eval_algorithms.R`, `R/eval_fuzzy.R`, `R/eval_dea.R`, `R/eval_inequality.R`, `R/eval_regional.R`, and `R/eval_system.R` form one evaluation stack. Shared preprocessing feeds the weighting and evaluation methods.
- `R/pred_grey.R` and `R/pred_markov.R` hold the older prediction families. `R/pred_ts.R` is the largest predictive module and combines transformation, diagnostics, modeling, forecasting, and plotting helpers in one file.
- `R/diffEq_models.R` bundles the reusable ODE solver with the built-in population and epidemic models, plus the epidemic visualizations and metrics.
- Keep changes as small as possible across families; these modules are grouped by model family intentionally, and interface stability matters more than reshaping the layout.
