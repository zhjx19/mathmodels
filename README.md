
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mathmodels:  A Powerful R Toolkit for Mathematical Modeling

<!-- badges: start -->
<!-- badges: end -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/mathmodels)](https://cran.r-project.org/package=mathmodels)
[![Online Manual](https://img.shields.io/badge/Online-Manual-blue)](https://zhjx19.github.io/mathmodels-book/)

An R package providing a versatile toolkit for mathematical modeling, developed as a companion to the book *Mathematical Modeling: Algorithms and Programming Implementation* (China Machine Press). It focuses on implementing rigorous algorithms in a user-friendly manner.

**Current Version (0.0.8)** adds ODE-based differential equation modeling, epidemic compartment models (SI/SIS/SIR/SEIR), Lotka–Volterra predator–prey dynamics, comprehensive epidemic visualization and metrics, and Markov chain prediction — on top of the established evaluation algorithms (AHP, Entropy, CRITIC, PCA, TOPSIS, Fuzzy, RSR, DEA), inequality measures (Gini, Theil), and grey prediction models (GM(1,1), GM(1,N), Verhulst).

## Key Features

*   **Comprehensive Modeling**: Covers evaluation methods (AHP, CRITIC, Entropy, TOPSIS, Fuzzy, RSR, DEA), differential equation models (SI/SIS/SIR/SEIR, Lotka–Volterra, Malthus, Logistic), Markov chain prediction, grey prediction (GM(1,1), GM(1,N), Verhulst), and inequality measures (Gini, Theil).
*   **ODE Solver & Epidemic Toolkit**: String-formula based `ode_solver()` for arbitrary ODE systems, ready-to-use epidemic compartment models, plus dedicated visualization and metrics functions.
*   **Tidyverse Integration**: Seamlessly works with `|>` and `tidyverse` tools for smooth data manipulation and batch processing.

## Installation

You can install the latest development version of `mathmodels` directly
from [GitHub](https://github.com/zhjx19/mathmodels) use:

``` r
remotes::install_github("zhjx19/mathmodels")
```

Or download to current path, unzip and use:

```r
# install.packages("lpSolveAPI")   # If necessary, install dependency packages first
install.packages("mathmodels-main", repos=NULL, type="source")
```

## Getting Started

```r
library(mathmodels)

# --- Epidemic compartment modeling ---
library(ggplot2)

result = model_sir(
  init   = c(S = 990, I = 10, R = 0),
  params = c(beta = 0.3, gamma = 0.1),
  times  = seq(0, 100, by = 0.5)
)

# Visualize
plot_compartments(result, compartments = c("S", "I", "R"))

# Compute epidemic metrics
metrics = epidemic_metrics(result, params = c(beta = 0.3, gamma = 0.1), N = 1000)
metrics$summary$R0         # basic reproduction number
metrics$summary$peak_time  # time of peak infection

# --- AHP (Analytic Hierarchy Process) ---
A = matrix(c(1,   1/2, 4, 3,   3,
             2,   1,   7, 5,   5,
             1/4, 1/7, 1, 1/2, 1/3,
             1/3, 1/5, 2, 1,   1,
             1/3, 1/5, 3, 1,   1), byrow = TRUE, nrow = 5)
AHP(A)
```

## Learning More

For detailed documentation, tutorials, and in-depth examples on using the `mathmodels` package, please refer to the comprehensive online manual:

📘 **[mathmodels Package Manual - Simplifying Mathematical Modeling (Online Book, In Chinese)](https://zhjx19.github.io/mathmodels-book/)**

This online book is the definitive guide to the package's functionalities. Currently implemented modules include:

- **Differential equation models**: Malthus, Logistic, SI, SIS, SIR, SEIR, Lotka–Volterra with `ode_solver()` and `model_*()` functions
- **Epidemic visualization**: `plot_compartments()`, `plot_incidence()`, `plot_Rt_estimate()`, etc.
- **Epidemic metrics**: `epidemic_metrics()` for R0, peak, attack rate, and trajectory analysis
- **Markov chain prediction**: `markov_chain()` and `GM11_markov()`
- Indicator data preprocessing
- AHP, Entropy weighting, CRITIC, PCA weighting
- Weight combination techniques
- TOPSIS, Grey Relational Analysis (GRA)
- Rank Sum Ratio (RSR), Fuzzy Comprehensive Evaluation (FCE)
- Data Envelopment Analysis (CCR/BCC/SBM, Malmquist)
- Inequality Measures (Gini, Theil Index)
- Regional Economics (LQ/HHI/EG Index)
- Coupling coordination degree and obstacle degree







