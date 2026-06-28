# mathmodels: A Powerful R Toolkit for Mathematical Modeling

An R package providing a versatile toolkit for mathematical modeling,
developed as a companion to the book *Mathematical Modeling: Algorithms
and Programming Implementation* (China Machine Press). It focuses on
implementing rigorous algorithms in a user-friendly manner.

**Current Version (0.0.9)** adds a complete time series toolkit (SARIMA,
ETS, GARCH, with integrated transformation and visualization) and
overhauls epidemic modeling with streamlined visualization
(`epi_plots.R`) — on top of differential equation models, grey
prediction, Markov chain models, and a rich suite of evaluation methods
(AHP, Entropy, CRITIC, PCA, TOPSIS, Fuzzy, RSR, DEA).

## Key Features

- **Evaluation Models** — AHP, Entropy weighting, CRITIC, PCA weighting,
  TOPSIS, Grey Relational Analysis (GRA), Rank Sum Ratio (RSR), Fuzzy
  Comprehensive Evaluation (FCE), Data Envelopment Analysis
  (CCR/BCC/SBM, Malmquist), plus inequality measures (Gini, Theil
  Index), coupling coordination degree, and obstacle degree models.
- **Prediction Models** — Grey prediction (GM(1,1), GM(1,N), Verhulst),
  Markov chain prediction
  ([`markov_chain()`](https://zhjx19.github.io/mathmodels/reference/markov.md),
  [`GM11_markov()`](https://zhjx19.github.io/mathmodels/reference/markov.md)),
  and a full **time series toolkit**:
  [`ts_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_transform.md)
  /
  [`ts_back_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_back_transform.md)
  for Box-Cox–model–forecast–back-transform workflows;
  [`ts_ets()`](https://zhjx19.github.io/mathmodels/reference/ts_ets.md),
  [`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md),
  [`ts_garch()`](https://zhjx19.github.io/mathmodels/reference/ts_garch.md),
  [`ts_sarima_garch()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima_garch.md)
  for modeling;
  [`ts_stl()`](https://zhjx19.github.io/mathmodels/reference/ts_stl.md)
  for decomposition;
  [`ts_test()`](https://zhjx19.github.io/mathmodels/reference/ts_test.md)
  for stationarity/ARCH tests;
  [`ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/ts_forecast.md)
  for unified forecasting; plus dedicated `plot_ts_*()` visualization
  functions.
- **Differential Equation Models** — String-formula
  [`ode_solver()`](https://zhjx19.github.io/mathmodels/reference/ode_solver.md)
  for arbitrary ODE systems; ready-to-use population models (Malthus,
  Logistic), epidemic compartment models (SI, SIS, SIR, SEIR), and
  Lotka–Volterra predator–prey model, all with a unified `init` +
  `params` interface. Epidemic visualization
  ([`plot_compartments()`](https://zhjx19.github.io/mathmodels/reference/plot_compartments.md),
  [`plot_incidence()`](https://zhjx19.github.io/mathmodels/reference/plot_incidence.md),
  [`plot_phase_si()`](https://zhjx19.github.io/mathmodels/reference/plot_phase_si.md),
  [`plot_Rt_estimate()`](https://zhjx19.github.io/mathmodels/reference/plot_Rt_estimate.md))
  and metrics
  ([`epi_metrics()`](https://zhjx19.github.io/mathmodels/reference/epi_metrics.md)
  for R0, peak, attack rate).
- **Tidyverse Integration** — Seamlessly works with `|>` and
  `dplyr`/`ggplot2` for smooth data manipulation and batch processing.

## Installation

You can install the latest development version of `mathmodels` directly
from [GitHub](https://github.com/zhjx19/mathmodels) use:

``` r

remotes::install_github("zhjx19/mathmodels")
```

Or download to current path, unzip and use:

``` r

# install.packages("lpSolveAPI")   # If necessary, install dependency packages first
install.packages("mathmodels-main", repos=NULL, type="source")
```

## Getting Started

``` r

library(mathmodels)

# --- AHP (Analytic Hierarchy Process) ---
A = matrix(c(1,   1/2, 4, 3,   3,
             2,   1,   7, 5,   5,
             1/4, 1/7, 1, 1/2, 1/3,
             1/3, 1/5, 2, 1,   1,
             1/3, 1/5, 3, 1,   1), byrow = TRUE, nrow = 5)
AHP(A)

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
metrics = epi_metrics(result, beta = 0.3, gamma = 0.1, N = 1000)
metrics$R0         # basic reproduction number
metrics$peak_time  # time of peak infection
```

## Learning More

For detailed documentation, tutorials, and in-depth examples on using
the `mathmodels` package, please refer to the comprehensive online
manual:

📘 **[mathmodels Package Manual - Simplifying Mathematical Modeling
(Online Book, In Chinese)](https://zhjx19.github.io/mathmodels-book/)**

This online book is the definitive guide to the package’s
functionalities. Currently implemented modules include:

- **Differential equation models**: Malthus, Logistic, SI, SIS, SIR,
  SEIR, Lotka–Volterra with
  [`ode_solver()`](https://zhjx19.github.io/mathmodels/reference/ode_solver.md)
  and `model_*()` functions; epidemic visualization
  ([`plot_compartments()`](https://zhjx19.github.io/mathmodels/reference/plot_compartments.md),
  [`plot_incidence()`](https://zhjx19.github.io/mathmodels/reference/plot_incidence.md),
  [`plot_phase_si()`](https://zhjx19.github.io/mathmodels/reference/plot_phase_si.md),
  [`plot_Rt_estimate()`](https://zhjx19.github.io/mathmodels/reference/plot_Rt_estimate.md))
  and metrics
  ([`epi_metrics()`](https://zhjx19.github.io/mathmodels/reference/epi_metrics.md))
- **Time series**:
  [`ts_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_transform.md),
  [`ts_back_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_back_transform.md),
  [`ts_ets()`](https://zhjx19.github.io/mathmodels/reference/ts_ets.md),
  [`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md),
  [`ts_garch()`](https://zhjx19.github.io/mathmodels/reference/ts_garch.md),
  [`ts_sarima_garch()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima_garch.md),
  [`ts_stl()`](https://zhjx19.github.io/mathmodels/reference/ts_stl.md),
  [`ts_test()`](https://zhjx19.github.io/mathmodels/reference/ts_test.md),
  [`ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/ts_forecast.md),
  `plot_ts_*()`
- **Markov chain prediction**:
  [`markov_chain()`](https://zhjx19.github.io/mathmodels/reference/markov.md)
  and
  [`GM11_markov()`](https://zhjx19.github.io/mathmodels/reference/markov.md)
- Indicator data preprocessing
- AHP, Entropy weighting, CRITIC, PCA weighting
- Weight combination techniques
- TOPSIS, Grey Relational Analysis (GRA)
- Rank Sum Ratio (RSR), Fuzzy Comprehensive Evaluation (FCE)
- Data Envelopment Analysis (CCR/BCC/SBM, Malmquist)
- Inequality Measures (Gini, Theil Index)
- Regional Economics (LQ/HHI/EG Index)
- Coupling coordination degree and obstacle degree
