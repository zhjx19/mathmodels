# mathmodels - A Powerful Toolkit for Mathematical Modeling

### Introduction

Welcome to the `mathmodels` package! This is a comprehensive and
user-friendly R toolkit designed for practitioners and researchers in
the field of mathematical modeling. At its core, `mathmodels` embodies
the philosophy of “simplifying the complex,” encapsulating sophisticated
mathematical algorithms into intuitive functions. This allows you to
focus on model construction and interpretation, rather than getting
bogged down in implementation details.

The `mathmodels` package is developed as a companion to the book
“[Mathematical Modeling: Algorithms and Programming
Implementation](https://item.jd.com/13272295.html)” (China Machine
Press). The current version (0.0.14) spans the full modeling workflow:
**data cleaning and multiple imputation**, **statistical inference**,
**multivariate statistics**, a rich suite of **evaluation algorithms**
(weighting, TOPSIS, fuzzy evaluation, DEA, inequality measures, …),
**prediction models** (grey, Markov, regression, time series),
**interpolation and curve fitting**, and **differential equation /
epidemic models**.

### Key Features

1.  **Data Cleaning** — missing-value diagnostics
    ([`na_summary()`](https://zhjx19.github.io/mathmodels/reference/na_summary.md),
    [`na_mcar_test()`](https://zhjx19.github.io/mathmodels/reference/na_mcar_test.md)
    for Little’s MCAR test) with heatmap / bar visuals; imputation via
    simple fills
    ([`impute()`](https://zhjx19.github.io/mathmodels/reference/impute.md):
    group mean / median / mode, linear & spline interpolation,
    constants), a unified model-based framework
    ([`impute_model()`](https://zhjx19.github.io/mathmodels/reference/impute_model.md):
    kNN, linear regression, decision trees), and **multiple imputation**
    ([`impute_multiple()`](https://zhjx19.github.io/mathmodels/reference/impute_multiple.md):
    FCS chained equations with pmm/norm for numeric columns and
    logreg/polyreg for categorical ones) with **Rubin’s-rules pooling**
    ([`rubin_pool()`](https://zhjx19.github.io/mathmodels/reference/rubin_pool.md)
    for scalar estimates,
    [`pool_fit()`](https://zhjx19.github.io/mathmodels/reference/pool_fit.md)
    for lm/glm fits, including the Barnard-Rubin small-sample df
    adjustment); outlier detection and handling
    ([`outlier()`](https://zhjx19.github.io/mathmodels/reference/outlier.md));
    Excel-style
    [`pivot_table()`](https://zhjx19.github.io/mathmodels/reference/pivot_table.md).
2.  **Statistical Inference** — tidy hypothesis tests: t-test /
    Wilcoxon, ANOVA / Kruskal-Wallis, ANCOVA, correlation, chi-square,
    normality and variance-homogeneity tests (`stat_*()`), each
    returning one tidy row per test with effect sizes and `.by` slice
    grouping; generic **bootstrap** and **permutation** frameworks;
    [`export_table()`](https://zhjx19.github.io/mathmodels/reference/export_table.md)
    renders results as Word/LaTeX three-line tables.
3.  **Multivariate Statistics** — PCA
    ([`mv_pca()`](https://zhjx19.github.io/mathmodels/reference/mv_pca.md)),
    EFA
    ([`mv_efa()`](https://zhjx19.github.io/mathmodels/reference/mv_efa.md)),
    hierarchical and k-means clustering with validity metrics
    ([`mv_hclust()`](https://zhjx19.github.io/mathmodels/reference/mv_hclust.md),
    [`mv_kmeans()`](https://zhjx19.github.io/mathmodels/reference/mv_kmeans.md),
    [`mv_cluster_metrics()`](https://zhjx19.github.io/mathmodels/reference/mv_cluster_metrics.md)),
    discriminant analysis
    ([`mv_discrim()`](https://zhjx19.github.io/mathmodels/reference/mv_discrim.md)),
    correspondence analysis
    ([`mv_corresp()`](https://zhjx19.github.io/mathmodels/reference/mv_corresp.md)),
    MDS
    ([`mv_mds()`](https://zhjx19.github.io/mathmodels/reference/mv_mds.md)),
    canonical and partial correlation
    ([`mv_cancor()`](https://zhjx19.github.io/mathmodels/reference/mv_cancor.md),
    [`mv_pcor()`](https://zhjx19.github.io/mathmodels/reference/mv_pcor.md)),
    all with tidy components and `plot_mv_*()` visuals.
4.  **Evaluation Models** — AHP, Entropy weighting, CRITIC, CV, PCA and
    factor-analysis weighting (`weight_*()` family), weight combination
    ([`combine_weights()`](https://zhjx19.github.io/mathmodels/reference/combine_weights.md)),
    TOPSIS, Grey Relational Analysis (GRA), Rank Sum Ratio (RSR), Fuzzy
    Comprehensive Evaluation (FCE) with the `mf_*` membership-function
    family, Data Envelopment Analysis (CCR/BCC/SBM, Malmquist),
    inequality measures (Gini, Theil index), regional economics
    (LQ/HHI/EG), and coupling coordination / obstacle degree models.
5.  **Prediction Models** — grey prediction (GM(1,1), GM(1,N), DGM(2,1),
    Verhulst), Markov chain prediction
    ([`markov_chain()`](https://zhjx19.github.io/mathmodels/reference/markov_chain.md),
    [`GM11_markov()`](https://zhjx19.github.io/mathmodels/reference/GM11_markov.md)),
    a regression prediction toolkit
    ([`reg_lm()`](https://zhjx19.github.io/mathmodels/reference/reg_lm.md),
    [`reg_logistic()`](https://zhjx19.github.io/mathmodels/reference/reg_logistic.md),
    [`reg_poisson()`](https://zhjx19.github.io/mathmodels/reference/reg_poisson.md),
    [`reg_negbin()`](https://zhjx19.github.io/mathmodels/reference/reg_negbin.md)
    with stepwise selection;
    [`reg_predict()`](https://zhjx19.github.io/mathmodels/reference/reg_predict.md),
    [`reg_diagnostics()`](https://zhjx19.github.io/mathmodels/reference/reg_diagnostics.md),
    `plot_reg_*()`), and a lightweight time-series workflow built around
    `ts_df`
    ([`as_ts_df()`](https://zhjx19.github.io/mathmodels/reference/as_ts_df.md),
    [`complete_ts_df()`](https://zhjx19.github.io/mathmodels/reference/complete_ts_df.md)
    /
    [`impute_ts_df()`](https://zhjx19.github.io/mathmodels/reference/impute_ts_df.md),
    [`ts_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_transform.md)
    /
    [`ts_back_transform()`](https://zhjx19.github.io/mathmodels/reference/ts_back_transform.md),
    [`ts_test()`](https://zhjx19.github.io/mathmodels/reference/ts_test.md),
    [`ts_stl()`](https://zhjx19.github.io/mathmodels/reference/ts_stl.md),
    [`ts_ets()`](https://zhjx19.github.io/mathmodels/reference/ts_ets.md),
    [`ts_sarima()`](https://zhjx19.github.io/mathmodels/reference/ts_sarima.md),
    [`ts_arimax()`](https://zhjx19.github.io/mathmodels/reference/ts_arimax.md),
    [`ts_forecast()`](https://zhjx19.github.io/mathmodels/reference/ts_forecast.md),
    plus `plot_ts*()`).
6.  **Interpolation and Curve Fitting** — piecewise linear, cubic
    spline, polynomial and Hermite interpolation
    ([`interp_linear()`](https://zhjx19.github.io/mathmodels/reference/interp_linear.md),
    [`interp_spline()`](https://zhjx19.github.io/mathmodels/reference/interp_spline.md),
    [`interp_poly()`](https://zhjx19.github.io/mathmodels/reference/interp_poly.md),
    [`interp_hermite()`](https://zhjx19.github.io/mathmodels/reference/interp_hermite.md)),
    plus curve fitting
    ([`poly_fit()`](https://zhjx19.github.io/mathmodels/reference/poly_fit.md),
    [`curve_fit()`](https://zhjx19.github.io/mathmodels/reference/curve_fit.md),
    [`growth_fit()`](https://zhjx19.github.io/mathmodels/reference/growth_fit.md))
    with automatic starting values for logistic, Gompertz and
    Michaelis-Menten nonlinear models.
7.  **Differential Equation Models** — string-formula
    [`ode_solver()`](https://zhjx19.github.io/mathmodels/reference/ode_solver.md)
    for arbitrary ODE systems; ready-to-use population models (Malthus,
    Logistic), epidemic compartment models (SI, SIS, SIR, SEIR) and the
    Lotka-Volterra predator-prey model, all with a unified `init` +
    `params` interface; epidemic visualization
    ([`plot_compartments()`](https://zhjx19.github.io/mathmodels/reference/plot_compartments.md),
    [`plot_incidence()`](https://zhjx19.github.io/mathmodels/reference/plot_incidence.md),
    [`plot_phase_si()`](https://zhjx19.github.io/mathmodels/reference/plot_phase_si.md),
    [`plot_Rt_estimate()`](https://zhjx19.github.io/mathmodels/reference/plot_Rt_estimate.md))
    and metrics
    ([`epi_metrics()`](https://zhjx19.github.io/mathmodels/reference/epi_metrics.md)
    for R0, peak, attack rate).
8.  **Elegant, Tidyverse-friendly Interface** — standardized argument
    names, `|>` pipes and `dplyr`/`tidyr`/`ggplot2` integration, with
    data preprocessing (normalization, scaling, indicator direction)
    handled internally.

### Quick Start

To start using `mathmodels`, first ensure it’s installed and loaded:

``` r

# Install from GitHub (recommended)
remotes::install_github("zhjx19/mathmodels")

# Or install from a local source package
install.packages("mathmodels.tar.gz", repos = NULL, type = "source")
```

``` r

library(dplyr)
library(mathmodels)
```

Once loaded, you can immediately use its functions. For example, perform
a quick entropy weighting TOPSIS analysis using the built-in
`water_quality` dataset:

``` r

# View example data
water_quality
```

    ## # A tibble: 20 × 5
    ##       ID    O2    PH  germ nutrient
    ##    <dbl> <dbl> <dbl> <dbl>    <dbl>
    ##  1     1  4.69  6.59    51    11.9 
    ##  2     2  2.03  7.86    19     6.46
    ##  3     3  9.11  6.31    46     8.91
    ##  4     4  8.61  7.05    46    26.4 
    ##  5     5  7.13  6.5     50    23.6 
    ##  6     6  2.39  6.77    38    24.6 
    ##  7     7  7.69  6.79    38     6.01
    ##  8     8  9.3   6.81    27    31.6 
    ##  9     9  5.45  7.62     5    18.5 
    ## 10    10  6.19  7.27    17     7.51
    ## 11    11  7.93  7.53     9     6.52
    ## 12    12  4.4   7.28    17    25.3 
    ## 13    13  7.46  8.24    23    14.4 
    ## 14    14  2.01  5.55    47    26.3 
    ## 15    15  2.04  6.4     23    17.9 
    ## 16    16  7.73  6.14    52    15.7 
    ## 17    17  6.35  7.58    25    29.5 
    ## 18    18  8.29  8.41    39    12.0 
    ## 19    19  3.54  7.27    54     3.16
    ## 20    20  7.44  6.26     8    28.4

``` r

# Preprocess data: handle centered and interval-type indicators
df = water_quality |>
  mutate(
    # Optimal PH is 7 -> Centered-type positive
    PH = rescale_middle(PH, 7),
    # Optimal nutrient is 10-20 -> Interval-type positive
    nutrient = rescale_interval(nutrient, 10, 20)
  )
```

Calculate weights using the Entropy Method (normalization is handled
internally). First, specify the indicator directions: since the
indicator columns PH and nutrient have been pre-normalized, their
directions are set to `NA` to avoid repeating the normalization process.

``` r

idx = c("+", NA, "-", NA)
res = weight_entropy(df[-1], idx)
res$w                         # View calculated weights
```

    ##        O2        PH      germ  nutrient 
    ## 0.3153202 0.1813149 0.3506929 0.1526721

When performing TOPSIS (which has built-in normalization), you need to
specify the indicator directions in order to correctly determine the
positive and negative ideal solutions.

``` r

idx = c("+","+","-","+")
topsis(df[2:5], res$w, idx)
```

    ##         1         2         3         4         5         6         7         8 
    ## 0.3841556 0.4852199 0.5021310 0.5084128 0.4226200 0.3945085 0.5502261 0.6258551 
    ##         9        10        11        12        13        14        15        16 
    ## 0.7260529 0.7138447 0.8091923 0.6159619 0.6007565 0.1593096 0.4974685 0.4250799 
    ##        17        18        19        20 
    ## 0.5565555 0.4815844 0.3070290 0.7098355

The obtained result is the relative closeness score.

The same tidy workflow extends to the other modules. For instance,
**multiple imputation** of missing values followed by **Rubin’s-rules
pooling** of a regression model:

``` r

set.seed(1)
d = data.frame(x = rnorm(40), z = rnorm(40))
d$y = 0.5 * d$x + rnorm(40)
d$x[c(3, 9)] = NA                      # introduce missingness

imps = impute_multiple(d, .cols = x, m = 5, maxit = 3, seed = 2)  # 5 imputations
fits = lapply(imps, \(dat) lm(y ~ x, dat))                        # fit each
pool_fit(fits)                         # pooled coefficients (Rubin's rules)
```

    ## Pooled estimates (Rubin's rules)
    ## # A tibble: 2 × 7
    ##   term        estimate std_error     lcl   ucl p_value    fmi
    ##   <chr>          <dbl>     <dbl>   <dbl> <dbl>   <dbl>  <dbl>
    ## 1 (Intercept)    0.134     0.133 -0.135  0.403  0.319  0.0520
    ## 2 x              0.277     0.153 -0.0325 0.587  0.0778 0.0533

### Getting Help and In-Depth Learning

This Vignette offers just a glimpse into the capabilities of
`mathmodels`. For a thorough understanding of all algorithms, their
detailed principles, function usage, and practical case studies, we
highly recommend consulting our comprehensive **online companion
manual**:

📘 **[mathmodels Package Manual - Simplifying Mathematical Modeling
(Online Book)](https://zhjx19.github.io/mathmodels-book/)**

This online manual serves as the definitive documentation for the
`mathmodels` package. It consolidates content from the package vignettes
and covers the evaluation algorithms, statistical inference,
multivariate statistics, data cleaning, prediction models and
differential equation models, with R implementations and extensive
examples, making it your best resource for learning and using this
package.

------------------------------------------------------------------------

We hope `mathmodels` proves to be a valuable and reliable companion on
your journey through mathematical modeling and scientific research.
