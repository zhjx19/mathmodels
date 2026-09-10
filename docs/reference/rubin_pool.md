# Pool Estimates via Rubin's Rules

Pools point estimates and standard errors from multiple imputations
according to Rubin's rules. Degrees of freedom follow the Barnard-Rubin
(1999) small-sample adjustment when complete-data degrees of freedom are
supplied via `edf`.

## Usage

``` r
rubin_pool(est, err, alpha = 0.05, edf = Inf)
```

## Arguments

- est:

  Numeric vector of point estimates, one per imputation.

- err:

  Numeric vector of corresponding standard errors.

- alpha:

  Significance level for the confidence interval. Default 0.05.

- edf:

  Complete-data degrees of freedom (positive scalar or `Inf`). When
  finite, the Barnard-Rubin adjusted df are used; default `Inf` applies
  the Rubin large-sample df.

## Value

A list of class `rubin_pool`: `estimate`, `std_error`, `m`,
`var_between`, `var_within`, `var_total`, `df`, `r`, `fmi`, `rel_eff`,
`lcl`, `ucl`, `t`, `p_value` (two-sided test against 0), `alpha`, `edf`.

## Details

With Q_i and U_i the estimate and variance from imputation i, the pooled
estimate is the mean of Q_i; total variance T = Ubar + (1 + 1/m) B where
Ubar is the mean variance and B the between-imputation variance. The
fraction of missing information uses the standard Rubin definition (r +
2/(df+3))/(r+1); the relative efficiency is (1 + lambda/m)^-1. The
fraction of missing information and relative efficiency use the final
(Barnard-Rubin adjusted) degrees of freedom when `edf` is finite.

## Examples

``` r
est = c(0.345, 0.303, 0.298, 0.378)
err = c(0.056, 0.101, 0.099, 0.034)
rubin_pool(est, err)
#> Rubin pooling (m = 4 )
#>   estimate = 0.331000, std_error = 0.088634, df = 58.3
#>   95.0% CI: [0.153597, 0.508403], t = 3.7345, p = 0.0004297, FMI = 0.2521
rubin_pool(est, err, edf = 120)
#> Rubin pooling (m = 4 )
#>   estimate = 0.331000, std_error = 0.088634, df = 35.6
#>   95.0% CI: [0.151165, 0.510835], t = 3.7345, p = 0.0006572, FMI = 0.2670
```
