# Bootstrap Estimation and Confidence Intervals

Generic bootstrap framework built on
[`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html). The
`statistic` follows the boot convention `function(data, indices, ...)` —
use `indices` to index the original data inside the statistic (the
framework never copies the data in the resampling loop). Confidence
intervals via
[`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html): BCa
(with automatic jackknife), percentile, or basic.

## Usage

``` r
stat_bootstrap(
  data,
  statistic,
  ...,
  B = 2000,
  type = c("bca", "percentile", "basic"),
  conf_level = 0.95,
  seed = NULL,
  .by = NULL,
  parallel = "no",
  ncpus = 1
)
```

## Arguments

- data:

  A data frame (or any object accepted by `statistic`).

- statistic:

  A function `function(data, indices, ...)` returning a scalar or a
  named numeric vector.

- ...:

  Additional arguments passed to `statistic` (and to
  [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html) where
  unambiguous).

- B:

  Number of bootstrap resamples.

- type:

  `"bca"` (default), `"percentile"`, or `"basic"`.

- conf_level:

  Confidence level.

- seed:

  Optional integer seed for reproducibility.

- .by:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Optional slice columns; the bootstrap is run separately within each
  slice.

- parallel:

  [`boot::boot`](https://rdrr.io/pkg/boot/man/boot.html) parallel
  option: `"no"` (default), `"multicore"` (not on Windows), or `"snow"`.

- ncpus:

  Number of workers for parallel boot.

## Value

A tibble with class `stat_infer`: `.by` identifiers, `term`, `estimate`,
`se`, `bias`, `conf_low`, `conf_high`, `type`, `B`.

## Examples

``` r
set.seed(1)
stat_bootstrap(mtcars, statistic = \(d, i) mean(d$mpg[i]))
#> Bootstrap (bca CI, B = 2000) 
#> # A tibble: 1 × 8
#>   term  estimate    se     bias conf_low conf_high type      B
#> * <chr>    <dbl> <dbl>    <dbl>    <dbl>     <dbl> <chr> <dbl>
#> 1 value     20.1  1.07 -0.00861     18.2      22.3 bca    2000
stat_bootstrap(mtcars,
  statistic = \(d, i) c(mean = mean(d$mpg[i]), sd = stats::sd(d$mpg[i])))
#> Bootstrap (bca CI, B = 2000) 
#> # A tibble: 2 × 8
#>   term  estimate    se     bias conf_low conf_high type      B
#> * <chr>    <dbl> <dbl>    <dbl>    <dbl>     <dbl> <chr> <dbl>
#> 1 mean     20.1  1.05   0.00224    18.2      22.4  bca    2000
#> 2 sd        6.03 0.722 -0.124       4.75      7.50 bca    2000
```
