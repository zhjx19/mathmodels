# Membership Functions for Fuzzy Logic

A collection of functions to compute membership values for various fuzzy
sets, including triangular, trapezoidal, Gaussian, generalized bell,
two-parameter Gaussian, sigmoid, difference of sigmoids, product of
sigmoids, Z-shaped, PI-shaped, and S-shaped membership functions.
Includes a function to visualize membership functions using ggplot2.
These are designed for evaluation models in mathematical modeling,
compatible with `fuzzy_eval` in the `mathmodels` package.

## Usage

``` r
mf_tri(x, params)

mf_trap(x, params)

mf_gauss(x, params)

mf_gbell(x, params)

mf_gauss2(x, params)

mf_sigmoid(x, params)

mf_dsigmoid(x, params)

mf_psigmoid(x, params)

mf_z(x, params)

mf_pi(x, params)

mf_s(x, params)

plot_mf(mf, xlim = c(0, 10), main = NULL)
```

## Arguments

- x:

  Numeric vector, input values for which to compute membership.

- params:

  Numeric vector, parameters defining the membership function:

  - For `mf_tri`: `c(a, b, c)`, where `a < b < c` (left base, peak,
    right base).

  - For `mf_trap`: `c(a, b, c, d)`, where `a < b <= c < d` (left base,
    left top, right top, right base); one-sided shoulders are supported
    with `a = b = -Inf` (no rise) or `c = d = Inf` (no fall).

  - For `mf_gauss`: `c(sigma, c)`, where `sigma > 0` (spread, center).

  - For `mf_gbell`: `c(a, b, c)`, where `a > 0`, `b > 0` (width, shape,
    center).

  - For `mf_gauss2`: `c(s1, c1, s2, c2)`, where `s1 > 0`, `s2 > 0` (left
    spread, left center, right spread, right center).

  - For `mf_sigmoid`: `c(a, b)`, where `a > 0` (slope, inflection
    point).

  - For `mf_dsigmoid`: `c(a1, c1, a2, c2)`, where `a1 > 0`, `a2 > 0`
    (slopes and inflection points for two sigmoids).

  - For `mf_psigmoid`: `c(a1, c1, a2, c2)`, where `a1 > 0`, `a2 > 0`
    (slopes and inflection points for two sigmoids).

  - For `mf_z`: `c(a, b)`, where `a < b` (left base, right base).

  - For `mf_pi`: `c(a, b, c, d)`, where `a < b <= c < d` (left base,
    left shoulder, right shoulder, right base).

  - For `mf_s`: `c(a, b)`, where `a < b` (left base, right base).

- mf:

  Function, a membership function with fixed parameters (e.g.,
  `function(x) mf_tri(x, c(2, 5, 8))`).

- xlim:

  Numeric vector of length 2, x-axis limits for plotting (default
  `c(0, 10)`).

- main:

  Character, plot title (default `NULL`, no title).

## Value

- For membership functions (`mf_tri`, `mf_trap`, `mf_gauss`, `mf_gbell`,
  `mf_gauss2`, `mf_sigmoid`, `mf_psigmoid`, `mf_z`, `mf_pi`, `mf_s`): A
  numeric vector of membership values in `[0, 1]`, same length as `x`.

- For `mf_dsigmoid`: A numeric vector holding the difference of two
  sigmoids. Values lie in \\\[-1, 1\]\\ and are **negative when
  `c1 > c2`**, so unlike the other membership functions it is *not*
  guaranteed to be in `[0, 1]`. Clip the result to `[0, 1]`, or use
  `mf_psigmoid()` instead, before passing it to
  [`fuzzy_eval()`](https://zhjx19.github.io/mathmodels/reference/fuzzy_eval.md).

- For `plot_mf`: A ggplot2 object, plotting the membership function.

## Details

These functions support evaluation models in mathematical modeling:

- `mf_tri`: Triangular membership, linear rise from `a` to `b` (peak)
  and fall to `c`.

- `mf_trap`: Trapezoidal membership, linear rise from `a` to `b`,
  plateau from `b` to `c`, fall to `d`.

- `mf_gauss`: Gaussian membership, bell-shaped curve centered at `c`
  with spread `sigma`.

- `mf_gbell`: Generalized bell membership, bell-shaped curve with width
  `a`, shape `b`, and center `c`.

- `mf_gauss2`: Two-parameter Gaussian membership, combining two
  Gaussians with spreads `s1`, `s2` and centers `c1`, `c2`.

- `mf_sigmoid`: Sigmoid membership, S-shaped curve with slope `a` and
  inflection point `b`.

- `mf_dsigmoid`: Difference of two sigmoids, combining slopes `a1`, `a2`
  and inflection points `c1`, `c2`.

- `mf_psigmoid`: Product of two sigmoids, combining slopes `a1`, `a2`
  and inflection points `c1`, `c2`.

- `mf_z`: Z-shaped membership, decreasing from 1 at `a` to 0 at `b`.

- `mf_pi`: PI-shaped membership, rising from `a` to `b`, plateau from
  `b` to `c`, falling to `d`.

- `mf_s`: S-shaped membership, increasing from 0 at `a` to 1 at `b`.

- `plot_mf`: Plots a membership function over `xlim` using ggplot2,
  suitable for tidyverse workflows.

Membership values can be used to construct fuzzy evaluation matrices for
`fuzzy_eval`. Implemented in base R, except `plot_mf`, which requires
ggplot2.

## Examples

``` r
# Define input values
x = 0:10

# Triangular membership
mf_tri(x, params = c(3, 6, 8))
#>  [1] 0.0000000 0.0000000 0.0000000 0.0000000 0.3333333 0.6666667 1.0000000
#>  [8] 0.5000000 0.0000000 0.0000000 0.0000000

# Trapezoidal membership
mf_trap(x, params = c(1, 5, 7, 8))
#>  [1] 0.00 0.00 0.25 0.50 0.75 1.00 1.00 1.00 0.00 0.00 0.00

# Gaussian membership
mf_gauss(x, params = c(2, 5))
#>  [1] 0.04393693 0.13533528 0.32465247 0.60653066 0.88249690 1.00000000
#>  [7] 0.88249690 0.60653066 0.32465247 0.13533528 0.04393693

# Generalized bell membership
mf_gbell(x, params = c(2, 4, 6))
#>  [1] 0.0001523926 0.0006549308 0.0038910506 0.0375531759 0.5000000000
#>  [6] 0.9961089494 1.0000000000 0.9961089494 0.5000000000 0.0375531759
#> [11] 0.0038910506

# Two-parameter Gaussian membership
mf_gauss2(x, params = c(1, 3, 3, 4))
#>  [1] 0.0111090 0.1353353 0.6065307 1.0000000 1.0000000 0.9459595 0.8007374
#>  [8] 0.6065307 0.4111123 0.2493522 0.1353353

# Sigmoid membership
mf_sigmoid(x, params = c(2, 4))
#>  [1] 0.0003353501 0.0024726232 0.0179862100 0.1192029220 0.5000000000
#>  [6] 0.8807970780 0.9820137900 0.9975273768 0.9996646499 0.9999546021
#> [11] 0.9999938558

# Difference of sigmoids membership
dsigmoid_y = mf_dsigmoid(x, params = c(5, 2, 5, 7))

# Product of sigmoids membership
mf_psigmoid(x, params = c(2, 3, -5, 8))
#>  [1] 2.472623e-03 1.798621e-02 1.192029e-01 5.000000e-01 8.807971e-01
#>  [6] 9.820135e-01 9.974821e-01 9.929740e-01 4.999773e-01 6.692810e-03
#> [11] 4.539783e-05

# Z-shaped membership
mf_z(x, params = c(3, 7))
#>  [1] 1.000 1.000 1.000 1.000 0.875 0.500 0.125 0.000 0.000 0.000 0.000

# PI-shaped membership
mf_pi(x, params = c(1, 4, 5, 10))
#>  [1] 0.0000000 0.0000000 0.2222222 0.7777778 1.0000000 1.0000000 0.9200000
#>  [8] 0.6800000 0.3200000 0.0800000 0.0000000

# S-shaped membership
mf_s(x, params = c(1, 8))
#>  [1] 0.00000000 0.00000000 0.04081633 0.16326531 0.36734694 0.63265306
#>  [7] 0.83673469 0.95918367 1.00000000 1.00000000 1.00000000

if (FALSE) { # \dontrun{
# Visualize membership functions
plot_mf(\(x) mf_tri(x, c(3, 6, 8)), main = "Triangular MF")
plot_mf(\(x) mf_trap(x, c(1, 5, 7, 8)), main = "Trapezoidal MF")
plot_mf(\(x) mf_gauss(x, c(2, 5)), main = "Gaussian MF")
plot_mf(\(x) mf_gbell(x, c(2, 4, 6)), main = "Generalized Bell MF")
plot_mf(\(x) mf_gauss2(x, c(1, 3, 3, 4)), main = "Two-Parameter Gaussian MF")
plot_mf(\(x) mf_sigmoid(x, c(2, 4)), main = "Sigmoid MF")
plot_mf(\(x) mf_dsigmoid(x, c(5, 2, 5, 7)), main = "Difference of Sigmoids MF")
plot_mf(\(x) mf_psigmoid(x, c(2, 3, -5, 8)), main = "Product of Sigmoids MF")
plot_mf(\(x) mf_z(x, c(3, 7)), main = "Z-Shaped MF")
plot_mf(\(x) mf_pi(x, c(1, 4, 5, 10)), main = "PI-Shaped MF")
plot_mf(\(x) mf_s(x, c(1, 8)), main = "S-Shaped MF")
} # }
```
