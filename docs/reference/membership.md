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
tri_mf(x, params)

trap_mf(x, params)

gauss_mf(x, params)

gbell_mf(x, params)

gauss2mf(x, params)

sigmoid_mf(x, params)

dsigmoid_mf(x, params)

psigmoid_mf(x, params)

z_mf(x, params)

pi_mf(x, params)

s_mf(x, params)

plot_mf(mf, xlim = c(0, 10), main = NULL)
```

## Arguments

- x:

  Numeric vector, input values for which to compute membership.

- params:

  Numeric vector, parameters defining the membership function:

  - For `tri_mf`: `c(a, b, c)`, where `a <= b <= c` (left base, peak,
    right base).

  - For `trap_mf`: `c(a, b, c, d)`, where `a <= b <= c <= d` (left base,
    left top, right top, right base).

  - For `gauss_mf`: `c(sigma, c)`, where `sigma > 0` (spread, center).

  - For `gbell_mf`: `c(a, b, c)`, where `a > 0`, `b > 0` (width, shape,
    center).

  - For `gauss2mf`: `c(s1, c1, s2, c2)`, where `s1 > 0`, `s2 > 0` (left
    spread, left center, right spread, right center).

  - For `sigmoid_mf`: `c(a, b)`, where `a > 0` (slope, inflection
    point).

  - For `dsigmoid_mf`: `c(a1, c1, a2, c2)`, where `a1 > 0`, `a2 > 0`
    (slopes and inflection points for two sigmoids).

  - For `psigmoid_mf`: `c(a1, c1, a2, c2)`, where `a1 > 0`, `a2 > 0`
    (slopes and inflection points for two sigmoids).

  - For `z_mf`: `c(a, b)`, where `a < b` (left base, right base).

  - For `pi_mf`: `c(a, b, c, d)`, where `a < b < c < d` (left base, left
    shoulder, right shoulder, right base).

  - For `s_mf`: `c(a, b)`, where `a < b` (left base, right base).

- mf:

  Function, a membership function with fixed parameters (e.g.,
  `function(x) tri_mf(x, c(2, 5, 8))`).

- xlim:

  Numeric vector of length 2, x-axis limits for plotting (default
  `c(0, 10)`).

- main:

  Character, plot title (default `NULL`, no title).

## Value

- For membership functions (`tri_mf`, `trap_mf`, `gauss_mf`, `gbell_mf`,
  `gauss2mf`, `sigmoid_mf`, `dsigmoid_mf`, `psigmoid_mf`, `z_mf`,
  `pi_mf`, `s_mf`): A numeric vector of membership values in `[0, 1]`,
  same length as `x`.

- For `plot_mf`: A ggplot2 object, plotting the membership function.

## Details

These functions support evaluation models in mathematical modeling:

- `tri_mf`: Triangular membership, linear rise from `a` to `b` (peak)
  and fall to `c`.

- `trap_mf`: Trapezoidal membership, linear rise from `a` to `b`,
  plateau from `b` to `c`, fall to `d`.

- `gauss_mf`: Gaussian membership, bell-shaped curve centered at `c`
  with spread `sigma`.

- `gbell_mf`: Generalized bell membership, bell-shaped curve with width
  `a`, shape `b`, and center `c`.

- `gauss2mf`: Two-parameter Gaussian membership, combining two Gaussians
  with spreads `s1`, `s2` and centers `c1`, `c2`.

- `sigmoid_mf`: Sigmoid membership, S-shaped curve with slope `a` and
  inflection point `b`.

- `dsigmoid_mf`: Difference of two sigmoids, combining slopes `a1`, `a2`
  and inflection points `c1`, `c2`.

- `psigmoid_mf`: Product of two sigmoids, combining slopes `a1`, `a2`
  and inflection points `c1`, `c2`.

- `z_mf`: Z-shaped membership, decreasing from 1 at `a` to 0 at `b`.

- `pi_mf`: PI-shaped membership, rising from `a` to `b`, plateau from
  `b` to `c`, falling to `d`.

- `s_mf`: S-shaped membership, increasing from 0 at `a` to 1 at `b`.

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
tri_mf(x, params = c(3, 6, 8))
#>  [1] 0.0000000 0.0000000 0.0000000 0.0000000 0.3333333 0.6666667 1.0000000
#>  [8] 0.5000000 0.0000000 0.0000000 0.0000000

# Trapezoidal membership
trap_mf(x, params = c(1, 5, 7, 8))
#>  [1] 0.00 0.00 0.25 0.50 0.75 1.00 1.00 1.00 0.00 0.00 0.00

# Gaussian membership
gauss_mf(x, params = c(2, 5))
#>  [1] 0.04393693 0.13533528 0.32465247 0.60653066 0.88249690 1.00000000
#>  [7] 0.88249690 0.60653066 0.32465247 0.13533528 0.04393693

# Generalized bell membership
gbell_mf(x, params = c(2, 4, 6))
#>  [1] 0.0001523926 0.0006549308 0.0038910506 0.0375531759 0.5000000000
#>  [6] 0.9961089494 1.0000000000 0.9961089494 0.5000000000 0.0375531759
#> [11] 0.0038910506

# Two-parameter Gaussian membership
gauss2mf(x, params = c(1, 3, 3, 4))
#>  [1] 0.0111090 0.1353353 0.6065307 1.0000000 1.0000000 0.9459595 0.8007374
#>  [8] 0.6065307 0.4111123 0.2493522 0.1353353

# Sigmoid membership
sigmoid_mf(x, params = c(2, 4))
#>  [1] 0.0003353501 0.0024726232 0.0179862100 0.1192029220 0.5000000000
#>  [6] 0.8807970780 0.9820137900 0.9975273768 0.9996646499 0.9999546021
#> [11] 0.9999938558

# Difference of sigmoids membership
dsigmoid_y = dsigmoid_mf(x, params = c(5, 2, 5, 7))

# Product of sigmoids membership
psigmoid_mf(x, params = c(2, 3, -5, 8))
#>  [1] 2.472623e-03 1.798621e-02 1.192029e-01 5.000000e-01 8.807971e-01
#>  [6] 9.820135e-01 9.974821e-01 9.929740e-01 4.999773e-01 6.692810e-03
#> [11] 4.539783e-05

# Z-shaped membership
z_mf(x, params = c(3, 7))
#>  [1] 1.000 1.000 1.000 1.000 0.875 0.500 0.125 0.000 0.000 0.000 0.000

# PI-shaped membership
pi_mf(x, params = c(1, 4, 5, 10))
#>  [1] 0.0000000 0.0000000 0.2222222 0.7777778 1.0000000 1.0000000 0.9200000
#>  [8] 0.6800000 0.3200000 0.0800000 0.0000000

# S-shaped membership
s_mf(x, params = c(1, 8))
#>  [1] 0.00000000 0.00000000 0.04081633 0.16326531 0.36734694 0.63265306
#>  [7] 0.83673469 0.95918367 1.00000000 1.00000000 1.00000000

if (FALSE) { # \dontrun{
# Visualize membership functions
plot_mf(\(x) tri_mf(x, c(3, 6, 8)), main = "Triangular MF")
plot_mf(\(x) trap_mf(x, c(1, 5, 7, 8)), main = "Trapezoidal MF")
plot_mf(\(x) gauss_mf(x, c(2, 5)), main = "Gaussian MF")
plot_mf(\(x) gbell_mf(x, c(2, 4, 6)), main = "Generalized Bell MF")
plot_mf(\(x) gauss2mf(x, c(1, 3, 3, 4)), main = "Two-Parameter Gaussian MF")
plot_mf(\(x) sigmoid_mf(x, c(2, 4)), main = "Sigmoid MF")
plot_mf(\(x) dsigmoid_mf(x, c(5, 2, 5, 7)), main = "Difference of Sigmoids MF")
plot_mf(\(x) psigmoid_mf(x, c(2, 3, -5, 8)), main = "Product of Sigmoids MF")
plot_mf(\(x) z_mf(x, c(3, 7)), main = "Z-Shaped MF")
plot_mf(\(x) pi_mf(x, c(1, 4, 5, 10)), main = "PI-Shaped MF")
plot_mf(\(x) s_mf(x, c(1, 8)), main = "S-Shaped MF")
} # }
```
