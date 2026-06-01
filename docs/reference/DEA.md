# Data Envelopment Analysis and Malmquist Productivity Index

Computes standard and super-efficiency Data Envelopment Analysis (DEA)
models (CCR, BCC, and slacks-based), including efficiency scores,
slacks, lambdas, targets, returns to scale, and references, with support
for undesirable outputs. Additionally, calculates the Malmquist
productivity index to measure productivity changes over time, including
efficiency change (EC), technical change (TC), and other decomposition
components.

## Usage

``` r
basic_DEA(
  data,
  inputs,
  outputs,
  ud_outputs = NULL,
  orientation = "io",
  rts = "vrs"
)

super_DEA(
  data,
  inputs,
  outputs,
  ud_outputs = NULL,
  orientation = "io",
  rts = "vrs"
)

basic_SBM(
  data,
  inputs,
  outputs,
  ud_outputs = NULL,
  orientation = "io",
  rts = "vrs"
)

super_SBM(data, inputs, outputs, orientation = "io", rts = "vrs")

malmquist(
  data,
  period,
  inputs,
  outputs,
  orientation = "oo",
  rts = "vrs",
  type1 = "glob",
  type2 = "rd"
)
```

## Arguments

- data:

  For DEA functions (`basic_DEA`, `super_DEA`, `basic_SBM`,
  `super_SBM`): A data frame where the first column contains DMU
  (Decision Making Unit) names/identifiers, and subsequent columns are
  input/output variables. For `malmquist`: A long-format data frame with
  a period column, DMU identifiers, and input/output variables.

- inputs:

  A numeric vector of column indices or a character vector of column
  names indicating input variables.

- outputs:

  A numeric vector of column indices or a character vector of column
  names indicating (desirable) output variables.

- ud_outputs:

  Optional. A numeric vector indicating the position of undesirable
  outputs within the `outputs` parameter. Defaults to `NULL`. Not
  applicable for `super_SBM` or `malmquist`.

- orientation:

  Character string. Model orientation: `"io"` (input-oriented, default
  for DEA functions) or `"oo"` (output-oriented, default for
  `malmquist`).

- rts:

  Character string. Returns to scale assumption: `"vrs"` (variable
  returns to scale, default) or `"crs"` (constant returns to scale).

- period:

  For `malmquist` only: A numeric or character index/name indicating the
  column in `data` containing time periods.

- type1:

  For `malmquist` only: Reference technology for Malmquist index:
  `"cont"` (contemporary), `"seq"` (sequential), or `"glob"` (global,
  default).

- type2:

  For `malmquist` only: Decomposition method for Malmquist index:
  `"fgnz"` (Färe et al., 1994), `"rd"` (Ray and Desli, 1997, default),
  `"gl"` (generalized), or `"bias"` (biased).

## Value

For DEA functions (`basic_DEA`, `super_DEA`, `basic_SBM`, `super_SBM`):
A list containing six elements:

- `efficiencies`: A named numeric vector of efficiency scores for each
  DMU, standardized to (0, 1\] for both input- and output-oriented
  models.

- `slacks`: A data frame or matrix containing slack values for inputs
  and outputs (including undesirable outputs, if specified).

- `lambdas`: A matrix or data frame of intensity variables
  (\\\lambda\\), representing the contribution of reference DMUs to the
  efficiency frontier (self-excluded in super-efficiency models).

- `targets`: A data frame or matrix of efficient target values for
  inputs and outputs, adjusted for undesirable outputs in DDF models.

- `returns`: A character vector indicating returns-to-scale status for
  each DMU: `"crs"` (constant), `"irs"` (increasing), or `"drs"`
  (decreasing).

- `references`: A matrix or data frame listing reference DMUs (peers)
  contributing to the efficiency frontier (\\\lambda \> 0\\).

For `malmquist`: A data frame containing:

- `Period`: Time period transitions (e.g., "t~t+1").

- `DMU`: Decision Making Unit identifiers.

- `mi`: Malmquist productivity index, measuring total productivity
  change.

- `ec`: Efficiency change (EC), computed directly or as `mi / tc` if not
  available from the model.

- `tc`: Technical change (TC), measuring frontier shift.

- `pech`: Pure efficiency change (if applicable, based on decomposition
  method).

- `sech`: Scale efficiency change (if applicable, based on decomposition
  method).

## Details

This package provides a unified interface for computing efficiency
scores and productivity changes using the deaR package. It includes five
functions: `basic_DEA`, `super_DEA`, `basic_SBM`, `super_SBM`, and
`malmquist`, each tailored to specific DEA or productivity analysis
models.

- **DEA Models**:

  - `basic_DEA`: Implements standard radial DEA models (CCR for CRS, BCC
    for VRS) as described by Charnes et al. (1978) and Banker et al.
    (1984), optimizing radial efficiency (input contraction or output
    expansion).

  - `super_DEA`: Computes super-efficiency radial DEA, excluding the
    evaluated DMU from the reference set to allow efficiency scores
    beyond 1 (output-oriented) or below 1 (input-oriented) for efficient
    DMUs (Andersen & Petersen, 1993).

  - `basic_SBM`: Implements standard Slacks-Based Measure (SBM) models
    (Tone, 2001), optimizing input and output slacks for a non-radial
    efficiency measure.

  - `super_SBM`: Combines SBM with super-efficiency properties,
    excluding the evaluated DMU from the reference set. Note:
    `super_SBM` does not support undesirable outputs.

- **Malmquist Productivity Index**:

  - `malmquist`: Calculates the Malmquist productivity index to measure
    productivity changes over time, decomposing it into efficiency
    change (EC) and technical change (TC), with optional further
    decomposition into pure efficiency change (PECH) and scale
    efficiency change (SECH) based on `type2` (Färe et al., 1994; Ray &
    Desli, 1997). If EC is unavailable (e.g., under `rts = "vrs"` and
    `type1 = "glob"`), it is computed as `ec = mi / tc` to ensure
    consistent output.

**Orientation**:

- Input-oriented (`"io"`): Minimizes inputs while maintaining outputs.
  Efficiency scores are in \\(0, 1\]\\ (\\\theta \leq 1\\ for radial
  models, \\\rho\\ or \\\delta \leq 1\\ for SBM models).

- Output-oriented (`"oo"`): Maximizes outputs for given inputs.
  Efficiency scores are in \\(0, 1\]\\. Radial models output \\\eta \geq
  1\\, converted to \\1/\eta\\; SBM models output \\1/\rho^\*\\ or
  \\1/\delta\\.

**Returns to Scale (RTS)**:

- CRS (`"crs"`): Assumes constant returns to scale, suitable for
  long-run analysis.

- VRS (`"vrs"`): Allows variable returns to scale, with increasing
  (`"irs"`) or decreasing (`"drs"`) returns determined by the sum of
  intensity variables (\\\lambda\\).

**Undesirable Outputs**:

- Supported in `basic_DEA`, `super_DEA`, and `basic_SBM` using
  directional distance functions (DDF) with direction vector \\(g_y,
  -g_b)\\ to increase desirable outputs and decrease undesirable outputs
  (Färe & Grosskopf, 2004). Not supported in `super_SBM` or `malmquist`.

**Malmquist-Specific Parameters**:

- `type1`: Defines the reference technology for the Malmquist index:

  - `"cont"`: Contemporary technology, using each period's frontier.

  - `"seq"`: Sequential technology, incorporating all prior periods.

  - `"glob"`: Global technology, using a single frontier across all
    periods.

- `type2`: Specifies the decomposition method:

  - `"fgnz"`: Färe et al. (1994) decomposition.

  - `"rd"`: Ray and Desli (1997) decomposition.

  - `"gl"`: Generalized decomposition.

  - `"bias"`: Bias-corrected decomposition.

**Handling NA Values**:

- Super-efficiency models (`super_DEA`, `super_SBM`) may return `NA` for
  efficient DMUs due to infeasible linear programming solutions,
  especially under VRS (Andersen & Petersen, 1993). Users can replace
  `NA` with standard efficiency scores or exclude affected DMUs.

- For `malmquist`, if `ec` is `NULL` (e.g., under `rts = "vrs"` and
  `type1 = "glob"`), it is computed as `mi / tc` to ensure a complete
  result.

The package leverages deaR for robust computation, handling zero values
internally and ensuring compatibility with input/output specifications.
Efficiency scores are standardized to \\(0, 1\]\\ for DEA models, and
Malmquist results are formatted as a data frame for easy analysis.

## Examples

``` r
# Sample DEA data
data_dea = data.frame(
  DMU = paste0("DMU", 1:5),
  input1 = c(10, 20, 15, 25, 30),
  input2 = c(5, 8, 7, 10, 12),
  output = c(100, 150, 120, 180, 200),
  ud_output = c(10, 15, 12, 20, 25)
)

# Standard DEA
result = basic_DEA(data_dea, inputs = 2:3, outputs = 4)
#> Warning: We note that we call 'targets' to the 'efficient projections' in the
#>           strongly efficient frontier.
result$efficiencies
#>      DMU1      DMU2      DMU3      DMU4      DMU5 
#> 1.0000000 1.0000000 0.9166667 1.0000000 1.0000000 

# DEA with undesirable outputs
result = basic_DEA(data_dea, inputs = 2:3, outputs = 4:5, ud_outputs = 2)
#> Warning: Undesirable (bad) outputs with no output-oriented model.
#> Warning: We note that we call 'targets' to the 'efficient projections' in the
#>           strongly efficient frontier.
result$efficiencies
#>      DMU1      DMU2      DMU3      DMU4      DMU5 
#> 1.0000000 1.0000000 0.9333333 1.0000000 1.0000000 

# Super-efficiency DEA
result = super_DEA(data_dea, inputs = 2:3, outputs = 4)
#> Warning: We note that we call 'targets' to the 'efficient projections' in the
#>           strongly efficient frontier.
result$efficiencies
#>      DMU1      DMU2      DMU3      DMU4      DMU5 
#> 1.5000000 1.0156250 0.9166667 1.0400000        NA 

# Standard SBM
result = basic_SBM(data_dea, inputs = 2:3, outputs = 4)
#> Warning: We note that we call 'targets' to the 'efficient projections' in the
#>           strongly efficient frontier.
result$efficiencies
#>      DMU1      DMU2      DMU3      DMU4      DMU5 
#> 1.0000000 1.0000000 0.9047619 1.0000000 1.0000000 

# Super-efficiency SBM
result = super_SBM(data_dea, inputs = 2:3, outputs = 4)
#> Warning: For oriented models and non constant returns to scale, feasibility is not
#>              guaranteed. Proceed with caution, some DMUs results may be missing!
#> Warning: We note that we call 'targets' to the 'efficient projections' in the
#>           strongly efficient frontier.
result$efficiencies
#>     DMU1     DMU2     DMU3     DMU4     DMU5 
#> 1.450000 1.007812 1.000000 1.040000       NA 

# Sample Malmquist data (long format)
data_malm = data.frame(
  DMU = rep(paste0("DMU", 1:5), 3),
  Period = rep(1:3, each = 5),
  input1 = c(10, 20, 15, 25, 30, 12, 22, 17, 27, 32, 14, 24, 19, 29, 34),
  input2 = c(5, 8, 7, 10, 12, 6, 9, 8, 11, 13, 7, 10, 9, 12, 14),
  output = c(100, 150, 120, 180, 200, 110, 160, 130, 190, 210, 120, 170, 140, 200, 220)
)
malmquist(data_malm, period = 2, inputs = 3:4, outputs = 5)
#> # A tibble: 10 × 7
#>    Period DMU      mi    ec    tc  pech  sech
#>    <chr>  <chr> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1 1~2    DMU1  0.917 0.922 0.994  1    0.922
#>  2 1~2    DMU2  0.948 0.972 0.976  1    0.972
#>  3 1~2    DMU3  0.948 0.953 0.995  1.00 0.949
#>  4 1~2    DMU4  0.949 0.949 1      1    0.949
#>  5 1~2    DMU5  0.923 0.923 1      1    0.923
#>  6 2~3    DMU1  0.936 0.940 0.995  1    0.940
#>  7 2~3    DMU2  0.976 0.979 0.998  1    0.979
#>  8 2~3    DMU3  0.958 0.962 0.996  1.00 0.958
#>  9 2~3    DMU4  0.955 0.955 1      1    0.955
#> 10 2~3    DMU5  0.929 0.929 1      1    0.929
```
