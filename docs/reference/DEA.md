# DEA efficiency analysis

DEA efficiency analysis

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

super_DEA(data, inputs, outputs, orientation = "io", rts = "vrs")

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

  A data frame. 1st column = DMU names.

- inputs:

  Column indices/names of input variables.

- outputs:

  Column indices/names of output variables.

- ud_outputs:

  Column indices (within `outputs`) of undesirable outputs, or NULL.

- orientation:

  `"io"` (input-oriented, default) or `"oo"` (output-oriented). All DEA
  functions return Farrell efficiency (0~1), where 1 = efficient. For
  output orientation: value = 1/\\\phi\\, where \\\phi \ge 1\\ is the
  Shephard distance.

- rts:

  `"vrs"` (variable returns to scale, default) or `"crs"` (constant).

- x0, y0:

  Used by Malmquist cross-period evaluation (not for direct calls).

## Value

A list with components: `efficiencies`, `lambdas`, `slacks`, `targets`,
`returns`, `model`, `orientation`, `dmu`.

## Functions

- `basic_DEA()`: Standard radial DEA (CCR/BCC)

- `super_DEA()`: Super-efficiency DEA (self excluded, radial only)

- `basic_SBM()`: Standard Slacks-Based Measure (SBM, Tone 2001)

- `super_SBM()`: Super-efficiency SBM (self excluded, no undesirable
  outputs)

- `malmquist()`: Malmquist productivity index

## Examples

``` r
df = data.frame(
  DMU = paste0("DMU", 1:7),
  x1  = c(20, 60, 40, 60, 70, 30, 50),
  x2  = c(151, 200, 120, 170, 250, 210, 90),
  y1  = c(100, 210, 150, 240, 220, 80, 200)
)
basic_DEA(df, inputs = 2:3, outputs = 4, rts = "crs")
#> $efficiencies
#>      DMU1      DMU2      DMU3      DMU4      DMU5      DMU6      DMU7 
#> 1.0000000 0.8203125 0.8910124 0.9570042 0.7295271 0.5437352 1.0000000 
#> 
#> $lambdas
#>           DMU1 DMU2 DMU3 DMU4 DMU5 DMU6       DMU7
#> DMU1 1.0000000    0    0    0    0    0 0.00000000
#> DMU2 0.6562500    0    0    0    0    0 0.72187500
#> DMU3 0.3719008    0    0    0    0    0 0.56404959
#> DMU4 0.5159501    0    0    0    0    0 0.94202497
#> DMU5 0.7866205    0    0    0    0    0 0.70668973
#> DMU6 0.7375887    0    0    0    0    0 0.03120567
#> DMU7 0.0000000    0    0    0    0    0 1.00000000
#> 
#> $slacks
#> $slacks$inputs
#>      x1 x2
#> DMU1  0  0
#> DMU2  0  0
#> DMU3  0  0
#> DMU4  0  0
#> DMU5  0  0
#> DMU6  0  0
#> DMU7  0  0
#> 
#> $slacks$outputs
#>                y1
#> DMU1 5.542233e-13
#> DMU2 7.673862e-13
#> DMU3 5.684342e-13
#> DMU4 1.080025e-12
#> DMU5 3.694822e-13
#> DMU6 3.836931e-13
#> DMU7 5.400125e-13
#> 
#> 
#> $targets
#> $targets$inputs
#>    [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> x1   20   60   40   60   70   30   50
#> x2  151  200  120  170  250  210   90
#> 
#> $targets$outputs
#>    [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> y1  100  210  150  240  220   80  200
#> 
#> 
#> $returns
#>  DMU1  DMU2  DMU3  DMU4  DMU5  DMU6  DMU7 
#> "CRS" "CRS" "CRS" "CRS" "CRS" "CRS" "CRS" 
#> 
#> $model
#> [1] "DEA"
#> 
#> $orientation
#> [1] "io"
#> 
#> $dmu
#> [1] "DMU1" "DMU2" "DMU3" "DMU4" "DMU5" "DMU6" "DMU7"
#> 
basic_DEA(df, inputs = 2:3, outputs = 4, rts = "vrs")
#> $efficiencies
#>      DMU1      DMU2      DMU3      DMU4      DMU5      DMU6      DMU7 
#> 1.0000000 0.8571429 0.9519868 1.0000000 0.7755102 0.7072571 1.0000000 
#> 
#> $lambdas
#>           DMU1 DMU2 DMU3      DMU4 DMU5 DMU6       DMU7
#> DMU1 1.0000000    0    0 0.0000000    0    0 0.00000000
#> DMU2 0.2142857    0    0 0.7857143    0    0 0.00000000
#> DMU3 0.3973510    0    0 0.0000000    0    0 0.60264901
#> DMU4 0.0000000    0    0 1.0000000    0    0 0.00000000
#> DMU5 0.1428571    0    0 0.8571429    0    0 0.00000000
#> DMU6 0.9594096    0    0 0.0000000    0    0 0.04059041
#> DMU7 0.0000000    0    0 0.0000000    0    0 1.00000000
#> 
#> $slacks
#> $slacks$inputs
#>                x1           x2
#> DMU1 3.552714e-15 2.842171e-14
#> DMU2 7.105427e-15 5.500000e+00
#> DMU3 0.000000e+00 1.421085e-14
#> DMU4 0.000000e+00 0.000000e+00
#> DMU5 7.105427e-15 2.659184e+01
#> DMU6 0.000000e+00 0.000000e+00
#> DMU7 0.000000e+00 0.000000e+00
#> 
#> $slacks$outputs
#>                y1
#> DMU1 0.000000e+00
#> DMU2 0.000000e+00
#> DMU3 1.026490e+01
#> DMU4 1.136868e-13
#> DMU5 0.000000e+00
#> DMU6 2.405904e+01
#> DMU7 0.000000e+00
#> 
#> 
#> $targets
#> $targets$inputs
#>    [,1]  [,2] [,3] [,4]     [,5] [,6] [,7]
#> x1   20  60.0   40   60  70.0000   30   50
#> x2  151 194.5  120  170 223.4082  210   90
#> 
#> $targets$outputs
#>    [,1] [,2]     [,3] [,4] [,5]    [,6] [,7]
#> y1  100  210 160.2649  240  220 104.059  200
#> 
#> 
#> $returns
#>  DMU1  DMU2  DMU3  DMU4  DMU5  DMU6  DMU7 
#> "VRS" "VRS" "VRS" "VRS" "VRS" "VRS" "VRS" 
#> 
#> $model
#> [1] "DEA"
#> 
#> $orientation
#> [1] "io"
#> 
#> $dmu
#> [1] "DMU1" "DMU2" "DMU3" "DMU4" "DMU5" "DMU6" "DMU7"
#> 
df = data.frame(
  DMU = paste0("DMU", 1:7),
  x1  = c(20, 60, 40, 60, 70, 30, 50),
  x2  = c(151, 200, 120, 170, 250, 210, 90),
  y1  = c(100, 210, 150, 240, 220, 80, 200)
)
super_DEA(df, inputs = 2:3, outputs = 4, rts = "crs")
#> $efficiencies
#>      DMU1      DMU2      DMU3      DMU4      DMU5      DMU6      DMU7 
#> 1.2500000 0.8203125 0.8910124 0.9570042 0.7295271 0.5437352 1.5740741 
#> 
#> $lambdas
#>           DMU1 DMU2 DMU3      DMU4 DMU5 DMU6       DMU7
#> DMU1        NA    0    0 0.4166667    0    0 0.00000000
#> DMU2 0.6562500   NA    0 0.0000000    0    0 0.72187500
#> DMU3 0.3719008    0   NA 0.0000000    0    0 0.56404959
#> DMU4 0.5159501    0    0        NA    0    0 0.94202497
#> DMU5 0.7866205    0    0 0.0000000   NA    0 0.70668973
#> DMU6 0.7375887    0    0 0.0000000    0   NA 0.03120567
#> DMU7 0.0000000    0    0 0.8333333    0    0         NA
#> 
#> $slacks
#> $slacks$inputs
#>                x1           x2
#> DMU1 3.552714e-15 1.179167e+02
#> DMU2 0.000000e+00 0.000000e+00
#> DMU3 0.000000e+00 0.000000e+00
#> DMU4 0.000000e+00 2.842171e-14
#> DMU5 0.000000e+00 0.000000e+00
#> DMU6 0.000000e+00 0.000000e+00
#> DMU7 2.870370e+01 4.547474e-13
#> 
#> $slacks$outputs
#>                y1
#> DMU1 0.000000e+00
#> DMU2 7.673862e-13
#> DMU3 5.684342e-13
#> DMU4 1.136868e-13
#> DMU5 3.694822e-13
#> DMU6 0.000000e+00
#> DMU7 0.000000e+00
#> 
#> 
#> $targets
#> $targets$inputs
#>        [,1] [,2] [,3] [,4] [,5] [,6]    [,7]
#> x1 20.00000   60   40   60   70   30 21.2963
#> x2 33.08333  200  120  170  250  210 90.0000
#> 
#> $targets$outputs
#>    [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> y1  100  210  150  240  220   80  200
#> 
#> 
#> $returns
#>  DMU1  DMU2  DMU3  DMU4  DMU5  DMU6  DMU7 
#> "CRS" "CRS" "CRS" "CRS" "CRS" "CRS" "CRS" 
#> 
#> $model
#> [1] "DEA"
#> 
#> $orientation
#> [1] "io"
#> 
#> $dmu
#> [1] "DMU1" "DMU2" "DMU3" "DMU4" "DMU5" "DMU6" "DMU7"
#> 
df = data.frame(
  DMU = paste0("DMU", 1:7),
  x1  = c(20, 60, 40, 60, 70, 30, 50),
  x2  = c(151, 200, 120, 170, 250, 210, 90),
  y1  = c(100, 210, 150, 240, 220, 80, 200)
)
basic_SBM(df, inputs = 2:3, outputs = 4, rts = "crs")
#> $efficiencies
#>      DMU1      DMU2      DMU3      DMU4      DMU5      DMU6      DMU7 
#> 1.0000000 0.6737500 0.7500000 0.8176471 0.5908571 0.4190476 1.0000000 
#> 
#> $lambdas
#>      DMU1 DMU2 DMU3 DMU4 DMU5 DMU6 DMU7
#> DMU1    1    0    0    0    0    0  0.0
#> DMU2    0    0    0    0    0    0  1.2
#> DMU3    0    0    0    0    0    0  0.8
#> DMU4    0    0    0    0    0    0  1.2
#> DMU5    0    0    0    0    0    0  1.4
#> DMU6    0    0    0    0    0    0  0.6
#> DMU7    0    0    0    0    0    0  1.0
#> 
#> $slacks
#> $slacks$inputs
#>      x1  x2
#> DMU1  0   0
#> DMU2  0  92
#> DMU3  0  48
#> DMU4  0  62
#> DMU5  0 124
#> DMU6  0 156
#> DMU7  0   0
#> 
#> $slacks$outputs
#>      y1
#> DMU1  0
#> DMU2 30
#> DMU3 10
#> DMU4  0
#> DMU5 60
#> DMU6 40
#> DMU7  0
#> 
#> 
#> $targets
#> $targets$inputs
#>    [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> x1   20   60   40   60   70   30   50
#> x2  151  108   72  108  126   54   90
#> 
#> $targets$outputs
#>    [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> y1  100  240  160  240  280  120  200
#> 
#> 
#> $returns
#>  DMU1  DMU2  DMU3  DMU4  DMU5  DMU6  DMU7 
#> "CRS" "CRS" "CRS" "CRS" "CRS" "CRS" "CRS" 
#> 
#> $model
#> [1] "DEA"
#> 
#> $orientation
#> [1] "io"
#> 
#> $dmu
#> [1] "DMU1" "DMU2" "DMU3" "DMU4" "DMU5" "DMU6" "DMU7"
#> 
df = data.frame(
  DMU = paste0("DMU", 1:7),
  x1  = c(20, 60, 40, 60, 70, 30, 50),
  x2  = c(151, 200, 120, 170, 250, 210, 90),
  y1  = c(100, 210, 150, 240, 220, 80, 200)
)
super_SBM(df, inputs = 2:3, outputs = 4, rts = "crs")
#> $efficiencies
#>      DMU1      DMU2      DMU3      DMU4      DMU5      DMU6      DMU7 
#>        NA 0.6737500 0.7500000 0.8176471 0.5908571 0.4190476        NA 
#> 
#> $lambdas
#>      DMU1 DMU2 DMU3 DMU4 DMU5 DMU6 DMU7
#> DMU1   NA   NA   NA   NA   NA   NA   NA
#> DMU2    0   NA    0    0    0    0  1.2
#> DMU3    0    0   NA    0    0    0  0.8
#> DMU4    0    0    0   NA    0    0  1.2
#> DMU5    0    0    0    0   NA    0  1.4
#> DMU6    0    0    0    0    0   NA  0.6
#> DMU7   NA   NA   NA   NA   NA   NA   NA
#> 
#> $slacks
#> $slacks$inputs
#>      x1  x2
#> DMU1 NA  NA
#> DMU2  0  92
#> DMU3  0  48
#> DMU4  0  62
#> DMU5  0 124
#> DMU6  0 156
#> DMU7 NA  NA
#> 
#> $slacks$outputs
#>      y1
#> DMU1 NA
#> DMU2 30
#> DMU3 10
#> DMU4  0
#> DMU5 60
#> DMU6 40
#> DMU7 NA
#> 
#> 
#> $targets
#> $targets$inputs
#>    [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> x1   NA   60   40   60   70   30   NA
#> x2   NA  108   72  108  126   54   NA
#> 
#> $targets$outputs
#>    [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> y1   NA  240  160  240  280  120   NA
#> 
#> 
#> $returns
#>  DMU1  DMU2  DMU3  DMU4  DMU5  DMU6  DMU7 
#> "CRS" "CRS" "CRS" "CRS" "CRS" "CRS" "CRS" 
#> 
#> $model
#> [1] "DEA"
#> 
#> $orientation
#> [1] "io"
#> 
#> $dmu
#> [1] "DMU1" "DMU2" "DMU3" "DMU4" "DMU5" "DMU6" "DMU7"
#> 
panel = data.frame(
  DMU    = rep(paste0("DMU", 1:5), 3),
  Period = rep(1:3, each = 5),
  x1     = c(10, 20, 15, 25, 30, 12, 22, 17, 27, 32, 14, 24, 19, 29, 34),
  y1     = c(100, 150, 120, 180, 200, 110, 160, 130, 190, 210, 120, 170, 140, 200, 220)
)
malmquist(panel, period = "Period", inputs = 3, outputs = 4,
  rts = "crs", type1 = "cont", type2 = "fgnz")
#> # A tibble: 10 × 7
#>    Period DMU      mi    ec    tc  pech  sech
#>    <chr>  <chr> <dbl> <dbl> <dbl> <dbl> <dbl>
#>  1 1~2    DMU1  0.917  1    0.917    NA    NA
#>  2 1~2    DMU2  0.970  1.06 0.917    NA    NA
#>  3 1~2    DMU3  0.956  1.04 0.917    NA    NA
#>  4 1~2    DMU4  0.977  1.07 0.917    NA    NA
#>  5 1~2    DMU5  0.984  1.07 0.917    NA    NA
#>  6 2~3    DMU1  0.935  1    0.935    NA    NA
#>  7 2~3    DMU2  0.974  1.04 0.935    NA    NA
#>  8 2~3    DMU3  0.964  1.03 0.935    NA    NA
#>  9 2~3    DMU4  0.980  1.05 0.935    NA    NA
#> 10 2~3    DMU5  0.986  1.05 0.935    NA    NA
```
