# Inequality Indices

Computes inequality indices for individual or grouped data: `gini0`
calculates the Gini coefficient for individual sample data. `gini`
calculates the Gini coefficient for grouped data using income and
population shares. `theil0` calculates the Theil index for individual
sample data. `theil` calculates the Theil index for grouped average
data. `theil0_g` calculates the Theil index and decomposition for
grouped sample data. `theil_g` calculates the Theil index and
decomposition for grouped average data. `theil_g2_cross` calculates the
Theil index and decomposition for two-level cross-grouped average data.
`theil_g2_nest` calculates the Theil index and decomposition for
two-level nested grouped average data.

For `theil`, `theil_g`, `theil_g2_cross`, `theil_g2_nest`: Name of
population variable (character).

## Usage

``` r
gini0(y)

gini(y, pop)

theil0(y)

theil(y, pop)

theil0_g(data, group, y)

theil_g(data, group, y, pop)

theil_g2_cross(data, group1, group2, y, pop)

theil_g2_nest(data, group1, group2, y, pop)
```

## Arguments

- y:

  For `gini0`, `gini`, `theil0`: Numeric vector of individual incomes.

- pop:

  For `gini`: Numeric vector of group populations or population shares.
  For `theil`, `theil0_g`, `theil_g`, `theil_g2_cross`, `theil_g2_nest`:
  Name of income variable (character).

- data:

  For `theil0_g`, `theil_g`, `theil_g2_cross`, `theil_g2_nest`: Data
  frame containing variables.

- group:

  For `theil0_g`, `theil_g`: Name of grouping variable (e.g., province).

- group1:

  For `theil_g2_cross`, `theil_g2_nest`: Name of first grouping variable
  (e.g., region or province).

- group2:

  For `theil_g2_cross`, `theil_g2_nest`: Name of second grouping
  variable (e.g., type or city).

## Value

For `gini0`, `gini`: Numeric Gini coefficient (0 to 1). For `theil0`,
`theil`: Numeric Theil index. For `theil0_g`, `theil_g`,
`theil_g2_cross`, `theil_g2_nest`: List with two vectors:

- `theil`: theil (Theil index and its decomposition),

- `ratio`: ratio (contribution rates of each component).

## Examples

``` r
# Sample data
income = c(10, 20, 30, 40, 100)
pop = c(100, 150, 200, 250, 300)

# Gini coefficient (individual data)
gini0(income)
#> [1] 0.4

# Gini coefficient (grouped data)
gini(income, pop)
#> [1] 0.225

# Theil index (individual sample)
data = data.frame(g = c("A","A",rep("B",10),rep("A",6)),
                  y = c(10,10,rep(8,4),rep(6,6),rep(4,4),2,2))
theil0(data$y)
#> [1] 0.07907822

# Theil index (grouped average)
data2 = data |> dplyr::count(g, y, name = "pop")
theil(data2$y, data2$pop)
#> [1] 0.07907822

# Theil index with grouping (sample data)
theil0_g(data, "g", "y")
#> $theil
#>          T         Tb         Tw          A          B 
#> 0.07907822 0.01127992 0.06779830 0.16568710 0.01021666 
#> 
#> $ratio
#>          T         Tb         Tw          A          B 
#> 1.00000000 0.14264257 0.85735743 0.77601127 0.08134615 
#> 

# Theil index with grouping (average data)
theil_g(data2, "g", "y", "pop")
#> $theil
#>          T         Tb         Tw          A          B 
#> 0.07907822 0.01127992 0.06779830 0.16568710 0.01021666 
#> 
#> $ratio
#>          T         Tb         Tw          A          B 
#> 1.00000000 0.14264257 0.85735743 0.77601127 0.08134615 
#> 

# Theil index with two-level cross-grouping
data3 = data.frame(
  industry = c("A", "A", "A", "A", "B", "B", "B", "B"),
  area = c("East", "East", "West", "West", "East", "East", "West", "West"),
  province = c("Shanghai", "Beijing", "Sichuan", "Yunnan", "Shanghai", "Beijing", "Sichuan", "Yunnan"),
  avg_wage = c(500, 400, 80, 60, 300, 250, 50, 40),
  emp_num = c(100, 80, 90, 70, 120, 100, 80, 60)
)
theil_g2_cross(data3, "industry", "area", "avg_wage", "emp_num")
#> $theil
#>           T          Tb          Tw        Tw_b        Tw_w           A 
#> 0.279831207 0.018158987 0.261672220 0.255126585 0.006545635 0.287178079 
#>           B 
#> 0.226327305 
#> 
#> $ratio
#>          T         Tb         Tw       Tw_b       Tw_w          A          B 
#> 1.00000000 0.06489265 0.93510735 0.91171598 0.02339137 0.59609568 0.33901168 
#> 

# Theil index with two-level nested grouping
data4 = data.frame(
  province = c("A", "A", "A", "A", "B", "B"),
  city = c("A1", "A1", "A2", "A2", "B1", "B1"),
  industry = c("Manu", "Serv", "Manu", "Serv", "Manu", "Serv"),
  y = c(50000, 45000, 60000, 55000, 70000, 65000),
  pop = c(10000, 8000, 15000, 12000, 10000, 8000)
)
theil_g2_nest(data4, "province", "city", "y", "pop")
#> $theil
#>            T        Tb_g1        Tb_g2           Tw            A            B 
#> 0.0095463309 0.0058045406 0.0027974591 0.0009443311 0.0042077485 0.0000000000 
#>         A_A1         A_A2         B_B1 
#> 0.0013579722 0.0009278223 0.0006738577 
#> 
#> $ratio
#>          T      Tb_g1      Tb_g2         Tw          A          B       A_A1 
#> 1.00000000 0.60803891 0.29304024 0.09892085 0.29304024 0.00000000 0.03360868 
#>       A_A2       B_B1          A          B 
#> 0.04165351 0.02365866 0.09892085 0.09892085 
#> 
```
