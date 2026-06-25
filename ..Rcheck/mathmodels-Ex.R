pkgname <- "mathmodels"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('mathmodels')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("AHP")
### * AHP

flush(stderr()); flush(stdout())

### Name: AHP
### Title: AHP: Analytic Hierarchy Process
### Aliases: AHP

### ** Examples

A = matrix(c(1,   1/2, 4, 3,   3,
             2,   1,   7, 5,   5,
             1/4, 1/7, 1, 1/2, 1/3,
             1/3, 1/5, 2, 1,   1,
             1/3, 1/5, 3, 1,   1), byrow = TRUE, nrow = 5)
AHP(A)



cleanEx()
nameEx("DEA")
### * DEA

flush(stderr()); flush(stdout())

### Name: DEA
### Title: DEA efficiency analysis
### Aliases: DEA basic_DEA super_DEA basic_SBM super_SBM malmquist

### ** Examples

df = data.frame(
  DMU = paste0("DMU", 1:7),
  x1  = c(20, 60, 40, 60, 70, 30, 50),
  x2  = c(151, 200, 120, 170, 250, 210, 90),
  y1  = c(100, 210, 150, 240, 220, 80, 200)
)
basic_DEA(df, inputs = 2:3, outputs = 4, rts = "crs")
basic_DEA(df, inputs = 2:3, outputs = 4, rts = "vrs")
df = data.frame(
  DMU = paste0("DMU", 1:7),
  x1  = c(20, 60, 40, 60, 70, 30, 50),
  x2  = c(151, 200, 120, 170, 250, 210, 90),
  y1  = c(100, 210, 150, 240, 220, 80, 200)
)
super_DEA(df, inputs = 2:3, outputs = 4, rts = "crs")
df = data.frame(
  DMU = paste0("DMU", 1:7),
  x1  = c(20, 60, 40, 60, 70, 30, 50),
  x2  = c(151, 200, 120, 170, 250, 210, 90),
  y1  = c(100, 210, 150, 240, 220, 80, 200)
)
basic_SBM(df, inputs = 2:3, outputs = 4, rts = "crs")
df = data.frame(
  DMU = paste0("DMU", 1:7),
  x1  = c(20, 60, 40, 60, 70, 30, 50),
  x2  = c(151, 200, 120, 170, 250, 210, 90),
  y1  = c(100, 210, 150, 240, 220, 80, 200)
)
super_SBM(df, inputs = 2:3, outputs = 4, rts = "crs")
panel = data.frame(
  DMU    = rep(paste0("DMU", 1:5), 3),
  Period = rep(1:3, each = 5),
  x1     = c(10, 20, 15, 25, 30, 12, 22, 17, 27, 32, 14, 24, 19, 29, 34),
  y1     = c(100, 150, 120, 180, 200, 110, 160, 130, 190, 210, 120, 170, 140, 200, 220)
)
malmquist(panel, period = "Period", inputs = 3, outputs = 4,
  rts = "crs", type1 = "cont", type2 = "fgnz")



cleanEx()
nameEx("combine_preds")
### * combine_preds

flush(stderr()); flush(stdout())

### Name: combine_preds
### Title: Combine Multiple Prediction Results
### Aliases: combine_preds

### ** Examples

# Example: Combine three prediction results
preds = c(100, 102, 98)  # E.g., from grey prediction, ARIMA, or ML models
combine_preds(preds)



cleanEx()
nameEx("combine_weights")
### * combine_weights

flush(stderr()); flush(stdout())

### Name: combine_weights
### Title: Combine Subjective and Objective Weights
### Aliases: combine_weights

### ** Examples

w_subj = c(0.4, 0.3, 0.2, 0.1)
w_obj = c(0.25, 0.2, 0.3, 0.25)
combine_weights(w_subj, w_obj, type = "linear", alpha = 0.6)
combine_weights(w_subj, w_obj, type = "multiplicative")
combine_weights(w_subj, w_obj, type = "game")
combine_weights(w_subj, w_obj, type = "game_linear")




cleanEx()
nameEx("compute_mf")
### * compute_mf

flush(stderr()); flush(stdout())

### Name: compute_mf
### Title: Compute fuzzy membership vector and return corresponding
###   membership functions.
### Aliases: compute_mf compute_mf_funs

### ** Examples

# Example: SO2 concentration = 0.07, thresholds = c(0.05, 0.15, 0.25, 0.5)
th = c(0.05, 0.15, 0.25, 0.5)
compute_mf(0.07, th)

## Not run: 
##D mfs = compute_mf_funs(th)
##D plots = lapply(mfs, \(x) plot_mf(x, xlim = c(0, 0.6)))
##D gridExtra::grid.arrange(grobs = plots, nrow = 2)
## End(Not run)





cleanEx()
nameEx("critic_weight")
### * critic_weight

flush(stderr()); flush(stdout())

### Name: critic_weight
### Title: CRITIC Weight Method
### Aliases: critic_weight

### ** Examples

# Example: Using CRITIC method on a simple dataset
X = data.frame(
  x1 = c(3, 5, 2, 7),
  x2 = c(10, 20, 15, 25)
)
index = c("+", "-")
critic_weight(X, index)
critic_weight(X, index, method = "entropy")



cleanEx()
nameEx("cv_weight")
### * cv_weight

flush(stderr()); flush(stdout())

### Name: cv_weight
### Title: Coefficient of Variation Weighting
### Aliases: cv_weight

### ** Examples

X = data.frame(x1 = c(10, 20, 15), x2 = c(5, 10, 8))
cv_weight(X)




cleanEx()
nameEx("defuzzify")
### * defuzzify

flush(stderr()); flush(stdout())

### Name: defuzzify
### Title: Defuzzification Methods for Fuzzy Comprehensive Evaluation
### Aliases: defuzzify

### ** Examples

# Example: Defuzzify fuzzy evaluation vectors for three schemes
mu = c(0.318, 0.351, 0.203, 0.128)
scores = c(30, 60, 75, 90)  # Scores for "Poor", "Fair", "Good", "Excellent"
defuzzify(mu, scores, method = "weighted_average")
defuzzify(mu, scores, method = "max_membership")
defuzzify(mu, scores, method = "centroid")



cleanEx()
nameEx("entropy_weight")
### * entropy_weight

flush(stderr()); flush(stdout())

### Name: entropy_weight
### Title: Entropy Weight Method
### Aliases: entropy_weight

### ** Examples

X = data.frame(
  x1 = c(3, 5, 2, 7),
  x2 = c(10, 20, 15, 25)
)
index = c("+", "-")
entropy_weight(X, index)




cleanEx()
nameEx("fuzzy_eval")
### * fuzzy_eval

flush(stderr()); flush(stdout())

### Name: fuzzy_eval
### Title: Fuzzy Comprehensive Evaluation
### Aliases: fuzzy_eval

### ** Examples

w = c(0.3, 0.3, 0.3, 0.1)  # weights (e.g., from AHP or entropy)

# fuzzy evaluation matrix (3 grades for 4 factors)
R = matrix(c(0.8, 0.7, 0.6, 0.7,
             0.1, 0.2, 0.2, 0.1,
             0.1, 0.1, 0.2, 0.2), nrow = 3, byrow = TRUE)
# Apply fuzzy comprehensive evaluation
fuzzy_eval(w, R, type = 3)  # Weighted sum




cleanEx()
nameEx("grey_analysis")
### * grey_analysis

flush(stderr()); flush(stdout())

### Name: grey_analysis
### Title: Grey Relational Analysis Functions
### Aliases: grey_analysis grey_corr grey_corr_topsis

### ** Examples

# Grey correlation degree
ref = c(0.9, 0.8, 0.7)
cmp = data.frame(
  x1 = c(0.9, 0.7, 0.8),
  x2 = c(0.8, 0.9, 0.7),
  x3 = c(0.7, 0.8, 0.9)
)
grey_corr(ref, cmp, rho = 0.5)

# Grey correlation evaluation
X = data.frame(x1 = c(8, 7, 6), x2 = c(150, 180, 200), x3 = c(60, 80, 100))
w = c(0.3, 0.4, 0.3)
idx = c("+", "+", "+")
grey_corr_topsis(X, w, idx, rho = 0.5)




cleanEx()
nameEx("grey_models")
### * grey_models

flush(stderr()); flush(stdout())

### Name: grey_models
### Title: Grey Prediction Models
### Aliases: grey_models GM11 GM1N DGM21 verhulst

### ** Examples

# Sample time series for GM11, DGM21, Verhulst
x = c(100, 120, 145, 175, 210)

# GM11
result = GM11(x)
result$fitted    # Fitted values
result$pnext     # Next prediction
result$f(6:8)    # Predict next 3 periods

# DGM21
x = c(2.874,3.278,3.39,3.679,3.77,3.8)
result = DGM21(x)
result$fitted    # Fitted values
result$pnext     # Next prediction
result$f(6:8)    # Predict next 3 periods

# Verhulst
x = c(4.93,2.33,3.87,4.35,6.63,7.15,5.37,6.39,7.81,8.35)
result = verhulst(x)
result$fitted    # Fitted values
result$pnext     # Next prediction
result$f(6:8)    # Predict next 3 periods

# Sample data for GM1N
data = data.frame(
  factor1 = c(50, 55, 60, 65, 70),
  factor2 = c(20, 22, 25, 28, 30),
  output = c(100, 120, 145, 175, 210)
)
result = GM1N(data)
result$fitted




cleanEx()
nameEx("inequality")
### * inequality

flush(stderr()); flush(stdout())

### Name: inequality
### Title: Inequality Indices
### Aliases: inequality gini0 gini theil0 theil theil0_g theil_g
###   theil_g2_cross theil_g2_nest

### ** Examples

# Sample data
income = c(10, 20, 30, 40, 100)
pop = c(100, 150, 200, 250, 300)

# Gini coefficient (individual data)
gini0(income)

# Gini coefficient (grouped data)
gini(income, pop)

# Theil index (individual sample)
data = data.frame(g = c("A","A",rep("B",10),rep("A",6)),
                  y = c(10,10,rep(8,4),rep(6,6),rep(4,4),2,2))
theil0(data$y)

# Theil index (grouped average)
data2 = data |> dplyr::count(g, y, name = "pop")
theil(data2$y, data2$pop)

# Theil index with grouping (sample data)
theil0_g(data, "g", "y")

# Theil index with grouping (average data)
theil_g(data2, "g", "y", "pop")

# Theil index with two-level cross-grouping
data3 = data.frame(
  industry = c("A", "A", "A", "A", "B", "B", "B", "B"),
  area = c("East", "East", "West", "West", "East", "East", "West", "West"),
  province = c("Shanghai", "Beijing", "Sichuan", "Yunnan", "Shanghai", "Beijing", "Sichuan", "Yunnan"),
  avg_wage = c(500, 400, 80, 60, 300, 250, 50, 40),
  emp_num = c(100, 80, 90, 70, 120, 100, 80, 60)
)
theil_g2_cross(data3, "industry", "area", "avg_wage", "emp_num")

# Theil index with two-level nested grouping
data4 = data.frame(
  province = c("A", "A", "A", "A", "B", "B"),
  city = c("A1", "A1", "A2", "A2", "B1", "B1"),
  industry = c("Manu", "Serv", "Manu", "Serv", "Manu", "Serv"),
  y = c(50000, 45000, 60000, 55000, 70000, 65000),
  pop = c(10000, 8000, 15000, 12000, 10000, 8000)
)
theil_g2_nest(data4, "province", "city", "y", "pop")




cleanEx()
nameEx("linear_sum")
### * linear_sum

flush(stderr()); flush(stdout())

### Name: linear_sum
### Title: Linear weighted synthesis
### Aliases: linear_sum

### ** Examples

data = data.frame(
  GDP = c(0.85, 0.72, 0.91, NA),
  Employment = c(0.78, 0.85, 0.67, 0.73),
  Environment = c(0.65, 0.72, NA, 0.81)
)
w = c(0.5, 0.3, 0.2)
linear_sum(data, w)



cleanEx()
nameEx("markov")
### * markov

flush(stderr()); flush(stdout())

### Name: markov
### Title: Markov Chain and Grey-Markov Prediction Models
### Aliases: markov markov_chain GM11_markov

### ** Examples

# --- Markov chain prediction ---
# Weather states: Rainy, Cloudy, Sunny
S = factor(c("Sunny", "Sunny", "Cloudy", "Rainy", "Sunny",
             "Cloudy", "Sunny", "Sunny", "Rainy", "Cloudy",
             "Sunny", "Cloudy", "Rainy", "Sunny", "Sunny",
             "Cloudy", "Sunny", "Rainy", "Cloudy", "Sunny"),
            levels = c("Rainy", "Cloudy", "Sunny"))
markov_chain(S, s0 = "Cloudy", n_steps = 3)

# --- Grey-Markov prediction ---
X = c(174, 179, 183, 189, 207, 234, 220.5, 256, 270, 285)
GM11_markov(X, n_ahead = 3)




cleanEx()
nameEx("membership")
### * membership

flush(stderr()); flush(stdout())

### Name: membership
### Title: Membership Functions for Fuzzy Logic
### Aliases: membership tri_mf trap_mf gauss_mf gbell_mf gauss2mf
###   sigmoid_mf dsigmoid_mf psigmoid_mf z_mf pi_mf s_mf plot_mf

### ** Examples

# Define input values
x = 0:10

# Triangular membership
tri_mf(x, params = c(3, 6, 8))

# Trapezoidal membership
trap_mf(x, params = c(1, 5, 7, 8))

# Gaussian membership
gauss_mf(x, params = c(2, 5))

# Generalized bell membership
gbell_mf(x, params = c(2, 4, 6))

# Two-parameter Gaussian membership
gauss2mf(x, params = c(1, 3, 3, 4))

# Sigmoid membership
sigmoid_mf(x, params = c(2, 4))

# Difference of sigmoids membership
dsigmoid_y = dsigmoid_mf(x, params = c(5, 2, 5, 7))

# Product of sigmoids membership
psigmoid_mf(x, params = c(2, 3, -5, 8))

# Z-shaped membership
z_mf(x, params = c(3, 7))

# PI-shaped membership
pi_mf(x, params = c(1, 4, 5, 10))

# S-shaped membership
s_mf(x, params = c(1, 8))

## Not run: 
##D # Visualize membership functions
##D plot_mf(\(x) tri_mf(x, c(3, 6, 8)), main = "Triangular MF")
##D plot_mf(\(x) trap_mf(x, c(1, 5, 7, 8)), main = "Trapezoidal MF")
##D plot_mf(\(x) gauss_mf(x, c(2, 5)), main = "Gaussian MF")
##D plot_mf(\(x) gbell_mf(x, c(2, 4, 6)), main = "Generalized Bell MF")
##D plot_mf(\(x) gauss2mf(x, c(1, 3, 3, 4)), main = "Two-Parameter Gaussian MF")
##D plot_mf(\(x) sigmoid_mf(x, c(2, 4)), main = "Sigmoid MF")
##D plot_mf(\(x) dsigmoid_mf(x, c(5, 2, 5, 7)), main = "Difference of Sigmoids MF")
##D plot_mf(\(x) psigmoid_mf(x, c(2, 3, -5, 8)), main = "Product of Sigmoids MF")
##D plot_mf(\(x) z_mf(x, c(3, 7)), main = "Z-Shaped MF")
##D plot_mf(\(x) pi_mf(x, c(1, 4, 5, 10)), main = "PI-Shaped MF")
##D plot_mf(\(x) s_mf(x, c(1, 8)), main = "S-Shaped MF")
## End(Not run)




cleanEx()
nameEx("pca_weight")
### * pca_weight

flush(stderr()); flush(stdout())

### Name: pca_weight
### Title: PCA-Based Weighting Method
### Aliases: pca_weight

### ** Examples

# Example: Using PCA to compute indicator weights
ind = c("+","+","-","-")
pca_weight(iris[1:10, 1:4], ind, nfs = 2)



cleanEx()
nameEx("preprocess")
### * preprocess

flush(stderr()); flush(stdout())

### Name: preprocess
### Title: Preprocessing Functions for Data Normalization and
###   Standardization
### Aliases: preprocess standardize normalize rescale rescale_middle
###   rescale_interval rescale_extreme rescale_initial rescale_mean
###   to_positive

### ** Examples

# Standardization
x = c(4, 1, NA, 5, 8)
standardize(x)

# L2 norm normalization
normalize(x)

# Min-Max normalization (positive direction)
rescale(x)                # Scale to \code{[0, 1]}
rescale(x, type = "-", a = 0.002, b = 0.996)  # Reverse scaling

# Negative-to-positive transformation
to_positive(x)                       # Min-max transformation
to_positive(x, type = "reciprocal")  # Reciprocal transformation

# Centered-type normalization
PH = 6:9
rescale_middle(PH, 7)

# Interval-type normalization
Temp = c(35.2, 35.8, 36.6, 37.1, 37.8, 38.4)
rescale_interval(Temp, 36, 37)

# Extreme-based normalization
rescale_extreme(x)         # min(x)/x
rescale_extreme(x, "-")    # x/max(x)

# Initial-based normalization
rescale_initial(x)

# Mean-based normalization
rescale_mean(x)




cleanEx()
nameEx("rank_sum_ratio")
### * rank_sum_ratio

flush(stderr()); flush(stdout())

### Name: rank_sum_ratio
### Title: Rank Sum Ratio (RSR) Evaluation
### Aliases: rank_sum_ratio

### ** Examples

# Example data
data = data.frame(ID = c("A", "B", "C"), X1 = c(10, 20, 15), X2 = c(5, 10, 8))
w = c(0.4, 0.6)
rank_sum_ratio(data, w, method = "int")




cleanEx()
nameEx("read_nbs")
### * read_nbs

flush(stderr()); flush(stdout())

### Name: read_nbs
### Title: Read and Combine National Bureau of Statistics XLS Files
### Aliases: read_nbs

### ** Examples

## Not run: 
##D paths = c("file1.xls", "file2.xls")
##D data = read_nbs(paths)
## End(Not run)



cleanEx()
nameEx("regional_economics")
### * regional_economics

flush(stderr()); flush(stdout())

### Name: regional_economics
### Title: Regional Economics Functions
### Aliases: regional_economics LQ HHI EG

### ** Examples

# Example data
data = data.frame(
  region = c("National", "Region_A", "Region_B"),
  total = c(10000, 4000, 6000),
  industry1 = c(2000, 1000, 1000),
  industry2 = c(6000, 2000, 4000)
)

# Calculate Location Quotient
LQ(data, region, total, starts_with("industry"))

# Calculate HHI
x = c(50, 30, 20)
# Calculate the raw HHI
HHI(x)
# Calculate the standard (scaled) HHI
HHI(x, scaled = TRUE)

# Example data for EG
eg_data = data.frame(
 region = c("R1", "R1", "R1", "R1", "R2", "R2", "R3", "R3", "R1", "R2", "R2", "R3"),
 industry = c("A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B"),
 employment = c(250, 200, 150, 100, 20, 15, 10, 5, 50, 200, 150, 50)
)
EG(eg_data, region, industry, y = employment)




cleanEx()
nameEx("system_evaluation")
### * system_evaluation

flush(stderr()); flush(stdout())

### Name: system_evaluation
### Title: System Evaluation Functions for Coupling and Obstacle Analysis
### Aliases: system_evaluation coupling_degree obstacle_degree

### ** Examples

# Sample normalized subsystem scores
df = data.frame(
  ID = LETTERS[1:6],
  s1 = c(0.0162, 0.1782, 0.5490, 0.6730, 0.0207, 0.9875),
  s2 = c(0.2720, 0.6824, 0.0593, 0.4812, 0.8891, 0.5573),
  s3 = c(0.2655, 0.3721, 0.5729, 0.9082, 0.2017, 0.8984)
)
# Coupling Degree Analysis
coupling_degree(df, id_cols = "ID")        # Equal weights
coupling_degree(df, c(0.4, 0.3, 0.3), id_cols = "ID",
                type = "adjusted")         # "adjusted" coupling degree
# Obstacle Degree Analysis
obstacle_degree(df, id_cols = "ID")        # Equal weights
obstacle_degree(df, c(0.4, 0.3, 0.3), id_cols = "ID")




cleanEx()
nameEx("topsis")
### * topsis

flush(stderr()); flush(stdout())

### Name: topsis
### Title: TOPSIS Method for Multi-Criteria Decision Making
### Aliases: topsis

### ** Examples

A = data.frame(
  X1 = c(2, 5, 3),  # "+"
  X2 = c(8, 1, 6)   # "-"
)
w = c(0.6, 0.4)
idx = c("+","-")
topsis(A, w, idx)



cleanEx()
nameEx("water_quality")
### * water_quality

flush(stderr()); flush(stdout())

### Name: water_quality
### Title: Water Quality Dataset
### Aliases: water_quality
### Keywords: datasets

### ** Examples

# Load the dataset
data(water_quality)

# Preview the data
head(water_quality)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
