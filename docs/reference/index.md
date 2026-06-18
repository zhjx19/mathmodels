# Package index

## Data Preprocessing

Functions for scaling, normalizing, and standardizing indicator data.

- [`standardize()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`normalize()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`rescale()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`rescale_middle()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`rescale_interval()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`rescale_extreme()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`rescale_initial()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`rescale_mean()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  [`to_positive()`](https://zhjx19.github.io/mathmodels/reference/preprocess.md)
  : Preprocessing Functions for Data Normalization and Standardization

## Subjective Weighting

Analytic Hierarchy Process (AHP).

- [`AHP()`](https://zhjx19.github.io/mathmodels/reference/AHP.md) : AHP:
  Analytic Hierarchy Process

## Objective Weighting

Data-driven methods for determining indicator weights.

- [`entropy_weight()`](https://zhjx19.github.io/mathmodels/reference/entropy_weight.md)
  : Entropy Weight Method
- [`critic_weight()`](https://zhjx19.github.io/mathmodels/reference/critic_weight.md)
  : CRITIC Weight Method
- [`pca_weight()`](https://zhjx19.github.io/mathmodels/reference/pca_weight.md)
  : PCA-Based Weighting Method
- [`cv_weight()`](https://zhjx19.github.io/mathmodels/reference/cv_weight.md)
  : Coefficient of Variation Weighting

## Weight Combination

Methods for combining subjective and objective weights.

- [`combine_weights()`](https://zhjx19.github.io/mathmodels/reference/combine_weights.md)
  : Combine Subjective and Objective Weights

## Comprehensive Evaluation

Multi-criteria decision making and evaluation methods.

- [`topsis()`](https://zhjx19.github.io/mathmodels/reference/topsis.md)
  : TOPSIS Method for Multi-Criteria Decision Making
- [`linear_sum()`](https://zhjx19.github.io/mathmodels/reference/linear_sum.md)
  : Linear weighted synthesis
- [`fuzzy_eval()`](https://zhjx19.github.io/mathmodels/reference/fuzzy_eval.md)
  : Fuzzy Comprehensive Evaluation
- [`rank_sum_ratio()`](https://zhjx19.github.io/mathmodels/reference/rank_sum_ratio.md)
  : Rank Sum Ratio (RSR) Evaluation

## Grey Relational Analysis

Grey system theory for relational analysis and TOPSIS.

- [`grey_corr()`](https://zhjx19.github.io/mathmodels/reference/grey_analysis.md)
  [`grey_corr_topsis()`](https://zhjx19.github.io/mathmodels/reference/grey_analysis.md)
  : Grey Relational Analysis Functions

## Fuzzy Sets & Membership Functions

Fuzzy set membership functions, computation, and defuzzification.

- [`compute_mf_funs()`](https://zhjx19.github.io/mathmodels/reference/compute_mf.md)
  [`compute_mf()`](https://zhjx19.github.io/mathmodels/reference/compute_mf.md)
  : Compute fuzzy membership vector and return corresponding membership
  functions.
- [`defuzzify()`](https://zhjx19.github.io/mathmodels/reference/defuzzify.md)
  : Defuzzification Methods for Fuzzy Comprehensive Evaluation
- [`tri_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`trap_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`gauss_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`gbell_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`gauss2mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`sigmoid_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`dsigmoid_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`psigmoid_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`z_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`pi_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`s_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  [`plot_mf()`](https://zhjx19.github.io/mathmodels/reference/membership.md)
  : Membership Functions for Fuzzy Logic

## Data Envelopment Analysis (DEA)

CCR, BCC, SBM, and Malmquist productivity index.

- [`basic_DEA()`](https://zhjx19.github.io/mathmodels/reference/DEA.md)
  [`super_DEA()`](https://zhjx19.github.io/mathmodels/reference/DEA.md)
  [`basic_SBM()`](https://zhjx19.github.io/mathmodels/reference/DEA.md)
  [`super_SBM()`](https://zhjx19.github.io/mathmodels/reference/DEA.md)
  [`malmquist()`](https://zhjx19.github.io/mathmodels/reference/DEA.md)
  : DEA efficiency analysis

## Inequality Measures

Gini coefficient, Theil index, and their decompositions.

- [`gini0()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  [`gini()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  [`theil0()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  [`theil()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  [`theil0_g()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  [`theil_g()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  [`theil_g2_cross()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  [`theil_g2_nest()`](https://zhjx19.github.io/mathmodels/reference/inequality.md)
  : Inequality Indices

## Regional Economics

Location Quotient, Herfindahl-Hirschman Index, and Ellison-Glaeser
Index.

- [`LQ()`](https://zhjx19.github.io/mathmodels/reference/regional_economics.md)
  [`HHI()`](https://zhjx19.github.io/mathmodels/reference/regional_economics.md)
  [`EG()`](https://zhjx19.github.io/mathmodels/reference/regional_economics.md)
  : Regional Economics Functions

## Grey Prediction Models

GM(1,1), GM(1,N), DGM(2,1), and Verhulst models.

- [`GM11()`](https://zhjx19.github.io/mathmodels/reference/grey_models.md)
  [`GM1N()`](https://zhjx19.github.io/mathmodels/reference/grey_models.md)
  [`DGM21()`](https://zhjx19.github.io/mathmodels/reference/grey_models.md)
  [`verhulst()`](https://zhjx19.github.io/mathmodels/reference/grey_models.md)
  : Grey Prediction Models

## System Evaluation

Coupling coordination degree and obstacle degree models.

- [`coupling_degree()`](https://zhjx19.github.io/mathmodels/reference/system_evaluation.md)
  [`obstacle_degree()`](https://zhjx19.github.io/mathmodels/reference/system_evaluation.md)
  : System Evaluation Functions for Coupling and Obstacle Analysis

## Datasets

Built-in example datasets for demonstration and testing.

- [`water_quality`](https://zhjx19.github.io/mathmodels/reference/water_quality.md)
  : Water Quality Dataset

## Utilities

Helper functions for predictions, data I/O, and other tasks.

- [`combine_preds()`](https://zhjx19.github.io/mathmodels/reference/combine_preds.md)
  : Combine Multiple Prediction Results
- [`read_nbs()`](https://zhjx19.github.io/mathmodels/reference/read_nbs.md)
  : Read and Combine National Bureau of Statistics XLS Files
