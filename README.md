
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mathmodels:  A Powerful R Toolkit for Mathematical Modeling

<!-- badges: start -->
<!-- badges: end -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/mathmodels)](https://cran.r-project.org/package=mathmodels)
[![Online Manual](https://img.shields.io/badge/Online-Manual-blue)](https://zhjx19.github.io/mathmodels-book/)

An R package providing a versatile toolkit for mathematical modeling, developed as a companion to the book *Mathematical Modeling: Algorithms and Programming Implementation* (China Machine Press). It focuses on implementing rigorous algorithms in a user-friendly manner.

**Current Version (0.0.5) Focus**: Evaluation algorithms, including data preprocessing, subjective (AHP) and objective (Entropy, CRITIC, PCA) weighting methods, weight combination techniques, comprehensive evaluation (TOPSIS, Fuzzy, RSR, DEA), inequality measures (Gini, Theil), and grey prediction models (GM(1,1), GM(1,N), Verhulst).

## Key Features

*   **Rigorous & Comprehensive**: Implements core algorithms for evaluation and analysis with attention to detail.
*   **User-Friendly Interface**: Simplifies complex workflows with intuitive function design and intelligent handling of data preprocessing (standardization, normalization, directionality).
*   **Tidyverse Integration**: Seamlessly works with `|>` and `tidyverse` tools for smooth data manipulation and batch processing.

## Installation

You can install the latest development version of `mathmodels` directly
from [GitHub](https://github.com/zhjx19/mathmodes) use:

``` r
remotes::install_github("zhjx19/mathmodels")
```

Or download to current path, unzip and use:

```r
# install.packages("deaR")    # Install the dependency packages first
install.packages("mathmodels-main", repos=NULL, type="source")
```

## Getting Started

```r
library(mathmodels)
A = matrix(c(1,   1/2, 4, 3,   3,
             2,   1,   7, 5,   5,
             1/4, 1/7, 1, 1/2, 1/3,
             1/3, 1/5, 2, 1,   1,
             1/3, 1/5, 3, 1,   1), byrow = TRUE, nrow = 5)
# Analytic Hierarchy Process
AHP(A)
## $w
## [1] 0.26360349 0.47583538 0.05381460 0.09806829 0.10867824
## 
## $CR
## [1] 0.01609027
## 
## $Lmax
## [1] 5.072084
## 
## $CI
## [1] 0.0180211
```

## Learning More

For detailed documentation, tutorials, and in-depth examples on using the `mathmodels` package, please refer to the comprehensive online manual:

📘 **[mathmodels Package Manual - Simplifying Mathematical Modeling (Online Book)](https://zhjx19.github.io/mathmodels-book/)**

This online book is the definitive guide to the package's functionalities. Currently, only evaluation algorithms are fully implemented, including:

- Indicator data preprocessing
- AHP
- Entropy weighting
- CRITIC weighting
- PCA weighting
- Combination of subjective and objective weighting
- TOPSIS
- Grey Relational Analysis (GRA)
- Rank Sum Ratio (RSR) method
- Fuzzy Comprehensive Evaluation (FCE)
- Data Envelopment Analysis (CCR/BCC/SBM, Malmquist)
- Inequality Measures (Gini, Theil Index)
- Regional Economics (LQ/HHI/EG Index)
- Coupling coordination degree and obstacle degree







