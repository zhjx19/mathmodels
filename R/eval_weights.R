# Evaluation Algorithms — Weighting Methods
#
# Includes: AHP, combine_weights, weight_critic, weight_cv,
# weight_entropy, linear_sum, and weight_pca
#


# ============================================================================
# AHP — 层次分析法
# ============================================================================

#' @title AHP: Analytic Hierarchy Process
#'
#' @description AHP is a multi-criteria decision
#' analysis method developed by Saaty, which can also be used to
#' determine indicator weights.
#'
#' @param A a numeric matrix, i.e. pairwise comparison matrix
#'
#' @return a list object that contains: w (Weight vector), CR (Consistency ratio),
#' Lmax (Maximum eigenvalue), CI (Consistency index)
#' @export
#'
#' @examples
#' A = matrix(c(1,   1/2, 4, 3,   3,
#'              2,   1,   7, 5,   5,
#'              1/4, 1/7, 1, 1/2, 1/3,
#'              1/3, 1/5, 2, 1,   1,
#'              1/3, 1/5, 3, 1,   1), byrow = TRUE, nrow = 5)
#' AHP(A)


AHP = function(A) {

  if(!is.matrix(A)) stop("A must be a matrix.")
  n = nrow(A)
  if(n != ncol(A)) stop("A must be a square matrix.")
  if(n > 15) stop("A must have at most 15 rows/columns.")
  if(any(A <= 0)) stop("A must contain only positive values.")

  for(i in 1:n) {
    for(j in i:n) {
      if(abs(A[i,j] * A[j,i] - 1) > 1.5e-8)
        warning(paste0("A[", i, ",", j, "] violate reciprocal rule!"))
    }
  }

  rlt = eigen(A)
  Lmax = Re(rlt$values[1])   # Maximum eigenvalue

  # Weight vector
  w = Re(rlt$vectors[,1]) / sum(Re(rlt$vectors[,1]))

  # Consistency ratio
  # Saaty's random Consistency indexes (extended to n=15)
  RI = c(0, 0, 0.58, 0.90, 1.12, 1.24, 1.32, 1.41, 1.45, 1.49,
         1.51, 1.54, 1.56, 1.57, 1.59)
  # For n <= 2 the comparison matrix cannot be inconsistent
  if (n <= 2) {
    CI = 0
    CR = 0
  } else {
    # Consistency index
    CI = (Lmax - n) / (n - 1)
    # Consistency ratio
    CR = CI / RI[n]
  }
  list(w = w, CR = CR, Lmax = Lmax, CI = CI)
}


# ============================================================================
# combine_weights — 组合主客观权重
# ============================================================================

#' @title Combine Subjective and Objective Weights
#' @description Combines subjective and objective weights using linear, multiplicative, or game theory-based methods (geometric mean or linear system).
#' @param w_subj Numeric vector of subjective weights.
#' @param w_obj Numeric vector of objective weights.
#' @param type Character string specifying the combination method: "linear", "multiplicative", "game", or "game_linear".
#' @param alpha Numeric value between 0 and 1, used only for the linear method to weight subjective weights. Defaults to 0.5.
#' @return A numeric vector of combined weights. For `multiplicative` and `game`,
#' the result is normalized to sum to 1. For `linear` and `game_linear`, the
#' theoretically-defined combination is returned as is: it sums to 1 when both
#' input weight vectors do (only the two normalized methods are rescaled).
#' @details
#' The function supports four methods:
#' - Linear: Combines weights as alpha * w_subj + (1 - alpha) * w_obj.
#' - Multiplicative: Combines weights as w_subj * w_obj, requiring positive weights, then normalizes to sum to 1.
#' - Game: Uses the geometric mean (sqrt(w_subj * w_obj)) to balance weights, then normalizes to sum to 1.
#' - Game_linear: Uses a game-theoretic approach by solving a linear system based on the cross-product of weights.
#' Inputs are expected to be weight vectors that each sum to 1.

#' @examples
#' w_subj = c(0.4, 0.3, 0.2, 0.1)
#' w_obj = c(0.25, 0.2, 0.3, 0.25)
#' combine_weights(w_subj, w_obj, type = "linear", alpha = 0.6)
#' combine_weights(w_subj, w_obj, type = "multiplicative")
#' combine_weights(w_subj, w_obj, type = "game")
#' combine_weights(w_subj, w_obj, type = "game_linear")
#'
#' @export
combine_weights = function(w_subj, w_obj, type = "linear", alpha = 0.5) {
  if(!is.numeric(w_subj) || is.matrix(w_subj))
    stop("w_subj must be a numeric vector.")
  if(!is.numeric(w_obj) || is.matrix(w_obj))
    stop("w_obj must be a numeric vector.")
  if(length(w_subj) != length(w_obj))
    stop("w_subj and w_obj must have the same length.")
  if (any(w_subj < 0) || any(w_obj < 0)) stop("Weights must be non-negative.")

  # Combine weights based on type
  w = switch(type,
             linear = {
               if (alpha < 0 || alpha > 1) stop("Alpha must be between 0 and 1.")
               alpha * w_subj + (1 - alpha) * w_obj
             },
             multiplicative = {
               if (any(w_subj <= 0) || any(w_obj <= 0)) stop("Weights must be positive for multiplicative method.")
               w_subj * w_obj
             },
             game = sqrt(w_subj * w_obj),
             game_linear = {
               # Solve linear system for game theory-based weights
               w_mat = cbind(w_subj, w_obj)
               S = crossprod(w_mat)
               alpha_vec = tryCatch(
                 solve(S, diag(S)),
                 error = \(e) MASS::ginv(S) %*% diag(S))
               alpha_vec = alpha_vec / sum(alpha_vec)
               as.vector(w_mat %*% alpha_vec)
             },
             stop('type must be one of "linear", "multiplicative", "game", "game_linear".')
  )
  # If needed, perform normalization
  if (type %in% c("multiplicative", "game")) {
    w = w / sum(w)
  }
  w
}


# ============================================================================
# weight_critic — CRITIC 权重法
# ============================================================================

#' @title CRITIC Weight Method
#'
#' @description Computes objective weights of indicators and scores of samples using the CRITIC method.
#' The method considers both the contrast intensity (e.g., standard deviation or entropy)
#' and conflict among indicators (based on correlation) to determine indicator importance.
#' This version supports different methods for contrast intensity and correlation types.
#'
#' @param X A numeric data frame or matrix where rows represent samples (observations)
#'          and columns represent indicators (variables).
#' @param index A character vector indicating the direction of each indicator:
#'              Use `"+"` for positive indicators (higher is better),
#'              `"-"` for negative indicators (lower is better),
#'              and `NA` for already normalized indicators (no rescaling will be applied).
#'
#'              If `index = NULL` (default), all indicators are treated as `NA`, meaning no normalization is performed.
#' @param method Character scalar; specifies the method used to compute contrast intensity.
#'               Options: `"std"` (standard deviation, default), or `"entropy"` (based on information redundancy).
#' @param cor_method Character scalar; specifies the method for computing correlations.
#'                   Options: `"pearson"` (default), `"spearman"`, or `"kendall"`.
#' @param epsilon A small constant used to replace exact 0s and 1s in the data to prevent log(0) errors.
#'                Default is 0.002. Only used when `method = "entropy"`.
#'
#' @return A list containing:
#' \item{w}{Numeric vector of weights for each indicator.}
#' \item{s}{Numeric vector of scores for each sample (row), scaled by 100.}
#'
#' @importFrom stats cor
#' @export
#' @examples
#' # Example: Using CRITIC method on a simple dataset
#' X = data.frame(
#'   x1 = c(3, 5, 2, 7),
#'   x2 = c(10, 20, 15, 25)
#' )
#' index = c("+", "-")
#' weight_critic(X, index)
#' weight_critic(X, index, method = "entropy")

weight_critic = function(X, index = NULL, method = "std",
                         cor_method = "pearson", epsilon = 0.002) {
  X = .check_matrix_like(X, "X must be a data frame or matrix.", coerce = TRUE)
  if (anyNA(X)) stop("X must not contain NA.")
  if(ncol(X) < 2)
    stop("X must have at least 2 columns.")
  if(nrow(X) < 3)
    stop("X must have at least 3 rows (correlations are degenerate with 2 rows).")
  index = .check_index(index, ncol(X), "index must have length equal to ncol(X), or be NULL.",
                       'index must contain only "+", "-", or NA', allow_na = TRUE)

  # Normalize the data
  X = .rescale_indicators(X, index)

  # Compute sum of conflicts (1 - correlation) for each indicator
  R = apply(1 - cor(X, method = cor_method), 2, sum)

  # Compute standard deviation or entropy of each indicator
  S = switch (method,
    "std" = apply(X, 2, stats::sd),
    "entropy" = {
      if (any(X < 0))
        stop("weight_critic() method = 'entropy' requires non-negative data.")
      X[X == 0] = epsilon

      X[X == 1] = 1 - epsilon
      P = data.frame(lapply(X, \(x) x / sum(x)))
      e = sapply(P, \(x) - sum(x * log(x)) /log(nrow(P)))
      1 - e
    })

  # Compute indicator importance (C)
  C = S * R
  # Compute weights
  w = C / sum(C)
  # Compute sample scores
  s = as.vector(100 * as.matrix(X) %*% w)

  list(w = w, s = s)
}


# ============================================================================
# weight_cv — 变异系数权重法
# ============================================================================

#' Coefficient of Variation Weighting
#'
#' @description
#' Computes weights for indicators using the Coefficient of Variation (CV) method.
#' Weights are derived by normalizing the CV (standard deviation divided by mean)
#' for each indicator.
#'
#' @param X Numeric matrix or data frame with positive indicator data.
#'
#' @return Numeric vector of weights for the indicators, summing to 1.
#'
#' @details
#' The `weight_cv` function calculates weights using the CV method.
#' For each column in `data`, the CV is computed as the standard deviation
#' divided by the mean. Weights are obtained by normalizing the CVs to sum to 1.
#' This lightweight implementation uses base R and assumes all columns are numeric
#' indicators.
#'
#' @examples
#' X = data.frame(x1 = c(10, 20, 15), x2 = c(5, 10, 8))
#' weight_cv(X)
#'
#' @export
weight_cv = function(X) {
  .check_matrix_like(X, "X must be a data frame or matrix.")
  if(ncol(X) < 2)
    stop("X must have at least 2 columns.")
  if (anyNA(X)) stop("X must not contain NA.")
  if(any(apply(X, 2, \(x) any(x <= 0, na.rm = TRUE))))
    stop("All values in X must be positive for CV weighting.")

  means = apply(X, 2, mean, na.rm = TRUE)
  sds = apply(X, 2, stats::sd, na.rm = TRUE)
  cv = sds / means     # Compute coefficient of variation
  cv / sum(cv)
}


# ============================================================================
# weight_entropy — 熵权法
# ============================================================================

#' @title Entropy Weight Method
#'
#' @description Computes the weights of indicators and scores of samples based on the entropy method.
#' This method objectively determines the importance of each indicator according to
#' the amount of information it contains.
#'
#' @param X A numeric data frame or matrix where rows represent samples (observations)
#'          and columns represent indicators (variables).
#'
#' @param index A character vector indicating the direction of each indicator.
#'              Use `"+"` for positive indicators (higher is better),
#'              `"-"` for negative indicators (lower is better),
#'              and `NA` for already normalized indicators (no rescaling will be applied,
#'              but minor adjustments will still be made to avoid log(0) errors).
#'              If `index = NULL` (default), all indicators are treated as `NA`,
#'              meaning no normalization or rescaling is performed,
#'              but a small adjustment is still applied to prevent log(0) errors.
#' @param epsilon A small constant used to replace exact 0s and 1s in the data to prevent log(0) errors.
#'                Default is 0.002.
#'
#' @return A list containing:
#' \item{w}{Numeric vector of weights for each indicator.}
#'
#' \item{s}{Numeric vector of scores for each sample (row), scaled by 100.}
#'
#' @export
#'
#' @examples
#' X = data.frame(
#'   x1 = c(3, 5, 2, 7),
#'   x2 = c(10, 20, 15, 25)
#' )
#' index = c("+", "-")
#' weight_entropy(X, index)
#'

weight_entropy = function(X, index = NULL, epsilon = 0.002) {
  X = .check_matrix_like(X, "X must be a data frame or matrix.", coerce = TRUE)
  if (anyNA(X)) stop("X must not contain NA.")
  if(ncol(X) < 2)
    stop("X must have at least 2 columns.")
  if(nrow(X) < 2)
    stop("X must have at least 2 rows.")
  index = .check_index(index, ncol(X), "index must have length equal to ncol(X), or be NULL.",
                       'index must contain only "+", "-", or NA', allow_na = TRUE)
  if(any(is.na(index))) {
    X_na = X[is.na(index)]
    X_na[X_na == 0] = epsilon
    X_na[X_na == 1] = 1- epsilon
    X[is.na(index)] = X_na
  }
  pos = which(index == "+")
  neg = which(index == "-")
  # Normalize the data (the epsilon fix above runs first; NA columns are
  # untouched by either step)
  X = .rescale_indicators(X, index, a = epsilon, b = 1 - epsilon)
  # Compute proportion p(i,j) of sample i in indicator j
  if (any(X < 0)) stop("weight_entropy() requires non-negative data.")
  # Mixed scales: + / - columns were rescaled to [0, 1] just above, so any
  # index = NA column is expected to be on that scale too. Only NA columns can
  # exceed 1 here; when one does, `s` combines incommensurable units.
  # (index = NULL rescales nothing, so it is left alone.)
  na_cols = is.na(index)
  if ((length(pos) + length(neg) > 0) && any(na_cols) &&
      any(X[, na_cols, drop = FALSE] > 1))
    warning("weight_entropy(): columns with index = NA are expected to be ",
            "already normalized to [0, 1]; found values > 1 while other ",
            "columns were rescaled, so the returned 100-scale scores mix ",
            "units and may be meaningless.", call. = FALSE)
  P = data.frame(lapply(X, \(x) x / sum(x)))
  # Compute entropy e(j) for each indicator j
  e = sapply(P, \(x) -sum(x * log(x)) / log(nrow(P)))
  d = 1 - e          # Compute redundancy degree
  w = d / sum(d)     # Compute weight vector
  # Compute sample scores
  s = as.vector(100 * as.matrix(X) %*% w)
  list(w = w, s = s)
}


# ============================================================================
# weight_fa — 因子分析权重法
# ============================================================================

#' @title Factor-Analysis-Based Weighting Method
#'
#' @description Computes indicator weights via exploratory factor analysis
#' (`stats::factanal`). The (optionally rotated) loading matrix is combined
#' with the variance share of each factor to derive objective weights:
#' \deqn{w_j \propto \sum_i \mathrm{share}_i \cdot |a_{ji}|}
#' where \eqn{a_{ji}} is the loading of indicator \eqn{j} on factor \eqn{i}
#' and \code{share_i} is the proportion of total variance explained by
#' factor \eqn{i}. Optionally handles positive/negative indicator directions
#' and supports pre-standardized data.
#'
#' @param X A numeric data frame or matrix where rows represent samples and
#'          columns represent indicators.
#' @param index A character vector indicating the direction of each indicator.
#'              Use `"+"` for positive indicators (higher is better),
#'              `"-"` for negative indicators (lower is better),
#'              and `NA` for already standardized indicators (no
#'              standardization will be applied).
#'              If `index = NULL` (default), all indicators are treated as
#'              `NA`, meaning no standardization is performed.
#' @param nfs Number of common factors to extract; must be at least 1 and
#'        smaller than the number of indicators (default 1).
#' @param rotation Rotation passed to `stats::factanal`: `"varimax"`
#'        (default), `"promax"`, `"none"`, or the name of any available
#'        rotation function.
#' @param method Weighting method: `"abs"` uses absolute loadings
#'        \eqn{|a_{ji}|} (default), `"squared"` uses \eqn{a_{ji}^2}.
#'
#' @return A list containing:
#' \item{w}{Numeric vector of normalized weights for each indicator.}
#' \item{s}{Numeric vector of scores for each sample.}
#' \item{loading}{The (rotated) loading matrix.}
#' \item{share}{Variance share of each retained factor.}
#' \item{uniqueness}{Uniquenesses of the indicators.}
#'
#' @examples
#' set.seed(1)
#' X = data.frame(x1 = rnorm(50), x2 = rnorm(50), x3 = rnorm(50))
#' weight_fa(X, nfs = 1)
#'
#' @export
weight_fa = function(X, index = NULL, nfs = 1, rotation = "varimax",
                     method = "abs") {
  X = .check_matrix_like(X, "X must be a data frame or matrix.", coerce = TRUE)
  if (anyNA(X))
    stop("X must not contain NA.")
  if (ncol(X) < 2)
    stop("X must have at least 2 columns.")
  if (nrow(X) < 2)
    stop("X must have at least 2 rows.")
  if (!all(vapply(X, is.numeric, logical(1))))
    stop("X must be all numeric.")
  index = .check_index(index, ncol(X), "index must have length equal to ncol(X), or be NULL.",
                       'index must contain only "+", "-", or NA', allow_na = TRUE)
  if (!is.numeric(nfs) || length(nfs) != 1 || nfs < 1 ||
      nfs != floor(nfs)) {
    stop("nfs must be a positive integer.", call. = FALSE)
  }
  if (nfs >= ncol(X)) {
    stop("nfs must be smaller than the number of indicators.",
         call. = FALSE)
  }
  if (nrow(X) <= ncol(X)) {
    stop("X must have more rows than columns for factor analysis.",
         call. = FALSE)
  }
  if (!method %in% c("abs", "squared")) {
    stop('method must be "abs" or "squared".', call. = FALSE)
  }
  if (!is.character(rotation) || length(rotation) != 1) {
    stop("rotation must be a single character value.", call. = FALSE)
  }
  if (!rotation %in% c("varimax", "promax", "none") &&
      !exists(rotation, mode = "function", envir = .GlobalEnv)) {
    stop("rotation must be 'varimax', 'promax', 'none', or the name of ",
         "an available rotation function.", call. = FALSE)
  }
  if (any(vapply(X, \(x) stats::sd(x) == 0, logical(1)))) {
    stop("X must not contain constant columns.", call. = FALSE)
  }

  # Standardize (directional) columns. X is a data frame here (guaranteed by the
  # guard above), so column subsetting is safe for a matrix argument too.
  idx_std = which(!is.na(index))
  X[idx_std] = lapply(X[idx_std], \(x) as.vector(scale(x)))
  neg = which(index == "-")
  if (length(neg) > 0) {
    X[, neg] = lapply(X[, neg, drop = FALSE], `*`, -1)
  }

  fit = stats::factanal(X, factors = nfs, rotation = rotation,
    scores = "none")
  A = unclass(fit$loadings)
  # Variance share of each retained factor (standardized total variance = p)
  share = colSums(A^2) / ncol(X)

  beta = switch(method,
    "abs" = abs(A) %*% share,
    "squared" = (A^2) %*% share)
  w = as.vector(beta / sum(beta))

  s = as.vector(as.matrix(X) %*% w)

  list(w = w, s = s, loading = A, share = as.vector(share),
       uniqueness = as.numeric(fit$uniquenesses))
}


# ============================================================================
# linear_sum — 线性加权综合
# ============================================================================

#' Linear weighted synthesis
#'
#' @param data Indicator matrix or data frame (rows = samples, columns = indicators)
#' @param w Weight vector (length must match ncol(data))
#'
#' @return Vector of comprehensive scores
#'
#' @export
#'
#' @examples
#' data = data.frame(
#'   GDP = c(0.85, 0.72, 0.91, NA),
#'   Employment = c(0.78, 0.85, 0.67, 0.73),
#'   Environment = c(0.65, 0.72, NA, 0.81)
#' )
#' w = c(0.5, 0.3, 0.2)
#' linear_sum(data, w)

linear_sum = function(data, w) {
  .check_matrix_like(data, "data must be a data frame or matrix.")
  if(!is.numeric(w) || is.matrix(w))
    stop("w must be a numeric vector.")
  if(length(w) != ncol(data))
    stop("Length of w must equal ncol(data).")
  as.matrix(data) %*% w |> as.vector()
}


# ============================================================================
# weight_pca — PCA 权重法
# ============================================================================

#' @title PCA-Based Weighting Method
#'
#' @description Computes indicator weights using Principal Component Analysis (PCA).
#' The method extracts principal components and uses their variance contribution
#' to derive objective weights for indicators. Optionally handles positive/negative
#' directions of indicators, and supports pre-standardized data.
#'
#' @param X A numeric data frame or matrix where rows represent samples and
#'          columns represent indicators.
#'
#' @param index A character vector indicating the direction of each indicator.
#'              Use `"+"` for positive indicators (higher is better),
#'              `"-"` for negative indicators (lower is better),
#'              and `NA` for already standardized indicators (no standardization will be applied).
#'
#'              If `index = NULL` (default), all indicators are treated as `NA`,
#'              meaning no standardization is performed.
#'
#' @param nfs Number of principal components to use; by default, all are used.
#'
#' @param varimax Whether to perform Varimax rotation, default is TRUE.
#'
#' @param method Weighting Method. "abs" uses absolute loading values \code{|a_{ji}|} (default), "squared" uses \code{a_{ji}^2}.
#'
#' @return A list containing:
#' \item{w}{Numeric vector of normalized weights for each indicator.}
#'
#' \item{s}{Numeric vector of scores for each sample.}
#'
#' \item{lambda}{Eigenvalues of principal components (explained variance).}
#'
#' \item{varP}{Proportion of variance explained by selected PCs.}
#'
#' @importFrom stats prcomp
#' @export
#' @examples
#' # Example: Using PCA to compute indicator weights
#' ind = c("+","+","-","-")
#' weight_pca(iris[1:10, 1:4], ind, nfs = 2)

weight_pca = function(X, index = NULL, nfs = NULL, varimax = TRUE,
                      method = "abs") {
  X = .check_matrix_like(X, "X must be a data frame or matrix.", coerce = TRUE)
  if(ncol(X) < 2)
    stop("X must have at least 2 columns.")
  if(nrow(X) < 2)
    stop("X must have at least 2 rows.")
  index = .check_index(index, ncol(X), "index must have length equal to ncol(X), or be NULL.",
                       'index must contain only "+", "-", or NA', allow_na = TRUE)

  neg = which(index == "-")

  # Standardize data
  idx = !is.na(index)
  X[idx] = scale(X[idx])
  X[,neg] = -X[,neg]

  # Perform PCA with all components and varimax rotation
  if(is.null(nfs)) nfs = ncol(X)

  pc = prcomp(X, center = FALSE, scale. = FALSE) # Already scaled above
  loadings = pc$rotation[, 1:nfs, drop = FALSE]  # Extract first nfs loadings

  # Varimax rotation
  if(varimax && nfs > 1) {
    rotated = varimax(loadings)
    A = unclass(rotated$loadings)
  } else {
    A = loadings
  }

  # Eigenvalues = variance explained by each PC
  eigenvalues = pc$sdev ^ 2
  lambda = eigenvalues[1:nfs]

  # Variance accounted by each component (proportion for selected components)
  varP = lambda / sum(eigenvalues)

  # Compute weighted contribution of each loading
  beta = switch(method,
                "abs" = abs(A) %*% varP,
                "squared" = (A^2) %*% varP)
  # Normalize weights
  w = beta / sum(beta)

  # Compute sample scores
  s = as.vector(as.matrix(X) %*% w)

  list(w = as.vector(w), s = s, cum_contrib = sum(varP), loading = A,
       lambda = lambda, varP = varP)
}
