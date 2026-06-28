# 评价模型 — 权重计算方法
# 
# 包含: AHP, combine_weights, critic_weight, cv_weight, 
#       entropy_weight, linear_sum, pca_weight
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

  # Consistency index
  CI = (Lmax - n) / (n - 1)

  # Consistency ratio
  # Saaty's random Consistency indexes (extended to n=15)
  RI = c(0, 0, 0.58, 0.90, 1.12, 1.24, 1.32, 1.41, 1.45, 1.49,
         1.51, 1.48, 1.56, 1.57, 1.59)
  CR = ifelse(n == 2, 0, CI / RI[n])
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
#' @return A numeric vector of combined weights, normalized to sum to 1.
#' @details
#' The function supports four methods:
#' - Linear: Combines weights as alpha * w_subj + (1 - alpha) * w_obj.
#' - Multiplicative: Combines weights as w_subj * w_obj, requiring positive weights.
#' - Game: Uses the geometric mean (sqrt(w_subj * w_obj)) to balance weights.
#' - Game_linear: Uses a game-theoretic approach by solving a linear system based on the cross-product of weights.
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
               alpha_vec = solve(S, diag(S))
               alpha_vec = alpha_vec / sum(alpha_vec)
               as.vector(w_mat %*% alpha_vec)
             }
  )
  # If needed, perform normalization
  if (type %in% c("multiplicative", "game")) {
    w = w / sum(w)
  }
  w
}


# ============================================================================
# critic_weight — CRITIC 权重法
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
#' critic_weight(X, index)
#' critic_weight(X, index, method = "entropy")

critic_weight = function(X, index = NULL, method = "std",
                         cor_method = "pearson", epsilon = 0.002) {
  if(!is.data.frame(X) && !is.matrix(X))
    stop("X must be a data frame or matrix.")
  if(ncol(X) < 2)
    stop("X must have at least 2 columns.")
  if(nrow(X) < 2)
    stop("X must have at least 2 rows.")
  if(!is.null(index) && length(index) != ncol(X))
    stop("index must have length equal to ncol(X), or be NULL.")

  if(is.null(index)) index = rep(NA, ncol(X))
  pos = which(index == "+")
  neg = which(index == "-")

  # Normalize the data
  X[,pos] = lapply(X[,pos, drop = FALSE], rescale)
  X[,neg] = lapply(X[,neg, drop = FALSE], rescale, type = "-")

  # Compute sum of conflicts (1 - correlation) for each indicator
  R = apply(1 - cor(X, method = cor_method), 2, sum)

  # Compute standard deviation or entropy of each indicator
  S = switch (method,
    "std" = apply(X, 2, sd),
    "entropy" = {
      X[X == 0] = epsilon
      X[X == 1] = 1 - epsilon
      P = data.frame(lapply(X, \(x) x / sum(x)))
      e = sapply(P, function(x) - sum(x * log(x)) /log(nrow(P)))
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
# cv_weight — 变异系数权重法
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
#' The `cv_weight` function calculates weights using the CV method.
#' For each column in `data`, the CV is computed as the standard deviation
#' divided by the mean. Weights are obtained by normalizing the CVs to sum to 1.
#' This lightweight implementation uses base R and assumes all columns are numeric
#' indicators.
#'
#' @examples
#' X = data.frame(x1 = c(10, 20, 15), x2 = c(5, 10, 8))
#' cv_weight(X)
#'
#' @export
cv_weight = function(X) {
  if(!is.data.frame(X) && !is.matrix(X))
    stop("X must be a data frame or matrix.")
  if(ncol(X) < 2)
    stop("X must have at least 2 columns.")
  if(any(apply(X, 2, function(x) any(x <= 0, na.rm = TRUE))))
    stop("All values in X must be positive for CV weighting.")

  means = apply(X, 2, mean, na.rm = TRUE)
  sds = apply(X, 2, sd, na.rm = TRUE)
  cv = sds / means     # Compute coefficient of variation
  cv / sum(cv)
}


# ============================================================================
# entropy_weight — 熵权法
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
#' entropy_weight(X, index)
#'

entropy_weight = function(X, index = NULL, epsilon = 0.002) {
  if(!is.data.frame(X) && !is.matrix(X))
    stop("X must be a data frame or matrix.")
  if(ncol(X) < 2)
    stop("X must have at least 2 columns.")
  if(nrow(X) < 2)
    stop("X must have at least 2 rows.")
  if(!is.null(index) && length(index) != ncol(X))
    stop("index must have length equal to ncol(X), or be NULL.")

  if(is.null(index)) index = rep(NA, ncol(X))
  if(any(is.na(index))) {
    X_na = X[is.na(index)]
    X_na[X_na == 0] = epsilon
    X_na[X_na == 1] = 1- epsilon
    X[is.na(index)] = X_na
  }
  pos = which(index == "+")
  neg = which(index == "-")
  # Normalize the data
  X[,pos] = lapply(X[,pos, drop = FALSE],
                    function(x) rescale(x, a = epsilon, b = 1-epsilon))
  X[,neg] = lapply(X[,neg, drop = FALSE],
                    function(x) rescale(x, type = "-", a = epsilon, b = 1-epsilon))
  # Compute proportion p(i,j) of sample i in indicator j
  P = data.frame(lapply(X, function(x) x / sum(x)))
  # Compute entropy e(j) for each indicator j
  e = sapply(P, \(x) -sum(x * log(x)) / log(nrow(P)))
  d = 1 - e          # Compute redundancy degree
  w = d / sum(d)     # Compute weight vector
  # Compute sample scores
  s = as.vector(100 * as.matrix(X) %*% w)
  list(w = w, s = s)
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
  if(!is.data.frame(data) && !is.matrix(data))
    stop("data must be a data frame or matrix.")
  if(!is.numeric(w) || is.matrix(w))
    stop("w must be a numeric vector.")
  if(length(w) != ncol(data))
    stop("Length of w must equal ncol(data).")
  as.matrix(data) %*% w |> as.vector()
}


# ============================================================================
# pca_weight — PCA 权重法
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
#' pca_weight(iris[1:10, 1:4], ind, nfs = 2)

pca_weight = function(X, index = NULL, nfs = NULL, varimax = TRUE,
                      method = "abs") {
  if(!is.data.frame(X) && !is.matrix(X))
    stop("X must be a data frame or matrix.")
  if(ncol(X) < 2)
    stop("X must have at least 2 columns.")
  if(nrow(X) < 2)
    stop("X must have at least 2 rows.")
  if(!is.null(index) && length(index) != ncol(X))
    stop("index must have length equal to ncol(X), or be NULL.")

  if(is.null(index)) index = rep(NA, ncol(X))
  pos = which(index == "+")
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
