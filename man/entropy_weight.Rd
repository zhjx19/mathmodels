% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entropy_weight.R
\name{entropy_weight}
\alias{entropy_weight}
\title{Entropy Weight Method}
\usage{
entropy_weight(X, index = NULL, epsilon = 0.002)
}
\arguments{
\item{X}{A numeric data frame or matrix where rows represent samples (observations)
and columns represent indicators (variables).}

\item{index}{A character vector indicating the direction of each indicator.
Use \code{"+"} for positive indicators (higher is better),
\code{"-"} for negative indicators (lower is better),
and \code{NA} for already normalized indicators (no rescaling will be applied,
but minor adjustments will still be made to avoid log(0) errors).
If \code{index = NULL} (default), all indicators are treated as \code{NA},
meaning no normalization or rescaling is performed,
but a small adjustment is still applied to prevent log(0) errors.}

\item{epsilon}{A small constant used to replace exact 0s and 1s in the data to prevent log(0) errors.
Default is 0.002.}
}
\value{
A list containing:
\item{w}{Numeric vector of weights for each indicator.}

\item{s}{Numeric vector of scores for each sample (row), scaled by 100.}
}
\description{
Computes the weights of indicators and scores of samples based on the entropy method.
This method objectively determines the importance of each indicator according to
the amount of information it contains.
}
\examples{
X = data.frame(
  x1 = c(3, 5, 2, 7),
  x2 = c(10, 20, 15, 25)
)
index = c("+", "-")
entropy_weight(X, index)

}
