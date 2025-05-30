% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entropy_weight.R
\name{entropy_weight}
\alias{entropy_weight}
\title{Entropy Weight Method for Indicator Weighting and Sample Scoring}
\usage{
entropy_weight(X, index = NULL)
}
\arguments{
\item{X}{A numeric data frame or matrix where rows represent samples (observations)
and columns represent indicators (variables).}

\item{index}{A character vector indicating the direction of each indicator.
Use \code{"+"} for positive indicators (higher is better) and \code{"-"} for negative
indicators (lower is better). If not provided, all indicators are assumed
to be positive.}
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
# Example: Using entropy weight method on a simple dataset
X = data.frame(
  x1 = c(3, 5, 2, 7),
  x2 = c(10, 20, 15, 25)
)
index = c("+", "-")
entropy_weight(X, index)
}
