% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topsis.R
\name{topsis}
\alias{topsis}
\title{TOPSIS Method for Multi-Criteria Decision Making}
\usage{
topsis(X, w = NULL, index = NULL)
}
\arguments{
\item{X}{A numeric matrix or data frame where rows represent alternatives and columns represent criteria.}

\item{w}{A numeric vector of weights for each criterion. Must be non-negative and sum to 1.
If not provided, equal weights are used.}

\item{index}{A character vector indicating the direction of each indicator:
Use \code{"+"} for positive indicators (higher is better),
\code{"-"} for negative indicators (lower is better).
If \code{index = NULL} (default), all indicators are treated as \code{"+"}.}
}
\value{
A named numeric vector of relative closeness scores (in \link{0, 1}) for each alternative.
Higher values indicate better alternatives.
Names are taken from \code{rownames(X)} or default to "Sample1", "Sample2", etc.
}
\description{
Implements the Technique for Order of Preference by Similarity to Ideal Solution (TOPSIS)
to rank alternatives based on multiple criteria. The function normalizes the decision matrix using Min-Max method,
applies weights, and computes relative closeness to the ideal solution.
}
\details{
The TOPSIS method ranks alternatives by:
\enumerate{
\item Normalizing the decision matrix using Min-Max normalization.
\item Applying weights to form a weighted normalized matrix.
\item Identifying positive and negative ideal solutions based on indicator directions.
\item Computing Euclidean distances to ideal solutions.
\item Calculating relative closeness as \code{S0 / (S0 + Sstar)}, where \code{S0}
is the distance to the negative ideal and \code{Sstar} is the distance to the positive ideal.
}
This implementation supports both positive and negative indicators via the \code{index} parameter.
}
\examples{
A = data.frame(
  X1 = c(2, 5, 3),  # "+"
  X2 = c(8, 1, 6)   # "-"
)
w = c(0.6, 0.4)
idx = c("+","-")
topsis(A, w, idx)
}
