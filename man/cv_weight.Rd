% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv_weight.R
\name{cv_weight}
\alias{cv_weight}
\title{Coefficient of Variation Weighting}
\usage{
cv_weight(X)
}
\arguments{
\item{data}{Numeric matrix or data frame with positive indicator data.}
}
\value{
Numeric vector of weights for the indicators, summing to 1.
}
\description{
Computes weights for indicators using the Coefficient of Variation (CV) method.
Weights are derived by normalizing the CV (standard deviation divided by mean)
for each indicator.
}
\details{
The \code{cv_weight} function calculates weights using the CV method.
For each column in \code{data}, the CV is computed as the standard deviation
divided by the mean. Weights are obtained by normalizing the CVs to sum to 1.
This lightweight implementation uses base R and assumes all columns are numeric
indicators.
}
\examples{
X = data.frame(x1 = c(10, 20, 15), x2 = c(5, 10, 8))
cv_weight(X)

}
