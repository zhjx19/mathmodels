#' Water Quality Dataset
#'
#' @title Water Quality Dataset
#'
#' @description
#' A dataset containing water quality evaluation metrics for 20 rivers, including
#' dissolved oxygen (O2, positive indicator), pH value (PH, centered indicator),
#' total bacteria count (germ, negative indicator), and plant nutrient content
#' (nutrient, interval indicator with optimal range 10-20). This dataset is suitable
#' for multi-criteria decision analysis, such as weight calculation and fuzzy
#' comprehensive evaluation in the \code{mathmodels} package.
#'
#' @format
#' A data frame with 20 rows and 5 columns:
#' \describe{
#'   \item{ID}{Numeric, unique identifier for each river (1 to 20).}
#'   \item{O2}{Numeric, dissolved oxygen content (mg/L), higher values are better (positive indicator).}
#'   \item{PH}{Numeric, pH value, values closer to 7 are optimal (centered indicator).}
#'   \item{germ}{Numeric, total bacteria count, lower values are better (negative indicator).}
#'   \item{nutrient}{Numeric, plant nutrient content (mg/L), optimal range is 10-20 (interval indicator).}
#' }
#'
#' @source
#' Simulated data for water quality evaluation, created for demonstration purposes.
#'
#' @examples
#' # Load the dataset
#' data(water_quality)
#'
#' # Preview the data
#' head(water_quality)
#'
#' @name water_quality
"water_quality"
