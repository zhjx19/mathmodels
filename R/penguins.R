#' @title Size measurements for adult foraging penguins near Palmer Station, Antarctica
#'
#' @description Includes measurements for penguin species, island in Palmer Archipelago,
#' size (flipper length, body mass, bill dimensions), and sex.
#'
#' @docType data
#'
#' @usage data(penguins)
#'
#' @format A tibble with 344 rows and 8 variables:
#' \describe{
#'   \item{species}{a factor denoting penguin species (Adélie, Chinstrap and Gentoo)}
#'   \item{island}{a factor denoting island in Palmer Archipelago, Antarctica (Biscoe, Dream or Torgersen)}
#'   \item{bill_length_mm}{a number denoting bill length (millimeters)}
#'   \item{bill_depth_mm}{a number denoting bill depth (millimeters)}
#'   \item{flipper_length_mm}{an integer denoting flipper length (millimeters)}
#'   \item{body_mass_g}{an integer denoting body mass (grams)}
#'   \item{sex}{a factor denoting penguin sex (female, male)}
#'   \item{year}{an integer denoting the study year (2007, 2008, or 2009)}
#' }
#'
#' @references This dataset referenced from the palmerpenguins package.
#' @keywords datasets
#' @examples
#'
#' data(penguins)
#' head(penguins)
#'
"penguins"
