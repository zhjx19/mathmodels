#' Read and Combine National Bureau of Statistics XLS Files
#'
#' This function reads multiple XLS files downloaded from the National Bureau of Statistics of China
#' website (\url{http://www.stats.gov.cn}) without requiring renaming. It extracts the variable name
#' from the "indicator: XXX" text in cell A1 of each file, skips the first 3 rows of headers, reads
#' the first 31 rows of data, and combines the data into a single data frame by stacking vertically
#' based on the indicator names. It also simplifies region names by removing suffixes such as "city",
#' "province", "autonomous region", "clan", or "Uyghur".
#'
#' @param paths A character vector containing the file paths to the XLS files.
#' @return A tibble containing the combined data with an additional column `indicator` representing
#'   the original indicator names and a `region` column with simplified region names.
#' @importFrom readxl read_excel
#' @importFrom dplyr pull mutate rename
#' @importFrom purrr map map_chr set_names list_rbind
#' @importFrom stringr str_remove str_remove_all
#'
#' @examples
#' \dontrun{
#' paths = c("file1.xls", "file2.xls")
#' data = read_nbs(paths)
#' }

#' @export
read_nbs = function(paths) {
  vars = purrr::map_chr(paths, \(path) {
    readxl::read_excel(path, range = "A1:A2") |>
      dplyr::pull(1) |> stringr::str_remove("\\u6307\\u6807\\uff1a")
  })
  purrr::map(paths, \(path) readxl::read_excel(path, skip = 3, n_max = 31)) |>
    rlang::set_names(vars) |>
    purrr::list_rbind(names_to = "indicator") |>
    dplyr::rename(region = "\u5730\u533a") |>
    dplyr::mutate(region = stringr::str_remove_all(region, "\\u5e02|\\u7701|\\u81ea\\u6cbb\\u533a|.\\u65cf|\\u7ef4\\u543e\\u5c14"))
}
