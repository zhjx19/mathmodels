# ============================================================================
# Three-line (booktabs-style) table export: Word via flextable, LaTeX via
# hand-built booktabs code. Accepts tidy data frames (descriptive statistics,
# stat_* results) and regression results from the reg_* module.
# ============================================================================

#' Three-line Table Export (Word / LaTeX)
#'
#' Formats a tidy data frame (descriptive statistics, `stat_*` results, or
#' any table-like data frame) or a `reg_*()` regression result as a
#' three-line (booktabs) table, and writes it to Word (via `flextable`,
#' mirroring `as_flextable(x) |> font(...) |> save_as_docx()`) or to LaTeX
#' booktabs code.
#'
#' @param x A data frame, or the list returned by `reg_lm()`,
#'   `reg_logistic()`, `reg_poisson()`, or `reg_negbin()`.
#' @param output `"docx"` (Word) or `"latex"` (booktabs code).
#' @param path File path to write. Required for `output = "docx"`;
#'   for `"latex"` the code is also returned invisibly when `path` is NULL.
#' @param font Font name for the Word table.
#' @param caption Table caption.
#' @param note Table footnote (rendered below the table).
#' @param digits Digits for numeric columns.
#' @param font_size Base font size (Word only).
#'
#' @return `output = "latex"`: the LaTeX code as a character vector
#'   (invisibly when `path` is given). `output = "docx"`: the `flextable`
#'   object invisibly.
#'
#' @examples
#' d = data.frame(Var = c("mean", "sd"), Value = c(20.09, 6.03))
#' code = export_table(d, output = "latex", caption = "Summary")
#' \dontrun{
#' export_table(d, output = "docx", path = "table1.docx")
#' }
#' @export
export_table = function(x, output = c("docx", "latex"), path = NULL,
                        font = "Times New Roman", caption = NULL,
                        note = NULL, digits = 3, font_size = 10) {
  UseMethod("export_table")
}

# Extract a plain table from a reg_* result list.
.as_coef_table = function(x, digits = 3) {
  coef = if (!is.list(x)) NULL else {
    if (!is.null(x$coefficients) && is.data.frame(x$coefficients)) {
      x$coefficients
    } else if (!is.null(x$coefficient) && is.data.frame(x$coefficient)) {
      x$coefficient
    } else {
      NULL
    }
  }
  if (is.null(coef)) {
    stop("export_table.list() requires a reg_* result with a ",
         "`coefficient(s)` data frame.", call. = FALSE)
  }
  keep = c("term", "estimate", "std.error", "statistic", "p.value")
  keep = keep[keep %in% names(coef)]
  tab = coef[, keep, drop = FALSE]
  if ("p.value" %in% names(tab)) {
    tab$p.value = format_p(tab$p.value)
    tab$signif = .sig_stars(suppressWarnings(as.numeric(
      coef[["p.value"]])))
  }
  note_extra = NULL
  mi = if (!is.null(x$model_info) && is.data.frame(x$model_info)) {
    x$model_info
  } else {
    NULL
  }
  if (!is.null(mi) && nrow(mi) > 0) {
    mcols = intersect(c("metric", "value", "estimate"), names(mi))
    if (length(mcols) >= 2) {
      note_extra = paste0(mi[[mcols[1]]], " = ",
        format(round(as.numeric(mi[[mcols[2]]]), digits)),
        collapse = "; ")
    }
  }
  attr(tab, "note_extra") = note_extra
  tab
}

#' @export
export_table.list = function(x, output = "docx", path = NULL,
                             font = "Times New Roman", caption = NULL,
                             note = NULL, digits = 3, font_size = 10) {
  output = match.arg(output, c("docx", "latex"))
  tab = .as_coef_table(x, digits = digits)
  extra = attr(tab, "note_extra")
  note = if (!is.null(note) && !is.null(extra)) {
    paste0(note, "; ", extra)
  } else if (!is.null(extra)) {
    extra
  } else {
    note
  }
  export_table.data.frame(tab, output = output, path = path, font = font,
    caption = caption, note = note, digits = digits, font_size = font_size)
}

# Format p-values: "<0.001" style with fixed digits.
format_p = function(p, digits = 3) {
  out = ifelse(is.na(p), "", ifelse(p < 10^(-digits),
    paste0("<", format(10^(-digits))), format(round(p, digits))))
  out
}

#' @export
export_table.data.frame = function(x, output = "docx", path = NULL,
                                   font = "Times New Roman", caption = NULL,
                                   note = NULL, digits = 3, font_size = 10) {
  output = match.arg(output, c("docx", "latex"))
  if (!is.data.frame(x)) stop("x must be a data frame.", call. = FALSE)
  tab = x

  # Format numeric columns; p-like columns get the <0.001 treatment
  for (nm in names(tab)) {
    if (is.numeric(tab[[nm]])) {
      if (grepl("p$", nm) && !grepl("prop|term", nm)) {
        tab[[nm]] = format_p(tab[[nm]], digits)
      } else {
        tab[[nm]] = format(round(tab[[nm]], digits), nsmall = digits)
      }
    }
  }
  tab[] = lapply(tab, as.character)

  if (output == "latex") {
    code = .latex_booktabs(tab, caption = caption, note = note)
    if (!is.null(path)) writeLines(code, path)
    invisible(code)
  } else {
    ft = flextable::flextable(tab)
    ft = flextable::theme_booktabs(ft)
    ft = flextable::font(ft, fontname = font, part = "all")
    ft = flextable::fontsize(ft, size = font_size, part = "all")
    if (!is.null(caption)) ft = flextable::set_caption(ft, caption)
    if (is.null(path)) {
      stop("path is required for output = 'docx'.", call. = FALSE)
    }
    flextable::save_as_docx(ft, path = path)
    invisible(ft)
  }
}

# Build LaTeX booktabs code from an already-formatted character data frame.
.latex_booktabs = function(tab, caption = NULL, note = NULL) {
  esc = function(s) {
    s = gsub("\\\\", "\\\\textbackslash{}", s)
    s = gsub("([_&%$#{}])", "\\\\\\1", s)
    s
  }
  ncol_ = ncol(tab)
  align = paste0(c("l", rep("r", ncol_ - 1)), collapse = "")
  header = paste0(esc(names(tab)), collapse = " & ")
  body = apply(tab, 1, function(r) paste0(esc(r), collapse = " & "))
  lines = c(
    "% Requires \\usepackage{booktabs}",
    "\\begin{table}[htbp]",
    "\\centering",
    if (!is.null(caption)) paste0("\\caption{", esc(caption), "}"),
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    header,
    "\\midrule",
    body,
    "\\bottomrule",
    "\\end{tabular}",
    if (!is.null(note)) c("\\par\\footnotesize{", esc(note), "}")
      else NULL,
    "\\end{table}"
  )
  lines
}
