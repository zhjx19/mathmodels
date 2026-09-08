# Three-line Table Export (Word / LaTeX)

Formats a tidy data frame (descriptive statistics, `stat_*` results, or
any table-like data frame) or a `reg_*()` regression result as a
three-line (booktabs) table, and writes it to Word (via `flextable`,
mirroring `as_flextable(x) |> font(...) |> save_as_docx()`) or to LaTeX
booktabs code.

## Usage

``` r
export_table(
  x,
  output = c("docx", "latex"),
  path = NULL,
  font = "Times New Roman",
  caption = NULL,
  note = NULL,
  digits = 3,
  font_size = 10
)
```

## Arguments

- x:

  A data frame, or the list returned by
  [`reg_lm()`](https://zhjx19.github.io/mathmodels/reference/reg_lm.md),
  [`reg_logistic()`](https://zhjx19.github.io/mathmodels/reference/reg_logistic.md),
  [`reg_poisson()`](https://zhjx19.github.io/mathmodels/reference/reg_poisson.md),
  or
  [`reg_negbin()`](https://zhjx19.github.io/mathmodels/reference/reg_negbin.md).

- output:

  `"docx"` (Word) or `"latex"` (booktabs code).

- path:

  File path to write. Required for `output = "docx"`; for `"latex"` the
  code is also returned invisibly when `path` is NULL.

- font:

  Font name for the Word table.

- caption:

  Table caption.

- note:

  Table footnote (rendered below the table).

- digits:

  Digits for numeric columns.

- font_size:

  Base font size (Word only).

## Value

`output = "latex"`: the LaTeX code as a character vector (invisibly when
`path` is given). `output = "docx"`: the `flextable` object invisibly.

## Examples

``` r
d = data.frame(Var = c("mean", "sd"), Value = c(20.09, 6.03))
code = export_table(d, output = "latex", caption = "Summary")
if (FALSE) { # \dontrun{
export_table(d, output = "docx", path = "table1.docx")
} # }
```
