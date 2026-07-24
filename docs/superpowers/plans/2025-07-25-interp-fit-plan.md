# interp_fit.R 模块——实现计划

> **针对代理工作器：** 推荐使用 superpowers:subagent-driven-development，或使用 superpowers:executing-plans 来逐步实现此计划。步骤使用复选框 (`- [ ]`) 语法进行跟踪。

**目标：** 构建 `R/interp_fit.R` 模块，包含 4 个插值函数和 3 个曲线拟合函数，所有拟合函数返回统一列表结构。

**架构：** 采用扁平分治方案 A——分散函数共享 5 个内部辅助函数。线性化曲线拟合通过变量变换 + `lm()` 实现；增长曲线拟合通过 `minpack.lm::nlsLM()` 实现，并自动生成初值。

**技术栈：** R, stats, tibble, minpack.lm (新增依赖)

## 全局约束

- 使用基础管道 `|>` 而非 `%>%`
- 单行匿名函数使用 `\() ...`；多行匿名函数使用 `function() { ... }`
- 所有导出函数必须包含 roxygen2 文档
- 测试文件：`tests/testthat/test-interp_fit.R`
- 完成后运行 `devtools::document()` 以重新生成 NAMESPACE 和 man/
- 添加到 `NEWS.md`

---

### 任务 1：添加 minpack.lm 依赖项

**文件：**
- 修改：`DESCRIPTION:16`（在 `dplyr` 和 `forecast` 之间添加 `minpack.lm`）

**依赖项：**
- 无

**产出：**
- `DESCRIPTION` 中现在包含 `minpack.lm`

- [ ] **第 1 步：编辑 DESCRIPTION**

在 `DESCRIPTION` 文件的 Imports 字段中添加 `minpack.lm`。插入位置在 `lpSolveAPI,` 之后。

编辑前：
```
Imports: 
    lpSolveAPI,
    dplyr (>= 1.0.0),
```

编辑后：
```
Imports: 
    lpSolveAPI,
    minpack.lm,
    dplyr (>= 1.0.0),
```

- [ ] **第 2 步：验证描述**

运行：`Rscript -e "read.dcf('DESCRIPTION', fields = 'Imports')"`
预期结果：输出中包含 `minpack.lm`。

- [ ] **第 3 步：提交**

```bash
git add DESCRIPTION
git commit -m "build: add minpack.lm dependency for interp_fit"
```

---

### 任务 2：创建测试文件骨架

**文件：**
- 创建：`tests/testthat/test-interp_fit.R`

**依赖项：**
- 无（独立）

**产出：**
- 可运行的测试文件，只包含用于占位的 "true" 测试

- [ ] **第 1 步：编写测试骨架**

```r
# =============================================================================
# Unit tests for the interpolation & curve fitting module interp_fit.R
# =============================================================================
```

将此内容写入 `tests/testthat/test-interp_fit.R`。

- [ ] **第 2 步：验证文件存在且结构正确**

运行：`Rscript -e "testthat::test_file('tests/testthat/test-interp_fit.R')"`
预期结果：一个空的测试运行，零个结果。

- [ ] **第 3 步：提交**

```bash
git add tests/testthat/test-interp_fit.R
git commit -m "test: add test skeleton for interp_fit"
```

---

### 任务 3：实现 `._validate_xy()` 内部辅助函数并编写测试

**文件：**
- 创建：`R/interp_fit.R`
- 修改：`tests/testthat/test-interp_fit.R`（追加）

**依赖项：**
- 无

**产出：**
- 实现 `._validate_xy(x, y, min_len = 2)` 内部函数
- 所有校验级别的测试均通过

- [ ] **第 1 步：编写失败的测试**

追加到 `tests/testthat/test-interp_fit.R`：

```r
# --- ._validate_xy ---

test_that("._validate_xy() passes for valid input", {
  expect_no_error(mathmodels:::`._validate_xy`(c(1,2,3), c(4,5,6)))
})

test_that("._validate_xy() rejects non-numeric x", {
  expect_error(mathmodels:::`._validate_xy`(c("a","b"), c(1,2)),
               "x and y must be numeric")
})

test_that("._validate_xy() rejects non-numeric y", {
  expect_error(mathmodels:::`._validate_xy`(c(1,2), c("a","b")),
               "x and y must be numeric")
})

test_that("._validate_xy() rejects unequal lengths", {
  expect_error(mathmodels:::`._validate_xy`(c(1,2,3), c(1,2)),
               "x and y must have the same length")
})

test_that("._validate_xy() rejects too-short input", {
  expect_error(mathmodels:::`._validate_xy`(c(1), c(2), min_len = 2),
               "x and y must have length >= 2")
})

test_that("._validate_xy() rejects duplicate x values", {
  expect_error(mathmodels:::`._validate_xy`(c(1,1,2), c(1,2,3)),
               "x must not contain duplicate values")
})

test_that("._validate_xy() rejects NA in x", {
  expect_error(mathmodels:::`._validate_xy`(c(1,NA,3), c(4,5,6)),
               "x and y must not contain NA")
})

test_that("._validate_xy() rejects Inf in y", {
  expect_error(mathmodels:::`._validate_xy`(c(1,2,3), c(1,Inf,3)),
               "x and y must not contain")
})
```

- [ ] **第 2 步：运行测试，预期失败**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`
预期结果：8 个错误，均与 `object '._validate_xy' not found` 相关。

- [ ] **第 3 步：创建 `R/interp_fit.R` 并实现 `._validate_xy`**

```r
# =============================================================================
# Interpolation and curve fitting for mathmodels
# Implements: linear/poly/spline/hermite interpolation,
#             polynomial fit, linearizable curve fit, nonlinear growth fit
# =============================================================================

# -----------------------------------------------------------------------------
# Internal Helpers
# -----------------------------------------------------------------------------

# Validate x, y input pair
._validate_xy = function(x, y, min_len = 2) {
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors.", call. = FALSE)
  }
  if (length(x) != length(y)) {
    stop("x and y must have the same length.", call. = FALSE)
  }
  if (length(x) < min_len) {
    stop(sprintf("x and y must have length >= %d.", min_len), call. = FALSE)
  }
  if (anyDuplicated(x) != 0) {
    stop("x must not contain duplicate values.", call. = FALSE)
  }
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("x and y must not contain NA, Inf, or NaN.", call. = FALSE)
  }
  invisible(TRUE)
}
```

- [ ] **第 4 步：运行测试——预期全部通过**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`
预期结果：所有 8 个测试通过。

- [ ] **第 5 步：提交**

```bash
git add R/interp_fit.R tests/testthat/test-interp_fit.R
git commit -m "feat: add ._validate_xy internal helper"
```

---

### 任务 4：实现 `interp_linear()` 并编写测试

**文件：**
- 修改：`R/interp_fit.R`（追加约 50 行）
- 修改：`tests/testthat/test-interp_fit.R`（追加）

**依赖项：**
- `._validate_xy`（任务 3）

**产出：**
- 导出 `interp_linear(x, y, xout)`，返回 `tibble(x, y)` + `attr("method") = "linear"`
- 测试：正确性、边缘情况、输入校验

- [ ] **第 1 步：编写失败的测试**

追加到 `tests/testthat/test-interp_fit.R`：

```r
# --- interp_linear ---

test_that("interp_linear() returns expected structure", {
  res = interp_linear(c(1,2,3), c(2,4,6), xout = c(1, 1.5, 2, 2.5, 3))
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("x", "y"))
  expect_identical(attr(res, "method"), "linear")
  expect_equal(nrow(res), 5)
})

test_that("interp_linear() passes through known points", {
  res = interp_linear(c(1,3,5), c(10,30,50), xout = c(1,3,5))
  expect_equal(res$y, c(10, 30, 50))
})

test_that("interp_linear() interpolates midpoints", {
  res = interp_linear(c(0,2), c(0,4), xout = 1)
  expect_equal(res$y, 2)
})

test_that("interp_linear() extrapolates outside range (constant)", {
  res = interp_linear(c(1,2), c(1,2), xout = 0)
  expect_true(is.na(res$y))
})

test_that("interp_linear() validates input via ._validate_xy", {
  expect_error(interp_linear(c(1,1), c(1,2), xout = c(1,2)), "duplicate")
})

test_that("interp_linear() handles single xout", {
  res = interp_linear(c(1,2,3), c(10,20,30), xout = 2)
  expect_equal(nrow(res), 1)
  expect_equal(res$x, 2)
  expect_equal(res$y, 20)
})
```

- [ ] **第 2 步：运行测试——预期 6 个失败**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`
预期结果：`interp_linear` 测试失败并显示 "could not find function"。

- [ ] **第 3 步：实现 `interp_linear`**

追加到 `R/interp_fit.R` 中，位于内部辅助函数之后：

```r
# -----------------------------------------------------------------------------
# Interpolation
# -----------------------------------------------------------------------------

#' Linear Interpolation
#'
#' Performs linear (piecewise) interpolation at specified points.
#'
#' @param x Numeric vector of observed x-values (no duplicates).
#' @param y Numeric vector of observed y-values, same length as `x`.
#' @param xout Numeric vector of x-values where interpolation is desired.
#'
#' @return A tibble with columns `x` and `y`, containing the interpolated
#'   values at `xout`. The attribute `method` is set to `"linear"`.
#'
#' @examples
#' x = c(1, 2, 3, 4, 5)
#' y = c(2, 4, 1, 5, 3)
#' interp_linear(x, y, xout = c(1.5, 2.5, 3.5))
#'
#' @export
interp_linear = function(x, y, xout) {
  ._validate_xy(x, y, min_len = 2)
  yout = stats::approx(x, y, xout = xout)$y
  res = tibble::tibble(x = xout, y = yout)
  attr(res, "method") = "linear"
  res
}
```

- [ ] **第 4 步：运行测试——预期的 14 个测试（8 个 validate 测试 + 6 个 linear 测试），全部通过**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`
预期结果：14 个测试，0 个失败。

- [ ] **第 5 步：提交**

```bash
git add R/interp_fit.R tests/testthat/test-interp_fit.R
git commit -m "feat: add interp_linear()"
```

---

### 任务 5：实现 `interp_spline()` 并编写测试

**文件：**
- 修改：`R/interp_fit.R`（追加约 40 行）
- 修改：`tests/testthat/test-interp_fit.R`（追加）

**依赖项：**
- `._validate_xy`（任务 3）

**产出：**
- 导出 `interp_spline(x, y, xout, spar = NULL)`，返回 `tibble(x, y)` + `attr("method")`

- [ ] **第 1 步：编写失败的测试**

追加到 `tests/testthat/test-interp_fit.R`：

```r
# --- interp_spline ---

test_that("interp_spline() returns expected structure", {
  res = interp_spline(c(1,2,3,4,5), c(2,4,1,5,3), xout = c(1,2,3))
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("x", "y"))
  expect_identical(attr(res, "method"), "spline")
})

test_that("interp_spline() passes through known points", {
  res = interp_spline(c(1,3,5,7), c(10,30,50,70), xout = c(1,3,5,7))
  expect_equal(res$y, c(10, 30, 50, 70))
})

test_that("interp_spline() interpolates between points", {
  res = interp_spline(c(0,2), c(0,4), xout = 1)
  expect_gt(res$y, 0)
  expect_lt(res$y, 4)
})

test_that("interp_spline() with spar uses smooth.spline", {
  res = interp_spline(c(1,2,3,4,5), c(2.1,3.9,5.8,8.2,9.9), xout = c(1,3,5), spar = 0.5)
  expect_identical(attr(res, "method"), "smooth_spline")
  expect_equal(nrow(res), 3)
})

test_that("interp_spline() requires at least 4 points", {
  expect_error(interp_spline(c(1,2,3), c(1,2,3), xout = c(1,2)),
               "length >= 4")
})
```

- [ ] **第 2 步：运行测试——预期 5 个失败**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`

- [ ] **第 3 步：实现 `interp_spline`**

追加到 `R/interp_fit.R` 中，位于 `interp_linear` 之后：

```r
#' Cubic Spline Interpolation
#'
#' Performs cubic spline interpolation. When `spar` is `NULL`, uses
#' `stats::spline()` with the FMM method. Otherwise uses
#' `stats::smooth.spline()` for smoothing.
#'
#' @param x Numeric vector of observed x-values (no duplicates, length >= 4).
#' @param y Numeric vector of observed y-values.
#' @param xout Numeric vector of x-values where interpolation is desired.
#' @param spar Smoothing parameter passed to `smooth.spline()`. `NULL` for
#'   standard (non-smoothing) spline.
#'
#' @return A tibble with columns `x` and `y`. The attribute `method` is set to
#'   `"spline"` or `"smooth_spline"` depending on the variant used.
#'
#' @examples
#' x = c(1, 2, 3, 4, 5)
#' y = c(2, 4, 1, 5, 3)
#' interp_spline(x, y, xout = c(1.5, 2.5, 3.5))
#' interp_spline(x, y, xout = c(1.5, 2.5, 3.5), spar = 0.8)
#'
#' @export
interp_spline = function(x, y, xout, spar = NULL) {
  ._validate_xy(x, y, min_len = 4)
  if (is.null(spar)) {
    yout = stats::spline(x, y, xout = xout, method = "fmm")$y
    method = "spline"
  } else {
    fit = stats::smooth.spline(x, y, spar = spar)
    yout = predict(fit, xout)$y[, 1]
    method = "smooth_spline"
  }
  res = tibble::tibble(x = xout, y = yout)
  attr(res, "method") = method
  res
}
```

- [ ] **第 4 步：运行测试——预期的 19 个测试，全部通过**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`

- [ ] **第 5 步：提交**

```bash
git add R/interp_fit.R tests/testthat/test-interp_fit.R
git commit -m "feat: add interp_spline()"
```

---

### 任务 6：实现 `interp_poly()` 和 `interp_hermite()` 并编写测试

**文件：**
- 修改：`R/interp_fit.R`（追加约 70 行）
- 修改：`tests/testthat/test-interp_fit.R`（追加）

**依赖项：**
- `._validate_xy`（任务 3）

**产出：**
- 导出 `interp_poly(x, y, xout, degree = 3)`
- 导出 `interp_hermite(x, y, xout)`
- 两个函数的测试

- [ ] **第 1 步：编写失败的测试**

追加到 `tests/testthat/test-interp_fit.R`：

```r
# --- interp_poly ---

test_that("interp_poly() returns expected structure", {
  res = interp_poly(c(1,2,3,4,5), c(2,4,1,5,3), xout = c(1,2,3), degree = 2)
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("x", "y"))
  expect_identical(attr(res, "method"), "poly")
})

test_that("interp_poly() degree=1 fits a straight line", {
  res = interp_poly(c(0,2,4), c(0,4,8), xout = c(0,2,4), degree = 1)
  expect_equal(res$y, c(0, 4, 8))
})

test_that("interp_poly() degree too high for data triggers lm error", {
  expect_error(interp_poly(c(1,2,3), c(1,2,3), xout = c(1,2,3), degree = 5))
})

test_that("interp_poly() validates min_len based on degree", {
  expect_error(interp_poly(c(1,2), c(1,2), xout = c(1,2), degree = 3),
               "length >= 4")
})

# --- interp_hermite ---

test_that("interp_hermite() returns expected structure", {
  res = interp_hermite(c(1,2,3,4,5), c(2,4,1,5,3), xout = c(1,2,3))
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("x", "y"))
  expect_identical(attr(res, "method"), "hermite")
})

test_that("interp_hermite() passes through known points", {
  res = interp_hermite(c(1,3,5,7), c(10,30,50,70), xout = c(1,3,5,7))
  expect_equal(res$y, c(10, 30, 50, 70))
})

test_that("interp_hermite() preserves monotonicity", {
  x = c(1, 2, 3, 4, 5)
  y = c(1, 2, 4, 7, 11)
  res = interp_hermite(x, y, xout = seq(1, 5, by = 0.1))
  expect_true(all(diff(res$y) >= -1e-10))
})

test_that("interp_hermite() validates input", {
  expect_error(interp_hermite(c(1,1,3), c(1,2,3), xout = c(1,2)), "duplicate")
})
```

- [ ] **第 2 步：运行测试——预期 8 个失败**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`

- [ ] **第 3 步：实现两个函数**

追加到 `R/interp_fit.R` 中，位于 `interp_spline` 之后：

```r
#' Polynomial Interpolation
#'
#' Fits a polynomial of specified degree via `lm()` and evaluates it at
#' `xout`.
#'
#' @param x Numeric vector of observed x-values (no duplicates).
#' @param y Numeric vector of observed y-values.
#' @param xout Numeric vector of x-values where interpolation is desired.
#' @param degree Integer, degree of the polynomial. Minimum 1.
#'
#' @return A tibble with columns `x` and `y`. The attribute `method` is set
#'   to `"poly"`.
#'
#' @examples
#' x = c(1, 2, 3, 4, 5)
#' y = c(2, 4, 1, 5, 3)
#' interp_poly(x, y, xout = c(1.5, 2.5, 3.5), degree = 3)
#'
#' @export
interp_poly = function(x, y, xout, degree = 3) {
  ._validate_xy(x, y, min_len = degree + 1)
  fit = lm(y ~ poly(x, degree, raw = TRUE))
  yout = predict(fit, newdata = data.frame(x = xout))
  res = tibble::tibble(x = xout, y = as.numeric(yout))
  attr(res, "method") = "poly"
  res
}

#' Piecewise Cubic Hermite Interpolation
#'
#' Shape-preserving (monotone) piecewise cubic Hermite interpolation using
#' the Fritsch-Carlson method via `splinefun(method = "monoH.FC")`.
#'
#' @param x Numeric vector of observed x-values (no duplicates, length >= 2).
#' @param y Numeric vector of observed y-values.
#' @param xout Numeric vector of x-values where interpolation is desired.
#'
#' @return A tibble with columns `x` and `y`. The attribute `method` is set
#'   to `"hermite"`.
#'
#' @examples
#' x = c(1, 2, 3, 4, 5)
#' y = c(1, 2, 4, 7, 11)
#' interp_hermite(x, y, xout = c(1.5, 2.5, 3.5))
#'
#' @export
interp_hermite = function(x, y, xout) {
  ._validate_xy(x, y, min_len = 2)
  sf = stats::splinefun(x, y, method = "monoH.FC")
  yout = sf(xout)
  res = tibble::tibble(x = xout, y = yout)
  attr(res, "method") = "hermite"
  res
}
```

- [ ] **第 4 步：运行测试——预期 27 个测试，全部通过**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`

- [ ] **第 5 步：提交**

```bash
git add R/interp_fit.R tests/testthat/test-interp_fit.R
git commit -m "feat: add interp_poly() and interp_hermite()"
```

---

### 任务 7：实现拟合辅助函数（`._compute_fit_stats`, `._tidy_fit_coefs`, `._build_fit_residuals`, `._assemble_return`）并编写测试

**文件：**
- 修改：`R/interp_fit.R`（追加约 70 行）
- 修改：`tests/testthat/test-interp_fit.R`（追加）

**依赖项：**
- `._validate_xy`（任务 3）

**产出：**
- 4 个内部辅助函数的实现，用于组装统一的拟合返回值
- 通过测试，验证其输出结构

- [ ] **第 1 步：编写测试**

追加到 `tests/testthat/test-interp_fit.R`：

```r
# --- fitting helpers ---

test_that("._compute_fit_stats() returns expected columns", {
  obs = c(1, 2, 3)
  fitted = c(1.1, 2.2, 2.9)
  info = mathmodels:::._compute_fit_stats(obs, fitted, df_residual = 1)
  expect_s3_class(info, "tbl_df")
  expect_named(info, c("r.squared", "adj.r.squared", "aic", "bic", "sigma"))
  expect_true(info$r.squared >= 0 && info$r.squared <= 1)
})

test_that("._compute_fit_stats() gives r.squared=1 for perfect fit", {
  obs = c(1, 2, 3)
  fitted = c(1, 2, 3)
  info = mathmodels:::._compute_fit_stats(obs, fitted, df_residual = 0)
  expect_equal(info$r.squared, 1)
  expect_equal(info$sigma, 0)
})

test_that("._tidy_fit_coefs() tidy from lm model", {
  fit = lm(y ~ x, data = data.frame(x = 1:10, y = 2 * (1:10) + rnorm(10, sd = 0.1)))
  coefs = mathmodels:::._tidy_fit_coefs(fit)
  expect_s3_class(coefs, "tbl_df")
  expect_named(coefs, c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high"))
  expect_equal(coefs$term[1], "(Intercept)")
})

test_that("._build_fit_residuals() returns expected structure", {
  resid_tbl = mathmodels:::._build_fit_residuals(c(1, 2, 3), c(1.1, 2.2, 2.9))
  expect_s3_class(resid_tbl, "tbl_df")
  expect_named(resid_tbl, c(".observed", ".fitted", ".residual"))
  expect_equal(resid_tbl$.observed, c(1, 2, 3))
  expect_equal(resid_tbl$.residual[1], 1 - 1.1)
})
```

- [ ] **第 2 步：运行测试——预期 4 个失败**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`

- [ ] **第 3 步：实现 4 个辅助函数**

追加到 `R/interp_fit.R`，位于验证部分且位于插值函数之前：

```r
# Compute fit stats on the original scale
._compute_fit_stats = function(observed, fitted, df_residual) {
  ss_res = sum((observed - fitted)^2)
  ss_tot = sum((observed - mean(observed))^2)
  r_sq   = if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  n      = length(observed)
  p      = n - df_residual
  adj_r_sq = if (ss_tot > 0 && df_residual > 0) {
    1 - (1 - r_sq) * (n - 1) / df_residual
  } else NA_real_
  sigma = if (df_residual > 0) sqrt(ss_res / df_residual) else 0
  nll   = -sum(stats::dnorm(observed, mean = fitted, sd = sigma, log = TRUE))
  aic   = 2 * p + 2 * nll
  bic   = log(n) * p + 2 * nll
  tibble::tibble(
    r.squared     = round(r_sq, 4),
    adj.r.squared = round(adj_r_sq, 4),
    aic           = round(aic, 4),
    bic           = round(bic, 4),
    sigma         = round(sigma, 4)
  )
}

# Build tidy coefficient tibble from lm or nls model
._tidy_fit_coefs = function(model, level = 0.95) {
  coefs = stats::coef(model)
  se    = sqrt(diag(stats::vcov(model)))
  df_residual = if (inherits(model, "lm")) {
    model$df.residual
  } else {
    length(stats::residuals(model)) - length(coefs)
  }
  t_crit = stats::qt((1 - level) / 2, df_residual, lower.tail = FALSE)
  tibble::tibble(
    term      = names(coefs),
    estimate  = as.numeric(coefs),
    std.error = as.numeric(se),
    statistic = as.numeric(coefs) / as.numeric(se),
    p.value   = 2 * stats::pt(abs(as.numeric(coefs) / as.numeric(se)),
                               df_residual, lower.tail = FALSE),
    conf.low  = as.numeric(coefs) - t_crit * as.numeric(se),
    conf.high = as.numeric(coefs) + t_crit * as.numeric(se)
  )
}

# Build residuals tibble
._build_fit_residuals = function(observed, fitted) {
  tibble::tibble(
    .observed = observed,
    .fitted   = fitted,
    .residual = observed - fitted
  )
}

# Assemble unified return value
._assemble_return = function(model, observed, fitted, formula, type, x, y) {
  coefs  = ._tidy_fit_coefs(model)
  df_residual = if (inherits(model, "lm")) {
    model$df.residual
  } else {
    length(fitted) - length(stats::coef(model))
  }
  info   = ._compute_fit_stats(observed, fitted, df_residual)
  resid  = ._build_fit_residuals(observed, fitted)
  list(
    model       = model,
    coefficient = coefs,
    model_info  = info,
    residuals   = resid,
    formula     = formula,
    type        = type,
    input       = tibble::tibble(x = x, y = y)
  )
}
```

- [ ] **第 4 步：运行测试——预期 31 个测试，全部通过**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`

- [ ] **第 5 步：提交**

```bash
git add R/interp_fit.R tests/testthat/test-interp_fit.R
git commit -m "feat: add fitting internal helpers"
```

---

### 任务 8：实现 `poly_fit()` 并编写测试

**文件：**
- 修改：`R/interp_fit.R`（追加约 45 行）
- 修改：`tests/testthat/test-interp_fit.R`（追加）

**依赖项：**
- 所有内部辅助函数（任务 3、任务 7）

**产出：**
- 导出 `poly_fit(x, y, degree = 1)`，返回统一的拟合列表
- 测试通过

- [ ] **第 1 步：编写测试**

追加到 `tests/testthat/test-interp_fit.R`：

```r
# --- poly_fit ---

test_that("poly_fit() returns expected list structure", {
  set.seed(42)
  x = 1:20
  y = 3 + 2*x + rnorm(20, sd = 1)
  res = poly_fit(x, y, degree = 1)
  expect_type(res, "list")
  expect_setequal(names(res), c("model", "coefficient", "model_info", "residuals", "formula", "type", "input"))
  expect_s3_class(res$model_info, "tbl_df")
  expect_s3_class(res$coefficient, "tbl_df")
  expect_s3_class(res$residuals, "tbl_df")
  expect_equal(res$type, "poly")
})

test_that("poly_fit() degree=2 recovers quadratic coefficients", {
  set.seed(42)
  x = seq(-5, 5, length.out = 30)
  y = 1 + 2*x + 0.5*x^2 + rnorm(30, sd = 0.5)
  res = poly_fit(x, y, degree = 2)
  expect_true(abs(res$coefficient$estimate[3] - 0.5) < 0.2)
})

test_that("poly_fit() validates input", {
  expect_error(poly_fit(c(1,2), c(1,2), degree = 3), "length >= 4")
})

test_that("poly_fit() residuals match data length", {
  x = 1:10
  y = 2 + 3*x + rnorm(10, sd = 0.5)
  res = poly_fit(x, y, degree = 1)
  expect_equal(nrow(res$residuals), 10)
  expect_equal(nrow(res$input), 10)
})
```

- [ ] **第 2 步：运行测试——预期 4 个失败**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`

- [ ] **第 3 步：实现**

追加到 `R/interp_fit.R` 中，位于插值函数之后：

```r
# -----------------------------------------------------------------------------
# Curve Fitting
# -----------------------------------------------------------------------------

#' Polynomial Fit
#'
#' Fits a polynomial regression model via `lm()`.
#'
#' @param x Numeric vector of observed predictor values.
#' @param y Numeric vector of observed response values.
#' @param degree Integer, degree of the polynomial. Default is 1 (linear).
#'
#' @return A named list with elements: `model` (lm object), `coefficient`
#'   (tibble with estimates, SE, t-value, p-value, 95% CI), `model_info`
#'   (tibble with r.squared, adj.r.squared, AIC, BIC, sigma), `residuals`
#'   (tibble of .observed, .fitted, .residual), `formula` (character),
#'   `type` (character), and `input` (tibble of x, y).
#'
#' @examples
#' x = 1:20
#' y = 3 + 2*x + rnorm(20, sd = 2)
#' poly_fit(x, y, degree = 1)
#' poly_fit(x, y, degree = 2)
#'
#' @export
poly_fit = function(x, y, degree = 1) {
  ._validate_xy(x, y, min_len = degree + 1)
  fit = lm(y ~ poly(x, degree, raw = TRUE))
  fitted_vals = as.numeric(stats::fitted(fit))
  formula_str = paste0("y ~ poly(x, ", degree, ", raw = TRUE)")
  ._assemble_return(fit, y, fitted_vals, formula_str, "poly", x, y)
}
```

- [ ] **第 4 步：运行测试——预期 35 个测试，全部通过**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`

- [ ] **第 5 步：提交**

```bash
git add R/interp_fit.R tests/testthat/test-interp_fit.R
git commit -m "feat: add poly_fit()"
```

---

### 任务 9：实现 `curve_fit()` 并编写测试

**文件：**
- 修改：`R/interp_fit.R`（追加约 90 行）
- 修改：`tests/testthat/test-interp_fit.R`（追加）

**依赖项：**
- 所有内部辅助函数（任务 3、任务 7）
- 由 `lm()` 驱动，无需额外包

**产出：**
- 导出 `curve_fit(x, y, type = c("exp", "power", "log", "hyperbolic"))`
- 4 种曲线类型的测试

- [ ] **第 1 步：编写测试**

追加到 `tests/testthat/test-interp_fit.R`：

```r
# --- curve_fit ---

test_that("curve_fit() returns expected list structure", {
  set.seed(42)
  x = 1:20
  y = 2 * exp(0.15 * x) * exp(rnorm(20, sd = 0.1))
  res = curve_fit(x, y, type = "exp")
  expect_type(res, "list")
  expect_setequal(names(res), c("model", "coefficient", "model_info", "residuals", "formula", "type", "input"))
  expect_equal(res$type, "exp")
})

test_that("curve_fit(type='exp') recovers parameters", {
  set.seed(42)
  x = seq(0, 5, length.out = 50)
  y = 2 * exp(0.5 * x) * exp(rnorm(50, sd = 0.05))
  res = curve_fit(x, y, type = "exp")
  expect_true(abs(res$coefficient$estimate[1] - log(2)) < 0.3)
  expect_true(abs(res$coefficient$estimate[2] - 0.5) < 0.1)
})

test_that("curve_fit(type='power') recovers parameters", {
  set.seed(42)
  x = seq(1, 10, length.out = 50)
  y = 3 * x^2 * exp(rnorm(50, sd = 0.1))
  res = curve_fit(x, y, type = "power")
  expect_true(abs(exp(res$coefficient$estimate[1]) - 3) < 1)
  expect_true(abs(res$coefficient$estimate[2] - 2) < 0.3)
})

test_that("curve_fit(type='log') recovers parameters", {
  set.seed(42)
  x = seq(1, 20, length.out = 50)
  y = 1 + 3 * log(x) + rnorm(50, sd = 0.3)
  res = curve_fit(x, y, type = "log")
  expect_true(abs(res$coefficient$estimate[1] - 1) < 0.5)
  expect_true(abs(res$coefficient$estimate[2] - 3) < 0.3)
})

test_that("curve_fit(type='hyperbolic') recovers parameters", {
  set.seed(42)
  x = seq(1, 20, length.out = 50)
  y = 5 + 10/x + rnorm(50, sd = 0.3)
  res = curve_fit(x, y, type = "hyperbolic")
  expect_true(abs(res$coefficient$estimate[1] - 5) < 1)
  expect_true(abs(res$coefficient$estimate[2] - 10) < 2)
})

test_that("curve_fit() rejects negative y for exp type", {
  expect_error(curve_fit(1:5, c(-1, 1, 2, 3, 4), type = "exp"),
               "y must be positive")
})

test_that("curve_fit() rejects non-positive x for power type", {
  expect_error(curve_fit(c(1, -2, 3), c(1, 2, 3), type = "power"),
               "x must be positive")
})

test_that("curve_fit() rejects invalid type", {
  expect_error(curve_fit(1:5, 1:5, type = "unknown"),
               "should be one of")
})
```

- [ ] **第 2 步：运行测试——预期 8 个失败**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`

- [ ] **第 3 步：实现**

追加到 `R/interp_fit.R`，位于 `poly_fit` 之后：

```r
#' Linearizable Curve Fit
#'
#' Fits common empirical curves via variable transformation and `lm()`.
#' Supported types: exponential, power-law, logarithmic, and hyperbolic.
#'
#' @param x Numeric vector of observed predictor values.
#' @param y Numeric vector of observed response values.
#' @param type Character, one of `"exp"`, `"power"`, `"log"`, `"hyperbolic"`.
#'
#' @details
#' Each curve is transformed to a linear form:
#' - `"exp"`: fits `log(y) ~ x`, back-transforms to `y = a * exp(b * x)`
#' - `"power"`: fits `log(y) ~ log(x)`, back-transforms to `y = a * x^b`
#' - `"log"`: fits `y ~ log(x)`, gives `y = a + b * log(x)`
#' - `"hyperbolic"`: fits `y ~ 1/x`, gives `y = a + b / x`
#'
#' Fit statistics are computed on the original scale.
#'
#' @return A named list, see [poly_fit()] for details.
#'
#' @examples
#' set.seed(42)
#' x = 1:20
#' y = 2 * exp(0.15 * x) * rlnorm(20, sdlog = 0.1)
#' curve_fit(x, y, type = "exp")
#'
#' @export
curve_fit = function(x, y, type = c("exp", "power", "log", "hyperbolic")) {
  type = match.arg(type)
  ._validate_xy(x, y, min_len = 3)

  if (type == "exp") {
    if (any(y <= 0)) stop("y must be positive for exponential fit.", call. = FALSE)
    ly = log(y)
    fit = lm(ly ~ x)
    coefs = stats::coef(fit)
    a_hat = exp(coefs[1])
    b_hat = coefs[2]
    fitted_vals = a_hat * exp(b_hat * x)
    formula_str = sprintf("y = %.4g * exp(%.4g * x)", a_hat, b_hat)
  } else if (type == "power") {
    if (any(x <= 0)) stop("x must be positive for power fit.", call. = FALSE)
    if (any(y <= 0)) stop("y must be positive for power fit.", call. = FALSE)
    lx = log(x); ly = log(y)
    fit = lm(ly ~ lx)
    coefs = stats::coef(fit)
    a_hat = exp(coefs[1])
    b_hat = coefs[2]
    fitted_vals = a_hat * x^b_hat
    formula_str = sprintf("y = %.4g * x^(%.4g)", a_hat, b_hat)
  } else if (type == "log") {
    if (any(x <= 0)) stop("x must be positive for logarithmic fit.", call. = FALSE)
    lx = log(x)
    fit = lm(y ~ lx)
    coefs = stats::coef(fit)
    a_hat = coefs[1]; b_hat = coefs[2]
    fitted_vals = a_hat + b_hat * log(x)
    formula_str = sprintf("y = %.4g + %.4g * log(x)", a_hat, b_hat)
  } else if (type == "hyperbolic") {
    if (any(x == 0)) stop("x must be non-zero for hyperbolic fit.", call. = FALSE)
    inv_x = 1/x
    fit = lm(y ~ inv_x)
    coefs = stats::coef(fit)
    a_hat = coefs[1]; b_hat = coefs[2]
    fitted_vals = a_hat + b_hat / x
    formula_str = sprintf("y = %.4g + %.4g / x", a_hat, b_hat)
  }

  ._assemble_return(fit, y, fitted_vals, formula_str, type, x, y)
}
```

- [ ] **第 4 步：运行测试——预期 43 个测试，全部通过**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`

- [ ] **第 5 步：提交**

```bash
git add R/interp_fit.R tests/testthat/test-interp_fit.R
git commit -m "feat: add curve_fit()"
```

---

### 任务 10：实现 `growth_fit()` 并编写测试

**文件：**
- 修改：`R/interp_fit.R`（追加约 120 行——最大的函数）
- 修改：`tests/testthat/test-interp_fit.R`（追加）

**依赖项：**
- 所有内部辅助函数（任务 3、任务 7）
- `minpack.lm::nlsLM()`（任务 1 中的依赖项）

**产出：**
- 导出 `growth_fit(x, y, type = c("logistic", "gompertz", "saturation", "mm"))`
- 4 种增长曲线类型的测试

- [ ] **第 1 步：编写测试**

追加到 `tests/testthat/test-interp_fit.R`：

```r
# --- growth_fit ---

test_that("growth_fit() returns expected list structure", {
  set.seed(42)
  x = seq(0, 10, length.out = 30)
  y = 100 / (1 + exp(-1.5 * (x - 5))) + rnorm(30, sd = 2)
  res = growth_fit(x, y, type = "logistic")
  expect_type(res, "list")
  expect_setequal(names(res), c("model", "coefficient", "model_info", "residuals", "formula", "type", "input"))
  expect_equal(res$type, "logistic")
  expect_equal(nrow(res$residuals), 30)
})

test_that("growth_fit(type='logistic') recovers parameters", {
  set.seed(42)
  x = seq(0, 10, length.out = 50)
  L_true = 100; k_true = 1.5; x0_true = 5
  y = L_true / (1 + exp(-k_true * (x - x0_true))) + rnorm(50, sd = 2)
  res = growth_fit(x, y, type = "logistic")
  est = res$coefficient$estimate
  expect_true(abs(est[1] - L_true) < 10)
  expect_true(abs(est[2] - k_true) < 0.5)
  expect_true(abs(est[3] - x0_true) < 1)
})

test_that("growth_fit(type='gompertz') recovers parameters", {
  set.seed(42)
  x = seq(0, 10, length.out = 50)
  L_true = 100; k_true = 0.5; x0_true = 4
  y = L_true * exp(-exp(-k_true * (x - x0_true))) + rnorm(50, sd = 2)
  res = growth_fit(x, y, type = "gompertz")
  est = res$coefficient$estimate
  expect_true(abs(est[1] - L_true) < 10)
  expect_true(abs(est[2] - k_true) < 0.3)
  expect_true(abs(est[3] - x0_true) < 1.5)
})

test_that("growth_fit(type='saturation') recovers parameters", {
  set.seed(42)
  x = seq(0, 10, length.out = 50)
  a_true = 50; b_true = 0.5
  y = a_true * (1 - exp(-b_true * x)) + rnorm(50, sd = 1)
  res = growth_fit(x, y, type = "saturation")
  est = res$coefficient$estimate
  expect_true(abs(est[1] - a_true) < 5)
  expect_true(abs(est[2] - b_true) < 0.2)
})

test_that("growth_fit(type='mm') recovers parameters", {
  set.seed(42)
  x = seq(0.5, 20, length.out = 50)
  a_true = 30; b_true = 3
  y = a_true * x / (b_true + x) + rnorm(50, sd = 0.5)
  res = growth_fit(x, y, type = "mm")
  est = res$coefficient$estimate
  expect_true(abs(est[1] - a_true) < 5)
  expect_true(abs(est[2] - b_true) < 2)
})

test_that("growth_fit() rejects invalid type", {
  expect_error(growth_fit(1:10, 1:10, type = "unknown"), "should be one of")
})

test_that("growth_fit() validates input", {
  expect_error(growth_fit(c(1,2), c(1,2), type = "logistic"), "length >= 3")
})
```

- [ ] **第 2 步：运行测试——预期 7 个失败**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`

- [ ] **第 3 步：实现**

追加到 `R/interp_fit.R`，位于 `curve_fit` 之后：

```r
#' Nonlinear Growth Curve Fit
#'
#' Fits nonlinear growth, saturation, and sigmoid curves via
#' `minpack.lm::nlsLM()` with automatic starting value generation.
#'
#' @param x Numeric vector of observed predictor values.
#' @param y Numeric vector of observed response values.
#' @param type Character, one of `"logistic"`, `"gompertz"`, `"saturation"`,
#'   or `"mm"`.
#'
#' @details
#' Models and starting value strategies:
#'
#' **Logistic**: `y = L / (1 + exp(-k * (x - x0)))`
#' - `L0 = max(y) * 1.05`
#' - `lm(log((L0 - y) / y) ~ x)` for `k, x0` initials
#'
#' **Gompertz**: `y = L * exp(-exp(-k * (x - x0)))`
#' - `L0 = max(y) * 1.05`
#' - `lm(log(-log(y / L0)) ~ x)` for `k, x0` initials
#'
#' **Saturation** (exponential approach): `y = a * (1 - exp(-b * x))`
#' - `a0 = max(y) * 1.1`
#' - `lm(log(a0 - y) ~ x)` for `b` initial
#'
#' **Michaelis-Menten**: `y = a * x / (b + x)`
#' - `a0 = max(y) * 1.1`
#' - Lineweaver-Burk `lm(1/y ~ 1/x)` for `a, b` initials
#'
#' Falls back to heuristic defaults if linearized estimation fails.
#'
#' @return A named list, see [poly_fit()] for details. `model` contains
#'   the `nls` object.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x = seq(0, 10, length.out = 30)
#' y = 100 / (1 + exp(-1.5 * (x - 5))) + rnorm(30, sd = 2)
#' growth_fit(x, y, type = "logistic")
#' }
#'
#' @export
growth_fit = function(x, y, type = c("logistic", "gompertz", "saturation", "mm")) {
  type = match.arg(type)
  ._validate_xy(x, y, min_len = 3)

  if (!requireNamespace("minpack.lm", quietly = TRUE)) {
    stop("Package 'minpack.lm' is required for growth_fit().", call. = FALSE)
  }

  nlsLM = minpack.lm::nlsLM

  if (type == "logistic") {
    L0  = max(y) * 1.05
    start = tryCatch({
      y_adj = pmax(y, 1e-8)
      y_adj = pmin(y_adj, L0 - 1e-8)
      ratio = log((L0 - y_adj) / y_adj)
      lin = stats::lm(ratio ~ x)
      co = stats::coef(lin)
      c(L = L0, k = -co[2], x0 = co[1] / (-co[2]))
    }, error = function(e) c(L = L0, k = 1/diff(range(x)), x0 = mean(x)))
    fit = nlsLM(y ~ L / (1 + exp(-k * (x - x0))),
                start = list(L = start[["L"]], k = start[["k"]], x0 = start[["x0"]]))
    fitted_vals = stats::fitted(fit)
    formula_str = sprintf("y = %.4g / (1 + exp(-%.4g * (x - %.4g)))",
                          stats::coef(fit)[["L"]], stats::coef(fit)[["k"]],
                          stats::coef(fit)[["x0"]])

  } else if (type == "gompertz") {
    L0  = max(y) * 1.05
    start = tryCatch({
      y_adj = pmax(y / L0, 1e-8)
      y_adj = pmin(y_adj, 1 - 1e-8)
      dublog = log(-log(y_adj))
      lin = stats::lm(dublog ~ x)
      co = stats::coef(lin)
      c(L = L0, k = -co[2], x0 = co[1] / (-co[2]))
    }, error = function(e) c(L = L0, k = 1/diff(range(x)), x0 = mean(x)))
    fit = nlsLM(y ~ L * exp(-exp(-k * (x - x0))),
                start = list(L = start[["L"]], k = start[["k"]], x0 = start[["x0"]]))
    fitted_vals = stats::fitted(fit)
    formula_str = sprintf("y = %.4g * exp(-exp(-%.4g * (x - %.4g)))",
                          stats::coef(fit)[["L"]], stats::coef(fit)[["k"]],
                          stats::coef(fit)[["x0"]])

  } else if (type == "saturation") {
    a0  = max(y) * 1.1
    start = tryCatch({
      diff_y = pmax(a0 - y, 1e-8)
      lin = stats::lm(log(diff_y) ~ x)
      co = stats::coef(lin)
      c(a = a0, b = -co[2])
    }, error = function(e) c(a = a0, b = 1/diff(range(x))))
    fit = nlsLM(y ~ a * (1 - exp(-b * x)),
                start = list(a = start[["a"]], b = start[["b"]]))
    fitted_vals = stats::fitted(fit)
    formula_str = sprintf("y = %.4g * (1 - exp(-%.4g * x))",
                          stats::coef(fit)[["a"]], stats::coef(fit)[["b"]])

  } else if (type == "mm") {
    a0  = max(y) * 1.1
    start = tryCatch({
      pos = x > 0 & y > 0
      lin = stats::lm(I(1/y[pos]) ~ I(1/x[pos]))
      co = stats::coef(lin)
      c(a = 1/co[1], b = co[2]/co[1])
    }, error = function(e) c(a = a0, b = stats::median(x)))
    fit = nlsLM(y ~ a * x / (b + x),
                start = list(a = start[["a"]], b = start[["b"]]))
    fitted_vals = stats::fitted(fit)
    formula_str = sprintf("y = %.4g * x / (%.4g + x)",
                          stats::coef(fit)[["a"]], stats::coef(fit)[["b"]])
  }

  ._assemble_return(fit, y, fitted_vals, formula_str, type, x, y)
}
```

- [ ] **第 4 步：运行测试——预期 50 个测试，全部通过**

运行：`Rscript -e "devtools::load_all(); devtools::test(filter = '^interp')"`

- [ ] **第 5 步：提交**

```bash
git add R/interp_fit.R tests/testthat/test-interp_fit.R
git commit -m "feat: add growth_fit() with auto-start-values"
```

---

### 任务 11：生成文档、包检查并更新 NEWS

**文件：**
- 修改：`NEWS.md`（前置新版本章节）
- 生成：`man/interp_linear.Rd`、`man/interp_spline.Rd`、`man/interp_poly.Rd`、`man/interp_hermite.Rd`
- 生成：`man/poly_fit.Rd`、`man/curve_fit.Rd`、`man/growth_fit.Rd`
- 修改：`NAMESPACE`（通过 roxygen2）

**依赖项：**
- 所有之前的任务

**产出：**
- 7 个帮助页面已构建
- NAMESPACE 已更新，包含 7 个新的导出项
- NEWS.md 已更新，包含此功能
- `devtools::check()` 通过，无错误或警告

- [ ] **第 1 步：通过 roxygen2 生成文档**

运行：`Rscript -e "devtools::document()"`
预期结果：无错误。NAMESPACE 现在包含 7 个新的 `export()` 条目；`man/` 包含 7 个新的 `.Rd` 文件。

- [ ] **第 2 步：检查生成的文档**

运行：`Rscript -e "devtools::document()"`
验证 `man/interp_linear.Rd`、`man/curve_fit.Rd`、`man/growth_fit.Rd` 等文件存在。

- [ ] **第 3 步：更新 NEWS.md**

将新版本章节前置到 `NEWS.md`：

```markdown
# mathmodels 0.0.11

## 插值与曲线拟合工具包（新增）

- **插值**（`interp_fit.R`）：
  - `interp_linear()`：分段线性插值。
  - `interp_spline()`：三次样条插值，可选通过 `spar` 进行平滑处理。
  - `interp_poly()`：给定次数的多项式插值。
  - `interp_hermite()`：保形分段三次 Hermite 插值。

- **曲线拟合**（`interp_fit.R`）：
  - `poly_fit()`：多项式回归。
  - `curve_fit()`：通过变量变换实现可线性化曲线的统一接口
    （指数、幂律、对数、双曲）。
  - `growth_fit()`：通过 `minpack.lm::nlsLM()` 进行非线性生长曲线拟合，
    自动生成初值（Logistic、Gompertz、指数饱和、Michaelis-Menten）。

- **依赖项**（新增）：`minpack.lm`，用于基于 Levenberg-Marquardt 的非线性最小二乘拟合。
```

- [ ] **第 4 步：运行完整的测试套件**

运行：`Rscript -e "devtools::test()"`
预期结果：所有现有测试 + 50 个新测试均通过，无失败。

- [ ] **第 5 步：运行包检查**

运行：`Rscript -e "devtools::check()"`
预期结果：0 个错误，0 个警告，0 个非 base 相关的注意事项。

- [ ] **第 6 步：格式化代码**

运行：`air format R/interp_fit.R tests/testthat/test-interp_fit.R`

- [ ] **第 7 步：提交**

```bash
git add man/*.Rd NAMESPACE NEWS.md R/interp_fit.R tests/testthat/test-interp_fit.R
git commit -m "docs: generate docs and update NEWS for interp_fit"
```

---

### 任务 12：最终清理与验证

**依赖项：**
- 所有之前的任务

**产出：**
- 干净的状态，所有内容已提交并通过

- [ ] **第 1 步：最终包检查**

运行：`Rscript -e "devtools::check()"`
预期结果：0 个错误，0 个警告。

- [ ] **第 2 步：检查 pkgdown 引用覆盖率**

运行：`Rscript -e "pkgdown::check_pkgdown()"`

- [ ] **第 3 步：如果 pkgdown 检查通过，则提交**

```bash
git add -A
git commit -m "chore: final review pass for interp_fit"
```
