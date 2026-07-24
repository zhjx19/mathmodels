# interp_fit.R — 插值与曲线拟合模块设计

**日期**: 2025-07-25  
**状态**: 待实现  
**文件**: `R/interp_fit.R`

---

## 概要

为 mathmodels 包新增插值和曲线拟合功能模块，所有代码集中在 `R/interp_fit.R`。采用方案 A（扁平分治）：内部共享少量 helper，对外暴露 7 个导出函数。

---

## 1. 依赖

**新增 Imports**（需加入 `DESCRIPTION`）:

- `minpack.lm` — Levenberg-Marquardt 非线性最小二乘 (`nlsLM()`)，比 `stats::nls()` 的 Gauss-Newton 算法更稳定，对自动生成初值误差容错更好

**已有依赖**（无需新增）:

- `stats` — `approx`, `spline`, `splinefun`, `lm`, `nls`, `AIC`, `BIC`, `coef`, `fitted`, `predict`, `residuals`, `vcov`
- `tibble` — 返回值
- `ggplot2`, `patchwork` — 画图（本次不含）

---

## 2. 插值函数（4 个导出函数）

### 公共规则

- **输入**: `x`, `y` 为数值向量，等长，`x` 无重复值
- **输出**: `tibble(x, y)`，`attr(return, "method")` 记录插值方法名
- **输入校验**: 检查 `NA`/`Inf`/`NaN`、`x` 最小长度

### 2.1 `interp_linear(x, y, xout)`

- **底层**: `stats::approx(x, y, xout)$y`
- **x 最小长度**: 2

### 2.2 `interp_spline(x, y, xout, spar = NULL)`

- **底层**: `spar` 为 `NULL` 时 `stats::spline(x, y, xout, method = "fmm")`；否则 `stats::smooth.spline(x, y, spar = spar)` 后 predict
- **x 最小长度**: 4
- **attr**: `"spline"` 或 `"smooth_spline"`

### 2.3 `interp_poly(x, y, xout, degree = 3)`

- **底层**: `lm(y ~ poly(x, degree, raw = TRUE))` 后用 `predict()` 在 `xout` 处求值
- **x 最小长度**: `degree + 1`
- **attr**: `"poly"`

### 2.4 `interp_hermite(x, y, xout)`

- **底层**: `stats::splinefun(x, y, method = "monoH.FC")(xout)`（Fritsch-Carlson 保形算法）
- **x 最小长度**: 2
- **attr**: `"hermite"`

---

## 3. 拟合函数（3 个导出函数）

### 公共返回值结构

所有拟合函数返回相同的 list 结构：

```r
list(
  model       = <lm | nls>          # 原始 R 模型对象
  coefficient = tibble(term, estimate, std.error, statistic, p.value, conf.low, conf.high),
  model_info  = tibble(r.squared, adj.r.squared, aic, bic, sigma)  # sigma = RSE
  residuals   = tibble(.observed, .fitted, .residual),
  formula     = <character>,        # 模型表达式字符串
  type        = <character>,        # 曲线类型标识
  input       = tibble(x, y)
)
```

`r.squared` 计算为 `1 - SS_res / SS_tot`（原始尺度上），`sigma` 为残差标准误 `sqrt(SS_res / df_residual)`。

### 3.1 `poly_fit(x, y, degree = 1)`

- **底层**: `lm(y ~ poly(x, degree, raw = TRUE))`
- **type**: `"poly"`
- **初值**: 不需要（线性）

### 3.2 `curve_fit(x, y, type = c("exp", "power", "log", "hyperbolic"))`

变量变换 + `lm()` 线性化，然后反变换回原参数。

| type | 模型 | 变换 | 反变换 | 前提 |
|------|------|------|--------|------|
| `"exp"` | `y = a·exp(b·x)` | `ln(y) ~ x` | `a=exp(intercept)`, `b=slope` | `y > 0` |
| `"power"` | `y = a·x^b` | `ln(y) ~ ln(x)` | `a=exp(intercept)`, `b=slope` | `x > 0, y > 0` |
| `"log"` | `y = a + b·ln(x)` | `y ~ ln(x)` | 直接 `lm` 结果 | `x > 0` |
| `"hyperbolic"` | `y = a + b/x` | `y ~ 1/x` | 直接 `lm` 结果 | `x != 0` |

- **注意**: `exp` 和 `power` 的 `a=exp(intercept)` 存在对数变换偏差，首版不做 smearing 校正（后续可加）。
- **model slot**: 存储变换后的 `lm` 对象。

### 3.3 `growth_fit(x, y, type = c("logistic", "gompertz", "saturation", "mm"))`

- **底层**: `minpack.lm::nlsLM()`（Levenberg-Marquardt 算法）
- **自动初值生成**（核心设计）：

| type | 模型 | 初值策略 |
|------|------|----------|
| `"logistic"` | `y = L / (1 + exp(-k·(x - x0)))` | `L0 = max(y) * 1.05`；`lm(log((L0-y)/y) ~ x)` 估计 `k, x0`；fallback: `c(L=max(y)*1.05, k=1/diff(range(x)), x0=mean(x))` |
| `"gompertz"` | `y = L · exp(-exp(-k·(x - x0)))` | `L0 = max(y) * 1.05`；`lm(log(-log(y/L0)) ~ x)` 估计 `k, x0`；fallback 同上 |
| `"saturation"` | `y = a · (1 - exp(-b·x))` | `a0 = max(y) * 1.1`；`lm(log(a0-y) ~ x)` 估计 `b`；fallback: `c(a=max(y)*1.1, b=1/diff(range(x)))` |
| `"mm"` | `y = a·x / (b + x)` | `a0 = max(y) * 1.1`；Lineweaver-Burk `lm(1/y ~ 1/x)` 估计 `a, b`；fallback: `c(a=max(y)*1.1, b=median(x))` |

初值容错：线性化估计失败时降级为 fallback，仍失败则 `stop()` 报错并建议用户手动提供初值。

---

## 4. 内部辅助函数

| 函数 | 用途 |
|------|------|
| `._validate_xy(x, y, min_len = 2)` | 输入校验：数值向量、等长、无 NA/Inf/NaN、`x` 无重复、长度 ≥ min_len |
| `._compute_fit_stats(observed, fitted, df_residual)` | 计算 `r.squared`, `adj.r.squared`, `sigma`, 用于返回值 `model_info` |
| `._tidy_fit_coefs(model)` | 从 lm/nls 提取系数 tibble（含 CI） |
| `._build_fit_residuals(observed, fitted)` | 构建残差 tibble |
| `._assemble_return(model, coefs, info, resid, formula, type, x, y)` | 统一组装返回值 list |

---

## 5. 输入校验

`._validate_xy()` 检查条目：

1. `x`, `y` 为数值向量（`is.numeric`, `is.vector`）
2. `length(x) == length(y)`
3. `length(x) >= min_len`
4. `x` 无重复值（`anyDuplicated(x) == 0`）
5. 无 `NA`, `Inf`, `NaN`（`all(is.finite())`）

---

## 6. 测试计划

新增 `tests/testthat/test-interp_fit.R`，覆盖：

- **插值**: 每种方法的正确性（已知点精确通过）、错误输入（NA、不等长、x 重复）、边界条件（xout 超出范围）
- **curve_fit**: 各类型变换正确性、前提条件验证（y ≤ 0 时 exp 报错等）
- **growth_fit**: 各类型初值自动生成、nlsLM 收敛、返回值结构、初值生成失败 fallback
- **poly_fit**: 不同 degree 的正确性、degree 过大时的行为

---

## 7. 不包含的内容（首版不实现）

- 画图函数（后续迭代）
- 画图函数（后续迭代）
- Richards、Hill 等模型（后续按需添加）
- 对数变换偏差的 smearing 校正
- 样条插值中 `spar` 的自动选择（cross-validation）

---

## 8. 命名空间

导出函数（加入 `NAMESPACE`）:

```
export(interp_linear)
export(interp_spline)
export(interp_poly)
export(interp_hermite)
export(poly_fit)
export(curve_fit)
export(growth_fit)
```

内部函数（`.` 前缀，不导出）:

```
._validate_xy
._compute_fit_stats
._tidy_fit_coefs
._build_fit_residuals
._assemble_return
```
