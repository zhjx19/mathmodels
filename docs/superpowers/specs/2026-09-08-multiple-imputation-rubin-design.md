# 多重插补 + Rubin 合并（v1）设计文档

日期：2026-09-08
状态：已获批准
目标文件：`R/data_cleaning.R` 主题扩展（新功能另起文件）

## 背景与决策记录

`data_cleaning.R`（772 行）已有缺失值诊断（`na_summary()`、`na_mcar_test()` 含自研
EM 算法）、简单插补（`impute()`：mean/median/mode/linear/spline/constant）、模型插补
（`impute_model()`：knn/lm/tree）、异常值处理（`outlier()`）与透视表。本次新增多重
插补（MI）与 Rubin 规则合并，打通"插补 → 逐套拟合 → 合并推断"全流程。

已确认的决策：

1. **自研轻量 FCS 链式方程引擎**（方案 2），不导入 mice。数值列严格 proper
   （贝叶斯线性回归抽参 + pmm/norm）；分类列近似 proper（参数扰动，与 mice 的
   polyreg/logreg 同级）。mice 进 Suggests 仅作测试对照（运行时零依赖）。
2. **API 方案 B**：`impute_multiple()` + `rubin_pool()`（标量核心）+
   `pool_fit()`（lm/glm 系数一条龙）。
3. **方法矩阵**：数值 → `pmm`（默认）/`norm`；二分类 → `logreg`；多分类 →
   `polyreg`。v1 不含 polr、cart、rf。
4. 数据范围：数值 + 分类（含科研用途）。
5. 用户全局编码铁律优先：`=` 赋值、`|>` 管道、`\(x)` 匿名函数。

## 文件布局与依赖

- 新文件 `R/mi_impute.R`（~450 行）：`impute_multiple()` + FCS 内部函数
  （`.fcs_*` 族）+ `print.mi_impute()`。
- 新文件 `R/mi_pool.R`（~200 行）：`rubin_pool()` + `pool_fit()` +
  print 方法。
- DESCRIPTION：`Imports` 加 `nnet`（随 R 发行，multinom 用）；`Suggests` 加
  `mice`（测试对照）。MASS 已有，无需动。
- 测试：`tests/testthat/test-mi_impute.R`、`tests/testthat/test-mi_pool.R`。

## `impute_multiple()` 设计

### 接口

```r
impute_multiple(data, .cols = dplyr::everything(),
                m = 5, maxit = 5, seed = NULL, methods = NULL)
```

- `.cols`：tidy-select，默认所有列；只对含缺失的列执行插补，无缺失列原样保留。
- `m`：插补套数（默认 5，要求 >= 2）。
- `maxit`：链式迭代轮数（默认 5）。
- `seed`：可复现种子；m 套共用一个 RNG 流。
- `methods`：命名向量/列表覆盖自动选择。自动规则：数值 → `pmm`；因子/字符
  且水平数 = 2 → `logreg`；水平数 > 2 → `polyreg`。逻辑列视作二分类。

### FCS 引擎（内部 `.fcs_*` 函数族）

1. **初始化**：数值列均值、分类列众数填充（仅链式迭代起步用）。
2. **迭代**：每轮按列序依次，目标列当前缺失位置用其余列作回归元重填：
   - `pmm`：贝叶斯 lm 抽参（β|σ² ~ N(β̂, σ²(X'X)⁻¹)，σ² ~ scale × Inv-χ²）
     → 线性预测 → 预测均值匹配（从观测值中取预测最接近的 3 个 donor 的
     均值；mice 默认是从近邻随机抽一个 donor，此处为确定性均值变体，
     有助于 seed 复现与测试稳定性）。
   - `norm`：同参数抽取 + N(0, σ²) 残差噪声。
   - `logreg`：glm(binomial) → β ~ N(β̂, V̂) → 逆 logit → rbinom 抽类别。
   - `polyreg`：nnet::multinom → β ~ N(β̂, V̂)（V̂ 经 tryCatch 奇异保护，
     降级点预测 + 一次性警告）→ 按线性预测概率多项抽样。
3. **m 套独立链**：各自独立初始化、迭代、RNG 抽取，产 m 个完整数据集。
4. 因子特征经 `model.matrix` 哑元化进回归；完全共线/完全分离/拟合失败时
   tryCatch 跳过该列该迭代（保留上次填充值），累计计数后统一警告。
5. 收敛诊断：每迭代每列记录插补单元格的均值/SD（`trace` 表）。

### 返回值

`"mi_impute"` 对象：m 个 data frame 的 list，可直接 `lapply(imps, ...)`。
属性：`methods`（每列实际方法）、`iterations`、`trace`（tibble：imputation、
iteration、variable、mean、sd）、`seed`、`n_imputed`（每列插补单元格数）。
`print.mi_impute()`：简报（m、各列方法、插补格数、迭代数）。

## `rubin_pool()` 设计

### 接口

```r
rubin_pool(est, err, alpha = 0.05, edf = Inf)
```

参考 `kstools::ks_mianalyze`（SAS PROC MIANALYZE 同源）简化：删 `v0` 兼容
别名、`theta0`（p 值对 0 检验）、knit_print、四表 summary。保留全部统计核心。

### 统计量（snake_case 组件名）

- `estimate`：Q̄ = mean(est)；`std_error`：√T
- 方差分解：`var_within` Ū、`var_between` B、`var_total` T = Ū + (1+1/m)B
- `df`：Rubin (1987) 大样本 ν_m；`edf` 有限时 Barnard-Rubin (1999) 调整
  ν_m* = (1/ν_m + 1/ν_obs)⁻¹，ν_obs = (1-γ)ν₀(ν₀+1)/(ν₀+3)，
  γ = (1+1/m)B/T
- 诊断：`r`（方差相对增加）、`fmi`（缺失信息分数）、`rel_eff`（相对效率
  (1+λ/m)⁻¹）
- 推断：`lcl`/`ucl`（t 分位 CI）、`t`、`p_value`
- 边界处理：m=1 → B=0、r=0、df=Inf；err 全 0 → 同上退化路径；est/err 长度
  不等、含 NA/Inf、err<0、alpha 越界 → cli 风格报错（包内现有 stop 风格）。

`print` 方法：单表（Estimate/SE/CI/df/FMI/p）。

## `pool_fit()` 设计

### 接口

```r
pool_fit(fits, alpha = 0.05)
```

- `fits`：lm/glm 拟合列表（通常 `lapply(imps, \(d) lm(y ~ ., d))`）。
- 从每个 fit 提取 `coef()` 与 `sqrt(diag(vcov()))` → 逐系数 `rubin_pool()`。
- **自动传 `edf = fits[[1]]$df.residual`**（Barnard-Rubin 小样本调整自动生效）。
- 系数名跨套不对齐 → 报错；NA 系数（秩亏）跳过并警告。
- 返回 tibble：`term, estimate, std_error, lcl, ucl, df, t, p_value, fmi,
  r, rel_eff`，配简洁 print。

## 测试与验证

### `test-mi_pool.R`

- m=2 小例子手算闭式核对（Q̄、T、df、CI）。
- Barnard-Rubin 调整：构造已知 ν₀/γ 验证 ν_m* 公式。
- 边界：m=1、err 全 0、NA/Inf 输入、alpha 越界报错。
- **mice 对照**：`skip_if_not_installed("mice")`，同输入与
  `mice::pool` 逐数对照（estimate/se/p）。

### `test-mi_impute.R`

- 形状不变、类型不变、因子水平不变；无缺失列原样保留。
- seed 可复现；m 套互异；`.cols` 选择生效。
- 方法自动匹配（数值/二分类/多分类）；`methods` 覆盖生效；非法方法名报错。
- 插补值填充缺失位且不改动观测位。
- **mice 数值接近对照**：同一数据上 pmm/norm 合并结果与 mice 数值接近
  （容差断言，非逐位相等）——自研 proper 性的可信度背书。

## v1 明确不做

- `polr`（有序分类按 polyreg 处理）、cart/randomForest 方法。
- `plot.mi_impute` 收敛图（trace 表已存数据，图留 v2）。
- vignette 章节（后续按需）。
