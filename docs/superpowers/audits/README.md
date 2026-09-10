# mathmodels 包级对照审核总报告

- **日期:** 2026-09-09
- **范围:** 全部自研算法（跳过纯封装现成函数的模块），对照权威参考包或独立重实现。
- **方法:** 同输入双跑对照（stats/forecast/ineq/naniar/Benchmarking/…）+ 教科书公式独立重实现 + 手算小例 + 直接似然优化验证。

## 审核结果总览

| 模块 | 函数数 | 方法 | 结果 | 严重 bug | 已修 |
|---|---|---|---|---|---|
| stat_infer | 13 | stats/coin/car 双跑 | 0 计算错误；7 Important 语义/契约 | 0 | 7 |
| pred_grey | 5 | 闭式解独立重实现 | 31/31 PASS | 0 | 0 |
| eval_inequality | 9 | ineq/DescTools + 分解恒等式 | 1 Critical | 1 | 1 |
| eval_weights | 8 | 独立公式 + Saaty 表 | 1 Important (RI[12]) | 0 | 1 |
| eval_fuzzy | 21 | 教科书形状/公式 | 38/38 PASS | 0 | 0 |
| eval_preprocess | 9 | 手算/公式 | 20/20 PASS；1 文档反转 | 0 | 1（文档） |
| eval_algorithms | 4 | 独立 TOPSIS/GRA/RSR | 匹配 | 0 | 0 |
| eval_system | 2 | 公式/不变量 | 匹配 | 0 | 0 |
| eval_regional | 3 | 独立 LQ/HHI/EG(1997) | 匹配 | 0 | 0 |
| data_cleaning（自研） | — | naniar + 直接 MLE | 匹配（na_mcar 0.9% 差异属实现细节） | 0 | 0 |
| interp_fit | 7 | 封装 stats/nlsLM | 抽查匹配 | 0 | 0 |
| eval_dea | basic_DEA | Benchmarking | 完全匹配（VRS/CRS） | 0 | 0 |
| diffEq_models | — | deSolve 封装 + 标准方程 | 守恒/闭式解通过 | 0 | 0 |
| mv_stats / pred_ts | — | 封装 stats/forecast | 按用户要求跳过 | — | — |

## 修复清单（commit）

| commit | 模块 | 问题 | 修复 |
|---|---|---|---|
| `c8d4686` | stat_infer | 7 Important：perm-test 自定义统计量语义、pairwise 相关、ANCOVA Type II、describe 逐变量删除、t 检验 estimate 符号等 | 全部修复 + 回归测试 |
| `f00454f` | eval_inequality | **Critical**：`theil_g2_cross` 组内权重用人口份额（应为收入份额）→ 分解总和不等于真实总量 | 改收入份额 + 恒等式回归测试 |
| `ead97af` | eval_weights | AHP RI 表 n=12 为 1.48（应为 Saaty 1.54）→ CR 错 | 改正 + 回归测试 |
| `26fa4fc` | eval_preprocess | `rescale_extreme` 文档与代码相反 | 修正文档 |

## 未决 / 建议

1. **eval_dea 的 SBM / super-efficiency / Malmquist** 未做独立数值对照（无现成参考在范围内）——若用于正式发表，建议补充验证。
2. **eval_system 的 "adjusted" 耦合协调度**（引 Wang et al. 2021）无法在无原文情况下独立确认精确常数，已用不变量（相等→1、值域）验证；建议对照原文。
3. **na_mcar_test 与 naniar 的 0.9% 数值差** 属 EM/统计量实现细节（包 EM 已用直接似然优化证明收敛到真 MLE），非 bug。

## 结论

自研算法的**核心数值全部正确**：所有能与参考包对照的统计量、评价结果、指数、DEA 效率值均精确或高精度一致。发现并修复 **1 个 Critical（Theil 交叉分解）**、**1 个 Important（AHP RI）**、**7 个 stat_infer 语义/契约问题**、1 处文档反转。纯封装模块（mv_stats/pred_ts/stat_infer 多数/DEA 的 LP）按用户要求未做重复审查。
