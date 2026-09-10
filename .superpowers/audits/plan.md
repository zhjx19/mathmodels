# 包级对照审核计划

日期：2026-09-09
授权：用户批准范围（A 类全审 + B 类公式核对，bug 直接修，全程留底）
交付：`docs/superpowers/audits/<模块>-audit.md` + 修复 commit + 本 ledger

## 参考包安装状态

已装：psych, pcor/ppcor, ineq, DescTools, FactoMineR, Benchmarking, markovchain,
imputeTS, coin, truncdist, mice, naniar, forecast, cluster, MASS, car, pracma。
未装：`grey`（CRAN 已下架）→ pred_grey 用 GM(1,1) 闭式解手算核对替代。

## 批次划分（每批 = 1 个对照 subagent + 1 个修复 subagent（如发现）+ 1 个复核）

| 批次 | 模块 | 参考基准 | 函数 |
|---|---|---|---|
| 1 | stat_infer.R | stats / coin / car | stat_describe, stat_t_test, stat_wilcox_test, stat_anova, stat_kruskal, stat_ancova, stat_cor, stat_chisq, stat_normality, stat_hov_test, stat_bootstrap, stat_perm_test |
| 2 | mv_stats.R | stats::princomp/factanal/cancor, psych, cluster, FactoMineR, ppcor | mv_pca, mv_efa, mv_discrim, mv_mds, mv_corresp, mv_cancor, mv_pcor, mv_hclust, mv_kmeans, mv_cluster_metrics |
| 3 | pred_ts.R | forecast, imputeTS | ts_stl, ts_test, ts_transform/back, ts_ets, ts_sarima, ts_arimax, ts_forecast, impute_ts_df |
| 4 | data_cleaning.R | naniar, mice, stats | na_mcar_test(Little), impute(线性/样条 vs pracma/splines), impute_model, outlier |
| 5 | pred_grey + pred_markov + eval_inequality | 闭式手算 / markovchain / ineq+DescTools | grey_pred, gm11 等, markov 预测, gini/theil/atkinson |
| 6 | interp_fit + eval_weights | pracma / splines, 公式 | interp_linear, interp_spline, lm/mono 拟合, 熵权/AHP/变异系数 |
| 7 | eval_preprocess + eval_fuzzy + eval_system + eval_dea | 公式核对 + Benchmarking(DEA) | 标准化/归一化, 隶属函数, AHP, dea |
| 8 | diffEq_models | deSolve 数值解 | SIR/SEIR 等 ODE |
| 9 | 其余（C 类行为审查） | 无 | eval_algorithms/regional, water_quality, export_table, read_nbs, pred_reg |

## 约定

- 对照脚本放 `.superpowers/audits/probe-<模块>.R`，可复跑；发现写 `docs/superpowers/audits/<模块>-audit.md`。
- bug 分级：Critical（错误结果）/ Important（崩溃或错误语义）/ Minor（文档/边角）。
- 修复走 subagent，配回归测试，逐批 commit；全部记录留 ledger。
- C 类（批次 9）只做边界行为审查，不做数值对照。
- 环境：工作区有用户在制 0.0.14 发行文件，staging 一律外科手术式。
Audit batch A (pred_grey + eval_inequality): DONE. pred_grey 31/31 PASS 0 bugs. eval_inequality: BUG theil_g2_cross within-weight used pop-share not income-share -> decomposition total wrong; fixed f00454f + regression tests; 21/21 PASS after.
