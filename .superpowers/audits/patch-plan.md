# mathmodels 0.0.14 补丁实施计划

来源：`D:\WorkBuddy\mathmodels技能\mathmodels-0.0.14-修复补丁建议.md`
已逐条实证（见 `.superpowers/audits/verify-patch-claims.R` / `.log`、`probe-dea-na.R`）。

| # | 文件 | 动作 | 状态 |
|---|---|---|---|
| 1 | `R/eval_fuzzy.R:606-612` | type=5 加 `r0 == 0` 分支；守卫改 `!is.finite()` | ✅ 完成 |
| 2 | `R/eval_inequality.R:412-415` | 末块改名 `.within`；补 `@return` 文档 | ✅ 完成 |
| 2b | `R/eval_inequality.R:401-405` | **新发现**：`.within` 值也错（每组都取全局合计），改为按组计算 | ✅ 完成 |
| 3 | `R/eval_fuzzy.R:37-41` | `mf_dsigmoid` 从 `[0,1]` 摘出单列 | ✅ 完成 |
| 4 | `R/pred_markov.R:59-62` | `dimnames(P)` 命名 | ✅ 完成 |
| 5 | `R/eval_weights.R:348` | 加 warning（限定混尺度场景，`index = NULL` 不告警） | ✅ 完成 |
| 6 | `R/eval_dea.R:318,367` | basic_DEA + super_DEA 收集失败 DMU 并 warning | ✅ 完成 |
| 7 | `R/eval_dea.R:389-393` | basic_SBM 守卫 —— 保持 | 不改 |
| 8 | `R/eval_regional.R` LQ | 契约已文档化 | 不改 |

## 回归测试

- [x] `test-fuzzy_more.R`：type 5 零等级容忍 + 全零报错 + `mf_dsigmoid` 值域锁定 + 负值报错文案
- [x] `test-inequality.R`：`theil_g2_nest` ratio 名唯一 + 每组取自身 + 分量和为 1
- [x] `test-markov.R`：`dimnames(trans_mat)` 有名 + `as_tibble()` 不再报错
- [x] `test-entropy.R`：混尺度触发 warning、归一化/`index = NULL` 保持安静
- [x] `test-DEA.R`：正常路径零 warning + mock 求解器验证失败告警路径

## 收尾

- [x] `roxygen2::roxygenise()`（8.1.0 顺带重排了 NAMESPACE 格式）
- [x] 全量 36 个测试文件：0 failed / 0 error（基线对比确认零新增 warning）
- [x] 更新 `NEWS.md`

## 风格约定

包源码用 `=` 赋值、原生 `|>`（与全局 AGENTS.md 一致，非 `<-`）。
