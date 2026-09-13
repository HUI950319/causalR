# 第 0 周 · 环境配置

> 本文件是**实际执行记录**，不是计划。装了什么、踩了什么坑，都写在这里。
> 计划见 [00-roadmap.md](00-roadmap.md) 末尾的「第 0 周 · 准备」清单。

## 环境

| 项 | 版本 / 路径 | 备注 |
|---|---|---|
| R |  |  |
| RStudio |  |  |
| Quarto |  |  |
| Git |  |  |
| renv |  | 首次 `renv::init()` 日期： |

## R 包安装记录

分阶段装，不要一次性全装——版本冲突会在错误的时间点爆发。

### 阶段一（W1–6）

```r
install.packages(c(
  "tidyverse", "ggdag", "dagitty",
  "MatchIt", "WeightIt", "cobalt", "halfmoon",
  "marginaleffects", "tipr", "causaldata"
))
```

### 阶段二（W7–10）

```r
install.packages(c("ipw", "gfoRmula", "survival", "adjustedCurves"))
```

### 阶段三（W11–16）

```r
install.packages(c("TrialEmulation", "data.table", "mice"))
```

### 阶段四（W17–20）

```r
install.packages(c("SuperLearner", "tmle", "ltmle", "grf", "policytree", "DoubleML", "mlr3"))
remotes::install_github("nt-williams/lmtp")
```

> tlverse 系列（sl3 / tmle3 / hal9001 / origami）只能从 GitHub 装，常撞 API 限流。
> 先 `usethis::create_github_token()` + `gitcreds::gitcreds_set()`。
> 见 [00-roadmap.md](00-roadmap.md) 的「常见卡点 6」。

## 数据申请

| 数据 | 状态 | 提交日期 | 获批日期 |
|---|---|---|---|
| NHEFS（随 What If 免费） | | | |
| MIMIC-IV / PhysioNet（审批 2–4 周） | | | |
| CITI 培训 | | | |
| 自有队列（伦理 / 脱敏） | | | |

## 踩坑记录

