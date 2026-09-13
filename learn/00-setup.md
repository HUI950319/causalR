# 第 0 周 · 环境配置

> 本文件是**实际执行记录**，不是计划。装了什么、踩了什么坑，都写在这里。
> 计划见 [00-roadmap.md](00-roadmap.md) 末尾的「第 0 周 · 准备」清单。

## 环境

| 项 | 版本 / 路径 | 备注 |
|---|---|---|
| R | 4.4.3 (2025-02-28 ucrt)｜`D:\software\R\R-4.4.3` | 学习周代码在 **Windows R** 跑；包构建 / 测试仍在 WSL R 4.4.3 |
| RStudio | 2026.01.0+392｜`D:\software\R\RStudio2026` | |
| Quarto | CLI 1.8.25｜R 包 `quarto` 1.5.1 | |
| Git | 2.51.1.windows.1 | `core.hooksPath = .githooks` 已启用 |
| renv |  1.1.5 | 首次 `renv::init()` 日期：**未执行** |
| 用户库 | `C:/Users/ouyan/AppData/Local/R/win-library/4.4` | 阶段一的包都装在这里，不进系统库 |
| CRAN 镜像 | `https://mirrors.tuna.tsinghua.edu.cn/CRAN/` | 清华源，`available.packages()` 0.7 s |

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

**实际安装｜2026-09-13｜Windows R 用户库**

| 包 | 版本 | |
|---|---|---|
| `ggdag` | 0.2.13 | 本次安装 |
| `dagitty` | 0.3-4 | 本次安装 |
| `halfmoon` | 0.2.0 | 本次安装 |
| `marginaleffects` | 0.32.0 | 本次安装 |
| `tipr` | 1.0.2 | 本次安装 |
| `causaldata` | 0.1.4 | 本次安装 |
| `tidyverse` `MatchIt` `WeightIt` `cobalt` | — | 此前已装 |

装的时候用了 `dependencies = TRUE`，把 `marginaleffects` 的 Suggests 树
（rstan / brms / lavaan / fixest 等）也一起拉了进来，共解包 125 个二进制包。
下次分阶段安装用默认的 `dependencies = NA` 就够，阶段四的 tlverse 尤其不要用 `TRUE`。

唯一告警：`dependency 'Icens' is not available`——Bioconductor 包，是 `interval`
（`marginaleffects` 的 Suggests）的依赖，不在阶段一需要的范围内，六个目标包
`requireNamespace()` 全部通过。

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
| NHEFS（随 What If 免费） | ✅ 不用下载 | — | — |
| MIMIC-IV / PhysioNet | ✅ 已获取 | | |
| CITI 培训 | | | |
| 自有队列（伦理 / 脱敏） | | | |

NHEFS 不必去 What If 官网下载：`causaldata::nhefs_complete` 直接带（1566 × 67），
W1–W6 全程够用。MIMIC-IV 放 `data-raw/raw/`，该目录整体不入库。

## 踩坑记录

### `.csv.gz` 绕过了两道数据守卫（2026-09-13 已修）

MIMIC-IV 从 PhysioNet 下下来是 `.csv.gz`，而原来的守卫都只认 `.csv`：

- `.gitignore` 的「双保险」段只有 `*.csv`，没有 `*.gz`
- `.githooks/pre-commit` 的正则是 `\.(csv|...)$`，`.csv.gz` 不匹配

`data-raw/` 下因为有 `data-raw/*` 通配所以安全，但 `weeks/`、`output/`、`learn/`
下的 `.csv.gz` 两道都漏——W16 的 EHR 管线恰好就在 `weeks/w16-ehr-pipeline/` 工作。

修法：`.gitignore` 增加 `*.gz` `*.bz2` `*.xz` `*.zst` `*.zip` `*.7z`；
pre-commit 正则加上可选压缩后缀，并单独拦裸压缩包。

用 `git check-ignore` 逐路径验证过，白名单（`data-raw/synthetic/`、`data/`、
`inst/extdata/`、源码、文档）不受影响。
