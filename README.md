# causalR

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Study notes](https://img.shields.io/badge/study_notes-online-1b4079.svg)](https://hui950319.github.io/causalR/)
[![Code site](https://img.shields.io/badge/code_site-cibookex--r-2b7489.svg)](https://hui950319.github.io/cibookex-r/)
<!-- badges: end -->

观察性研究的因果推断与**目标试验模拟（Target Trial Emulation, TTE）**工具包，
面向临床与流行病学研究者，R 路线。

仓库同时是一份 **24 周学习记录**：`learn/` 是路线与笔记，`weeks/` 是每周的分析代码，
`R/` 里的函数是从这些代码中被反复用到、最终沉淀下来的那部分。

### 📖 在线读书笔记

**<https://hui950319.github.io/causalR/>**

每周一章，结构是「全章地图 → 逐节精读 → 关键陷阱 → 自检问题 → 交付说明」。
站点由 Quarto book 生成（源文件在 [`notes/`](notes/)，产物在 `docs/`），支持全文搜索与明暗主题。

| 已上线 | 主题 |
|---|---|
| [W1 · 第 1 章](https://hui950319.github.io/causalR/notes/w01-potential-outcomes.html) | 潜在结果框架与三大识别假设 |
| [W1 · 第 2 章](https://hui950319.github.io/causalR/notes/w01-randomized-experiments.html) | 随机试验：可交换性、标准化与 IP 加权 |

### 💻 配套代码站点

**<https://hui950319.github.io/cibookex-r/>**

原书**第二部分（第 11–17 章）**的 R 与 Stata 代码，已渲染成书（另含 PDF 与 EPUB）。
fork 自 Tom Palmer 的 [cibookex-r](https://github.com/remlapmot/cibookex-r)，
其中 R 代码原作者为 Joy Shi 与 Sean McGrath，Stata 代码原作者为 Eleanor Murray 与 Roger Logan。
上面读书笔记中第 11–17 章的页顶，都有指向对应代码页的链接。

> **该站点按 GPL-3 发布，与本仓库的 MIT 许可相互独立。**

---

## ⚠️ 数据声明

**本仓库为公开仓库，不包含任何患者数据。**

- MIMIC-IV 受 [PhysioNet DUA](https://physionet.org/content/mimiciv/) 约束，禁止再分发
- 院内 HIS / 自有队列含患者信息，任何形式都不入库
- `data-raw/` 采用「默认拒绝 + 白名单」策略，只有 `data-raw/synthetic/` 下的
  小体量合成数据会被跟踪
- `.githooks/pre-commit` 会拦截误暂存的数据文件

复现分析需要你自己按 [`data-raw/README.md`](data-raw/README.md) 申请对应数据。

---

## 仓库结构

```
causalR/
├── R/              包函数（跨周复用的代码沉淀于此）
├── man/  tests/    文档与测试
│
├── learn/          学习资料
│   ├── 00-roadmap.md      ★ 24 周路线，唯一的计划来源
│   ├── 00-setup.md          环境配置与踩坑记录
│   ├── notes/               每周笔记
│   ├── protocol/          ★ 目标试验协议、TARGET 自查、DAG（会进论文）
│   └── refs/                references.bib（PDF 不入库）
│
├── weeks/          每周代码，w01 … w20
├── milestones/     M1–M4 交付物（Quarto 报告）
├── data-raw/       数据（除 synthetic/ 外均不入库）
├── output/         最终图表（cache/ 不入库）
└── paper/          W21–24 投稿产出
```

`learn/`、`weeks/`、`milestones/`、`data-raw/`、`output/`、`paper/` 都在
`.Rbuildignore` 里，不进入构建产物——`R CMD build` 出来的仍是一个干净的 R 包。

---

## 安装

```r
# install.packages("remotes")
remotes::install_github("HUI950319/causalR")
```

当前版本 `0.0.0.9000`，处于实验阶段，已提供 6 个导出函数：

| 分析 | 计算函数 | 绘图函数 |
|---|---|---|
| 倾向评分加权与平衡诊断 | `get_PSW()` | `plt_PSW()` |
| 倾向评分匹配与匹配样本诊断 | `get_PSM()` | `plt_PSM()` |
| 未测量混杂敏感性分析（线性模型、Cox、DML、IV） | `get_sens()` | `plt_sens()` |

计算函数返回权重或敏感性分析结果，并附诊断与模型对象；绘图函数接收对应结果。参数及示例见
`?get_PSW`、`?get_PSM`、`?get_sens` 和相应绘图函数的帮助页。

---

## 包函数的生长规则

不预先设计 API。函数只在**被真实用到两次以上**时才存在：

| 阶段 | 触发 | 动作 |
|---|---|---|
| 1 | 某段代码只在一周里用到 | 留在 `weeks/wNN/` |
| 2 | 第二周又用到同一段逻辑 | 抽进 `R/`，写 roxygen + 测试 |
| 3 | `R/` 攒到一组成体系的函数 | 打 tag，写 vignette |

后续候选函数（按实际复用情况决定）：

- `ccw_*()` —— clone-censor-weight 的克隆与删失构造（W14 手写模板）
- `diag_weights()` —— 权重诊断：ESS、极值比、SMD
- `plt_love()` / `plt_forest_sens()` / `plt_trt_timeline()` —— W22 的三张图
- `pool_rubin()` —— 多重插补下多估计量结果的 Rubin 合并（W16 管线）

---

## 开发约定

### 依赖

| 包 | 位置 | 说明 |
|---|---|---|
| `WeightIt`、`MatchIt` | `Imports` | 倾向评分建模与匹配 |
| `halfmoon` | `Suggests` | 协变量平衡诊断及倾向评分分布图 |
| `sensemakr`、`survival`、`tipr`、`dml.sensemakr`、`iv.sensemakr` | `Suggests` | 按所选敏感性分析方法使用相应后端 |
| `ggplotify` | `Suggests` | 将部分后端的基础图形转换为 ggplot |
| `RegR` | `Suggests` | 非空 `save` 参数通过 `RegR::save_plt()` 保存 PDF |

可选依赖在调用对应功能时检查。三个绘图函数的 `save = list()` 或 `NULL`
均不保存文件；传入非空保存参数时需要安装 `RegR`。完整依赖见 `DESCRIPTION`。

依赖方向只能向下，不引入反向硬依赖。

### 换行符

`.gitattributes` 已设 `eol=lf`。Windows 下编辑、WSL 下跑测试不会出现整文件假修改。

### 提交前

```bash
git config core.hooksPath .githooks   # 首次克隆后执行一次，启用数据守卫
```

---

## 学习进度

| 阶段 | 周次 | 主题 | 里程碑 | 状态 |
|---|---|---|---|---|
| 〇 | W0 | 环境配置、数据申请 | — | ⬜ |
| 一 | W1–6 | 估计量全景：PO / DAG / PS / g-formula / DR / 敏感性 | M1 | ⬜ |
| 二 | W7–10 | g-methods：时变混杂 / MSM / g-formula / 竞争风险 | M2 | ⬜ |
| 三 | W11–16 | 目标试验模拟：协议 / TARGET / 序贯试验 / CCW / EHR 管线 | M3 | ⬜ |
| 四 | W17–20 | 因果机器学习：SuperLearner / TMLE / lmtp / CATE | M4 | ⬜ |
| 五 | W21–24 | 产出：完整分析 / 可视化 / 按 TARGET 撰写 / 投稿 | — | ⬜ |

细化到周的进度表和自检问题在 [`learn/00-roadmap.md`](learn/00-roadmap.md)。

> 硬性节点：**里程碑达不到标准不要进入下一阶段。**
> 尤其是 M3——time zero 定义错了，第四阶段学再多 ML 也只是把错的问题算得更精确。

---

## 主要参考

- Hernán & Robins, *Causal Inference: What If* —— https://miguelhernan.org/whatifbook
- TARGET 声明（Cashin et al. 2025, JAMA）—— https://www.target-guideline.org
- *Causal Inference in R* —— https://www.r-causal.org/
- TMLE 医学入门（Karim）—— https://ehsanx.github.io/TMLEworkshop/

完整资源索引见 [`learn/00-roadmap.md`](learn/00-roadmap.md#资源总索引)。

---

## License

代码 [MIT](LICENSE.md)。`learn/` 下的笔记为个人学习记录；其中改编自
[causal-inference-visual-guides](https://github.com/kathoffman/causal-inference-visual-guides)（CC-BY）
的图示已在对应位置注明出处。
