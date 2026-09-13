# causalR

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

观察性研究的因果推断与**目标试验模拟（Target Trial Emulation, TTE）**工具包，
面向临床与流行病学研究者，R 路线。

仓库同时是一份 **24 周学习记录**：`learn/` 是路线与笔记，`weeks/` 是每周的分析代码，
`R/` 里的函数是从这些代码中被反复用到、最终沉淀下来的那部分。

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

> 当前版本 `0.0.0.9000`，尚无导出函数。包的内容随学习进度生长，见下。

---

## 包函数的生长规则

不预先设计 API。函数只在**被真实用到两次以上**时才存在：

| 阶段 | 触发 | 动作 |
|---|---|---|
| 1 | 某段代码只在一周里用到 | 留在 `weeks/wNN/` |
| 2 | 第二周又用到同一段逻辑 | 抽进 `R/`，写 roxygen + 测试 |
| 3 | `R/` 攒到一组成体系的函数 | 打 tag，写 vignette |

预计最先沉淀下来的（约 W14 前后）：

- `ccw_*()` —— clone-censor-weight 的克隆与删失构造（W14 手写模板）
- `diag_weights()` —— 权重诊断：ESS、极值比、SMD
- `plt_love()` / `plt_forest_sens()` / `plt_trt_timeline()` —— W22 的三张图
- `pool_rubin()` —— 多重插补下多估计量结果的 Rubin 合并（W16 管线）

---

## 开发约定

### 依赖

| 包 | 位置 | 说明 |
|---|---|---|
| `UtilsR` | 可用 | https://github.com/HUI950319/UtilsR 公开，可以 `Imports` |
| `RegR` | **待决** | 仓库当前私有。公开的 causalR 不能硬依赖它，否则他人装不上 |

**待决项（推迟到 W22 第一个需要存图的函数出现时再定）**：
`RegR::save_plt()` / `save_tb()` 的复用在公开包里走不通，三个选择——
① 把 RegR 转公开；② `RegR` 放 `Suggests`，缺失时降级为直接 `ggsave()`；
③ causalR 自己实现。在做出选择前，不要在 `R/` 里引入任何 RegR 调用。

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
