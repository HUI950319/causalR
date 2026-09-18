# Beyond the ATE：lmtp 工作坊中译

本目录是 [**Beyond the ATE: Estimating the causal effects of binary, categorical,
continuous, and multivariate exposures in R**](https://www.beyondtheate.com/)
工作坊讲义的逐章中文翻译。

原作者：Nick Williams、Kara Rudolph、Iván Díaz。
原始仓库：<https://github.com/nt-williams/lmtp-workshop>，采用 **GNU GPL-3.0** 许可。

GPL-3.0 明确允许修改与再分发（翻译属于「修改后的版本」），条件是保留许可证、
标明作者并标注改动。本目录据此发布：每页顶部保留来源、作者与许可声明，
并注明「本页为中文翻译，非原文」。

## 翻译约定

- 数学记号、公式和 LaTeX 宏（`\dd`、`\P`、`\E` 等）按原文保留，宏定义见
  [`macros.qmd`](macros.qmd)。
- R 代码、变量名、包名和函数参数一律不译；代码内的注释译为中文。
- 原站的交互式代码块由 [webR](https://docs.r-wasm.org/webr/latest/) 在浏览器中
  执行。本站没有安装 webR 扩展，这些代码块改为普通静态代码块呈现，
  **不在本站重新执行**，因此不附运行结果。需要实际运行请回原站。
- 引用键与参考文献沿用原始 `references.bib`（本目录内为 `lmtp-refs.bib`）。
- 原文插图直接复制到 [`images/`](images/)，未做改动。

## 章节对照

| 本目录 | 原站页面 | 标题 |
|---|---|---|
| `00-welcome.qmd` | `index.html` | Welcome! |
| `00-instructors.qmd` | `00_instructors.html` | Instructors |
| `01-introduction.qmd` | `01_info_introduction.html` | The Causal Model and Notation |
| `02-defining-interventions.qmd` | `02_info_d.html` | Defining Interventions |
| `03-estimators.qmd` | `03_info_estimators.html` | Estimators |
| `04-lmtp-package.qmd` | `04_lmtp.html` | The lmtp package |
| `05-static.qmd` | `05_R_static.html` | Static effects and causal contrasts |
| `06-dtr.qmd` | `06_R_dtr.html` | Dynamic Treatment Regimes |
| `07-mtp.qmd` | `07_R_mtp.html` | Modified Treatment Policies |
| `08-ipsi.qmd` | `08_R_ipsi.html` | IPSI |
| `09-survival.qmd` | `09_R_survival.html` | Survival Analysis |
| `10-multivariate.qmd` | `10_R_multivariate.html` | Multivariate Exposures |
| `11-conclusion.qmd` | `11_info_conclusion.html` | Final Remarks |
