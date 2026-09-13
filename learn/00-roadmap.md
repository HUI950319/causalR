# 因果推断与目标试验模拟 · 24 周学习路线（R 路线）

> **适用对象**：有 R 基础的医学 / 临床流行病学研究者
> **主线**：潜在结果框架 → g-methods → 目标试验模拟（TTE） → 因果机器学习 → 可投稿产出
> **总投入**：约 24 周，每周 8–14 小时

**关于本文件**
本文件是 wolai 在线笔记（顶级页面 https://www.wolai.com/oNRPztAGPU8ASFLbCs4QrT ）对应内容的本地离线版，
覆盖「阶段一～阶段五 + 里程碑 + 卡点 + 资源索引 + 进度表」。
在线版另有两节未包含在此文件中：`00 · 总览、自测与每周记录模板` 与 `01 · 第 0 周：环境配置`。
第 0 周的待办清单可在本文件末尾的「进度追踪表」中找到简版。

---

## 目录

- [阶段一 · 估计量全景（W1–6）](#阶段一--估计量全景w16)
- [阶段二 · g-methods 与时变混杂（W7–10）](#阶段二--g-methods-与时变混杂w710)
- [阶段三 · 目标试验模拟（W11–16）](#阶段三--目标试验模拟w1116)
- [阶段四 · 因果机器学习（W17–20）](#阶段四--因果机器学习w1720)
- [阶段五 · 产出（W21–24，可选）](#阶段五--产出w2124可选)
- [里程碑验收标准](#里程碑验收标准)
- [常见卡点与应对](#常见卡点与应对)
- [资源总索引](#资源总索引)
- [进度追踪表](#进度追踪表)
- [附录 · 内容来源与准确性说明](#附录--内容来源与准确性说明)

---

## 阶段一 · 估计量全景（W1–6）

> 🎯 **阶段目标**：搞清楚各种估计量之间的关系，能在单时点暴露的场景独立完成一个规范的因果分析。本阶段末交付里程碑 M1。

### 第 1 周 · 潜在结果框架与三大识别假设｜8–10 h

| 环节 | 内容 |
|---|---|
| 阅读 | Hernán & Robins《What If》第 1–3 章，官方免费下载（含数据与 R/Stata/SAS/Python/Julia 代码）：https://miguelhernan.org/whatifbook |
| 跑码 | migariane 仓库的 R 代码，只跑「回归调整」和「G-formula」两节：https://github.com/migariane/TutorialCausalInferenceEstimators |
| 产出 | 手绘「估计量演进图」：回归调整 → G-formula → 倾向评分 → IPW → 双稳健 → TMLE，每个箭头上写清楚「解决了前者的什么问题」 |
| 自检 | ① 为什么个体因果效应不可识别而平均因果效应可以？② 一致性假设在「手术 vs 保守治疗」上具体意味着什么？③ 正性假设违背时 IPW 权重会发生什么？ |

migariane 这个仓库明确面向流行病学家，按历史顺序讲各估计量如何解决前一个的缺陷，并提供 R / Stata / Python 三套带注释代码。

> ⭐ 如果只允许 clone 一个仓库开始，就是它。它把「这些方法到底什么关系」一次性理顺。

### 第 2 周 · DAG、混杂与碰撞偏倚｜10–12 h

| 环节 | 内容 |
|---|---|
| 阅读 | 《What If》第 6–8 章；edX 上 HarvardX 的 Causal Diagrams: Draw Your Assumptions Before Your Conclusions（免费旁听，可一周刷完） |
| 跑码 | ggdag 文档 https://r-causal.github.io/ggdag/ ；dagitty 官网 https://www.dagitty.net/ ；用 ggdag 画出你自己研究问题的 DAG |
| 产出 | 你自己研究问题的 DAG 图（会一直用到第 24 周，也会进论文和基金本子）；用 dagitty::adjustmentSets() 列出最小充分调整集 |
| 自检 | ① DAG 里有没有中介变量被误当成混杂调整了？② 有没有变量是暴露与结局的共同结果（碰撞）？③ 关键混杂未测量时打算怎么处理？ |

> ⭐ **本周是整个计划里性价比最高的一周。** 医学研究里绝大多数「调整了一堆协变量还是有偏」的问题，根源都在这里。

### 第 3 周 · 倾向评分：匹配与加权｜10–12 h

| 环节 | 内容 |
|---|---|
| 阅读 | 《Causal Inference in R》倾向评分章节，在线版 https://www.r-causal.org/ ，源码 https://github.com/r-causal/causal-inference-in-R |
| 跑码 | MatchIt + cobalt 做匹配与平衡诊断（lalonde）；WeightIt 做 IPTW / ATT / overlap 权重（NHEFS）；halfmoon 做平衡可视化 https://github.com/r-causal/halfmoon |
| 产出 | 匹配 vs IPTW 结果对比表；SMD 的 Love plot；权重分布图 + 极端权重的识别与处理记录 |
| 自检 | ① ATE / ATT / ATO 分别是什么人群的效应？你的临床问题需要哪一个？② SMD < 0.1 就算平衡好了吗？还需要看什么？③ 权重最大值到了 80，你会怎么做？截尾的代价是什么？ |

包文档：MatchIt https://kosukeimai.github.io/MatchIt/ ｜WeightIt https://ngreifer.github.io/WeightIt/ ｜cobalt https://ngreifer.github.io/cobalt/

### 第 4 周 · 标准化与 g-formula（单时点）｜10–12 h

| 环节 | 内容 |
|---|---|
| 阅读 | 《What If》第 13 章；Kat Hoffman 视觉指南的 G-Computation.pdf https://github.com/kathoffman/causal-inference-visual-guides |
| 跑码 | 手写 g-computation（不用包）：拟合结局模型 → 全体设为暴露 → 预测 → 取均值 → 对未暴露重复 → 相减；再用 marginaleffects::avg_comparisons() 对比 https://marginaleffects.com/ ；Bootstrap 求置信区间 |
| 产出 | 手写版与包版结果一致的证明；一段能给同事讲明白的「g-computation 在做什么」的口语解释 |
| 自检 | ① g-formula 和「在回归模型里放一个暴露变量看系数」的区别是什么？② 结局模型有交互项时，回归系数还能解释成因果效应吗？③ 为什么 g-computation 的标准误不能直接用回归的标准误？ |

Kat Hoffman 这套图解为 CC-BY 授权，组会汇报、开题、基金本子的方法学示意图可以直接改用（注明出处）。医学评审专家看不懂公式，看得懂这个。作者博客：https://khstats.com/

### 第 5 周 · 双稳健估计入门｜10–12 h

| 环节 | 内容 |
|---|---|
| 阅读 | migariane 仓库的双稳健与 TMLE 章节；视觉指南的 TMLE.pdf |
| 跑码 | AIPW 手写实现；把结局模型故意设错看结果，再把暴露模型故意设错看结果 |
| 产出 | 一张 2×2 表：结局模型对/错 × 暴露模型对/错，四种组合下 g-computation、IPW、AIPW 三种估计量的偏倚对比 |
| 自检 | ① 「两个模型只要有一个对就行」——为什么？直觉上怎么理解？② 两个都错会怎样？③ 双稳健能解决未测量混杂吗？ |

> 💡 **故意把模型设错，是理解「双稳健」最快的方式。** 这一周的价值全在这个实验里。

### 第 6 周 · 敏感性分析与偏倚量化｜10–12 h

| 环节 | 内容 |
|---|---|
| 阅读 | VanderWeele & Ding 2017, Ann Intern Med（E-value 原始文献）doi:10.7326/M16-2607 ；阴性对照结局的方法学综述 |
| 跑码 | tipr 包做未测量混杂的敏感性分析 https://r-causal.github.io/tipr/ ；在自己的 DAG 上设计一个阴性对照结局 |
| 产出 | **里程碑 M1**：NHEFS 数据的完整因果分析报告 |
| 自检 | ① E-value = 1.8 意味着什么？怎么向临床读者解释？② 你的研究里能找到什么阴性对照结局？③ 定量偏倚分析和「讨论部分写一句局限性」的区别在哪？ |

E-value 在线计算器：https://www.evalue-calculator.com/

### 阶段一结束时应该具备的能力

- [ ] 能独立画 DAG 并论证调整集
- [ ] 能解释匹配、加权、标准化三条路线各自估计的参数
- [ ] 能手写 g-computation 和 AIPW，不依赖包
- [ ] 能读懂并批评一篇用倾向评分的临床论文

---

## 阶段二 · g-methods 与时变混杂（W7–10）

> 🎯 **阶段目标**：理解时变混杂为什么让传统方法全部失效，掌握三种 g-methods。本阶段末交付里程碑 M2。

> ⚠️ **这是 TTE 的真正地基。** 跳过这四周直接学 TTE 的人，通常会在「per-protocol 分析里怎么处理治疗转换」这一步彻底卡住，或者错误地只用基线协变量调整时变混杂。

### 第 7 周 · 时变混杂的本质｜10–12 h

| 环节 | 内容 |
|---|---|
| 阅读 | 《What If》第 19–20 章 https://miguelhernan.org/whatifbook ；Daniel et al. 2013, Stat Med, "Methods for dealing with time-dependent confounding" |
| 跑码 | 自己写生成过程模拟一个时变混杂数据集；分别用「不调整 / 基线调整 / 时变调整（普通回归）」三种方法估计，观察偏倚方向 |
| 产出 | 模拟结果表 + 为什么「把时变混杂放进回归」会引入新偏倚的解释 |
| 自检 | ① 什么变量既是前一时点治疗的结果、又是后一时点治疗的原因、还影响结局？举一个你专业里的例子。② 把这种变量放进回归会发生什么？为什么叫「过度调整」？③ 为什么 g-methods 能同时处理它？ |

**临床例子提示**：抗生素治疗中的 CRP、机械通气中的氧合指数、化疗中的中性粒细胞计数——这些都是典型的时变混杂兼中介。

### 第 8 周 · 边际结构模型与 IPTW｜12 h

| 环节 | 内容 |
|---|---|
| 阅读 | 《What If》第 12 章、第 21 章；Robins, Hernán & Brumback 2000, Epidemiology（MSM 原始文献） |
| 跑码 | ipw 包拟合 MSM https://cran.r-project.org/package=ipw ；手写稳定化权重；权重诊断：均值应接近 1、极值分布、随时间的累积 |
| 产出 | MSM 分析全流程脚本 + 权重诊断报告 |
| 自检 | ① 稳定化权重的分子放什么？为什么放它？② 权重均值 = 1.4 说明什么问题？③ MSM 里的「边际」是相对什么而言的？ |

**权重会随时间累乘**，随访时点越多，极端权重出现得越早越猛。这是纵向分析和横断面分析最大的实操差别。

### 第 9 周 · 参数化 g-formula｜12 h

| 环节 | 内容 |
|---|---|
| 阅读 | 《What If》第 21 章；gfoRmula 包的方法学论文 |
| 跑码 | gfoRmula 跑通自带示例 https://github.com/CausalInference/gfoRmula ；对比同一问题上 IPTW-MSM 与 g-formula 的结果 |
| 产出 | 两种方法的结果对比 + 差异来源分析；干预策略的设定代码（static / dynamic 两种） |
| 自检 | ① g-formula 需要拟合几类模型？分别是什么？② 什么是 g-null 悖论？③ 动态治疗策略在 g-formula 里怎么表达？ |

g-formula 相比 IPTW 的优势是**能直接表达复杂的动态策略**（如「CD4 降到某阈值就启动治疗」），代价是需要正确设定更多模型。

### 第 10 周 · 生存结局与竞争风险｜12 h

| 环节 | 内容 |
|---|---|
| 阅读 | 《What If》第 17 章；Young et al. 关于竞争风险下因果估计量的论文 |
| 跑码 | IPCW（逆概率删失加权）；加权 Kaplan-Meier；adjustedCurves 包画调整后生存曲线 https://cran.r-project.org/package=adjustedCurves |
| 产出 | **里程碑 M2**：时变暴露 + 生存结局的完整分析 |
| 自检 | ① 竞争风险下，「总效应」和「直接效应」分别对应什么临床问题？② 把竞争事件当删失处理，估计的是什么？合理吗？③ IPCW 的删失模型该放哪些变量？ |

> 🚨 肿瘤、老年医学、ICU 研究里竞争风险几乎必然存在。**「把死亡当删失」在多数临床问题里估计的是一个没有实际意义的假想量**，这一点要能向审稿人解释清楚。

### 阶段二结束时应该具备的能力

- [ ] 能识别自己数据里哪些变量构成时变混杂
- [ ] 能手写稳定化权重并做权重诊断
- [ ] 能用 g-formula 表达一个动态治疗策略
- [ ] 能解释为什么普通 Cox 模型在时变暴露下会出错

---

## 阶段三 · 目标试验模拟（W11–16）

> 🎯 **阶段目标**：能独立设计和执行一个符合 TARGET 规范的 TTE 研究。本阶段末交付里程碑 M3。

### 第 11 周 · TTE 框架与协议七要素｜12 h

**本周阅读量最大，优先保证。本周不写代码，专心写协议。**

#### 阅读（按顺序）

1. Hernán & Robins 2016, Am J Epidemiol, "Using Big Data to Emulate a Target Trial When a Randomized Trial Is Not Available" — doi:10.1093/aje/kwv254
2. Hernán, Wang & Leaf 2022, JAMA, "Target Trial Emulation: A Framework for Causal Inference From Observational Data" — 简明版，适合给临床同事看
3. Hernán, Dahabreh, Dickerman & Swanson 2025, Ann Intern Med — "why and when is it helpful"，**讲清楚什么时候不该用 TTE**
4. Matthews et al. 2022, BMJ — 实操 tutorial

#### 产出：目标试验协议 v1

七要素逐条填写：

| # | 要素 | 你的填写 |
|---|---|---|
| 1 | 入组标准 Eligibility criteria | |
| 2 | 治疗策略 Treatment strategies | |
| 3 | 分组方式 Assignment procedures | |
| 4 | 随访期 Follow-up period | |
| 5 | 结局 Outcome | |
| 6 | 因果对照 Causal contrast（ITT / per-protocol） | |
| 7 | 分析计划 Analysis plan | |

再做一张**映射表**：目标试验的每一要素 ↔ 你的观察数据里用什么变量实现 ↔ 存在什么妥协。

#### 自检

① 你的 time zero 是什么时刻？这个时刻的定义有没有用到之后才知道的信息？
② 入组标准的评估时点和 time zero 是同一时刻吗？如果不是，会引入什么？
③ 你要估计 ITT 还是 per-protocol？临床上哪个更有意义？

> ⚠️ **TTE 最难的从来不是代码，是 time zero 和入组标准的定义。** 这一周花的时间越多，后面越顺。

### 第 12 周 · TARGET 声明逐条精读｜10 h

#### 阅读

- Cashin, Hansford, Hernán, Swanson, Lee, Jones 等, 2025, JAMA，**TARGET 声明**（21 条清单，分 6 个部分：摘要、前言、方法、结果、讨论、其他信息）
- 官网下载清单与 Explanation & Elaboration 文档：https://www.target-guideline.org
- Hansford, McAuley & Cashin 2025, PLoS Medicine（TARGET 导读）

#### TARGET 的四条核心要求

1. 明确标识研究为对目标试验的观察性模拟
2. 概述因果问题及为何采用目标试验模拟
3. 清楚说明目标试验协议——因果估计量、识别假设、数据分析计划——以及如何映射到观察数据
4. 报告每个因果估计量的估计值、精度，以及评估结果对假设和设计／分析选择敏感性的附加分析结果

> 🚨 **这不是可选项**：PLOS Medicine 已宣布所有依赖 TTE（或声称如此）的投稿必须提交完整的 TARGET 清单，逐条对应到稿件相应章节，且这是在 STROBE 等其他适用指南之外的额外要求。其他期刊正在跟进。

#### 产出

**TARGET 21 条逐条自查表**，对照第 11 周写的协议，标出哪些条目目前无法满足、原因是什么。

#### 自检

① 你的研究属于 TARGET 的适用范围吗？（当前范围是「模拟平行组、个体随机化、且对基线混杂进行调整」的目标试验）
② 哪几条你现在肯定写不出来？缺什么？
③ 如果用了工具变量或时变治疗策略，超出了主声明范围，你打算怎么报告？

### 第 13 周 · 序贯试验路线：TrialEmulation 包｜12–14 h

| 环节 | 内容 |
|---|---|
| 阅读 | Su, Rezvani, Seaman, Starr & Gravestock, arXiv:2402.12083 https://arxiv.org/abs/2402.12083 ；Getting Started vignette https://causal-lda.github.io/TrialEmulation/ |
| 跑码 | TrialEmulation 自带数据（data_censored 等）跑通全流程。GitHub https://github.com/Causal-LDA/TrialEmulation ｜CRAN https://cran.r-project.org/package=TrialEmulation |
| 产出 | 序贯试验分析全流程脚本；ITT 与 per-protocol 结果对比 + 差异解释 |
| 自检 | ① 什么是「序贯试验」？为什么一个人可以进入多个试验？② 扩展后的数据集为什么会变得非常大？病例对照抽样怎么解决？③ per-protocol 分析里治疗转换的人怎么处理？为什么需要删失权重？ |

由剑桥 MRC 生物统计所与 Roche 合作维护。提供：序贯试验的数据准备与扩展、处理治疗转换与依赖性删失的逆概率治疗权重与删失权重计算、时间-事件结局的边际结构模型拟合、针对用户指定目标人群的边际 ITT 与 per-protocol 效应估计与推断。

> 💡 **对国内医院数据特别重要的一点**：它能通过分块处理数据和病例对照抽样，在 R 的内存限制内处理 EHR 级别的大数据集。

### 第 14 周 · Clone-censor-weight 路线｜14 h

| 环节 | 内容 |
|---|---|
| 阅读 | Maringe et al. 2020, Int J Epidemiol, "Reflection on modern methods: trial emulation in the presence of immortal-time bias" — doi:10.1093/ije/dyaa057（**R 和 Stata 代码在补充材料里**）；Zhao, Lyu & Yoshida 2021, IJE |
| 跑码 | 手写 CCW 全流程（**这一周必须手写，不能用包**）；对比不加权的克隆分析 vs 加权后的分析 |
| 产出 | CCW 完整脚本（可复用模板）；克隆前后的样本量、事件数、协变量分布变化表 |
| 自检 | ① 克隆为什么能消除不朽时间偏倚？② 克隆引入的组内相关性怎么处理？（提示：稳健方差 / 聚类）③ 什么时候该用 CCW，什么时候该用序贯试验？ |

#### Maringe 五步法

1. 明确目标试验与入组标准
2. 克隆患者
3. 定义删失时间与生存时间
4. 估计权重以处理设计引入的信息性删失
5. 分析

**核心机制**：在 time zero 把每个符合条件的个体克隆到每一个治疗策略组，当其实际治疗与被分配的策略不再相容时对该克隆进行删失，再用逆概率删失权重校正由此引入的选择偏倚。

CCW 主要覆盖三类场景：宽限期（grace period）、时间相关的静态策略、动态策略，以及三者的组合。

> 🚨 **关于不朽时间偏倚的严重性**：近期基于 EHR 的 HPV 疫苗接种 TTE 研究显示，不朽时间偏倚可以直接**反转估计关联的方向**。这不是「结果偏保守」的问题，是结论完全相反的问题。

### 第 15 周 · TTE + 双稳健：从设计到 ML 估计｜14 h

| 环节 | 内容 |
|---|---|
| 阅读 | Hoffman et al. 2022, JAMA Network Open（激素治疗 COVID-19 死亡率的 TTE）；lmtp 方法学论文（Williams & Díaz） |
| 跑码 | 跑通 https://github.com/kathoffman/steroids-trial-emulation |
| 产出 | 跑通全流程 + 用它的可视化脚本画出你自己数据的治疗时间线图；一段文字：换成你的临床场景，Markov 阶数该设几？依据是什么？ |
| 自检 | ① 什么是修正治疗策略（modified treatment policy）？和静态、动态策略的区别？② SDR 和 TMLE 的区别是什么？为什么两者都比 IPW 和 g-formula 更被推荐？③ SuperLearner 库该放哪些学习器？只放 LASSO 够吗？ |

#### 仓库内容

- `analysis.R` —— 500 人模拟数据的精简版分析
- `report_results.R` —— 结果整理
- `trt_timeline_viz.R` —— 患者治疗时间线图
- `forest_plot_viz.R` —— 森林图

主分析用开源 R 包 lmtp https://github.com/nt-williams/lmtp ，SuperLearner 库通过 sl3 构建（演示代码里除 LASSO 和 mean 外的学习器被注释掉以加速，交叉验证折数从论文的 10 折降到 5 折）。

#### 重点关注它的一个设定

作者用了 **Markov 假设 = 2**，即认为前两个时间窗（每窗 48 小时）的时变混杂因素足以捕捉下一时点机制的混杂——这个决定**基于临床知识而非统计准则**。

> ⭐ 这正是 TTE 里最考验临床医生的地方，也是你相对纯统计背景研究者的优势所在。把这个决策过程想透，比会调包重要得多。

### 第 16 周 · 真实 EHR 管线与自有数据映射｜14 h

| 环节 | 内容 |
|---|---|
| 阅读 | OHDSI / RCT-DUPLICATE 相关的基准研究文献（了解 TTE 在真实世界的验证结果） |
| 跑码 | https://github.com/Zhengxian-Fan/target-trial-emulation 的 demo 生成器跑通 |
| 产出 | **里程碑 M3**：自有数据的 TTE 方案书；自有数据结构整理完成（长格式纵向表、时变协变量、删失指示、结局时间） |
| 自检 | ① 你的数据里基线协变量缺失率多少？缺失机制是 MCAR / MAR / MNAR？② 插补应该在克隆／扩展之前还是之后做？为什么？③ 三种估计量结果不一致时，你怎么判断该信哪个？ |

#### 管线结构（基于 CPRD）

PySpark 抽数 → R 做多重插补 → Python 做因果估计。

- **Step 1**：抽取数据、构建序贯试验、应用入组标准、生成合并队列
- **Step 2**：对基线协变量做 MICE 多重插补
- **Step 3**：在插补数据上跑 PSM / IPTW / TMLE 三种估计量，用 Rubin 规则合并结果

带一个**纯合成队列生成器，不需要真实数据授权就能立刻跑通全流程**；另有基于真实协变量、模拟结局的半合成数据集用于验证估计量。

虽然是 Python 主导，但**「多重插补 + 多估计量对照 + Rubin 合并」这个三件套的架构值得完整照抄到你的 HIS 数据上**——国内医院数据缺失严重，这是必需品。

MIMIC-IV 申请：https://physionet.org/content/mimiciv/ ｜配套抽数代码：https://github.com/MIT-LCP/mimic-code

### 阶段三结束时应该具备的能力

- [ ] 能写出一份经得起审稿的目标试验协议
- [ ] 能说清楚自己的 time zero 并论证它没有使用未来信息
- [ ] 序贯试验和 CCW 两条路线都能手动实现
- [ ] 能按 TARGET 21 条逐项自查

---

## 阶段四 · 因果机器学习（W17–20）

> 🎯 **阶段目标**：理解为什么需要 ML、什么时候 ML 会帮倒忙，掌握 TMLE 和 CATE 的正确用法。本阶段末交付里程碑 M4。

> ⚠️ **进入这一阶段前的一个警告**：机器学习不解决混杂，只解决「混杂调整模型设错了」。未测量混杂用什么算法都救不回来。

### 第 17 周 · SuperLearner 与交叉拟合｜12 h

| 环节 | 内容 |
|---|---|
| 阅读 | TMLE workshop 的 SuperLearner 章节 https://ehsanx.github.io/TMLEworkshop/ （明确说明不要求深厚的因果推断或高等统计基础）；视觉指南的 Superlearner.pdf |
| 跑码 | SuperLearner 包：构建学习器库、看交叉验证风险、看各学习器权重 https://cran.r-project.org/package=SuperLearner ；对比单一 GLM vs 单一随机森林 vs SuperLearner 集成 |
| 产出 | 学习器库设定的决策记录：为什么选这几个、样本量多大时该放几个 |
| 自检 | ① 为什么直接把随机森林的预测值当因果效应是错的？（关键词：正则化偏倚）② 交叉拟合（cross-fitting）解决什么问题？③ 学习器库里放 20 个模型一定比放 5 个好吗？ |

> 💡 **经验法则**：样本量 < 500 时，学习器库放 4–6 个（mean、glm、glmnet、ranger、earth）通常足够；盲目堆学习器会显著增加计算时间而收益递减。

### 第 18 周 · TMLE｜12–14 h

| 环节 | 内容 |
|---|---|
| 阅读 | TMLE workshop 全书 https://ehsanx.github.io/TMLEworkshop/ ；Advanced Epidemiological Methods 的 ML in causal inference 章节 https://ehsanx.github.io/EpiMethods/ |
| 跑码 | RHC 数据上跑 TMLE（二分类结局 + 连续结局各一次）https://cran.r-project.org/package=tmle ；与前面学过的 g-computation、IPW、AIPW 结果对比 |
| 产出 | 四方法对照表（估计值 + 95% CI + 计算耗时） |
| 自检 | ① TMLE 的 targeting 这一步在做什么？为什么需要它？② TMLE 和 AIPW 都是双稳健，区别在哪？③ 有效样本量（ESS）掉到原样本的 30%，说明什么？ |

Karim 这套教程覆盖 TMLE 全流程：从构建初始结局模型与暴露模型，到通过倾向评分做定向调整，再到治疗效应估计；强调为暴露模型和结局模型设定多样化的 SuperLearner 库、确定有效样本量、以及候选学习器的选择；用 tmle 包演示应用，并与默认 SuperLearner 库和传统回归做完整对比。连续结局的 TMLE（含变量转换与结果解释）单独有一章。

**进阶补充**：Mondol & Karim 2024, Am J Epidemiol，双交叉拟合 TMLE（DC-TMLE）实操指南，用公开临床数据集演示，解决 TMLE 结合复杂算法时可能出现的偏倚与置信区间覆盖不足问题。

### 第 19 周 · 纵向 TMLE 与 lmtp｜14 h

| 环节 | 内容 |
|---|---|
| 阅读 | lmtp 包文档 https://github.com/nt-williams/lmtp + 方法学论文 |
| 跑码 | lmtp 自带示例（点暴露生存、纵向）；回到第 15 周的 steroids 仓库，这次逐行读懂而不只是跑通 |
| 产出 | 纵向 TMLE 脚本模板；与第 8–9 周的 IPTW-MSM、g-formula 结果对比 |
| 自检 | ① ltmle 和 lmtp 的适用场景差别？② 什么是「修正治疗策略」？它比传统的确定性干预好在哪？③ 时点数增加时，估计的稳定性会怎样？ |

lmtp 实现四种估计器：TMLE（`lmtp_tmle`）、序贯双稳健（`lmtp_sdr`）、参数化 g-formula、IPW。**作者强烈建议只用前两个**，因为它们的理论性质允许在使用机器学习的同时保持有效置信区间和 p 值的计算能力。

ltmle 包：https://cran.r-project.org/package=ltmle

### 第 20 周 · CATE、异质性与 DML 概览｜12–14 h

| 环节 | 内容 |
|---|---|
| 阅读 | Athey、Spiess、Wager（Stanford）的 Machine Learning and Causal Inference（ECON 293 / MGTECON 634），19 讲 × 30 min，覆盖 ATE 与 CATE。课程材料与软件教程：https://web.stanford.edu/~swager/teaching.html ；grf 文档 https://grf-labs.github.io/grf/ |
| 跑码 | grf::causal_forest() 在 NHEFS 或自己的数据上 https://github.com/grf-labs/grf ；RATE 检验；policytree 做治疗规则学习（选做）；DoubleML for R 跑一个示例 https://github.com/DoubleML/doubleml-for-r |
| 产出 | **里程碑 M4**：多估计量对照报告 |
| 自检 | ① 因果森林的变量重要性能用来「筛选亚组」吗？为什么不能？② RATE 检验在检验什么？如果不显著说明什么？③ CATE 估计不稳定的主要来源是什么？怎么验证？ |

DoubleML 实现 Chernozhukov 等 2018 的框架，三要素是 Neyman 正交性、高质量的机器学习估计、样本分割；nuisance 部分的估计基于 mlr3 生态，可用其中任何学习器。文档：https://docs.doubleml.org/

> 🚨 **CATE 是这个领域最容易被误用的工具。** 医学论文里大量「机器学习识别出获益亚组」的结论经不起独立数据验证。用它可以，但必须配 RATE 检验或外部验证，且在讨论部分明确说明是探索性的。

### 三条技术线的优先级（医学研究者）

| 优先级 | 技术线 | 主要包 | 适用 |
|---|---|---|---|
| 1 | TMLE / Targeted Learning | tmle、ltmle、lmtp、SuperLearner | 临床流行病学主流，与 TTE 天然衔接 |
| 2 | CATE / 异质性 | grf、policytree | 精准医学、亚组探索 |
| 3 | 双重机器学习 DML | DoubleML、mlr3 | 概念必懂，医学场景实际用得少 |

**中介分析补充**：CMAverse，因果中介 + 交互分解，医学论文用得很多 https://github.com/BS1125/CMAverse

### 阶段四结束时应该具备的能力

- [ ] 能解释为什么 ML 预测好 ≠ 因果估计准
- [ ] 能独立配置 SuperLearner 库并说明理由
- [ ] 能在同一问题上跑通四种估计量并解释差异
- [ ] 能识别并批评论文中 CATE 的误用

---

## 阶段五 · 产出（W21–24，可选）

> 🎯 **阶段目标**：把前 20 周的学习转成可投稿的成果。这一阶段是可选的，但强烈建议做——不产出的学习会很快遗忘。

### 第 21 周 · 完整分析执行

在自有数据上执行第 16 周的方案书，完成主分析 + 全部敏感性分析。

#### 执行顺序

1. 数据锁定（此后不再改动原始数据，任何修改走脚本）
2. 入组与排除流程图（CONSORT 风格，记录每一步剔除人数和原因）
3. 基线特征表
4. 权重估计与诊断
5. 主分析（ITT）
6. 主分析（per-protocol）
7. 敏感性分析全套
8. 阴性对照分析

#### 分析前必须先注册

即使不是 RCT，也建议在 OSF https://osf.io/ 或 ClinicalTrials.gov https://clinicaltrials.gov/ 预注册分析计划。

> ⭐ **TTE 研究的可信度很大程度上取决于「你是不是先定协议再看结果」**，预注册是最有力的证明。

### 第 22 周 · 结果整理与可视化

| 图表 | 工具 | 说明 |
|---|---|---|
| 基线特征表 | gtsummary https://www.danieldsjoberg.com/gtsummary/ | 按治疗策略分组，报 SMD 而非 p 值 |
| 入组流程图 | DiagrammeR 或手绘 | 每一步剔除人数与原因 |
| 权重诊断图 | halfmoon / cobalt | 权重分布 + 极值标注 |
| Love plot | cobalt | 加权前后 SMD 对比 |
| 调整后生存曲线 | adjustedCurves | 配风险人数表 |
| 森林图 | 参考 steroids 仓库的 forest_plot_viz.R | 主分析 + 各敏感性分析并列 |
| 治疗时间线图 | 参考同仓库 trt_timeline_viz.R | 展示时变暴露模式，审稿人很吃这套 |

参考仓库：https://github.com/kathoffman/steroids-trial-emulation

> 💡 把森林图做成「主分析 + 所有敏感性分析」并列的形式。这比在正文里写一段「结果稳健」有说服力得多。

### 第 23 周 · 按 TARGET 清单撰写

#### 写作顺序

**方法部分先写，前言最后写。** 方法写完你才知道这篇文章真正的贡献是什么，前言才能写准。

#### 逐条对照

边写边填 TARGET 21 条清单，每条标注对应的稿件位置（章节 + 段落）。投稿时这份清单要作为补充材料提交。清单下载：https://www.target-guideline.org

#### 几个容易被挑刺的点

- **time zero 的定义必须在方法部分第一段就说清楚**，不要埋在中间
- 识别假设要逐条论证，不能只写一句「假设无未测量混杂」
- 每一处「协议 ↔ 数据」的妥协都要明写，藏着反而更容易被发现
- per-protocol 分析的删失规则要具体到变量和阈值

### 第 24 周 · 内审与投稿准备

#### 内审

请一位统计背景同事按你的 TARGET 自查表逐条挑刺。**让他专门找「说不清楚的地方」，而不是找错别字。**

#### 补充材料清单

- [ ] 完整分析代码（建议放 GitHub 或 OSF，给 DOI）
- [ ] 目标试验协议全文
- [ ] 协议 ↔ 观察数据映射表
- [ ] TARGET 21 条清单（逐条标注位置）
- [ ] 敏感性分析完整结果
- [ ] 预注册链接

#### 期刊选择提示

PLOS Medicine 已强制要求 TARGET 清单，其他期刊正在跟进。**投 TTE 研究时优先选已明确采纳 TARGET 的期刊**——审稿人懂这套方法，沟通成本低得多。

### 延伸产出

| 产出 | 来源 |
|---|---|
| 基金申请书方法学章节 | 第 11 周的协议 + 第 2 周的 DAG |
| 科室方法学讲座 | 第 3 周的视觉指南 + 第 5 周的 2×2 对照实验 |
| 院内可复用分析模板 | 第 14 周的 CCW 脚本 + 第 16 周的管线架构 |
| 研究生教学材料 | 全套周记录 + 自检问题 |

---

## 里程碑验收标准

> ⚠️ 四个里程碑是这个计划的硬性节点。**达不到标准不要进入下一阶段**——后面的内容会建立在你没掌握的基础上。

### M1（第 6 周末）· NHEFS 完整因果分析报告

#### 必须包含

- [ ] DAG 图 + 最小充分调整集的说明
- [ ] 倾向评分模型 + 平衡诊断（SMD 表 + Love plot）
- [ ] 权重分布图 + 极端权重处理说明
- [ ] 至少三种估计量的结果（g-computation / IPW / AIPW）
- [ ] Bootstrap 或稳健方差的置信区间
- [ ] E-value 敏感性分析
- [ ] 一段面向临床读者的结果解释（不含统计术语）

#### 不合格信号

只报了点估计和 p 值，没有诊断图。

### M2（第 10 周末）· 时变暴露 + 生存结局分析

#### 必须包含

- [ ] 时变混杂的识别与论证（哪个变量、为什么）
- [ ] 稳定化权重的计算 + 权重随时间的诊断
- [ ] MSM 结果
- [ ] g-formula 结果
- [ ] 两种方法差异的解释
- [ ] 调整后生存曲线
- [ ] 竞争风险的处理说明（如适用）

#### 不合格信号

用基线协变量调整了时变混杂还不自知。

### M3（第 16 周末）· 自有数据 TTE 方案书

#### 必须包含

- [ ] 目标试验协议七要素完整填写
- [ ] 协议 ↔ 观察数据的映射表（含每一处妥协的说明）
- [ ] time zero 的定义 + 论证它没有使用未来信息
- [ ] 识别假设逐条论证（可交换性 / 正性 / 一致性），每条说明在你的数据里是否可信
- [ ] 估计量选择的理由（序贯试验 vs CCW）
- [ ] 缺失数据处理计划
- [ ] 敏感性分析清单
- [ ] TARGET 21 条自查表（标出当前无法满足的条目）

#### 不合格信号

说不清 time zero 到底是哪一刻。

> 🚨 这是四个里程碑里最重要的一个。**M3 通不过，第四阶段学再多 ML 也没用**——因为你会把一个定义错的因果问题算得越来越精确。

### M4（第 20 周末）· 多估计量对照报告

同一个问题，用以下方法各算一遍并解释差异：

- [ ] IPTW-MSM
- [ ] 参数化 g-formula
- [ ] TMLE（或 lmtp 的 SDR）
- [ ] （可选）DML
- [ ] 结果差异的来源分析：是模型设定差异、正性违背、还是样本量不足？
- [ ] 计算成本对比（对实际项目排期有意义）

#### 为什么这一步重要

**这一步的收获比读十篇方法学论文都大。** 当四种方法给出接近的结果时，你才有底气说结论稳健；当它们差异很大时，差异本身就告诉你数据的哪个假设出了问题。

### 自评量表

每个里程碑完成后，用下面三个问题给自己打分：

| 问题 | 打分 |
|---|---|
| 我能不看代码，向同事口头讲清楚这个分析在做什么吗？ | ⬜ 能 ⬜ 勉强 ⬜ 不能 |
| 如果审稿人问「为什么不用另一种方法」，我能答上来吗？ | ⬜ 能 ⬜ 勉强 ⬜ 不能 |
| 换一个数据集，我能独立重做一遍吗？ | ⬜ 能 ⬜ 勉强 ⬜ 不能 |

**出现两个「不能」，就回去重做这个阶段。**

---

## 常见卡点与应对

### 速查表

| # | 卡点 | 症状 | 应对 |
|---|---|---|---|
| 1 | 跳过 g-methods | 到第 13 周发现完全看不懂 per-protocol 分析 | 立刻回到第 7 周，不要硬撑 |
| 2 | time zero 定义错 | 结果好得不真实（HR 0.3 之类） | 检查是否用了 time zero 之后才知道的信息 |
| 3 | 正性假设违背 | 权重极值 > 100，有效样本量暴跌 | 看权重分布；考虑改变目标人群（ATO 而非 ATE）；截尾必须报告 |
| 4 | 不做交叉拟合 | 置信区间过窄，覆盖率不足 | ML 估计器必须交叉拟合 |
| 5 | 把因果森林当筛变量工具 | 「机器学习发现 XX 亚组获益」 | 预测重要性 ≠ 因果重要性，必须配 RATE 或外部验证 |
| 6 | tlverse 装不上 | HTTP error 403 API rate limit | 设 GitHub PAT；或用 Binder 云端环境；正式项目用 CRAN 包 |
| 7 | 数据量跑不动 | 序贯试验扩展后内存爆 | TrialEmulation 的分块处理 + 病例对照抽样；或改用 data.table |
| 8 | 只做主分析 | 审稿人第一条意见就是「缺敏感性分析」 | 阴性对照结局 + E-value + 关键设计选择的替代方案，三样起步 |

### 1 · 跳过 g-methods

最常见也最致命。症状是到第 13 周学序贯试验时，完全无法理解「为什么需要删失权重」「per-protocol 和 ITT 为什么要分开算」。

**根因**：没有建立「时变混杂 → 传统方法失效 → g-methods」这条逻辑链。

**应对**：回到第 7 周，重做那个模拟实验（不调整 / 基线调整 / 时变调整三种方法的偏倚对比）。这个实验做透了，后面自然通。

### 2 · time zero 定义错误

**症状**：效应量大得不可思议。手术组 vs 非手术组的 HR 做到 0.3，多半不是手术神奇，是不朽时间偏倚。

**自查三问**：

- time zero 的定义有没有用到之后才知道的信息？（如「接受了手术的患者」——手术发生在 time zero 之后）
- 入组标准的评估时点和 time zero 是不是同一时刻？
- 分组依据是在 time zero 就能确定的吗？

> 🚨 近期 EHR 研究显示，不朽时间偏倚可以直接**反转估计关联的方向**。这不是「保守一点」的问题。

### 3 · 正性假设违背

EHR 数据里极其常见——某些协变量组合下几乎所有人都接受了同一种治疗。

**诊断**：

- 权重最大值 / 中位数的比值
- 有效样本量（ESS）相对原样本的比例
- 倾向评分分布在两组间的重叠区域

**应对优先级**：

1. 改变目标人群（用 ATO / overlap 权重，只估计有重叠区域的效应）
2. 收紧入组标准，把结构性无重叠的人排除掉
3. 截尾（必须在方法部分报告截尾阈值和截尾比例）

**不能做的**：偷偷截尾不报告。

### 4 · 不做交叉拟合

用机器学习估计 nuisance 模型时，如果不做样本分割，会因为过拟合导致置信区间过窄、覆盖率不足。

tmle、lmtp 默认支持；手动用 SuperLearner 时要自己确保训练和估计用的是不同折的数据。

### 5 · 因果森林的误用

**最常见的错误说法**：「因果森林的变量重要性显示 A、B、C 是效应修饰因子。」

**为什么错**：变量重要性反映的是这个变量对 CATE 预测的贡献，在高维、弱信号、样本量有限时极不稳定，换个随机种子结果就变。

**正确做法**：

- 用 RATE 检验整体是否存在可检测的异质性
- 如果 RATE 不显著，就不要报告亚组结论
- 有条件的话做外部数据验证
- 无论如何，讨论部分要写明这是探索性的

### 6 · tlverse 装不上

tlverse 系列（sl3、tmle3、hal9001、origami）至今未上 CRAN，只能从 GitHub 装，常撞 API 限流。

```r
usethis::create_github_token()   # 浏览器生成 token
gitcreds::gitcreds_set()         # 粘贴 token
```

**替代方案**：tlverse 的 workshop 仓库虽已归档，但都配了 Binder，点一下就有预装好全部材料的云端 RStudio，不用在本地折腾依赖。组织主页：https://github.com/tlverse

**生产环境建议**：正式分析用 CRAN 上的 tmle / ltmle / lmtp，这几个包维护活跃。

### 7 · 数据量跑不动

序贯试验的数据扩展会让数据集膨胀数十倍（每个人进入多个试验）。

**应对**：

- TrialEmulation 内置的分块处理和病例对照抽样
- 全程用 data.table 而非 dplyr
- 权重估计和结局模型分开跑，中间结果落盘

### 8 · 只做主分析

**审稿人必问的三件事**：

1. 未测量混杂怎么办？（→ E-value 或定量偏倚分析）
2. 怎么知道没有系统性偏倚？（→ 阴性对照结局）
3. 换个分析选择结果还一样吗？（→ 关键设计选择的替代方案）

这三样是起步配置，不是加分项。

### 遇到没列在这里的问题

先问自己：**这是代码问题，还是因果问题？**

- 代码问题 → 查包文档、GitHub issues
- 因果问题 → 回到 DAG 和目标试验协议

> 💡 **八成的「代码跑不出来」其实是因果问题没想清楚。**

---

## 资源总索引

### 核心仓库

| 仓库 | 用途 | 阶段 | 地址 |
|---|---|---|---|
| TutorialCausalInferenceEstimators | 估计量演进全景，面向流行病学家，R/Stata/Python 三语言 | 一 | https://github.com/migariane/TutorialCausalInferenceEstimators |
| causal-inference-visual-guides | 图解 g-computation / IPW / TMLE / SuperLearner，CC-BY 可改用 | 一、四 | https://github.com/kathoffman/causal-inference-visual-guides |
| causal_inference_notebook | 《What If》第 2 部分 R 代码 | 一 | https://github.com/malcolmbarrett/causal_inference_notebook |
| cibookex-r | 《What If》全书习题的 R 与 Stata 代码，已渲染成书 | 一、二 | https://remlapmot.github.io/cibookex-r/ |
| causal-inference-in-R | tidyverse 风格因果推断教材（Quarto 源码全开） | 一 | https://github.com/r-causal/causal-inference-in-R |
| gfoRmula | 参数化 g-formula | 二 | https://github.com/CausalInference/gfoRmula |
| TrialEmulation | 序贯试验 TTE 工业级实现（剑桥 MRC + Roche） | 三 | https://github.com/Causal-LDA/TrialEmulation |
| steroids-trial-emulation | TTE + 双稳健的最佳过渡样本 | 三 | https://github.com/kathoffman/steroids-trial-emulation |
| target-trial-emulation | 真实 EHR 完整管线（含合成数据生成器） | 三 | https://github.com/Zhengxian-Fan/target-trial-emulation |
| TMLEworkshop | 医学背景最友好的 TMLE 入口 | 四 | https://github.com/ehsanx/TMLEworkshop |
| EpiMethods | Advanced Epidemiological Methods 全书 | 四 | https://github.com/ehsanx/EpiMethods |
| lmtp | 纵向修正治疗策略，TMLE + SDR | 四 | https://github.com/nt-williams/lmtp |
| grf | 因果森林、CATE | 四 | https://github.com/grf-labs/grf |
| doubleml-for-r | 双重机器学习 | 四 | https://github.com/DoubleML/doubleml-for-r |
| tlverse-handbook | Targeted Learning 系统教材 | 四 | https://github.com/tlverse/tlverse-handbook |
| CMAverse | 因果中介 + 交互分解 | 选修 | https://github.com/BS1125/CMAverse |
| awesome-causal-inference | 导航（偏计量，按需取用） | 全程 | https://github.com/matteocourthoud/awesome-causal-inference |

### 浏览入口

比 awesome list 更新更快的方式——直接刷 GitHub topic 页面：

- https://github.com/topics/tmle
- https://github.com/topics/g-computation
- https://github.com/topics/inverse-probability-weights
- https://github.com/topics/targeted-learning

这些页面下面全是真实的分析代码仓库，包括 MSM、纵向数据、生存分析的各种实现。

### 官方文档站

| 资源 | 地址 |
|---|---|
| Causal Inference: What If（含数据与多语言代码） | https://miguelhernan.org/whatifbook |
| TARGET 指南官网 | https://www.target-guideline.org |
| TrialEmulation 文档站 | https://causal-lda.github.io/TrialEmulation/ |
| Causal Inference in R | https://www.r-causal.org/ |
| Targeted Learning 手册 | https://tlverse.org/tlverse-handbook/ |
| Advanced Epidemiological Methods | https://ehsanx.github.io/EpiMethods/ |
| TMLE 医学入门 | https://ehsanx.github.io/TMLEworkshop/ |
| grf 文档 | https://grf-labs.github.io/grf/ |
| DoubleML 文档 | https://docs.doubleml.org/ |
| marginaleffects | https://marginaleffects.com/ |
| MatchIt | https://kosukeimai.github.io/MatchIt/ |
| WeightIt | https://ngreifer.github.io/WeightIt/ |
| cobalt | https://ngreifer.github.io/cobalt/ |
| ggdag | https://r-causal.github.io/ggdag/ |
| dagitty | https://www.dagitty.net/ |
| tipr | https://r-causal.github.io/tipr/ |
| E-value 在线计算器 | https://www.evalue-calculator.com/ |
| Kat Hoffman 博客 | https://khstats.com/ |

### 视频课程（只推荐两门）

| 课程 | 作者 | 体量 | 地址 |
|---|---|---|---|
| Machine Learning and Causal Inference（ECON 293 / MGTECON 634） | Athey、Spiess、Wager（Stanford, 2022） | 19 × 30 min | https://web.stanford.edu/~swager/teaching.html |
| Introduction to Causal Inference | Petersen（Berkeley）、Balzer（UMass） | 完整课件 | https://ctml.berkeley.edu/introduction-causal-inference |

Berkeley 这门课引入 Causal Roadmap 七步框架：明确研究问题 → 定义因果模型与目标效应 → 界定观察数据 → 评估可识别性 → 设定统计估计问题 → 选择并实施估计量 → 恰当解释结果；统计方法涵盖 G-computation、IPW 与 TMLE。曾获美国统计学会 2014 年 Causality in Statistics Education Award。

**补充课程**：

- Hernán 的 Causal Diagrams（HarvardX，edX 平台搜索 "Causal Diagrams: Draw Your Assumptions Before Your Conclusions"）——第 2 周用，医学方向的正统入口
- Kosuke Imai（Harvard, 2022）11 × 50 min，覆盖 PO、ATE、IV、RD、匹配、IPW、固定效应、DiD、CATE，见 awesome list 的 courses.md
- Brady Neal, Introduction to Causal Inference https://www.bradyneal.com/causal-inference-course
- Schuler & van der Laan, Introduction to Modern Causal Inference https://alejandroschuler.github.io/mci/
- Peng Ding, A First Course in Causal Inference https://arxiv.org/abs/2305.18793
- Stefan Wager, STATS 361 讲义与草稿教材 https://web.stanford.edu/~swager/teaching.html

课程总索引（awesome list 的 courses.md）：https://github.com/matteocourthoud/awesome-causal-inference/blob/main/src/courses.md

### 必读文献清单（按阅读顺序）

#### 基础理论

1. Hernán & Robins,《Causal Inference: What If》https://miguelhernan.org/whatifbook
2. VanderWeele & Ding 2017, Ann Intern Med（E-value）doi:10.7326/M16-2607

#### g-methods

3. Robins, Hernán & Brumback 2000, Epidemiology（MSM 原始文献）
4. Daniel et al. 2013, Stat Med（时变混杂方法综述）

#### 目标试验模拟

5. Hernán & Robins 2016, Am J Epidemiol（框架原文）doi:10.1093/aje/kwv254
6. Hernán, Wang & Leaf 2022, JAMA（简明版）
7. Hernán, Dahabreh, Dickerman & Swanson 2025, Ann Intern Med（何时该用／不该用）
8. Matthews et al. 2022, BMJ（实操 tutorial）
9. Maringe et al. 2020, Int J Epidemiol（CCW，附 R/Stata 代码）doi:10.1093/ije/dyaa057
10. **Cashin, Hansford, Hernán 等 2025, JAMA（TARGET 声明）** ← 必须精读
11. Hansford, McAuley & Cashin 2025, PLoS Medicine（TARGET 导读）
12. Hoffman et al. 2022, JAMA Netw Open（激素 TTE 实例）

#### 因果机器学习

13. Su, Rezvani, Seaman, Starr & Gravestock, arXiv:2402.12083（TrialEmulation 包）https://arxiv.org/abs/2402.12083
14. Williams & Díaz（lmtp 包）
15. Chernozhukov et al. 2018, Econom J（DML）
16. Mondol & Karim 2024, Am J Epidemiol（双交叉拟合 TMLE 的 R 实操指南）

### 数据源

| 数据 | 地址 |
|---|---|
| NHEFS（What If 配套） | https://miguelhernan.org/whatifbook |
| MIMIC-IV | https://physionet.org/content/mimiciv/ |
| MIMIC 抽数代码 | https://github.com/MIT-LCP/mimic-code |
| causaldata R 包 | https://cran.r-project.org/package=causaldata |

### 中文资源现状

基本空白。唯一像样的是《Causal Inference for the Brave and True》中文翻译版 https://github.com/xieliaing/CausalInferenceIntro ，但**全部代码基于 Python，面向计量经济学、量化社会学与策略评估**，与医学 TTE 不是同一条技术路线，只能当中文概念扫盲。

Datawhale 社区有人提过做因果推断中文教程，但那是基于图模型（Pearl 路线）的方向，跟 Hernán 的潜在结果路线又是另一支。

> ⚠️ **这条路上没有中文捷径。** 好在 Hernán 的书和 r-causal 的书写得都不难。

### GitHub 资源生态的一个判断

这个领域的 GitHub 资源是「计量经济学重、临床流行病学轻」的。热门的 awesome list 基本都是 DiD、IV、断点回归那一套，对医学 TTE 帮助有限。

**目前没有一个像样的 awesome-target-trial-emulation 列表**——TTE 在 GitHub 上是散的，得自己拼。上面这份清单就是拼出来的结果。

---

## 进度追踪表

> 💡 每周勾一次。**连续两周没勾，就该重新评估时间投入了。**

### 第 0 周 · 准备

- [ ] R / RStudio / Quarto / Git 环境就绪
- [ ] 第一批 R 包安装完成
- [ ] 仓库全部 clone 完成
- [ ] MIMIC-IV / PhysioNet 申请已提交（审批 2–4 周，务必第 0 周提交）
- [ ] CITI 培训完成
- [ ] Zotero collection 建好（四个子文件夹）
- [ ] 自有数据的脱敏与结构整理已排期

### 阶段一 · 估计量全景（W1–6）

- [ ] W1 潜在结果与识别假设 → 估计量演进图
- [ ] W2 DAG 与混杂 → 自己研究问题的 DAG
- [ ] W3 倾向评分 → 匹配 vs IPTW 对比
- [ ] W4 g-formula → 手写版与包版结果一致
- [ ] W5 双稳健 → 2×2 模型设错对照表
- [ ] W6 敏感性分析 → M1 交付

### 阶段二 · g-methods（W7–10）

- [ ] W7 时变混杂本质 → 三方法模拟对比
- [ ] W8 MSM 与 IPTW → 权重诊断报告
- [ ] W9 参数化 g-formula → 两法结果对比
- [ ] W10 生存与竞争风险 → M2 交付

### 阶段三 · 目标试验模拟（W11–16）

- [ ] W11 协议七要素 → 自己问题的目标试验协议 v1
- [ ] W12 TARGET 精读 → 21 条自查表
- [ ] W13 TrialEmulation → 序贯试验脚本
- [ ] W14 Clone-censor-weight → 手写 CCW 模板
- [ ] W15 TTE + 双稳健 → steroids 仓库跑通
- [ ] W16 真实 EHR 管线 → M3 交付

### 阶段四 · 因果机器学习（W17–20）

- [ ] W17 SuperLearner 与交叉拟合
- [ ] W18 TMLE → 四方法对照表
- [ ] W19 纵向 TMLE / lmtp
- [ ] W20 CATE 与 DML → M4 交付

### 阶段五 · 产出（W21–24，可选）

- [ ] W21 自有数据完整分析
- [ ] W22 结果整理与可视化
- [ ] W23 按 TARGET 撰写
- [ ] W24 内审与投稿准备

### 投入记录

| 周次 | 日期 | 计划时长 | 实际时长 | 完成度 | 备注 |
|---|---|---|---|---|---|
| W0 | | 8 h | | | |
| W1 | | 9 h | | | |
| W2 | | 11 h | | | |
| W3 | | 11 h | | | |
| W4 | | 11 h | | | |
| W5 | | 11 h | | | |
| W6 | | 11 h | | | |
| W7 | | 11 h | | | |
| W8 | | 12 h | | | |
| W9 | | 12 h | | | |
| W10 | | 12 h | | | |
| W11 | | 12 h | | | |
| W12 | | 10 h | | | |
| W13 | | 13 h | | | |
| W14 | | 14 h | | | |
| W15 | | 14 h | | | |
| W16 | | 14 h | | | |
| W17 | | 12 h | | | |
| W18 | | 13 h | | | |
| W19 | | 14 h | | | |
| W20 | | 13 h | | | |

### 里程碑状态

| 里程碑 | 目标周 | 实际完成 | 自评 |
|---|---|---|---|
| M1 · NHEFS 完整分析 | W6 | | ⬜ 通过 ⬜ 需返工 |
| M2 · 时变暴露 + 生存 | W10 | | ⬜ 通过 ⬜ 需返工 |
| M3 · 自有数据 TTE 方案书 | W16 | | ⬜ 通过 ⬜ 需返工 |
| M4 · 多估计量对照 | W20 | | ⬜ 通过 ⬜ 需返工 |

### 进度落后时的处理

| 落后程度 | 处理 |
|---|---|
| 1 周内 | 正常，下周补 |
| 2–3 周 | 把后续每周时长下调，拉长总周期，**不要压缩内容** |
| 超过 4 周 | 重新评估：是时间不够，还是某一周卡住了？如果是卡住，回到自检问题定位缺口 |

> ⚠️ **唯一不能做的事**：为了赶进度跳过阶段一到阶段三的任何一周。

---

## 附录 · 内容来源与准确性说明

- 所有 GitHub 仓库地址、包文档站地址均已核实可访问
- 论文给出 DOI 的三条已核实；其余仅给出题名、作者、期刊，请通过 PubMed 或 Zotero 按题名检索获取准确 DOI，**不要凭记忆补全 DOI**
- edX 课程未给直链，因平台 URL 结构变动频繁，按课程名搜索更可靠
- R 包生态更新较快，若某个包 API 变化，以其官方文档站为准
