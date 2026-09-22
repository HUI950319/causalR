# 异质性治疗效应（HTE）模拟：详细中文注释版 ----
# 原作者：Toshiaki Komura；原文日期：2025-05-27。
# 来源：https://github.com/Toshi934/HTE_simulation/blob/main/Simulation.Rmd
# 固定版本：https://github.com/Toshi934/HTE_simulation/blob/11ff26bcdc26db1d1efa9c304b4df4d510085d08/Simulation.Rmd
# 获取日期：2026-09-23；原始文件保存在同目录 Simulation.Rmd。
# 转换原则：按原顺序提取全部 R 代码块，仅添加注释；不修改可执行表达式、
# 参数、随机种子、变量名或统计方法。下文“原文注意”指出需理解的原始实现。
#
# 学习目标与统一符号：
#   W / treatment：是否接受治疗，0 为对照，1 为治疗。
#   Y / outcome：二分类结局；这里没有规定 1 是有益事件还是不良事件。
#   X：治疗前协变量矩阵，一行对应一名模拟个体。
#   mu0(x)、mu1(x)：E[Y | X=x, W=0]、E[Y | X=x, W=1]。
#   e(x)：P(W=1 | X=x)，即倾向性评分。
#   m(x)：E[Y | X=x]，对治疗状态边缘化后的结局均值。
#   tau(x)：E[Y(1)-Y(0) | X=x]，即条件平均治疗效应（CATE）。
#   ATE：在目标人群中对 tau(X) 取平均；HTE：tau(x) 随 x 改变。
#   对本例的二分类结局，tau 位于风险差尺度，不是 OR、RR 或 HR。
#   例如 tau=0.10 表示结局发生概率增加 10 个百分点，不能直接称为“获益”。
#   每人的 tau.hat 是按其协变量预测的条件平均效应，并非可观测的个体因果差。
#
# 因果解释依赖一致性、可交换性和正值性等条件。本例随机生成治疗，因此
# 真实倾向性评分为 0.5；用于观察性数据时，机器学习不能自动消除未测混杂。
# 原文演示一次模拟数据上的拟合、诊断和解释，没有重复模拟、独立测试集，
# 也没有计算各方法相对真实 CATE 的偏倚、RMSE 或区间覆盖率。
#
# 使用方式：在 RStudio 中按第 1—10 步依次运行，章节依赖前面的对象。
# 如需整文件执行并显示未显式 print 的图表，可在项目根目录运行：
# source("learn/Simulation.R", encoding = "UTF-8", echo = TRUE, print.eval = TRUE)
# 所需包均列在第 1 步；library() 仅加载，不会安装。需要时可自行运行：
# install.packages(c("tidyverse", "grf", "kernelshap", "shapviz", "lmtest",
#                    "sandwich", "broom", "ggplot2", "gtsummary", "ranger", "bcf"))
# 完整流程包含 10000 人、因果森林调参、全样本 SHAP 和 MCMC，运行可能较久。
# bcf 默认可能向当前工作目录写入树样本和日志，详见第 10 步。
#
# 结果对象索引：
#   cf、data$tau.hat          因果森林及其训练样本袋外 CATE 预测；
#   forest.ate / calibration.plot  按 CATE 排序分组后的 AIPW 均值、标准误及图；
#   var.imp / shap_values / sv     分裂重要性、SHAP 分解与可视化对象；
#   data.s / data.t / data.x       S / T / X-learner 的 tau.hat 列；
#   data.dr / data.r / data.bcf    DR / R-learner / BCF 的 tau.hat 列。
# 注意：data.s$treatment 最终被统一改成 1；原始分组保留在 data0$treatment。
#
# 注释核对资料（接口及解释以实际安装版本为准）：
#   https://grf-labs.github.io/grf/reference/causal_forest.html
#   https://grf-labs.github.io/grf/reference/predict.causal_forest.html
#   https://grf-labs.github.io/grf/reference/test_calibration.html
#   https://grf-labs.github.io/grf/reference/variable_importance.html
#   https://grf-labs.github.io/grf/REFERENCE.html
#   https://imbs-hl.github.io/ranger/reference/ranger.html
#   https://cran.r-project.org/web/packages/bcf/bcf.pdf
# 第 5—9 步的原文参考：Salditt M, et al.
# A Tutorial Introduction to Heterogeneous Treatment Effect Estimation with Meta-learners.
# https://doi.org/10.1007/s10488-023-01303-9


# 阅读路线：原文说明 ----
### This file introduces HTE analysis using machine learning algorithms.
### This file is structured as follows:
### Step1: load packages
### Step2: test data generation
### Step3: estimation of CATE using causal forest
### Step4: calibration tests of causal forest
### Step5: S-learner
### Step6: T-learner
### Step7: X-learner
### Step8: DR-learner
### Step9: R-learner
### Step10: Bayesian causal forest

### Codes for Step5-9 were adopted from: Salditt M, et al. A Tutorial Introduction to Heterogeneous Treatment Effect Estimation with Meta-learners.(Adm Policy Ment Health. doi:10.1007/s10488-023-01303-9)

# 第 1 步：加载 R 包 ----
### Install packages
# tidyverse 提供 %>% 管道、dplyr 数据整理等工具。
# grf 用于广义随机森林，此处调用 causal_forest() 估计条件治疗效应。
# kernelshap 计算预测函数的 SHAP 分解；shapviz 负责展示结果。
# lmtest::coeftest() 与 sandwich::vcovHC() 联合给出异方差稳健标准误。
# broom 整理结果；ggplot2 绘图；gtsummary 生成分组描述表。
# ranger 拟合普通随机森林，是第 5—10 步的基础学习器；bcf 拟合贝叶斯因果森林。
# 原文标题写 Install packages，但下列代码实际上只是加载已安装的包。
library(tidyverse)
library(grf)
library(kernelshap)
library(shapviz)
library(lmtest)
library(sandwich)
library(broom)
library(ggplot2)
library(gtsummary)
library(ranger)
library(bcf)

# 第 2 步：生成随机治疗、协变量与二分类结局 ----
# Sample size
# 生成 10000 名个体。变量名虽然来自临床，但分布和系数仅用于模拟教学。
N <- 10000

# Treatment
# 设定随机种子，固定后续这段模拟的随机数序列。
# rbinom(N, 1, 0.5) 为每人独立生成一次成功概率 0.5 的伯努利试验。
# 治疗先独立生成，后续协变量不参与分配，故这里模拟的是随机治疗场景。
set.seed(1)
treatment <- rbinom(N, 1, 0.5)

# Covariates
# rnorm(n, mean, sd) 的第三个参数是标准差，不是方差。
# 年龄均值 40、标准差 10；收缩压均值 120、标准差 10；HbA1c 均值 5、标准差 1。
# 各协变量使用独立随机抽样；没有模拟真实临床变量之间的相关结构。
age <- rnorm(N, 40, 10)
systolic.blood.pressure <- rnorm(N, 120, 10)
hba1c <- rnorm(N, 5, 1)
# eGFR 先按均值 60、标准差 20 抽样，再把超界值压到 [10,120]。
# 这是截尾压界（会在边界堆积数值），不是从截断正态分布重新抽样。
eGFR <- rnorm(N, 60, 20)
eGFR <- pmin(pmax(eGFR, 10), 120)
# medication 是另一项 0/1 协变量，与本研究 treatment 不是同一变量。
medication <- rbinom(N, 1, 0.5)

# Parameters for HTE
# b 是“概率压界之前”的治疗系数；它随 HbA1c、eGFR 和 medication 增大。
# 除以 20 是原作者选择的效应缩放；年龄没有进入结局生成公式。
# 关键区别：后面会把结局概率压到 [0,1]，所以最终真实 CATE 未必等于 b。
b <- (0.1 + 0.05*hba1c + 0.05*eGFR + 0.3*medication) / 20

# Outcome
# 令 eta0 = -0.1 + 0.0005*(SBP-样本均值) + 0.01*HbA1c
#             + 0.001*(eGFR-样本均值) + 0.01*medication。
# 本行计算 eta0 + b*W，是概率尺度的加性表达式，没有使用 logit 链接。
# 两个减均值项按当前模拟样本中心化；regY1 实际用了每人的真实 W，
# 尽管名字含 Y1，它并不是所有人“均接受治疗”时的潜在结局概率。
regY1 <- -0.1 + b*treatment + 0.0005*(systolic.blood.pressure - mean(systolic.blood.pressure)) + 0.01*hba1c + 0.001*(eGFR- mean(eGFR)) + 0.01*medication
# 把概率限制到 [0,1]，避免 rbinom() 收到非法概率。
# 给定当前样本的中心化常数，真实风险差应是：
#   tau_true = clip(eta0+b) - clip(eta0)，clip(u)=min(max(u,0),1)。
# 仅在两种治疗状态的概率均未被压界时，tau_true 才等于 b。
# 压界还会使原本仅影响基线风险的收缩压影响最终风险差。
regY1 <- pmin(pmax(regY1, 0), 1)
# 按每人的事件概率产生一个 0/1 观测结局；并未同时观测 Y(0) 和 Y(1)。
outcome <- rbinom(N, 1, regY1)

# Create a data frame
# 整理成每行一人的数据框。后面增加 tau.hat 和 ranking 时，
# 会明确选取协变量列，避免把结局、治疗或估计效应混入 X。
data <- data.frame(treatment, outcome, age, systolic.blood.pressure, hba1c,  eGFR, medication)

# 第 3 步：拟合因果森林并预测 CATE ----
### Set up for cross-fitting
# num.rankings=5：按预测 CATE 分成五组；Q1 最低，Q5 最高。
# num.folds=10：构造十个分组标签，供森林 clusters 和组内分位数计算使用。
num.rankings <- 5
num.folds <- 10

# n 为样本量。%% 是取余；先取余、再排序、再加 1，得到 1—10 的折号。
# 这里没有随机打乱行号；本例每折 1000 人，按行顺序形成连续块。
# 原文称 cross-fitting，但实际只拟合一个指定 clusters=folds 的森林，
# 利用簇外/袋外预测；并非显式循环训练十个“留一折”的独立森林。
# 若换成按时间或中心排序的真实数据，需要重新考虑分折设计。
n = data %>% nrow()
folds <- sort(seq(n) %% num.folds) + 1

### Run causal forest ###
# grf 的输入分为：Y 结局向量、W 数值型治疗向量、X 数值协变量矩阵。
# 列顺序必须在训练、预测与变量名映射时保持一致。
Y <-  data$outcome
W <-  data$treatment
X <- data[, c("age", "systolic.blood.pressure", "hba1c", "eGFR", "medication")]
X <- as.matrix(X)

# 重新设置种子以固定该拟合步骤。2000 指主森林的树数；
# tune.parameters="all" 还会进行调参，实际工作量超过只长 2000 棵树。
# 未提供 Y.hat、W.hat 时，grf 会另估计 m(x) 和 e(x)。
# 原文也未显式设置 honesty；其默认开启，把树的分裂与叶内估计分开。
# clusters=folds 使抽样按折号作为簇处理，并影响 OOB 预测与后续推断。
# 这十个簇是算法构造的分组，不代表十个真实研究中心。
set.seed(1)
cf <- causal_forest(X = X, # covariate matrix
                           Y = Y, # outcome vector
                           W = W, # exposure vector
                           num.trees = 2000, # grow 2000 trees
                           tune.parameters = "all", # tune parameters
                           clusters = folds)

# Obtain estimated CATE
# predict(cf) 未传 newdata，返回训练样本的袋外（OOB）预测。
# $predictions 提取 CATE 向量；estimate.variance=F 关闭逐人方差估计。
# 原文的 F 在未被重新赋值时等价于 FALSE。
# 第一行存入 data，第二行另存同一类预测以便分组；不是两种不同估计量。
data$tau.hat <- predict(cf, estimate.variance = F)$predictions # predict CATEs
tau.hat = predict(cf)$predictions

# 第 4 步：校准、效应分组、变量重要性和 SHAP ----
### 1) BLP analysis ###
# 4.1 森林预测的最佳线性校准诊断（BLP）。
# mean.forest.prediction 的系数接近 1，支持平均效应尺度校准；
# differential.forest.prediction 的系数接近 1，支持异质性幅度校准。
# 后者显著大于 0 可作为存在异质性的证据；输出的单侧 P 值检验正信号，
# 不能把“小 P 值”直接理解成通过了“系数等于 1”的检验。
test_calibration(cf)



### 2) Calibration plot ###
# Obrain CATE quintile ranking
# 4.2 在每一折内部按 OOB CATE 的五分位数分组，再合并所有折。
# 先用长度 n 的空向量占位，循环只回填当前折对应的位置。
ranking <- rep(NA, n)
for (fold in seq(num.folds)){
# seq(0,1,by=1/5) 取 0%、20%、40%、60%、80%、100% 六个边界。
# cut() 据此产生 1—5 的有序区间标签；include.lowest=TRUE 包含最小值。
# 当前写法将 factor 的内部整数码填入 ranking，得到数值 1—5。
# 原文注意：预测值有大量并列时，分位数边界可能重复，使 cut() 报错。
# 组别是折内相对排名，合并后的同一 Q 不保证具有完全相同的数值边界。
  tau.hat.quintiles <- quantile(tau.hat[folds == fold], probs = seq(0, 1, by=1/num.rankings))
  ranking[folds == fold] <- cut(tau.hat[folds == fold], tau.hat.quintiles, include.lowest=TRUE,labels=seq(num.rankings))
}

# Computing AIPW scores
# 提取森林估计的三个辅助量：tau.hat、e.hat 和 m.hat。
# 它们分别对应条件效应、治疗概率和不区分治疗状态的结局均值。
tau.hat <- data$tau.hat
e.hat <- cf$W.hat # P [W=1|X]
m.hat <- cf$Y.hat # E [Y|X]
# Estimating mu.hat(X, 1) and mu. hat (X, 0) for observations
# 由 m=e*mu1+(1-e)*mu0 以及 tau=mu1-mu0 可得：
#   mu0 = m-e*tau；mu1 = m+(1-e)*tau。
# 原文 mu.hat.O 的后缀是大写字母 O，不是数字 0；仍表示对照结局均值。
# 这些代数重构值没有再压到 [0,1]，有限样本下可能超出概率范围。
mu.hat.O <- m.hat - e.hat * tau.hat
mu.hat.1 <- m.hat + (1 - e.hat) * tau.hat
# AIPW scores
# AIPW（增广逆概率加权）伪结局：
#   Gamma = tau + W/e*(Y-mu1) - (1-W)/(1-e)*(Y-mu0)。
# 治疗组用治疗结局残差校正，对照组用对照结局残差校正；校正项带逆概率权重。
# 在适当的识别和估计条件下，Gamma 的条件均值对应 CATE。
# 本步没有对 e.hat 截断；若接近 0 或 1，校正项可能很不稳定。
# 单个 Gamma 可能超出 [-1,1]，它是伪结局而非一个人的已知真实风险差。
aipw.scores <- tau.hat + data[,"treatment"]/(e.hat)*(data[,"outcome"] - mu.hat.1) - (1-data[,"treatment"])/(1-e.hat)*(data[,"outcome"] - mu.hat.O)
# 无截距回归 ~0+factor(ranking) 为每个 Q 设置一个指示变量。
# 因此每个回归系数就是该组 AIPW 分数的均值，即按预测效应分组的效应估计。
ols <- lm(aipw.scores ~ 0 + factor(ranking))

# vcovHC(...,"HC2") 给出异方差稳健协方差；coeftest()[,1:2]
# 取系数和标准误。组名与方法名一起存成 forest.ate。
# 这里使用的是普通 HC2，没有把 folds 传给聚类稳健协方差函数；
# 它也不与 test_calibration(cf) 的推断设置自动保持一致。
forest.ate <- data.frame ("AIPW", paste0("Q", seq(length(unique(ranking)))), coeftest(ols, vcov = vcovHC(ols, "HC2")) [,1:2])
colnames(forest.ate) <- c("method", "ranking", "estimate", "std.err")
rownames(forest.ate) <- NULL #forest.ate

# Plot estimated ATEs by CATE quintile
# 横轴 Q1—Q5，纵轴为各组 AIPW 均值，误差线为估计值 ±2 个标准误。
# 在正态近似下可近似看作 95% 区间，严格 95% 正态区间常用 1.96。
# 这幅图展示按效应排序的组间差异，没有把组内平均预测值与 AIPW 均值
# 同时画出，因此不是“预测值对观测值加 45 度线”的传统校准图。
# 原文注意：ylab 写 percentage point，但数据没有乘以 100；
# 图上 0.1 实际表示 10 个百分点。此处保留原标签，避免悄悄修改原代码。
calibration.plot <- ggplot(forest.ate) +
  aes (x = ranking, y = estimate, group = method, color = method) +
  geom_point(position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = estimate -2 * std.err, ymax = estimate +2 * std.err), width = .2, position = position_dodge (0.2)) +
  ylab("ATE (percentage point)") + xlab("CATE Ranking") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())

# 显式 print() 在 source() 中也会画出此图；其余裸图表表达式需要
# 交互式逐段运行，或 source(..., print.eval=TRUE) 才会自动打印。
print(calibration.plot)



### 3) Variable importance ###
# 4.3 分裂变量重要性。名称顺序与 X 的五列完全一致。
covariate.list <- c("age", "systolic.blood.pressure", "hba1c", "eGFR", "medication")

# Obtain variable importance
# variable_importance() 汇总变量在森林分裂中使用的频率并按树深度加权。
# 这不是 SHAP、不是置换重要性，也不是某变量存在治疗交互作用的 P 值。
# tidy()、data.frame() 整理输出；x 是原文所用的重要性数值列名。
# mutate() 添加名称，arrange(desc(x)) 按重要性从高到低排列。
# 若更换 broom/grf 版本后列结构不同，应先检查结果结构再调整映射。
# 本次核验中 tidy.numeric() 提示已弃用，但仍返回原文需要的 x 列；
# 此处保留原调用以保持代码一致，未来升级时需要留意该接口。
var.imp <- variable_importance(cf) %>%
  tidy() %>%
  data.frame() %>%
  mutate(varname = covariate.list) %>%
  arrange(desc(x))

# Order the factor levels by value descending
# 因子水平按数值升序设置，配合 coord_flip() 后较大值显示在图的上方。
var.imp$varname <- factor(var.imp$varname, levels = var.imp$varname[order(var.imp$x, decreasing = FALSE)])

# Barplot of variable importance
# stat="identity" 直接使用给定的重要性数值作为柱高，不是计数柱状图。
# coord_flip() 交换横纵轴；后续 theme 设置背景及坐标文字颜色。
ggplot(var.imp, aes(x=varname, y=x)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ylab("Importance") +
  xlab("Variables") +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  )



### 4) SHAP Value ###
# Obtain background data by random sampling 100 rows from X
# 4.4 SHAP：从 X 中不放回抽取 100 人作为背景，重新设种子便于复现。
# 背景样本决定“相对于什么基线人群”解释预测，改变背景会改变归因。
set.seed(1)
X_bg <- X[sample(nrow(X), 100),]

# Estimate SHAP values.
# kernelshap 用背景数据替换部分特征，调用指定预测函数，估计 SHAP。
# 这里解释的是森林的 CATE 预测，而不是结局风险；X=X 表示解释全部 10000 人。
# 计算量可能较大。对每人，基线预测加各变量 SHAP 之和近似还原该人的预测。
# 原文注意：pred_fun 显式传入 X，使用的是新数据预测接口，而非
# predict(cf) 的训练样本 OOB 路径；解释值不应直接当作 OOB 校准结果。
# SHAP 描述模型如何使用特征，不能据此断言干预某协变量会改变治疗效应。
shap_values <- kernelshap(
  object = cf, # the causal forest model from above.
  X = X, # the feature matrix
  bg_X = X_bg, # the background data
  pred_fun = function(object, X) predict(object, X)$predictions # a custom prediction function to compute CATEs
)

# Convert the output to a shapviz object for plotting
# 把 kernelshap 结果转换成 shapviz 对象，保留 SHAP 矩阵和特征值。
sv <- shapviz(shap_values)

# Plot SHAP values
# sv_importance() 默认柱图汇总平均绝对 SHAP 值，衡量预测归因的大小。
# 绝对值汇总不保留正负方向；正负作用需查看个体值、蜂群图或依赖图。
sv_importance(sv, bar_width = 0.6) +
  ggtitle("SHAP Values") +
  theme_minimal()



### 5) Comparison of CATE quintile ###
# Label CATE ranking
# 4.5 把 CATE 分位组写回数据，用于比较不同预测效应人群的特征。
data$ranking <- ranking

# Summarize each CATE subgroup
# tbl_summary() 按 ranking 分列，连续变量展示均值（标准差），保留两位小数。
# medication 等二分类变量按 gtsummary 的自动类型规则展示；bold_labels() 加粗标签。
# 这里只做描述性汇总，没有组间检验，也不能把组间特征差异等同于效应修饰证据。
data[, c("ranking", "age", "systolic.blood.pressure", "hba1c", "eGFR", "medication")] %>%
  tbl_summary(statistic = list(all_continuous() ~ "{mean} ({sd})"),
              digits = all_continuous() ~ 2,
              by = ranking) %>%
  bold_labels()

# 第 5 步：S-learner——一个结局模型，切换治疗状态 ----
# Create sample for this section
# 从 data 重新选取原始七列，排除刚增加的 tau.hat 和 ranking。
# data0 是后续各学习器共用的数据；covariates 只含五个治疗前特征。
# 共同约定：ranger 的数值型 0/1 outcome 默认走回归森林，估计条件均值。
# 对象的 $predictions 是训练样本 OOB 预测；predict(model,newdata) 则用保存的森林。
# keep.inbag=TRUE 保存每棵树的入袋次数，不等于显式实施多折交叉拟合。
# seed=1 控制该次森林拟合，num.threads=10 要求使用十个线程。
data0 <- data[c("treatment", "outcome", "age", "systolic.blood.pressure", "hba1c", "eGFR", "medication")]
covariates <- c("age", "systolic.blood.pressure", "hba1c", "eGFR", "medication")

# S 指 single：用一个模型学习 mu(x,w)，把 W 作为普通特征与 X 一起输入。
# data.s 是副本，下面修改它的治疗列不会改变 data0。
data.s <- data0 # Replicate data

# Build outcome model with random forest
# 训练一个包含治疗特征的结局回归森林。
# 当治疗特征很少被树选中时，S-learner 可能把估计效应压向零。
s_learner <- ranger(y = data.s$outcome, x = data.s[, c("treatment", covariates)], keep.inbag = TRUE, seed = 1, num.threads = 10)

# Predict outcome when untreated
# 构造“所有人均未治疗”的预测数据。对真实对照者，直接取其 OOB 预测；
# 对真实治疗者，使用把 W 改为 0 后的预测，插补未治疗时的结局均值。
# 后者来自含该训练个体的完整森林，不应一概称为严格样本外预测。
# 原始治疗状态始终从 data0$treatment 读取，不能使用已被覆盖的 data.s$treatment。
data.s$treatment <- 0
mu0.hat.s <- rep(0, nrow(data.s))
mu0.hat.s[data0$treatment == 0] <- s_learner$predictions[data0$treatment == 0]
mu0.hat.s[data0$treatment == 1] <- predict(s_learner, data.s)$predictions[data0$treatment == 1]

# Predict outcome when treated
# 同理构造“所有人均接受治疗”的数据：真实治疗者取 OOB 预测，
# 真实对照者取 W 改为 1 后的反事实预测。
# rep(0,nrow(...)) 先分配长度 n 的向量，再按原始分组回填。
data.s$treatment <- 1
mu1.hat.s <- rep(0,  nrow(data.s))
mu1.hat.s[data0$treatment == 1] <- s_learner$predictions[data0$treatment == 1]
mu1.hat.s[data0$treatment == 0] <- predict(s_learner, data.s)$predictions[data0$treatment == 0]

# Calculate CATE
# 同一人的两种预测结局相减：tau_S(x)=mu_hat(x,1)-mu_hat(x,0)。
# 此时 data.s$treatment 全为 1，这是预测数据的构造结果，不是原始分组。
data.s$tau.hat <- mu1.hat.s - mu0.hat.s

# 第 6 步：T-learner——两组分别拟合结局模型 ----
# T 指 two：对照组与治疗组分别拟合 mu0(x) 和 mu1(x)。
# 两个模型只使用 X，不再把 treatment 当特征；各模型只能使用本组样本。
data.t <- data0 # Replicate data
data.t_0 <- data.t[data.t$treatment ==0,] # Control group
data.t_1 <- data.t[data.t$treatment ==1,] # Treatment group

# Build outcome model with random forest to predict mu0
# 先在对照组学习 mu0(x)：本组个体用 OOB 均值，治疗组用外推预测。
# 逻辑索引按原始行顺序回填，使 mu0.hat.t 与 data.t 的每行对应。
t_learner_m0 <- ranger(y = data.t_0$outcome, x = data.t_0[, covariates], keep.inbag = TRUE, seed = 1, num.threads = 10)
mu0.hat.t <- rep(0, nrow(data.t))
mu0.hat.t[data.t$treatment == 0] <- t_learner_m0$predictions
mu0.hat.t[data.t$treatment == 1] <- predict(t_learner_m0, data.t_1)$predictions

# Build outcome model with random forest to predict mu1
# 再在治疗组学习 mu1(x)：本组用 OOB 均值，对照组用外推预测。
# 两组模型独立，样本较少的组可能估计不稳定；本例分配约为 1:1。
t_learner_m1 <- ranger(y = data.t_1$outcome, x = data.t_1[, covariates], keep.inbag = TRUE, seed = 1, num.threads = 10)
mu1.hat.t <- rep(0, nrow(data.t))
mu1.hat.t[data.t$treatment == 1] <- t_learner_m1$predictions
mu1.hat.t[data.t$treatment == 0] <- predict(t_learner_m1, data.t_0)$predictions

# Calculate CATE
# 逐人取两个结局模型预测之差，得到 T-learner 的条件风险差。
data.t$tau.hat <- mu1.hat.t - mu0.hat.t

# 第 7 步：X-learner——插补效应，再分组建模并加权 ----
# X-learner 分三层：先拟合两组结局，再构造效应伪结局并分组建模，
# 最后按倾向性评分组合两种效应预测。它比 T-learner 多了效应回归这一层。
data.x <- data0 # Replicate data
data.x_0 <- data.x[data.x$treatment ==0,] # Control group
data.x_1 <- data.x[data.x$treatment ==1,] # Treatment group

# Build outcome model with random forest to predict mu0
# 第一层：在对照组拟合 mu0(x)，按“本组 OOB、另一组外推”回填全样本。
x_learner_m0 <- ranger(y = data.x_0$outcome, x = data.x_0[, covariates], keep.inbag = TRUE, seed = 1, num.threads = 10)
mu0.hat.x <- rep(0, nrow(data.x))
mu0.hat.x[data.x$treatment == 0] <- x_learner_m0$predictions
mu0.hat.x[data.x$treatment == 1] <- predict(x_learner_m0, data.x_1)$predictions

# Build outcome model with random forest to predict mu1
# 在治疗组拟合 mu1(x)，同样保留按原始行顺序排列的全样本预测。
# mu0.hat.x、mu1.hat.x 在原文中被计算和保留，但下一步直接重新调用
# predict() 构造伪结局，没有使用这两个完整向量。
x_learner_m1 <- ranger(y = data.x_1$outcome, x = data.x_1[, covariates], keep.inbag = TRUE, seed = 1, num.threads = 10)
mu1.hat.x <- rep(0, nrow(data.x))
mu1.hat.x[data.x$treatment == 1] <- x_learner_m1$predictions
mu1.hat.x[data.x$treatment == 0] <- predict(x_learner_m1, data.x_0)$predictions

# Compute the pseudo-outcome via the estimated outcome models
# 第二层的效应伪结局：
#   对照者 D0 = mu1_hat(X)-Y，用插补的治疗结局减去实际对照结局；
#   治疗者 D1 = Y-mu0_hat(X)，用实际治疗结局减去插补的对照结局。
# 两个方向都对应“治疗减对照”，但个体观测噪声仍在，D0/D1 不是真实个体效应。
psi.x0 <- predict(x_learner_m1, data.x_0)$predictions - data.x_0$outcome
psi.x1 <- data.x_1$outcome - predict(x_learner_m0, data.x_1)$predictions

# Fit models for the pseudo-outcome for each treatment group
# 分别在两组回归 D0 和 D1，得到 tau0(x) 与 tau1(x)。
# 此处的 0/1 指效应伪结局来自哪一组，不是两个潜在结局模型。
x_learner_tau0 <- ranger(y = psi.x0, x = data.x_0[, covariates], keep.inbag = TRUE, seed = 1, num.threads = 10)
x_learner_tau1 <- ranger(y = psi.x1, x = data.x_1[, covariates], keep.inbag = TRUE, seed = 1, num.threads = 10)

# Compute treatment effect for each treatment group
# 对全样本获得两种效应预测：建模组用各自 OOB，另一组用完整森林预测。
# 原文索引使用全局 treatment，它与 data.x$treatment 在当前流程中相同；
# 若先筛选或重排 data.x，必须同步索引，否则会错配个体。
tau0.hat.x <- rep(0, nrow(data.x))
tau0.hat.x[treatment == 0] <- x_learner_tau0$predictions
tau0.hat.x[treatment == 1] <- predict(x_learner_tau0, data.x_1)$predictions
tau1.hat.x <- rep(0, nrow(data.x))
tau1.hat.x[treatment == 1] <- x_learner_tau1$predictions
tau1.hat.x[treatment == 0] <- predict(x_learner_tau1, data.x_0)$predictions

# Calculate propensity score
# 第三层：用概率森林估计 e(x)=P(W=1|X=x)。probability=TRUE 返回各类概率。
# 本例 0/1 标签对应预测矩阵的列 "0"、"1"，[,2] 因而取治疗概率；
# 换用其他治疗编码时不能假定第二列始终代表治疗组。
# 这里没有显式 seed；结果依赖执行到此处的随机数状态，单独重跑未必相同。
# 本模拟的真实 e(x)=0.5；原文仍拟合倾向性评分模型以演示通用流程。
ps.x <- ranger(y = data.x$treatment, x = data.x[, covariates], probability = TRUE)
ps.hat.x <- ps.x$predictions[,2]
# Ensure positivity
# 将估计倾向性评分限制到 [0.01,0.99]。
# 截断可以缓解数值极端，但不能让真实数据中缺乏重叠的区域变得可识别。
epsilon <- 0.01
ps.hat.x <- ifelse(ps.hat.x < epsilon, epsilon, ifelse(ps.hat.x > 1 - epsilon, 1 - epsilon, ps.hat.x))

# Compute CATE as propensity score-weighted combination of the group-specific estimates
# 组合公式 tau_X(x)=e_hat(x)*tau0_hat(x)+(1-e_hat(x))*tau1_hat(x)。
# 注意权重对应：治疗概率乘“从对照组伪结局学到的”tau0，而非 tau1。
data.x$tau.hat <- ps.hat.x * tau0.hat.x + (1 - ps.hat.x) * tau1.hat.x

# 第 8 步：DR-learner——对双重稳健伪结局进行回归 ----
# DR 指 doubly robust。先估计两种结局均值和治疗概率，再构造
# AIPW 型伪结局，最后让随机森林学习该伪结局如何随 X 变化。
data.dr <- data0 # Replicate data
data.dr_0 <- data.dr[data.dr$treatment ==0,] # Control group
data.dr_1 <- data.dr[data.dr$treatment ==1,] # Treatment group

# Build outcome model with random forest to predict mu0
# 对照结局模型 mu0：对照者取 OOB 预测，治疗者取对照模型外推值。
dr_learner_m0 <- ranger(y = data.dr_0$outcome, x = data.dr_0[, covariates], keep.inbag = TRUE, seed = 1, num.threads = 10)
mu0.hat.dr <- rep(0, nrow(data.dr))
mu0.hat.dr[data.dr$treatment == 0] <- dr_learner_m0$predictions
mu0.hat.dr[data.dr$treatment == 1] <- predict(dr_learner_m0, data.dr_1)$predictions

# Build outcome model with random forest to predict mu1
# 治疗结局模型 mu1：治疗者取 OOB 预测，对照者取治疗模型外推值。
dr_learner_m1 <- ranger(y = data.dr_1$outcome, x = data.dr_1[, covariates], keep.inbag = TRUE, seed = 1, num.threads = 10)
mu1.hat.dr <- rep(0, nrow(data.dr))
mu1.hat.dr[data.dr$treatment == 1] <- dr_learner_m1$predictions
mu1.hat.dr[data.dr$treatment == 0] <- predict(dr_learner_m1, data.dr_0)$predictions

# Calculate propensity score
# 概率森林估计 e(x)，取预测矩阵的治疗类别列；和 X-learner 一样
# 依赖 0/1 编码与列顺序，且该调用没有显式指定 seed。
ps.dr <- ranger(y = data.dr$treatment, x = data.dr[, covariates], probability = TRUE)
ps.hat.dr <- ps.dr$predictions[,2]
# Ensure positivity
# 将 e_hat 限于 [0.01,0.99]，避免后续逆概率权重的分母趋近于零。
epsilon <- 0.01
ps.hat.dr <- ifelse(ps.hat.dr < epsilon, epsilon, ifelse(ps.hat.dr > 1 - epsilon, 1 - epsilon, ps.hat.dr))

# Compute the pseudo-outcome of the DR-learner
# 增广残差项为 W/e*(Y-mu1) - (1-W)/(1-e)*(Y-mu0)。
# 治疗者只贡献第一项，对照者只贡献第二项。
augmented.term <- 1/ps.hat.dr * (data.dr$treatment * (data.dr$outcome - mu1.hat.dr)) -
  1/(1 - ps.hat.dr) * ((1 - data.dr$treatment) * (data.dr$outcome - mu0.hat.dr))
# psi_DR = mu1_hat-mu0_hat + 增广残差项。
# 双重稳健指在相关识别、正则与估计条件下，倾向性评分正确或两组结局均值
# 均正确时，伪结局可获得正确的条件均值；并非任何一个模型随意拟合就保证准确。
# 最后的效应回归也要能学习这个条件均值。原文使用 OOB 辅助预测，
# 没有展示独立分折的完整多阶段交叉拟合，不能直接据此宣称推断已有效。
psi.dr <- mu1.hat.dr - mu0.hat.dr + augmented.term

# Fit a random forest model to the pseudo-outcome
# 把连续的伪结局作为 y，只用协变量 X 回归；这里不再加入 treatment。
# 伪结局可能超出风险差范围，回归森林本身不强制 CATE 落在 [-1,1]。
tau.dr <- ranger(y = psi.dr, x = data.dr[, covariates], keep.inbag = TRUE, seed = 1, num.threads = 10)

#Compute the CATE as the predictions from the pseudo-outcome regression
# 取最终效应森林的训练样本 OOB 预测，作为 DR-learner 的 tau.hat。
data.dr$tau.hat <- tau.dr$predictions

# 第 9 步：R-learner——对治疗与结局残差进行建模 ----
# R 指 Robinson 残差化思路。目标是学习 tau(x)，使
#   Y-m_hat(X) 约等于 (W-e_hat(X))*tau(X)。
# 与 T/DR 不同，这里先估计总体结局均值 m(x)，没有分别拟合 mu0 和 mu1。
data.r <- data0 # Replicate data

# Build outcome model with random forest to predict mu
# 结局模型只输入 X，不输入 W，因此学习 m(x)=E[Y|X=x]。
# $predictions 提供 OOB 预测。这里重用了 m.hat 名称，会覆盖第 4 步的同名对象。
r_learner <- ranger(y = data.r$outcome, x = data.r[, covariates], keep.inbag = TRUE, seed = 1, num.threads = 10)
m.hat <- r_learner$predictions

# Calculate propensity score
# 估计治疗概率并取治疗类别列，随后压到 [0.01,0.99]。
# 真实模拟倾向性评分仍是 0.5；这里按原文使用估计值。
ps.r <- ranger(y = data.r$treatment, x = data.r[, covariates], probability = TRUE)
ps.hat.r <- ps.r$predictions[,2]
# Ensure positivity
epsilon <- 0.01
ps.hat.r <- ifelse(ps.hat.r < epsilon, epsilon, ifelse(ps.hat.r > 1 - epsilon, 1 - epsilon, ps.hat.r))

# Compute the pseudo-outcome
# 构造治疗残差 W-e_hat(X) 和结局残差 Y-m_hat(X)。
# 二者相除形成 psi_R；若治疗残差很小，伪结局的绝对值可能很大。
# 对当前 0/1 治疗和截断概率，治疗残差不会恰好为零。
resid.treat <- data.r$treatment - ps.hat.r
resid.out <- data.r$outcome - m.hat
psi.r <- resid.out / resid.treat

# Compute weight
# 理论权重 w=(W-e_hat)^2，用于降低近零治疗残差对应的大伪结局的影响。
# 加权平方损失 sum[w*(psi_R-tau(X))^2] 在代数上等于
# sum[(Y-m_hat-(W-e_hat)*tau(X))^2]，这就是残差化目标的由来。
w <- resid.treat ^ 2

# Regress pseudo-outcome on covariates using weights w
# 原文把理论权重传入 ranger 的 case.weights。
# 重要接口区别：ranger 的 case.weights 控制训练样本被抽入每棵树的概率，
# 不是逐个节点直接优化上述精确加权平方损失的开关。
# 因此这里应理解为用加权抽样森林实现的近似做法，不能仅凭参数名称断言
# 它与严格最小化加权 R-loss 的实现完全等价。
# 加权抽样还会提高高权重个体的入袋频率；树太少时，这些个体可能没有
# 可用的袋外树，导致 OOB 预测缺失。因此缩小示例时不宜把树数降得过低。
tau.r <- ranger(y = psi.r, x = data.r[, covariates], case.weights = w, keep.inbag = TRUE, seed = 1, num.threads = 10)

#Compute the CATE as the predictions from the weighted pseudo-outcome regression
# 取效应森林的 OOB 预测，写入 data.r$tau.hat。
data.r$tau.hat <- tau.r$predictions

# 第 10 步：贝叶斯因果森林（BCF） ----
# BCF 同时学习基线预后函数和治疗效应函数，采用不同的树集合与先验。
# 其基本模型为 Y=mu(X,e_hat(X))+tau(X)*W+误差，通常使用正态误差。
# 原文把二分类 outcome 直接传入该连续结局模型，并未指定二项分布或 logit 链接；
# 应将其视为原文的均值建模演示，不能称为严格的二分类概率模型。
data.bcf <- data0 # Replicate data

# Calculate propensity score
# 先用概率森林估计 e(x)，BCF 将其作为输入，有助于处理预后与治疗分配的关系。
# 本例治疗随机分配，所以真实评分恒为 0.5；截断流程与前面一致。
ps.bcf <- ranger(y = data.bcf$treatment, x = data.bcf[, covariates], probability = TRUE)
ps.hat.bcf <- ps.bcf$predictions[,2]
# Ensure positivity
epsilon <- 0.01
ps.hat.bcf <- ifelse(ps.hat.bcf < epsilon, epsilon, ifelse(ps.hat.bcf > 1 - epsilon, 1 - epsilon, ps.hat.bcf))

# Run Bayesian causal forest
# y 是观测结局；z 是 0/1 治疗。
# x_control 是基线预后函数的特征，x_moderate 是效应函数的特征，
# 二者这里使用同样五列，且必须转为数值矩阵。
# pihat 是估计治疗概率；默认 include_pi="control" 将其加入预后部分。
# nburn=100 为丢弃的预热迭代；nsim=100 为保留的后验抽样次数。
# 这些次数很少，不能仅凭这段代码声称链已收敛或区间估计稳定。
# random_seed 固定该步随机数，n_threads 指定线程数。
# 原文注意：bcf 默认可在工作目录保存树样本/日志；执行前应了解输出位置。
# bcf <- bcf(...) 把结果命名为 bcf，可用 bcf::bcf() 明确引用包内函数。
bcf <- bcf(y = data.bcf$outcome, # Outcome
           z = data.bcf$treatment, # Treatment
           x_control = as.matrix(data.bcf[, covariates]), # Covariates for prognostic score function
           x_moderate = as.matrix(data.bcf[, covariates]), # Covariates for tau function
           pihat = ps.hat.bcf, # Propensity score
           nburn = 100, # The number of burn-in iterations
           nsim = 100, # The number of iterations used for CATE estimation
           random_seed = 1,
           n_threads = 10)

# bcf$tau 的每列对应一个人，每行对应一次保留的后验抽样。
# colMeans() 对每人的后验样本取均值，得到训练个体的 CATE 后验均值。
# 这不是 OOB 预测，与前述森林的 OOB 评估口径不同。
# 本段未计算后验可信区间，也未绘制轨迹图或评估收敛；后验均值不能替代这些诊断。
data.bcf$tau.hat <- colMeans(bcf$tau)
