###########################################################
#            倾向得分四种常用用法 —— R 全流程演示
#    Matching / Stratification / Adjustment / IPTW
#
# 研究问题（示例）：
#   暴露：PA_CAT（身体活动是否为高：1=高，0=非高/低）
#   结局：lungfun（二分类结局）
#   目标：估计“高身体活动 vs 非高/低身体活动”对 lungfun 的平均处理效应（ATE）
#
# 整体流程：
#  0. 数据读取与变量类型处理
#  1. 单因素筛选协变量（教学演示用），确定进入 PS 模型的变量
#  2. Logistic 倾向得分模型 + PS 分布图（检查重叠/支持性）
#  3.0 未调整 ATE（基线对照）
#  3.1 PS 匹配：最近邻 + 卡钳 + 平衡性 + ATE
#  3.2 PS 分层：4 分位分层 + 平衡性 + 分层加权 ATE
#  3.3 PS 校正：结局模型中调整 PS
#  3.4 IPTW：生成 ATE 权重 + 平衡性 + ATE
#  4. 汇总：未调整 + 4 种 PS 方法的 ATE / CI / P 值
###########################################################

##---------------------------------------------------------
## -1. 依赖包：如未安装则先安装（一次性）
##---------------------------------------------------------
pkgs_needed <- c(
  "tableone",  # 基线表 / SMD
  "ggplot2",   # 作图
  "MatchIt",   # PS 匹配
  "cobalt",    # 平衡性诊断（SMD、Love plot）
  "WeightIt",  # PS 分层权重 / IPTW 权重
  "survey"     # 加权回归（svyglm）
)

pkgs_to_install <- pkgs_needed[!vapply(pkgs_needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(pkgs_to_install) > 0) {
  install.packages(pkgs_to_install, dependencies = TRUE)
}

###########################################################
## 0. 读取数据 & 变量预处理
###########################################################

# 设置当前工作目录（按自己的电脑路径调整）
setwd("C:/Users/77489/Desktop/倾向得分培训班+训练营_25.12.23/录课")

# 读取数据（data.csv 需放在上述工作目录下）
data <- read.csv("data.csv")

# 候选协变量（潜在混杂因素）
# 说明：这些变量将用于构建倾向得分模型（PS = P(A=1|X)）
names(data)
covariates <- c("sex","hypertension","diabetes","heart_disease","residence","marital",
                "smoking","age","education_level","cardio_metabolic",
                "alcohol","bmi")

# 需要按分类变量处理的协变量（转为 factor，便于后续自动生成哑变量/分层比较）
categorical_vars <- c("sex","hypertension","diabetes","heart_disease","residence","marital",
                      "smoking","education_level","cardio_metabolic",
                      "alcohol")

# 统一把分类协变量转为 factor（否则可能被当作连续变量进入模型）
data[categorical_vars] <- lapply(data[categorical_vars], factor)

# 暴露变量检查（应为 0/1 或二分类因子）
table(data$PA_CAT)


###########################################################
# 1. 单因素分析 —— 筛选进入 PS 模型的协变量（教学演示）
#
# 目的（教学演示用）：
#   仅为了演示“变量筛选”的流程；真实研究中更推荐基于领域知识/因果图(DAG)
#   来确定混杂因素，而不是完全依赖 P 值筛选。
#
# 做法：
#   以 lungfun 为因变量，对每个协变量做单因素回归；
#   P < 0.05 的协变量进入后续 PS 模型。
#
# 说明：
#   - 分类变量：查看所有虚拟变量（不含截距）的 P 值，只要有一个水平显著则纳入
#   - 连续变量：查看该变量的系数 P 值
###########################################################

# 保存“通过筛选”的协变量名
sig_vars <- c()

# 保存单因素结果（用于展示：变量类型、P 值、是否显著）
result_table <- data.frame()

# 循环：逐个协变量做单因素回归
# v=covariates[1]
for (v in covariates) {
  # 单因素回归：lungfun ~ 协变量
  # 用 data[[v]] 便于在循环中按列名取变量
  fit <- glm(lungfun ~ data[[v]], data = data, family = "binomial")
  
  # 回归系数表（估计值、标准误、统计量、P 值）
  coef_table <- summary(fit)$coefficients
  
  # 分类变量：看所有 dummy 项的 P 值（不含截距）
  if (v %in% categorical_vars) {
    pvals <- coef_table[-1, 4]
    sig <- ifelse(any(pvals < 0.05), "Yes", "No")
    
    result_table <- rbind(
      result_table,
      data.frame(variable   = v,
                 type       = "categorical",
                 p_values   = paste(round(pvals,4), collapse=", "),
                 significant= sig)
    )
    if (sig=="Yes") sig_vars <- c(sig_vars, v)
    
  } else {
    # 连续变量：只有一个系数（第二行）
    pval <- coef_table[2,4]
    sig  <- ifelse(pval < 0.05, "Yes", "No")
    
    result_table <- rbind(
      result_table,
      data.frame(variable   = v,
                 type       = "continuous",
                 p_values   = round(pval,4),
                 significant= sig)
    )
    if (sig=="Yes") sig_vars <- c(sig_vars, v)
  }
}

# 输出单因素筛选结果
cat("\n==== 单因素分析结果 ====\n")
print(result_table)

cat("\n进入 PS 模型的协变量：\n")
print(sig_vars)


##---------------------------------------------------------
## 原始基线表（未做任何 PS 处理前）
## 目的：展示 PA_CAT 两组在协变量上的初始不平衡（可看 SMD/检验）
##---------------------------------------------------------
library(tableone)

# sig_vars 中属于分类变量的部分（CreateTableOne 建议显式指定）
categorical_vars_sig <- intersect(sig_vars, categorical_vars) 

# Table 1：按暴露 PA_CAT 分组比较协变量分布
table1 <- CreateTableOne(
  vars       = sig_vars,
  strata     = "PA_CAT",
  data       = data,
  factorVars = categorical_vars_sig
)

print(table1, showAllLevels = TRUE, test = TRUE)



###########################################################
# 2. 倾向得分模型（Logistic 回归）
#
# 目标：
#   估计每个个体“高身体活动”的概率
#   PS = P(PA_CAT=1 | X)
#
# 模型：
#   PA_CAT ~ sig_vars
#
# 关键检查：
#   画 PS 分布图，看两组是否有足够重叠（支持性/可比性/positivity）
###########################################################

library(ggplot2)

# 动态拼接 PS 模型公式
formula_ps <- as.formula(
  paste("PA_CAT ~", paste(sig_vars, collapse=" + "))
)

# Logistic 回归拟合 PS
ps_model <- glm(formula_ps, data=data, family=binomial())

# 预测倾向得分（概率尺度）
data$pscore <- predict(ps_model, type="response")

# PS 分布图：观察两组的重叠（重叠越好，支持性越强）
ggplot(data, aes(pscore, color=factor(PA_CAT), fill=factor(PA_CAT))) +
  geom_density(alpha=0.3, size=1.2) +
  scale_color_manual(values=c("#1f77b4","#ff7f0e"),
                     labels = c("0 = 非高/低", "1 = 高")) +
  scale_fill_manual(values=c("#1f77b4","#ff7f0e"),
                    labels = c("0 = 非高/低", "1 = 高")) +
  labs(title="倾向得分分布（原始）", x="倾向得分", y="密度") +
  theme_minimal(base_size = 14)+
  theme(legend.title = element_blank())



###########################################################
# 3.0 未调整的平均处理效应（原始数据）
#
# 作为对照（粗略关联，不控制混杂）：
#   ATE = E(lungfun|A=1) - E(lungfun|A=0)
#
# 说明：
#   t.test(y ~ A) 的默认“差值方向”与因子水平顺序有关；
#   这里统一输出为（A=1 - A=0）
###########################################################

# 均值差（方向：高身体活动 - 非高/低）
ATE_unadj <- with(data,
                  mean(lungfun[PA_CAT == 1]) - mean(lungfun[PA_CAT == 0]))

# t 检验：给出均值差的 CI 和 P 值
t_unadj <- t.test(lungfun ~ PA_CAT, data = data)

cat("\n==== 未调整 ATE ====\n")
cat("ATE =", round(ATE_unadj, 3), "\n")
# t.test 默认输出为（组0 - 组1），这里统一成（组1 - 组0）
cat("95% CI = [", round(-t_unadj$conf.int[2], 3), ",",
    round(-t_unadj$conf.int[1], 3), "]\n")
cat("P值 =", round(t_unadj$p.value, 4), "\n")



###########################################################
# 3. 四种倾向得分方法
#    3.1 匹配 / 3.2 分层 / 3.3 校正 / 3.4 IPTW
#
# 平衡性诊断常用口径：
#   - SMD（标准化均值差）绝对值 < 0.1 作为经验阈值
#   - Love plot 用于可视化“处理前 vs 处理后”的整体改善情况
###########################################################

###########################################################
# 3.1 方法①：倾向得分匹配（PS Matching）
#
# 设置要点：
#   method="nearest"  : 最近邻匹配
#   caliper=0.2       : 卡钳宽度（限制可匹配的距离，避免“硬配”）
#   replace=FALSE     : 不放回（每个对照最多匹配一次）
#   distance="logit"  : 匹配距离基于 logit(PS)
###########################################################

library(MatchIt)
sig_vars
m.out <- matchit(
  formula_ps,
  data = data,
  method = "nearest",
  distance = "logit",
  caliper = 0.2,
  replace = FALSE
)
matched_data <- match.data(m.out)
head(matched_data)

## ---- ① 匹配后 PS 分布（用于直观看重叠） ----
ggplot(matched_data, aes(x = pscore, 
                         color = factor(PA_CAT), 
                         fill = factor(PA_CAT))) +
  geom_density(alpha = 0.3, size = 1.2) +
  labs(title = "匹配后两组的倾向得分分布（1:1卡钳匹配）",
       x = "倾向得分",
       y = "密度",
       color = "PA_CAT",
       fill = "PA_CAT") +
  scale_color_manual(values = c("#1f77b4","#ff7f0e"),
                     labels = c("0 = 非高/低","1 = 高")) +
  scale_fill_manual(values = c("#1f77b4","#ff7f0e"),
                    labels = c("0 = 非高/低","1 = 高")) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())

## ---- ② 匹配前后协变量分布（TableOne）----
## 说明：可同时查看均值/比例以及（可选的）SMD 指标
tab_before <- CreateTableOne(vars = sig_vars, strata = "PA_CAT", data = data)
tab_after  <- CreateTableOne(vars = sig_vars, strata = "PA_CAT", data = matched_data)
tab_after

## ---- ③ 平衡性诊断（cobalt + Love plot） ----
library(cobalt)

# bal.tab：汇总匹配前/后 SMD
# 关键参数：
#   un=TRUE           : 同时给出未调整（匹配前）的不平衡
#   s.d.denom="pooled": SMD 的标准化分母用合并标准差
#   binary="std"      : 二分类变量用标准化口径计算 SMD
bal_match <- bal.tab(
  m.out,
  un = TRUE,
  s.d.denom = "pooled",
  binary    = "std"
)

print(bal_match)

# love.plot：统一可视化 SMD（常用阈值 0.1）
love.plot(
  bal_match,
  stats     = "mean.diffs",
  drop.distance = TRUE,
  abs       = TRUE,
  threshold = 0.1,
  var.order = "unadjusted",
  colors    = c("red", "blue")
)

## ---- ④ 匹配后的 ATE（matched_data 上两组均值差） ----
fit_match <- t.test(lungfun ~ PA_CAT, data = matched_data)

ATE_match <- -(fit_match$estimate[[1]] - fit_match$estimate[[2]])
CI_match  <- -rev(fit_match$conf.int)
p_match   <- fit_match$p.value

cat("\n==== PS Matching: ATE ====\n")
cat("ATE =", round(ATE_match, 3), "\n")
cat("95% CI = [", round(CI_match[1], 3), ",", round(CI_match[2], 3), "]\n")
cat("P值 =", round(p_match, 4), "\n")



###########################################################
# 3.2 方法②：PS 分层（Stratification）—— WeightIt 的 MMWS
#
# 设置要点：
#   subclass=4    : 按 PS 分位数形成 4 层（Q1~Q4）
#   stabilize=TRUE: 稳定化权重（减少极端权重，提高估计稳定性）
#
# 直观理解：
#   先“分层让组内更可比”，再用分层权重把每层的比较结果合并为总体 ATE
###########################################################
library(WeightIt)
library(cobalt)
library(survey)

# WeightIt：估计 PS + 指定 4 层
w_strat_obj <- weightit(
  formula_ps,
  data      = data,
  method    = "glm",
  estimand  = "ATE",
  stabilize = TRUE,
  subclass  = 4
)

# 从 PS 生成 MMWS 权重（并带 subclass 属性）
w_mmws <- get_w_from_ps(
  ps        = w_strat_obj$ps,
  treat     = data$PA_CAT,
  estimand  = "ATE",
  subclass  = 4,
  stabilize = TRUE
)

# 写回：PS（分层用）、分层标签、分层权重
# 注：这里保留 data$pscore（第2步 glm 的 PS），便于 3.3 的 PS 校正复用
data$pscore_strat <- w_strat_obj$ps
data$ps_strata <- factor(attr(w_mmws, "subclass")[,1],
                         levels = 1:4, labels = paste0("Q", 1:4))
data$w_strat   <- w_mmws[1:nrow(data)]  


###########################################################
# 3.2.1 分层平衡性诊断（在估计 ATE 之前）
#
# A) 分层内：每一层内部对比（更强调“可比性”）
# B) 整体：用分层权重后整体 SMD 是否改善（更强调“合并后的总体平衡”）
###########################################################
library(tableone)

cat("\n================ PS 分层后的协变量平衡性（Table1 + SMD）================\n")

strata_levels <- levels(data$ps_strata)
smd_list <- list()

for (s in strata_levels) {
  cat("\n---- 分层", s, "内的协变量平衡情况 ----\n")
  
  tab_s <- CreateTableOne(
    vars       = sig_vars,
    strata     = "PA_CAT",
    data       = subset(data, ps_strata == s),
    factorVars = intersect(sig_vars, categorical_vars)
  )
  
  print(tab_s, showAllLevels = TRUE, smd = TRUE)
  smd_list[[s]] <- ExtractSmd(tab_s)
}

smd_strata <- do.call(cbind, smd_list)
colnames(smd_strata) <- strata_levels

cat("\n==== 各分层的 SMD 表 ====\n")
print(round(smd_strata, 3))


cat("\n================ MMWS 分层权重：加权前后整体平衡性（cobalt）================\n")

bal_mmws <- bal.tab(
  formula_ps,
  data      = data,
  weights   = data$w_strat,
  method    = "weighting",
  un        = TRUE,
  s.d.denom = "pooled",
  binary    = "std"
)

print(bal_mmws)

love.plot(
  bal_mmws,
  stats     = "mean.diffs",
  abs       = TRUE,
  threshold = 0.1,
  var.order = "unadjusted",
  colors    = c("red", "blue")
)

###########################################################
# 3.2.2 分层加权 ATE（survey）
#
# 加权后拟合：lungfun ~ PA_CAT
# 解释：PA_CAT 的回归系数即（A=1 相对 A=0 的）ATE 估计
###########################################################
design_strat <- svydesign(ids = ~1, weights = ~w_strat, data = data)
fit_strat    <- svyglm(lungfun ~ PA_CAT, design = design_strat)

coef_tab <- summary(fit_strat)$coef
rn <- rownames(coef_tab)
coef_name <- rn[grepl("^PA_CAT", rn)][1]

ATE_strat <- coef_tab[coef_name, "Estimate"]
SE_strat  <- coef_tab[coef_name, "Std. Error"]
CI_strat  <- c(ATE_strat - 1.96*SE_strat, ATE_strat + 1.96*SE_strat)
p_strat   <- coef_tab[coef_name, "Pr(>|t|)"]



###########################################################
# 3.3 方法③：PS 校正（PS Adjustment）
#
# 基本思路：
#   在结局回归模型中同时调整暴露和 PS：
#   lungfun ~ PA_CAT + pscore
#
# 说明：
#   这是一种“模型校正”策略，依赖线性模型形式假设（含线性项的充分性）
###########################################################
fit_adj <- lm(lungfun ~ PA_CAT + pscore, data=data)

est_adj <- summary(fit_adj)$coef["PA_CAT", ]

ATE_adj <- est_adj["Estimate"]
SE_adj  <- est_adj["Std. Error"]
CI_adj  <- c(ATE_adj - 1.96*SE_adj,
             ATE_adj + 1.96*SE_adj)
p_adj   <- est_adj["Pr(>|t|)"]



###########################################################
# 3.4 方法④：IPTW 加权（ATE 权重）—— WeightIt + 平衡性 + 结局模型
#
# 核心步骤：
#   1）估计 PS，并生成 ATE 目标下的 IPTW 权重
#   2）检查加权前后协变量 SMD（是否达到“加权后可比”）
#   3）在加权样本上拟合结局模型：lungfun ~ PA_CAT
#
# 关键参数：
#   estimand="ATE"：目标是总体平均处理效应（而非 ATT 等）
###########################################################
library(WeightIt)
library(cobalt)

w_ipw_obj <- weightit(
  formula_ps,
  data      = data,
  method    = "glm",
  estimand  = "ATE"
)

data$pscore_ipw <- w_ipw_obj$ps 
data$w_ate      <- w_ipw_obj$weights


cat("\n==== IPTW 加权后的协变量平衡性（cobalt）====\n")

bal_iptw <- bal.tab(
  w_ipw_obj,
  un        = TRUE,
  s.d.denom = "pooled",
  binary    = "std"
)

print(bal_iptw)

love.plot(
  bal_iptw,
  stats     = "mean.diffs",
  abs       = TRUE,
  threshold = 0.1,
  var.order = "unadjusted",
  colors    = c("red", "blue")
)

# 加权结局模型（WeightIt 提供的便捷接口）
fit_ipw <- lm_weightit(
  lungfun ~ PA_CAT,
  data     = data,
  weightit = w_ipw_obj
)

coef_tab <- coef(summary(fit_ipw))
rn <- rownames(coef_tab)
coef_name <- rn[grepl("^PA_CAT", rn)][1]

ATE_ipw <- coef_tab[coef_name, "Estimate"]
SE_ipw  <- coef_tab[coef_name, "Std. Error"]
CI_ipw  <- confint(fit_ipw, parm = coef_name, level = 0.95)
p_ipw   <- coef_tab[coef_name, "Pr(>|z|)"]



###########################################################
# 4. ATE 汇总表（未调整 + 4 种 PS 方法）
#
# 输出口径统一：
#   ATE：A=1 相对 A=0
#   CI ：95% 置信区间
#   P  ：对应检验的 P 值
###########################################################

result_summary <- data.frame(
  Method = c("未调整", "PS 匹配", "PS 分层", "PS 校正", "IPTW"),
  
  ATE = c(
    ATE_unadj,
    ATE_match,
    ATE_strat,
    ATE_adj,
    ATE_ipw
  ),
  
  CI_lower = c(
    -t_unadj$conf.int[2],
    CI_match[1],
    CI_strat[1],
    CI_adj[1],
    CI_ipw[1]
  ),
  
  CI_upper = c(
    -t_unadj$conf.int[1],
    CI_match[2],
    CI_strat[2],
    CI_adj[2],
    CI_ipw[2]
  ),
  
  P_value = c(
    t_unadj$p.value,
    p_match,
    p_strat,
    p_adj,
    p_ipw
  )
)

# 结果展示：数值统一保留 4 位小数
result_summary_round <- data.frame(
  Method = result_summary$Method,
  round(result_summary[ , -1], 4)
)

cat("\n\n==================== 五种方法 ATE 汇总表 ====================\n")
print(result_summary_round)
