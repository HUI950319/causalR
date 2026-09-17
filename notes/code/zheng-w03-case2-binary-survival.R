###########################################################
#  二分类暴露（PA_CAT）对生存结局（lungfun, lungfun_time）的因果效应估计 —— HR
#  Matching / Stratification / Adjustment / IPTW
#
#  暴露：PA_CAT（1=高, 0=低/非高）
#  结局：lungfun（1=事件发生，0=未发生或删失）
#        lungfun_time（随访时间，月）
#  目标：估计 ATE 的风险比 HR（高 vs 低）
#
#  总体流程：
#   -0. 依赖包检查与安装
#    0. 数据读取与变量预处理（类型、取值检查）
#    1. 单因素 Cox（教学演示）筛选进入 PS 模型的协变量
#    2. 倾向得分模型（Logistic）+ PS 分布图（检查支持性/重叠）
#    3.0 未调整 Cox：HR（基线对照）
#    3.1 PS 匹配 + 平衡性 + 匹配后 Cox（cluster 稳健方差）
#    3.2 PS 分层（MMWS）+ 平衡性 + 加权 Cox（svycoxph）
#    3.3 PS 校正：Cox 中调整 PS（样条以放松线性假设）
#    3.4 IPTW（ATE）+ 平衡性 + 加权 Cox（svycoxph）
#    4. 汇总：五种方法 HR / 95%CI / P 值
###########################################################

##---------------------------------------------------------
## -1. 依赖包：如未安装则先安装（一次性）
##---------------------------------------------------------
pkgs_needed <- c(
  "survival", # Cox 模型 / Surv()
  "ggplot2",  # 作图
  "MatchIt",  # PS 匹配
  "cobalt",   # 平衡性诊断（SMD、Love plot）
  "WeightIt", # PS 分层权重 / IPTW 权重
  "survey",   # 加权 Cox（svycoxph）
  "splines"   # ns() 样条函数（base 推荐包，通常无需安装）
)

pkgs_to_install <- pkgs_needed[!vapply(pkgs_needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(pkgs_to_install) > 0) {
  install.packages(pkgs_to_install, dependencies = TRUE)
}

###########################################################
## 0. 读取数据 & 变量预处理
###########################################################
rm(list = ls())
setwd("C:/Users/77489/Desktop/倾向得分培训班+训练营_25.12.23/录课")

data <- read.csv("data.csv")

# 候选协变量（潜在混杂因素）
names(data)
covariates <- c("sex","hypertension","diabetes","heart_disease","residence","marital",
                "smoking","age","education_level","cardio_metabolic",
                "alcohol","bmi")

# 需要按分类变量处理的协变量（转为 factor，便于哑变量展开与分层比较）
categorical_vars <- c("sex","hypertension","diabetes","heart_disease","residence","marital",
                      "smoking","education_level","cardio_metabolic",
                      "alcohol")

# 统一把分类协变量转为 factor（避免被误当作连续变量）
data[categorical_vars] <- lapply(data[categorical_vars], factor)

# 暴露变量检查（应为 0/1）
table(data$PA_CAT)

# 生存结局：确保类型正确（Surv(time, event) 要求 time 数值、event 0/1）
data$lungfun <- as.integer(data$lungfun)
data$lungfun_time <- as.numeric(data$lungfun_time)


###########################################################
# 1. 单因素分析 —— 筛选进入 PS 模型的协变量（教学演示）
#
# 做法：
#   用生存结局做单因素 Cox：Surv(time, event) ~ X
#   P < 0.05 的协变量进入后续 PS 模型（仅用于教学演示）
#
# 说明：
#   - 分类变量：可能产生多行系数（各水平 vs 参照水平）
#              任一水平 P < 0.05 即纳入
#   - 连续变量：单行系数
###########################################################
library(survival)
sig_vars <- c()
result_table <- data.frame()

for (v in covariates) {
  df_tmp <- data.frame(time = data$lungfun_time,
                       event = data$lungfun,
                       x = data[[v]])
  
  fit <- coxph(Surv(time, event) ~ x, data = df_tmp)
  coef_tab <- summary(fit)$coefficients
  
  if (v %in% categorical_vars) {
    pvals <- coef_tab[, "Pr(>|z|)"]
    sig <- ifelse(any(pvals < 0.05), "Yes", "No")
    result_table <- rbind(
      result_table,
      data.frame(variable=v, type="categorical",
                 p_values=paste(round(pvals, 4), collapse=", "),
                 significant=sig)
    )
    if (sig == "Yes") sig_vars <- c(sig_vars, v)
  } else {
    pval <- coef_tab[1, "Pr(>|z|)"]
    sig <- ifelse(pval < 0.05, "Yes", "No")
    result_table <- rbind(
      result_table,
      data.frame(variable=v, type="continuous",
                 p_values=round(pval, 4),
                 significant=sig)
    )
    if (sig == "Yes") sig_vars <- c(sig_vars, v)
  }
}

cat("\n==== 单因素 Cox 结果（用于筛选 PS 变量）====\n")
print(result_table)
cat("\n进入 PS 模型的协变量：\n")
print(sig_vars)


###########################################################
# 2. 倾向得分模型（Logistic）
#
# 目标：
#   PS = P(PA_CAT=1 | X)
#
# 说明：
#   画 PS 分布密度图，检查两组重叠（支持性/positivity）
###########################################################
library(ggplot2)

formula_ps <- as.formula(paste("PA_CAT ~", paste(sig_vars, collapse = " + ")))
ps_model <- glm(formula_ps, data = data, family = binomial())
data$pscore <- predict(ps_model, type = "response")

ggplot(data, aes(pscore, color=factor(PA_CAT), fill=factor(PA_CAT))) +
  geom_density(alpha=0.3, linewidth=1.0) +
  labs(title="倾向得分分布（原始）", x="PS = P(PA_CAT=1|X)", y="密度") +
  theme_minimal(base_size = 13) +
  theme(legend.title = element_blank())

###########################################################
# 3.0 未调整效应（Cox）：HR
#
# 模型：
#   Surv(time, event) ~ PA_CAT
#
# 输出：
#   HR = exp(beta)，并给出 95% CI 与 P 值
###########################################################
fit_unadj <- coxph(Surv(lungfun_time, lungfun) ~ PA_CAT, data = data)
sum_unadj <- summary(fit_unadj)

HR_unadj <- exp(coef(fit_unadj)["PA_CAT"])
CI_unadj <- exp(confint(fit_unadj)["PA_CAT", ])
p_unadj  <- sum_unadj$coefficients["PA_CAT", "Pr(>|z|)"]

###########################################################
# 3.1 方法①：PS Matching（MatchIt） + Cox（cluster=subclass）
#
# 匹配设置要点：
#   method="nearest"  : 最近邻匹配
#   distance="logit"  : 距离基于 logit(PS)
#   caliper=0.2       : 卡钳宽度（限制可匹配的距离，避免“硬配”）
#   replace=FALSE     : 不放回
#
# 匹配后 Cox：
#   - 使用 cluster(subclass) 计算稳健方差（考虑匹配对子类相关）
#   - weights=weights 使用 MatchIt 生成的匹配权重
###########################################################
library(MatchIt)
library(cobalt)

m.out <- tryCatch(
  matchit(formula_ps, data = data,
          method = "nearest", distance = "logit",
          caliper = 0.2, replace = FALSE,
          estimand = "ATE"),
  error = function(e) {
    message("提示：你的 MatchIt 版本可能不支持 estimand='ATE'，将按默认设置运行（matching 通常更接近 ATT/重叠子人群）。")
    matchit(formula_ps, data = data,
            method = "nearest", distance = "logit",
            caliper = 0.2, replace = FALSE)
  }
)

matched_data <- match.data(m.out)

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

## ---- ② 匹配前后协变量分布（TableOne） ----
library(tableone)
tab_before <- CreateTableOne(vars = sig_vars, strata = "PA_CAT", data = data)
tab_after  <- CreateTableOne(vars = sig_vars, strata = "PA_CAT", data = matched_data)
tab_after

# 平衡性诊断：匹配前后 SMD（常用阈值 |SMD|<0.1）
bal_match <- bal.tab(m.out, un = TRUE, s.d.denom = "pooled", binary = "std")
print(bal_match)
love.plot(bal_match, stats="mean.diffs", abs=TRUE, threshold=0.1,
          var.order="unadjusted", drop.distance=TRUE)

fit_match <- coxph(Surv(lungfun_time, lungfun) ~ PA_CAT + cluster(subclass),
                   data = matched_data, weights = weights)

sum_match <- summary(fit_match)
HR_match <- exp(coef(fit_match)["PA_CAT"])
CI_match <- exp(confint(fit_match)["PA_CAT", ])
p_match  <- sum_match$coefficients["PA_CAT", "Pr(>|z|)"]

###########################################################
# 3.2 方法②：PS 分层（MMWS） + 加权 Cox（svycoxph）
#
# 分层设置：
#   subclass=4     : 将样本按 PS 分位数分为 4 层（Q1~Q4）
#   stabilize=TRUE : 稳定化权重，减少极端权重带来的方差膨胀
#
# 分层权重（MMWS）：
#   在“分层 × 暴露组”内为常数，使加权后协变量总体更平衡
###########################################################
library(WeightIt)
library(survey)

w_strat_obj <- weightit(
  formula_ps,
  data      = data,
  method    = "glm",
  estimand  = "ATE",
  stabilize = TRUE,
  subclass  = 4
)

w_mmws <- get_w_from_ps(
  ps        = w_strat_obj$ps,
  treat     = data$PA_CAT,
  estimand  = "ATE",
  subclass  = 4,
  stabilize = TRUE
)

data$w_strat <- w_mmws[1:nrow(data)]

# 平衡性：分层权重加权前后 SMD
bal_mmws <- bal.tab(formula_ps, data = data, weights = data$w_strat,
                    method="weighting", un=TRUE, s.d.denom="pooled", binary="std")
print(bal_mmws)
love.plot(bal_mmws, stats="mean.diffs", abs=TRUE, threshold=0.1,
          var.order="unadjusted", drop.distance=TRUE)

# 加权 Cox（ATE）
design_strat <- svydesign(ids = ~1, weights = ~w_strat, data = data)
fit_strat <- svycoxph(Surv(lungfun_time, lungfun) ~ PA_CAT, design = design_strat)

coef_strat <- coef(summary(fit_strat))
HR_strat <- exp(coef(fit_strat)["PA_CAT"])
SE_strat <- coef_strat["PA_CAT", "se(coef)"]
CI_strat <- exp(c(log(HR_strat) - 1.96*SE_strat, log(HR_strat) + 1.96*SE_strat))
p_strat  <- coef_strat["PA_CAT", "Pr(>|z|)"]

###########################################################
# 3.3 方法③：PS 校正（Cox 调整 PS）
#
# 模型：
#   Surv(time, event) ~ PA_CAT + f(PS)
#
# 说明：
#   这里用自然样条 ns(pscore, df=3) 放松“PS 线性”假设
#   df 为样条自由度：越大越灵活，但也更易过拟合（教学中 3 常用）
###########################################################
library(splines)

fit_adj <- coxph(Surv(lungfun_time, lungfun) ~ PA_CAT + ns(pscore, df = 3), data = data)
sum_adj <- summary(fit_adj)

HR_adj <- exp(coef(fit_adj)["PA_CAT"])
CI_adj <- exp(confint(fit_adj)["PA_CAT", ])
p_adj  <- sum_adj$coefficients["PA_CAT", "Pr(>|z|)"]

###########################################################
# 3.4 方法④：IPTW（ATE） + 加权 Cox（svycoxph）
#
# IPTW 权重（ATE）：
#   目标是让“加权后的总体”在协变量分布上类似随机试验
#
# 关键检查：
#   - Love plot 看加权后 SMD 是否明显下降
#   - 如出现极端权重，可考虑 stabilize、截尾/修剪（此处仅演示）
###########################################################
w_ipw_obj <- weightit(
  formula_ps,
  data     = data,
  method   = "glm",
  estimand = "ATE"
)

data$w_ate <- w_ipw_obj$weights

# 平衡性：IPTW 加权前后 SMD
bal_iptw <- bal.tab(w_ipw_obj, un=TRUE, s.d.denom="pooled", binary="std")
print(bal_iptw)
love.plot(bal_iptw, stats="mean.diffs", abs=TRUE, threshold=0.1,
          var.order="unadjusted", drop.distance=TRUE)

# 加权 Cox（ATE）
design_ipw <- svydesign(ids = ~1, weights = ~w_ate, data = data)
fit_ipw <- svycoxph(Surv(lungfun_time, lungfun) ~ PA_CAT, design = design_ipw)

coef_ipw <- coef(summary(fit_ipw))
HR_ipw <- exp(coef(fit_ipw)["PA_CAT"])
SE_ipw <- coef_ipw["PA_CAT", "se(coef)"]
CI_ipw <- exp(c(log(HR_ipw) - 1.96*SE_ipw, log(HR_ipw) + 1.96*SE_ipw))
p_ipw  <- coef_ipw["PA_CAT", "Pr(>|z|)"]

###########################################################
# 4. 汇总表：HR / 95% CI / P
#
# 输出口径：
#   HR > 1 ：提示 PA_CAT=1 组风险更高
#   HR < 1 ：提示 PA_CAT=1 组风险更低
###########################################################
result_summary <- data.frame(
  Method   = c("未调整", "PS 匹配", "PS 分层(MMWS)", "PS 校正", "IPTW(ATE)"),
  HR       = c(HR_unadj, HR_match, HR_strat, HR_adj, HR_ipw),
  CI_lower = c(CI_unadj[1], CI_match[1], CI_strat[1], CI_adj[1], CI_ipw[1]),
  CI_upper = c(CI_unadj[2], CI_match[2], CI_strat[2], CI_adj[2], CI_ipw[2]),
  P_value  = c(p_unadj, p_match, p_strat, p_adj, p_ipw)
)

result_summary_round <- data.frame(
  Method = result_summary$Method,
  round(result_summary[, -1], 4)
)

cat("\n\n==================== 五种方法 HR 结果汇总 ====================\n")
print(result_summary_round)
