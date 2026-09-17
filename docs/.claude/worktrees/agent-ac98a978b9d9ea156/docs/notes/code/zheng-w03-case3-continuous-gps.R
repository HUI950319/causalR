###########################################################
# 连续型处理（PA）对二分类结局（lungfun 0/1）的因果效应估计 —— RD(概率差)口径
# GPS-Matching / GPS-Stratification / GPS-Adjustment / GPS-Weighting
#
# 目标 estimand：
#   ATE(RD) = E[Y(a+1) - Y(a)]，即 PA 每 +1 单位，lungfun 发病概率的平均变化
#
# 口径说明：
#   - 统一用“线性概率模型”(Gaussian) 输出 RD，便于解释与不同方法间可比
#   - RD 可理解为“风险差/概率差”，不是 OR/HR
#
# 平衡性诊断说明（连续型暴露）：
#   - 连续暴露下通常用“协变量与 PA 的相关系数”作为平衡性指标（越接近 0 越好）
#   - 这里不能用二分类暴露常见的 SMD（标准化均值差）口径
#   - cobalt::bal.tab() 可以输出 correlations，但目前不支持“连续暴露 + subclass(分层)”
#     来直接给出总体 love plot；因此“总体（分层后）love plot”用 ggplot 自行实现
###########################################################

##---------------------------------------------------------
## -1. 依赖包：如未安装则先安装（一次性）
##---------------------------------------------------------
pkgs_needed <- c(
  "CausalGPS",     # 连续暴露的 GPS 估计与 counter weight
  "SuperLearner",  # estimate_gps 的 sl_lib（此处用轻量 SL.glm）
  "ggplot2",       # 作图
  "cobalt",        # 平衡性诊断（相关性 + love plot）
  "survey",        # svydesign / svyglm
  "WeightIt"       # 连续暴露加权（ps 方法）
)

pkgs_to_install <- pkgs_needed[!vapply(pkgs_needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(pkgs_to_install) > 0) {
  install.packages(pkgs_to_install, dependencies = TRUE)
}

library(CausalGPS)
library(SuperLearner)
library(ggplot2)
library(cobalt)
library(survey)
library(WeightIt)

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

# 需要按分类变量处理的协变量（转为 factor，便于哑变量展开/相关性计算更合理）
categorical_vars <- c("sex","hypertension","diabetes","heart_disease","residence","marital",
                      "smoking","education_level","cardio_metabolic",
                      "alcohol")

# 统一把分类协变量转为 factor（避免被误当作连续变量）
data[categorical_vars] <- lapply(data[categorical_vars], factor)

# 连续处理：PA
# 结局：lungfun (0/1)
data$id <- seq_len(nrow(data))

###########################################################
## 1. 单因素分析 —— 筛选进入 GPS 模型的协变量（教学演示）
##    做法：lungfun ~ 协变量 的单因素线性概率模型；P < 0.05 进入 GPS
##
## 说明：
##   - 这是教学演示用的“筛选流程”，真实研究更推荐基于领域知识/DAG 决定混杂因素
##   - 因子变量：任一水平显著即纳入
##   - 连续变量：看对应系数的 P 值
###########################################################
sig_vars <- c()
result_table <- data.frame()

for (v in covariates) {
  df_tmp <- data.frame(lungfun = data$lungfun, x = data[[v]])
  fit <- glm(lungfun ~ x, data = df_tmp, family = "binomial")
  coef_table <- summary(fit)$coefficients
  
  if (v %in% categorical_vars) {
    # 因子变量：看各水平（相对参考组）的 p 值，任一<0.05则入模
    pvals <- coef_table[-1, 4]
    sig <- ifelse(any(pvals < 0.05), "Yes", "No")
    
    result_table <- rbind(
      result_table,
      data.frame(variable=v, type="categorical",
                 p_values=paste(round(pvals,4), collapse=", "),
                 significant=sig)
    )
    if (sig == "Yes") sig_vars <- c(sig_vars, v)
    
  } else {
    # 连续变量
    pval <- coef_table[2, 4]
    sig  <- ifelse(pval < 0.05, "Yes", "No")
    
    result_table <- rbind(
      result_table,
      data.frame(variable=v, type="continuous",
                 p_values=round(pval,4), significant=sig)
    )
    if (sig == "Yes") sig_vars <- c(sig_vars, v)
  }
}

cat("\n==== 单因素分析结果 ====\n")
print(result_table)

cat("\n进入 GPS 模型的协变量：\n")
print(sig_vars)

###########################################################
## 2. GPS 模型（连续处理的“倾向得分”：条件密度 f(PA|X)）
##
## 说明：
##   - 连续暴露的 GPS 可理解为“给定 X 的 PA 条件密度”
##   - estimate_gps() 输出 gps（conditional density），后续用于匹配型 counter weights
###########################################################
formula_gps <- as.formula(paste("PA ~", paste(sig_vars, collapse=" + ")))

gps_obj <- estimate_gps(
  .data       = data,
  .formula    = formula_gps,
  gps_density = "kernel",
  sl_lib      = c("SL.glm")
)

data$gps <- gps_obj$.data$gps

# 一个更常用的“balancing score”：e_hat = E[PA|X],用于后续GPS分层法
data$e_pa_hat <-  gps_obj$.data$e_gps_pred

###########################################################
## 3.0 未调整效应（RD口径）：lungfun ~ PA （线性概率模型）
##
## 解释：
##   - PA 的回归系数即 PA 每增加 1 单位，发病概率的平均变化（RD）
###########################################################
fit_unadj <- lm(lungfun ~ PA, data = data)
est_unadj <- coef(summary(fit_unadj))["PA", ]

ATE_unadj <- est_unadj["Estimate"]
SE_unadj  <- est_unadj["Std. Error"]
CI_unadj  <- c(ATE_unadj - 1.96*SE_unadj, ATE_unadj + 1.96*SE_unadj)
p_unadj   <- est_unadj["Pr(>|t|)"]

###########################################################
## 3.1 方法①：GPS 匹配（CausalGPS matching）
##     输出 RD：在匹配计数权重(counter_weight)下拟合线性概率模型
##
## 关键参数（compute_counter_weight）：
##   - ci_appr="matching" : 使用匹配型 counter weights
##   - delta_n            : 匹配邻域/平滑强度（越大通常越“宽松”，权重更平滑）
##   - dist_measure="l1"  : 距离度量（L1 更直观）
##   - scale              : 距离缩放（影响匹配邻域大小，越大通常权重更平滑）
##
## 平衡性诊断（连续暴露）：
##   - 用“协变量与 PA 的相关系数”而非 SMD（SMD 主要用于二分类/多分类处理）
###########################################################
cw_match <- compute_counter_weight(
  gps_obj  = gps_obj,
  ci_appr  = "matching",
  delta_n  = 4,
  dist_measure = "l1",
  scale    = 1
)

data$w_match <- cw_match$.data$counter_weight

# 权重分布检查：关注极端权重（可提示是否需要调整 delta_n/scale 或做截尾）
#     - 这些阈值不是硬标准，只是帮助你快速判断是否存在严重的极端权重
#      - 一般来说，P(w>10) 非常小（接近 0）会更放心；如果已经明显>0，建议进一步处理
#         - 常见操作：delta_n 从 2.5 → 4 → 6（逐步加）
#         - 常见操作：scale 从 0.5 → 1 → 1.5（逐步加）
#      调参策略：
#         - 每次只改一个参数（例如先调 delta_n），观察权重分布 + 平衡性是否同时改善
#         - 目标是：权重不极端 + 平衡性可接受（|cor| 接近 0，且多数变量<0.1）
summary(data$w_match)
quantile(data$w_match, probs = c(.9,.95,.99,.995,.999), na.rm=TRUE)
cat("max w_match =", max(data$w_match, na.rm=TRUE), "\n")
cat("P(w_match>10) =", mean(data$w_match>10), "\n")
cat("P(w_match>20) =", mean(data$w_match>20), "\n")
cat("P(w_match>50) =", mean(data$w_match>50), "\n")

# 平衡性：协变量与 PA 的相关性（越接近0越好）
bal_match <- bal.tab(
  PA ~ .,
  data    = data[, c("PA", sig_vars), drop = FALSE],
  weights = data$w_match,
  method  = "weighting",
  un      = TRUE,
  stats   = "correlations"
)
print(bal_match)

love.plot(
  bal_match,
  stats     = "correlations",
  abs       = TRUE,
  threshold = 0.1,
  var.order = "unadjusted"
)

# RD 结局模型（Gaussian -> 线性概率模型）
design_match <- svydesign(ids = ~1, weights = ~w_match, data = data)
fit_match    <- svyglm(lungfun ~ PA, design = design_match, family = gaussian())

est_match <- coef(summary(fit_match))["PA", ]
ATE_match <- est_match["Estimate"]
SE_match  <- est_match["Std. Error"]
CI_match  <- c(ATE_match - 1.96*SE_match, ATE_match + 1.96*SE_match)
p_match   <- est_match["Pr(>|t|)"]

###########################################################
## 3.2 方法②：GPS 分层（Stratification）
##     用 e_pa_hat 分 4 层；层内 RD（lm），再按样本占比加权合并
##
## 分层的目的：
##   - 让层内的 PA 与协变量关系更接近“随机/弱相关”，提高可比性
##
## 平衡性诊断（连续暴露）：
##   - 依然使用“相关系数”作为指标（不使用 SMD）
##   - cobalt 的 bal.tab() 可以做“每一层内部”的相关性平衡诊断
##   - 但 cobalt 目前不支持：连续暴露 + subclass(分层) 来直接输出“分层后的总体 love plot”
##     因此：总体（分层后）相关性用 ggplot 自行实现
##
## 总体相关性的原理（教学口径）：
##   - 用“层内相关性”的总体汇总：等价于控制 strata 后的相关性
##   - 实现方式：对 PA 与每个协变量分别对 strata 残差化，再计算残差相关
###########################################################
# gps 是 密度值 f(PA∣X)。它不是“暴露大小”的排序量，受带宽、分布形状影响很大，
# 而且不同个体的 gps 大小不一定对应“暴露高/低”或“可比性强/弱”。
# 用它分层，层的含义不直观，还可能把“密度高（更常见）/密度低（更罕见）”的人混在一起
s_breaks <- quantile(data$e_pa_hat, probs = seq(0, 1, 0.25), na.rm = TRUE)
s_breaks[1] <- s_breaks[1] - 1e-8

data$strata <- cut(
  data$e_pa_hat,
  breaks = s_breaks,
  include.lowest = TRUE,
  labels = paste0("Q", 1:4)
)
table(data$strata)

##-----------------------------
## 3.2.1 分层内平衡性（逐层相关性 + love plot）
##-----------------------------
strata_levels <- levels(data$strata)
bal_list <- list()

cat("\n================ 分层内平衡性（correlations）================\n")

for (s in strata_levels) {
  cat("\n---- 分层", s, "----\n")
  ds <- droplevels(subset(data, strata == s))
  
  # 只在该层内做平衡性：剔除“该层内无变异”的变量（否则无法计算相关性/模型矩阵）
  vars_s <- sig_vars
  one_level_factors <- vars_s[sapply(vars_s, function(v) {
    is.factor(ds[[v]]) && nlevels(ds[[v]]) < 2
  })]
  
  if (length(one_level_factors) > 0) {
    cat("提示：该层以下因子仅 1 个水平，已从 bal.tab 中剔除：",
        paste(one_level_factors, collapse = ", "), "\n")
    vars_s <- setdiff(vars_s, one_level_factors)
  }
  
  bal_s <- bal.tab(
    PA ~ .,
    data  = ds[, c("PA", vars_s), drop = FALSE],
    un    = TRUE,
    stats = "correlations"
  ) # cobalt 目前明确限制：continuous treatments 不兼容 subclasses,因此后面用ggplot来画
  print(bal_s)
  bal_list[[s]] <- bal_s
}


##-----------------------------
## 3.2.2 总体（分层后）相关性 love plot：ggplot 实现
##-----------------------------
# 说明：
#   - 目标是得到“层内相关性”的总体汇总指标
#   - 做法：对 PA 与协变量分别对 strata 做残差化（去掉层间差异），再相关
#   - 这个相关系数可视为控制 strata 后的相关/层内 pooled correlation

X <- model.matrix(
  as.formula(paste("~", paste(sig_vars, collapse = " + "))),
  data = data
)[, -1, drop = FALSE]

A <- data$PA
S <- data$strata

# 未分层：cor(PA, X_j)
cor_un <- apply(X, 2, function(x) cor(A, x, use = "pairwise.complete.obs"))

# 分层后总体：cor(resid(PA|strata), resid(X_j|strata))
A_res <- resid(lm(A ~ S))
cor_strat <- apply(X, 2, function(x) {
  x_res <- resid(lm(x ~ S))
  cor(A_res, x_res, use = "pairwise.complete.obs")
})

df_plot2 <- rbind(
  data.frame(var = names(cor_un), cor = cor_un, group = "Overall"),
  data.frame(var = names(cor_strat), cor = cor_strat, group = "Stratified (within-strata)")
)
df_plot2$abs_cor <- abs(df_plot2$cor)

ggplot(df_plot2, aes(x = abs_cor, y = reorder(var, abs_cor), color = group)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0.1, linetype = 2) +
  labs(x = "|Correlation with PA|", y = NULL, color = NULL) +
  theme_bw()

##-----------------------------
## 3.2.3 分层合并效应（教学版）
##-----------------------------
beta_h <- se_h <- w_h <- numeric(length(strata_levels))
N_total <- nrow(data)

for (j in seq_along(strata_levels)) {
  s  <- strata_levels[j]
  ds <- subset(data, strata == s)
  
  fit_s <- lm(lungfun ~ PA, data = ds)
  est_s <- coef(summary(fit_s))["PA", ]
  
  beta_h[j] <- est_s["Estimate"]
  se_h[j]   <- est_s["Std. Error"]
  w_h[j]    <- nrow(ds) / N_total
}

ATE_strat <- sum(w_h * beta_h)
SE_strat  <- sqrt(sum((w_h^2) * (se_h^2)))  # 教学版合并SE（近似）
CI_strat  <- c(ATE_strat - 1.96*SE_strat, ATE_strat + 1.96*SE_strat)
p_strat   <- 2 * (1 - pnorm(abs(ATE_strat / SE_strat)))

###########################################################
## 3.3 方法③：GPS 校正（Adjustment）—— RD口径
##     常用做法：lungfun ~ PA + gps + gps^2 + PA:gps （加入非线性与交互）
##
## 说明：
##   - 这里把 gps 作为调整项，并加入二次项与交互，提升函数形式灵活性（教学示例）
###########################################################
fit_adj <- lm(lungfun ~ PA + gps + I(gps^2) + PA:gps, data = data)
est_adj <- coef(summary(fit_adj))["PA", ]

ATE_adj <- est_adj["Estimate"]
SE_adj  <- est_adj["Std. Error"]
CI_adj  <- c(ATE_adj - 1.96*SE_adj, ATE_adj + 1.96*SE_adj)
p_adj   <- est_adj["Pr(>|t|)"]

###########################################################
## 3.4 方法④：GPS 加权（IPTW 类比）—— RD口径
##     用 WeightIt 生成连续处理权重；加权后用 svyglm(Gaussian)估计 RD
##
## 说明：
##   - 连续暴露下的加权思想与 IPTW 类似：加权后让 PA 与协变量更接近“独立/弱相关”
##   - 平衡性同样用 correlations 口径
###########################################################
w_wt_obj <- weightit(
  formula_gps,
  data     = data[, c("PA", sig_vars), drop = FALSE],
  method   = "ps",
  estimand = "ATE"
)

data$w_ipw <- w_wt_obj$weights

# 平衡性诊断：correlations
bal_ipw <- bal.tab(
  PA ~ .,
  data    = data[, c("PA", sig_vars), drop = FALSE],
  weights = data$w_ipw,
  method  = "weighting",
  un      = TRUE,
  stats   = "correlations"
)
print(bal_ipw)

love.plot(
  bal_ipw,
  stats     = "correlations",
  abs       = TRUE,
  threshold = 0.1,
  var.order = "unadjusted"
)

# RD 加权结局模型（Gaussian -> 线性概率模型）
design_ipw <- svydesign(ids = ~1, weights = ~w_ipw, data = data)
fit_ipw    <- svyglm(lungfun ~ PA, design = design_ipw, family = gaussian())

est_ipw <- coef(summary(fit_ipw))["PA", ]
ATE_ipw <- est_ipw["Estimate"]
SE_ipw  <- est_ipw["Std. Error"]
CI_ipw  <- c(ATE_ipw - 1.96*SE_ipw, ATE_ipw + 1.96*SE_ipw)
p_ipw   <- est_ipw["Pr(>|t|)"]

###########################################################
## 4. 汇总：未调整 + 4 种 GPS 方法（统一为 RD）
###########################################################
result_summary <- data.frame(
  Method  = c("未调整(RD)", "GPS 匹配(RD)", "GPS 分层(RD)", "GPS 校正(RD)", "GPS 加权(RD)"),
  ATE     = c(ATE_unadj, ATE_match, ATE_strat, ATE_adj, ATE_ipw),
  CI_lower= c(CI_unadj[1], CI_match[1], CI_strat[1], CI_adj[1], CI_ipw[1]),
  CI_upper= c(CI_unadj[2], CI_match[2], CI_strat[2], CI_adj[2], CI_ipw[2]),
  P_value = c(p_unadj, p_match, p_strat, p_adj, p_ipw)
)

result_summary_round <- data.frame(
  Method = result_summary$Method,
  round(result_summary[, -1], 4)
)

cat("\n\n==================== 五种方法结果汇总（统一RD口径）====================\n")
print(result_summary_round)
