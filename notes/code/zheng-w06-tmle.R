# ==========================================================
# R包准备
# ==========================================================
# 
# .libPaths(c("D:/Program Files (x86)/R-4.5.2"))
# .libPaths()

# 安装必要的包
# install.packages("AIPW")
# install.packages("SuperLearner")
# install.packages("tmle")
devtools::install_github("benkeser/survtmle")
#加载必要R包
library(AIPW)
library(SuperLearner)
library(tmle)
library(readr)
library(ctmle)
library(survtmle)
library(ggplot2)

# ==========================================================
# 数据准备
# ==========================================================
set.seed(123)
setwd("D:/继教班2025/20260114TMLE及其衍生方法")
data <- read_csv("data.csv")
colnames(data)

Y <- data$death1  #结局变量
A <- data$RHC     #暴露变量
covars <- c("age1","sex1","edu1","race1","income1"
            ,"Cardiovascular","Psychiatric","Pulmonary")  #协变量
W <- data[, covars]
covariates <- as.matrix(W)
colnames(W)
colSums(is.na(cbind(Y, A, W)))
# ==========================================================
# 方法一 AIPTW
# ==========================================================
#设置SuperLearner库
SL.library <- c("SL.glm", "SL.ranger", "SL.xgboost", "SL.mean")
AIPW_SL <- AIPW$new(
  Y = Y,
  A = A,
  W = covariates,
  Q.SL.library = c("SL.mean", "SL.glm"),  # 结果模型使用的算法
  g.SL.library = c("SL.mean", "SL.glm"),   # 倾向得分模型使用的算法
  k_split = 3,
  verbose = TRUE   # 先开着，看得清楚
)$fit()$summary(g.bound = 0.025) # 
# 假设 g.bound = 0.025  默认值
# 估计的倾向得分设置上下限，防止极端权重。原始倾向得分会被修剪到 [0.025, 0.975] 范围内。

# 打印结果
print(AIPW_SL$result, digits = 2)


# 效应图示

data_AIPW_PLOT<-data.frame(AIPW_SL$result)
data_AIPW_PLOT$Group<-rownames(data_AIPW_PLOT)


# 对于效应值
ggplot(data_AIPW_PLOT, aes(x = Group, y = Estimate, ymin = X95..LCL, ymax =X95..UCL)) +
  geom_pointrange(size = 1, color = "blue") +
  geom_errorbar(width = 0.2) +
  labs(title = "Estimates with 95% Confidence Intervals",
       x = "Group", y = "Risk Estimate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ==========================================================
# 方法二  TMLE
# ==========================================================
SL.library <- c("SL.glm", "SL.ranger", "SL.xgboost", "SL.mean")
tmle_fit <- tmle(Y = as.vector(Y), A = as.vector(A),W = W,
                   Q.SL.library=SL.library,
                   g.SL.library=SL.library,
                   family="gaussian")
# 查看结果
tmle_fit


# 创建整齐格式的数据框
extract_tmle_estimates <- function(tmle_fit) {
  effects <- c("EY0", "EY1", "RR", "OR", "ATE", "ATC", "ATT")
  result_list <- list()
  
  for (effect in effects) {
    if (!is.null(tmle_fit$estimates[[effect]])) {
      ci <- tmle_fit$estimates[[effect]]$CI
      psi <- tmle_fit$estimates[[effect]]$psi
      
      result_list[[effect]] <- data.frame(
        Effect = effect,
        Estimate = psi,
        CI_lower = ci[1],
        CI_upper = ci[2],
        SE = ifelse(!is.null(tmle_fit$estimates[[effect]]$var.psi), 
                    sqrt(tmle_fit$estimates[[effect]]$var.psi), NA)
      )
    }
  }
  
  do.call(rbind, result_list)
}

# 提取并整理结果
tmle_summary <- extract_tmle_estimates(tmle_fit)
print(tmle_summary)


#  效应绘图 

ggplot(tmle_summary, aes(x = Effect, y = Estimate, ymin = CI_lower, ymax =CI_upper)) +
  geom_pointrange(size = 1, color = "blue") +
  geom_errorbar(width = 0.2) +
  labs(title = "Estimates with 95% Confidence Intervals",
       x = "Effects", y = "Estimate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# ==========================================================
# 方法三 DML
# ==========================================================
library(DoubleML)
library(mlr3)
library(mlr3learners)
library(readr)        # 如果还没加载

# 方法1：最推荐、最干净（强烈建议）
library(data.table)

# 把 tibble 直接转成 data.table（最常用写法）
setDT(data)   # ← 这一行就够了！会直接修改原对象

# 之后就可以正常创建了
dml_data <- DoubleMLData$new(
  data   = data,
  y_col  = "death1",     # 建议用列名字符串
  d_cols = "RHC",
  x_cols = covars
)

# 7. 设置机器学习模型（推荐做法）
# 结果模型（g） - 回归问题
# ml_g <- lrn("regr.ranger", 
#            num.trees = 500, 
#            min.node.size = 10,
#            max.depth = 15)

# 倾向得分模型（m） - 二分类问题
# ml_m <- lrn("classif.ranger", 
#            num.trees = 500,
#            min.node.size = 10,
#            max.depth = 15,
#            predict_type = "prob")   # 重要！需要概率预测

# 8. 建立 PLR（部分线性回归）模型 - 最常用的因果推断设定
set.seed(2026)
dml_irm <- DoubleMLIRM$new(
  data = dml_data,
  ml_g = lrn("classif.ranger", predict_type = "prob", num.trees = 500),
  ml_m = lrn("classif.ranger", predict_type = "prob", num.trees = 500),
  score = "ATE",           # 最常见的平均处理效应
  n_folds = 5,
  n_rep = 3
)$fit()$summary()

dml_irm


# ==========================================================
# 方法四 C-TMLE
# ==========================================================
SL.library <- c("SL.glm", "SL.glmnet", "SL.ranger", "SL.mean")
ctmle_fit <- ctmleDiscrete(
  Y = Y,
  A = A,
  W = W,
  family = "binomial",
  SL.library = SL.library
)
summary(ctmle_fit)


# ==========================================================
# 方法五 CV-TMLE
# ==========================================================

# library(tmle)
# library(SuperLearner)

SL.library <- c("SL.glm", "SL.glmnet", "SL.ranger", "SL.mean")
set.seed(2026)
cv_tmle_sl <- tmle(
  Y = Y,                    # 结局变量（连续型结局，对应 family = "binomial"）
  
  A = A,                    # 暴露 / 处理变量（通常为 0/1 的二分类变量）
  
  W = W,                    # 基线协变量矩阵，用于控制混杂（Q 模型和 g 模型）
  
  family = "binomial",      # 指定结局分布：gaussian 表示连续结局
  
  Q.SL.library = SL.library,# 结局模型 Q(Y | A, W) 使用 Super Learner
                            # 通过交叉验证组合多种机器学习算法
  
  g.SL.library = SL.library,# 倾向得分模型 g(A | W) 使用 Super Learner
                            # 用于估计处理分配机制，保证双重稳健性
  
  cvQinit = TRUE,            # 启用 Cross-Validated TMLE：
                            # 对初始结局回归 Q 使用交叉验证（cross-fitting），
                            # 减少机器学习过拟合带来的偏倚

   V.Q=  3 ,                  # Number of cross-validation folds for super learner estimation of Q

   V.g=  3                   # Number of cross-validation folds for super learner estimation of g

)


cv_tmle_sl



# ==========================================================
# 方法六 Survival-TMLE 
# ==========================================================
ftime <- data$Survival_time180
ftype <- data$death1       # 1=死亡, 0=删失
trt   <- data$RHC

W <- data[, ..covars]

set.seed(2026)
fit_surv_tmle <- survtmle(
  ftime = ftime,
  ftype = ftype,
  trt   = trt,
  adjustVars = W,
  SL.trt = "SL.glm",
  
  # 事件风险模型 P(T=t | A, W)
  glm.ftime = "trt + age1 + sex1 + edu1 + race1 + income1 +
               Cardiovascular + Psychiatric + Pulmonary + t",
  
  # 删失模型 P(C=t | W)
  glm.ctime = "age1 + sex1 + edu1 + race1 + income1 +
               Cardiovascular + Psychiatric + Pulmonary + t",
  
  method = "hazard",
  t0 = 180,
  
  # 建議加上這行，方便之後檢查模型
  returnModels = TRUE
)

tp <- timepoints(fit_surv_tmle, times = seq(1, 180, by = 10))
plot(tp)

calculate_effects <- function(fit_surv_tmle) {
  # 提取估计值
  surv0 <- fit_surv_tmle$est[1, 1]
  surv1 <- fit_surv_tmle$est[2, 1]
  
  # 计算风险
  risk0 <- 1 - surv0
  risk1 <- 1 - surv1
  
  # 风险比
  rr <- risk1 / risk0
  log_rr <- log(rr)
  
  # 比值比
  or_risk <- (risk1/surv1) / (risk0/surv0)
  log_or <- log(or_risk)
  
  # 提取方差矩阵
  var_matrix <- fit_surv_tmle$var
  
  # 计算 RR 的标准误（Delta方法）
  var_log_rr <- var_matrix[2,2] / (risk1^2) + 
    var_matrix[1,1] / (risk0^2) - 
    2 * var_matrix[1,2] / (risk1 * risk0)
  se_log_rr <- sqrt(var_log_rr)
  
  # 计算 OR 的标准误（近似）
  se_log_or <- sqrt(1/(risk1*surv1) + 1/(risk0*surv0)) * 
    sqrt(var_matrix[2,2] + var_matrix[1,1])
  
  # 置信区间
  rr_ci <- exp(log_rr + c(-1.96, 1.96) * se_log_rr)
  or_ci <- exp(log_or + c(-1.96, 1.96) * se_log_or)
  
  # 返回结果
  results <- list(
    survival = c(trt0 = surv0, trt1 = surv1),
    risk = c(trt0 = risk0, trt1 = risk1),
    risk_ratio = c(
      estimate = rr,
      lower = rr_ci[1],
      upper = rr_ci[2]
    ),
    odds_ratio = c(
      estimate = or_risk,
      lower = or_ci[1],
      upper = or_ci[2]
    )
  )
  
  return(results)
}

# 使用函数
effects <- calculate_effects(fit_surv_tmle)

cat("=== 生存分析结果 ===\n")
cat(sprintf("对照组生存概率: %.3f\n", effects$survival["trt0.0 1"]))
cat(sprintf("治疗组生存概率: %.3f\n", effects$survival["trt1.1 1"]))
cat(sprintf("\n风险比 (RR): %.3f (95%% CI: %.3f-%.3f)\n",
            effects$risk_ratio["estimate.1 1"],
            effects$risk_ratio["lower"],
            effects$risk_ratio["upper"]))
cat(sprintf("比值比 (OR): %.3f (95%% CI: %.3f-%.3f)\n",
            effects$odds_ratio["estimate.1 1"],
            effects$odds_ratio["lower"],
            effects$odds_ratio["upper"]))
# 解释结果：
# 风险比 (RR): 0.846 (治疗组风险是对照组的84.6%)
# 比值比 (OR): 约 0.723(取决于计算方法)











