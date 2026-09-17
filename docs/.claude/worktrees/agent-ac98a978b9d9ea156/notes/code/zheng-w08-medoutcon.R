# ============================================================================
# ============================================================================
# 课程演示：因果中介分析 (Causal Mediation Analysis) 
# 研究问题：身体活动量(PA_CAT)是否通过抑郁症状(CESD_CAT)影响慢性肺病(Lung)？
# 变量说明：
#   - 结局变量(Y): lungfun      (慢性肺病      0：未患慢性肺病；1：患有慢性肺病)#   
#   - 暴露变量(X): PA_CAT       (身体活动量     0: 低           1: 正常， )
#   - 中介变量(M): CESD_CAT     (抑郁症状评分   0：未抑郁       1：抑郁) 
#   - 协变量(C): sex, age       (性别 0: 女性 1: 男性;  年龄(岁))
# ============================================================================
#1.安装包，加载包
#安装包
remotes::install_github("nhejazi/medoutcon")
remotes::install_github("tlverse/sl3")

#加载包
library(medoutcon)#用于中介分析和因果推断的包，基于目标最大似然估计（TMLE）方法
library(sl3)#现代机器学习框架
library(data.table)#数据操作和转换包
library(dplyr)#高性能数据处理包
library(speedglm)#快速广义线性模型
#2.数据准备
# 查看当前工作空间
getwd()

# 设置工作空间
setwd("C:/机器学习因果推断课程中介分析王老师")

#导入数据
data <- readr::read_csv("data.csv")

df <- data            #将原始数据框赋值给新的对象df
W <- c("sex", "age")  #定义了一组基准协变量
A <- "PA_CAT"       #定义处理/暴露变量
M <- "CESD_CAT"   #定义中介变量
Y <- "lungfun"    #定义结局变量


# 用纯GLM learner，也可以修改为其它的机器学习模型
glm_learner <- Lrnr_glm$new() #创建一个广义线性模型（GLM）的学习器


#计算自然直接效应
set.seed(123)    # 设置随机种子，确保结果可重复
nde_glm <- medoutcon(
  W = df[, W], A = df[[A]], Z = NULL, M = df[[M]], Y = df[[Y]],
  effect = "direct",      #估计自然直接效应
  # 以下四个参数都使用GLM进行估计
   g_learners = glm_learner,# 处理机制模型估计 P(A|W) - 给定协变量下处理的概率
  h_learners = glm_learner,  #中介机制模型估计 P(M|A,W) - 给定处理和协变量下中介变量的分布
  b_learners = glm_learner, # 结果机制模型估计 E[Y|A,M,W] - 给定处理、中介和协变量下的结果期望
  q_learners = glm_learner # 中介密度比模型，用于计算中介变量的密度比，是 TMLE 的关键部分
)

#计算自然间接效应
set.seed(123)
nie_glm <- medoutcon(
  W = df[, W], A = df[[A]], Z = NULL,#没有其他中介变量
  M = df[[M]], Y = df[[Y]],
  effect = "indirect", # 估计自然间接效应
  g_learners = glm_learner,
  h_learners = glm_learner,
  b_learners = glm_learner,
  q_learners = glm_learner
)
# 计算中介百分比
set.seed(123)

pm_glm <- medoutcon(
  W = df[, W], A = df[[A]], Z = NULL, M = df[[M]], Y = df[[Y]],
  effect = "pm",
  g_learners = glm_learner,
  h_learners = glm_learner,
  b_learners = glm_learner,
  q_learners = glm_learner
)

cat("GLM结果:\n")
cat("NDE:", nde_glm$theta, "\n")
cat("NIE:", nie_glm$theta, "\n")
cat("TE:", nde_glm$theta + nie_glm$theta, "\n")
cat("PM (手动):", nie_glm$theta / (nde_glm$theta + nie_glm$theta), "\n")
cat("PM (medoutcon):", pm_glm$theta, "\n")


results <- function(nde_glm,nie_glm){
  NDE <- nde_glm %>% summary %>% as.data.frame() %>% 
    dplyr::mutate(
      效应类型 = "自然直接效应 (NDE)",
      标准误 = nde_glm$var %>% sqrt, #对方差进行开方
      p值 = 2 * pnorm(-abs(nde_glm$theta / sqrt(nde_glm$var)))
    ) %>% 
    dplyr::select(效应类型,估计值 = param_est,标准误,CI下限 = lwr_ci,CI上限 = upr_ci,p值)
  NIE <- nie_glm %>% summary %>% as.data.frame() %>% 
    dplyr::mutate(
      效应类型 = "自然间接效应 (NIE)",
      标准误 = nie_glm$var %>% sqrt, #对方差进行开方
      p值 = 2 * pnorm(-abs(nie_glm$theta / sqrt(nie_glm$var)))
    ) %>% 
    dplyr::select(效应类型,估计值 = param_est,标准误,CI下限 = lwr_ci,CI上限 = upr_ci,p值)
  rbind(NDE,NIE) %>% as.data.frame()
}
results(nde_glm,nie_glm)


results <- function(nde_glm, nie_glm,pm_glm=NULL){
  # NDE结果
  NDE <- nde_glm %>% summary %>% as.data.frame() %>% 
    dplyr::mutate(
      效应类型 = "自然直接效应 (NDE)",
      标准误 = nde_glm$var %>% sqrt,
      p值 = 2 * pnorm(-abs(nde_glm$theta / sqrt(nde_glm$var)))
    ) %>% 
    dplyr::select(效应类型, 估计值 = param_est, 标准误, CI下限 = lwr_ci, CI上限 = upr_ci, p值)
  # NIE结果
  NIE <- nie_glm %>% summary %>% as.data.frame() %>% 
    dplyr::mutate(
      效应类型 = "自然间接效应 (NIE)",
      标准误 = nie_glm$var %>% sqrt,
      p值 = 2 * pnorm(-abs(nie_glm$theta / sqrt(nie_glm$var)))
    ) %>% 
    dplyr::select(效应类型, 估计值 = param_est, 标准误, CI下限 = lwr_ci, CI上限 = upr_ci, p值)
  # 总效应 (TE = NDE + NIE)
  te_est <- nde_glm$theta + nie_glm$theta
  # 方差：假设NDE和NIE独立，var(TE) = var(NDE) + var(NIE)
  # 注意：这是简化假设，实际上它们可能相关
  te_var <- nde_glm$var + nie_glm$var
  te_se <- sqrt(te_var)
  te_ci_lower <- te_est - 1.96 * te_se
  te_ci_upper <- te_est + 1.96 * te_se
  te_pval <- 2 * pnorm(-abs(te_est / te_se))
  TE <- data.frame(
    效应类型 = "总效应 (TE)",
    估计值 = te_est,
    标准误 = te_se,
    CI下限 = te_ci_lower,
    CI上限 = te_ci_upper,
    p值 = te_pval
  )
  # 合并结果
  result_table <- rbind(NDE, NIE, TE)
  
  if (!is.null(pm_glm)) {
    cat("\n============== 中介比例 ==============\n")
    cat("方法1 (手算 NIE/TE):", round((nie_glm$theta / te_est) * 100, 2), "%\n")
    cat("方法2 (medoutcon pm):", round(pm_glm$theta * 100, 2), "%\n")
    cat("  标准误:", round(sqrt(pm_glm$var), 3), "\n")
    cat("  95% CI: [", round(pm_glm$theta - 1.96*sqrt(pm_glm$var), 3), 
        ",", round(pm_glm$theta + 1.96*sqrt(pm_glm$var), 3), "]\n")
    cat("=====================================\n\n")
  } else {
    # 计算中介比例
    mediation_prop <- (nie_glm$theta / te_est) * 100
    cat("\n中介比例 (NIE/TE):", round(mediation_prop, 2), "%\n\n")
  }
  return(result_table)
}

# 调用函数--手搓中介百分比
final_results <- results(nde_glm, nie_glm)
print(final_results)



# 调用函数--采用程序计算的中介百分比
final_results <- results(nde_glm, nie_glm,pm_glm)
print(final_results)