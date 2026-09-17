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
# ----------------------------
# 1.安装和加载必要的R包
# ============================================================================
# 安装包（如果尚未安装）
install.packages("mediation")  # 核心中介分析包
install.packages("tidyverse")  # 数据整理和可视化工具集
install.packages("broom")      # 整理模型结果
install.packages("EValue")     # E-value敏感性分析

# 加载所有需要的包
library(mediation)     # 用于中介效应分析
library(tidyverse)     # 包含dplyr(数据处理)、ggplot2(绘图)等
library(broom)         # 将模型结果转为整洁的数据框
library(EValue)        # 进行敏感性分析的E-value计算（备选方法）


#2.数据准备
# 查看当前工作空间
getwd()

# 设置工作空间
setwd("C:/机器学习因果推断课程中介分析王老师")

#3.导入数据
data<- read.csv("data.csv")
# 查看数据结构
str(data)
head(data)
# 3. 建立中介模型
# ----------------------------
# 注意：mediation包需要两个回归模型

# 模型1: 中介变量模型 (M ~ X + covariates)
# 使用逻辑回归，因为M是二分类变量
med_model <- glm(
  formula = CESD_CAT ~ PA_CAT + sex + age,          # 中介变量受自变量和协变量影响
  family = binomial(link = "logit"),  # 二分类变量使用logit连接
  data = data
)

# 查看中介模型的系数
summary(med_model)
tidy(med_model)  # 使用broom包整理结果

# 模型2: 因变量模型 (Y ~ X + M + covariates)
# 使用逻辑回归，因为Y是二分类变量
out_model <- glm(
  formula = lungfun ~ PA_CAT+ CESD_CAT + sex + age,      # 因变量受自变量、中介变量和协变量影响
  family = binomial(link = "logit"),  # 二分类变量使用logit连接
  data = data
)

# 查看因变量模型的系数
summary(out_model)
tidy(out_model)  # 使用broom包整理结果

# ----------------------------
# 4. 进行中介效应分析
# ----------------------------
# 使用mediate()函数计算中介效应

med_result <- mediate(
  model.m = med_model,    # 中介变量模型
  model.y = out_model,    # 结果变量模型
  treat = "PA_CAT",            # 自变量的名称
  mediator = "CESD_CAT",         # 中介变量的名称
  sims = 100,            # Bootstrap重复次数，建议1000次以上
  boot = TRUE,            # 使用非参数bootstrap
  robustSE = TRUE,        # 使用稳健标准误（对二分类变量推荐）
  treat.value = 1,        # 处理组值（用于计算效应）
  control.value = 0,      # 对照组值（用于计算效应）
  boot.ci.type = "perc"   # Bootstrap置信区间类型："perc"百分位法
)

# 查看详细的中介分析结果
summary(med_result)


# 解释输出：
# ACME (Average Causal Mediation Effect): 平均中介效应（间接效应）
# ADE (Average Direct Effect): 平均直接效应
# Total Effect: 总效应 = 间接效应 + 直接效应
# Prop. Mediated: 中介效应占总效应的比例

# ----------------------------
# 5. 结果可视化
# ----------------------------
# 绘制中介效应图

#自定义绘图（使用ggplot2）
# 提取Bootstrap分布结果
if(!is.null(med_result$d0.sims)) {
  effect_data <- data.frame(
    ACME = med_result$d0.sims,  # 中介效应分布
    ADE = med_result$z0.sims    # 直接效应分布
  )
  
  # 绘制中介效应和直接效应的分布
  p <- effect_data %>%
    pivot_longer(cols = everything(), names_to = "Effect", values_to = "Value") %>%
    ggplot(aes(x = Value, fill = Effect)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ Effect, scales = "free") +
    labs(title = "中介效应和直接效应的Bootstrap分布",
         x = "效应大小", y = "密度") +
    theme_minimal()
  print(p)
}

# ----------------------------
# 6. 敏感性分析
# ----------------------------
# 对于二分类中介变量，mediation包的medsens()函数不适用

# 使用E-value进行敏感性分析
# E-value表示需要多强的未观测混淆才能解释观察到的效应
# 注意：这需要将效应大小转换为风险比(risk ratio)

# 首先，我们需要获得效应估计

# 计算总效应、间接效应、直接效应E-value
total_effect <- med_result$tau.coef  # 总效应估计值
indirect_effect<- med_result$d0#间接效应估计值
direct_effect<- med_result$z0#直接效应估计值
# 将log-odds转换为近似的风险比（近似公式）
# 注意：这是近似计算，严格来说需要得到风险比
# 计算总效应的E-value
evalue_result_total <- EValue::evalues.RR(
  est = exp(total_effect),  # 估计的风险比
  lo = exp(med_result$tau.ci[1]),  # 置信区间下限
  hi = exp(med_result$tau.ci[2])   # 置信区间上限
)
# 计算间接效应的E-value
evalue_result_indirect <- EValue::evalues.RR(
  est = exp(indirect_effect),  # 估计的风险比
  lo = exp(med_result$d0.ci[1]),  # 置信区间下限
  hi = exp(med_result$d0.ci[2])   # 置信区间上限
)
# 计算直接效应的E-value
evalue_result_direct <- EValue::evalues.RR(
  est = exp(direct_effect),  # 估计的风险比
  lo = exp(med_result$z0.ci[1]),  # 置信区间下限
  hi = exp(med_result$z0.ci[2])   # 置信区间上限
)

# E-value敏感性分析结果:
evalue_result_total
evalue_result_indirect
evalue_result_direct

