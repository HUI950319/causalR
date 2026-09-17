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

# 安装r包
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("devtools")) install.packages("devtools")
#安装CMAverse包
devtools::install_local("C:/机器学习因果推断课程中介分析王老师/CMAverse") #本地安装
devtools::install_github("BS1125/CMAverse")     # 2.在线安装（上面安装不成功，备选）

# 加载r包
library(dplyr)
library(ggplot2)
library(devtools)
library(CMAverse)

#2.数据准备
# 查看当前工作空间
getwd()

# 设置工作空间
setwd("C:/机器学习因果推断课程中介分析王老师")

#导入数据
data<- read.csv("data.csv")
# ----------------------------

#########二、CMAverse 包的因果中介分析########

## 自然效应模型 ne（Natural Effect Model）
set.seed(2024)
cm_ne <- cmest(
  data = data,
  model = "ne",             # 指定模型为自然效应模型
  outcome = "lungfun",      # Y：结局变量（二分类）
  exposure = "PA_CAT",      # X：暴露变量（二分类）
  mediator = "CESD_CAT",    # M：中介变量（二分类）
  basec = c("sex", "age"),  # C：协变量（混杂因素）
  yreg = "logistic",        # Y 的回归模型（logistic 回归）
  mreg = list("logistic"),  # M 的回归模型（logistic 回归）
  EMint = FALSE,            # 是否指定暴露–中介交互（此处无）
  estimation = "imputation",# 使用反事实插补进行效应估计
  inference  = "bootstrap", # 使用 bootstrap 进行推断
  nboot      = 100,         # bootstrap 次数（建议 ≥1000）
  yval = 1,                 # Y=1 的效应（肺功能异常）
  mval = list(1)            # M=1 的中介类别（抑郁）
)
summary(cm_ne)  


# 从 cm_ne 中提取关键效应指标 
df_ne <- data.frame(
  effect = c("NDE", "NIE", "TE"), # 三种效应类型NDE/NIE/TE
  OR = c(                        
    cm_ne$effect.pe["Rpnde"],     # 点估计值（OR值）
    cm_ne$effect.pe["Rpnie"],
    cm_ne$effect.pe["Rte"]
  ),
  lower = c(
    cm_ne$effect.ci.low["Rpnde"],  # CI下限
    cm_ne$effect.ci.low["Rpnie"],
    cm_ne$effect.ci.low["Rte"]
  ),
  upper = c(
    cm_ne$effect.ci.high["Rpnde"], # CI上限
    cm_ne$effect.ci.high["Rpnie"],
    cm_ne$effect.ci.high["Rte"]
  )
)
df_ne

#绘制带误差线的森林图
ggplot(df_ne, aes(x = effect, y = OR, color = effect)) +
  geom_point(size = 4) +    # 点估计值
  geom_errorbar(aes(ymin = lower, ymax = upper),  # 置信区间
                width = 0.15, size = 1) +        
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") + # 参考线(OR=1)
  coord_flip() +   # 翻转坐标轴
  theme_bw() +     # 黑白主题
  labs(
    title = "Natural Effect Model",
    x = "",
    y = "Odds Ratio (OR)"
  )

# 敏感性分析
sens_ne <- cmsens(cm_ne, sens="uc") #指定敏感性分析的类型为uc（未测量混杂）
print(sens_ne$evalues)              # 输出结果    