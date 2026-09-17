
# ============================================================================
# 课程演示：G-computation 与 倾向性评分(IPTW) 的实现
# —— 在中间量估计中引入 Super Learner
#
# 数据：rhc.csv（RHC 研究数据）
#
# 研究问题（示例）：
#   实施 RHC (A=1) 相比不实施 RHC (A=0)，是否影响死亡风险 (Y)？
#
# 变量说明：
#   - 处理/干预 A：RHC
#       1 = 行 RHC
#       0 = 未行 RHC
#
#   - 结局 Y：death
#       Yes / No （分析中转换为 1 / 0）
#
#   - 混杂变量 W（基线协变量）：
#       Cardiovascular, Congestive.HF, Age, Edu, DASIndex,
#       APACHE.score, blood.pressure, WBC, Heart.rate, DNR.status
#       （这些变量同时影响 RHC 的实施与死亡风险）
#
# 本脚本实现两条经典因果估计路径：
#
#   (1) G-computation（G 计算）：
#       - 目标：估计 E[Y(1)] 与 E[Y(0)]
#       - 做法：使用 Super Learner 拟合结果模型 E(Y | A, W)，
#               并基于该模型模拟反事实结局
#
#   (2) IPTW（逆概率加权）：
#       - 目标：构造一个“处理近似随机分配”的加权总体
#       - 做法：使用 Super Learner 估计倾向性评分
#               e(W) = P(A = 1 | W)，并据此生成权重
#
# 重要说明：
#   - G-computation 与 IPTW 是因果推断的估计框架
#   - Super Learner 仅作为灵活的机器学习工具，
#     用于估计上述方法所需的中间量
#   - 因果解释依赖于无未测混杂、重叠性和一致性等假设
# ============================================================================


# ===============================
# 0. 安装并加载需要的包
# ===============================

# install.packages() 只需要第一次安装时运行一次
install.packages(c(
  "readr",        # 读 csv
  "dplyr",        # 数据整理
  "SuperLearner", # Super Learner 框架
  "cobalt",       # 平衡性检查 (love plot / SMD)
  "WeightIt",     # 倾向性评分加权工具（可调用 super learner）
  "survey",       # 加权后推断（标准误、置信区间）
  "boot"          # bootstrap（用于 G-computation 的 CI）
))

# 每次打开 R 都需要 library() 加载包 
library(readr)
library(dplyr)
library(SuperLearner)
library(cobalt)
library(WeightIt)
library(survey)
library(boot)

# ===============================
#### 1. 读取数据 ####
# ===============================
data <- read_csv("C:/Users/zhangyan/Desktop/ML CR/rhc.csv")

# glimpse() 可以快速显示变量类型、缺失情况、前几行
glimpse(data)

# ===============================
#### 2. 定义变量+数据预处理 ####
# ===============================

## 2.1 定义变量

# 干预变量名: RHC (RHC干预:1=行RHC，0=未行RHC)
A <- "RHC"

# 结局变量名: death (是否死亡:No，Yes)
Y <- "death"

# 基线协变量 W（混杂变量）：影响 A 和 Y
W <- c("Cardiovascular","Congestive.HF","Age","Edu","DASIndex", 
       "APACHE.score","blood.pressure","WBC","Heart.rate","DNR.status")

# 只保留分析需要的列（避免其他列干扰建模/出错）
data2 <- data %>% select(all_of(c(A, Y, W)))


## 2.2 检查 A 和 Y 的原始类型 
# 很多数据集里A/Y是factor或字符型,需转换成0/1数值才能用于binomial模型

str(data2[[A]])
str(data2[[Y]])

# 将 Y 转为 0/1 (# death是 "Yes"/"No",把 Yes->1, No->0)
data2[[Y]] <- ifelse(data2[[Y]] == "Yes", 1, 0)

# 查看 A 和 Y 的分布
table(data2[[A]])
table(data2[[Y]])

# ===============================
#### 3. G-computation + Super Learner #### 
#    目标：估计 E[Y(1)]与E[Y(0)], ATE = E[Y(1)]-E[Y(0)] (风险差)
# ===============================

## 3.1 用 Super Learner 拟合结局模型

Y_vec <- data2[[Y]]  # 结局向量
X_outcome <- data2 %>% 
          select(all_of(c(A, W)))  # 自变量 = 干预(A) + 混杂(w)

# 使用 Super Learner结局模型
set.seed(2026)   # 结果可复现
fit_outcome <- SuperLearner(
  Y = Y_vec,   
  X = X_outcome,    # A + w
  family = binomial(),  # 结局是0/1->二项分布(logit 链接)
  SL.library = c(
    "SL.mean",   # 基线：只用均值
    "SL.glm",    # Logistic 回归
    "SL.glmnet"  # LASSO
    )
  )    # SL.library 中可添加："SL.xgboost", "SL.randomforest"

# 输出 Super Learner 结果(Risk:预测误差, Coef:权重)
fit_outcome

# 查看各学习器在集成中的权重（系数）
fit_outcome$coef


## 3.2 构造反事实数据：把所有人的 A 人为设为 1 / 0

# 反事实世界1：所有人都做 RHC (A=1)
# 反事实世界0：所有人都不做 RHC (A=0)

X1 <- X_outcome; X1[[A]] <- 1   # 所有人都 RHC
X0 <- X_outcome; X0[[A]] <- 0   # 所有人都不 RHC

# 用拟合好的Q 模型预测每个人在两种世界下的死亡概率（反事实风险）
p1 <- predict(fit_outcome, newdata = X1)$pred
p0 <- predict(fit_outcome, newdata = X0)$pred


## 3.3 计算边际风险与 ATE（风险差 RD）
# E[Y(1)]：把所有人的反事实概率取平均
# E[Y(0)]：同理

risk1_g <- mean(p1)         # 如果所有人都做RHC的死亡风险
risk0_g <- mean(p0)         # 如果所有人都不做RHC的死亡风险
RD_g <- risk1_g - risk0_g   # 风险差 ATE

res_g <- c(risk_if_A1 = risk1_g, risk_if_A0 = risk0_g, RD = RD_g)
print(res_g)


## 3.4 使用bootstrap 计算G-computation 的 percentile CI 
# bootstrap：重复“抽样 -> 拟合结局回归函数 -> 做反事实预测 -> 求 ATE”
# 得到 ATE 的经验分布，从而构造置信区间

gcomp_ate <- function(data, indices) {
  # data：原始数据；indices：bootstrap 重抽样行号
  d <- data[indices, ]
  
  # (1) 拟合结局模型 Q(A,W)
  sl_fit_boot <- SuperLearner(
    Y = d[[Y]],
    X = d[, c(A, W)],
    SL.library = c("SL.mean","SL.glm","SL.glmnet"),
    family = binomial()
  )
  
  # (2) 构造反事实数据（只需要 A+W 列即可）
  d_A1 <- d_A0 <- d
  d_A1[[A]] <- 1
  d_A0[[A]] <- 0
  
  # (3) 预测反事实风险（onlySL=TRUE 表示使用集成器的预测）
  pred_A1 <- predict(sl_fit_boot, newdata = d_A1[, c(A, W)], onlySL = TRUE)$pred
  pred_A0 <- predict(sl_fit_boot, newdata = d_A0[, c(A, W)], onlySL = TRUE)$pred
  
  # (4) 返回该次 bootstrap 的 ATE（风险差）
  mean(pred_A1) - mean(pred_A0)
}

# 执行bootstrap（考虑到计算时间，这里用5次，实际应用中建议1000+次）
set.seed(2026)
boot_results <- boot(data2, gcomp_ate, R = 5)

# 计算 percentile CI
CI_g <- boot.ci(boot_results, type = "perc")

# CI_g$percent 是一个向量：其中第4和第5个元素分别是 2.5% 和 97.5% 分位点
CI_low_g  <- CI_g $ percent[4]
CI_high_g <- CI_g $ percent[5]

# G 计算的结果:点估计（用原样本 RD_g）+ 置信区间（bootstrap percentile）
c(ATE = RD_g, CI_low = CI_low_g, CI_high = CI_high_g)


# ===============================
#### 4. IPTW（倾向性评分加权） + Super Learner ####
#   目标：用 e(W)=P(A=1|W) 构造权重，使加权后A与W近似独立（像随机试验）
# ===============================

## 4.1 用 WeightIt(method="super") 估计 PS 并生成权重

set.seed(2026)
w.out <- weightit(
  formula   = as.formula(paste(A, "~", paste(W, collapse = " + "))),
              # 等价于：A ~ W1 + W2 + ... + Wp
  data      = data2,
  estimand  = "ATE",      # 目标是平均处理效应ATE
  method    = "super",    # Super Learner 来拟合 e(W)
  SL.library = c(
      "SL.mean",      # 用总体均值预测 A（基线模型）
      "SL.glm",       # Logistic回归
      "SL.glmnet"     # LASSO
    )
  )

# 查看权重概况
summary(w.out)

# 查看权重分布
w <- w.out$weights 
summary(w)


## 4.2 平衡性检查:加权后两组是否更像随机？
# 看加权前后标准化均值差（SMD）
# 常用经验标准：|SMD| < 0.1 认为平衡不错

bal.tab(
  w.out,                 # weightit() 生成的倾向性评分权重对象
  un = TRUE,             # un = TRUE：显示“未加权”和“加权后”的结果
  thresholds = c(m = 0.1) # 设置平衡性阈值：|SMD| < 0.1 
)

# love.plot：可视化加权前后 SMD
love.plot(
  w.out,                 # 同样使用 weightit() 输出对象
  binary = "std",        # 对二分类变量使用“标准化均值差(SMD)”作为平衡性度量
  thresholds = c(m = 0.1),# 在图中画出|SMD| = 0.1 的阈值线
  abs = TRUE,            # 使用“绝对值 SMD”
  var.order = "unadjusted", # 按“未加权时的不平衡程度”排序变量
  line = TRUE            # 用线连接未加权和加权后
)

## 4.3 用 survey 在加权伪总体中估计风险差 + 95%CI
# survey 包把权重当作“抽样权重”，从而可计算标准误和置信区间

# 建立加权设计对象
design.w <- svydesign(ids = ~1,   # 每一行数据是一个独立个体
                      weights = ~w,  # 抽样权重
                      data = data2)

# 用线性概率模型(gaussian)来直接估计风险差（RD）与 CI
fit_rd <- svyglm(
  as.formula(paste(Y, "~", A)),  # 模型公式：Y ~ A
  design = design.w,             # 指定 survey 设计对象
  family = gaussian()
  )
# 提取 RD（RHC 的回归系数）以及置信区间
ATE_w <- coef(fit_rd)[A]    # ATE（风险差 RD）
CI_w  <- confint(fit_rd)[A, ]   # 95% CI
# CI_w 是长度为2的向量：下限、上限
CI_low_w <- CI_w [1]
CI_high_w <- CI_w [2]

# 输出倾向性评分加权的结果
c(ATE = ATE_w, CI_low = CI_low_w, CI_high = CI_high_w)




