###安装包####
###一次性安装下列包######
install.packages("AIPW")
install.packages("SuperLearner")
install.packages("tableone")
install.packages("grf")
install.packages("DoubleML")
install.packages("mlr3")
install.packages("mlr3learners")
install.packages("xgboost")

install.packages("randomForest")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("hstats")
install.packages("tidyr")
install.packages("mlr3misc")
install.packages("boot")
install.packages("mlr3tuning")

install.packages("patchwork")
install.packages("pdp")
install.packages("tmle")
install.packages("EpiForsk")
install.packages("hstats")
install.packages("DeepLearningCausal")

# pkgs_needed <- c(
# 	"AIPW",  
# 	"SuperLearner",   # AIPW依赖此包进行机器学习建模
# 	"tableone",  # 基线表 / SMD
# 	"SuperLearner",      # 用于S/T/X-learner的基础学习器
# 	"grf",               # 用于因果森林
# 	"DoubleML",          # 用于DML
# 	"mlr3",              # 机器学习框架
# 	"mlr3learners",      # 机器学习学习器
# 	"xgboost",           # XGBoost
# 	"randomForest",      # 随机森林
# 	"ggplot2",           # 绘图
# 	"dplyr",             # 数据处理
# 	"tidyr",             # 数据整理
# 	"hstats",            # 因果森林
# 	"boot",
# 	"mlr3misc",
# 	"mlr3",
# 	"mlr3learners",
# 	"mlr3tuning",
# 	"DoubleML",
# 	"patchwork",
# 	"pdp",               # 用于PDP图
# 	"tmle",              # tmle
# 	"EpiForsk"
# )

# # 当DoubleML安装不顺利的替代方法
# # 关闭RStudio，重新打开（以管理员身份）
# # 然后运行：
# remove.packages(c("mlr3misc", "mlr3", "mlr3learners", "DoubleML", "digest", "cli", "checkmate"))
# 
# # 安装remotes（如果没有）
# if (!require(remotes)) install.packages("remotes")
# 
# # 安装mlr3生态系统（按顺序）
# install.packages("mlr3misc")
# install.packages("mlr3")
# install.packages("mlr3learners")
# install.packages("mlr3tuning")  # DoubleML需要这个
# # 使用官方文档中的命令
# remotes::install_github("DoubleML/doubleml-for-r")
# # 检查版本
# packageVersion("DoubleML")

# pkgs_to_install <- pkgs_needed[!vapply(pkgs_needed, requireNamespace, logical(1), quietly = TRUE)]
# if (length(pkgs_to_install) > 0) {
# 	install.packages(pkgs_to_install, dependencies = TRUE)
# }


###1.数据导入######
# 设置当前工作目录（按自己的电脑路径调整）
setwd("D:/因果训练营/异质性分析课程")


# 读取数据（data.csv 需放在上述工作目录下）
data <- read.csv("data1.csv")

# 协变量
names(data)
covariates <- c("sex","hypertension","diabetes","heart_disease","residence","marital",
								"smoking","age","education_level","cardio_metabolic",
								"alcohol","bmi")

# 需要按分类变量处理的协变量（转为 factor，便于后续自动生成哑变量/分层比较）
categorical_vars <- c("sex","hypertension","diabetes","heart_disease","residence","marital",
											"smoking","education_level","cardio_metabolic",
											"alcohol")

# 统一把分类协变量转为 factor
data[categorical_vars] <- lapply(data[categorical_vars], factor)

# 暴露变量检查（应为 0/1 或二分类因子）
table(data$PA_CAT)
str(data$PA_CAT)
####均衡性检验####
## 目的：展示 PA_CAT 两组在协变量上的初始不平衡（可看 SMD/检验）
library(tableone)

# 根据实际数据调整变量列表
# 定义想要比较的变量
vars <- c("age", "bmi", "sex", "hypertension", "diabetes", 
					"heart_disease", "residence", "marital", "smoking", 
					"education_level", "cardio_metabolic", "alcohol")

# 分类变量列表
categorical_vars <- c("sex", "hypertension", "diabetes", "heart_disease", 
											"residence", "marital", "smoking", "education_level", 
											"cardio_metabolic", "alcohol")

# 创建表1
table_one <- CreateTableOne(vars = vars, 
														strata = "PA_CAT",  # 暴露分组变量
														data = data, 
														factorVars = categorical_vars,  # 指定分类变量
														test = TRUE)  # 进行统计检验

# 打印结果，显示SMD（标准化均数差）
print(table_one, smd = TRUE)

###2.疗效评估（计算ATE）#####
#####2.1 AIPW方法######
library(AIPW)
library(SuperLearner)
# ①准备数据
# 首先，需要确定结局变量（Y）
# - lungfun: 肺功能（二分类）
Y <- data$lungfun  

# ②治疗变量（暴露变量）
# PA_CAT: 0=不活动, 1=活动
A <- data$PA_CAT

# ③协变量矩阵（混杂因素）
# 注意：SuperLearner要求所有协变量为数值型
# 方法1：使用您已定义的协变量列表
covariates <- c("sex", "hypertension", "diabetes", "heart_disease", 
								"residence", "marital", "smoking", "age", 
								"education_level", "cardio_metabolic", "alcohol", "bmi")

# 创建协变量矩阵，将因子变量转换为数值型
W <- data[, covariates]

# 将因子变量转换为数值型（SuperLearner要求）
# 对于二分类因子，转换为0/1
# 对于多分类因子，创建哑变量
for(var in names(W)) {
	if(is.factor(W[[var]])) {
		if(nlevels(W[[var]]) == 2) {
			# 二分类变量：转换为0/1
			W[[var]] <- as.numeric(W[[var]]) - 1
		} else {
			# 多分类变量：创建哑变量
			dummies <- model.matrix(~ W[[var]] - 1)
			colnames(dummies) <- paste0(var, "_", levels(W[[var]]))
			# 移除原始变量，添加哑变量
			W <- cbind(W, dummies)
			W[[var]] <- NULL
		}
	}
}

# 确保所有变量都是数值型
W <- as.data.frame(lapply(W, as.numeric))

# ④检查数据完整性，当然也可以采用缺失数据填补的方法，简单一些的就是随机森林的方法
# 移除缺失值（如果有）
complete_cases <- complete.cases(Y, A, W)
if(sum(!complete_cases) > 0) {
	cat("移除", sum(!complete_cases), "个缺失值样本\n")
	Y <- Y[complete_cases]
	A <- A[complete_cases]
	W <- W[complete_cases, ]
}

# ⑤设定机器学习算法库
# 根据您的数据特点选择合适的算法(为快速运行，这里仅选用其中一种运行)
sl_lib <- c("SL.glm"           # 线性/逻辑回归（基础模型）
						#"SL.gam",           # 广义可加模型（捕捉非线性）
						#"SL.randomForest",  # 随机森林（捕捉复杂交互）
						#"SL.glmnet",      # 带正则化的回归（防止过拟合）
)           #"SL.xgboost"       


# ⑥运行AIPW算法
set.seed(456)  # 设置随机种子，确保结果可重复

aipw_result <- AIPW$new(
	Y = Y, 
	A = A, 
	W = W,
	Q.SL.library = sl_lib,      # 结果模型的算法库
	g.SL.library = sl_lib,      # 倾向得分模型的算法库
	k_split = 5,                 # 5折交叉验证
	verbose = TRUE               # 显示进度信息
)$fit()

# ⑦ 查看结果
cat("\n========== AIPW分析结果 ==========\n")
aipw_result$summary()
print(aipw_result$result)

# ⑧提取关键结果
# 风险差（Risk Difference）就是平均处理效应（ATE）
ATE_est <- aipw_result$result[3, "Estimate"]  # 第3行是Risk Difference
ATE_se <- aipw_result$result[3, "SE"]
ATE_ci_lower <- aipw_result$result[3, "95% LCL"]
ATE_ci_upper <- aipw_result$result[3, "95% UCL"]
# 计算z统计量和p值
z_score <- ATE_est / ATE_se
ATE_p <- 2 * (1 - pnorm(abs(z_score)))###因AIPW本身不能用于计算P值，需通过95%Cl进行计算

cat("\n========== 平均处理效应 (ATE) ==========\n")
cat(sprintf("估计值: %.3f\n", ATE_est))
cat(sprintf("标准误: %.3f\n", ATE_se))
cat(sprintf("95%% CI: [%.3f, %.3f]\n", ATE_ci_lower, ATE_ci_upper))
cat(sprintf("P值: %.4f\n", ATE_p))

####2.3.异质性分析方法--疾病风险评分法（Disease Risk Score, DRS）#####
### 3.1 数据准备 
library(dplyr)

setwd("D:/因果训练营/异质性分析课程")
data <- read.csv("data1.csv")

covariates <- c("sex", "hypertension", "diabetes", "heart_disease", 
								"residence", "marital", "smoking", "age", 
								"education_level", "cardio_metabolic", "alcohol", "bmi")

categorical_vars <- c("sex", "hypertension", "diabetes", "heart_disease", 
											"residence", "marital", "smoking", "education_level", 
											"cardio_metabolic", "alcohol")
data[categorical_vars] <- lapply(data[categorical_vars], factor)

# 处理缺失值
complete_idx <- complete.cases(data[, covariates])
data <- data[complete_idx, ]

### 3.2 拟合风险模型 

risk_model <- glm(lungfun ~ sex + hypertension + diabetes + heart_disease + 
										residence + marital + smoking + age + education_level + 
										cardio_metabolic + alcohol + bmi,
									data = data,
									family = binomial())

### 再对照组中拟合模型
# control_data <- data[data$PA_CAT == 0, ]
# risk_model <- glm(lungfun ~ sex + hypertension + diabetes + heart_disease + 
#                     residence + marital + smoking + age + education_level + 
#                     cardio_metabolic + alcohol + bmi,
#                   data = control_data,
#                   family = binomial())

### 3.3 计算风险评分 
data$risk_score <- predict(risk_model, newdata = data, type = "response")


### 3.4 对风险评分划分4组
data$risk_quartile <- cut(data$risk_score,
													breaks = quantile(data$risk_score, probs = seq(0, 1, 0.25)),
													include.lowest = TRUE,
													labels = c("Q1", "Q2", "Q3", "Q4"))

####3.5 在各亚组计算ATE，传统回归方法
# 交互效应模型
model_interaction <- glm(lungfun ~ PA_CAT * risk_score,
												 data = data,
												 family = binomial())
summary(model_interaction)

# 在各组计算ATE，传统回归方法，以Q1组为例，此时没有混杂因素，单因素分析即可。
data11<-data[data$risk_quartile=="Q1",]
table(data11$risk_quartile)
model_main <- glm(lungfun ~ PA_CAT,
									data = data11,
									family = binomial())
summary(model_main)


# ##### 4. 效应评分法亚组分析 - 各层ATE ####

library(DeepLearningCausal)
library(SuperLearner)
# 设置随机种子
set.seed(456)

# 2. 准备数据
# 设置当前工作目录（按自己的电脑路径调整）
setwd("D:/因果训练营/异质性分析课程")

# 读取数据（data.csv 需放在上述工作目录下）
data <- read.csv("data1.csv")

# 协变量（潜在混杂因素）
names(data)
# 定义协变量
covariates <- c("sex","hypertension","diabetes","heart_disease","residence","marital",
								"smoking","age","education_level","cardio_metabolic",
								"alcohol","bmi")

# 需要转为 factor 的分类变量
categorical_vars <- c("sex","hypertension","diabetes","heart_disease","residence","marital",
											"smoking","education_level","cardio_metabolic","alcohol")

# 转换分类变量为 factor
data[categorical_vars] <- lapply(data[categorical_vars], factor)

# 删除缺失值（只针对协变量）
complete_idx <- complete.cases(data[, covariates])
data_clean <- data[complete_idx, ]

# 提取变量（虽然 metalearner_ensemble 直接使用 data_clean 即可，但这里保留以备后用）
Y <- data_clean$lungfun
A <- data_clean$PA_CAT

# 定义 Super Learner 基础学习器库 #######
# 可根据需要补充，这里仅使用其中一种用于演示
SL_library <- c(
	# "SL.glmnet",        # 弹性网络
	# "SL.xgboost",       # XGBoost
	# "SL.randomForest",  # 随机森林
	# "SL.nnet",          # 神经网络,一般不用
	"SL.glm"            # 广义线性模型（基线）
)

# 定义公式（可根据需要修改，纳入所需要的协变量）
formula_str <- "lungfun ~ age + bmi + smoking + education_level + sex"

#### 3.1 S-Learner#####################
# S-Learner with Super Learner ensemble
result_S <- metalearner_ensemble(
	cov.formula = as.formula(formula_str),
	data = data_clean,
	treat.var = "PA_CAT",
	meta.learner.type = "S.Learner",
	SL.learners = SL_library,
	nfolds = 5,
	family = binomial(),
	binary.preds = FALSE
)

cate_S <- result_S$CATEs
ate_S <- mean(cate_S)
cat(sprintf("S-Learner ATE = %.6f\n", ate_S))
cat("CATE 分位数 (10%, 25%, 50%, 75%, 90%):\n")
print(quantile(cate_S, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)))

#### 3.2 T-Learner###########
# T-Learner with Super Learner ensemble
result_T <- metalearner_ensemble(
	cov.formula = as.formula(formula_str),
	data = data_clean,
	treat.var = "PA_CAT",
	meta.learner.type = "T.Learner",
	SL.learners = SL_library,
	nfolds = 5,
	family = binomial(),
	binary.preds = FALSE
)

cate_T <- result_T$CATEs
ate_T <- mean(cate_T)
cat(sprintf("T-Learner ATE = %.6f\n", ate_T))
cat("CATE 分位数 (10%, 25%, 50%, 75%, 90%):\n")
print(quantile(cate_T, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)))

#### 3.3 X-Learner#################
# X-Learner with Super Learner ensemble
result_X <- metalearner_ensemble(
	cov.formula = as.formula(formula_str),
	data = data_clean,
	treat.var = "PA_CAT",
	meta.learner.type = "X.Learner",
	SL.learners = SL_library,
	nfolds = 5,
	family = binomial(),
	binary.preds = FALSE
)

cate_X <- result_X$CATEs
ate_X <- mean(cate_X)
cat(sprintf("X-Learner ATE = %.6f\n", ate_X))
cat("CATE 分位数 (10%, 25%, 50%, 75%, 90%):\n")
print(quantile(cate_X, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)))

#### 3.4 DML (Double Machine Learning)--R-learn#################
# R-Learner with Super Learner ensemble
result_R <- metalearner_ensemble(
	cov.formula = as.formula(formula_str),
	data = data_clean,
	treat.var = "PA_CAT",
	meta.learner.type = "R.Learner",
	SL.learners = SL_library,
	nfolds = 5,
	family = binomial(),
	binary.preds = FALSE
)

cate_R <- result_R$CATEs
ate_R <- mean(cate_R)
cat(sprintf("R-Learner ATE = %.6f\n", ate_R))
cat("CATE 分位数 (10%, 25%, 50%, 75%, 90%):\n")
print(quantile(cate_R, probs = c(0.1, 0.25, 0.5, 0.75, 0.9)))
##### 结果汇总#####
# 构建 ATE 比较表格
results_ATE <- data.frame(
	Method = c("S-Learner", "T-Learner", "X-Learner", "R-Learner"),
	ATE = c(ate_S, ate_T, ate_X, ate_R)
)

print(results_ATE)

###4.因果森林########
###### 4.1 构建因果森林######
# 1. 设置工作目录和加载包
setwd("D:/因果训练营/异质性分析课程")

library(grf)
library(ggplot2)
library(dplyr)
library(patchwork)
library(pdp)  # 用于PDP图

# 2. 读取数据
data <- read.csv("data1.csv")

# 3. 数据预处理
covariates <- c("sex", "hypertension", "diabetes", "heart_disease", 
								"residence", "marital", "smoking", "age", 
								"education_level", "cardio_metabolic", "alcohol", "bmi")

categorical_vars <- c("sex", "hypertension", "diabetes", "heart_disease", 
											"residence", "marital", "smoking", "education_level", 
											"cardio_metabolic", "alcohol")

# 转换分类变量
for(var in categorical_vars) {
	data[[var]] <- as.numeric(as.character(data[[var]]))
}

# 定义处理变量和结局变量
W <- data$PA_CAT
Y <- data$lungfun

# 准备协变量矩阵
X <- as.matrix(data[, covariates])

# 删除缺失值
complete_idx <- complete.cases(X, W, Y)
X <- X[complete_idx, ]
W <- W[complete_idx]
Y <- Y[complete_idx]
data_clean <- data[complete_idx, ]

cat("样本量:", nrow(data_clean), "\n")
cat("处理组:", sum(W==1), "\n")
cat("对照组:", sum(W==0), "\n\n")

###### 4.1.1 拟合因果森林######
set.seed(42)
fit <- causal_forest(
	X = X,
	Y = Y,
	W = W,
	num.trees = 2000,
	mtry = 4,
	sample.fraction = 0.5,
	honesty = TRUE
)

######4.1.2 计算CATE#################
# 计算CATE
cate <- predict(fit)$predictions
data_clean$CATE <- cate

# CATE描述性统计
cat("CATE范围: [", round(min(cate), 4), ", ", round(max(cate), 4), "]\n")
cat("CATE均值:", round(mean(cate), 4), "\n")
cat("CATE中位数:", round(median(cate), 4), "\n")
cat("CATE标准差:", round(sd(cate), 4), "\n")
cat("CATE<0(获益)比例:", round(mean(cate < 0) * 100, 2), "%\n\n")

# CATE分布图
#简化版
hist(data_clean$CATE ) 
#ggplot法
p1 <- ggplot(data_clean, aes(x = CATE)) +
	geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
	geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
	labs(title = "个体处理效应(CATE)分布",
			 subtitle = "负值表示高身体活动降低肺功能异常风险",
			 x = "CATE (风险差)",
			 y = "频数") +
	theme_minimal()
print(p1)

###### 4.1.3 特征重要性#############
importance <- variable_importance(fit)
importance_df <- data.frame(
	Feature = covariates,
	Importance = importance
) %>% arrange(desc(Importance))

cat("特征重要性排名:\n")
print(importance_df)

# 特征重要性图
p2 <- ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
	geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
	coord_flip() +
	labs(title = "特征重要性 (决定疗效异质性的关键因素)",
			 x = "特征",
			 y = "重要性得分") +
	theme_minimal()
print(p2)

#####4.1.4 部分依赖图(PDP) 
library(hstats)
# 定义预测函数
pred_fun <- function(object, newdata) {
	predict(object, newdata)$predictions
}

# 第一种方法，利用hstats包绘制
# 绘制单个变量的PDP
pdp_age <- plot(partial_dep(fit, v = "age", X = X, pred_fun = pred_fun))

# 第二种方法，带置信区间#####
# 定义预测函数
pred_fun <- function(object, newdata) {
	predict(object, newdata)$predictions
}

# 选择最重要的两个特征绘制PDP
top_features <- importance_df$Feature[1:2]
cat("最重要的两个特征:", top_features[1], "和", top_features[2], "\n\n")

# 方法：手动计算带置信区间的PDP
# 函数：计算带置信区间的PDP
compute_pdp_with_ci <- function(fit, X_data, feature_name, grid_points = 50) {
	# 获取特征的最小值和最大值
	feature_min <- min(X_data[, feature_name])
	feature_max <- max(X_data[, feature_name])
	
	# 创建网格点
	grid_vals <- seq(feature_min, feature_max, length.out = grid_points)
	
	# 存储结果
	results <- data.frame(
		x = grid_vals,
		yhat_mean = NA,
		yhat_lower = NA,
		yhat_upper = NA
	)
	
	# 对每个网格点计算CATE的分布
	for(i in 1:length(grid_vals)) {
		# 创建新数据，固定当前特征为特定值
		X_new <- X_data
		X_new[, feature_name] <- grid_vals[i]
		
		# 预测CATE
		cate_pred <- predict(fit, X_new)$predictions
		
		# 计算均值、2.5%和97.5%分位数
		results$yhat_mean[i] <- mean(cate_pred)
		results$yhat_lower[i] <- quantile(cate_pred, 0.025)
		results$yhat_upper[i] <- quantile(cate_pred, 0.975)
	}
	
	return(results)
}

# 为第一个特征计算PDP（带置信区间）
cat("计算", top_features[1], "的PDP（带置信区间）...\n")
pdp_1_ci <- compute_pdp_with_ci(fit, X, top_features[1], grid_points = 50)

# 绘制第一个特征的PDP图（带置信区间）
p3 <- ggplot(pdp_1_ci, aes(x = x, y = yhat_mean)) +
	geom_line(color = "steelblue", size = 1.2) +
	geom_ribbon(aes(ymin = yhat_lower, ymax = yhat_upper), 
							alpha = 0.3, fill = "steelblue") +
	geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.8) +
	labs(title = paste("部分依赖图:", top_features[1]),
			 subtitle = paste(top_features[1], "对处理效应的影响（阴影区域为95%置信区间）"),
			 x = top_features[1],
			 y = "条件平均处理效应 (CATE)") +
	theme_minimal() +
	theme(plot.title = element_text(hjust = 0.5),
				plot.subtitle = element_text(hjust = 0.5))
print(p3)

# 为第二个特征计算PDP（带置信区间）
cat("\n计算", top_features[2], "的PDP（带置信区间）...\n")
pdp_2_ci <- compute_pdp_with_ci(fit, X, top_features[2], grid_points = 50)

# 绘制第二个特征的PDP图（带置信区间）
p4 <- ggplot(pdp_2_ci, aes(x = x, y = yhat_mean)) +
	geom_line(color = "steelblue", size = 1.2) +
	geom_ribbon(aes(ymin = yhat_lower, ymax = yhat_upper), 
							alpha = 0.3, fill = "steelblue") +
	geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.8) +
	labs(title = paste("部分依赖图:", top_features[2]),
			 subtitle = paste(top_features[2], "对处理效应的影响（阴影区域为95%置信区间）"),
			 x = top_features[2],
			 y = "条件平均处理效应 (CATE)") +
	theme_minimal() +
	theme(plot.title = element_text(hjust = 0.5),
				plot.subtitle = element_text(hjust = 0.5))
print(p4)

## 4.1.5 结合shap法######
# ####3提示，运行时间较长,不要轻易尝试#####
# 加载包
library(kernelshap)
library(shapviz)

# 定义预测函数：返回 CATE
pred_fun <- function(object, newdata) {
	predict(object, newdata)$predictions
}

# 可选：结合 SHAP （背景数据可采样以加速，这里取小样本 100 例进行演示）
# kernelshap 需要调用预测函数数千次（每行每个特征），所需时间非常长
set.seed(123)
bg_X <- X[sample(nrow(X), 100), ]

system.time({
	ks <- kernelshap(fit, X = X, bg_X = bg_X, pred_fun = pred_fun)
})

# 转为 shapviz 对象
shp <- shapviz(ks)

# 绘图
sv_importance(shp)               # 变量重要性条形图
sv_importance(shp, kind = "bee") # 蜂群图
sv_dependence(shp, v = "age")    # 依赖图（可换变量名）

### ##5. 因果森林的结合###########################
#5.1 利用因果森林的包grf 整合aipw和tmle方法#####
# 研究假设: 高身体活动是否能降低肺功能异常风险？
# 原假设H0: 高身体活动对肺功能无影响 (ATE = 0)
# 备择假设H1: 高身体活动对肺功能有影响 (ATE ≠ 0)

# 总的效应
# 平均处理效应（使用target.sample="control"避免警告）
ate <- average_treatment_effect(fit, target.sample = "all",method = c("AIPW"))

cat("=== 平均处理效应(ATE)检验结果 ===\n")
cat("ATE估计值:", round(ate[1], 4), "\n")
cat("标准误:", round(ate[2], 4), "\n")
cat("95%置信区间: [", round(ate[1] - 1.96*ate[2], 4), ", ", 
		round(ate[1] + 1.96*ate[2], 4), "]\n")
cat("Z统计量:", round(ate[1]/ate[2], 4), "\n")
cat("P值:", format.pval(2 * pnorm(-abs(ate[1]/ate[2])), eps = 0.001), "\n\n")

###5.2 亚组分析######
#### 5.2.1 CATE五等分亚组分析###################
data_clean$CATE_group <- cut(data_clean$CATE,
														 breaks = quantile(data_clean$CATE, seq(0, 1, 0.2)),
														 labels = c("Q1(疗效最差)", "Q2", "Q3", "Q4", "Q5(疗效最好)"),
														 include.lowest = TRUE)

cat("各分位数组样本量:\n")
print(table(data_clean$CATE_group))

# 计算各组的ATE
cat("\n各分位数组的ATE结果:\n")
for(group in levels(data_clean$CATE_group)) {
	idx <- which(data_clean$CATE_group == group)
	if(length(idx) > 0) {
		ate_group <- average_treatment_effect(fit, subset = idx, target.sample = "all",method = c("AIPW"))
		p_val <- 2 * pnorm(-abs(ate_group[1]/ate_group[2]))
		
		cat("\n", group, ":\n")
		cat("  样本量:", length(idx), "\n")
		cat("  平均CATE:", round(mean(data_clean$CATE[idx]), 4), "\n")
		cat("  ATE:", round(ate_group[1], 4), "\n")
		cat("  95% CI: [", round(ate_group[1] - 1.96*ate_group[2], 4), ", ",
				round(ate_group[1] + 1.96*ate_group[2], 4), "]\n")
		cat("  P值:", format.pval(p_val, eps = 0.001), 
				ifelse(p_val < 0.05, " (显著)", " (不显著)"), "\n")
	}
}


#### 5.2.3 基于重要特征的亚组分析######

top1 <- importance_df$Feature[1]  # age
top2 <- importance_df$Feature[2]  # bmi

# 使用中位数分割
cut1 <- median(data_clean[[top1]])
cut2 <- median(data_clean[[top2]])

data_clean$group1 <- ifelse(data_clean[[top1]] <= cut1, 
														paste0(top1, "低(≤", round(cut1,1), ")"), 
														paste0(top1, "高(>", round(cut1,1), ")"))
data_clean$group2 <- ifelse(data_clean[[top2]] <= cut2, 
														paste0(top2, "低(≤", round(cut2,1), ")"), 
														paste0(top2, "高(>", round(cut2,1), ")"))
data_clean$subgroup <- paste(data_clean$group1, data_clean$group2, sep = " & ")

cat("亚组划分依据:\n")
cat("  ", top1, "中位数:", round(cut1, 2), "\n")
cat("  ", top2, "中位数:", round(cut2, 2), "\n\n")

for(subg in unique(data_clean$subgroup)) {
	idx <- which(data_clean$subgroup == subg)
	if(length(idx) > 0) {
		ate_sub <- average_treatment_effect(fit, subset = idx,target.sample = "all",method = c("AIPW"))
		p_val <- 2 * pnorm(-abs(ate_sub[1]/ate_sub[2]))
		
		cat("\n", subg, ":\n")
		cat("  样本量:", length(idx), "\n")
		cat("  平均CATE:", round(mean(data_clean$CATE[idx]), 4), "\n")
		cat("  ATE:", round(ate_sub[1], 4), "\n")
		cat("  95% CI: [", round(ate_sub[1] - 1.96*ate_sub[2], 4), ", ",
				round(ate_sub[1] + 1.96*ate_sub[2], 4), "]\n")
		cat("  P值:", format.pval(p_val, eps = 0.001),
				ifelse(p_val < 0.05, " (显著)", " (不显著)"), "\n")
	}
}


# 保存结果
# write.csv(data_clean, "analysis_results.csv", row.names = FALSE)
# write.csv(importance_df, "feature_importance.csv", row.names = FALSE)

### 5.3 因果森林 + TMLE 组合###############
library(grf)
library(tmle)
library(SuperLearner)
library(dplyr)

# 设置工作目录
setwd("D:/因果训练营/异质性分析课程")

# 读取数据
data <- read.csv("data1.csv")

# 定义变量
covariates <- c("sex","hypertension","diabetes","heart_disease","residence","marital",
								"smoking","age","education_level","cardio_metabolic",
								"alcohol","bmi")

categorical_vars <- c("sex","hypertension","diabetes","heart_disease","residence","marital",
											"smoking","education_level","cardio_metabolic","alcohol")

# 转换分类变量为因子
data[categorical_vars] <- lapply(data[categorical_vars], factor)

# 结局变量
outcome_var <- "lungfun"

# 处理缺失值
data_complete <- data[complete.cases(data[, c(outcome_var, "PA_CAT", covariates)]), ]

# 准备数值型数据（GRF需要）
X <- as.data.frame(lapply(data_complete[, covariates], function(x) {
	if(is.factor(x)) as.numeric(x) else x
}))
W <- data_complete$PA_CAT
Y <- data_complete[[outcome_var]]

# 确保Y是0/1数值型
Y <- as.numeric(Y)
if(!all(Y %in% c(0, 1))) {
	Y <- ifelse(Y == min(Y), 0, 1)
}

# 步骤1：因果森林识别重要变量 
set.seed(123)
cf <- causal_forest(X, Y, W)
var_imp <- variable_importance(cf)
names(var_imp) <- covariates
top_features <- names(sort(var_imp, decreasing = TRUE))[1:2]
cat("\n重要效应修饰因子：", top_features, "\n")

# 步骤2：预测 CATE（用于后续亚组划分）
cate_all <- predict(cf)$predictions
data_clean <- data_complete
data_clean$CATE <- cate_all

# 步骤3：定义提取 TMLE 6 个指标的辅助函数
extract_tmle_6 <- function(tmle_fit) {
	effects <- c("EY0", "EY1", "ATE", "ATT", "ATC", "RR", "OR")
	res <- data.frame()
	for (e in effects) {
		if (!is.null(tmle_fit$estimates[[e]])) {
			ci <- tmle_fit$estimates[[e]]$CI
			res <- rbind(res, data.frame(
				Effect = e,
				Estimate = tmle_fit$estimates[[e]]$psi,
				CI_lower = ci[1],
				CI_upper = ci[2],
				pvalue = tmle_fit$estimates[[e]]$pvalue
			))
		}
	}
	return(res)
}


# 5.3.1 全数据集上的 TMLE（完整 6 个指标）############
# 为减少运行时间，这里仅运行一种
SL_library <- c("SL.glm")#, "SL.glmnet", "SL.randomForest", "SL.xgboost", "SL.gam"

tmle_full <- tmle(Y = Y, A = W, W = X,
									Q.SL.library = SL_library,
									g.SL.library = SL_library,
									family = "binomial", V.Q = 5, V.g = 5, verbose = FALSE)

full_results <- extract_tmle_6(tmle_full)
full_results$Subgroup <- "全数据集"
full_results$N <- length(Y)
print(full_results[, c("Subgroup", "N", "Effect", "Estimate", "CI_lower", "CI_upper", "pvalue")], row.names = FALSE)


# 5.3.2 CATE 五等分亚组 TMLE############
# 创建五等分分组
data_clean$CATE_group <- cut(data_clean$CATE,
														 breaks = quantile(data_clean$CATE, seq(0,1,0.2)),
														 labels = c("Q1(疗效最差)", "Q2", "Q3", "Q4", "Q5(疗效最好)"),
														 include.lowest = TRUE)

cate_groups <- levels(data_clean$CATE_group)
all_cate_results <- lapply(cate_groups, function(grp) {
	idx <- which(data_clean$CATE_group == grp)
	if(length(idx) == 0) return(NULL)
	Y_sub <- Y[idx]; A_sub <- W[idx]; X_sub <- X[idx, , drop = FALSE]
	
	tmle_sub <- tryCatch({
		tmle(Y = Y_sub, A = A_sub, W = X_sub,
				 Q.SL.library = SL_library, g.SL.library = SL_library,
				 family = "binomial", V.Q = 5, V.g = 5, verbose = FALSE)
	}, error = function(e) { cat("Error in", grp, ":", e$message, "\n"); return(NULL) })
	
	if(is.null(tmle_sub)) return(NULL)
	res <- extract_tmle_6(tmle_sub)
	res$Subgroup <- grp
	res$N <- length(idx)
	return(res)
})
all_cate_results <- do.call(rbind, all_cate_results[!sapply(all_cate_results, is.null)])
print(all_cate_results[, c("Subgroup", "N", "Effect", "Estimate", "CI_lower", "CI_upper", "pvalue")], row.names = FALSE)

# 5.3.3年龄亚组 + BMI 亚组的 TMLE 分析（分别独立）################
library(tmle)
library(SuperLearner)


# 定义提取 7 个指标的函数
extract_tmle_7 <- function(tmle_fit) {
	effects <- c("EY0", "EY1", "ATE", "ATT", "ATC", "RR", "OR")
	res <- data.frame()
	for (e in effects) {
		if (!is.null(tmle_fit$estimates[[e]])) {
			ci <- tmle_fit$estimates[[e]]$CI
			res <- rbind(res, data.frame(
				Effect = e,
				Estimate = tmle_fit$estimates[[e]]$psi,
				CI_lower = ci[1],
				CI_upper = ci[2],
				pvalue = tmle_fit$estimates[[e]]$pvalue
			))
		}
	}
	return(res)
}

# 定义运行单个亚组 TMLE 的函数
run_tmle_subgroup <- function(idx, subgroup_name) {
	if(length(idx) == 0) {
		cat("亚组", subgroup_name, "样本量为0，跳过\n")
		return(NULL)
	}
	Y_sub <- Y[idx]
	A_sub <- W[idx]
	X_sub <- X[idx, , drop = FALSE]
	
	tmle_fit <- tryCatch({
		tmle(Y = Y_sub, A = A_sub, W = X_sub,
				 Q.SL.library = SL_library,
				 g.SL.library = SL_library,
				 family = "binomial",
				 V.Q = 5, V.g = 5, verbose = FALSE)
	}, error = function(e) {
		cat("TMLE 失败 (", subgroup_name, "): ", e$message, "\n")
		return(NULL)
	})
	if(is.null(tmle_fit)) return(NULL)
	res <- extract_tmle_7(tmle_fit)
	res$Subgroup <- subgroup_name
	res$N <- length(idx)
	return(res)
}

# ===================== 1. 年龄亚组 =====================
cat("\n========== 年龄亚组 TMLE 结果 ==========\n")
age_median <- median(data_clean$age, na.rm = TRUE)
idx_age_low <- which(data_clean$age <= age_median)
idx_age_high <- which(data_clean$age > age_median)

age_results <- lapply(list(
	list(idx = idx_age_low, name = paste0("age低 (≤", round(age_median,1), ")")),
	list(idx = idx_age_high, name = paste0("age高 (>", round(age_median,1), ")"))
), function(x) run_tmle_subgroup(x$idx, x$name))

age_results <- do.call(rbind, age_results[!sapply(age_results, is.null)])
if(!is.null(age_results) && nrow(age_results) > 0) {
	print(age_results[, c("Subgroup", "N", "Effect", "Estimate", "CI_lower", "CI_upper", "pvalue")], row.names = FALSE)
} else {
	cat("年龄亚组无有效结果\n")
}

# ===================== 2. BMI 亚组 =====================
cat("\n========== BMI 亚组 TMLE 结果 ==========\n")
bmi_median <- median(data_clean$bmi, na.rm = TRUE)
idx_bmi_low <- which(data_clean$bmi <= bmi_median)
idx_bmi_high <- which(data_clean$bmi > bmi_median)

bmi_results <- lapply(list(
	list(idx = idx_bmi_low, name = paste0("bmi低 (≤", round(bmi_median,1), ")")),
	list(idx = idx_bmi_high, name = paste0("bmi高 (>", round(bmi_median,1), ")"))
), function(x) run_tmle_subgroup(x$idx, x$name))

bmi_results <- do.call(rbind, bmi_results[!sapply(bmi_results, is.null)])
if(!is.null(bmi_results) && nrow(bmi_results) > 0) {
	print(bmi_results[, c("Subgroup", "N", "Effect", "Estimate", "CI_lower", "CI_upper", "pvalue")], row.names = FALSE)
} else {
	cat("BMI亚组无有效结果\n")
}


### 5.4 因果森林 + CausalForestDynamicSubgroups 组合，结合AIPW算法#####
library(grf)
library(EpiForsk)
library(dplyr)
library(ggplot2)
# 设置工作目录
setwd("D:/因果训练营/异质性分析课程")

# 读取数据
data <- read.csv("data1.csv")

# 定义变量
covariates <- c("sex","hypertension","diabetes","heart_disease","residence","marital",
								"smoking","age","education_level","cardio_metabolic",
								"alcohol","bmi")

categorical_vars <- c("sex","hypertension","diabetes","heart_disease","residence","marital",
											"smoking","education_level","cardio_metabolic","alcohol")

# 转换分类变量为因子
data[categorical_vars] <- lapply(data[categorical_vars], factor)

# 结局变量
outcome_var <- "lungfun"

# 处理缺失值
data_complete <- data[complete.cases(data[, c(outcome_var, "PA_CAT", covariates)]), ]

# 准备数值型数据（GRF需要）
X <- as.data.frame(lapply(data_complete[, covariates], function(x) {
	if(is.factor(x)) as.numeric(x) else x
}))
W <- data_complete$PA_CAT
Y <- data_complete[[outcome_var]]

# 确保Y是0/1数值型
Y <- as.numeric(Y)
if(!all(Y %in% c(0, 1))) {
	Y <- ifelse(Y == min(Y), 0, 1)
}

# 步骤1：拟合因果森林
set.seed(123)
cf <- causal_forest(X, Y, W)

# 步骤2：CausalForestDynamicSubgroups 
cf_ds <- CausalForestDynamicSubgroups(
	forest = cf,
	n_rankings = 3,    # 分为3组（高、中、低效应）
	n_folds = 5        # 5折交叉验证
)

# 输出结果 
cat("\n========== 各亚组ATE估计 ==========\n")
print(cf_ds$forest_rank_ate)

cat("\n========== 亚组间差异检验 ==========\n")
print(cf_ds$forest_rank_diff_test)

cat("\n========== 个体CATE及亚组归属（前20行） ==========\n")
print(head(cf_ds$forest_subgroups, 20))

# 可视化
# 森林图
print(cf_ds$forest_rank_ate_plot)

# 热力图（协变量分布）
print(cf_ds$heatmap)


