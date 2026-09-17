
##第一部分  准备工作

################################################################
# 一、加载所需要的包
################################################################

#1.下载需要的包
if(!require("dplyr")) install.packages("dplyr")
if(!require("survival")) install.packages("survival")
if(!require("autoReg")) install.packages("autoReg")
if(!require("rrtable")) install.packages("rrtable")#导出word
if(!require("boot")) install.packages("boot")
if(!require("stdReg2")) install.packages("stdReg2")
if(!require("RISCA")) install.packages("RISCA")
if(!require("AIPW")) install.packages("AIPW")
if(!require("SuperLearner")) install.packages("SuperLearner")

#2.加载需要的包
library(dplyr)
library(survival)
library(autoReg)
library(rrtable)#导出word
library(boot)
library(stdReg2)
library(RISCA)
library(AIPW)
library(SuperLearner)
library(ggplot2)

################################################################
# 二、数据导入
################################################################

# 1.获取工作空间
getwd()
# 2.设置工作空间
# setwd("D:/BaiduSyncdisk/讲座/2025年讲座/因果训练营")

# 3.导入数据data,该数据在之前已经被整理成为data
data<- read.csv("dataAIPW.csv")
# 4.查看数据
str(data)

# 第二部分  AIPW方法
################################################################
#一、普通AIPW方法（二分类结局）
################################################################
#1. 设置协变量
covariates_set <- c("年龄", "性别", "教育年份", "种族", "收入",
                    "体重","DNR状态", "医疗保险", "癌症", "心血管疾病史", 
                    "充血性心力衰竭", "痴呆", "精神病史", "肺病", "肾病",
                    "肝病", "上消化道出血","肿瘤", "免疫抑制", "转院", 
                    "确定性心肌梗死", "DASI指数", "APACHE评分", "格拉斯哥昏迷评分","平均血压",
                    "白细胞计数", "心率", "呼吸频率","体温","氧合指数", 
                    "白蛋白", "红细胞比容","胆红素", "肌酐","钠", 
                    "钾", "动脉血二氧化碳分压", "血清PH值")

vars<-c("性别", "种族", "收入","DNR状态",
        "医疗保险", "癌症", "心血管疾病史", 
                    "充血性心力衰竭", "痴呆", "精神病史", "肺病", "肾病",
                    "肝病", "上消化道出血","肿瘤", "免疫抑制", "转院", 
                    "确定性心肌梗死")

data[vars]<-lapply(data[vars],as.factor)
str(data)

#2.构建普通AIPW模型
aipw_simple <- AIPW$new(Y = data$死亡,
                        A = data$右心导管术,
                        W = data[, covariates_set], # 这会创建一个只含协变量的数据框
                        Q.SL.library = "SL.glm",
                        g.SL.library = "SL.glm",
                        k_split = 1,
                        verbose = TRUE)
# 3. 拟合模型
aipw_simple$fit()

# 4. 查看结果
aipw_simple$summary(g.bound = 0.025) 

#5.查看倾向得分
aipw_simple$plot.p_score()
aipw_simple$plot.ip_weights()


################################################################
#二、AIPW+SL方法（二分类结局）
################################################################
#1. 同样设置协变量
covariates_set <- c("年龄", "性别", "教育年份", "种族", "收入",
                    "体重","DNR状态", "医疗保险", "癌症", "心血管疾病史", 
                    "充血性心力衰竭", "痴呆", "精神病史", "肺病", "肾病",
                    "肝病", "上消化道出血","肿瘤", "免疫抑制", "转院", 
                    "确定性心肌梗死", "DASI指数", "APACHE评分", "格拉斯哥昏迷评分","平均血压",
                    "白细胞计数", "心率", "呼吸频率","体温","氧合指数", 
                    "白蛋白", "红细胞比容","胆红素", "肌酐","钠", 
                    "钾", "动脉血二氧化碳分压", "血清PH值")


#2.构建超级学习者算法
# 创建AIPW对象并拟合
aipw_super <- AIPW$new(Y = data$死亡,
                       A = data$右心导管术,
                       W = data[, covariates_set],# 这会创建一个只含协变量的数据框
                       # 超级学习器库：集成多种算法预测结局Y
                       # 注意：结局是二分类，所有算法需兼容binomial
                       Q.SL.library = c("SL.glm",        # 广义线性模型（逻辑回归）
                                        "SL.glmnet",     # 弹性网络（带正则化）
                                        "SL.ranger",     # 随机森林
                                        "SL.mean" ),
                                             # 均值模型（作为基准）
                       # 超级学习器库：集成多种算法预测暴露A
                       g.SL.library = c("SL.glm",
                                        "SL.glmnet",
                                        "SL.ranger",
                                        "SL.mean" ) ,
                                      
                       # 关键：使用5折交叉拟合以消除过拟合偏误
                       k_split = 5,
                       # 是否在控制台显示详细拟合过程（第一次运行建议设为TRUE）
                       verbose = TRUE)

# 拟合模型（这可能需要一些时间，取决于数据量和算法复杂度）
aipw_super$fit()


# 5. 查看结果

aipw_super$summary(g.bound = 0.025) 
print(aipw_super$result)

aipw_super$plot.p_score()
aipw_super$plot.ip_weights()



#  第三部分 双重机器学习
################################################################
#一、双重机器学习（二分类结局）
################################################################
if(!require("DoubleML")) install.packages("DoubleML")
if(!require("mlr3")) install.packages("mlr3")

library(DoubleML)
library(mlr3)

covariates_set <- c("年龄", "性别", "教育年份", "种族", "收入",
                    "体重","DNR状态", "医疗保险", "癌症", "心血管疾病史", 
                    "充血性心力衰竭", "痴呆", "精神病史", "肺病", "肾病",
                    "肝病", "上消化道出血","肿瘤", "免疫抑制", "转院", 
                    "确定性心肌梗死", "DASI指数", "APACHE评分", "格拉斯哥昏迷评分","平均血压",
                    "白细胞计数", "心率", "呼吸频率","体温","氧合指数", 
                    "白蛋白", "红细胞比容","胆红素", "肌酐","钠", 
                    "钾", "动脉血二氧化碳分压", "血清PH值")


#定义数据
dml_data <- DoubleMLData$new(
  data = data,
  y_col = "死亡",
  d_cols = "右心导管术",
  x_cols = covariates_set  # 这里使用你正面定义的协变量列表
)

#定义分类器（以随机森林为例）

learner_classif <- lrn("classif.ranger", num.trees = 500, predict_type = "prob")
learner_regr <- lrn("regr.ranger", num.trees = 500)


dml_plr <- DoubleMLIRM$new(
  data = dml_data,
  ml_g = learner_regr,   # 用于E[Y|X]，回归学习器
  ml_m = learner_classif, # 用于E[D|X]，分类学习器 (必须)
  n_folds = 5,           # 交叉拟合折数
  n_rep = 3,             # 重复交叉拟合次数，增加稳定性
  score = "ATE"          # 估计平均处理效应
)

dml_plr$fit()

dml_plr$summary()

print(dml_plr$confint())

###############################################################

