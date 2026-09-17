
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

#2.加载需要的包
library(dplyr)
library(survival)
library(autoReg)
library(rrtable)#导出word
library(boot)
library(stdReg2)
library(RISCA)

################################################################
# 二、数据导入
################################################################

# 1.获取工作空间
getwd()
# 2.设置工作空间
# setwd("C:/Users/fengbaotongji/Desktop/因果推断训练营")
setwd("E:/Rpackage/causalR/learn/courses/zheng-w02-regression-gformula")

# 3.导入数据
Obsdata<- read.csv("ObsFinal.csv")
# 4.查看数据
str(Obsdata)


################################################################
# 三、数据预处理和变量编码
################################################################
# 1.结局变量
Obsdata$死亡_180天 <- ifelse(Obsdata$死亡_180天 == "Yes", 1, 0)  #赋值需转换为0,1
Obsdata$死亡_180天 <- as.numeric(Obsdata$死亡_180天)
Obsdata$生存时间_180天 <- as.numeric(Obsdata$生存时间_180天)

# 2.处理变量
Obsdata$右心导管术 <- as.numeric(Obsdata$右心导管术)

# 3.分类变量因子化并创建标签
Obsdata$癌症 <- factor(Obsdata$癌症,levels = c("No","Yes","Metastatic"),
                  labels = c("无", "有", "转移性"))

Obsdata$性别 <- factor(Obsdata$性别,levels = c("Male","Female"),
                  labels = c("男","女" ))
Obsdata$种族 <- factor(Obsdata$种族,levels = c("white","black","other"),
                  labels = c("白人", "黑人", "其他"))
Obsdata$DNR状态 <- factor(Obsdata$DNR状态,levels = c("No", "Yes"),
                     labels = c("否", "是"))
Obsdata$医疗保险 <- factor(Obsdata$医疗保险,
                    levels = c("No insurance", "Private", "Medicare",
                               "Medicaid","Private & Medicare",
                               "Medicare & Medicaid"),
                    labels = c("没有保险", "私人", "医疗保险",
                               "医疗补助", "私人与医疗保险",
                               "医疗保险与医疗补助"))
Obsdata$收入 <- factor(Obsdata$收入,
                  levels = c("Under $11k","$11-$25k","$25-$50k","> $50k"),
                  labels = c("不到1.1万美元", "1.1-2.5万美元",
                             "2.5-5万美元", ">5万美元"))

Obsdata$心血管疾病史 <- factor(Obsdata$心血管疾病史,levels = c(0, 1),
                      labels = c("否", "是"))

Obsdata$充血性心力衰竭 <- factor(Obsdata$充血性心力衰竭,levels = c(0, 1),
                      labels = c("否", "是"))

Obsdata$痴呆 <- factor(Obsdata$痴呆,levels = c(0, 1),
                      labels = c("否", "是"))

Obsdata$精神病史 <- factor(Obsdata$精神病史,levels = c(0, 1),
                      labels = c("否", "是"))

Obsdata$肺病 <- factor(Obsdata$肺病,levels = c(0, 1),
                      labels = c("否", "是"))

Obsdata$肾病 <- factor(Obsdata$肾病,levels = c(0, 1),
                      labels = c("否", "是"))

Obsdata$肝病 <- factor(Obsdata$肝病,levels = c(0, 1),
                      labels = c("否", "是"))

Obsdata$上消化道出血 <- factor(Obsdata$上消化道出血,levels = c(0, 1),
                      labels = c("否", "是"))

Obsdata$肿瘤 <- factor(Obsdata$肿瘤,levels = c(0, 1),
                      labels = c("否", "是"))

Obsdata$免疫抑制 <- factor(Obsdata$免疫抑制,levels = c(0, 1),
                      labels = c("否", "是"))

Obsdata$转院 <- factor(Obsdata$转院,levels = c(0, 1),
                      labels = c("否", "是"))

Obsdata$确定性心肌梗死 <- factor(Obsdata$确定性心肌梗死,levels = c(0, 1),
                      labels = c("否", "是"))


# 4.转换变量
Obsdata$年龄_分类 <- cut(Obsdata$年龄,breaks=c(-Inf, 20, 40, 60, 80, Inf),label=F)
str(Obsdata)

# 5.定量变量转换为数值型
#产生变量字符串
quantitative_vars <- c("年龄", "教育年份", "DASI指数", "APACHE评分", "格拉斯哥昏迷评分",
                       "平均血压", "白细胞计数", "心率", "呼吸频率", "体温", "氧合指数",
                       "白蛋白", "红细胞比容", "胆红素", "肌酐", "钠", "钾",
                       "动脉血二氧化碳分压", "血清PH值", "体重")
Obsdata[quantitative_vars] <- lapply(Obsdata[quantitative_vars], as.numeric)
str(Obsdata)

# 6.提取数据集用于后续分析
data0 <- dplyr::select(Obsdata,
                          c("死亡_180天","生存时间_180天","右心导管术",
                            "年龄", "性别", "教育年份", "种族", "收入",
                            "体重","DNR状态", "医疗保险", "癌症", "心血管疾病史",
                            "充血性心力衰竭", "痴呆", "精神病史", "肺病", "肾病",
                            "肝病", "上消化道出血","肿瘤", "免疫抑制", "转院",
                            "确定性心肌梗死", "DASI指数", "APACHE评分", "格拉斯哥昏迷评分","平均血压",
                            "白细胞计数", "心率", "呼吸频率","体温","氧合指数",
                            "白蛋白", "红细胞比容","胆红素", "肌酐","钠",
                            "钾", "动脉血二氧化碳分压", "血清PH值"))

data<- dplyr::rename(data0, 死亡 = 死亡_180天, 生存时间 = 生存时间_180天) #重命名
str(data)


# 第二部分  回归方法控制混杂

################################################################
#一、单因素+多因素logistic回归控制混杂（二分类结局）
################################################################

#1.构建logistic回归模型
death.log <- glm(死亡 ~ 右心导管术+年龄+性别+教育年份+种族+收入+
                   体重+DNR状态+医疗保险+癌症+心血管疾病史+
                   充血性心力衰竭+痴呆+精神病史+肺病+肾病+
                   肝病+上消化道出血+肿瘤+免疫抑制+转院+
                   确定性心肌梗死+DASI指数+APACHE评分+格拉斯哥昏迷评分+平均血压+
                   白细胞计数+心率+呼吸频率+体温+氧合指数+
                   白蛋白+红细胞比容+胆红素+肌酐+钠+
                   钾+动脉血二氧化碳分压+血清PH值,
                 data=data,
                 family=binomial(link = "logit"))

#2.查看回归结果
summary(death.log)

#3.autoReg包的使用
logmodel<-autoReg(death.log,uni=TRUE,multi=TRUE,threshold=0.05)
logmodel
#uni为TRUE指输出单因素模型结果，multi为TRUE输出多因素模型结果，threshold纳入条件

#4.导出到docx，可编辑数据
table2docx(logmodel,"单因素+多因素logistic回归控制混杂")



################################################################
#二、单因素+多因素cox回归控制混杂（生存资料结局）
################################################################

#1.构建cox风险评估模型
death.cox<-coxph(Surv(data$生存时间, data$死亡) ~
                   右心导管术+年龄+性别+教育年份+种族+收入+
                   体重+DNR状态+医疗保险+癌症+心血管疾病史+
                   充血性心力衰竭+痴呆+精神病史+肺病+肾病+
                   肝病+上消化道出血+肿瘤+免疫抑制+转院+
                   确定性心肌梗死+DASI指数+APACHE评分+格拉斯哥昏迷评分+平均血压+
                   白细胞计数+心率+呼吸频率+体温+氧合指数+
                   白蛋白+红细胞比容+胆红素+肌酐+钠+
                   钾+动脉血二氧化碳分压+血清PH值,
                 data = data)
#2.查看回归结果
summary(death.cox)

#3.autoReg包的使用
coxmodel<-autoReg(death.cox,uni=TRUE,multi=TRUE,threshold=0.05)
coxmodel

#4.导出到docx，可编辑数据
table2docx(coxmodel,"单因素+多因素cox控制混杂")



################################################################
###三、logistic回归多模型策略控制混杂（二分类结局）
################################################################
#建立多个模型，逐步调整协变量

#1.模型1，不调整协变量
death.log1 <- glm(死亡 ~ 右心导管术,data=data,family=binomial(link = "logit"))
summary(death.log1)
logmodel1<-autoReg(death.log1,uni=TRUE,multi=FALSE,threshold=1)#只显示单因素。
logmodel1
table2docx(logmodel1,"多模型策略logistic回归model1")

#2.模型2，调整年龄和性别变量
death.log2 <- glm(死亡 ~ 右心导管术+年龄+性别,data=data,family=binomial(link = "logit"))
summary(death.log2)
logmodel2<-autoReg(death.log2,uni=FALSE,multi=TRUE,threshold=1)
logmodel2
table2docx(logmodel2,"多模型策略logistic回归model2")

#3.模型3，调整所有变量
death.log3 <- glm(死亡 ~ 右心导管术+年龄+性别+教育年份+种族+收入+
                    体重+DNR状态+医疗保险+癌症+心血管疾病史+
                    充血性心力衰竭+痴呆+精神病史+肺病+肾病+
                    肝病+上消化道出血+肿瘤+免疫抑制+转院+
                    确定性心肌梗死+DASI指数+APACHE评分+格拉斯哥昏迷评分+平均血压+
                    白细胞计数+心率+呼吸频率+体温+氧合指数+
                    白蛋白+红细胞比容+胆红素+肌酐+钠+
                    钾+动脉血二氧化碳分压+血清PH值,
                  data=data,
                  family=binomial(link = "logit"))
summary(death.log3)
logmodel3<-autoReg(death.log3,uni=FALSE,multi=TRUE,threshold=1)
logmodel3
table2docx(logmodel3,"多模型策略logistic回归model3")


################################################################
#四、cox回归多模型策略控制混杂（生存资料结局）
################################################################
#建立多个模型，逐步调整协变量

#1.模型1，不调整协变量
death.cox1 <- coxph(Surv(data$生存时间, data$死亡) ~ 右心导管术,data = data)
summary(death.cox1)
coxmodel1<-autoReg(death.cox1,uni=TRUE,multi=FALSE,threshold=1)#只显示单因素。
coxmodel1
table2docx(coxmodel1,"多模型策略cox回归model1")

#2.模型2，调整年龄和性别变量
death.cox2 <- coxph(Surv(data$生存时间, data$死亡) ~ 右心导管术+年龄+性别,data = data)
summary(death.cox2)
coxmodel2<-autoReg(death.cox2,uni=FALSE,multi=TRUE,threshold=1)
coxmodel2
table2docx(coxmodel2,"多模型策略cox回归model2")

#3.模型3，调整所有变量
death.cox3 <- coxph(Surv(data$生存时间, data$死亡) ~
                      右心导管术+年龄+性别+教育年份+种族+收入+
                      体重+DNR状态+医疗保险+癌症+心血管疾病史+
                      充血性心力衰竭+痴呆+精神病史+肺病+肾病+
                      肝病+上消化道出血+肿瘤+免疫抑制+转院+
                      确定性心肌梗死+DASI指数+APACHE评分+格拉斯哥昏迷评分+平均血压+
                      白细胞计数+心率+呼吸频率+体温+氧合指数+
                      白蛋白+红细胞比容+胆红素+肌酐+钠+
                      钾+动脉血二氧化碳分压+血清PH值,
                    data = data)
summary(death.cox3)
coxmodel3<-autoReg(death.cox3,uni=FALSE,multi=TRUE,threshold=1)
coxmodel3
table2docx(coxmodel3,"多模型策略cox回归model3")




#### 第三部分：G计算因果推断

################################################################
#一、G计算开展因果推断——（手工5步法，以二分类结局为例）
################################################################
# 0.G计算协变量筛选：DAG+单因素分析 —— 筛选进入G计算的协变量
# 首先，定义协变量集合
covariates_set <- c("年龄", "性别", "教育年份", "种族", "收入",
                    "体重","DNR状态", "医疗保险", "癌症", "心血管疾病史",
                    "充血性心力衰竭", "痴呆", "精神病史", "肺病", "肾病",
                    "肝病", "上消化道出血","肿瘤", "免疫抑制", "转院",
                    "确定性心肌梗死", "DASI指数", "APACHE评分", "格拉斯哥昏迷评分","平均血压",
                    "白细胞计数", "心率", "呼吸频率","体温","氧合指数",
                    "白蛋白", "红细胞比容","胆红素", "肌酐","钠",
                    "钾", "动脉血二氧化碳分压", "血清PH值")
# 构建公式
formula_str <- paste("死亡 ~", paste(covariates_set, collapse = " + "))
# 展示公式
cat("使用的公式:\n", formula_str, "\n")

#单因素logistic回归筛选变量
death.log <- glm(as.formula(formula_str), data=data,family=binomial(link = "logit"))
summary(death.log)
logmodel<-autoReg(death.log,uni=TRUE,multi=T,threshold=0.05)#只显示单因素。
logmodel

# 1.用暴露变量、协变量和结局拟合回归模型
# 定义变量集合，进纳入单因素有意义的协变量
covariates_g<- c("年龄", "收入","体重","DNR状态", "医疗保险",
                 "癌症", "心血管疾病史", "充血性心力衰竭", "痴呆", "精神病史",
                 "肝病", "上消化道出血","肿瘤", "免疫抑制","DASI指数",
                 "APACHE评分", "格拉斯哥昏迷评分","平均血压","白细胞计数", "体温",
                 "白蛋白", "红细胞比容","胆红素", "肌酐","钾",
                 "动脉血二氧化碳分压")

# 构建公式
formula_g <- paste("死亡 ~ 右心导管术 +", paste(covariates_g, collapse = " + "))
# 展示公式
cat("使用的公式:\n", formula_g, "\n")

# 拟合回归模型
mod <- glm(as.formula(formula_g),
           family=binomial(link = "logit"),
           data=data)

# 2.创建反事实人群
# 所有受试者右心导管术=1（即接受右心导管术）
newdata1 = data.frame( 右心导管术 = 1, # 强制所有观测的处理变量 右心导管术 = 1（反事实设定）
                      # 保持原始数据中除 右心导管术 外的所有协变量不变
                      # 使用 dplyr::select() 选择 ObsData 中除右心导管术列外的所有列
                      dplyr::select(data, !右心导管术))


# 所有受试者右心导管术=0（即未接受右心导管术）
newdata0 = data.frame( 右心导管术 = 0,  # 强制所有观测的处理变量 右心导管术 = 0（反事实设定）
                      # 保持原始数据中除 右心导管术 外的所有协变量不变
                      # 使用 dplyr::select() 选择 ObsData 中除右心导管术列外的所有列
                      dplyr::select(data, !右心导管术))


# 3.以第一步中拟合回归模型mod为基础，在反事实人群中，预测疗效结局
# 在强制设置条件为所有受试者右心导管术=1（即接受右心导管术），预测疗效结局（Pred.Y1）
data$Pred.Y1 <- predict(mod,
                           newdata = newdata1,
                           type ="response")
Y.1 <- data$Pred.Y1

# 在强制设置条件为所有受试者右心导管术=0（即未接受右心导管术），预测疗效结局（Pred.Y0）
data$Pred.Y0 <- predict(mod,
                        newdata = newdata0,
                        type ="response")
Y.0 <- data$Pred.Y0


# 4.因果效应估计ATE
# 计算风险差 RD
RD_ATE <- mean((Y.1) - (Y.0), na.rm = TRUE)
RD_ATE

# 计算RR (相对危险度度)
RR_ATE <- mean(Y.1, na.rm = TRUE) / mean(Y.0, na.rm = TRUE)
RR_ATE

# 5. 使用R包boot，估计g计算的Bootstrap置信区间
# 估计ATE (风险差 RD)95%置信区间
gcomp.boot.RD <-function(formula = mod, data = data, indices) # 定义估算ATE的函数
  {
  boot_sample <- data[indices, ]
  fit.boot <- glm(formula, data = boot_sample)
  Pred.Y1 <- predict(fit.boot,
                     newdata = newdata1,
                     type ="response")
  Pred.Y0 <- predict(fit.boot,
                     newdata = newdata0,
                     type ="response")
  Pred.TE <- mean(Pred.Y1, na.rm = TRUE) - mean(Pred.Y0, na.rm = TRUE)
  return(Pred.TE)
}

set.seed(123)     # 设置随机种子以确保结果可重现

gcomp.RD <- boot(data=data,
                  statistic=gcomp.boot.RD,
                  R=250,
                  formula=mod)# 设置参数：R=250,生成250个bootstrap样本
plot(gcomp.RD)

RD.CI1 <- boot.ci(gcomp.RD, type="norm",conf=0.95)  # 基于正态近似的95%置信区间
RD.CI1
RD.CI2 <- boot.ci(gcomp.RD, type="perc",conf=0.95)  # 基于百分位数的95%置信区间
RD.CI2

# 估计ATE (风险差 RR)95%置信区间
gcomp.boot.RR <-function(formula = mod, data = data, indices) # 定义估算ATE的函数
{
  boot_sample <- data[indices, ]
  fit.boot <- glm(formula, data = boot_sample)
  Pred.Y1 <- predict(fit.boot,
                     newdata = newdata1,
                     type ="response")
  Pred.Y0 <- predict(fit.boot,
                     newdata = newdata0,
                     type ="response")
  Pred.TE <- mean(Pred.Y1, na.rm = TRUE) / mean(Pred.Y0, na.rm = TRUE)

}
set.seed(123)     # 设置随机种子以确保结果可重现

gcomp.RR <- boot(data=data,
                 statistic=gcomp.boot.RR,
                 R=250,
                 formula=mod)# 设置参数：R=250,生成250个bootstrap样本
plot(gcomp.RR)

RR.CI <- boot.ci(gcomp.RR, type="perc",conf=0.95)  # 基于百分位数的95%置信区间
RR.CI


################################################################
#二、使用R包stdReg2进行g公式标准化估计——二分类结局，计算RR值
################################################################

#1.定义变量集合
#除协变量外，需加入处理变量
covariates_std<- c("右心导管术","年龄", "收入","体重","DNR状态", "医疗保险",
                   "癌症", "心血管疾病史", "充血性心力衰竭", "痴呆", "精神病史",
                   "肝病", "上消化道出血","肿瘤", "免疫抑制","DASI指数",
                   "APACHE评分", "格拉斯哥昏迷评分","平均血压","白细胞计数", "体温",
                   "白蛋白", "红细胞比容","胆红素", "肌酐","钾",
                   "动脉血二氧化碳分压")

# 构建公式字符串
formula_std <- paste("死亡 ~", paste(covariates_std, collapse = " + "))
cat("使用的公式:\n", formula_std, "\n")

#2.使用R包stdReg2进行标准化估计
result <- standardize_glm(
  as.formula(formula_std),
  family = "binomial",
  data = data,
  values = list(右心导管术 = c(0, 1)),
  contrasts = c("difference", "ratio"),
  reference = 0,
  ci_level = 0.95,
  ci_type = "plain",
  transforms = NULL)

#3.结果解释

#查看标准化估计的汇总结果
print(result)

#获取整洁格式的结果表格
tidy_result <- tidy(result)
print(tidy_result)

################################################################
#三、使用R包RISCA 包估计ATE——二分类结局，计算OR值
################################################################
#使用R包RISCA 包中的 gc.logistic 函数估计二分类结局的 ATE
#1. 定义变量集合，需包含处理变量和至少一个协变量
covariates_std<- c("右心导管术","年龄", "收入","体重","DNR状态", "医疗保险",
                   "癌症", "心血管疾病史", "充血性心力衰竭", "痴呆", "精神病史",
                   "肝病", "上消化道出血","肿瘤", "免疫抑制","DASI指数",
                   "APACHE评分", "格拉斯哥昏迷评分","平均血压","白细胞计数", "体温",
                   "白蛋白", "红细胞比容","胆红素", "肌酐","钾",
                   "动脉血二氧化碳分压")

# 构建公式字符串
formula_std <- paste("死亡 ~", paste(covariates_std, collapse = " + "))
cat("使用的公式:\n", formula_std, "\n")

# 2.拟合多变量逻辑回归模型
model <- glm(
  as.formula(formula_std),  # 包含处理和协变量
  data = data,
  family = binomial(link = "logit")# 注意：必须使用 family = binomial(link = "logit")
)

#3. 使用gc.logistic 函数，通过 G计算估计 ATE

ate_result <- gc.logistic(
  glm.obj = model,
  data = data,
  group = "右心导管术",     #处理变量名称
  effect = "ATE",           #边际效应类型，"ATE"表示总体平均处理效应
  var.method = "bootstrap", #方差估计方法，"simulations"（模拟）或"bootstrap"（自助法）
  iterations = 250,         #模拟/自助法迭代次数（建议至少1000）
  n.cluster = 1             #并行计算使用的核心数
)


#4. 查看和分析结果

# 打印所有结果
print(ate_result)

#提取logOR
marginal_logOR <- ate_result$logOR  # 对数OR及其统计量

#转换为边际OR
marginal_OR <- data.frame(
  OR = exp(marginal_logOR[, "estimate"]),
  CI_lower = exp(marginal_logOR[, "ci.lower"]),
  CI_upper = exp(marginal_logOR[, "ci.upper"])
)
print(marginal_OR)



################################################################
#四、使用R包RISCA 包中的 gc.survival 函数估计生存结局的 ATE
################################################################
#1.构建Cox比例风险模型
# 定义模型变量集合
covariates_g<- c("右心导管术","年龄", "收入","体重","DNR状态", "医疗保险",
                 "癌症", "心血管疾病史", "充血性心力衰竭", "痴呆", "精神病史",
                 "肺病","肝病", "上消化道出血","肿瘤", "DASI指数",
                 "APACHE评分", "格拉斯哥昏迷评分","平均血压","白细胞计数", "体温",
                 "白蛋白", "红细胞比容","胆红素", "肌酐","钾",
                 "动脉血二氧化碳分压","血清PH值")

# 构建公式字符串
formula_g <- paste("Surv(生存时间, 死亡) ~", paste(covariates_g, collapse = " + "))
cat("使用的公式:\n", formula_g, "\n")

cox.cdt <- coxph(
  as.formula(formula_g),
  data = data,
  x = TRUE # 必须设置为TRUE！
)

# 2.使用G-computation计算边际效应

set.seed(123)  # 设置随机种子以保证结果可重现
# 运行gc.survival函数
gc.ate <- gc.survival(
  object = cox.cdt,                # 上面创建的校正Cox模型
  data = data,                     # 使用的数据集
  group = "右心导管术",            # 治疗/暴露变量名称
  times = "生存时间",        # 生存时间变量名称
  failures = "死亡",         # 事件状态变量名称（1=事件，0=删失）
  max.time = max(data$生存时间), # 计算RMST的最大时间（通常设为最大随访时间）
  effect = "ATE",                  # 边际效应类型：ATE=总体平均效应
                                   # 其他选项："ATT"=治疗组效应, "ATU"=未治疗组效应
  iterations = 250,                # Bootstrap重复次数（示例用250，实际建议≥1000）
  n.cluster = 1                    # 并行计算使用的核心数（1=不并行）
                                   # 如果数据量大，可增加核心数加速计算
)

# 3.查看和分析结果

# 3.1 查看完整输出
print(gc.ate)

# 3.2 提取边际HR的详细信息
marginal_logHR <- gc.ate$logHR  # 对数HR及其统计量

# 转换为普通边际HR
marginal_HR <- data.frame(
  HR = exp(marginal_logHR[, "estimate"]),
  CI_lower = exp(marginal_logHR[, "ci.lower"]),
  CI_upper = exp(marginal_logHR[, "ci.upper"]),
  p_value = marginal_logHR[, "p.value"]
)
print(marginal_HR)

# 3.3 查看边际生存概率表
# 包含每个时间点两组的生存概率和风险人数
surv_table <- gc.ate$table.surv
head(surv_table)

# 4.可视化结果
# 4.1 绘制边际生存曲线
# gc.survival输出的对象可以直接用plot()绘图
plot(gc.ate,
     main = "G-computation校正的生存曲线",
     ylab = "生存概率",
     xlab = "随访时间（天）",
     col = c("black", "red"),      # 曲线颜色（0组和1组）
     lty = c(1, 2),                # 线型
     lwd = 2)                      # 线宽

# 添加图例
legend("topright",
       legend = c("未治疗组", "治疗组"),
       col = c("black", "red"),
       lty = c(1, 2),
       lwd = 2,
       bty = "n")

# 4.2 也可以提取数据后使用ggplot2绘图（更美观）
library(ggplot2)

# 从table.surv中提取数据
surv_data <- gc.ate$table.surv
surv_data$variable <- factor(surv_data$variable,
                             levels = c(0, 1),
                             labels = c("未治疗", "治疗"))

ggplot(surv_data, aes(x = times, y = survival, color = variable)) +
  geom_line(size = 1.5) +
  labs(title = "G-computation校正的边际生存曲线",
       x = "随访时间（天）",
       y = "生存概率",
       color = "治疗组") +
  scale_color_manual(values = c("未治疗" = "black", "治疗" = "red")) +
  theme_minimal() +
  theme(legend.position = "top")


