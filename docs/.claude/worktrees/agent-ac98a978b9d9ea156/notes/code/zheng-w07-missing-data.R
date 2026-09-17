

#设置工作空间
setwd("D:/数据")


#读取数据
data<-read.csv("data_missing.csv")
colnames(data)


####缺失数据分析####

#检查缺失值情况
install.packages("dplyr")
library(dplyr)
datana<-data

# 缺失模式
install.packages("mice")
library(mice)
datana1<-datana[,c("height" , "weight" , "PA","PA_CAT", "CESD","CESD_CAT")]
md.pattern(datana)


#各变量缺失的比例
var <- names(datana1)
tab <- list()
total <- nrow(datana1)
for (i in 1:length(var)){
  tab[i]=data.frame(Variables=var[i],
                    Total=total,
                    Freq=total-sum(is.na(datana1[,var[i]])),
                    Missing=sum(is.na(datana1[,var[i]])),
                    miss_p=sprintf("%0.2f",sum(is.na(datana1[,var[i]]))/total*100)) %>% list()}
#缺失值合并
na<-do.call(rbind, tab)
View(na)


#缺失比例可视化
install.packages("VIM")
library(VIM)
mice_plot <- aggr(datana1, col=c('navyblue','yellow'),
                  numbers=TRUE, sortVars=TRUE,
                  labels=names(datana1), cex.axis=.7,
                  gap=5, ylab=c("Missing data","Pattern"),)


#填补
#填补前处理
str(datana)
vars<-c("CESD_CAT","sex","hypertension","diabetes","heart_disease","residence","smoking","education_level","cardio_metabolic", "alcohol")
datana[vars] <- lapply(datana[vars], factor)


#多重插补
# m为模型数量，通常为5，maxit为模型迭代次数，通常为50，此处为节省时间进行缩减
install.packages("mice")
library(mice)
multi_imputed<-mice(datana, m=5, seed = 123)
meth<-multi_imputed$method 

#根据不同类型数据选择不同插补方法
meth["height"]="pmm"           #定量数据
meth[c("CESD_CAT","sex","hypertension","diabetes","heart_disease","residence","smoking","cardio_metabolic", "alcohol")]="logreg"    #二分类
meth["education_level"]="polyreg"                 #无序多分类

multi_imputed1<-mice(datana, m=5,method=meth, seed = 123)
multi_imputed1$method 

imp1<-complete(multi_imputed,1)
imp2<-complete(multi_imputed,2)

names(datana)

#填补后回归分析
fit<-with(multi_imputed1,{
  glm(lungfun~PA+CESD_CAT+sex+age+diabetes,family = binomial(link ="logit"))
})

pooled<-pool(fit)
summary(pooled)  


#倾向性匹配
install.packages("MatchThem")
library(MatchThem)

#读取数据查看数据
datana<- read.csv("data_missing.csv")
summary(datana)

#查看数据类型
str(datana)

#使用 mice 包对肺炎数据集进行多重填补
library(mice)
multi_imputed<-mice(datana, m=5, seed = 123) #生成 5 个填补数据集

names(datana)

# 使用matchthem( )来匹配多个填补数据集
matched.datasets <- matchthem(
  lungfun ~ CESD_CAT+sex+age+hypertension+diabetes+heart_disease+smoking,
  datasets = multi_imputed,
  approach = 'within', #内部匹配
  method = 'nearest', #邻近匹配
  caliper = 0.2, #倾向得分标准差 5%
  ratio = 2) #2:1 的匹配比


# 评价匹配效果
library(cobalt)
bal.tab(matched.datasets,stats = c("m","ks"),imp.fun= "max") #m: 绝对标准化均值差，ks: Kolmogorov-Smirnov (KS) 统计量

# 匹配效果可视化
love.plot(bal.tab(matched.datasets, m.threshold=0.1),                            
          stat = "mean.diffs", grid=TRUE, stars="raw", abs = F)
names(datana)

#匹配后分析
# 对每一份匹配后的填补数据集进行二元logistic回归分析
install.packages("survey")
library(survey)
matched.models <- with(matched.datasets,
                       svyglm(lungfun ~ PA, family = quasibinomial()),
                       cluster = TRUE)

# 利用MICE包将每份填补数据集估计的效应值进行综合
match.results <- pool(matched.models)
summary(match.results)


#倾向性得加加权
# 使用 weightthem( ) 来加权多个填补数据集
weighted.datasets <- weightthem(lungfun ~ CESD_CAT+sex+age+hypertension+diabetes+heart_disease+smoking,
                                datasets = multi_imputed,
                                approach = 'across', #交叉加权和逻辑回归倾向得分加权
                                method = 'ps',  
                                estimand = "ATE") #估计值加权样本的平均治疗效果（ATE）


# 评价加权效果
library(cobalt)
bal.tab(weighted.datasets, stats = c('m', 'ks'), imp.fun = 'max') #m: 绝对标准化均值差，ks: Kolmogorov-Smirnov (KS) 统计量

# 加权效果可视化
love.plot(bal.tab(weighted.datasets, m.threshold=0.1),                         
          stat = "mean.diffs", grid=TRUE, stars="raw",abs = F)

# 对每一份加权后的填补数据集进行二元logistic回归分析
library(survey)
weighted.models <- with(weighted.datasets,
                        svyglm( lungfun ~ PA, family = quasibinomial()))

# 将每份填补数据集估计的效应值进行综合
weighted.results <- pool(weighted.models)
summary(weighted.results, conf.int = TRUE)



####多重插补的双重稳健法TMLE
# 加载必要的包
#install.packages("tmle")
library(mice)
library(tmle)

#设置工作空间
setwd("D:/数据")


#读取数据
data<-read.csv("data_missing.csv")
datana<-data

#多重插补
set.seed(123)
multi_imputed<-mice(datana, m=5, seed = 123)

# 提取分析数据
#协变量，要是数据框的形式
W <- c("height" , "weight" , "CESD")  

#暴露
A<-"PA_CAT"

#结局
Y<-"lungfun"

#提取m
m<-multi_imputed$m

library(dplyr)

for (i in 1:m) {
  dat_i <- complete(multi_imputed, action = i)
  # 取出 Y, A, W
  Y_dat <- dat_i[,Y]
  A_dat <- dat_i[,A]
  W_dat <- dat_i %>% select(all_of(W)) %>% as.data.frame()

  #构建TMLE模型
  tmle_fit <- tmle(Y = Y_dat , A = A_dat, W = W_dat,
                   Q.SL.library = c("SL.glm", "SL.mean"),  # 简单SL库
                   g.SL.library = c("SL.glm", "SL.mean"),
                   family = "binomial")
  #ATE结果合并
  tab <- data.frame(psi = tmle_fit$estimates$ATE$psi,
                    var.psi = tmle_fit$estimates$ATE$var.psi)
  if (i==1){
    table <- tab
  }
  else {
    table <- rbind(table,tab)
  }
}


#rubin合并
rubin_pool_scalar <- function(df, est_col = "psi", var_col = "var.psi",
                              conf_level = 0.95, tiny_B = 1e-12) {
  stopifnot(is.data.frame(df), est_col %in% names(df), var_col %in% names(df))
  
  Q <- as.numeric(df[[est_col]])      # Q_i
  U <- as.numeric(df[[var_col]])      # U_i
  ok <- is.finite(Q) & is.finite(U)
  Q <- Q[ok]; U <- U[ok]
  m <- length(Q)
  if (m < 2) stop("至少需要 2 个插补数据集的结果才能做 Rubin 合并。")
  
  Qbar <- mean(Q)          # 合并点估计
  Ubar <- mean(U)          # 插补内方差均值
  B    <- stats::var(Q)    # 插补间方差（样本方差，分母 m-1）
  
  # 总方差
  if (is.na(B) || B < tiny_B) {
    Tvar <- Ubar
    df_br <- Inf
  } else {
    Tvar <- Ubar + (1 + 1/m) * B
    df_br <- (m - 1) * (1 + Ubar / ((1 + 1/m) * B))^2   # Barnard–Rubin df
  }
  
  se <- sqrt(Tvar)
  alpha <- 1 - conf_level
  crit <- if (is.finite(df_br)) stats::qt(1 - alpha/2, df = df_br) else stats::qnorm(1 - alpha/2)
  
  ci_low  <- Qbar - crit * se
  ci_high <- Qbar + crit * se
  
  t_stat <- Qbar / se
  p_val <- if (is.finite(df_br)) {
    2 * stats::pt(abs(t_stat), df = df_br, lower.tail = FALSE)
  } else {
    2 * stats::pnorm(abs(t_stat), lower.tail = FALSE)
  }
  
  # 评价插补效果：RIV / FMI
  RIV <- if (Tvar > 0) ((1 + 1/m) * B) / Ubar else NA_real_
  FMI <- if (Tvar > 0) ((1 + 1/m) * B) / Tvar else NA_real_
  
  out <- data.frame(
    m = m,
    est = Qbar,
    se = se,
    df = df_br,
    lcl = ci_low,
    ucl = ci_high,
    t = t_stat,
    p_value = p_val,
    Ubar = Ubar,
    B = B,
    Tvar = Tvar,
    RIV = RIV,
    FMI = FMI
  )
  return(out)
}

#rubin合并的结果
pooled <- rubin_pool_scalar(table, est_col = "psi", var_col = "var.psi")
pooled
