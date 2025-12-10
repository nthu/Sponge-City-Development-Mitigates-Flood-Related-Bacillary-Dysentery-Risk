# 创建新的输出目录
output_dir <- "E:\\dlnm\\china"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("已创建输出目录:", output_dir, "\n")
} else {
  cat("输出目录已存在:", output_dir, "\n")
}

################时空分层的病例交叉设计#####
setwd("E:\\毕业论文\\1011\\dlnm")
data <- data.table::fread("1012fina_dlnm.csv")

library(dplyr)
library(splines)
library(gnm)
library(lubridate)
library(nlme)
library(mgcv)
library(tsModel)
library(ggplot2)
library(eoffice) 
library(doParallel)
Sys.setlocale("LC_TIME", "English")
data$date <-as.Date(data$date)
data$dow<-as.factor(weekdays(data$date)) 
data$month<-as.factor(months(data$date)) #
data$year<-as.factor(format(data$date,format="%Y")) 
data$code<-as.factor(data$code) ##
data$stratum<-as.factor (data$code:data$year:data$month:data$dow)#

temlag02<-runMean(data$Tmean,0:2,group = data$code)
rhlag02<-runMean(data$RH,0:2,group = data$code)
prelag02<-runMean(data$pre,0:2,group = data$code)

#########洪涝对儿童细菌性痢疾影响总#####
##############flood洪涝#####
num_cores <- 16
if (num_cores < 1) {
  num_cores <- 1
}
cl <- makeCluster(num_cores)

# 注册并行后端，使 foreach 知道要使用这个集群
registerDoParallel(cl)

# 打印出使用的核心数
cat("正在使用", num_cores, "个核心进行并行计算...\n")

# ***************************************************************
# 定义结果存储矩阵
# ***************************************************************
Flood <- matrix(NA, 28 + 1, 7,
                dimnames = list(paste("lag", 0:28),
                                c("coef", "se", "P-value", "AIC", "RR", "RR_lower", "RR_upper")))

# ***************************************************************
# 使用 foreach 循环代替 for 循环
# ***************************************************************
# .combine = "rbind" 告诉 foreach 将每个核心返回的结果（向量）按行绑定成一个矩阵
# .packages = c(...) 确保每个核心都加载了所需的包
# .export = c(...) 确保每个核心都能够访问主会话中的变量（如 data）
results <- foreach(i = 0:28, .combine = "rbind",
                   .packages = c("MASS", "gnm", "dlnm", "splines", "tsModel"),
                   .export = "data") %dopar% {
                     
                     # Floodlag 的计算部分
                     Floodlag <- runMean(data$flood, i, group = data$code)
                     
                     # theta 估计模型，这里我们移除了 stratum 以节省内存
                     model_theta_est <- glm.nb(case ~ Floodlag + ns(temlag02, 6) + as.factor(holiday) + 
                                                 ns(rhlag02, 3),
                                               data = data)
                     
                     # 提取 theta 值
                     estimated_theta <- model_theta_est$theta
                     
                     # gnm 模型拟合
                     model1 <- gnm(case ~ Floodlag + ns(temlag02, 6) + as.factor(holiday) +
                                     ns(rhlag02, 3),
                                   family = negative.binomial(estimated_theta), 
                                   eliminate = factor(stratum), 
                                   data = data)
                     
                     # 提取系数和标准误
                     summary_model <- summary(model1)
                     coef <- summary_model[["coefficients"]][1, 1]
                     se <- summary_model[["coefficients"]][1, 2]
                     p_value <- summary_model[["coefficients"]][1, 4]
                     aic <- summary_model[["aic"]]
                     
                     # 计算RR及其95%置信区间
                     rr <- exp(coef)
                     rr_lower <- exp(coef - 1.96 * se)
                     rr_upper <- exp(coef + 1.96 * se)
                     
                     # 返回一个包含所有结果的向量。foreach 会将这些向量按行合并。
                     return(c(coef, se, p_value, aic, rr, rr_lower, rr_upper))
                   }

# ***************************************************************
# 停止集群并清理
# ***************************************************************
stopCluster(cl)

# ***************************************************************
# 将并行计算的结果赋值给 Flood 矩阵并保存
# ***************************************************************
Flood <- results
write.csv(Flood, file = file.path(output_dir, "chinaflood.csv"), row.names = F)
# 打印最终结果
print(Flood)


##严重洪涝SF#####
num_cores <- 16
if (num_cores < 1) {
  num_cores <- 1
}
cl <- makeCluster(num_cores)

# 注册并行后端，使 foreach 知道要使用这个集群
registerDoParallel(cl)

# 打印出使用的核心数
cat("正在使用", num_cores, "个核心进行并行计算...\n")

# ***************************************************************
# 定义结果存储矩阵
# ***************************************************************
Flood <- matrix(NA, 28 + 1, 7,
                dimnames = list(paste("lag", 0:28),
                                c("coef", "se", "P-value", "AIC", "RR", "RR_lower", "RR_upper")))

# ***************************************************************
# 使用 foreach 循环代替 for 循环
# ***************************************************************
# .combine = "rbind" 告诉 foreach 将每个核心返回的结果（向量）按行绑定成一个矩阵
# .packages = c(...) 确保每个核心都加载了所需的包
# .export = c(...) 确保每个核心都能够访问主会话中的变量（如 data）
results <- foreach(i = 0:28, .combine = "rbind",
                   .packages = c("MASS", "gnm", "dlnm", "splines", "tsModel"),
                   .export = "data") %dopar% {
                     
                     # Floodlag 的计算部分
                     Floodlag <- runMean(data$severe_flood, i, group = data$code)
                     
                     # theta 估计模型，这里我们移除了 stratum 以节省内存
                     model_theta_est <- glm.nb(case ~ Floodlag + ns(temlag02, 6) + as.factor(holiday) + 
                                                 ns(rhlag02, 3),
                                               data = data)
                     
                     # 提取 theta 值
                     estimated_theta <- model_theta_est$theta
                     
                     # gnm 模型拟合
                     model1 <- gnm(case ~ Floodlag + ns(temlag02, 6) + as.factor(holiday) +
                                     ns(rhlag02, 3),
                                   family = negative.binomial(estimated_theta), 
                                   eliminate = factor(stratum), 
                                   data = data)
                     
                     # 提取系数和标准误
                     summary_model <- summary(model1)
                     coef <- summary_model[["coefficients"]][1, 1]
                     se <- summary_model[["coefficients"]][1, 2]
                     p_value <- summary_model[["coefficients"]][1, 4]
                     aic <- summary_model[["aic"]]
                     
                     # 计算RR及其95%置信区间
                     rr <- exp(coef)
                     rr_lower <- exp(coef - 1.96 * se)
                     rr_upper <- exp(coef + 1.96 * se)
                     
                     # 返回一个包含所有结果的向量。foreach 会将这些向量按行合并。
                     return(c(coef, se, p_value, aic, rr, rr_lower, rr_upper))
                   }

# ***************************************************************
# 停止集群并清理
# ***************************************************************
stopCluster(cl)

# ***************************************************************
# 将并行计算的结果赋值给 Flood 矩阵并保存
# ***************************************************************
Flood <- results
Flood <- as.data.frame(Flood)
names(Flood) <- c("coef", "se", "p_value", "aic", "rr", "rr_lower", "rr_upper")
write.csv(Flood, file = file.path(output_dir, "chinasf.csv"))
# 打印最终结果
print(Flood)



################Pre_school#####
##############flood洪涝#####
num_cores <- 16
if (num_cores < 1) {
  num_cores <- 1
}
cl <- makeCluster(num_cores)

# 注册并行后端，使 foreach 知道要使用这个集群
registerDoParallel(cl)

# 打印出使用的核心数
cat("正在使用", num_cores, "个核心进行并行计算...\n")

# ***************************************************************
# 定义结果存储矩阵
# ***************************************************************
Flood <- matrix(NA, 28 + 1, 7,
                dimnames = list(paste("lag", 0:28),
                                c("coef", "se", "P-value", "AIC", "RR", "RR_lower", "RR_upper")))

# ***************************************************************
# 使用 foreach 循环代替 for 循环
# ***************************************************************
# .combine = "rbind" 告诉 foreach 将每个核心返回的结果（向量）按行绑定成一个矩阵
# .packages = c(...) 确保每个核心都加载了所需的包
# .export = c(...) 确保每个核心都能够访问主会话中的变量（如 data）
results <- foreach(i = 0:28, .combine = "rbind",
                   .packages = c("MASS", "gnm", "dlnm", "splines", "tsModel"),
                   .export = "data") %dopar% {
                     
                     # Floodlag 的计算部分
                     Floodlag <- runMean(data$flood, i, group = data$code)
                     
                     model_theta_est <- glm.nb(Pre_school ~ Floodlag + ns(temlag02, 6) + as.factor(holiday) + 
                                                 ns(rhlag02, 3),
                                               data = data)
                     
                     # 提取 theta 值
                     estimated_theta <- model_theta_est$theta
                     
                     # gnm 模型拟合
                     model1 <- gnm(Pre_school ~ Floodlag + ns(temlag02, 6) + as.factor(holiday) +
                                     ns(rhlag02, 3),
                                   family = negative.binomial(estimated_theta), 
                                   eliminate = factor(stratum), 
                                   data = data)
                     
                     # 提取系数和标准误
                     summary_model <- summary(model1)
                     coef <- summary_model[["coefficients"]][1, 1]
                     se <- summary_model[["coefficients"]][1, 2]
                     p_value <- summary_model[["coefficients"]][1, 4]
                     aic <- summary_model[["aic"]]
                     
                     # 计算RR及其95%置信区间
                     rr <- exp(coef)
                     rr_lower <- exp(coef - 1.96 * se)
                     rr_upper <- exp(coef + 1.96 * se)
                     
                     # 返回一个包含所有结果的向量。foreach 会将这些向量按行合并。
                     return(c(coef, se, p_value, aic, rr, rr_lower, rr_upper))
                   }

# ***************************************************************
# 停止集群并清理
# ***************************************************************
stopCluster(cl)

# ***************************************************************
# 将并行计算的结果赋值给 Flood 矩阵并保存
# ***************************************************************
Flood <- results
Flood <- as.data.frame(Flood)
names(Flood) <- c("coef", "se", "p_value", "aic", "rr", "rr_lower", "rr_upper")
write.csv(Flood, file = file.path(output_dir, "Pre_schoolflood.csv"), row.names = F)
# 打印最终结果
print(Flood)


##严重洪涝SF#####
num_cores <- 16
if (num_cores < 1) {
  num_cores <- 1
}
cl <- makeCluster(num_cores)

# 注册并行后端，使 foreach 知道要使用这个集群
registerDoParallel(cl)

# 打印出使用的核心数
cat("正在使用", num_cores, "个核心进行并行计算...\n")

# ***************************************************************
# 定义结果存储矩阵
# ***************************************************************
Flood <- matrix(NA, 28 + 1, 7,
                dimnames = list(paste("lag", 0:28),
                                c("coef", "se", "P-value", "AIC", "RR", "RR_lower", "RR_upper")))

# ***************************************************************
# 使用 foreach 循环代替 for 循环
# ***************************************************************
# .combine = "rbind" 告诉 foreach 将每个核心返回的结果（向量）按行绑定成一个矩阵
# .packages = c(...) 确保每个核心都加载了所需的包
# .export = c(...) 确保每个核心都能够访问主会话中的变量（如 data）
results <- foreach(i = 0:28, .combine = "rbind",
                   .packages = c("MASS", "gnm", "dlnm", "splines", "tsModel"),
                   .export = "data") %dopar% {
                     
                     # Floodlag 的计算部分
                     Floodlag <- runMean(data$severe_flood, i, group = data$code)
                     
                     model_theta_est <- glm.nb(Pre_school ~ Floodlag + ns(temlag02, 6) + as.factor(holiday) + 
                                                 ns(rhlag02, 3),
                                               data = data)
                     
                     # 提取 theta 值
                     estimated_theta <- model_theta_est$theta
                     
                     # gnm 模型拟合
                     model1 <- gnm(Pre_school ~ Floodlag + ns(temlag02, 6) + as.factor(holiday) +
                                     ns(rhlag02, 3),
                                   family = negative.binomial(estimated_theta), 
                                   eliminate = factor(stratum), 
                                   data = data)
                     
                     # 提取系数和标准误
                     summary_model <- summary(model1)
                     coef <- summary_model[["coefficients"]][1, 1]
                     se <- summary_model[["coefficients"]][1, 2]
                     p_value <- summary_model[["coefficients"]][1, 4]
                     aic <- summary_model[["aic"]]
                     
                     # 计算RR及其95%置信区间
                     rr <- exp(coef)
                     rr_lower <- exp(coef - 1.96 * se)
                     rr_upper <- exp(coef + 1.96 * se)
                     
                     # 返回一个包含所有结果的向量。foreach 会将这些向量按行合并。
                     return(c(coef, se, p_value, aic, rr, rr_lower, rr_upper))
                   }

# ***************************************************************
# 停止集群并清理
# ***************************************************************
stopCluster(cl)

# ***************************************************************
# 将并行计算的结果赋值给 Flood 矩阵并保存
# ***************************************************************
Flood <- results
Flood <- as.data.frame(Flood)
names(Flood) <- c("coef", "se", "p_value", "aic", "rr", "rr_lower", "rr_upper")
write.csv(Flood, file = file.path(output_dir, "Pre_schoolsf.csv"), row.names = F)
# 打印最终结果
print(Flood)


nb_function <- function(num_cores = 16, outcome, variale, n_covariate, covariate, ns_df, output){
  num_cores <- num_cores
  if (num_cores < 1) {
    num_cores <- 1
  }
  cl <- makeCluster(num_cores)
  
  # 注册并行后端，使 foreach 知道要使用这个集群
  registerDoParallel(cl)
  
  # 打印出使用的核心数
  cat("正在使用", num_cores, "个核心进行并行计算...\n")
  
  # ***************************************************************
  # 定义结果存储矩阵
  # ***************************************************************
  Flood <- matrix(NA, 28 + 1, 7,
                  dimnames = list(paste("lag", 0:28),
                                  c("coef", "se", "P-value", "AIC", "RR", "RR_lower", "RR_upper")))
  
  # ***************************************************************
  # 使用 foreach 循环代替 for 循环
  # ***************************************************************
  # .combine = "rbind" 告诉 foreach 将每个核心返回的结果（向量）按行绑定成一个矩阵
  # .packages = c(...) 确保每个核心都加载了所需的包
  # .export = c(...) 确保每个核心都能够访问主会话中的变量（如 data）
  results <- foreach(i = 0:28, .combine = "rbind",
                     .packages = c("MASS", "gnm", "dlnm", "splines", "tsModel"),
                     .export = "data") %dopar% {
                       
                       # Floodlag 的计算部分
                       Floodlag <- runMean(data[[variale]], i, group = data$code)
                       if(n_covariate == 1){
                         formula <- as.formula(paste0(outcome, "~", "Floodlag + ns(", covariate[1], ",", ns_df[1], ")", "+ as.factor(holiday)"))
                       } else{
                         formula <- as.formula(paste0(outcome, "~", "Floodlag + ns(", covariate[1], ",", ns_df[1], ")", "+ as.factor(holiday) + ns(",  covariate[2], ", ", ns_df[2], ")"))
                       }
                       
                       model_theta_est <- glm.nb(formula,
                                                 data = data)
                       
                       # 提取 theta 值
                       estimated_theta <- model_theta_est$theta
                       
                       # gnm 模型拟合
                       model1 <- gnm(formula,
                                     family = negative.binomial(estimated_theta), 
                                     eliminate = factor(stratum), 
                                     data = data)
                       
                       # 提取系数和标准误
                       summary_model <- summary(model1)
                       coef <- summary_model[["coefficients"]][1, 1]
                       se <- summary_model[["coefficients"]][1, 2]
                       p_value <- summary_model[["coefficients"]][1, 4]
                       aic <- summary_model[["aic"]]
                       
                       # 计算RR及其95%置信区间
                       rr <- exp(coef)
                       rr_lower <- exp(coef - 1.96 * se)
                       rr_upper <- exp(coef + 1.96 * se)
                       
                       # 返回一个包含所有结果的向量。foreach 会将这些向量按行合并。
                       return(c(coef, se, p_value, aic, rr, rr_lower, rr_upper))
                     }
  
  # ***************************************************************
  # 停止集群并清理
  # ***************************************************************
  stopCluster(cl)
  
  # ***************************************************************
  # 将并行计算的结果赋值给 Flood 矩阵并保存
  # ***************************************************************
  Flood <- results
  Flood <- as.data.frame(Flood)
  names(Flood) <- c("coef", "se", "p_value", "aic", "rr", "rr_lower", "rr_upper")
  write.csv(Flood, file = file.path(output_dir, output), row.names = F)
  if(file.exists(file.path(output_dir, output))){
    cat("成功输出文件", file.path(output_dir, output))
  }
  # 打印最终结果
  print(Flood)
  return(Flood)
}


################Elders#####
##############flood洪涝#####
Flood <- nb_function(num_cores = 16, 
                     outcome = "Elders", 
                     variale = "flood", 
                     n_covariate = 2, 
                     covariate = c("temlag02", "rhlag02"), 
                     ns_df = c(6, 3), 
                     output = "Eldersflood.csv")


##严重洪涝SF#####
Flood <- nb_function(num_cores = 16, 
                     outcome = "Elders", 
                     variale = "severe_flood", 
                     n_covariate = 2, 
                     covariate = c("temlag02", "rhlag02"), 
                     ns_df = c(6, 3), 
                     output = "Elderssf.csv")
