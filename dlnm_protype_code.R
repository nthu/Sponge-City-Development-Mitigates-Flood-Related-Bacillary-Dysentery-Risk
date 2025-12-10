################################ initial #######################################
graphics.off()
rm(list = ls(all = TRUE))
gc(verbose=TRUE)
cat("\014")

library(rstudioapi)
library(this.path)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
print(this.path::this.dir())

# 创建新的输出目录
output_dir <- "E:\\dlnm\\pro_type"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("已创建输出目录:", output_dir, "\n")
} else {
  cat("输出目录已存在:", output_dir, "\n")
}


################时空分层的病例交叉设计#####
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
library(foreach)
library(MASS)
library(dlnm)

################时空分层的病例交叉设计#####
setwd("E:\\毕业论文\\1011\\dlnm")
data <- data.table::fread("1012fina_dlnm.csv")

# 定义要分析的省份类型列表
pro_types <- c("直辖市", "一般省份", "边疆省区", "经济发达省")

# 定义事件类型列表
event_types <- c("flood", "severe_flood")

# 定义年龄组列表
age_groups <- c("case", "Pre_school", "School", "Adults", "Elders")

# 创建通用分析函数
analyze_protype_event <- function(pro_type_name, event_type, age_group = "case", num_cores = 16) {
  
  # 参数验证
  if(!pro_type_name %in% pro_types) {
    stop("无效的省份类型: ", pro_type_name)
  }
  if(!event_type %in% event_types) {
    stop("无效的事件类型: ", event_type)
  }
  if(!age_group %in% age_groups) {
    stop("无效的年龄组: ", age_group)
  }
  
  # 筛选省份类型数据 
  data_protype <- subset(data, pro_type == pro_type_name)
  
  # 使用局部变量避免全局污染
  current_data <- data_protype
  
  # 设置时区
  Sys.setlocale("LC_TIME", "English")
  current_data$date <- as.Date(current_data$date)
  current_data$dow <- as.factor(weekdays(current_data$date)) 
  current_data$month <- as.factor(months(current_data$date)) 
  current_data$year <- as.factor(format(current_data$date, format = "%Y")) 
  current_data$code <- as.factor(current_data$code)
  current_data$stratum <- as.factor(current_data$code:current_data$year:current_data$month:current_data$dow)
  
  # 计算气象变量 - 使用current_data
  temlag02 <- runMean(current_data$Tmean, 0:2, group = current_data$code)
  rhlag02 <- runMean(current_data$RH, 0:2, group = current_data$code)
  prelag02 <- runMean(current_data$pre, 0:2, group = current_data$code)
  
  # 设置并行
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  cat("正在处理省份类型:", pro_type_name, "事件:", event_type, "年龄组:", age_group, "\n")
  cat("正在使用", num_cores, "个核心进行并行计算...\n")
  
  # 使用 foreach 并行计算
  results <- foreach(i = 0:28, .combine = "rbind",
                     .packages = c("MASS", "gnm", "dlnm", "splines", "tsModel"),
                     .export = c("current_data", "temlag02", "rhlag02", "prelag02", 
                                 "event_type", "age_group")) %dopar% {
                                   
                                   tryCatch({
                                     # 计算暴露变量 - 使用current_data
                                     Explag <- runMean(current_data[[event_type]], i, group = current_data$code)
                                     
                                     # 根据事件类型设置模型公式
                                     if(event_type %in% c("flood", "severe_flood")) {
                                       # 洪涝和严重洪涝使用温度和湿度作为协变量
                                       formula <- as.formula(paste0(age_group, " ~ Explag + ns(temlag02, 6) + as.factor(holiday) + ns(rhlag02, 3)"))
                                     } 
                                     # 估计 theta - 使用current_data
                                     model_theta_est <- glm.nb(formula, data = current_data)
                                     estimated_theta <- model_theta_est$theta
                                     
                                     # gnm 模型拟合 - 使用current_data
                                     model1 <- gnm(formula,
                                                   family = negative.binomial(estimated_theta), 
                                                   eliminate = factor(stratum), 
                                                   data = current_data)
                                     
                                     # 模型收敛检查
                                     if(!model1$converged) {
                                       warning("模型未收敛: ", pro_type_name, " - ", event_type, " - ", age_group, " - lag ", i)
                                     }
                                     
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
                                     
                                     # 返回结果
                                     return(c(coef, se, p_value, aic, rr, rr_lower, rr_upper))
                                     
                                   }, error = function(e) {
                                     # 错误处理
                                     cat("错误在 ", pro_type_name, "-", event_type, "-", age_group, "-lag", i, ":", e$message, "\n")
                                     return(rep(NA, 7))
                                   })
                                 }
  
  # 停止集群和清理
  stopCluster(cl)
  gc(verbose = FALSE)
  
  # 将结果转换为数据框
  results_df <- as.data.frame(results)
  names(results_df) <- c("coef", "se", "p_value", "aic", "rr", "rr_lower", "rr_upper")
  
  # 添加省份类型和事件信息
  results_df$pro_type <- pro_type_name
  results_df$event_type <- event_type
  results_df$age_group <- age_group
  results_df$lag <- 0:28
  
  # 创建输出文件名 - 修改为使用output_dir
  output_file <- file.path(output_dir, paste0(pro_type_name, "_", age_group, "_", event_type, ".csv"))
  write.csv(results_df, file = output_file, row.names = FALSE)
  
  cat("成功输出文件:", output_file, "\n")
  
  return(results_df)
}

# 主循环 - 处理所有省份类型、所有事件类型和所有年龄组
all_results <- list()

# 进度跟踪
total_iterations <- length(pro_types) * length(event_types) * length(age_groups)
current_iteration <- 0

for (pro_type in pro_types) {
  protype_results <- list()
  
  for (event in event_types) {
    for (age_group in age_groups) {
      current_iteration <- current_iteration + 1
      cat(sprintf("进度: %d/%d - 省份类型: %s, 事件: %s, 年龄组: %s\n", 
                  current_iteration, total_iterations, pro_type, event, age_group))
      
      # 运行分析
      result <- analyze_protype_event(pro_type, event, age_group)
      
      # 存储结果
      protype_results[[paste(event, age_group, sep = "_")]] <- result
    }
  }
  
  # 存储当前省份类型的所有结果
  all_results[[pro_type]] <- protype_results
  
  # 保存当前省份类型的所有结果到一个RData文件 - 修改为使用output_dir
  saveRDS(protype_results, file = file.path(output_dir, paste0(pro_type, "_all_results.rds")))
  
  # 清理内存
  rm(protype_results)
  gc(verbose = FALSE)
}

# 合并所有结果
final_results <- do.call(rbind, lapply(all_results, function(protype) {
  do.call(rbind, protype)
}))

# 保存最终结果 
write.csv(final_results, file = file.path(output_dir, "all_protypes_events_agegroups_results.csv"), row.names = FALSE)

# 创建可视化函数
create_plots <- function(results_df, pro_type_name, event_type, age_group) {
  # 过滤数据
  plot_data <- results_df %>%
    filter(pro_type == pro_type_name, event_type == event_type, age_group == age_group)
  
  # 检查是否有有效数据
  if(nrow(plot_data) == 0 || all(is.na(plot_data$rr))) {
    cat("无有效数据用于绘图:", pro_type_name, event_type, age_group, "\n")
    return(NULL)
  }
  
  # 创建基础图形
  p <- ggplot(plot_data, aes(x = lag, y = rr)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    geom_line(color = "blue") +
    geom_ribbon(aes(ymin = rr_lower, ymax = rr_upper), alpha = 0.2, fill = "blue") +
    geom_point(color = "blue", size = 1.5) +
    labs(title = paste("省份类型:", pro_type_name, "事件:", event_type, "年龄组:", age_group),
         x = "Lag (days)", y = "RR") +
    theme_minimal()
  
  return(p)
}

# 为每个省份类型、事件和年龄组创建图形
for (pro_type in pro_types) {
  for (event in event_types) {
    for (age_group in age_groups) {
      # 创建图形
      p <- create_plots(final_results, pro_type, event, age_group)
      
      if(!is.null(p)) {
        # 保存图形 - 修改为使用output_dir
        file_name <- file.path(output_dir, paste0(pro_type, "_", event, "_", age_group, "_plot.png"))
        ggsave(file_name, p, width = 8, height = 6)
        cat("已保存图形:", file_name, "\n")
      }
    }
  }
}

# 创建汇总图形（所有省份类型的事件效应）
summary_plots <- list()

for (event in event_types) {
  for (age_group in age_groups) {
    # 过滤数据
    plot_data <- final_results %>%
      filter(event_type == event, age_group == age_group)
    
    # 检查是否有有效数据
    if(nrow(plot_data) == 0 || all(is.na(plot_data$rr))) {
      cat("无有效数据用于汇总绘图:", event, age_group, "\n")
      next
    }
    
    # 创建图形
    p <- ggplot(plot_data, aes(x = lag, y = rr, color = pro_type)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      geom_line() +
      labs(title = paste("事件:", event, "年龄组:", age_group),
           x = "Lag (days)", y = "RR",
           color = "省份类型") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # 保存图形 - 修改为使用output_dir
    file_name <- file.path(output_dir, paste0("summary_", event, "_", age_group, "_plot.png"))
    ggsave(file_name, p, width = 10, height = 7)
    
    # 存储图形
    summary_plots[[paste(event, age_group, sep = "_")]] <- p
    
    cat("已保存汇总图形:", file_name, "\n")
  }
}

cat("所有分析已完成！结果已保存到CSV文件和PNG图形中。\n")

# 保存当前全局环境所有对象 
save.image(file.path(output_dir, "workspace-pro-type_2025-10-12.RData"))
cat("工作空间已保存到:", file.path(output_dir, "workspace-pro-type_2025-10-12.RData"), "\n")