# 创建新的输出目录
output_dir <- "E:\\dlnm\\region"
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

# 定义要分析的地区列表
regions <- c("东北地区", "华北地区", "华东地区", "华南地区", "华中地区", "西北地区", "西南地区")

# 定义事件类型列表
event_types <- c("flood", "severe_flood")

# 定义年龄组列表
age_groups <- c("case", "Pre_school", "School", "Adults", "Elders")

# 创建通用分析函数（修改后版本）
analyze_region_event <- function(region_name, event_type, age_group = "case", num_cores = 16) {
  # 筛选地区数据
  data_region <- subset(data, region == region_name)
  analysis_data <- data_region  # 使用新的变量名，避免覆盖全局data
  
  # 设置时区
  Sys.setlocale("LC_TIME", "English")
  analysis_data$date <- as.Date(analysis_data$date)
  analysis_data$dow <- as.factor(weekdays(analysis_data$date)) 
  analysis_data$month <- as.factor(months(analysis_data$date)) 
  analysis_data$year <- as.factor(format(analysis_data$date, format = "%Y")) 
  analysis_data$code <- as.factor(analysis_data$code)
  analysis_data$stratum <- as.factor(analysis_data$code:analysis_data$year:analysis_data$month:analysis_data$dow)
  
  # 计算气象变量
  temlag02 <- runMean(analysis_data$Tmean, 0:2, group = analysis_data$code)
  rhlag02 <- runMean(analysis_data$RH, 0:2, group = analysis_data$code)
  prelag02 <- runMean(analysis_data$pre, 0:2, group = analysis_data$code)
  
  # 设置并行
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  cat("正在处理地区:", region_name, "事件:", event_type, "年龄组:", age_group, "\n")
  cat("正在使用", num_cores, "个核心进行并行计算...\n")
  
  # 使用 foreach 并行计算
  results <- foreach(i = 0:28, .combine = "rbind",
                     .packages = c("MASS", "gnm", "dlnm", "splines", "tsModel"),
                     .export = c("analysis_data", "temlag02", "rhlag02", "prelag02")) %dopar% {
                       
                       tryCatch({
                         # 计算暴露变量
                         Explag <- runMean(analysis_data[[event_type]], i, group = analysis_data$code)
                         
                         # 根据事件类型设置模型公式
                         if(event_type %in% c("flood", "severe_flood")) {
                           # 洪涝和严重洪涝使用温度和湿度作为协变量
                           formula <- as.formula(paste0(age_group, " ~ Explag + ns(temlag02, 6) + as.factor(holiday) + ns(rhlag02, 3)"))
                         } 
                         
                         # 估计 theta
                         model_theta_est <- glm.nb(formula, data = analysis_data)
                         estimated_theta <- model_theta_est$theta
                         
                         # gnm 模型拟合
                         model1 <- gnm(formula,
                                       family = negative.binomial(estimated_theta), 
                                       eliminate = factor(stratum), 
                                       data = analysis_data)
                         
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
                         return(rep(NA, 7))
                       })
                     }
  
  # 停止集群
  stopCluster(cl)
  
  # 将结果转换为数据框
  results_df <- as.data.frame(results)
  names(results_df) <- c("coef", "se", "p_value", "aic", "rr", "rr_lower", "rr_upper")
  
  # 添加地区和事件信息
  results_df$region <- region_name
  results_df$event_type <- event_type
  results_df$age_group <- age_group
  results_df$lag <- 0:28
  
  # 创建输出文件名（使用output_dir）
  output_file <- file.path(output_dir, paste0(region_name, "_", age_group, "_", event_type, ".csv"))
  write.csv(results_df, file = output_file, row.names = FALSE)
  
  cat("成功输出文件:", output_file, "\n")
  
  return(results_df)
}

# 主循环 - 处理所有地区、所有事件类型和所有年龄组
all_results <- list()

for (region in regions) {
  region_results <- list()
  
  for (event in event_types) {
    for (age_group in age_groups) {
      # 运行分析
      result <- analyze_region_event(region, event, age_group)
      
      # 存储结果
      region_results[[paste(event, age_group, sep = "_")]] <- result
    }
  }
  
  # 存储当前地区的所有结果
  all_results[[region]] <- region_results
  
  # 保存当前地区的所有结果到一个RData文件（建议也使用output_dir）
  saveRDS(region_results, file = file.path(output_dir, paste0(region, "_all_results.rds")))
}

# 合并所有结果
final_results <- do.call(rbind, lapply(all_results, function(region) {
  do.call(rbind, region)
}))

# 保存最终结果到output_dir
write.csv(final_results, file = file.path(output_dir, "all_regions_events_agegroups_results.csv"), row.names = FALSE)

# 创建可视化函数
create_plots <- function(results_df, region_name, event_type, age_group) {
  # 过滤数据
  plot_data <- results_df %>%
    filter(region == region_name, event_type == event_type, age_group == age_group)
  
  # 创建基础图形
  p <- ggplot(plot_data, aes(x = lag, y = rr)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    geom_line(color = "blue") +
    geom_ribbon(aes(ymin = rr_lower, ymax = rr_upper), alpha = 0.2, fill = "blue") +
    geom_point(color = "blue", size = 1.5) +
    labs(title = paste("地区:", region_name, "事件:", event_type, "年龄组:", age_group),
         x = "Lag (days)", y = "RR") +
    theme_minimal()
  
  return(p)
}

# 为每个地区、事件和年龄组创建图形
for (region in regions) {
  for (event in event_types) {
    for (age_group in age_groups) {
      # 创建图形
      p <- create_plots(final_results, region, event, age_group)
      
      # 保存图形（使用output_dir）
      file_name <- file.path(output_dir, paste0(region, "_", event, "_", age_group, "_plot.png"))
      ggsave(file_name, p, width = 8, height = 6)
    }
  }
}

# 创建汇总图形（所有地区的事件效应）
summary_plots <- list()

for (event in event_types) {
  for (age_group in age_groups) {
    # 过滤数据
    plot_data <- final_results %>%
      filter(event_type == event, age_group == age_group)
    
    # 创建图形
    p <- ggplot(plot_data, aes(x = lag, y = rr, color = region)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      geom_line() +
      labs(title = paste("事件:", event, "年龄组:", age_group),
           x = "Lag (days)", y = "RR",
           color = "地区") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # 保存图形（使用output_dir）
    file_name <- file.path(output_dir, paste0("summary_", event, "_", age_group, "_plot.png"))
    ggsave(file_name, p, width = 10, height = 7)
    
    # 存储图形
    summary_plots[[paste(event, age_group, sep = "_")]] <- p
  }
}

# 保存当前全局环境所有对象（使用output_dir）
save.image(file.path(output_dir, "workspace_2025-10-12.RData"))

cat("所有分析已完成！结果已保存到目录:", output_dir, "\n")