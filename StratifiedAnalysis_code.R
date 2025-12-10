# 读取数据
library(tidyverse)
library(data.table)
library(MASS)
library(splines)
library(doParallel)
library(zoo)
setwd("F:\\flood\\did")
data <- data.table::fread("study_dataset.csv")

# ==================== 多期DID数据预处理 ====================
data_processed_corrected <- data %>%
  mutate(
    date = as.Date(date),
    
    # 处理组标识
    treated = ifelse(group == "sponge", 1, 0),
    
    # 明确区分两个队列
    cohort_2015 = ifelse(policy_year == "Exposed_2015", 1, 0),
    cohort_2016 = ifelse(policy_year == "Exposed_2016", 1, 0),
    control_2015 = ifelse(policy_year == "Non_Exposed_2015", 1, 0),
    control_2016 = ifelse(policy_year == "Non_Exposed_2016", 1, 0),
    
    # 为每个队列创建独立的相对年份
    relative_year_2015 = as.numeric(format(date, "%Y")) - 2015,
    relative_year_2016 = as.numeric(format(date, "%Y")) - 2016,
    
    # 队列特定的政策后变量
    post_2015 = ifelse((cohort_2015 == 1 | control_2015 == 1) & 
                         date >= as.Date("2015-04-02"), 1, 0),
    post_2016 = ifelse((cohort_2016 == 1 | control_2016 == 1) & 
                         date >= as.Date("2016-04-25"), 1, 0)
  ) %>%
  # 分别过滤各队列的数据
  filter(
    (cohort_2015 == 1 & relative_year_2015 >= -4 & relative_year_2015 <= 4) |
      (control_2015 == 1 & relative_year_2015 >= -4 & relative_year_2015 <= 4) |
      (cohort_2016 == 1 & relative_year_2016 >= -4 & relative_year_2016 <= 4) |
      (control_2016 == 1 & relative_year_2016 >= -4 & relative_year_2016 <= 4)
  ) %>%
  mutate(
    # 时间固定效应
    dow = as.factor(weekdays(date)),
    month = as.factor(month(date)),
    year = as.factor(year(date)),
    code = as.factor(code)
  )

# ==================== 滞后计算函数 ====================
calculate_lag_mean <- function(x, lag, group) {
  require(data.table)
  dt <- data.table(x = x, group = group)
  dt[, lag_mean := frollmean(x, n = lag + 1, align = "right", na.rm = TRUE), by = group]
  return(dt$lag_mean)
}

# 计算气象变量滞后
data_processed_corrected$temlag02 <- calculate_lag_mean(data_processed_corrected$Tmean, 2, data_processed_corrected$code)
data_processed_corrected$rhlag02 <- calculate_lag_mean(data_processed_corrected$RH, 2, data_processed_corrected$code)

# ==================== 分层分析函数 - 整合DID框架 ====================
# ==================== 分层分析：全部数据中海绵城市 vs 非海绵城市 ====================

stratified_city_analysis <- function(outcome = "case", 
                                     variable = "flood", 
                                     city_type = "sponge",  # "sponge" 或 "non_sponge"
                                     output_dir = "City_Stratified_Results",
                                     max_lag = 28) {
  
  cat("\n=== 开始城市分层DID分析 ===\n")
  cat("结局变量:", outcome, "\n")
  cat("暴露变量:", variable, "\n")
  cat("城市类型:", city_type, "\n")
  
  # 构建全部数据（all队列）的数据集
  analysis_data_all <- bind_rows(
    # 2015队列
    data_processed_corrected %>% 
      filter(cohort_2015 == 1 | control_2015 == 1) %>%
      mutate(
        post = post_2015,
        relative_year = relative_year_2015,
        cohort_treated = ifelse(cohort_2015 == 1, 1, 0),
        cohort_id = "2015"
      ),
    # 2016队列
    data_processed_corrected %>% 
      filter(cohort_2016 == 1 | control_2016 == 1) %>%
      mutate(
        post = post_2016,
        relative_year = relative_year_2016,
        cohort_treated = ifelse(cohort_2016 == 1, 1, 0),
        cohort_id = "2016"
      )
  ) %>%
    # 创建唯一的城市-队列标识
    mutate(
      code_cohort = as.factor(paste(code, cohort_id, sep = "_"))
    )
  
  # 根据城市类型筛选数据
  if (city_type == "sponge") {
    analysis_data <- analysis_data_all %>% filter(treated == 1)
  } else if (city_type == "non_sponge") {
    analysis_data <- analysis_data_all %>% filter(treated == 0)
  } else {
    stop("城市类型必须是 'sponge' 或 'non_sponge'")
  }
  
  # 创建相对年份因子，安全地设置参考水平
  analysis_data$rel_year_factor <- factor(analysis_data$relative_year)
  
  # 检查-1是否存在于相对年份中
  if("-1" %in% levels(analysis_data$rel_year_factor)) {
    analysis_data$rel_year_factor <- relevel(analysis_data$rel_year_factor, ref = "-1")
    cat("使用相对年份-1作为参考水平\n")
  } else {
    # 如果没有-1，使用最小的相对年份作为参考
    min_year <- min(analysis_data$relative_year, na.rm = TRUE)
    analysis_data$rel_year_factor <- relevel(analysis_data$rel_year_factor, ref = as.character(min_year))
    cat("警告: 相对年份-1不存在，使用", min_year, "作为参考水平\n")
  }
  
  # 初始化结果数据框（只关注DID效应）
  results_df <- data.frame(
    lag = 0:max_lag,
    city_type = city_type,
    outcome = outcome,
    variable = variable,
    coef_did = NA,
    se_did = NA,
    p_did = NA,
    rr_did = NA,
    rr_did_lower = NA,
    rr_did_upper = NA,
    aic = NA,
    n_obs = NA,
    n_cities = NA
  )
  
  # 逐个lag分析
  for(lag in 0:max_lag) {
    
    cat("\r处理 Lag", lag, "/", max_lag, "...")
    
    tryCatch({
      # 计算滞后暴露
      lag_exposure <- calculate_lag_mean(analysis_data[[variable]], lag, analysis_data$code)
      
      # 创建分析数据
      model_data <- analysis_data %>%
        mutate(
          exposure_lag = lag_exposure,
          exposure_binary = ifelse(exposure_lag > 0, 1, 0)
        ) %>%
        filter(!is.na(exposure_lag))
      
      # 使用事件研究法模型
      formula_es <- as.formula(paste0(
        outcome, " ~ exposure_binary * rel_year_factor + ",
        "ns(temlag02, df = 6) + ns(rhlag02, df = 3) + ",
        "as.factor(holiday) + factor(year) + factor(month) + factor(dow) + ",
        "factor(code_cohort)"  # 城市-队列固定效应
      ))
      
      # 拟合模型
      model <- glm.nb(formula_es, data = model_data)
      
      # 提取DID效应（政策后第一年的暴露效应）
      did_term <- "exposure_binary:rel_year_factor0"  # 政策后第一年
      
      if(did_term %in% rownames(summary(model)$coefficients)) {
        coef_summary <- summary(model)$coefficients
        results_df$coef_did[lag + 1] <- coef_summary[did_term, "Estimate"]
        results_df$se_did[lag + 1] <- coef_summary[did_term, "Std. Error"]
        results_df$p_did[lag + 1] <- coef_summary[did_term, "Pr(>|z|)"]
        results_df$rr_did[lag + 1] <- exp(results_df$coef_did[lag + 1])
        results_df$rr_did_lower[lag + 1] <- exp(results_df$coef_did[lag + 1] - 1.96 * results_df$se_did[lag + 1])
        results_df$rr_did_upper[lag + 1] <- exp(results_df$coef_did[lag + 1] + 1.96 * results_df$se_did[lag + 1])
      } else {
        # 备选：如果没有交互项，使用政策后的暴露效应
        alt_term <- "exposure_binary"
        if(alt_term %in% rownames(summary(model)$coefficients)) {
          coef_summary <- summary(model)$coefficients
          results_df$coef_did[lag + 1] <- coef_summary[alt_term, "Estimate"]
          results_df$se_did[lag + 1] <- coef_summary[alt_term, "Std. Error"]
          results_df$p_did[lag + 1] <- coef_summary[alt_term, "Pr(>|z|)"]
          results_df$rr_did[lag + 1] <- exp(results_df$coef_did[lag + 1])
          results_df$rr_did_lower[lag + 1] <- exp(results_df$coef_did[lag + 1] - 1.96 * results_df$se_did[lag + 1])
          results_df$rr_did_upper[lag + 1] <- exp(results_df$coef_did[lag + 1] + 1.96 * results_df$se_did[lag + 1])
        }
      }
      
      # 模型拟合统计
      results_df$aic[lag + 1] <- AIC(model)
      results_df$n_obs[lag + 1] <- nrow(model$model)
      results_df$n_cities[lag + 1] <- length(unique(model_data$code))
      
    }, error = function(e) {
      cat("\nLag", lag, "出错:", as.character(e), "\n")
    })
  }
  
  cat("\n城市分层分析完成！\n")
  
  # 创建输出目录
  if(!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # 保存结果
  output_file <- paste0(output_dir, "/City_Stratified_", city_type, "_", outcome, "_", variable, ".csv")
  write.csv(results_df, output_file, row.names = FALSE)
  cat("结果已保存至:", output_file, "\n")
  
  return(results_df)
}

# ==================== 批量执行城市分层分析 ====================

# 定义分析参数
outcomes <- c("case")
variables <- c("flood", "severe_flood")
city_types <- c("sponge", "non_sponge")

# 执行批量城市分层分析
city_stratified_results <- list()

for(city_type in city_types) {
  for(outcome in outcomes) {
    for(variable in variables) {
      
      cat("\n", paste0(rep("=", 60), collapse = ""), "\n")
      cat("开始城市分层分析: city_type =", city_type, ", outcome =", outcome, ", variable =", variable, "\n")
      
      result <- stratified_city_analysis(
        outcome = outcome,
        variable = variable,
        city_type = city_type,
        output_dir = "City_Stratified_Results",
        max_lag = 28
      )
      
      # 存储结果
      key <- paste(city_type, outcome, variable, sep = "_")
      city_stratified_results[[key]] <- result
      
      # 打印简要结果
      valid_results <- sum(!is.na(result$coef_did))
      cat("分析完成. 有效DID结果数量:", valid_results, "/", nrow(result), "\n")
      
      if(valid_results > 0) {
        # 显示前5个lag的DID效应
        cat("前5个lag的DID效应:\n")
        print(result[1:min(5, nrow(result)), c("lag", "rr_did", "p_did")])
      }
    }
  }
}

#####不同年龄#####

# 定义分析参数
outcomes <- c("Pre_school", "School", "Adults", "Elders")
variables <- c("flood", "severe_flood")
city_types <- c("sponge", "non_sponge")

# 执行批量城市分层分析
city_stratified_results <- list()

for(city_type in city_types) {
  for(outcome in outcomes) {
    for(variable in variables) {
      
      cat("\n", paste0(rep("=", 60), collapse = ""), "\n")
      cat("开始城市分层分析: city_type =", city_type, ", outcome =", outcome, ", variable =", variable, "\n")
      
      result <- stratified_city_analysis(
        outcome = outcome,
        variable = variable,
        city_type = city_type,
        output_dir = "City_Stratified_Results",
        max_lag = 28
      )
      
      # 存储结果
      key <- paste(city_type, outcome, variable, sep = "_")
      city_stratified_results[[key]] <- result
      
      # 打印简要结果
      valid_results <- sum(!is.na(result$coef_did))
      cat("分析完成. 有效DID结果数量:", valid_results, "/", nrow(result), "\n")
      
      if(valid_results > 0) {
        # 显示前5个lag的DID效应
        cat("前5个lag的DID效应:\n")
        print(result[1:min(5, nrow(result)), c("lag", "rr_did", "p_did")])
      }
    }
  }
}




# ==================== 可视化比较函数 ====================
setwd("F:\\Flood\\did\\City_Stratified_Results")
data <- data.table::fread("All_City_Stratified_sponge_case_flood.csv")
library(ggplot2)
library(dplyr)
library(tidyr)

# Nature风格颜色
nature_colors <- c("#1F77B4", "#FF7F0E", "#2CA02C")

# 专业主题设置
theme_nature <- function() {
  theme_minimal(base_family = "Arial") +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.title = element_text(face = "bold"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2)
    )
}

# 数据准备函数
prepare_city_data <- function(data) {
  data %>%
    filter(!is.na(rr_did)) %>%
    mutate(
      city_type = factor(city_type, levels = c("sponge", "non_sponge")),
      flood_type = ifelse(variable == "flood", "General Flood", "Severe Flood"),
      # 基于置信区间判断显著性
      significance = case_when(
        rr_did_lower > 1 | rr_did_upper < 1 ~ "*",
        TRUE ~ "NS"
      )
    )
}

# 美化后的分面显示函数
plot_rr_by_lag_nature <- function(data) {
  plot_data <- prepare_city_data(data)
  
  p <- ggplot(plot_data, aes(x = lag, y = rr_did, color = city_type, linetype = city_type)) +
    # 添加参考线
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
    
    # 添加置信区间
    geom_ribbon(aes(ymin = rr_did_lower, ymax = rr_did_upper, fill = city_type),
                alpha = 0.15, color = NA) +
    
    # 添加线条和点
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = significance), size = 1.5, alpha = 0.8) +
    
    # 分面显示两个变量
    facet_wrap(~ flood_type, ncol = 2, scales = "fixed") +
    
    # 设置颜色和样式
    scale_color_manual(
      values = c("sponge" = "#1F77B4", "non_sponge" = "#FF7F0E"),
      labels = c("sponge" = "Sponge City", "non_sponge" = "Non-Sponge City")
    ) +
    scale_fill_manual(
      values = c("sponge" = "#1F77B4", "non_sponge" = "#FF7F0E"),
      labels = c("sponge" = "Sponge City", "non_sponge" = "Non-Sponge City")
    ) +
    scale_linetype_manual(
      values = c("sponge" = "solid", "non_sponge" = "dashed"),
      labels = c("sponge" = "Sponge City", "non_sponge" = "Non-Sponge City")
    ) +
    scale_shape_manual(values = c("*" = 16, "NS" = 1)) +
    scale_x_continuous(breaks = seq(0, 28, 7)) +
    
    # 设置坐标轴和标题
    labs(
      x = "Lag (Days)",
      y = "Risk Ratio",
      color = "City Type",
      fill = "City Type", 
      linetype = "City Type",
      shape = "Significance"
    ) +
    
    # 应用Nature主题
    theme_nature()
  
  return(p)
}

# 生成美化后的图形
plot_nature <- plot_rr_by_lag_nature(data)
print(plot_nature)
plot_nature 
topptx(filename = "F:\\Flood\\did\\DID_Results\\stri_interaction.pptx",width =8, height = 4)


##########age####
# 设置工作目录
setwd("H:/毕业论文/1011/Flood/did/City_Stratified_Results/age")

# 获取所有CSV文件
csv_files <- list.files(pattern = "\\.csv$", full.names = TRUE)

# 检查是否找到文件
if (length(csv_files) == 0) {
  stop("在指定路径下未找到CSV文件")
}

# 读取并合并所有CSV文件
merged_data <- do.call(rbind, lapply(csv_files, function(file) {
  df <- read.csv(file)
  df$source_file <- basename(file)  # 添加源文件名
  return(df)
}))

# 保存合并后的文件
write.csv(merged_data, "result.csv", row.names = FALSE)

data <- merged_data
library(ggplot2)
library(dplyr)
library(tidyr)
library(eoffice)  # 用于导出PPT

# Nature风格颜色
nature_colors <- c("#1F77B4", "#FF7F0E", "#2CA02C")

# 专业主题设置
theme_nature <- function() {
  theme_minimal(base_family = "Arial") +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.title = element_text(face = "bold"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2)
    )
}

# 数据准备函数 - 针对年龄组数据
prepare_age_data <- function(data) {
  data %>%
    filter(!is.na(rr_did)) %>%
    mutate(
      city_type = factor(city_type, levels = c("sponge", "non_sponge")),
      flood_type = ifelse(variable == "flood", "General Flood", "Severe Flood"),
      # 基于置信区间判断显著性
      significance = case_when(
        rr_did_lower > 1 | rr_did_upper < 1 ~ "*",
        TRUE ~ "NS"
      ),
      # 美化年龄组标签
      age_group = case_when(
        outcome == "Pre_school" ~ "Pre-school children",
        outcome == "School" ~ "School-aged children", 
        outcome == "Adults" ~ "Adults",
        outcome == "Elders" ~ "Elderly",
        TRUE ~ outcome
      )
    ) %>%
    # 设置年龄组顺序
    mutate(age_group = factor(age_group, 
                              levels = c("Pre-school children", "School-aged children", "Adults", "Elderly")))
}

# 分年龄组可视化函数 - 网格布局
plot_rr_by_age_nature <- function(data) {
  plot_data <- prepare_age_data(data)
  
  p <- ggplot(plot_data, aes(x = lag, y = rr_did, color = city_type, linetype = city_type)) +
    # 添加参考线
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
    
    # 添加置信区间
    geom_ribbon(aes(ymin = rr_did_lower, ymax = rr_did_upper, fill = city_type),
                alpha = 0.15, color = NA) +
    
    # 添加线条和点
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = significance), size = 1.5, alpha = 0.8) +
    
    # 分面显示：年龄组为行，洪水类型为列
    facet_grid(age_group ~ flood_type, scales = "free_y") +
    
    # 设置颜色和样式
    scale_color_manual(
      values = c("sponge" = "#1F77B4", "non_sponge" = "#FF7F0E"),
      labels = c("sponge" = "Sponge City", "non_sponge" = "Non-Sponge City")
    ) +
    scale_fill_manual(
      values = c("sponge" = "#1F77B4", "non_sponge" = "#FF7F0E"),
      labels = c("sponge" = "Sponge City", "non_sponge" = "Non-Sponge City")
    ) +
    scale_linetype_manual(
      values = c("sponge" = "solid", "non_sponge" = "dashed"),
      labels = c("sponge" = "Sponge City", "non_sponge" = "Non-Sponge City")
    ) +
    scale_shape_manual(values = c("*" = 16, "NS" = 1)) +
    scale_x_continuous(breaks = seq(0, 28, 7)) +
    
    # 设置坐标轴和标题
    labs(
      x = "Lag (Days)",
      y = "Risk Ratio",
      color = "City Type",
      fill = "City Type", 
      linetype = "City Type",
      shape = "Significance",
      #title = "Risk Ratios by Age Group and Flood Type"
    ) +
    
    # 应用Nature主题
    theme_nature() +
    theme(
      strip.text = element_text(size = 10),
      axis.text = element_text(size = 9)
    )
  
  return(p)
}
# 1. 生成所有年龄组的综合图形（网格布局）
plot_all_ages <- plot_rr_by_age_nature(data)
print(plot_all_ages)

# 导出综合图形
topptx(plot_all_ages, 
       filename = "H:\\毕业论文\\1011\\Flood\\did\\City_Stratified_Results\\age\\All_Age_Groups_Grid.pptx",
       width = 8, height = 8)
# 分年龄组可视化函数 - 并排布局（替代方案）
plot_rr_by_age_nature_wrap <- function(data) {
  plot_data <- prepare_age_data(data)
  
  p <- ggplot(plot_data, aes(x = lag, y = rr_did, color = city_type, linetype = city_type)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_ribbon(aes(ymin = rr_did_lower, ymax = rr_did_upper, fill = city_type),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = significance), size = 1.5, alpha = 0.8) +
    
    # 使用facet_wrap，每个年龄组一个图形
    facet_wrap(~ age_group + flood_type, ncol = 2, scales = "free_y") +
    
    scale_color_manual(
      values = c("sponge" = "#1F77B4", "non_sponge" = "#FF7F0E"),
      labels = c("sponge" = "Sponge City", "non_sponge" = "Non-Sponge City")
    ) +
    scale_fill_manual(
      values = c("sponge" = "#1F77B4", "non_sponge" = "#FF7F0E"),
      labels = c("sponge" = "Sponge City", "non_sponge" = "Non-Sponge City")
    ) +
    scale_linetype_manual(
      values = c("sponge" = "solid", "non_sponge" = "dashed"),
      labels = c("sponge" = "Sponge City", "non_sponge" = "Non-Sponge City")
    ) +
    scale_shape_manual(values = c("*" = 16, "NS" = 1)) +
    scale_x_continuous(breaks = seq(0, 28, 7)) +
    
    labs(
      x = "Lag (Days)",
      y = "Risk Ratio",
      color = "City Type",
      fill = "City Type", 
      linetype = "City Type",
      shape = "Significance"
    ) +
    
    theme_nature()
  
  return(p)
}

# 单一年龄组可视化函数
plot_single_age_nature <- function(data, age_group_name) {
  plot_data <- prepare_age_data(data) %>%
    filter(age_group == age_group_name)
  
  p <- ggplot(plot_data, aes(x = lag, y = rr_did, color = city_type, linetype = city_type)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_ribbon(aes(ymin = rr_did_lower, ymax = rr_did_upper, fill = city_type),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = significance), size = 1.5, alpha = 0.8) +
    
    facet_wrap(~ flood_type, ncol = 2, scales = "fixed") +
    
    scale_color_manual(
      values = c("sponge" = "#1F77B4", "non_sponge" = "#FF7F0E"),
      labels = c("sponge" = "Sponge City", "non_sponge" = "Non-Sponge City")
    ) +
    scale_fill_manual(
      values = c("sponge" = "#1F77B4", "non_sponge" = "#FF7F0E"),
      labels = c("sponge" = "Sponge City", "non_sponge" = "Non-Sponge City")
    ) +
    scale_linetype_manual(
      values = c("sponge" = "solid", "non_sponge" = "dashed"),
      labels = c("sponge" = "Sponge City", "non_sponge" = "Non-Sponge City")
    ) +
    scale_shape_manual(values = c("*" = 16, "NS" = 1)) +
    scale_x_continuous(breaks = seq(0, 28, 7)) +
    
    labs(
      x = "Lag (Days)",
      y = "Risk Ratio",
      title = paste("Risk Ratios for", age_group_name),
      color = "City Type",
      fill = "City Type", 
      linetype = "City Type",
      shape = "Significance"
    ) +
    
    theme_nature()
  
  return(p)
}

# 使用示例：
# 假设您的数据已经加载并命名为data，且包含age_group列




# 2. 生成所有年龄组的综合图形（并排布局）
plot_all_ages_wrap <- plot_rr_by_age_nature_wrap(data)
print(plot_all_ages_wrap)

# 3. 分别生成每个年龄组的图形
age_groups <- c("Pre-school", "School Age", "Adults", "Elderly")

for(age in age_groups) {
  plot_single <- plot_single_age_nature(data, age)
  print(plot_single)
  
  # 导出单个年龄组图形
  filename <- paste0("F:\\毕业论文\\1011\\Flood\\did\\City_Stratified_Results\\age\\", 
                     gsub(" ", "_", age), "_Results.pptx")
  topptx(plot_single, filename = filename, width = 8, height = 4)
}

# 检查数据
print("数据概览:")
print(head(data))
print("年龄组分布:")
if("age_group" %in% names(data)) {
  print(table(data$age_group))
} else {
  print("数据中没有age_group列，请检查数据格式")
}