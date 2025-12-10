############## 事件研究法：海绵城市政策前后的变化（前后4年）#####
# 创建相对时间变量（允许政策前后4年）

setwd("F:\\flood\\did")
data <- data.table::fread("study_dataset.csv")

# 使用与多期DID相同的预处理
data_processed <- data %>%
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
    # 2015队列：处理组+控制组，相对年份在±4年内
    (cohort_2015 == 1 & relative_year_2015 >= -4 & relative_year_2015 <= 4) |
      (control_2015 == 1 & relative_year_2015 >= -4 & relative_year_2015 <= 4) |
      # 2016队列：处理组+控制组，相对年份在±4年内
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

# 使用相同的滞后计算函数
calculate_lag_mean <- function(x, lag, group) {
  require(data.table)
  dt <- data.table(x = x, group = group)
  dt[, lag_mean := frollmean(x, n = lag + 1, align = "right", na.rm = TRUE), by = group]
  return(dt$lag_mean)
}

# 计算气象变量滞后
data_processed$temlag02 <- calculate_lag_mean(data_processed$Tmean, 2, data_processed$code)
data_processed$rhlag02 <- calculate_lag_mean(data_processed$RH, 2, data_processed$code)

# 修正后的事件研究法函数（基于多期DID框架）
event_study_analysis <- function(cohort = "2015", outcome = "case", variable = "flood", 
                                 lag_value = 28, output_file) {
  
  cat("\n=== 开始事件研究法分析 ===\n")
  cat("队列:", cohort, "\n")
  cat("结局变量:", outcome, "\n")
  cat("暴露变量:", variable, "\n")
  cat("滞后天数:", lag_value, "\n")
  
  # 选择对应队列的数据
  if (cohort == "2015") {
    analysis_data <- data_processed %>% 
      filter(cohort_2015 == 1 | control_2015 == 1) %>%
      mutate(
        post = post_2015,
        relative_year = relative_year_2015,
        cohort_treated = ifelse(cohort_2015 == 1, 1, 0)
      )
  } else if (cohort == "2016") {
    analysis_data <- data_processed %>% 
      filter(cohort_2016 == 1 | control_2016 == 1) %>%
      mutate(
        post = post_2016,
        relative_year = relative_year_2016,
        cohort_treated = ifelse(cohort_2016 == 1, 1, 0)
      )
  } else if (cohort == "all") {
    # 合并两个队列的数据
    analysis_data <- bind_rows(
      # 2015队列
      data_processed %>% 
        filter(cohort_2015 == 1 | control_2015 == 1) %>%
        mutate(
          post = post_2015,
          relative_year = relative_year_2015,
          cohort_treated = ifelse(cohort_2015 == 1, 1, 0),
          cohort_id = "2015"
        ),
      # 2016队列
      data_processed %>% 
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
        code_cohort = as.factor(paste(code, cohort_id, sep = "_")),
        # 统一的相对年份（基于各自队列）
        relative_year = ifelse(cohort_id == "2015", relative_year_2015, relative_year_2016)
      )
  } else {
    stop("必须指定具体的队列: '2015', '2016' 或 'all'")
  }
  
  # 计算滞后暴露
  if (cohort == "all") {
    # 对于合并队列，按城市-队列分组计算滞后
    lag_exposure <- calculate_lag_mean(analysis_data[[variable]], lag_value, analysis_data$code_cohort)
  } else {
    # 对于单个队列，按城市分组计算滞后
    lag_exposure <- calculate_lag_mean(analysis_data[[variable]], lag_value, analysis_data$code)
  }
  
  # 创建分析数据
  model_data <- analysis_data %>%
    mutate(
      exposure_lag = lag_exposure,
      exposure_binary = ifelse(exposure_lag > 0, 1, 0),
      # 创建相对年份因子，用于事件研究法
      rel_year_factor = factor(relative_year)
    ) %>%
    filter(!is.na(exposure_lag))
  
  # 重新设置基准年（政策前一年，t-1）
  model_data$rel_year_factor <- relevel(model_data$rel_year_factor, ref = "-1")
  
  # 初始化结果数据框
  time_windows <- -4:4
  results_df <- data.frame(
    relative_year = time_windows,
    cohort = cohort,
    outcome = outcome,
    variable = variable,
    lag = lag_value,
    coef = NA,
    se = NA,
    p_value = NA,
    rr = NA,
    rr_lower = NA,
    rr_upper = NA,
    n_obs = NA
  )
  
  # 使用三重交互项的事件研究法模型
  if (cohort == "all") {
    # 对于合并队列，使用城市-队列固定效应
    formula_es <- as.formula(paste0(
      outcome, " ~ exposure_binary * cohort_treated * rel_year_factor + ",
      "ns(temlag02, df = 6) + ns(rhlag02, df = 3) + ",
      "as.factor(holiday) + factor(code_cohort) + factor(year) + factor(month) + factor(dow)"
    ))
  } else {
    # 对于单个队列，使用城市固定效应
    formula_es <- as.formula(paste0(
      outcome, " ~ exposure_binary * cohort_treated * rel_year_factor + ",
      "ns(temlag02, df = 6) + ns(rhlag02, df = 3) + ",
      "as.factor(holiday) + factor(code) + factor(year) + factor(month) + factor(dow)"
    ))
  }
  
  # 拟合模型
  model <- glm.nb(formula_es, data = model_data)
  coef_summary <- summary(model)$coefficients
  
  # 提取每个相对年份的DID效应
  for (i in seq_along(time_windows)) {
    time_win <- time_windows[i]
    
    # 跳过基准年（t-1）
    if (time_win == -1) {
      results_df$coef[i] <- 0  # 基准年设为0
      results_df$se[i] <- NA
      results_df$p_value[i] <- NA
      results_df$rr[i] <- 1
      results_df$rr_lower[i] <- NA
      results_df$rr_upper[i] <- NA
      results_df$n_obs[i] <- nrow(model_data)
      next
    }
    
    # 构建交互项名称
    did_term <- paste0("exposure_binary:cohort_treated:rel_year_factor", time_win)
    
    if (did_term %in% rownames(coef_summary)) {
      results_df$coef[i] <- coef_summary[did_term, "Estimate"]
      results_df$se[i] <- coef_summary[did_term, "Std. Error"]
      results_df$p_value[i] <- coef_summary[did_term, "Pr(>|z|)"]
      results_df$rr[i] <- exp(results_df$coef[i])
      results_df$rr_lower[i] <- exp(results_df$coef[i] - 1.96 * results_df$se[i])
      results_df$rr_upper[i] <- exp(results_df$coef[i] + 1.96 * results_df$se[i])
      results_df$n_obs[i] <- nrow(model_data)
    } else {
      cat("未找到交互项:", did_term, "\n")
    }
  }
  
  # 保存结果
  write.csv(results_df, output_file, row.names = FALSE)
  cat("结果已保存至:", output_file, "\n")
  
  return(results_df)
}

# 运行事件研究法分析

# 批量分析函数
batch_event_study_analysis <- function() {
  outcomes <- c("case", "Pre_school", "School", "Adults", "Elders")
  variables <- c("flood", "severe_flood")
  cohorts <- c("2015", "2016", "all")
  lag_value <- 28
  
  all_results <- list()
  
  for(cohort in cohorts) {
    for(outcome in outcomes) {
      for(variable in variables) {
        output_file <- paste0("event_study_", cohort, "_", outcome, "_", variable, ".csv")
        
        cat("\n分析: cohort =", cohort, ", outcome =", outcome, ", variable =", variable, "\n")
        
        result <- event_study_analysis(
          cohort = cohort,
          outcome = outcome,
          variable = variable,
          lag_value = lag_value,
          output_file = output_file
        )
        
        key <- paste(cohort, outcome, variable, sep = "_")
        all_results[[key]] <- result
      }
    }
  }
  
  return(all_results)
}

# 执行批量分析
all_event_results <- batch_event_study_analysis()

###########事件研究法结果可视化####
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)
library(eoffice)  # 用于topptx函数

# Nature Water期刊配色
nature_colors <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")

# 专业主题设置
theme_nature <- function() {
  theme_minimal(base_family = "Arial") +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(size = rel(0.9)),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
}

# 1. 从批量分析结果中提取数据并合并
#setwd("F:\\毕业论文\\1011\\flood\\did\\Event_Study_Results\\case")
#all_event_results  <- data.table::fread("allcase.csv")

# 假设all_event_results是批量分析的结果列表
extract_event_results <- function(all_event_results) {
  # 提取case结局的结果
  case_results <- list()
  
  # 定义需要提取的组合
  combinations <- list(
    c("2015", "flood", "General Flood"),
    c("2016", "flood", "General Flood"), 
    c("all", "flood", "General Flood"),
    c("2015", "severe_flood", "Severe Flood"),
    c("2016", "severe_flood", "Severe Flood"),
    c("all", "severe_flood", "Severe Flood")
  )
  
  for(comb in combinations) {
    cohort <- comb[1]
    variable <- comb[2]
    flood_type <- comb[3]
    
    key <- paste(cohort, "case", variable, sep = "_")
    if(key %in% names(all_event_results)) {
      result_df <- all_event_results[[key]] %>%
        mutate(flood_type = flood_type)
      case_results[[key]] <- result_df
    }
  }
  
  # 合并所有结果
  all_results <- bind_rows(case_results) %>%
    mutate(
      cohort = factor(cohort, levels = c("all","2015", "2016")),
      significance = case_when(
        p_value < 0.05 ~ "*",
        TRUE ~ "NS"
      ),
      # 创建时间标签
      time_label = case_when(
        relative_year == -4 ~ "t-4",
        relative_year == -3 ~ "t-3", 
        relative_year == -2 ~ "t-2",
        relative_year == -1 ~ "t-1",
        relative_year == 0 ~ "t0",
        relative_year == 1 ~ "t+1",
        relative_year == 2 ~ "t+2",
        relative_year == 3 ~ "t+3",
        relative_year == 4 ~ "t+4"
      ),
      time_label = factor(time_label, levels = c("t-4", "t-3", "t-2", "t-1", "t0", "t+1", "t+2", "t+3", "t+4"))
    )
  
  return(all_results)
}

# 从批量分析结果中提取数据
all_results <- extract_event_results(all_event_results)

# 2. 创建事件研究轨迹图（三个队列）
create_event_study_plot <- function(data, flood_type_sel) {
  plot_data <- data %>% filter(flood_type == flood_type_sel)
  
  ggplot(plot_data, aes(x = time_label, y = rr, color = cohort, group = cohort)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
    # 置信区间
    geom_errorbar(aes(ymin = rr_lower, ymax = rr_upper), 
                  width = 0.2, alpha = 0.6, position = position_dodge(0.4)) +
    # 效应轨迹线
    geom_line(linewidth = 0.8, position = position_dodge(0.4)) +
    # 点估计
    geom_point(aes(shape = significance), size = 2.5, position = position_dodge(0.4)) +
    # 美学设置
    scale_color_manual(
      values = c("all" = nature_colors[3],"2015" = nature_colors[1], 
                 "2016" = nature_colors[2]
                 ),
      labels = c("all" = "All","2015" = "2015", 
                 "2016" = "2016"
                 )
    ) +
    scale_shape_manual(values = c("*" = 16, "NS" = 1)) +
    labs(
      title = paste(flood_type_sel),
      x = "Relative Year",
      y = "Risk Ratio (DID)",
      color = "Implementation Cohort",
      shape = "Significance"
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}
p1_general <- create_event_study_plot(all_results, "General Flood")
p1_general 
p1_severe <- create_event_study_plot(all_results, "Severe Flood")
p1_severe
# 组合事件研究图
combined_events <- (p1_general + p1_severe) + 
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = 'A')
combined_events
topptx(combined_events, filename = "F:\\Flood\\did\\DID_Results\\Events_combined.pptx", width = 8, height = 4)


# 扩展的提取函数，支持所有年龄组
extract_age_group_results <- function(all_event_results) {
  # 定义所有需要分析的年龄组
  age_groups <- c("Pre_school", "School", "Adults", "Elders")
  
  all_age_results <- list()
  
  for(age_group in age_groups) {
    age_results <- list()
    
    # 定义需要提取的组合
    combinations <- list(
      c("2015", "flood", "General Flood"),
      c("2016", "flood", "General Flood"), 
      c("all", "flood", "General Flood"),
      c("2015", "severe_flood", "Severe Flood"),
      c("2016", "severe_flood", "Severe Flood"),
      c("all", "severe_flood", "Severe Flood")
    )
    
    for(comb in combinations) {
      cohort <- comb[1]
      variable <- comb[2]
      flood_type <- comb[3]
      
      key <- paste(cohort, age_group, variable, sep = "_")
      if(key %in% names(all_event_results)) {
        result_df <- all_event_results[[key]] %>%
          mutate(
            flood_type = flood_type,
            age_group = age_group
          )
        age_results[[key]] <- result_df
      }
    }
    
    # 合并当前年龄组的结果
    if(length(age_results) > 0) {
      age_group_df <- bind_rows(age_results) %>%
        mutate(
          cohort = factor(cohort, levels = c("all", "2015", "2016")),
          significance = case_when(
            p_value < 0.05 ~ "*",
            TRUE ~ "NS"
          ),
          # 创建时间标签
          time_label = case_when(
            relative_year == -4 ~ "t-4",
            relative_year == -3 ~ "t-3", 
            relative_year == -2 ~ "t-2",
            relative_year == -1 ~ "t-1",
            relative_year == 0 ~ "t0",
            relative_year == 1 ~ "t+1",
            relative_year == 2 ~ "t+2",
            relative_year == 3 ~ "t+3",
            relative_year == 4 ~ "t+4"
          ),
          time_label = factor(time_label, levels = c("t-4", "t-3", "t-2", "t-1", "t0", "t+1", "t+2", "t+3", "t+4"))
        )
      
      all_age_results[[age_group]] <- age_group_df
    }
  }
  
  # 合并所有年龄组的结果
  final_results <- bind_rows(all_age_results) %>%
    mutate(
      age_group = factor(age_group, 
                         levels = c("Pre_school", "School", "Adults", "Elders"),
                         labels = c("Pre-school", "School-age", "Adults", "Elderly"))
    )
  
  return(final_results)
}

# 从批量分析结果中提取所有年龄组数据
all_age_results <- extract_age_group_results(all_event_results)

# 3. 创建年龄组特定的事件研究轨迹图
create_age_group_event_plot <- function(data, age_group_sel, flood_type_sel) {
  plot_data <- data %>% 
    filter(age_group == age_group_sel, flood_type == flood_type_sel)
  
  # 检查数据是否存在
  if(nrow(plot_data) == 0) {
    return(ggplot() + 
             geom_text(aes(x = 0.5, y = 0.5, label = "No data available"), 
                       size = 6) + 
             theme_void() +
             labs(title = paste(age_group_sel, "-", flood_type_sel)))
  }
  
  ggplot(plot_data, aes(x = time_label, y = rr, color = cohort, group = cohort)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
    # 置信区间
    geom_errorbar(aes(ymin = rr_lower, ymax = rr_upper), 
                  width = 0.2, alpha = 0.6, position = position_dodge(0.4)) +
    # 效应轨迹线
    geom_line(linewidth = 0.8, position = position_dodge(0.4)) +
    # 点估计
    geom_point(aes(shape = significance), size = 2, position = position_dodge(0.4)) +
    # 美学设置
    scale_color_manual(
      values = c("all" = nature_colors[3], "2015" = nature_colors[1], 
                 "2016" = nature_colors[2]),
      labels = c("all" = "All", "2015" = "2015", 
                 "2016" = "2016")
    ) +
    scale_shape_manual(values = c("*" = 16, "NS" = 1)) +
    labs(
      title = paste(age_group_sel, "-", flood_type_sel),
      x = "Relative Year",
      y = "Risk Ratio (DID)",
      color = "Implementation Cohort",
      shape = "Significance"
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title = element_text(size = 10)
    ) +
    # 调整y轴范围以适应不同年龄组的数据范围
    coord_cartesian(ylim = c(0, 2.5))  # 可根据实际数据调整
}

# 4. 创建所有年龄组的组合图
create_all_age_groups_plot <- function(data) {
  # 定义年龄组和洪水类型
  age_groups <- c("Pre-school", "School-age", "Adults", "Elderly")
  flood_types <- c("General Flood", "Severe Flood")
  
  # 创建所有子图
  plot_list <- list()
  
  for(age_group in age_groups) {
    for(flood_type in flood_types) {
      plot_name <- paste(age_group, flood_type, sep = "_")
      plot_list[[plot_name]] <- create_age_group_event_plot(data, age_group, flood_type)
    }
  }
  
  # 组合所有子图
  combined_plot <- wrap_plots(plot_list, ncol = 2, nrow = 4) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  return(combined_plot)
}
# 方法1：所有年龄组的大组合图（4行×2列）
all_age_groups_combined <- create_all_age_groups_plot(all_age_results)
all_age_groups_combined
topptx(all_age_groups_combined, filename = "F:\\Flood\\did\\DID_Results\\ageEvents_combined.pptx",
       width = 8, height =8)


