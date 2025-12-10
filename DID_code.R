# 读取数据
library(tidyverse)
library(data.table)
library(MASS)
library(splines)
library(doParallel)
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
    control_2015 = ifelse(policy_year == "Non_Exposed_2015", 1, 0),  # 控制组可用于两个队列
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
# ==================== 简化的滞后计算函数 ====================
calculate_lag_mean <- function(x, lag, group) {
  require(data.table)
  dt <- data.table(x = x, group = group)
  dt[, lag_mean := frollmean(x, n = lag + 1, align = "right", na.rm = TRUE), by = group]
  return(dt$lag_mean)
}

# 计算气象变量滞后
data_processed_corrected$temlag02 <- calculate_lag_mean(data_processed_corrected$Tmean, 2, data_processed_corrected$code)
data_processed_corrected$rhlag02 <- calculate_lag_mean(data_processed_corrected$RH, 2, data_processed_corrected$code)

# ==================== 多期DID分析函数 ====================
corrected_did_analysis <- function(outcome, variable, cohort = "2015", output_file, max_lag = 28) {
  
  cat("\n=== 开始多期DID分析 ===\n")
  cat("结局变量:", outcome, "\n")
  cat("暴露变量:", variable, "\n")
  cat("分析队列:", cohort, "\n")
  
  # 选择对应队列的数据
  if (cohort == "2015") {
    analysis_data <- data_processed_corrected %>% 
      filter(cohort_2015 == 1 | control_2015 == 1) %>%
      mutate(
        post = post_2015,
        relative_year = relative_year_2015,
        cohort_treated = ifelse(cohort_2015 == 1, 1, 0)  # 该队列的处理组
      )
  } else if (cohort == "2016") {
    analysis_data <- data_processed_corrected %>% 
      filter(cohort_2016 == 1 | control_2016 == 1) %>%
      mutate(
        post = post_2016,
        relative_year = relative_year_2016,
        cohort_treated = ifelse(cohort_2016 == 1, 1, 0)  # 该队列的处理组
      )
  } else {
    stop("必须指定具体的队列: '2015' 或 '2016'")
  }
  
  # 创建相对年份因子，用于事件研究法
  analysis_data$rel_year_factor <- factor(analysis_data$relative_year)
  analysis_data$rel_year_factor <- relevel(analysis_data$rel_year_factor, ref = "-1")
  
  # 初始化结果数据框
  results_df <- data.frame(
    lag = 0:max_lag,
    cohort = cohort,
    outcome = outcome,
    variable = variable,
    coef_main = NA,
    se_main = NA,
    p_main = NA,
    rr_main = NA,
    rr_main_lower = NA,
    rr_main_upper = NA,
    coef_did = NA,
    se_did = NA,
    p_did = NA,
    rr_did = NA,
    rr_did_lower = NA,
    rr_did_upper = NA,
    aic = NA,
    n_obs = NA
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
      
      # 使用队列特定的处理组变量
      # 标准三重交互项模型
      formula_did <- as.formula(paste0(
        outcome, " ~ exposure_binary * cohort_treated * post + ",
        "ns(temlag02, df = 6) + ns(rhlag02, df = 3) + ",
        "as.factor(holiday) + factor(code) + factor(year) + factor(month) + factor(dow)"
      ))
      
      # 事件研究法模型（推荐）
      formula_es <- as.formula(paste0(
        outcome, " ~ exposure_binary * cohort_treated * rel_year_factor + ",
        "ns(temlag02, df = 6) + ns(rhlag02, df = 3) + ",
        "as.factor(holiday) + factor(code) + factor(year) + factor(month) + factor(dow)"
      ))
      
      # 拟合模型（使用事件研究法）
      model <- glm.nb(formula_es, data = model_data)
      
      # 提取系数
      coef_summary <- summary(model)$coefficients
      
      # 主效应
      if("exposure_binary" %in% rownames(coef_summary)) {
        results_df$coef_main[lag + 1] <- coef_summary["exposure_binary", "Estimate"]
        results_df$se_main[lag + 1] <- coef_summary["exposure_binary", "Std. Error"]
        results_df$p_main[lag + 1] <- coef_summary["exposure_binary", "Pr(>|z|)"]
        results_df$rr_main[lag + 1] <- exp(results_df$coef_main[lag + 1])
        results_df$rr_main_lower[lag + 1] <- exp(results_df$coef_main[lag + 1] - 1.96 * results_df$se_main[lag + 1])
        results_df$rr_main_upper[lag + 1] <- exp(results_df$coef_main[lag + 1] + 1.96 * results_df$se_main[lag + 1])
      }
      
      # 提取DID效应（政策后第一年的三重交互项）
      did_term <- "exposure_binary:cohort_treated:rel_year_factor0"  # 政策后第一年
      if(did_term %in% rownames(coef_summary)) {
        results_df$coef_did[lag + 1] <- coef_summary[did_term, "Estimate"]
        results_df$se_did[lag + 1] <- coef_summary[did_term, "Std. Error"]
        results_df$p_did[lag + 1] <- coef_summary[did_term, "Pr(>|z|)"]
        results_df$rr_did[lag + 1] <- exp(results_df$coef_did[lag + 1])
        results_df$rr_did_lower[lag + 1] <- exp(results_df$coef_did[lag + 1] - 1.96 * results_df$se_did[lag + 1])
        results_df$rr_did_upper[lag + 1] <- exp(results_df$coef_did[lag + 1] + 1.96 * results_df$se_did[lag + 1])
      } else {
        # 备选：使用标准三重交互项
        did_term_alt <- "exposure_binary:cohort_treated:post"
        if(did_term_alt %in% rownames(coef_summary)) {
          results_df$coef_did[lag + 1] <- coef_summary[did_term_alt, "Estimate"]
          results_df$se_did[lag + 1] <- coef_summary[did_term_alt, "Std. Error"]
          results_df$p_did[lag + 1] <- coef_summary[did_term_alt, "Pr(>|z|)"]
          results_df$rr_did[lag + 1] <- exp(results_df$coef_did[lag + 1])
          results_df$rr_did_lower[lag + 1] <- exp(results_df$coef_did[lag + 1] - 1.96 * results_df$se_did[lag + 1])
          results_df$rr_did_upper[lag + 1] <- exp(results_df$coef_did[lag + 1] + 1.96 * results_df$se_did[lag + 1])
        }
      }
      
      # 模型拟合统计
      results_df$aic[lag + 1] <- AIC(model)
      results_df$n_obs[lag + 1] <- nrow(model$model)
      
    }, error = function(e) {
      cat("\nLag", lag, "出错:", as.character(e), "\n")
    })
  }
  
  cat("\n分析完成！\n")
  
  # 保存结果
  write.csv(results_df, output_file, row.names = FALSE)
  cat("结果已保存至:", output_file, "\n")
  
  return(results_df)
}

# ==================== 全部数据分析函数 ====================

analyze_all_cohorts_direct <- function(outcome, variable, output_file, max_lag = 28) {
  
  cat("\n=== 开始直接分析全部数据（2015+2016队列合并） ===\n")
  cat("结局变量:", outcome, "\n")
  cat("暴露变量:", variable, "\n")
  cat("分析队列: all\n")
  
  # 合并两个队列的数据
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
  
  # 创建相对年份因子
  analysis_data_all$rel_year_factor <- factor(analysis_data_all$relative_year)
  analysis_data_all$rel_year_factor <- relevel(analysis_data_all$rel_year_factor, ref = "-1")
  
  # 初始化结果数据框
  results_df <- data.frame(
    lag = 0:max_lag,
    cohort = "all",
    outcome = outcome,
    variable = variable,
    coef_main = NA,
    se_main = NA,
    p_main = NA,
    rr_main = NA,
    rr_main_lower = NA,
    rr_main_upper = NA,
    coef_did = NA,
    se_did = NA,
    p_did = NA,
    rr_did = NA,
    rr_did_lower = NA,
    rr_did_upper = NA,
    aic = NA,
    n_obs = NA,
    n_treated = NA,
    n_control = NA
  )
  
  # 逐个lag分析
  for(lag in 0:max_lag) {
    
    cat("\r处理 Lag", lag, "/", max_lag, "...")
    
    tryCatch({
      # 计算滞后暴露
      lag_exposure <- calculate_lag_mean(analysis_data_all[[variable]], lag, analysis_data_all$code)
      
      # 创建分析数据
      model_data <- analysis_data_all %>%
        mutate(
          exposure_lag = lag_exposure,
          exposure_binary = ifelse(exposure_lag > 0, 1, 0)
        ) %>%
        filter(!is.na(exposure_lag))
      
      # 使用Stacked DID模型
      formula_es <- as.formula(paste0(
        outcome, " ~ exposure_binary * cohort_treated * rel_year_factor + ",
        "ns(temlag02, df = 6) + ns(rhlag02, df = 3) + ",
        "as.factor(holiday) + factor(year) + factor(month) + factor(dow) + ",
        "factor(code_cohort)"  # 城市-队列固定效应
      ))
      
      # 拟合模型
      model <- glm.nb(formula_es, data = model_data)
      
      # 提取系数
      coef_summary <- summary(model)$coefficients
      
      # 主效应
      if("exposure_binary" %in% rownames(coef_summary)) {
        results_df$coef_main[lag + 1] <- coef_summary["exposure_binary", "Estimate"]
        results_df$se_main[lag + 1] <- coef_summary["exposure_binary", "Std. Error"]
        results_df$p_main[lag + 1] <- coef_summary["exposure_binary", "Pr(>|z|)"]
        results_df$rr_main[lag + 1] <- exp(results_df$coef_main[lag + 1])
        results_df$rr_main_lower[lag + 1] <- exp(results_df$coef_main[lag + 1] - 1.96 * results_df$se_main[lag + 1])
        results_df$rr_main_upper[lag + 1] <- exp(results_df$coef_main[lag + 1] + 1.96 * results_df$se_main[lag + 1])
      }
      
      # 提取DID效应（政策后第一年的三重交互项）
      did_term <- "exposure_binary:cohort_treated:rel_year_factor0"  # 政策后第一年
      if(did_term %in% rownames(coef_summary)) {
        results_df$coef_did[lag + 1] <- coef_summary[did_term, "Estimate"]
        results_df$se_did[lag + 1] <- coef_summary[did_term, "Std. Error"]
        results_df$p_did[lag + 1] <- coef_summary[did_term, "Pr(>|z|)"]
        results_df$rr_did[lag + 1] <- exp(results_df$coef_did[lag + 1])
        results_df$rr_did_lower[lag + 1] <- exp(results_df$coef_did[lag + 1] - 1.96 * results_df$se_did[lag + 1])
        results_df$rr_did_upper[lag + 1] <- exp(results_df$coef_did[lag + 1] + 1.96 * results_df$se_did[lag + 1])
      } else {
        # 备选：使用标准三重交互项
        did_term_alt <- "exposure_binary:cohort_treated:post"
        if(did_term_alt %in% rownames(coef_summary)) {
          results_df$coef_did[lag + 1] <- coef_summary[did_term_alt, "Estimate"]
          results_df$se_did[lag + 1] <- coef_summary[did_term_alt, "Std. Error"]
          results_df$p_did[lag + 1] <- coef_summary[did_term_alt, "Pr(>|z|)"]
          results_df$rr_did[lag + 1] <- exp(results_df$coef_did[lag + 1])
          results_df$rr_did_lower[lag + 1] <- exp(results_df$coef_did[lag + 1] - 1.96 * results_df$se_did[lag + 1])
          results_df$rr_did_upper[lag + 1] <- exp(results_df$coef_did[lag + 1] + 1.96 * results_df$se_did[lag + 1])
        }
      }
      
      # 模型拟合统计
      results_df$aic[lag + 1] <- AIC(model)
      results_df$n_obs[lag + 1] <- nrow(model$model)
      results_df$n_treated[lag + 1] <- sum(model_data$cohort_treated)
      results_df$n_control[lag + 1] <- sum(1 - model_data$cohort_treated)
      
    }, error = function(e) {
      cat("\nLag", lag, "出错:", as.character(e), "\n")
    })
  }
  
  cat("\n全部数据分析完成！\n")
  
  # 保存结果
  write.csv(results_df, output_file, row.names = FALSE)
  cat("结果已保存至:", output_file, "\n")
  
  return(results_df)
}

# ==================== 批量执行所有组合的分析 ====================

# 定义分析参数
outcomes <- c("case", "Pre_school", "School", "Adults", "Elders")
variables <- c("flood", "severe_flood")
cohorts <- c("2015", "2016", "all")  # 包含全部队列

# 创建结果目录
if(!dir.exists("DID_Results")) {
  dir.create("DID_Results")
}

# 执行批量分析
all_results <- list()

for(cohort in cohorts) {
  for(outcome in outcomes) {
    for(variable in variables) {
      
      output_file <- paste0("DID_Results/DID_", cohort, "_", outcome, "_", variable, ".csv")
      
      cat("\n", paste0(rep("=", 60), collapse = ""), "\n")
      cat("开始分析: cohort =", cohort, ", outcome =", outcome, ", variable =", variable, "\n")
      
      if (cohort == "all") {
        # 使用直接合并方法分析全部数据
        result <- analyze_all_cohorts_direct(
          outcome = outcome,
          variable = variable,
          output_file = output_file,
          max_lag = 28
        )
      } else {
        # 使用原有函数分析单个队列
        result <- corrected_did_analysis(
          outcome = outcome,
          variable = variable,
          cohort = cohort,
          output_file = output_file,
          max_lag = 28
        )
      }
      
      # 存储结果
      key <- paste(cohort, outcome, variable, sep = "_")
      all_results[[key]] <- result
      
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

# ==================== 生成汇总报告 ====================

generate_summary_report <- function(all_results) {
  
  summary_df <- data.frame()
  
  for(key in names(all_results)) {
    result <- all_results[[key]]
    
    # 计算统计量
    n_valid <- sum(!is.na(result$coef_did))
    n_total <- nrow(result)
    sig_count <- sum(result$p_did < 0.05, na.rm = TRUE)
    
    # 找到最显著的滞后
    if(n_valid > 0) {
      min_p_idx <- which.min(result$p_did)
      best_lag <- result$lag[min_p_idx]
      best_rr <- result$rr_did[min_p_idx]
      best_p <- result$p_did[min_p_idx]
    } else {
      best_lag <- best_rr <- best_p <- NA
    }
    
    # 解析key
    parts <- strsplit(key, "_")[[1]]
    cohort <- parts[1]
    outcome <- parts[2]
    variable <- parts[3]
    
    summary_df <- rbind(summary_df, data.frame(
      cohort = cohort,
      outcome = outcome,
      variable = variable,
      n_valid = n_valid,
      n_total = n_total,
      completion_rate = round(n_valid/n_total * 100, 1),
      significant_count = sig_count,
      best_lag = best_lag,
      best_rr = round(best_rr, 3),
      best_p = round(best_p, 4),
      stringsAsFactors = FALSE
    ))
  }
  
  # 保存汇总报告
  write.csv(summary_df, "DID_Results/analysis_summary_report.csv", row.names = FALSE)
  
  cat("\n", paste0(rep("=", 60), collapse = ""), "\n")
  cat("分析汇总报告:\n")
  print(summary_df)
  
  # 按队列统计
  cat("\n按队列统计:\n")
  cohort_summary <- summary_df %>%
    group_by(cohort) %>%
    summarise(
      n_analyses = n(),
      avg_completion = round(mean(completion_rate), 1),
      n_significant = sum(significant_count > 0),
      avg_best_rr = round(mean(best_rr, na.rm = TRUE), 3)
    )
  print(cohort_summary)
  
  return(summary_df)
}

# 生成汇总报告
summary_report <- generate_summary_report(all_results)

