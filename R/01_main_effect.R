# =============================================================================
# 01_main_effect.R
# Main short-term association analysis: flood -> bacillary dysentery (BD).
# Time-stratified case-crossover design + distributed-lag non-linear model
# (DLNM), fit as a conditional negative-binomial model via gnm (stratum =
# city:year:month:day-of-week). Reproduces manuscript equation (1) and its
# sensitivity (no-precipitation, lag windows) and subgroup (region, climate
# zone, age group) analyses.
#
# Input : data/daily_panel.csv  (daily city panel)
# Output: results/<run_id>/...       (single-lag & cumulative RR tables)
# Run from the repository root:  Rscript R/01_main_effect.R
# =============================================================================

if (!require("pacman")) install.packages("pacman", repos = "https://cloud.r-project.org")
pacman::p_load(data.table, foreach, doParallel, gnm, MASS, dlnm, splines, parallel)

resolve_project_root <- function(default = ".") {
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(cwd, "data", "daily_panel.csv"))) {
    return(cwd)
  }
  
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = TRUE)
    for (candidate in c(dirname(script_path), file.path(dirname(script_path), ".."))) {
      candidate <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
      if (file.exists(file.path(candidate, "data", "daily_panel.csv"))) {
        return(candidate)
      }
    }
  }
  
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(active_path)) {
      for (candidate in c(dirname(active_path), file.path(dirname(active_path), ".."))) {
        candidate <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
        if (file.exists(file.path(candidate, "data", "daily_panel.csv"))) {
          return(candidate)
        }
      }
    }
  }
  
  normalizePath(default, winslash = "/", mustWork = TRUE)
}

# ===========================  USER CONFIGURATION  ============================
project_root <- resolve_project_root()
input_file <- file.path("data", "daily_panel.csv")
output_root <- "results"
# run_mode (analysis "tier"), set via the RUN_MODE environment variable; controls
# how much of the main-effect stage runs:
#   "smoke"  — load packages and read/preprocess the data only (no model fitting)
#   "pilot"  — fast: main analysis, total cases x two flood definitions (seconds)
#   "middle" — all outcomes (total + 4 age groups) x all 8 flood definitions, plus
#              the no-precipitation and lag-window (7/14/21/28) sensitivity analyses
#   "full"   — everything in "middle" PLUS the region- and climate-zone subgroup
#              analyses (~15 min on the simulated data) -- the default
# Default "full" so run_all.R reproduces the complete analysis; set RUN_MODE=pilot
# (or middle) for a faster demo.
run_mode <- Sys.getenv("RUN_MODE", "full")
requested_cores <- 4L   # parallel workers; raise on a multi-core machine
run_label <- "main-effect"
# ================================================

# ============================  UTILITY FUNCTIONS  ============================
# Write a data.table to CSV as UTF-8 with a byte-order mark, so spreadsheet
# software detects the encoding when results files are opened directly.
write_csv_utf8 <- function(x, file) {
  fwrite(as.data.table(x), file = file, bom = TRUE)
}

# Coerce to numeric, tolerating factor/character columns (unparseable -> NA).
to_numeric_safe <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(as.character(x)))
}

# Quantile of x ignoring NAs; returns NA for an all-missing vector.
safe_quantile <- function(x, p) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

# Within-group moving average of x over the given lags (e.g. 0:2 = 3-day mean),
# computed separately for each city so it never averages across city boundaries.
runMean <- function(x, lags, group) {
  dt <- data.table(x = x, group = group)
  dt[, avg := {
    lag_sum <- Reduce(`+`, lapply(lags, function(lag) shift(x, n = lag, fill = NA)))
    lag_sum / length(lags)
  }, by = group]
  dt$avg
}

# =========================  MODEL-FITTING FUNCTIONS  =========================
# Start a PSOCK worker cluster (portable across operating systems) and pin each
# worker to a single BLAS / data.table thread, so nested parallelism does not
# oversubscribe the cores.
prepare_cluster_simple <- function(num_cores) {
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  invisible(clusterCall(cl, function() {
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      BLIS_NUM_THREADS = "1",
      VECLIB_MAXIMUM_THREADS = "1"
    )
    if (requireNamespace("data.table", quietly = TRUE)) {
      data.table::setDTthreads(1L)
    }
    NULL
  }))
  
  invisible(clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(data.table)
      library(foreach)
      library(doParallel)
      library(gnm)
      library(MASS)
      library(dlnm)
      library(splines)
      library(parallel)
    })
    NULL
  }))
  
  cl
}

# Maximum-likelihood estimate of the negative-binomial dispersion (theta) at the
# given fitted means; returns NA if the estimate does not converge.
estimate_theta_ml <- function(y, mu, weights = NULL, limit = 100L, eps = 1e-08) {
  if (is.null(weights)) weights <- rep(1, length(y))
  mu <- pmax(as.numeric(mu), eps)
  theta_est <- tryCatch(
    suppressWarnings(
      theta.ml(
        y = as.numeric(y),
        mu = mu,
        n = sum(weights),
        weights = weights,
        limit = as.integer(limit),
        trace = FALSE
      )
    ),
    error = function(e) NA_real_
  )
  if (!is.finite(theta_est) || theta_est <= 0) return(NA_real_)
  as.numeric(theta_est)
}

# Subset a DLNM cross-basis to selected rows while preserving its attributes
# (df, lag, basis arguments) that plain matrix indexing would otherwise drop.
subset_crossbasis_rows <- function(cb, keep_rows) {
  cb_subset <- cb[keep_rows, , drop = FALSE]
  structure(
    cb_subset,
    df = attr(cb, "df"),
    range = attr(cb, "range"),
    lag = attr(cb, "lag"),
    argvar = attr(cb, "argvar"),
    arglag = attr(cb, "arglag"),
    class = class(cb)
  )
}

# Cumulative-lag windows to report: the candidate windows that fit within max_lag,
# always including max_lag itself.
resolve_cumulative_windows <- function(max_lag, candidate_windows = c(7L, 14L, 21L, 28L)) {
  max_lag <- as.integer(max_lag)
  windows <- as.integer(candidate_windows[candidate_windows <= max_lag])
  if (!max_lag %in% windows) windows <- c(windows, max_lag)
  sort(unique(windows))
}

# First row of a crosspred matrix as a plain numeric vector (RR by lag at the
# unit exposure contrast).
extract_pred_row <- function(x) {
  if (is.null(dim(x))) return(as.numeric(x))
  as.numeric(x[1, , drop = TRUE])
}

# Fit the stratified conditional count model: a Poisson gnm first, then iteratively
# re-estimate the negative-binomial dispersion (theta) to convergence, refitting at
# each step; fall back to Poisson if the negative-binomial fit fails.
fit_nb_gnm_model <- function(dt, fm, outcome, model_maxit = 500L, theta_maxit = 8L, theta_tol = 1e-06) {
  y <- dt[[outcome]]
  weights <- rep(1, length(y))
  
  fit_poisson <- function() {
    tryCatch(
      gnm(
        fm,
        family = poisson,
        eliminate = factor(stratum),
        data = dt,
        control = gnmControl(maxit = model_maxit, trace = FALSE)
      ),
      error = function(e) NULL
    )
  }
  
  fit_negbin <- function(theta_value) {
    tryCatch(
      gnm(
        fm,
        family = negative.binomial(theta_value),
        eliminate = factor(stratum),
        data = dt,
        control = gnmControl(maxit = model_maxit, trace = FALSE)
      ),
      error = function(e) NULL
    )
  }
  
  poisson_fit <- fit_poisson()
  if (is.null(poisson_fit)) {
    return(list(
      fit = NULL,
      family_used = "poisson",
      theta_est = NA_real_,
      theta_converged = FALSE,
      theta_iterations = 0L
    ))
  }
  
  theta_current <- tryCatch(
    suppressWarnings(as.numeric(glm.nb(fm, data = dt)$theta)),
    error = function(e) NA_real_
  )
  if (!is.finite(theta_current) || theta_current <= 0) {
    theta_current <- estimate_theta_ml(y = y, mu = fitted(poisson_fit), weights = weights)
  }
  if (!is.finite(theta_current) || theta_current <= 0) {
    return(list(
      fit = poisson_fit,
      family_used = "poisson",
      theta_est = NA_real_,
      theta_converged = FALSE,
      theta_iterations = 0L
    ))
  }
  
  theta_converged <- FALSE
  theta_iterations <- 0L
  nb_fit <- NULL
  
  for (iter in seq_len(theta_maxit)) {
    theta_iterations <- iter
    nb_fit <- fit_negbin(theta_current)
    if (is.null(nb_fit)) break
    theta_next <- estimate_theta_ml(y = y, mu = fitted(nb_fit), weights = weights)
    if (!is.finite(theta_next) || theta_next <= 0) break
    if (abs(theta_next - theta_current) <= theta_tol * max(1, abs(theta_current))) {
      theta_current <- theta_next
      theta_converged <- TRUE
      break
    }
    theta_current <- theta_next
  }
  
  if (is.null(nb_fit)) {
    return(list(
      fit = poisson_fit,
      family_used = "poisson",
      theta_est = NA_real_,
      theta_converged = FALSE,
      theta_iterations = 0L
    ))
  }
  
  final_nb_fit <- fit_negbin(theta_current)
  if (is.null(final_nb_fit)) final_nb_fit <- nb_fit
  
  list(
    fit = final_nb_fit,
    family_used = "negbin",
    theta_est = theta_current,
    theta_converged = theta_converged,
    theta_iterations = theta_iterations
  )
}

# Core estimator for one outcome x exposure: build the distributed-lag cross-basis
# for flood, keep complete cases, check the exposure is frequent enough to identify
# a model, fit the conditional NB-DLNM, and return single-lag and cumulative RR by lag.
fit_dlm_crossbasis <- function(dt, outcome, variable,
                               include_pre = TRUE,
                               max_lag = 28L,
                               lag_df = 5L,
                               min_exposed_days = 10L,
                               min_exposed_strata = 10L,
                               min_exposed_positive_rows = 5L,
                               min_exposed_outcome_sum = 10L) {
  needed_cols <- c(outcome, variable, "holiday", "temlag02", "rhlag02", "stratum", "code")
  if (include_pre && "prelag02" %in% names(dt) && !all(is.na(dt$prelag02))) {
    needed_cols <- c(needed_cols, "prelag02")
  }
  missing_cols <- setdiff(needed_cols, names(dt))
  if (length(missing_cols) > 0L) {
    return(list(single = NULL, cumulative = NULL, error = paste("missing columns:", paste(missing_cols, collapse = ", "))))
  }
  
  lag_df_eff <- max(1L, min(as.integer(lag_df), as.integer(max_lag)))
  cb_error <- NULL
  cb_full <- tryCatch(
    crossbasis(
      x = dt[[variable]],
      lag = max_lag,
      argvar = list(fun = "lin"),
      arglag = list(fun = "ns", df = lag_df_eff),
      group = droplevels(factor(dt$code))
    ),
    error = function(e) {
      cb_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(cb_full)) {
    return(list(single = NULL, cumulative = NULL, error = paste0("crossbasis failed: ", cb_error)))
  }
  
  complete_mask <- complete.cases(dt[, ..needed_cols]) & complete.cases(cb_full)
  dt_cc <- dt[complete_mask]
  cb <- subset_crossbasis_rows(cb_full, complete_mask)
  if (nrow(dt_cc) == 0L) {
    return(list(single = NULL, cumulative = NULL, error = "No complete cases after filtering."))
  }
  
  exposed_days <- dt_cc[[variable]] == 1
  n_exposed_days <- sum(exposed_days, na.rm = TRUE)
  n_exposed_strata <- uniqueN(dt_cc$stratum[exposed_days])
  n_exposed_positive_rows <- sum(exposed_days & dt_cc[[outcome]] > 0, na.rm = TRUE)
  exposed_outcome_sum <- sum(dt_cc[[outcome]][exposed_days], na.rm = TRUE)
  
  diag_stats <- list(
    n_complete = nrow(dt_cc),
    n_exposed_days = n_exposed_days,
    n_exposed_strata = n_exposed_strata,
    n_exposed_positive_rows = n_exposed_positive_rows,
    exposed_outcome_sum = exposed_outcome_sum
  )
  
  if (n_exposed_days < min_exposed_days ||
      n_exposed_strata < min_exposed_strata ||
      n_exposed_positive_rows < min_exposed_positive_rows ||
      exposed_outcome_sum < min_exposed_outcome_sum) {
    return(list(single = NULL, cumulative = NULL, error = "Insufficient data (exposed days/strata/cases).", diag = diag_stats))
  }
  
  base_vars <- c("ns(temlag02, 6)", "as.factor(holiday)", "ns(rhlag02, 3)")
  if (include_pre && "prelag02" %in% names(dt_cc) && !all(is.na(dt_cc$prelag02))) {
    base_vars <- c(base_vars, "ns(prelag02, 3)")
  }
  fm <- as.formula(paste(outcome, "~ cb +", paste(base_vars, collapse = " + ")))
  
  model_fit <- fit_nb_gnm_model(dt = dt_cc, fm = fm, outcome = outcome)
  final_fit <- model_fit$fit
  if (is.null(final_fit)) {
    return(list(single = NULL, cumulative = NULL, error = "Model fitting failed.", diag = diag_stats))
  }
  
  pred_error <- NULL
  pred <- tryCatch(
    crosspred(cb, final_fit, at = 1, cen = 0, cumul = TRUE, model.link = "log"),
    error = function(e) {
      pred_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(pred)) {
    return(list(single = NULL, cumulative = NULL, error = paste0("crosspred failed: ", pred_error), diag = diag_stats))
  }
  
  lag_index <- 0:max_lag
  single_res <- data.table(
    lag = lag_index,
    rr = extract_pred_row(pred$matRRfit),
    rr_lower = extract_pred_row(pred$matRRlow),
    rr_upper = extract_pred_row(pred$matRRhigh),
    coef = extract_pred_row(pred$matfit),
    se = extract_pred_row(pred$matse)
  )
  
  if (nrow(single_res) != length(lag_index)) {
    return(list(single = NULL, cumulative = NULL, error = "Unexpected single-lag prediction dimensions.", diag = diag_stats))
  }
  
  cum_windows <- resolve_cumulative_windows(max_lag)
  cum_idx <- match(cum_windows, lag_index)
  cumulative_res <- data.table(
    interval = paste0("lag0-", cum_windows),
    rr = extract_pred_row(pred$cumRRfit)[cum_idx],
    rr_lower = extract_pred_row(pred$cumRRlow)[cum_idx],
    rr_upper = extract_pred_row(pred$cumRRhigh)[cum_idx],
    coef = extract_pred_row(pred$cumfit)[cum_idx],
    se = extract_pred_row(pred$cumse)[cum_idx]
  )
  
  single_res[, `:=`(
    n_complete = diag_stats$n_complete,
    n_exposed_days = diag_stats$n_exposed_days,
    n_exposed_strata = diag_stats$n_exposed_strata,
    n_exposed_positive_rows = diag_stats$n_exposed_positive_rows,
    exposed_outcome_sum = diag_stats$exposed_outcome_sum,
    family_used = model_fit$family_used,
    theta_est = model_fit$theta_est,
    theta_converged = model_fit$theta_converged,
    theta_iterations = model_fit$theta_iterations,
    error_msg = NA_character_
  )]
  
  cumulative_res[, `:=`(
    n_complete = diag_stats$n_complete,
    n_exposed_days = diag_stats$n_exposed_days,
    n_exposed_strata = diag_stats$n_exposed_strata,
    n_exposed_positive_rows = diag_stats$n_exposed_positive_rows,
    exposed_outcome_sum = diag_stats$exposed_outcome_sum,
    family_used = model_fit$family_used,
    theta_est = model_fit$theta_est,
    theta_converged = model_fit$theta_converged,
    theta_iterations = model_fit$theta_iterations,
    error_msg = NA_character_
  )]
  
  list(single = single_res, cumulative = cumulative_res, error = NA_character_)
}

# Run one analysis job (one row of the job table): optionally restrict to a region
# or climate subgroup, fit the model, write its RR tables, and return a one-row
# status record.
run_one_job <- function(job, data,
                        min_exposed_days,
                        min_exposed_strata,
                        min_exposed_positive_rows,
                        min_exposed_outcome_sum,
                        lag_df) {
  job <- as.list(job)
  dt_job <- copy(data)
  
  if (!is.na(job$subgroup_type) && !is.na(job$subgroup_name)) {
    if (identical(job$subgroup_type, "region")) {
      dt_job <- dt_job[region == as.character(job$subgroup_name)]
    }
    if (identical(job$subgroup_type, "climate")) {
      dt_job <- dt_job[climate_code == as.integer(job$subgroup_name)]
    }
  }
  
  status <- data.table(
    analysis_scope = job$analysis_scope,
    subgroup_type = if (is.na(job$subgroup_type)) NA_character_ else as.character(job$subgroup_type),
    subgroup_name = if (is.na(job$subgroup_name)) NA_character_ else as.character(job$subgroup_name),
    outcome = as.character(job$outcome),
    variable = as.character(job$variable),
    pre_control = if (isTRUE(job$include_pre)) "with_pre" else "no_pre",
    lag_window = as.integer(job$lag_window),
    single_path = as.character(job$single_path),
    cum_path = as.character(job$cum_path),
    status = NA_character_,
    runtime_sec = NA_real_,
    n_complete = NA_integer_,
    n_exposed_days = NA_integer_,
    family_used = NA_character_,
    theta_est = NA_real_,
    error_msg = NA_character_
  )
  
  t0 <- Sys.time()
  
  if (nrow(dt_job) == 0L) {
    status[, `:=`(status = "empty_subset", runtime_sec = 0, error_msg = "No rows in subgroup subset.")]
    return(status)
  }
  
  res <- fit_dlm_crossbasis(
    dt = dt_job,
    outcome = job$outcome,
    variable = job$variable,
    include_pre = isTRUE(job$include_pre),
    max_lag = as.integer(job$lag_window),
    lag_df = min(lag_df, as.integer(job$lag_window)),
    min_exposed_days = min_exposed_days,
    min_exposed_strata = min_exposed_strata,
    min_exposed_positive_rows = min_exposed_positive_rows,
    min_exposed_outcome_sum = min_exposed_outcome_sum
  )
  
  if (!is.null(res$single) && !is.null(res$cumulative)) {
    res$single[, `:=`(
      analysis_scope = job$analysis_scope,
      subgroup_type = if (is.na(job$subgroup_type)) NA_character_ else job$subgroup_type,
      subgroup_name = if (is.na(job$subgroup_name)) NA_character_ else job$subgroup_name,
      outcome = job$outcome,
      variable = job$variable,
      exposure_family = job$exposure_family,
      definition_type = job$definition_type,
      exposure_label = job$exposure_label,
      pre_control = if (isTRUE(job$include_pre)) "with_pre" else "no_pre",
      lag_window = job$lag_window,
      output_file = basename(job$single_path)
    )]
    
    res$cumulative[, `:=`(
      analysis_scope = job$analysis_scope,
      subgroup_type = if (is.na(job$subgroup_type)) NA_character_ else job$subgroup_type,
      subgroup_name = if (is.na(job$subgroup_name)) NA_character_ else job$subgroup_name,
      outcome = job$outcome,
      variable = job$variable,
      exposure_family = job$exposure_family,
      definition_type = job$definition_type,
      exposure_label = job$exposure_label,
      pre_control = if (isTRUE(job$include_pre)) "with_pre" else "no_pre",
      lag_window = job$lag_window,
      output_file = basename(job$cum_path)
    )]
    
    dir.create(dirname(job$single_path), recursive = TRUE, showWarnings = FALSE)
    write_csv_utf8(res$single, job$single_path)
    write_csv_utf8(res$cumulative, job$cum_path)
    
    status[, `:=`(
      status = "success",
      runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
      n_complete = res$single$n_complete[1L],
      n_exposed_days = res$single$n_exposed_days[1L],
      family_used = res$single$family_used[1L],
      theta_est = res$single$theta_est[1L],
      error_msg = NA_character_
    )]
  } else {
    diag_info <- res$diag
    status[, `:=`(
      status = "fit_failed",
      runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
      n_complete = if (!is.null(diag_info)) diag_info$n_complete else NA_integer_,
      n_exposed_days = if (!is.null(diag_info)) diag_info$n_exposed_days else NA_integer_,
      error_msg = ifelse(length(res$error) == 0L, "Unknown error.", as.character(res$error))
    )]
  }
  
  status
}

# Run a whole table of jobs -- in parallel across workers when a cluster is
# supplied, otherwise sequentially.
run_job_table <- function(job_table, data,
                          min_exposed_days,
                          min_exposed_strata,
                          min_exposed_positive_rows,
                          min_exposed_outcome_sum,
                          lag_df,
                          cl = NULL) {
  if (nrow(job_table) == 0L) return(data.table())
  
  if (!is.null(cl) && nrow(job_table) > 1L) {
    results <- foreach(
      i = seq_len(nrow(job_table)),
      .inorder = FALSE,
      .packages = c("data.table", "foreach", "doParallel", "gnm", "MASS", "dlnm", "splines", "parallel")
    ) %dopar% {
      run_one_job(
        job = job_table[i],
        data = data,
        min_exposed_days = min_exposed_days,
        min_exposed_strata = min_exposed_strata,
        min_exposed_positive_rows = min_exposed_positive_rows,
        min_exposed_outcome_sum = min_exposed_outcome_sum,
        lag_df = lag_df
      )
    }
  } else {
    results <- lapply(seq_len(nrow(job_table)), function(i) {
      run_one_job(
        job = job_table[i],
        data = data,
        min_exposed_days = min_exposed_days,
        min_exposed_strata = min_exposed_strata,
        min_exposed_positive_rows = min_exposed_positive_rows,
        min_exposed_outcome_sum = min_exposed_outcome_sum,
        lag_df = lag_df
      )
    })
  }
  
  rbindlist(results, fill = TRUE)
}

if (!run_mode %in% c("full", "middle", "pilot", "smoke")) {
  stop("run_mode must be smoke, pilot, middle, or full.")
}

requested_cores <- as.integer(requested_cores)
if (length(requested_cores) != 1L || is.na(requested_cores) || requested_cores < 1L) {
  stop("requested_cores must be an integer >= 1.")
}

pilot_mode <- identical(run_mode, "pilot")
middle_mode <- identical(run_mode, "middle")
full_mode <- identical(run_mode, "full")

# ===========================  ANALYSIS PARAMETERS  ===========================
outcomes <- c("case", "Pre_school", "School", "Adults", "Elders")
# Eight flood-exposure definitions used by the main and sensitivity analyses:
# severity (general-only / severe / any), two hydrologic alternatives (high soil
# moisture, high runoff) for general and severe floods, and a >=50 mm rainstorm proxy.
exposure_dict <- data.table(
  variable = c("GF_general_only", "SF", "GF", "GF_1", "GF_2", "SF_1", "SF_2", "GF_3"),
  file_suffix = c("GFgeneralOnly", "SF", "GF", "GF1", "GF2", "SF1", "SF2", "GF3"),
  exposure_family = c(
    "severity_exclusive", "severity_nested", "severity_nested",
    "hydrologic_general", "hydrologic_general",
    "hydrologic_severe", "hydrologic_severe",
    "rainfall_proxy"
  ),
  definition_type = c(
    "exclusive_general_flood", "severe_flood", "nested_any_flood",
    "alternative_swvl1", "alternative_runoff",
    "alternative_swvl1_severe", "alternative_runoff_severe",
    "rainstorm_50p"
  ),
  exposure_label = c(
    "general flood (excluding severe)", "severe flood", "general flood or worse",
    "general flood + high soil moisture", "general flood + high runoff",
    "severe flood + high soil moisture", "severe flood + high runoff",
    "rainstorm or above (daily precip >=50mm)"
  )
)
selected_exposures <- exposure_dict$variable
pilot_exposures <- c("GF_general_only", "SF")
max_lag_main <- 28L
lag_windows <- c(7L, 14L, 21L, 28L)
lag_df <- 5L
# Minimum exposed days / strata / positive rows / case sum required before a model
# is attempted -- skips ill-identified outcome x exposure x subgroup cells.
min_exposed_days <- 10L
min_exposed_strata <- 10L
min_exposed_positive_rows <- 5L
min_exposed_outcome_sum <- 10L

if (pilot_mode) {
  selected_exposures <- pilot_exposures
  outcomes <- "case"
}

# ==========================  OUTPUT PATHS & RUN ID  ==========================
project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
setwd(project_root)
if (!file.exists(file.path(project_root, input_file))) {
  stop("Input file not found: ", file.path(project_root, input_file), ". Please check project_root and input_file settings.")
}
output_root_abs <- normalizePath(file.path(project_root, output_root), winslash = "/", mustWork = FALSE)
dir.create(output_root_abs, recursive = TRUE, showWarnings = FALSE)

run_id <- paste(
  "Flood_Dysentery_main_effect",
  format(Sys.time(), "%y-%m-%d-%H-%M-%S"),
  gsub("[^A-Za-z0-9_-]+", "-", trimws(run_label)),
  sep = "_"
)
output_dir <- file.path(output_root_abs, run_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dir_main <- file.path(output_dir, "main")
dir_nopre <- file.path(output_dir, "nopre")
dir_lag <- file.path(output_dir, "lag")
dir_region <- file.path(output_dir, "region")
dir_climate <- file.path(output_dir, "climate")
invisible(lapply(c(dir_main, dir_nopre, dir_lag, dir_region, dir_climate), dir.create, recursive = TRUE, showWarnings = FALSE))

cat("==========  RUN CONFIGURATION  ==========\n")
cat("Project root: ", project_root, "\n", sep = "")
cat("Input file: ", input_file, "\n", sep = "")
cat("Run mode: ", run_mode, "\n", sep = "")
cat("Cores: ", requested_cores, "\n", sep = "")
cat("Run label: ", run_label, "\n", sep = "")
cat("Output dir: ", output_dir, "\n", sep = "")
cat("Exposures: ", paste(selected_exposures, collapse = ", "), "\n", sep = "")
cat(
  "Analysis flow: ",
  if (identical(run_mode, "smoke")) {
    "smoke test only"
  } else if (pilot_mode) {
    "pilot (main analysis, reduced outcomes/exposures)"
  } else if (middle_mode) {
    "middle (all outcomes/exposures + no-precipitation + lag-window sensitivity)"
  } else {
    "full (middle + region and climate-zone subgroup analyses)"
  },
  "\n",
  sep = ""
)
cat("===========================================\n")

# ======================  DATA LOADING & PREPROCESSING  =======================
data <- fread(file.path(project_root, input_file))
setDT(data)

required_cols <- c("code", "date", "GF", "SF", "tem", "RH", "holiday", "region", "climate_code")
missing_required <- setdiff(required_cols, names(data))
if (length(missing_required) > 0L) {
  stop("Missing required columns: ", paste(missing_required, collapse = ", "))
}

data[, date := as.Date(date)]
data[, region := as.character(region)]
numeric_cols <- intersect(c("case", "Pre_school", "School", "Adults", "Elders", "GF", "SF", "ro", "swvl1", "tem", "RH", "pre", "holiday", "climate_code"), names(data))
data[, (numeric_cols) := lapply(.SD, to_numeric_safe), .SDcols = numeric_cols]
data[, climate_code := as.integer(climate_code)]
setorderv(data, c("code", "date"))

# General flood excluding severe (the severity-exclusive definition).
if (!"GF_general_only" %in% names(data)) {
  data[, GF_general_only := as.integer(GF == 1L & SF == 0L)]
}

# Hydrologic alternatives: flood days that also show high soil moisture
# (swvl1 > city-specific p75) or high runoff (ro > city-specific p95).
if (all(c("ro", "swvl1") %in% names(data))) {
  data[, swvl1_thresh75 := safe_quantile(swvl1, 0.75), by = code]
  data[, ro_thresh95 := safe_quantile(ro, 0.95), by = code]
  data[, GF_1 := as.integer(GF == 1L & swvl1 > swvl1_thresh75)]
  data[, SF_1 := as.integer(SF == 1L & swvl1 > swvl1_thresh75)]
  data[, GF_2 := as.integer(GF == 1L & ro > ro_thresh95)]
  data[, SF_2 := as.integer(SF == 1L & ro > ro_thresh95)]
  data[, c("swvl1_thresh75", "ro_thresh95") := NULL]
} else if (!all(c("GF_1", "SF_1", "GF_2", "SF_2") %in% names(data))) {
  stop("Missing ro/swvl1, and GF_1 / SF_1 / GF_2 / SF_2 not found in data.")
}

# Rainfall proxy: a rainstorm day is daily precipitation >= 50 mm.
if ("pre" %in% names(data)) {
  data[, GF_3 := as.integer(pre >= 50)]
} else if (!"GF_3" %in% names(data)) {
  stop("Missing pre, and GF_3 not found in data.")
}

data[, `:=`(
  code = factor(code),
  year = factor(format(date, "%Y")),
  month = factor(format(date, "%m")),
  dow = factor(format(date, "%u"))
)]
# Case-crossover stratum = city x year x month x day-of-week; the conditional model
# contrasts flood vs non-flood days within each stratum, removing all within-stratum
# confounding (season, long-term trend, day-of-week).
data[, stratum := factor(paste(code, year, month, dow, sep = ":"))]
# 0-2 day moving means of temperature, humidity and precipitation (confounder
# controls entered as splines in the DLNM).
data[, temlag02 := runMean(tem, 0:2, group = code)]
data[, rhlag02 := runMean(RH, 0:2, group = code)]
data[, prelag02 := if ("pre" %in% names(data)) runMean(pre, 0:2, group = code) else NA_real_]

exposure_summary <- rbindlist(lapply(exposure_dict$variable, function(v) {
  meta <- exposure_dict[variable == v]
  data[, .(
    variable = v,
    exposure_family = meta$exposure_family[1L],
    definition_type = meta$definition_type[1L],
    exposure_label = meta$exposure_label[1L],
    exposed_days = sum(get(v) == 1L, na.rm = TRUE),
    exposed_codes = uniqueN(code[get(v) == 1L]),
    exposed_strata = uniqueN(stratum[get(v) == 1L]),
    exposed_positive_rows = sum(get(v) == 1L & case > 0, na.rm = TRUE),
    exposed_case_sum = sum(case[get(v) == 1L], na.rm = TRUE)
  )]
}), fill = TRUE)
write_csv_utf8(exposure_summary, file.path(output_dir, "exposure_summary.csv"))

if (identical(run_mode, "smoke")) {
  cat("\n========== Smoke test completed ==========\n")
  cat("Packages loaded, data read and preprocessed. Model fitting not started.\n")
  quit(save = "no", status = 0L)
}

# ============================  PARALLEL BACKEND  =============================
num_cores <- max(1L, requested_cores)
cl <- NULL
if (num_cores > 1L) {
  cl <- prepare_cluster_simple(num_cores)
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  clusterExport(
    cl,
    varlist = c(
      "data", "write_csv_utf8", "estimate_theta_ml", "subset_crossbasis_rows",
      "resolve_cumulative_windows", "extract_pred_row", "fit_nb_gnm_model",
      "fit_dlm_crossbasis", "run_one_job", "lag_df",
      "min_exposed_days", "min_exposed_strata", "min_exposed_positive_rows",
      "min_exposed_outcome_sum"
    ),
    envir = environment()
  )
}

# ============================  ANALYSIS DRIVERS  =============================
# Build a job table = outcomes x exposures (x subgroup levels), attaching each
# exposure's metadata and the output file paths.
build_jobs <- function(outcomes, variables, analysis_scope, include_pre, lag_window,
                       subgroup_type = NA_character_, subgroup_values = NA_character_,
                       path_builder) {
  if (length(subgroup_values) == 1L && is.na(subgroup_values)) {
    jobs <- CJ(outcome = outcomes, variable = variables, unique = TRUE)
    jobs[, subgroup_name := NA_character_]
    subgroup_type_value <- NA_character_
  } else {
    jobs <- CJ(outcome = outcomes, variable = variables, subgroup_name = as.character(subgroup_values), unique = TRUE)
    subgroup_type_value <- subgroup_type
  }
  
  jobs <- merge(jobs, exposure_dict, by = "variable", all.x = TRUE, sort = FALSE)
  jobs[, `:=`(
    analysis_scope = analysis_scope,
    subgroup_type = subgroup_type_value,
    include_pre = include_pre,
    lag_window = as.integer(lag_window)
  )]
  
  paths <- path_builder(jobs)
  jobs[, `:=`(single_path = paths$single_path, cum_path = paths$cum_path)]
  jobs[]
}

# Run and print a labelled block of jobs.
run_block <- function(title, jobs) {
  cat("\n========== ", title, " ==========\n", sep = "")
  if (nrow(jobs) == 0L) {
    return(data.table())
  }
  
  run_job_table(
    job_table = jobs,
    data = data,
    min_exposed_days = min_exposed_days,
    min_exposed_strata = min_exposed_strata,
    min_exposed_positive_rows = min_exposed_positive_rows,
    min_exposed_outcome_sum = min_exposed_outcome_sum,
    lag_df = lag_df,
    cl = cl
  )
}

# ==============================  MAIN ANALYSIS  ==============================
main_jobs <- build_jobs(
  outcomes = outcomes,
  variables = selected_exposures,
  analysis_scope = "main_pre",
  include_pre = TRUE,
  lag_window = max_lag_main,
  path_builder = function(jobs) {
    jobs[, .(
      single_path = file.path(dir_main, paste0(outcome, "_", file_suffix, "_single.csv")),
      cum_path = file.path(dir_main, paste0(outcome, "_", file_suffix, "_cumulative.csv"))
    )]
  }
)

all_status <- list(run_block("Main analysis", main_jobs))

# ==========================  SENSITIVITY ANALYSES  ===========================
if (middle_mode || full_mode) {
  nopre_jobs <- build_jobs(
    outcomes = outcomes,
    variables = selected_exposures,
    analysis_scope = "sensitivity_nopre",
    include_pre = FALSE,
    lag_window = max_lag_main,
    path_builder = function(jobs) {
      jobs[, .(
        single_path = file.path(dir_nopre, paste0(outcome, "_", file_suffix, "_single_nopre.csv")),
        cum_path = file.path(dir_nopre, paste0(outcome, "_", file_suffix, "_cumulative_nopre.csv"))
      )]
    }
  )
  all_status[[length(all_status) + 1L]] <- run_block("Sensitivity analysis without precipitation control", nopre_jobs)
  
  lag_status_all <- list()
  for (lag_w in lag_windows) {
    dir_lag_w <- file.path(dir_lag, paste0("lag", lag_w))
    dir.create(dir_lag_w, recursive = TRUE, showWarnings = FALSE)
    lag_jobs <- build_jobs(
      outcomes = outcomes,
      variables = selected_exposures,
      analysis_scope = "sensitivity_lag",
      include_pre = TRUE,
      lag_window = lag_w,
      path_builder = function(jobs) {
        jobs[, .(
          single_path = file.path(dir_lag_w, paste0(outcome, "_", file_suffix, "_lag", lag_w, "_single.csv")),
          cum_path = file.path(dir_lag_w, paste0(outcome, "_", file_suffix, "_lag", lag_w, "_cumulative.csv"))
        )]
      }
    )
    lag_status_all[[length(lag_status_all) + 1L]] <- run_block(paste0("Lag window sensitivity analysis (lag", lag_w, ")"), lag_jobs)
  }
  all_status[[length(all_status) + 1L]] <- rbindlist(lag_status_all, fill = TRUE)
}

# Region- and climate-zone subgroup analyses (full tier only)
if (full_mode) {
  region_levels <- sort(unique(na.omit(as.character(data$region))))
  region_jobs <- build_jobs(
    outcomes = outcomes,
    variables = selected_exposures,
    analysis_scope = "subgroup_region",
    include_pre = TRUE,
    lag_window = max_lag_main,
    subgroup_type = "region",
    subgroup_values = region_levels,
    path_builder = function(jobs) {
      jobs[, prefix := gsub("\\s+", "_", subgroup_name)]
      jobs[, .(
        single_path = file.path(dir_region, paste0(prefix, "_", outcome, "_", file_suffix, "_single.csv")),
        cum_path = file.path(dir_region, paste0(prefix, "_", outcome, "_", file_suffix, "_cumulative.csv"))
      )]
    }
  )
  all_status[[length(all_status) + 1L]] <- run_block("Region subgroup analysis", region_jobs)
  
  climate_levels <- sort(unique(na.omit(as.integer(data$climate_code))))
  climate_jobs <- build_jobs(
    outcomes = outcomes,
    variables = selected_exposures,
    analysis_scope = "subgroup_climate",
    include_pre = TRUE,
    lag_window = max_lag_main,
    subgroup_type = "climate",
    subgroup_values = climate_levels,
    path_builder = function(jobs) {
      jobs[, .(
        single_path = file.path(dir_climate, paste0("climate", subgroup_name, "_", outcome, "_", file_suffix, "_single.csv")),
        cum_path = file.path(dir_climate, paste0("climate", subgroup_name, "_", outcome, "_", file_suffix, "_cumulative.csv"))
      )]
    }
  )
  all_status[[length(all_status) + 1L]] <- run_block("Climate zone subgroup analysis", climate_jobs)
}

# =======================  RESULTS: AGGREGATE & EXPORT  =======================
all_status <- Filter(function(x) nrow(x) > 0L, all_status)

task_status <- if (length(all_status) > 0L) rbindlist(all_status, fill = TRUE) else data.table()
task_status_export <- task_status[, .(
  analysis_scope, subgroup_type, subgroup_name, outcome, variable,
  pre_control, lag_window, status, runtime_sec, n_complete,
  n_exposed_days, family_used, theta_est, error_msg
)]
write_csv_utf8(task_status_export, file.path(output_dir, "task_status.csv"))

success_status <- task_status[status == "success"]
if (nrow(success_status) > 0L) {
  # read theta_est etc. as double (a near-Poisson cell can produce a very large
  # dispersion that fread would otherwise type as integer64, mismatching other files)
  read_result <- function(f) fread(f, integer64 = "double")
  all_single <- rbindlist(lapply(success_status$single_path[file.exists(success_status$single_path)], read_result), fill = TRUE)
  all_cumulative <- rbindlist(lapply(success_status$cum_path[file.exists(success_status$cum_path)], read_result), fill = TRUE)
  
  write_csv_utf8(all_single, file.path(output_dir, "all_results_single.csv"))
  write_csv_utf8(all_cumulative, file.path(output_dir, "all_results_cumulative.csv"))
}

cat("\n========== Analysis completed ==========\n")
cat("Results directory: ", output_dir, "\n", sep = "")
cat("Status summary:\n")
if (nrow(task_status) > 0L) {
  print(task_status[, .N, by = status])
} else {
  cat("No task results generated.\n")
}



