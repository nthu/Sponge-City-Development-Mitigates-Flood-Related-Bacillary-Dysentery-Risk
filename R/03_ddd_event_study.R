# =============================================================================
# 03_ddd_event_study.R
#
# Difference-in-difference-in-differences (DDD) event study of the Sponge City
# Development policy effect on the flood -> bacillary-dysentery (BD) association
# (manuscript equation 2). A distributed-lag non-linear model (DLNM) cross-basis
# for flood exposure is interacted with treated-city event time:
#
#   log(mu_ict) = a_i + g_y + d_m + sum_l beta_l Flood_{t-l}
#               + sum_tau sum_l gamma_{tau,l} 1[treated & RelYear=tau] Flood_{t-l}
#               + ns(weather) + X_it
#
# RelYear is missing for controls, so it is excluded from complete-case
# deletion and each event-time term cb_lag * 1[treated & RelYear=tau] is
# 0 for controls. Controls thus identify the baseline flood effect.
# The standardised DDD RR for year tau is exp(gamma_{tau,.}); tau = -1 is the
# reference. Fit as a conditional negative-binomial model with city, year and
# month fixed effects and city-clustered SEs, for total BD cases and the four
# age groups, plus the parallel-trend and placebo tests.
# =============================================================================

message("== DDD event study ==")

dt <- fread(DDD_PANEL, encoding = "UTF-8")
dt[, date := as.Date(date)]
setorder(dt, code, date)   # order rows by city and date for the cross-basis lag computation

dt[, flood   := as.integer(GF == 1L)]
dt[, treated := as.integer(sponge == 1L)]
dt[, start_date := fifelse(policy_year == 2015L, START_2015,
                    fifelse(policy_year == 2016L, START_2016, as.Date(NA)))]
# relative year is defined for treated cities; missing for controls
dt[, rel_year_exact := fifelse(treated == 1L, floor(as.numeric(date - start_date) / 365.25), NA_real_)]
dt[, `:=`(year_fe = factor(year(date)), month_fe = factor(month(date)))]
dt[, `:=`(temlag02 = calc_lag_mean(tem, 2L, code),
          rhlag02  = calc_lag_mean(RH,  2L, code),
          prelag02 = calc_lag_mean(pre, 2L, code))]
# annual mean temperature covariate, derived from daily temperature
dt[, tem_annual := mean(tem, na.rm = TRUE), by = .(code, year(date))]

# Fit the DDD event-study model for one outcome. Returns single-lag and cumulative RR tables
# (RR = exp(gamma_tau); reference year -1 = 1 at every lag).
fit_ddd <- function(data, outcome, saturated = FALSE) {
  d <- copy(data)
  window <- (d$treated == 0L) | (d$treated == 1L & !is.na(d$rel_year_exact) &
              d$rel_year_exact >= min(REL_WINDOW) & d$rel_year_exact <= max(REL_WINDOW))
  d <- d[window]
  rel_levels <- sort(unique(d[treated == 1L, rel_year_exact])); rel_levels <- rel_levels[!is.na(rel_levels)]

  cb <- dlnm::crossbasis(d$flood, lag = MAX_LAG, argvar = list(fun = "lin"),
                         arglag = list(fun = "ns", df = LAG_DF), group = d$code)
  cb_cols <- paste0("cb_", seq_len(ncol(cb)))
  cb_mat  <- as.data.table(cb); setnames(cb_mat, cb_cols)
  d <- cbind(d, cb_mat); keep <- complete.cases(cb_mat); d <- d[keep]
  cb <- subset_crossbasis(cb, keep); colnames(cb) <- cb_cols

  needed <- c(outcome, cb_cols, "temlag02", "rhlag02", "prelag02", DDD_COVARIATES,
              "treated", "rel_year_exact", "code", "year_fe", "month_fe")
  # RelYear is missing for controls, so it is excluded from complete-case deletion
  cc_cols <- setdiff(needed, "rel_year_exact")
  keep2 <- complete.cases(d[, ..cc_cols]); d <- d[keep2]; cb <- subset_crossbasis(cb, keep2); colnames(cb) <- cb_cols
  md <- copy(d[, ..needed]); md[, treated := as.integer(treated)]

  rel_label <- function(x) ifelse(x < 0L, paste0("m", abs(x)), paste0("p", x))
  triple_map <- list(); triple_all <- character(0); d_main <- character(0)
  for (r in setdiff(rel_levels, REF_REL_YEAR)) {
    lab <- rel_label(r); dcol <- paste0("D_", lab)
    md[, (dcol) := as.integer(treated == 1L & !is.na(rel_year_exact) & rel_year_exact == r)]
    d_main <- c(d_main, dcol)
    cols_r <- character(0)
    for (cv in cb_cols) { tc <- paste0(cv, "_D_", lab); md[, (tc) := get(cv) * get(dcol)]
      triple_all <- c(triple_all, tc); cols_r <- c(cols_r, tc) }
    triple_map[[as.character(r)]] <- cols_r
  }
  controls <- c("ns(temlag02, 3)", "ns(rhlag02, 2)", "ns(prelag02, 2)", DDD_COVARIATES)
  # saturated = TRUE adds the bare Treated x 1[RelYear=tau] main effects (D_tau), so
  # the flood interaction is identified net of treated-specific event-time level shifts.
  rhs <- c(cb_cols, triple_all, if (saturated) d_main, controls)
  fm <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + "),
                         "| code + year_fe + month_fe"))
  model <- fixest::fenegbin(fm, data = md, cluster = ~code, fixef.rm = "perfect")

  ca <- coef(model); va <- vcov(model)
  ref_s <- data.table(rel_year = REF_REL_YEAR, lag = 0:MAX_LAG, rr = 1, lower = 1, upper = 1, term = "reference")
  ref_c <- data.table(rel_year = REF_REL_YEAR, lag = 0:MAX_LAG, cum_rr = 1, cum_lower = 1, cum_upper = 1, term = "reference")
  eff <- list(ref_s); cum <- list(ref_c)
  for (r in setdiff(rel_levels, REF_REL_YEAR)) {
    cols <- triple_map[[as.character(r)]]; if (!all(cols %in% names(ca))) next
    g <- ca[cols]; Vg <- va[cols, cols, drop = FALSE]
    names(g) <- cb_cols; rownames(Vg) <- cb_cols; colnames(Vg) <- cb_cols
    pr <- dlnm::crosspred(cb, coef = g, vcov = Vg, at = 1, cen = 0, cumul = TRUE, model.link = "log")
    eff[[as.character(r)]] <- data.table(rel_year = as.numeric(r), lag = 0:MAX_LAG,
      rr = pred_row(pr$matRRfit), lower = pred_row(pr$matRRlow), upper = pred_row(pr$matRRhigh), term = "ddd")
    cum[[as.character(r)]] <- data.table(rel_year = as.numeric(r), lag = 0:MAX_LAG,
      cum_rr = pred_row(pr$cumRRfit), cum_lower = pred_row(pr$cumRRlow), cum_upper = pred_row(pr$cumRRhigh), term = "ddd")
  }
  list(model = model, effects = rbindlist(eff, fill = TRUE)[order(rel_year, lag)],
       cumulative = rbindlist(cum, fill = TRUE)[order(rel_year, lag)], triple_map = triple_map)
}

out_dir <- file.path(RESULTS_DIR, "ddd_event_study")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
results <- list()
for (oc in DDD_OUTCOMES) {
  message("  fitting outcome: ", oc)
  res <- tryCatch(fit_ddd(dt, oc), error = function(e) { warning(oc, ": ", conditionMessage(e)); NULL })
  if (is.null(res)) next
  results[[oc]] <- res
  fwrite(res$effects, file.path(out_dir, paste0(oc, "_single_lag_RR.csv")))
  fwrite(res$cumulative, file.path(out_dir, paste0(oc, "_cumulative_RR.csv")))
  cat(sprintf("    %s year-3 lag0 RR = %.4f\n", oc, res$effects[rel_year == 3 & lag == 0, rr]))
}

main_res <- results[["case"]]
if (!is.null(main_res)) {
  pre_terms <- unlist(main_res$triple_map[intersect(c("-4","-3","-2"), names(main_res$triple_map))], use.names = FALSE)
  pt <- wald_chisq(main_res$model, pre_terms)
  cat("  parallel-trend Wald test (yrs -4,-3,-2):  P =", round(pt$p_value, 4), "\n")
  fwrite(pt, file.path(out_dir, "parallel_trend_test.csv"))

  # placebo: shift treated policy start two years earlier
  dtp <- copy(dt)
  dtp[treated == 1L, start_date := start_date %m-% years(2)]
  dtp[, rel_year_exact := fifelse(treated == 1L, floor(as.numeric(date - start_date) / 365.25), NA_real_)]
  pres <- tryCatch(fit_ddd(dtp, "case"), error = function(e) NULL)
  if (!is.null(pres)) {
    pl_terms <- unlist(pres$triple_map[intersect(c("0","1"), names(pres$triple_map))], use.names = FALSE)
    pl <- wald_chisq(pres$model, pl_terms)
    cat("  placebo (lead-2y) Wald test (fake yrs 0,1):  P =", round(pl$p_value, 4), "\n")
    fwrite(pl, file.path(out_dir, "placebo_lead2_test.csv"))
    fwrite(pres$effects, file.path(out_dir, "placebo_lead2_single_lag_RR.csv"))
    fwrite(pres$cumulative, file.path(out_dir, "placebo_lead2_cumulative_RR.csv"))
  }

  # saturated-model sensitivity: refit with the bare Treated x 1[RelYear=tau] main
  # effects added, so the flood interaction is identified net of any treated-specific
  # event-time level shift. Compare the lag-0 RR by event year to the primary model.
  sat <- tryCatch(fit_ddd(dt, "case", saturated = TRUE),
                  error = function(e) { warning("saturated model: ", conditionMessage(e)); NULL })
  if (!is.null(sat)) {
    cmp <- merge(main_res$effects[term == "ddd" & lag == 0L, .(rel_year, rr_primary = rr)],
                 sat$effects[term == "ddd" & lag == 0L, .(rel_year, rr_saturated = rr)], by = "rel_year")
    cat("  saturated-model sensitivity (lag0 RR by event year):\n")
    print(cmp[, .(rel_year, rr_primary = round(rr_primary, 4), rr_saturated = round(rr_saturated, 4))])
    fwrite(cmp, file.path(out_dir, "saturated_sensitivity_lag0.csv"))
  }
}
# ---- effect modification by seven pre-policy city characteristics ---------------
# Median-split the matched cities by each modifier's pre-policy (PRE_YEARS) mean,
# then within each half fit a Poisson event study (fepois; lag spline df = 2;
# city and month fixed effects; covariates gdp/water/str/den/doc) and report the
# post-policy (years 0-4) average lag-0 RR (geometric mean; CI from the across-year
# variability of the yearly log-RRs).
# muffle benign fixest/fepois notes on sparse capacity-subgroups (results unaffected;
# fit errors are still caught per subgroup)
hetero_quiet <- function(expr) withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"))
run_hetero_ddd <- function(data, lag_df = 2L) {
  d <- copy(data); d[, month := factor(month(date))]
  d[treated == 0L, rel_year_exact := REF_REL_YEAR]
  d <- d[treated == 0L | (treated == 1L & rel_year_exact >= min(REL_WINDOW) & rel_year_exact <= max(REL_WINDOW))]
  rl <- sort(unique(d[treated == 1L, rel_year_exact])); if (!length(rl)) return(NULL)
  d[, rel_year_factor := relevel(factor(rel_year_exact, levels = rl), ref = as.character(REF_REL_YEAR))]
  cb <- dlnm::crossbasis(d$flood, lag = MAX_LAG, argvar = list(fun = "lin"), arglag = list(fun = "ns", df = lag_df), group = d$code)
  cb_cols <- paste0("cb_", seq_len(ncol(cb))); cm <- as.data.table(cb); setnames(cm, cb_cols)
  d <- cbind(d, cm); keep <- complete.cases(cm); d <- d[keep]; cb <- subset_crossbasis(cb, keep); colnames(cb) <- cb_cols
  needed <- c("case", cb_cols, "temlag02", "rhlag02", "prelag02", DDD_COVARIATES, "treated", "rel_year_factor", "code", "month")
  md <- na.omit(d[, ..needed])
  it <- vapply(cb_cols, function(v) sprintf("i(rel_year_factor, var = %s, ref = '-1')", v), character(1))
  ctrl <- c("ns(temlag02, 2)", "ns(rhlag02, 2)", "ns(prelag02, 2)", DDD_COVARIATES)
  fm <- as.formula(paste("case ~", paste(c(cb_cols, "i(rel_year_factor, ref = '-1')", it, ctrl), collapse = " + "), "| code + month"))
  model <- hetero_quiet(fixest::fepois(fm, data = md, cluster = ~code, fixef.rm = "perfect"))
  ca <- coef(model); va <- hetero_quiet(vcov(model)); cn <- names(ca); mpos <- match(cb_cols, cn); if (anyNA(mpos)) return(NULL)
  cbm <- ca[mpos]; Vm <- va[mpos, mpos]
  base0 <- pred_row(dlnm::crosspred(cb, coef = cbm, vcov = Vm, at = 1, cen = 0, model.link = "log")$matRRfit)[1]
  rels <- setdiff(unique(gsub("rel_year_factor::(.*):cb_.*", "\\1", grep("rel_year_factor::.*:cb_", cn, value = TRUE))), "-1")
  out <- list()
  for (r in rels) {
    pat <- paste0("rel_year_factor::", r, ":cb_"); idx <- grep(pat, cn); if (length(idx) != length(cb_cols)) next
    idx <- idx[order(as.numeric(gsub(pat, "", cn[idx])))]; tot <- cbm + ca[idx]
    Vt <- Vm + va[idx, idx] + va[mpos, idx] + t(va[mpos, idx])
    rr0 <- pred_row(dlnm::crosspred(cb, coef = tot, vcov = Vt, at = 1, cen = 0, model.link = "log")$matRRfit)[1]
    out[[r]] <- data.table(rel_year = as.numeric(r), rr = rr0 / base0)
  }
  rbindlist(out)
}
if (!"SE" %in% names(dt)) dt[, SE := st / expenditure]
if (!"Stru" %in% names(dt)) dt[, Stru := `3_gdp` / `2_gdp`]
hetero_vars <- c("invest", "pt", "bed", "wspd", "dpd", "SE", "Stru")
hetero_vars <- hetero_vars[hetero_vars %in% names(dt)]
em_hetero <- rbindlist(lapply(hetero_vars, function(v) {
  pm <- dt[year(date) %in% PRE_YEARS & is.finite(get(v)), .(pv = mean(get(v), na.rm = TRUE)), by = code]
  med <- median(pm$pv, na.rm = TRUE)
  rbindlist(lapply(c("high", "low"), function(grp) {
    codes <- pm$code[if (grp == "high") pm$pv > med else pm$pv <= med]
    res <- tryCatch(run_hetero_ddd(dt[code %in% codes]), error = function(e) NULL); if (is.null(res)) return(NULL)
    y <- res[rel_year >= 0, log(rr)]; y <- y[is.finite(y)]; if (length(y) < 2L) return(NULL)
    m <- mean(y); s <- stats::sd(y) / sqrt(length(y))
    data.table(variable = v, group = grp, n_cities = length(codes),
               rr = round(exp(m), 4), lower = round(exp(m - 1.96 * s), 4), upper = round(exp(m + 1.96 * s), 4), n_years = length(y))
  }))
}))
if (nrow(em_hetero)) {
  cat("  effect modification (post-policy avg lag0 RR by subgroup):\n"); print(em_hetero)
  fwrite(em_hetero, file.path(out_dir, "effect_modification_summary.csv"))
}

message("== DDD event study complete -> ", out_dir, " ==")
