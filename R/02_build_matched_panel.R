# =============================================================================
# 02_build_matched_panel.R  —  build the matched daily panel for the DDD analysis
#
# Pipeline:
#   raw annual + daily inputs -> annual city panel -> pre-policy covariate means
#   -> Mahalanobis 1:1 nearest-neighbour matching per cohort -> matched annual
#   panel (interpolated) -> expand to daily -> join daily climate/case variables
#   -> write results/derived/matched_daily_panel.csv
# Run after sourcing R/00_common.R (run_all.R does this).
# =============================================================================

message("== building matched daily panel ==")

# ---- 1. read inputs ----------------------------------------------------------
infra <- fread(file.path(DATA_DIR, "annual_infrastructure.csv"), encoding = "UTF-8")   # annual infrastructure
reg   <- fread(file.path(DATA_DIR, "city_lookup.csv"),  encoding = "UTF-8")    # city -> region lookup
econ  <- fread(file.path(DATA_DIR, "annual_economy.csv"), encoding = "UTF-8")   # annual economic
daily <- fread(file.path(DATA_DIR, "daily_panel.csv"), encoding = "UTF-8") # daily climate / cases

# annual_economy.csv carries three leading index columns (province/city names).
econ <- econ[, -(1:3)]

# keep only the annual covariates used downstream
infra <- infra[, .SD, .SDcols = intersect(c("code","year","water","wspd","dpd","str","gcr","gdr","pt","invest"), names(infra))]
econ_keep <- intersect(c("code","year","gdp","2_gdp","3_gdp","expenditure","st","bed","doc","areas","den"), names(econ))
econ <- econ[, .SD, .SDcols = econ_keep]

# ---- 2. annual city panel ----------------------------------------------------
panel <- CJ(code = unique(reg$code), year = YEARS)
panel <- merge(panel, infra, by = c("code","year"), all.x = TRUE)
panel <- merge(panel, econ,  by = c("code","year"), all.x = TRUE)

# annual mean temperature / total precipitation, and region/climate, from daily
daily[, date := as.Date(date)]
annual_clim <- daily[, .(tem = mean(tem, na.rm = TRUE), pre = sum(pre, na.rm = TRUE)),
                     by = .(code, year = year(date))]
city_info <- unique(daily[, .(code, region, climate_code)], by = "code")
panel <- merge(panel, annual_clim, by = c("code","year"), all.x = TRUE)
panel <- merge(panel, city_info, by = "code", all.x = TRUE)

# treatment status + cohort policy year
panel[, sponge := as.integer(code %in% ALL_TREAT)]
panel[, policy_year := fifelse(code %in% TREAT_2015, 2015L, fifelse(code %in% TREAT_2016, 2016L, NA_integer_))]

# ---- 3. pre-policy covariate means + matching --------------------------------
# keep cities with >= 3 pre-policy observations
have <- panel[year %in% PRE_YEARS, .(n = .N), by = code][n >= 3L, code]
panel <- panel[code %in% have]
setorder(panel, code, year)

# Matching covariate means use PRE-POLICY information only: gaps are interpolated
# within PRE_YEARS so no post-policy value can leak into a matching covariate, and
# cities are screened on the matching covariates themselves (not on variables that
# never enter the matching formula). The time-varying socioeconomic covariates for
# the downstream daily panel are interpolated separately in step 4.
match_src <- c("tem", "water", "str", "den")
pre_panel <- copy(panel[year %in% PRE_YEARS, .SD, .SDcols = c("code", "year", match_src)])
setorder(pre_panel, code, year)
for (v in match_src) pre_panel[, (v) := safe_interpolate(get(v), year), by = code]
pre_panel <- pre_panel[!code %in% pre_panel[, .(bad = any(vapply(.SD, function(x) sum(!is.na(x)) < 2L, logical(1)))),
                                            by = code, .SDcols = match_src][bad == TRUE, code]]
city_means <- pre_panel[, .(tem_mean = mean(tem, na.rm = TRUE), water_mean = mean(water, na.rm = TRUE),
                            str_mean = mean(str, na.rm = TRUE), den_mean = mean(den, na.rm = TRUE)), by = code]
city_means <- merge(city_means, unique(panel[, .(code, sponge, policy_year)]), by = "code")
city_means <- city_means[complete.cases(city_means[, ..MATCH_COVARIATES])]

match_cohort <- function(treat_codes, pool) {
  treat <- city_means[code %in% treat_codes & sponge == 1L]
  if (nrow(treat) == 0L || nrow(pool) == 0L) return(NULL)
  fm <- as.formula(paste("sponge ~", paste(MATCH_COVARIATES, collapse = " + ")))
  m  <- MatchIt::matchit(fm, data = rbind(treat, pool), method = "nearest",
                         distance = "mahalanobis", replace = FALSE, ratio = 1)
  md <- as.data.table(MatchIt::match.data(m))
  md[, .(treat_code = code[sponge == 1L], control_code = code[sponge == 0L]), by = subclass]
}
set.seed(123)   # fixed seed for reproducible matching
pool   <- city_means[sponge == 0L]
p2015  <- match_cohort(TREAT_2015, pool)
pool   <- pool[!code %in% p2015$control_code]
p2016  <- match_cohort(TREAT_2016, pool)
all_pairs <- rbind(p2015[, policy_year := 2015L], p2016[, policy_year := 2016L])
fwrite(all_pairs[, .(policy_year, treat_code, control_code)], file.path(DERIVED_DIR, "matched_pairs.csv"))
message(sprintf("  matched %d treated-control pairs", nrow(all_pairs)))

# controls inherit the policy year of the treated city they were matched to
control_policy <- unique(all_pairs[, .(code = control_code, policy_year)])
matched_cities <- rbind(
  unique(city_means[sponge == 1L, .(code, sponge, policy_year)]),
  control_policy[, .(code, sponge = 0L, policy_year)]
)

# ---- 4. matched annual panel (interpolate socioeconomic covariates) ----------
interp_vars <- intersect(c("water","wspd","dpd","str","gcr","gdr","pt","invest",
                           "gdp","2_gdp","3_gdp","expenditure","st","bed","doc","areas","den"), names(panel))
matched_annual <- CJ(code = matched_cities$code, year = YEARS) |>
  merge(matched_cities, by = "code") |>
  merge(panel[, .SD, .SDcols = c("code","year","region","climate_code", interp_vars)],
        by = c("code","year"), all.x = TRUE)
setorder(matched_annual, code, year)
for (v in interp_vars) matched_annual[, (v) := locf_fill(safe_interpolate(get(v), year)), by = code]
matched_annual[, region := zoo::na.locf(zoo::na.locf(region, na.rm = FALSE), fromLast = TRUE), by = code]
matched_annual[, climate_code := climate_code[which(!is.na(climate_code))[1]], by = code]

# ---- 5. expand to daily and join daily climate / cases -----------------------
daily_grid <- CJ(code = matched_cities$code,
                 date = seq(min(daily$date), max(daily$date), by = "day"))
daily_grid[, year := year(date)]
clim_cols <- intersect(c("code","date","tem","pre","RH","case","Pre_school","School","Adults","Elders","GF","SF"), names(daily))
matched_daily <- daily_grid |>
  merge(matched_annual, by = c("code","year"), all.x = TRUE) |>
  merge(daily[, ..clim_cols], by = c("code","date"), all.x = TRUE)

fwrite(matched_daily, file.path(DERIVED_DIR, "matched_daily_panel.csv"))
message(sprintf("  wrote matched daily panel: %d rows x %d cols (%d cities)",
                nrow(matched_daily), ncol(matched_daily), uniqueN(matched_daily$code)))
