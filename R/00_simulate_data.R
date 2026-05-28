# =============================================================================
# 00_simulate_data.R  —  generate a SIMULATED dataset (public demo)
#
# The real surveillance data are restricted and cannot be published. This script
# deterministically generates synthetic raw inputs with the SAME schemas as the
# real files, so the full analysis pipeline (02 -> 03 -> 04, and 01) runs
# end-to-end. The simulated data reproduce the STRUCTURE and qualitative direction
# of the findings (a short-term flood->BD association; an attenuation in sponge
# cities after the policy) but not the point estimates of the real study.
#
# Writes into data/:  daily_panel.csv, city_lookup.csv, annual_infrastructure.csv,
# annual_economy.csv  — matching the real raw-file schemas.
#
# by R/02_build_matched_panel.R.
# Run from the project root:  Rscript R/00_simulate_data.R
# =============================================================================

suppressWarnings(suppressMessages(library(data.table)))
set.seed(105)
dir.create("data", recursive = TRUE, showWarnings = FALSE)

# 2011-2020 keeps the matching pre-period (2011-2014) and the full -4..+4 DDD
# event window around the 2015/2016 policy cohorts, while keeping the demo light.
years    <- 2011:2020
date_seq <- seq(as.Date("2011-01-01"), as.Date("2020-12-31"), by = "day")
n_days   <- length(date_seq)

treat2015 <- c(2208, 5109, 3417, 3701, 4501, 3211, 3603, 4106, 3304, 4201, 4307, 3502, 5000)
treat2016 <- c(3100, 5304, 1100, 1200, 6404, 3702, 4403, 4404, 3302, 4602, 6210, 3501, 2102, 6301)
all_treat <- c(treat2015, treat2016)
control_pool <- 7001:7030   # control pool to match the 27 treated cities against
codes  <- c(all_treat, control_pool)
n_city <- length(codes)
regions <- c("North", "East", "South", "West")

meta <- data.table(
  code = codes, city = paste0("City_", codes),
  region = sample(regions, n_city, replace = TRUE),
  climate_code = sample(1:4, n_city, replace = TRUE),
  is_treated = as.integer(codes %in% all_treat),
  policy_year = fifelse(codes %in% treat2015, 2015L, fifelse(codes %in% treat2016, 2016L, NA_integer_)),
  base_temp = runif(n_city, 8, 24), base_lograte = runif(n_city, -0.5, 1.2))

doy <- as.integer(format(date_seq, "%j")); yr <- as.integer(format(date_seq, "%Y")); mon <- as.integer(format(date_seq, "%m"))

# ---- daily climate + flood + BD counts (daily_panel.csv, 42-col schema) --
build_city <- function(i) {
  m <- meta[i]
  seas_t <- 12 * sin(2 * pi * (doy - 120) / 365)
  tem <- m$base_temp + seas_t + rnorm(n_days, 0, 3); d2m <- tem - runif(n_days, 2, 9)
  es <- 611.21 * exp(17.502 * tem / (240.97 + tem)); ea <- 611.20 * exp(17.502 * d2m / (240.97 + d2m))
  RH <- pmin(100, pmax(10, 100 * ea / es))
  wet <- 0.5 + 0.5 * pmax(0, sin(2 * pi * (doy - 120) / 365))
  pre <- rgamma(n_days, shape = 0.4 * wet + 0.05, scale = 8)
  swvl1 <- pmin(0.6, pmax(0.02, 0.18 + 0.08 * wet + 0.0008 * pre + rnorm(n_days, 0, 0.03)))
  ro <- pmax(0, 0.0005 * pre + rgamma(n_days, shape = 0.3, scale = 0.01))
  GF <- rbinom(n_days, 1, plogis(-4.4 + 1.1 * wet + 0.02 * pmax(0, pre - 15)))
  SF <- GF * rbinom(n_days, 1, 0.15); fl <- GF == 1
  pre[fl] <- pre[fl] + runif(sum(fl), 55, 140) + 60 * SF[fl]
  swvl1[fl] <- pmin(0.6, swvl1[fl] + runif(sum(fl), 0.08, 0.2)); ro[fl] <- ro[fl] + runif(sum(fl), 0.02, 0.12)
  # Flood "burden" = exponentially-weighted flood history (~10-day memory). General
  # and severe floods enter as separate burdens, so a severe flood carries an extra
  # effect on top of a general one.
  w <- exp(-(0:10) / 4)
  burden_g <- as.numeric(stats::filter(GF, w, sides = 1)); burden_g[is.na(burden_g)] <- 0
  burden_s <- as.numeric(stats::filter(SF, w, sides = 1)); burden_s[is.na(burden_s)] <- 0
  # Sponge-city policy: in treated cities the flood->BD effect is attenuated after the
  # policy start, ramping up to full size by event-year 3 (the year-3 peak).
  start <- if (!is.na(m$policy_year)) as.Date(paste0(m$policy_year, "-04-15")) else as.Date(NA)
  yrs_since <- if (!is.na(start)) pmax(0, as.numeric(date_seq - start) / 365.25) else rep(0, n_days)
  post <- if (!is.na(start)) as.integer(date_seq >= start) else rep(0L, n_days)
  protect <- (m$is_treated == 1L) * post * pmin(yrs_since / 3, 1)   # policy dose: 0 -> 1 by year 3
  # Demonstration effect sizes: a baseline flood effect (beta_g for general, beta_s
  # for severe) plus a per-event policy attenuation (rho) in treated cities after
  # the policy start.
  beta_g <- 0.080; beta_s <- 0.060; rho <- 0.170
  flood_effect  <-  beta_g * burden_g + beta_s * burden_s          # baseline flood->BD effect
  policy_effect <- -rho * burden_g * protect                       # treated, post-policy only
  seas_bd <- 0.45 * sin(2 * pi * (doy - 150) / 365)
  dow_e <- c(0,0,0,0,0.02,-0.05,-0.08)[as.integer(format(date_seq, "%u"))]
  hol <- as.integer(format(date_seq, "%m-%d") %in% c("01-01","05-01","10-01","10-02","10-03"))
  base <- m$base_lograte + seas_bd + dow_e - 0.1 * hol + policy_effect
  # age groups differ in flood vulnerability (fsens); policy protection is uniform.
  mk <- function(off, fsens, th) rnbinom(n_days, mu = pmax(exp(base + off + fsens * flood_effect), 1e-3), size = th)
  ps <- mk(-0.2, 1.20, 3); sc <- mk(-0.9, 1.30, 3); ad <- mk(0.3, 0.80, 4); el <- mk(-0.4, 1.00, 3)
  data.table(code=m$code, date=date_seq, tem=round(tem,3), tmax=round(tem+runif(n_days,1,6),3),
    tmin=round(tem-runif(n_days,1,6),3), d2m=round(d2m,3), RH=round(RH,2), swvl1=round(swvl1,4),
    swvl2=round(swvl1+rnorm(n_days,0,0.02),4), pre=round(pre,3), ro=round(ro,5), sro=round(ro*runif(n_days,.4,.9),5),
    case=ps+sc+ad+el, group1=0L,group2=0L,group3=0L,group4=0L,group5=0L,group6=0L,group7=0L,group8=0L,group9=0L,
    region=m$region, year=yr, month=mon,
    season=c("winter","winter","spring","spring","spring","summer","summer","summer","autumn","autumn","autumn","winter")[mon],
    dow=format(date_seq,"%u"), holiday=hol, policy=0L, Pre_school=ps, Adults=ad, School=sc, Elders=el,
    climate_code=m$climate_code, GF=GF, SF=SF, event_type="", GF_1=0L, GF_2=0L, SF_1=0L, SF_2=0L, GF_3=0L)
}
cat("Simulating", n_city, "cities x", n_days, "days ...\n")
daily <- rbindlist(lapply(seq_len(n_city), build_city))
daily[, swvl1_p75 := quantile(swvl1, .75, na.rm=TRUE), by=code][, ro_p95 := quantile(ro, .95, na.rm=TRUE), by=code]
daily[, `:=`(GF_1=as.integer(GF==1L & swvl1>swvl1_p75), GF_2=as.integer(GF==1L & ro>ro_p95),
             SF_1=as.integer(SF==1L & swvl1>swvl1_p75), SF_2=as.integer(SF==1L & ro>ro_p95), GF_3=as.integer(pre>=50))]
daily[, c("swvl1_p75","ro_p95") := NULL]
fwrite(daily, "data/daily_panel.csv", bom = TRUE)

# ---- code -> region lookup (city_lookup.csv) --------------------------------
fwrite(meta[, .(code, regname = city, region, code_1 = code)], "data/city_lookup.csv", bom = TRUE)

# ---- annual panels (real schemas) --------------------------------------------
grid <- CJ(code = codes, year = years); grid <- merge(grid, meta[, .(code, city, is_treated, policy_year)], by = "code")
ny <- nrow(grid); trend <- (grid$year - 2011) / 9

# annual_infrastructure.csv : city,code,sponge,water,wspd,dpd,str,gcr,gdr,pt,invest,year
infra <- data.table(city = grid$city, code = grid$code,
  sponge = fifelse(grid$is_treated == 1L, "sponge", "non_sponge"),
  water = round(pmin(100, 60 + 25*trend + 6*grid$is_treated + rnorm(ny,0,4)), 2),
  wspd  = round(pmax(0, 6 + 6*trend + rnorm(ny,0,1)), 2),
  dpd   = round(pmax(0, 5 + 7*trend + rnorm(ny,0,1)), 2),
  str   = round(pmin(100, 45 + 35*trend + 8*grid$is_treated + rnorm(ny,0,5)), 2),  # sewage-treatment rate
  gcr   = round(pmin(100, 30 + 20*trend + rnorm(ny,0,4)), 2),
  gdr   = round(pmin(100, 50 + 30*trend + rnorm(ny,0,5)), 2),
  pt    = round(pmax(0, 3 + 5*runif(ny) + 3*trend), 2),
  invest= round(pmax(1, 80 + 300*trend + rnorm(ny,0,30)), 1), year = grid$year)
fwrite(infra, "data/annual_infrastructure.csv", bom = TRUE)

# annual_economy.csv : <3 leading index cols> + code,year,gdp,1_gdp,2_gdp,3_gdp,expenditure,st,bed,doc,areas,den
sec <- round(pmax(1, 200 + 600*trend + 200*runif(ny) + rnorm(ny,0,40)), 1)
ter <- round(pmax(1, 150 + 800*trend + 300*runif(ny) + rnorm(ny,0,50)), 1)
econ <- data.table(prov = "Prov", provname = "ProvName", cityname = grid$city,
  code = grid$code, year = grid$year, gdp = round(sec + ter + pmax(1, 50 + 100*trend), 1),
  `1_gdp` = round(pmax(1, 50 + 50*runif(ny)), 1), `2_gdp` = sec, `3_gdp` = ter,
  expenditure = round(pmax(1, 50 + 200*trend + rnorm(ny,0,20)), 1),
  st = round(pmax(0, 2 + 8*trend*runif(ny) + rnorm(ny,0,1)), 2),
  bed = round(pmax(0, 30 + 40*trend + rnorm(ny,0,6)), 1),
  doc = round(pmax(0, 15 + 25*trend + rnorm(ny,0,4)), 1),
  areas = round(pmax(50, 500 + 1000*runif(ny)), 1),
  den = round(pmax(50, 300 + 1500*runif(ny) + 300*grid$is_treated + rnorm(ny,0,50)), 1))
fwrite(econ, "data/annual_economy.csv", bom = TRUE)

cat(sprintf("Done. Flood-day rate %.2f%%, mean daily BD %.2f. Files written to data/.\n",
            mean(daily$GF) * 100, mean(daily$case)))
