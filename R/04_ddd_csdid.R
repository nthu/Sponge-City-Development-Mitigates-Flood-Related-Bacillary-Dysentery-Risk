# =============================================================================
# 04_ddd_csdid.R  —  Callaway & Sant'Anna (2021) group-time ATT (CSDID)
#
# Complementary, robust-to-staggered-adoption estimator of the sponge-city
# policy effect on ANNUAL total bacillary-dysentery counts.
#
# Run after 02 (needs results/derived/matched_daily_panel.csv).
# =============================================================================

message("== CSDID (Callaway & Sant'Anna) ==")
if (!requireNamespace("did", quietly = TRUE)) install.packages("did", repos = "https://cloud.r-project.org")
library(did)

panel <- fread(file.path(DERIVED_DIR, "matched_daily_panel.csv"), encoding = "UTF-8")
panel[, date := as.Date(date)]

# annual total cases + annual covariate means per matched city
annual <- panel[, .(total_cases = sum(case, na.rm = TRUE),
                    water = mean(water, na.rm = TRUE), str = mean(str, na.rm = TRUE),
                    den = mean(den, na.rm = TRUE), tem = mean(tem, na.rm = TRUE)),
                by = .(code, year)]
annual[, sponge := as.integer(code %in% ALL_TREAT)]
# first-treatment year g: policy year for treated, 0 for never-treated controls.
# g is numeric (did encodes never-treated as Inf internally, requiring a double type).
annual[, g := as.numeric(fifelse(code %in% TREAT_2015, 2015L, fifelse(code %in% TREAT_2016, 2016L, 0L)))]
annual[, log_cases := log1p(total_cases)]
annual <- annual[complete.cases(annual[, .(log_cases, water, str, den, tem)])]
annual[, `:=`(code = as.integer(code), year = as.integer(year))]

out_dir <- file.path(RESULTS_DIR, "ddd_csdid")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

att <- att_gt(yname = "log_cases", tname = "year", idname = "code", gname = "g",
              xformla = ~ water + str + den + tem, data = as.data.frame(annual),
              panel = TRUE, allow_unbalanced_panel = TRUE,
              control_group = "nevertreated", est_method = "dr", clustervars = "code")

agg_simple  <- aggte(att, type = "simple",  na.rm = TRUE)
agg_dynamic <- aggte(att, type = "dynamic", na.rm = TRUE)
cat(sprintf("  overall ATT (log cases) = %.3f (SE %.3f)\n", agg_simple$overall.att, agg_simple$overall.se))

fwrite(data.table(att = agg_simple$overall.att, se = agg_simple$overall.se,
                  ci_lo = agg_simple$overall.att - 1.96 * agg_simple$overall.se,
                  ci_hi = agg_simple$overall.att + 1.96 * agg_simple$overall.se),
       file.path(out_dir, "csdid_overall_att.csv"))
fwrite(data.table(event_time = agg_dynamic$egt, att = agg_dynamic$att.egt, se = agg_dynamic$se.egt),
       file.path(out_dir, "csdid_event_study.csv"))
message("== CSDID complete -> ", out_dir, " ==")
