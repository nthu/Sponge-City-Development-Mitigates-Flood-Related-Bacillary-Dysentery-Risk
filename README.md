# Sponge City Development & Flood-Related Bacillary Dysentery — analysis code

Analysis code for the study of (1) the short-term association between floods and
bacillary dysentery (BD) and (2) the effect of China's Sponge City Development
(SCD) policy on that association. The policy analysis is a
difference-in-difference-in-differences (DDD) event study. This repository
runs the full workflow on a **simulated** dataset.

> **Data note.** The real surveillance data are not public (see
> [Data availability](#data-availability)). `R/00_simulate_data.R` deterministically
> generates a synthetic dataset with the same schema as the real inputs. The
> simulated data reproduce the *structure* and *qualitative direction* of the
> findings, not the published point estimates. Running the same code on the real
> data reproduces the study's quantitative results.

## What the code does

- **Main effect** (`R/01_main_effect.R`): time-stratified case-crossover design +
  distributed-lag non-linear model (DLNM), conditional negative binomial.
- **Policy evaluation** (`R/02`–`R/04`): builds a matched city panel (Mahalanobis
  1:1 matching), then a DDD event study combined with a DLNM. Covariates are GDP,
  water-supply rate, sewage-treatment rate, population density, and doctors per
  capita. Reports per-year and cumulative RRs for the total and four age groups,
  parallel-trend and placebo tests, a **saturated-model supplementary
  specification** (refit with the bare `Treated × 1[RelYear=τ]` main effects), an
  effect-modification analysis (post-policy RR within the high/low half of seven
  pre-policy city characteristics), and a complementary Callaway & Sant'Anna
  (CSDID) estimator.

## Repository structure

```
.
├── R/
│   ├── 00_simulate_data.R       # generate the simulated dataset -> data/
│   ├── 00_common.R              # configuration + shared utilities
│   ├── 01_main_effect.R         # main flood→BD DLNM
│   ├── 02_build_matched_panel.R # matching -> matched daily panel
│   ├── 03_ddd_event_study.R     # DDD event study (DLNM); case + age groups; saturated + effect modification
│   ├── 04_ddd_csdid.R           # Callaway–Sant'Anna estimator
│   └── run_all.R                # run the whole pipeline
├── data/        # simulated inputs (generated; git-ignored)
├── results/     # outputs (generated; git-ignored)
├── LICENSE      # MIT
└── README.md
```

## 1. System requirements

- **OS:** developed and tested on **macOS** (Darwin 25.5). The code uses only
  portable constructs — relative paths, PSOCK parallel clusters (not Unix-only
  `fork`/`mclapply`), CRAN-only packages, and ASCII-only simulated data — so it is
  expected to run unmodified on Linux and Windows, though it has not been tested
  there.
- **R:** ≥ 4.1; developed/tested on **R 4.5.1**.
- **R packages** (CRAN; tested versions): `data.table` 1.17.8, `dplyr` 1.1.4,
  `tidyr` 1.3.1, `lubridate` 1.9.4, `zoo` 1.8.14, `MatchIt` 4.7.2, `cobalt` 4.6.2,
  `fixest` 0.13.2, `dlnm` 2.4.10, `splines` (base), `MASS` 7.3.65, `did` 2.3.0,
  `ggplot2` 4.0.0, `scales` 1.4.0, `foreach`, `doParallel`, `gnm` 1.1.5.
- **Hardware:** no non-standard hardware; ≥ 8 GB RAM is sufficient for the demo.

## 2. Installation

```bash
git clone <repository-url> && cd Sponge-City-Development-Mitigates-Flood-Related-Bacillary-Dysentery-Risk
```
Packages auto-install on first run via `pacman`. To pre-install:
```r
install.packages("pacman")
pacman::p_load(data.table, dplyr, tidyr, lubridate, zoo, MatchIt, cobalt, fixest,
               dlnm, splines, MASS, did, ggplot2, scales, foreach, doParallel, gnm)
```
**Typical install time:** ~5–15 min on a normal desktop (seconds if already installed).

## 3. Demo

From the repository root:
```bash
Rscript R/run_all.R
```
This generates the simulated data and runs the main-effect, DDD, and CSDID
analyses. **Expected run time:** ~1–2 minutes total on a normal desktop.

**Expected output** (under `results/`):
- `results/<run_id>/` (main effect): single-lag & cumulative RR tables; the
  flood→BD lag-0 RR is > 1 and decays over the lag window.
- `results/ddd_event_study/`: per-year single-lag & cumulative DDD RRs for the total
  and four age groups; `parallel_trend_test.csv`; `placebo_lead2_*`; and
  `saturated_sensitivity_lag0.csv` (per-year RRs from the saturated supplementary
  specification, alongside the main model); and
  `effect_modification_summary.csv` (post-policy average lag-0 RR for the high/low
  half of seven pre-policy city characteristics). On the simulated data: flat
  pre-trends and a post-policy protective effect, similar across subgroups.
- `results/ddd_csdid/`: CSDID overall ATT and event-study table.

## 4. Instructions for use

To run on other data, place files with the same schema as the simulated ones in
`data/` (see `R/00_simulate_data.R` for the exact columns) and run `R/run_all.R`
(skipping `00_simulate_data.R`). `R/01_main_effect.R` defaults to the full analysis;
for a fast subset, run with `RUN_MODE=pilot` (or temporarily edit the `run_mode` line).

## Data availability

The surveillance data are available from the Chinese Center for Disease Control
and Prevention under restrictions (licensed use) and are not publicly available;
they are available from the authors on reasonable request with CDC permission.
Meteorological data are from ERA5-Land.

## License

MIT — see [LICENSE](LICENSE).
