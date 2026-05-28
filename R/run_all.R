# =============================================================================
# run_all.R  —  run the full analysis pipeline on the simulated dataset.
# Run from this project's root:  Rscript R/run_all.R
# Generates the simulated dataset, then runs the main-effect, DDD and CSDID
# analyses on the simulated dataset.
# =============================================================================
source("R/00_common.R")

message("\n########## 0. simulate data ##########")
source("R/00_simulate_data.R")

message("\n########## 1. main-effect DLNM ##########")
# self-contained (own packages + parallel cluster); run in a fresh process
if (!identical(system2("Rscript", "R/01_main_effect.R"), 0L)) stop("01_main_effect.R failed; aborting pipeline.", call. = FALSE)

message("\n########## 2. build matched panel ##########")
source("R/02_build_matched_panel.R")

message("\n########## 3. DDD event study ##########")
source("R/03_ddd_event_study.R")

message("\n########## 4. CSDID ##########")
source("R/04_ddd_csdid.R")

message("\nAll steps complete. See results/.")
