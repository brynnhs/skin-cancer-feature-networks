# Master script to run all phases

setwd(here::here())

source("R/01_load_validate.R")
source("R/02_preprocess.R")
source("R/03_split_conditions.R")
source("R/04_mgm_baseline.R")
source("R/05_stars_selection.R")
source("R/06_mgm_bootstrap.R")
source("R/07_network_metrics.R")
source("R/08_pcstable_baseline.R")
source("R/09_causal_bootstrap.R")
