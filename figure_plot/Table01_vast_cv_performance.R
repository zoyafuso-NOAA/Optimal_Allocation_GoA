###############################################################################
## Project:       Table 2
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   
###############################################################################
rm(list = ls())

##################################################
####   Set up directories based on whether the optimization is being conducted
####        on a multi-species or single-species level
##################################################
github_dir <- getwd()

##################################################
####   Load Data
####   Load Population CVs for use in the thresholds
##################################################
load("data/processed/optimization_data.RData")
load("data/processed/prednll_VAST_models.RData")

pred_jnll[match( common_names_all[c(spp_idx_opt, spp_idx_eval)], pred_jnll$spp_name), ]

rrmse[match( common_names_all[c(spp_idx_opt, spp_idx_eval)], rrmse$spp_name), ]
