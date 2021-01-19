###############################################################################
## Project:         Plot Multispecies Optimization Solutions
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     
###############################################################################
rm(list = ls())

##################################################
####  Install a forked version of the SamplingStrata Package from 
####  zoyafuso-NOAA's Github page
####
####  Import other required packages
##################################################
library(sp)
library(RColorBrewer)
library(raster)

##################################################
####   Set up directories based on whether the optimization is being conducted
####        on a multi-species or single-species level
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[1]

github_dir <- paste0(c("/Users/zackoyafuso/Documents", 
                       "C:/Users/Zack Oyafuso/Documents",
                       "C:/Users/zack.oyafuso/Work")[which_machine],
                     "/GitHub/Optimal_Allocation_GoA/")

output_dir <- paste0(c("/Users/zackoyafuso/Google Drive/"),
                     "MS_Optimizations/TechMemo/figures/")

##################################################
####   Load Data
####   Load Population CVs for use in the thresholds
##################################################
load(paste0(github_dir, "data/optimization_data.RData"))
load(paste0(github_dir, "data/Extrapolation_depths.RData"))
load(paste0(github_dir, "results/MS_optimization_knitted_results.RData"))

cells_by_district <- list()

for(sol_idx in paste0("sol_", settings$id)) {
 cells_by_district[[sol_idx]] <- do.call(rbind, 
                                         tapply(X = district_vals, 
                                                INDEX = res_df[, sol_idx], 
                                                FUN = table) )
 cells_by_district[[sol_idx]] <- 
  round(sweep(x = cells_by_district[[sol_idx]],
              MARGIN = 1,
              STATS = rowSums(cells_by_district[[sol_idx]]),
              FUN = "/"), 2)
}


