###############################################################################
## Project:       Knitting Result for univariate STRS optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   
###############################################################################
rm(list = ls())

###############################
## Set up directories
###############################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[3]

github_dir <- paste0(c("/Users/zackoyafuso/Documents", 
                       "C:/Users/Zack Oyafuso/Documents",
                       "C:/Users/zack.oyafuso/Work")[which_machine],
                     "/GitHub/Optimal_Allocation_GoA/")

###############################
## Import required packages
###############################
library(sp)
library(RColorBrewer)
library(raster)

###########################
## Load Data
###########################
load(paste0(github_dir, "/data/optimization_data.RData"))
load(paste0(github_dir, "/data/Extrapolation_depths.RData"))

###########################
## Empty Result objects
###########################
master_res_df <- data.frame(id = 1:N)
master_settings <- data.frame()
master_strata_list <- master_strata_stats_list <- list()
master_tradeoff <- list()

##########################
##########################

for (ispp in 1:ns_opt) {
  for (iboat in 1:nboats) {
    
    ## For a given species and boat scenario, collect all runs
    runs = dir(paste0(github_dir,
                      "results/Single_Species_Optimization/",
                      gsub(x = sci_names_opt[ispp], 
                           pattern = " ", 
                           replacement = "_"),
                      "/boat", iboat), 
               full.names = T)
    nruns = length(runs)
    
    ## For each run
    if (nruns > 0) {
      for (irun in 1:nruns) {
        
        #Load run
        temp_file <- paste0(github_dir, 
                            "results/Single_Species_Optimization/",
                            gsub(x = sci_names_opt[ispp], 
                                 pattern = " ", 
                                 replacement = "_"),
                            "/boat", iboat,
                            "/Run", irun, "/result_list.RData")
        
        if (file.exists(temp_file)) {
          
          load(temp_file)
          
          #master_settings: result of optimization (CV, sample size)
          master_settings <- rbind(
            master_settings,
            data.frame(irun = irun,
                       iboat = iboat,
                       ispp = ispp,
                       n = with(result_list$sum_stats, tapply(Allocation, Domain, sum)),
                       cv = as.numeric(result_list[[3]]))
          )
          
          #master_res_df: solution (which cell belongs to which stratum?)
          master_res_df <- cbind(master_res_df,
                                 result_list[[1]]$indices$X1)
          
          #master_strata_list: stratum-level details of solution
          master_strata_list <- c(master_strata_list, 
                                  list(result_list[[2]]))
          
          #master_strata_stats_list: stratum-level means and variances
          master_strata_stats_list <- c(master_strata_stats_list, 
                                        list(result_list$solution$aggr_strata))
        }
      }
    }
  }
}

####################################
## Subset those solutions that correspond to 1, 2, and 3 boats
####################################
master_settings$id = 1:nrow(master_settings)
master_settings_agg <- aggregate(n ~ irun + iboat + ispp, 
                                 data = master_settings, sum)

sol_idx <- c()

for (ispp in sort(unique(master_settings$ispp)) ) {
  for (isample in samples) {
    #Find solution closet to isample, append to sol_idx
    sol_idx <- c(sol_idx, 
                 with(master_settings[master_settings$ispp == ispp, ],
                      id[which.min(abs(n - isample))])
    )
  }
} 

settings <- master_settings[sol_idx, 1:4]
# res_df <- master_res_df[, 1 + sol_idx]
# strata_list <- master_strata_list[sol_idx]
# strata_stats_list <- master_strata_stats_list[sol_idx]
# 
# ####################################
# ## Save
# ####################################
# save(list = c("res_df", "settings", "strata_list", "strata_stats_list"),
#      file = paste0(github_dir, 
#                    "results/Single_Species_Optimization/", 
#                    "optimization_knitted_results.RData"))
