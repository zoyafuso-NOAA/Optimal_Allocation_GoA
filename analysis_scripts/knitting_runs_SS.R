###############################################################################
## Project:       Knitting Result for univariate STRS optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   
###############################################################################
rm(list = ls())

###############################
## Set up directories
###############################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[1]

github_dir <- paste0(c("/Users/zackoyafuso/Documents", 
                       "C:/Users/Zack Oyafuso/Documents",
                       "C:/Users/zack.oyafuso/Work")[which_machine],
                     "/GitHub/Optimal_Allocation_GoA/")

###############################
## Import required packages
###############################
library(devtools)
devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata")
library(SamplingStrata)

###########################
## Load Data
###########################
load(paste0(github_dir, "/data/optimization_data.RData"))

for(idom in c("full_domain", "district")) {
  
  frame <- switch( idom,
                   "full_domain" = frame_all,
                   "district" = frame_district)
  
  n_dom <- length(unique(frame$domainvalue))
  
  ###########################
  ## Empty Result objects
  ###########################
  master_res_df <- data.frame(id = 1:n_cells)
  master_settings <- master_settings_district <- data.frame()
  master_strata_list <- master_strata_stats_list <- list()
  master_tradeoff <- list()
  
  ##########################
  ##########################
  id <- 0
  for (ispp in 1:ns_all) {
    
    ## For a given species and boat scenario, collect all runs
    runs = dir(paste0(github_dir, "results/", idom, 
                      "/Single_Species_Optimization/",
                      gsub(x = sci_names_all[ispp], 
                           pattern = " ", 
                           replacement = "_")), 
               full.names = T)
    
    ## For each run
    if (length(runs) > 0) {
      for (irun in runs) {
        
        #Load run
        temp_file <- paste0(irun, "/result_list.RData")
        
        if (file.exists(temp_file)) {
          id <- id + 1
          
          load(temp_file)
          
          ## master_settings_district: result of optimization (CV, sample size) 
          ## by district (aka domain)
          master_settings_district <- rbind(
            master_settings_district,
            data.frame(id = id,
                       domain = 1:n_dom,
                       spp = ispp,
                       n = tapply(X = result_list$sum_stats$Allocation, 
                                  INDEX = result_list$sum_stats$Domain, 
                                  FUN = sum),
                       cv = as.numeric(result_list[[3]]))
          )
          
          ## master_settings: result of optimization (CV, sample size) 
          ## distict-aggregated 
          agg_strata <- result_list$solution$aggr_strata
          agg_strata$STRATO <- 1:nrow(agg_strata)
          agg_strata$DOM1 <- 1
          
          master_settings <- rbind(
            master_settings,
            data.frame(id = id,
                       spp = ispp,
                       n = result_list$n,
                       cv = as.numeric(SamplingStrata::expected_CV(
                         strata = agg_strata) ) )
          )
          
          #master_res_df: solution (which cell belongs to which stratum?)
          ## Solution: which strata is assigned to each extrapolation cell
          solution <- 
            switch(idom,
                   "full_domain" = result_list$solution$indices$X1,
                   "district" = as.factor(paste(
                     result_list$solution$framenew$DOMAINVALUE,
                     result_list$solution$framenew$STRATO))
                   
            )
          
          solution <- as.integer(solution)
          
          master_res_df <- cbind(master_res_df,
                                 solution)
          
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
  
  
  ####################################
  ## Subset those solutions that correspond to 1, 2, and 3 boats
  ####################################
  sol_idx <- c()
  
  for (ispp in 1:ns_all) {
    for (isample in samples) {
      #Find solution closet to isample, append to sol_idx
      sol_idx <- 
        c(sol_idx, with(master_settings[master_settings$spp == ispp, ],
                        id[which.min(abs(n - isample))])
        )
    }
  } 
  
  settings_agg <- master_settings[sol_idx,]
  settings_agg$iboat <- rep(1:n_boats, times = ns_all)
  res_df <- master_res_df[, 1 + sol_idx]
  strata_list <- master_strata_list[sol_idx]
  strata_stats_list <- master_strata_stats_list[sol_idx]
  
  vars_to_save <- c("settings_agg", "res_df", "strata_list", "strata_stats_list")
  
  for (ivar in vars_to_save) {
    assign(x = paste0(ivar, "_", idom),
           value = get(ivar))
  }
  
  if (idom == "district") {
    temp_settings_district <- subset(x = master_settings_district, 
                                     subset = id %in% sol_idx )
    
    temp_settings_district <- 
      tidyr::spread(data = subset(x = temp_settings_district,
                                  select = -n),
                    value = cv,
                    key = domain)
    
    settings_district <- data.frame()
    
    for (idx in sol_idx) {
      settings_district <- rbind(settings_district,
                                 subset(temp_settings_district, 
                                        subset = id %in% idx))  
    }
    settings_district$iboat = rep(1:n_boats, times = ns_all)
    vars_to_save <- c(vars_to_save, "settings")
  }
  
  ####################################
  ## Save
  ####################################
  save(list = paste0(vars_to_save, "_", idom),
       file = paste0(github_dir, 
                     "results/", idom, 
                     "/Single_Species_Optimization/",
                     "optimization_knitted_results.RData"))
  
}


