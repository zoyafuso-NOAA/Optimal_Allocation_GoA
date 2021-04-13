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
library(devtools)
devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata")
library(SamplingStrata)

###########################
## Load Data
###########################
load(paste0(github_dir, "/data/optimization_data.RData"))

for(idom in c("full_domain", "district")[1]) {
  
  n_dom <- ifelse(idom == "full_domain", 1, 5)
  
  frame <- switch( idom,
                   "full_domain" = frame_all,
                   "district" = frame_district)
  
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
                      common_names_all[ispp], "/"), 
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
                       species = common_names_all[ispp],
                       domain = 1:n_dom,
                       n = tapply(X = result_list$sum_stats$Allocation, 
                                  INDEX = result_list$sum_stats$Domain, 
                                  FUN = sum),
                       cv = as.numeric(result_list[[3]]))
          )
          
          ## master_settings: result of optimization (CV, sample size) 
          ## distict-aggregated 
          agg_strata <- result_list$solution$aggr_strata
          
          master_settings <- rbind(
            master_settings,
            data.frame(species = common_names_all[ispp],
                       id = id,
                       n = result_list$n,
                       cv = as.numeric(SamplingStrata::expected_CV(
                         strata = agg_strata) ) )
          )
          
          #master_res_df: solution (which cell belongs to which stratum?)
          ## Solution: which strata is assigned to each extrapolation cell
          solution <-paste0("DOM", result_list$solution$framenew$DOMAINVALUE,
                            " STR", result_list$solution$framenew$STRATO)
          
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
      temp_idx <- master_settings$species == common_names_all[ispp]
      
      #Find solution closet to isample, append to sol_idx
      sol_idx <- c(sol_idx, 
                   with(master_settings[temp_idx, ],
                        id[which.min(abs(n - isample))])
      )
    }
  } 
  
  settings_agg <- master_settings[sol_idx, -2]
  settings_agg$boat <- rep(1:n_boats, 
                           times = length(unique(settings_agg$species)))
  res_df <- master_res_df[, 1 + sol_idx]
  strata_list <- master_strata_list[sol_idx]
  strata_stats_list <- master_strata_stats_list[sol_idx]
  
  temp_settings_district <- subset(x = master_settings_district, 
                                   subset = id %in% sol_idx)
  temp_settings_district <-tidyr::spread(data = temp_settings_district,
                                         value = cv,
                                         key = domain)
  
  settings_district <- data.frame()
  for (idx in sol_idx) {
    settings_district <- rbind(settings_district,
                               subset(temp_settings_district, 
                                      subset = id %in% idx))  
  }
  settings_district$boat = rep(1:n_boats, 
                               times = length(unique(settings_district$species)))
  settings_district <- subset(x = settings_district,
                              select = -id)
  
  old_strata_stats_list <- strata_stats_list
  ####################################
  ## Adjust CVs to match sample sizes
  ####################################
  for (irow in 1:nrow(settings_agg)) {
    
    ## Temporary constants
    temp_boat <- settings_agg[irow, "boat"]
    temp_cv <- switch(idom,
                      "full_domain" = settings_agg[irow, "cv"],
                      "district" = as.numeric(settings_district[irow, paste(1:5)]))
    temp_n <- settings_agg[irow, "n"]
    
    temp_strata_stats_list <- subset(x = strata_stats_list[[irow]],
                                     select = -SOLUZ)
    n_dom <- ifelse(idom == "full_domain", 1, 5)
    
    ## Temporary variables
    temp_n_by_strata <- actual_n_by_strata <- 
      vector(length = nrow(temp_strata_stats_list))
    names(temp_n_by_strata) <- names(actual_n_by_strata) <- 
      temp_strata_stats_list$DOM1
    
    ## aiming for an optimal sample size within 2 stations of target
    outside_sample_bounds <- !(temp_n >= (samples[temp_boat] - 1) & 
                                 temp_n <= (samples[temp_boat] + 1)) 
    
    ## if outside that optimal window...
    if (outside_sample_bounds) {
      while (outside_sample_bounds) {
        
        ## Adjust CV slightly depending on current sample size
        over_or_under <- temp_n > samples[temp_boat]
        
        temp_cv <- temp_cv * (1 + ifelse(test = over_or_under == T, 
                                         yes = 0.0001,
                                         no = -0.0001))
        temp_n_by_reg <- vector(length = n_dom)
        
        ## Run bethel algorithm under new temp_cv for each domain
        for (i in 1:n_dom) {
          error_df <- data.frame("DOM" = "DOM1",
                                 "CV1" = temp_cv[i],
                                 "domainvalue"  = 1)
          
          temp_bethel <- SamplingStrata::bethel(
            errors = error_df,
            stratif = subset(x = temp_strata_stats_list, subset = DOM1 == i), 
            realAllocation = T, 
            printa = T)
          
          #Record bethel algorithm run
          temp_n_by_reg[i] <- sum(as.numeric(ceiling(temp_bethel)))
          temp_n_by_strata[names(temp_n_by_strata) == paste(i)] <- 
            ceiling(temp_bethel)
          actual_n_by_strata[names(actual_n_by_strata) == paste(i)] <- 
            temp_bethel
        }
        
        ## Update temp_n, repeat process until temp_n is within the optimal
        ## total sample size window
        temp_n <- sum(temp_n_by_reg)
        print(temp_n)
        
        outside_sample_bounds <- !(temp_n >= (samples[temp_boat] - 1) & 
                                     temp_n <= (samples[temp_boat] + 1)) 
      }
      
      ####################################
      ## Record new results
      ####################################
      temp_cv_agg <- SamplingStrata::expected_CV(
        strata = cbind(subset(temp_strata_stats_list,
                              select = -DOM1),
                       DOM1 = 1,
                       SOLUZ = actual_n_by_strata))
      
      settings_agg[irow, "cv"] <- temp_cv_agg
      settings_agg[irow, "n"] <- temp_n
      strata_stats_list[[irow]]$SOLUZ <- actual_n_by_strata 
      strata_list[[irow]]$Allocation <- temp_n_by_strata
      settings_district[irow, paste(1:n_dom)] <- temp_cv
    }
  }
  
  ###################################
  ## Assign result variables with domain name attached at the end
  ####################################
  vars_to_save <- c("settings_agg", "res_df", "settings_district",
                    "strata_list", "strata_stats_list")
  
  for (ivar in vars_to_save) assign(x = paste0(ivar, "_", idom),
                                    value = get(ivar))
  
  ####################################
  ## Save
  ####################################
  save(list = paste0(vars_to_save, "_", idom),
       file = paste0(github_dir, 
                     "results/", idom, 
                     "/Single_Species_Optimization/",
                     "optimization_knitted_results.RData"))
  
}


