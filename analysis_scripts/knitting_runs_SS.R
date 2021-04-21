###############################################################################
## Project:       Knitting Result for univariate STRS optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   
###############################################################################
rm(list = ls())

###############################
## Set up directories
###############################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[2]

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
  
  n_dom <- ifelse(idom == "full_domain", 1, 5)
  
  frame <- switch( idom,
                   "full_domain" = frame_all,
                   "district" = frame_district)
  
  ###########################
  ## Empty Result objects
  ###########################
  master_res_df <- data.frame(id = 1:n_cells)
  master_settings <- data.frame()
  master_strata_list <- master_strata_stats_list <- list()
  
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
          
          master_settings <- rbind(
            master_settings,
            data.frame(id = id,
                       species = common_names_all[ispp],
                       domain = 1:n_dom,
                       n = tapply(X = result_list$sum_stats$Allocation,
                                  INDEX = result_list$sum_stats$Domain,
                                  FUN = sum),
                       cv = as.numeric(result_list[[3]]) )
          )
          
          #master_res_df: solution (which cell belongs to which stratum?)
          ## Solution: which strata is assigned to each extrapolation cell
          solution <- result_list$sol_by_cell
          
          master_res_df <- cbind(master_res_df,
                                 solution)
          
          #master_strata_list: stratum-level details of solution
          master_strata_list <- c(master_strata_list, 
                                  list(result_list$sum_stats))
          
          #master_strata_stats_list: stratum-level means and variances
          master_strata_stats_list <- c(master_strata_stats_list, 
                                        list(result_list$solution$aggr_strata))
        }
      }
    }
  }
  
  settings_n_by_district <- tidyr::spread(subset(master_settings, 
                                                 select = -cv),
                                          value = n,
                                          key = domain)
  names(settings_n_by_district)[-c(1:2)] <- paste0("n_domain_", 1:n_dom)
  
  settings_total_n <- aggregate(n ~ id + species,
                                data = master_settings,
                                FUN = sum)
  
  settings_cv_by_district <- tidyr::spread(subset(master_settings, 
                                                  select = -n),
                                           value = cv,
                                           key = domain)
  names(settings_cv_by_district)[-c(1:2)] <- paste0("cv_domain_", 1:n_dom)
  
  ####################################
  ## Subset those solutions that correspond to 1, 2, and 3 boats
  ####################################
  sol_idx <- c()
  settings <- data.frame()
  
  for (ispp in 1:ns_all) {
    for (isample in 1:n_boats) {
      temp_idx <- settings_total_n$species == common_names_all[ispp]
      
      if(sum(temp_idx) != 0) {
        #Find solution closet to isample, append to sol_idx
        best_idx <- with(settings_total_n[temp_idx, ], 
                         id[which.min(abs(n - samples[isample]))])
        sol_idx <- c(sol_idx, best_idx)
        
        settings <- rbind(
          settings,
          cbind(
            settings_total_n[best_idx, ],
            boat = isample,
            settings_n_by_district[best_idx, paste0("n_domain_", 1:n_dom)],
            settings_cv_by_district[best_idx, paste0("cv_domain_", 1:n_dom)]))
      }
    }
  } 
  names(settings) <- c("id", "species", "n", "boat", 
                       paste0("n_domain_", 1:n_dom),
                       paste0("cv_domain_", 1:n_dom))
  res_df <- master_res_df[, 1 + sol_idx]
  strata_list <- master_strata_list[sol_idx]
  strata_stats_list <- master_strata_stats_list[sol_idx]
  
  ####################################
  ## Adjust CVs to match sample sizes
  ####################################
  for (irow in 1:nrow(settings)) {
    
    ## Temporary constants
    temp_boat <- settings[irow, "boat"]
    temp_cv <- as.numeric(settings[irow, paste0("cv_domain_", 1:n_dom)])
    temp_n <- settings[irow, "n"]
    
    temp_strata_stats_list <- subset(x = strata_stats_list[[irow]],
                                     select = -SOLUZ)
    
    ## Temporary variables
    temp_n_by_strata <- strata_list[[irow]]$Allocation
    actual_n_by_strata <- strata_stats_list[[irow]]$SOLUZ
    
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
    }
    
    ####################################
    ## Record new results
    ####################################
    temp_cv_agg <- SamplingStrata::expected_CV(
      strata = cbind(subset(temp_strata_stats_list,
                            select = -DOM1),
                     DOM1 = 1,
                     SOLUZ = actual_n_by_strata))
    
    settings[irow, "cv"] <- as.numeric(temp_cv_agg)
    settings[irow, "n"] <- temp_n
    settings[irow, paste0("cv_domain_", 1:n_dom)] <- round(temp_cv, 3)
    settings[irow, paste0("n_domain_", 1:n_dom)] <- temp_n_by_reg
    strata_stats_list[[irow]]$SOLUZ <- actual_n_by_strata 
    strata_list[[irow]]$Allocation <- temp_n_by_strata
  }
  
  settings <- settings[, c("species", "boat", "n", "cv", 
                           paste0("cv_domain_", 1:n_dom),
                           paste0("n_domain_", 1:n_dom))]
  rownames(settings) <- NULL
  
  ###################################
  ## Assign result variables with domain name attached at the end
  ####################################
  vars_to_save <- c("settings", "res_df", "strata_list", "strata_stats_list")
  
  ####################################
  ## Save
  ####################################
  save(list = vars_to_save,
       file = paste0(github_dir, 
                     "results/", idom, 
                     "/Single_Species_Optimization/",
                     "optimization_knitted_results.RData"))
  
}


