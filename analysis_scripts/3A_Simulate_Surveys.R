###############################################################################
## Project:       Simulate Surveys
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Simulate a Stratified Random Survey of the Gulf of 
##                Alaska Groundfish Survey Based on Current Stratifications
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
VAST_sim_data_dir <- "D:/VAST_Runs/"
github_dir <- getwd()

##################################################
####   Result directories
##################################################
if(!dir.exists("results/survey_simulations/"))
  dir.create("results/survey_simulations/")

##################################################
####   Set up cores
##################################################
library(foreach)
library(parallel)
library(doParallel)

cl <- parallel::makeCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)

##################################################
####   Optimization scenarios
##################################################
scen <- data.frame(survey_type = c("cur", rep("opt", 4) ),
                   strata = c("cur", 3, 5, 10, 15),
                   domain = c("full_domain", 
                              rep(c("district", "full_domain"), each = 2)))

foreach(irow = nrow(scen):1 ) %dopar% { ## loop over scenario type -- start
  
  isurvey <- scen$survey_type[irow]
  
  ##################################################
  ####   Import Libraries
  ##################################################
  library(readxl)
  library(sp)
  
  ##################################
  ## Import Strata Allocations and spatial grid and predicted density
  ##################################
  load('data/processed/optimization_data.RData')
  load('data/processed/grid_goa.RData')
  load("results/MS_optimization_knitted_results.RData")
  load("data/processed/prednll_VAST_models.RData")
  
  GOA_allocations <- 
    readxl::read_xlsx(path = "data/GOA 2019 stations by stratum.xlsx")
  GOA_allocations3 <- 
    readxl::read_xlsx(path = 'data/GOA2019_ 3 boat_825_RNDM_stations.xlsx')
  
  ##################################
  ## Specify Management Districts
  ##################################
  new_strata_labels = 1:length(unique(grid_goa$stratum))
  names(new_strata_labels) <- sort(unique(grid_goa$stratum))
  
  ##################################
  ## Rename Current Stratum labels
  ##################################
  grid_goa$stratum_new_label <- new_strata_labels[paste(grid_goa$stratum)]
  
  ##################################################
  ####   Create dataframe of effort allocations across boats
  ##################################################
  allocations <- data.frame(Stratum = sort(unique(GOA_allocations3$stratum)),
                            boat3 = aggregate(id ~ stratum, 
                                              data = GOA_allocations3, 
                                              FUN = length)$id,
                            boat2 = c(GOA_allocations$`Number stations`,
                                      rep(0, 5)))
  allocations$boat1 <- ceiling(allocations$boat2 / 2)
  
  allocations$boat1 <- ifelse(allocations$boat1 == 0, 0, 
                              ifelse(allocations$boat1 == 1, 2, 
                                     allocations$boat1))
  
  allocations <- rbind(data.frame(Stratum = 0, boat3 = 0, boat2 = 0, boat1 = 0),
                       allocations)
  
  allocations$Stratum <- 1:nrow(allocations)
  
  ##################################################
  ####   Load simulate survey function
  ##################################################
  source("modified_functions/sim_fns.R") 
  
  ##################################################
  ####   
  ##################################################
  scen_name <- paste0(scen$domain[irow], "_STR_", scen$strata[irow], "_")
  
  ##################################################
  ####   Result Objects
  ##################################################
  STRS_sim_index <- STRS_sim_cv <- STRS_rel_bias_est <- 
    array(dim = c(n_years, ns_all, n_boats, n_iters),
          dimnames = list(paste0("year_", 1:n_years),
                          common_names_all,
                          paste0("boat_", 1:n_boats),
                          NULL ))
  
  STRS_rel_bias_index_district <- STRS_log_bias_index_district <-
    STRS_index_district <-
    array(dim = c(n_years, ns_all, n_boats, n_districts, n_iters),
          dimnames = list(paste0("year_", 1:n_years),
                          common_names_all,
                          paste0("boat_", 1:n_boats),
                          paste0("district_", 1:n_districts),
                          NULL ))
  
  STRS_true_cv_array <- STRS_rrmse_cv_array <-
    array(dim = c(n_years, ns_all, n_boats),
          dimnames = list(paste0("year_", 1:n_years),
                          common_names_all,
                          paste0("boat_", 1:n_boats)))
  
  ##################################################
  ####   Simulate Survey
  ##################################################
  for (ispp in common_names_all) { ## loop over species--start
    
    ## Load species specific simulated data
    depth_in_model <- c(F, T)[which.min(pred_jnll[pred_jnll$spp_name == ispp, 
                                                  2:3])]
    if(ispp == "shortraker rockfish") depth_in_model <- FALSE
    
    load(paste0(VAST_sim_data_dir, ispp, 
                ifelse(test = depth_in_model, yes = "_depth", no = ""),
                "/simulated_data.RData"))
    
    for (iter in 1:n_iters) { ## loop over survey replicates--start
      set.seed(n_iters + iter)
      for (iboat in 2) { ## loop over boat-effort--start
        
        ## Temporary input objects for the survey simulation
        temp_density <- sweep(x = sim_data[, , iter] / 0.001,
                              MARGIN = 1,
                              STATS = grid_goa$Area_km2,
                              FUN = "/")
        
        temp_solution <- switch(
          isurvey,
          "cur" = grid_goa$stratum_new_label,
          "opt" = {
            idx <- which(settings$strata == as.numeric(scen$strata[irow]) &
                           settings$boat == iboat &
                           settings$domain == scen$domain[irow])
            res_df[, idx]
          })
        
        temp_allocation <- switch(
          isurvey,
          "cur" = allocations[, paste0("boat", iboat)],
          "opt" = {
            idx <- which(settings$strata == as.numeric(scen$strata[irow]) &
                           settings$boat == iboat &
                           settings$domain == scen$domain[irow])
            
            strata_list[[idx]]$Allocation
          })
        
        ## Subset the iter-th simulated dataset and turn into densities
        sim_survey <- do_STRS(input = list(
          "density" = temp_density,
          "cell_areas" = grid_goa$Area_km2,
          "solution" = temp_solution,
          "allocation" = temp_allocation,
          "true_index_district" = true_index_district[ispp, , ],
          "post_strata" = district_vals
        ))
        
        ## Record results
        STRS_sim_index[, ispp, iboat, iter] = sim_survey$strs_index
        STRS_sim_cv[, ispp, iboat, iter] = sim_survey$cv
        STRS_rel_bias_est[, ispp, iboat, iter] = sim_survey$rel_bias
        
        STRS_rel_bias_index_district[, ispp, iboat, , iter] <-
          sim_survey$bias_index_district
        STRS_log_bias_index_district[, ispp, iboat, , iter]  <-
          sim_survey$log_bias_index_district
        
      }  ## loop over boat-effort--end
      
      if(iter%%100 == 0) { ## update results--start
        
        print(paste0("Finished: Iteration ", iter))
        
        ##################################
        ## Calculate Performance Metrics
        ##################################
        for (iboat in 2) { ## subloop over boat-effort--start
          for(iyear in 1:n_years) { ## subloop over years--start
            
            ## Calculate True CV
            STRS_true_cv_array[iyear, ispp, iboat] <- temp_true_cv <-
              sd(STRS_sim_index[iyear, ispp, iboat, ], na.rm = T) / 
              (true_index[ispp, iyear])
            
            temp_sim_cv <- STRS_sim_cv[iyear, ispp, iboat, ]
            
            ## Calculate RRMSE of CV
            STRS_rrmse_cv_array[iyear, ispp, iboat] <-
              sqrt(mean((temp_sim_cv - temp_true_cv)^2, na.rm = T)) /
              mean(temp_sim_cv, na.rm = T)
            
            #Update progress to file
            update_file  <- file(paste0(github_dir, 
                                        "/results/survey_simulations/",
                                        scen_name, "progress.txt"))
            writeLines(c(scen_name,
                         paste0("Finished: Iteration ", iter, ", ",
                                ispp, ", ", " scenario, ",
                                scen$strata[irow], " strata")), update_file)
            close(update_file)
            
          } ## subloop over years--end
        } ## subloop over boat-effort--end
        
        
        ##################################
        ## Save results
        ##################################
        assign(value = STRS_sim_index,
               x = paste0(scen_name, "sim_index"))
        assign(value = STRS_sim_cv,
               x = paste0(scen_name, "sim_cv"))
        
        assign(value = STRS_rel_bias_est,
               x = paste0(scen_name, "rb_agg"))
        
        assign(value = STRS_rel_bias_index_district,
               x = paste0(scen_name, "rb_district"))
        assign(value = STRS_log_bias_index_district,
               x = paste0(scen_name, "log_rb_district"))
        
        assign(value = STRS_true_cv_array,
               x = paste0(scen_name, "true_cv"))
        assign(value = STRS_rrmse_cv_array,
               x = paste0(scen_name, "rrmse_cv"))
        
        save(list = paste0(scen_name,
                           c("sim_index", "sim_cv", "rb_agg",
                             "true_cv", "rrmse_cv",
                             "rb_district", "log_rb_district") ),
             file = paste0(github_dir, "/results/survey_simulations/",
                           scen_name, "simulation_results.RData"))
        
      } ## update results--end
    } ## loop over survey replicates--end
  } ## loop over species--end
  ## loop over data OM type -- end
} ## loop over scenario type -- end
