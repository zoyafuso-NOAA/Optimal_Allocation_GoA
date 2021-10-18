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
VAST_sim_data_dir <- "E:/VAST_Runs/"
github_dir <- getwd()

##################################################
####   Set up cores
##################################################
library(foreach)
library(parallel)
library(doParallel)

cl <- parallel::makeCluster(parallel::detectCores() - 1)
doParallel::registerDoParallel(cl)

##################################################
####   Load data used by all cores
##################################################
load('data/processed/optimization_data.RData')
load('data/processed/dens_vals_fixed_random.RData')
load('data/processed/dens_vals_measurement.RData')
load('data/processed/dens_vals_MLE.RData')

load('data/processed/grid_goa.RData')
load("data/processed/prednll_VAST_models.RData")
load("data/processed/VAST_fit_I_gct.RData")

##################################
#### Import current STRS design allocations and add to scenarios
##################################
GOA_allocations <- 
  readxl::read_xlsx(path = "data/GOA 2019 stations by stratum.xlsx")
GOA_allocations3 <- 
  readxl::read_xlsx(path = 'data/GOA2019_ 3 boat_825_RNDM_stations.xlsx')

## Specify Management Districts
new_strata_labels = 1:length(unique(grid_goa$stratum))
names(new_strata_labels) <- sort(unique(grid_goa$stratum))

## Rename Current Stratum labels
grid_goa$stratum_new_label <- new_strata_labels[paste(grid_goa$stratum)]

## Create dataframe of effort allocations across boats
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


scenarios <- data.frame(scen_name = c(rep(x = LETTERS[1:11], each = 2), 
                                      c("L", "M")),
                        scale_opt = c(rep(x = "full_domain", times = 2),
                                      rep(x = "district", times = 10 * 2),
                                      rep(x = "current", times = 2)),
                        strata = c(10, 15, 
                                   rep(x = c(3, 5), times = 10),
                                   "current", "current"),
                        opt_type = c(rep(x = "opt", times = 22), 
                                     c("current", "current")),
                        max_depth = c(rep(x = c(rep(x = 1000, 9*2), 
                                                rep(700, 2*2), c(1000, 700)))))

## Create result directories for the existing designs
curr_design_dir <- paste0(github_dir, "/results/scenario_", 
                          scenarios$scen_name[23:24], 
                          "/Multispecies_Optimization/Str_",
                          scenarios$strata[23:24], "/")
for(idir in curr_design_dir) 
  if(!dir.exists(idir)) 
    dir.create(idir, recursive = TRUE)

##################################################
####   Load simulate survey function
##################################################
source("modified_functions/sim_fns.R") 

foreach(iscen = 23:24) %dopar% {
  ##################################################
  ####   Load scenario results
  ##################################################
  if (scenarios$opt_type[iscen] == "opt") 
    load(paste0(github_dir, "/results/scenario_", scenarios$scen_name[iscen], 
                "/Multispecies_Optimization/Str_", scenarios$strata[iscen], 
                "/result_list.RData"))
  
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
  
  cell_idx <- rep(x = TRUE, times = n_cells)
  if(scenarios$max_depth[iscen] == 700) cell_idx[grid_goa$DEPTH_EFH > 700] <- F
  
  temp_solution <- switch(
    scenarios$opt_type[iscen],
    "current" = grid_goa$stratum_new_label[cell_idx],
    "opt" = result_list$sol_by_cell[cell_idx])
  
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
    
    ## Calculate true index based on the survey footprint
    temp_true_index_district <- 
      t(apply(X = I_gct[cell_idx, ispp, ],
              MARGIN = 2,
              FUN = function(x) tapply(X = x, 
                                       INDEX = district_vals[cell_idx],
                                       FUN = sum)))
    for (iter in 1:n_iters) { ## loop over survey replicates--start
      set.seed(n_iters + iter)
      for (iboat in 1:n_boats) { ## loop over boat-effort--start
        
        temp_allocation <- switch(
          scenarios$opt_type[iscen],
          "current" = allocations[, paste0("boat", iboat)],
          "opt" = result_list$sample_allocations[, paste0("boat_", iboat)])
        
        ## Temporary input objects for the survey simulation
        temp_density <- sweep(x = sim_data[cell_idx, , iter] / 0.001,
                              MARGIN = 1,
                              STATS = grid_goa$Area_km2[cell_idx],
                              FUN = "/")
        
        ## Subset the iter-th simulated dataset and turn into densities
        sim_survey <- do_STRS(input = list(
          "density" = temp_density,
          "cell_areas" = grid_goa$Area_km2[cell_idx],
          "solution" = temp_solution,
          "allocation" = temp_allocation,
          "true_index_district" = temp_true_index_district,
          "post_strata" = district_vals[cell_idx]
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
        for (iboat in 1:n_boats) { ## subloop over boat-effort--start
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
            update_file  <- file(paste0(github_dir, "/results/scenario_",
                                        scenarios$scen_name[iscen], 
                                        "/Multispecies_Optimization/Str_",
                                        scenarios$strata[iscen], 
                                        "/progress.txt"))
            writeLines(c(paste0("Finished: Iteration ", iter, 
                                ", Species: ", ispp)), update_file)
            close(update_file)
            
          } ## subloop over years--end
        } ## subloop over boat-effort--end
        
        
        ##################################
        ## Save results
        ##################################
        save(list = c("STRS_sim_index", "STRS_sim_cv",
                      "STRS_rel_bias_est", "STRS_rel_bias_index_district",
                      "STRS_log_bias_index_district",
                      "STRS_true_cv_array", "STRS_rrmse_cv_array" ),
             file = paste0(github_dir, "/results/scenario_", 
                           scenarios$scen_name[iscen], 
                           "/Multispecies_Optimization/Str_",
                           scenarios$strata[iscen], 
                           "/simulation_results.RData"))
        
      } ## update results--end
    } ## loop over survey replicates--end
  } ## loop over species--end  
}
