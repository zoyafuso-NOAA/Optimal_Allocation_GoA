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
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/results/")

VAST_sim_data_dir <- "C:/Users/Zack Oyafuso/Desktop/VAST_Runs/Simulate_Data/"

##################################################
####   Result directories
##################################################
for (idata in c("pred_dens", "sim_dens")) {
  if(!dir.exists(paste0(github_dir, idata, "_surveys/"))) {
    dir.create(paste0(github_dir, idata, "_surveys/"))
  }
}
rm(idata)

##################################################
####   Set up cores
##################################################
library(foreach)
library(parallel)
library(doParallel)

cl <- parallel::makeCluster(parallel::detectCores() - 1)
doParallel::registerDoParallel(cl)

##################################################
####   Optimization scenarios
##################################################
scen <- data.frame(survey_type = c("cur", rep("opt", 4) ),
                   strata = c("cur", 3, 5, 10, 15),
                   domain = c("full_domain", 
                              rep(c("district", "full_domain"), each = 2)))

foreach(irow = nrow(scen):1 ) %dopar% {
  
  ##################################################
  ####   Import Libraries
  ##################################################
  library(readxl)
  library(sp)
  
  ##################################################
  ####   Predicted densities
  ##################################################
  load(paste0(dirname(github_dir), "/data/fit_density.RData"))
  dimnames(D_gct)[[2]][24] <- "Pacific spiny dogfish"
  
  ##################################
  ## Import Strata Allocations and spatial grid and predicted density
  ##################################
  load(paste0(dirname(github_dir), '/data/optimization_data.RData'), verbose = T)
  load(paste0(dirname(github_dir), '/data/Extrapolation_depths.RData'))
  load(paste0(github_dir, "MS_optimization_knitted_results.RData"))
  
  GOA_allocations <- readxl::read_xlsx(
    path = paste0(dirname(github_dir), 
                  '/data/GOA 2019 stations by stratum.xlsx')
  )
  GOA_allocations3 <- readxl::read_xlsx(
    path = paste0(dirname(github_dir), 
                  '/data/GOA2019_ 3 boat_825_RNDM_stations.xlsx')
  ) 
  
  ##################################
  ## Specify Management Districts
  ##################################
  new_strata_labels = 1:length(unique(Extrapolation_depths$stratum))
  names(new_strata_labels) <- sort(unique(Extrapolation_depths$stratum))
  
  ##################################
  ## Rename Current Stratum labels
  ##################################
  Extrapolation_depths$stratum_new_label <- 
    new_strata_labels[paste(Extrapolation_depths$stratum)]
  
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
  
  n_dom = 5
  
  ##################################################
  ####   Load simulation functions
  ##################################################
  source( paste0(dirname(github_dir), "/modified_functions/sim_fns.R") )
  
  for(idata in c("pred_dens_surveys", "sim_dens_surveys")[]) {
    isurvey <- scen$survey_type[irow]
    
    scen_name <- paste0("SUR_", isurvey, "_", scen$domain[irow], "_STR_", 
                        scen$strata[irow], "_")
    
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
      array(dim = c(n_years, ns_all, n_boats, n_dom, n_iters),
            dimnames = list(paste0("year_", 1:n_years),
                            common_names_all,
                            paste0("boat_", 1:n_boats),
                            paste0("district_", 1:n_dom),
                            NULL ))
    
    STRS_true_cv_array <- STRS_rrmse_cv_array <-
      array(dim = c(n_years, ns_all, n_boats),
            dimnames = list(paste0("year_", 1:n_years),
                            common_names_all,
                            paste0("boat_", 1:n_boats)))
    
    ##################################################
    ####   Simulate Survey
    ##################################################
    for (ispp in common_names_all) {
      
      # ## Load species specific simulated data
      if(idata ==  "sim_dens_surveys")
        load(paste0(VAST_sim_data_dir, ispp, "/simulated_data.RData"))
      
      
      for (iter in 1:n_iters) {
        set.seed(n_iters + iter)
        for (iboat in 1:n_boats) {
          
          ## Subset the iter-th simulated dataset and turn into densities
          sim_survey <- do_STRS(input = list(
            "density" = switch(idata,
                               "sim_dens_surveys" = 
                                 sweep(x = sim_data[, , 1],
                                       MARGIN = 1,
                                       STATS = Extrapolation_depths$Area_km2,
                                       FUN = "/"),
                               "pred_dens_surveys" = D_gct[1:n_cells,
                                                           ispp,
                                                           1:n_years] ),
            "cell_areas" = Extrapolation_depths$Area_km2,
            
            "solution" = switch(
              isurvey,
              "cur" = Extrapolation_depths$stratum_new_label,
              "opt" = {
                idx <- which(settings$strata == as.numeric(scen$strata[irow]) &
                               settings$boat == iboat &
                               settings$domain == scen$domain[irow])
                
                res_df[, idx]
              }),
            
            "allocation" = switch(
              isurvey,
              "cur" = allocations[, paste0("boat", iboat)],
              "opt" = {
                idx <- which(settings$strata == as.numeric(scen$strata[irow]) &
                               settings$boat == iboat &
                               settings$domain == scen$domain[irow])
                
                strata_list[[idx]]$Allocation
              }),
            
            "true_density" = true_mean[ispp, ],
            
            # "true_index_district" = true_index_district[ispp, , ],
            "true_index_district" = switch(
              idata,
              "sim_dens_surveys" = 
                t(apply(sim_data[, , 1], 
                        MARGIN = 2,
                        FUN = function(x) 0.001 * tapply(X = x, 
                                                         INDEX = district_vals,
                                                         FUN = sum))),
              "pred_dens_surveys" = true_index_district[ispp, , ]),
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
        }
        
        
        if(iter%%50 == 0) {
          ##################################
          ## Calculate Performance Metric
          ##################################
          for (iboat in 1:n_boats) {
            for(iyear in 1:n_years) {
              
              ## Calculate True CV
              STRS_true_cv_array[iyear, ispp, iboat] <- temp_true_cv <-
                sd(STRS_sim_index[iyear, ispp, iboat, ], na.rm = T) /
                switch(idata,
                       "pred_dens_surveys" = true_index[ispp, iyear],
                       "sim_dens_surveys" = 0.001 * sum(sim_data[, iyear, 1]) )
              
              
              temp_sim_cv <- STRS_sim_cv[iyear, ispp, iboat, ]
              
              ## Calculate RRMSE of CV
              STRS_rrmse_cv_array[iyear, ispp, iboat] <-
                sqrt(mean((temp_sim_cv - temp_true_cv)^2, na.rm = T)) /
                mean(temp_sim_cv, na.rm = T)
              
              #Update progress to file
              update_file  <- file(paste0(github_dir, idata, "/",
                                          scen_name, "progress.txt"))
              writeLines(c(scen_name,
                           paste0("Finished: Iteration ", iter, ", ",
                                  ispp, ", ",
                                  isurvey, " scenario, ",
                                  scen$strata[irow], " strata")), update_file)
              close(update_file)
              
            }
          }
          
          
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
               file = paste0(github_dir,  idata, "/", scen_name,
                             "simulation_results.RData"))
        }
      }
      
      
    }
    
  }
}
