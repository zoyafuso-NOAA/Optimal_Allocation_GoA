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
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[3]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/results/")

##################################################
####   Import Libraries
##################################################
library(readxl)
library(sp)

##################################################
####   Load simulation functions
##################################################
source( paste0(dirname(github_dir), "/modified_functions/sim_fns.R") )

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load(paste0(dirname(github_dir), '/data/optimization_data.RData'))
load(paste0(dirname(github_dir), '/data/Extrapolation_depths.RData'))
load(paste0(dirname(github_dir), '/data/fit_density.RData'))
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

scen <- data.frame(survey_type = c("cur", rep("opt", 6) ),
                   strata = c("cur", 3, 5, 10, 10, 15, 20),
                   domain = c("full_domain", 
                              rep(c("district", "full_domain"), each = 3)))

n_iters = 500

for (irow in 7) {
  ##################################################
  ####   Result Objects
  ##################################################
  STRS_sim_mean <- STRS_sim_cv <- STRS_rel_bias_est <- STRS_log_bias_est <-
    array(dim = c(n_obs_cv, n_years, ns_all, n_boats, n_iters), 
          dimnames = list(paste0("obsCV=", obs_cv), 
                          paste0("year_", 1:n_years), 
                          sci_names_all, 
                          paste0("boat_", 1:n_boats), 
                          NULL ))
  
  STRS_rel_bias_index_district <- STRS_log_bias_index_district <- 
    STRS_index_district <- 
    array(dim = c(n_obs_cv, n_years, ns_all, n_boats, n_dom, n_iters),
          dimnames = list(paste0("obsCV=", obs_cv),
                          paste0("year_", 1:n_years),
                          sci_names_all,
                          paste0("boat_", 1:n_boats),
                          paste0("district_", 1:n_dom),
                          NULL ))
  
  STRS_true_cv_array <- STRS_rrmse_cv_array <-  
    array(dim = c(n_obs_cv, n_years, ns_all, n_boats), 
          dimnames = list(paste0("obsCV=", obs_cv), 
                          paste0("year_", 1:n_years), 
                          sci_names_all, 
                          paste0("boat_", 1:n_boats)))
  
  ##################################################
  ####   Simulate Survey
  ##################################################
  
  isurvey <- scen$survey_type[irow]
  
  for (iter in 1:n_iters) {
    # for (iter in 1:200) {
    set.seed(1000 + iter)
    for (ierror in 1:n_obs_cv) {
      for (iboat in 1:3) {
        
        sim_survey <- 
          do_STRS( 
            input = list(
              "density" = D_gct[, , years_included],
              
              "cell_areas" = Extrapolation_depths$Area_km2,
              
              "obs_CV" = obs_cv[ierror],
              
              "solution" = switch(
                isurvey,
                "cur" = Extrapolation_depths$stratum_new_label,
                "opt" = {
                  idx <- subset(
                    x = settings, 
                    subset = strata == as.numeric(scen$strata[irow]) &
                      boat == iboat &
                      domain == scen$domain[irow])$id
                  
                  res_df[, paste0("sol_", idx)]
                }),
              
              "allocation" = switch(
                isurvey,
                "cur" = allocations[, paste0("boat", iboat)],
                "opt" = {
                  idx <- subset(
                    x = settings, 
                    subset = strata == as.numeric(scen$strata[irow]) &
                      boat == iboat &
                      domain == scen$domain[irow])$id
                  
                  strata_list[[paste0("sol_", idx)]]$Allocation
                }),
              
              "true_density" = true_mean,
              
              "true_index_district" = true_index_district,
              
              "post_strata" = district_vals
            )
          )
        
        STRS_sim_mean[ierror, , , iboat, iter] = sim_survey$mean_denisty
        STRS_sim_cv[ierror, , , iboat, iter] = sim_survey$cv
        STRS_rel_bias_est[ierror, , , iboat, iter] = sim_survey$rel_bias
        STRS_log_bias_est[ierror, , , iboat, iter] = sim_survey$rel_log_bias
        STRS_rel_bias_index_district[ierror, , , iboat, , iter] <- 
          sim_survey$bias_index_district
        STRS_log_bias_index_district[ierror, , , iboat, , iter] <- 
          sim_survey$log_bias_index_district
      }
    }
    
    if(iter%%10 == 0) {
      ##################################
      ## Calculate Performance Metric
      ##################################
      for (ierror in 1:n_obs_cv) {
        for (iboat in 1:n_boats) {
          for (ispp in 1:ns_all) {
            for(iyear in 1:n_years) {
              
              ## Calculate True CV
              STRS_true_cv_array[ierror, iyear, ispp, iboat] <- temp_true_cv <- 
                sd(STRS_sim_mean[ierror, iyear, ispp, iboat, ], na.rm = T) / 
                true_mean[ispp, iyear]
              
              temp_sim_cv <- STRS_sim_cv[ierror, iyear, ispp, iboat, ]
              
              ## Calculate RRMSE of CV
              STRS_rrmse_cv_array[ierror, iyear, ispp, iboat] <-
                sqrt(mean((temp_sim_cv - temp_true_cv)^2, na.rm = T)) /  
                mean(temp_sim_cv, na.rm = T)
            }
          }
        }
      }
      
      ##################################
      ## Save results
      ##################################
      assign(value = STRS_sim_mean,
             x = paste0("SUR_", isurvey, "_", scen$domain[irow], "_STR_", 
                        scen$strata[irow], "_sim_mean"))
      assign(value = STRS_sim_cv,
             x = paste0("SUR_", isurvey, "_", scen$domain[irow], "_STR_", 
                        scen$strata[irow], "_sim_cv"))
      
      assign(value = STRS_rel_bias_est,
             x = paste0("SUR_", isurvey, "_", scen$domain[irow], "_STR_", 
                        scen$strata[irow], "_rb_agg"))
      assign(value = STRS_log_bias_est,
             x = paste0("SUR_", isurvey, "_", scen$domain[irow], "_STR_", 
                        scen$strata[irow], "_log_rb_agg"))
      
      assign(value = STRS_rel_bias_index_district,
             x = paste0("SUR_", isurvey, "_", scen$domain[irow], "_STR_", 
                        scen$strata[irow], "_rb_district"))
      assign(value = STRS_log_bias_index_district,
             x = paste0("SUR_", isurvey, "_", scen$domain[irow], "_STR_", 
                        scen$strata[irow], "_log_rb_district"))
      
      assign(value = STRS_true_cv_array,
             x = paste0("SUR_", isurvey, "_", scen$domain[irow], "_STR_", 
                        scen$strata[irow], "_true_cv"))
      assign(value = STRS_rrmse_cv_array,
             x = paste0("SUR_", isurvey, "_", scen$domain[irow], "_STR_", 
                        scen$strata[irow], "_rrmse_cv"))
      
      save(list = paste0(paste0("SUR_", isurvey, "_", scen$domain[irow],
                                "_STR_", scen$strata[irow]), 
                         c("_sim_mean", "_sim_cv", 
                           "_rb_agg", "_log_rb_agg",
                           "_rb_district", "_log_rb_district", 
                           "_true_cv", "_rrmse_cv") ),
           file = paste0(github_dir, "SUR_", isurvey, "_", 
                         scen$domain[irow], "_STR_", scen$strata[irow], 
                         "_simulation_results.RData"))
      
      print(paste("Finished with Iteration", iter))
    }
  }
}

