###############################################################################
## Project:       Simulate Surveys
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Simulate a Stratified Random Survey of the Gulf of 
##                Alaska Groundfish Survey Based on Current Stratifications
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(readxl)
library(sp)

##################################################
####   Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[3]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/results/")
VAST_dir <- "G:/Oyafuso/VAST_Runs_EFH/Single_Species/"

##################################################
####   Load simulation functions
##################################################
source( paste0(dirname(github_dir), "/modified_functions/sim_fns.R") )

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load(paste0(dirname(github_dir), '/data/optimization_data.RData'))
load(paste0(dirname(github_dir), '/data/RMSE_VAST_models.RData'))
load(paste0(dirname(github_dir), '/data/Extrapolation_depths.RData'))
load(paste0(dirname(github_dir), '/data/fit_density.RData'))
load(paste0(github_dir, "Spatiotemporal_Optimization",
            "/optimization_knitted_results.RData"))

GOA_allocations <- readxl::read_xlsx(
  path = paste0(dirname(github_dir), 
                '/data/GOA 2019 stations by stratum.xlsx')
)
GOA_allocations3 <- readxl::read_xlsx(
  path = paste0(dirname(github_dir), 
                '/data/GOA2019_ 3 boat_825_RNDM_stations.xlsx')
) 

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

##################################################
####   Subset 15 strata solutions
##################################################
sol_idx <- which(settings$strata == 15)
settings <- settings[sol_idx,]
res_df <- res_df[, 1 + sol_idx]
strata_list <- strata_list[sol_idx]

##################################################
####   Result Objects
##################################################
Current_sim_mean <- Current_sim_cv <- Current_rel_bias_est <- 
  STRS_sim_mean <- STRS_sim_cv <- STRS_rel_bias_est <-
  array(dim = c(nobs_CV, NTime, ns_all, nboats, Niters), 
        dimnames = list(paste0("obsCV=", obs_CV), 
                        paste0("year_", 1:NTime), 
                        sci_names_all, 
                        paste0("boat_", 1:nboats), 
                        NULL ))

Current_true_cv_array <- Current_rrmse_cv_array <- 
  STRS_true_cv_array <- STRS_rrmse_cv_array <-  
  array(dim = c(nobs_CV, NTime, ns_all, nboats), 
        dimnames = list(paste0("obsCV=", obs_CV), 
                        paste0("year_", 1:NTime), 
                        sci_names_all, 
                        paste0("boat_", 1:nboats)))

##################################################
####   Simulate Survey
##################################################

for (iter in 1:Niters) {
  
  set.seed(1000 + iter)
  
  for (ierror in 1:nobs_CV) {
    for (iboat in 1:3) {
      for (isurvey in c("Current", "STRS")) {
        
        sim_survey <- 
          do_STRS( input = list(
            "density" = D_gct[, , Years2Include],
            
            "obs_CV" = obs_CV[ierror],
            
            "solution" = switch(
              isurvey,
              "Current" = Extrapolation_depths$stratum,
              "STRS" = res_df[, iboat]),
            
            "allocation" = switch( 
              isurvey,
              "Current" = allocations[, paste0("boat", iboat)],
              "STRS" = strata_list[[iboat]]$Allocation),
            
            "true_density" = true_mean) )
        
        stmt <- paste0(isurvey, "_sim_mean",  
                       "[ierror, , , iboat, iter] = sim_survey$mean_denisty")
        eval(parse(text = stmt))
        
        stmt <- paste0(isurvey, "_sim_cv",  
                       "[ierror, , , iboat, iter] = sim_survey$cv")
        eval(parse(text = stmt))
        
        stmt <- paste0(isurvey, "_rel_bias_est",  
                       "[ierror, , , iboat, iter] = sim_survey$rel_bias")
        eval(parse(text = stmt))
      }
    }
  }
  if(iter%%10 == 0) print(paste("Finished with Iteration", iter))
}

##################################
## Calculate Performance Metric
##################################
for (ierror in 1:nobs_CV) {
  for (iboat in 1:3) {
    for (isurvey in c("Current", "STRS")) {
      for (ispp in 1:ns_all) {
        for(iyear in 1:NTime) {
          
          stmt <- paste0(isurvey, "_true_cv_array[ierror, iyear, ispp, iboat]",
                         " <- temp_true_cv <- sd(", isurvey, 
                         "_sim_mean[ierror, iyear, ispp, iboat,], na.rm = T) / ",
                         "true_mean[ispp, iyear]")
          eval(parse(text = stmt))
          
          temp_sim_cv <- get(paste0(isurvey, 
                                    "_sim_cv"))[ierror, iyear, ispp, iboat,]
          
          stmt <- paste0(isurvey, "_rrmse_cv_array[ierror, iyear, ispp, iboat]",
                         " <- sqrt(mean((temp_sim_cv - temp_true_cv)^2, na.rm = T)) / ", 
                         "mean(temp_sim_cv, na.rm = T)")
          eval(parse(text = stmt))
          
        }
      }
    }
  }
}

##################################
## Save results
##################################
save(list = c("Current_sim_mean", "STRS_sim_mean", "Current_sim_cv",        
              "STRS_sim_cv", "Current_rel_bias_est", "STRS_rel_bias_est",     
              "Current_true_cv_array", "STRS_true_cv_array", 
              "Current_rrmse_cv_array", "STRS_rrmse_cv_array"),
     file = paste0(github_dir, "simulation_result.RData"))
