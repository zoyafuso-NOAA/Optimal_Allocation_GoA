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
library(spatialEco)
library(sp)
library(VAST)

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
load(paste0(dirname(VAST_dir), "/sim_density.RData"))

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
                                    rep(0,5)))
allocations$boat1 <- ceiling(allocations$boat2 / 2)

allocations$boat1 <- ifelse(allocations$boat1 == 0, 0, 
                            ifelse(allocations$boat1 == 1, 2, 
                                   allocations$boat1))

allocations <- rbind(data.frame(Stratum = 0, boat3 = 0, boat2 = 0, boat1 = 0),
                     allocations)

##################################################
####   Result Objects
##################################################
Current_sim_mean <- Current_sim_cv <- Current_rel_bias_est <- 
  STRS_sim_mean <- STRS_sim_cv <- STRS_rel_bias_est <-
  array(dim = c(2, NTime, ns, nboats, Niters), 
        dimnames = list(c("pred_density", "plus_obs_error"), 
                        paste0("year_", 1:NTime), 
                        sci_names, 
                        paste0("boat_", 1:3), 
                        NULL ))

Current_true_cv_array <- Current_rrmse_cv_array <- 
  STRS_true_cv_array <- STRS_rrmse_cv_array <-  
  array(dim = c(2, NTime, ns, nboats), 
        dimnames = list(c("pred_density", "plus_obs_error"), 
                        paste0("year_", 1:NTime), 
                        sci_names, 
                        paste0("boat_", 1:3)))

##################################################
####   Simulate Survey
##################################################

for (iter in 1:1000) {
  
  set.seed(1000 + iter)
  
  for (isim in c("pred_density", "plus_obs_error")) {
    truth <- true_mean
      
      # switch(
      #   isim,
      #   "pred_density" = true_mean,
      #   "plus_obs_error" = t(apply(X = pred_density$plus_obs_error[ceiling(iter / 100), , , ],
      #                              MARGIN = 2:3, 
      #                              FUN = mean)))
    
    for (iboat in 1:3) {
      for (isurvey in c("Current", "STRS")) { #Current or Optimized Survey
        
        if(isurvey == "STRS") {
          #Load optimization data, only focusing on 15 strata for now
          sub_settings = subset(settings, strata == 15)
          idx <- which.min(abs(sub_settings$n - samples[iboat]))
        }
        
        sim_survey <- 
          do_STRS(
            input <- list(
              "density" = switch(
                isim,
                "pred_density" = pred_density$pred_density,
                "plus_obs_error" = pred_density$plus_obs_error[ceiling(iter / 100), , , ]),
              
              "solution" = switch(
                isurvey,
                "Current" = Extrapolation_depths$stratum,
                "STRS" = res_df[, 1 + idx]),
              
              "allocation" = switch( 
                isurvey,
                "Current" = allocations[, paste0("boat", iboat)],
                "STRS" = strata_list[[idx]]$Allocation),
              
              "true_density" = truth )
          )
        
        stmt <- paste0(isurvey, "_sim_mean",  
                       "[isim, , , iboat, iter] = sim_survey$mean_denisty")
        eval(parse(text = stmt))
        
        stmt <- paste0(isurvey, "_sim_cv",  
                       "[isim, , , iboat, iter] = sim_survey$cv")
        eval(parse(text = stmt))
        
        stmt <- paste0(isurvey, "_rel_bias_est",  
                       "[isim, , , iboat, iter] = sim_survey$rel_bias")
        eval(parse(text = stmt))
      }
    }
  }
  if(iter%%10 == 0) print(paste("Finished with Iteration", iter))
}

##################################
## Calculate Performance Metric
##################################
for (isim in c("pred_density", "plus_obs_error")) {
  truth <- switch(
    isim,
    "pred_density" = true_mean,
    "plus_obs_error" = t(apply(X = pred_density$plus_obs_error[ceiling(1 / 100), , , ],
                               MARGIN = 2:3, 
                               FUN = mean)))
  
  for (iboat in 1:3) {
    for (isurvey in c("Current", "STRS")) {
      for (ispp in 1:ns) {
        for(iyear in 1:NTime) {
          
          stmt <- paste0(isurvey, "_true_cv_array[isim, iyear, ispp, iboat]",
                         " <- temp_true_cv <- sd(", 
                         isurvey, "_sim_mean[isim, iyear, ispp, iboat,]) / ",
                         "truth[iyear, ispp]")
          eval(parse(text = stmt))
        
          temp_sim_cv <- get(paste0(isurvey, 
                                    "_sim_cv"))[isim, iyear, ispp, iboat,]
          
          stmt <- paste0(isurvey, "_rrmse_cv_array[isim, iyear, ispp, iboat]",
                         " <- sqrt(mean((temp_sim_cv - temp_true_cv)^2)) / ", 
                         "mean(temp_sim_cv)")
          eval(parse(text = stmt))
            
        }
      }
    }
  }
}


save(list = c("Current_sim_mean", "STRS_sim_mean", "Current_sim_cv",        
              "STRS_sim_cv", "Current_rel_bias_est", "STRS_rel_bias_est",     
              "Current_true_cv_array", "STRS_true_cv_array", 
              "Current_rrmse_cv_array", "STRS_rrmse_cv_array"),
     file = paste0(github_dir, "simulation_result.RData"))