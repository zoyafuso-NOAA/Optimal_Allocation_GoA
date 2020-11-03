###############################################################################
## Project:         Simulate Multispecies GoA Surveys
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     Simulate surveys based on the optimized survey designs
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[3]
VAST_model <- "11" 
which_domain <- c("full_domain", "trawlable")[1]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/",
                     which_domain, "/")

if(!dir.exists(paste0(github_dir, "Survey_Comparison_Simulations/"))) {
  dir.create(paste0(github_dir, "Survey_Comparison_Simulations/"))
}

result_dir <- paste0(github_dir, "Spatiotemporal_Optimization/")

##################################################
####   Load simulation functions
##################################################
source( paste0(dirname(dirname(github_dir)), "/modified_functions/sim_fns.R") )

##################################################
####    Load predicted density and optimization results
##################################################
# load(paste0(dirname(dirname(github_dir)), '/data/Extrapolation_depths.RData'))
load(paste0(github_dir, "optimization_data.RData"))
load(paste0(result_dir, "optimization_knitted_results.RData"))

##################################################
####   Result Objects
##################################################
# stratas = c(5, 10, 15, 20, 30, 60)
# NStratas = length(stratas)

STRS_sim_mean <- STRS_sim_cv <- STRS_rel_bias_est <-
  array(dim = c(NTime, ns, nboats, NStrata, Niters), 
        dimnames = list(NULL, sci_names, NULL))

STRS_true_cv_array <- STRS_rrmse_cv_array <-  
  array(dim = c(NTime, ns, nboats, NStrata), 
        dimnames = list(NULL, sci_names, NULL, NULL ))

##################################################
####   Simulating surveys from each optimized solution
##################################################
# for (istrata in c(1:NStratas)) {
for (istrata in 3) {
  for (iboat in 2) {
    
    #Load optimization data
    sub_settings = subset(settings, strata == stratas[istrata])
    
    idx <- which.min(abs(sub_settings$n - samples[iboat]))
    
    
    for (iter in 1:Niters) {
      
      temp_sim <-  do_STRS(density = frame_raw[, paste0("Y", 1:ns)],
                           cell_idx = (1:N),
                           strata = res_df[, 1 + idx],
                           strata_to_use = 1:nrow(strata_list[[idx]]), 
                           allocation = strata_list[[idx]]$Allocation,
                           true_density = true_mean,
                           time = frame_raw$year)
      
      STRS_sim_mean[, , iboat, istrata, iter] <- temp_sim$mean_denisty
      STRS_sim_cv[, , iboat, istrata, iter] <- temp_sim$cv
      STRS_rel_bias_est[, , iboat, istrata, iter] <- temp_sim$rel_bias
      
      if(iter%%50 == 0) print(paste("Done with iter", iter, ",", iboat, "boat")) 
    }
    
    STRS_true_cv_array[, , iboat, istrata] <- 
      as.matrix(
        apply(X = STRS_sim_mean[, , iboat, istrata, ],
              MARGIN = 1:2,
              FUN = sd) / true_mean) 
    
    STRS_rrmse_cv_array[, , iboat, istrata] <- 
      sqrt(apply(X = sweep(x = STRS_sim_cv[, , iboat, istrata, ], 
                           MARGIN = 1:2,
                           STATS = STRS_true_cv_array[, , iboat, istrata],
                           FUN = "-")^2,
                 MARGIN = 1:2,
                 FUN = mean)) / apply(STRS_sim_cv[, , iboat, istrata, ], 
                                      MARGIN = 1:2, mean)
    
  }
}


##################################################
####   Save output
##################################################
if(which_domain == "full_domain") {
  save(list = c("STRS_sim_mean", "STRS_sim_cv", "STRS_rel_bias_est",
                "STRS_true_cv_array", "STRS_rrmse_cv_array"),
       file = paste0(github_dir, 
                     "Spatiotemporal_Optimization/",
                     "STRS_Sim_Res_spatiotemporal.RData"))
}

if (which_domain == "trawlable") {
  for(ivar in c("STRS_sim_mean", "STRS_sim_cv", "STRS_rel_bias_est",
                "STRS_true_cv_array", "STRS_rrmse_cv_array")) {
    assign(x = paste0(ivar, "_trawl"),
           value = get(ivar))
  }
  save(list = paste0(c("STRS_sim_mean", "STRS_sim_cv", "STRS_rel_bias_est",
                       "STRS_true_cv_array", "STRS_rrmse_cv_array"), "_trawl"),
       file = paste0(github_dir, 
                     "Spatiotemporal_Optimization/",
                     "STRS_Sim_Res_spatiotemporal.RData"))
}

