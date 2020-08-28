###############################################################################
## Project:       Simple Random Sample Survey Simulation
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Simulate SRS survey and calculate performance metrics  
##                for a given operating model specified on line 13
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
which_machine <- c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[2]
VAST_model <- "10b" 

github_dir <- paste0(c('/Users/zackoyafuso/Documents/', 
                       'C:/Users/Zack Oyafuso/Documents/',
                       'C:/Users/zack.oyafuso/Work/')[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/")

##################################################
####   Load predicted density and optimization results
##################################################
load(paste0(github_dir, 'optimization_data.RData'))

##################################################
####   Result Objects
##################################################
sim_mean <- sim_cv <- array(dim = c(NTime, ns, nboats, Niters), 
                            dimnames = list(NULL, sci_names, NULL, NULL))

true_cv_array <- rrmse_cv_array <- rel_bias_est <- rel_bias_cv <- 
  array(dim = c(NTime, ns, nboats), dimnames = list(NULL, sci_names, NULL))

##################################################
####   Simulate SRS Survey
##################################################
set.seed(233)
for (iyear in 1:NTime) {
  for (isample in 1:nboats) {
    for (iter in 1:Niters) {
      
      #Take a random sample based on number of boats
      samplesize <- samples[isample]
      sample_vec <- sample(x = 1:N, size = samplesize)
      sample_df <- subset(frame_raw, year == iyear)[sample_vec,]
      
      #Calculate Mean and standard error
      temp_sim_mean <- colMeans(sample_df[, paste0('Y', 1:ns)])
      temp_sim_var <- apply(sample_df[, paste0('Y', 1:ns)], 
                            MARGIN = 2, 
                            FUN = var)
      temp_sim_se <- sqrt(temp_sim_var / sample_size)
      
      #Save Mean and CV
      sim_mean[iyear, ,isample, iter] <- temp_sim_mean
      sim_cv[iyear, ,isample, iter] <- temp_sim_se / temp_sim_mean
    }
  }
  print(paste0('Done with year', iyear))
}

############################
## Simulation Metrics
############################
for (iyear in 1:NTime) {
  for (isample in 1:nboats) {
    for (ispp in 1:ns) {
      iter_mean <- sim_mean[iyear, ispp, isample, ]
      iter_cv <- sim_cv[iyear, ispp, isample, ]
      
      temp_true_cv <- sd(iter_mean)/true_mean[iyear,ispp]
      true_cv_array[iyear, ispp, isample] <- temp_true_cv
      
      abs_bias <- iter_mean - true_mean[iyear,ispp]
      rel_bias_est[iyear, ispp, isample] <- 
        100 * mean(abs_bias / true_mean[iyear,ispp])
      
      abs_bias <- iter_cv - temp_true_cv
      rel_bias_cv[iyear, ispp, isample] <- 
        100 * mean(abs_bias / temp_true_cv)
      
      rrmse_cv_array[iyear, ispp, isample] <- 
        sqrt(mean((iter_cv - temp_true_cv)^2)) / mean(iter_cv)
    }
  }
}

#######################
## Save results
#######################
for (ivar in c('rrmse_cv_array', 'true_cv_array', 
               'sim_mean', 'sim_cv', 
               'rel_bias_est', 'rel_bias_cv')) {
  assign(x = paste0('SRS_', ivar), value = get(ivar))
}

save(file = paste0(github_dir, 'Survey_Comparison_Simulations/',
                   'Simple_RS_Simulation_Results.RData'),
     list = c(paste0('SRS_', c('rrmse_cv_array', 
                               'true_cv_array', 'sim_mean', 'sim_cv',
                               'rel_bias_est', 'rel_bias_cv')) ))
