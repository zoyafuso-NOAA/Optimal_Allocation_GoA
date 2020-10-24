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
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' =2, 'Zack_GI_PC' = 3)[1]
VAST_model <- "11" 

github_dir <- paste0(c('/Users/zackoyafuso/Documents/', 
                       'C:/Users/Zack Oyafuso/Documents/',
                       'C:/Users/zack.oyafuso/Work/')[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/")

##################################################
####   Load predicted density and optimization results
##################################################
load(paste0(github_dir, 'optimization_data.RData'))
load(paste0(dirname(github_dir), '/data/Extrapolation_depths.RData'))

##################################################
####   Create indices for trawlable and shallow cells
##################################################
trawl_idx <- Extrapolation_depths$Id %in% cells_trawlable
shallow_idx <- Extrapolation_depths$Id %in% cells_shallower_than_700m
trawl_shallow_idx <- apply(X = cbind(trawl_idx, shallow_idx),
                           MARGIN = 1,
                           FUN = all)

##################################################
####   Result Objects
##################################################
sim_mean <- sim_cv <- sim_mean_trawl <- sim_cv_trawl <- 
  array(dim = c(NTime, ns, nboats, Niters), 
        dimnames = list(NULL, sci_names, NULL, NULL))

true_cv_array <- rrmse_cv_array <- rel_bias_est <- true_cv_array_trawl <- rrmse_cv_array_trawl <- rel_bias_est_trawl <- 
  array(dim = c(NTime, ns, nboats), dimnames = list(NULL, sci_names, NULL))

##################################################
####   Simulate SRS Survey
##################################################
set.seed(233)
for (iyear in 1:NTime) {
  for (isample in 1:nboats) {
    for (iter in 1:Niters) {
      
      ##Full Domain
      #Take a random sample based on number of boats
      samplesize <- samples[isample]
      sample_vec <- sample(x = 1:N, size = samplesize)
      sample_df <- subset(frame_raw, year == iyear)[sample_vec,]
      
      #Calculate Mean and standard error
      temp_sim_mean <- colMeans(sample_df[, paste0('Y', 1:ns)])
      temp_sim_var <- apply(sample_df[, paste0('Y', 1:ns)], 
                            MARGIN = 2, 
                            FUN = var)
      temp_sim_se <- sqrt(temp_sim_var / samplesize)
      
      #Save Mean and CV
      sim_mean[iyear, ,isample, iter] <- temp_sim_mean
      sim_cv[iyear, ,isample, iter] <- temp_sim_se / temp_sim_mean
      
      
      ##SUbsetted Spatial Domain
      #Take a random sample based on number of boats
      samplesize <- samples[isample]
      sample_vec <- sample(x = which(trawl_shallow_idx == T), 
                           size = samplesize)
      sample_df <- subset(frame_raw, year == iyear)[sample_vec,]
      
      #Calculate Mean and standard error
      temp_sim_mean <- colMeans(sample_df[, paste0('Y', 1:ns)])
      temp_sim_var <- apply(sample_df[, paste0('Y', 1:ns)], 
                            MARGIN = 2, 
                            FUN = var)
      temp_sim_se <- sqrt(temp_sim_var / samplesize)
      
      #Save Mean and CV
      sim_mean_trawl[iyear, ,isample, iter] <- temp_sim_mean
      sim_cv_trawl[iyear, ,isample, iter] <- temp_sim_se / temp_sim_mean
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
      
      #Simulated Mean Estimates
      iter_mean <- sim_mean[iyear, ispp, isample, ]
      iter_mean_trawl <- sim_mean_trawl[iyear, ispp, isample, ]
      
      #Simulated CVs
      iter_cv <- sim_cv[iyear, ispp, isample, ]
      iter_cv_trawl <- sim_cv_trawl[iyear, ispp, isample, ]
      
      #True CV
      temp_true_cv <- sd(iter_mean)/true_mean[iyear,ispp]
      true_cv_array[iyear, ispp, isample] <- temp_true_cv
      
      temp_true_cv_trawl <- sd(iter_mean_trawl)/true_mean[iyear,ispp]
      true_cv_array_trawl[iyear, ispp, isample] <- temp_true_cv_trawl
      
      #Relative Bias of Estimate
      abs_bias <- iter_mean - true_mean[iyear,ispp]
      rel_bias_est[iyear, ispp, isample] <- 
        100 * mean(abs_bias / true_mean[iyear,ispp])
      
      abs_bias_trawl <- iter_mean_trawl - true_mean[iyear,ispp]
      rel_bias_est_trawl[iyear, ispp, isample] <- 
        100 * mean(abs_bias_trawl / true_mean[iyear,ispp])
      
      #RRMSE of CV
      rrmse_cv_array[iyear, ispp, isample] <- 
        sqrt(mean((iter_cv - temp_true_cv)^2)) / mean(iter_cv)
      
      rrmse_cv_array_trawl[iyear, ispp, isample] <- 
        sqrt(mean((iter_cv_trawl - temp_true_cv_trawl)^2)) / 
        mean(iter_cv_trawl)
    }
  }
}

#######################
## Save results
#######################
for (ivar in c('rrmse_cv_array', 'true_cv_array', 
               'sim_mean', 'sim_cv', 
               'rel_bias_est',
               
               'rrmse_cv_array_trawl', 'true_cv_array_trawl', 
               'sim_mean_trawl', 'sim_cv_trawl', 
               'rel_bias_est_trawl')) {
  assign(x = paste0('SRS_', ivar), value = get(ivar))
}

save(file = paste0(github_dir, 'Survey_Comparison_Simulations/',
                   'Simple_RS_Simulation_Results.RData'),
     list = c(paste0('SRS_', c('rrmse_cv_array', 'true_cv_array', 
                               'sim_mean', 'sim_cv', 
                               'rel_bias_est',
                               
                               'rrmse_cv_array_trawl', 'true_cv_array_trawl', 
                               'sim_mean_trawl', 'sim_cv_trawl', 
                               'rel_bias_est_trawl') ) ))
