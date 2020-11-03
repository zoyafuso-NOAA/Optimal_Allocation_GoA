###############################################################################
## Project:       Simulate Various Surveys
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Simulate SRS, STRS surveys
###############################################################################
rm(list = ls())

##################################################
####  Set up directories
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[2]
VAST_model <- "11" 
which_domain <- c("full_domain", "trawlable")[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/",
                     which_domain, "/")

if(!dir.exists(paste0(github_dir, "Survey_Comparison_Simulations/"))) {
  dir.create(paste0(github_dir, "Survey_Comparison_Simulations/"))
}

##################################################
####   Load simulation functions
##################################################
source( paste0(dirname(dirname(github_dir)), "/modified_functions/sim_fns.R") )

##################################################
####   Load Data
##################################################
load(paste0(github_dir, "optimization_data.RData"))
load(paste0(dirname(dirname(github_dir)), '/data/Extrapolation_depths.RData'))

##################################################
####   Create indices for trawlable and shallow cells
##################################################
trawl_idx <- Extrapolation_depths$Id %in% cells_trawlable
shallow_idx <- Extrapolation_depths$Id %in% cells_shallower_than_700m
trawl_shallow_idx <- apply(X = cbind(trawl_idx, shallow_idx),
                           MARGIN = 1,
                           FUN = all)

if(which_domain == "full_domain") domain_idx <- 1:N
if(which_domain == "trawlable") domain_idx <- trawl_shallow_idx

##################################################
####   Result Objects
##################################################
Niters = 1000

SRS_sim_mean <-  SRS_sim_cv <- SRS_rel_bias_est <- 
  array(dim = c(NTime, ns, nboats, Niters), 
        dimnames = list(NULL, sci_names, NULL, NULL))

SRS_true_cv_array <- SRS_rrmse_cv_array <-
  array(dim = c(NTime, ns, nboats), 
        dimnames = list(NULL, sci_names, NULL))


##################################################
####   Simulate SRS Survey
##################################################
set.seed(234)

for (iboat in 1:nboats) {
  for (iter in 1:Niters) {
    temp_sim <- do_SRS(density = frame_raw[, paste0("Y", 1:ns)],
                       true_density = true_mean,
                       time = frame_raw$year,
                       n = c(280, 550, 820)[iboat],
                       cell_idx = 1:N)
    
    SRS_sim_mean[, , iboat, iter] <- temp_sim$mean_denisty
    SRS_sim_cv[, , iboat, iter] <- temp_sim$cv
    SRS_rel_bias_est[, , iboat, iter] <- temp_sim$rel_bias
    
    if(iter%%100 == 0) print(paste("Done with iter", iter, ",", iboat, "boat")) 
  }
  
  SRS_true_cv_array[, , iboat] <- 
    as.matrix(
      apply(X = SRS_sim_mean[, , iboat, ],
            MARGIN = 1:2,
            FUN = sd) / true_mean) 
  
  SRS_rrmse_cv_array[, , iboat] <- 
    sqrt(apply(X = sweep(x = SRS_sim_cv[, , iboat, ], 
                         MARGIN = 1:2,
                         STATS = SRS_true_cv_array[, , iboat],
                         FUN = "-")^2,
               MARGIN = 1:2,
               FUN = mean)) / apply(SRS_sim_cv[, , iboat, ], MARGIN = 1:2, mean)
  
}

##################################################
####   Save output
##################################################
if(which_domain == "full_domain") {
  save(list = c("SRS_sim_mean", "SRS_sim_cv", "SRS_rel_bias_est",
                "SRS_true_cv_array", "SRS_rrmse_cv_array"),
       file = paste0(github_dir, 
                     "Survey_Comparison_Simulations/",
                     "Simple_RS_Simulation_Results.RData"))
}

if (which_domain == "trawlable") {
  for(ivar in c("SRS_sim_mean", "SRS_sim_cv", "SRS_rel_bias_est",
                "SRS_true_cv_array", "SRS_rrmse_cv_array")) {
    assign(x = paste0(ivar, "_trawl"),
           value = get(ivar))
  }
  
  save(list = paste0(c("SRS_sim_mean", "SRS_sim_cv", "SRS_rel_bias_est",
                       "SRS_true_cv_array", "SRS_rrmse_cv_array"), "_trawl"),
       file = paste0(github_dir, 
                     "Survey_Comparison_Simulations/",
                     "Simple_RS_Simulation_Results.RData"))
}
