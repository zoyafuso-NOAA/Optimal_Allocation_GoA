###############################################################################
## Project:       Simulate Current Stratified Random Survey
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

##################################################
####   Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[2]
VAST_model <- "11" 
which_domain <- c("full_domain", "trawlable")[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/",

if(!dir.exists(paste0(github_dir, "Survey_Comparison_Simulations/"))) {
  dir.create(paste0(github_dir, "Survey_Comparison_Simulations/"))
}

##################################################
####   Load simulation functions
##################################################
source( paste0(dirname(dirname(github_dir)), "/modified_functions/sim_fns.R") )

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load(paste0(github_dir, 'optimization_data.RData'))
load(paste0(dirname(dirname(github_dir)), '/data/Extrapolation_depths.RData'))

GOA_allocations <- readxl::read_xlsx(
  path = paste0(dirname(dirname(github_dir)), 
                '/data/GOA 2019 stations by stratum.xlsx')
)
GOA_allocations3 <- readxl::read_xlsx(
  path = paste0(dirname(dirname(github_dir)), 
                '/data/GOA2019_ 3 boat_825_RNDM_stations.xlsx')
) 

##################################################
####   Create indices for trawlable and shallow cells
##################################################
trawl_idx <- Extrapolation_depths$Id %in% cells_trawlable
shallow_idx <- Extrapolation_depths$Id %in% cells_shallower_than_700m
trawl_shallow_idx <- apply(X = cbind(trawl_idx, shallow_idx),
                           MARGIN = 1,
                           FUN = all)

if(which_domain == "full_domain") domain_idx <- rep(TRUE, N)
if(which_domain == "trawlable") domain_idx <- trawl_shallow_idx

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
  array(dim = c(NTime, ns, nboats, Niters), 
        dimnames = list(NULL, sci_names, NULL, NULL ))

Current_true_cv_array <- Current_rrmse_cv_array <- 
  array(dim = c(NTime, ns, nboats), 
        dimnames = list(NULL, sci_names, NULL))

##################################################
####   Simulate STRS based on current stratifications and allocations
##################################################
set.seed(23234)

for (iboat in 1:nboats) {
  for (iter in 1:Niters) {
    temp_sim <- do_STRS(density = frame_raw[, paste0("Y", 1:ns)],
                        cell_idx = which(domain_idx == T),
                        strata = Extrapolation_depths$stratum[domain_idx],
                        strata_to_use = allocations$Stratum, 
                        allocation = allocations[, paste0("boat", iboat)],
                        true_density = true_mean,
                        time = frame_raw$year)
    
    Current_sim_mean[, , iboat, iter] <- temp_sim$mean_denisty
    Current_sim_cv[, , iboat, iter] <- temp_sim$cv
    Current_rel_bias_est[, , iboat, iter] <- temp_sim$rel_bias
    
    if(iter%%50 == 0) print(paste("Done with iter", iter, ",", iboat, "boat")) 
  }
  
  Current_true_cv_array[, , iboat] <- 
    as.matrix(
      apply(X = Current_sim_mean[, , iboat, ],
            MARGIN = 1:2,
            FUN = sd) / true_mean) 
  
  Current_rrmse_cv_array[, , iboat] <- 
    sqrt(apply(X = sweep(x = Current_sim_cv[, , iboat, ], 
                         MARGIN = 1:2,
                         STATS = Current_true_cv_array[, , iboat],
                         FUN = "-")^2,
               MARGIN = 1:2,
               FUN = mean)) / apply(Current_sim_cv[, , iboat, ], MARGIN = 1:2, mean)
  
}

##################################################
####   Save output
##################################################
if(which_domain == "full_domain") {
  save(list = c("Current_sim_mean", "Current_sim_cv", "Current_rel_bias_est",
                "Current_true_cv_array", "Current_rrmse_cv_array"),
       file = paste0(github_dir, 
                     "Survey_Comparison_Simulations/",
                     "Survey_Simulation_Results.RData"))
}

if (which_domain == "trawlable") {
  for(ivar in c("Current_sim_mean", "Current_sim_cv", "Current_rel_bias_est",
                "Current_true_cv_array", "Current_rrmse_cv_array")) {
    assign(x = paste0(ivar, "_trawl"),
           value = get(ivar))
  }
  save(list = paste0(c("Current_sim_mean", "Current_sim_cv", "Current_rel_bias_est",
                       "Current_true_cv_array", "Current_rrmse_cv_array"), "_trawl"),
       file = paste0(github_dir, 
                     "Survey_Comparison_Simulations/",
                     "Survey_Simulation_Results.RData"))
}
