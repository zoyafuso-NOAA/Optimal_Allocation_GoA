###############################################################################
## Project:       Sensitivity Analysis for Survey Optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Simulate densities using the FishStatsUtils::simulate_data()
##                using type = 3 (simulate both random and fixed effects)
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
VAST_sim_data_dir <- "E:/VAST_Runs/"
github_dir <- getwd()

##################################################
####   Load packages
##################################################
library(VAST)
library(SamplingStrata)
library(raster)
library(sp)
library(RColorBrewer)

##################################################
####   Load data
##################################################
load("data/processed/optimization_data.RData")
load("data/processed/prednll_VAST_models.RData")
load("data/processed/grid_goa.RData")

##################################################
####   Result object
##################################################
if(!dir.exists("results/sensitivity_analysis_type3/")){
  dir.create("results/sensitivity_analysis_type3/")
}

##################################################
####   Loop over species and simulate data with the four different options
####   in the FishStatsUtils::simulate_data() function
##################################################

for (isim in 1:10) { ## Loop over replicates -- start
  
  ## Temporary result object
  dens_vals <- array(dim = c(n_cells, ns_opt, n_years),
                     dimnames = list(NULL, common_names_opt, NULL))
  
  for (ispp in common_names_opt) {  ## Loop over species -- start
    
    pred_jnll_score <- subset(x = pred_jnll, 
                              spp_name == ispp, 
                              select = c("FALSE", "TRUE"))
    depth_in_model <- as.logical(names(which.min(pred_jnll_score)))
    
    temp_VAST_dir <- paste0(VAST_sim_data_dir, ispp, 
                            ifelse(depth_in_model == TRUE,
                                   yes = "_depth", no = ""), "/")
    
    load(paste0(temp_VAST_dir, "fit_sim.RData"))
    
    dyn.load(paste0(temp_VAST_dir, "VAST_v12_0_0.dll"))
    
    sim_val <- FishStatsUtils::simulate_data(fit = fit_sim,
                                             type = 3,
                                             random_seed = 123423 + isim)
    
    pred_TF <- (1:length(sim_val$b_i))[-c(1:7900)]
    
    temp_density <- matrix(data = sim_val$b_i[pred_TF],
                           nrow = n_cells,
                           ncol = n_years)
    
    ## Temporary input objects for the survey simulation
    temp_density <- sweep(x = temp_density,
                          MARGIN = 1,
                          STATS = grid_goa$Area_km2,
                          FUN = "/")
    
    dens_vals[, ispp, ] <- temp_density
    
    print(paste("Finished with ", ispp, "Simulation", isim))
    dyn.unload(paste0(temp_VAST_dir, "VAST_v12_0_0.dll"))
  } ## Loop over species -- end
  
  ## Save simulated densities
  save(list = "dens_vals", 
       file = paste0("results/sensitivity_analysis_type3/dens_vals_sim", 
                     isim, ".RData"))
} ## Loop over replicates -- end
