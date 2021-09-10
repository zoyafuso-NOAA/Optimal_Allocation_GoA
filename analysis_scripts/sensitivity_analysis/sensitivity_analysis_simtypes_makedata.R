###############################################################################
## Project:       Sensitivity Analysis for Survey Optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
VAST_sim_data_dir <- "D:/VAST_Runs/"
github_dir <- getwd()

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
result_dir <- paste0(github_dir, "/results/sensitivity_analysis/")
if(!dir.exists(result_dir)) dir.create(result_dir)

dens_vals <- array(dim = c(n_cells, ns_opt, n_years),
                   dimnames = list(NULL, common_names_opt, NULL))

##################################################
####   Loop over species and simulate data with the four different options
####   in the FishStatsUtils::simulate_data() function
##################################################
for (itype in c(4)) {
  for (ispp in common_names_opt) {
    
    pred_jnll_score <- subset(x = pred_jnll, 
                              spp_name == ispp, 
                              select = c("FALSE", "TRUE"))
    depth_in_model <- as.logical(names(which.min(pred_jnll_score)))
    
    temp_VAST_dir <- paste0(VAST_sim_data_dir, ispp, 
                            ifelse(depth_in_model == TRUE,
                                   yes = "_depth", no = ""), "/")
    
    load(paste0(temp_VAST_dir, "fit_sim.RData"))
    
    dyn.load(paste0(temp_VAST_dir, "VAST_v12_0_0.dll"))
    
    if(itype == 0) {
      dens_vals[, ispp, ] <- 
        fit_sim$Report$D_gct[, 1, years_included]
    }
    
    if(itype != 0) {
      sim_val <- FishStatsUtils::simulate_data(fit = fit_sim,
                                               type = itype,
                                               random_seed = 123423)
      
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
    }
    
    dyn.unload(paste0(temp_VAST_dir, "VAST_v12_0_0.dll"))
    print(paste("Finished with Type", itype, ispp))
  }

  #Save

  save(dens_vals, 
       file = paste0(github_dir, 
                     "/results/sensitivity_analysis/dens_vals_type", itype,
                     ".RData"))
}

