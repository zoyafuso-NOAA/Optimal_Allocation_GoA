###############################################################################
## Project:       Sensitivity Analysis for Survey Optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
VAST_sim_data_dir <- "F:/VAST_Runs/"
github_dir <- getwd()

library(VAST)

##################################################
####   Load data
##################################################
load("data/processed/optimization_data.RData")
load("data/processed/prednll_VAST_models.RData")
load("data/processed/grid_goa.RData")

##################################################
####   Result object
##################################################
dens_vals <- array(dim = c(ns_opt, n_cells, n_years, 5),
                   dimnames = list(common_names_opt, NULL, NULL, 
                                   c(paste0("Type ", 0:4))))

##################################################
####   Loop over species and simulate data with the four different options
####   in the FishStatsUtils::simulate_data() function
##################################################

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
  
  dens_vals[ispp, , , "Type 0"] <- fit_sim$Report$D_gct[, 1, years_included]
  print("Finished with Type 0")
  
  for (itype in 1:4) {
    sim_val <- FishStatsUtils::simulate_data(fit = fit_sim, 
                                             type = itype, 
                                             random_seed = 123423)
    
    pred_TF <- (1:length(sim_val$b_i))[-c(1:7900)]
    
    dens_vals[ispp, , ,paste0("Type ", itype)] <- 
      matrix(data = sim_val$b_i[pred_TF] * 0.001, 
             nrow = n_cells, 
             ncol = n_years)
    
    print(paste("Finished with Type", itype))
  }
  
  dyn.unload(paste0(temp_VAST_dir, "VAST_v12_0_0.dll"))
  
}

##################################################
####   Sensitivty Scenarios
##################################################
scen <- expand.grid(type = 0:4,
                    strata_vars = 1:2,
                    deep_stations = c(TRUE, FALSE))
irow = 1


frame <- frame_district[, c("domainvalue", "id", 
                            paste0("X", scen$strata_vars[irow]),
                            # "X1", 
                            # "X2",
                            "WEIGHT",
                            paste0("Y", spp_idx_opt), 
                            paste0("Y", spp_idx_opt,
                                   "_SQ_SUM"))]
