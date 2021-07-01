###############################################################################
## Project:       Synthesize the Cross Validation 
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Loop through every CV fold and extract performance metrics
##                Maximum Gradient, RRMSE, Predictive Joint Negative LogLike.
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(tidyr)

##################################################
####   Set VAST_dir to the external drive that the VAST runs are stored
####   Import the VAST data input
##################################################
VAST_dir <- "D:/VAST_Runs/"
VAST_bias_corrected_dir <- "D:/VAST_Runs_Bias_Corrected/"
data_df <- read.csv("data/processed/goa_vast_data_input.csv")

##################################################
####   Result object
##################################################
cross_val_results <- data.frame()

##################################################
####   Loop through each model (species, depth, CV fold) and calculate RRMSE
####   synthesize RRMSE, pred_jnll, maximum gradient
##################################################

spp_names <- sort(unique(data_df$COMMON_NAME))
for (ispp in spp_names){ ## Loop through species -- start
  
  mean_obs_cpue <- with(subset(data_df, COMMON_NAME == ispp), 
                        mean(WEIGHT / EFFORT))
  
  for (depth_in_model in c("", 
                           "_depth")) { ## Loop through depth models -- start
    for (icv in 1:10) { ## Loop through cv folds -- start
      
      ## Only go through CV folds where a model run was successful 
      fit_file <- paste0(VAST_dir, ispp, depth_in_model, 
                         "/CV_", icv, "/", "fit.RData")
      
      if(file.exists(fit_file)) {
        
        ## Load performance metrics
        load(paste0(VAST_dir, ispp, depth_in_model, 
                    "/CV_", icv, "/", "crossval_fit_performance.RData"))
        load(paste0(VAST_dir, ispp, depth_in_model,
                    "/CV_", icv, "/", "parameter_estimates.RData"))
        
        ## rbind the rrmse, max gradient, and pred_jnll
        rmse <- with(cv_performance$cpues, 
                     sqrt(mean((pred_cpue - obs_cpue)^2, na.rm = TRUE)) )
        rrmse <- rmse / mean_obs_cpue
        
        cross_val_results <- 
          rbind(cross_val_results, 
                data.frame(spp_name = ispp,
                           depth_in_model = ifelse(depth_in_model == "", F, T),
                           cv_fold = icv,
                           max_gradient = parameter_estimates$max_gradient,
                           rrmse = rrmse,
                           pred_jnll = cv_performance$prednll )) 
        
      }  ## Loop through depth models -- start
    } ## Loop through cv folds -- end
  }
  
} ## Loop through species --end

##################################################
####   Synthesize Cross Validation Results
##################################################
converged <- spread(data = aggregate(max_gradient ~ depth_in_model + spp_name,
                                     FUN = function(x) sum(x < 1e4),
                                     data = cross_val_results),
                    key = depth_in_model,
                    value = max_gradient)

rrmse <- spread(data = aggregate(rrmse ~ depth_in_model + spp_name,
                                 FUN = function(x) round(mean(x, na.rm = T), 
                                                         digits = 3),
                                 data = cross_val_results,
                                 subset = max_gradient < 1e-4),
                key = depth_in_model,
                value = rrmse)

pred_jnll <- spread(data = aggregate(pred_jnll ~ depth_in_model + spp_name,
                                     FUN = function(x) round(mean(x, na.rm = T)),
                                     data = cross_val_results,
                                     subset = max_gradient < 1e-4),
                    key = depth_in_model,
                    value = pred_jnll)

##################################################
####   Save
##################################################
save(list = c("converged", "cross_val_results", "pred_jnll", "rrmse"),
     file = "data/processed/prednll_VAST_models.RData")

##################################################
####   Synthesize the densities and abundance indices for the best models
##################################################
load(paste0(VAST_dir, "arrowtooth flounder/fit.RData"))
year_idx <- 1 + as.integer(names(table(fit$data_list$t_i)))
n_years <- length(year_idx)
n_cells <- dim(fit$Report$D_gct)[1]
n_spp <- nrow(pred_jnll)
spp_names <- pred_jnll$spp_name

D_gct <- I_gct <- array(dim = c(n_cells, n_spp, n_years), 
                        dimnames = list(NULL, spp_names, NULL))

I_gc_bias_corrected <- array(dim = c(n_spp, n_years), 
                             dimnames = list(spp_names, NULL))

index <- index_bias_corrected <- data.frame()

for(irow in 1:nrow(pred_jnll)) { ## Loop over species -- start
  
  ## Extract file name of best model
  ispp <- pred_jnll$spp_name[irow]
  depth_in_model <- c(FALSE, TRUE)[which.min(pred_jnll[irow, 2:3])]
  filename <- paste0(VAST_dir, ispp, ifelse(test = depth_in_model, 
                                            yes = "_depth/", 
                                            no = "/"),
                     "fit.RData")
  
  ## Load data
  load(filename)
  
  ## Extract grid specific density and abundance
  D_gct[, irow, ] <- fit$Report$D_gct[, 1, year_idx]
  I_gct[, irow, ] <- fit$Report$Index_gctl[, 1, year_idx, 1]
  
  ## Rbind abundance index table
  temp_index <- read.csv(paste0(VAST_dir, ispp, 
                                ifelse(test = depth_in_model, 
                                       yes = "_depth/", 
                                       no = "/"), 
                                "diagnostics/Table_for_SS3.csv" ))
  
  index <- rbind(index, data.frame(Species_Name = ispp, 
                                   temp_index[year_idx, ]))
  
  ## Rbind bias-indexed abundance index table
  temp_index_BC <- read.csv(paste0(VAST_bias_corrected_dir, ispp, 
                                   ifelse(test = depth_in_model, 
                                          yes = "_depth/", 
                                          no = "/"),
                                   "diagnostics/Table_for_SS3.csv"))
  
  index_bias_corrected <- rbind(index_bias_corrected, 
                                data.frame(Species_Name = ispp, 
                                           temp_index_BC[year_idx, ]))
  
  
  print(paste0(ispp, ifelse(depth_in_model, " with Depth", " without Depth")))
}  ## Loop over species -- end

I_ct <- apply(I_gct, MARGIN = c(2:3), FUN = sum)

##################################################
####   Save
##################################################
save(list = "D_gct", file = "data/processed/VAST_fit_D_gct.RData")
save(list = "I_gct", file = "data/processed/VAST_fit_I_gct.RData")
save(list = "index", file = "data/processed/index.RData")
save(list = "index_bias_corrected", 
     file = "data/processed/index_bias_corrected.RData")
save(list = "I_ct", file = "data/processed/I_ct.RData")
