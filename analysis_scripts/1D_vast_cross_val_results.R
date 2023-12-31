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
VAST_dir <- "temp/"
data_df <- read.csv("data/processed/goa_data_geostat.csv")

##################################################
####   Result object
##################################################
cross_val_results <- data.frame()

##################################################
####   Loop through each model (species, depth, CV fold) and calculate RRMSE
####   synthesize RRMSE, pred_jnll, maximum gradient
##################################################

spp_names <- sort(x = unique(x = data_df$Species))
for (ispp in spp_names){ ## Loop through species -- start
  
  mean_obs_cpue <- with(subset(data_df, Species == ispp), 
                        mean(Catch_KG / AreaSwept_km2))
  
  for (depth_in_model in c("", 
                           "_depth")) { ## Loop through depth models -- start
    for (icv in 1:10) { ## Loop through cv folds -- start
      
      ## Only go through CV folds where a model run was successful 
      fit_file <- paste0(VAST_dir, ispp, depth_in_model, 
                         "/CV_", icv, "/", "fit.RDS")
      
      if(file.exists(fit_file)) {
        
        ## Load performance metrics
        cv_performance <- 
          readRDS(paste0(VAST_dir, ispp, depth_in_model, 
                         "/CV_", icv, "/", "crossval_fit_performance.RDS"))
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
pred_jnll <- spread(data = aggregate(pred_jnll ~ depth_in_model + spp_name,
                                     FUN = function(x) round(sum(x, na.rm = T)),
                                     data = cross_val_results),
                    key = depth_in_model,
                    value = pred_jnll)

##################################################
####   Save
##################################################
save(list = c("cross_val_results", "pred_jnll"),
     file = "data/processed/prednll_VAST_models.RData")

##################################################
####   Synthesize the densities and abundance indices for the best models
##################################################
fit <- readRDS(file = paste0(VAST_dir, "arrowtooth flounder/fit.RDS"))
year_idx <- 1 + as.integer(x = names(table(fit$data_list$t_i)))
n_years <- length(x = year_idx)
n_cells <- dim(x = fit$Report$D_gct)[1]
n_spp <- nrow(x = pred_jnll)
spp_names <- pred_jnll$spp_name

D_gct <- array(dim = c(n_cells, n_spp, n_years), 
               dimnames = list(NULL, spp_names, NULL))

index <- data.frame()

for (irow in 1:nrow(x = pred_jnll)) { ## Loop over species -- start
  
  ## Extract file name of best model
  ispp <- pred_jnll$spp_name[irow]
  depth_in_model <- c(FALSE, TRUE)[which.min(x = pred_jnll[irow, 2:3])]
  
  filename <- paste0(VAST_dir, ispp, ifelse(test = depth_in_model, 
                                            yes = "_depth/", 
                                            no = "/"),
                     "fit.RDS")
  
  ## Load data
  fit <- readRDS(file = filename)
  
  ## Extract grid specific density and abundance
  D_gct[, irow, ] <- fit$Report$D_gct[, 1, year_idx]
  
  ## Rbind abundance index table
  temp_index <- read.csv(paste0(VAST_dir, ispp, 
                                ifelse(test = depth_in_model, 
                                       yes = "_depth/", 
                                       no = "/"), 
                                "diagnostics/Index.csv" ))
  
  index <- rbind(index, data.frame(Species_Name = ispp, 
                                   temp_index[year_idx, ]))

  print(paste0(ispp, ifelse(depth_in_model, " with Depth", " without Depth")))
}  ## Loop over species -- end

##################################################
####   Save
##################################################
saveRDS(object = D_gct, file = "data/processed/VAST_fit_D_gct.RDS")
saveRDS(object = index, file = "data/processed/index.RDS")
