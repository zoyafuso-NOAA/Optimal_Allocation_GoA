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
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[1]
VAST_model <- "11" 

github_dir <- paste0(c('/Users/zackoyafuso/Documents/', 
                       'C:/Users/Zack Oyafuso/Documents/',
                       'C:/Users/zack.oyafuso/Work/')[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/")

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load(paste0(github_dir, 'optimization_data.RData'))
load(paste0(dirname(github_dir), '/data/Extrapolation_depths.RData'))

GOA_allocations <- readxl::read_xlsx(
  path = paste0(dirname(github_dir), 
                '/data/GOA 2019 stations by stratum.xlsx')
)
GOA_allocations3 <- readxl::read_xlsx(
  path = paste0(dirname(github_dir), 
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


##################################################
####   Create dataframe of effort allocations across boats
##################################################
allocations <- data.frame(Stratum = sort(unique(GOA_allocations3$stratum)),
                          boat3 = aggregate(id ~ stratum, 
                                            data = GOA_allocations3, 
                                            FUN = length)$id,
                          boat2 = c(GOA_allocations$`Number stations`,
                                    rep(0,5)))
allocations$boat1 = ceiling(allocations$boat2 / 2)

allocations$boat1 = ifelse(allocations$boat1 == 0, 0, 
                           ifelse(allocations$boat1 == 1, 2, 
                                  allocations$boat1))

##################################################
####   Attribute each grid point to the current stratification
##################################################
goa_grid <- sp::SpatialPointsDataFrame(
  coords = Extrapolation_depths[, c("E_km", "N_km")],
  proj4string=CRS("+proj=utm +zone=5N +units=km"),
  data = Extrapolation_depths)

##################################################
####   Result Objects
##################################################
sim_mean <- sim_cv <- sim_mean_trawl <- sim_cv_trawl <- 
  array(dim = c(NTime, ns, Niters, nboats), 
        dimnames = list(NULL, sci_names, NULL, NULL ))

true_cv_array <- rrmse_cv_array <- rel_bias_est <- true_cv_array_trawl <- rrmse_cv_array_trawl <- rel_bias_est_trawl <- 
  array(dim = c(NTime, ns, nboats), 
        dimnames = list(NULL, sci_names, NULL))

##################################################
####   Simulate STRS based on current stratifications and allocations
##################################################
set.seed(23234)

for (isample in 1:nboats) {
  
  #Adjust sample size proportionally
  nh <- allocations[, paste0('boat', isample)]
  sampled_strata <- allocations$Stratum[nh > 0]
  nstrata <- length(sampled_strata)
  
  nh <- nh[allocations$Stratum %in% sampled_strata]
  
  #strata constraints
  stratano <- rep(x = sampled_strata, 
                  times = nh)
  
  Nh <- table(goa_grid$stratum)[as.character(sampled_strata)]
  Wh <- Nh / N
  wh <- nh / Nh
  
  for (iyear in 1:NTime) {
    
    #Subset densities
    sub_df <- subset(frame_raw, 
                     year == iyear)
    
    for (iter in 1:Niters) {
      
      #Sample grid cells on full domain by stratum allocations
      sample_idx = c()
      for (istrata in 1:nstrata) {
        temp_nh <- nh[istrata]
        temp_strata <- sampled_strata[istrata]
        str_idx <- which(goa_grid$stratum == temp_strata)
        sample_idx <- c(sample_idx, 
                        sample(x = str_idx, 
                               size = temp_nh))
      }
      
      sample_df <- sub_df[sample_idx, ]
      
      #Sample grid cells on only trawlable cells by stratum allocations
      sample_idx_trawl = c()
      for (istrata in 1:nstrata) {
        temp_nh <- nh[istrata]
        temp_strata <- sampled_strata[istrata]
        str_idx <- which(goa_grid$stratum == temp_strata)
        str_idx <- str_idx[str_idx %in% which(trawl_shallow_idx == T)]
        
        sample_idx_trawl <- c(sample_idx_trawl, 
                              sample(x = str_idx, 
                                     size = temp_nh))
      }
      
      sample_df_trawl <- sub_df[sample_idx_trawl, ]
      
      #Calculate Stratum Mean Density and Variance
      stmt <- paste0('aggregate(cbind(',
                     paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                     ") ~ stratano, data = sample_df, FUN = mean)")
      sample_mean <- eval(parse(text = stmt))[, -1]
      stmt <- paste0('aggregate(cbind(',
                     paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                     ") ~ stratano, data = sample_df_trawl, FUN = mean)")
      sample_mean_trawl <- eval(parse(text = stmt))[, -1]
      
      stmt <- paste0('aggregate(cbind(',
                     paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                     ") ~ stratano, data = sample_df, FUN = var)")
      sample_var <- eval(parse(text = stmt))[, -1]
      stmt <- paste0('aggregate(cbind(',
                     paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                     ") ~ stratano, data = sample_df_trawl, FUN = var)")
      sample_var_trawl <- eval(parse(text = stmt))[, -1]
      
      #Calculate Total Mean, Variance, CV
      SRS_var <- colSums(sweep(x = sample_var, 
                               MARGIN = 1, 
                               STATS = Wh^2 * (1 - wh) / nh,
                               FUN = '*'))
      SRS_var_trawl <- colSums(sweep(x = sample_var_trawl, 
                                     MARGIN = 1, 
                                     STATS = Wh^2 * (1 - wh) / nh,
                                     FUN = '*'))
      
      SRS_mean <- colSums(sweep(x = sample_mean, 
                                MARGIN = 1, 
                                STATS = Wh,
                                FUN = '*'))
      SRS_mean_trawl <- colSums(sweep(x = sample_mean_trawl, 
                                      MARGIN = 1, 
                                      STATS = Wh,
                                      FUN = '*'))
      
      strata_cv <- sqrt(SRS_var) / SRS_mean 
      strata_cv_trawl <- sqrt(SRS_var_trawl) / SRS_mean_trawl
      
      #Record mean and CV values
      sim_mean[iyear, , iter, isample] <- SRS_mean
      sim_cv[iyear, , iter, isample] <- strata_cv
      
      sim_mean_trawl[iyear, , iter, isample] <- SRS_mean_trawl
      sim_cv_trawl[iyear, , iter, isample] <- strata_cv_trawl
      
    }
    print(paste('Finished with n =', samples[isample], 'Year', iyear))
  }
}

#################################
## Simulation Metrics
#################################
for (isample in 1:nboats) {
  for (iyear in 1:NTime) {
    for (ispp in sci_names) {
      
      #Simulated Mean Estimates
      iter_est <- sim_mean[iyear, ispp, , isample]
      iter_est_trawl <- sim_mean_trawl[iyear, ispp, , isample]
      
      #Simulated CVs
      iter_cv <- sim_cv[iyear, ispp, , isample]
      iter_cv_trawl <- sim_cv_trawl[iyear, ispp, , isample]
      
      #True CV
      true_cv <- sd(iter_est) / true_mean[iyear, ispp]
      true_cv_array[iyear, ispp, isample] <- true_cv
      
      true_cv_trawl <- sd(iter_est_trawl) / true_mean[iyear, ispp]
      true_cv_array_trawl[iyear, ispp, isample] <- true_cv_trawl
      
      #RRMSE of CV
      rrmse_cv_array[iyear, ispp, isample] <-
        sqrt(mean((iter_cv - true_cv)^2)) / mean(iter_cv)
      rrmse_cv_array_trawl[iyear, ispp, isample] <-
        sqrt(mean((iter_cv_trawl - true_cv_trawl)^2)) / mean(iter_cv_trawl)
      
      #Relative Bias of Mean Estimate
      abs_bias <- iter_est - true_mean[iyear, ispp]
      rel_bias_est[iyear, ispp, isample] <- 
        100 * mean(abs_bias / true_mean[iyear, ispp])
      abs_bias_trawl <- iter_est_trawl - true_mean[iyear, ispp]
      rel_bias_est_trawl[iyear, ispp, isample] <- 
        100 * mean(abs_bias_trawl / true_mean[iyear, ispp])

    }
  }
}

##################################
## Save
##################################
for(ivar in  c('rrmse_cv_array', 'true_cv_array',
               'sim_mean', 'sim_cv', 'rel_bias_est',
               
               'rrmse_cv_array_trawl', 'true_cv_array_trawl', 
               'sim_mean_trawl', 'sim_cv_trawl', 'rel_bias_est_trawl')){
  assign(x = paste0('Survey_', ivar), value = get(ivar))
}

save(file = paste0(github_dir, 'Survey_Comparison_Simulations/',
                   'Survey_Simulation_Results.RData'),
     list = c(paste0('Survey_', c('sim_mean', 'sim_cv',
                                  'rrmse_cv_array', 'true_cv_array', 
                                  'rel_bias_est',
                                  
                                  'rrmse_cv_array_trawl', 'true_cv_array_trawl', 
                                  'sim_mean_trawl', 'sim_cv_trawl', 
                                  'rel_bias_est_trawl'))))

