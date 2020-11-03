###############################################################################
## Project:         Simulate Single Sspecies GoA Survey
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     Simulate surveys based on the optimized survey designs
###############################################################################
rm(list = ls())

##################################################
####  Set up directories  
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[2]
VAST_model <- "11" 

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/")

result_dir <- paste0(github_dir, "Single_Species_Optimization/")

##################################################
####    Load predicted density and optimization results
##################################################
load(paste0(github_dir, "optimization_data.RData"))
load(paste0(result_dir, "optimization_knitted_results.RData"))

##################################################
####   Result Objects
##################################################
sim_mean <- sim_cv <- array(dim = c(NTime, ns, nboats, Niters), 
                            dimnames = list(NULL, sci_names, NULL, NULL))

true_cv_array <- rrmse_cv_array <- rel_bias_est <- rel_bias_cv <- 
  array(dim = c(NTime, ns, nboats), 
        dimnames = list(NULL, sci_names, NULL))

##################################################
####   Simulating surveys from each optimized solution
##################################################
for (ispp in 1:ns) {
  for (isample in 1:nboats) {
    
    #Load optimization data
    sub_settings = settings[settings$ispp == ispp, ]
    # sub_settings = settings
    
    which_run <- which.min(abs(sub_settings$n - samples[isample]))
    temp_run = sub_settings$id[which_run]
    
    strata_allocation <- strata_list[[temp_run]]$Allocation
    stratapop <- strata_list[[temp_run]]$Population
    stratanos <- res_df[, 1+temp_run]
    
    #Remove strata with only 1 sample allocated
    str_idx <- strata_allocation > 1
    
    for (iyear in 1:NTime) {
      for (iter in 1:Niters) {
        
        #Sample based on the stratification allocations
        sample_vec <- c()
        for (i in which(str_idx == T)) {
          available_cells <- which(stratanos == i)
          sample_cells <- sample(x = available_cells, 
                                 size = strata_allocation[i], 
                                 replace = F)
          sample_vec <- c(sample_vec, sample_cells)
        }
        
        #Organize sample set and total number of samples
        sample_vec <- sort(sample_vec)
        n <- length(sample_vec)
        stratano_samp <- stratanos[sample_vec]
        sample_df <- subset(frame_raw, year == iyear)[sample_vec, ]
        
        #Calculate Stratum Mean Density and Variance
        stmt <- paste0("aggregate(Y", ispp, " ~ stratano_samp, ",
                       "data = sample_df, FUN = mean)")
        sample_mean <- eval(parse(text = stmt))[, -1]
        stmt <- paste0("aggregate(Y", ispp, " ~ stratano_samp, ",
                       "data = sample_df, FUN = var)")
        sample_var <- eval(parse(text = stmt))[, -1]
        
        #How many samples are allocated in each strata
        #How many sampling units are in each strata
        Wh <- (stratapop / N)[str_idx]
        wh <- (strata_allocation / stratapop)[str_idx]
        
        #Calculate Total Abundance and Variance, calculate CV
        SRS_var <- 
          sum(sample_var * Wh^2 * (1 - wh) / strata_allocation[str_idx])
        
        SRS_mean <- sum(sample_mean * Wh)
        
        strata_cv <- sqrt(SRS_var) / SRS_mean 
        
        #Record mean and CV values
        sim_mean[iyear, ispp, isample, iter] <- SRS_mean
        sim_cv[iyear, ispp, isample, iter] <- strata_cv
        
        if (iter%%100 == 0) {
          print(paste0("Finished with: Iteration ", iter, ", ", "Year ", iyear,
                       ", ", sci_names[ispp], ", and ", isample, " Boat"))
        }
      }
    }
  }
}


##################################################
####   Simulation Metrics
##################################################
for(iyear in 1:NTime){
  for(isample in 1:3){
    for(ispp in sci_names){
      iter_est <- sim_mean[iyear, ispp, isample, ]
      iter_cv <- sim_cv[iyear, ispp, isample, ]
      true_cv <- sd(iter_est) / true_mean[iyear, ispp]
      
      true_cv_array[iyear, ispp, isample] <- true_cv
      
      rrmse_cv_array[iyear, ispp, isample] <- 
        sqrt(mean((iter_cv - true_cv)^2)) / mean(iter_cv)
      
      abs_bias <- iter_est - true_mean[iyear, ispp]
      rel_bias_est[iyear, ispp, isample] <- 
        mean(100* abs_bias / true_mean[iyear, ispp])
      
      abs_bias <- iter_cv - true_cv
      rel_bias_cv[iyear, ispp, isample] <- 
        mean(100 * abs_bias / true_cv)
    }
  }
}

##################################################
####   Save results
##################################################
for(ivar in  c("rrmse_cv_array", "true_cv_array", 
               "sim_mean", "sim_cv", "rel_bias_est", "rel_bias_cv")){
  assign(x = paste0("SS_", ivar), value = get(ivar))
}

save(file = paste0(github_dir, "Single_Species_Optimization/",
                   "SS_Sim_Res_spatiotemporal.RData"),
     list = c(paste0("SS_", 
                     c("rrmse_cv_array", "true_cv_array", "sim_mean", "sim_cv", 
                       "rel_bias_est", "rel_bias_cv"))))

