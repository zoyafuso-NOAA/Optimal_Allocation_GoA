###############################################################################
## Project:       Synthesize Optimization Results
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Synthesize all optimization results
##                sample sizes, expected CVs, solutions, allocations, etc
###############################################################################
rm(list = ls())

##################################################
####  Set up directories
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[2]
VAST_model <- "11" 
which_domain <- c("full_domain", "trawlable")[1]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/",
                     which_domain, "/")

##################################################
####   Load Data
##################################################
load(paste0(github_dir, "optimization_data.RData"))

##################################################
####   Define which optimization settings is being worked on
####
####   which_variance:
####   Spatial: spatial variance for stratum variance
####   Spatiotemporal: spatiotemporal variance for stratum variance
####
####   which_constraint: 
####   one_CV: One CV constraint applied to all species
####   "": species specific CV constraints, assumed to be the default
##################################################
# result_dir <- paste0(github_dir, which_variance)

##################################################
####   Empty Result Objects
##################################################
res_df <- data.frame(id = 1:N)
settings <- data.frame()
strata_stats_list <- strata_list <- list()
stratas <- c(15)
NStrata <- length(stratas)

##################################################
####   Collect optimization results from each strata
##################################################
for (istrata in 1:NStrata) {
   
   #Temporary result objects
   temp_strata <- stratas[istrata]
   temp_res_df <- data.frame(id = 1:N)
   temp_settings <- data.frame()
   temp_strata_stats_list <- temp_strata_list <- list()
   
   for (iboat in 1:nboats) {
      result_dir <- paste0(github_dir, "Spatiotemporal_Optimization/",
                           "boat", iboat, "/")
      
      runs <- grep(x = dir(result_dir), 
                   pattern = paste0("Str", temp_strata, "Run"), 
                   value = T )
      
      for (irun in runs) {
         temp_dir <- paste0(result_dir,  irun, "/result_list.RData")
         
         if (file.exists(temp_dir)) {
            load(temp_dir)
            
            #Solution: which strata is assigned to each extrapolation cell
            temp_res_df <- cbind(temp_res_df, 
                                 result_list[[1]]$indices$X1 )
            
            #Strata characteristics: sample size, population, sampling rate, 
            # strata variable cuts
            temp_strata_list <- c(temp_strata_list, 
                                  list(result_list[[2]]))
            
            #Strata statistics (mean and variance)
            temp_strata_stats_list <- c(temp_strata_stats_list, 
                                        list(result_list[[1]]$aggr_strata))
            
            #High-level settings: total sample size and expected CV across species
            species_cv <- result_list[[3]]
            attributes(species_cv)$dimnames[[1]] <- ""
            attributes(species_cv)$dimnames[[2]] <- paste0("CV_", 1:ns)
            cv <- max(as.numeric(species_cv))
            n <- result_list$n
            
            temp_settings <- rbind(temp_settings, 
                                   data.frame(strata = temp_strata, 
                                              boat = iboat,
                                              n, cv, species_cv))
         }
      }
   }
   
   ##################################################
   ####   Subset to only the solutions closest to 280, 550, and 820
   ##################################################
   idx <- sapply(X = samples, 
                 FUN = function(x) which.min(abs(temp_settings$n - x)) )
   
   settings <- rbind(settings, temp_settings[idx,])
   res_df <- cbind(res_df, temp_res_df[, 1 + idx])
   strata_list <- c(strata_list, temp_strata_list[idx])
   strata_stats_list <- c(strata_stats_list, temp_strata_stats_list[idx])
   
}

settings$id <- 1:nrow(settings)
names(res_df)[-1] <- paste0("sol_", 1:(ncol(res_df)-1))

##################################################
####   Save Objects
##################################################
save(list = c("res_df", "settings", "strata_list", "strata_stats_list"),
     file = paste0(dirname(result_dir), "/optimization_knitted_results.RData"))
