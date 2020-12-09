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
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[3]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")

##################################################
####   Load Data
##################################################
load(paste0(github_dir, "data/optimization_data.RData"))

##################################################
####   Empty Result Objects
##################################################
res_df <- data.frame(id = 1:N)
settings <- data.frame()
strata_stats_list <- strata_list <- list()

##################################################
####   Collect optimization results from each strata
##################################################
iid <- 1
#Temporary result objects
temp_res_df <- data.frame(id = 1:N)
temp_settings <- temp_settings_by_dom <- data.frame()
temp_strata_stats_list <- temp_strata_list <- list()

for (iboat in 1:nboats) {
   result_dir <- paste0(github_dir, "results/Spatiotemporal_Optimization/",
                        "boat", iboat, "/")
   
   runs <- grep(x = dir(result_dir), 
                pattern = "Run", 
                value = T )
   
   for (irun in runs) {
      temp_dir <- paste0(result_dir,  irun, "/result_list.RData")
      
      if (file.exists(temp_dir)) {
         load(temp_dir)
         
         #Solution: which strata is assigned to each extrapolation cell
         solution <- as.factor(paste(result_list$solution$framenew$DOMAINVALUE,
                                          result_list$solution$framenew$STRATO))
         solution <- as.integer(solution)
         
         temp_res_df <- cbind(temp_res_df, 
                              solution )
         
         #Strata characteristics: sample size, population, sampling rate, 
         # strata variable cuts
         temp_strata_list <- c(temp_strata_list, 
                               list(result_list[[2]]))
         
         #Strata statistics (mean and variance)
         temp_strata_stats_list <- c(temp_strata_stats_list, 
                                     list(result_list[[1]]$aggr_strata))
         
         #High-level settings: total sample size and expected CV across species
         species_cv <- result_list[[3]]
         # attributes(species_cv)$dimnames[[1]] <- ""
         attributes(species_cv)$dimnames[[2]] <- paste0("CV_", 1:ns_opt)
         n <- result_list$n
         
         temp_settings <- rbind(temp_settings, 
                                data.frame(id = iid,
                                           boat = iboat,
                                           n))
         
         temp_settings_by_dom <- rbind(temp_settings_by_dom, 
                                       data.frame(id = iid,
                                                  dom = 1:ndom,
                                                  boat = iboat,
                                                  species_cv))
         iid <- iid + 1
      }
   }
}

##################################################
####   Subset to only the solutions closest to 280, 550, and 820
##################################################
idx <- sapply(X = samples, 
              FUN = function(x) which.min(abs(temp_settings$n - x)) )

settings <- temp_settings[idx,]
settings_by_dom <- subset(x = temp_settings_by_dom, subset = id %in% idx)
res_df <- cbind(id = 1:N, temp_res_df[, 1 + idx])
names(res_df)[-1] <- paste0("sol_", 1:(ncol(res_df)-1))
strata_list <- temp_strata_list[idx]
strata_stats_list <- temp_strata_stats_list[idx]

##################################################
####   Save Objects
##################################################
save(list = c("res_df", "settings", "settings_by_dom",
              "strata_list", "strata_stats_list"),
     file = paste0(github_dir,
                   "results/Spatiotemporal_Optimization",
                   "/optimization_knitted_results.RData"))
