###############################################################################
## Project:       Knitting Result for univariate STRS optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   
###############################################################################
rm(list = ls())

###############################
## Import required packages
###############################
library(sp)
library(RColorBrewer)
library(raster)

###############################
## Set up directories
###############################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[2]

SamplingStrata_dir <- paste0(c("/Users/zackoyafuso/",
                               "C:/Users/Zack Oyafuso/",
                               "C:/Users/zack.oyafuso/")[which_machine],
                             "Downloads/SamplingStrata-master/R")

VAST_model <- "11" 
github_dir <- paste0(c("/Users/zackoyafuso/Documents", 
                       "C:/Users/Zack Oyafuso/Documents",
                       "C:/Users/zack.oyafuso/Work")[which_machine],
                     "/GitHub/Optimal_Allocation_GoA/model_", 
                     VAST_model, "/Single_Species_Optimization/")

###########################
## Load Data
###########################
load(paste0(dirname(github_dir), "/optimization_data.RData"))
load(paste0(dirname(dirname(github_dir)), "/data/Extrapolation_depths.RData"))

###########################
## Empty Result objects
###########################
master_res_df <- data.frame(id = 1:N)
master_settings <- data.frame()
master_strata_list <- list()
master_tradeoff <- list()

istrata <- 10

##########################
##########################

for (ispp in 1:ns) {
  runs = dir(paste0(github_dir, gsub(x = sci_names[ispp], 
                                     pattern = " ", 
                                     replacement = "_")), 
             full.names = T)
  
  if (length(runs) > 0) {
    temp_sample_size <- temp_cvs <- c()
    
    nruns = length(runs)
    
    for (irun in 1:nruns) {
      temp_file <- paste0(github_dir, 
                          gsub(x = sci_names[ispp], 
                               pattern = " ", 
                               replacement = "_"),
                          "/Str", istrata, "Run", irun, "/result_list.RData")
      
      if (file.exists(temp_file)) {
        load(temp_file)
        temp_sample_size <- c(temp_sample_size, result_list$n )
        temp_cvs <- c(temp_cvs, result_list[[3]])
        
        master_tradeoff[[ispp]] <-
          data.frame(cv = temp_cvs[order(temp_cvs, decreasing = T)], 
                     n = temp_sample_size[order(temp_cvs, decreasing = T)] )
      }
    }
    
    for (isample in 1:nboats) {
      irun <- which.min(abs(temp_sample_size - samples[isample]))
      temp_file <- paste0(github_dir, 
                          gsub(x = sci_names[ispp], 
                               pattern = " ", 
                               replacement = "_"),
                          "/Str", istrata, "Run", irun, "/result_list.RData")
      load(temp_file)
      
      master_settings <- rbind(master_settings,
                              data.frame(isample = isample,
                                         ispp = ispp,
                                         n = result_list$n,
                                         cv = as.numeric(result_list[[3]]) ))
      
      master_res_df <- cbind(master_res_df,
                            result_list[[1]]$indices$X1)
      
      master_strata_list <- c(master_strata_list, list(result_list[[2]]))
    }
  }
  
}

####################################
## Save
####################################
settings <- master_settings
settings$id = 1:nrow(settings)
res_df <- master_res_df
strata_list <- master_strata_list

save(list = c("res_df", "settings", "strata_list"),
     file = paste0(github_dir, "optimization_knitted_results.RData"))

