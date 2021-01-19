###############################################################################
## Project:         Extract table of strata boundaries
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     Extract the longitude and depth boundaries for 
##                  recommended optimization solutions
###############################################################################
rm(list = ls())

##################################################
####  Set up directories  
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[1]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")

##################################################
####  Load Data
##################################################
load(paste0(github_dir, 
            "data/Extrapolation_depths.RData"))
load(paste0(github_dir,
            "results/MS_optimization_knitted_results.RData"))

## Index for district-level, 2 boat, 3 and 5 districts per strata solutions
idx <- c(2, 5)

## Subset stratum settings here
strata_list <- strata_list[idx]

## Rescale UTM eastings back to UTM 5N
strata_list <- 
 lapply(X = strata_list,
        FUN = function(x) {
         x[, c("Lower_X1", "Upper_X1")] <- 
          x[, c("Lower_X1", "Upper_X1")] + min(Extrapolation_depths$E_km)
         
         x$Lower_X1 <- 
          Extrapolation_depths$Lon[match(x = round(x$Lower_X1), 
                                         table = round(Extrapolation_depths$E_km))]
         x$Upper_X1 <- 
          Extrapolation_depths$Lon[match(x = round(x$Upper_X1), 
                                         table = round(Extrapolation_depths$E_km))]
         
         x[, c("Lower_X1", "Upper_X1", "Lower_X2", "Upper_X2")] <- 
          round(x[, c("Lower_X1", "Upper_X1", "Lower_X2", "Upper_X2")], 1)
         
         return(x)
        })

