###############################################################################
## Project:         Plot Single Species Optimization Solutions
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     
###############################################################################
rm(list = ls())

##################################################
####  Import Libraries  
##################################################
library(RColorBrewer)

##################################################
####  Set up directories  ----
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[1]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")

output_dir <- paste0(c("/Users/zackoyafuso/", 
                       "C:/Users/Zack Oyafuso/")[which_machine], 
                     "Google Drive/MS_Optimizations/TechMemo/figures/")

##################################################
####    Load predicted density and optimization results 
##################################################
load(paste0(github_dir, "data/optimization_data.RData"))
load(paste0(github_dir, "/data/Extrapolation_depths.RData"))
load(paste0(github_dir, "results/Spatiotemporal_Optimization/",
            "optimization_knitted_results.RData"))

##################################################
####    
##################################################
par(mfcol = c(3, 3),
    mar = c(3,3,1,1))
for (istrata in 2:4) { 
 for (iboat in 1:nboats) { 
  
  #Which index to plot
  idx = settings$id[settings$boat == iboat & 
                     settings$strata == stratas[istrata]]
  
  temp_strata_list <- strata_list[[idx]]
  temp_nstrata <- nrow(temp_strata_list)
  
  plot(1,
       las = 1,
       ann = F,
       type = "n",
       xlim = range(Extrapolation_depths$E_km),
       ylim = range(Extrapolation_depths$DEPTH_EFH))
  
  for (temp_istrata in temp_nstrata:1) {
   rect(xleft = temp_strata_list$Lower_X2[temp_istrata] + 
         min(Extrapolation_depths$E_km),
        xright = temp_strata_list$Upper_X2[temp_istrata] + 
         min(Extrapolation_depths$E_km),
        ybottom = temp_strata_list$Lower_X1[temp_istrata],
        ytop = temp_strata_list$Upper_X1[temp_istrata],
        col = c(brewer.pal(n = 12, name = 'Paired'),
                brewer.pal(n = 12, name = 'Paired'))[temp_istrata])
  }
 }
}
