###############################################################################
## Project:         Plot Single Species Optimization Solutions
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     
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
load(paste0(dirname(github_dir), "/data/Extrapolation_depths.RData"))
load(paste0(result_dir, "optimization_knitted_results.RData"))

isample <- 1 #1, 2, or 3 boat solution

par(mfrow = c(3, 5), 
    mar = c(0, 0, 0, 0))

for (ispp in 1:ns) {
  
  #Which index to plot
  idx = settings$id[settings$isample == isample & settings$ispp == ispp]
  
  #Plot Solution
  goa <- SpatialPointsDataFrame(
    coords = Extrapolation_depths[,c("E_km", "N_km")],
    data = data.frame(Str_no = res_df[, 1 + idx]) )
  
  goa_ras <- raster(goa, 
                    resolution = 5)
  goa_ras <- rasterize(x = goa, 
                       y = goa_ras, 
                       field = "Str_no")
  
  image(goa_ras, 
        axes = F,
        ann = F,
        asp = 1,
        col = brewer.pal(n = 10, name = 'Paired'))
  box()
  
  #Simulate a sample solution
  temp_samples <- c()
  temp_strata <- nrow(strata_list[[idx]])
  temp_solution <- res_df[, idx + 1] 
  temp_allocation <- strata_list[[idx]]$Allocation
  
  for (istrata in 1:temp_strata) {
    temp_samples = c(temp_samples, 
                     sample(x = which(temp_solution == istrata),
                            size = temp_allocation[istrata])
    )
  }
  
  points(Extrapolation_depths[temp_samples, c("E_km", "N_km")], 
         pch = 16,
         cex = 0.5)
  
  legend("bottomright", 
         legend = sci_names[ispp],
         bty = "n",
         text.font = 3)
  
}


