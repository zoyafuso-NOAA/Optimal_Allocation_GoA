###############################################################################
## Project:         Plot Single Species Optimization Solutions
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     
###############################################################################
rm(list = ls())

##################################################
####  Set up output directories  
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[2]
output_dir <- paste0(c("/Users/zackoyafuso/", 
                       "C:/Users/Zack Oyafuso/")[which_machine], 
                     "Google Drive/MS_Optimizations/TechMemo/figures/")

##################################################
####  Import Libraries  
##################################################
library(sp)
library(raster)
library(RColorBrewer)

##################################################
####    Load predicted density and optimization results 
##################################################
load("data/processed/optimization_data.RData")
load("data/processed/grid_goa.RData")

##################################################
####   Set Constants
##################################################
isample <- 2 #1, 2, or 3 boat solution

##################################################
####   Plot
##################################################
{
  ## Open Device
  png(filename = paste0(output_dir, "Figure02_SS_solutions.png"),
      height = 190, width = 150, units = "mm",
      res = 500)
  
  ## Plot setup
  layout(mat = matrix(c(2, 1), nrow = 1))
  par(mar = c(0, 0, 0, 0))
  
  for (ispp in 1:ns_all) {
    
    ## Load Species Optimization Data
    load(paste0("results/scenario_A/Single_Species_Optimization/", 
                c(rev(common_names_eval), rev(common_names_opt))[ispp], 
                "/boat_", isample, "/result_list.RData"))
    
    if (ispp %in% c(ns_all/2 + 1, 1)) { 
      
      ## Base Plot
      plot(x = 1, y = 1, 
           type = "n", axes = F, ann = F,
           xlim = range(grid_goa$E_km),
           ylim = with(grid_goa, 
                       c(min(N_km), max(N_km) + 7 * y_range)),
           asp = 1)
      
      ## Since we are plotting multiple maps on a plot, we shift the plots
      ## using an offset value which is initialized at 0
      offset_val <- 0
    }
    
    ## Set up spatial object to plot solution
    goa <- sp::SpatialPointsDataFrame(
      coords = grid_goa[, c("E_km", "N_km")],
      data = data.frame("Str_no" = result_list$sol_by_cell,
                        stringsAsFactors = TRUE))
    
    goa_ras <- raster(goa,
                      resolution = 5)
    goa_ras <- rasterize(x = goa,
                         y = goa_ras,
                         field = "Str_no")
    
    ## Plot Solution
    offset_y <- 0.6 * y_range * offset_val
    goa_ras <- raster::shift(goa_ras, dy = offset_y)
    offset_val <- offset_val + 1
    
    n_strata <- length(unique(result_list$sol_by_cell))
    
    image(x = goa_ras,
          axes = F, ann = F,
          asp = 1,
          col = colorRampPalette(
            colors = brewer.pal(n = 10, name = 'Paired'))(n_strata),
          add = T)
    
    ## Species label
    text(x = min(grid_goa$E_km) + x_range * 0.72,
         y = min(grid_goa$N_km) + offset_y + y_range * 0.6,
         labels =  c(rev(common_names_eval), rev(common_names_opt))[ispp],
         col = ifelse(c(rev(common_names_eval),
                        rev(common_names_opt))[ispp] %in% common_names_opt, 
                      "black", 
                      "darkgrey"),
         cex = 0.5)
    
    ## Simulate station locations from a survey
    temp_samples <- c()
    temp_solution <- result_list$sol_by_cell
    temp_allocation <- result_list$sample_allocations[, isample]
    
    for (istrata in 1:n_strata) {
      temp_samples = c(temp_samples,
                       sample(x = which(temp_solution == istrata),
                              size = temp_allocation[istrata])
      )
    }
    
    ## Plot station locations
    temp_loc <- grid_goa[temp_samples, c("E_km", "N_km")]
    temp_loc$N_km <- temp_loc$N_km + offset_y
    points(x = temp_loc,
           pch = 16,
           cex = 0.2)
    
  }
  
  ## Close Device
  dev.off()
}
