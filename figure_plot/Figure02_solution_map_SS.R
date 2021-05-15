###############################################################################
## Project:         Plot Single Species Optimization Solutions
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     
###############################################################################
rm(list = ls())

##################################################
####  Set up directories  
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")

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
load(paste0(github_dir, "data/optimization_data.RData"))
load(paste0(github_dir, "/data/Extrapolation_depths.RData"))
load(paste0(github_dir, "results/full_domain/Single_Species_Optimization/", 
            "optimization_knitted_results.RData"))

##################################################
####   Set Constants
##################################################
idomain <- "full_domain"
isample <- 2 #1, 2, or 3 boat solution
x_range <- diff(range(Extrapolation_depths$E_km))
y_range <- diff(range(Extrapolation_depths$N_km))

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
    if (ispp %in% c(ns_all/2 + 1, 1)) { 
      
      ## Base Plot
      plot(x = 1, y = 1, 
           type = "n", axes = F, ann = F,
           xlim = range(Extrapolation_depths$E_km),
           ylim = with(Extrapolation_depths, 
                       c(min(N_km), max(N_km) + 7 * y_range)),
           asp = 1)
      
      ## Since we are plotting multiple maps on a plot, we shift the plots
      ## using an offset value which is initialized at 0
      offset_val <- 0
    }
    
    ## Which index to plot
    idx = which(settings$boat == isample &
                  settings$species == c(rev(common_names_eval),
                                        rev(common_names_opt))[ispp])
    
    ## Set up spatial object to plot solution
    goa <- SpatialPointsDataFrame(
      coords = Extrapolation_depths[, c("E_km", "N_km")],
      data = data.frame("Str_no" = res_df[, idx],
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
    
    n_strata <- length(unique(res_df[, idx]))
    
    image(x = goa_ras,
          axes = F, ann = F,
          asp = 1,
          col = colorRampPalette(
            colors = brewer.pal(n = 10, name = 'Paired'))(n_strata),
          add = T)
    
    ## Species label
    text(x = min(Extrapolation_depths$E_km) + x_range * 0.72,
         y = min(Extrapolation_depths$N_km) + offset_y + y_range * 0.6,
         labels =  c(rev(common_names_eval), rev(common_names_opt))[ispp],
         col = ifelse(c(rev(common_names_eval),
                        rev(common_names_opt))[ispp] %in% common_names_opt, 
                      "black", 
                      "darkgrey"),
         cex = 0.5)
    
    ## Simulate station locations from a survey
    temp_samples <- c()
    temp_strata <- nrow(strata_list[[idx]])
    temp_solution <- res_df[, idx ]
    temp_allocation <- strata_list[[idx]]$Allocation
    
    for (istrata in 1:temp_strata) {
      temp_samples = c(temp_samples,
                       sample(x = which(temp_solution == istrata),
                              size = temp_allocation[istrata])
      )
    }
    
    ## Plot station locations
    temp_loc <- Extrapolation_depths[temp_samples, c("E_km", "N_km")]
    temp_loc$N_km <- temp_loc$N_km + offset_y
    points(x = temp_loc,
           pch = 16,
           cex = 0.2)
    
  }
  
  ## Close Device
  dev.off()
}
