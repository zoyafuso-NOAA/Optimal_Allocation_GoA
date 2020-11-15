###############################################################################
## Project:         Plot Single Species Optimization Solutions
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     
###############################################################################
rm(list = ls())

##################################################
####  Import Libraries  
##################################################
library(sp)
library(raster)
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
####    Load predicted density and optimization results ----
##################################################
load(paste0(github_dir, "data/optimization_data.RData"))
common_names <- gsub(common_names, pattern = " ", replacement = "\n")
common_names[11] <- "blackspotted/\nrougheye\nrockfishes"
load(paste0(github_dir, "/data/Extrapolation_depths.RData"))
load(paste0(github_dir, "results/Single_Species_Optimization/optimization_knitted_results.RData"))

isample <- 2 #1, 2, or 3 boat solution

par(mfrow = c(1, 1), 
    mar = c(0, 0, 0, 0))

x_range <- diff(range(Extrapolation_depths$E_km))
y_range <- diff(range(Extrapolation_depths$N_km))

{png(filename = paste0(output_dir, "SS_solutions.png"),
     height = 160,
     width = 150,
     units = "mm",
     res = 500)
  
  par(mfrow = c(1,2),
      mar = c(0,0,0,0))
  
  for (ispp in 1:ns) {
    
    if (ispp%%8 == 1 ) plot(1, 
                            type = "n",
                            xlim = range(Extrapolation_depths$E_km),
                            ylim = with(Extrapolation_depths, 
                                        c(min(N_km) + 0.0*y_range, 
                                          max(N_km) + 5.5*y_range)), 
                            axes = F,
                            ann = F,
                            asp = 1)
    
    #Which index to plot
    idx = settings$id[settings$iboat == isample & settings$ispp == ispp]
    
    #Plot Solution
    goa <- SpatialPointsDataFrame(
      coords = Extrapolation_depths[,c("E_km", "N_km")],
      data = data.frame(Str_no = res_df[, 1 + idx]) )
    
    goa_ras <- raster(goa, 
                      resolution = 5)
    goa_ras <- rasterize(x = goa, 
                         y = goa_ras, 
                         field = "Str_no")
    offset_y <- 0.8 * y_range * (ispp-1)%%8
    goa_ras <- raster::shift(goa_ras, dy = offset_y )
    
    image(goa_ras, 
          axes = F,
          ann = F,
          asp = 1,
          col = brewer.pal(n = 10, name = 'Paired'),
          add = T)
    
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
    
    temp_loc <- Extrapolation_depths[temp_samples, c("E_km", "N_km")]
    temp_loc$N_km <- temp_loc$N_km + offset_y
    
    points(temp_loc, 
           pch = 16,
           cex = 0.3)
    
    text(x = min(Extrapolation_depths$E_km) + x_range * 0.725,
         y = min(Extrapolation_depths$N_km) + offset_y + y_range*0.55,
         labels = common_names[ispp],
         cex = 0.75)
  }
  
  dev.off()}
