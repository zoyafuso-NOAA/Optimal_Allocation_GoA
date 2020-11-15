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
####  Set up directories  
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
####    Set up constants
##################################################
load(paste0(github_dir, "data/optimization_data.RData"))
load(paste0(github_dir, "/data/Extrapolation_depths.RData"))
load(paste0(github_dir, "results/Spatiotemporal_Optimization/",
            "optimization_knitted_results.RData"))

common_names <- gsub(common_names, pattern = " ", replacement = "\n")
common_names[11] <- "blackspotted/\nrougheye\nrockfishes"

x_range <- diff(range(Extrapolation_depths$E_km))
y_range <- diff(range(Extrapolation_depths$N_km))


##################################################
####    Plot
##################################################
{
  png(filename = paste0(output_dir, "MS_solutions.png"),
      width = 190,
      height = 200,
      units = "mm",
      res = 500)
  
  par(mfcol = c(3, 2),
      mar = c(0,0,0,0))
  
  for (sample_survey in c(FALSE, TRUE)) { 
    for (istrata in 2:4) {
      
      #Base Plot layer
      plot(1, 
           type = "n",
           xlim = range(Extrapolation_depths$E_km),
           ylim = with(Extrapolation_depths, 
                       c(min(N_km), 
                         max(N_km) + 1.5 * y_range)), 
           axes = F,
           ann = F,
           asp = 1)
      box()
      
      #Strata label
      mtext(side = 1,
            line = -1.5,
            cex = 1.25,
            font = 2,
            text = paste(stratas[istrata], "Strata"))
      
      for (iboat in 1:3) {
        #Which index to plot
        idx = settings$id[settings$boat == iboat & 
                            settings$strata == stratas[istrata]]
        
        #Plot Solution
        goa <- SpatialPointsDataFrame(
          coords = Extrapolation_depths[, c("E_km", "N_km")],
          data = data.frame(Str_no = res_df[, 1 + idx]) )
        
        goa_ras <- raster(goa, 
                          resolution = 5)
        goa_ras <- rasterize(x = goa, 
                             y = goa_ras, 
                             field = "Str_no")
        offset_y <- 0.70 * y_range * (iboat-1)
        goa_ras <- raster::shift(goa_ras, dy = offset_y )
        
        image(goa_ras, 
              axes = F,
              ann = F,
              asp = 1,
              col = c(brewer.pal(n = 12, 
                                 name = 'Paired'),
                      brewer.pal(n = 12, 
                                 name = 'Paired'))[1:nrow(strata_list[[idx]])],
              add = T)
        
        #Sample Size label
        text(x = min(Extrapolation_depths$E_km) + x_range*0.7,
             y = min(Extrapolation_depths$N_km) + offset_y + y_range*0.6,
             label = paste("n =", samples[iboat]),
             cex = 1.5)
        
        if (sample_survey) {
          #Simulate a sample solution
          temp_samples <- c()
          temp_strata <- nrow(strata_list[[idx]])
          temp_solution <- res_df[, idx + 1]
          temp_allocation <- strata_list[[idx]]$Allocation
          
          for (temp_istrata in 1:temp_strata) {
            temp_samples = c(temp_samples,
                             sample(x = which(temp_solution == temp_istrata),
                                    size = temp_allocation[temp_istrata]) )
          }
          
          temp_loc <- Extrapolation_depths[temp_samples, c("E_km", "N_km")]
          temp_loc$N_km <- temp_loc$N_km + offset_y
          
          points(temp_loc,
                 pch = 16,
                 cex = 0.3)
        }
      }
    }
  }
  dev.off()
}
