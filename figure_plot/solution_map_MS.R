###############################################################################
## Project:         Plot Multispecies Optimization Solutions
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     
###############################################################################
rm(list = ls())

##################################################
####  Install a forked version of the SamplingStrata Package from 
####  zoyafuso-NOAA's Github page
####
####  Import other required packages
##################################################
library(sp)
library(RColorBrewer)
library(raster)

##################################################
####   Set up directories based on whether the optimization is being conducted
####        on a multi-species or single-species level
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[1]

github_dir <- paste0(c("/Users/zackoyafuso/Documents", 
                       "C:/Users/Zack Oyafuso/Documents")[which_machine],
                     "/GitHub/Optimal_Allocation_GoA/")

output_dir <- paste0(c("/Users/zackoyafuso/Google Drive/",
                       "C:/Users/Zack Oyafuso/Google Drive/")[which_machine],
                     "MS_Optimizations/TechMemo/figures/")

##################################################
####   Load Data
####   Load Population CVs for use in the thresholds
##################################################
load(paste0(github_dir, "data/optimization_data.RData"))
load(paste0(github_dir, "data/Extrapolation_depths.RData"))
load(paste0(github_dir, "results/MS_optimization_knitted_results.RData"))

districts$E_UTM[4] <- 1240
districts$W_UTM[5] <- 1240

##################################################
####   Constants
##################################################
x_range <- diff(range(Extrapolation_depths$E_km))
y_range <- diff(range(Extrapolation_depths$N_km))

##################################################
####    Plot
##################################################
{
  
  for(plot_stations in c(T, F)) {
    png(filename = paste0(output_dir,
                          "MS_solutions",
                          ifelse(plot_stations, "_with_stations", ""),
                          ".png"),
        width = 150,
        height = 200,
        units = "mm",
        res = 500)
    
    layout(mat = matrix(data = c(1, 4,
                                 2, 5,
                                 3, 6,
                                 7, 7), ncol = 2, byrow = T),
           heights = c(1,1,1,0.4))
    
    par(mar = c(0, 0, 0, 0))
    
    for (idomain in c("district", "full_domain")) {
      istrata <- list("district" = c(3, 5, 10),
                      "full_domain" = c(10, 15, 20))[[idomain]][1]
      
      for(istrata in list("district" = c(3, 5, 10),
                          "full_domain" = c(10, 15, 20))[[idomain]]) {
        #Base Plot layer
        plot(1, 
             type = "n",
             xlim = range(Extrapolation_depths$E_km),
             ylim = with(Extrapolation_depths, 
                         c(min(N_km), 
                           max(N_km) + 2.1 * y_range)), 
             axes = F,
             ann = F,
             asp = 1)
        
        
        #Strata label
        mtext(side = 3,
              line = -2,
              cex = 1.25,
              font = 2,
              text = paste(istrata, 
                           c("full_domain" = "Strata", 
                             "district" = "Strata Per District")[idomain] ))
        
        box()
        
        for (iboat in 1:n_boats) {
          sol_idx <- paste0("sol_", (subset(x = settings, 
                                            select = id,
                                            subset = domain == idomain & 
                                              strata == istrata & 
                                              boat == iboat)))
          
          goa <- sp::SpatialPointsDataFrame(
            coords = Extrapolation_depths[, c("E_km", "N_km")],
            data = data.frame(Str_no = res_df[, sol_idx]) )
          goa_ras <- raster::raster(x = goa, 
                                    resolution = 10)
          goa_ras <- raster::rasterize(x = goa, 
                                       y = goa_ras, 
                                       field = "Str_no")
          offset_y <- 0.9 * y_range * (iboat - 1)
          goa_ras <- raster::shift(goa_ras, dy = offset_y )
          
          n_strata <- nrow(strata_list[[sol_idx]])
          strata_colors <- 
            colorRampPalette(brewer.pal(n = 11, 
                                        name = "Paired"))(n_strata)
          image(goa_ras,
                add = TRUE,
                col =  sample(strata_colors))
          
          #Sample Size label
          text(x = min(Extrapolation_depths$E_km) + x_range*0.15,
               y = min(Extrapolation_depths$N_km) + offset_y + y_range*0.55,
               label = paste("n =", samples[iboat]),
               cex = 1.5,
               font = 2)
          
          rect(xleft = districts$W_UTM,
               xright = districts$E_UTM,
               ybottom = tapply(X = Extrapolation_depths$N_km,
                                INDEX = district_vals,
                                FUN = min) + offset_y,
               ytop = tapply(X = Extrapolation_depths$N_km,
                             INDEX = district_vals,
                             FUN = max) + offset_y)
          
          if (plot_stations) {
            #Simulate a sample solution
            temp_samples <- c()
            temp_strata <- nrow(strata_list[[sol_idx]])
            temp_solution <- res_df[, sol_idx]
            temp_allocation <- strata_list[[sol_idx]]$Allocation
            
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
    
    ## Districts Legend
    goa <- sp::SpatialPointsDataFrame(
      coords = Extrapolation_depths[, c("E_km", "N_km")],
      data = data.frame(Str_no = res_df[, sol_idx]) )
    goa_ras <- raster::raster(x = goa, 
                              resolution = 10)
    goa_ras <- raster::rasterize(x = goa, 
                                 y = goa_ras, 
                                 field = "Str_no")
    
    par(mar = c(0,0,0.5,0))
    image(goa_ras,
          col = "lightgrey",
          asp = 1,
          axes = F,
          ann = F)
    
    rect(xleft = districts$W_UTM,
         xright = districts$E_UTM,
         ybottom = tapply(X = Extrapolation_depths$N_km,
                          INDEX = district_vals,
                          FUN = min),
         ytop = tapply(X = Extrapolation_depths$N_km,
                       INDEX = district_vals,
                       FUN = max) )
    
    text(x = tapply(X = Extrapolation_depths$E_km, 
                    INDEX = district_vals,
                    FUN = mean),
         y = tapply(X = Extrapolation_depths$N_km, 
                    INDEX = district_vals,
                    FUN = mean),
         labels = districts$district,
         font = 2,
         cex = 0.8)
    
    dev.off()
  }
  
}
