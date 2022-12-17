###############################################################################
## Project:         Plot Multispecies Optimization Solutions
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     
###############################################################################
rm(list = ls())

##################################################
####  Import required packages
##################################################
library(sp)
library(RColorBrewer)
library(raster)

##################################################
####   Set up directories based on whether the optimization is being conducted
####        on a multi-species or single-species level
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[2]
output_dir <- paste0(c("/Users/zackoyafuso/Google Drive/",
                       "C:/Users/Zack Oyafuso/Google Drive/")[which_machine],
                     "MS_Optimizations/TechMemo/figures/")

##################################################
####   Load Data
####   Load Population CVs for use in the thresholds
##################################################
load("data/processed/optimization_data.RData")
load("data/processed/grid_goa.RData")

seed <- 234233
districts$E_UTM[4] <- 1240
districts$W_UTM[5] <- 1240

##################################################
####    Plot
##################################################
{
  
  ## Open device
  png(filename = paste0(output_dir, "Figure04_MS_solutions.png"),
      width = 170, height = 170, units = "mm", res = 500)

  ## Panel Layout
  layout(mat = matrix(data = c(1, 2,
                               3, 4,
                               5, 5), ncol = 2, byrow = T),
         heights = c(1, 1, 0.4))
  par(mar = c(0, 0, 0, 0))
  
  for (iscen in c("A", "B")) {
    
    for(istrata in list("B" = c(3, 5),
                        "A" = c(10, 15))[[iscen]]) {
      ## Base Plot layer
      plot(x = 1, y = 1, asp = 1,
           type = "n", axes = F, ann = F,
           xlim = range(grid_goa$E_km),
           ylim = with(grid_goa, 
                       c(min(N_km), 
                         max(N_km) + 2.1 * y_range)))
      box()
      
      #Strata label
      mtext(side = 3,
            line = -2,
            cex = 1.25,
            font = 2,
            text = paste(c("3" = "C)", "5" = "D)", 
                           "10" = "A)", "15" = "B)")[paste(istrata)],
                         istrata, 
                         c("A" = "Strata", 
                           "B" = "Strata Per Area")[iscen] ))
      
      for (iboat in 1:n_boats) {
        
        ## Load Solution
        load(paste0("results/scenario_", iscen, "/Multispecies_Optimization/", 
                    "Str_", istrata, "/boat_", iboat, "/result_list.RData"))
        
        ## Set up spatial object to plot solution
        goa <- sp::SpatialPointsDataFrame(
          coords = grid_goa[, c("E_km", "N_km")],
          data = data.frame(Str_no = result_list$sol_by_cell) )
        goa_ras <- raster::raster(x = goa, 
                                  resolution = 10)
        goa_ras <- raster::rasterize(x = goa, 
                                     y = goa_ras, 
                                     field = "Str_no")
        offset_y <- 0.9 * y_range * (iboat - 1)
        goa_ras <- raster::shift(x = goa_ras, dy = offset_y )
        
        ## Set strata colors
        n_strata <- length(unique(result_list$sol_by_cell))
        strata_cols <- colorRampPalette(
          c(brewer.pal(n = 11, name = "Paired"),
            brewer.pal(n = 11, name = "Spectral")))(n_strata)
        
        set.seed(seed)
        strata_cols <- sample(x = strata_cols, replace = F)
        
        ## Plot solution
        image(x = goa_ras,
              add = TRUE,
              col = strata_cols)
        
        ## Sample Size label
        text(x = min(grid_goa$E_km) + x_range * 0.15,
             y = min(grid_goa$N_km) + offset_y + y_range * 0.55,
             label = paste("n =", samples[iboat]),
             cex = 1.5,
             font = 2)
        
        ## District boxes
        rect(xleft = districts$W_UTM,
             xright = districts$E_UTM,
             ybottom = tapply(X = grid_goa$N_km,
                              INDEX = district_vals,
                              FUN = min) + offset_y,
             ytop = tapply(X = grid_goa$N_km,
                           INDEX = district_vals,
                           FUN = max) + offset_y)
        
        
        ## Simulate a sample solution
        temp_samples <- c()
        temp_solution <- result_list$sol_by_cell
        temp_allocation <- result_list$sample_allocations[, iboat]
        
        for (temp_istrata in 1:n_strata) {
          temp_samples = c(temp_samples,
                           sample(x = which(temp_solution == temp_istrata),
                                  size = temp_allocation[temp_istrata]) )
        }
        
        ## Plot stations
        temp_loc <- grid_goa[temp_samples, c("E_km", "N_km")]
        temp_loc$N_km <- temp_loc$N_km + offset_y
        points(x = temp_loc,
               pch = 16,
               cex = 0.3)
        
      }
    }
  }
  
  ## Districts Legend
  goa <- sp::SpatialPointsDataFrame(
    coords = grid_goa[, c("E_km", "N_km")],
    data = data.frame(Str_no =result_list$sol_by_cell) )
  goa_ras <- raster::raster(x = goa, 
                            resolution = 10)
  goa_ras <- raster::rasterize(x = goa, 
                               y = goa_ras, 
                               field = "Str_no")
  
  par(mar = c(0,0,0.5,0))
  image(x = goa_ras,
        col = "lightgrey",
        asp = 1,
        axes = F,
        ann = F)
  
  rect(xleft = districts$W_UTM,
       xright = districts$E_UTM,
       ybottom = tapply(X = grid_goa$N_km,
                        INDEX = district_vals,
                        FUN = min),
       ytop = tapply(X = grid_goa$N_km,
                     INDEX = district_vals,
                     FUN = max) )
  
  text(x = tapply(X = grid_goa$E_km, 
                  INDEX = district_vals,
                  FUN = mean),
       y = tapply(X = grid_goa$N_km, 
                  INDEX = district_vals,
                  FUN = mean),
       labels = districts$district,
       font = 2,
       cex = 0.8)
  
  ## Close Device
  dev.off()
}
