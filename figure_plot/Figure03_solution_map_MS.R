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

seed <- 234233
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
  
  ## Open device
  png(filename = paste0(output_dir, "Figure03_MS_solutions.png"),
      width = 170, height = 170, units = "mm",
      res = 500)
  
  layout(mat = matrix(data = c(1, 2,
                               3, 4,
                               5, 5), ncol = 2, byrow = T),
         heights = c(1, 1, 0.4))
  par(mar = c(0, 0, 0, 0))
  
  for (idomain in c("district", "full_domain")) {
    
    istrata <- list("district" = c(3, 5),
                    "full_domain" = c(10, 15))[[idomain]][1]
    
    for(istrata in list("district" = c(3, 5),
                        "full_domain" = c(10, 15))[[idomain]]) {
      ## Base Plot layer
      plot(x = 1, y = 1, asp = 1,
           type = "n", axes = F, ann = F,
           xlim = range(Extrapolation_depths$E_km),
           ylim = with(Extrapolation_depths, 
                       c(min(N_km), 
                         max(N_km) + 2.1 * y_range)))
      box()
      
      #Strata label
      mtext(side = 3,
            line = -2,
            cex = 1.25,
            font = 2,
            text = paste(c("3" = "A)", "5" = "B)", 
                           "10" = "C)", "15" = "D)")[paste(istrata)],
                         istrata, 
                         c("full_domain" = "Strata", 
                           "district" = "Strata Per District")[idomain] ))
      
      for (iboat in 1:n_boats) {
        
        ## Find index that corresponds to the solution scenario
        sol_idx <- with(settings, 
                        which(domain == idomain & 
                                strata == istrata & 
                                boat == iboat))
        
        ## Set up spatial object to plot solution
        goa <- sp::SpatialPointsDataFrame(
          coords = Extrapolation_depths[, c("E_km", "N_km")],
          data = data.frame(Str_no = res_df[, sol_idx]) )
        goa_ras <- raster::raster(x = goa, 
                                  resolution = 10)
        goa_ras <- raster::rasterize(x = goa, 
                                     y = goa_ras, 
                                     field = "Str_no")
        offset_y <- 0.9 * y_range * (iboat - 1)
        goa_ras <- raster::shift(x = goa_ras, dy = offset_y )
        
        ## Set strata colors
        n_strata <- nrow(strata_list[[sol_idx]])
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
        text(x = min(Extrapolation_depths$E_km) + x_range * 0.15,
             y = min(Extrapolation_depths$N_km) + offset_y + y_range * 0.55,
             label = paste("n =", samples[iboat]),
             cex = 1.5,
             font = 2)
        
        ## District boxes
        rect(xleft = districts$W_UTM,
             xright = districts$E_UTM,
             ybottom = tapply(X = Extrapolation_depths$N_km,
                              INDEX = district_vals,
                              FUN = min) + offset_y,
             ytop = tapply(X = Extrapolation_depths$N_km,
                           INDEX = district_vals,
                           FUN = max) + offset_y)
        
        
        ## Simulate a sample solution
        temp_samples <- c()
        temp_strata <- nrow(strata_list[[sol_idx]])
        temp_solution <- res_df[, sol_idx]
        temp_allocation <- strata_list[[sol_idx]]$Allocation
        
        for (temp_istrata in 1:temp_strata) {
          temp_samples = c(temp_samples,
                           sample(x = which(temp_solution == temp_istrata),
                                  size = temp_allocation[temp_istrata]) )
        }
        
        ## Plot stations
        temp_loc <- Extrapolation_depths[temp_samples, c("E_km", "N_km")]
        temp_loc$N_km <- temp_loc$N_km + offset_y
        points(x = temp_loc,
               pch = 16,
               cex = 0.3)
        
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
  image(x = goa_ras,
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
