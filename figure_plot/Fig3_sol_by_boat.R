###############################################################################
## Project:      What do the stratifications look like?
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:  Figure 3: Representative examples of strata boundary maps 
##               arising from solutions for the species-specific CV constraint 
##               optimization for five, ten, and fifteen strata across the 
##               three effort scenarios. The colors represent distinct strata. 
###############################################################################
rm(list = ls())

##################################################
####   Import required packages
##################################################
library(sp)
library(raster)
library(RColorBrewer)

############################
## Set up directories
#############################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[2]

VAST_model <- "6g"
github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/")

figure_dir <- paste0(c("/Users/zackoyafuso/Google Drive/", 
                       "C:/Users/Zack Oyafuso/Google Drive/")[which_machine],
                     "MS_Optimizations/figure_plot/")

########################
## Load Data
########################
load(paste0(dirname(github_dir), "/data/Extrapolation_depths.RData" ))
load(paste0(github_dir, "optimization_data.RData"))
load(paste0(github_dir, "Spatiotemporal_Optimization_Scheme2/",
            "spatiotemporal_Flexible_optimization_results.RData"))

########################
## Constants
########################
samples <- c(820, 550, 280)
stratas <- c(5, 10, 15)
yrange <- diff(range(Extrapolation_depths[, "N_km"]))
plot_random_sample <- F

{
  png(file = paste0(figure_dir, "Fig3_sol_by_boat",
                    ifelse(plot_random_sample == T, "_withsamples", ""),
                    ".png"),
      width = 240, 
      height = 100, 
      units = "mm", 
      res = 1000)
  
  ####################################
  ## Plot Settings
  ####################################
  par(mfrow = c(1, 3), 
      mar = c(0, 0, 3, 0))
  
  for(istrata in stratas){
    
    ##################################
    ## Empty Plot
    ##################################
    plot(1, 
         type = "n", 
         axes = F, 
         ann = F,
         xlim = range(Extrapolation_depths[, c("E_km")]),
         ylim = c(min(Extrapolation_depths[, c("N_km")]) - 1.5 * yrange,
                  max(Extrapolation_depths[, c("N_km")])))
    
    mtext(side = 3, 
          paste(istrata, "Strata"), 
          font = 2) 
    
    offset <- 0
    for(isample in samples){
      sub_settings <- subset(settings, 
                             nstrata == istrata)
      
      isol <- sub_settings$id[which.min(abs(sub_settings$n - isample))]
      
      #Plot Solution
      goa <- sp::SpatialPointsDataFrame(
        coords = Extrapolation_depths[, c("E_km", "N_km")],
        data = data.frame(stratum = res_df[, isol+1]) )
      goa_ras <- raster::raster(x = goa, 
                                resolution = 5)
      goa_ras <- raster::rasterize(x = goa, 
                                   y = goa_ras, 
                                   field = "stratum")
      
      goa_ras <- raster::shift(goa_ras, 
                               dy = offset)
      
      
      if(istrata != 15) sol_col = brewer.pal(n = istrata, 
                                             name = "Paired")
      if(istrata == 15) sol_col = c(brewer.pal(n = 12, 
                                               name = "Paired"),
                                    rev(brewer.pal(n = 11, 
                                               name = "Spectral"))[11-istrata] )
      
      image(goa_ras, 
            add = T,
            asp = 1, 
            axes = F, 
            col = sol_col)
      
      #Plot sample size label
      xrange <- diff(goa_ras@extent[1:2])
      yrange <- diff(goa_ras@extent[3:4])
      text(x = goa_ras@extent[1] + xrange * 0.70, 
           y = goa_ras@extent[3] + yrange * 0.60,
           labels = paste0("n = ", isample, "\n"), 
           cex = 1.5)
      
      if(plot_random_sample) {
        #Sample based on the stratification allocations
        sample_vec <- c()
        for (i in 1:istrata ){
          available_cells <- which(res_df[, isol+1] == i)
          
          if(length(available_cells) > 0){
            sample_cells <- sample(x = available_cells, 
                                   size = strata_list[[isol]]$Allocation[i], 
                                   replace = F)
            sample_vec <- c(sample_vec, sample_cells)
          }
          
        }
        
        #Organize sample set and total number of samples
        sample_vec <- sort(sample_vec)
        sample_pts <- Extrapolation_depths[sample_vec, c("E_km", "N_km")]
        sample_pts[, "N_km"] <- sample_pts[, "N_km"] + offset
        
        points(sample_pts, 
               pch = 16, 
               cex = 0.2)
      }
      
      #Update offset
      offset <- offset - yrange*.75
    }
  }
  dev.off()
}
