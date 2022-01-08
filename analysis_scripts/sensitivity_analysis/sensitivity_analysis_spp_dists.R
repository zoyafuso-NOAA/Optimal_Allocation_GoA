###############################################################################
## Project:       Sensitivity Analysis for Survey Optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
VAST_sim_data_dir <- "F:/VAST_Runs/"
github_dir <- getwd()

library(VAST)
library(SamplingStrata)
library(raster)
library(sp)
library(RColorBrewer)

##################################################
####   Load data
##################################################
load("data/processed/optimization_data.RData")
load("data/processed/prednll_VAST_models.RData")
load("data/processed/grid_goa.RData")

##################################################
####   Result object
##################################################
dens_vals <- array(dim = c(n_cells, ns_opt, n_years, 5),
                   dimnames = list(NULL, common_names_opt, NULL, 
                                   c(paste0("Type ", 0:4))))

##################################################
####   Loop over species and simulate data with the four different options
####   in the FishStatsUtils::simulate_data() function
##################################################

for (ispp in common_names_opt) {
  
  pred_jnll_score <- subset(x = pred_jnll, 
                            spp_name == ispp, 
                            select = c("FALSE", "TRUE"))
  depth_in_model <- as.logical(names(which.min(pred_jnll_score)))
  
  temp_VAST_dir <- paste0(VAST_sim_data_dir, ispp, 
                          ifelse(depth_in_model == TRUE,
                                 yes = "_depth", no = ""), "/")
  
  load(paste0(temp_VAST_dir, "fit_sim.RData"))
  
  dyn.load(paste0(temp_VAST_dir, "VAST_v12_0_0.dll"))
  
  dens_vals[, ispp, , "Type 0"] <- fit_sim$Report$D_gct[, 1, years_included]
  print("Finished with Type 0")
  
  for (itype in 1:4) {
    sim_val <- FishStatsUtils::simulate_data(fit = fit_sim,
                                             type = itype,
                                             random_seed = 123423)
    
    pred_TF <- (1:length(sim_val$b_i))[-c(1:7900)]
    
    temp_density <- matrix(data = sim_val$b_i[pred_TF],
                           nrow = n_cells,
                           ncol = n_years)
    
    ## Temporary input objects for the survey simulation
    temp_density <- sweep(x = temp_density,
                          MARGIN = 1,
                          STATS = grid_goa$Area_km2,
                          FUN = "/")
    
    dens_vals[, ispp, , paste0("Type ", itype)] <- temp_density
    
    print(paste("Finished with Type", itype))
  }
  
  dyn.unload(paste0(temp_VAST_dir, "VAST_v12_0_0.dll"))
  
}

########################################
## Predicted Density
########################################
yrange_diff <- diff(range(grid_goa$N_km))

{
  pdf("C:/Users/Zack Oyafuso/Desktop/sensitivity_sim_data.pdf", 
      width = 8, height = 11, onefile = TRUE)
  
  par(mar = c(0, 0, 0, 0), mfcol = c(5, 3), oma = c(0, 2, 2, 0))
  for(ispp in common_names_opt) {
    for (itype in 0:4){ 
      ## Base Layer
      plot(1, type = "n", asp = 1,
           xlim = range(grid_goa$E_km),
           ylim = with(grid_goa,
                       c(min(N_km) + 0.0 * yrange_diff,
                         max(N_km) + 1.25 * yrange_diff)),
           axes = F, ann = F)
      box()
      
      ## Subtitle
      if(itype == 0) {
        mtext(side = 3,
              text = ispp,
              line = 0,
              cex = 1)
      }
      if(ispp %in% common_names_opt[c(1,4,7,10,13)] ) {
        mtext(side = 2, text = paste("Type", itype))
      }

      
      ## Calculate quantiles of the density distribution
      vals  = dens_vals[, ispp, , paste0("Type ", itype)]
      val_cuts = c(0,quantile(x = vals[vals > 10],
                              probs = seq(0, 1, length = 9) ))
      
      #Add legend
      val_cuts_legend = round(val_cuts[-1])
      colors = c("lightgrey", brewer.pal(n = 7, name = "Oranges"), "black")
      legend(x = x_range[1],
             y = y_range[1] - yrange_diff * 0.15,
             fill = colors,
             bty = "n",
             ncol = 3,
             cex = 0.7,
             legend = c("< 10", paste0("10-", val_cuts_legend[2]),
                        paste0(val_cuts_legend[2:(length(val_cuts_legend)-1)], "-",
                               val_cuts_legend[3:length(val_cuts_legend)])) )
      
      ## Loop over years and plot spatial distributions
      for (iyear in c(1, 5, 11)) {
        
        #Extract density values for a species in a year,
        vals <- dens_vals[, ispp, iyear, paste0("Type ", itype)]
        # vals  = report$D_gct[, 1, years_included[iyear]]
        
        #plot density
        goa = sp::SpatialPointsDataFrame(
          coords = grid_goa[, c("E_km", "N_km")],
          data = data.frame(density = vals) )
        goa_ras = raster::raster(x = goa,
                                 resolution = 10)
        goa_ras = raster::rasterize(x = goa,
                                    y = goa_ras,
                                    field = "density")
        
        #Discretize into quantiles
        values(goa_ras) = cut(x = values(goa_ras),
                              breaks = val_cuts)
        
        offset_y <- 0.6 * yrange_diff * (c("1" = 0, "5" = 1, "11" = 2)[paste0(iyear)] )
        goa_ras <- raster::shift(goa_ras,
                                 dy = offset_y )
        
        #lay image
        image(x = goa_ras,
              asp = 1,
              axes = F,
              ann = F,
              add = T,
              col = colors)
        
        #Year label
        text(x = goa_ras@extent[1] + 0.7 * diff(goa_ras@extent[1:2]),
             y = goa_ras@extent[3]+ 0.7 * diff(goa_ras@extent[3:4]),
             labels = year_set[years_included[iyear]],
             cex = 1)
      }
    }
    
  }
  
  dev.off()
}
