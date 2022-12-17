###############################################################################
## Project:       Sensitivity: Observation Error
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Compare Relative Bias across districts across species
##                for current and optimized STRS surveys
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(raster)

##################################################
####  Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[2]
output_dir <- paste0(c("/Users/zackoyafuso/",
                       "C:/Users/Zack Oyafuso/")[which_machine],
                     "Google Drive/MS_Optimizations/TechMemo/figures/")

##################################
## Import plotting percentile time series function
##################################
source("modified_functions/plot_percentiles.R")

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load('data/processed/optimization_data.RData')
load('data/processed/grid_goa.RData')

##################################
## Load Survey Simulations for the current optimization (scenario L),
## full-scale optimization with 15 strata (scenario A) and 
## area-level optimization with 5 strata per area (scenario B)
##################################
main_rb <- array(dim = c(3, n_years, ns_all, n_districts, n_iters),
                 dimnames = list(c("current", "full_domain", "district"),
                                 paste0("year_", 1:n_years),
                                 common_names_all, 
                                 paste0("district_", 1:n_districts), 
                                 NULL) )

load(paste0("results/scenario_L/Multispecies_Optimization/Str_current/",
            "simulation_results.RData"))
main_rb["current", , , , ] <- STRS_rel_bias_est[, , "boat_2", ]

load(paste0("results/scenario_A/Multispecies_Optimization/Str_15/",
            "simulation_results.RData"))
main_rb["full_domain", , , , ] <- STRS_rel_bias_est[, , "boat_2", ]

load(paste0("results/scenario_B/Multispecies_Optimization/Str_5/",
            "simulation_results.RData"))
main_rb["district", , , , ] <- STRS_rel_bias_est[, , "boat_2", ]

##################################
## General layout of plots
##################################
gen_layout <- rbind(matrix(rep(17, 6), nrow = 1),
                    cbind(matrix(data = 1:15,
                                 ncol = 5,
                                 byrow = T),
                          matrix(rep(16, 3), ncol = 1)))

legend_layout <- matrix(c(7, 1, 1, 1,
                          7, 1, 1, 1,
                          7, 2, 2, 2,
                          7, 3:5,
                          7, 3:5,
                          7, 6, 6, 6), nrow = 4)

##################################################
####   This plot is split into four pages. Set the species indices for
####   each page.
##################################################
spp_idx_1 <- spp_idx_opt[1:7]
spp_idx_2 <- spp_idx_opt[8:14]
spp_idx_3 <- c(spp_idx_opt[15], spp_idx_eval[1:6])
spp_idx_4 <- spp_idx_eval[7:11]

for (which_spp in paste(1:4)) {
  
  ## Open png file
  # png(filename = paste0(output_dir, "Figure07_RB_district_", which_spp, ".png"),
  #     width = 190,
  #     height = 220,
  #     units = "mm",
  #     res = 500)
  
  ## Set up plot layout
  par(mar = c(0.25, 0, 0.25, 0), oma = c(1.25, 5, 0, 0))
  plot_layout <- rbind( cbind(gen_layout + 17*0, gen_layout + 17*1),
                        cbind(gen_layout + 17*2, gen_layout + 17*3),
                        switch(paste0(which_spp %in% paste0(1:3)),
                               "TRUE" = rbind(cbind(gen_layout + 17*4, 
                                                    gen_layout + 17*5),
                                              cbind(gen_layout + 17*6, 
                                                    legend_layout + 17*7)),
                               "FALSE" = cbind(gen_layout + 17*4, 
                                               legend_layout + 17*5))
                        
  )
  
  layout(mat =  plot_layout, 
         widths = c(rep(1, 11), 0.1),
         heights = c(0.75, 1, 1, 1,
                     0.75, 1, 1, 1,
                     0.75, 1, 1, 1,
                     0.75, 1, 1, 1))
  
  ## Loop over species
  for (ispp in get(paste0("spp_idx_", which_spp)) ) {
    
    ## Calculate y-max for the plot: for each survey type (MARGIN = 1),
    ## calculate the 95% percentile for each year and calculate the max
    ## absolute value of that output
    y_max <- max(abs(apply(main_rb[, , ispp, ,],
                           MARGIN = 1, 
                           FUN = function(x) 
                             plot_percentiles(values = x, plot = F))))
    
    y_max <- ifelse(y_max * 1.5 < 0.35, 0.35, y_max * 1.5)
    
    for (itype in dimnames(main_rb)[[1]] ) {  ## Loop over design -- start
      for (idistrict in 1:5) { ## Loop over district -- start
        
        ## Base plot
        plot(1,
             type = "n",
             xlim = c(0, 12),
             ylim = c(-y_max, y_max),
             axes = F,
             ann = F)
        box()
        
        ## Color of the background determines the survey type
        rect(xleft = par("usr")[1], 
             xright = par("usr")[2], 
             ybottom = par("usr")[3], 
             ytop = par("usr")[4], 
             col = switch(itype,
                          "current" = "white",
                          "full_domain" = "grey90",
                          "district" = "grey50"))
        

        ## Time axis
        axis(side = 1, 
             at = 1:11, 
             labels = NA,
             tck = -0.05)
        
        if (itype == "district") {
          if(idistrict %in% c(1, 3, 5) ) {
            axis(side = 1, 
                 at = c(1, 11), 
                 labels = c("Yr 1", "Yr 11"),
                 lwd = F, 
                 tick = F, 
                 line = -0.5, 
                 cex.axis = 0.75)
          }
        }
        
        
        ## Plot time series
        plot_this <- main_rb[itype,
                             ,
                             ispp,
                             idistrict,
                             ]
        
        plot_percentiles(values = plot_this,
                         xs = 1:11, 
                         pt.cex = 0.25,
                         pt.colors = "black")
        
        if(idistrict == 1) {
          axis(side = 2,
               las = 1,
               at = pretty(c(-y_max, y_max), n = 3),
               cex.axis = 0.7,
               tck = -0.1)
        }
        
        ## District Labels
        if(itype == "current") mtext(side = 3, 
                                     line = -0.75,
                                     text = districts$district[idistrict],
                                     cex = 0.5)
        
        abline(h = 0, lwd = 0.5, lty = "dotted")
      } ## Loop over district -- end
    } ## Loop over design scenario -- end
    
    ## Species Label
    plot(1,type = "n", axes = F, ann = F)  
    plot(1,type = "n", axes = F, ann = F, xlim = c(0, 1), ylim = c(0, 1))
    text(x = 0.4, 
         y = 0.35, 
         labels = common_names_all[ispp], 
         cex = 1.4, 
         font = 2)
  }
  
  par(mar = c(0.25, 0, 0.25, 0))
  goa <- sp::SpatialPointsDataFrame(
    coords = grid_goa[, c("E_km", "N_km")],
    data = data.frame(Str_no = as.integer(district_vals) ) )
  goa_ras <- raster::raster(x = goa, 
                            resolution = 10)
  goa_ras <- raster::rasterize(x = goa, 
                               y = goa_ras, 
                               field = "Str_no")
  
  raster::image(goa_ras, asp = 1, axes = F, col = palette()[-1])
  segments(x0 = districts$W_UTM,
           x1 = districts$E_UTM,
           y0 = tapply(X = grid_goa$N_km,
                       INDEX = as.integer(district_vals),
                       FUN = min) + c(-50, 500, -100, 200, -100),
           xpd = NA,
           lwd = 2)
  
  text(x = rowMeans(cbind(districts$W_UTM, districts$E_UTM)),
       y = tapply(X = grid_goa$N_km,
                  INDEX = district_vals,
                  FUN = min) + c(-400, 450, -450, 150, -450),
       cex = 0.90,
       labels = districts$district,
       pos = 3,
       xpd = NA)
  
  plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = F)
  # box()
  rect(xleft = 0.1, xright = 0.3,
       ybottom = 0.05, ytop = 0.95,
       col = "dodgerblue",
       border = F)
  rect(xleft = 0.1, xright = 0.3,
       ybottom = 0.25, ytop = 0.75,
       col = "cadetblue1",
       border = F)
  points(0.2, 0.5, pch = 15)
  text(x = c(0.7),
       y = c(0.05, 0.25, 0.5, 0.75, 0.95),
       labels = paste0(c(5, 25, 50, 75, 95), "%"))
  
  plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = F)
  box()
  text(0.5, 0.5, "Existing\nSTRS Design")
  
  plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = F)
  rect(xleft = -2, xright = 2,
       ybottom = -2, ytop = 2,
       col = "grey90")
  box()
  text(0.5, 0.5, "Gulf-Wide\n(15 Strata)\nProposed STRS Design")
  
  plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = F)
  rect(xleft = -2, xright = 2,
       ybottom = -2, ytop = 2,
       col = "grey50")
  box()
  text(0.5, 0.5, "Area-Level\n(5 Strata per Area)\nProposed STRS Design")
  
  mtext(side = 2, 
        outer = T, 
        text = "Bias Ratio (Simulated Value / True Value)", 
        line = 3, 
        font = 2)
  
  # dev.off()
}
