###############################################################################
## Project:       Sensitivity: Observation Error
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Compare Relative Bias across districts across species
##                for current and optimized STRS surveys
###############################################################################
rm(list = ls())

library(raster)

##################################################
####  Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/",
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine],
                     "GitHub/Optimal_Allocation_GoA/")

output_dir <- paste0(c("/Users/zackoyafuso/",
                       "C:/Users/Zack Oyafuso/")[which_machine],
                     "Google Drive/MS_Optimizations/TechMemo/figures/")

##################################
## Import plotting percentile time series function
##################################
source(paste0(github_dir, "modified_functions/plot_percentiles.R"))

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load(paste0(github_dir, '/data/optimization_data.RData'))
load(paste0(github_dir, '/data/Extrapolation_depths.RData'))

scen <- data.frame(survey_type = c("cur", rep("opt", 4) ),
                   strata = c("cur", 3, 5, 10, 15),
                   domain = c("full_domain",
                              rep(c("district", "full_domain"), each = 2)))

for (irow in 1:nrow(scen)) {
  scen_name <- paste0("SUR_", scen$survey_type[irow], "_",
                      scen$domain[irow], "_STR_", scen$strata[irow], "_")
  file_name <- paste0(github_dir, "results/sim_dens_surveys/",
                      scen_name, "simulation_results.RData")
  
  load(file_name)
}

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

# layout(legend_layout)
plot_legend <- function() {
  goa <- sp::SpatialPointsDataFrame(
    coords = Extrapolation_depths[, c("E_km", "N_km")],
    data = data.frame(Str_no = as.integer(district_vals) ) )
  goa_ras <- raster::raster(x = goa, 
                            resolution = 10)
  goa_ras <- raster::rasterize(x = goa, 
                               y = goa_ras, 
                               field = "Str_no")
  
  raster::image(goa_ras, asp = 1, axes = F, col = palette()[-1])
  segments(x0 = districts$W_UTM,
           x1 = districts$E_UTM,
           y0 = tapply(X = Extrapolation_depths$N_km,
                       INDEX = as.integer(district_vals),
                       FUN = min) + c(-50, 500, -100, 200, -100),
           xpd = NA,
           lwd = 2)
  
  text(x = rowMeans(cbind(districts$W_UTM, districts$E_UTM)),
       y = tapply(X = Extrapolation_depths$N_km,
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
  text(0.5, 0.5, "Current\nSTRS Design")
  
  plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = F)
  rect(xleft = -2, xright = 2,
       ybottom = -2, ytop = 2,
       col = "grey90")
  box()
  text(0.5, 0.5, "Gulf-Wide\n(10 Strata)\nOpt. STRS Design")
  
  plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = F)
  rect(xleft = -2, xright = 2,
       ybottom = -2, ytop = 2,
       col = "grey50")
  box()
  text(0.5, 0.5, "District-Level\n(3 Strata per District)\nOpt. STRS Design")
  # plot(1, type = "n",  axes = F)

}



spp_idx_1 <- spp_idx_opt[1:7]
spp_idx_2 <- spp_idx_opt[8:14]
spp_idx_3 <- c(spp_idx_opt[15], spp_idx_eval[1:6])
spp_idx_4 <- spp_idx_eval[7:11]

for (which_spp in paste(1:4)) {
  
  ## Open png file
  png(filename = paste0(output_dir, "Figure07_RB_district_", which_spp, ".png"),
      width = 190,
      height = 220,
      units = "mm",
      res = 500)
  
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
    ## Loop over three survey scenarios
    ## 1) Current
    ## 4) Gulf-wide optiization, 10 strata
    ## 2) District-level optimiztion, 3 strata per district
    
    ## Calculate what should be the ylimits for each set of species plots
    plot_this <- list()
    for (irow in c(1, 4, 2) ) { 
      ## subset result object based on survey scenario
      scen_name <- paste0("SUR_", scen$survey_type[irow], "_",
                          scen$domain[irow], "_STR_", scen$strata[irow], "_")
      # rb_district <- get(paste0(scen_name, "log_rb_district"))
      rb_district <- get(paste0(scen_name, "log_rb_district"))
      
      for (idistrict in 1:5) {
        plot_this <- c(plot_this, 
                       list(plot_percentiles(values = rb_district[1:n_years, 
                                                                  ispp, 
                                                                  "boat_2",
                                                                  idistrict, 
                                                                  1:n_iters],
                                             plot = F) ))
      }
      
    }
    
    ymax_ <- max(unlist(plot_this), 1, na.rm = T) 
    ymin_ <- min(unlist(plot_this), -1, na.rm = T)
    
    for (irow in c(1, 4, 2) ) { ## Loop over design scenario -- start
      ## subset result object based on survey scenario
      scen_name <- paste0("SUR_", scen$survey_type[irow], "_",
                          scen$domain[irow], "_STR_", scen$strata[irow], "_")
      rb_district <- get(paste0(scen_name, "log_rb_district"))
      
      for (idistrict in 1:5) { ## Loop over district -- start
        
        ## Base plot
        plot(1,
             type = "n",
             xlim = c(0, 12),
             ylim = c(ymin_, ymax_),
             axes = F,
             ann = F)
        box()
        
        ## Color of the background determines the survey type
        rect(xleft = -5, 
             xright = 15, 
             ybottom = -2, 
             ytop = 5, 
             col = ifelse(irow %in% 1, 
                          "white",
                          ifelse(irow %in% 4, "grey90", 
                                 "grey50")))
        
        ## District Labels
        if(irow == 1) mtext(side = 3, 
                            line = -0.75,
                            text = districts$district[idistrict],
                            cex = 0.5)
        
        ## Time axis
        axis(side = 1, 
             at = 1:11, 
             labels = NA,
             tck = -0.05)
        
        if (irow == 2) {
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
        plot_this <- rb_district[,
                                 ispp,
                                 "boat_2" ,
                                 idistrict,
                                 ]
        
        plot_percentiles(values = plot_this,
                         xs = 1:11, 
                         pt.cex = 0.25,
                         pt.colors = "black")
        
        ## 
        if(idistrict == 1) {
          axis(side = 2, 
               las = 1, 
               at = log10(c(0.1, 1, 10, 100, 1000)),
               labels = NA,
               cex.axis = 0.7,
               tck = -0.1)
          
          axis(side = 2, 
               las = 1, 
               at = log10(c(0.1, 1, 10, 100, 1000)),
               labels = c(0.1, 1, 10, 100, 1000),
               cex.axis = 0.75,
               lwd = 0,
               line = -0.25)
          
          axis(side = 2, 
               las = 1, 
               at = log10(c(0.05, 0.5, 5, 50, 500)),
               labels = NA,
               cex.axis = 0.7,
               tck = -0.05)
        }
        
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
  
  par(mar = c(0.25, 0, 0.25, 1))
  plot_legend()
  
  mtext(side = 2, 
        outer = T, 
        text = "Bias Ratio (Simulated Value / True Value)", 
        line = 3, 
        font = 2)
  
  dev.off()
}
