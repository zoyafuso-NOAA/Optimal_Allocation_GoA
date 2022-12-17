##############################################################################
# Project:       
# Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
# Description:   Compare Relative Bias across species
#                for current and optimized STRS surveys
##############################################################################
rm(list = ls())

##################################################
####  Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2)[2]

output_dir <- paste0(c("/Users/zackoyafuso/",
                       "C:/Users/Zack Oyafuso/")[which_machine],
                     "Google Drive/MS_Optimizations/TechMemo/")

##################################
## Import plotting percentile time series function
##################################
source("modified_functions/plot_percentiles.R")

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load("data/processed/optimization_data.RData")
load("data/processed/grid_goa.RData")

##################################
## Load Survey Simulations for the current optimization (scenario L),
## full-scale optimization with 15 strata (scenario A) and 
## area-level optimization with 5 strata per area (scenario B)
##################################
main_rb <- array(dim = c(3, n_years, ns_all, n_boats, n_iters),
                 dimnames = list(c("current", "full_domain", "district"),
                                 paste0("year_", 1:n_years),
                                 common_names_all, 
                                 paste0("boat_", 1:n_boats), 
                                 NULL) )

load(paste0("results/scenario_L/Multispecies_Optimization/Str_current/",
            "simulation_results.RData"))
main_rb["current", , , , ] <- STRS_rel_bias_est

load(paste0("results/scenario_A/Multispecies_Optimization/Str_15/",
            "simulation_results.RData"))
main_rb["full_domain", , , , ] <- STRS_rel_bias_est

load(paste0("results/scenario_K/Multispecies_Optimization/Str_5/",
            "simulation_results.RData"))
main_rb["district", , , , ] <- STRS_rel_bias_est

figs <- data.frame(boat = 1:3,
                   filename = c("appendix/Appendix C plots/Appendix C5",
                                "figures/Figure06",
                                "appendix/Appendix C plots/Appendix C6"))


##################################
## General layout of plots
##################################
gen_layout <- matrix(c(5, 1,2,3, 4,4,4,4), ncol = 2)

for (irow in 1:nrow(figs)) {
  for (spp_group in 1:2) { ## loop over species groupings (opt or eval) -- start
    
    ## Open Device
    png(filename = paste0(output_dir, figs$filename[irow], "_RB_full_domain_",
                          "boat_", figs$boat[irow],
                          c("_opt", "_eval")[spp_group], "_spp.png"),
        width = 190, height = c(220, 170)[spp_group], units = "mm",
        res = 500)
    
    ## Plot Layout
    par(mar = c(.5, 0, 0.25, 0), oma = c(1, 5, 0, 0))
    plot_layout <- rbind(
      cbind(gen_layout + 5 * 0, gen_layout + 5 * 1, 
            gen_layout + 5 * 2, gen_layout + 5 * 3),
      cbind(gen_layout + 5 * 4, gen_layout + 5 * 5, 
            gen_layout + 5 * 6, gen_layout + 5 * 7), 
      cbind(gen_layout + 5 * 8, gen_layout + 5 * 9, 
            gen_layout + 5 * 10, gen_layout + 5 * 11),
      switch(paste0(spp_group),
             "1" = cbind(gen_layout + 5 * 12, gen_layout + 5 * 13, 
                         gen_layout + 5 * 14, gen_layout + 5 * 15),
             "2" = NULL)
    )
    layout(mat =  plot_layout, 
           widths = c(1,0.4, 1,0.4, 1,0.4, 1,0.1),
           heights = c(0.5, 1, 1, 1))
    
    ## Loop over species
    for (ispp in list(spp_idx_opt, spp_idx_eval)[[spp_group]] ) {
      
      ## Calculate y-max for the plot: for each survey type (MARGIN = 1),
      ## calculate the 95% percentile for each year and calculate the max
      ## absolute value of that output
      y_max <- max(abs(apply(main_rb[, , ispp, paste0("boat_", figs$boat[irow]), ],
                             MARGIN = 1, 
                             FUN = function(x) 
                               plot_percentiles(values = x, plot = F))))
      
      y_max <- max(y_max, 25)
      
      for (itype in dimnames(main_rb)[[1]] ) {
        
        ## Base plot
        plot(x = 1, y = 1, type = "n", axes = F, ann = F,
             xlim = c(0, 12), ylim = c(-y_max, y_max))
        
        ## Color of the background corresponds to the  the survey type
        rect(xleft = -5, xright = 15, ybottom = -500, ytop = 500, 
             col = switch(itype,
                          "current" = "white",
                          "full_domain" = "grey90",
                          "district" = "grey50") )
        axis(side = 2,
             las = 1,
             cex.axis = 0.75,
             at = pretty(x = c(-y_max, y_max), n = 3) )
        box()
        
        ## Time axis
        axis(side = 1,
             at = 1:11,
             labels = NA,
             tck = -0.05)
        
        if(itype == "district") {
          axis(side = 1, 
               at = c(1, 11), 
               labels = c("Yr 1", "Yr 11"),
               lwd = F, 
               tick = F, 
               line = -0.5, 
               cex.axis = 1)
        }
        
        ## Plot bias distributions
        plot_percentiles(values = main_rb[itype, , ispp, paste0("boat_", figs$boat[irow]), ],
                         xs = 1:11, 
                         pt.cex = 0.5,
                         pt.colors = "black")
        abline(h = 0, lwd = 0.5, lty = "dotted")
      }
      
      plot(1,type = "n", axes = F, ann = F)  
      plot(1,type = "n", axes = F, ann = F, xlim = c(0, 1), ylim = c(0, 1))
      text(x = 0.5, y = 0.3, 
           labels = common_names_all[ispp], 
           cex = 1.25, 
           font = 2, 
           xpd = NA)
    } ## Loop over three survey scenarios -- end
    
    ## Plot Legend
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
          text = "Percent Bias (100% (Sim - True) / True)", 
          line = 3, 
          font = 2)
    
    ## Close Device
    dev.off()
  } ## loop over species groupings (opt or eval) -- end
}
