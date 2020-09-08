###############################################################################
## Project:      Sample Size Versus Precision
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:  Figure 4: Total sample size required to meet each value of 
##               upper coefficient of variation (CV) constraint, accounting 
##               only for spatial variability (top) or both spatial and 
##               temporal variability (bottom). This example uses a five-strata 
##               optimization, but qualitative results were consistent with 
##               solutions given additional strata (Supplementary SX). Both 
##               optimizations were conducted under the One-CV constraint 
##               approach where all species have the same upper CV constraint 
##               in the optimization. Horizontal dotted grey lines indicate the
##               sampling levels for one, two, and three boats.
##
##               SFigure 
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
which_machine = c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[2]

VAST_model <- "6g"
github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/")

figure_dir <- paste0(c("/Users/zackoyafuso/", 
                       "C:/Users/Zack Oyafuso/")[which_machine],
                     "Google Drive/MS_Optimizations/figure_plot/")

##################################################
####   Plot Settings
##################################################
samples <- c(280, 550, 820)
plot_settings <- 
  data.frame(type = c("Spatial", "Spatiotemporal"),
             xmin = c(0.09, 0.15),
             xmax = c(0.24, 0.30),
             xlabel = c(0.19, 0.2775))

##################################################
####   Figure 4
##################################################

{ 
  png(filename = paste0(figure_dir, "Fig4_N_CV_Tradeoff.png"),
      width = 90, 
      height = 150, 
      res = 500, 
      units = "mm")
  
  par(mfrow = c(2, 1), 
      mar = c(2, 0, 1, 0.5), 
      oma = c(2, 4, 0, 0))
  
  for (irow in 1:2) {
    
    ##########################
    ## Load Data
    ##########################
    load(paste0(github_dir, plot_settings$type[irow], 
                "_Optimization/optimization_knitted_results.RData"))
    
    #############################
    ## Main Plot
    #############################
    subsettings = subset(settings, 
                         strata == 5)
    subsettings = subsettings[order(subsettings$cv),]
    
    plot(n ~ cv, 
         data = subsettings, 
         las = 1, 
         pch = 16,
         ylim = c(0, 1100),
         ann = F, 
         cex.axis = 0.85)
    
    lines(n ~ cv, 
          data = subsettings)
    
    #1, 2, and 3-boat sample size lines
    abline(h = samples, 
           col = "grey", 
           lty = "dashed")
    
    #Boat labels
    text(x = plot_settings$xlabel[irow], 
         y = samples, 
         labels = paste(1:3, "Boat"), 
         pos = 1)
    
    #Plot subtitle
    legend("top", 
           legend = paste(plot_settings$type[irow], "Optimization"), 
           bty = "n")
    
  }
  
  #Plot Axes Names
  mtext(side = 1, 
        text = "Upper CV Constraint", 
        outer = T, 
        line = 0.5)
  
  mtext(side = 2, 
        text = "Total Optimized Sample Size", 
        outer = T, 
        line = 3)
  
  dev.off()
  }


#####################################################
## Supplemental Figure 4
#####################################################
plot_settings <- 
  data.frame(type = c("Spatial", "Spatiotemporal"),
             xmin = c(0.09, 0.15),
             xmax = c(0.24, 0.30),
             xlabel = c(0.22, 0.2775))
{
  png(paste0(figure_dir, "Supplemental_Figures/SFig4_N_CV_Tradeoff.png"),
      width = 140, 
      height = 180, 
      res = 500, 
      units = "mm")
  
  par(mfcol = c(6, 2), 
      mar = c(0, 0, 0, 0), 
      oma = c(4, 6, 3, 0.5))
  
  for(irow in 1:2){
    
    ##########################
    ## Load Data
    ##########################
    load(paste0(github_dir, plot_settings$type[irow], 
                "_Optimization/optimization_knitted_results.RData"))
    
    for (istrata in c(5, 10, 15, 20, 30, 60)) {
      
      #############################
      ## Plot
      #############################
      if (irow == 1) {
        sub_settings = subset(settings, 
                              strata == istrata)
        sub_settings = sub_settings[order(sub_settings$cv, decreasing = T), ] 
        
        plot(n ~ cv, 
             data = sub_settings, 
             las = 1, 
             pch = 16,
             ylim = c(0, 1100), 
             xlim = c(plot_settings$xmin[irow], 
                      plot_settings$xmax[irow]),
             ann = F, 
             axes = F)
        box()
        axis(side = 2,
             las = 1)
        if (istrata == 60) axis(side = 1)
        
        lines(n ~ cv, 
              data = sub_settings, 
              subset = strata == istrata)
      }
      
      if (irow == 2) {
        plot(n ~ cv, 
             data = settings, 
             subset = strata == istrata, 
             las = 1, 
             pch = 16,
             ylim = c(0, 1100), 
             xlim = c(plot_settings$xmin[irow], 
                      plot_settings$xmax[irow]),
             ann = F,
             axes = F)
        box()
        if (istrata == 60) axis(side = 1)
        
        lines(n ~ cv, 
              data = settings, 
              subset = strata == istrata)
        }
      
      abline(h = samples, 
             col = "grey", 
             lty = "dashed")
      
      text(x = plot_settings$xlabel[irow], 
           y = samples, 
           labels = paste(1:3, "Boat"), 
           pos = 1)
      
      text(x = plot_settings$xlabel[irow],
           y = plot_settings$ymax[irow],
           labels = paste(istrata, "Strata"), 
           font = 2, 
           pos = 1, 
           cex = 1.25)
      
      if (istrata == 5) mtext(side = 3, 
                              text = paste(plot_settings$type[irow], 
                                           "Optimization"),
                              line = 1)
    }
  }
  
  mtext(side = 1, 
        text = "Upper CV Constraint", 
        outer = T, 
        line = 3)
  mtext(side = 2, 
        text = "Total Sample Size", 
        outer = T, 
        line = 4)
  dev.off()
}
