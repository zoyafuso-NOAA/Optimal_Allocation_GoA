###############################################################################
## Project:       Sensitivity: Observation Error
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Compare Relative Bias, True CV, and RRMSE of CV across species
##                for current and optimized STRS surveys
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2)[1]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")
result_dir <- paste0(c("/Users/zackoyafuso/", 
                       "C:/Users/Zack Oyafuso/")[which_machine],
                     "Google Drive/MS_Optimizations/TechMemo/figures")

##################################################
####   Import Libraries
##################################################
library(RColorBrewer)

##################################################
####   Load Data
##################################################
load(paste0(github_dir, "data/optimization_data.RData"))
load(paste0(github_dir, "results/simulation_result.RData"))
source(paste0(github_dir, "modified_functions/plot_percentiles.R"))

##################################################
####   Plot
##################################################
for (spp_set in  c("opt", "eval")) {
  
  ## Set up png file
  png(filename = paste0(result_dir, 
                        "/Comparison_",
                        spp_set,
                        ".png"),
      width = 190,
      height = c("opt" = 220, "eval" = 120)[spp_set],
      units = "mm",
      res = 500)
  
  ## Set up plot layout based on which species are being plotted
  layout(mat = matrix(data = list(
    "opt" = c(1, 2, 3,       25, 26, 27,  
              4, 5, 6,       28, 29, 30,
              7, 8, 9,       31, 32, 33,
              10, 11, 12,    34, 35, 36,
              13, 14, 15,    37, 38, 39,
              16, 17, 18,    40, 41, 42,
              19, 20, 21,    43, 44, 45,
              22, 23, 24,    46, 46, 46),
    "eval" = c(1, 2, 3,       13, 14, 15,  
               4, 5, 6,       16, 17, 18,   
               7, 8, 9,       19, 20, 21,    
               10, 11, 12,    22, 22, 22))[[spp_set]],
    ncol = 6,
    byrow = T),
    widths = rep(x = c(1,0.5, 0.5), times = 3) )
  
  par(mar = c(0.25, 1.75, 0.25, 1.75), 
      oma = c(0, 1, 3, 1))
  
  ## Some constants based on which species to plot
  ns = get(paste0("ns_", spp_set))
  spp_idx <- get(paste0("spp_idx_", spp_set))
  common_names <- get(paste0("common_names_", spp_set))
  isample = 2 # 2-boat scenario 
  
  for (ispp in 1:ns) {
    
    ## Relative Bias
    par(mar = c(0.25, 1.75, 0.25, 1.75))
    
    ylim_ <- max(abs(apply(
      X = cbind(STRS_rel_bias_est["obsCV=0", 
                                  1:NTime, 
                                  spp_idx[ispp], 
                                  isample,
                                  1:Niters], 
                Current_rel_bias_est["obsCV=0", 
                                     1:NTime, 
                                     spp_idx[ispp], 
                                     isample,
                                     1:Niters]),
      MARGIN = 1,
      FUN = quantile,
      probs = c(0.025, 0.975),
      na.rm = T)))
  
    plot(1,  
         type = "n",
         ylim = c(-ylim_, ylim_),
         xlim = c(1, 23),
         axes = F,
         ann = F)
    box()
    
    plot_percentiles(
      values = STRS_rel_bias_est["obsCV=0",
                                 1:NTime,
                                 spp_idx[ispp],
                                 isample,
                                 1:Niters],
      xs = 1:11, 
      inner_color = "cadetblue1",
      outer_color = "dodgerblue")
    
    plot_percentiles(
      values = Current_rel_bias_est["obsCV=0",
                                    1:NTime,
                                    spp_idx[ispp],
                                    isample,
                                    1:Niters],
      xs = 13:23, 
      inner_color = "tomato2",
      outer_color =  "firebrick" )
    
    if (ispp %in% c(1, c("opt" = 9, "eval" = 5)[spp_set])) {
      mtext(side = 3, 
            text = "Percent\nBias", 
            font = 2)
    }
    abline(h = 0)
    axis(side = 2, las = 1)
    
    ## True CV
    par(mar = c(2, 1.75, 0.25, 1.5))
    ylim_ <- max(abs(
      cbind(STRS_true_cv_array["obsCV=0",
                               1:NTime,
                               spp_idx[ispp], 
                               isample],
            Current_true_cv_array["obsCV=0",
                                  1:NTime,
                                  spp_idx[ispp],
                                  isample])))
    
    boxplot(cbind(STRS_true_cv_array["obsCV=0",
                                     1:NTime,
                                     spp_idx[ispp],
                                     isample],
                  Current_true_cv_array["obsCV=0",
                                        1:NTime,
                                        spp_idx[ispp],
                                        isample]),
            ylim = c(0, 1.25 * ylim_),
            las = 1,
            axes = F,
            pch = 16,
            col  = c("dodgerblue", "firebrick"))
    box()
    axis(side = 2, 
         las = 1)
    if (ispp %in% c(1, c("opt" = 9, "eval" = 5)[spp_set])) {
      mtext(side = 3, 
            text = "True\nCV", 
            font = 2)
    }
    
    ## RRMSE of CV
    ylim_ <- max(abs(cbind(
      STRS_rrmse_cv_array["obsCV=0", 
                          1:NTime,
                          spp_idx[ispp],
                          isample],
      Current_rrmse_cv_array["obsCV=0", 
                             1:NTime,
                             spp_idx[ispp], 
                             isample])))
    
    boxplot(cbind(
      STRS_rrmse_cv_array["obsCV=0", 
                          1:NTime, 
                          spp_idx[ispp], 
                          isample],
      Current_rrmse_cv_array["obsCV=0",
                             1:NTime,
                             spp_idx[ispp],
                             isample]),
      ylim = c(0, 1.25 * ylim_),
      las = 1,
      names = NA,
      axes = F,
      pch = 16,
      col = c("dodgerblue", "firebrick"))
    box()
    axis(side = 2, 
         las = 1)
    if (ispp %in% c(1, c("opt" = 9, "eval" = 5)[spp_set])) {
      mtext(side = 3, 
            text = "RRMSE\nCV", 
            font = 2)
    }
    
    #Species Label
    text(x = -1, 
         y = par("usr")[3] - diff(par("usr")[3:4]) * 0.15,
         labels = common_names[ispp],
         xpd = NA,
         cex = 1.25)
  }
  
  ##Legend
  par(mar = c(0, 1, 0, 2))
  plot(1,
       xlim = c(0, 1),
       ylim = c(0, 1),
       type = "n",
       axes = F,
       ann = F)
  
  #Percentile legend
  rect(xleft =   c(0.1, 0.1, 0.4, 0.4), 
       xright =  c(0.2, 0.2, 0.5, 0.5),
       ybottom = c(0.0, 0.3, 0.0, 0.3), 
       ytop =    c(0.8, 0.5, 0.8, 0.5),
       col =     c("dodgerblue", "cadetblue1", "firebrick", "tomato2"),
       border = F)
  points(x = c(0.15, 0.45), 
         y = c(0.4, 0.4), 
         pch = 15)
  segments(x0 = c(0.15, 0.45),
           y0 = 0.0,
           y1 = 0.8,
           lty = "dashed")
  text(x = 0.3,
       y = c(0.025, 0.275, 0.4, 0.525, 0.75),
       labels = paste0(c(5,25,50,75,95), "%"))
  text(x = c(0.15, 0.45),
       y = c(0.925, 0.925),
       labels = c("Opt.\nSurvey", "Curr.\nSurvey") )
  
  #Years Legend
  points(x = seq(from = 0.6, to = 1, length = 11),
         y = rep(x = 0.5, times = 11),
         col = rev(grey.colors(n = 11, start = 0, end = 0.9)),
         pch = 15,
         cex = 2)
  text(x = c(0.6, 1),
       y = 0.65,
       labels = c(1996, 2019) )
  
  dev.off()
}

