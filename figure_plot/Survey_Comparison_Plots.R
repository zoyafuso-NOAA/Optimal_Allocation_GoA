###############################################################################
## Project:       Survey Comparison Plots
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Bias in Estimate, True CV, and RRMSE of CV
###############################################################################
rm(list = ls())

#######################################
## Set up directories
#######################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[2]
VAST_model <- "11" 

github_dir = paste0(c("/Users/zackoyafuso/Documents/", 
                      "C:/Users/Zack Oyafuso/Documents/",
                      "C:/Users/zack.oyafuso/Work/")[which_machine], 
                    "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/")

figure_dir = paste0(c("/Users/zackoyafuso/Google Drive/", 
                      "C:/Users/Zack Oyafuso/Google Drive/")[which_machine],
                    "MS_Optimizations/TechMemo/figures/")

##################################################
####   Load Results
##################################################
load( paste0(github_dir, "optimization_data.RData") )
load( paste0(github_dir, "Survey_Comparison_Simulations/",
             "Survey_Simulation_Results.RData") )
load( paste0(github_dir, "Survey_Comparison_Simulations/",
             "Simple_RS_Simulation_Results.RData") )
load( paste0(github_dir, "Spatiotemporal_Optimization/",
             "STRS_Sim_Res_spatiotemporal.RData") )

##################################################
####   Constants
##################################################
istrata <- 2 #10 strata index

common_names <- c("arrowtooth flounder", "walleye pollock", "Pacific cod",
                  "rex sole", "flathead sole", "Pacific halibut",
                  "southern rock sole", "northern rock sole", 
                  "Dover sole", "Pacific ocean perch",
                  "blackspotted/rougheye\nrockfishes",
                  "silvergray rockfish", "northern rockfish", "dusky rockfish",
                  "shortspine thornyhead")

spp_order = c(1,  3,  5,
              6,  7,  8,
              13, 14, 2,
              4,  9,  10,
              12, 11, 15)
# spp_order <- c(8,  3,  5, 7, 6,
#                1, 2, 14, 4, 13, 
#                9, 10, 11, 12, 15)

# spp_order <- order(
#   apply(X = Survey_rel_bias_est[,,1], 
#         MARGIN = 2, 
#         FUN = median),
#   decreasing = T)

##################################################
####   Bias in Estimate
##################################################
{
  # png(filename = paste0(figure_dir, "Bias_Est.png"),
  #     units = "in", 
  #     width = 11, 
  #     height = 5, 
  #     res = 500)
  
  par(mfrow = c(3, 5), 
      mar = c(0.5, 4, 0.5, 0), 
      oma = c(2, 1, 2, 0.5))
  
  # png(filename = paste0(figure_dir, 'Bias_Est.png'),
  # units = 'mm', width = 200, height = 150, res = 500)
  
  for (ispp in spp_order) {
    
    ymax <- 1.1 * max(abs(Survey_rel_bias_est[,ispp,]), 
                      na.rm  = T)
    ymax <- max(ymax, 5)
    
    plot(1, 
         type = "n", 
         xlim = c(0, 12),
         ylim = c(-ymax,ymax),  
         las = 1, 
         axes = F, 
         ann = F)
    
    
    if(ispp == spp_order[3]) 
      legend(x = -3, 
             y = 9, 
             col = c("red", "blue", "black"), 
             text.col = c("red", "blue", "black"),
             legend = paste(1:nboats, "Boat"),
             xpd = NA, 
             pch = 0,  
             horiz = T, 
             bty = "n", 
             cex = 1.5, 
             lty = 1, 
             x.intersp = 0.25)
    
    abline(h = 0, 
           lty = "dotted")
    box()
    axis(side = 2, 
         las = 1)
    
    legend("topright", 
           legend = common_names[ispp], 
           bty = "n")
    
    
    #Survey Type X-axis label for bottom row
    if (ispp %in% spp_order[ns:(ns-4)])  
      axis(side = 1, 
           at = c(2, 6, 10), 
           labels = c("Current", "Optimized", "SRS"))
    
    
    
    #Plot Relative Biases
    boxplot(Survey_rel_bias_est[, ispp, ], 
            add = T, 
            at = 1:3, 
            axes = F, 
            border = c("red", "blue", "black"), 
            pch = 16)
    
    boxplot(STRS_rel_bias_est[,ispp , , istrata], 
            add = T, 
            at = 5:7, 
            axes = F, 
            border = c("red", "blue", "black"), 
            pch = 16)
    
    boxplot(SRS_rel_bias_est[,ispp , ], 
            add = T, 
            at = 9:11, 
            axes = F, 
            border = c("red", "blue", "black"), 
            pch = 16)
  }
  mtext(side = 2, 
        text = "Relative Percent Bias", 
        outer = T, 
        line = -0.5)
  # dev.off()
}

####################################
## Bias in CV
####################################

# {
#   # png(filename = paste0(figure_dir, "Bias_CV.png"),
#   #     units = "mm", width = 200, height = 150, res = 500)
#   par(mfrow = c(5, 3), mar = c(0.5, 4, 0.5, 0), oma = c(2, 1, 2, 0.5))
#   
#   for (ispp in c(1,  3,  5,
#                  6,  7,  8,
#                  13, 14, 2,
#                  4,  9,  10,
#                  12, 11, 15)){
#   
#     ymax = max(abs(Survey_rel_bias_cv[,ispp,]), na.rm  = T)
#     ymax = max(ymax, 5)
#     
#     plot(1, 
#          type = "n", 
#          ylim = c(-ymax,ymax), xlim = c(0,8), 
#          las = 1,
#          axes = F, 
#          ann = F)
#     if(ispp == 3) legend(x = -2, y = 10, 
#                          legend = paste(1:3, "Boat"),
#                          col = c("red","blue","black"), 
#                          text.col = c("red","blue","black"),
#                          xpd = NA, 
#                          pch = 0, 
#                          horiz = T, 
#                          bty = "n", 
#                          cex = 1.5, 
#                          lty = 1,
#                          x.intersp = 0.25)
#     abline(h=0, lty = "dotted")
#     box()
#     axis(side = 2, las = 1)
#     legend("top", sci_names[ispp], bty = "n")
#     if(ispp %in% c(12,11,15))
#       axis(side = 1, at = c(2,6),
#            labels = c("Current", "Optimized"),
#            cex.axis = 1)
#     
#     boxplot(Survey_rel_bias_cv[,ispp,], add = T, at = 1:3, axes = F,
#             border = c("red", "blue", "black"), pch = 16)
#     boxplot(STRS_rel_bias_cv[,ispp,,istrata], add = T, at = 5:7, axes = F,
#             border = c("red", "blue", "black"), pch = 16)
#   }
#   mtext(side = 2, "Relative Percent Bias", outer = T, line = -.5)
#   # dev.off()
# }

####################################
## True CV
####################################
{
  # png(filename = paste0(figure_dir, "True_CV.png"),
  #     units = "in", 
  #     width = 11, 
  #     height = 5, 
  #     res = 500)
  
  par(mfrow = c(3, 5), 
      mar = c(0.5, 4, 0.5, 0), 
      oma = c(2, 1, 2, 0.5))
  
  for (ispp in spp_order) {
    
    ymax = 1.1 * max( c(Survey_true_cv_array[,ispp,], 
                        STRS_true_cv_array[,ispp,,]), 
                      na.rm = T) 
    
    #Base Plot
    plot(1, 
         type = "n", 
         xlim = c(0, 12),
         ylim = c(0, ymax),  
         las = 1,
         axes = F, 
         ann = F)
    
    if (ispp == spp_order[3]) {
      legend(x = -3, 
             y = 0.115, 
             col = c("red", "blue", "black"),
             legend = paste(1:3, "Boat"), 
             pch = 0, 
             horiz = T,
             bty = "n", 
             cex = 1.5, 
             lty = 1, 
             x.intersp = 0.25,
             text.col = c("red", "blue", "black"), 
             xpd = NA)
    }
    
    axis(side = 2, 
         las = 1)
    
    if (ispp %in% spp_order[ns:(ns-4)]) 
      axis(side = 1, 
           at = c(2, 6, 10), 
           labels = c("Current", "Optimized", "SRS"),
           cex.axis =1)
    
    legend("topright", 
           legend = common_names[ispp], 
           bty = "n", 
           text.font = 3)
    box()
    axis(side = 2, 
         las = 1)
    
    boxplot(Survey_true_cv_array[, ispp, ], 
            add = T, 
            at = 1:3,
            axes = F, 
            border = c("red", "blue", "black"), 
            pch = 16)
    
    boxplot(STRS_true_cv_array[, ispp, , istrata ], 
            add = T, 
            at = 5:7, 
            axes = F, 
            border = c("red", "blue", "black"), 
            pch = 16)
    
    boxplot(SRS_true_cv_array[, ispp, ], 
            add = T, 
            at = 9:11, 
            axes = F, 
            border = c("red", "blue", "black"), 
            pch = 16)
    
  }
  mtext(side = 2, 
        text = "True CV", 
        outer = T, 
        line = -0.5)
  # dev.off()
}

####################################
## RRMSE of CV
####################################
{
  # png(filename = paste0(figure_dir, "RRMSE_CV.png"),
  #     units = "in", 
  #     width = 11, 
  #     height = 5, 
  #     res = 500)
  
  par(mfrow = c(3, 5), 
      mar = c(0.5, 4, 0.5, 0), 
      oma = c(2, 1, 2, 0.5))
  
  for (ispp in spp_order) {
    
    ymax = 1.1 * max( c(Survey_rrmse_cv_array[,ispp,], 
                        STRS_rrmse_cv_array[,ispp,,]), 
                      na.rm = T) 
    
    #Base Plot
    plot(1, 
         type = "n", 
         xlim = c(0, 12),
         ylim = c(0, ymax),  
         las = 1,
         axes = F, 
         ann = F)
    
    if (ispp == spp_order[3]) {
      legend(x = -3, 
             y = 0.505, 
             col = c("red", "blue", "black"),
             legend = paste(1:3, "Boat"), 
             pch = 0, 
             horiz = T,
             bty = "n", 
             cex = 1.5, 
             lty = 1, 
             x.intersp = 0.25,
             text.col = c("red", "blue", "black"), 
             xpd = NA)
    }
    
    axis(side = 2, 
         las = 1)
    
    if (ispp %in% spp_order[ns:(ns-4)]) 
      axis(side = 1, 
           at = c(2, 6, 10), 
           labels = c("Current", "Optimized", "SRS"),
           cex.axis =1)
    
    legend("topright", 
           legend = common_names[ispp], 
           bty = "n", 
           text.font = 3)
    box()
    axis(side = 2, 
         las = 1)
    
    boxplot(Survey_rrmse_cv_array[, ispp, ], 
            add = T, 
            at = 1:3,
            axes = F, 
            border = c("red", "blue", "black"), 
            pch = 16)
    
    boxplot(STRS_rrmse_cv_array[, ispp, , istrata ], 
            add = T, 
            at = 5:7, 
            axes = F, 
            border = c("red", "blue", "black"), 
            pch = 16)
    
    boxplot(SRS_rrmse_cv_array[, ispp,  ], 
            add = T, 
            at = 9:11, 
            axes = F, 
            border = c("red", "blue", "black"), 
            pch = 16)
    
  }
  mtext(side = 2, 
        text = "RRMSE of CV", 
        outer = T, 
        line = -0.5)
  # dev.off()
}

load( paste0(github_dir, "Spatiotemporal_Optimization/",
             "optimization_knitted_results.RData") )
apply(STRS_true_cv_array[,,,2], MARGIN = 2:3, median)
settings[8,]
