###############################################################################
## Project:      PowerPoint Plots for SAFS Talk 
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:  
###############################################################################
rm(list = ls()) 

##################################################
####   Load packages 
##################################################
library(sp)
library(raster)
library(RColorBrewer)

##################################################
####   Set up Directories
##################################################
github_dir <- "C:/Users/Zack Oyafuso/Documents/GitHub/Optimal_Allocation_GoA/results/full_domain/Single_Species_Optimization/Pacific cod/"
output_dir <- "C:/Users/Zack Oyafuso/Documents/GitHub/Optimal_Allocation_GoA/videos/Survey_Optimization_SS/"

##################################################
####   Base Plots
##################################################
base_plot <- function(abline_h = NULL, abline_col = NULL,
                      plot_sol = NULL, past_sol = NULL) {
  par(mar = c(5,5,1,3))
  plot(x = 1, y = 1, type = "n",
       ylim = c(0, 900), xlim = c(0, .15),
       las = 1,
       xlab = "Expected CV", ylab = "Total sample size")
  
}


n_runs <- length(dir(github_dir))
result_pts <- data.frame()
for (irun in 1:n_runs) {
  load(paste0(github_dir, "Run", irun, "/result_list.RData"))
  
  result_pts <- rbind(result_pts,
                      c(result_list$CV_constraints, result_list$n))
  
}

png(filename = paste0(output_dir, "plot0.png"), 
    width = 5.5, height = 5.5, units = "in", res = 500)
base_plot()
dev.off()

png(filename = paste0(output_dir, "plot1.png"), 
    width = 5.5, height = 5.5, units = "in", res = 500)
base_plot()
points(result_pts[1, ], pch = 16)
dev.off()

png(filename = paste0(output_dir, "plot2.png"), 
    width = 5.5, height = 5.5, units = "in", res = 500)
base_plot()
points(result_pts[1, ], pch = 16, col = 'grey')
abline(v = result_pts[2, 1], lty = "dotted", col = "grey")
dev.off()

png(filename = paste0(output_dir, "plot3.png"), 
    width = 5.5, height = 5.5, units = "in", res = 500)
base_plot()
points(result_pts[1, ], pch = 16, col = 'grey')
abline(v = result_pts[2, 1], lty = "dotted", col = "grey")
points(result_pts[2, ], pch = 16)
dev.off()

png(filename = paste0(output_dir, "plot4.png"), 
    width = 5.5, height = 5.5, units = "in", res = 500)
base_plot()
points(result_pts, pch = 16)
lines(result_pts, pch = 16)
dev.off()

png(filename = paste0(output_dir, "plot5.png"), 
    width = 5.5, height = 5.5, units = "in", res = 500)
base_plot()
points(result_pts, pch = 16)
lines(result_pts, pch = 16)
abline(h = c(292, 550, 825), lty = "dotted", col = "grey")
text(x = 0.13, y = c(292, 550, 825), labels = paste(1:3, "Boats"), pos = 3)

sol_idx <- sapply(X = c(292, 550, 825),
                  FUN = function(x) which.min(abs(x - result_pts[, 2])))
points(result_pts[sol_idx, ], col = "blue", pch = 16, cex = 2)
segments(x0 = result_pts[sol_idx, 1],
         y0 = -100,
         y1 = result_pts[sol_idx, 2],
         lty = "dotted", col = "grey")
dev.off()
