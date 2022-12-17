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
## Import plotting percentile time series function and constants
##################################
source("modified_functions/plot_percentiles.R")
load("data/processed/optimization_data.RData")

##################################################
####   Result Object
##################################################
index <- array(dim = c(2, n_years, ns_all, n_boats, n_iters),
               dimnames = list(c("opt", "current"),
                               paste0("year_", 1:n_years),
                               common_names_all,
                               paste0("boat_", 1:n_boats),
                               NULL))

##################################################
####   Load data
##################################################
load("results/scenario_J/Multispecies_Optimization/Str_5/simulation_results.RData")
index["opt", , , ,] <- STRS_sim_index

load("results/scenario_M/Multispecies_Optimization/Str_current/simulation_results.RData")
index["current", , , ,] <- STRS_sim_index

##################################################
####   Plot
##################################################
png(filename = paste0(output_dir, 
                      "/appendix/Appendix C plots/Appendix C7_Index_Ratio.png"),
    width = 190, height = 220, units = "mm", res = 500)

par(mfrow = c(7, 4), mar = c(0, 2, 2, 1), oma = c(4, 3, 0, 0))
for (ispp in c(common_names_opt, common_names_eval)){
  ratio <- index["current", , ispp, "boat_2", ] / 
    index["opt", , ispp, "boat_2", ]
  
  y_max <- max(max(plot_percentiles(values = ratio, plot = F)) - 1, 0.25)
  
  ## Base plot
  plot(x = 1, y = 1, type = "n", axes = F, ann = F,
       xlim = c(0, 12), 
       ylim = c(max(1 - y_max, 0), 1 + y_max))
  mtext(side = 3, text = ispp)
  
  plot_percentiles(values = ratio, plot = T, 
                   xs = 1:n_years, 
                   pt.colors = "black",
                   inner_color = "chartreuse",
                   outer_color = "darkgreen")
  axis(side = 2, las = 1)
  axis(side = 1, labels = F, tcl = 0.25, at = 1:n_years)
  
  if(ispp %in% common_names_eval[ns_eval:(ns_eval-3)]) axis(side = 1,
                                                            at = 1:11,
                                                            tick = FALSE)
  abline(h = 1)
  box()
  
}
mtext(side = 1, text = "Survey Year", outer = TRUE, line = 2.5, font = 2)
mtext(side = 2, text = "Index Ratio (Proposed Design / Exisiting Design)",
      outer = TRUE, font = 2, line = 1)
dev.off()
