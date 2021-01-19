###############################################################################
## Project:       Sensitivity: Observation Error
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Compare True CV, and RRMSE of CV across species
##                for surveys with varying added observation error for 
##                two-boat optimized surveys
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
####   Load Data
##################################################
load(paste0(github_dir, "data/optimization_data.RData"))
load(paste0(github_dir, "results/MS_optimization_knitted_results.RData"))
source(paste0(github_dir, "modified_functions/plot_percentiles.R"))

scen <- data.frame(survey_type = c("cur", "opt", "opt"),
                   strata = c("cur", 3, 10),
                   domain = c("full_domain", "district", "full_domain"))

for (irow in 1:nrow(scen)) {
 scen_name <- paste0("SUR_", scen$survey_type[irow], "_",
                     scen$domain[irow], "_STR_", scen$strata[irow], "_")
 file_name <- paste0(github_dir, "results/",
                     scen_name, "simulation_results.RData")
 
 load(file_name)
}

true_cv <- lapply(X = 1:nrow(scen),
                  FUN = function(x) 
                   with(scen[x, ], 
                        get(paste0("SUR_", survey_type, 
                                   "_", domain, 
                                   "_STR_", strata, 
                                   "_true_cv")))[,,,"boat_2"])

rrmse_cv <- lapply(X = 1:nrow(scen),
                   FUN = function(x) 
                    with(scen[x, ], 
                         get(paste0("SUR_", survey_type, 
                                    "_", domain, 
                                    "_STR_", strata, 
                                    "_rrmse_cv")))[,,,"boat_2"])

for (imetric in c("rrmse_cv", "true_cv")) {
 
 # Set up png file
 png(filename = paste0(result_dir,
                       "/Sensitivity_ObsError_", imetric,
                       ".png"),
     width = 190,
     height = 200,
     units = "mm",
     res = 500)
 
 par(mfrow = c(6, 4), mar = c(2, 3, 2, 1), oma = c(0, 2, 0, 0))
 for (ispp in c(spp_idx_opt, spp_idx_eval)) {
  plot_this <- lapply(X = get(imetric), FUN = function(x) x[, , ispp])
  plot_this <- lapply(X = obs_cv, 
                      FUN = function(x) 
                       lapply(X = plot_this, 
                              FUN = function(y) y[paste0("obsCV=", x),]))
  plot_this <- lapply(plot_this, 
                      FUN = function(x) do.call(rbind, x))
  plot_this <- t(do.call(rbind, plot_this))
  
  boxplot( plot_this,
           at = c(1:3, 5:7, 9:11, 13:15),
           axes = F,
           ann = F,
           ylim = c(0, max(plot_this)),
           col = c("white", "green", "blue"),
           pch = 16)
  
  mtext(side = 3, text = common_names_all[ispp], cex = 0.75)
  box()
  axis(side = 2, las = 1)
  axis(side = 1, at = c(2, 6, 10, 14), 
       labels = paste0(obs_cv*100, "%"),
       cex.axis = 0.9,
       line = -0.5,
       tick = F)
  axis(side = 1, at = c(2, 6, 10, 14), 
       labels = NA, tck = -0.05)
 }
 
 plot(1, 
      type = "n",
      axes = F, 
      ann = F, 
      xlim = c(0, 1), 
      ylim = c(0, 1))
 legend(x = 0, y = 1, 
        legend = c("Current Survey", 
                   "District-Level Optimization", 
                   "Gulf-Wide Optimization"),
        fill = c("white", "green", "blue"),
        bty = "n",
        cex = 1.5,
        xpd = NA)
 
 mtext(side = 2, 
       outer = T, 
       text = switch(imetric, 
                     "rrmse_cv" = "RRMSE of CV",
                     "true_cv" = "True CV"), 
       line = 0.5, 
       font = 2)
 
 dev.off()
}



