###############################################################################
## Project:       Sensitivity: Observation Error
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Compare relative bias of total abundance index across species
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

rb <- lapply(X = 1:nrow(scen),
             FUN = function(x) 
              with(scen[x, ], 
                   get(paste0("SUR_", survey_type, 
                              "_", domain, 
                              "_STR_", strata, 
                              "_rb_agg")))[,,,"boat_2",])

# Set up png file
# png(filename = paste0(result_dir,
#                       "/Sensitivity_ObsError_", imetric,
#                       ".png"),
#     width = 190,
#     height = 200,
#     units = "mm",
#     res = 500)

par(mfrow = c(8, 3), mar = c(1, 3, 1, 1), oma = c(0, 2, 0, 0))

for (ispp in c(spp_idx_opt, spp_idx_eval)) {
 plot_this <- lapply(X = rb, FUN = function(x) x[, , ispp,])
 plot_this <- lapply(X = obs_cv, 
                     FUN = function(x) 
                      lapply(X = plot_this, 
                             FUN = function(y) y[paste0("obsCV=", x), ,]))
 # plot_this <- lapply(plot_this, 
 #                     FUN = function(x) do.call(rbind, x))
 # plot_this <- t(do.call(rbind, plot_this))
 
 opt_ylim <- max(abs(unlist(lapply(X = plot_this,
                                   FUN = function(x) 
                                    lapply(X = x, FUN = quantile,
                                           probs = c(0.025, 0.975),
                                           na.rm = T)))))
 
 plot(1,  
      type = "n",
      ylim = c(-opt_ylim, opt_ylim),
      xlim = c(1, 155),
      axes = F,
      ann = F)
 
 mtext(side = 3, text = common_names_all[ispp], cex = 0.75)
 
 # if (ispp %in% c(1, c("opt" = 9, "eval" = 5)[spp_set])) {
 #   mtext(side = 3, 
 #         text = "Percent\nBias", 
 #         font = 2)
 # }
 
 x_ids <- 1
 for (icv in 1:n_obs_cv) {
  for (isurvey in 1:3) {
   plot_percentiles(values = plot_this[[icv]][[isurvey]],
                    xs = x_ids:(x_ids + 10),
                    inner_color = c("white", "green", "dodgerblue"),
                    outer_color = "dodgerblue",
                     pt.cex = 0.5)
   
   x_ids <- x_ids + 12
  }
  x_ids <- x_ids + 5
 }
 
 abline(h = 0)
 axis(side = 2, las = 1)
 box()
 
}

plot(1, 
     type = "n",
     axes = F, 
     ann = F, 
     xlim = c(0, 1), 
     ylim = c(0, 1))
legend("center", 
       legend = c("Current Survey", 
                  "District-Level Optimization", 
                  "Gulf-Wide Optimization"),
       fill = c("white", "green", "blue"),
       bty = "n",
       cex = 1.5,
       xpd = NA)

mtext(side = 2, 
      outer = T, 
      text = "Relative Percent Bias (100% (Sim - True) / True)", 
      line = 0, 
      font = 2)
