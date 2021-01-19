###############################################################################
## Project:       True and RRMSE of CV, Survey Comparison
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Compare True CV, and RRMSE of CV across species
##                for current and optimized STRS surveys
###############################################################################
rm(list = ls())

##################################################
####  Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[1]

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

scen <- data.frame(survey_type = c("cur", rep("opt", 6) ),
                   strata = c("cur", 3, 5, 10, 10, 15, 20),
                   domain = c("full_domain", 
                              rep(c("district", "full_domain"), each = 3)))

for (irow in 1:7) {
 scen_name <- paste0("SUR_", scen$survey_type[irow], "_", 
                     scen$domain[irow], "_STR_", scen$strata[irow], "_")
 file_name <- paste0(github_dir, "results/", 
                     scen_name, "simulation_results.RData")
 
 load(file_name)
 
}

for (imetric in c("true_cv", "rrmse_cv")) {
 
 png(filename = paste0(output_dir, imetric, ".png"), 
     width = 190, 
     height = 200,
     units = "mm", 
     res = 500)
 
 layout(mat = matrix(c(1:ns_all, ns_all+1, ns_all+1), ncol = 4, byrow = TRUE))
 par(mar = c(0.5, 3, 2.5, 1),
     oma = c(0, 2.5, 0, 0))
 
 for (ispp in c(spp_idx_opt, spp_idx_eval)) {
  merged_metric <- 
   lapply(X = with(scen, paste0("SUR_", survey_type, "_", 
                                domain, "_STR_", strata, "_", imetric )), 
          FUN = function(x) get(x)["obsCV=0",
                                   1:n_years,
                                   ispp,
                                   paste0("boat_2")])
  ylim_ <- max(unlist(merged_metric))
  
  boxplot(merged_metric,
          ylim = c(0, 1.25 * ylim_),
          las = 1,
          axes = F,
          pch = 16,
          col  = c("white", "cyan", "cornflowerblue", "blue4",
                   "darkolivegreen1", "chartreuse1", "darkgreen"),
          cex = 0.5,
          at = c(1, 3:5, 7:9))
  
  abline(v = c(2, 6),
         lty = "dotted", 
         col = "black")
  
  box()
  axis(side = 2, 
       las = 1)
  mtext(side = 3, text = common_names_all[ispp], 
        col = ifelse(ispp %in% spp_idx_opt, "black", "darkgrey"),
        cex = 0.8)
  
 }
 
 par(mar = c(0,0,0,0))
 plot(1,
      xlim = c(0,1),
      ylim = c(0,1),
      type = "n",
      axes = F,
      ann = F)
 
 legend(x = -0.05,
        y = 1,
        legend = c("Current Survey",
                   "District-Level Optimized Survey, 3 Strata per District",
                   "District-Level Optimized Survey, 5 Strata per District",
                   "District-Level Optimized Survey, 10 Strata per District",
                   "Gulf-Wide Optimized Survey, 10 Total Strata",
                   "Gulf-Wide Optimized Survey, 15 Total Strata",
                   "Gulf-Wide Optimized Survey, 20 Total Strata"),
        fill = c("white", 
                 "cyan", "cornflowerblue", "blue4",
                 "darkolivegreen1", "chartreuse1", "darkgreen"),
        xpd = NA,
        cex = 1.25,
        bty = "n")
 
 mtext(side = 2, 
       outer = T, 
       line = 0.5,
       font = 2,
       text = switch(imetric, "rrmse_cv" = "RRMSE of CV",
                                   "true_cv" = "True CV") )
 
 dev.off()
}

