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
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2)[2]
output_dir <- paste0(c("/Users/zackoyafuso/",
                       "C:/Users/Zack Oyafuso/")[which_machine],
                     "Google Drive/MS_Optimizations/TechMemo/")

##################################################
####   Set up appendix output directories
##################################################
if(!dir.exists(paste0(output_dir, "/appendix/Appendix C plots/" )))
  dir.create(paste0(output_dir, "/appendix/Appendix C plots/" ))

##################################
## Import plotting percentile time series function
##################################
source("modified_functions/plot_percentiles.R")
library(SamplingStrata)
library(readxl)

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load("data/processed/optimization_data.RData")
load("data/processed/grid_goa.RData")

##################################
## Load Survey Simulations for all scenarios
##################################
true_cv <- rrmse_cv <- array(dim = c(13, n_years, ns_all, n_boats),
                             dimnames = list(LETTERS[c(1:13)],
                                             paste0("year_", 1:n_years),
                                             common_names_all, 
                                             paste0("boat_", 1:n_boats)) )

load(paste0("results/scenario_A/Multispecies_Optimization/Str_15/",
            "simulation_results.RData"))
true_cv["A", , ,] <- STRS_true_cv_array
rrmse_cv["A", , ,] <- STRS_rrmse_cv_array

for (iscen in LETTERS[2:11]) {
  load(paste0("results/scenario_", iscen, "/Multispecies_Optimization/Str_5/",
              "simulation_results.RData"))
  
  true_cv[iscen, , , ] <- STRS_true_cv_array
  rrmse_cv[iscen, , , ] <- STRS_rrmse_cv_array
}

load(paste0("results/scenario_L/Multispecies_Optimization/Str_current/",
            "simulation_results.RData"))
true_cv["L", , , ] <- STRS_true_cv_array
rrmse_cv["L", , , ] <- STRS_rrmse_cv_array

load(paste0("results/scenario_M/Multispecies_Optimization/Str_current/",
            "simulation_results.RData"))
true_cv["M", , , ] <- STRS_true_cv_array
rrmse_cv["M", , , ] <- STRS_rrmse_cv_array

##################################################
####   Plot
##################################################
figs <- cbind( expand.grid(metric = c("true_cv", "rrmse_cv"),
                           boat = 1:3,
                           stringsAsFactors = FALSE),
               filename = c("appendix/Appendix C plots/Appendix C1",
                            "appendix/Appendix C plots/Appendix C2",
                            "figures/Figure08_true_cv",
                            "figures/Figure09_true_cv",
                            "appendix/Appendix C plots/Appendix C3",
                            "appendix/Appendix C plots/Appendix C4"))

for (irow in 1:nrow(figs)) { # Loop over metric -- start
  for (spp_group in 1:2) { # Loop over species groups -- start
    
    ## Open device
    png(filename = paste0(output_dir, figs$filename[irow], "_",
                          figs$metric[irow], "_boat_", figs$boat[irow], "_",
                          c("opt", "eval")[spp_group], "_spp.png"),
        width = 190, height = 210, units = "mm", res = 500)
    
    ## Format panels based on species group
    par(mfrow = switch(spp_group, "1" = c(5, 3), "2" = c(4, 3)), 
        mar = c(3, 2.5, 1, 0.5), oma = c(1, 2.5, 0, 0))
    
    ## Loop over species and 
    for (ispp in list(spp_idx_opt, 
                      spp_idx_eval)[[spp_group]] ) { # Loop over spp -- start
      
      ymax <- max(t(get(figs$metric[irow])[, , ispp, figs$boat[irow]]))
      
      boxplot(t(get(figs$metric[irow])[, , ispp, figs$boat[irow]]), 
              main = common_names_all[ispp],
              ylim = c(0, ymax),
              col = c(rep("white", 11), "grey", "grey"), pch = 16,
              at = c(1:11, 12:13), las = 1, cex.axis = 0.75)

    } # Loop over spp -- end
    
    ## Axis Titles
    mtext(side = 2, line = 1, outer = TRUE,
          text = switch(figs$metric[irow], 
                        "true_cv" = "True CV", 
                        "rrmse_cv" = "RRMSE of CV"),
          font = 2)
    
    mtext(side = 1, line = -0.5, outer = TRUE,
          text = "Design Scenario",
          font = 2)
    
    ## Close Device
    dev.off()
  } # Loop over species groups -- start
} # Loop over metric -- end

# par(mfrow = c(5, 3), mar = c(3, 3, 1, 1))
# for (ispp in spp_idx_opt) {
#   boxplot(t(true_cv[, , ispp]), 
#           main = common_names_all[ispp],
#           # ylim = c(0, max(t(true_cv[, , ispp]))), 
#           at = c(1:11, 13:14), las = 1); abline(v = 12)
# }
# 
# par(mfrow = c(4, 3), mar = c(3, 3, 1, 1))
# for (ispp in spp_idx_eval) {
#   boxplot(t(true_cv[, , ispp]), 
#           main = common_names_all[ispp],
#           # ylim = c(0, max(t(true_cv[, , ispp]))), 
#           at = c(1:11, 13:14), las = 1); abline(v = 12)
# }
# 
# par(mfrow = c(5, 3), mar = c(3, 3, 1, 1))
# for (ispp in spp_idx_opt) {
#   boxplot(t(rrmse_cv[, , ispp]), 
#           main = common_names_all[ispp],
#           # ylim = c(0, max(t(true_cv[, , ispp]))), 
#           at = c(1:11, 13:14), las = 1); abline(v = 12)
# }
# 
# par(mfrow = c(4, 3), mar = c(3, 3, 1, 1))
# for (ispp in spp_idx_eval) {
#   boxplot(t(rrmse_cv[, , ispp]), 
#           main = common_names_all[ispp],
#           # ylim = c(0, max(t(true_cv[, , ispp]))), 
#           at = c(1:11, 13:14), las = 1); abline(v = 12)
# }
# 
# ##################################################
# ####   Plot True CV and RRMSE
# ##################################################
# for(irow in 1:nrow(plot_scen)){
#   plot_filename <- paste0(output_dir, plot_scen$fig_name[irow])
#   imetric <- plot_scen$metric[irow]
#   iboat <- plot_scen$boat[irow]
#   
#   ## Open Device
#   png(filename = plot_filename,
#       width = 170, height = 200, units = "mm", res = 500)
#   
#   
#   ## Plot layout
#   layout(mat = matrix(c(1:ns_all, ns_all + 1, ns_all + 1), 
#                       ncol = 4, 
#                       byrow = TRUE))
#   par(mar = c(0.5, 3, 1.5, 1),
#       oma = c(0, 2.5, 0, 0))
#   
#   for (ispp in c(spp_idx_opt, spp_idx_eval)) { ## Loop over species -- start
#     
#     ## Calculate max y for plotting
#     merged_metric <- 
#       lapply(X = with(scen, paste0(domain, "_STR_", strata, "_", imetric )), 
#              FUN = function(x) get(x)[1:n_years,
#                                       common_names_all[ispp],
#                                       paste0("boat_", iboat) ])
#     ymax_ <- max(unlist(merged_metric)) * 1.1
#     
#     ## Plot metric
#     boxplot(merged_metric,
#             ylim = c(0, max(0.10, ymax_) ),
#             las = 1,
#             axes = F,
#             pch = 16,
#             col  = c("white", "cyan", "blue4", "chartreuse1", "darkgreen"),
#             cex = 0.5,
#             at = c(1, 3:4, 6:7))
#     box()
#     axis(side = 2, las = 1)
#     
#     ## Separate scenarios
#     abline(v = c(2, 5),
#            lty = "dotted", 
#            col = "black")
#     
#     ## Species label
#     mtext(side = 3, text = common_names_all[ispp], 
#           col = ifelse(ispp %in% spp_idx_opt, "black", "darkgrey"),
#           cex = 0.8)
#     
#   } ## Loop over species -- start
#   
#   ## Plot legend
#   par(mar = c(0, 0, 0, 0))
#   plot(x = 1, y = 1,
#        xlim = c(0,1), ylim = c(0,1),
#        type = "n", axes = F, ann = F)
#   
#   legend(x = 0,
#          y = 0.9,
#          legend = c("Exisiting STRS Design",
#                     "Area-Level Optimized STRS Design, 3 Strata per Area",
#                     "Area-Level Optimized STRS Design, 5 Strata per Area",
#                     "Gulf-Wide Optimized STRS Design, 10 Total Strata",
#                     "Gulf-Wide Optimized STRS Design, 15 Total Strata"),
#          fill = c("white", 
#                   "cyan",  "blue4",
#                   "chartreuse1", "darkgreen"),
#          xpd = NA,
#          cex = 1,
#          bty = "n")
#   
#   ## y-axis label
#   mtext(side = 2, 
#         outer = T, 
#         line = 0.5,
#         font = 2,
#         text = switch(imetric, 
#                       "rrmse_cv" = "RRMSE of CV",
#                       "true_cv" = "True CV") )
#   
#   ## Close Device
#   dev.off()
# }

