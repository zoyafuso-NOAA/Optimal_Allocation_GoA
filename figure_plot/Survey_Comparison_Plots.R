###############################################################################
## Project:       Synthesize Optimization Results
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Synthesize all optimization results
##                sample sizes, expected CVs, solutions, allocations, etc
###############################################################################
rm(list = ls())

##################################################
####  Set up directories
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[2]
VAST_model <- "11" 
which_domain <- c("full_domain", "trawlable")[1]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/",
                     which_domain, "/")

result_dir <- paste0(c("/Users/zackoyafuso/", 
                       "C:/Users/Zack Oyafuso/")[which_machine],
                     "Google Drive/MS_Optimizations/GF_Talk/figures/")

common_names <- c("arrowtooth flounder", "", "",
                  "rex sole", "", "",
                  "", "", "dover sole",
                  "Pacific Ocean perch", "rougheye/blackspotted RFs", 
                  "",
                  "", "", "shortspine thornyhead")

##################################################
####   Load Results
##################################################
load(paste0(github_dir, 
            "Spatiotemporal_Optimization/STRS_Sim_Res_spatiotemporal.RData"))
load(paste0(github_dir, 
            "Survey_Comparison_Simulations/Survey_Simulation_Results.RData"))

##################################################
####   True CV
##################################################
{png(filename = paste0(result_dir, "True_CV.png"),
     width = 10,
     height = 4,
     units = "in",
     res = 500)
  
  par(mar = c(1,4,1,1), 
      mfrow = c(2,3),
      oma = c(0,2,0,0))
  
  # for (ispp in c(1, 4, 9, 11, 15)) {
  for (ispp in 1:15) {
    temp_plot <- cbind(Curr = Current_true_cv_array[,ispp,2],
                       Opt = STRS_true_cv_array[,ispp,2,3])
    boxplot(temp_plot,
            ylim = c(0, max(temp_plot)),
            col = c("cadetblue1", "tomato2"),
            las = 1,
            names = NA,
            cex.axis = 1.25)
    mtext(side = 1, 
          line = -2,
          text = common_names[ispp])
    
  }
  
  plot(1,
       xlim = c(0, 1),
       ylim = c(0,1),
       type = "n",
       axes = F,
       ann = F)
  
  legend("center", 
         legend = c("Current Survey", "Optimized Survey"),
         fill = c("cadetblue1", "tomato2"),
         cex = 2)
  
  mtext(side = 2, 
        outer = T, 
        text = "True CV", 
        font = 2, 
        cex = 1.25,
        line = 0)
  
  dev.off()}

##################################################
####   RRMSE
##################################################
{png(filename = paste0(result_dir, "RRMSE_CV.png"),
     width = 10,
     height = 4,
     units = "in",
     res = 500)
  par(mar = c(1,4,1,1), 
      mfrow = c(2,3),
      oma = c(0,2,0,0))
  
  for (ispp in c(1, 4, 9, 11, 15)) {
  # for (ispp in 1:15) {
    temp_plot <- cbind(Curr = Current_rrmse_cv_array[,ispp,2],
                       Opt = STRS_rrmse_cv_array[,ispp,2,3])
    boxplot(temp_plot,
            ylim = c(0, max(temp_plot)),
            col = c("cadetblue1", "tomato2"),
            las = 1,
            names = NA,
            cex.axis = 1.25)
    mtext(side = 1, 
          line = -2,
          text = common_names[ispp])
    
  }
  
  plot(1,
       xlim = c(0, 1),
       ylim = c(0,1),
       type = "n",
       axes = F,
       ann = F)
  
  legend("center", 
         legend = c("Current Survey", "Optimized Survey"),
         fill = c("cadetblue1", "tomato2"),
         cex = 2)
  
  mtext(side = 2, 
        outer = T, 
        text = "True CV", 
        font = 2, 
        cex = 1.25,
        line = 0)
  
  dev.off()}
