###############################################################################
## Project:
## Author:
## Description:
###############################################################################
rm(list = ls())

##################################################
####  Set up directories
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")

output_dir <- paste0(c("/Users/zackoyafuso/", 
                       "C:/Users/Zack Oyafuso/")[which_machine],
                     "Google Drive/MS_Optimizations/TechMemo/figures")

##################################################
####   Load Data
##################################################
load(paste0(github_dir, "data/optimization_data.RData"))

##################################################
####   Constants
##################################################
common_names[11] <- "blackspotted/rougheye\nrockfishes"

##################################################
####   Knit all results
##################################################
settings <- data.frame()

for (iboat in 1:nboats) {
  result_dir <- paste0(github_dir, "results/Spatiotemporal_Optimization/",
                       "boat", iboat, "/")
  
  runs <- grep(x = dir(result_dir), 
               pattern = paste0("Str15Run"), 
               value = T )
  
  for (irun in runs) {
    temp_dir <- paste0(result_dir,  irun, "/result_list.RData")
    
    if (file.exists(temp_dir)) {
      load(temp_dir)
      
      #High-level settings: total sample size and expected CV across species
      species_cv <- result_list[[3]]
      attributes(species_cv)$dimnames[[1]] <- ""
      attributes(species_cv)$dimnames[[2]] <- paste0("CV_", 1:ns)
      cv <- max(as.numeric(species_cv))
      n <- result_list$n
      
      settings <- rbind(settings, 
                        data.frame(boat = iboat,
                                   n, cv, species_cv))
    }
  }
}

##################################################
####   
##################################################

{
  png(filename = paste0(output_dir, "/CV_N_bySpp.png"),
      width = 190,
      height = 190,
      units = "mm",
      res = 500)
  
  par(mfrow = c(5, 3),
      mar = c(0, 4, 0, 0),
      oma = c(4, 2, 1, 1))
  
  spp_order <- order(settings[1, paste0("CV_", 1:ns)])
  
  for (ispp in spp_order) {
    
    plot(x = 1,
         y = 1,
         type = "n",
         xlim = c(0, 900), 
         xlab = "",
         ylim = c(0, 1.25 * max(settings[, paste0("CV_", ispp)])), 
         ylab = "",
         las = 1,
         axes = F)
    box()
    abline(v = c(280, 550, 820),
           lty = "dashed",
           col = "grey")
    
    axis(side = 2, las = 1)
    if(ispp %in% spp_order[(ns-2):ns]) axis(side = 1)
    
    points(x = settings$n,
           y = settings[, paste0("CV_", ispp)],
           pch = 16)
    lines(x = settings$n,
          y = settings[, paste0("CV_", ispp)])
    legend("bottomright", 
           legend = common_names[ispp],
           bty = "n",
           cex = 1.25)
  }
  
  mtext(side = 1, "Optimized Sample Size", line = 2.5, outer = T, font = 2)
  mtext(side = 2, "Upper CV Limit", line = 0, outer = T, font = 2)
  
  dev.off()
}
