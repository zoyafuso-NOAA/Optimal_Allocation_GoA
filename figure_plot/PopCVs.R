###############################################################################
## Project:         Plot Population CVs of different sampling designs
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     
###############################################################################
rm(list = ls())

##################################################
####  Set up directories  
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[1]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")
github_dir2 <- paste0(c("/Users/zackoyafuso/Documents/", 
                        "C:/Users/Zack Oyafuso/Documents/")[which_machine], 
                      "GitHub/MS_OM_GoA/data/")

output_dir <- paste0(c("/Users/zackoyafuso/", 
                       "C:/Users/Zack Oyafuso/")[which_machine], 
                     "Google Drive/MS_Optimizations/TechMemo/figures/")

##################################################
####  Load Data
##################################################
load(paste0(github_dir, 
            "data/optimization_data.RData"))
load(paste0(github_dir, 
            "data/fit_density.RData"))
load(paste0(github_dir, 
            "results/Spatiotemporal_Optimization/",
            "optimization_knitted_results.RData"))

load(paste0(github_dir, 
            "results/Population_Variances.RData"))

##################################################
####  Plot
##################################################

png(filename = paste0(output_dir, "PopCVs.png"),
    width = 190,
    height = 150,
    units = "mm",
    res = 500)

par(mar = c(4,11,1,1))
plot(1,
     type = "n",
     pch = 16,
     cex = 1.5,
     ylim = c(1, ns_opt),
     xlim = c(0, 0.4),
     axes = F,
     ann = F)

box()
abline(h = 1:ns_opt, 
       col = "lightgrey", 
       lty = "dashed")
mtext(side = 1, 
      text = "Population CV",
      line = 2.5,
      cex = 1.25,
      font = 2)

matpoints(y = 1:ns_opt,
          x = cbind(SRS_Pop_CV[spp_idx_opt, 2], 
                    SS_STRS_Pop_CV[, 2],
                    unlist(settings[settings$strata == 15 &
                                      settings$boat == 2,
                                    paste0("CV_", 1:ns_opt)]),
                    Current_STRS_Pop_CV[spp_idx_opt, 2]),
          pch = 16,
          col = c("black", "red", "blue", "green"),
          cex = 1)

axis(side = 1, at = seq(from = 0, to = 0.8, by = 0.05))
axis(side = 2, 
     labels = common_names_opt, 
     at = 1:ns_opt,
     las = 1)

legend("bottomright",
       pch = c(rep(16, 4), NA, NA),
       lty = c(rep(NA, 4), 1, 1), 
       lwd = 2,
       col = c("red", "green", "blue", "black", "black", "brown"),
       legend = c("SS STRS", "Current STRS", "Optimized STRS", "SRS", 
                  "Est. VAST CVs",
                  "Est. DBE CVs"),
       cex = 1)

dev.off()