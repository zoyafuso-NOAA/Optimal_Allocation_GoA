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
common_names[11] <- "blackspotted/rougheye\nrockfishes" 

load(paste0(github_dir, 
            "results/Spatiotemporal_Optimization/",
            "optimization_knitted_results.RData"))

load(paste0(github_dir, 
            "results/Population_Variances.RData"))

##################################################
####  Import Observed DBEs, calculate ranges of sample CVs
##################################################
DBE <- readRDS(paste0(github_dir2, 
                      "GOA_biomass_indices_wnames.rds"))
DBE <- subset(DBE, 
              SPECIES_NAME %in% sci_names &
                YEAR %in% (1996:2019)[Years2Include])

DBE$BIOMASS_CV <- sqrt(DBE$VAR_WGT_CPUE) / DBE$MEAN_WGT_CPUE

DBE_CV <- aggregate(BIOMASS_CV ~ SPECIES_NAME, 
                    data = DBE,
                    FUN = range)$BIOMASS_CV

load(paste0(github_dir, 
            "data/vast_index.RData"))

vast_ranges <- aggregate(cv ~ spp, 
                         data = vast_index,
                         FUN = range)$cv

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
     ylim = c(1, ns),
     xlim = c(0, 0.75),
     axes = F,
     ann = F)

box()
abline(h = 1:ns, 
       col = "lightgrey", 
       lty = "dashed")
mtext(side = 1, 
      text = "Population CV",
      line = 2.5,
      cex = 1.25,
      font = 2)

segments(x0 = vast_ranges[, 1],
         x1 = vast_ranges[, 2],
         y0 = (1:ns) + 0.2 , 
         lwd = 2)

segments(x0 = DBE_CV[, 1],
         x1 = DBE_CV[, 2],
         y0 = (1:ns) - 0.2 , 
         lwd = 2,
         col = "brown")

matpoints(y = 1:ns,
          x = cbind(SRS_Pop_CV[, 2], 
                    SS_STRS_Pop_CV[, 2],
                    unlist(settings[settings$strata == 15 &
                                      settings$boat == 2, 
                                    paste0("CV_", 1:ns)]),
                    Current_STRS_Pop_CV[, 2]),
          pch = 16,
          col = c("black", "red", "blue", "green"),
          cex = 1)

axis(side = 1, at = seq(from = 0, to = 0.8, by = 0.05))
axis(side = 2, 
     labels = common_names, 
     at = 1:ns,
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