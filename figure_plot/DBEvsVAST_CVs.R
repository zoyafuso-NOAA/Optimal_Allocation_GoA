###############################################################################
## Project:         Model vs design-based index estimates
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

# load(paste0(github_dir, 
#             "results/Spatiotemporal_Optimization/",
#             "optimization_knitted_results.RData"))

# load(paste0(github_dir, 
#             "data/Population_Variances.RData"))

##################################################
####  Import Observed DBEs, calculate ranges of sample CVs
##################################################
DBE <- readRDS(paste0(github_dir2, 
                      "GOA_biomass_indices_wnames.rds"))
DBE <- subset(DBE, 
              SPECIES_NAME %in% sci_names_all &
                YEAR %in% (1996:2019)[years_included])

DBE$BIOMASS_CV <- sqrt(DBE$VAR_WGT_CPUE) / DBE$MEAN_WGT_CPUE
DBE <- DBE[order(DBE$YEAR),]

load(paste0(github_dir, 
            "data/vast_index.RData"))
vast_index <- subset(vast_index, 
                     spp %in% sci_names_all)

##################################################
####  Plot
##################################################

{
  png(filename = paste0(output_dir, "DBEvsVAST_CVs.png"),
      width = 190,
      height = 210,
      units = "mm",
      res = 500)
  
  par(mfrow = c(7, 3),
      mar = c(0, 2.5, 1, 1),
      oma = c(4, 3.5, 0.5, 0))
  
  for (ispp in c(spp_idx_opt[-11], spp_idx_eval[-c(2:3)]) ) {
    
    #Knit the VAST and DBE output into one list
    temp_list <- list(
      DBE = subset(DBE, SPECIES_NAME == sci_names_all[ispp], 
                   select = c(YEAR, TOTAL_BIOMASS, BIOMASS_VAR)),
      VAST = subset(vast_index, 
                    subset = spp == sci_names_all[ispp],
                    select = c(year, est, cv))
    )
    temp_list$DBE$BIOMASS_SD <- sqrt(temp_list$DBE$BIOMASS_VAR)
    temp_list$VAST$sd <- with(temp_list$VAST, est*cv)
    
    #Base plot
    plot(1,
         type = "n",
         las = 1,
         ann = F,
         ylim = c(0, 1.25 * max(cbind(with(temp_list$DBE, 
                                           BIOMASS_SD/TOTAL_BIOMASS),
                                      temp_list$VAST$cv))),
         xlim = c(1995, 2020),
         axes = F)
    
    # Year label for the last row of plots
    if(ispp %in% c(16, 20, 22)) axis(side = 1, 
                                         at = seq(1995, 2020, by = 5))
    axis(side = 2, las = 1)
    box()
    
    matpoints(x = temp_list$DBE$YEAR,
              y = cbind(with(temp_list$DBE, BIOMASS_SD/TOTAL_BIOMASS),
                        temp_list$VAST$cv) ,
              pch = c(16, 16),
              col = c("red", "black"))
    matlines(x = temp_list$DBE$YEAR,
             y = cbind(with(temp_list$DBE, BIOMASS_SD/TOTAL_BIOMASS),
                       temp_list$VAST$cv) ,
             lty = 1,
             col = c("red", "black"))
    # boxplot(cbind(with(temp_list$DBE, BIOMASS_SD/TOTAL_BIOMASS),
    #               temp_list$VAST$cv),
    #         add = T,
    #         at = c(2022, 2024),
    #         col = c("red", "white"),
    #         axes = F,
    #         pch = 16)
    
    mtext(side = 3, 
          text = common_names_all[ispp], 
          line = -1.5,
          col = ifelse(ispp %in% spp_idx_eval, "darkgrey", "black"),
          cex = 0.75)
    
  }
  
  mtext(side = 1, 
        outer = T, 
        text = "Year", 
        line = 3,
        font = 2)
  mtext(side = 2, 
        outer = T, 
        text = "Sample CV of Abundance Index", 
        line = 1.5,
        font = 2)
  
  #Plot Legend
  par(mar = c(0,2,3,2))
  plot(1, 
       type = "n", 
       axes = F, 
       ann = F)
  legend("bottom", 
         legend = c("DBE", "VAST"), 
         col = c("red", "black"), 
         pch = 16, 
         lty = 1, 
         cex = 1.5)
  
  dev.off()
}
