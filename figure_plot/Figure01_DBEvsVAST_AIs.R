###############################################################################
## Project:         VASt Model vs design-based index estimates
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     Plot abundance indices and SDs for each species included
##                  in the operating model
###############################################################################
rm(list = ls())

##################################################
####  Set up directories  
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[2]

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

##################################################
####  Import Observed DBEs, calculate ranges of sample CVs
##################################################
load(paste0(github_dir, 
            "data/VAST_fit_I_ct.RData"))
vast_index$species[vast_index$species == "spiny dogfish"] <-
  "Pacific spiny dogfish"

DBE <- readRDS(paste0(github_dir2, 
                      "GOA_biomass_indices_wnames.rds"))
DBE$COMMON_NAME[DBE$COMMON_NAME == "rougheye and blackspotted rockfish"] <-
  "BS and RE rockfishes"
DBE$COMMON_NAME[DBE$COMMON_NAME == "spiny dogfish"] <- "Pacific spiny dogfish"

DBE <- subset(DBE, 
              COMMON_NAME %in% common_names_all &
                YEAR %in% (1996:2019)[years_included])
DBE <- DBE[order(DBE$YEAR),]

##################################################
####  Plot
##################################################

{
  png(filename = paste0(output_dir, "Figure_01_DBEvsVAST_AIs.png"),
      width = 190,
      height = 210,
      units = "mm",
      res = 500)
  
  par(mfrow = c(7, 4),
      mar = c(0, 2.5, 1, 1),
      oma = c(3, 2.5, 0.5, 0))
  
  for (ispp in c(common_names_opt, common_names_eval) ) {
    
    #Knit the VAST and DBE output into one list
    temp_list <- list(
      DBE = subset(DBE, COMMON_NAME == ispp, 
                   select = c(YEAR, TOTAL_BIOMASS, BIOMASS_VAR)),
      VAST = subset(vast_index, 
                    subset = species == ispp,
                    select = c(Year, Estimate_metric_tons, SD_mt))
    )
    
    temp_list$DBE$BIOMASS_SD <- sqrt(temp_list$DBE$BIOMASS_VAR)
    
    #Base plot
    plot(1,
         type = "n",
         las = 1,
         ann = F,
         ylim = c(0, 1.5 * 0.001 * max(cbind(temp_list$DBE$TOTAL_BIOMASS,
                                             temp_list$VAST$Estimate_metric_tons))),
         xlim = c(1995, 2020),
         axes = F)
    
    # Year label for the last row of plots
    if(ispp %in% common_names_eval[8:11]) 
      axis(side = 1, at = seq(1995, 2020, by = 5))
    axis(side = 2, las = 1)
    box()
    
    matpoints(x = temp_list$DBE$YEAR,
              y = 0.001 * cbind(temp_list$DBE$TOTAL_BIOMASS,
                                temp_list$VAST$Estimate_metric_tons) ,
              pch = c(16, 16),
              col = c("red", "black"))
    
    matlines(x = temp_list$DBE$YEAR,
             y = 0.001 * cbind(temp_list$DBE$TOTAL_BIOMASS,
                               temp_list$VAST$Estimate_metric_tons) ,
             lty = 1,
             col = c("red", "black"))
    
    segments(x0 = temp_list$DBE$YEAR,
             x1= temp_list$DBE$YEAR,
             y0 = 0.001 * with(temp_list$DBE, TOTAL_BIOMASS - BIOMASS_SD),
             y1 = 0.001 * with(temp_list$DBE, TOTAL_BIOMASS + BIOMASS_SD),
             col = "red")
    
    segments(x0 = temp_list$VAST$Year,
             x1= temp_list$VAST$Year,
             y0 = 0.001 * with(temp_list$VAST, Estimate_metric_tons - SD_mt),
             y1 = 0.001 * with(temp_list$VAST, Estimate_metric_tons + SD_mt),
             col = "black")
    
    
    
    mtext(side = 3, 
          text = ispp, 
          line = -1.5,
          col = ifelse(ispp %in% common_names_eval, "darkgrey", "black"),
          cex = 0.75)
    
  }
  
  mtext(side = 1, 
        outer = T, 
        text = "Year", 
        line = 2,
        font = 2)
  mtext(side = 2, 
        outer = T, 
        text = "Abundance Index (+/- 1 SD, units: thousand metric tons)", 
        line = 1,
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
