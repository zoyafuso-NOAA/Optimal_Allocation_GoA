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
DBE <- DBE[order(DBE$YEAR),]

load(paste0(github_dir, 
            "data/vast_index.RData"))

##################################################
####  Plot
##################################################

{
 png(filename = paste0(output_dir, "DBEvsVAST_CVs.png"),
     width = 190,
     height = 190,
     units = "mm",
     res = 500)
 
 par(mfrow = c(5, 3),
     mar = c(0, 3, 1, 1),
     oma = c(4, 3, 0.5, 0))
 
 for (ispp in (1:ns)[-11] ) {
  temp_list <- list(DBE = subset(DBE, SPECIES_NAME == sci_names[ispp], 
                                 select = c(YEAR, TOTAL_BIOMASS, BIOMASS_VAR)),
                    VAST = subset(vast_index, 
                                  subset = spp == sci_names[ispp],
                                  select = c(year, est, cv))
  )
  temp_list$DBE$BIOMASS_SD <- sqrt(temp_list$DBE$BIOMASS_VAR)
  temp_list$VAST$sd <- with(temp_list$VAST, est*cv)
  plot(1,
       type = "n",
       las = 1,
       ann = F,
       ylim = c(0, 1.1 * max(cbind(with(temp_list$DBE, 
                                       BIOMASS_SD/TOTAL_BIOMASS),
                                  temp_list$VAST$cv))),
       xlim = c(1995, 2025),
       axes = F)
  
  if(ispp %in% 13:15) axis(side = 1, 
                           at = seq(1995, 2020, by = 5), 
                           cex.axis = 0.75)
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
  boxplot(cbind(with(temp_list$DBE, BIOMASS_SD/TOTAL_BIOMASS),
                temp_list$VAST$cv),
          add = T,
          at = c(2022, 2024),
          col = c("red", "white"),
          # border = c("black", "darkgrey"),
          axes = F,
          pch = 16)
  mtext(side = 1, text = common_names[ispp], line = -1.5)
  
 }
 
 mtext(side = 1, 
       outer = T, 
       text = "Year", 
       line = 3,
       font = 2)
 mtext(side = 2, 
       outer = T, 
       text = "Sample CV of Abundance Index", 
       line = 1,
       font = 2)
 
 #Plot Legend
 par(mar = c(1,2,3,2))
 plot(1, 
      type = "n", 
      axes = F, 
      ann = F)
 legend("bottom", 
        legend = c("DBE", "VAST"), 
        col = c("red", "black"), 
        pch = 16, 
        lty = 1, 
        cex = 2)

 dev.off()
}
