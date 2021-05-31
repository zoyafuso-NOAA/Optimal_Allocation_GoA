###############################################################################
## Project:       Survey Feasibility (total distance travelled)
## Author:        Megsie Siple (margaret.siple@noaa.gov)
##                Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Simulate several 2-boat surveys station locations
##                Calculate shortest distance path across two boats
##                Calculate total distance travelled
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(TSP)

##################################################
####   Set up directories
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents", 
                       "C:/Users/Zack Oyafuso/Documents",
                       "C:/Users/zack.oyafuso/Work")[which_machine],
                     "/GitHub/Optimal_Allocation_GoA/")

output_dir <- paste0(c("/Users/zackoyafuso/",
                       "C:/Users/Zack Oyafuso/")[which_machine],
                     "Google Drive/MS_Optimizations/TechMemo/figures/")

##################################################
####   Import result objects
##################################################
load(file = paste0(github_dir, '/data/Extrapolation_depths.RData'))
load(file = paste0(github_dir, "/results/historical_surveys_distances.RData"))
load(file = paste0(github_dir, '/results/MS_optimization_knitted_results.RData'))

##################################################
####   Import functions to conduct TSP
##################################################
source(file = paste0(github_dir, '/modified_functions/do_TSP.R'))
source(file = paste0(github_dir, '/modified_functions/sample_stations.R'))

##################################################
####   Draw stations under a district-level 3-strata STRS (idx = 2)
####   Calculate distance matrix of each drawn station for the first boat
##################################################
idx <- 2
settings[idx, ]

set.seed(54662)
temp_locs <- sample_stations(strata_labels = strata_list[[idx]]$stratum_id,
                             nh = strata_list[[idx]]$Allocation,
                             grid_strata = res_df[, idx ])
dist_km_y_b <- as.matrix(dist(temp_locs[, c("E_km", "N_km")]))

##################################################
####   Randomly assign stations to boats (2 boats in total)
####   Calculate the shortest path that visits each station once using TSP
##################################################
boat_idx <- sample(x = 1:2, size = nrow(temp_locs), replace = TRUE)
temp_path <- do_tsp(locs = temp_locs[boat_idx == 1, ])

##################################################
####   Using the station order, calculate 
##################################################
opt_nearest_dists <- NULL
take_out_idx <- which(boat_idx == 2)

for(istation in temp_path$path) {
  
  opt_nearest_dists <- rbind(opt_nearest_dists, 
                             sort(dist_km_y_b[istation, -take_out_idx])[2:3]  )
  
  ## Query stations that have already been sampled as well as those not 
  ## assigned to boat 1
  take_out_idx <- c(take_out_idx, istation)
}

##################################################
####   Plot
##################################################

{
  ## Open Device
  png(filename = paste0(output_dir, "Figure11_closest_station_dists.png"),
      width = 6.5, height = 3, units = "in", res = 500)
  
  ## Plot Layout
  par(mfrow = c(1, 2), 
      mar = c(2, 2.5, 1, 1), 
      oma = c(2, 2, 1, 1))
  
  for (idist in 1:2) { ## Loop through closest and 2nd closest dists -- start
    
    ## Base Plot
    plot(x = 1, y = 1, type = "n", axes = F, ann = F,
         xlim = c(1995, 2030), ylim = c(0, 300), las = 1)
    box()
    axis(side = 1, at = c(seq(from = 1995, to = 2020, by = 5), 2028 ), 
         labels = c(seq(from = 1995, to = 2020, by = 5), "Opt. Survey"))
    axis(side = 2, las = 1)
    mtext(side = 3, 
          text = c("Distance to Closest Station",
                   "Distance to Second-Closest Station",
                   font = 2)[idist],
          line = 0.5, font = 2)
    
    ## Loop through years and create a boxplot of nearest distances
    for (iyear in c(1996, 1999, 2003, 2005, 
                    2007, 2009, 2011, 2013, 
                    2015, 2017, 2019)) { ## Loop years -- start
      plot_this <- nearest_dists$actual_path[[paste0("year_", iyear)]]
      boxplot(plot_this$boat_1[, idist],
              add = T,
              at = iyear,
              axes = F, width = 1,
              pch = 16, cex = 0.25)
    } ## Loop years -- end
    
    ## Plot distribution of the optimal solution
    plot_this <- opt_nearest_dists[, idist]
    boxplot(plot_this,
            add = T,
            at = 2028,
            axes = F, width = 1,
            pch = 16, cex = 0.25)
    
  } ## Loop through closest and second closest dists -- end
  
  ## Outer labels
  mtext(side = 1, text = "Year", outer = TRUE, font = 2, line = 0.5)
  mtext(side = 2, text = "Distance (km)", outer = TRUE, font = 2, line = 0.5)
  
  ## Close Device
  dev.off()
}
