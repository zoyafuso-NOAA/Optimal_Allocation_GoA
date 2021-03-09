###############################################################################
## Project:         Plot Single Species Optimization Solutions
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     
###############################################################################
rm(list = ls())

##################################################
####  Set up directories  
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")

output_dir <- paste0(c("/Users/zackoyafuso/", 
                       "C:/Users/Zack Oyafuso/")[which_machine], 
                     "Google Drive/MS_Optimizations/TechMemo/figures/")

##################################################
####  Import Libraries  
##################################################
library(sp)
library(raster)
library(RColorBrewer)

##################################################
####    Load predicted density and optimization results 
##################################################
load(paste0(github_dir, 
            "data/optimization_data.RData"))
load(paste0(github_dir, 
            "/data/Extrapolation_depths.RData"))
load(paste0(github_dir, "results/full_domain/Single_Species_Optimization/", 
            "optimization_knitted_results.RData"))

settings <- settings_agg_full_domain
res_df <- res_df_full_domain
strata_list <- strata_list_full_domain

##################################################
####   Set Constants
##################################################
common_names_opt <- c("arrowtooth\nflounder", "walleye pollock", "Pacific cod", 
                      "rex sole", "flathead sole", "Pacific halibut", 
                      "southern\nrock sole", "northern\nrock sole", "Dover sole",
                      "Pacific ocean\nperch", "blackspotted/\nrougheye\nrockfishes",
                      "silvergrey\nrockfish", "northern\nrockfish", 
                      "dusky\nrockfish",  "shortspine\nthornyhead")
common_names_eval_labels <- c("sablefish", "skates spp.", "Octopus spp.", 
                              "Atka mackerel", "shortraker\nrockfish", 
                              "harlequin\nrockfish", "Pacific\nspiny dogfish")

idomain <- "full_domain"

isample <- 2 #1, 2, or 3 boat solution
x_range <- diff(range(Extrapolation_depths$E_km))
y_range <- diff(range(Extrapolation_depths$N_km))

##################################################
####   Plot
##################################################
{png(filename = paste0(output_dir, "SS_solutions.png"),
     height = 170,
     width = 150,
     units = "mm",
     res = 500)
  
  ## Plot setup
  par(mfrow = c(1, 2),
      mar = c(0, 0, 0, 0))
  
  for (ispp in 1:ns_all) {
    if (ispp%%11 == 1) { 
      
      ## Base Plot
      plot(1, 
           type = "n",
           xlim = range(Extrapolation_depths$E_km),
           ylim = with(Extrapolation_depths, 
                       c(min(N_km) + 0*y_range, max(N_km) + 5.75*y_range)), 
           axes = F, 
           ann = F,
           asp = 1)
      # box()
      
      offset_val <- 1
    }
    
    ## Which index to plot
    idx = which(settings$boat == isample &
                  settings$spp == c(spp_idx_opt, spp_idx_eval)[ispp])
    
    #Plot Solution
    goa <- SpatialPointsDataFrame(
      coords = Extrapolation_depths[, c("E_km", "N_km")],
      data = data.frame(Str_no = res_df[, idx]) )
    
    goa_ras <- raster(goa,
                      resolution = 5)
    goa_ras <- rasterize(x = goa,
                         y = goa_ras,
                         field = "Str_no")
    offset_y <- 0.6 * y_range * (offset_val-1)
    goa_ras <- raster::shift(goa_ras, dy = offset_y )
    offset_val <- offset_val + 1
    
    n_strata <- length(unique(res_df[, idx]))
    
    image(goa_ras,
          axes = F,
          ann = F0,
          asp = 1,
          col = colorRampPalette(
            colors = brewer.pal(n = 10, name = 'Paired'))(n_strata) ,
          add = T)
    
    text(x = min(Extrapolation_depths$E_km) + x_range * 0.72,
         y = min(Extrapolation_depths$N_km) + offset_y + y_range*0.6,
         labels = c(common_names_opt, common_names_eval_labels)[ispp],
         cex = 0.5)
    
    #Simulate a sample solution
    temp_samples <- c()
    temp_strata <- nrow(strata_list[[idx]])
    temp_solution <- res_df[, idx ]
    temp_allocation <- strata_list[[idx]]$Allocation
    
    for (istrata in 1:temp_strata) {
      temp_samples = c(temp_samples,
                       sample(x = which(temp_solution == istrata),
                              size = temp_allocation[istrata])
      )
    }
    
    temp_loc <- Extrapolation_depths[temp_samples, c("E_km", "N_km")]
    temp_loc$N_km <- temp_loc$N_km + offset_y
    
    points(temp_loc,
           pch = 16,
           cex = 0.25)
    
  }
  dev.off()
}
