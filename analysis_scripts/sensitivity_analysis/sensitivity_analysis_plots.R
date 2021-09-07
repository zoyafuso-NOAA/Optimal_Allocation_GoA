setwd("C:/Users/Zack Oyafuso/Documents/GitHub/Optimal_Allocation_GoA/")

library(raster); library(sp)
library(RColorBrewer)

load("data/processed/optimization_data.RData")
load("data/processed/grid_goa.RData")

scen <- expand.grid(strata = c(3, 5),
                    type = 0:4,
                    strata_vars = 1:2,
                    deep_stations = c(TRUE, FALSE))

{
  pdf("C:/Users/Zack Oyafuso/Google Drive/MS_Optimizations/Update Documents/Updates 7 September 2021/sensitivity.pdf", width = 8, height = 11, onefile = TRUE)
  
  for (strata_vars in 2:1) {
    par(mar = c(0, 0, 0, 0), mfrow = c(2, 2))
    
    for (istrata in c(3, 5)) {
      for (deep_stations in c(T, F)) {
        plot(1, type = "n", ann = F, axes = F, asp = 1,
             xlim = range(grid_goa$E_km),
             ylim = c(min(grid_goa$N_km), max(grid_goa$N_km) + y_range*2.25) )
        box()
        mtext(side = 3, line = -3, text = paste0(istrata, " strata per area\n",
                                                 ifelse(deep_stations == TRUE,
                                                        "Depths up to 1000 m",
                                                        "Depths > 300 m -> 300 m")))
        
        y_offset = 0
        
        for(itype in c(3, 4, 1, 0)) {
          load(paste0("results/sensitivity_analysis/",
                      "Str_", istrata, 
                      "/Type", itype, 
                      "_stratavars", strata_vars, "_deep", 
                      deep_stations, "/result_list.RData"))
          
          
          
          ##Save a plot of the solution
          goa <- sp::SpatialPointsDataFrame(
            coords = grid_goa[, c("E_km", "N_km")],
            data = data.frame(Str_no = result_list$sol_by_cell) )
          goa_ras <- raster::raster(x = goa, 
                                    resolution = 5)
          goa_ras <- raster::rasterize(x = goa, 
                                       y = goa_ras, 
                                       field = "Str_no")
          goa_ras <- raster::shift(x = goa_ras, dy = y_range * 0.7 * y_offset)
          
          cols <- rep(rev(brewer.pal(n = 9, name = "Paired")), 3)
          
          image(goa_ras, asp = 1, axes = FALSE, ann = F, add = T, col = cols)
          
          text(x = min(grid_goa$E_km) + x_range*0.7,
               y = min(grid_goa$N_km) + y_range*0.55 + y_range * 0.7 * y_offset, 
               labels = switch(paste0(itype),
                               "0" = "Predicted\nDensity",
                               "1" = "Simulate\nMeasurement Error",
                               "3" = "Simulate Random\nand Fixed Effects",
                               "4" = "Simulate\nRandom Effects"),
               cex = 0.9)
          
          y_offset = y_offset + 1
          
        }
        
        
      }
    }
  }
  
  
  dev.off()
}
