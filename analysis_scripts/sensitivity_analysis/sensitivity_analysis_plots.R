library(raster); 
library(sp)
library(RColorBrewer)

load("data/processed/optimization_data.RData")
load("data/processed/grid_goa.RData")

{
  for (strata_vars in 2) {
    
    png(paste0("C:/Users/Zack Oyafuso/Google Drive/MS_Optimizations/", 
               "past_talks/PlanTeam_2021/sensitivity_by_type_", 
               strata_vars, "vars.png"), 
        width = 8, height = 7, units = "in", res = 500)
    
    par(mar = c(0, 0, 0, 0), mfrow = c(1, 1))
    
    plot(1, type = "n", ann = F, axes = F, asp = 1,
         xlim = range(grid_goa$E_km),
         ylim = c(min(grid_goa$N_km), max(grid_goa$N_km) + y_range*1.5) )
    box()
    
    y_offset = 0
    for(itype in c(3, 1, 0)) {
      load(paste0("results/sensitivity_analysis/",
                  "Str_3/Type", itype, 
                  "_stratavars", strata_vars, "_deepTRUE/result_list.RData"))
      
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
           y = min(grid_goa$N_km) + y_range*0.6 + y_range * 0.7 * y_offset, 
           labels = switch(paste0(itype),
                           "0" = "Predicted\nDensity",
                           "1" = "Simulate\nMeasurement Error",
                           "3" = "Simulate Random\nand Fixed Effects"),
           cex = 1.25)
      
      y_offset = y_offset + 1
      
    }
    dev.off()
  }
}
