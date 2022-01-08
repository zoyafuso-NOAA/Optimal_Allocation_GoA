library(raster); 
library(sp)
library(RColorBrewer)

load("data/processed/optimization_data.RData")
load("data/processed/grid_goa.RData")

{
  for (strata_vars in 2) {
    
    png(paste0("C:/Users/Zack Oyafuso/Google Drive/MS_Optimizations/", 
               "past_talks/PlanTeam_2021/sensitivity_type3_", 
               strata_vars, "vars.png"), 
        width = 7, height = 12, units = "in", res = 500)
    
    par(mar = c(0, 0, 0, 0), mfrow = c(1, 1))
    
    plot(1, type = "n", ann = F, axes = F, asp = 1,
         xlim = range(grid_goa$E_km),
         ylim = c(min(grid_goa$N_km), max(grid_goa$N_km) + y_range*5) )
    # box()
    
    y_offset = 0
    for(isim in 1:10) {
      load(paste0("results/sensitivity_analysis_type3/Str_3/sim", isim, 
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
      goa_ras <- raster::shift(x = goa_ras, dy = y_range * 0.55 * y_offset)
      
      cols <- rep(rev(brewer.pal(n = 9, name = "Paired")), 3)
      image(goa_ras, asp = 1, axes = FALSE, ann = F, add = T, col = cols)

      y_offset = y_offset + 1
      
    }
    dev.off()
  }
}
