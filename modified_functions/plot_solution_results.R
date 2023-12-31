###############################################################################
## Project:       Plot function 
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Create function to output basic solution maps
###############################################################################

plot_solution_results <- function( file_name,
                                   grid_object,
                                   districts_object,
                                   district_values,
                                   sol_by_cell,
                                   draw_stations,
                                   allocations){
  
  ## Import libraries
  # library(sp); library(raster)
  library(terra)
  
  ## Set up spatial object
  goa <- terra::vect(x = cbind(grid_object,
                               Sol = sol_by_cell), 
                     geom = c("Eastings", "Northings"))
  goa_ras <- terra::rast(x = goa, resolution = 10000)
  goa_ras <- terra::rasterize(x = goa, y = goa_ras, field = "Sol")
  
  ## Set up plot
  png(filename = file_name, width = 6, height = 3, units = "in", res = 500)
  
  ## Set up panel layout
  par(mfrow = c(1, 1), mar = c(1, 1, 1, 1))
  
  ## Plot solution
  strata_colors <- colorRampPalette(
    brewer.pal(n = 11,
               name = "Paired"))(length(x = unique(x = plot_solution)) )
  
  # if (any(temp_ids == 0)) strata_colors <- c("black", strata_colors)
  
  plot(goa_ras, axes = F, asp = 1, col = sample(strata_colors) )
  
  ## Draw boundaries of management areas
  # rect(xleft = districts_object$W_UTM,
  #      xright = districts_object$E_UTM,
  #      ybottom = tapply(X = grid_object$N_km,
  #                       INDEX = district_values,
  #                       FUN = min),
  #      ytop = tapply(X = grid_object$N_km,
  #                    INDEX = district_values,
  #                    FUN = max))
  # 
  # text(x = rowMeans(districts_object[, c("W_UTM", "E_UTM")]),
  #      y = tapply(X = grid_object$N_km,
  #                 INDEX = district_values,
  #                 FUN = max),
  #      labels = districts_object$district,
  #      pos = 3)
  # box()
  
  ## Draw samples if needed
  # if (draw_stations) {
  #   
  #   ## Take a random sample based on the allocation and stratum
  #   sample_vec <- c()
  #   for(istrata in 1:length(allocations)) {
  #     sample_vec <- c(sample_vec,
  #                     sample(x = which(sol_by_cell == istrata),
  #                            size = allocations[istrata]) )
  #   }
  #   
  #   points(grid_object[sample_vec, c("E_km", "N_km")],
  #          pch = 16, cex = 0.5)
  # }
  
  ## Close device
  dev.off()
}
