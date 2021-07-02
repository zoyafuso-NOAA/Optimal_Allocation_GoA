###############################################################################
## Project:         MS Solutions with strata boundary breakdowns
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     Plot MS solutions with where the boundaries are with 
##                  respect to bathymetry and longitude
###############################################################################
rm(list = ls())

##################################################
####  Set up directories  
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[2]
output_dir <- paste0(c("/Users/zackoyafuso/Google Drive/",
                       "C:/Users/Zack Oyafuso/Google Drive/")[which_machine],
                     "MS_Optimizations/TechMemo/figures/")

##################################################
####  Imported Libraries
##################################################
library(RColorBrewer)

##################################################
####  Transparent colors fn from Mark Gardener 2015 www.dataanalytics.org.uk
##################################################
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

##################################################
####  Load Data
##################################################
load("data/processed/optimization_data.RData")
load("data/processed/grid_goa.RData")
load("results/MS_optimization_knitted_results.RData")

##################################################
####  Random seed for consistent colors
##################################################
seed <- 234233

##################################################
## Rescale UTM eastings back to UTM 5N
##################################################
strata_list <- 
  lapply(X = strata_list,
         FUN = function(x) {
           x[, c("Lower_X1", "Upper_X1")] <- 
             x[, c("Lower_X1", "Upper_X1")] + min(grid_goa$E_km)
           
           temp_idx <- match(x = round(x$Lower_X1), 
                             table = round(grid_goa$E_km))
           x$Lower_X1 <- grid_goa$Lon[temp_idx]
           
           temp_idx <- match(x = round(x$Upper_X1), 
                             table = round(grid_goa$E_km))
           x$Upper_X1 <- grid_goa$Lon[temp_idx]
           
           x[, c("Lower_X1", "Upper_X1", "Lower_X2", "Upper_X2")] <- 
             round(x[, c("Lower_X1", "Upper_X1", "Lower_X2", "Upper_X2")], 1)
           
           return(x)
         })

##################################################
## Function to plot the strata boundary plot
##################################################
plot_strata_boundary <- function(
  strata = strata_list[[idx]],
  sol = res_df[, idx],
  extrapolation_grid = grid_goa,
  seed_value = seed,
  y_axis_pos = 2) {
  
  ## Calculate number of strata
  n_strata <- nrow(strata)
  
  ## Create strata colors 
  strata_cols <- colorRampPalette(
    c(brewer.pal(n = 11, name = "Paired"),
      brewer.pal(n = 11, name = "Spectral")))(n_strata)
  
  set.seed(seed_value)
  strata_cols <- sample(strata_cols, replace = F)
  
  ## Calculate ECDF of depth (used to scale y-axis)
  ecdf_depth <- ecdf(extrapolation_grid$DEPTH_EFH)
  
  ## Base plot
  plot(1, 
       type = "n", 
       xlim = range(extrapolation_grid$Lon), 
       ylim = c(1, 0),
       axes = F, ann = F,
       xaxs = "i", yaxs = "i")
  box()
  
  ## y-axis for depth
  axis(side = y_axis_pos,
       at = ecdf_depth(c(0, 50, 100, 150, 200, 300, 1000)),
       labels =        c(0, 50, 100, 150, 200, 300, 1000),
       las = 1)
  
  ## Loop over strata
  for (istrata in 1:n_strata) {
    
    ## Calculate a version of the stratum with more transparency
    temp_col <- t_col(strata_cols[istrata], percent = 50)
    
    ## add extrapolation grid locations
    with(extrapolation_grid[sol == istrata, ],
         points(y = ecdf_depth(DEPTH_EFH), 
                x = Lon,
                pch = 16,
                col = temp_col,
                cex = 0.5)
    )
  }
  
  ## Add labels for stratum labels
  text(x = tapply(extrapolation_grid$Lon, 
                  res_df[, idx],
                  FUN = median),
       y = ecdf_depth(tapply(extrapolation_grid$DEPTH_EFH, 
                             sol,
                             FUN = median)),
       labels = 1:n_strata,
       font = 2,
       cex = 0.75,
       col = "gray28")
}

##################################################
## Plot
##################################################
{
  
  ## Open device
  png(filename = paste0(output_dir, "Figure04_strata_boundaries.png"),
      height = 190, width = 170, units = "mm",
      res = 500)
  
  ## Plot layout
  layout(mat = matrix(data = c(1, 2, 5, 6, 3, 4, 7, 8), ncol = 2),
         heights = c(1, 0.75))
  par(oma = c(2.5, 1, 0, 0))
  
  ## Loop over district level (3 and 5 strata per district) and the gulf-wide
  ## (10 and 15 strata total) solutions. idx refers to the index of the solution
  ## in the settings dataframe
  
  for (idx in c(2, 5, 8, 11)) { ## Loop over scenarios -- start
    
    ## Plot solution map
    par(mar = c(0.5, 3, 0, 0))
    
    ## Set up spatial object
    goa <- sp::SpatialPointsDataFrame(
      coords = grid_goa[, c("E_km", "N_km")],
      data = data.frame(Str_no = res_df[, idx]) )
    goa_ras <- raster::raster(x = goa, 
                              resolution = 10)
    goa_ras <- raster::rasterize(x = goa, 
                                 y = goa_ras, 
                                 field = "Str_no")
    
    ## set up colors
    strata_cols <- colorRampPalette(
      c(brewer.pal(n = 11, name = "Paired"),
        brewer.pal(n = 11, name = "Spectral")))(nrow(strata_list[[idx]]))
    set.seed(seed)
    strata_cols <- sample(strata_cols, replace = F)
    
    plot(x = 1, y = 1, 
         type = "n", axes = F, ann = F, asp = 1,
         xlim = goa_ras@extent[1:2],
         ylim = goa_ras@extent[3:4] + c(0, 0.7 * y_range ))
    
    ## Plot image
    raster::image(x = goa_ras, add = TRUE,
                  axes = F, ann = F,
                  asp = 1,
                  col = strata_cols)
    
    
    ## District labels
    segments(x0 = districts$W_UTM,
             x1 = districts$E_UTM,
             y0 = tapply(X = grid_goa$N_km,
                         INDEX = district_vals,
                         FUN = min) - 100,
             xpd = NA,
             lwd = 2)
    
    ## Simulate a sample solution
    temp_samples <- c()
    temp_strata <- nrow(strata_list[[idx]])
    temp_solution <- res_df[, idx]
    temp_allocation <- strata_list[[idx]]$Allocation
    
    for (temp_istrata in 1:temp_strata) {
      temp_samples = c(temp_samples,
                       sample(x = which(temp_solution == temp_istrata),
                              size = temp_allocation[temp_istrata]) )
    }
    
    ## Plot stations
    temp_loc <- grid_goa[temp_samples, c("E_km", "N_km")]
    temp_loc$N_km <- temp_loc$N_km
    points(x = temp_loc,
           pch = 16,
           cex = 0.3)
    
    ## District labels
    text(x = rowMeans(cbind(districts$W_UTM, districts$E_UTM)),
         y = tapply(X = grid_goa$N_km,
                    INDEX = district_vals,
                    FUN = min) - 150,
         cex = 0.90,
         labels = districts$district,
         pos = 3)
    
    ## Plot solution
    raster::image(raster::shift(x = goa_ras,
                                dy = 0.65 * y_range),
                  asp = 1,
                  col = strata_cols,
                  axes = F,
                  ann = F,
                  add = TRUE)
    
    ## Stratum labels
    text(x = tapply(grid_goa$E_km, 
                    res_df[, idx],
                    median),
         y = tapply(grid_goa$N_km, 
                    res_df[, idx],
                    median) + 0.65 * y_range,
         labels = 1:nrow(strata_list[[idx]]),
         font = 2,
         cex = 0.5)
    
    ## Solution label
    text(x = min(grid_goa$E_km) + x_range*0.15,
         y = min(grid_goa$N_km) + y_range*1.5,
         font = 2, cex = 1.25, xpd = NA,
         labels = paste0(c("2" = "A) ", "5" = "B) ", 
                           "8" = "C) ", "11" = "D) ")[paste0(idx)],
                         ifelse(test = settings$domain[idx] == "district",
                                yes = "District-Level\n",
                                no = "Gulf-Wide\n"),
                         ifelse(test = settings$domain[idx] == "district",
                                yes = paste0(settings$strata[idx],
                                             " Strata\nper District"),
                                no = paste0(settings$strata[idx],
                                            " Strata"))))
    
    ## Solution boundary map
    par(mar = c(2, 4.24, 0, 1))
    plot_strata_boundary(strata = strata_list[[idx]],
                         sol = res_df[, idx],
                         extrapolation_grid = grid_goa,
                         y_axis_pos = 2)
    axis(side = 1)
  }  ## Loop over scenarios -- end
  
  ## Outer axes labels
  mtext(side = 1, text = "Longitude", outer = T, line = 1.25, font = 2)
  mtext(side = 2, text = "Bottom Depth (m)", outer = T, line = -0.25, font = 2)
  
  ## Close device
  dev.off()
}
