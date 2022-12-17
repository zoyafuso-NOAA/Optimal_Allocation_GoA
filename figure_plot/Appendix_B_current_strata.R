###############################################################################
## Project:         Overlay every current GoA stratum on optimized destigns
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
###############################################################################
rm(list = ls())

##################################################
####  Set up directories  
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[2]

output_dir <- paste0(c("/Users/zackoyafuso/Google Drive/",
                       "C:/Users/Zack Oyafuso/Google Drive/")[which_machine],
                     "MS_Optimizations/TechMemo/appendix/Appendix B plots/")

if (!dir.exists(output_dir)) dir.create(output_dir)

##################################################
####  Libraries
##################################################
library(rgdal)
library(rgeos)
library(RColorBrewer)

##################################################
####  Import Data
####  AK_land: alaska land polygons
####  goa_strata: Current GoA strata polygons
####  MS_optimization_knitted_results: optimal solutions
####  grid_goa.RData: goa_grid that corresponds to the solutions
####  optimization_data.RData: constants 
##################################################
AK_land <- rgdal::readOGR("data/shapefiles/AKland.shp")
AK_land <- sp::spTransform(x = AK_land,
                           CRSobj = "+proj=utm +units=km +zone=5")

goa_strata <- rgdal::readOGR("data/shapefiles/goa_strata.shp")
goa_strata <- sp::spTransform(x = goa_strata,
                              CRSobj = "+proj=utm +units=km +zone=5")

load("data/processed/grid_goa.RData")
load("data/processed/optimization_data.RData")

n_districts <- 5
n_depth_zones <- 4

##################################################
####  Current GoA Strata names
##################################################
strata_groupings <- list(
  "West" = 
    list("Depth Zone 1-100 m" = c("Fox Islands" = 10,
                                  "Davidson Bank" = 11,
                                  "Lower Alaska Peninsula" = 12,
                                  "Shumagin Bank" = 13),
         "Depth Zone 101-200 m" = c("Sanak Gully" = 110,
                                    "Shumagin Outer Shelf" = 111,
                                    "W. Shumagin Gully" = 112),
         "Depth Zone 201-300 m" = c("Shumagin Slope" = 210),
         "Depth Zone 301-1000 m" = c("Shumagin Slope" = 310,
                                     "Shumagin Slope" = 410,
                                     "Shumagin Slope" = 510)),
  
  "Chirikof" = 
    list("Depth Zone 1-100m" = c("Upper AK Peninsula" = 20,
                                 "Semidi Bank" = 21,
                                 "Chirikof Bank" = 22),
         "Depth Zone 101-200 m" = c("E Shumagin Gully" = 120,
                                    "Shelikof Edge" = 121,
                                    "Chirikof Outer Shelf" = 122),
         "Depth Zone 201-300 m" = c("Lower Chirikof Gully" = 220,
                                    "Chirikof Slope" = 221,
                                    "Upper Shelikof Gully" = 232),
         "Depth Zone 301-1000 m" = c("Chirikof Slope" = 320,
                                     "Chirikof Slope" = 420,
                                     "Chirikof Slope" = 520)),
  
  "Kodiak" = 
    list("Depth Zone 1-100 m" = c("Albatross Shallows" = 30,
                                  "Albatross Banks" = 31,
                                  "Lower Cook Inlet" = 32,
                                  "Kenai Peninsula" = 33,
                                  "N. Kodiak Shallows" = 35),
         "Depth Zone 101-200 m" = c("Albatross Gullies" = 130,
                                    "Portlock Flats" = 131,
                                    "Barren Islands" = 132,
                                    "Kenai Flats" = 133,
                                    "Kodiak Outer Shelf" = 134),
         "Depth Zone 201-300 m" = c("Kenai Gullies" = 230,
                                    "Kodiak Slope" = 231),
         "Depth Zone 301-1000 m" = c("Kodiak Slope" = 330,
                                     "Kodiak Slope" = 430,
                                     "Kodiak Slope" = 530)),
  
  "Yakutat" = 
    list("Depth Zone 1-100 m" = c("Yakutat Shallows" = 40,
                                  "Middleton Shallows" = 41),
         "Depth Zone 101-200 m" = c("Middleton Shelf" = 140, 
                                    "Yakataga Shelf" = 141,
                                    "Yakutat Flats" = 142),
         "Depth Zone 201-300 m" = c("Yakutat Gullies" = 240, 
                                    "Yakutat Slope" = 241),
         "Depth Zone 301-1000 m" = c("Yakutat Gullies" = 340, 
                                     "Yakutat Slope" = 341,
                                     "Yakutat Slope" = 440, 
                                     "Yakutat Slope" = 540)),
  "Southeast" = 
    list("Depth Zone 1-100 m" = c("SE Shallows" = 50),
         "Depth Zone 101-200 m" = c("Fairweather Shelf" = 143,
                                    "Baranof-Chichagof Shelf" = 150, 
                                    "Prince of Wales Shelf" = 151),
         "Depth Zone 201-300 m" = c("Baranof-Chichagof Slope" = 250, 
                                    "Prince of Wales Slope/Gullies" = 251),
         "Depth Zone 301-1000 m" = c("SE Deep Gullies" = 350, 
                                     "SE Slope" = 351,
                                     "SE Slope" = 450, 
                                     "SE Slope" = 550)) )

##################################################
####  Plot settings detials
##################################################
plot_settings <- data.frame(district = districts$district,
                            x_max_offset = c(0, 3.55, 3.6, 0, 3.5),
                            x_offset = c(0, 400, 500, 0, 700),
                            y_min_offset = c(3.5, 0.15, 0.3, 5.5, 0.2),
                            y_offset = c(325, 0, 0, 250, 0),
                            x_legend = c(0, -0.1, 0, 0.15, 0.1),
                            y_legend = c(1.05, 1.25, 1.5, 1.75, 1.35),
                            legend_cex = c(0.44, 0.5, 0.5, 0.45, 0.44),
                            y_inter = c(0.9, 1, 1, 1, 1))

##################################################
####  Plot
##################################################
for (iscen in scenarios$scen_name[-1]) {
  
  set.seed(23433)
  
  ## Set up plot
  png(filename = paste0(output_dir, "Appendix B scenario ", iscen, ".png"),
      units = "in", width = 6, height = 6.5, res = 500, family = "serif")
  
  ## Set up plot layout
  par(mar = c(0, 0, 0, 0), 
      oma = c(0.5, 0.5, 0.5, 0.5))
  layout(matrix(data = c(1, 1, 1, 1, 1, 4, 4, 4, 4, 4,
                         6, 6, 2, 2, 2, 3, 3, 3, 5, 5), ncol = 2),
         widths = c(1, 2))
  
  load(paste0("results/scenario_", iscen, "/Multispecies_Optimization/",
              "Str_5/boat_2/result_list.RData"))
  
  for(idistrict_idx in 1:n_districts) { ## For each district--start loop
    
    ## Subset those cells that are in the given district 
    sub_res_df <- result_list$sol_by_cell[district_vals == paste(idistrict_idx)]
    
    ## Set up spatial object
    goa <- sp::SpatialPointsDataFrame(
      coords = grid_goa[district_vals == idistrict_idx, 
                        c("E_km", "N_km")],
      data = data.frame(Str_no = sub_res_df) )
    goa_ras <- raster::raster(x = goa, 
                              resolution = 3.8)
    goa_ras <- raster::rasterize(x = goa, 
                                 y = goa_ras, 
                                 field = "Str_no")
    
    x_extent = diff(goa_ras@extent[1:2])
    y_extent = diff(goa_ras@extent[3:4])
    district_name <- districts$district[idistrict_idx]
    
    ## Base plot and subtitle
    plot(1, 
         type = "n", 
         xlim = as.numeric(districts[idistrict_idx, 
                                     c("W_UTM", "E_UTM")]) +
           c(0, x_extent * plot_settings$x_max_offset[idistrict_idx]),
         ylim = goa_ras@extent[3:4] +
           c(0, y_extent * plot_settings$y_min_offset[idistrict_idx]), 
         axes = F, ann = F, asp = 1)
    box()
    mtext(side = 1, 
          text = paste0(LETTERS[idistrict_idx], ") ", district_name), 
          line = -1.25, font = 2)
    
    ## Plot current strata by depth zone
    for(istratas in 1:n_depth_zones) { ## For each depth zone--start loop
      
      set.seed(idistrict_idx)
      
      ## Plot land
      plot( raster::shift(
        x = raster::intersect(x = AK_land,
                              sp::bbox(goa_ras)),
        dx = plot_settings$x_offset[idistrict_idx] * (istratas - 1),
        dy = plot_settings$y_offset[idistrict_idx] * (istratas - 1)),
        add = TRUE, col = "tan", border = F)
      
      depth_zone <- attributes(strata_groupings[[district_name]][istratas])$names
      temp_strata <- strata_groupings[[district_name]][[istratas]]
      
      raster::image( 
        raster::shift( x = goa_ras, 
                       dx = plot_settings$x_offset[idistrict_idx] * (istratas - 1),
                       dy = plot_settings$y_offset[idistrict_idx] * (istratas - 1)), 
        col = sample(brewer.pal(n = 9, name = "Set3"), size = 5),
        asp = 1,
        add = TRUE)
      
      ## Subset current strata within the district and dissolve if there are 
      ## multiple records to a particular stratum
      temp_goa_strata <- subset(goa_strata, STRATUM %in% temp_strata)
      temp_goa_strata <- temp_goa_strata[order(temp_goa_strata$STRATUM), ]
      temp_goa_strata <- rgeos::gUnaryUnion(spgeom = temp_goa_strata, 
                                            id = temp_goa_strata@data$STRATUM)
      
      ## Since we are plotting multple rasters on one plot, offset the raster
      temp_goa_strata <- 
        raster::shift(x = temp_goa_strata, 
                      dx = plot_settings$x_offset[idistrict_idx] * (istratas - 1),
                      dy = plot_settings$y_offset[idistrict_idx] * (istratas - 1) )
      
      ## Plot current strata on top of solution
      plot(temp_goa_strata, add = TRUE, lwd = 0.5,
           border = c("black", "purple", 
                      "blue", "red", "green")[1:length(temp_strata)])
      
      
      
      legend(x = districts$W_UTM[idistrict_idx] + 
               x_extent * plot_settings$x_legend[idistrict_idx] + 
               plot_settings$x_offset[idistrict_idx] * (istratas - 1),
             y = goa_ras@extent[3] + 
               y_extent * plot_settings$y_legend[idistrict_idx] + 
               plot_settings$y_offset[idistrict_idx] * (istratas - 1) ,
             legend = paste0("Str ", temp_strata, ": ", names(temp_strata)),
             lty = 1,
             col = c("black", "purple", 
                     "blue", "red", "green")[1:length(temp_strata)],
             title = depth_zone,
             cex = plot_settings$legend_cex[idistrict_idx],
             y.intersp = plot_settings$y_inter[idistrict_idx], 
             bg = "white")
      
    } ## For each depth zone--end loop
  }  ## For each district--end loop
  
  ## Districts Legend: first plot spatial domain
  goa <- sp::SpatialPointsDataFrame(
    coords = grid_goa[, c("E_km", "N_km")],
    data = data.frame(Str_no = 1:nrow(grid_goa))) 
  
  goa_ras <- raster::raster(x = goa, 
                            resolution = 10)
  goa_ras <- raster::rasterize(x = goa, 
                               y = goa_ras, 
                               field = "Str_no")
  
  ## Plot
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  raster::image(goa_ras,
                col = "grey",
                asp = 1,
                axes = F,
                ann = F)
  
  ## Plot land
  plot( AK_land,
        add = TRUE, col = "tan", border = F)
  
  ## Plot and label district boxes
  rect(xleft = districts$W_UTM,
       xright = districts$E_UTM,
       ybottom = tapply(X = grid_goa$N_km,
                        INDEX = district_vals,
                        FUN = min),
       ytop = tapply(X = grid_goa$N_km,
                     INDEX = district_vals,
                     FUN = max) )
  
  text(x = tapply(X = grid_goa$E_km, 
                  INDEX = district_vals,
                  FUN = mean),
       y = tapply(X = grid_goa$N_km, 
                  INDEX = district_vals,
                  FUN = mean),
       labels = LETTERS[1:5],
       font = 2)
  
  box(which = "figure")
  
  dev.off()
  
}
