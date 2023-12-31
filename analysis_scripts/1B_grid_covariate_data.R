###############################################################################
## Project:     VAST covariates across goa grid
## Author:      Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description: Assign bathymetry value from the EFH data lyaer to each grid in
##              the Gulf of Alaska grid.
###############################################################################
rm(list = ls())

##################################################
####   Import packages
##################################################
library(terra)

##################################################
####   CRSs used
##################################################
lonlat_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

##################################################
####   Merge together bathymetry rasters
##################################################
goa_bathy <-
  # terra::rast("//AKC0SS-n086/AKC_PubliC/Dropbox/Zimm/GEBCO/GOA/goa_bathy/")
terra::rast("C:/Users/zack.oyafuso/Desktop/goa_bathy/")

##################################################
####   Import Current Strata (current_survey_strata)
####   Import spatial domain outline mask (current_survey_mask)
####   Import Extrapolation grid (goa_grid)
####   Import CPUE data (data)
####   Untrawlable areas (goa_grid_nountrawl)
##################################################
current_survey_strata <- terra::vect(x = "data/shapefiles/goa_strata.shp")
current_survey_strata <- 
  current_survey_strata[current_survey_strata$STRATUM != 0]
current_survey_strata <- 
  terra::project(x = current_survey_strata, terra::crs(goa_bathy))

goa_grid <- read.csv("data/extrapolation_grid/GOAThorsonGrid_Less700m.csv")
goa_grid <- goa_grid[, c("Id", "Shape_Area", "Longitude", "Latitude")]
goa_grid$Shape_Area <- goa_grid$Shape_Area / 1000 / 1000 #Convert to km2 
names(goa_grid) <- c( "Id", "Area_km2", "Lon", "Lat")

goa_data_geostat = read.csv("data/processed/goa_data_geostat.csv")

##################################################
####   Transform extrapolation grid to aea, extract bathymetry values onto grid
##################################################
grid_shape = terra::vect(x = goa_grid[, c("Lon", "Lat", "Area_km2")], 
                         geom = c("Lon", "Lat"), keepgeom = TRUE,
                         crs = lonlat_crs)
grid_shape_aea <- terra::project(x = grid_shape, terra::crs(x = goa_bathy))
grid_shape_aea[, c("Eastings", "Northings")] <- terra::crds(x = grid_shape_aea)

grid_shape_aea$Depth_m <-
  terra::extract(x = goa_bathy, y = grid_shape_aea, ID = FALSE)$goa_bathy

##################################################
####   Remove cells not already in the current GOA strata
##################################################
grid_shape_aea <- terra::intersect(x = grid_shape_aea,
                                    y = current_survey_strata[, "STRATUM"])

##################################################
####   Remove cells with depths outside the observed range to the range
##################################################
grid_shape_aea <- grid_shape_aea[
  (grid_shape_aea$Depth_m >= min(x = goa_data_geostat$Depth_m) & 
     grid_shape_aea$Depth_m <= 700),
]

##################################################
####   scale grid bathymetry values to standard normal, using the mean and sd
####   of the BTS data
##################################################
BTS_mean <- mean(log10(x = goa_data_geostat$Depth_m))
BTS_sd   <- sd(log10(x = goa_data_geostat$Depth_m))

grid_shape_aea$LOG10_DEPTH_M <- log10(x = grid_shape_aea$Depth_m)
grid_shape_aea$LOG10_DEPTH_M_CEN <- 
  (grid_shape_aea$LOG10_DEPTH_M - BTS_mean) / BTS_sd

##################################################
####   Save
##################################################
if(!dir.exists(paths = "data/processed/")) dir.create(path = "data/processed/")
write.csv(as.data.frame(x = grid_shape_aea), 
          row.names = F, 
          file = "data/processed/goa_interpolation_grid.csv")

saveRDS(object = as.data.frame(x = grid_shape_aea), 
        file = "data/processed/goa_interpolation_grid.RDS")
