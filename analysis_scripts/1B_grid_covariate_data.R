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
library(rgdal)
library(raster)
library(SpaDES)

##################################################
####   CRSs used
##################################################
lonlat_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
utm_crs <- "+proj=utm +zone=5N +units=km"

##################################################
####   Merge together bathymetry rasters
##################################################
split_bathy <- list()
n_split_rasters <- length(dir("data/raw_data/split_goa_bathy_ras/")) / 2
for (i in 1:n_split_rasters) {
  split_bathy[[i]] <- raster::raster(paste0("data/raw_data/",
                                            "split_goa_bathy_ras",
                                            "/goa_bathy_processed_", 
                                            i, ".grd"))
}

bathy <- SpaDES.tools::mergeRaster(split_bathy)
rm(split_bathy)

##################################################
####   Import Current Strata (current_survey_strata)
####   Import spatial domain outline mask (current_survey_mask)
####   Import Extrapolation grid (goa_grid)
####   Import CPUE data (data)
####   Untrawlable areas (goa_grid_nountrawl)
##################################################
current_survey_strata <- rgdal::readOGR("data/shapefiles/goa_strata.shp")
current_survey_mask <- rgdal::readOGR("data/shapefiles/goagrid_polygon.shp")
current_survey_mask <- sp::spTransform(x = current_survey_mask,
                                       CRSobj = crs(bathy))
current_survey_mask <- rgeos::gUnaryUnion(spgeom = current_survey_mask)

goa_grid <- read.csv("data/extrapolation_grid/GOAThorsonGrid.csv")

goa_grid <- goa_grid[, c("Id", "Shape_Area", "Longitude", "Latitude")]
goa_grid$Shape_Area <- goa_grid$Shape_Area / 1000 / 1000 #Convert to km2 
names(goa_grid) <- c( "Id", "Area_km2", "Lon", "Lat")

data = read.csv("data/processed/goa_vast_data_input.csv")

goa_grid_nountrawl <- read.csv("data/extrapolation_grid/GOA_ALL_nountrawl.csv")

##################################################
####   Transform extrapolation grid to aea, extract bathymetry values onto grid
##################################################
grid_shape = sp::SpatialPointsDataFrame(
  coords = goa_grid[, c("Lon", "Lat")],
  data = goa_grid,
  proj4string = CRS(lonlat_crs))

grid_shape_aea = sp::spTransform(x = grid_shape,
                                 CRSobj = crs(bathy))
grid_shape_aea@data$DEPTH_EFH =  raster::extract(x = bathy,
                                                 y = grid_shape_aea,
                                                 method = "simple")

##################################################
####   Remove cells not already in the goa stratification
##################################################
grid_shape_aea <- raster::intersect(x = grid_shape_aea,
                                    y = current_survey_mask)

##################################################
####   Remove cells with depths outside the observed range to the range
##################################################
grid_shape_aea <- subset(grid_shape_aea,
                         DEPTH_EFH >= min(data$DEPTH_EFH) & 
                           DEPTH_EFH <= max(data$DEPTH_EFH))

##################################################
####   Plot depth covariate of the extrapolation grid
##################################################
# spplot(grid_shape_aea[, "DEPTH_EFH"], 
#        col.regions = rev(terrain.colors(1000)),
#        pch = 16, cex = 0.1,
#        cuts = 100, colorkey = T)

grid_goa <- grid_shape_aea@data
grid_goa[, c("E_km", "N_km")] <- 
  project(xy = coordinates(grid_goa[, c("Lon", "Lat")]), 
          proj = utm_crs )

##################################################
####   Create indices to easily subset <700 m cells and untrawlable cells
##################################################
grid_goa$shallower_than_700m <- cells_shallower_than_700m <- 
  grid_goa$DEPTH_EFH <= 700
grid_goa$trawlable <- cells_trawlable <- grid_goa$Id %in% goa_grid_nountrawl$Id
grid_goa$shallow_trawlable <- 
  rowSums(grid_goa[, c("shallower_than_700m", "trawlable")]) == 2

##################################################
####   scale grid bathymetry values to standard normal, using the mean and sd
####   of the BTS data
##################################################
BTS_mean <- mean(log(data$DEPTH_EFH))
BTS_sd   <- sd(log(data$DEPTH_EFH))

grid_goa$LOG_DEPTH_EFH <- log(grid_goa$DEPTH_EFH)
grid_goa$LOG_DEPTH_EFH_CEN <- (grid_goa$LOG_DEPTH_EFH - BTS_mean) / BTS_sd
grid_goa$LOG_DEPTH_EFH_CEN_SQ <- grid_goa$LOG_DEPTH_EFH_CEN ^ 2

##################################################
####   Add current strata labels to each grid cell
##################################################
grid_goa$stratum <- raster::extract(x = current_survey_strata, 
                                    y = grid_shape_aea)$STRATUM
grid_goa$stratum[is.na(grid_goa$stratum)] <- 0

##################################################
####   Save
##################################################
if(!dir.exists("data/processed/")) dir.create("data/processed/")
save(list = "grid_goa", file = paste0("data/processed/grid_goa.RData"))
write.csv(grid_goa, row.names = F, file = "data/processed/grid_goa.csv")

