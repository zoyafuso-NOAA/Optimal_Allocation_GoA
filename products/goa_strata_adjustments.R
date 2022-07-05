##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       GOA Restratification Adjustments
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Workflow to take the new stratum boundaries and create new
##                strata polygons and intersect its boundaries with the 5 km
##                GOA grid. Then within each 5km grid cell, calculate the 
##                total area and perimeter of the stratum component. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Packages ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
library(rgdal)
library(raster)
library(sp)
library(rgeos)
library(stars)
library(SpaDES)
library(RColorBrewer)
library(spatialEco)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load depth modifications ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
depth_mods <- read.csv("products/depth_modifications.csv")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Merge bathy rasters ----
##   These rasters are really dense, so they are stored in parts and then 
##   "puzzled" togethered using the SpaDES.tools::mergeRaster() function.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
split_bathy <- list()
n_split_rasters <- length(dir("data/raw_data/split_goa_bathy_ras/")) / 2
for (i in 1:n_split_rasters) { ## Loop over subraster -- start
  ## Import raster
  temp_raster <- raster::raster(paste0("data/raw_data/",
                                       "split_goa_bathy_ras",
                                       "/goa_bathy_processed_", 
                                       i, ".grd"))
  
  ## Append raster to the split_bathy list
  split_bathy[[i]] <- temp_raster
} ## Loop over subraster -- end

## Merge rasters together
bathy <- SpaDES.tools::mergeRaster(split_bathy)
rm(split_bathy, i, n_split_rasters, temp_raster)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Shapefiles ---- 
##   goa_domain is a mask of the survey footprint (in lat-lon here)
##   ak_land is AK land
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_domain <- rgdal::readOGR("data/shapefiles/GOAdissolved.shp")
ak_land <- rgdal::readOGR(dsn = "data/shapefiles/alaska_dcw.shp")
goa_grid <- rgdal::readOGR("data/shapefiles/goagrid.shp")
goa_grid_untrawl <- rgdal::readOGR("data/shapefiles/goagrid2019_landuntrawlsndmn.shp")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Management areas ----
##   Create masks for each management area in aea projection and transform
##   goa_domain to aea projection. 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
Shumagin_shape <- raster::crop(x = goa_domain, 
                               y = extent(c(-176, -159, 50, 65)))
Shumagin_shape <- sp::spTransform(x = Shumagin_shape, CRSobj = crs(ak_land))

Chirikof_shape <- raster::crop(x = goa_domain, 
                               y = extent(c(-159, -154, 50, 65)))
Chirikof_shape <- sp::spTransform(x = Chirikof_shape, CRSobj = crs(ak_land))

Kodiak_shape <- raster::crop(x = goa_domain, 
                             y = extent(c(-154, -147, 50, 65)))
Kodiak_shape <- sp::spTransform(x = Kodiak_shape, CRSobj = crs(ak_land))

Yakutat_shape <- raster::crop(x = goa_domain, 
                              y = extent(c(-147, -140, 50, 65)))
Yakutat_shape <- sp::spTransform(x = Yakutat_shape, CRSobj = crs(ak_land))

Southeast_shape <- raster::crop(x = goa_domain, 
                                y = extent(c(-140, -132, 50, 65)))
Southeast_shape <- sp::spTransform(x = Southeast_shape, CRSobj = crs(ak_land))

goa_domain <- sp::spTransform(x = goa_domain, CRSobj = crs(ak_land))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create strata polygons ----
##   For each management area, create new strata based on depth specifications
##   and append to strata_list
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
strata_list <- list()
for (idistrict in unique(depth_mods$manage_area)) {
  
  ## Crop bathymetry raster to just the management area
  district_outline <- get(paste0(idistrict, "_shape"))
  district_bathy <- raster::crop(x = bathy, 
                                 y = district_outline)
  
  ## Define modified stratum depth boundaries
  depth_boundary <- subset(x = depth_mods, 
                           subset = manage_area == idistrict, 
                           select = c("lower_depth_m", "upper_depth_m"))
  
  
  ## Define each raster cell based on the stratum depth boundaries
  values(district_bathy) <- 
    as.integer(as.character(cut(x = values(district_bathy), 
                                breaks = c(0, depth_boundary$upper_depth_m), 
                                labels = 1:nrow(depth_boundary)) ))
  
  ## crop out any part of district_bathy that is outside the survey footprint
  district_bathy <- raster::mask(x = district_bathy, mask = goa_domain)
  
  ## Convert raster to polygon
  strata_poly <- stars::st_as_stars(district_bathy) %>% 
    sf::st_as_sf(merge = TRUE) ## this is the raster to polygons part
  
  ## Convert back to spatial object
  strata_poly <- sf::as_Spatial(strata_poly)
  
  ## Create dataframe of stratum information
  strata_poly@data <- 
    data.frame(manage_area = idistrict, 
               stratum = paste0(idistrict, "_", strata_poly@data$dblbnd))
  
  ## Append to strata_list
  strata_list <- c(strata_list, list(strata_poly))
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Merge strata into one object
##   Calculate the Area and perimeter of each stratum shape
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
strata_list <- raster::bind(strata_list)
strata_list$AREA_KM2 <- raster::area(strata_list) / 1e6
strata_list$PER_KM <- spatialEco::polyPerimeter(x = strata_list) / 1e3

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Take the current goa grid and merge all stratum polygons within a 5x5 km
##   cell. This essentially resets the stratum data in each grid cell.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
goa_full_grid <- raster::aggregate(
  x = goa_grid,
  by = "ID",
  sums = list(list(sum, "AREA_KM2")))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   New stations ----
##   A "station" is defined as the intersection of the grid cells and the
##   stratum polygons. A grid cell can consist of more than one strata. 
##   
##   Intersect the new strata polygons with the 5x5 km grid. This function
##   takes a few minutes. Then calculate the area and perimeter of each station
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stations <- raster::intersect(goa_full_grid, strata_list)
stations <- raster::aggregate(x = stations,
                              by = c("ID", "manage_area", "stratum"))
stations$PER_KM <- spatialEco::polyPerimeter(x = stations) / 1e3
stations$AREA_KM2 <- rgeos::gArea(spgeom = stations, byid = T) / 1e6

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Assign untrawlable stations based on data up to 2019
##   Intersect the untrawlable areas with the new stations and calculate
##   the area of the untrawlable area in the station
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
overlap_with_trawl_polygon <- raster::intersect(stations, goa_grid_untrawl)
overlap_with_trawl_polygon@data <- overlap_with_trawl_polygon@data[, 1:3]
names(overlap_with_trawl_polygon) <- c("ID", "manage_area", "stratum")
overlap_with_trawl_polygon <- 
  raster::aggregate(x = overlap_with_trawl_polygon,
                    by = c("ID", "manage_area", "stratum"))

overlap_with_trawl_polygon$AREA_KM2 <- 
  raster::area(overlap_with_trawl_polygon) / 1e6

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
stations$trawlable <- TRUE
stations$untrawl_area_km2 <- 0

for (irow in 1:length(overlap_with_trawl_polygon)) {
  temp_id <- overlap_with_trawl_polygon$ID[irow]
  temp_str <- overlap_with_trawl_polygon$stratum[irow]
  
  temp_idx <- which(stations$ID == temp_id & stations$stratum == temp_str)
  temp_untrawl <- subset(x = overlap_with_trawl_polygon, 
                         subset = ID == temp_id & stratum == temp_str)
  
  stations$trawlable[temp_idx] <- FALSE
  stations$untrawl_area_km2[temp_idx] <- temp_untrawl$AREA_KM2
  
  rm(temp_idx, irow)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create export directory and export
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
if(!dir.exists(paste0(getwd(), "/products/updated_goa_strata")) ) {
  dir.create(paste0(getwd(), "/products/updated_goa_strata"))
}

names(strata_list) <- c("MGT_AREA", "STRATUM", "AREA_KM2", "PER_KM")
writeOGR(obj = strata_list, 
         dsn = paste0(getwd(), "/products/",
                      "updated_goa_strata/updated_goa_strata.shp"), 
         layer = "updated_goa_strata", 
         driver = "ESRI Shapefile")

names(stations) <- c("ID", "MGT_AREA", "STRATUM", "PER_KM" ,  
                     "AREA_KM2", "TRAWL", "UT_AR_KM2")

writeOGR(obj = stations,
         dsn = paste0(getwd(), "/products/",
                      "updated_goa_strata/updated_stations.shp"), 
         layer = "updated_stations", 
         driver = "ESRI Shapefile")

names(overlap_with_trawl_polygon) <- c("ID", "MGT_AREA", "STRATUM", "AREA_KM2")

writeOGR(obj = overlap_with_trawl_polygon,
         dsn = paste0(getwd(), "/products/",
                      "updated_goa_strata/UT_areas.shp"), 
         layer = "overlap_with_trawl_polygon", 
         driver = "ESRI Shapefile")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plots
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

pdf(file = "products/updated_goa_strata/updated_strata.pdf", 
    width = 6, height = 6, onefile = TRUE)
for (iarea in unique(depth_mods$manage_area)) {
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  n_strata <- with(depth_mods, table(manage_area))[iarea] 
  temp_area <- subset(stations, MGT_AREA == iarea)
  temp_area <- temp_area[order(temp_area$STRATUM), ]
  
  plot(temp_area, axes = F, col = "white", border = F)
  
  for (istratum in 1:n_strata) {
    plot(subset(stations, 
                STRATUM == unique(temp_area$STRATUM)[istratum]), 
         col =   c(RColorBrewer::brewer.pal("Spectral", 
                                            n = n_strata - 1), 
                   "gray" )[istratum], 
         border = F, add = TRUE)
  }
  
  plot(overlap_with_trawl_polygon, col = "black", add = TRUE, border = FALSE)
  plot(ak_land, add = TRUE, col = "tan", border = F)
}
dev.off()
