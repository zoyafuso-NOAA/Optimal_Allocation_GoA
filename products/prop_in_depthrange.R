###############################################################################
## Project:       
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Calculate how often raster cells in a stratum fall within the
##                depth range of the stratum
###############################################################################
rm(list = ls())

##################################################
####   Import packages
##################################################
library(rgdal)
library(raster)
library(SpaDES)
library(RColorBrewer)

##################################################
####   Merge together bathymetry rasters
##################################################
split_bathy <- list()
n_split_rasters <- length(dir("data/raw_data/split_goa_bathy_ras/")) / 2
for (i in 1:n_split_rasters) {
  temp_raster <- raster::raster(paste0("data/raw_data/",
                                       "split_goa_bathy_ras",
                                       "/goa_bathy_processed_", 
                                       i, ".grd"))
  split_bathy[[i]] <- temp_raster
}

bathy <- SpaDES.tools::mergeRaster(split_bathy)
rm(split_bathy)

###############################################################################
####   Import strata shapefiles
###############################################################################
goa_strata <- rgdal::readOGR(dsn = "data/shapefiles/goa_strata.shp")
goa_strata <- subset(goa_strata, STRATUM > 0)
strata <- sort(unique(goa_strata@data$STRATUM))
strata[strata < 100] <- paste0("0", strata[strata < 100])
strata <- strata[order(substr(x = strata, start = 2, stop = 2))]
management_area <- 
  c("Shumagin", "Chirikof", "Kodiak", 
    "Yakutat", "SE")[as.integer(substr(x = strata, 
                                       start = 2, 
                                       stop = 2))]

prop_in_depth <- data.frame()

for (iarea in c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "SE")) {
  par(mfrow = c(4, 4), mar = c(3, 5, 3, 0))
  for (istratum in as.integer(strata[management_area == iarea])) {
    min_depth <- c(0, 101, 201, 301, 501, 701)[ceiling(istratum / 100)]
    max_depth <- c(100, 200, 300, 500, 700, 1000)[ceiling(istratum / 100)]
    
    extracted_vals <- 
      na.omit(unlist(raster::extract(x = bathy, 
                                     y = subset(goa_strata, 
                                                STRATUM == istratum))))
    correct_depths <- sum(extracted_vals >= min_depth & 
                            extracted_vals < max_depth,
                          na.rm = TRUE)
    
    Nh <- length(extracted_vals)
    deeper <- sum(extracted_vals > max_depth, na.rm = TRUE)
    shallower <- sum(extracted_vals < min_depth, na.rm = TRUE)
    
    prop_shallower <- round(shallower / Nh, 2)
    prop_correct <- round(correct_depths /Nh, 2)
    prop_deeper <- round(deeper /Nh, 2)
    
    prop_in_depth <- rbind(prop_in_depth,
                           data.frame(management_area = iarea, 
                                      stratum = istratum,
                                      min_depth_m = min_depth,
                                      max_depth_m = max_depth,
                                      prop_shallower = prop_shallower,
                                      prop_within_depth = prop_correct,
                                      prop_deeper = prop_deeper))
    
    print(prop_in_depth[nrow(prop_in_depth), ])
  }
}

###############################################################################
####   Save
###############################################################################  
save(prop_in_depth, 
     file = "products/depth_distributions_strata/prop_in_depthrange.csv")
