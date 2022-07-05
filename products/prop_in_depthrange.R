###############################################################################
## Project:       How well are stations witin their existing depth stratum
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Calculate how often raster cells in a stratum fall within the
##                depth range of the stratum
##                Current depth strata are (in meters):
##                 0-100, 101-200, 201-300, 301-500, 501-700, and 701-1000
###############################################################################
rm(list = ls())

##################################################
####   Import packages
##################################################
library(rgdal)

###############################################################################
####   Import strata shapefiles
####   Import haul information
###############################################################################
goa_strata <- rgdal::readOGR(dsn = "data/shapefiles/goa_strata.shp")
goa_strata <- subset(goa_strata, STRATUM > 0)
strata <- sort(unique(goa_strata@data$STRATUM))
strata[strata < 100] <- paste0("0", strata[strata < 100])
strata <- strata[order(substr(x = strata, start = 2, stop = 2))]
management_area <- c("Shumagin", "Chirikof", "Kodiak", 
                     "Yakutat", "SE")[as.integer(substr(x = strata, 
                                                        start = 2, 
                                                        stop = 2))]

haul <- read.csv("data/raw_data/haul.csv")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Empty dataframe
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
prop_in_depth <- data.frame()

for (iarea in c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "SE")) {
  
  ## Loop over the strata in a given iarea
  for (istratum in as.integer(strata[management_area == iarea])) {
    
    ## Extract depth range of istratum given the code of the istratum
    min_depth <- c(0, 101, 201, 301, 501, 701)[ceiling(istratum / 100)]
    max_depth <- c(100, 200, 300, 500, 700, 1000)[ceiling(istratum / 100)]
    
    ## Subset those hauls that occurred in istratum 
    temp_hauls <- subset(haul, subset = STRATUM == istratum)
    temp_hauls <- temp_hauls[!is.na(temp_hauls$BOTTOM_DEPTH), ]
    
    ## Calculate the proportion of the hauls that occurred within, shallower, 
    ## and deeper than the specified depth range
    prop_correct <- sum(with(temp_hauls, 
                             BOTTOM_DEPTH <= max_depth & 
                               BOTTOM_DEPTH >= min_depth)) / nrow(temp_hauls)
    prop_shallower <- sum(with(temp_hauls, 
                               BOTTOM_DEPTH < min_depth)) / nrow(temp_hauls)
    prop_deeper <- sum(with(temp_hauls, 
                            BOTTOM_DEPTH > max_depth)) / nrow(temp_hauls)
    
    ## Append to result dataframe
    prop_in_depth <- 
      rbind(prop_in_depth,
            data.frame(management_area = iarea, 
                       stratum = istratum,
                       min_depth_m = min_depth,
                       max_depth_m = max_depth,
                       prop_shallower = round(prop_shallower, 2),
                       prop_within_depth = round(prop_correct, 2),
                       prop_deeper = round(prop_deeper, 2)))
  }
}


###############################################################################
####   Save
###############################################################################  
write.csv(x = prop_in_depth, 
          file = "products/prop_in_depthrange.csv",
          row.names = F)
