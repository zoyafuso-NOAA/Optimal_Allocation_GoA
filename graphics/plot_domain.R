###############################################################################
## Project:       Plot Spatial Domain
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Plot spatial domains used in analysis
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(raster)
library(sp)
library(rnaturalearth)
library(ocedata)

##################################################
####   Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[2]
VAST_model <- "11" 
which_domain <- c("full_domain", "trawlable")[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/",
                     which_domain, "/")

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load(paste0(github_dir, 'optimization_data.RData'))
load(paste0(dirname(dirname(github_dir)), '/data/Extrapolation_depths.RData'))

##################################################
####   Import Land objects, clip to save object space
##################################################
AK = rgdal::readOGR(paste0(dirname(dirname(github_dir)), 
                           "/data/shapefiles/AKland.shp")) 
AK = sp::spTransform(AK, 
                     CRS = paste0("+proj=utm +zone=5 +lat_1=55 +lat_2=65",
                                  " +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0",
                                  " +ellps=GRS80 +datum=NAD83 +units=km",
                                  " +no_defs"))

CA = rgdal::readOGR(paste0(dirname(dirname(github_dir)), 
                           "/data/shapefiles/canada_dcw.shp")) 
CA = sp::spTransform(CA, 
                     CRS = paste0("+proj=utm +zone=5 +lat_1=55 +lat_2=65",
                                  " +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0",
                                  " +ellps=GRS80 +datum=NAD83 +units=km",
                                  " +no_defs"))
{
png(filename = paste0(dirname(dirname(github_dir)), "/graphics/domain.png"),
    width = 6,
    height = 4,
    units = "in",
    res = 500)
par(mar = c(0,0,0,0), mfrow = c(2,1))
plot(N_km ~ E_km, 
     data = Extrapolation_depths, 
     pch = 15,
     cex = 0.2, 
     asp = 1, 
     axes = F,
     ann = F)
legend("bottomright", 
       legend = "Full domain (black)", 
       bty = "n")
plot(AK, add = T, col = "tan", border = F)
plot(CA, add = T, col = "tan", border = F)
box()

plot(N_km ~ E_km, 
     data = Extrapolation_depths[trawl_shallow_idx, ], 
     pch = 15,
     cex = 0.2, 
     asp = 1, 
     axes = F,
     ann = F)
legend("bottomright", 
       legend = "Untrawlable areas (red) removed", 
       bty = "n")
points(N_km ~ E_km, 
       data = Extrapolation_depths[!trawl_shallow_idx, ], 
       pch = 15,
       cex = 0.2, 
       col = "red")
plot(AK, add = T, col = "tan", border = F)
plot(CA, add = T, col = "tan", border = F)
box()

dev.off()
}
