###############################################################################
## Project:       Optimized Strata Boundaries
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Create excel worksheet 
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(raster)
library(sp)
library(rgdal)
library(RColorBrewer)
library(xlsx); library(openxlsx)

##################################################
####   Import Data
##################################################
ak_land <- rgdal::readOGR("data/shapefiles/AKland.shp")
load('data/processed/optimization_data.RData')
load('data/processed/grid_goa.RData')

##################################################
####   Scenarios to pull
##################################################
scen <- c("D", "E", "J", "K")
scenarios <- subset(x = scenarios, 
                    subset = scen_name %in% scen, 
                    select = -c(scale_opt, data_type, max_depth))
names(scenarios) <- c("scenario", "stratum variables", "Notes")
scenarios$`stratum variables` <- 
  sapply(X = scenarios$`stratum variables`,
         FUN = function(x) switch(x,
                                  "depth_lon" = "depth and longitude",
                                  "depth" = "depth"))

scenarios$`Notes` <- 
  sapply(X = paste0(scenarios$`Notes`),
         FUN = function(x) switch(x,
                                  "300" = "depths > 300 m  are set to 300 m in optimization",
                                  "1000" = ""))

##################################################
####   Result Objects
##################################################
xlsx::write.xlsx(x = scenarios, 
                 file = "products/STRS_Solutions.xlsx",
                 sheetName = "scenario overview", 
                 row.names = FALSE)

strs_sols <- sp::SpatialPointsDataFrame(
  coords = grid_goa[, c("Lon", "Lat")], 
  proj4string = crs("+proj=longlat +datum=WGS84"), 
  data = grid_goa[, c("Lon", "Lat")])
# strs_sols_list <- list()

##################################################
####   Loop over scenarios, load, and pull out solution s
##################################################
for (iscen in scen) {
  load(paste0("results/scenario_", iscen, 
              "/Multispecies_Optimization/Str_5/boat_2/result_list.RData"))
  
  result_list$sol_by_cell[!grid_goa$shallower_than_700m] <- 0 
  
  strs_sols@data[, paste0("scenario_", iscen)] <- result_list$sol_by_cell
  stratum_area <- round(tapply(X = grid_goa$Area_km2, 
                               INDEX = result_list$sol_by_cell, 
                               FUN = sum))
  
  temp <- result_list$sum_stats[, c("Domain", "Allocation",
                                    "Lower_X1", "Upper_X1") ]
  temp$Domain <- districts$district[temp$Domain]
  names(temp)[3:4] <- c("lower_depth", "upper_depth")
  
  # if (iscen %in% c("D", "J")) {
  #   temp <- cbind(temp,
  #                 result_list$sum_stats[, c("Lower_X2", "Upper_X2") ])
  #   names(temp)[5:6] <-  c("lower_lon", "upper_lon")
  # }
  
  temp <- cbind(temp, 
                total_area_km2 = stratum_area[names(stratum_area) != "0"])
  
  xlsx::write.xlsx(x = temp, 
                   file = "products/STRS_Solutions.xlsx",
                   sheetName = paste0("scenario_", iscen), 
                   row.names = FALSE, 
                   append = TRUE)
  
  ## Plot
  png(filename = paste0("products/scenario_", iscen, ".png"),
      res = 500, width = 8, height = 4, units = "in")
  par(mar = c(0, 0, 0, 0))
  no_str <- length(unique(result_list$sol_by_cell))
  strata_cols <- c("black",
                   rep(x = RColorBrewer::brewer.pal(n = 10, name = "Paired"),
                       3))[1:no_str]
  plot(strs_sols,
       col = strata_cols[1 + strs_sols@data[, paste0("scenario_", iscen)]],
       cex = 0.18, pch = 15)
  
  legend("bottom", legend = "**black indicates areas > 700 m", bty = "n")
  dev.off()
  
  wb <- loadWorkbook(file = "products/STRS_Solutions.xlsx")
  img <- system.file("extdata", 
                     paste0("products/scenario_", iscen, ".png"), 
                     package = "openxlsx")
  insertImage(wb = wb, 
              sheet = paste0("scenario_", iscen), 
              file = paste0("products/scenario_", iscen, ".png"), 
              width = 8, height = 4, startCol = 7, startRow = 3)
  saveWorkbook(wb, "products/STRS_Solutions.xlsx", overwrite = TRUE)
  
}

strs_sols <- sp::spTransform(x = strs_sols,
                             CRSobj = crs(ak_land))
