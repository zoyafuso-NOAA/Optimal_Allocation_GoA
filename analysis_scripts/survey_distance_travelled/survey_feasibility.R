###############################################################################
## Project:       Survey Feasibility (total distance travelled)
## Author:        Megsie Siple (margaret.siple@noaa.gov)
##                Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Simulate several 2-boat surveys station locations
##                Calculate shortest distance path across two boats
##                Calculate total distance travelled
###############################################################################
rm(list = ls())

library(TSP); library(sp); library(tidyverse)

##################################################
####   Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")

source(file = paste0(github_dir, '/modified_functions/sample_stations.R'))
source(file = paste0(github_dir, '/modified_functions/do_TSP.R'))

load(file = paste0(github_dir, '/data/Extrapolation_depths.RData'))
load(file = paste0(github_dir, '/results/MS_optimization_knitted_results.RData'))
load(file = paste0(github_dir, "data/processed_haul_data.RData"))

GOA_allocations <- readxl::read_xlsx(
  path = paste0(github_dir, '/data/GOA 2019 stations by stratum.xlsx'))
GOA_allocations3 <- readxl::read_xlsx(
  path = paste0(github_dir, '/data/GOA2019_ 3 boat_825_RNDM_stations.xlsx')) 

##################################
## Specify Management Districts
##################################
new_strata_labels = 1:length(unique(Extrapolation_depths$stratum))
names(new_strata_labels) <- sort(unique(Extrapolation_depths$stratum))

##################################
## Rename Current Stratum labels
##################################
Extrapolation_depths$stratum_new_label <- 
  new_strata_labels[paste(Extrapolation_depths$stratum)]

##################################################
####   Create dataframe of effort allocations across boats
##################################################
allocations <- data.frame(Stratum = sort(unique(GOA_allocations3$stratum)),
                          boat3 = aggregate(id ~ stratum, 
                                            data = GOA_allocations3, 
                                            FUN = length)$id,
                          boat2 = c(GOA_allocations$`Number stations`,
                                    rep(0, 5)))

allocations <- rbind(data.frame(Stratum = 0, boat3 = 0, boat2 = 0),
                     allocations)

allocations$Stratum <- 1:nrow(allocations)

##################################################
####   
##################################################

##################################################
####   Simulate survey using the current STRS Design
####   Calculate the shortest path for each boat
##################################################
curr_tour_length <- matrix(nrow = 1000, ncol = 2)
for (iter in 1:nrow(curr_tour_length)) {
  set.seed(iter)
  temp_locs <- sample_stations(strata_labels = allocations$Stratum,
                               nh = allocations$boat2)
  boat_idx <- sample(x = 1:2, size = nrow(temp_locs), replace = TRUE)
  
  for (iboat in 1:2) {
    temp_path <- do_tsp(locs = temp_locs[boat_idx == iboat, ])
    curr_tour_length[iter, iboat] <- temp_path$tour_len 
  }
  
  if(iter%%10 == 0) print(iter)
}

##################################################
####   Simulate survey using an optimized STRS Design
####   Calculate the shortest path for each boat
##################################################
idx <- 2
sol_tour_length <- matrix(nrow = 1000, ncol = 2)
for (iter in 1:nrow(sol_tour_length)) {
  set.seed(iter)
  temp_locs <- sample_stations(strata_labels = strata_list[[idx]]$stratum_id,
                               nh = strata_list[[idx]]$Allocation,
                               grid_strata = res_df[, idx ])
  
  # depth_quants <- median(temp_locs$DEPTH_EFH)
  boat_idx <- sample(x = 1:2, size = nrow(temp_locs), replace = TRUE)
  
  for (iboat in 1:2) {
    temp_path <- do_tsp(locs = temp_locs[boat_idx == iboat, ])
    sol_tour_length[iter, iboat] <- temp_path$tour_len 
  }
  if(iter%%10 == 0) print(iter)
}

##################################################
####   
##################################################
par(mfrow = c(1, 1), mar = c(3, 6, 1, 1), oma = c(0, 0, 0, 0))
boxplot(cbind(rowSums(sol_tour_length),
              rowSums(curr_tour_length)),
        las = 1,
        ylim = c(23000, 26000),
        names = c("Optimized Survey", "Current Survey"))
mtext(side = 2, text = "Total Distance Travelled (km)", line = 4)

shapiro.test(rowSums(sol_tour_length))
shapiro.test(rowSums(curr_tour_length))
var.test(x = rowSums(sol_tour_length),
         y = rowSums(curr_tour_length))
t.test(x = rowSums(sol_tour_length),
       y = rowSums(curr_tour_length))

##################################################
####   Plot an example of a two-boat survey and pathway
##################################################
par(mfrow = c(3, 3), mar = c(0.5, 1, 0.5, 1))

locs <- haul %>%
  filter(YEAR == 1996) %>%
  distinct(STATIONID, lat, lon, BOTTOM_DEPTH, END_LATITUDE, END_LONGITUDE, 
           HAUL, DATE, HAULJOIN, VESSEL) %>%
  arrange(DATE, HAULJOIN) %>%
  rowid_to_column("order")

## How many boats were used in the survey
vessel_ids <- sort(unique(locs$VESSEL))
for(iboat in 1:2){ ## Loop across boats -- start
  
  boat_lab <- paste0("boat_", iboat)
  ##################################################
  ####   Calculation One: Actual Survey Path
  ##################################################
  ## Subset stations by year and the boat actually used in the survey
  surv_pts_y_b <- locs %>% 
    dplyr::filter(VESSEL == vessel_ids[iboat])
  ## reset order of actual stations
  surv_pts_y_b$order <- 1:nrow(surv_pts_y_b)
  locs_y_b <- sf::st_as_sf(x = surv_pts_y_b,
                           coords = c("lon", "lat"),
                           crs = 4326, agr = "constant")
  
  locs_utm <- sp::spTransform(
    x = SpatialPoints(coords = surv_pts_y_b[, c("lon", "lat")],
                      proj4string =  CRS("+proj=longlat +datum=WGS84")) , 
    CRSobj = CRS("+proj=utm +zone=5 +datum=NAD83 +units=km") )
  dist_km_y_b <- 
    matrix(as.numeric(sf::st_distance(locs_y_b) / 1000),
           nrow = nrow(locs_y_b),
           dimnames = list(locs_y_b$order, locs_y_b$order))
  
  ## Calculate cumulative distance travelled by the actual survey
  current_surey_distance <- 0
  for(istation in 2:nrow(locs_y_b)) {
    current_surey_distance[istation] <- dist_km_y_b[istation - 1, istation]
  }
  
  ## Plot actual survey path travelled by the survey boat
  plot(locs_utm@coords,
       type = "l", asp = 1,
       axes = F, ann = F); box()
  points(locs_utm@coords,
         pch = 16,
         cex = 0.5)
  
  tour_length <- round(sum(current_surey_distance))
  legend("topleft",
         paste0("Actual Survey Path, 2019 (Boat ", iboat, ")\n", 
                nrow(locs_y_b), " Stations\n",
                "Station Length: ", round(sum(current_surey_distance)), " km"),
         bty = "n",
         cex = 0.75)
}


set.seed(1)
temp_locs <- sample_stations(strata_labels = allocations$Stratum,
                             nh = allocations$boat2)
boat_idx <- sample(x = 1:2, size = nrow(temp_locs), replace = TRUE)

for (iboat in 1:2) {
  temp_path <- do_tsp(locs = temp_locs[boat_idx == iboat, ])
  locs_b <- temp_locs[boat_idx == iboat, c("E_km", "N_km")]
  
  plot(locs_b[temp_path$path, ], type = "l", asp = 1, 
       axes = F, ann = F); box()
  points(locs_b[temp_path$path, ], pch = 16, cex = 1/2)
  legend("topleft",
         legend = paste0("Current STRS Design (Boat ", iboat, ")\n", 
                         nrow(locs_b), " Stations\n",
                         "Station Length: ", round(temp_path$tour_len), " km"),
         bty = "n",
         cex = 0.75)
}

set.seed(1)
idx <- 2
temp_locs <- sample_stations(strata_labels = strata_list[[idx]]$stratum_id,
                             nh = strata_list[[idx]]$Allocation,
                             grid_strata = res_df[, idx ])

# depth_quants <- median(temp_locs$DEPTH_EFH)
boat_idx <- sample(x = 1:2, size = nrow(temp_locs), replace = TRUE)

for (iboat in 1:2) {
  temp_path <- do_tsp(locs = temp_locs[boat_idx == iboat, ])
  locs_b <- temp_locs[boat_idx == iboat, c("E_km", "N_km")]
  
  plot(locs_b[temp_path$path, ], type = "l", asp = 1,
       axes = F, ann = F); box()
  points(locs_b[temp_path$path, ], pch = 16, cex = 1/2)
  legend("topleft",
         legend = paste0("Optimized STRS Design (Boat ", iboat, ")\n", 
                         nrow(locs_b), " Stations\n",
                         "Station Length: ", round(temp_path$tour_len), " km"),
         bty = "n",
         cex = 0.75)
}

##################################################
####   Save results
##################################################
save(list = c("curr_tour_length", "sol_tour_length"), 
     file = paste0(github_dir, "/results/sim_tours.RData"))

apply(curr_tour_length, MARGIN = 2, FUN = mean)
