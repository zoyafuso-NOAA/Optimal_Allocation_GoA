###############################################################################
## Project:       Survey Feasibility (total distance travelled)
## Author:        Megsie Siple (margaret.siple@noaa.gov)
##                Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Simulate several 2-boat surveys station locations
##                Calculate shortest distance path across two boats
##                Calculate total distance travelled
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(TSP)

##################################################
####   Set up directories
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents", 
                       "C:/Users/Zack Oyafuso/Documents",
                       "C:/Users/zack.oyafuso/Work")[which_machine],
                     "/GitHub/Optimal_Allocation_GoA/")

##################################################
####   Load constants, spatial domain grid, and optimization solutions
##################################################
load(file = paste0(github_dir, "data/Extrapolation_depths.RData"))
load(file = paste0(github_dir, "data/processed_haul_data.RData"))
load(file = paste0(github_dir, "results/MS_optimization_knitted_results.RData"))

years_vec <- c(1996, 1999, 2003, 2005, 2007, 
               2009, 2011, 2013, 2015, 2017, 2019)
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
####   Source functions to calculate nearest station
##################################################
# source(file = paste0(github_dir, "modified_functions/get_next_station_1.R"))

##################################################
####   Calculate cumulative distance of current stations with ordering of 
####   stations based on the station choice algorithm
##################################################
nearest_dists <- list()
tour_dists <- array(dim = c(length(years_vec), 2, 3),
                    dimnames = list(paste0("year_", years_vec),
                                    c("actual", "tsp"),
                                    paste0("boat_", 1:3)) )

{
  pdf(paste0("C:/Users/Zack Oyafuso/Google Drive/MS_Optimizations/",
             "Update Documents/Updates 25 May 2021/path_trajectories.pdf"),
      width = 7, height = 8, onefile = T  )
for (iyear in years_vec[]) { ## Loop over years -- start
  
  year_lab <- paste0("year_", iyear)
  ## Subset haul data for a given iyear and order based on the actual station
  ## order observed in the survey
  locs <- haul %>%
    filter(YEAR == iyear) %>%
    distinct(STATIONID, lat, lon, BOTTOM_DEPTH, END_LATITUDE, END_LONGITUDE, 
             HAUL, DATE, HAULJOIN, VESSEL) %>%
    arrange(DATE, HAULJOIN) %>%
    rowid_to_column("order")
  
  ## How many boats were used in the survey
  vessel_ids <- sort(unique(locs$VESSEL))
  nboats <- length(vessel_ids)
  
  ## Depth quantiles used to allocate number of stations across boats for the 
  ## travelling salesman problem
  depth_quants <- quantile(x = locs$BOTTOM_DEPTH)
  locs$whichboat <- sample(x = 1:nboats, size = nrow(locs), replace = TRUE)
  
  par(mfrow = c(3, 2), mar = c(3, 3, 1, 1), oma = c(1, 1, 3.5, 0))
  for(b in 1:nboats){ ## Loop across boats -- start
    
    boat_lab <- paste0("boat_", b)
    ##################################################
    ####   Calculation One: Actual Survey Path
    ##################################################
    ## Subset stations by year and the boat actually used in the survey
    surv_pts_y_b <- locs %>% 
      dplyr::filter(VESSEL == vessel_ids[b])
    ## reset order of actual stations
    surv_pts_y_b$order <- 1:nrow(surv_pts_y_b)
    locs_y_b <- sf::st_as_sf(x = surv_pts_y_b,
                             coords = c("lon", "lat"),
                             crs = 4326, agr = "constant")
    
    dist_km_y_b <- 
      matrix(as.numeric(sf::st_distance(locs_y_b) / 1000),
             nrow = nrow(locs_y_b),
             dimnames = list(locs_y_b$order, locs_y_b$order))
    
    
    ## Calculate cumulative distance travelled by the actual survey
    current_surey_distance <- 0
    temp_var <- nearest_dists$actual_path[[year_lab]][[boat_lab]]
    
    for(istation in 2:nrow(locs_y_b)) {
      current_surey_distance[istation] <- dist_km_y_b[istation - 1, istation]
      
      take_out_idx <- 1:(istation - 1)
      temp_var <- rbind(temp_var, 
                        sort(dist_km_y_b[istation - 1, -take_out_idx])[1:2]  )
    }
    
    dimnames(temp_var) <- list(NULL, c("closest", "second_closest"))   
    nearest_dists$actual_path[[year_lab]][[boat_lab]] <- temp_var
    
    ## Plot actual survey path travelled by the survey boat
    plot(x = 1, y = 1, type = "n",
         xlim = c(-170, -130),
         ylim = c(52, 60), las = 1,
         xlab = "Lon", ylab = "Lat")
    legend("topleft", legend = paste("Boat", b), bty = "n", text.font = 2)
    lines(lat ~ lon,
          data = surv_pts_y_b)
    points(lat ~ lon,
           data = surv_pts_y_b,
           pch = 16,
           cex = 0.5)
    
    tour_length <- round(sum(current_surey_distance))
    legend("bottomright",
           legend = c( paste("Total Length:", tour_length, "km"),
                       paste(nrow(locs_y_b), "Stations")),
           bty = "n")
    
    if(b == 1) mtext(side = 3, text = "Actual Survey Path")
    
    tour_dists[paste0("year_", iyear), 1, b] <- tour_length
    
    ##################################################
    ####   Calculation Two: Calculate the path as calculated using the 
    ####   station decision heuristic and boats separated by depth
    ##################################################
    ## Subset stations by the boat assigned by depths
    # surv_pts_y_b <- locs %>%
    #   dplyr::filter(whichboat == b)
    # surv_pts_y_b$Id <- 1:nrow(surv_pts_y_b)
    # surv_pts_y_b$DEPTH_EFH <- surv_pts_y_b$BOTTOM_DEPTH
    # surv_pts_y_b$Lon <- surv_pts_y_b$lon
    # 
    # locs_y_b <- sf::st_as_sf(x = surv_pts_y_b,
    #                          coords = c("lon", "lat"),
    #                          crs = 4326, agr = "constant")
    # 
    # dist_km_y_b <- matrix(as.numeric(sf::st_distance(locs_y_b) / 1000),
    #                       nrow = nrow(locs_y_b),
    #                       dimnames = list(locs_y_b$Id, locs_y_b$Id))
    # 
    # ## Initialize the station order (boat_plan) with the most western station
    # western_end <- surv_pts_y_b %>%
    #   dplyr::arrange(lon) %>%
    #   dplyr::slice(1)
    # 
    # boat_plan <- boat_distance <- closest_dist <- second_closest_dist <- 
    #   rep(x = NA, times = nrow(surv_pts_y_b))
    # boat_plan[1] <- western_end$Id
    # boat_distance[1] <- 0
    # west_id <- western_end$Id
    # closest_dist[1] <- sort(dist_km_y_b[west_id, ])[2]
    # second_closest_dist[1] <- sort(dist_km_y_b[west_id, ])[3]
    # 
    # # Loop through stations, assigning a "next station" for each one
    # for (istation in 2:nrow(surv_pts_y_b)) { ## loop over stations -- start
    #   calc_next_station <- get_next_station_1(
    #     stationId = boat_plan[istation - 1],
    #     already_sampled = boat_plan[1:(istation - 1)],
    #     distances = dist_km_y_b[boat_plan[istation - 1], ],
    #     survey_locations = surv_pts_y_b
    #   )
    #   
    #   boat_plan[istation] <- calc_next_station$selection
    #   boat_distance[istation] <- calc_next_station$distance
    #   closest_dist[istation] <- calc_next_station$closest_dist
    #   second_closest_dist[istation] <- calc_next_station$second_closest_dist
    # } ## loop over stations -- end
    # 
    # nearest_dists$nn[[year_lab]][[boat_lab]] <- 
    #   cbind(closest_dist, second_closest_dist)
    # 
    # plot(x = 1, y = 1, type = "n",
    #      xlim = c(-170, -130),
    #      ylim = c(52, 60), las = 1,
    #      xlab = "Lon", ylab = "Lat")
    # lines(lat ~ lon,
    #       data = surv_pts_y_b[as.integer(boat_plan), ])
    # points(lat ~ lon,
    #        data = surv_pts_y_b[as.integer(boat_plan), ],
    #        pch = 16,
    #        cex = 0.5)
    # tour_length <- round(sum(boat_distance))
    # 
    # legend("bottomright",
    #        legend = c(paste("Total Length:", tour_length, "km"),
    #                   paste(nrow(surv_pts_y_b), "Stations")),
    #        bty = "n")
    # if(b == 1) mtext(side = 3, 
    #                  text = "Nearest-Neighbor path\nstations allocated by depth")
    
    ##################################################
    ####   Calculation Two: Calculate the shortest path using the TSP package
    ####   Convert the distance matrix to a TSP object
    ####   Insert a "dummy" location to find the shortest hamiltonian path
    ####   (i.e., a path which visits each node in the graph exactly once)
    ####   Solve the TSP and remove the "dummy" location from the tour
    ##################################################
  
    tsp_data <- TSP::TSP(x = dist_km_y_b)
    tsp_data <- TSP::insert_dummy(tsp_data, label = "cut")
    solution <- TSP::solve_TSP(tsp_data, method = "nearest_insertion")
    path <- as.integer(names(TSP::cut_tour(solution, "cut")))
    
    temp_var <- nearest_dists$tsp[[year_lab]][[boat_lab]]
    take_out_idx <- c()
    
    for(istation in 1:length(path)) {
      take_out_idx <- c(take_out_idx, path[1:istation])
      temp_var <- rbind(temp_var, 
                        sort(dist_km_y_b[path[istation], -take_out_idx])[1:2])
    }
    
    dimnames(temp_var) <- list(NULL, c("closest", "second_closest"))   
    nearest_dists$tsp[[year_lab]][[boat_lab]] <- temp_var
    
    ## Plot shortest survey path solved using the TSP
    plot(x = 1, y = 1, type = "n",
         xlim = c(-170, -130),
         ylim = c(52, 60), las = 1,
         xlab = "Lon", ylab = "Lat")
    lines(lat ~ lon,
          data = surv_pts_y_b[path, ])
    points(lat ~ lon,
           data = surv_pts_y_b[path, ],
           pch = 16,
           cex = 0.5)
    
    tour_length <- round(attributes(solution)$tour_length)
    
    tour_dists[paste0("year_", iyear), 2, b] <- tour_length
    
    legend("bottomright",
           legend = c(paste("Total Length:", tour_length, "km"),
                      paste(nrow(surv_pts_y_b), "Stations")),
           bty = "n")
    if(b == 1) mtext(side = 3, 
                     text = "TSP path\nstations allocated by depth")
  }
  
  mtext(side = 3, 
        line = 2,
        font = 2,
        text = iyear, 
        outer = TRUE)
  
}  ## Loop over years -- end

dev.off()
}

##################################################
####   Save
##################################################
save(list = c("tour_dists", "nearest_dists"),
     file = paste0(github_dir, "/results/historical_surveys_distances.RData") )
