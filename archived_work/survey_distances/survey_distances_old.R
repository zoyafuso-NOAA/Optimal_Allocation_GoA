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
####   Set up cores
##################################################
library(foreach)
library(parallel)
library(doParallel)

cl <- parallel::makeCluster(parallel::detectCores() - 1)
doParallel::registerDoParallel(cl)

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

##################################################
####   Source functions to calculate nearest station
##################################################
source(file = paste0(github_dir, "modified_functions/get_next_station_1.R"))

##################################################
####   Calculate cumulative distance of current stations with ordering of 
####   stations based on the station choice algorithm
##################################################
df_list <- NULL

for (iyear in years_vec[1]) { ## Loop over years -- start
  
  ## Subset for a give year iyear
  locs <- haul %>%
    filter(YEAR == iyear) %>%
    distinct(STATIONID, lat, lon, BOTTOM_DEPTH, END_LATITUDE, END_LONGITUDE, 
             HAUL, DATE, HAULJOIN, VESSEL) %>%
    arrange(DATE, HAULJOIN) %>%
    rowid_to_column("order")
  
  ## How many boats were used in the survey
  nboats <- length(unique(locs$VESSEL))
  
  ## Depth quantiles used to allocate number of stations across boats
  depth_quants <- quantile(x = locs$BOTTOM_DEPTH)
  locs$whichboat <- dplyr::case_when(
    ## If two boats, split stations by the median
    nboats == 2 ~ as.integer(cut(x = locs$BOTTOM_DEPTH, 
                                 breaks = depth_quants[paste0(c(0, 50, 100), "%")], 
                                 labels = c(1, 2), 
                                 include.lowest =  T)),
    
    ## If three boats, split boats by the 25th and 75th perceniles
    ## Note: boat two will have more stations allocated
    nboats == 3 ~ as.integer(cut(x = locs$BOTTOM_DEPTH, 
                                 breaks = depth_quants[paste0(c(0, 25, 75, 100), "%")], 
                                 labels = c(1, 2, 3), 
                                 include.lowest =  T)))
  
  for(b in 1:nboats){ ## Loop across boats -- start
    ## Subset stations for boat b
    surv_pts_boat <- locs %>% 
      dplyr::filter(whichboat == b)
    
    surv_pts_boat$DEPTH_EFH <- surv_pts_boat$BOTTOM_DEPTH
    surv_pts_boat$Lon <- surv_pts_boat$lon
    surv_pts_boat$Id <- 1:nrow(surv_pts_boat)
    
    ## Subset most western station for boat b
    western_end <- surv_pts_boat %>%
      dplyr::arrange(lon) %>%
      dplyr::slice(1)
    
    ## Subset total sample size for boat b
    sample_size <- nrow(surv_pts_boat)
    
    
    ##
    locs_sf <- sf::st_as_sf(x = surv_pts_boat,
                            coords = c("lon", "lat"),
                            crs = 4326, agr = "constant")
    
    distance_matrix_km <- matrix(as.numeric(sf::st_distance(locs_sf) / 1000),
                                 nrow = nrow(locs_sf),
                                 dimnames = list(locs_sf$Id, locs_sf$Id))
    
    ## Initialize the station order (boat_plan) with the most western station
    boat_plan <- boat_distance <- closest_dist <- second_closest_dist <- 
      rep(x = NA, times = sample_size)
    boat_plan[1] <- western_end$Id
    boat_distance[1] <- 0
    west_id <- western_end$Id
    closest_dist[1] <- sort(distance_matrix_km[west_id, ])[2]
    second_closest_dist[1] <- sort(distance_matrix_km[west_id, ])[3]
    
    
    # Loop through stations, assigning a "next station" for each one
    for (i in 2:sample_size) { ## loop over stations -- start
      calc_next_station <- get_next_station_1(
        stationId = boat_plan[i - 1],
        already_sampled = boat_plan[1:(i-1)],
        distances = distance_matrix_km[boat_plan[i - 1], ],
        survey_locations = surv_pts_boat
      )
      
      boat_plan[i] <- calc_next_station$selection
      boat_distance[i] <- calc_next_station$distance
      closest_dist[i] <- calc_next_station$closest_dist
      second_closest_dist[i] <- calc_next_station$second_closest_dist
    } ## loop over stations -- end
    
    df_list <- rbind(df_list,
                     data.frame(year = iyear,
                                boatno = b,
                                lat = surv_pts_boat$lat[as.integer(boat_plan)],
                                lon = surv_pts_boat$lon[as.integer(boat_plan)],
                                Id = boat_plan,
                                nwd_order = 1:sample_size,
                                distance = boat_distance,
                                closest_dist,
                                second_closest_dist))
      
  }  ## Loop across boats -- end
  
  print(paste("Done with", iyear))
}  ## Loop over years -- end

dplyr::spread(data = aggregate(distance ~ boadno + year,
                               data = df_list,
                               FUN = function(x) round(sum(x))),
              value = distance,
              key = boadno)

##################################################
####  Pick solution and get survey information 
####  Query which solution to use based on 2-boat, district-level solution 
####  with three strata per district
##################################################
idx <- which(settings$strata == 3 & 
               settings$boat == 2 & 
               settings$domain == "district")

##################################################
####  Extract characteristics of the solution
##################################################
strata_no <- strata_list[[idx]]$stratum_id # stratum "ID"
nh <- strata_list[[idx]]$Allocation # stratum effort (n of locations to visit)
nstrata <- length(x = nh)
solution <- res_df[, idx] # Optimized solution

##################################################
####   Result object
##################################################
sims_list <- list()
nboats <- 2  ## Only considering two-boat solutions

##################################################
####   Simulation
##################################################
sims_list <- foreach (sim = 1:1000, 
                      .combine = rbind, 
                      .packages = "tidyverse") %dopar% 
  { # loop over replicates -- start
    
    ## Set seed
    set.seed(2323 + sim)
    
    ## Randomly draw samples from the optimized survey solution
    sample_vec <- c()
    for (istrata in 1:nstrata) {
      sample_vec <- c(sample_vec, 
                      sample(x = which(solution == istrata),
                             size = nh[istrata]))
    }
    
    ## These are the stations to visit
    survey_pts <- Extrapolation_depths[sample_vec, c("Lon", "Lat", "E_km", 
                                                     "N_km", "Id", "stratum",
                                                     "trawlable", "DEPTH_EFH")]
    
    ## Convert survey location points to a sf class
    survey_sf <- sf::st_as_sf( x = survey_pts,
                               coords = c("Lon", "Lat"),
                               crs = 4326,
                               agr = "constant")
    
    ## Pairwise distances between survey points in km
    ## st_distance() provides distances in m
    distance_matrix_km <-
      matrix(data = as.numeric(sf::st_distance(x = survey_sf) / 1000),
             nrow = nrow(x = survey_sf),
             dimnames = list(survey_sf$Id, survey_sf$Id))
    
    ## Assign total number of stations for each boat. Assume n = 2 boats for now
    nperboat <- nrow(survey_pts) / 2
    if (!is.integer(x = nperboat)) {
      n1 <- floor(x = nperboat)
      n2 <- ceiling(x = nperboat)
      n3 <- 0
    } else {
      (n1 <- n2 <- nperboat)
      n3 <- 0
    }
    
    ## Assign shallow stations (i.e., stations less than median depth) to boat 1
    ## and deeper tations to boat 2
    depth_quants <- quantile(x = survey_pts$DEPTH_EFH)
    survey_pts <- survey_pts %>%
      #assign stations to boats by depth:
      dplyr::mutate(whichboat = ifelse(test = DEPTH_EFH < depth_quants["50%"], 
                                       yes = 1, 
                                       no = 2))
    
    ## Setup survey plan for multiple boats
    df_list <- list()
    
    for(b in 1:nboats){
      ## Subset stations for boat b
      surv_pts_boat <- survey_pts %>% 
        dplyr::filter(whichboat == b)
      
      ## Subset most western station for boat b
      western_end <- surv_pts_boat %>%
        dplyr::arrange(Lon) %>%
        dplyr::slice(1)
      
      ## Subset total sample size for boat b
      sample_size = c(n1, n2, n3)[b]
      
      ## Initialize the station order (boat_plan) with the most western station
      boat_plan <- boat_distance <- closest_dist <- second_closest_dist <- 
        rep(x = NA, times = sample_size)
      boat_plan[1] <- western_end$Id
      boat_distance[1] <- 0
      west_id <- paste(western_end$Id)
      closest_dist[1] <- sort(distance_matrix_km[west_id, ])[2]
      second_closest_dist[1] <- sort(distance_matrix_km[west_id, ])[3]
      
      # Loop through stations, assigning a "next station" for each one
      for (i in 2:length(boat_plan)) { ## loop over stations -- start
        calc_next_station <- get_next_station_1(
          stationId = boat_plan[i - 1],
          already_sampled = boat_plan[1:i],
          distances = distance_matrix_km[paste(boat_plan[i - 1]), ]
        )
        
        boat_plan[i] <- calc_next_station$selection
        boat_distance[i] <- calc_next_station$distance
        closest_dist[i] <- calc_next_station$closest_dist
        second_closest_dist[i] <- calc_next_station$second_closest_dist
      } ## loop over stations -- end
      
      df_list[[b]] <- data.frame(boadno = b,
                                 Id = boat_plan,
                                 nwd_order = 1:sample_size,
                                 distance = boat_distance,
                                 closest_dist,
                                 second_closest_dist,
                                 whichsim = sim)
    }
    
    df_list
    
  }  ## Loop over replicate -- end


##################################################
####   Stop core cluster
##################################################
parallel::stopCluster(cl)

##################################################
####   Summarize simulation results
##################################################
cum_dists <- with(do.call(rbind.data.frame, sims_list),
                  tapply(X = distance, 
                         INDEX = list(whichsim, boadno), 
                         FUN = sum))

cum_dists_sum <- round(apply(X = cum_dists, 
                             MARGIN = 2, 
                             FUN = function(x) cbind(mean(x), sd(x))))

##################################################
####   Save Results
##################################################
save(list = c("cum_dists", "cum_dists_sum", "sims_list"),
     file = paste0(github_dir, "/results/survey_distances.RData"))
