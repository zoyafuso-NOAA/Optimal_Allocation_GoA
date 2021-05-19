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
nboats <- 2  ## Only considering two-boat solutions
load(file = paste0(github_dir, "data/Extrapolation_depths.RData"))
load(file = paste0(github_dir, "results/MS_optimization_knitted_results.RData"))

##################################################
####   Source functions to calculate nearest station
##################################################
source(file = paste0(github_dir, "modified_functions/get_next_station_1.R"))

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
      closest_dist[1] <- sort(distance_matrix_km[paste(western_end$Id), ])[2]
      second_closest_dist[1] <- sort(distance_matrix_km[paste(western_end$Id), ])[3]
      
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
