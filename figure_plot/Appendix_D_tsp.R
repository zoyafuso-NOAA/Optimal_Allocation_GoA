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
library(readxl)
library(tidyverse)

##################################################
####   Set up directories
##################################################
result_dir <- paste0("C:/Users/Zack Oyafuso/Google Drive/MS_Optimizations/",
                     "TechMemo/appendix/Appendix D plots/")

# result_dir <- paste0("C:/Users/Zack Oyafuso/Desktop/")

if(!dir.exists(result_dir)) dir.create(result_dir)

##################################################
####   Load constants, spatial domain grid, and optimization solutions
##################################################
load(file = "data/processed/grid_goa.RData")
load(file = "data/processed/processed_haul_data.RData")

years_vec <- c(1996, 1999, 2003, 2005, 2007, 
               2009, 2011, 2013, 2015, 2017, 2019)
GOA_allocations <- 
  readxl::read_xlsx(path = 'data/GOA 2019 stations by stratum.xlsx')
GOA_allocations3 <- 
  readxl::read_xlsx(path = 'data/GOA2019_ 3 boat_825_RNDM_stations.xlsx')

##################################
## Specify Management Districts
##################################
new_strata_labels <- 1:length(unique(grid_goa$stratum))
names(new_strata_labels) <- sort(unique(grid_goa$stratum))

##################################
## Rename Current Stratum labels
##################################
grid_goa$stratum_new_label <- new_strata_labels[paste(grid_goa$stratum)]

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
####   Calculate cumulative distance of current stations with ordering of 
####   stations based on the station choice algorithm
##################################################
nearest_dists <- list()
tour_dists <- array(dim = c(length(years_vec), 2, 3),
                    dimnames = list(paste0("year_", years_vec),
                                    c("actual", "tsp"),
                                    paste0("boat_", 1:3)) )

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
  
  png(paste0(result_dir, "Appendix D", which(iyear == years_vec), ".png"),
      width = 7, height = 9,
      units = "in", res = 500)
  
  par(mfrow = c(3, 2), mar = c(3, 3, 1, 1), oma = c(1.5, 1.5, 10, 0))
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
    
    if(b == 1) mtext(side = 3, line = 0.5, text = "Actual Survey Path")
    
    tour_dists[paste0("year_", iyear), 1, b] <- tour_length
    
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
    if(b == 1) {
      
      mtext(side = 3, line = 0.5,
            text = "TSP path\nstations allocated by depth")
      
      space_indent <- ifelse(which(iyear == years_vec) %in% 10:11,
                             "\n                                         ",
                             "\n                                       ") 
      
      text(x = -230.5, 
           y = 65, xpd = NA, adj = 0, family = "serif",
           cex  = 1.5,
           labels = paste0("Appendix Figure D-", which(iyear == years_vec), 
                           ". -- Survey paths observed by the survey (left) ",
                           "for each boat (rows) along with the",
                           space_indent,
                           "optimal shortest path from solving the ",
                           "Travelling Salesperson Problem (right)",
                           space_indent,
                           "using the same stations allocated to each ",
                           "boat in ", iyear, ". Surveys start at the ",
                           space_indent,
                           "most western station and traverse eastward." ))
    }
    
    
  } ## Loop over boats -- end
  
  mtext(side = 1, text = "Longitude", outer = TRUE, font = 2,
        line = ifelse(test = nboats == 3, yes = 0, no = -18))
  mtext(side = 2, text = "Latitude", outer = TRUE, font = 2)
  dev.off()
}  ## Loop over years -- end

##################################################
####   Save
##################################################
# save(list = c("tour_dists", "nearest_dists"),
#      file = paste0(github_dir, "/results/historical_surveys_distances.RData") )
