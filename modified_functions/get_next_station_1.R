# get_next_station_1
# Survey design 1: Pick stations based on proximity, depth, west-to-east --------
# After completing each survey station, three decision rules:
# 1) which station is closest and unsampled?
# 2) which station is the furthest west unsampled station? (i.e., don't skip over sites going west to east)
# 3) is one of the above 2 deeper than the other? If so, pick the deepest station (to prioritize deeper, longer trawls first per Ned & Wayne)


get_next_station_1 <- function(stationId = boat_plan[i - 1],
                               already_sampled = boat_plan[1:i],
                               distances = distance_matrix_km[paste(boat_plan[i - 1]), ],
                               survey_locations = surv_pts_boat){
  
  depths <- subset(survey_locations, select = c(Id, DEPTH_EFH))
  longs <- subset(survey_locations, select = c(Id, Lon))
  
  names(depths) <- c("Id","depth")
  names(longs) <- c("Id", "lon")
  
  closest <- names(which.min(distances[!names(distances) %in% already_sampled]))
  
  furthest_w_unsampled <- longs %>%
    filter(Id != stationId) %>%
    filter(!Id %in% already_sampled) %>%
    slice_min(lon) %>%
    dplyr::select(Id) %>%
    as.character()
  
  if(closest == furthest_w_unsampled){
    selection = closest} else{
      depth1 <- depths %>% 
        filter(Id == closest) %>% 
        dplyr::select(depth)
      depth2 <- depths %>% 
        filter(Id == furthest_w_unsampled) %>% 
        dplyr::select(depth)
      
      ind <- which.max(c(depth1, depth2))
      selection <- c(closest, furthest_w_unsampled)[ind]
    }
  
  ## Calculate distance between the next station the previous station
  distance_travelled <- distances[selection]
  
  return(list(selection = selection, distance = distance_travelled))
}
