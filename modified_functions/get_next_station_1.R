###############################################################################
## Function:      Calculate next station to go to  based on proximity, depth, 
##                and west-to-east direction
## Author:        Megsie Siple (margaret.siple@noaa.gov)
##                modified by Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Three decision rules to choose the next station
##                  1) which station is closest and unsampled?
##                  2) which station is the furthest west unsampled station? 
##                     (i.e., don't skip over sites going west to east)
##                  3) is one of the above 2 deeper than the other? If so, pick 
##                     the deepest station (to prioritize deeper, longer trawls
##                     first per Ned & Wayne)
##
## Arguments:
##                 stationId: current stationId to find the next station for
##                 already_sampled: character vector of previously sampled 
##                                  stations
##                 distances: numeric vector of distances from stationId to 
##                            the other stations
##                 survey_locations: depths and locations of all stations
###############################################################################

get_next_station_1 <- function(
  stationId = boat_plan[i - 1],
  already_sampled = boat_plan[1:i],
  distances = distance_matrix_km[paste(boat_plan[i - 1]), ],
  survey_locations = surv_pts_boat)
{
  
  ## Isolate depths and longitudes of stations
  depths <- subset(survey_locations, select = c(Id, DEPTH_EFH))
  names(depths) <- c("Id","depth")
  longs <- subset(survey_locations, select = c(Id, Lon))
  names(longs) <- c("Id", "lon")
  
  ## Select closest stationId to stationId that has not already been sampled
  closest <- names(which.min(distances[!names(distances) %in% already_sampled]))
  
  ## Calculate the most western unsampled station 
  furthest_w_unsampled <- longs %>%
    filter(Id != stationId) %>%
    filter(!Id %in% already_sampled) %>%
    slice_min(lon) %>%
    dplyr::select(Id) %>%
    as.character()
  
  ## Ideally, furthest_w_unsampled should be the closest but if not,
  ## choose the deeper station
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
  
  ## Calculate the distance travelled between the current stationId and the 
  ## selected next stationId
  distance_travelled <- distances[selection]
  
  return(list(selection = selection, 
              distance = distance_travelled))
}
