###############################################################################
## Function:     sample_stations
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Arguments:    strata_labels: unique labels for each stratum
##               nh: sample allocations across strata
##               grid_strata: stratum labels of each grid cell
##               grid_locs: eastings, westings, and depth of each gridd cell
## Description:  sample cells across strata according to sample allocations 
## Output:       location, stratum id, and depth of sampled stations
###############################################################################

sample_stations <- function(
  strata_labels, 
  nh, 
  grid_strata = Extrapolation_depths$stratum_new_label,
  grid_locs = Extrapolation_depths[, c("E_km", "N_km", "DEPTH_EFH")]) {
  
  if(length(strata_labels) != length(nh)) 
    stop("strata_labels and nh are not the same length")
  
  chosen_stations <- c()
  for (h in 1:length(nh)) {
    chosen_stations <- c(chosen_stations, 
                         sample(x = which(grid_strata == strata_labels[h]), 
                                size = nh[h]))
  }
  
  return_val <- cbind(stratum = rep(strata_labels, nh),
                      grid_locs[chosen_stations, ]) 
  rownames(return_val) <- 1:nrow(return_val)
  return( return_val )
}
