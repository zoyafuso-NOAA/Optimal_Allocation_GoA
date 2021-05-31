###############################################################################
## Function:     do_tsp
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Arguments:    locs: dataframe with fields "E_km", and "N_km" for eastings
##                     and northings
## Description:  Turns locs into a dist class to be used in the TSP solver
##               in the TSP package. Add dummy city to convert the TSP into 
##               a hamiltonian path problem (i.e., remove the roundtrip part
##               of th TSP). 
## Output:       list of the station order and the tour length. 
###############################################################################
do_tsp <- function(locs) {
  
  ## Turn locs into a dist obj. to be used in the TSP solver, create TSP obj.
  dist_km_y_b <- dist(x = locs[, c("E_km", "N_km")])
  tsp_data <- TSP::TSP(x = dist_km_y_b)
  
  ## Insert dummy city to hack the solver to calculate the hamiltonian path
  tsp_data <- TSP::insert_dummy(tsp_data, label = "cut")
  
  ## Solve
  solution <- TSP::solve_TSP(tsp_data, method = "nearest_insertion")
  path <- TSP::cut_tour(solution, "cut")

  return(list(path = as.integer(names(path)), 
              tour_len = attributes(solution)$tour_length))
}