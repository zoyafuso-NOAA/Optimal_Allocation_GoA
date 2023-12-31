###############################################################################
## Project:       Data synthesis for stratified survey optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Create dataset used for all optimization runs based on a 
##                Gulf of Alaska groundfish VAST spatiotemporal 
##                operating single-species models 
##
##                Calculate true mean density across years for each species
##
##                Set up other constants used in downstream processes
###############################################################################
rm(list = ls())

##################################################
####   Load the true density, true index, and spatial domain dataset
##################################################
load(file = "data/processed/prednll_VAST_models.RData")
grid_goa <- readRDS(file = "data/processed/goa_interpolation_grid.RDS")

##################################################
####   Constants used throughout all scripts
##################################################

## Years to use
year_set <- 1996:2023
years_included <- c(1, 4, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28)
n_years <- length(x = years_included)

## Number of sampling grids
n_cells <- nrow(x = grid_goa)

## Scientific and common names used in optimization
common_names <- pred_jnll$spp_name
ns <- length(x = common_names)

## Specify Management Districts
districts <- data.frame("reg_area" = c("WRA", "CRA", "CRA", "ERA", "ERA"),
                        "district" = c("West", "Chirikof", "Kodiak", 
                                       "Yakutat", "Southeast"),
                        "domainvalue" = 1:5,
                        "W_lon" = c(-170, -159, -154, -147, -140),
                        "E_lon" = c(-159, -154, -147, -140, -132))

district_vals <- cut(x = grid_goa$Lon, 
                     breaks = c(-170, -159, -154, -147, -140, -132), 
                     labels = 1:5)
districts[, c("W_UTM", "E_UTM")] <-
  do.call(rbind,tapply(X = grid_goa$Eastings,
                       INDEX = district_vals,
                       FUN = range) )

n_districts <- nrow(x = districts)

## ranges of the spatial domain for plotting
x_range <- diff(x = range(grid_goa$Eastings))
y_range <- diff(x = range(grid_goa$Northings))

## Number of times to simulate survey
n_iters <- 1000

depth_discrete_cutoff <- c(300, 700)
strata <- 4:5
target_n <- 520

##################################################
####   Save Data, Species densities separately
##################################################
save(list = c("districts", "district_vals", "n_districts", 
              "depth_discrete_cutoff", "strata", "target_n",
              "x_range", "y_range",
              "ns",  "common_names", "n_cells", "n_iters",
              "year_set", "years_included", "n_years"),
     file = "data/processed/optimization_data.RData")

