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
load("data/processed/VAST_fit_D_gct.RData")
load( "data/processed/VAST_fit_I_gct.RData")
load("data/processed/prednll_VAST_models.RData")
load("data/processed/grid_goa.RData")

##################################################
####   Constants used throughout all scripts
##################################################

## Years to use
year_set <- 1996:2019
years_included <- c(1, 4, 8, 10, 12, 14, 16, 18, 20, 22, 24)
n_years <- dim(D_gct)[3]

## Number of sampling grids
n_cells <- nrow(grid_goa)

## Scientific and common names used in optimization
common_names_all <- pred_jnll$spp_name

ns_all <- length(common_names_all)

spp_idx_opt <- c(25, 14, #cods 
                 1, 7, 18, 12, 24, 5, 15, #flatfishes
                 16, 4, 23, 6, 13, 22 #rockfish types
                 )
common_names_opt <- common_names_all[spp_idx_opt]
ns_opt <- length(common_names_opt)

## Scientific and common names not used in optimization, but evaluated
## when simulating surveys
spp_idx_eval <- (1:ns_all)[-spp_idx_opt]
common_names_eval <- common_names_all[spp_idx_eval]
ns_eval <- length(common_names_eval)

## Sample sizes across 1, 2, and 3 boats
samples <- c(292, 550, 825)
n_boats <- length(samples)

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
  do.call(rbind,tapply(X = grid_goa$E_km, 
                       INDEX = district_vals, 
                       FUN = range) )

n_districts <- nrow(districts)

## International North Pacific Fisheries Commission statistical areas
inpfc_vals_current <- district_vals
inpfc_vals_current[grid_goa$stratum %in% 
                        c(10:13, 110:112, 210, 310, 410, 510)] <- 1
inpfc_vals_current[grid_goa$stratum %in% 
                        c(20:22, 120:122, 220:221, 320, 420, 520)] <- 2
inpfc_vals_current[grid_goa$stratum %in% 
                        c(30:33, 35, 130:134, 230:232, 330, 430, 530)] <- 3
inpfc_vals_current[grid_goa$stratum %in% 
                        c(40:41, 140:143, 240:241, 340:341, 440, 540)] <- 4
inpfc_vals_current[grid_goa$stratum %in% 
                        c(50, 150:151, 250:251, 350:351, 450, 550)] <- 5

## ranges of the spatial domain for plotting
x_range <- diff(range(grid_goa$E_km))
y_range <- diff(range(grid_goa$N_km))

## Number of times to simulate survey
n_iters <- 1000

##################################################
####   Our df will have fields for:
####   domain: only one domain so the value is just 1
####   id: unique ID for each sampling cell
####   X1: strata variable 1: longitude in eastings (km). Because the 
####           optimization does not read in negative values, I shift the 
####           values so that the lowest value is 0
####   X2: strata variable 2: depth of cell (m)
####
####   Variables used to more efficiently calcualte stratum variance 
####
####   WEIGHT: number of observed years 
####   Y1, Y2, ... : density for a given cell summed across observed years
####   Y1_SQ_SUM, Y2_SQ_SUM, ... : density-squared for a given cell, 
####           summed across observed years
##################################################
frame_all <- cbind(
  data.frame(domainvalue = 1,
             id = 1:n_cells,
             X1 = with(grid_goa, E_km - min(E_km)),
             X2 = grid_goa$DEPTH_EFH,
             WEIGHT = n_years),
  
  matrix(data = apply(X = D_gct,
                      MARGIN = c(1, 2), 
                      FUN = sum),
         ncol = ns_all,
         dimnames = list(NULL, paste0("Y", 1:ns_all))),
  
  matrix(data = apply(X = D_gct,
                      MARGIN = c(1, 2), 
                      FUN = function(x) sum(x^2)),
         ncol = ns_all,
         dimnames = list(NULL, paste0("Y", 1:ns_all, "_SQ_SUM")))
)

frame_district <- cbind(data.frame(
  domainvalue = cut(x = grid_goa$Lon, 
                    breaks = c(-170, -159, -154, -147, -140, -132), 
                    labels = 1:5),
  id = 1:n_cells,
  X1 = with(grid_goa, E_km - min(E_km)),
  X2 = grid_goa$DEPTH_EFH,
  WEIGHT = n_years),
  
  matrix(data = apply(X = D_gct,
                      MARGIN = c(1, 2), 
                      FUN = sum),
         ncol = ns_all,
         dimnames = list(NULL, paste0("Y", 1:ns_all))),
  
  matrix(data = apply(X = D_gct,
                      MARGIN = c(1, 2), 
                      FUN = function(x) sum(x^2)),
         ncol = ns_all,
         dimnames = list(NULL, paste0("Y", 1:ns_all, "_SQ_SUM")))
)

##################################################
####   Calculate true mean density and true abundance index along with
####   the true abundance index within districts 
##################################################
true_mean <- apply(X = D_gct, 
                   MARGIN = 2:3,
                   FUN = mean)

true_index <- apply(X = I_gct, 
                    MARGIN = 2:3,
                    FUN = sum)

true_index_district <- apply(X = I_gct, 
                             MARGIN = 2:3,
                             FUN = function(x) tapply(x, 
                                                      INDEX = district_vals, 
                                                      FUN = sum))
true_index_district <- aperm(a = true_index_district, 
                             perm = c(2,3,1))

dimnames(true_index)[[1]] <- dimnames(true_mean)[[1]] <- 
  dimnames(true_index_district)[[1]] <- common_names_all


##################################################
####   Save Data
##################################################
save(list = c("frame_all", "frame_district",
              "districts", "district_vals", "n_districts", "inpfc_vals_current",
              "true_mean", "true_index", "true_index_district",
              "ns_all", "ns_eval", "ns_opt", 
              "x_range", "y_range",
              "common_names_all", "common_names_eval", "common_names_opt",
              "spp_idx_eval", "spp_idx_opt",
              "year_set", "years_included", "n_years", 
              "n_cells", "samples", "n_boats", "n_iters"),
     file = "data/processed/optimization_data.RData")
