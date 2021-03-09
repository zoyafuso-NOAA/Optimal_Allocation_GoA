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
####    Set up directories here first 
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents",
                       "C:/Users/Zack Oyafuso/Documents",
                       "C:/Users/zack.oyafuso/Work")[which_machine],
                     "/GitHub/Optimal_Allocation_GoA/")

##################################################
####   Load the true density, true index, and spatial domain dataset
##################################################
load(paste0(github_dir,  "data/fit_density.RData"))
load(paste0(github_dir,  "data/fit_index.RData"))
load(paste0(github_dir, "/data/Extrapolation_depths.RData"))

##################################################
####   Constants used throughout all scripts
##################################################

## Years to use
year_set <- 1996:2019
years_included <- c(1, 4, 8, 10, 12, 14, 16, 18, 20, 22, 24)
n_years <- length(years_included)

## Number of sampling grids
n_cells <- nrow(Extrapolation_depths)

## Scientific and common names used in optimization
sci_names_opt <- c("Atheresthes stomias", "Gadus chalcogrammus",
                   "Gadus macrocephalus", "Glyptocephalus zachirus",
                   "Hippoglossoides elassodon", "Hippoglossus stenolepis", 
                   "Lepidopsetta bilineata", "Lepidopsetta polyxystra",
                   "Microstomus pacificus", "Sebastes alutus", "Sebastes B_R",
                   "Sebastes brevispinis", "Sebastes polyspinis", 
                   "Sebastes variabilis", "Sebastolobus alascanus" )

common_names_opt <- c("arrowtooth flounder", "walleye pollock", "Pacific cod",
                      "rex sole", "flathead sole", "Pacific halibut", 
                      "southern rock sole", "northern rock sole", 
                      "Dover sole", "Pacific ocean perch", 
                      "BS and RE rockfishes", 
                      "silvergrey rockfish", "northern rockfish", 
                      "dusky rockfish", "shortspine thornyhead")

ns_opt <- length(sci_names_opt)

## Scientific and common names not used in optimization, but evaluated
## when simulating surveys

sci_names_eval <- c("Anoplopoma fimbria", "Beringraja spp", "Octopus spp",
                    "Pleurogrammus monopterygius", "Sebastes borealis",
                    "Sebastes variegatus", "Squalus suckleyi")

common_names_eval <- c("sablefish", "skates spp", "Octopus spp", 
                       "Atka mackerel", "shortraker rockfish",
                       "harlequin rockfish", "Pacific spiny dogfish")
common_names_eval_labels <- c("sablefish", "skates spp.", "Octopus spp.", 
                              "Atka mackerel", "shortraker rockfish",
                              "harlequin rockfish", "Pacific spiny dogfish")

ns_eval <- length(sci_names_eval)

## In case we need it, all species names together
sci_names_all <- sort(c(sci_names_opt, sci_names_eval))
common_names_all <- c(common_names_opt, common_names_eval)[order(c(sci_names_opt, sci_names_eval))]
common_names_all_labels <- c(common_names_opt, common_names_eval_labels)[order(c(sci_names_opt, sci_names_eval))]

ns_all <- ns_opt + ns_eval

spp_idx_opt <- which(sci_names_all %in% sci_names_opt)
spp_idx_eval <- which(sci_names_all %in% sci_names_eval)

## Sample sizes across 1, 2, and 3 boats
samples <- c(292, 550, 825)
n_boats <- length(samples)

## Number of strata to input into optimization
stratas <- c(10, 15)
n_strata <- length(stratas)

## Specify Management Districts
districts <- data.frame("reg_area" = c("WRA", "CRA", "CRA", "ERA", "ERA"),
                        "district" = c("West", "Chirikof", "Kodiak", 
                                       "Yakutat", "Southeast"),
                        "domainvalue" = 1:5,
                        "W_lon" = c(-170, -159, -154, -147, -140),
                        "E_lon" = c(-159, -154, -147, -140, -132))

n_dom <- nrow(districts)

district_vals <- cut(x = Extrapolation_depths$Lon, 
                     breaks = c(-170, -159, -154, -147, -140, -132), 
                     labels = 1:5)
districts[, c("W_UTM", "E_UTM")] <- 
  do.call(rbind,tapply(X = Extrapolation_depths$E_km, 
                       INDEX = district_vals, 
                       FUN = range) )

## International North Pacific Fisheries Commission statistical areas
inpfc_vals_current <- district_vals
inpfc_vals_current[Extrapolation_depths$stratum %in% 
                        c(10:13, 110:112, 210, 310, 410, 510)] <- 1
inpfc_vals_current[Extrapolation_depths$stratum %in% 
                        c(20:22, 120:122, 220:221, 320, 420, 520)] <- 2
inpfc_vals_current[Extrapolation_depths$stratum %in% 
                        c(30:33, 35, 130:134, 230:232, 330, 430, 530)] <- 3
inpfc_vals_current[Extrapolation_depths$stratum %in% 
                        c(40:41, 140:143, 240:241, 340:341, 440, 540)] <- 4
inpfc_vals_current[Extrapolation_depths$stratum %in% 
                        c(50, 150:151, 250:251, 350:351, 450, 550)] <- 5

## Number of times to simulate survey
n_iters <- 1000
obs_cv <- c(0, 0.25, 0.5, 1) #low to high sampling CVs
n_obs_cv <- length(obs_cv)

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
             X1 = with(Extrapolation_depths, E_km - min(E_km)),
             X2 = Extrapolation_depths$DEPTH_EFH,
             WEIGHT = n_years),
  
  matrix(data = apply(X = D_gct[, , years_included],
                      MARGIN = c(1, 2), 
                      FUN = sum),
         ncol = ns_all,
         dimnames = list(NULL, paste0("Y", 1:ns_all))),
  
  matrix(data = apply(X = D_gct[, , years_included],
                      MARGIN = c(1, 2), 
                      FUN = function(x) sum(x^2)),
         ncol = ns_all,
         dimnames = list(NULL, paste0("Y", 1:ns_all, "_SQ_SUM")))
)

frame_district <- cbind(data.frame(
  domainvalue = cut(x = Extrapolation_depths$Lon, 
                    breaks = c(-170, -159, -154, -147, -140, -132), 
                    labels = 1:5),
  id = 1:n_cells,
  X1 = with(Extrapolation_depths, E_km - min(E_km)),
  X2 = Extrapolation_depths$DEPTH_EFH,
  WEIGHT = n_years),
  
  matrix(data = apply(X = D_gct[, , years_included],
                      MARGIN = c(1, 2), 
                      FUN = sum),
         ncol = ns_all,
         dimnames = list(NULL, paste0("Y", 1:ns_all))),
  
  matrix(data = apply(X = D_gct[, , years_included],
                      MARGIN = c(1, 2), 
                      FUN = function(x) sum(x^2)),
         ncol = ns_all,
         dimnames = list(NULL, paste0("Y", 1:ns_all, "_SQ_SUM")))
)

##################################################
####   Calculate true mean density and true abundance index along with
####   the true abundance index within districts 
##################################################
true_mean <- apply(X = D_gct[, , years_included], 
                   MARGIN = 2:3,
                   FUN = mean)

true_index <- apply(X = Index[, , years_included], 
                    MARGIN = 2:3,
                    FUN = sum)

true_index_district <- apply(X = Index[,, years_included], 
                             MARGIN = 2:3,
                             FUN = function(x) tapply(x, 
                                                      INDEX = district_vals, 
                                                      FUN = sum))
true_index_district <- aperm(a = true_index_district, 
                             perm = c(2,3,1))

##################################################
####   Save Data
##################################################
save(list = c("frame_all", "frame_district",
              "districts", "district_vals", "n_dom", "inpfc_vals_current",
              "true_mean", "true_index", "true_index_district",
              "ns_all", "ns_eval", "ns_opt", 
              "common_names_all", "common_names_eval", "common_names_opt",
              "common_names_eval_labels", "common_names_all_labels",
              "sci_names_all", "sci_names_eval", "sci_names_opt",
              "spp_idx_eval", "spp_idx_opt",
              "year_set", "years_included", "n_years", 
              "n_cells", "samples", "n_boats", "n_iters", 
              "obs_cv", "n_obs_cv",
              "stratas", "n_strata"),
     file = paste0(github_dir, "data/optimization_data.RData"))
