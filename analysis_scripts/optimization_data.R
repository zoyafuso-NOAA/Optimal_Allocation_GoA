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
load(paste0(github_dir,  "data/fit_Index.RData"))
load(paste0(github_dir, "/data/Extrapolation_depths.RData"))

##################################################
####   Constants used throughout all scripts
##################################################

## Years to use
Year_Set <- 1996:2019
Years2Include <- c(1, 4, 8, 10, 12, 14, 16, 18, 20, 22, 24)
NTime <- length(Years2Include)

## Number of sampling grids
N <- nrow(Extrapolation_depths)

## Scientific and common names used in optimization
sci_names_opt <- c("Atheresthes stomias", "Gadus chalcogrammus",
                   "Gadus macrocephalus", "Glyptocephalus zachirus",
                   "Hippoglossoides elassodon", "Hippoglossus stenolepis", 
                   "Lepidopsetta bilineata", "Lepidopsetta polyxystra",
                   "Microstomus pacificus", "Sebastes alutus", "Sebastes B_R",
                   "Sebastes brevispinis", "Sebastes polyspinis", 
                   "Sebastes variabilis", "Sebastolobus alascanus" )

common_names_opt <- c("arrowtooth flounder", "Alaska pollock", "Pacific cod",
                      "rex sole", "flathead sole", "Pacific halibut", 
                      "southern rock sole", "northern rock sole", 
                      "Pacific Dover sole", "Pacific ocean perch", 
                      "blackspotted/rougheye\nrockfishes", 
                      "silvergrey rockfish", "northern rockfish", 
                      "dusky rockfish", "shortspine thornyhead")

ns_opt <- length(sci_names_opt)

## Scientific and common names not used in optimization, but evaluated
## when simulating surveys

sci_names_eval <- c(  "Anoplopoma fimbria", "Beringraja spp.", 
                      "Octopus spp.", "Pleurogrammus monopterygius",
                      "Sebastes borealis",
                      # "Sebastes ruberrimus",
                      "Sebastes variegatus", "Squalus suckleyi")

common_names_eval <- c("sablefish", "skates spp.", "Octopus spp.", 
                       "Atka mackerel", "shortraker rockfish",
                       # "yelloweye rockfish",
                       "harlequin rockfish", "spiny dogfish")

ns_eval <- length(sci_names_eval)

## In case we need it, all species names together
sci_names_all <- sort(c(sci_names_opt, sci_names_eval))
common_names_all <- c(common_names_opt, 
                      common_names_eval)[order(c(sci_names_opt, 
                                                 sci_names_eval))]
ns_all <- ns_opt + ns_eval

spp_idx_opt <- which(sci_names_all %in% sci_names_opt)
spp_idx_eval <- which(sci_names_all %in% sci_names_eval)

## Sample sizes across 1, 2, and 3 boats
samples <- c(280, 550, 820)
nboats <- length(samples)

## Number of strata to input into optimization
stratas <- c(10, 15, 20)
NStrata <- length(stratas)


## Number of times to simulate survey
Niters <- 1000
obs_CV <- c(0, 0.1, 0.25, 0.5, 1) #low to high sampling CVs
nobs_CV <- length(obs_CV)

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
frame <- cbind(data.frame(domainvalue = 1,
                          id = 1:N,
                          X1 = with(Extrapolation_depths, E_km - min(E_km)),
                          X2 = Extrapolation_depths$DEPTH_EFH,
                          WEIGHT = NTime),
               
               matrix(data = apply(X = D_gct[, , Years2Include],
                                   MARGIN = c(1, 2), 
                                   FUN = sum),
                      ncol = ns_all,
                      dimnames = list(NULL, paste0("Y", 1:ns_all))),
               
               matrix(data = apply(X = D_gct[, , Years2Include],
                                   MARGIN = c(1, 2), 
                                   FUN = function(x) sum(x^2)),
                      ncol = ns_all,
                      dimnames = list(NULL, paste0("Y", 1:ns_all, "_SQ_SUM")))
)

##################################################
####   Calculate "true" mean density and "true" abundance index
##################################################
true_mean <- apply(X = D_gct[,,Years2Include], 
                   MARGIN = 2:3,
                   FUN = mean)

true_index <- apply(X = Index[,, Years2Include], 
                    MARGIN = 2:3,
                    FUN = sum)

##################################################
####   Save Data
##################################################
save(list = c("frame", 
              "true_mean", "true_index", 
              "ns_all", "ns_eval", "ns_opt", 
              "common_names_all", "common_names_eval", "common_names_opt",
              "sci_names_all", "sci_names_eval", "sci_names_opt",
              "spp_idx_eval", "spp_idx_opt",
              "Year_Set", "Years2Include", "NTime", 
              "N", "samples", "nboats", "Niters", 
              "obs_CV", "nobs_CV",
              "stratas", "NStrata"),
     file = paste0(github_dir, "data/optimization_data.RData"))
