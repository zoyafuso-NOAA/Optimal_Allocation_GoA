###############################################################################
## Project:       Simple Optimization Example
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Run a very simple and short optimization on Gulf of Alaska
##                groundfishes. This is a reduced form of the actual 
##                multispecies optimization.
###############################################################################
rm(list = ls())

##################################################
####   Setup a directories using either an Rproject or by using a variable
####   as shown below
##################################################
github_dir <- paste0("C:/Users/Zack Oyafuso/Documents/GitHub",
                     "/Optimal_Allocation_GoA/")

##################################################
####  Install a forked version of the SamplingStrata Package from 
####  zoyafuso-NOAA's Github page
####
####  Import other required packages
##################################################
library(devtools)
devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata")
library(SamplingStrata)

library(sp)
library(RColorBrewer)
library(raster)

##################################################
####   Load input data
##################################################
load(paste0(github_dir, "data/optimization_data.RData"), verbose = T)
load(paste0(github_dir, "/data/Extrapolation_depths.RData"), verbose = T)

##################################################
####   For speed, we'll subset just three species:
####   Atheresthes stomias, Gadus macrocephalus, Sebastolobus alascanus   
####
####   The frame df contains information for each cell in the spatial domain
####   The X variabces are the habitat variables 
####   X1: scaled UTM eastings
####   X2: depth (m)
####
####   The next three quantities are used to speed up the population stratum 
####   variance caluclation, which includes both spatial and temporal variation
####   across all cells within a stratum:
####
####   YS: VAST-predicted density of species S summed over observed years
####   YS_SQ_SUM: VAST-predicted squared density of species s summed 
####   over observed years
####   WEIGHT: number of observed years
##################################################
which_spp <- c(2, 5, 21)
n_spp <- length(which_spp)

sci_names_all[which_spp]

frame <- frame_all[, c("id", "X1", "X2", "WEIGHT",
                       paste0("Y", which_spp),
                       paste0("Y", which_spp, "_SQ_SUM"),
                       "domainvalue")]
names(frame)[names(frame) %in% paste0("Y", which_spp)] <- paste0("Y", 1:n_spp)
names(frame)[names(frame) %in% paste0("Y", which_spp, "_SQ_SUM")] <- 
  paste0("Y", 1:n_spp, "_SQ_SUM")

head(frame)

##################################################
####   Specify the population CV to optimize over: the optimization will 
###    minimize total sample size subject to these population CV constraints.  
###    Specify one CV (upper-constraint) value for each species of interest
##################################################
CV_constraints = rep(0.1, n_spp)
#Create CV dataframe
cv <- list()
for (spp in 1:n_spp) cv[[paste0("CV", spp)]] <- as.numeric(CV_constraints[spp])
cv[["DOM"]] <- 1
cv[["domainvalue"]] <- 1
(cv <- as.data.frame(cv))

##################################################
####   Run optimization
####   If you want to save the output, first setwd() to the directory you want
####   the output saved to, then turn the writeFiles argumenet to TRUE
####   Iterations are set to 50 for speed but in practice should be in the 
####   hundreds. Population size is 10 for speed but in practice should be 
####   higher (e.g., 30 or 50). See ?SamplingStrata::optimStrata for 
####   descriptions of the other arguments
##################################################
num_of_strata = 10

solution <- SamplingStrata::optimStrata(method = "continuous",
                                        errors = cv, 
                                        framesamp = frame,
                                        iter = 50,
                                        pops = 10,
                                        elitism_rate = 0.1,
                                        mut_chance = 1 / (num_of_strata + 1),
                                        nStrata = num_of_strata,
                                        showPlot = T,
                                        parallel = F,
                                        writeFiles = F)

##################################################
####   Save result objects:
####   Stratum Characteristcs (sum_stats), 
####   optimized CV (opt_CV): may be slightly lower than specified constraints,
####   and which cells belong to each stratum (solution_by_strata)
##################################################
(sum_stats <- summaryStrata(solution$framenew,
                           solution$aggr_strata,
                           progress=FALSE) )

(opt_CV <- expected_CV(strata = solution$aggr_strata))
solution_by_strata <- solution$framenew$STRATO

##################################################
####   Plot coarse map of the solution
##################################################
goa <- sp::SpatialPointsDataFrame(
  coords = Extrapolation_depths[, c("E_km", "N_km")],
  data = data.frame(Str_no = solution_by_strata) )
goa_ras <- raster::raster(goa, 
                          resolution = 5)
goa_ras <- raster::rasterize(x = goa, 
                             y = goa_ras, 
                             field = "Str_no")

plot(goa_ras, 
     axes = F, 
     col = colorRampPalette(brewer.pal(name = "Paired", n = 12))(num_of_strata))
