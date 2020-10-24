###############################################################################
## Project:       Spatiotemporal Survey Optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Conduct SamplingStrata R package multispecies stratified
##                survey optimization
###############################################################################
rm(list = ls())

##################################################
####    Import required packages
##################################################
library(sp)
library(RColorBrewer)
library(raster)

##################################################
####   Set up directories
####
####   Set up some constants of the optimization
####   Multispeceis: Spatiotemporal Variance, species specific CV constraints
####   Single_Species: Spatiotemporal Variance, univariate optimization, 
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[3]
VAST_model <- "11" 

SamplingStrata_dir <- paste0(c("/Users/zackoyafuso/",
                               "C:/Users/Zack Oyafuso/",
                               "C:/Users/zack.oyafuso/")[which_machine],
                             "Downloads/SamplingStrata-master/R")

which_method = c("Multi_Species" = 1,
                 "Single_Species" = 2)[1]

github_dir <- paste0(c("/Users/zackoyafuso/Documents", 
                       "C:/Users/Zack Oyafuso/Documents",
                       "C:/Users/zack.oyafuso/Work")[which_machine],
                     "/GitHub/Optimal_Allocation_GoA/model_", 
                     VAST_model, "/", 
                     c("Spatiotemporal_Optimization/", 
                       "Single_Species_Optimization/")[which_method])

##################################################
####   Load functions from SamplingStrata packages into global environment
####   Load modified buildStrataDF function to incorporate spatiotemporal 
####   stratum variance instead of spatial variance
##################################################
for (ifile in dir(SamplingStrata_dir, full.names = T)) source(ifile)
source(paste0(dirname(dirname(github_dir)), 
              "/modified_functions/buildStrataDF_Zack.R"))

##################################################
####   Load Data
####   Load Population CVs for use in the thresholds
##################################################
load(paste0(dirname(github_dir), "/optimization_data.RData"))
load(paste0(dirname(dirname(github_dir)), "/data/Extrapolation_depths.RData"))
load(paste0(dirname(github_dir), "/Population_Variances.RData"))

##################################################
####   Some Constants
##################################################
stratas <- c(5,10,15,20,30,60)
ns <- c(15, 1)[which_method]

##################################################
####   If Single_Species: subset just the one species
##################################################
if (which_method == 2) {
  SS_which_species <- 13 #which species are we doing?
  
  frame <- frame[,c("id", "X1", "X2", paste0("Y", SS_which_species),
                    "domainvalue")]
  
  frame_raw <- frame_raw[,c("id", "X1", "X2", 
                            paste0("Y", SS_which_species),
                            "domainvalue", "year")]
  
  names(frame)[4] <- names(frame_raw)[4] <- "Y1"
  
  github_dir = paste0(github_dir, gsub(x = sci_names[SS_which_species], 
                                       pattern = ' ', 
                                       replacement = '_'), '/')
  if(!dir.exists(github_dir)) dir.create(github_dir)
  
  # Lower CV threshold is not needed for a single-species analysis
  threshold <- matrix(data = 0,
                      nrow = ns,
                      ncol = 3)
}

##################################################
####   Run optimization
##################################################
par(mfrow = c(6,6), 
    mar = c(2,2,0,0))
isample <- 3

for (istrata in 3) {
  
  temp_strata <- stratas[istrata]
  
  ##Initial Condition
  Run <- 1
  current_n <- 0
  CV_constraints <- SRS_Pop_CV[, isample] 
  
  #Create CV dataframe
  cv <- list()
  for (spp in 1:ns) cv[[paste0("CV", spp)]] <- as.numeric(CV_constraints[spp])
  cv[["DOM"]] <- 1
  cv[["domainvalue"]] <- 1
  cv <- as.data.frame(cv)
  
  while (current_n <= c(280, 550, 820)[isample] ) {
    
    #Set wd for output files, create a directory if it doesn"t exist yet
    temp_dir = paste0(github_dir, "boat", isample, "/Str", temp_strata, 
                      "Run", Run)
    if(!dir.exists(temp_dir)) dir.create(temp_dir, recursive = T)
    
    setwd(temp_dir)
    
    #Run optimization
    solution <- optimStrata(method = "continuous",
                            errors = cv, 
                            framesamp = frame,
                            iter = 300,
                            pops = 50,
                            elitism_rate = 0.1,
                            mut_chance = 1 / (temp_strata + 1),
                            nStrata = temp_strata,
                            showPlot = T,
                            writeFiles = T)
    
    sum_stats <- summaryStrata(solution$framenew,
                               solution$aggr_strata,
                               progress=FALSE) 
    
    #Plot Solution
    goa <- sp::SpatialPointsDataFrame(
      coords = Extrapolation_depths[,c("E_km", "N_km")],
      data = data.frame(Str_no = solution$framenew$STRATO) )
    goa_ras <- raster::raster(x = goa, 
                              resolution = 5)
    goa_ras <- raster::rasterize(x = goa, 
                                 y = goa_ras, 
                                 field = "Str_no")
    
    png(filename = "solution.png", 
        width = 5, 
        height = 5, 
        units = "in", 
        res = 500)
    plot(goa_ras, axes = F, 
         col = terrain.colors(temp_strata)[sample(temp_strata)])
    dev.off()
    
    #Save Output
    CV_constraints <- expected_CV(strata = solution$aggr_strata)
    current_n <- sum(sum_stats$Allocation)

    result_list <- list(solution = solution, 
                        sum_stats = sum_stats, 
                        CV_constraints = CV_constraints, 
                        n = current_n)
    save(list = "result_list", file = "result_list.RData")
    
    #Set up next run by changing upper CV constraints
    Run <- Run + 1
    
    CV_constraints <- 0.9*CV_constraints + 0.1*(SS_STRS_Pop_CV[, isample]) 
    
    #Create CV dataframe in the formmat of SamplingStrata
    cv <- list()
    for (spp in 1:ns) cv[[paste0("CV", spp)]] <- as.numeric(CV_constraints[spp])
    cv[["DOM"]] <- 1
    cv[["domainvalue"]] <- 1
    cv <- as.data.frame(cv)
  }
}

