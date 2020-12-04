###############################################################################
## Project:       Spatiotemporal Survey Optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Conduct SamplingStrata R package multispecies stratified
##                survey optimization
###############################################################################
rm(list = ls())

##################################################
####   Set up directories based on whether the optimization is being conducted
####        on a multi-species or single-species level
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[3]

github_dir <- paste0(c("/Users/zackoyafuso/Documents", 
                       "C:/Users/Zack Oyafuso/Documents",
                       "C:/Users/zack.oyafuso/Work")[which_machine],
                     "/GitHub/Optimal_Allocation_GoA/results/",
                     "Single_Species_Optimization/")

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
####   Load Data
####   Load Population CVs for use in the thresholds
##################################################
load(paste0(dirname(dirname(github_dir)), "/data/optimization_data.RData"))
load(paste0(dirname(dirname(github_dir)), "/data/Extrapolation_depths.RData"))
load(paste0(dirname(github_dir), "/Population_Variances.RData"))

##################################################
####   Some Constants
##################################################
#Choose a boat level

which_species <- 3
frame <- frame[, c("domainvalue", "id", "X1", "X2", "WEIGHT",
                   paste0("Y", which_species), 
                   paste0("Y", which_species, "_SQ_SUM"))]

names(frame)[6:7] <- paste0("Y", c("1", "1_SQ_SUM") )

github_dir = paste0(github_dir, 
                    gsub(x = sci_names_opt[which_species], 
                         pattern = ' ', 
                         replacement = '_'), '/')
if(!dir.exists(github_dir)) dir.create(github_dir, recursive = T)


##################################################
####   Run optimization
##################################################
par(mfrow = c(6,6), 
    mar = c(2,2,0,0))

for (isample in 1:nboats) {
  #Create CV dataframe
  cv <- list()
  for (spp in 1:1) 
    cv[[paste0("CV", spp)]] <- SRS_Pop_CV[[isample]][, which_species] * 0.90^4
  cv[["DOM"]] <- 1:5
  cv[["domainvalue"]] <- 1:5
  cv <- as.data.frame(cv)
  
  ##Initial Condition
  Run <- 1
  current_n <- 0
  
  while (current_n <= c(280, 550, 820)[isample] ) {
    
    #Set wd for output files, create a directory if it doesn"t exist yet
    temp_dir = paste0(github_dir, "boat", isample, "/Run", Run)
    if(!dir.exists(temp_dir)) dir.create(temp_dir, recursive = T)
    
    setwd(temp_dir)
    
    #Run optimization
    solution <- optimStrata(method = "continuous",
                            errors = cv, 
                            framesamp = frame,
                            iter = 300,
                            pops = 300,
                            elitism_rate = 0.1,
                            mut_chance = 1 / (c(5,5,5,5,5) + 1),
                            nStrata = c(5,5,5,5,5),
                            showPlot = T,
                            writeFiles = T)
    
    sum_stats <- summaryStrata(solution$framenew,
                               solution$aggr_strata,
                               progress=FALSE) 
    
    plot_solution <- as.factor(paste(solution$framenew$DOMAINVALUE,
                                     solution$framenew$STRATO))
    plot_solution <- as.integer(plot_solution)
    
    ##Save a plot of the olution
    goa <- sp::SpatialPointsDataFrame(
      coords = Extrapolation_depths[, c("Lon", "Lat")],
      data = data.frame(Str_no = plot_solution) )
    goa_ras <- raster::raster(x = goa, 
                              resolution = 0.075)
    goa_ras <- raster::rasterize(x = goa, 
                                 y = goa_ras, 
                                 field = "Str_no")
    
    png(filename = "solution.png",
        width = 5,
        height = 5,
        units = "in",
        res = 500)
    
    plot(goa_ras, 
         axes = F, 
         col = rep(brewer.pal(n = 10, name = "Spectral"), times = 3) )
    abline(v = districts$E_lon)
    
    dev.off()
    
    ## Save Output
    CV_constraints <- expected_CV(strata = solution$aggr_strata)
    current_n <- sum(sum_stats$Allocation)
    
    result_list <- list(solution = solution, 
                        sum_stats = sum_stats, 
                        CV_constraints = CV_constraints, 
                        n = current_n)
    save(list = "result_list", file = "result_list.RData")
    print(current_n)
    
    ## Set up next run by changing upper CV constraints
    Run <- Run + 1
    CV_constraints <- 0.9*CV_constraints
    
    #Create CV dataframe in the formmat of SamplingStrata
    cv <- list()
    cv[[paste0("CV1")]] <- as.numeric(CV_constraints)
    cv[["DOM"]] <- 1:5
    cv[["domainvalue"]] <- 1:5
    cv <- as.data.frame(cv)
    
  }
}


