###############################################################################
## Project:       Spatiotemporal Survey Optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Conduct SamplingStrata R package multispecies stratified
##                survey optimization
###############################################################################
rm(list = ls())

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
####   Set up directories based on whether the optimization is being conducted
####        on a multi-species or single-species level
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[3]

github_dir <- paste0(c("/Users/zackoyafuso/Documents", 
                       "C:/Users/Zack Oyafuso/Documents",
                       "C:/Users/zack.oyafuso/Work")[which_machine],
                     "/GitHub/Optimal_Allocation_GoA/")

##################################################
####   Load Data
####   Load Population CVs for use in the thresholds
##################################################
load(paste0(github_dir, "/data/optimization_data.RData"))
load(paste0(github_dir, "/data/Extrapolation_depths.RData"))

##################################################
####   Constants to specify before doing optimization
##################################################
which_domain <- c("full_domain", "district")[2]
district_vals <- switch(which_domain,
                        "full_domain" = rep(1, n_cells), 
                        "district" = district_vals)
n_dom <- length(unique(district_vals))

frame <- switch( which_domain,
                 "full_domain" = frame_all,
                 "district" = frame_district)[, c("domainvalue", "id", 
                                                  "X1", "X2", "WEIGHT",
                                                  paste0("Y", spp_idx_opt), 
                                                  paste0("Y", spp_idx_opt,
                                                         "_SQ_SUM"))]

names(frame)[names(frame) %in% paste0("Y", spp_idx_opt)] <- 
  paste0("Y", 1:ns_opt)
names(frame)[names(frame) %in% paste0("Y", spp_idx_opt, "_SQ_SUM")] <- 
  paste0("Y", 1:ns_opt, "_SQ_SUM")

##################################################
####   Run optimization
##################################################
par(mfrow = c(6,6), 
    mar = c(2,2,0,0))

#Choose a boat level
isample <- 1
istrata <- 1

for (istrata in 1) {
  
  temp_strata <- rep(stratas[istrata], n_dom)
  
  ##Initial Condition
  Run <- 1
  current_n <- 0
  # CV_constraints <- SRS_Pop_CV[, isample] 
  
  ## Initiate CVs to be those calculated under SRS
  srs_stats <- SamplingStrata::buildStrataDF( 
    dataset = cbind( subset(frame, select = -c(X1, X2)),
                     X1 = 1))
  
  srs_n <- as.numeric(samples[isample] * table(district_vals) / n_cells)
  
  srs_var <- srs_stats[, paste0("S", 1:ns_opt)]^2 * (1 - srs_n / n_cells) / srs_n
  srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:ns_opt)]
  # names(srs_cv) <- paste0("S", 1:ns_opt)
  
  cv <- list()
  for (spp in 1:ns_opt) 
    cv[[paste0("CV", spp)]] <- as.numeric(srs_cv[, spp])
  cv[["DOM"]] <- 1:n_dom
  cv[["domainvalue"]] <- 1:n_dom
  cv <- as.data.frame(cv)
  
  
  while (current_n <= c(280, 550, 820)[isample] ) {
    
    #Set wd for output files, create a directory if it doesn"t exist yet
    temp_dir = paste0(github_dir, 
                      "results/", which_domain, "/Multi_Species_Optimization",
                      "/boat", isample, "/Str", temp_strata, "Run", Run)
    if(!dir.exists(temp_dir)) dir.create(temp_dir, recursive = T)
    
    setwd(temp_dir)
    
    #Run optimization
    solution <- optimStrata(method = "continuous",
                            errors = cv, 
                            framesamp = frame[],
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
    
    CV_constraints <- 0.95*CV_constraints + 0.05*(SS_STRS_Pop_CV[, isample]) 
    
    #Create CV dataframe in the formmat of SamplingStrata
    cv <- list()
    for (spp in 1:ns_opt) 
      cv[[paste0("CV", spp)]] <- as.numeric(CV_constraints[spp])
    cv[["DOM"]] <- 1
    cv[["domainvalue"]] <- 1
    cv <- as.data.frame(cv)
  }
}

