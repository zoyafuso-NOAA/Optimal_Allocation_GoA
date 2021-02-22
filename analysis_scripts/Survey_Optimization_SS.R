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
which_domain <- c("full_domain", "district")[1]

for (which_species in 1) {
  
  ##################################################
  ####   Constants to set up based on which_domain and which_species
  ##################################################
  frame <- switch( which_domain,
                   "full_domain" = frame_all,
                   "district" = frame_district)[, c("domainvalue", "id", 
                                                    "X1", "X2", "WEIGHT",
                                                    paste0("Y", which_species), 
                                                    paste0("Y", which_species,
                                                           "_SQ_SUM"))]
  names(frame)[6:7] <- paste0("Y", c("1", "1_SQ_SUM") )
  n_dom <- length(unique(frame$domainvalue))
  no_strata <- switch(which_domain,
                      "full_domain" = 10,
                      "district" = rep(5, n_dom))
  
  result_dir = paste0(github_dir, 
                      "results/", which_domain, 
                      "/Single_Species_Optimization/", 
                      gsub(x = sci_names_all[which_species], 
                           pattern = ' ', 
                           replacement = '_'), '/')
  if(!dir.exists(result_dir)) dir.create(path = result_dir, recursive = T)
  
  ##################################################
  ####   Run optimization
  ##################################################
  
  ##Initial Conditions
  run <- 1
  current_n <- 0
  
  ## Initiate CVs to be those calculated under SRS
  srs_stats <- SamplingStrata::buildStrataDF( 
    dataset = cbind( subset(frame, select = -c(X1, X2)),
                     X1 = 1))
  
  srs_n <- as.numeric(280 * table(frame$domainvalue) / n_cells)
  
  srs_var <- srs_stats$S1^2 * (1 - srs_n / n_cells) / srs_n
  srs_cv <- sqrt(srs_var) / srs_stats$M1
  
  cv <- list()
  cv[["CV1"]] <- srs_cv
  cv[["DOM"]] <- 1:n_dom
  cv[["domainvalue"]] <- 1:n_dom
  cv <- as.data.frame(cv)
  
  while (current_n <= 820 ) {
    
    #Set wd for output files, create a directory if it doesn"t exist yet
    temp_dir = paste0(result_dir, "Run", run)
    if(!dir.exists(temp_dir)) dir.create(temp_dir, recursive = T)
    
    setwd(temp_dir)
    
    #Run optimization
    par(mfrow = c(6,6), 
        mar = c(2,2,0,0))
    solution <- optimStrata(method = "continuous",
                            errors = cv, 
                            framesamp = frame,
                            iter = 300,
                            pops = 50,
                            elitism_rate = 0.1,
                            mut_chance = 1 / (no_strata[1] + 1),
                            nStrata = no_strata,
                            showPlot = T,
                            writeFiles = T)
    
    sum_stats <- summaryStrata(solution$framenew,
                               solution$aggr_strata,
                               progress=FALSE) 
    
    plot_solution <- as.factor(paste(solution$framenew$DOMAINVALUE,
                                     solution$framenew$STRATO))
    plot_solution <- as.integer(plot_solution)
    
    ##Save a plot of the solution
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
    
    par(mfrow = c(1, 1), 
        mar = c(1, 1, 1, 1))
    plot(goa_ras, 
         axes = F, 
         asp = 1,
         col = colorRampPalette(
           brewer.pal(n = 11, 
                      name = "Paired"))(length(unique(plot_solution)) ) )
    
    rect(xleft = districts$W_lon,
         xright = districts$E_lon,
         ybottom = tapply(X = Extrapolation_depths$Lat, 
                          INDEX = district_vals,
                          FUN = min), 
         ytop = tapply(X = Extrapolation_depths$Lat, 
                       INDEX = district_vals,
                       FUN = max))
    
    text(x = rowMeans(districts[, c("W_lon", "E_lon")]),
         y = tapply(X = Extrapolation_depths$Lat, 
                    INDEX = district_vals,
                    FUN = max),
         labels = districts$district,
         pos = 3)
    
    box()
    dev.off()
    
    ## Save Output
    CV_constraints <- expected_CV(strata = solution$aggr_strata)
    current_n <- sum(sum_stats$Allocation)
    
    result_list <- list(solution = solution, 
                        sum_stats = sum_stats, 
                        CV_constraints = CV_constraints, 
                        n = current_n)
    save(list = "result_list", file = "result_list.RData")
    
    ## Set up next run by changing slightly reducing the CV constraints
    ## CVs are reduced proportionally, based on the effort level
    run <- run + 1
    effort_level <- as.integer(cut(x = current_n, 
                                   breaks = c(0, 200, samples, 1000), 
                                   labels = 1:5))
    CV_constraints <- CV_constraints * c(0.80, 0.90, 0.95, 0.975)[effort_level]
    
    #Create CV dataframe in the format of SamplingStrata
    cv <- list()
    cv[["CV1"]] <- as.numeric(CV_constraints)
    cv[["DOM"]] <- 1:n_dom
    cv[["domainvalue"]] <- 1:n_dom
    cv <- as.data.frame(cv)
  }
}
