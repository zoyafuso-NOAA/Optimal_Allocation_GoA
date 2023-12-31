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
# library(sp)
library(RColorBrewer)
# library(raster)

github_dir <- getwd()

##################################################
####   Load Data
####   Load Population CVs for use in the thresholds
####   Source plotting function
##################################################
load("data/processed/optimization_data.RData")
D_gct <- readRDS("data/processed/VAST_fit_D_gct.RDS")
grid_goa <- readRDS("data/processed/goa_interpolation_grid.RDS")
source("modified_functions/plot_solution_results.R")

##################################################
####   Collect optimization results from each strata
##################################################

for (iscen in 1:2){
  
  depth_input <- as.numeric(x = as.character(
    x = cut(x = grid_goa$Depth_m, 
            breaks = seq(from = 0, to = 1000, by = 25),
            labels = seq(from = 0, to = 975, by = 25),
            right = F)
  ))
  
  ## Change depths > 300 m if needed
  if (depth_discrete_cutoff[iscen] == 300) 
    depth_input[depth_input > 300] <- 300
  
  ##################################################
  ####   Our df will have fields for:
  ####   domain: only one domain so the value is just 1
  ####   id: unique ID for each sampling cell
  ####   X1: strata variable 2: depth of cell (m) 
  ####   X2: strata variable 1: longitude in eastings (km). Because the 
  ####       optimization does not read in negative values, I shift the 
  ####       values so that the lowest value is 0
  ####
  ####   Variables used to more efficiently calcualte stratum variance 
  ####
  ####   WEIGHT: number of observed years 
  ####   Y1, Y2, ... : density for a given cell summed across observed years
  ####   Y1_SQ_SUM, Y2_SQ_SUM, ... : density-squared for a given cell, 
  ####           summed across observed years
  ##################################################
  frame <- cbind(
    data.frame(domainvalue = as.numeric(district_vals),
               id = 1:n_cells,
               X1 = depth_input,
               WEIGHT = n_years),
    
    matrix(data = apply(X = D_gct,
                        MARGIN = 1:2,
                        FUN = sum),
           ncol = ns,
           dimnames = list(NULL, paste0("Y", 1:ns))),
    
    matrix(data = apply(X = D_gct,
                        MARGIN = 1:2,
                        FUN = function(x) sum(x^2)),
           ncol = ns,
           dimnames = list(NULL, paste0("Y", 1:ns, "_SQ_SUM")))
  )
  
  for (strata_per_district in strata) { ## Start loop over strata per area
    ## Compile Single-Species CVs to establish a lower limit for a given boat
    ss_cvs <- vector(length = ns); names(ss_cvs) <- common_names
    for (ispp in common_names) {
      load( paste0(github_dir, "/results/depthcut_", 
                   depth_discrete_cutoff[iscen], "m/Str_", strata_per_district, 
                   "/Single_Species_Optimization/", ispp, "/result_list.RData"))
      ss_cvs[ispp] <- result_list$cvs["actual_cv"]
    }
    
    ## Initiate CVs to be those calculated under SRS, assign to a variable 
    ## named cv_constraints
    ## buildStrataDF calculates the stratum means and variances, X1 = 1 
    ##     means to calculate those statics on the whole domain
    srs_stats <- SamplingStrata::buildStrataDF(
      dataset = cbind( frame[, -grep(x = names(frame), pattern = "X")],       
                       X1 = 1))
    
    srs_n <- as.numeric(target_n * table(frame$domainvalue) / nrow(frame))
    srs_var <- as.matrix(srs_stats[, paste0("S", 1:ns)])^2
    
    srs_var <- sweep(x = srs_var, 
                     MARGIN = 1, 
                     STATS = (1 - srs_n / n_cells) / srs_n, 
                     FUN = "*")
    
    srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:ns)]
    
    cv_constraints <- srs_cv
    
    ## Create CV constraint df
    cv <- list()
    for (spp in 1:ns) 
      cv[[paste0("CV", spp)]] <- as.numeric(cv_constraints[, spp])
    cv[["DOM"]] <- 1:n_districts
    cv[["domainvalue"]] <- 1:n_districts
    cv <- as.data.frame(cv)
    
    ## Set the result directory to which optimization outputs will save 
    result_dir = paste0(github_dir, "/results/depthcut_", 
                        depth_discrete_cutoff[iscen], "m/",
                        "Str_", strata_per_district, "/",
                        "Multi_Species_Optimization/")
    if(!dir.exists(result_dir)) dir.create(path = result_dir, recursive = T)
    setwd(result_dir)
    
    solution <- SamplingStrata::optimStrata(
      method = "continuous",
      errors = cv, 
      framesamp = frame,
      iter = 300,
      pops = 100,
      elitism_rate = 0.1,
      mut_chance = 1 / (rep(strata_per_district, n_districts) + 1),
      nStrata = rep(strata_per_district, n_districts),
      showPlot = T,
      writeFiles = T)
    
    
    ## Organize result outputs
    solution$aggr_strata$STRATO <- as.integer(solution$aggr_strata$STRATO)
    solution$aggr_strata <- 
      solution$aggr_strata[order(solution$aggr_strata$DOM1,
                                 solution$aggr_strata$STRATO), ]
    
    sum_stats <- summaryStrata(solution$framenew,
                               solution$aggr_strata,
                               progress=FALSE)
    sum_stats$stratum_id <- 1:nrow(sum_stats)
    sum_stats$Population <- sum_stats$Population / n_years
    sum_stats$wh <- sum_stats$Allocation / sum_stats$Population
    sum_stats$Wh <- sum_stats$Population / n_cells
    sum_stats <- cbind(sum_stats,
                       subset(x = solution$aggr_strata,
                              select = -c(STRATO, N, COST, CENS, DOM1, X1)))
    
    ##Plot Solution
    plot_solution <- as.factor(paste0(
      "DOM", solution$framenew$DOMAINVALUE,
      " STR", solution$framenew$STRATO))
    plot_solution <- as.integer(plot_solution)
    
    ##################################################
    ####   Tune CV to hit 1, 2 or 3 boats (292, 550, 825 stations) 
    ####   Assume CVs at the gulf-scale regardless on whether scale of the 
    ####   optimization
    ##################################################
    temp_frame <- frame
    temp_frame$domainvalue <- 1
    srs_stats <- SamplingStrata::buildStrataDF(
      dataset = cbind( temp_frame[, -grep(x = names(temp_frame), 
                                          pattern = "X")],
                       X1 = 1))
    
    srs_n <- sum(as.numeric(target_n * table(frame$domainvalue) / nrow(x = frame)))
    srs_var <- as.matrix(srs_stats[, paste0("S", 1:ns)])^2
    
    ## SRS statistics
    srs_var <- sweep(x = srs_var, 
                     MARGIN = 1, 
                     STATS = (1 - srs_n / nrow(frame)) / srs_n, 
                     FUN = "*")
    srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:ns)]
    
    sample_allocations <- c()
    cv_by_boat <- list()
    
    error_df <- data.frame("DOM" = "DOM1",
                           srs_cv,
                           "domainvalue"  = 1)
    names(error_df)[2:(1 + ns)] <- paste0("CV", 1:ns)
    
    temp_stratif <- solution$aggr_strata
    temp_stratif$N <- temp_stratif$N / n_years
    temp_stratif$DOM1 <- 1
    
    temp_bethel <- SamplingStrata::bethel(
      errors = error_df,
      stratif = temp_stratif, 
      realAllocation = T, 
      printa = T)
    temp_n <- sum(ceiling(temp_bethel))
    
    iter = 1
    while (temp_n != target_n & iter < 10000){
      over_under <- temp_n > target_n
      CV_adj <- ifelse(over_under == TRUE, 
                       yes = 1.001,
                       no = 0.999)
      
      updated_cv_constraint <- 
        as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]) * (CV_adj) + 
        ss_cvs * (1  - CV_adj)
      
      error_df[, paste0("CV", 1:ns)] <- as.numeric(updated_cv_constraint)
      
      temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                            errors = error_df, 
                                            printa = TRUE)
      
      temp_n <- sum(as.numeric(temp_bethel))
      
      print(paste0("n = ", temp_n) )
      iter = iter + 1
    }
    
    ##################################################
    ####   Updated nh
    ##################################################
    sample_allocations <- as.numeric(temp_bethel)
    cvs <-
      data.frame(species = common_names,
                 cv_constraint = as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]),
                 actual_cv =     as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"]))
    
    ##################################################
    ####   Save output
    ##################################################
    result_list <- list(solution = solution,
                        sum_stats = sum_stats,
                        cvs = cvs,
                        sample_allocations = sample_allocations,
                        sol_by_cell = plot_solution)
    save(list = "result_list", file = "result_list.RData")
    
  } ## Start loop over strata per area
}
