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
github_dir <- getwd()

##################################################
####   Load Data
####   Load Population CVs for use in the thresholds
####   Source plotting function
##################################################
load("data/processed/optimization_data.RData")
load('data/processed/dens_vals_fixed_random.RData')
load('data/processed/dens_vals_measurement.RData')
load('data/processed/dens_vals_MLE.RData')

load("data/processed/grid_goa.RData")
source("modified_functions/plot_solution_results.R")

##################################################
####   Constants to specify before doing optimization
##################################################
for (iscen in 1:5) { ## Start loop over scenarios
  
  ## Domain is the term used in the SamplingStrata package, is used to 
  ## distinguish whether the optimization is done gulf-wide (n_dom == 1) or 
  ## at each of the five management districts (n_dom == n_districts)
  which_domain <- scenarios$scale_opt[iscen]
  n_dom <- ifelse(test = which_domain == "full_domain", 
                  yes = 1, 
                  no = n_districts)
  
  
  ## depth input
  depth_input <- grid_goa$DEPTH_EFH
  
  ## For the gulf-wide optimization, use 10 strata
  ## For the district-level optimization, use 5 strata per district
  no_strata <- switch(which_domain,
                      "full_domain" = 10,
                      "district" = rep(5, n_dom))
  
  ## Change depths > 300 m if needed
  if (scenarios$depth_dis[iscen] == 300) 
    depth_input[depth_input > 300] <- 300
  
  ## Subset depths < 700 m if needed
  cell_idx <- rep(x = TRUE, times = n_cells)
  if (scenarios$max_depth[iscen] == 700) 
    cell_idx[grid_goa$DEPTH_EFH > 700] <- FALSE
  
  depth_input <- grid_goa$DEPTH_EFH[cell_idx]
  lon_input <- with(grid_goa[cell_idx, ], E_km - min(E_km))
  
  
  domain_input <- switch(which_domain,
                         "full_domain" = rep(1, sum(cell_idx)),
                         "district" = district_vals[cell_idx])
  
  ## Subset strata variables
  stratum_var_input <- switch(scenarios$stratum_vars[iscen], 
                              "depth_lon" = data.frame(X1 = depth_input,
                                                       X2 = lon_input),
                              "depth" = data.frame(X1 = depth_input))
  
  for (which_species in c(spp_idx_opt, 
                          spp_idx_eval)) { ## Start loop over species
    
    ## Which density values are we using?
    density_input <- 
      get(x = paste0("dens_vals_", 
                     scenarios$data_type[iscen]))[cell_idx, which_species, ]
    
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
      data.frame(domainvalue = domain_input,
                 id = (1:n_cells)[cell_idx],
                 stratum_var_input,
                 WEIGHT = n_years),
      
      matrix(data = apply(X = density_input,
                          MARGIN = 1,
                          FUN = sum),
             ncol = 1,
             dimnames = list(NULL, paste0("Y1"))),
      
      matrix(data = apply(X = density_input,
                          MARGIN = 1,
                          FUN = function(x) sum(x^2)),
             ncol = 1,
             dimnames = list(NULL, paste0("Y", 1, "_SQ_SUM")))
    )
    
    
    for (iboat in 1:n_boats) { ## Start loop over boats
      
      ## Set the result directory to which optimization outputs will save 
      result_dir = paste0(github_dir, "/results/scenario_", 
                          scenarios$scen_name[iscen], 
                          "/Single_Species_Optimization/", 
                          common_names_all[which_species], 
                          '/boat_', iboat, "/")
      if(!dir.exists(result_dir)) dir.create(path = result_dir, recursive = T)
      
      ##################################################
      ####   Run optimization at SRS CV constraints
      ##################################################
      ## Initiate CVs to be those calculated under simple random sampling (SRS)
      srs_stats <- SamplingStrata::buildStrataDF(
        dataset = cbind( frame[, -grep(x = names(frame), pattern = "X")],       
                         X1 = 1))
      
      srs_n <- as.numeric(samples[iboat] * table(frame$domainvalue) / n_cells)
      
      ## SRS statistics
      srs_var <- srs_stats$S1^2 * (1 - srs_n / n_cells) / srs_n
      srs_cv <- sqrt(srs_var) / srs_stats$M1
      
      ## cv is a data input to the SamplingStrata package, assign the initial 
      ## cv constraints 
      cv <- list()
      cv[["CV1"]] <- srs_cv
      cv[["DOM"]] <- 1:n_dom
      cv[["domainvalue"]] <- 1:n_dom
      cv <- as.data.frame(cv)
      
      ## Set wd for output files, create a directory if it doesn"t exist yet
      setwd(result_dir)
      
      #Run optimization, set up a plot layout to show optimization updates
      par(mfrow = c(6, 6), mar = c(2, 2, 0, 0))
      
      solution <- optimStrata(method = "continuous",
                              errors = cv, 
                              framesamp = frame,
                              iter = 300,
                              pops = 100,
                              elitism_rate = 0.1,
                              mut_chance = 1 / (no_strata[1] + 1),
                              nStrata = no_strata,
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
      
      plot_solution <- 
        switch(which_domain,
               "full_domain" = solution$indices$X1,
               "district" = as.factor(paste0(
                 "DOM", solution$framenew$DOMAINVALUE,
                 " STR", solution$framenew$STRATO))
        )
      
      plot_solution <- as.integer(plot_solution)
      
      ##################################################
      ####   Save a plot of the solution
      ##################################################
      temp_ids <- rep(0, n_cells)
      temp_ids[cell_idx] <- plot_solution
      
      plot_solution_results(file_name = paste0("solution.png"),
                            grid_object =  grid_goa,
                            districts_object = districts,
                            district_values = district_vals,
                            sol_by_cell = temp_ids,
                            draw_stations = FALSE)
      
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
      srs_n <- as.numeric(samples[iboat] * table(temp_frame$domainvalue) / n_cells)
      
      ## SRS statistics
      srs_var <- srs_stats$S1^2 * (1 - srs_n / n_cells) / srs_n
      srs_cv <- sqrt(srs_var) / srs_stats$M1
      
      sample_allocations <- matrix(nrow = nrow(sum_stats), 
                                   ncol = 3,
                                   dimnames = list(NULL, paste0("boat_", 1:3)))
      cv_by_boat <- data.frame(total_n = vector(length = 3),
                               cv_constraint = vector(length = 3),
                               actual_cv = vector(length = 3))
      
      error_df <- data.frame("DOM" = "DOM1",
                             "CV1" = srs_cv,
                             "domainvalue"  = 1)
      
      temp_stratif <- solution$aggr_strata
      temp_stratif$N <- temp_stratif$N / n_years
      temp_stratif$DOM1 <- 1
      
      temp_bethel <- SamplingStrata::bethel(
        errors = error_df,
        stratif = temp_stratif, 
        realAllocation = T, 
        printa = T)
      temp_n <- sum(ceiling(temp_bethel))
      
      while (temp_n != samples[iboat]){
        over_under <- temp_n > samples[iboat]
        CV_adj <- ifelse(over_under == TRUE, 
                         yes = 1.001,
                         no = 0.999)
        
        error_df$CV1 <- 
          as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]) * CV_adj
        temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                              errors = error_df, 
                                              printa = TRUE)
        
        temp_n <- sum(as.numeric(temp_bethel))
        
        print(paste0("n = ", temp_n, ", CV = ", 
                     as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"])) )
      }
      
      ##################################################
      ####   Update sample_allocations with optimal allocation
      ##################################################
      sample_allocations[, paste0("boat_", isample)] <- as.numeric(temp_bethel)
      cv_by_boat[isample, "cv_constraint"] <- 
        as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])
      cv_by_boat[isample, "actual_cv"] <- 
        as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"])
      cv_by_boat[isample, "total_n"] <- temp_n
      
      ##################################################
      ####   Plot solution with a random draw of the design
      ##################################################
      plot_solution_results(file_name = paste0("solution_with_stations_boat_",
                                               isample, ".png"),
                            grid_object =  grid_goa,
                            districts_object = districts,
                            district_values = district_vals,
                            sol_by_cell = temp_ids,
                            draw_stations = TRUE, 
                            allocations = sample_allocations[, isample])
      
      ##################################################
      ####   Save output
      ##################################################
      result_list <- list(solution = solution,
                          sum_stats = sum_stats,
                          cvs = cv_by_boat,
                          sample_allocations = sample_allocations,
                          sol_by_cell = temp_ids)
      save(list = "result_list", file = "result_list.RData")
      
    } ## End loop over boats
  } ## End loop over species
} ## End loop over scenarios

