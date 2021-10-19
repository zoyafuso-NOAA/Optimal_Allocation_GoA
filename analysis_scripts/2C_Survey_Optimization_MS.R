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
####   Collect optimization results from each strata
##################################################
for (iscen in c(11)){
  which_domain <- scenarios$scale_opt[iscen]
  n_dom <- ifelse(test = which_domain == "full_domain", 
                  yes = 1, 
                  no = n_districts)
  
  ## Subset depths < 700 m if needed
  cell_idx <- rep(x = TRUE, times = n_cells)
  if (scenarios$max_depth[iscen] == 700) cell_idx[grid_goa$DEPTH_EFH > 700] <- FALSE
  
  domain_input <- switch(which_domain,
                         "full_domain" = rep(1, sum(cell_idx)),
                         "district" = district_vals[cell_idx])
  
  ## Change depths > 300 m if needed
  depth_input <- grid_goa$DEPTH_EFH[cell_idx]
  if (scenarios$depth_dis[iscen] == 300) depth_input[depth_input > 300] <- 300
  
  lon_input <- with(grid_goa[cell_idx, ], E_km - min(E_km))
  
  ## Subset strata variables
  
  stratum_var_input <- switch(scenarios$stratum_vars[iscen], 
                              "depth_lon" = data.frame(X1 = depth_input,
                                                       X2 = lon_input),
                              "depth" = data.frame(X1 = depth_input))
  
  ## For the gulf-wide optimization, use 10 or 15 strata
  ## For the district-level optimization, use 3 or 5 strata per district
  no_strata <- switch(which_domain,
                      "full_domain" = c(10, 15),
                      "district" = c(3, 5))
  
  ## Which density values are we using?
  density_input <- 
    get(x = paste0("dens_vals_", 
                   scenarios$data_type[iscen]))[cell_idx, spp_idx_opt, ]
  
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
                        MARGIN = 1:2,
                        FUN = sum),
           ncol = ns_opt,
           dimnames = list(NULL, paste0("Y", 1:ns_opt))),
    
    matrix(data = apply(X = density_input,
                        MARGIN = 1:2,
                        FUN = function(x) sum(x^2)),
           ncol = ns_opt,
           dimnames = list(NULL, paste0("Y", 1:ns_opt, "_SQ_SUM")))
  )
  
  
  ## Compile Single-Species CVs to establish a lower limit
  ss_cvs <- vector(length = ns_opt); names(ss_cvs) <- common_names_opt
  for (ispp in common_names_opt){
    load( paste0(github_dir, "/results/scenario_", scenarios$scen_name[iscen], 
                 "/Single_Species_Optimization/", ispp, "/result_list.RData"))
    ss_cvs[ispp] <- result_list$cvs[1, "actual_cv"]
  }
  
  for (isample in 1:3) {
    
    ## Initiate CVs to be those calculated under SRS, assign to a variable 
    ## named cv_constraints
    ## buildStrataDF calculates the stratum means and variances, X1 = 1 
    ##     means to calculate those statics on the whole domain
    srs_stats <- SamplingStrata::buildStrataDF(
      dataset = cbind( frame[, -grep(x = names(frame), pattern = "X")],       
                       X1 = 1))
    
    srs_n <- as.numeric(samples[isample] * table(frame$domainvalue) / n_cells)
    srs_var <- as.matrix(srs_stats[, paste0("S", 1:ns_opt)])^2
    
    srs_var <- sweep(x = srs_var, 
                     MARGIN = 1, 
                     STATS = (1 - srs_n / n_cells) / srs_n, 
                     FUN = "*")
    
    srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:ns_opt)]
    
    cv_constraints <- srs_cv
    
    ## Create CV constraint df
    cv <- list()
    for (spp in 1:ns_opt) 
      cv[[paste0("CV", spp)]] <- 
      as.numeric(switch(which_domain, 
                        "district" = cv_constraints[, spp],
                        "full_domain" = cv_constraints[spp]))
    cv[["DOM"]] <- 1:n_dom
    cv[["domainvalue"]] <- 1:n_dom
    cv <- as.data.frame(cv)
    
    for(temp_strata in no_strata) {
      
      ## Set the result directory to which optimization outputs will save 
      result_dir = paste0(github_dir, "/results/scenario_", 
                          scenarios$scen_name[iscen], 
                          "/Multispecies_Optimization/Str_", temp_strata, 
                          '/boat_', iboat, "/")
      if(!dir.exists(result_dir)) dir.create(path = result_dir, recursive = T)
      setwd(result_dir)
      
      #Run optimization
      if(which_domain == "full_domain") par(mfrow = c(6,6), 
                                            mar = c(2,2,0,0))
      
      solution <- optimStrata(method = "continuous",
                              errors = cv, 
                              framesamp = frame,
                              iter = 300,
                              pops = 100,
                              elitism_rate = 0.1,
                              mut_chance = 1 / (rep(temp_strata, n_dom) + 1),
                              nStrata = rep(temp_strata, n_dom),
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
      
      srs_n <- as.numeric(samples[isample] * table(temp_frame$domainvalue) / n_cells)
      srs_var <- as.matrix(srs_stats[, paste0("S", 1:ns_opt)])^2
      
      ## SRS statistics
      srs_var <- sweep(x = srs_var, 
                       MARGIN = 1, 
                       STATS = (1 - srs_n / n_cells) / srs_n, 
                       FUN = "*")
      srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:ns_opt)]
      
      sample_allocations <- matrix(nrow = nrow(sum_stats), 
                                   ncol = 3,
                                   dimnames = list(NULL, paste0("boat_", 1:3)))
      cv_by_boat <- list()
      
      error_df <- data.frame("DOM" = "DOM1",
                             srs_cv,
                             "domainvalue"  = 1)
      names(error_df)[2:(1 + ns_opt)] <- paste0("CV", 1:ns_opt)
      
      temp_stratif <- solution$aggr_strata
      temp_stratif$N <- temp_stratif$N / n_years
      temp_stratif$DOM1 <- 1
      
      temp_bethel <- SamplingStrata::bethel(
        errors = error_df,
        stratif = temp_stratif, 
        realAllocation = T, 
        printa = T)
      temp_n <- sum(ceiling(temp_bethel))
      
      
      while (temp_n != c(292, 550, 825)[isample]){
        over_under <- temp_n > c(292, 550, 825)[isample]
        CV_adj <- ifelse(over_under == TRUE, 
                         yes = 1.001,
                         no = 0.999)
        
        updated_cv_constraint <- 
          as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]) * (CV_adj) + 
          ss_cvs * (1  - CV_adj)
        
        error_df[, paste0("CV", 1:ns_opt)] <- as.numeric(updated_cv_constraint)
        
        temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                              errors = error_df, 
                                              printa = TRUE)
        
        temp_n <- sum(as.numeric(temp_bethel))
        
        print(paste0("n = ", temp_n) )
      }
      
      ##################################################
      ####   Updated nh
      ##################################################
      sample_allocations[, paste0("boat_", isample)] <- as.numeric(temp_bethel)
      cv_by_boat[[paste0("boat_", isample)]] <-
        data.frame(species = common_names_opt,
                   cv_constraint = as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]),
                   actual_cv =     as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"]))
      
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
    }
    
    ##################################################
    ####   Save output
    ##################################################
    result_list <- list(solution = solution,
                        sum_stats = sum_stats,
                        cvs = cv_by_boat,
                        sample_allocations = sample_allocations,
                        sol_by_cell = temp_ids)
    save(list = "result_list", file = "result_list.RData")
  }
}
