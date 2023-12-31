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
library(terra)
library(RColorBrewer)

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
D_gct <- readRDS("data/processed/VAST_fit_D_gct.RDS")
grid_goa <- readRDS("data/processed/goa_interpolation_grid.RDS")
source("modified_functions/plot_solution_results.R")

##################################################
####   Constants to specify before doing optimization
##################################################
# iscen = 2
for (iscen in 1) { ## Start loop over depth cutoff scenarios
  
  ## Domain is the term used in the SamplingStrata package, is used to 
  ## distinguish whether the optimization is done gulf-wide (n_dom == 1) or 
  ## at each of the five management districts (n_dom == n_districts)
  
  depth_input <- as.numeric(x = as.character(
    x = cut(x = grid_goa$Depth_m, 
            breaks = seq(from = 0, to = 1000, by = 25),
            labels = seq(from = 0, to = 975, by = 25),
            right = F)
  ))
  
  ## Change depths > 300 m if needed
  if (depth_discrete_cutoff[iscen] == 300) 
    depth_input[depth_input > 300] <- 300
  
  for (which_species in 1:ns) { ## Start loop over species
    
    for (strata_per_district in strata) { ## Start loop over strata per area
      
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
      frame <- cbind(data.frame(domainvalue = as.numeric(district_vals),
                                id = 1:n_cells,
                                X1 = depth_input,
                                WEIGHT = n_years),
                     
                     matrix(data = apply(X = D_gct[, which_species, ],
                                         MARGIN = 1,
                                         FUN = sum),
                            ncol = 1,
                            dimnames = list(NULL, paste0("Y1"))),
                     
                     matrix(data = apply(X = D_gct[, which_species, ],
                                         MARGIN = 1,
                                         FUN = function(x) sum(x^2)),
                            ncol = 1,
                            dimnames = list(NULL, paste0("Y", 1, "_SQ_SUM"))))
      
      ## Set the result directory to which optimization outputs will save 
      result_dir = paste0(github_dir, "/results/depthcut_", 
                          depth_discrete_cutoff[iscen], "m/",
                          "Str_", strata_per_district, "/",
                          "Single_Species_Optimization/", 
                          common_names[which_species], "/")
      if(!dir.exists(paths = result_dir)) 
        dir.create(path = result_dir, recursive = T)
      
      ##################################################
      ####   Run optimization at SRS CV constraints
      ##################################################
      ## Initiate CVs to be those calculated under simple random sampling (SRS)
      srs_stats <- SamplingStrata::buildStrataDF(
        dataset = cbind( frame[, -grep(x = names(frame), pattern = "X")],       
                         X1 = 1))
      
      srs_n <- as.numeric(target_n * table(frame$domainvalue) / n_cells)
      
      ## SRS statistics
      srs_var <- srs_stats$S1^2 * (1 - srs_n / n_cells) / srs_n
      srs_cv <- sqrt(srs_var) / srs_stats$M1
      
      ## cv is a data input to the SamplingStrata package, assign the initial 
      ## cv constraints 
      cv <- list()
      cv[["CV1"]] <- srs_cv
      cv[["DOM"]] <- 1:n_districts
      cv[["domainvalue"]] <- 1:n_districts
      cv <- as.data.frame(x = cv)
      
      ## Set wd for output files, create a directory if it doesn"t exist yet
      setwd(result_dir)
      
      solution <- optimStrata(method = "continuous",
                              errors = cv, 
                              framesamp = frame,
                              # iter = 10,
                              # pops = 10,
                              iter = 300,
                              pops = 100,
                              elitism_rate = 0.1,
                              mut_chance = 1 / (strata_per_district + 1),
                              nStrata = rep(strata_per_district, n_districts),
                              showPlot = F,
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
      srs_n <- as.numeric(target_n * table(temp_frame$domainvalue) / n_cells)
      
      ## SRS statistics
      srs_var <- srs_stats$S1^2 * (1 - srs_n / n_cells) / srs_n
      srs_cv <- sqrt(srs_var) / srs_stats$M1
      
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
      
      iter = 1
      while (temp_n != target_n & iter < 10000){
        over_under <- temp_n > target_n
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
        iter = iter + 1
      }
      
      ##################################################
      ####   Update sample_allocations with optimal allocation
      ##################################################
      cvs <- c("planned_cv" = NA, "actual_cv" = NA)
      sample_allocations <- as.numeric(temp_bethel)
      cvs["planned_cv"] <- 
        as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])
      cvs["actual_cv"] <- 
        as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"])
      
      ##################################################
      ####   Save output
      ##################################################
      result_list <- list(solution = solution,
                          sum_stats = sum_stats,
                          cvs = cvs,
                          sample_allocations = sample_allocations,
                          sol_by_cell = plot_solution)
      save(list = "result_list", file = "result_list.RData")
      
    } ## End loop over strata per area
  } ## End loop over species
} ## End loop over depth cutoff scenarios
