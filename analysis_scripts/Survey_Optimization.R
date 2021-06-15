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
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[2]

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
####   Create optimization scenarios
##################################################
scen <- data.frame(nstrata = c(3,5, 10,15),
                   which_domain = rep(c("district", "full_domain"), each = 2))

##################################################
####   Collect optimization results from each strata
##################################################
for (irow in 1:nrow(scen)) { ## Loop through scen dataframe -- start
  for(isample in 1:n_boats) {
    
    ##  Constants to specify before doing optimization
    ##  The optimization at the full_domain (gulf-wide) uses a different
    ##      data input than the optimization at the district-level
    which_domain <- scen$which_domain[irow]
    
    frame <- switch( which_domain,
                     "full_domain" = frame_all,
                     "district" = frame_district)[, c("domainvalue", "id", 
                                                      "X1", "X2", "WEIGHT",
                                                      paste0("Y", spp_idx_opt), 
                                                      paste0("Y", spp_idx_opt,
                                                             "_SQ_SUM"))]

    ## Only subset the species included in the optimization, indexed by 
    ## spp_idx_opt, then reorganize field names
    names(frame)[names(frame) %in% paste0("Y", spp_idx_opt)] <- 
      paste0("Y", 1:ns_opt)
    names(frame)[names(frame) %in% paste0("Y", spp_idx_opt, "_SQ_SUM")] <- 
      paste0("Y", 1:ns_opt, "_SQ_SUM")
    
    ## n_dom: domain is the term used in the optimization
    ##        1 domain if conducting a gulf-wide optimization
    ##        5 domains if conduct a district-level optimization (i.e., 
    ##                 five managment districts)
    ##
    ## Assign the number of strata in each domain
    n_dom <- length(unique(frame$domainvalue))
    temp_strata <- rep(x = scen$nstrata[irow], times = n_dom)
    
    ## Initial Condition
    run <- 1 
    current_n <- 0
    
    ## Initiate CVs to be those calculated under SRS, assign to a variable 
    ## named cv_constraints
    ## buildStrataDF calculates the stratum means and variances, X1 = 1 
    ##     means to calculate those statics on the whole domain
    srs_stats <- SamplingStrata::buildStrataDF( 
      dataset = cbind( subset(frame, select = -c(X1, X2)),
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
    
    ## Load the single species optimized CVs
    load(paste0(github_dir, "results/", which_domain, 
                "/Single_Species_Optimization/",
                "optimization_knitted_results.RData"))
    
    ss_strs_pop_cv <- subset(x = settings,
                             subset = boat == isample,
                             select = c("species", paste0("cv_domain_", 
                                                          1:n_dom)))
    
    ss_strs_pop_cv <- ss_strs_pop_cv[match(common_names_opt, 
                                           ss_strs_pop_cv$species), ]
    ss_strs_pop_cv <- t(ss_strs_pop_cv[, -1])
    colnames(ss_strs_pop_cv) <- common_names_opt
    
    ## Run optimization
    while (current_n <= samples[isample] ) { 
      
      #Set wd for output files, create a directory if it doesn"t exist yet
      temp_dir = paste0(github_dir, "results/", which_domain, 
                        "/Multi_Species_Optimization/boat", isample,
                        "/Str", temp_strata[1], "/Run", run)
      
      if(!dir.exists(temp_dir)) dir.create(temp_dir, recursive = T)
      
      setwd(temp_dir)
      
      #Run optimization
      if(which_domain == "full_domain") par(mfrow = c(6,6), 
                                            mar = c(2,2,0,0))
      
      solution <- optimStrata(method = "continuous",
                              errors = cv, 
                              framesamp = frame,
                              iter = 300,
                              pops = 50,
                              elitism_rate = 0.1,
                              mut_chance = 1 / (temp_strata[1] + 1),
                              nStrata = temp_strata,
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
      sum_stats <- sum_stats[, c(10, 1:4, 15, 11:12, 6:9)]
      
      plot_solution <- 
        switch(which_domain,
               "full_domain" = solution$indices$X1,
               "district" = as.factor(paste0(
                 "DOM", solution$framenew$DOMAINVALUE,
                 " STR", solution$framenew$STRATO))
        )
      
      plot_solution <- as.integer(plot_solution)
      
      ## Save Output
      CV_constraints <- expected_CV(strata = solution$aggr_strata)
      current_n <- sum(sum_stats$Allocation)
      result_list <- list(solution = solution, 
                          sum_stats = sum_stats, 
                          CV_constraints = CV_constraints, 
                          n = current_n,
                          sol_by_cell = plot_solution)
      save(list = "result_list", file = "result_list.RData")
      
      ##Save a plot of the solution
      goa <- sp::SpatialPointsDataFrame(
        coords = Extrapolation_depths[, c("E_km", "N_km")],
        data = data.frame(Str_no = plot_solution) )
      goa_ras <- raster::raster(x = goa, 
                                resolution = 5)
      goa_ras <- raster::rasterize(x = goa, 
                                   y = goa_ras, 
                                   field = "Str_no")
      
      png(filename = "solution.png",
          width = 6,
          height = 3,
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
      
      rect(xleft = districts$W_UTM,
           xright = districts$E_UTM,
           ybottom = tapply(X = Extrapolation_depths$N_km, 
                            INDEX = district_vals,
                            FUN = min), 
           ytop = tapply(X = Extrapolation_depths$N_km, 
                         INDEX = district_vals,
                         FUN = max))
      
      text(x = rowMeans(districts[, c("W_UTM", "E_UTM")]),
           y = tapply(X = Extrapolation_depths$N_km, 
                      INDEX = district_vals,
                      FUN = max),
           labels = districts$district,
           pos = 3)
      box()
      dev.off()
      
      
      png(filename = "solution_with_stations.png",
          width = 6,
          height = 3,
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
      
      rect(xleft = districts$W_UTM,
           xright = districts$E_UTM,
           ybottom = tapply(X = Extrapolation_depths$N_km, 
                            INDEX = district_vals,
                            FUN = min), 
           ytop = tapply(X = Extrapolation_depths$N_km, 
                         INDEX = district_vals,
                         FUN = max))
      
      text(x = rowMeans(districts[, c("W_UTM", "E_UTM")]),
           y = tapply(X = Extrapolation_depths$N_km, 
                      INDEX = district_vals,
                      FUN = max),
           labels = districts$district,
           pos = 3)
      box()
      
      #Take a random sample based on the allocation and stratum
      sample_vec <- c()
      for(istrata in 1:nrow(sum_stats)) {
        sample_vec <- c(sample_vec,
                        sample(x = which(plot_solution == istrata),
                               size = sum_stats$Allocation[istrata]) )
      }
      
      points(Extrapolation_depths[sample_vec, c("E_km", "N_km")],
             pch = 16, cex = 0.5)
      
      dev.off()
      
      
      #Set up next run by changing upper CV constraints
      run <- run + 1
      
      cv_constraints <- 0.95 * cv_constraints + 0.05 * ss_strs_pop_cv
      
      #Create CV dataframe in the formmat of SamplingStrata
      cv <- list()
      for (spp in 1:ns_opt) 
        cv[[paste0("CV", spp)]] <- as.numeric(cv_constraints[, spp])
      cv[["DOM"]] <- 1:n_dom
      cv[["domainvalue"]] <- 1:n_dom
      cv <- as.data.frame(cv)
    }
  }
}  ## Loop through scen dataframe -- start
