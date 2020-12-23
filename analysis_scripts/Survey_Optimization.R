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
scen <- data.frame(nstrata = c(3,5,10, 10,15,20),
                   which_domain = rep(c("district", "full_domain"), each = 3))

##################################################
####   Collect optimization results from each strata
##################################################
for (irow in 5) {
  for(isample in 1:n_boats) {
    
    ##################################################
    ####   Constants to specify before doing optimization
    ##################################################
    which_domain <- scen$which_domain[irow]
    
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
    
    n_dom <- length(unique(frame$domainvalue))
    
    temp_strata <- rep(x = scen$nstrata[irow], times = n_dom)
    
    ##Initial Condition
    run <- 1
    current_n <- 0
    
    ## Load SRS information to initialize the starting points for the CVs
    load(paste0(github_dir, "results/", which_domain, "/srs_pop_cv.RData"))
    
    cv_constraints <- get(paste0("srs_pop_cv_", which_domain))[, , isample] 
    cv_constraints <- switch(which_domain, 
                             "district" = cv_constraints[spp_idx_opt, ],
                             "full_domain" = cv_constraints[spp_idx_opt])
    
    cv <- list()
    for (spp in 1:ns_opt) 
      cv[[paste0("CV", spp)]] <- 
      as.numeric(switch(which_domain, 
                        "district" = cv_constraints[spp, ],
                        "full_domain" = cv_constraints[spp]))
    cv[["DOM"]] <- 1:n_dom
    cv[["domainvalue"]] <- 1:n_dom
    cv <- as.data.frame(cv)
    
    load(paste0(github_dir, "results/", which_domain, 
                "/Single_Species_Optimization/",
                "optimization_knitted_results.RData"))
    
    ss_strs_pop_cv <- 
      switch(which_domain,
             "district" = t(subset(x = settings_district,
                                   subset = iboat == isample & 
                                     spp %in% spp_idx_opt,
                                   select = paste(1:5))),
             "full_domain" = unlist(subset(x = settings_agg_full_domain,
                                           subset = iboat == isample & 
                                             spp %in% spp_idx_opt,
                                           select = cv)))
    
    ##################################################
    ####   Run optimization
    ##################################################
    while (current_n <= c(280, 550, 820)[isample] ) {
      
      #Set wd for output files, create a directory if it doesn"t exist yet
      temp_dir = paste0(github_dir, "results/", which_domain, 
                        "/Multi_Species_Optimization/boat", isample,
                        "/Str", temp_strata[1], "/Run", run)
      
      if(!dir.exists(temp_dir)) dir.create(temp_dir, recursive = T)
      
      setwd(temp_dir)
      
      #Run optimization
      if(which_domain == "full_domain") par(mfrow = c(6,6), mar = c(2,2,0,0))
      
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
      
      sum_stats <- summaryStrata(solution$framenew,
                                 solution$aggr_strata,
                                 progress=FALSE) 
      
      #Plot Solution
      plot_solution <- as.factor(paste(solution$framenew$DOMAINVALUE,
                                       solution$framenew$STRATO))
      plot_solution <- as.integer(plot_solution)
      
      goa <- sp::SpatialPointsDataFrame(
        coords = Extrapolation_depths[,c("Lon", "Lat")],
        data = data.frame(Str_no = plot_solution) )
      goa_ras <- raster::raster(x = goa, 
                                resolution = 0.075)
      goa_ras <- raster::rasterize(x = goa, 
                                   y = goa_ras, 
                                   field = "Str_no")
      
      png(filename = "solution.png", 
          width = 5, 
          height = 3, 
          units = "in", 
          res = 500)
      
      par(mfrow = c(1, 1), 
          mar = c(1, 1, 1, 1))
      plot( goa_ras, 
            axes = F, 
            asp = 1,
            col = colorRampPalette(
              brewer.pal(n = 11, 
                         name = "Paired"))(sum(temp_strata))[sample(1:sum(temp_strata))] ) 
      
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
      
      png(filename = "solution_with_simulated_survey.png", 
          width = 5, 
          height = 3, 
          units = "in", 
          res = 500)
      
      par(mfrow = c(1, 1), 
          mar = c(1, 1, 1, 1))
      plot( goa_ras, 
            axes = F, 
            asp = 1,
            col = colorRampPalette(
              brewer.pal(n = 11, 
                         name = "Paired"))(sum(temp_strata))[sample(1:sum(temp_strata))] ) 
      
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
      
      #Simulate a sample solution
      temp_samples <- c()
      temp_strata <- nrow(sum_stats)
      temp_solution <- solution$framenew$STRATO
      temp_allocation <- sum_stats$Allocation
      
      for (temp_istrata in 1:temp_strata) {
        temp_samples = c(temp_samples,
                         sample(x = which(temp_solution == temp_istrata),
                                size = temp_allocation[temp_istrata]) )
      }
      
      temp_loc <- Extrapolation_depths[temp_samples, c("Lon", "Lat")]
      
      points(temp_loc,
             pch = 16,
             cex = 0.25)
      
      
      dev.off()
      
      #Save Output
      cv_constraints <- expected_CV(strata = solution$aggr_strata)
      current_n <- sum(sum_stats$Allocation)
      
      result_list <- list(solution = solution, 
                          sum_stats = sum_stats, 
                          cv_constraints = cv_constraints, 
                          n = current_n)
      save(list = "result_list", file = "result_list.RData")
      
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
}
