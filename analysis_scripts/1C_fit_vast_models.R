###############################################################################
## Project:      Univariate VAST model runs
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Contributors: Lewis Barnett (lewis.barnett@noaa.gov)
##               Jim Thorson"s VAST wiki example
##           (https://github.com/James-Thorson-NOAA/VAST/wiki/Crossvalidation)
## Description:  Run single-species VAST models wit and without depth
##               as a covariate. Run 10-fold Cross Validation for each Model
##
##               Software versions:
##               R version 4.0.3
##               VAST version 3.6.1
##               FishStatsUtils 2.8.0
###############################################################################
rm(list = ls())

##################################################
####   Model output is saved outside of the repository due to the immense 
####   size of the result outputs
##################################################
VAST_dir <- c("E:/VAST_Runs/")
if(!dir.exists(VAST_dir)) dir.create(VAST_dir, recursive = T)

temp_res <- paste0(getwd(), "/results/temp_res/")
if(!dir.exists(temp_res)) dir.create(temp_res, recursive = T)

##################################################
####   Load packages, make sure versions are consistent
##################################################
library(VAST)

{
  switch(EXPR = R.version$version.string == "R version 4.0.2 (2020-06-22)",
         "TRUE" = print("R version is consistent with R version 4.0.2"),
         "FALSE" = print(
           paste0("R version ", 
                  R.version$version.string, 
                  " is not consistent with R version 4.0.2 (2020-06-22).",
                  "Update R version to be consistent")))
  
  switch(EXPR = packageVersion("VAST") == "3.6.1",
         "TRUE" = print("R package VAST is consistent with version 3.6.1"),
         "FALSE" = print(
           paste0("R version ", packageVersion("VAST"), 
                  " is not consistent with version 3.6.1.",
                  "Update R package VAST to be consistent")))
  
  switch(EXPR = packageVersion( "FishStatsUtils") == "2.8.0",
         "TRUE" = "R package FishStatsUtils is consistent with version 2.8.0",
         "FALSE" = print(
           paste0("R version ", packageVersion( "FishStatsUtils"), 
                  " is not consistent with version 2.8.0.",
                  "Update R package FishStatsUtils to be consistent")))
}

vast_cpp_version <- "VAST_v12_0_0"

##################################################
####   Import CPUE dataset, species set spreadsheet
##################################################
master_data <- read.csv(file = "data/processed/goa_vast_data_input.csv" )

#################################################
## Loop over species to fit models with and without depth covariates
#################################################
spp_names <- sort(unique(master_data$COMMON_NAME))

for (ispp in spp_names[1:6]) {
  for (depth_in_model in c(F, T)) {
    
    ##################################################
    ## Create directory to store model results
    ##################################################
    result_dir <- paste0(VAST_dir, 
                         ispp, 
                         ifelse(test = depth_in_model,  
                                yes = "_depth", 
                                no = ""), 
                         "/")
    
    if (!dir.exists(result_dir)) dir.create(result_dir)  
    
    ##################################################
    ####   Subset species
    ##################################################
    data <- subset(master_data, 
                   COMMON_NAME == ispp)
    
    ##################################################
    ####   Prepare the dataframe for catch-rate data in the VAST format
    ##################################################
    data_geostat <- data.frame( spp = data$COMMON_NAME,
                                Year = data$YEAR,
                                Catch_KG = data$WEIGHT,
                                AreaSwept_km2 = data$EFFORT,
                                Lat = data$LATITUDE,
                                Lon = data$LONGITUDE, 
                                stringsAsFactors = T)
    
    data_geostat[, c("LOG_DEPTH", "LOG_DEPTH2") ] <-
      data[, c("LOG_DEPTH_EFH_CEN", "LOG_DEPTH_EFH_CEN_SQ")]
    
    ##################################################
    ####   Assign 10 fold partitions of the data
    ##################################################
    n_fold <- 10
    years <- paste0(unique(data_geostat$Year))
    NTime <- length(unique(data_geostat$Year))
    
    #Create unique stationID from the latlon. To make sure the ids are unique,
    #we use the table function to make sure there are 7900 records (as of 2019)=
    data_geostat$latlon <- paste0(data_geostat$Lat, data_geostat$Lon)
    table(table(data_geostat$latlon))
    
    #split data_geostat by year, then on each year-split, randomly assign 
    #fold numbers to the each unique station
    set.seed(2342)
    foldno <- lapply(
      #Split data_geostat by Year
      X = split.data.frame(data_geostat, 
                           f = data_geostat$Year),
      
      #For each year split, randomly assign fold numbers so that each year is 
      #equally split into n_folds folds
      FUN = function(x) {
        unique_loc <- unique(x$latlon)
        fold_no <- sample(x = 1:n_fold, 
                          size = length(unique_loc), 
                          replace = T)
        return(split(unique_loc, fold_no))
      })
    
    #Attach fold number to the data_geostat
    for (iyear in years) {
      for (ifold in paste(1:n_fold)) {
        data_geostat[data_geostat$latlon %in% foldno[[iyear]][[ifold]] , 
                     "fold"] = as.integer(ifold) 
      }
    }
    
    ##################################################
    ####   Spatial settings: The following settings define the spatial resolution 
    ####   for the model, and whether to use a grid or mesh approximation
    ####   Stratification for results
    ##################################################
    settings <- FishStatsUtils::make_settings( 
      Version = vast_cpp_version,
      n_x = 500,   # Number of knots
      Region = "User", #User inputted extrapolation grid
      purpose = "index2",
      bias.correct = FALSE,
      FieldConfig = c(
        "Omega1" = 1,   #Spatial random effect on occurence 
        "Epsilon1" = 1, #Spatiotemporal random effect on occurence 
        "Omega2" = 1,   #Spatial random effect on positive response 
        "Epsilon2" = 1  #Spatiotemporal random effect on positive response
      ), 
      "Options" = c("Calculate_Range" = F, 
                    "Calculate_effective_area" = F),
      ObsModel = c(2, 1),
      max_cells = Inf,
      use_anisotropy = T)
    
    ##################################################
    ####   Import "true" and not interpolated covariate 
    ####   data if using depth covariates
    ##################################################
    load("data/processed/grid_goa.RData")
    
    n_g <- nrow(grid_goa) #number of grid cells
    n_t <- diff(range(data_geostat$Year)) + 1 #Number of total years
    n_p <- 2 #two density covariates
    
    X_gtp <- array(dim = c(n_g, n_t, n_p) )
    for (i in 1:n_t) {
      X_gtp[, i, ] <- as.matrix(grid_goa[, c("LOG_DEPTH_EFH_CEN", 
                                             "LOG_DEPTH_EFH_CEN_SQ")])
    }
    
    ##################################################
    ####   Fit the model and save output
    ##################################################
    fit <- switch(paste0(depth_in_model),
                  "FALSE" = FishStatsUtils::fit_model( 
                    "settings" = settings,
                    "working_dir" = temp_res,
                    "Lat_i" = data_geostat[, "Lat"],
                    "Lon_i" = data_geostat[, "Lon"],
                    "t_i" = data_geostat[, "Year"],
                    "c_i" = as.numeric(data_geostat[, "spp"]) - 1,
                    "b_i" = data_geostat[, "Catch_KG"],
                    "a_i" = data_geostat[, "AreaSwept_km2"],
                    "getJointPrecision" = TRUE,
                    "newtonsteps" = 1,
                    "test_fit" = F,
                    "input_grid" = grid_goa),
                  
                  "TRUE" = FishStatsUtils::fit_model( 
                    "settings" = settings,
                    "working_dir" = temp_res,
                    "Lat_i" = data_geostat[, "Lat"],
                    "Lon_i" = data_geostat[, "Lon"],
                    "t_i" = data_geostat[, "Year"],
                    "c_i" = as.numeric(data_geostat[, "spp"]) - 1,
                    "b_i" = data_geostat[, "Catch_KG"],
                    "a_i" = data_geostat[, "AreaSwept_km2"],
                    "getJointPrecision" = TRUE,
                    "newtonsteps" = 1,
                    "test_fit" = F,
                    "input_grid" = grid_goa[, c("Area_km2", "Lon", "Lat")],
                    
                    ##Additional arguments for covariates
                    "X1_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
                    "X2_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
                    "covariate_data" = cbind(data_geostat[, c("Lat",
                                                              "Lon",
                                                              "LOG_DEPTH",
                                                              "LOG_DEPTH2",
                                                              "Catch_KG")],
                                             Year = NA),
                    "X_gtp" = X_gtp))
    
    ##################################################
    ####   Diagnostics plots
    ##################################################
    if(!dir.exists(paste0(result_dir, "/diagnostics"))) {
      dir.create(paste0(result_dir, "/diagnostics"))
    }
    
    plot(x = fit,
         working_dir = paste0(result_dir, "diagnostics/"))
    
    ##################################################
    ####   Save original model fit and copy output from the temporary res dir
    ##################################################
    save(list = "fit", file = paste0(result_dir, "/fit.RData"))
    file.copy(from = dir(temp_res, full.names = T), to = result_dir)
    
    ##################################################
    ####   10-fold Cross Validation
    ##################################################
    n_fold <- 10
    for (fI in 1:n_fold) { 
      if (!dir.exists(paste0(result_dir, "CV_", fI))) {
        dir.create(paste0(result_dir, "CV_", fI))
        
        file.copy(from = paste0(result_dir, get_latest_version(), 
                                c(".cpp", ".dll", ".o")),
                  to = paste0(result_dir, "CV_", fI, "/", 
                              get_latest_version(), 
                              c(".cpp", ".dll", ".o")))
        
      }
    } 
    
    # Loop through partitions, refitting each time with a different PredTF_i
    for (fI in 1:n_fold) {
      PredTF_i <- ifelse( test = data_geostat$fold == fI, 
                          yes = TRUE, 
                          no = FALSE )
      
      fit_CV <- switch(paste0(depth_in_model),
                       "FALSE" = FishStatsUtils::fit_model( 
                         "settings" = settings,
                         "working_dir" = temp_res,
                         "Lat_i" = data_geostat[, "Lat"],
                         "Lon_i" = data_geostat[, "Lon"],
                         "t_i" = data_geostat[, "Year"],
                         "c_i" = as.numeric(data_geostat[, "spp"]) - 1,
                         "b_i" = data_geostat[, "Catch_KG"],
                         "a_i" = data_geostat[, "AreaSwept_km2"],
                         "getJointPrecision" = TRUE,
                         "newtonsteps" = 1,
                         "test_fit" = F,
                         "input_grid" = grid_goa,
                         "PredTF_i" = PredTF_i,
                         "Parameters" = fit$ParHat),
                       
                       "TRUE" = FishStatsUtils::fit_model( 
                         "settings" = settings,
                         "working_dir" = temp_res,
                         "Lat_i" = data_geostat[, "Lat"],
                         "Lon_i" = data_geostat[, "Lon"],
                         "t_i" = data_geostat[, "Year"],
                         "c_i" = as.numeric(data_geostat[, "spp"]) - 1,
                         "b_i" = data_geostat[, "Catch_KG"],
                         "a_i" = data_geostat[, "AreaSwept_km2"],
                         "getJointPrecision" = TRUE,
                         "newtonsteps" = 1,
                         "test_fit" = F,
                         "input_grid" = grid_goa[, c("Area_km2", "Lon", "Lat")],
                         "PredTF_i" = PredTF_i,
                         
                         ##Additional arguments for covariates
                         "X1_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
                         "X2_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
                         "covariate_data" = cbind(data_geostat[, c("Lat",
                                                                   "Lon",
                                                                   "LOG_DEPTH",
                                                                   "LOG_DEPTH2",
                                                                   "Catch_KG")],
                                                  Year = NA),
                         "X_gtp" = X_gtp,
                         "Parameters" = fit$ParHat))
      
      ## Save fit
      save(list = "fit_CV",  
           file = paste0(result_dir, "CV_", fI, "/fit.RData"))
      
      ## Save predicted and observed CPUEs
      obs_cpue <- with(data_geostat[PredTF_i, ], Catch_KG / AreaSwept_km2)
      pred_cpue <- fit_CV$Report$D_i[PredTF_i]
      
      cv_performance <- list(cpues = data.frame(cv_fold = fI, 
                                                obs_cpue, 
                                                pred_cpue),
                             prednll = fit_CV$Report$pred_jnll)
      
      save(cv_performance,
           file = paste0(result_dir, "CV_", fI, 
                         "/crossval_fit_performance.RData"))
      
      ## Copy other outputs to the CV-fold directory
      file.copy(from = dir(temp_res, full.names = T), 
                to = paste0(result_dir, "CV_", fI, "/"))
      
    }
    
    ##################################################
    #### Prediction Grid: df of the grid to simulate data onto
    ##################################################
    grid_df <- data.frame()
    for(itime in unique(data$YEAR)) {
      grid_df <- rbind(grid_df,
                       data.frame(spp = ispp,
                                  Year = rep(itime, nrow(grid_goa)),
                                  Catch_KG = mean(data$WEIGHT),
                                  AreaSwept_km2 = grid_goa$Area_km2,
                                  Lat = grid_goa$Lat,
                                  Lon = grid_goa$Lon,
                                  LOG_DEPTH = grid_goa$LOG_DEPTH_EFH_CEN,
                                  LOG_DEPTH2 = grid_goa$LOG_DEPTH_EFH_CEN_SQ,
                                  stringsAsFactors = T)
      )
    }
    
    ###################################################
    ## Add New Points: set catch to NAs?
    ###################################################
    data_geostat_with_grid <- rbind(data_geostat[, names(grid_df)],
                                    grid_df)
    
    ##################################################
    ####   Fit the model and save output
    ##################################################
    pred_TF <- rep(1, nrow(data_geostat_with_grid))
    pred_TF[1:nrow(data)] <- 0
    
    fit_sim <-
      switch(paste0(depth_in_model),
             "FALSE" = FishStatsUtils::fit_model( 
               "settings" = settings,
               "working_dir" = temp_res,
               "Lat_i" = data_geostat_with_grid[, "Lat"],
               "Lon_i" = data_geostat_with_grid[, "Lon"],
               "t_i" = data_geostat_with_grid[, "Year"],
               "c_i" = as.numeric(data_geostat_with_grid[, "spp"]) - 1,
               "b_i" = data_geostat_with_grid[, "Catch_KG"],
               "a_i" = data_geostat_with_grid[, "AreaSwept_km2"],
               "getJointPrecision" = TRUE,
               "newtonsteps" = 1,
               "test_fit" = F,
               "input_grid" = grid_goa, 
               "PredTF_i" = pred_TF,
               "Parameters" = fit$ParHat),
             
             "TRUE" = FishStatsUtils::fit_model( 
               "settings" = settings,
               "working_dir" = temp_res,
               "Lat_i" = data_geostat_with_grid[, "Lat"],
               "Lon_i" = data_geostat_with_grid[, "Lon"],
               "t_i" = data_geostat_with_grid[, "Year"],
               "c_i" = as.numeric(data_geostat_with_grid[, "spp"]) - 1,
               "b_i" = data_geostat_with_grid[, "Catch_KG"],
               "a_i" = data_geostat_with_grid[, "AreaSwept_km2"],
               "getJointPrecision" = TRUE,
               "newtonsteps" = 1,
               "test_fit" = F,
               "input_grid" = grid_goa[, c("Area_km2", "Lon", "Lat")],
               "Parameters" = fit$ParHat,
               
               ##Additional arguments for covariates
               "X1_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
               "X2_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
               "covariate_data" = cbind(data_geostat[, c("Lat",
                                                         "Lon",
                                                         "LOG_DEPTH",
                                                         "LOG_DEPTH2",
                                                         "Catch_KG")],
                                        Year = NA),
               "X_gtp" = X_gtp, 
               "PredTF_i" = pred_TF))
    
    ##################################################
    ####   Save model fit
    ##################################################
    save(list = "fit_sim", file = paste0(result_dir, "/fit_sim.RData"))
    
    ##################################################
    ####   Simulate 1000 iterations of data
    ##################################################
    sim_data <- array(data = NA, dim = c(nrow(grid_goa),
                                         length(unique(data$YEAR)),
                                         1000))
    
    for (isim in 1:1000) {
      Sim1 <- FishStatsUtils::simulate_data(fit = fit_sim, 
                                            type = 1, 
                                            random_seed = isim)
      sim_data[, , isim] <- matrix(data = Sim1$b_i[pred_TF == 1] * 0.001, 
                                   nrow = nrow(grid_goa), 
                                   ncol = length(unique(data$YEAR)))
      if(isim%%100 == 0) print(paste("Done with", ispp, "Iteration", isim))
    }
    
    save(sim_data, file = paste0(result_dir, "/simulated_data.RData"))
    
  }
}
