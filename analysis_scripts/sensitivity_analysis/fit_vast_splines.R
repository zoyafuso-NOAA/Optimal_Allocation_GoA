###############################################################################
## Project:      Univariate VAST model runs with basis splines
## Author:       Zack Oyafuso (zack.oyafuso@noaa.gov)
## Contributors: Lewis Barnett (lewis.barnett@noaa.gov)
##               Jim Thorson"s VAST wiki example
##           (https://github.com/James-Thorson-NOAA/VAST/wiki/Crossvalidation)
## Description:  Run single-species VAST models with depth as a density 
##               covariate modelled as a basis splinecovariate. Run 10-fold 
##               Cross Validation
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
VAST_dir <- c("F:/VAST_Runs_splines/")
if(!dir.exists(VAST_dir)) dir.create(VAST_dir, recursive = T)

temp_res <- paste0(getwd(), "/results/temp_res/")
if(!dir.exists(temp_res)) dir.create(temp_res, recursive = T)

##################################################
####   Load packages, make sure versions are consistent
##################################################
library(VAST)
library(splines)

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

for (ispp in spp_names[4]) {
  
  ##################################################
  ## Create directory to store model results
  ##################################################
  result_dir <- paste0(VAST_dir, ispp, "_depth/")
  
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
  fit <- FishStatsUtils::fit_model( 
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
    "X1_formula" = ~ splines::bs(LOG_DEPTH, degree = 3, intercept = FALSE),
    "X2_formula" = ~ splines::bs(LOG_DEPTH, degree = 3, intercept = FALSE),
    
    # "X1_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
    # "X2_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
    "covariate_data" = cbind(data_geostat[, c("Lat",
                                              "Lon",
                                              "LOG_DEPTH",
                                              # "LOG_DEPTH2",
                                              "Catch_KG")],
                             Year = NA),
    "X_gtp" = X_gtp)
  
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
    
    fit_CV <- FishStatsUtils::fit_model( 
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
      "X1_formula" = ~ splines::bs(LOG_DEPTH, degree = 3, intercept = FALSE),
      "X2_formula" = ~ splines::bs(LOG_DEPTH, degree = 3, intercept = FALSE),
      
      # "X1_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
      # "X2_formula" =  "Catch_KG ~ LOG_DEPTH + LOG_DEPTH2",
      "covariate_data" = cbind(data_geostat[, c("Lat",
                                                "Lon",
                                                "LOG_DEPTH",
                                                # "LOG_DEPTH2",
                                                "Catch_KG")],
                               Year = NA),
      "X_gtp" = X_gtp,
      "Parameters" = fit$ParHat)
    
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
  
  
}
