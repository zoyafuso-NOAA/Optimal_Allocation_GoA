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
# VAST_dir <- c("E:/VAST_Runs/")
# if(!dir.exists(VAST_dir)) dir.create(VAST_dir, recursive = T)
# 
# temp_res <- paste0(getwd(), "/results/temp_res/")
# if(!dir.exists(temp_res)) dir.create(temp_res, recursive = T)

##################################################
####   Load packages, make sure versions are consistent
##################################################
library(VAST)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Package version checks ----
##   Updated every year
##   2023 TOR is in this google doc: 
##   https://docs.google.com/document/d/18CeXcHhHK48hrtkiC6zygXlHj6YVrWEd/edit
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
current_year <- 2023
vast_cpp_version <- "VAST_v14_0_1"
pck_version <- c("VAST" = "3.10.0",
                 "FishStatsUtils" = "2.12.0",
                 "Matrix" = "1.5-3",
                 "TMB" = "1.9.2",
                 "DHARMa" = "0.4.6")

for (pck in 1:length(pck_version)) {
  temp_version <- packageVersion(pkg = names(pck_version)[pck])
  
  if(temp_version == pck_version[pck])
    message(paste0("The version of the '", names(pck_version)[pck], 
                   "' package (", temp_version, ") is consistent",
                   " with the ", current_year, " TOR."))
  
  if(!temp_version == pck_version[pck])
    message(paste0("WARNING: ", 
                   "The version of the '", names(pck_version)[pck], 
                   "' package (", temp_version, ") is NOT consistent",
                   " with the ", current_year, " TOR. Please update the '", 
                   names(pck_version)[pck], "' package to ", 
                   pck_version[pck]))
  
  rm(pck, temp_version)
}
##################################################
####   Import CPUE dataset, species set spreadsheet
##################################################
goa_data_geostat <- read.csv(file = "data/processed/goa_data_geostat.csv" )
goa_interpolation_grid <- 
  read.csv(file = "data/processed/goa_interpolation_grid.csv")

#################################################
## Loop over species to fit models with and without depth covariates
#################################################
spp_names <- sort(x = unique(x = goa_data_geostat$Species))

for (depth_in_model in c(F, T)[2]) {
  for (ispp in spp_names[7]) {
    
    ##################################################
    ## Create directory to store model results
    ##################################################
    result_dir <- paste0(getwd(), "/temp/", ispp, ifelse(test = depth_in_model,  
                                                         yes = "_depth/", 
                                                         no = "/"))
    
    if (!dir.exists(paths = result_dir)) 
      dir.create(path = result_dir, recursive = TRUE)  
    
    ##################################################
    ####   Subset species
    ##################################################
    temp_df <- subset(x = goa_data_geostat, subset = Species == ispp)
    
    ##################################################
    ####   Assign 10 fold partitions of the data
    ##################################################
    n_fold <- 10
    years <- paste0(sort(x = unique(x = temp_df$Year)))
    NTime <- length(x = unique(x = temp_df$Year))
    
    #split data_geostat by year, then on each year-split, randomly assign 
    #fold numbers to the each unique station (hauljoin)
    set.seed(2342)
    foldno <- lapply(
      #Split data_geostat by Year
      X = split.data.frame(x = temp_df, f = temp_df$Year),
      #For each year split, randomly assign fold numbers so that each year is 
      #equally split into n_folds folds
      FUN = function(x) {
        unique_loc <- unique(x = x$Hauljoin)
        fold_no <- sample(x = 1:n_fold, 
                          size = length(unique_loc), 
                          replace = T)
        return(split(unique_loc, fold_no))
      })
    
    #Attach fold number to the data_geostat
    for (iyear in years) {
      for (ifold in paste(1:n_fold)) {
        temp_df[temp_df$Hauljoin %in% foldno[[iyear]][[ifold]] , 
                "fold"] = as.integer(x = ifold) 
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
      bias.correct = TRUE,
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
    n_g <- nrow(x = goa_interpolation_grid) #number of grid cells
    n_t <- length(x = unique(x = temp_df$Year)) #Number of total years
    n_p <- 1 #two density covariates
    
    X_gtp <- array(dim = c(n_g, n_t, n_p) )
    for (i in 1:n_t) X_gtp[, i, ] <- 
      as.matrix(x = goa_interpolation_grid[, c("LOG10_DEPTH_M_CEN")])
    
    ##################################################
    ####   Fit the model and save output
    ##################################################
    ####  Set general arguments to FishStatsUtils::fit_model()
    vast_arguments <- list(
      
      ## Input settings
      "working_dir" = result_dir,
      "settings" = settings,
      
      ## Interpolation grid locations and total areas
      "input_grid" = goa_interpolation_grid[, c("Area_km2", "Lon", "Lat")],
      
      ## Data inputs
      "Lat_i" = temp_df[, "Lat"],
      "Lon_i" = temp_df[, "Lon"],
      "t_i" = temp_df[, "Year"],
      "c_i" = rep(0, nrow(x = temp_df)),
      "b_i" = units::as_units(x = temp_df$Catch_KG, "kg"),
      "a_i" = units::as_units(x = temp_df$AreaSwept_km2, "km2"),
      
      ## Output settings
      "getJointPrecision" = TRUE,
      
      ## Model tuning
      "newtonsteps" = 1,
      "test_fit" = F,
      
      ## Covariate data: bathy from grids come from the EFH bathymetry layer
      ## and bathy from data_geostat are those collected from the survey
      "covariate_data" = cbind(
        rbind(temp_df[, c("Lat", "Lon", "LOG10_DEPTH_M")],
              goa_interpolation_grid[, c("Lat", "Lon", "LOG10_DEPTH_M")]),
        Year = NA))
    
    if (depth_in_model)
      vast_arguments[["X1_formula"]] <- vast_arguments[["X2_formula"]] <-
      ~ splines::bs(LOG10_DEPTH_M,
                    degree = 2,
                    intercept = FALSE)
    
    ## Initial Fit
    fit <- do.call(what = FishStatsUtils::fit_model, args = vast_arguments)
    
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
    ## Save original model fit
    saveRDS(object =  fit, file = paste0(result_dir, "/fit.RDS"))
    
    ## Plot depth effects
    if (depth_in_model) {
      
      # Must add data-frames to global environment (hope to fix in future)
      covariate_data_full = fit$effects$covariate_data_full
      catchability_data_full = fit$effects$catchability_data_full
      depth_effect_pred1 <-
        VAST::Effect.fit_model( fit,
                                focal.predictors = c("LOG10_DEPTH_M"),
                                which_formula = "X1",
                                xlevels = 100,
                                transformation = list(link=identity,
                                                      inverse=identity) )
      depth_effect_pred2 <-
        VAST::Effect.fit_model( fit,
                                focal.predictors = c("LOG10_DEPTH_M"),
                                which_formula = "X2",
                                xlevels = 100,
                                transformation = list(link=identity,
                                                      inverse=identity) )
      
      ## Save a reduced form of the predicted effects
      depth_effect_pred1 <-
        with(depth_effect_pred1,
             data.frame(x = x, pred = fit, lower = lower, upper = upper))
      
      depth_effect_pred2 <-
        with(depth_effect_pred2,
             data.frame(x = x, pred = fit, lower = lower, upper = upper))
      
      pdf(file = paste0(result_dir, "depth_effects.pdf"),
          width = 5, height = 5, onefile = TRUE )
      for (ipred in 1:2) {
        temp <- get(paste0("depth_effect_pred", ipred))
        plot(x = temp$LOG10_DEPTH_M, y = temp$pred, type = "n", xaxt = "n",
             ylim = range(temp[, c("lower", "upper")]),
             las = 1, ylab = "Marginal Effect", xlab = "Depth (m)",
             main = paste0("Depth Effect on ", c("1st", "2nd")[ipred],
                           " Predictor\n", ispp))
        with(temp, polygon(x = c(LOG10_DEPTH_M, rev(LOG10_DEPTH_M)),
                           y = c(lower, rev(upper)),
                           col = "cornflowerblue", border = F))
        with(temp, lines(x = LOG10_DEPTH_M, y = pred, lwd = 2))
        axis(side = 1, labels = NA, tck = -0.015,
             at = log10(c(seq(10, 90, 10), seq(100, 1000, 100))))
        axis(side = 1, labels = NA, tck = -0.03,
             at = log10(c(10,20,50,100,200,500,1000)))
        axis(side = 1, at = log10(c(10,20,50,100,200,500,1000)),
             labels = c(10,20,50,100,200,500,1000))
      }
      dev.off()
      
      save(list = c("depth_effect_pred1", "depth_effect_pred2"),
           file = paste0(result_dir, "depth_effect.RData"))
    }
    
    ##################################################
    ####   10-fold Cross Validation
    ####   First, save MLE parameter values to the VAST arguments to be used
    ####   as starting values for the subsequent CV runs
    ##################################################
    vast_arguments$Parameters <- fit$ParHat
    
    ## Loop through partitions, refitting each time with a different PredTF_i
    for (fI in 1:n_fold) { ## Loop over folds -- start
      
      ## Create directory for CV run
      if (!dir.exists(paths = paste0(result_dir, "CV_", fI))){
        dir.create(path = paste0(result_dir, "CV_", fI))
        file.copy(from = paste0(result_dir, "Kmeans_knots-",
                                settings$n_x, ".RData"),
                  to = paste0(result_dir, "CV_", fI, "/"))
      }
      
      ## Update which indices are withheld for prediction
      PredTF_i <- ifelse( test = temp_df$fold == fI,
                          yes = TRUE,
                          no = FALSE )
      vast_arguments$PredTF_i <- PredTF_i
      vast_arguments$working_dir <- paste0(result_dir, "CV_", fI, "/")
      
      ## Turn bias correction off for Cross Validation runs  
      vast_arguments$settings$bias.correct <- FALSE
      
      ## Fit CV run
      fit_CV <- do.call(what = FishStatsUtils::fit_model, args = vast_arguments)
      
      ## Save predicted and observed CPUEs
      obs_cpue <- with(temp_df[PredTF_i, ], Catch_KG / AreaSwept_km2)
      pred_cpue <- fit_CV$Report$D_i[PredTF_i]
      
      cv_performance <-
        list(cpues = data.frame(cv_fold = fI,
                                obs_cpue,
                                pred_cpue),
             prednll = fit_CV$Report$pred_jnll,
             max_gradient = fit_CV$parameter_estimates$max_gradient)
      
      ## Save fit
      saveRDS(object =  fit_CV,
              file = paste0(result_dir, "CV_", fI, "/fit.RDS"))
      saveRDS(object = cv_performance,
              file = paste0(result_dir, "CV_", fI,
                            "/crossval_fit_performance.RDS"))
      
      if (depth_in_model) {
        # Must add data-frames to global environment (hope to fix in future)
        covariate_data_full = fit$effects$covariate_data_full
        catchability_data_full = fit$effects$catchability_data_full
        depth_effect_pred1 <-
          VAST::Effect.fit_model(fit,
                                 focal.predictors = c("LOG10_DEPTH_M"),
                                 which_formula = "X1",
                                 xlevels = 100,
                                 transformation = list(link=identity,
                                                       inverse=identity) )
        depth_effect_pred2 <-
          VAST::Effect.fit_model(fit,
                                 focal.predictors = c("LOG10_DEPTH_M"),
                                 which_formula = "X2",
                                 xlevels = 100,
                                 transformation = list(link=identity,
                                                       inverse=identity) )
        
        ## Save a reduced form of the predicted effects
        depth_effect_pred1 <-
          with(depth_effect_pred1,
               data.frame(x = x, pred = fit, lower = lower, upper = upper))
        
        depth_effect_pred2 <-
          with(depth_effect_pred2,
               data.frame(x = x, pred = fit, lower = lower, upper = upper))
        
        save(list = c("depth_effect_pred1", "depth_effect_pred2"),
             file = paste0(result_dir, "CV_", fI, "/depth_effect.RData"))
      }
      
    } ## Loop over folds -- end
    
    ##################################################
    #### Prediction Grid: df of the grid to simulate data onto
    ##################################################
    grid_df <- data.frame()
    for (itime in sort(x = unique(x = temp_df$Year))) {
      grid_df <- 
        rbind(grid_df,
              data.frame(
                Species = ispp,
                Year = rep(itime, nrow(x = goa_interpolation_grid)),
                Catch_KG = mean(temp_df$Catch_KG),
                AreaSwept_km2 = goa_interpolation_grid$Area_km2,
                Lat = goa_interpolation_grid$Lat,
                Lon = goa_interpolation_grid$Lon,
                LOG10_DEPTH_M_CEN = goa_interpolation_grid$LOG10_DEPTH_M_CEN,
                stringsAsFactors = T)
        )
    }
    
    ###################################################
    ## Add New Points: set catch to NAs?
    ###################################################
    data_geostat_with_grid <- rbind(temp_df[, names(x = grid_df)],
                                    grid_df)
    
    ##################################################
    ####   Fit the model and save output
    ##################################################
    pred_TF <- rep(1, nrow(x = data_geostat_with_grid))
    pred_TF[1:nrow(x = temp_df)] <- 0
    
    vast_arguments$Lat_i <- data_geostat_with_grid[, "Lat"]
    vast_arguments$Lon_i <- data_geostat_with_grid[, "Lon"]
    vast_arguments$t_i <- data_geostat_with_grid[, "Year"]
    vast_arguments$a_i <- data_geostat_with_grid[, "AreaSwept_km2"]
    vast_arguments$b_i <- data_geostat_with_grid[, "Catch_KG"]
    vast_arguments$c_i <- rep(x = 0, nrow(x = data_geostat_with_grid))
    vast_arguments$PredTF_i <- pred_TF
      
    ## Fit CV run and save model fit
    vast_arguments$settings$bias.correct <- TRUE
    fit_sim <- do.call(what = FishStatsUtils::fit_model, args = vast_arguments)
    saveRDS(object = "fit_sim", 
            file = paste0(result_dir, "/fit_sim.RDS"))
    
    ##################################################
    ####   Simulate 1000 iterations of data
    ##################################################
    sim_data <- array(data = NA, dim = c(n_g, n_t, 1000))
    
    for (isim in 1:1000) { ## Loop over replicates -- start
      Sim1 <- FishStatsUtils::simulate_data(fit = fit_sim, 
                                            type = 1, 
                                            random_seed = isim)
      sim_data[, , isim] <- matrix(data = Sim1$b_i[pred_TF == 1] * 0.001, 
                                   nrow = n_g, 
                                   ncol = n_t)
      if(isim%%100 == 0) print(paste("Done with", ispp, "Iteration", isim))
    } ## Loop over replicates -- end
    
    ## Save simulated data
    saveRDS(object = sim_data, 
            file = paste0(result_dir, "/simulated_data.RDS"))
    
  }
}
