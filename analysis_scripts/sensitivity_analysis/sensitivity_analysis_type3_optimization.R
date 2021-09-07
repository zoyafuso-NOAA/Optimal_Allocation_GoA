###############################################################################
## Project:       Sensitivity Analysis for Survey Optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
VAST_sim_data_dir <- "E:/VAST_Runs/"
github_dir <- getwd()

library(VAST)
library(SamplingStrata)
library(raster)
library(sp)
library(RColorBrewer)

##################################################
####   Load data
##################################################
load("data/processed/optimization_data.RData")
load("data/processed/prednll_VAST_models.RData")
load("data/processed/grid_goa.RData")

##################################################
####   Sensitivty Scenarios
##################################################
scen <- expand.grid(strata_vars = 1:2,
                    deep_stations = c(TRUE, FALSE),
                    isim = 1:10)

for (irow in 1:12) {
  frame <- frame_district[, c("domainvalue", "id", 
                              "WEIGHT",
                              paste0("Y", spp_idx_opt), 
                              paste0("Y", spp_idx_opt,
                                     "_SQ_SUM"))]
  
  ## Only subset the species included in the optimization, indexed by 
  ## spp_idx_opt, then reorganize field names
  names(frame)[names(frame) %in% paste0("Y", spp_idx_opt)] <- 
    paste0("Y", 1:ns_opt)
  names(frame)[names(frame) %in% paste0("Y", spp_idx_opt, "_SQ_SUM")] <- 
    paste0("Y", 1:ns_opt, "_SQ_SUM")
  
  ##################################################
  ####   Import simulated dataset
  ##################################################
  load(paste0(github_dir, "/results/sensitivity_analysis_type3/dens_vals_sim", 
         scen$isim[irow],".RData"))
  
  ##################################################
  ####   Input simulated density type
  ##################################################
  frame[, c(paste0("Y", 1:ns_opt), paste0("Y", 1:ns_opt, "_SQ_SUM") )] <-
    cbind(matrix(data = apply(X = dens_vals,
                              MARGIN = c(1, 2), 
                              FUN = sum),
                 ncol = ns_opt,
                 dimnames = list(NULL, paste0("Y", 1:ns_opt))),
          
          matrix(data = apply(X = dens_vals,
                              MARGIN = c(1, 2), 
                              FUN = function(x) sum(x^2)),
                 ncol = ns_opt,
                 dimnames = list(NULL, paste0("Y", 1:ns_opt, "_SQ_SUM"))) )
  
  ##################################################
  ####   Input strata variables
  ##################################################
  frame[, paste0("X", 1:scen$strata_vars[irow])] <- 
    frame_district[, switch(paste0(scen$strata_vars[irow]),
                            "1" = c("X2"), 
                            "2" = c("X1", "X2"))]
  
  ##################################################
  ####   Discretize deep stations if needed
  ##################################################
  if (scen$deep_stations[irow] == FALSE) {
    if (scen$strata_vars[irow] == 1) frame$X1 <- ifelse(test = frame$X1 > 300, 
                                                        yes = 1000, 
                                                        no = frame$X1) 
    if (scen$strata_vars[irow]  == 2) frame$X2 <- ifelse(test = frame$X2 > 300, 
                                                         yes = 1000, 
                                                         no = frame$X2) 
  }
  
  ## Initiate CVs to be those calculated under SRS, assign to a variable 
  ## named cv_constraints
  ## buildStrataDF calculates the stratum means and variances, X1 = 1 
  ##     means to calculate those statics on the whole domain
  
  temp_vars <- names(frame) %in% paste0("X", 1:scen$strata_vars[irow])
  srs_stats <- SamplingStrata::buildStrataDF( 
    dataset = cbind( subset(frame, select = names(frame)[!temp_vars]),
                     X1 = 1))
  
  srs_n <- as.numeric(samples[2] * table(frame$domainvalue) / n_cells)
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
    cv[[paste0("CV", spp)]] <- cv_constraints[, spp]
  cv[["DOM"]] <- 1:n_districts
  cv[["domainvalue"]] <- 1:n_districts
  cv <- as.data.frame(cv)
  
  #Set wd for output files, create a directory if it doesn"t exist yet
  temp_dir = paste0(github_dir, "/results/sensitivity_analysis_type3/", 
                    "sim", scen$isim[irow], "_stratavars", scen$strata_vars[irow],
                    "_deep", scen$deep_stations[irow], "/")
  
  if(!dir.exists(temp_dir)) dir.create(temp_dir, recursive = T)
  
  setwd(temp_dir)
  
  solution <- optimStrata(method = "continuous",
                          errors = cv, 
                          framesamp = frame,
                          iter = 300,
                          pops = 50,
                          # iter = 10,
                          # pops = 10,
                          elitism_rate = 0.1,
                          mut_chance = 1 / 6,
                          nStrata = rep(5, n_districts),
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
  
  plot_solution <- as.factor(paste0(
    "DOM", solution$framenew$DOMAINVALUE,
    " STR", solution$framenew$STRATO))
  
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
    coords = grid_goa[, c("E_km", "N_km")],
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
       ybottom = tapply(X = grid_goa$N_km, 
                        INDEX = district_vals,
                        FUN = min), 
       ytop = tapply(X = grid_goa$N_km, 
                     INDEX = district_vals,
                     FUN = max))
  
  text(x = rowMeans(districts[, c("W_UTM", "E_UTM")]),
       y = tapply(X = grid_goa$N_km, 
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
       ybottom = tapply(X = grid_goa$N_km, 
                        INDEX = district_vals,
                        FUN = min), 
       ytop = tapply(X = grid_goa$N_km, 
                     INDEX = district_vals,
                     FUN = max))
  
  text(x = rowMeans(districts[, c("W_UTM", "E_UTM")]),
       y = tapply(X = grid_goa$N_km, 
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
  
  points(grid_goa[sample_vec, c("E_km", "N_km")],
         pch = 16, cex = 0.5)
  
  dev.off()
}

