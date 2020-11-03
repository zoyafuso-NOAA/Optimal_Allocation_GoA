###############################################################################
## Project:       Simulate Surveys
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Simulate a Stratified Random Survey of the Gulf of 
##                Alaska Groundfish Survey Based on Current Stratifications
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(readxl)
library(spatialEco)
library(sp)
library(VAST)

##################################################
####   Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[3]
VAST_model <- "11" 

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/")
VAST_dir <- "G:/Oyafuso/VAST_Runs_EFH/Single_Species/"

##################################################
####   Load simulation functions
##################################################
source( paste0(dirname(github_dir), "/modified_functions/sim_fns.R") )

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load(paste0(github_dir, 'full_domain/optimization_data.RData'))

load(paste0(github_dir, '/RMSE_VAST_models.RData'))
load(paste0(dirname(github_dir), '/data/Extrapolation_depths.RData'))
# Extrapolation_depths <- subset(Extrapolation_depths, stratum != 0)

GOA_allocations <- readxl::read_xlsx(
  path = paste0(dirname(github_dir), 
                '/data/GOA 2019 stations by stratum.xlsx')
)
GOA_allocations3 <- readxl::read_xlsx(
  path = paste0(dirname(github_dir), 
                '/data/GOA2019_ 3 boat_825_RNDM_stations.xlsx')
) 

##################################################
####   Create dataframe of effort allocations across boats
##################################################
allocations <- data.frame(Stratum = sort(unique(GOA_allocations3$stratum)),
                          boat3 = aggregate(id ~ stratum, 
                                            data = GOA_allocations3, 
                                            FUN = length)$id,
                          boat2 = c(GOA_allocations$`Number stations`,
                                    rep(0,5)))
allocations$boat1 <- ceiling(allocations$boat2 / 2)

allocations$boat1 <- ifelse(allocations$boat1 == 0, 0, 
                            ifelse(allocations$boat1 == 1, 2, 
                                   allocations$boat1))

allocations <- rbind(data.frame(Stratum = 0, boat3 = 0, boat2 = 0, boat1 = 0),
                     allocations)

##################################################
####   Result Objects
##################################################
Current_sim_mean <- Current_sim_cv <- Current_rel_bias_est <- 
  Current_sim_mean_trawl <- Current_sim_cv_trawl <- Current_rel_bias_est_trawl <- 
  STRS_sim_mean <- STRS_sim_cv <- STRS_rel_bias_est <-
  STRS_sim_mean_trawl <- STRS_sim_cv_trawl <- STRS_rel_bias_est_trawl <-
  array(dim = c(3, NTime, ns, nboats, Niters), 
        dimnames = list(NULL, 
                        paste0("year_", 1:NTime), 
                        sci_names, 
                        paste0("boat_", 1:3), 
                        NULL ))

Current_true_cv_array <- Current_rrmse_cv_array <- 
  Current_true_cv_array_trawl <- Current_rrmse_cv_array_trawl <- 
  STRS_true_cv_array <- STRS_rrmse_cv_array <-  
  STRS_true_cv_array_trawl <- STRS_rrmse_cv_array_trawl <-  
  array(dim = c(3, NTime, ns, nboats), 
        dimnames = list(NULL, 
                        paste0("year_", 1:NTime), 
                        sci_names, 
                        paste0("boat_", 1:3)))

res_obj_names <- apply(X = expand.grid(c("Current_", "STRS_"),
                                       c("sim_mean", "sim_cv", "rel_bias_est", 
                                         "true_cv_array", "rrmse_cv_array"),
                                       c("", "_trawl")),
                       MARGIN = 1,
                       FUN = paste0,
                       collapse = "")

##################################################
####   
##################################################
# load("~/GitHub/Optimal_Allocation_GoA/model_11/simulation_result.RData")
#Load optimized solutions calculated on the full domain
# load("~/GitHub/Optimal_Allocation_GoA/model_11/fit_density.RData")
load(paste0(github_dir, "full_domain/Spatiotemporal_Optimization",
            "/optimization_knitted_results.RData"))

for (iter in 1:1000) {
  
  set.seed(1000 + iter)
  
  #Temporary matrix to hold simulated values
  #3-D array that holds the 
  #  1) Predicted density
  #  2) simulated data with measurement error
  #  3) simulated data with simulated random and fixed effects
  temp_sim <- array(dim = c(3,
                            N,
                            ns,
                            NTime),
                    dimnames = list(NULL, NULL, paste0("Y", 1:ns), NULL))
  
  #Populate temp_sim with single-species VAST models
  for (ispp in 1:ns) {
    chosen_model = paste0(VAST_dir,
                          sci_names[ispp],
                          ifelse(RMSE$depth_in_model[ispp],
                                 "_depth",
                                 ""), "/")
  
  #Load fitted model
  load(paste0(chosen_model, "fit.RData"))
  
  # #DYnamic load the dll file
  dyn.load(paste0(chosen_model, "VAST_v12_0_0"))
  
  #Predicted density
  temp_sim[1,,ispp,] <- fit$Report$D_gct[, , Years2Include]
  # temp_sim[1,,,] <- D_gct[, , Years2Include]
  
  #simulated data with measurement error
  temp_sim[2,,ispp,] <-
    FishStatsUtils::simulate_data(fit = fit,
                                  type = 1)$D_gct[ , ,Years2Include]
  
  #simulated data with simulated random and fixed effects
  temp_sim[3,,ispp,] <-
    FishStatsUtils::simulate_data(fit = fit,
                                  type = 3)$D_gct[ , ,Years2Include]
  
  #Unload the library (I don't know whether this step is necessary...)
  dyn.unload(paste0(chosen_model, "VAST_v12_0_0"))
  }
  
  #Loop through each simulation type
  for (idomain in c("full_domain", "trawlable")) {
    
    # load(paste0(github_dir, idomain,"/Spatiotemporal_Optimization",
    #             "/optimization_knitted_results.RData"))
    
    for (isurvey in c("Current", "STRS")) { #Current or Optimized Survey
      for (isim in 1:3) {
        for (iboat in 1:3) {
          
          if(isurvey == "STRS") {
            #Load optimization data, only focusing on 15 strata for now
            sub_settings = subset(settings, strata == 15)
            idx <- which.min(abs(sub_settings$n - samples[iboat]))
          }
          
          sim_survey <- 
            do_STRS(
              input <- list(
                "density" = switch(
                  idomain,
                  "full_domain" = temp_sim[isim, , ,],
                  "trawlable" = temp_sim[isim, , ,]),
                
                "solution" = switch(
                  paste(isurvey, idomain),
                  "Current full_domain" = with(Extrapolation_depths, 
                                               stratum),
                  "Current trawlable" = with(Extrapolation_depths, 
                                             stratum),
                  "STRS full_domain" = res_df[, 1 + idx],
                  "STRS trawlable" = res_df[, 1 + idx]),
                
                "allocation" = switch( 
                  paste(isurvey, idomain),
                  "Current full_domain" = allocations[, paste0("boat", iboat)],
                  "Current trawlable" = allocations[, paste0("boat", iboat)],
                  "STRS full_domain" = strata_list[[idx]]$Allocation,
                  "STRS trawlable" = strata_list[[idx]]$Allocation),
                
                "true_density" = true_mean,
                
                "domain_idx" = switch(
                  idomain,
                  "full_domain" = rep(T, N),
                  "trawlable" = Extrapolation_depths$shallow_trawlable
                ))
              )
          
          survey_name <- switch(
            isurvey, 
            "Current" = "Current_", 
            "STRS" = "STRS_")
          domain_name <- switch(
            idomain, 
            "full_domain" = "", 
            "trawlable" = "_trawl")
          
          stmt <- paste0(survey_name, "sim_mean", domain_name, 
                         "[isim, , , iboat, iter] = sim_survey$mean_denisty")
          eval(parse(text = stmt))
          
          stmt <- paste0(survey_name, "sim_cv", domain_name, 
                         "[isim, , , iboat, iter] = sim_survey$cv")
          eval(parse(text = stmt))
          
          stmt <- paste0(survey_name, "rel_bias_est", domain_name, 
                         "[isim, , , iboat, iter] = sim_survey$rel_bias")
          eval(parse(text = stmt))
          
          if (iter%%5 == 0) {
            stmt <- paste0(survey_name, "true_cv_array", domain_name, 
                           "[isim, , , iboat] <- ", "as.matrix(apply(X = ", 
                           survey_name, "sim_mean", domain_name, 
                           "[isim, , , iboat, ], MARGIN = 1:2, FUN = sd, ",
                           "na.rm = T) / true_mean)")
            eval(parse(text = stmt))
            
            stmt <- paste0(survey_name, "rrmse_cv_array", domain_name,
                           "[isim, , , iboat] <- sqrt(apply(X = sweep(x = ",
                           survey_name, "sim_cv", domain_name, 
                           "[isim, , , iboat,], STATS = ", survey_name, 
                           "true_cv_array", domain_name, "[isim, , , iboat],",
                           " MARGIN = 1:2, FUN = '-')^2,  
                           MARGIN = 1:2, FUN = mean, na.rm = T)) / apply(", 
                           survey_name, "sim_cv", domain_name, 
                           "[isim, , , iboat, ], MARGIN = 1:2, ",
                           "FUN = mean, na.rm = T)")
            eval(parse(text = stmt))
          } 
        }
      }
    }
  }
  

  
  if(iter%%10 == 0) {
    
    print(paste("Finished with Iteration", iter))
    save(list = res_obj_names,
                        file = paste0(github_dir, "simulation_result.RData"))
    }
}


