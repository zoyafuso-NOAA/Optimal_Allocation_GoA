###############################################################################
## Project:       Calculate Population Variances
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Calculate population variances for different sampling schemes
##                Simple Random Sampling
##                Current Stratified Random Sampling Design
##                These population variances are more comparable to the 
##                optimization, as the optimization is 
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[3]

github_dir <- paste0(c('/Users/zackoyafuso/Documents/', 
                       'C:/Users/Zack Oyafuso/Documents/',
                       'C:/Users/zack.oyafuso/Work/')[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")

##################################################
####  Install a forked version of the SamplingStrata Package from 
####  zoyafuso-NOAA's Github page
####
####  Import other required packages
##################################################
library(devtools)
devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata")
library(SamplingStrata)
library(readxl)
library(tidyr)

##################################
## Import Operating Model
##################################
load(paste0(github_dir, "data/Extrapolation_depths.RData") )
load(paste0(github_dir, 'data/fit_density.RData'))
load(paste0(github_dir, 'data/optimization_data.RData'))
load(paste0(github_dir, "results/Single_Species_Optimization/",
            "optimization_knitted_results.RData"))

frame <- cbind(data.frame(domainvalue = 1,
                          id = 1:N,
                          X1 = with(Extrapolation_depths, E_km - min(E_km)),
                          X2 = Extrapolation_depths$DEPTH_EFH,
                          WEIGHT = NTime),
               
               matrix(data = apply(X = D_gct[, , Years2Include],
                                   MARGIN = c(1, 2), 
                                   FUN = sum),
                      ncol = ns_all,
                      dimnames = list(NULL, paste0("Y", 1:ns_all))),
               
               matrix(data = apply(X = D_gct[, , Years2Include],
                                   MARGIN = c(1, 2), 
                                   FUN = function(x) sum(x^2)),
                      ncol = ns_all,
                      dimnames = list(NULL, paste0("Y", 1:ns_all, "_SQ_SUM")))
)

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
GOA_allocations <- readxl::read_xlsx(
  path = paste0(github_dir, '/data/GOA 2019 stations by stratum.xlsx')
)
GOA_allocations3 <- readxl::read_xlsx(
  path = paste0(github_dir, '/data/GOA2019_ 3 boat_825_RNDM_stations.xlsx')
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
allocations$boat1 = ceiling(allocations$boat2 / 2)

allocations$boat1 = ifelse(allocations$boat1 == 0, 0, 
                           ifelse(allocations$boat1 == 1, 2, 
                                  allocations$boat1))

allocations <- rbind(c("Stratum" = 0, "boat3" = 0, "boat2" = 0, "boat1" = 0),
                     allocations)

##################################
## Calculate Population CVs under Simple Random Sampling
##################################
frame_SRS <- subset(frame, 
                    select = c("domainvalue", "id", "WEIGHT", 
                               paste0("Y", 1:ns_all), 
                               paste0("Y", 1:ns_all, "_SQ_SUM")))
frame_SRS$X1 = 1

SRS_mean_sds <- buildStrataDF(dataset = frame_SRS)
SRS_Pop_CV <- sapply(X = samples,
                     FUN = function(x) {
                       temp_mean <- SRS_mean_sds[, paste0("M", 1:ns_all)]
                       temp_sd <- SRS_mean_sds[, paste0("S", 1:ns_all)]
                       return(temp_sd / sqrt(x) / temp_mean)
                     })
rownames(SRS_Pop_CV) <- sci_names_all

##################################
## Calculate Population CVs under current STRS sampling
##################################
Current_STRS_Pop_CV <- matrix(nrow = ns_all, 
                              ncol = nboats,
                              dimnames = list(sci_names_all, NULL))

for (iboat in 1:nboats) {
  #Adjust sample size proportionally
  nh <- allocations[, paste0('boat', iboat)]
  
  #strata constraints
  # strata_labels <- paste(allocations$Stratum[sampled_strata])
  Nh <- table(Extrapolation_depths$stratum)
  Wh <- Nh / N
  wh <- nh / Nh
  
  #Calculate Strata means and sds (calculated over time as well)
  frame_current <- subset(frame, 
                          select = c("domainvalue", "id", "WEIGHT", 
                                     paste0("Y", 1:ns_all), 
                                     paste0("Y", 1:ns_all, "_SQ_SUM")))
  frame_current$X1 = Extrapolation_depths$stratum
  frame_current$X1[is.na(frame_current$X1)] <- 0
  
  STRS_mean_sds <- buildStrataDF(dataset = frame_current)
  
  STRS_mean <- colSums(sweep(x = STRS_mean_sds[, paste0("M", 1:ns_all)], 
                             MARGIN = 1, 
                             STATS = Wh,
                             FUN = '*'))
  
  STRS_var_temp <- sweep(x = STRS_mean_sds[, paste0("S", 1:ns_all)]^2, 
                         MARGIN = 1, 
                         STATS = Wh^2 * (1 - wh) / nh,
                         FUN = '*')
  
  ## For those strata with zero effort allocated
  STRS_var_temp <- apply(X = STRS_var_temp, 
                         MARGIN = 1:2,
                         FUN = function(x) ifelse(is.finite(x), x, 0))
  
  STRS_var <- colSums(STRS_var_temp)
  
  strata_cv <- sqrt(STRS_var) / STRS_mean 
  
  Current_STRS_Pop_CV[, iboat ] <- strata_cv
}


SS_STRS_Pop_CV <- tidyr::spread(data = settings[,c("iboat", "ispp", "cv")], 
                                value = cv, 
                                key = iboat)[, -1]
rownames(SS_STRS_Pop_CV) = sci_names_opt

##################################
## Save
##################################
save(list = c("SRS_Pop_CV", "Current_STRS_Pop_CV", "SS_STRS_Pop_CV"), 
     file = paste0(github_dir, "results/Population_Variances.RData"))
