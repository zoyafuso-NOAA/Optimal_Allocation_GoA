###############################################################################
## Project:       Data synthesis for stratified survey optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Create dataset used for all optimization runs based on a 
##                Gulf of Alaska Multispecies Groundfish VAST Spatiotemporal 
##                Operating Model specified on line 24
##
##                Calculate true mean density across years for each species
##
##                Set up other constants used in downstream processes
###############################################################################
rm(list = ls())

##################################################
####  Import required packages   
##################################################
library(SamplingStrata)

##################################################
####    Set up directories
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[3]
VAST_model <- "11" 

github_dir <- paste0(c("/Users/zackoyafuso/Documents", 
                       "C:/Users/Zack Oyafuso/Documents",
                       "C:/Users/zack.oyafuso/Work")[which_machine],
                     "/GitHub/Optimal_Allocation_GoA/")

VAST_dir <- paste0(c("/Users/zackoyafuso/Google Drive/GOA_", 
                     "C:/Users/Zack Oyafuso/Google Drive/GOA_", 
                     "G:/Oyafuso/")[which_machine],
                   "VAST_Runs_EFH/VAST_output", VAST_model, "/")

##################################################
####   Set up Result Directories
##################################################
result_dir <- paste0(github_dir, "model_", VAST_model, "/")

if(!dir.exists(result_dir)){
  dir.create(result_dir)
  dir.create(paste0(result_dir, "Survey_Comparison_Simulations/"))
  dir.create(paste0(result_dir, "Spatiotemporal_Optimization/"))
  dir.create(paste0(result_dir, "Single_Species_Optimization/"))
}

##################################################
####   Load VAST output and bathymetry data
##################################################
load(paste0(VAST_dir, "/fit.RData"))

load(paste0(github_dir, "data/Extrapolation_depths.RData"))
spp_df <- read.csv(file = paste0(github_dir, "data/spp_df.csv"), 
                   check.names = F, 
                   header = T, 
                   row.names = "modelno")

##################################################
####   Constants
##################################################

Year_Set <- 1996:2019
Years2Include <- c(1,  4,  8, 10, 12, 14, 16, 18, 20, 22, 24)
NTime <- length(Years2Include)

#Number of sampling grids
N <- nrow(Extrapolation_depths)

#Species names
which_spp_idx <- unlist(spp_df[VAST_model,])
sci_names <- sort( names(spp_df)[which_spp_idx] )

ns <- length(sci_names)

#Sample sizes
samples <- c(280, 550, 820)
nboats <- length(samples)

#Number of times to simulate survey
Niters <- 1000

stratas <- c(5, 10, 15, 20, 30, 60)
NStrata <- length(stratas)

##################################################
####   Create indices for trawlable and shallow cells
##################################################
trawl_idx <- Extrapolation_depths$Id %in% cells_trawlable
shallow_idx <- Extrapolation_depths$Id %in% cells_shallower_than_700m
trawl_shallow_idx <- apply(X = cbind(trawl_idx, shallow_idx),
                           MARGIN = 1,
                           FUN = all)

##################################################
####   Create the data inputs to SamplingStrata
##################################################
df <- df_trawl <- df_raw <- df_raw_trawl <- NULL 

##################################################
####   Mean density across years, full domain
##################################################
df <- cbind(
  data.frame(Domain = 1,
             x = 1:N,
             lon = Extrapolation_depths$E_km - min(Extrapolation_depths$E_km),
             depth = Extrapolation_depths$DEPTH_EFH),
  
  #Mean Density across years
  apply(X = fit$Report$D_gct[,,Years2Include], 
        MARGIN = 1:2, 
        FUN = mean )
)
names(df)[-(1:4)] <- gsub(x = sci_names, pattern = " ", replacement = "_")

frame <- SamplingStrata::buildFrameDF(df = df,
                                      id = "x",
                                      X = c("depth", "lon"),
                                      Y = gsub(x = sci_names, 
                                               pattern = " ", 
                                               replacement = "_"),
                                      domainvalue = "Domain")

#################################################
####   Mean density across years, trawlable areas
##################################################
df_trawl <- cbind(
  data.frame(Domain = 1,
             x = (1:N)[trawl_shallow_idx],
             lon = Extrapolation_depths$E_km[trawl_shallow_idx] - 
               min(Extrapolation_depths$E_km[trawl_shallow_idx]),
             depth = Extrapolation_depths$DEPTH_EFH[trawl_shallow_idx]),
  
  #Mean Density across years
  apply(X = fit$Report$D_gct[trawl_shallow_idx,,Years2Include], 
        MARGIN = 1:2, 
        FUN = mean )
)
names(df_trawl)[-(1:4)] <- gsub(x = sci_names, pattern = " ", replacement = "_")

frame_trawl <- SamplingStrata::buildFrameDF(df = df_trawl,
                                            id = "x",
                                            X = c("depth", "lon"),
                                            Y = gsub(x = sci_names, 
                                                     pattern = " ", 
                                                     replacement = "_"),
                                            domainvalue = "Domain")

##################################################
####   Predicted density for each observed year and cell, Full domain
##################################################
for (iT in 1:NTime) {
  df_raw <- rbind(df_raw, cbind(
    data.frame(Domain = 1,
               x = 1:N,
               year = iT,
               lon = Extrapolation_depths$E_km - min(Extrapolation_depths$E_km),
               depth = Extrapolation_depths$DEPTH_EFH),
    fit$Report$D_gct[,,Years2Include[iT]] )
  )
}
names(df_raw)[-(1:5)] <- gsub(x = sci_names, 
                              pattern = " ",
                              replacement = "_")

frame_raw <- SamplingStrata::buildFrameDF(df = df_raw,
                                          id = "x",
                                          X = c("depth", "lon"),
                                          Y = gsub(x = sci_names, 
                                                   pattern = " ", 
                                                   replacement = "_"),
                                          domainvalue = "Domain")

##################################################
####   Predicted density for each observed year and cell, trawlable
##################################################
for (iT in 1:NTime) {
  df_raw_trawl <- rbind(df_raw_trawl, cbind(
    data.frame(Domain = 1,
               x = (1:N)[trawl_shallow_idx],
               year = iT,
               lon = Extrapolation_depths$E_km[trawl_shallow_idx] - 
                 min(Extrapolation_depths$E_km[trawl_shallow_idx]),
               depth = Extrapolation_depths$DEPTH_EFH[trawl_shallow_idx]),
    fit$Report$D_gct[trawl_shallow_idx,,Years2Include[iT]] )
  )
}
names(df_raw_trawl)[-(1:5)] <- gsub(x = sci_names, 
                                    pattern = " ", 
                                    replacement = "_")

frame_raw_trawl <- SamplingStrata::buildFrameDF(df = df_raw_trawl,
                                          id = "x",
                                          X = c("depth", "lon"),
                                          Y = gsub(x = sci_names, 
                                                   pattern = " ", 
                                                   replacement = "_"),
                                          domainvalue = "Domain")

##################################################
####   Calculate "true" mean density and "true" abundance index
##################################################
frame_raw$year <- rep(x = 1:NTime, 
                      each = N)
frame_raw_trawl$year <- rep(x = 1:NTime, 
                            each = sum(trawl_shallow_idx))

stmt <- paste0("aggregate(cbind(",
               paste0("Y", 1:(ns-1), sep = ",", collapse = ""), "Y",ns, 
               ") ~ year, data = frame_raw, FUN = mean)")
true_mean <- eval(parse(text = stmt))[,-1]
colnames(true_mean) <- sci_names

true_index <- t(apply(X = fit$Report$Index[,, Years2Include], 
                      MARGIN = 2:3,
                      FUN = sum))

##################################################
####   Save Data
##################################################
save(list = c("frame", "frame_trawl", "frame_raw", "frame_raw", "true_mean", 
              "true_index", "ns", "Years2Include","NTime", "N", "sci_names", 
              "samples", "nboats", "Niters", "stratas", "NStrata"),
     file = paste0(github_dir, "model_", VAST_model, 
                   "/optimization_data.RData"))
