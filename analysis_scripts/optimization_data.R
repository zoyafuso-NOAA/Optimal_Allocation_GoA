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
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents", 
                       "C:/Users/Zack Oyafuso/Documents",
                       "C:/Users/zack.oyafuso/Work")[which_machine],
                     "/GitHub/Optimal_Allocation_GoA/")

##################################################
####   Load VAST output and bathymetry data
##################################################
load(paste0(github_dir,  "data/fit_density.RData"))
load(paste0(github_dir,  "data/fit_Index.RData"))
load(paste0(github_dir, "/data/Extrapolation_depths.RData"))

##################################################
####   Constants
##################################################
Year_Set <- 1996:2019
Years2Include <- c(1,  4,  8, 10, 12, 14, 16, 18, 20, 22, 24)
NTime <- length(Years2Include)

#Number of sampling grids
N <- nrow(Extrapolation_depths)

#Species names
sci_names <- c("Atheresthes stomias", "Gadus chalcogrammus",
               "Gadus macrocephalus", "Glyptocephalus zachirus",
               "Hippoglossoides elassodon", "Hippoglossus stenolepis", 
               "Lepidopsetta bilineata", "Lepidopsetta polyxystra",
               "Microstomus pacificus", "Sebastes alutus", "Sebastes B_R",
               "Sebastes brevispinis", "Sebastes polyspinis", 
               "Sebastes variabilis", "Sebastolobus alascanus" )

common_names <- c("arrowtooth flounder", "Alaska pollock", "Pacific cod",
                  "rex sole", "flathead sole", "Pacific halibut", 
                  "southern rock sole", "northern rock sole", 
                  "Pacific Dover sole", "Pacific ocean perch", 
                  "blackspotted/rougheye rockfishes", "silvergrey rockfish",
                  "northern rockfish", "dusky rockfish",
                  "shortspine thornyhead")
  
ns <- length(sci_names)

#Sample sizes
samples <- c(280, 550, 820)
nboats <- length(samples)

#Number of times to simulate survey
Niters <- 1000

stratas <- c(5, 10, 15, 20, 30, 60)
NStrata <- length(stratas)

##################################################
####   Empty result objects
##################################################
df <- df_raw <- NULL

##################################################
####   Mean density across years, full domain
##################################################
df <- cbind(
  data.frame(Domain = 1,
             x = 1:N,
             lon = with(Extrapolation_depths, E_km - min(E_km)),
             depth = Extrapolation_depths$DEPTH_EFH),
  
  #Mean Density across years
  apply(X = D_gct[,,Years2Include], 
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

##################################################
####   Predicted density for each observed year and cell, Full domain
##################################################
for (iT in 1:NTime) {
  df_raw <- rbind(df_raw, cbind(
    data.frame(Domain = 1,
               x = 1:N,
               year = iT,
               lon = with(Extrapolation_depths, E_km - min(E_km)),
               depth = Extrapolation_depths$DEPTH_EFH),
    D_gct[,,Years2Include[iT]] )
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

frame_raw$year <- rep(x = 1:NTime, each = N)

##################################################
####   Calculate "true" mean density and "true" abundance index
##################################################
stmt <- paste0("aggregate(cbind(",
               paste0("Y", 1:(ns-1), sep = ",", collapse = ""), "Y",ns, 
               ") ~ year, data = frame_raw, FUN = mean)")
true_mean <- eval(parse(text = stmt))[,-1]
colnames(true_mean) <- sci_names

true_index <- t(apply(X = Index[,, Years2Include], 
                      MARGIN = 2:3,
                      FUN = sum))

##################################################
####   Save Data
##################################################
save(list = c("frame", "frame_raw", "true_mean", "true_index", "ns", 
              "Years2Include","NTime", "N", "sci_names", "common_names",
              "samples", "nboats", "Niters", "stratas", "NStrata"),
     file = paste0(github_dir, "data/optimization_data.RData"))
