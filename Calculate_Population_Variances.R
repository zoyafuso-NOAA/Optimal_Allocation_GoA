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
####   Import Libraries
##################################################
library(VAST)
library(readxl)
library(rgdal)
library(spatialEco)

##################################################
####   Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[3]
VAST_model <- "11" 

github_dir <- paste0(c('/Users/zackoyafuso/Documents/', 
                       'C:/Users/Zack Oyafuso/Documents/',
                       'C:/Users/zack.oyafuso/Work/')[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/")

##################################
## Import Operating Model
##################################
load("G:/Oyafuso/VAST_Runs_EFH/VAST_output11/fit.RData")
load(paste0(dirname(github_dir), "/data/Extrapolation_depths.RData") )
load(paste0(github_dir, 'optimization_data.RData'))

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
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
allocations$boat1 = ceiling(allocations$boat2 / 2)

allocations$boat1 = ifelse(allocations$boat1 == 0, 0, 
                           ifelse(allocations$boat1 == 1, 2, 
                                  allocations$boat1))

##################################################
####   Attribute each grid point to the current stratification
##################################################
goa_grid <- rgdal::readOGR(
  "C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/data/shapefiles/goa_strata.shp")
goa_grid = sp::spTransform(x = goa_grid, 
                           CRSobj = "+proj=utm +zone=5N +units=km")
temp <- spatialEco::point.in.poly(
  x = SpatialPoints(coords = Extrapolation_depths[, c("E_km", "N_km")],
                    proj4string=CRS("+proj=utm +zone=5N +units=km")),
  y = goa_grid)

plot(subset(goa_grid, STRATUM != 0), 
     col = terrain.colors(100)[sample(60)],
     border = F)
plot(subset(goa_grid, STRATUM == 0), col = "brown", add = T)

points(temp, 
       pch = 16,
       cex = 0.5)

points(subset(temp, STRATUM == 0), col = 'red')

##################################
## Calculate Population CVs under Simple Random Sampling
##################################
SRS_Pop_CV <- Current_STRS_Pop_CV <- matrix(nrow = ns, 
                                            ncol = nboats,
                                            dimnames = list(sci_names, NULL))

for (ispp in 1:15) {
  SRS_var = var(as.vector(fit$Report$D_gct[, ispp, Years2Include])) 
  SRS_mean = mean(as.vector(fit$Report$D_gct[, ispp, Years2Include]))
  
  SRS_CV = sqrt(SRS_var / samples) / SRS_mean
  
  SRS_Pop_CV[ispp, ] <- SRS_CV
}

##################################
## Calculate Population CVs under current STRS sampling
##################################

for (isample in 1:nboats) { 
  #Adjust sample size proportionally
  nh <- allocations[, paste0('boat', isample)]
  sampled_strata <- nh > 0
  nstrata <- sum(sampled_strata)
  
  nh <- nh[sampled_strata]
  
  #strata constraints
  Nh <- table(temp$STRATUM)[paste(allocations$Stratum[sampled_strata])]
  Wh <- Nh / N
  wh <- nh / Nh
  
  #Calculate Total Mean, Variance, CV
  stratano = rep(temp@data$STRATUM, 11)
  
  stmt <- paste0('aggregate(cbind(',
                 paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                 ") ~ stratano, data = frame_raw, FUN = mean)")
  strata_mean <- eval(parse(text = stmt))
  strata_mean <- subset(strata_mean, stratano %in% allocations$Stratum)[,-1]
  strata_mean <- strata_mean[sampled_strata,]
  
  stmt <- paste0('aggregate(cbind(',
                 paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                 ") ~ stratano, data = frame_raw, FUN = var)")
  strata_var <- eval(parse(text = stmt))
  strata_var <- subset(strata_var, stratano %in% allocations$Stratum)[,-1]
  strata_var <- strata_var[sampled_strata,]
  
  #Calculate Total Mean, Variance, CV
  SRS_var <- colSums(sweep(x = strata_var, 
                           MARGIN = 1, 
                           STATS = Wh^2 * (1 - wh) / nh,
                           FUN = '*'))
  
  SRS_mean <- colSums(sweep(x = strata_mean, 
                            MARGIN = 1, 
                            STATS = Wh,
                            FUN = '*'))
  
  strata_cv <- sqrt(SRS_var) / SRS_mean 
  
  Current_STRS_Pop_CV[, isample ] <- strata_cv
}

##################################
## Save
##################################
save(list = c("SRS_Pop_CV", "Current_STRS_Pop_CV"), 
     file = paste0(github_dir, "Population_Variances.RData"))
