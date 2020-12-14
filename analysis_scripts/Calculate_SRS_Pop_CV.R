###############################################################################
## Project:       Calculate Population Variances
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Calculate population variances under Simple Random Sampling
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
# load(paste0(github_dir, "data/Extrapolation_depths.RData") )
load(paste0(github_dir, 'data/optimization_data.RData'))

for (idom in c("full_domain", "district") ) {
  
  ## specify data inputs and constants depending on domain type
  frame <- switch( idom,
                   "full_domain" = frame_all,
                   "district" = frame_district)[, c("domainvalue", "id", 
                                                    "WEIGHT",
                                                    paste0("Y", 1:ns_all), 
                                                    paste0("Y", 1:ns_all,
                                                           "_SQ_SUM"))]
  
  n_dom <- length(unique(frame$domainvalue))
  
  ##################################
  ## Calculate Population CVs under Simple Random Sampling
  ##################################
  frame_srs <- subset(frame,
                      select = c("domainvalue", "id", "WEIGHT", 
                                 paste0("Y", 1:ns_all), 
                                 paste0("Y", 1:ns_all, "_SQ_SUM")))
  frame_srs$X1 = frame$domainvalue
  
  srs_mean_sds <- buildStrataDF(dataset = frame_srs)
  
  temp_srs_pop_cv <- array(dim = c(ns_all, n_dom, n_boats),
                           dimnames = list(sci_names_all, NULL, NULL))
  
  for (ispp in 1:ns_all) {
    for (iboat in 1:n_boats) {
      srs_n <- samples[iboat] * table(frame$domainvalue) / n_cells
      srs_var <- srs_mean_sds[, paste0("S", ispp)]^2 * 
        (1 - srs_n / n_cells) / srs_n
      srs_mean <- srs_mean_sds[, paste0("M", ispp)]
      temp_srs_pop_cv[ispp, , iboat] <- sqrt(srs_var) / srs_mean
    }
  }
  
  assign(x = paste0("srs_pop_cv_", idom), value = temp_srs_pop_cv) 
  
  ##################################
  ## Save
  ##################################
  save(list = paste0("srs_pop_cv_", idom), 
       file = paste0(github_dir, "results/", idom, "/srs_pop_cv.RData") )
  
}

