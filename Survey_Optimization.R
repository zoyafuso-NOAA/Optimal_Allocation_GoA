###############################################################################
## Project:       Spatiotemporal Survey Optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Conduct SamplingStrata R package multispecies stratified
##                survey optimization
###############################################################################

##################################################
####    Import required packages
##################################################
library(sp); library(RColorBrewer); library(raster)

##################################################
####   Set up directories
##################################################
rm(list = ls())
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[3]
VAST_model = "11" 

SamplingStrata_dir = paste0(c('/Users/zackoyafuso/',
                              'C:/Users/Zack Oyafuso/',
                              'C:/Users/zack.oyafuso/')[which_machine],
                            'Downloads/SamplingStrata-master/R')

github_dir = paste0(c('/Users/zack.oyafuso/Documents', 
                      'C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/Optimal_Allocation_GoA/model_', VAST_model, '/')

##################################################
####   Load functions from SamplingStrata packages into global environment
####   Load modified buildStrataDF function if using spatiotemporal 
####   stratum variance
##################################################
for(ifile in dir(SamplingStrata_dir, full.names = T)) source(ifile)
source(paste0(dirname(github_dir), '/modified_functions/buildStrataDF_Zack.R'))

##################################################
####   Load Data
##################################################
load(paste0(github_dir, 'optimization_data.RData'))
load(paste0(dirname(github_dir), '/data/Extrapolation_depths.RData'))
load(paste0(github_dir, 'Survey_Comparison_Simulations/',
            'Survey_Simulation_Results.RData'))

##################################################
####   Set up some constants of the optimization
####   Constrained: One CV constraint applied to all species
####   Flexible: species specific CV constraints
####   Flexible_wo_RockF: Flexible without Sebastses polyspinnis and 
####                      Sebastes variabilis
####   Single_Species: univariate optimization, one CV constraint
##################################################
which_method = c('Flexible' = 1,
                 'Flexible_wo_RockF' = 2,
                 'Single_Species' = 3)[1]

stratas = c(5,10,15,20,30,60)

ns = c(15, 13, 1)[which_method]

dirname = paste0(github_dir, 
                 c('Spatiotemporal_Optimization',
                   'Spatiotemporal_Optimization_noRFs',
                   'Single_Species_Optimization')[which_method], '/')
if(!dir.exists(dirname)) dir.create(dirname)

##################################################
####   If Flexible_wo_RockF: take the two species out of the dataframes   
####   If Single_Species: subset just the one species
##################################################
SS_which_species = 1 #which species are we doing?
if(which_method == 3){
  SS_which_species = 1 #which species are we doing?
  frame = frame[,c('id', 'X1', 'X2', paste0('Y', SS_which_species),
                   'domainvalue')]
  
  frame_raw = frame_raw[,c('id', 'X1', 'X2', paste0('Y', SS_which_species),
                           'domainvalue', 'year')]
  
  names(frame)[4] = names(frame_raw)[4] = 'Y1'
}

##################################################
####  lower CV threshold
##################################################
threshold = list(apply(Survey_true_cv_array, 
                       MARGIN = 2:3, FUN = median),
                 apply(Survey_true_cv_array, 
                       MARGIN = 2:3, FUN = median),
                 matrix(0, nrow = ns, ncol = 3))[[which_method]]


##################################################
####   Run optimization
##################################################
par(mfrow = c(6,6), mar = c(2,2,0,0))

for(istrata in 2){
  
  temp_strata = stratas[istrata]
  
  ##Initial Condition
  Run = 1
  isample = 1
  current_n = 0
  
  ##Initial Upper CV constraints
  CV_constraints = list( rep(c(.4, 0.3, 0.2)[isample], ns),
                         rep(c(.4, 0.3, 0.2)[isample], ns),
                         rep(c(.4, 0.3, 0.2)[isample], ns),
                         c(0.09, 0.20, 0.10, 
                           0.09, 0.15, 0.07,
                           0.05, 0.09, 0.15,
                           0.09, 0.20, 0.06,
                           0.30, 0.20, 0.06)[SS_which_species])[[which_method]]
  
  creep_rate = c(0.1, 0.05, 0.025)[isample]
  
  #Create CV dataframe
  cv = list()
  for(spp in 1:ns) 
    cv[[paste0('CV',spp)]] = as.numeric(CV_constraints[spp])
  cv[['DOM']] = 1
  cv[['domainvalue']] = 1
  cv <- as.data.frame(cv)
  
  while(current_n <= 820){ #Run until you reach 820 samples
    
    #Set wd for output files, create a directory if it doesn't exist yet
    temp_dir = paste0(dirname, 'Str', temp_strata, 'Run',Run)
    if(!dir.exists(temp_dir)) dir.create(temp_dir)
    
    setwd(temp_dir)
    
    #Run optimization
    solution <- optimStrata(method = "continuous",
                            errors = cv, 
                            framesamp = frame,
                            iter = 300,
                            pops = 30,
                            elitism_rate = 0.1,
                            mut_chance = 1 / (temp_strata + 1),
                            nStrata = temp_strata,
                            showPlot = T,
                            parallel = F,
                            writeFiles = T)
    
    sum_stats = summaryStrata(solution$framenew,
                              solution$aggr_strata,
                              progress=FALSE) 
    
    #Plot Solution
    goa = SpatialPointsDataFrame(
      coords = Extrapolation_depths[,c('E_km', 'N_km')],
      data = data.frame(Str_no = solution$framenew$STRATO) )
    goa_ras = raster(goa, resolution = 5)
    goa_ras =rasterize(x = goa, y = goa_ras, field = 'Str_no')
    
    png(filename = 'solution.png', width = 5, height = 5, units = 'in', 
        res = 500)
    plot(goa_ras, axes = F, 
         col = terrain.colors(temp_strata)[sample(temp_strata)])
    dev.off()
    
    #Save Output
    CV_constraints = expected_CV(strata = solution$aggr_strata)
    current_n = sum(sum_stats$Allocation)
    isample = ifelse(current_n < 280, 1, #1 boat
                     ifelse(current_n < 550, 2, #2 boat
                            3)) #3 boat
    result_list = list(solution, sum_stats, CV_constraints, n = current_n)
    save(list = 'result_list', file = 'result_list.RData')
    
    #Set up next run by changing upper CV constraints
    Run = Run + 1
    creep_rate = c(0.1, 0.05, 0.05)[isample]
    CV_constraints = CV_constraints * (1 - creep_rate) 
    for(ispp in 1:ns){
      CV_constraints[ispp]=ifelse(CV_constraints[ispp]<threshold[ispp,isample],
                                  threshold[ispp, isample],
                                  CV_constraints[ispp])
    }
    
    #Create CV dataframe in the formmat of SamplingStrata
    cv = list()
    for(spp in 1:ns) 
      cv[[paste0('CV',spp)]] = as.numeric(CV_constraints[spp])
    cv[['DOM']] = 1
    cv[['domainvalue']] = 1
    cv <- as.data.frame(cv)
  }
}
