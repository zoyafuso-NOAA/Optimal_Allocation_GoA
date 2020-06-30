################################################
## Spatiotemporal Constrained CV Survey Optimization
## Single Species
################################################
rm(list = ls())

###############################
## Import required packages
###############################
library(sp); library(RColorBrewer); library(raster)

###############################
## Set up directories
###############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[2]

SamplingStrata_dir = paste0(c('/Users/zackoyafuso/',
                              'C:/Users/Zack Oyafuso/',
                              'C:/Users/zack.oyafuso/')[which_machine],
                            'Downloads/SamplingStrata-master/R')

github_dir = paste0(c('/Users/zack.oyafuso/Documents', 
                      'C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/Optimal_Allocation_GoA/')

#########################
## Load functions from SamplingStrata packages into global environment
## Load modified buildStrataDF function if using spatiotemporal modification
#########################
for(ifile in dir(SamplingStrata_dir, full.names = T)) source(ifile)
source(paste0(github_dir, 'modified_functions/buildStrataDF_Zack.R'))

###########################
## Load Data
###########################
load(paste0(github_dir, 'data/optimization_data.RData'))
load(paste0(github_dir, 'data/Extrapolation_depths.RData'))

###########################
## Load Current CV Simulation
###########################
# stratas = c(5,10,15,20,30,60)
ns = 15
creep_rate = 0.05
threshold = 0.05

master_frame = frame
master_frame_raw = frame_raw

############################
## Optimizer
############################
ispp = 7
frame = master_frame[,c('id', 'X1', 'X2', paste0('Y',ispp), 'domainvalue')]
frame_raw = master_frame_raw[,c('id', 'X1', 'X2', 
                                paste0('Y',ispp), 'domainvalue')]

names(frame) = names(frame_raw) = c('id', 'X1', 'X2', 'Y1', 'domainvalue')
#Initial Condition
Run = 1
CV_constraints = 0.03
current_n = 0

#Create CV dataframe
cv = list()
cv[[paste0('CV',1)]] = CV_constraints
cv[['DOM']] = 1
cv[['domainvalue']] = 1
cv <- as.data.frame(cv)

par(mfrow = c(6,6), mar = c(2,2,0,0))
while(current_n <= 820){
  
  #Set wd for output files
  temp_dir = paste0(github_dir, 'Single_Species_Optimization/', 
                    'Spp', ispp,'Run',Run)
  if(!dir.exists(temp_dir)) dir.create(temp_dir)
  setwd(temp_dir)
  
  #Run optimization
  solution <- optimStrata(method = "continuous",
                          errors = cv, 
                          framesamp = frame,
                          iter = 50,
                          pops = 20,
                          elitism_rate = 0.1,
                          mut_chance = 1 / (10 + 1),
                          nStrata = 10,
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
  plot(goa_ras, axes = F, 
       col = terrain.colors(10)[sample(10)])
  
  #Save Output
  CV_constraints = expected_CV(strata = solution$aggr_strata)
  current_n = sum(sum_stats$Allocation)
  result_list = list(solution, sum_stats, CV_constraints, n = current_n)
  save(list = 'result_list', file = 'result_list.RData')
  
  #Set up next run
  Run = Run + 1
  CV_constraints =  CV_constraints * (1 - creep_rate)
  
  #Create CV dataframe
  cv = list()
  for(spp in 1:1) 
    cv[[paste0('CV',1)]] = as.numeric(CV_constraints)
  cv[['DOM']] = 1
  cv[['domainvalue']] = 1
  cv <- as.data.frame(cv)
}

