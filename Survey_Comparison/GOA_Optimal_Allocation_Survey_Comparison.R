################################################
## Optimization: Lowest CV for a given sample size
################################################
rm(list = ls())

###############################
## Import required packages
###############################
library(sp); library(RColorBrewer); library(raster)

###############################
## Set up directories
###############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[3]
optimization_type = c('_spatial', '_spatiotemporal')[2]
modelno = "6g"

SamplingStrata_dir = paste0(c('/Users/zackoyafuso/',
                              'C:/Users/Zack Oyafuso/',
                              'C:/Users/zack.oyafuso/')[which_machine],
                            'Downloads/SamplingStrata-master/R')

github_dir = paste0(c('/Users/zack.oyafuso/Documents', 
                      'C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/Optimum_Allocation/')

output_wd = paste0(c('/Users/zackoyafuso/Documents/', 
                     'C:/Users/Zack Oyafuso/Documents/',
                     'C:/Users/zack.oyafuso/Work/', 
                     'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                   "GitHub/MS_OM_GoA/Optimum_Allocation/model_", modelno,
                   optimization_type)

if(!dir.exists(output_wd)) dir.create(output_wd)

#########################
## Load functions from SamplingStrata packages into global environment
## Load modified buildStrataDF function if using spatiotemporal modification
#########################
if(optimization_type == '_spatiotemporal'){
  for(ifile in dir(SamplingStrata_dir, full.names = T)) source(ifile)
  source(paste0(github_dir, '/buildStrataDF_Zack.R'))
}

if(optimization_type == '_spatial') library(SamplingStrata)

###########################
## Load Data
###########################
load(paste0(output_wd, '/optimization_data_model_', 
            modelno, '.RData'))
if(optimization_type == '_spatial') rm(frame_raw)

###########################
## Load Current CV Simulation
###########################
load( paste0(output_wd, '/Survey_Simulation_Results.RData') )

CV_constraints = apply(Survey_true_cv_array, MARGIN = 2:3, mean)
CV_constraints = ifelse(CV_constraints < 0.1, 0.1, CV_constraints)
stratas = c(5,10,15,20,25,30,40,50,60)
ns = 15

############################
## Optimizer
############################
par(mfrow = c(6,6), mar = c(2,2,0,0))

for(isample in 3){
  #Create CV dataframe
  cv = list()
  for(spp in 1:ns) 
    cv[[paste0('CV',spp)]] = as.numeric(CV_constraints[spp,isample])
  cv[['DOM']] = 1
  cv[['domainvalue']] = 1
  cv <- as.data.frame(cv)
  
  for(istrata in 1:length(stratas)){   
    
    #Set wd for output files
    temp_dir = paste0(output_wd, '/Survey_Comparison/Str_', 
                      stratas[istrata], 'Boat_', isample, '/')
    if(!dir.exists(temp_dir)) dir.create(temp_dir)
    setwd(temp_dir)
    
    #Run optimization
    solution <- optimStrata(method = "continuous",
                            errors = cv, 
                            framesamp = frame,
                            iter = 200,#ifelse(stratas[istrata] <= 20, 100, 150),
                            pops = 30,
                            elitism_rate = 0.1,
                            mut_chance = 1 / (stratas[istrata] + 1),
                            nStrata = stratas[istrata],
                            showPlot = T,
                            parallel = F,
                            writeFiles = T)
    
    sum_stats = summaryStrata(solution$framenew,
                              solution$aggr_strata,
                              progress=FALSE) 
    
    #Plot Solution
    goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')],
                                 data = cbind(solution$framenew[,paste0('Y',1:ns)],
                                              Str_no = solution$framenew$STRATO,
                                              depth = solution$framenew$X1,
                                              lon = solution$framenew$X2) )
    goa_ras = raster(goa, resolution = 5)
    goa_ras =rasterize(x = goa, y = goa_ras, field = 'Str_no')
    plot(goa_ras, col = terrain.colors(10)[-10], axes = F)
    
    #Save Output
    result_list = list(solution, sum_stats)
    save(list = 'result_list', file = 'result_list.RData')
  }
}
