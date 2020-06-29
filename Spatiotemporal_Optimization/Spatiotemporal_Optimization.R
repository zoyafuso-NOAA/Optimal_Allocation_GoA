################################################
## Spatiotemporal Constrained CV Survey Optimization
################################################
rm(list = ls())

###############################
## Import required packages
###############################
library(sp); library(RColorBrewer); library(raster)

###############################
## Set up directories
###############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[1]

SamplingStrata_dir = paste0(c('/Users/zackoyafuso/',
                              'C:/Users/Zack Oyafuso/',
                              'C:/Users/zack.oyafuso/')[which_machine],
                            'Downloads/SamplingStrata-master/R')

github_dir = paste0(c('/Users/zackoyafuso/Documents', 
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

############################
## Settings for optimizer
############################
nstrata = c(5,10,15,20,25,30,40,50,60)
settings = data.frame()
res_df = data.frame(id = 1:N)
strata_list = list()

par(mfrow = c(4,3), mar = c(0,0,0,0))
for(istrata in nstrata){
  
  #Initial Condition
  current_n = 10000
  current_CV = 0.08
  
  while(current_n > 280){
    
    #Create CV dataframe
    cv = list()
    for(spp in 1:ns) cv[[paste0('CV', spp)]] = current_CV
    cv[['DOM']] = 1
    cv[['domainvalue']] = 1
    cv <- as.data.frame(cv)
    
    #Run optimization
    solution <- optimStrata(method = "continuous",
                            errors = cv, 
                            framesamp = frame,
                            iter = 200,
                            pops = 50,
                            elitism_rate = 0.1,
                            mut_chance = 1 / (istrata + 1),
                            nStrata = istrata,
                            showPlot = F,
                            parallel = F)
    
    sum_stats = summaryStrata(solution$framenew,
                              solution$aggr_strata,
                              progress=FALSE) 
    
    #Update settings, res_df, and strata_list
    settings = rbind(settings, data.frame(nstrata = istrata,
                                          cv = current_CV,
                                          n = sum(sum_stats$Allocation)))
    strata_list = c(strata_list, sum_stats)
    res_df = cbind(res_df, solution$indices$X1)
    
    #Calculate total sample size, used to determine whether to stop 
    #optimization or repeat the optimization with a higher CV constraint
    current_n = sum(sum_stats$Allocation)
    
    #Output the results of the optimzation to the console
    plot(1, type = 'n', xlim=c(0,5), ylim = c(0,5), axes = F, ann = F); box()
    text(x = 2.5, y = 2.5, paste0("Just Saved:\n", istrata, ' strata,\n',
                                  current_CV*100, '% CV,\n', current_n, 
                                  ' sample size'))
    
    #Update the next CV level
    current_CV = current_CV + 0.01
    
    #Plot Solution
    goa = SpatialPointsDataFrame(
      coords = Extrapolation_depths[,c('E_km', 'N_km')],
      data = data.frame(Str_no = solution$framenew$STRATO) )
    goa_ras = raster(goa, resolution = 5)
    goa_ras =rasterize(x = goa, y = goa_ras, field = 'Str_no')
    plot(goa_ras, col = terrain.colors(10)[-10], axes = F)
    
    #Save Output
    save(list = c('res_df', 'strata_list', 'settings'),
         file = paste0(github_dir, 'Spatiotemporal_Optimization/optimization_',
                       istrata,'_strata.RData'))
  }
}