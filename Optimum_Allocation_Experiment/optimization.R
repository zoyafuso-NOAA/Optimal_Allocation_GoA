######################################
## Survey Optimization
## Compare optimal surveys for each of the four scenarios
## Use lat and lon as the stratum variables
## One Species
## CV = 5%,  5, 10, 20 strata
######################################
rm(list = ls())

library(sp); library(raster)

######################################
## Import OM
######################################
which_machine = c('Zack_Mac' = 1)[1]
output_wd = paste0(c('/Users/zackoyafuso/Documents/')[which_machine],
                   'GitHub/MS_OM_GoA/Optimum_Allocation_Experiment/')
load(paste0(output_wd, '/OM.RData'))

SamplingStrata_dir = paste0(c('/Users/zackoyafuso', 
                              'C:/Users/Zack Oyafuso',
                              'C:/Users/zack.oyafuso',
                              'C:/Users/zack.oyafuso')[which_machine],
                            '/Downloads/SamplingStrata-master/R')
github_dir = paste0(c('/Users/zackoyafuso/Documents', 
                      'C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/Optimum_Allocation/')

#########################
## Load functions from SamplingStrata packages into global environment
## Load modified buildStrataDF function
#########################
for(ifile in dir(SamplingStrata_dir, full.names = T)) source(ifile)
source(paste0(github_dir, '/buildStrataDF_Zack.R'))

########################
## Settings for survey optimization
########################
settings = expand.grid(nstrata = c(5),
                       iscen = 1:nrow(ST_settings),
                       cv = 0.05)
res_df = data.frame(id = 1:nrow(frame))
strata_list = list()

par(mfrow = c(4,2))
for(irow in 1:nrow(settings)){
 ##indices
 iscen = settings$iscen[irow]
 istrata = settings$nstrata[irow]
 current_CV = settings$cv[irow]
 
 ##DFs for the optimizer
 frame = data.frame(id = 1:nrow(OM[[iscen]]$loc_xy),
                    X1 = OM[[iscen]]$loc_xy$x,
                    X2 = OM[[iscen]]$loc_xy$y,
                    Y1 = rowMeans(OM[[iscen]]$log_d_rt),
                    domainvalue = 1)
 
 frame_raw = data.frame(id = rep(1:nrow(OM[[iscen]]$loc_xy), times = 10),
                        X1 = rep(OM[[iscen]]$loc_xy$x, times = 10),
                        X2 = rep(OM[[iscen]]$loc_xy$y, times = 10),
                        Y1 = as.vector(OM[[iscen]]$log_d_rt),
                        domainvalue = 1)
 
 #Create CV dataframe
 cv = list()
 cv[['CV1']] = current_CV
 cv[['DOM']] = 1
 cv[['domainvalue']] = 1
 cv <- as.data.frame(cv)
 
 #Run optimization
 solution <- optimStrata(method = "continuous",
                         errors = cv, 
                         framesamp = frame,
                         iter = 50,
                         pops = 30,
                         elitism_rate = 0.1,
                         mut_chance = 1 / (istrata + 1),
                         nStrata = istrata,
                         showPlot = F,
                         parallel = F)
 
 sum_stats = list(summaryStrata(solution$framenew,
                           solution$aggr_strata,
                           progress=FALSE) )
 
 #Plot Solution
 domain = SpatialPointsDataFrame(coords = OM[[iscen]]$loc_xy,
                                 data = data.frame(Str_no = solution$framenew$STRATO) )
 domain_ras = raster(domain, resolution = 1/100)
 goa_ras =rasterize(x = domain, y = domain_ras, field = 'Str_no')
 plot(goa_ras, col = terrain.colors(10)[-10], axes = F)
 
 #Update settings, res_df, and strata_list
 settings$n[irow] = sum(sum_stats[[1]]$Allocation)
 strata_list = c(strata_list, sum_stats)
 res_df = cbind(res_df, solution$indices$X1)
}

save(list = c('res_df', 'settings', 'strata_list'), 
     file = paste0(output_wd, 'optimization.RData'))
