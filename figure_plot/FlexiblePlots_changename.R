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

github_dir = paste0(c('/Users/zackoyafuso/Documents', 
                      'C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/Optimum_Allocation/')

output_wd = paste0(c('/Users/zackoyafuso/Documents/', 
                     'C:/Users/Zack Oyafuso/Documents/',
                     'C:/Users/zack.oyafuso/Work/', 
                     'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                   "GitHub/MS_OM_GoA/Optimum_Allocation/model_", modelno,
                   optimization_type, '/Flexible_Optimization/')

###########################
## Load Data
###########################
load(paste0(github_dir, "model_", modelno,
            optimization_type,'/optimization_data_model_', 
            modelno, '.RData'))
if(optimization_type == '_spatial') rm(frame_raw)

stratas = c(5,10,15,20,25,30,40,50,60)
NStrata = length(stratas)
ns = 15
spp_cv = samplesizes = list()

for(istrata in c(1:3, 9)){
 temp_strata = paste0('Str_', stratas[istrata])
 runs = grep(x = dir(output_wd, full.names = T), 
             pattern = paste0('Thres10Str', stratas[istrata]),
             value = T)
 
 nruns = length(runs)
 for(irun in 1:nruns){
  load(paste0(output_wd, 'Thres10Str', stratas[istrata], 'Run', 
              irun, '/result_list.RData')
  )
  samplesizes[[temp_strata]]$n = c(samplesizes[[temp_strata]]$n, result_list$n) 
  spp_cv[[temp_strata]]$cv = rbind(spp_cv[[temp_strata]]$cv, result_list[[3]])
 }
}

istrata = 'Str_10'
nruns = length(samplesizes[[istrata]]$n)
spp_order = order(spp_cv[[istrata]]$cv[nruns,])
run_order = order(samplesizes[[istrata]]$n)
par(mar = c(5,5,2,1), mfrow = c(1,2))

matplot( t(spp_cv[[istrata]]$cv[,spp_order]), 
         type = 'l', lty = 1, las = 1, xlab = 'Species', ylim = c(0,0.3),
         ylab = 'Expected Spatiotemporal CV',
         col = 'black')
abline(h = 0.1, col = 'darkgrey', lty = 'dashed')
text(x = 1:ns, y = t(spp_cv[[istrata]]$cv[,spp_order]), 
     rep(paste(1:nruns),each = ns ))
         
plot(samplesizes[[istrata]]$n[run_order], type = 'l', cex = 2,
     xlab = 'Run Number', ylab = 'Total Sample Size', las = 1, ylim = c(0,850))
abline(h = c(280, 550, 820), col = 'darkgrey', lty = 'dashed')
text(1:nruns, samplesizes[[istrata]]$n[run_order], paste(1:nruns))
#mtext(side = 3, istrata, outer = T, line = -2)

sort( sapply(spp_cv[[istrata]]$cv[nruns,spp_order], 
             FUN = function(x) max(0.95*x, 0.1)) )


#Plot Solution
par(mfrow = c(1,1), mar = c(0,0,0,0))
istrata = 1
temp_n = samplesizes[[istrata]]$n
idx = which.min(abs(temp_n - 280))

load(paste0(output_wd, 'Thres10Str', stratas[istrata], 'Run', 
            idx, '/result_list.RData'))

goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')],
                             data = data.frame(Str_no = result_list[[1]]$indices$X1) )
goa_ras = raster(goa, resolution = 5)
goa_ras =rasterize(x = goa, y = goa_ras, field = 'Str_no')
plot(goa_ras, col = brewer.pal(n = 10, name = 'Paired'), axes = F)

