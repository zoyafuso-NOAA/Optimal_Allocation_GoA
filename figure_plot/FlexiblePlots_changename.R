rm(list = ls())

###############################
## Import required packages
###############################
library(sp); library(RColorBrewer); library(raster)

###############################
## Set up directories
###############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[1]

github_dir = paste0(c('/Users/zackoyafuso/Documents', 
                      'C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/Optimal_Allocation_GoA/')

###########################
## Load Data
###########################
load(paste0(github_dir, 'data/optimization_data.RData'))
stratas = c(5,10,15,20,30,40,50,60)
NStrata = length(stratas)
ns = 15
spp_cv = samplesizes = list()

for(istrata in c(1:4, 8)){
 temp_strata = paste0('Str_', stratas[istrata])
 runs = grep(x = dir(paste0(github_dir, 'Spatiotemporal_Optimization_Scheme2/'), 
                     full.names = T), 
             pattern = paste0('Thres10Str', stratas[istrata]),
             value = T)
 
 nruns = length(runs)
 for(irun in 1:nruns){
  load(paste0(github_dir, 'Spatiotemporal_Optimization_Scheme2/Thres10Str',
              stratas[istrata], 'Run', irun, '/result_list.RData') )
  samplesizes[[temp_strata]]$n = c(samplesizes[[temp_strata]]$n, result_list$n) 
  spp_cv[[temp_strata]]$cv = rbind(spp_cv[[temp_strata]]$cv, result_list[[3]])
 }
}

istrata = 'Str_20'
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

