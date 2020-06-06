#################################
## Allocation Tables
#################################
#############################
## Show Simulation Metrics for Simulated Survey Strata,
## Simple Random Sampling, and Stratified Random Sampling
#############################
rm(list = ls())

library(sp); library(raster)

############################
## Set up directories
#############################
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2)[1]
modelno = '6g'
github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/')[which_machine],
                    'GitHub/MS_OM_GoA/Optimum_Allocation/', 
                    'model_', modelno, '/')

output_dir = paste0(c('/Users/zackoyafuso/', 
                      'C:/Users/Zack Oyafuso/')[which_machine],
                    'Google Drive/MS_Optimizations/figure_plot/')

load(paste0(github_dir, 'Survey_Simulation_Results.RData'))
load(paste0(github_dir, 'Simple_RS_Simulation_Results.RData'))
load(paste0(github_dir, 'Stratified_RS_Simulation_Results.RData'))
load(paste0(github_dir, 'optimization_results.RData'))

settings$id = 1:nrow(settings)


for(isample in c(800, 550, 280)){
 sub_settings = as.data.frame(t(sapply(X = split(settings, f = settings$nstrata),
                                       FUN = function(x) x[which.min(abs(x$n - isample)),] )))
 
 assign(x = paste0('Allocation_', isample),
        value = t(sapply(sub_settings$id,
          FUN = function(x) quantile(strata_list[[x]]$Allocation))))
 
 assign(x = paste0('Population_', isample),
        value = t(sapply(sub_settings$id,
                         FUN = function(x) quantile(strata_list[[x]]$Population))))
 
 assign(x = paste0('SamplingRate_', isample),
        value = t(sapply(sub_settings$id,
                         FUN = function(x) quantile(strata_list[[x]]$SamplingRate))))
} 

save(list = paste0(c('Allocation_', 'Population_', 'SamplingRate_'), rep(c(800, 550, 280), each = 3)), file=paste0(output_dir, 'allocation_tables.RData') )

