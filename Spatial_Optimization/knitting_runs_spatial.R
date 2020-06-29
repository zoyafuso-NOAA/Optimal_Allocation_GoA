#########################
## Knitting together optimization results
## If method == 'spatial'
#########################
rm(list = ls())

############################
## Set up directories
#############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[1]

output_wd = paste0(c('/Users/zackoyafuso/Documents/', 
                     'C:/Users/Zack Oyafuso/Documents/',
                     'C:/Users/zack.oyafuso/Work/', 
                     'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                   "GitHub/Optimal_Allocation_GoA/")

load(paste0(output_wd, 'data/optimization_data.RData'))
load(paste0(output_wd, 'Spatial_Optimization/optimization_20_strata.RData'))

###################
## Constants
###################
nruns = nrow(settings)

####################
## Knit together strata_list
####################
master_strata_list = list()
for(icol in 1:nruns ){
  master_strata_list[[icol]] = data.frame(strata_list[1:9 + 9*(icol -1)])
}

strata_list = master_strata_list

##########################
## Save 
##########################
save(list = c('res_df', 'settings', 'strata_list'),
     file = paste0(output_wd, 'Spatial_Optimization/',
                   'spatial_only_optimization_results.RData'))

