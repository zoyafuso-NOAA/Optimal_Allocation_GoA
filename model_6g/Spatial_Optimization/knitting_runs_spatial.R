#########################
## Knitting together optimization results
## If method == 'spatial'
#########################
rm(list = ls())

############################
## Set up directories
#############################
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[1]
VAST_model = "6g"
output_wd = paste0(c('/Users/zackoyafuso/Documents/', 
                     'C:/Users/Zack Oyafuso/Documents/',
                     'C:/Users/zack.oyafuso/Work/', 
                     'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                   "GitHub/Optimal_Allocation_GoA/model_", VAST_model, '/')

load(paste0(output_wd, 'optimization_data.RData'))

####################
## Result lIst
####################
master_strata_list = list()
master_settings = data.frame()
master_res_df = data.frame(id = 1:N)

for(istrata in c(5,10,15,20,30,60)){
  load(paste0(output_wd, 'Spatial_Optimization/optimization_',
              istrata, '_strata.RData'))
  
  master_settings = rbind(master_settings, settings)
  
  nruns = nrow(settings)
  temp_strata_list = list()
  
  for(icol in 1:nruns ){
    temp_strata_list[[icol]] = data.frame(strata_list[1:9 + 9*(icol -1)])
  }
  master_strata_list = c(master_strata_list, temp_strata_list)
  
  master_res_df = cbind(master_res_df, res_df[,-1])
}

strata_list = master_strata_list
res_df = master_res_df
settings = master_settings
names[settings]
##########################
## Save 
##########################
save(list = c('res_df', 'settings', 'strata_list'),
     file = paste0(output_wd, 'Spatial_Optimization/',
                   'spatial_only_optimization_results.RData'))

