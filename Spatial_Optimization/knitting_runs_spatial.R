#########################
## Knitting together optimization results
## If method == 'spatial'
#########################
rm(list = ls())

############################
## Set up directories
#############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[2]
modelno = "6g"
optimization_type = c('_spatial', '_spatiotemporal')[1]

output_wd = paste0(c('/Users/zackoyafuso/Documents/', 
                     'C:/Users/Zack Oyafuso/Documents/',
                     'C:/Users/zack.oyafuso/Work/', 
                     'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                   "GitHub/MS_OM_GoA/Optimum_Allocation/model_", modelno,
                   optimization_type)

load(paste0(output_wd, '/optimization_data_model_', modelno, '.RData'))
load(paste0(output_wd, '/optimization_60_strata.RData'))

####################
## Knit together strata_list
####################
master_strata_list = list()
for(icol in 1:(ncol(res_df)-1) ){
  master_strata_list[[icol]] = data.frame(strata_list[1:9 + 9*(icol -1)])
}

strata_list = master_strata_list

save(list = c('res_df', 'settings', 'strata_list'),
     file = paste0(output_wd, '/optimization_results.RData'))

