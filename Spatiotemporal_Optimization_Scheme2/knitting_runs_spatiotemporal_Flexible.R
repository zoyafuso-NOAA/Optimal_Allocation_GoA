#########################
## Knitting together optimization results for the Flexible (Scheme 2)
## Spatiotemporal Optimization
#########################

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

######
stratas = c(5,10,15,20,25,30,40,50,60)
NStrata = length(stratas)

settings = data.frame()
strata_list = list()
res_df = data.frame(id = 1:N)

for(istrata in c(1:4, 9)){
 temp_strata = paste0('Str_', stratas[istrata])
 runs = grep(x = dir(paste0(github_dir, 'Spatiotemporal_Optimization_Scheme2'),
                     full.names = T), 
             pattern = paste0('Thres10Str', stratas[istrata]),
             value = T)
 
 nruns = length(runs)
 for(irun in 1:nruns){
  temp_file = paste0(github_dir, 
                     'Spatiotemporal_Optimization_Scheme2/Thres10Str',
                     stratas[istrata], 'Run', irun, '/result_list.RData')
  
  if(file.exists(temp_file)){
   load(temp_file)
   
   spp_cvs = result_list[[3]]
   colnames(spp_cvs) = paste0('CV_', 1:ns)
   rownames(spp_cvs) = NULL
   
   settings = rbind(settings, 
                    cbind(data.frame(nstrata = stratas[istrata], 
                                     n = result_list$n), 
                          spp_cvs))
   
   strata_list = c(strata_list, list(result_list[[2]]))
   
   res_df = cbind(res_df, result_list[[1]]$indices$X1)
  }
 }
}

settings$id = 1:nrow(settings)
names(res_df)[-1] = paste0('sol_', 1:nrow(settings))

########################
## Save
########################
save(list = c('res_df', 'settings', 'strata_list'),
     file = paste0(github_dir, 'Spatiotemporal_Optimization_Scheme2/',
                   'spatiotemporal_Flexible_optimization_results.RData'))
