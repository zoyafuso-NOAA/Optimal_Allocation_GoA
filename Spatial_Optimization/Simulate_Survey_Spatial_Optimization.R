######################################
## Simulate Stratified Random Sampling according to the optimized
## stratifications for hte Spatial Only Optimization 
######################################
rm(list = ls())

###############################
## Set up directories
###############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[3]

github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                     'C:/Users/Zack Oyafuso/Documents/',
                     'C:/Users/zack.oyafuso/Work/', 
                     'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                   "GitHub/Optimal_Allocation_GoA/")

#########################
## Load predicted density and optimization results
#########################
load(paste0(github_dir, 'Spatial_Optimization/spatial_only_optimization_results.RData'))
load(paste0(github_dir, 'data/optimization_data.RData'))

#Constants
ids = as.numeric(rownames(frame))
stratas = c(5,10,15,20,30,60)
Nstrata = length(stratas)
Niters = 1000

###########################
## Result Objects
###########################
sim_mean = sim_cv = array(dim = c(NTime, ns, Nstrata, 3, Niters), 
                          dimnames = list(paste0('Year_', 1:NTime),
                                          sci_names, 
                                          NULL, 
                                          NULL))
##########################
## Simulating each optimization
##########################
settings$id = 1:nrow(settings)

for(istrata in 1:Nstrata){
 for(isample in 1:3) {
  
  #Load optimization data
  sub_settings = subset(settings, nstrata == stratas[istrata])
  
  temp_run = sub_settings$id[which.min(abs(sub_settings$n-
                                            c(280,550,820)[isample]))]
  
  strata_allocation = strata_list[[temp_run]]$Allocation
  stratapop = strata_list[[temp_run]]$Population
  stratanos = res_df[,1+temp_run]
  
  #Remove strata with only 1 sample allocated
  str_idx = strata_allocation > 1
  
  for(iyear in 1:NTime){
   for(iter in 1:Niters){
    #Sample based on the stratification allocations
    sample_vec = c()
    for(i in which(str_idx == T) ){
     available_cells = which(stratanos == i)
     sample_cells = sample(x = available_cells, 
                           size = strata_allocation[i], 
                           replace = F)
     sample_vec = c(sample_vec, sample_cells)
    }
    
    #Organize sample set and total number of samples
    sample_vec = sort(sample_vec)
    n = length(sample_vec)
    stratano_samp =  stratanos[sample_vec]
    sample_df = subset(frame_raw, year == iyear)[sample_vec,]
    
    #Calculate Stratum Mean Density and Variance
    stmt = paste0('aggregate(cbind(',
                  paste0('Y', 1:(ns-1), sep=',', collapse=''),'Y',ns, 
                  ") ~ stratano_samp, data = sample_df, FUN = mean)")
    sample_mean = eval(parse(text = stmt))[,-1]
    stmt = paste0('aggregate(cbind(',
                  paste0('Y', 1:(ns-1), sep=',', collapse=''), 'Y',ns, 
                  ") ~ stratano_samp, data = sample_df, FUN = var)")
    sample_var = eval(parse(text = stmt))[,-1]
    
    #How many samples are allocated in each strata
    #How many sampling units are in each strata
    Wh = (stratapop/N)[str_idx]
    wh = (strata_allocation/stratapop)[str_idx]
    
    #Calculate Total Abundance and Variance, calculate CV
    SRS_var = colSums(sweep(x = sample_var, MARGIN = 1, 
                            STATS = (Wh)^2*(1-wh)/
                             strata_allocation[str_idx],
                            FUN = '*'))
    
    SRS_mean = colSums(sweep(x = sample_mean, MARGIN = 1, 
                             STATS = Wh,
                             FUN = '*'))
    
    strata_cv = sqrt(SRS_var) / SRS_mean 
    
    #Record mean and CV values
    sim_mean[paste0('Year_',iyear),,istrata,isample,iter] = SRS_mean
    sim_cv[paste0('Year_',iyear),,istrata,isample,iter] = strata_cv
    
    if(iter%%100 == 0){
     print(paste0('Finished with: Iteration ', iter, ', ', 'Year ', iyear, 
                  ', ', stratas[istrata], ' Strata and ', isample, ' Boat'))
    }
    
   }
  }
  
 }
}

#################################
## Simulation Metrics
#################################
#True CV, Cv of Cv, Rrmse of Cv
true_cv_array = rrmse_est_array = rrmse_cv_array = 
 array(dim = c(NTime, ns, Nstrata, 3), 
       dimnames = list(paste0('Year_', 1:NTime), sci_names, NULL, NULL ))

for(iyear in 1:NTime){
 for(istrata in 1:Nstrata){
  for(isample in 1:3){
   for(spp in sci_names){
    
    iter_est = sim_mean[paste0('Year_', iyear),spp,istrata,isample,]
    iter_cv = sim_cv[paste0('Year_', iyear), spp, istrata,isample, ]
    true_cv = sd(iter_est) / true_mean[iyear, spp]
    
    true_cv_array[paste0('Year_', iyear),spp,istrata,isample] = true_cv
    
    rrmse_cv_array[paste0('Year_', iyear), spp, istrata,isample] = 
     sqrt(mean((iter_cv-true_cv)^2)) / mean(iter_cv)
    
    rrmse_est_array[paste0('Year_', iyear), spp, istrata,isample] = 
     sqrt(mean((iter_est-true_mean[iyear,spp])^2))/
     true_mean[iyear,spp]
   }
  }
 }
 
}

#######################
## Save results
#######################
for(ivar in  c('rrmse_est_array', 'rrmse_cv_array', 'true_cv_array', 
               'sim_mean', 'sim_cv')){
 assign(x=paste0('STRS_', ivar), value = get(ivar))
}

save(file=paste0(github_dir,'Spatial_Optimization/STRS_Sim_Res_Spatial.RData'),
     list = c(paste0('STRS_', c('rrmse_est_array', 'rrmse_cv_array', 
                                'true_cv_array', 'sim_mean', 'sim_cv')),
              'Niters'))

