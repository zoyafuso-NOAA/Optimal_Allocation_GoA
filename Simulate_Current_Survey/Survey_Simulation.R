##################################
## Simulate a Stratified Random Survey of the Gulf of Alaska Groundfish Survey 
## Based on Current Stratifications
##################################
rm(list = ls())

##################################
## Set up directories
##################################
library(VAST)
library(readxl)

##################################
## Set up directories
##################################
which_machine = c('Zack_PC' =1, 'Zack_GI_PC'=2)[1]

github_dir = paste0(c('C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/Optimal_Allocation_GoA/')

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load(paste0(github_dir, 'data/optimization_data.RData'))

GOA_allocations = read_xlsx(paste0(github_dir, 'data/',
                                   'GOA 2019 stations by stratum.xlsx'))

GOA_grid = make_extrapolation_info(Region = 'Gulf_of_Alaska')$Data_Extrap
rm(gulf_of_alaska_grid)

##################################
## Constants
##################################
Niter = 1000
nstrata = nrow(GOA_allocations)

##################################
## Simulate Survey
##################################
sim_mean = sim_cv = array(dim = c(NTime, ns, Niter, 3), 
                          dimnames = list(paste0('Year_', 1:NTime),
                                          sci_names,
                                          NULL,
                                          paste0('n=', c(280,550,820)) ))

set.seed(23234)

for(isample in 1:3){
  
  #Adjust sample size proportionally
  nh = floor(GOA_allocations$`Number stations`*c(0.5,1,1.5)[isample])
  #Makes sure min number of samples within a strata is 2
  nh = sapply(nh, FUN = function(x) max(x, 2))
  
  #strata constraints
  stratano = rep(GOA_allocations$Stratum, nh)
  Wh = (table(GOA_grid$GOA_STRATUM)/N)[paste(GOA_allocations$Stratum)]
  wh = vector(length = length(Wh)); names(wh) = names(Wh)
  
  for(istrata in 1:nstrata){
    temp_allocation = nh[istrata]
    temp_pop =  table(GOA_grid$GOA_STRATUM)[paste(names(Wh))[istrata]]
    wh[istrata] = temp_allocation / temp_pop
  }
  
  for(iyear in 1:NTime){
    
    #Subset densities
    sub_df = subset(frame_raw, year == iyear)
    
    for(iter in 1:Niter){
      
      #Sample grid cells by stratum allocations
      sample_idx = c()
      for(istrata in 1:nstrata){
        temp_nh = nh[istrata]
        str_idx = which(GOA_grid$GOA_STRATUM == GOA_allocations$Stratum[istrata])
        sample_idx = c(sample_idx, sample(x = str_idx, size = temp_nh))
      }
      
      sample_df = frame_raw[sample_idx,]
      
      #Calculate Stratum Mean Density and Variance
      stmt = paste0('aggregate(cbind(',
                    paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                    ") ~ stratano, data = sample_df, FUN = mean)")
      sample_mean = eval(parse(text = stmt))[,-1]
      stmt = paste0('aggregate(cbind(',
                    paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                    ") ~ stratano, data = sample_df, FUN = var)")
      sample_var = eval(parse(text = stmt))[,-1]
      
      #Calculate Total Mean, Variance, CV
      SRS_var = colSums(
        sweep(x = sample_var, MARGIN = 1, 
              STATS = (Wh)^2*(1-wh)/GOA_allocations$`Number stations`,
              FUN = '*'))
      
      SRS_mean = colSums(sweep(x = sample_mean, MARGIN = 1, 
                               STATS = Wh,
                               FUN = '*'))
      
      strata_cv = sqrt(SRS_var) / SRS_mean 
      
      #Record mean and CV values
      sim_mean[paste0('Year_', iyear), , iter, isample] = SRS_mean
      sim_cv[paste0('Year_', iyear), , iter, isample] = strata_cv
    }
  }
  print(paste('Finished with n =', c(280,550, 820)[isample]))
}


#################################
## Simulation Metrics
#################################
#True CV, Cv of Cv, Rrmse of Cv
true_cv_array = rrmse_cv_array = rrmse_est_array = 
  array(dim = c(NTime, ns, 3), 
        dimnames = list(paste0('Year_', 1:NTime), sci_names, NULL ))

for(isample in 1:3){
  for(iyear in 1:NTime){
    for(spp in sci_names){
      
      iter_est = sim_mean[paste0('Year_', iyear), spp, , isample]
      iter_cv = sim_cv[paste0('Year_', iyear), spp, , isample ]
      true_cv = sd(iter_est) / true_mean[iyear, spp]
      temp_true_cv = true_mean[iyear,spp]
      true_cv_array[paste0('Year_', iyear), spp, isample ] = true_cv
      
      rrmse_cv_array[paste0('Year_', iyear), spp, isample] =
        sqrt(mean((iter_cv-true_cv)^2)) / mean(iter_cv)
      
      rrmse_est_array[iyear,spp,isample] = 
        sqrt(mean(((iter_est - temp_true_cv )^2)))/temp_true_cv
    }
  }
}


##################################
## Save
##################################
for(ivar in  c('rrmse_cv_array', 'true_cv_array', 'rrmse_est_array',
               'sim_mean', 'sim_cv')){
  assign(x=paste0('Survey_', ivar), value = get(ivar))
}

save(file = paste0(github_dir, 'Simulate_Current_Survey/',
                   'Survey_Simulation_Results.RData'),
     list = c(paste0('Survey_', c('rrmse_cv_array', 'rrmse_est_array',
                                  'true_cv_array', 'sim_mean', 'sim_cv'))))

