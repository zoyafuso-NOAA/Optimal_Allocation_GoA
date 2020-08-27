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
GOA_allocations3 = read_xlsx(paste0(github_dir, 'data/',
                                   'GOA2019_ 3 boat_825_RNDM_stations.xlsx')) 

allocations = data.frame(Stratum = sort(unique(GOA_allocations3$stratum)),
                         boat3 = aggregate(id ~ stratum, 
                                           data = GOA_allocations3, 
                                           FUN = length)$id,
                         boat2 = c(GOA_allocations$`Number stations`,
                                   rep(0,5)))
allocations$boat1 = ceiling(allocations$boat2 / 2)

allocations$boat1 = ifelse(allocations$boat1 == 0, 0, 
                           ifelse(allocations$boat1 == 1, 2, allocations$boat1))

GOA_grid = make_extrapolation_info(Region = 'Gulf_of_Alaska')$Data_Extrap
rm(gulf_of_alaska_grid)

##################################
## Constants
##################################
Niter = 100

##################################
## Simulate Survey
##################################
sim_mean = sim_cv = array(dim = c(NTime, ns, Niter, 3), 
                          dimnames = list(paste0('Year_', 1:NTime),
                                          sci_names,
                                          NULL,
                                          paste0('n=', c(280,550,820)) ))

set.seed(23234)

options(warn = 2)

for(isample in 1:3){
  
  #Adjust sample size proportionally
  nh = allocations[,paste0('boat', isample)]
  sampled_strata = nh > 0
  nstrata = sum(sampled_strata)
  
  nh = nh[sampled_strata]
  
  #strata constraints
  stratano = rep(allocations$Stratum[sampled_strata], nh)
  
  Nh = table(GOA_grid$GOA_STRATUM)[sampled_strata]
  Wh = Nh / N
  Wh = Wh
  # wh = vector(length = length(Wh)); names(wh) = names(Wh)
  wh = nh / Nh
  
  for(iyear in 1:NTime){
    
    #Subset densities
    sub_df = subset(frame_raw, year == iyear)
    
    for(iter in 1:Niter){
      
      #Sample grid cells by stratum allocations
      sample_idx = c()
      for(istrata in 1:nstrata){
        temp_nh = nh[istrata]
        str_idx = which(GOA_grid$GOA_STRATUM == allocations$Stratum[istrata])
        sample_idx = c(sample_idx, sample(x = str_idx, size = temp_nh))
      }
      
      sample_df = sub_df[sample_idx,]
      
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
              STATS = (Wh)^2*(1-wh)/nh,
              FUN = '*'))
      
      SRS_mean = colSums(sweep(x = sample_mean, MARGIN = 1, 
                               STATS = Wh,
                               FUN = '*'))
      
      strata_cv = sqrt(SRS_var) / SRS_mean 
      
      #Record mean and CV values
      sim_mean[paste0('Year_', iyear), , iter, isample] = SRS_mean
      sim_cv[paste0('Year_', iyear), , iter, isample] = strata_cv
    }
    
    print(paste('Finished with n =', c(280,550,820)[isample], 'Year', iyear))
  }
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
      
      true_cv_array[paste0('Year_', iyear), spp, isample ] = true_cv
      
      rrmse_cv_array[paste0('Year_', iyear), spp, isample] =
        sqrt(mean((iter_cv-true_cv)^2)) / mean(iter_cv)
      
      temp_true_cv = true_mean[iyear,spp]
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


par(mar = c(1,3,1,1), mfrow = c(5,3))
for(ispp in 1:ns){
  boxplot(true_cv_array[,ispp,], las = 1)
}

