####################################
## Simulation of Survey
####################################

rm(list = ls())

############################¬
## Import Libraries
############################
library(rgdal); library(raster); library(rgeos); library(tidyr)

############################
## Set up directories
#############################
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2)[1]
modelno = '6g'
github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/')[which_machine], 'GitHub/MS_OM_GoA/')
VAST_dir = paste0(c('/Users/zackoyafuso/Google Drive/VAST_Runs/', 
                    'C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/'),
                  'VAST_output', modelno)[which_machine]

############################
## Load Data
############################
load(paste0(VAST_dir, '/Spatial_Settings.RData'))
load(paste0(github_dir, '/Optimum_Allocation/model_', modelno, 
            '/optimization_data_model_', modelno, '.RData'))
sci_names = levels(Data_Geostat$spp)
frame_raw$year = rep(sort(unique(Data_Geostat$Year)), each = nrow(frame))
frame_raw$stratum = rep(gulf_of_alaska_grid$GOA_STRATUM, times = 11)

survey_data = subset(read.csv(paste0(github_dir, '/data/data/',
                                     'GOA_multspp_with_strata.csv')),
                     SPECIES_NAME == 'Sebastes alutus', 
                     select = c(YEAR, STRATUM))

#How many samples/extrapolation cells in each stratum
samples_by_str = table(survey_data)
cells_by_str = table(gulf_of_alaska_grid$GOA_STRATUM)

############################
## Simulate survey with same stratifications and allocations as currently done
############################
Niters = 100
sim_mean = sim_cv = array(dim = c(NTime, ns, Niters),
                          dimnames = list(rownames(samples_by_str), sci_names, NULL))

set.seed(234234)

for(iyear in rownames(samples_by_str)){
 
 temp_frame = subset(frame_raw, year == iyear)
 
 for(iter in 1:Niters){
  temp_samples = data.frame()
  
  for(istrata in colnames(samples_by_str)){
   temp_stratum = subset(temp_frame, stratum == istrata)
   which_samples = sample(x = 1:nrow(temp_stratum),
                          size = samples_by_str[iyear,istrata])
   temp_samples = rbind(temp_samples, temp_stratum[which_samples,])
  }
  
  stmt = paste0('aggregate(cbind(',
                paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                ") ~ stratum, data = temp_samples, FUN = mean)")
  sample_mean = eval(parse(text = stmt))[,-1]
  stmt = paste0('aggregate(cbind(',
                paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                ") ~ stratum, data = temp_samples, FUN = var)")
  sample_var = eval(parse(text = stmt))[,-1]
  
  temp_strata_allocation = samples_by_str[iyear,]
  temp_strata_allocation = temp_strata_allocation[temp_strata_allocation > 0]
  temp_stratapop = cells_by_str[names(temp_strata_allocation)]
  N = nrow(frame)
  
  SRS_var = colSums(sweep(x = sample_var, MARGIN = 1, 
                          STATS = (temp_stratapop/N)^2*(1 - temp_strata_allocation/temp_stratapop)/temp_strata_allocation,
                          FUN = '*'), na.rm = T)
  
  SRS_mean = colSums(sweep(x = sample_mean, MARGIN = 1, 
                           STATS = temp_stratapop / N,
                           FUN = '*'))
  
  strata_cv = sqrt(SRS_var) / SRS_mean 
  
  sim_mean[iyear, , iter] = SRS_mean
  sim_cv[iyear, , iter] = strata_cv
 }
 print(paste('Done with Year', iyear))
}

#########################
## Simulaton Metrics
#########################
#True density
stmt = paste0('aggregate(cbind(',
              paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
              ") ~ year, data = frame_raw, FUN = mean)")
true_mean = eval(parse(text = stmt))[,-1]
colnames(true_mean) = sci_names
rownames(true_mean) = rownames(samples_by_str)

#Result Object
true_cv_array = cv_cv_array = rrmse_cv_array = 
 array(dim = c(NTime, ns), 
       dimnames = list(rownames(samples_by_str), sci_names ))

#Calculate Metrics
for(iyear in rownames(samples_by_str)){
  for(spp in sci_names){
   true_cv_array[iyear, spp] = sd(sim_mean[iyear, spp,]) / true_mean[iyear,spp]
   
   cv_cv_array[iyear, spp] =  sd(sim_cv[iyear, spp,])/mean(sim_cv[iyear, spp,])
   
   rrmse_cv_array[iyear, spp] = 
    (sum((sim_cv[iyear, spp, ] - true_cv_array[iyear, spp])^2) / Niters)^0.5 / 
    mean(sim_cv[iyear, spp, ])
 }
}

#######################
## Save results
#######################
for(ivar in  c('cv_cv_array', 'rrmse_cv_array', 'true_cv_array', 
               'sim_mean', 'sim_cv')){
 assign(x=paste0('survey_', ivar), value = get(ivar))
}
 
save(file = paste0(github_dir, 'Optimum_Allocation/model_', 
                   modelno, '/Survey_Simulation_Results.RData'),
     list = c(paste0('survey_', c('cv_cv_array', 'rrmse_cv_array', 
                                  'true_cv_array', 'sim_mean', 'sim_cv')),
              'true_mean', 'sci_names', 'NTime', 'ns', 'Niters', 'N'))

