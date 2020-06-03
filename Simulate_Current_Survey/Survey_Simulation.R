##################################
## Simulate Survey Based on Current Stratifications
##################################
rm(list = ls())

library(VAST)
library(readxl)

##################################
## Set up directories
##################################
which_machine = c('Zack_PC' =1, 'Zack_GI_PC'=2)[1]
modelno = '6g'

github_dir = paste0(c('C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/')
VAST_dir = paste0(c('C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output', modelno, '/')

##################################
## Import Strata Allocations and Predicted Densities
##################################
load(paste0(github_dir, 'Optimum_Allocation/model_', modelno, 
            '_spatiotemporal/optimization_data_model_', modelno, '.RData'))

GOA_allocations = read_xlsx(paste0(github_dir, 'data/data/',
                                   'GOA 2019 stations by stratum.xlsx'))
GOA_grid = make_extrapolation_info(Region = 'Gulf_of_Alaska')$Data_Extrap
rm(gulf_of_alaska_grid)
load(paste0(VAST_dir, 'VAST_MS_GoA_Run.RData'))
D_gcy = Save$Report$Index_gcyl[,,,1]

##################################
## Constants
##################################
ns = 15
Niter = 1000
N = nrow(GOA_grid)
sci_names = c("Atheresthes stomias", "Gadus chalcogrammus", 
              "Gadus macrocephalus", "Glyptocephalus zachirus" , 
              "Hippoglossoides elassodon", "Hippoglossus stenolepis", 
              "Lepidopsetta bilineata", "Lepidopsetta polyxystra", 
              "Limanda aspera", "Microstomus pacificus",
              "Sebastes alutus", "Sebastes B_R", "Sebastes polyspinis", 
              "Sebastes variabilis", "Sebastolobus alascanus" )

Years2Include = 1 + as.vector(sort(unique(Save$TmbData$t_iz)))
NTime = length(Years2Include)
nstrata = nrow(GOA_allocations)

##################################
## Simulate Survey
##################################
sim_mean = sim_cv = array(dim = c(NTime, ns, Niter, 3), 
                          dimnames = list(paste0('Year_', 1:NTime),
                                          sci_names,
                                          NULL,
                                          paste0('n=', c(280,550,820)) ))

for(isamp in 1:3){
  
  #Adjust sample size proportionally
  nh = floor(GOA_allocations$`Number stations`*c(0.5,1,1.5)[isamp])
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
  
  for(iter in 1:Niter){
    for(iyear in 1:NTime){
      
      #Sample grid cells by stratum allocations
      sample_idx = c()
      for(istrata in 1:nstrata){
        temp_nh = nh[istrata]
        str_idx = which(GOA_grid$GOA_STRATUM == GOA_allocations$Stratum[istrata])
        sample_idx = c(sample_idx, sample(x = str_idx, size = temp_nh))
      }
      
      sample_df = as.data.frame(D_gcy[sample_idx,,Years2Include[iyear]])
      
      #Calculate Stratum Mean Density and Variance
      stmt = paste0('aggregate(cbind(',
                    paste0('V', 1:(ns-1), sep = ',', collapse = ''), 'V',ns, 
                    ") ~ stratano, data = sample_df, FUN = mean)")
      sample_mean = eval(parse(text = stmt))[,-1]
      stmt = paste0('aggregate(cbind(',
                    paste0('V', 1:(ns-1), sep = ',', collapse = ''), 'V',ns, 
                    ") ~ stratano, data = sample_df, FUN = var)")
      sample_var = eval(parse(text = stmt))[,-1]
      
      #Calculate Total Abundance and Variance, calculate CV
      SRS_var = colSums(
        sweep(x = sample_var, MARGIN = 1, 
              STATS = (Wh)^2*(1-wh)/GOA_allocations$`Number stations`,
              FUN = '*'))
      
      SRS_mean = colSums(sweep(x = sample_mean, MARGIN = 1, 
                               STATS = Wh,
                               FUN = '*'))
      
      strata_cv = sqrt(SRS_var) / SRS_mean 
      
      #Record mean and CV values
      sim_mean[paste0('Year_', iyear), , iter, isamp] = SRS_mean
      sim_cv[paste0('Year_', iyear), , iter, isamp] = strata_cv
    }
  }
  print(paste('Finished with n =', c(280,550, 820)[isamp]))
}


#################
# true density
#################
frame_raw$year = rep(1:NTime, each = N)
stmt = paste0('aggregate(cbind(',
              paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
              ") ~ year, data = frame_raw, FUN = mean)")
true_mean = eval(parse(text = stmt))[,-1]
colnames(true_mean) = sci_names

#################################
## Simulation Metrics
#################################
#True CV, Cv of Cv, Rrmse of Cv
true_cv_array = rrmse_cv_array = rrmse_est_array = 
  array(dim = c(NTime, ns, 3), 
        dimnames = list(paste0('Year_', 1:NTime), sci_names, NULL ))

for(isamp in 1:3){
  for(iyear in 1:NTime){
    for(spp in sci_names){
      
      iter_est = sim_mean[paste0('Year_', iyear), spp, , isamp]
      iter_cv = sim_cv[paste0('Year_', iyear), spp, , isamp ]
      true_cv = sd(iter_est) / true_mean[iyear, spp]
      temp_true_cv = true_mean[iyear,spp]
      true_cv_array[paste0('Year_', iyear), spp, isamp ] = true_cv

      rrmse_cv_array[paste0('Year_', iyear), spp, isamp] =
        sqrt(mean((iter_cv-true_cv)^2)) / mean(iter_cv)
      
      rrmse_est_array[iyear,spp,isamp] = 
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

save(file = paste0(github_dir, 'Optimum_Allocation/model_', modelno, 
                   '_spatiotemporal/Survey_Simulation_Results.RData'),
     list = c(paste0('Survey_', c('rrmse_cv_array', 'rrmse_est_array',
                                  'true_cv_array', 'sim_mean', 'sim_cv'))))

##################################
## 
##################################