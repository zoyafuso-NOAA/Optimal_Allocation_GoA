#############################
## Simple Random Sampling 550 and 800 samples
#############################

rm(list = ls())

###############################
## Set up directories
###############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[1]
VAST_model = "6g"
github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/Optimum_Allocation/model_', VAST_model)

VAST_dir = paste0(c('/Users/zackoyafuso/Google Drive/', 
                    'C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/',
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output', VAST_model)

output_wd = c(paste0('/Users/zackoyafuso/Documents/GitHub/MS_OM_GoA/',
                     'Optimum_Allocation/model_', VAST_model),
              paste0("C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/",
                     "Optimum_Allocation/model_", VAST_model),
              paste0("C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/",
                     "Optimum_Allocation/model_", VAST_model),
              paste0("C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/",
                     "Optimum_Allocation/model_", VAST_model))[which_machine]


#########################
## Load data
#########################
load(paste0(github_dir, '/optimization_data_model_', VAST_model, '.RData'))

N = nrow(frame)
nsamples = seq(350, 1000, by = 50)
Niters = 100
sci_names = c("Atheresthes stomias", "Gadus chalcogrammus", "Gadus macrocephalus", 
              "Glyptocephalus zachirus" , "Hippoglossoides elassodon", 
              "Hippoglossus stenolepis", "Lepidopsetta bilineata", 
              "Lepidopsetta polyxystra", "Limanda aspera", "Microstomus pacificus",
              "Sebastes alutus", "Sebastes B_R", "Sebastes polyspinis", 
              "Sebastes variabilis", "Sebastolobus alascanus" )

frame_raw$year = rep(1:NTime, each = N)

#################
# true density
#################
stmt = paste0('aggregate(cbind(',
              paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
              ") ~ year, data = frame_raw, FUN = mean)")
true_mean = eval(parse(text = stmt))[,-1]
colnames(true_mean) = sci_names

sim_mean = sim_cv = array(dim = c(NTime, ns, length(nsamples), Niters), 
                          dimnames = list(paste0('Year_', 1:NTime),
                                          sci_names, 
                                          paste0('size_', nsamples), 
                                          NULL))

#Do Runs
set.seed(233)
for(iyear in 1:NTime){
 for(isample in nsamples){
  for(iter in 1:Niters){
   sample_vec = sample(1:N, size = isample)
   sample_df = subset(frame_raw, year == iyear)[sample_vec,]
   
   test_sim_mean = colMeans(sample_df[,paste0('Y', 1:ns)])
   test_sim_se= apply(sample_df[,paste0('Y', 1:ns)], MARGIN = 2, sd)/sqrt(isample)
   
   sim_mean[paste0('Year_', iyear), , 
            paste0('size_', isample), iter] = test_sim_mean
   sim_cv[paste0('Year_', iyear), , 
            paste0('size_', isample), iter] = test_sim_se / test_sim_mean
  }
 }
 print(paste0('Done with year', iyear))
}

############################
## Simulation Metrics
############################
#True CV
true_cv_array = cv_cv_array = rrmse_cv_array = 
 array(dim = c(NTime, ns, length(nsamples)), 
       dimnames = list(paste0('Year_', 1:NTime),
                       sci_names, 
                       paste0('size_', nsamples)))

for(iyear in 1:NTime){
 for(isample in nsamples){
  for(spp in sci_names){
   
   true_cv_array[paste0('Year_', iyear), spp, 
                 paste0('size_', isample)] = 
    sd(sim_mean[paste0('Year_', iyear), spp, 
                paste0('size_', isample),]) / true_mean[iyear,spp]
   
   cv_cv_array[paste0('Year_', iyear), spp, 
               paste0('size_', isample)] = 
    sd(sim_cv[paste0('Year_', iyear), spp, 
              paste0('size_', isample), ] ) / 
    mean(sim_cv[paste0('Year_', iyear), spp, 
                paste0('size_', isample), ] )
   
   rrmse_cv_array[paste0('Year_', iyear), spp, 
                  paste0('size_', isample)] = 
    (sum((sim_cv[paste0('Year_', iyear), spp, 
                 paste0('size_', isample),] - 
           true_cv_array[paste0('Year_', iyear), spp, 
                         paste0('size_', isample)])^2) / 
      Niters)^0.5 / mean(sim_cv[paste0('Year_', iyear), spp, 
                                paste0('size_', isample) ,])
  }
    
 }
}

#######################
## Save results
#######################
for(ivar in  c('cv_cv_array', 'rrmse_cv_array', 'true_cv_array', 
               'sim_mean', 'sim_cv')){
   assign(x=paste0('SRS_', ivar), value = get(ivar))
}

save(file = paste0(github_dir, '/Simple_RS_Simulation_Results.RData'),
     list = c(paste0('SRS_', c('cv_cv_array', 'rrmse_cv_array', 
                                  'true_cv_array', 'sim_mean', 'sim_cv')),
              'true_mean', 'sci_names', 'NTime', 'ns', 'Niters', 'N', 
              'nsamples'))

