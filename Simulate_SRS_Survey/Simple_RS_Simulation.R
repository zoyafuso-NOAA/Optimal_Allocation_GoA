#############################
## Simple Random Sampling 550 and 800 samples
#############################
rm(list = ls())

###############################
## Set up directories
###############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[1]

github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/',
                      'C:/Users/zack.oyafuso/Work/', 
                      'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                    "GitHub/Optimal_Allocation_GoA/")

#########################
## Load predicted density and optimization results
#########################
load(paste0(github_dir, 'data/optimization_data.RData'))

#Constants
ids = as.numeric(rownames(frame))
Niters = 100

###########################
## Result Objects
###########################
sim_mean = sim_cv = array(dim = c(NTime, ns, 3, Niters), 
                          dimnames = list(paste0('Year_', 1:NTime),
                                          sci_names, 
                                          NULL))
#Do Runs
set.seed(233)
for(iyear in 1:NTime){
   for(isample in 1:3){
      for(iter in 1:Niters){
         samplesize = c(280,550,820)[isample]
         sample_vec = sample(1:N, size = samplesize)
         sample_df = subset(frame_raw, year == iyear)[sample_vec,]
         
         test_sim_mean = colMeans(sample_df[,paste0('Y', 1:ns)])
         test_sim_se= apply(sample_df[,paste0('Y', 1:ns)], MARGIN = 2, sd)/sqrt(samplesize)
         
         sim_mean[iyear,,isample,iter] = test_sim_mean
         sim_cv[iyear,,isample,iter] = test_sim_se / test_sim_mean
      }
   }
   print(paste0('Done with year', iyear))
}

############################
## Simulation Metrics
############################
#True CV
true_cv_array = cv_cv_array = rrmse_cv_array = 
   array(dim = c(NTime, ns, 3), 
         dimnames = list(NULL,
                         sci_names, 
                         NULL))

for(iyear in 1:NTime){
   for(isample in 1:3){
      for(spp in 1:ns){
         
         true_cv_array[iyear, spp, isample] = 
            sd(sim_mean[iyear,spp,isample,]) / true_mean[iyear,spp]
         
         rrmse_cv_array[iyear, spp, isample] = 
            (sum((sim_cv[iyear,spp,isample,] -  true_cv_array[iyear, spp, isample])^2) / 
                Niters)^0.5 / mean(sim_cv[iyear, spp, isample,])
         
         
      }
      
   }
}

#######################
## Save results
#######################
for(ivar in  c('rrmse_cv_array', 'true_cv_array', 
               'sim_mean', 'sim_cv')){
   assign(x=paste0('SRS_', ivar), value = get(ivar))
}

save(file = paste0(github_dir, 'Simulate_SRS_Survey/Simple_RS_Simulation_Results.RData'),
     list = c(paste0('SRS_', c('rrmse_cv_array', 
                               'true_cv_array', 'sim_mean', 'sim_cv')),
              "Niters"))

