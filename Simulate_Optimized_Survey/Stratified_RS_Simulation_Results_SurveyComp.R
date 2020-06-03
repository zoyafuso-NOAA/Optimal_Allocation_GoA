######################################
## Simulate Random Sampling according to the optimized
## stratifications
######################################
rm(list = ls())

###############################
## Set up directories
###############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[2]
VAST_model = "6g"
optimization_type = c('_spatial', '_spatiotemporal')[2]

VAST_dir = paste0(c('/Users/zackoyafuso/Google Drive/', 
                    'C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/',
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output', VAST_model)


output_wd = paste0(c('/Users/zackoyafuso/Documents/', 
                     'C:/Users/Zack Oyafuso/Documents/',
                     'C:/Users/zack.oyafuso/Work/', 
                     'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                   "GitHub/MS_OM_GoA/Optimum_Allocation/model_", VAST_model,
                   optimization_type)

#########################
## Load data
#########################
load(paste0(output_wd, '/optimization_data_model_', VAST_model, '.RData'))


#Constants
ids = 1:nrow(frame)
N = length(ids)
strata = c(5,10,15,20,25,30,40,50,60)
Nstrata = length(strata)
Niters = 1000
sci_names = c("Atheresthes stomias", "Gadus chalcogrammus", 
              "Gadus macrocephalus", "Glyptocephalus zachirus" , 
              "Hippoglossoides elassodon", "Hippoglossus stenolepis", 
              "Lepidopsetta bilineata", "Lepidopsetta polyxystra", 
              "Limanda aspera", "Microstomus pacificus",
              "Sebastes alutus", "Sebastes B_R", "Sebastes polyspinis", 
              "Sebastes variabilis", "Sebastolobus alascanus" )

#Add year column to the raw dataframe, modify res_df 
frame_raw$year = rep(1:NTime, each = N)

#################
## Function to reproduce unique seed
#################
# getseed = function(temp_time, temp_nstrata, temp_CV, temp_Niter){
#  
#  NTime = 11 #11 observed Years
#  NStrata = 60 #60 possible strata, only a subset are considered
#  NCV = length(seq(0, 0.5, 0.005)) #101 possible CV values
#  NIter = 1000 #1000 iterations, only 100 explored
#  
#  temp_CV2 = which(as.integer(temp_CV*1000)==as.integer(seq(0,0.5,0.005)*1000))
#  
#  return((NTime*NStrata*NCV) * (temp_Niter-1) +
#          (NTime * NStrata) * (temp_CV2-1) +
#          (NTime) * (temp_nstrata - 1) +
#          temp_time )
# }

#################
# true density
#################
stmt = paste0('aggregate(cbind(',
              paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
              ") ~ year, data = frame_raw, FUN = mean)")
true_mean = eval(parse(text = stmt))[,-1]
colnames(true_mean) = sci_names

sim_mean = sim_cv = array(dim = c(NTime, ns, Nstrata, 3, Niters), 
                          dimnames = list(paste0('Year_', 1:NTime),
                                          sci_names, 
                                          NULL,
                                          NULL, 
                                          NULL))

##########################
## Simulating each optimization
##########################
#for(irow in 1:nrow(settings)) {

for(isample in 1:3){
  for(istrata in 1:Nstrata){
    
    temp_dir = paste0(output_wd, '/Survey_Comparison/Str_', strata[istrata],
                      'Boat_', isample)
    
    if(dir.exists(temp_dir)){
      print(temp_dir)
      load(paste0(temp_dir, '/result_list.RData') )
      
      strata_allocation = result_list[[2]]$Allocation
      stratapop = result_list[[2]]$Population
      
      #Good strata (with > 1 cells)
      temp_strata = (1:length(strata_allocation))[strata_allocation >= 2]
      
      for(iyear in 1:NTime){
        for(iter in 1:Niters){
          
          #Set Unique Seed
          # temp_seed =  getseed(temp_time = iyear, 
          #                      temp_nstrata = settings$nstrata[irow], 
          #                      temp_CV = settings$cv[irow], 
          #                      temp_Niter = iter)
          # set.seed(temp_seed)
          
          #Sample based on the stratification allocations
          sample_vec = c()
          sol = result_list[[1]]$indices$X1
          
          for(i in temp_strata ){
            available_cells = which(sol == i)
            sample_cells = sample(x = available_cells, 
                                  size = strata_allocation[i], 
                                  replace = F)
            sample_vec = c(sample_vec, sample_cells)
          }
          
          #Organize sample set and total number of samples
          sample_vec = sort(sample_vec)
          n = length(sample_vec)
          stratano =  sol[sample_vec]
          sample_df = subset(frame_raw, year == iyear)[sample_vec,]
          
          #Calculate Stratum Mean Density and Variance
          stmt = paste0('aggregate(cbind(',
                        paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                        ") ~ stratano, data = sample_df, FUN = mean)")
          sample_mean = eval(parse(text = stmt))[,-1]
          stmt = paste0('aggregate(cbind(',
                        paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                        ") ~ stratano, data = sample_df, FUN = var)")
          sample_var = eval(parse(text = stmt))[,-1]
          
          #How many samples are allocated in each strata
          #How many sampling units are in each strata
          temp_strata_allocation = strata_allocation[temp_strata] #n_h
          temp_stratapop = stratapop[temp_strata] #N_h
          
          Wh = temp_stratapop/N
          wh = temp_strata_allocation/temp_stratapop
          
          #Calculate Total Abundance and Variance, calculate CV
          SRS_var = colSums(sweep(x = sample_var, MARGIN = 1, 
                                  STATS = (Wh)^2*(1-wh)/temp_strata_allocation,
                                  FUN = '*'))
          
          SRS_mean = colSums(sweep(x = sample_mean, MARGIN = 1, 
                                   STATS = temp_stratapop / N,
                                   FUN = '*'))
          
          strata_cv = sqrt(SRS_var) / SRS_mean 
          
          #Record mean and CV values
          sim_mean[paste0('Year_', iyear), , istrata,isample,iter] = SRS_mean
          sim_cv[paste0('Year_', iyear), , istrata,isample,iter] = strata_cv
          
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
  for(isample in 1:3){
    for(istrata in 1:Nstrata){
      for(spp in sci_names){
        
        iter_est = sim_mean[paste0('Year_', iyear), spp, istrata,isample, ]
        iter_cv = sim_cv[paste0('Year_', iyear), spp, istrata,isample, ]
        true_cv = sd(iter_est) / true_mean[iyear, spp]
        
        true_cv_array[paste0('Year_', iyear), spp, istrata,isample] = true_cv
        
        rrmse_cv_array[paste0('Year_', iyear), spp, istrata,isample] = 
          sqrt(mean((iter_cv-true_cv)^2)) / mean(iter_cv)
        
        rrmse_est_array[paste0('Year_', iyear), spp, istrata,isample]  = 
          sqrt(mean((iter_est - true_mean[iyear,spp])^2)) / true_mean[iyear,spp]
      }
    }
  }
}

#######################
## Save results
#######################
for(ivar in c('rrmse_cv_array', 'rrmse_est_array', 
              'true_cv_array', 'sim_mean', 'sim_cv')){
  assign(x=paste0('STRS_', ivar), value = get(ivar))
}

save(file = paste0(output_wd, '/Stratified_RS_Simulation_Results_SurveyComp.RData'),
     list = c(paste0('STRS_', c('rrmse_cv_array', 'rrmse_est_array', 
                                'true_cv_array', 'sim_mean', 'sim_cv')),
              'true_mean', 'sci_names', 'NTime', 'ns', 'Niters', 'N'))
