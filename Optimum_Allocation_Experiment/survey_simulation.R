######################################
## Simulate Random Sampling according to the optimized
## stratifications
######################################
rm(list = ls())

######################################
## Set up directories
######################################
which_machine = c('Zack_Mac' = 1)[1]
output_wd = paste0(c('/Users/zackoyafuso/Documents/')[which_machine],
                   'GitHub/MS_OM_GoA/Optimum_Allocation_Experiment/')

######################################
## Import OM
######################################
load(paste0(output_wd, '/OM.RData'))
load(paste0(output_wd, '/optimization.RData'))
res_df = res_df[,-1]

#####################################
## Constants
#####################################
Niters = 100
NTime = 10
NEM = NOM = 3
N = nrow(res_df)

#####################################
## Result Objects
#####################################

sim_mean = sim_cv = array(dim = c(NTime, NOM, NEM, Niters),
                          dimnames = list(paste0('Year', 1:NTime),
                                          paste0('OM', 1:NOM),
                                          paste0('EM', 1:NEM),
                                          NULL))

for(iEM in 1:nrow(ST_settings)) {
   
   #Set Estimation Model
   strata_allocation = strata_list[[iEM]]$Allocation
   stratapop = strata_list[[iEM]]$Population
   
   for(iOM in 1:nrow(ST_settings)){
      ##DFs for the optimizer
      frame = data.frame(id = 1:nrow(OM[[iOM]]$loc_xy),
                         X1 = OM[[iOM]]$loc_xy$x,
                         X2 = OM[[iOM]]$loc_xy$y,
                         Y1 = rowMeans(OM[[iOM]]$log_d_rt),
                         domainvalue = 1)
      
      frame_raw = data.frame(id = rep(1:nrow(OM[[iOM]]$loc_xy), times = 10),
                             X1 = rep(OM[[iOM]]$loc_xy$x, times = 10),
                             X2 = rep(OM[[iOM]]$loc_xy$y, times = 10),
                             Y1 = as.vector(OM[[iOM]]$log_d_rt),
                             year = rep(1:NTime, each = nrow(frame)),
                             domainvalue = 1)
      
      for(iyear in 1:NTime){
         for(iter in 1:Niters){
            
            #Set Unique Seed
            # temp_seed =  getseed(temp_time = iyear, 
            #                      temp_nstrata = settings$nstrata[irow], 
            #                      temp_CV = settings$cv[irow], 
            #                      temp_Niter = iter)
            # set.seed(temp_seed)
            
            #Sample based on the stratification allocations
            sample_vec = c()
            for(i in 1:length(strata_allocation) ){
               available_cells = which(res_df[,iEM] == i)
               sample_cells = sample(x = available_cells, 
                                     size = strata_allocation[i], 
                                     replace = F)
               sample_vec = c(sample_vec, sample_cells)
            }
            
            #Organize sample set and total number of samples
            sample_vec = sort(sample_vec)
            n = length(sample_vec)
            stratano =  res_df[sample_vec,iEM]
            sample_df = subset(frame_raw, year == iyear)[sample_vec,]
            
            #Calculate Stratum Mean Density and Variance
            sample_mean = aggregate(Y1 ~ stratano, 
                                    data = sample_df, FUN = mean)[,-1]
            sample_var = aggregate(Y1 ~ stratano, 
                                   data = sample_df, FUN = var)[,-1]
            
            #How many samples are allocated in each strata
            #How many sampling units are in each strata
            #temp_strata_allocation = strata_allocation[temp_strata] #n_h
            #temp_stratapop = stratapop[temp_strata] #N_h
            temp_strata_allocation = strata_allocation #n_h
            temp_stratapop = stratapop #N_h
            
            #Calculate Total Abundance and Variance, calculate CV
            SRS_var = sum(sample_var * (temp_stratapop/N)^2*(1 - temp_strata_allocation/temp_stratapop)/temp_strata_allocation, na.rm = T)
            
            SRS_mean = sum( sample_mean * temp_stratapop / N )
            
            strata_cv = sqrt(SRS_var) / SRS_mean 
            
            #Record mean and CV values
            sim_mean[paste0('Year', iyear), iOM, iEM, iter] = SRS_mean
            sim_cv[paste0('Year', iyear), iOM, iEM, iter] = strata_cv
            
         }
      }
      print(paste0('OM:', iOM, ', EM: ', iEM))
   }
}

boxplot( t(sim_cv[,paste0('OM', 1), paste0('EM', 3),]), ylim = c(0,0.1) )
abline(h = .05, col = 'red', lty = 'dashed')

#################################
## Simulation Metrics
#################################
true_mean = sapply(X = OM, FUN = function(x) colMeans(x$log_d_rt) )

#True CV, Cv of Cv, Rrmse of Cv
true_cv_array = rrmse_cv_array = array(dim = c(NTime, NOM, NEM), 
                                       dimnames = list(paste0('Year',1:NTime),
                                                       paste0('OM', 1:NOM),
                                                       paste0('EM', 1:NEM)) )

for(iyear in 1:NTime){
   for(iOM in 1:NOM ){
      for(iEM in 1:NEM){
         iter_ests = sim_mean[paste0('Year',iyear),iOM,iEM, ]
         iter_cvs = sim_cv[paste0('Year',iyear),iOM,iEM, ]
         true_cv = sd(iter_ests)/true_mean[iyear,iOM]
         RMSE = (sum((iter_cvs - true_cv)^2) / Niters)^0.5
         
         true_cv_array[paste0('Year', iyear), iOM, iEM] = true_cv
         rrmse_cv_array[paste0('Year',iyear),iOM,iEM] =  RMSE / mean(iter_cvs)
      }
   }
}

#######################
## Save 
#######################
save(list = c('true_mean', 'sim_cv', 'sim_mean', 
              'true_cv_array', 'rrmse_cv_array'), 
     file = paste0(output_wd, 'survey_simulation.RData'))
