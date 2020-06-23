#############################
## Show Simulation Metrics for Simulated Survey Strata,
## across Strata
## Supplemental Figures
#############################
rm(list = ls())

library(SamplingStrata); library(RColorBrewer)

############################
## Set up directories
#############################
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2)[2]
modelno = '6g'
github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/')[which_machine],
                    'GitHub/MS_OM_GoA/Optimum_Allocation/', 
                    'model_', modelno, '/')

output_dir = paste0(c('/Users/zackoyafuso/', 
                      'C:/Users/Zack Oyafuso/')[which_machine],
                    'Google Drive/MS_Optimizations/figure_plot/')

load(paste0(github_dir, 'optimization_results.RData'))
load(paste0(github_dir, "optimization_data_model_", modelno ,".RData"  ))

settings$id = 1:nrow(settings)

contingency = list()

samples = c(280, 550, 820)
strata = c(5,10,15)

for(istrata in strata){
  for(isample in samples){
    ##########################
    ## Create Input for bethel algorithm
    ##########################
    
    #Subset runs for a given strata number and sample level
    sub_settings = subset(settings, nstrata == istrata)
    idx = sub_settings$id[which.min(abs(sub_settings$n - isample))]
    
    #Extract stratification information:
    #number of strata, stratas, total units within strata
    stratano = rep(res_df[,idx+1], times= 11)
    nstrata = length(unique(stratano))
    stratapop = strata_list[[idx]]$Population
  
    #Stratum Means and Variances
    stmt = paste0('aggregate(cbind(', paste0('Y', 1:(ns-1), ',', collapse = ''),
                  'Y', ns, ') ~ stratano, FUN=mean, data=frame_raw, drop = F)')
    sample_mean = eval(parse(text = stmt))[,-1]
    names(sample_mean) = paste0('M', 1:ns)
    
    stmt = paste0('aggregate(cbind(', paste0('Y', 1:(ns-1), ',', collapse = ''),
                  'Y', ns, ') ~ stratano, FUN=var, data = frame_raw, drop = F)')
    sample_var = eval(parse(text = stmt))[,-1]
    names(sample_var) = paste0('S', 1:ns)
    
    sample_mean[is.na(sample_var)] = 0
    sample_var[is.na(sample_var)] = 0
    
    
    df = cbind(data.frame(stratum = 1:nstrata,
                          N = stratapop),
               sample_mean, sqrt(sample_var),
               data.frame(cens = 0,
                          cost = 1,
                          DOM1 = 'tot'))
    
    #############################
    ## Result Object
    #############################
    allocation_tab = matrix(ncol = istrata, nrow = isample-(istrata*5))
    cvs = c()
   
    start_CV = sub_settings$cv[which.min(abs(sub_settings$n - isample))]
    total_n = sub_settings$n[which.min(abs(sub_settings$n - isample))]
    current_n = total_n
    
    change_in_n = total_n - current_n
    icv = start_CV
    
    for(i in 1:(isample-(istrata*5))){
      while(change_in_n < i){
        stmt = paste0('cbind(',paste0("CV",1:ns,'=',icv,collapse=','),')')
        CV = eval(parse(text = stmt))
        errors = cbind(data.frame(DOM = 'DOM1'), CV, domainvalue = 1)
        
        n = bethel(stratif = df, errors = errors, printa=TRUE, 
                   epsilon = 1e-11, maxiter = 200)
        total_n = sum(n)
        change_in_n = current_n - total_n
        
        cv_inc = ifelse(i > 10, -mean(diff(cvs[(i-1):(i-11)]) ), 0.0001 )
        
        icv = icv + cv_inc
      }
      
      cvs = c(cvs, icv - cv_inc)
      allocation_tab[i,1:length(n)] = n
      if((i%%50) == 0) print(paste(istrata, 'Strata,', isample, 'Samples,', 
                                   'Run', i, 'CVinc =', round(cv_inc,4)))
    }
    
    contingency[[paste0("Strata_", istrata)]][[paste0('Samples_', isample)]] =
      list(cvs = cvs, allocation_tab = allocation_tab)
    save(list = 'contingency', file = paste0(github_dir, 'contingency.RData'))
  }
}

dev.off()
istrata = 5; isample = 820

par(mar = c(5,5,1,5))
allocation_tab = contingency[[paste0('Strata_', istrata)]][[paste0('Samples_', isample)]]$allocation_tab
cvs = contingency[[paste0('Strata_', istrata)]][[paste0('Samples_', isample)]]$cvs

matplot(rowSums(allocation_tab), allocation_tab, pch = 16, cex = 1.5,
        col= brewer.pal(n = istrata, name = 'Paired'),
        xlab = 'Total Sample Size', ylab = 'Sample Size', las = 1)
par(new = T)
plot(na.omit(rowSums(allocation_tab)), cvs, type = 'l', axes = F, ann = F, )
axis(side = 4, las = 1)
mtext(side = 4, 'Spatiotemporal CV Upper Bound', line = 3)
abline(v = seq(50,250,50), h = seq(0.4,1.4,0.2),
       lty = 'dashed', col = 'lightgrey')
legend('top', legend = paste('Stratum ', 1:istrata), 
       fill = brewer.pal(n = istrata, name = 'Paired'), ncol = 2)

