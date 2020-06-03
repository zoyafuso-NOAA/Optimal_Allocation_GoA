######################################
## Plot Performance Metrics
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
#load(paste0(output_wd, '/OM.RData'))
#load(paste0(output_wd, '/optimization.RData'))
load(paste0(output_wd, '/survey_simulation.RData'))

par(mfcol = c(2,2), mar = c(0,0,0,0), oma = c(4,5,2,1))
for(iOM in 1:2){
 boxplot( true_cv_array[,paste0('OM', iOM),1:2],
          xlab = 'EM type', ylab = 'True CV', 
          ylim = c(0,0.12), axes = F)
 abline(h = 0.05, lty = 'dashed')
 box()
 
 colmed = apply(true_cv_array[,paste0('OM', iOM),1:2], MARGIN = 2, median )
 print(100 * diff(colmed) / colmed[2])
 
 if(iOM == 1) {
  mtext(side = 2, 'True CV', line = 3.5)
  axis(side = 2, las = 1)
 }
 mtext(side = 3, c('Population A',
                   'Population B')[iOM], line = 0.5)
 
 boxplot( rrmse_cv_array[,paste0('OM', iOM),1:2],
          xlab = 'EM type', ylab = 'RRMSE of CV',
          ylim = c(0,0.28), axes = F)
 box()
 if(iOM == 1) {
  mtext(side = 2, 'RRMSE of CV', line = 3.5)
  axis(side = 2, las = 1)
 }
 colmed = apply(rrmse_cv_array[,paste0('OM', iOM),1:2], MARGIN = 2, median )
 print(100 * diff(colmed) / colmed[2])
 axis(side = 1, at = 1:2, c('Population A', 'Population B'))
}
mtext(side = 1, 'Survey Optimization Type', outer = T, line = 2.5)
