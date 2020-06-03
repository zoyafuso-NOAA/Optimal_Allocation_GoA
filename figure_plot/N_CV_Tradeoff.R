#############################
## Show Simulation Metrics for Simulated Survey Strata,
## Simple Random Sampling, and Stratified Random Sampling
#############################
rm(list = ls())

############################
## Set up directories
#############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[2]
optimization_type = c('_spatial', '_spatiotemporal')[2]
VAST_model = "6g"

output_wd = paste0(c('/Users/zackoyafuso/Documents/', 
                     'C:/Users/Zack Oyafuso/Documents/',
                     'C:/Users/zack.oyafuso/Work/', 
                     'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                   "GitHub/MS_OM_GoA/Optimum_Allocation/model_", VAST_model,
                   optimization_type)
PP_dir = paste0(c('/Users/zackoyafuso/', 
                  'C:/Users/Zack Oyafuso/')[which_machine],
                'Google Drive/MS_Optimizations/powerpoint_plot/')

load(paste0(output_wd, '/Stratified_RS_Simulation_Results.RData'))
load(paste0(output_wd, '/optimization_results.RData'))

samples = c(280, 550, 820)

#############################
## Plot
#############################
{
  png(paste0(PP_dir, 'N_CV_Example_', optimization_type, '.png'), 
      width = 5, height = 5, res = 500, units = 'in')
 par(mfrow = c(1,1), mar = c(4,4,1,1))
 plot(n ~ cv, data = settings, subset = nstrata == 5, las = 1, pch = 16,
      ylim = c(0,1100), col = 'darkgrey',
      xlab = 'Upper Spatiotemporal CV Constraint', ylab = 'Total Sample Size')
 lines(n ~ cv, data = settings, subset = nstrata == 5, col = 'grey')
 
 abline(h = samples, col = 'grey', lty = 'dashed')
 text(x = 0.28, y = samples, c('1 Boat', '2 Boats', '3 Boats'), 
      pos = 1, col = 'darkgrey')
 
 sub_settings = subset(settings, nstrata == 5)
 for(isample in samples){
  points(sub_settings[which.min(abs(sub_settings$n - isample)), c('cv', 'n')], 
         pch = 16, cex = 2, col = 'black')
 }
 dev.off()
 }
