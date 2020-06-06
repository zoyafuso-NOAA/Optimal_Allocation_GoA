#############################
## Show Simulation Metrics for Simulated Survey Strata,
## Simple Random Sampling, and Stratified Random Sampling
#############################
rm(list = ls())

############################
## Set up directories
#############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[1]

github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/',
                      'C:/Users/zack.oyafuso/Work/', 
                      'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                    "GitHub/Optimal_Allocation_GoA/")

PP_dir = paste0(c('/Users/zackoyafuso/', 
                  'C:/Users/Zack Oyafuso/')[which_machine],
                'Google Drive/MS_Optimizations/powerpoint_plot/')

##########################
## Settings for Spatial or Spatiotemporal Plot
##########################
samples = c(280, 550, 820)
plot_settings = data.frame(type = c('Spatial', 'Spatiotemporal'),
                           data_filename = c('spatial_only_', 
                                             'spatiotemporal_'),
                           ymax = c(1300, 1100),
                           xlabel = c(0.18, 0.28))

{ 
  png(paste0(PP_dir, 'N_CV_Tradeoff.png'),
      width = 5, height = 7, res = 500, units = 'in')
  par(mfrow = c(2,1), mar = c(2,3,1,1), oma = c(2,2,0,0))
  for(irow in 1:2){
    ##########################
    ## Load Data
    ##########################
    load(paste0(github_dir, plot_settings$type[irow], '_Optimization/',
                plot_settings$data_filename[irow], 'optimization_results.RData'))
    
    #############################
    ## Plot
    #############################
    plot(n ~ cv, data = settings, subset = nstrata == 5, las = 1, pch = 16,
         ylim = c(0, plot_settings$ymax[irow]),
         ann = F)
    #xlab = 'Upper Spatiotemporal CV Constraint', ylab = 'Total Sample Size')
    lines(n ~ cv, data = settings, subset = nstrata == 5)
    
    abline(h = samples, col = 'grey', lty = 'dashed')
    text(x = plot_settings$xlabel[irow], y = samples, 
         paste(1:3, 'Boats'), pos = 1)
    legend('top', legend = paste(plot_settings$type[irow], 'Optimization'), 
           bty = 'n')
    
  }
  
  
  mtext(side = 1, 'Upper CV Constraint', outer = T, line = 0.5)
  mtext(side = 2, 'Total Optimized Sample Size', outer = T, line = 0.5)
  
  dev.off()
}
