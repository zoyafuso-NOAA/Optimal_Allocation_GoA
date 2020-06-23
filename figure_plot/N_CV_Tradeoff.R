#############################
## Show Simulation Metrics for Simulated Survey Strata,
## Simple Random Sampling, and Stratified Random Sampling
#############################
rm(list = ls())

############################
## Set up directories
#############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[2]

github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/',
                      'C:/Users/zack.oyafuso/Work/', 
                      'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                    "GitHub/Optimal_Allocation_GoA/")

PP_dir = paste0(c('/Users/zackoyafuso/', 
                  'C:/Users/Zack Oyafuso/')[which_machine],
                'Google Drive/MS_Optimizations/powerpoint_plot/')

figure_dir = paste0(c('/Users/zackoyafuso/', 
                  'C:/Users/Zack Oyafuso/')[which_machine],
                'Google Drive/MS_Optimizations/figure_plot/')

##########################
## Settings for Spatial or Spatiotemporal Plot
##########################
samples = c(280, 550, 820)
plot_settings = data.frame(type = c('Spatial', 'Spatiotemporal'),
                           data_filename = c('spatial_only_', 
                                             'spatiotemporal_'),
                           ymax = c(1300, 1100),
                           xmin = c(0.08, 0.15),
                           xmax = c(0.20, 0.30),
                           xlabel = c(0.18, 0.28))

{ 
  png(paste0(figure_dir, 'Fig3_N_CV_Tradeoff.png'),
      width = 90, height = 150, res = 500, units = 'mm')
  par(mfrow = c(2,1), mar = c(2,0,1,0.5), oma = c(2,4,0,0))
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
         ann = F, cex.axis = 0.85)
    lines(n ~ cv, data = settings, subset = nstrata == 5)
    
    abline(h = samples, col = 'grey', lty = 'dashed')
    text(x = plot_settings$xlabel[irow], y = samples, 
         paste(1:3, 'Boat'), pos = 1)
    legend('top', legend = paste(plot_settings$type[irow], 'Optimization'), 
           bty = 'n')
    
  }
  
  
  mtext(side = 1, 'Upper CV Constraint', outer = T, line = 0.5)
  mtext(side = 2, 'Total Optimized Sample Size', outer = T, line = 2.5)
  
  dev.off()
}


#####################################################
## Supplemental Figure
#####################################################
{png(paste0(figure_dir, 'Supplemental_Figures/SFig2_N_CV_Tradeoff.png'),
    width = 140, height = 180, res = 500, units = 'mm')
par(mfcol = c(3,2), mar = c(2,3,1,1), oma = c(3,3,2,0))
for(irow in 1:2){
  ##########################
  ## Load Data
  ##########################
  load(paste0(github_dir, plot_settings$type[irow], '_Optimization/',
              plot_settings$data_filename[irow], 'optimization_results.RData'))
  for(istrata in c(5,10,15)){
    
    #############################
    ## Plot
    #############################
    plot(n ~ cv, data = settings, subset = nstrata == istrata, las=1, pch = 16,
         ylim = c(0, plot_settings$ymax[irow]), 
         xlim = c(plot_settings$xmin[irow], plot_settings$xmax[irow]),
         ann = F)
    #xlab = 'Upper Spatiotemporal CV Constraint', ylab = 'Total Sample Size')
    lines(n ~ cv, data = settings, subset = nstrata == istrata)
    
    abline(h = samples, col = 'grey', lty = 'dashed')
    text(x = plot_settings$xlabel[irow], y = samples, 
         paste(1:3, 'Boat'), pos = 1)
    
    text(x = plot_settings$xlabel[irow],
         y = plot_settings$ymax[irow],
         paste(istrata, 'Strata'), font = 2, pos = 1, cex = 1.25)
    
    if(istrata == 5) mtext(side = 3, 
                           paste(plot_settings$type[irow], 'Optimization'),
                           line = 1)
  }
}

mtext(side = 1, 'Upper CV Constraint', outer = T, line = 1)
mtext(side = 2, 'Total Sample Size', outer = T, line = 1)
dev.off()
}
