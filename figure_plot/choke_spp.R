#############################
## How does the spatiotemporal CV relate to true CV
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


figure_dir = paste0(c('/Users/zackoyafuso/', 
                      'C:/Users/Zack Oyafuso/')[which_machine],
                    'Google Drive/MS_Optimizations/figure_plot/')

load(paste0(github_dir, 'data/optimization_data.RData'))

################
## Plot Settings
################
plot_settings = data.frame(
  type = c('Spatial', 'Spatiotemporal'),
  data_filename = c('spatial_only_', 
                    'spatiotemporal_'),
  subtitle = c('Spatial Only\n(Constrained Optimization)', 
               'Spatiotemporal\n(Constrained Optimization)')
)

rockfish_cod_idx = c(2,3,11:15)
flatfish_idx = (1:ns)[-rockfish_cod_idx]

####################
{png(filename = paste0(figure_dir, 'Fig5_choke_spp.png'),
    width = 190, height = 120, units = 'mm', res =500)
  
  par(mfrow = c(1,3), mar = c(3,0,3,0), oma = c(1,11,0,1))
  
  for(irow in 1:2){
    load(paste0(github_dir, plot_settings$type[irow], 
                '_Optimization/STRS_Sim_Res_',
                plot_settings$type[irow], '.RData'))
    load(paste0(github_dir,plot_settings$type[irow],'_Optimization/',
                plot_settings$data_filename[irow], 
                'optimization_results.RData'))
    
    sub_settings = subset(settings, nstrata == 10)
    sample_idx = which.min(abs(sub_settings$n - 550))
    
    abs_diff = STRS_true_cv_array[,,2,2] - sub_settings$cv[sample_idx]
    rel_diff = 100* abs_diff / sub_settings$cv[sample_idx]
    
    boxplot(rel_diff, horizontal = TRUE, add = F, axes = F,
            pch = 16, cex = 0.5, ylim = c(-110,160),
            main = plot_settings$subtitle[irow])
    box()
    abline(v = 0, col = 'darkgrey', lty = 'dashed')
    axis(side = 1)
    if(irow == 1) axis(side = 2, sci_names, las = 1, font = 3, at = 1:ns)
    
  }
  mtext(side = 1, paste0('Percent Difference of the True CV ',
                         'Relative to the Upper CV Constraint'), 
        outer = T, line = 0)
  
  #######################
  ## Flexible Scheme
  #######################
  load(paste0(github_dir, 'Spatiotemporal_Optimization_Scheme2/',
              'spatiotemporal_Flexible_optimization_results.RData'))
  load(paste0(github_dir, 'Spatiotemporal_Optimization_Scheme2/',
              'STRS_Sim_Res_Spatiotemporal_Flexible.RData'))
  settings$id = 1:nrow(settings)
  
  sub_settings = subset(settings, nstrata == 10)
  idx = which.min(abs(sub_settings$n - 550))
  expected_cv = unlist(sub_settings[idx-1,paste0('CV_',1:ns)] * 0.95)
  expected_cv = sapply(expected_cv, function(x) max(x, 0.1))
  abs_diff = sweep(x = STRS_true_cv_array[,,2,2], MARGIN = 2, 
                   STATS = expected_cv, FUN = '-' )
  rel_diff = 100*sweep(x = abs_diff, MARGIN = 2, 
                       STATS = expected_cv, FUN = '/' )
  
  boxplot(rel_diff, horizontal = TRUE, add = F, axes = F,
          pch = 16, cex = 0.5, ylim = c(-110,160),
          main =  'Spatiotemporal\n(Flexible Optimization)',
          border = ifelse(expected_cv == 0.1, 'blue', 'black' ))
  box()
  abline(v = 0, col = 'darkgrey', lty = 'dashed')
  axis(side = 1)
  
  dev.off()
  
}

############################################
## Supplementary Figure
############################################
{png(filename = paste0(figure_dir, 'Supplemental_Figures/SFig3_choke_spp.png'),
     width = 190, height = 200, units = 'mm', res =500)
  
  par(mfrow = c(3,3), mar = c(3,0,3,0), oma = c(1,12,0,1))
  
  for(istrata in c(1:3)){
    for(irow in 1:2){
      load(paste0(github_dir, plot_settings$type[irow], 
                  '_Optimization/STRS_Sim_Res_',
                  plot_settings$type[irow], '.RData'))
      load(paste0(github_dir,plot_settings$type[irow],'_Optimization/',
                  plot_settings$data_filename[irow], 
                  'optimization_results.RData'))
      
      sub_settings = subset(settings, nstrata == c(5,10,15)[istrata])
      sample_idx = which.min(abs(sub_settings$n - 550))
      
      abs_diff = STRS_true_cv_array[,,istrata,2] - sub_settings$cv[sample_idx]
      rel_diff = 100* abs_diff / sub_settings$cv[sample_idx]
      
      boxplot(rel_diff, horizontal = TRUE, add = F, axes = F,
              pch = 16, cex = 0.5, ylim = c(-110,160),
              main = plot_settings$subtitle[irow])
      box()
      abline(v = 0, col = 'darkgrey', lty = 'dashed')
      axis(side = 1)
      if(irow == 1) axis(side = 2, sci_names, las = 1, font = 3, at = 1:ns)
      
    }
    mtext(side = 1, paste0('Percent Difference of the True CV ',
                           'Relative to the Upper CV Constraint'), 
          outer = T, line = 0)
    
    #######################
    ## Flexible Scheme
    #######################
    load(paste0(github_dir, 'Spatiotemporal_Optimization_Scheme2/',
                'spatiotemporal_Flexible_optimization_results.RData'))
    load(paste0(github_dir, 'Spatiotemporal_Optimization_Scheme2/',
                'STRS_Sim_Res_Spatiotemporal_Flexible.RData'))
    settings$id = 1:nrow(settings)
    
    sub_settings = subset(settings, nstrata == c(5,10,15)[istrata])
    idx = which.min(abs(sub_settings$n - 550))
    expected_cv = unlist(sub_settings[idx-1,paste0('CV_',1:ns)] * 0.95)
    expected_cv = sapply(expected_cv, function(x) max(x, 0.1))
    abs_diff = sweep(x = STRS_true_cv_array[,,istrata,2], MARGIN = 2, 
                     STATS = expected_cv, FUN = '-' )
    rel_diff = 100*sweep(x = abs_diff, MARGIN = 2, 
                         STATS = expected_cv, FUN = '/' )
    
    boxplot(rel_diff, horizontal = TRUE, add = F, axes = F,
            pch = 16, cex = 0.5, ylim = c(-110,160),
            main =  'Spatiotemporal\n(Flexible Optimization)',
            border = ifelse(expected_cv == 0.1, 'blue', 'black' ))
    box()
    abline(v = 0, col = 'darkgrey', lty = 'dashed')
    axis(side = 1)
    
    legend('bottomright', legend = paste(c(5,10,15)[istrata], 'Strata'),
           bty = 'n', cex = 1.75)
  }
  
  
  dev.off()
  
}

