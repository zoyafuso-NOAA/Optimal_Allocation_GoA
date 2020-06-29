#######################################
## Figures for Survey Comparison Analysis
#######################################
rm(list = ls())

#######################################
## Set up directories
#######################################

which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2)[2]

github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/',
                      'C:/Users/zack.oyafuso/Work/')[which_machine],
                    'GitHub/Optimal_Allocation_GoA/')

# PP_dir = paste0(c('/Users/zackoyafuso/Google Drive/', 
#                   'C:/Users/Zack Oyafuso/Google Drive/')[which_machine],
#                 'MS_Optimizations/powerpoint_plot/')

figure_dir = paste0(c('/Users/zackoyafuso/Google Drive/', 
                      'C:/Users/Zack Oyafuso/Google Drive/')[which_machine],
                    'MS_Optimizations/TechMemo/figures/')


##################################
## Load Results
##################################
load( paste0(github_dir, 'data/optimization_data.RData') )
load( paste0(github_dir, 'Simulate_Current_Survey/',
             'Survey_Simulation_Results.RData') )
load( paste0(github_dir, 'Spatiotemporal_Optimization_Scheme2/',
             'STRS_Sim_Res_Spatiotemporal_Flexible.RData') )

####################################
## Set to 10 stratas
####################################
istrata = 2

####################################
## Bias in Estimate
####################################

temp_bias = array(dim = c(11,15,3,2))
for(ispp in 1:15){
  # plot(1, type = 'n', ylim = c(-100,100), xlim = c(0,8), las = 1)
  for(isample in 1:3){
    
    #Survey
    abs_bias = sweep(x = Survey_sim_mean[,ispp,,isample], MARGIN = 1,
                     STATS = true_mean[,ispp], FUN = '-')
    rel_bias = 100*sweep(x = abs_bias, MARGIN = 1,
                         STATS = true_mean[,ispp], FUN = '/')
    
    temp_bias[,ispp,isample,1] = apply(rel_bias, MARGIN = 1, FUN = mean)
    
    #Optimized Survey
    abs_bias = sweep(x = STRS_sim_mean[,ispp,istrata,isample,], MARGIN = 1,
                     STATS = true_mean[,ispp], 
                     FUN = '-')
    rel_bias = 100*sweep(x = abs_bias, MARGIN = 1,
                         STATS = true_mean[,ispp], 
                         FUN = '/')
    temp_bias[,ispp,isample,2] = apply(rel_bias, MARGIN = 1, FUN = mean)
    
  }
  
}

{
  png(filename = paste0(figure_dir, 'Bias_Est.png'),
      units = 'mm', width = 190, height = 150, res = 500)
  par(mfrow = c(5,3), mar = c( 0.5,4,0.5,0), oma = c(2,1,2,0.5))
  for(ispp in 1:15){
    ymax = max(abs(temp_bias[,ispp,,]))
    plot(1, type = 'n', ylim = c(-ymax,ymax), xlim = c(0,8), las = 1, 
         axes = F, ann = F)
    if(ispp == 2) legend(x = -2, y = 8, col = c('red','blue','black'),xpd = NA,
                         pch = 0, legend = paste(1:3, 'Boat'), horiz = T, 
                         bty = 'n', cex = 1.5, lty = 1, x.intersp = 0.25,
                         text.col = c('red','blue','black'))
    abline(h=0, lty = 'dotted')
    box()
    axis(side = 2, las = 1)
    legend('topleft', sci_names[ispp], text.font = 3, bty = 'n')
    if(ispp %in% 13:15) 
      axis(side = 1, at = c(2,6), labels = c('Current Design', 
                                             'Optimized Design'), 
           cex.axis = 0.95)
    
    boxplot(temp_bias[,ispp,,1], add = T, at = 1:3, axes = F, 
            border = c('red', 'blue', 'black'))
    boxplot(temp_bias[,ispp,,2], add = T, at = 5:7, axes = F,
            border = c('red', 'blue', 'black'))
  }
  mtext(side = 2, 'Relative Percent Bias', outer = T, line = -.5)
  dev.off()
}

####################################
## Bias in CV
####################################
temp_bias = array(dim = c(11,15,3,2))

for(ispp in 1:15){
  # plot(1, type = 'n', ylim = c(-100,100), xlim = c(0,8), las = 1)
  for(isample in 1:3){
    
    #Survey
    abs_bias = sweep(x = Survey_sim_cv[,ispp,,isample], MARGIN = 1,
                     STATS = Survey_true_cv_array[,ispp,isample], FUN = '-')
    rel_bias = 100*sweep(x = abs_bias, MARGIN = 1,
                         STATS = Survey_true_cv_array[,ispp,isample], FUN = '/')
    
    temp_bias[,ispp,isample,1] = apply(rel_bias, MARGIN = 1, FUN = mean)
    
    #Optimized Survey
    abs_bias = sweep(x = STRS_sim_cv[,ispp,istrata,isample,], MARGIN = 1,
                     STATS = STRS_true_cv_array[,ispp,istrata,isample], 
                     FUN = '-')
    rel_bias = 100*sweep(x = abs_bias, MARGIN = 1,
                         STATS =STRS_true_cv_array[,ispp,istrata,isample], 
                         FUN = '/')
    temp_bias[,ispp,isample,2] = apply(rel_bias, MARGIN = 1, FUN = mean)
    
  }
  
}

{
  png(filename = paste0(figure_dir, 'Bias_CV.png'),
      units = 'mm', width = 190, height = 150, res = 500)
  par(mfrow = c(5,3), mar = c( 0.5,4,0.5,0), oma = c(2,1,2,0.5))
  for(ispp in 1:15){
    ymax = max(abs(temp_bias[,ispp,,]))
    plot(1, type = 'n', ylim = c(-ymax,ymax), xlim = c(0,8), las = 1, 
         axes = F, ann = F)
    if(ispp == 2) legend(x = -2, y = 7, col = c('red','blue','black'), xpd = NA,
                         pch = 0, legend = paste(1:3, 'Boat'), horiz = T, 
                         bty = 'n', cex = 1.5, lty = 1, x.intersp = 0.25,
                         text.col = c('red', 'blue', 'black'))
    abline(h=0, lty = 'dotted')
    box()
    axis(side = 2, las = 1)
    legend('bottom', sci_names[ispp], text.font = 3, bty = 'n')
    if(ispp %in% 13:15) 
      axis(side = 1, at = c(2,6), labels = c('Current Design', 
                                             'Optimized Design'), 
           cex.axis = 0.95)
    
    boxplot(temp_bias[,ispp,,1], add = T, at = 1:3, axes = F, 
            border = c('red', 'blue', 'black'))
    boxplot(temp_bias[,ispp,,2], add = T, at = 5:7, axes = F,
            border = c('red', 'blue', 'black'))
  }
  mtext(side = 2, 'Relative Percent Bias', outer = T, line = -.5)
  dev.off()
}

####################################
## True CV
####################################
{
  png(filename = paste0(figure_dir, 'True_CV.png'),
      units = 'mm', width = 190, height = 150, res = 500)
  par(mfrow = c(5,3), mar = c( 0.5,4,0.5,0), oma = c(2,1,2,0.5))
  for(ispp in 1:15){
    ymax = max(c(Survey_true_cv_array[,ispp,],
                 STRS_true_cv_array[,ispp,istrata,]))
    
    plot(1, type = 'n', ylim = c(0,ymax), xlim = c(0,8), las = 1,
         axes = F, ann = F)
    if(ispp == 2) legend(x = -2, y = 0.35, col = c('red','blue','black'), xpd = NA,
                         pch = 0, legend = paste(1:3, 'Boat'), horiz = T, 
                         bty = 'n', cex = 1.5, lty = 1, x.intersp = 0.25,
                         text.col = c('red', 'blue', 'black'))
    
    axis(side = 2, las = 1)
    if(ispp %in% 13:15) 
      axis(side = 1, at = c(2,6), labels = c('Current Design', 
                                             'Optimized Design'), 
           cex.axis = 0.95)
    legend('bottom', sci_names[ispp], bty = 'n', text.font = 3)
    box()
    axis(side = 2, las = 1)
    
    for(isample in 1:3){
      
      boxplot(Survey_true_cv_array[,ispp,isample], add = T, at = isample,
              axes = F, border = c('red', 'blue', 'black')[isample])
      
      boxplot(STRS_true_cv_array[,ispp,istrata,isample], add = T, 
              at = 4+isample, axes = F, 
              border = c('red', 'blue', 'black')[isample])
    }
  }
  mtext(side = 2, 'True CV', outer = T, line = -.5)
  dev.off()
}

####################################
## RRMSE of CV
####################################

{
  png(filename = paste0(figure_dir, 'RRMSE_CV.png'),
      units = 'mm', width = 190, height = 150, res = 500)
  par(mfrow = c(5,3), mar = c( 0.5,4,0.5,0), oma = c(2,1,2,0.5))
  for(ispp in 1:15){
    ymax = max(c(Survey_rrmse_cv_array[,ispp,],
                 STRS_rrmse_cv_array[,ispp,istrata,]))
    
    plot(1, type = 'n', ylim = c(0,ymax), xlim = c(0,8), las = 1,
         axes = F, ann = F)
    if(ispp == 2) legend(x = -2, y = 1.2, col = c('red','blue','black'), 
                         pch = 0, legend = paste(1:3, 'Boat'), horiz = T, 
                         bty = 'n', cex = 1.5, lty = 1, x.intersp = 0.25,
                         text.col = c('red', 'blue', 'black'), xpd = NA)
    
    axis(side = 2, las = 1)
    if(ispp %in% 13:15) 
      axis(side = 1, at = c(2,6), labels = c('Current Design', 
                                             'Optimized Design'), 
           cex.axis = 0.95)
    legend('topright', sci_names[ispp], bty = 'n', text.font = 3)
    box()
    axis(side = 2, las = 1)
    
    for(isample in 1:3){
      
      boxplot(Survey_rrmse_cv_array[,ispp,isample], add = T, at = isample,
              axes = F, border = c('red', 'blue', 'black')[isample])
      boxplot(STRS_rrmse_cv_array[,ispp,istrata,isample], add = T, 
              at = 4+isample, axes = F, 
              border = c('red', 'blue', 'black')[isample])
    }
  }
  mtext(side = 2, 'RRMSE of CV', outer = T, line = -.5)
  dev.off()
}
