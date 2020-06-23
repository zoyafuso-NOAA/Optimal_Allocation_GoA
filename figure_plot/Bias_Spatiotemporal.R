##################################
## Bias of Estimates
##################################
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

################################
## Load Data
################################
load(paste0(github_dir, 'data/optimization_data.RData'))

{
  png(file = paste0(PP_dir, paste0('Bias_Est.png')),
      width = 190, height = 150, units = 'mm', res = 1000)
  
  par(mar = c(0,4,0,0), oma = c(4,1,2,0.5), mfrow = c(5,3))
  
  for(itype in 2){
    load(paste0(github_dir, 
                c('Spatiotemporal_Optimization/', 
                  'Spatiotemporal_Optimization_Scheme2/')[itype],
                c('STRS_Sim_Res_Spatiotemporal.RData', 
                  'STRS_Sim_Res_Spatiotemporal_Flexible.RData')[itype] ))
    
    if(itype == 1) stratas = c(1,2,3,4,7,8)
    if(itype == 2) stratas = c(1:6)
    
    for(ispp in 1:ns){
      
      ylim_ = list(c(-60,60), c(-80,110), c(-60,70), 
                   c(-70,70), c(-70,70), c(-30,30), 
                   c(-80,110),c(-80,110),c(-130,130),
                   c(-110,160),c(-110,210),c(-110,130),
                   c(-110,160),c(-110,210),c(-110,110))[[ispp]]
      
      plot(1, type = 'n', xlim = c(-0.5,13.5), axes = F, ann = F, ylim=ylim_ )
      box(); 
      abline(v = seq(from = 1.75, by = 2.5, length = 5), lty = 'dashed', 
             col = 'lightgrey')
      abline(h = 0)
      axis(side = 2, las = 1, at = pretty(ylim_, 4) )
      if(ispp == 2) legend(x = -3, y = 190, legend = paste(1:3, 'Boat'),
                           fill = c('red', 'blue', 'white'), x.intersp = .5,
                           horiz = T, xpd = NA, cex = 1.5, bty = 'n')
      
      if(itype == 2) legend('bottom', sci_names[ispp], bty = 'n', 
                            text.font = 3 )
      if(ispp %in% c(13:15)) axis(side = 1, labels = c(5,10,15,20,30,60), 
                                  at = seq(from=0.5, by=2.5, length=6))
      offset = 0
      for(istrata in stratas){
        for(isample in 1:3){
          abs_bias = sweep(x = STRS_sim_mean[,ispp,istrata,isample,], MARGIN=1,
                           STATS = true_mean[,ispp],
                           FUN = '-')
          
          rel_bias = 100*sweep(x = abs_bias, MARGIN = 1,
                               STATS = true_mean[,ispp],
                               FUN = '/')
          
          boxplot( as.vector(rel_bias), add = T, axes = F, at = offset, 
                   col = c('red', 'blue', 'white')[isample],
                   pch = 16, cex = 0.25)
          offset = offset + 0.5
        }
        offset = offset + 1
      }
    }
  }
  
  # plot(1, type = 'n', axes = F, ann = F)
  # 
  # plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), axes = F, ann = F)
  # legend('bottom', legend = paste0(1:3, ' Boat'), horiz = T, cex = 1.5,
  #        fill = c('red', 'blue', 'white'))
  
  mtext(side = 1, 'Number of Strata', outer = T, line = 2.5)
  mtext(side = 2, 'Percent Relative Bias', outer = T, line = -1)
  dev.off()
}

##############################
## What proportion of surveys > 50% bias
##############################
rel_bias = array(dim = c(ns, 6, 3))
for(itype in 2){
  load(paste0(github_dir, 
              c('Spatiotemporal_Optimization/', 
                'Spatiotemporal_Optimization_Scheme2/')[itype],
              c('STRS_Sim_Res_Spatiotemporal.RData', 
                'STRS_Sim_Res_Spatiotemporal_Flexible.RData')[itype] ))
  
  if(itype == 1) stratas = c(1,2,3,4,7,8)
  if(itype == 2) stratas = c(1:6)
  
  for(ispp in 1:ns){
    for(istrata in stratas){
      for(isample in 1:3){
        abs_bias = sweep(x = STRS_sim_mean[,ispp,istrata,isample,], MARGIN = 1,
                         STATS = true_mean[,ispp],
                         FUN = '-')
        
        temp_rel_bias = 100*sweep(x = abs_bias, MARGIN = 1,
                                  STATS = true_mean[,ispp],
                                  FUN = '/')
        
        rel_bias[ispp,istrata,isample] = sum(as.vector(abs(temp_rel_bias))>50)
        
      }
    }
  }
}
rel_bias = round(rel_bias / (11000) * 100, 2)

#########################

{
  png(file = paste0(PP_dir, paste0('Bias_CV.png')),
      width = 190, height = 150, units = 'mm', res = 1000)
  
  layout(mat = matrix(c(1:4, 
                        8:11,
                        rep(16, 4),
                        5:7,15,
                        12:15), nrow = 4), 
         widths = c(1,1,0.25,1,1))
  par(mar = c(0,0,0,0), oma = c(4,4,3,1))
  
  for(itype in 1:2){
    load(paste0(github_dir,
                c('Spatiotemporal_Optimization/',
                  'Spatiotemporal_Optimization_Scheme2/')[itype],
                c('STRS_Sim_Res_Spatiotemporal.RData',
                  'STRS_Sim_Res_Spatiotemporal_Flexible.RData')[itype] ))
    
    if(itype == 1) stratas = c(1,2,3,4,7,8)
    if(itype == 2) stratas = c(1:6)
    
    for(ispp in c(1,3,11,13, 8,12,15) ){
      
      ylim_ = list(c(-80,150), c(-100,100), c(-80,120), 
                   c(-100,100), c(-100,100), c(-100,100), 
                   c(-100,100), c(-80,120), c(-100,100), 
                   c(-100,100), c(-100,120), c(-110,120),
                   c(-100,140), c(-100,100), c(-100,160) )[[ispp]]
      
      plot(1, type = 'n', xlim = c(-0.5,13.5), axes = F, ann = F, ylim=ylim_ )
      box(); abline(h = 0, lty = 'dashed')
      abline(v = seq(from = 1.75, by = 2.5, length = 5), lty = 'dashed', 
             col = 'darkgrey')
      if(itype == 1) axis(side = 2, las = 1, at = pretty(ylim_, 5) )
      
      if(itype == 2) legend('bottom', sci_names[ispp], bty = 'n', 
                            text.font = 3  )
      if(ispp %in% c(13,15)) axis(side = 1, at = seq(0.5, by=2.5, length=6),
                                  labels = c(5,10,15,20,30,60))
      if(ispp %in% c(1,8)) mtext(side = 3,c('Constrained\nSpatiotemporal',
                                            'Flexible\nSpatiotemporal')[itype])
      
      offset = 0
      for(istrata in stratas){
        for(isample in 1:3){
          abs_bias = sweep(x = STRS_sim_cv[,ispp,istrata,isample,], MARGIN = 1,
                           STATS = STRS_true_cv_array[,ispp,istrata,isample],
                           FUN = '-')
          
          rel_bias = sweep(x = abs_bias, MARGIN = 1,
                           STATS=STRS_true_cv_array[,ispp,istrata,isample],
                           FUN = '/') * 100
          
          boxplot( as.vector(rel_bias), add = T, axes = F, at = offset, 
                   col = c('red', 'blue', 'white')[isample],
                   pch = 16, cex = 0.25)
          offset = offset + 0.5
        }
        offset = offset + 1
      }
    }
  }
  
  plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), axes = F, ann = F)
  legend('bottom', legend = paste0(1:3, ' Boat'), horiz = T, cex = 1.5,
         fill = c('red', 'blue', 'white'))
  plot(1, type = 'n', axes = F, ann = F)
  
  mtext(side = 1, 'Number of Strata', outer = T, line = 2.5)
  mtext(side = 2, 'Percent Relative Bias', outer = T, line = 2.5)
  dev.off()
}

###############################################3
## SUpplementary Figuers
##############################################

{
  png(file = paste0(figure_dir, 'Supplemental_Figures/SFig6_Bias_Est.png'),
      width = 190, height = 190, units = 'mm', res = 1000)
  
  par(mar = c(0,4,0,0), oma = c(4,1,3,0.5), mfrow = c(5,3))
  
  for(itype in 1){
    load(paste0(github_dir, 
                c('Spatiotemporal_Optimization/', 
                  'Spatiotemporal_Optimization_Scheme2/')[itype],
                c('STRS_Sim_Res_Spatiotemporal.RData', 
                  'STRS_Sim_Res_Spatiotemporal_Flexible.RData')[itype] ))
    
    if(itype == 1) stratas = c(1,2,3,4,7,8)
    if(itype == 2) stratas = c(1:6)
    
    for(ispp in 1:ns){
      
      ylim_ = list(c(-60,60), c(-80,110), c(-60,70), 
                   c(-70,70), c(-70,70), c(-30,30), 
                   c(-80,110),c(-80,110),c(-130,130),
                   c(-110,160),c(-110,210),c(-110,130),
                   c(-110,160),c(-110,210),c(-110,110))[[ispp]]
      
      plot(1, type = 'n', xlim = c(-0.5,13.5), axes = F, ann = F, ylim=ylim_ )
      box(); 
      abline(v = seq(from = 1.75, by = 2.5, length = 5), lty = 'dashed', 
             col = 'lightgrey')
      abline(h = 0)
      axis(side = 2, las = 1, at = pretty(ylim_, 4) )
      if(ispp == 2) legend(x = -3, y = 190, legend = paste(1:3, 'Boat'),
                           fill = c('red', 'blue', 'white'), x.intersp = .5,
                           horiz = T, xpd = NA, cex = 1.5, bty = 'n')
      
      if(itype == 1) legend('bottom', sci_names[ispp], bty = 'n', 
                            text.font = 3 )
      if(ispp %in% c(13:15)) axis(side = 1, labels = c(5,10,15,20,30,60), 
                                  at = seq(from=0.5, by=2.5, length=6))
      offset = 0
      for(istrata in stratas){
        for(isample in 1:3){
          abs_bias = sweep(x = STRS_sim_mean[,ispp,istrata,isample,], MARGIN=1,
                           STATS = true_mean[,ispp],
                           FUN = '-')
          
          rel_bias = 100*sweep(x = abs_bias, MARGIN = 1,
                               STATS = true_mean[,ispp],
                               FUN = '/')
          
          boxplot( as.vector(rel_bias), add = T, axes = F, at = offset, 
                   col = c('red', 'blue', 'white')[isample],
                   pch = 16, cex = 0.25)
          offset = offset + 0.5
        }
        offset = offset + 1
      }
    }
  }
  
  # plot(1, type = 'n', axes = F, ann = F)
  # 
  # plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), axes = F, ann = F)
  # legend('bottom', legend = paste0(1:3, ' Boat'), horiz = T, cex = 1.5,
  #        fill = c('red', 'blue', 'white'))
  
  mtext(side = 1, 'Number of Strata', outer = T, line = 2.5)
  mtext(side = 2, 'Percent Relative Bias', outer = T, line = -1)
  dev.off()
}


#########################

{
  png(file = paste0(figure_dir,  'Supplemental_Figures/SFig7_Bias_CV.png'),
      width = 190, height = 200, units = 'mm', res = 1000)
  
  layout(mat = matrix(c(1:8, 
                        16:23,
                        rep(32, 8),
                        9:15, 31,
                        24:30,31), nrow = 8), 
         widths = c(1,1,0.25,1,1))
  par(mar = c(0,0,0,0), oma = c(4,4,3,1))
  
  for(itype in 1:2){
    load(paste0(github_dir,
                c('Spatiotemporal_Optimization/',
                  'Spatiotemporal_Optimization_Scheme2/')[itype],
                c('STRS_Sim_Res_Spatiotemporal.RData',
                  'STRS_Sim_Res_Spatiotemporal_Flexible.RData')[itype] ))
    
    if(itype == 1) stratas = c(1,2,3,4,7,8)
    if(itype == 2) stratas = c(1:6)
    
    for(ispp in 1:15 ){
      
      ylim_ = list(c(-80,150), c(-110,110), c(-80,120), 
                   c(-110,170), c(-110,110), c(-110,110), 
                   c(-110,110), c(-80,120), c(-110,110), 
                   c(-110,170), c(-110,120), c(-110,170),
                   c(-110,170), c(-110,110), c(-110,170) )[[ispp]]
      
      plot(1, type = 'n', xlim = c(-0.5,13.5), axes = F, ann = F, ylim=ylim_ )
      box(); abline(h = 0, lty = 'dashed')
      abline(v = seq(from = 1.75, by = 2.5, length = 5), lty = 'dashed', 
             col = 'darkgrey')
      if(itype == 1) axis(side = 2, las = 1, at = pretty(ylim_, 5) )
      
      if(itype == 2) legend('bottom', sci_names[ispp], bty = 'n', 
                            text.font = 3  )
      if(ispp %in% c(8,15)) axis(side = 1, at = seq(0.5, by=2.5, length=6),
                                 labels = c(5,10,15,20,30,60))
      if(ispp %in% c(1,9)) mtext(side = 3,c('Constrained\nSpatiotemporal',
                                            'Flexible\nSpatiotemporal')[itype])
      
      offset = 0
      for(istrata in stratas){
        for(isample in 1:3){
          abs_bias = sweep(x = STRS_sim_cv[,ispp,istrata,isample,], MARGIN = 1,
                           STATS = STRS_true_cv_array[,ispp,istrata,isample],
                           FUN = '-')
          
          rel_bias = sweep(x = abs_bias, MARGIN = 1,
                           STATS=STRS_true_cv_array[,ispp,istrata,isample],
                           FUN = '/') * 100
          
          boxplot( as.vector(rel_bias), add = T, axes = F, at = offset, 
                   col = c('red', 'blue', 'white')[isample],
                   pch = 16, cex = 0.25)
          offset = offset + 0.5
        }
        offset = offset + 1
      }
    }
  }
  
  plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), axes = F, ann = F)
  legend('bottom', legend = paste0(1:3, ' Boat'), horiz = T, cex = 1.5,
         fill = c('red', 'blue', 'white'))
  plot(1, type = 'n', axes = F, ann = F)
  
  mtext(side = 1, 'Number of Strata', outer = T, line = 2.5)
  mtext(side = 2, 'Percent Relative Bias', outer = T, line = 2.5)
  dev.off()
}


