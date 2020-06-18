##################################
## True CV and RRMSE of CV across strata for each species for 
## both optimization types
##################################
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

################################
## Load Data
################################
load(paste0(github_dir, 'data/optimization_data.RData'))

################################
## True CV across years (boxplot), sample size (color), strata (x-axis),
## and species (plot panel) for the Flexible Optimization
################################
{
  png(file = paste0(PP_dir, paste0('True_CV.png')),
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
      
      max_val = c(0.12, 0.32, 0.12, 
                  0.14, 0.17, 0.07, 
                  0.23, 0.22, 0.35, 
                  0.32, 0.44, 0.27,
                  0.35, 0.475, 0.24)[ispp]
      
      plot(1, type = 'n', xlim = c(-0.5,13.5), axes = F, ann = F, 
           ylim = c(0, max_val) )
      box(); 
      abline(v = seq(from = 1.75, by = 2.5, length = 5), lty = 'dashed', 
             col = 'lightgrey')
      axis(side = 2, las = 1, at = pretty(c(0,max_val), 3) )
      
      if(ispp == 2) legend(x = -3, y = 0.45, legend = paste(1:3, 'Boat'),
                           fill = c('red', 'blue', 'white'), x.intersp = .5,
                           horiz = T, xpd = NA, cex = 1.5, bty = 'n')
      
      if(itype == 2) legend('bottom', sci_names[ispp], bty = 'n', 
                            text.font = 3 )
      if(ispp %in% c(13:15)) axis(side = 1, labels = c(5,10,15,20,30,60), 
                                  at = seq(from=0.5, by=2.5, length=6))
    
      offset = 0
      for(istrata in stratas){
        for(isample in 1:3){
          boxplot( STRS_true_cv_array[,ispp,istrata,isample], add = T,
                   axes = F, at = offset, pch = 16, cex = 0.5, 
                   col = c('red', 'blue', 'white')[isample] )
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
  mtext(side = 2, 'True CV', outer = T, line = -1)
  dev.off()
}

##############################

{
  png(file = paste0(PP_dir, paste0('RRMSE_CV.png')),
      width = 190, height = 120, units = 'mm', res = 1000)
  
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
    
    for(ispp in c(1,3,11,13, 8,12,15)){
      
      max_val = c(0.6, 0.75, 0.5, 0.55,
                  0.25, 0.50, 0.5, 0.43,
                  0.35, 0.80, 0.75, 0.75,
                  0.55, 1.10, 0.55)[ispp]
      
      plot(1, type = 'n', xlim = c(-0.5,13.5), axes = F, ann = F, 
           ylim = c(0, max_val) )
      
      box()
      abline(v = seq(from = 1.75, by = 2.5, length = 5), lty = 'dashed', 
             col = 'darkgrey')
      if(itype == 1) axis(side = 2, las = 1, at = pretty(c(0, max_val), 3))
      
      if(itype == 2) legend('topright', sci_names[ispp], bty = 'n', 
                            text.font = 3 )
      if(ispp %in% c(13,15)) axis(side = 1, at = seq(0.5, by = 2.5, length = 6),
                                 labels = c(5,10,15,20,30,60))
      if(ispp %in% c(1,8)) mtext(side = 3, c('Constrained\nSpatiotemporal',
                                             'Flexible\nSpatiotemporal')[itype])
      
      offset = 0
      for(istrata in stratas){
        for(isample in 1:3){
          boxplot( STRS_rrmse_cv_array[,ispp,istrata,isample], add = T,
                   axes = F, at = offset, pch = 16, cex = 0.5,
                   col = c('red', 'blue', 'white')[isample] )
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
  mtext(side = 2, 'RRMSE of CV', outer = T, line = 2.5)
  dev.off()
}
