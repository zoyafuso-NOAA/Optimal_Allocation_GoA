##################################
## Bias of Estimates
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

istrata = 2
for(itype in 1:2){ #1: spatiotemporal, #2: flexible spatiotemporal
  
  load(paste0(github_dir, 
              c('Spatiotemporal_Optimization/', 
                'Spatiotemporal_Optimization_Scheme2/')[itype],
              c('STRS_Sim_Res_Spatiotemporal.RData', 
                'STRS_Sim_Res_Spatiotemporal_Flexible.RData')[itype] ))
  
  png(file = paste0(PP_dir, paste0('Bias_', 
                                   c('Spatiotemporal', 
                                     'Spatiotemporal_Felxible')[itype], 
                                   '.png')),
      width = 190, height = 220, units = 'mm', res = 1000)
  
  layout(mat = matrix(1:(15*7), nrow = 15, byrow = T),
         widths = c(1,1,1,0.4,1,1,1))
  par(mar = c(0,0,0,0), oma = c(2,5,2,0.5))
  
  for(ispp in 1:ns){
    for(itype in c('Est', 'CV')){
      for(isample in 1:3){
        
        if(itype == 'CV'){
          abs_bias = sweep(x = STRS_sim_cv[,ispp,istrata, isample,], MARGIN = 1,
                           STATS = STRS_true_cv_array[,ispp,istrata,isample],
                           FUN = '-')
          rel_bias = 100*sweep(x = abs_bias, MARGIN = 1,
                               STATS = STRS_true_cv_array[,ispp,istrata,isample],
                               FUN = '/')
        }
        
        if(itype == 'Est'){
          abs_bias = sweep(x = STRS_sim_mean[,ispp,istrata, isample,], MARGIN = 1,
                           STATS = true_mean[,ispp], FUN = '-')
          rel_bias = 100*sweep(x = abs_bias, MARGIN = 1,
                               STATS = true_mean[,ispp], FUN = '/')
        }
        
        bias_box = apply(rel_bias, MARGIN = 1,
                         FUN = function(x) quantile(x, 
                                                    probs = c(0.05, 0.25, 0.5,
                                                              0.75, 0.95)))
        #ymax = bias_box
        plot(1, type = 'n', ann = F, axes = F, 
             ylim = c(-1.5,1.5)*max(max(bias_box),25), xlim = c(0,1))
        box()
        abline(h = 0, col = 'darkgrey', lty = 'dashed')
        if(ispp == 1) mtext(side = 3, paste(isample, 'Boat'), line = 0.5,
                            font = 2, cex = 0.75)
        if(isample == 1) axis(side = 2, las = 1, at = seq(-50,50,by=25))
        if(isample == 2) mtext(side = 3, sci_names[ispp], line = -1,
                               font =3, cex = 0.5)
        if(ispp %in% 15){
          if(isample %in% c(1,3)) axis(side = 1, at = seq(0,1,length=11),
                                       label = c('Year 1',rep(NA,9),'Year 11'))
        }
        lines(y = bias_box['50%',], x = seq(0,1,length = 11))
        points(y = bias_box['50%',], x = seq(0,1,length = 11), pch = 16)
        lines(y = bias_box['25%',], x = seq(0,1,length = 11))
        lines(y = bias_box['75%',], x = seq(0,1,length = 11))
        lines(y = bias_box['5%',], x = seq(0,1,length = 11), lty = 'dotted')
        lines(y = bias_box['95%',], x = seq(0,1,length = 11), lty = 'dotted')
      }
      if(itype == 'Est') plot(1, type = 'n', axes = F, ann = F)
    }
  }
  mtext(side = 2, 'Relative Bias (%)', outer = T, line = 3)
  dev.off()
}
