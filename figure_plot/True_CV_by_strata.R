##################################
## Compare Optimized versus Current Survey Designs
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
load(paste0(github_dir, 'Spatiotemporal_Optimization/',
            'STRS_Sim_Res_Spatiotemporal.RData'))
load(paste0(github_dir, 'Spatiotemporal_Optimization/',
            'spatiotemporal_optimization_results.RData'))

#########################################
## Plot 
#########################################

{
  png(file = paste0(PP_dir, 'True_CV_Spatiotemporal.png'),
      width = 190, height = 190, units = 'mm', res = 1000)
  
  #Set up panel layout
  layout(mat = matrix(1:56, nrow = 8, byrow = T),
         widths = c(1,1,1,0.4,1,1,1))
  
  par(mar = c(0,0,0,0), oma = c(5,5,2,0.5))
  for(ispp in 1:ns){
    for(isample in 1:3){
      max_val = max(STRS_true_cv_array[,ispp,,])
      boxplot(STRS_true_cv_array[,ispp,,isample], 
              xlim = c(0,9), ylim = c(0, 1.1*max_val), axes = F)
      box()
      if(isample == 1) axis(side = 2, las = 1, 
                            at = pretty(c(0,max_val), 2))
      
      if(ispp %in% 14:15) {
        if(isample %in% c(1,3)){
          axis(side = 1, at = 1:8, labels = c(5,NA,15,NA,25,NA,40,NA))
        }
        if(isample ==2){
          axis(side = 1, at = 1:8, labels = c(NA,10,NA,20,NA,30,NA,50))
        }
      }
      
      if(isample == 2)  mtext(side = 3, sci_names[ispp], line = -1,
                              font =3, cex = 0.5)
      
      if(ispp %in% 1:2) mtext(side = 3, paste(isample, 'Boat'), line = 0.5,
                              font = 2)
    }
    
    if(ispp%%2 == 1) plot(1, type = 'n', axes = F, ann = F)
  }
  mtext(side = 2, 'True CV', outer = T, line = 3)
  mtext(side = 1, 'Number of Strata', outer = T, line = 3)
  dev.off()
}
