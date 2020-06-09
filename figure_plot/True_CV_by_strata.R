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


stratas = c(5,10,15,20,30,40,50,60)
Nstratas = length(stratas)

#########################################
## Plot 
#########################################

for(itype in c('Spatiotemporal', 'Spatiotemporal_Flexible')){
  
  if(itype == 'Spatiotemporal')  {
    load(paste0(github_dir, 
                'Spatiotemporal_Optimization/',
                'STRS_Sim_Res_Spatiotemporal.RData'))
  }
  if(itype == 'Spatiotemporal_Flexible') {
    load(paste0(github_dir, 
                'Spatiotemporal_Optimization_Scheme2/',
                'STRS_Sim_Res_Spatiotemporal_Flexible.RData'))
  }
  
  png(file = paste0(PP_dir, 'True_CV_', itype, '.png'),
      width = 190, height = 190, units = 'mm', res = 1000)
  
  #Set up panel layout
  layout(mat = matrix(1:(5*5), nrow = 5, byrow = T),
         widths = c(1,0.25,1,0.25,1))
  
  par(mar = c(0,0,0,0), oma = c(5,5,2,0.5))
  for(ispp in 1:ns){
    max_val = max(0.11, 1.1*max(STRS_true_cv_array[,ispp,,], na.rm = T))
    plot(1, type = 'n', ylim = c(0,max_val), xlim = c(3,30), axes = F, ann = F)
    box()
    abline(v = seq(from = 6.5, by = 5, length = 6), col = 'darkgrey', lty = 'dashed')
    mtext(side = 3, sci_names[ispp], line = -1, font =3, cex = 0.75)
    
    axis(side = 2, las = 1, at = seq(0.1,1,0.1) )
    if(max_val < 0.2) axis(side = 2, las = 1, at = seq(0.025,0.075,0.025) )
    
    if(ispp %in% 13:15) axis(side = 1, at = seq(from = 4, by = 5, length = 6), 
                             labels = c(5,10,15,20,30,60))
    
    
    offset = 1
    
    if(itype == 'Spatiotemporal_Flexible') stratas = 1:6
    if(itype == 'Spatiotemporal') stratas = c(1:5,8)
    
    for(istrata in stratas ){
      boxplot(STRS_true_cv_array[,ispp,istrata,], add = T, axes = F, pch = 16,
              cex = 0.75, 
              at = (offset*5-2):(offset*5), col = c('red', 'blue', 'white'))
      
      offset = offset + 1
      
    }
    
    if(ispp%%3 != 0) plot(1, type = 'n', axes = F, ann = F)
  }
  mtext(side = 2, 'True CV', outer = T, line = 3)
  mtext(side = 1, 'Number of Strata', outer = T, line = 3)
  dev.off()
}
