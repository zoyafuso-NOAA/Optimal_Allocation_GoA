##################################
## Compare Optimized versus Current Survey Designs
##################################
rm(list = ls())

library(VAST); library(sp); library(raster)

##################################
## Set up directories
##################################
which_machine = c('Zack_PC' =1, 'Zack_GI_PC'=2)[1]
modelno = '6g'

github_dir = paste0(c('C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/')
VAST_dir = paste0(c('C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output', modelno, '/')
PP_dir = paste0(c('C:/Users/Zack Oyafuso/')[which_machine],
                'Google Drive/MS_Optimizations/powerpoint_plot/')

################################
## Load Data
################################
load(paste0(github_dir, 'Optimum_Allocation/model_', modelno, 
            '_spatiotemporal/Survey_Simulation_Results.RData'))
load(paste0(github_dir, 'Optimum_Allocation/model_', modelno, 
            '_spatiotemporal/Stratified_RS_Simulation_Results.RData'))
load(paste0(github_dir, 'Optimum_Allocation/model_', modelno, 
            '_spatiotemporal/optimization_results.RData'))

settings$id = 1:nrow(settings)
sub_settings = subset(settings, nstrata ==10)

#########################################
## Plot 
#########################################

{
  png(file = paste0(PP_dir, 'Survey_Comparison.png'),
      width = 10, height = 6, units = 'in', res = 1000)
  
  #Set up panel layout
  layout(mat = matrix(c(rep(1,6), 2:31), nrow = 6, byrow = T),
         heights = c(0.5,1,1,1,1,1))
  
  #Legend
  par(oma = c(4,2,0,2), mar = c(2,0,0,0))
  plot(1, xlim = c(0,10), ylim = c(0,1), axes = F)
  legend('center', legend = paste(1:3, c('Boat', 'Boats', 'Boats')),
         fill = c('red', 'blue', 'white'), horiz = T, cex = 2, bty = 'n')
  
  for(ispp in 1:ns){
    
    #Plot the True CV
    par(mar = c(0.5,4,0,0))
    plot(1, type = 'n', xlim = c(0,5), axes = F, ann = F, 
         ylim = c(0, c(0.1, 0.31, 0.2,
                       0.15, 0.16, 0.11, 
                       0.22, 0.32, 0.55,
                       0.32, 0.42, 0.24,
                       0.55, 0.65, 0.28)[ispp]))
    box()
    axis(side = 2, las = 1, at = seq(0,2,0.1))
    
    
    for(isamp in 1:3){
      isol=sub_settings$id[which.min(abs(sub_settings$n-c(280,550,820)[isamp]))]
      
      boxplot(STRS_true_cv_array[,ispp,isol],  at = c(0.5,1,1.5)[isamp], 
              add = T, axes = F, col = c('red', 'blue', 'white')[isamp])
      boxplot(Survey_true_cv_array[,ispp,isamp], at = c(3.5,4,4.5)[isamp], 
              add = T, axes = F, col = c('red', 'blue', 'white')[isamp])
    }
    if(ispp %in% 1:3) mtext(side = 3, 'True CV', cex = 0.9)
    
    if(ispp %in% c(13:15)) {
      axis(side = 1, at = c(1, 4), cex.axis = 1, line = 1.5, tick = F,
           labels = c('Optimized\nSurvey\nDesign', 'Current\nSurvey\nDesign'))
      axis(side = 1, at = c(1, 4), labels = NA, tick = T)
    }
    
    #Plot the RRMSE of the CV
    par(mar = c(0.5,0,0,4))
    plot(1, type = 'n', xlim = c(0,5), axes = F, ann = F, 
         ylim = c(0, c(0.6, 1.2, 0.8,
                       0.8, 0.9, 0.8, 
                       0.65, 0.8, 0.8,
                       0.8, 0.8, 0.6,
                       1.6, 1.8, 0.45)[ispp]))
    box()
    axis(side = 4, las = 1, at = seq(0,2,0.25))
    
    for(isamp in 1:3){  
      isol=sub_settings$id[which.min(abs(sub_settings$n-c(280,550,820)[isamp]))]
      
      boxplot(STRS_rrmse_cv_array[,ispp,isol],  at =c(0.5,1,1.5)[isamp], 
              add = T, axes = F, col = c('red', 'blue', 'white')[isamp])
      boxplot(Survey_rrmse_cv_array[,ispp,isamp], at = c(3.5,4,4.5)[isamp], 
              add = T, axes = F, col = c('red', 'blue', 'white')[isamp])
    }
    if(ispp %in% 1:3) mtext(side = 3, 'RRMSE of CV', cex = 0.9)
    
    #Species Label filler
    legend(x = 0, cex = 0.9, box.col='white', 
           y = 1.03*c(0.6, 1.2, 0.8,
                       0.8, 0.9, 0.8, 
                       0.65, 0.8, 0.8,
                       0.8, 0.8, 0.6,
                       1.6, 1.8, 0.45)[ispp], 
           legend = sci_names[ispp], xpd = NA, text.font = 3, xjust = 0.5)
    
    if(ispp %in% c(13:15)) {
      axis(side = 1, at = c(1, 4), cex.axis = 1, line = 1.5, tick = F,
           labels = c('Optimized\nSurvey\nDesign', 'Current\nSurvey\nDesign'))
      axis(side = 1, at = c(1, 4), labels = NA, tick = T)
    }
    
  }
  mtext(side = 2, 'True CV', outer = T)
  mtext(side = 4, 'RRMSE of CV', outer = T)
  
  dev.off()
}
