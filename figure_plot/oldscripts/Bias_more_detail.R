##################################
## Distribution of Bias of Estimates across Simulated Surveys
## 
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
                    "GitHub/Optimal_Allocation_GoA/Spatiotemporal_Optimization_Scheme2/")

PP_dir = paste0(c('/Users/zackoyafuso/', 
                  'C:/Users/Zack Oyafuso/')[which_machine],
                'Google Drive/MS_Optimizations/powerpoint_plot/')

###########################
## 
###########################
load(paste0(github_dir, 'STRS_Sim_Res_Spatiotemporal_Flexible.RData'))
load( paste0(dirname(github_dir), '/data/optimization_data.RData') )

pdf(file = paste0(PP_dir, 'Bias_CV_distribution.pdf'), width = 8, height = 10, onefile = T)
par(mfrow = c(6,3), mar = c(0,0,0,0), oma = c(5,1,6,1))
for(ispp in 1:ns){
 for(istrata in 1:6){
  for(isample in 1:3){
   plot(1, type = 'n', ylim = c(0,1), xlim = c(-75,100), ann = F, axes = F)
   abline(v = 0, lty = 'dashed')
   box(); 
   
   if(istrata == 1) {
    axis(side = 3)
    mtext(side = 3, paste0(isample, ' Boat'), line = 2, cex = 1.5)
   }
   if(istrata == 6) axis(side = 1)
   
   if(isample == 2) {
    legend('topright', paste0(c(5,10,15,20,30,60)[istrata], ' Strata'), 
           bty = 'n', cex = 2)
    if(istrata == 1) mtext(side = 3, sci_names[ispp], outer = T, line = 4, 
                           font = 3, cex = 1.5)
   }
   
   for(iyear in 1:11){
    temp_response = 100*(STRS_sim_cv[iyear, ispp, istrata, isample, ] - 
                      STRS_true_cv_array[iyear, ispp, istrata, isample]) / STRS_true_cv_array[iyear, ispp, istrata, isample]
    temp_density = density(temp_response)
    temp_density$y = temp_density$y / max(temp_density$y)
    lines(temp_density, col = heat.colors(11)[iyear])
   }
   
   temp_response = sweep(x = STRS_sim_cv[, ispp, istrata, isample, ],  
                         STATS = STRS_true_cv_array[, ispp, istrata, isample],
                         MARGIN = 1, FUN = '-') 
   temp_response = 100*sweep(x = temp_response,  
                         STATS = STRS_true_cv_array[, ispp, istrata, isample],
                         MARGIN = 1, FUN = '/') 
   temp_response = apply(temp_response, MARGIN = 1, FUN = median)
   boxplot(temp_response, horizontal = T, add = T, axes = F, at = 0.25, xpd = NA)
   
  }
 }
 mtext(side = 1, 'Relative Percent Bias', outer = T, line = 3)
}

dev.off()

######################
pdf(file = paste0(PP_dir, 'Bias_Est_distribution.pdf'), width = 8, height = 10, onefile = T)
par(mfrow = c(6,3), mar = c(0,0,0,0), oma = c(5,1,6,1))
for(ispp in 1:ns){
 for(istrata in 1:6){
  for(isample in 1:3){
   plot(1, type = 'n', ylim = c(0,1), xlim = c(-75,100), ann = F, axes = F)
   abline(v = 0, lty = 'dashed')
   box(); 
   
   if(istrata == 1) {
    axis(side = 3)
    mtext(side = 3, paste0(isample, ' Boat'), line = 2, cex = 1.5)
   }
   if(istrata == 6) axis(side = 1)
   
   if(isample == 2) {
    legend('topright', paste0(c(5,10,15,20,30,60)[istrata], ' Strata'), 
           bty = 'n', cex = 2)
    if(istrata == 1) mtext(side = 3, sci_names[ispp], outer = T, line = 4, 
                           font = 3, cex = 1.5)
   }
   
   for(iyear in 1:11){
    temp_response = 100*(STRS_sim_mean[iyear, ispp, istrata, isample, ] - 
                      true_mean[iyear, ispp]) / true_mean[iyear, ispp]
    temp_density = density(temp_response)
    temp_density$y = temp_density$y / max(temp_density$y)
    lines(temp_density, col = heat.colors(11)[iyear])
   }
   
   temp_response = sweep(x = STRS_sim_mean[, ispp, istrata, isample, ],  
                         STATS = true_mean[, ispp],
                         MARGIN = 1, FUN = '-') 
   temp_response = 100*sweep(x = temp_response,  
                         STATS = true_mean[, ispp],
                         MARGIN = 1, FUN = '/') 
   temp_response = apply(temp_response, MARGIN = 1, FUN = median)
   boxplot(temp_response, horizontal = T, add = T, axes = F, at = 0.25, xpd = NA)
   
  }
 }
 mtext(side = 1, 'Relative Percent Bias', outer = T, line = 3)
}
dev.off()
