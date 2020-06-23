##########################
## Simulation Results
## Simple Random Sampling
##########################

rm(list = ls())

###############################
## Set up directories
###############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[1]
VAST_model = "6g"
github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/Optimum_Allocation/model_', VAST_model, '/')

load(paste0(github_dir, "Simple_RS_Simulation_Results.RData" ))

spp = 15

par(mfrow = c(2,1), mar = c(0,5,0,0), oma = c(3,0,2,1))
for(spp in sci_names){
   plot(1, type = 'n', xlim = c(30, 100), ylim = c(0,max(true_cv_array[,spp,])), 
        axes = F, ylab = 'True CV across years', xlab = '')
   axis(side = 2, las = 1)
   boxplot(true_cv_array[,spp,], at = nsamples/10, labels = NULL, add = T, axes = F)
   box()
   
   plot(1, type = 'n', xlim = c(30, 100), ylim = c(0,max(rrmse_cv_array[,spp,])), 
        axes = F, ylab = 'RRMSE CV across years', xlab = 'Simple Random Sample Size')
   axis(side = 2, las = 1); axis(side = 1, at = seq(30,100,10), labels = seq(300,1000,100))
   boxplot(rrmse_cv_array[,spp,], at = nsamples/10, labels = NULL, add = T, axes = F)
   box()
   mtext(side = 3, spp, outer = T, line = 0.5, cex = 1.5, font = 3)
}


