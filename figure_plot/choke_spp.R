#############################
## How does the spatiotemporal CV relate to true CV
#############################
rm(list = ls())

############################
## Set up directories
#############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[2]
optimization_type = c('_spatial', '_spatiotemporal')[1]
VAST_model = '6g'

output_wd = paste0(c('/Users/zackoyafuso/Documents/', 
                     'C:/Users/Zack Oyafuso/Documents/',
                     'C:/Users/zack.oyafuso/Work/', 
                     'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                   "GitHub/MS_OM_GoA/Optimum_Allocation/model_", VAST_model,
                   optimization_type)

paper_dir = paste0(c('/Users/zackoyafuso/', 
                     'C:/Users/Zack Oyafuso/')[which_machine],
                   'Google Drive/MS_Optimizations/figure_plot/')
PP_dir = paste0(c('/Users/zackoyafuso/', 
                  'C:/Users/Zack Oyafuso/')[which_machine],
                'Google Drive/MS_Optimizations/powerpoint_plot/')

load(paste0(output_wd, '/Stratified_RS_Simulation_Results.RData'))
load(paste0(output_wd, '/optimization_results.RData'))
rockfish_cod_idx = c(2,3,11:15)
flatfish_idx = (1:ns)[-rockfish_cod_idx]

settings$id = 1:nrow(settings)

{png(filename = paste0(PP_dir, 'choke_spp', optimization_type,'.png'), 
     width = 10, height = 5, units = 'in', res =500)
par(mfrow = c(1,2), mar = c(0,0,0,0), oma = c(5,5,1,1))
for(istrata in c(5)){
 sub_settings = subset(settings, nstrata == istrata)
 
 for(idx in list(flatfish_idx, rockfish_cod_idx)){
  mean_truecv = apply( STRS_true_cv_array[,idx,sub_settings$id], 
                       MARGIN = 3:2, mean)
  
  if(optimization_type == '_spatial') xlim_ = c(0.1,0.19)
  if(optimization_type == '_spatiotemporal') xlim_ = c(0.15,0.31)
  
  matplot(x = sub_settings$cv[(round((sub_settings$cv * 1000)) %% 10) == 0], 
          y = mean_truecv[(round((sub_settings$cv * 1000)) %% 10) == 0,], 
          lty = 1, type = 'b',  pch = paste0(1:length(idx)), col = 'grey',
          xlim = xlim_, ylim = c(0, 0.4), las = 1,
          axes = F, ann = F)
  box()
  if(identical(idx, flatfish_idx)) axis(side = 2, las = 1)
  if(istrata == 5) axis(side = 1)
  abline(a = 0, b = 1)
  legend('topleft', legend = sci_names[idx], ncol = 2, pch = paste0(1:length(idx)),
         col = 'black', bty = 'n', cex = 0.75)
 }
}

mtext(side = 1, 'Spatiotemporal Upper CV Constraint', outer = T, line = 3)
mtext(side = 2, 'True CV Average across Years', outer = T, line = 3)
dev.off()}
