#############################
## Show Simulation Metrics for Simulated Survey Strata,
## Simple Random Sampling, and Stratified Random Sampling
#############################
rm(list = ls())

############################
## Set up directories
#############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[2]
optimization_type = c('_spatial', '_spatiotemporal')[2]
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


plot_set = data.frame(spp = sci_names,
                      true_ymax = c(0.15, 0.30, 0.12, 
                                    0.15, 0.15, 0.08, 
                                    0.25, 0.25, 0.35, 
                                    0.35, 0.40, 0.35, 
                                    0.35, 0.40, 0.40),
                      rrmse_ymax =c(0.45, 0.50, 0.40, 
                                    0.50, 0.30, 0.50, 
                                    0.45, 0.50, 0.30, 
                                    0.66, 0.60, 0.60, 
                                    0.50, 0.80, 0.60))
par(mfcol = c(2,3), mar = c(1,4,2,0), oma = c(2,2,1,1))
# for(spp in c(1,3,14, 4,10,12)){
for(spp in 1:ns){
  # if(spp %in% c(1,4)){
  #   png(filename = paste0(PP_dir, 'Sim_Metrics_',
  #                         c('1' = 'NoChange', '4' = 'RRMSEChange')[paste(spp)],
  #                         '.png'),
  #       width = 9, height = 3.75, units = 'in', res = 500)
  #   par(mfcol = c(2,3), mar = c(1,4,2,0), oma = c(2,2,1,1))
  # }
  
  for(imetric in c('true', 'rrmse')){
    #survey_result = get(paste0('survey_', imetric, '_cv_array'))
    #SRS_result = get(paste0('SRS_', imetric, '_cv_array'))
    STRS_result = get(paste0('STRS_', imetric, '_cv_array'))
    
    #Empty plot for True CV
    plot(1, type = 'n', xlim = c(2,22), 
         ylim = c(0, plot_set[spp, paste0(imetric, '_ymax')] ), 
         ann = F, axes = F)
    legend('bottom', legend = paste('n =', c(820, 550, 280)), 
           fill = c('white', 'blue', 'red'), ncol = 3, bty = 'n')
    box()
    
    if(imetric == 'true') mtext(side = 3, sci_names[spp], font = 3)
    if(spp%in%c(1,4)) mtext(side = 2, 
                            ifelse(imetric == 'true',
                                   'True CV\nAcross Years',
                                   'RRMSE of\nCV Across Years'), 
                            line = 3)
    
    axis(side = 1, at = c(4,8,12,16,20), labels = NA)
    axis(side = 1, at = c(4,8,12,16,20), 
         labels = paste(c(5,10,20,40,60), '\nstrata'),
         line = 0.5, tick=F)
    axis(side = 2, las = 1)
    
    #Simulated Survey Strata
    #boxplot( survey_result[,spp], add = T, width = 1, axes = F, at = 1)
    
    #Simulated Simple Random Sampling, 350 and 550 samples
    # boxplot( SRS_result[,spp,c('size_550', 'size_350')], add = T, 
    #          axes = F, at = 3:4)
    
    #Simulated Stratified Random Sampling, 5 strata
    which_strata = which(settings$nstrata == 5)
    row_idx = which_strata[which.min(abs(settings[which_strata,'n']-800) )]
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1,
             axes = F, at = 3)
    
    row_idx = which_strata[which.min(abs(settings[which_strata,'n']-550) )]
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1,
             axes = F, at = 4, col = 'blue')
    
    row_idx = which_strata[which.min(abs(settings[which_strata,'n']-280) )]
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1,
             axes = F, at = 5, col = 'red')
    
    #Simulated Stratified Random Sampling, 10 strata
    which_strata = which(settings$nstrata == 10)
    row_idx = which_strata[which.min(abs(settings[which_strata,'n']-800) )]
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1,
             axes = F, at = 7)
    row_idx = which_strata[which.min(abs(settings[which_strata,'n']-550) )]
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1, col = 'blue',
             axes = F, at = 8)
    row_idx = which_strata[which.min(abs(settings[which_strata,'n']-280) )]
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1,
             axes = F, at = 9, col = 'red')
    
    #Simulated Stratified Random Sampling, 20 strata
    which_strata = which(settings$nstrata == 20)
    row_idx = which_strata[which.min(abs(settings[which_strata,'n']-800) )]
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1,
             axes = F, at = 11)
    row_idx = which_strata[which.min(abs(settings[which_strata,'n']-550) )]
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1,
             axes = F, at = 12, col = 'blue')
    row_idx = which_strata[which.min(abs(settings[which_strata,'n']-280) )]
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1,
             axes = F, at = 13, col = 'red')
    
    #Simulated Stratified Random Sampling, 40 strata
    which_strata = which(settings$nstrata == 40)
    row_idx = which_strata[which.min(abs(settings[which_strata,'n']-800) )]
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1,
             axes = F, at = 15)
    row_idx = which_strata[which.min(abs(settings[which_strata,'n']-550) )]
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1,
             axes = F, at = 16, col = 'blue')
    row_idx = which_strata[which.min(abs(settings[which_strata,'n']-280) )]
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1,
             axes = F, at = 17, col = 'red')
    
    #Simulated Stratified Random Sampling, 60 strata
    which_strata = which(settings$nstrata == 60)
    row_idx = which_strata[which.min(abs(settings[which_strata,'n']-800) )]
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1,
             axes = F, at = 19)
    
    row_idx = which_strata[which.min(abs(settings[which_strata,'n']-550) )]
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1, col = 'blue',
             axes = F, at = 20)
    row_idx = which_strata[which.min(abs(settings[which_strata,'n']-280) )]
    boxplot( STRS_result[,spp,row_idx], add = T, width = 1, col = 'red',
             axes = F, at = 21)
  }
  
  # if(spp%in% c(12,14)) dev.off() 
}

