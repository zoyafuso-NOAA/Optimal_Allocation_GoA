#########################
## CV Isobars
#########################
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
settings$id = 1:nrow(settings)

samples = c(820, 550, 280)
stratas = c(5,10,15)

sort_spp = names(sort(apply(STRS_true_cv_array, MARGIN = c(2), mean)))



{png(filename = paste0(PP_dir, 'CV_by_Spp.png'), 
    width = 10, height = 6, units = 'in', res = 500)
par(mar = c(3,7,3,1), oma = c(1,1,0,0), mfcol = c(1,3))
for(istrata in stratas){
 plot(1, type = 'n', ylim = c(0,57), xlim = c(0, 0.45), 
      axes = F, ann = F)
 abline(#v = seq(0.1, 0.4, 0.1), 
        h = seq(from = 3, by = 4, length = ns),
        col = 'lightgrey', lty = 'dashed')
 #legend('topleft', legend = paste(istrata, 'Strata'), bty = 'n', cex = 2)
 axis(side = 1, las = 1); axis(side = 3, las = 1)
 axis(side = 2, at = seq(from = 1, by = 4, length = ns), labels = NA)
 axis(side = 2, at = seq(from = 1, by = 4, length = ns), 
      labels = gsub(sort_spp, pattern = ' ', replacement = '\n'), las = 1,
      tick = F, line = 1, cex.axis = 0.9)
 box()
 legend('bottomright', legend = paste('n =', samples), fill = c('white', 'blue', 'red'),
        title = paste(istrata, 'Strata'), cex = 1.5)
 
 offset = 0
 for(spp in 1:ns){
  for(isample in samples){
   sub_settings = subset(settings, nstrata == istrata)
   
   isol = sub_settings$id[which.min(abs(sub_settings$n - isample))]
   
   boxplot(STRS_true_cv_array[,sort_spp[spp],isol], add = T, axes = F, horizontal = T,
           at = c('820' = 0,
                  '550' = 1,
                  '280' = 2)[paste(isample)] + 3*(spp-1) + offset,
           col = c('820' = 'white',
                   '550' = 'blue',
                   '280' = 'red')[paste(isample)])
  }
  offset = offset + 1
 }
}
mtext(side = 1, 'True CV', line = 0, outer = T, cex = 1)
dev.off()}
