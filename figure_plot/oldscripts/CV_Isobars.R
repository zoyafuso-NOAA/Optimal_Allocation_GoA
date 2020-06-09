#########################
## CV Isobars
#########################
rm(list = ls())

############################
## Set up directories
#############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[2]
VAST_model = "6g"
optimization_type = c('_spatial', '_spatiotemporal')[1]

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

load(paste0(output_wd, '/optimization_results.RData'))
settings$cv = as.integer(settings$cv * 1000)

#############################
## CV Isobars
#############################

for(figure_type in c('paper_dir', 'PP_dir')){
   # png(filename = paste0(get(figure_type), 'CV_Isobars.png'), res = 500,
   #     units = c('paper_dir' = 'in', 'PP_dir' = 'in')[figure_type],
   #     width = c('paper_dir' = 6, 'PP_dir' = 6)[figure_type], 
   #     height = c('paper_dir' = 5, 'PP_dir' = 5)[figure_type])
   
   par(mar = c(5,5,1,1), mfrow = c(1,1))
   plot(n ~ nstrata, data = settings, type = 'n', las = 1, 
        xlim = c(5,max(settings$nstrata) + 15), axes = F, ann = F )
   axis(side = 1, at = c(seq(5,25,5), seq(30,60,10)))
   axis(side = 2, las =1 )
   mtext(side = 1, 'Number of Strata', line = 3); mtext(side = 2, 'Total Sample Size', line = 3)
   box()
   abline(h = c(280, 550, 800), col='grey', lwd = 2, lty = 'dotted')
   
   for(icv in c(10:20)*10 ) {

      lines(n ~ nstrata, data = settings, subset = cv == icv)
      points(n ~ nstrata, data = settings, subset = cv == icv, 
             pch = 16, cex = 0.75)
      text(x = max(settings$nstrata[settings$cv == icv]),
           y = settings$n[settings$cv == icv][which.max(settings$nstrata[settings$cv == icv])],
           paste0(icv*.1, '% CV'),
           pos = 4, cex = 0.75)
   }

   text(x = 72.5, y = c(300, 570, 820), 
        c('1 Boat', '2 Boats', '3 Boats'), col = 'grey')
   # dev.off()
}
