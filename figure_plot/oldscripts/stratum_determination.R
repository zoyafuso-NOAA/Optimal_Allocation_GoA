###########################
## Parameter Tuning
## Reference: http://cran.r-project.org/web/packages/segmented/index.html
############################
library(segmented)
rm(list = ls())
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2)[1]
result_wd = paste0(c('/Users/zackoyafuso/Documents/', 
                     'C:/Users/Zack Oyafuso/Documents/')[which_machine], 
                   'GitHub/MS_OM_GoA/Optimum_Allocation/')

setwd(result_wd)
load('../Extrapolation_depths.RData')
VAST_model = '6g'

load(paste0('model_', VAST_model, '/optimization_ST_master.RData'))

temp_df = cbind(settings, n = sapply(strata_list, FUN = function(x) sum(x$Allocation)) )

########################
## Parameter choice
########################
perm_check = lm(n ~ as.factor(mut_change) 
                + as.factor(elitism_rate) 
                + nstata 
                + as.factor(cv), 
                data = temp_df, subset = nstata %in% c(5,7,10))

summary(perm_check)

########################
## Piecewise Regression
########################
{
 tiff('../figure_plot/strata_determination.tiff', res = 200,
      width = 90, height = 120, units = 'mm', compression = 'lzw')
 par(mfrow = c(1,1), mar = c(0,0,0,0), oma = c(3,5,1,1), family = 'serif')
 plot(n ~ nstata, data = temp_df, ann = F, axes = F, type = 'n',
      ylim = c(500,1100))
 
 for(CV in c(0.15, 0.20)){
  
  sub_df = subset(temp_df, 
                  subset = (cv==CV) & (mut_change==0.1) & (elitism_rate==0.1))
  my.lm = lm(n ~ nstata, data = sub_df)
  
  my.seg <- segmented(my.lm, 
                      seg.Z = ~ nstata, 
                      psi = 13)
  
  print(summary(my.seg))
  print(slope(my.seg))
  
  pred_stata = seq(min(temp_df$nstata), max(temp_df$nstata), length = 100)
  pred_val = predict(my.seg, newdata = data.frame(nstata = pred_stata))
  
  points(n ~ nstata, data = sub_df, pch = 16, cex = 0.5)
  lines(pred_stata, pred_val, lwd = 2, col = 'darkgrey')
  
  points(x = my.seg$psi[2],  
         y = predict(my.seg, newdata = data.frame(nstata = my.seg$psi[2])),
         pch = 15, cex = 1, col = 'darkgrey' )
  
  axis(side = 1, at = seq(5,20,5))
  axis(side = 2, las = 1)
  legend(x = 8.5, y = c('0.15'=1050, '0.2' = 650)[paste(CV)],
         paste0(CV*100, '% CV Constraint'), 
         bty = 'n', cex = 1)
  
 }
 
 box()
 mtext(side = 1, 'Number of Strata', outer = T, line = 2)
 mtext(side = 2, 'Total Sample Size', outer = T, line = 3)
 dev.off()
}
