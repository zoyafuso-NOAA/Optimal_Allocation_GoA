######################################
## Plot Performance Metrics
######################################
rm(list = ls())
library(sp); library(raster)

######################################
## Set up directories
######################################
which_machine = c('Zack_Mac' = 1)[1]
output_wd = paste0(c('/Users/zackoyafuso/Documents/')[which_machine],
                   'GitHub/MS_OM_GoA/Optimum_Allocation_Experiment/')

######################################
## Import OM
######################################
load(paste0(output_wd, '/OM.RData'))
load(paste0(output_wd, '/optimization.RData'))
#load(paste0(output_wd, '/survey_simulation.RData'))

par(mfrow = c(1,2), mar = c(1,1,1,1))
for(iscen in 1:2){
 #Plot Solution
 stratano = res_df[,1+iscen]
 domain = SpatialPointsDataFrame(coords = OM[[iscen]]$loc_xy,
                                 data = data.frame(Str_no = stratano) )
 domain_ras = raster(domain, resolution = 1/100)
 goa_ras =rasterize(x = domain, y = domain_ras, field = 'Str_no')
 plot(goa_ras, col = terrain.colors(6)[-6], axes = F, 
      main = paste0('Population ', c('A', 'B')[iscen]))
 legend('topleft', legend = paste0('Stratum ', 1:5, ': ', strata_list[[iscen]]$Allocation, ' samples'), fill = terrain.colors(6)[-6], bty = 'n')
}

