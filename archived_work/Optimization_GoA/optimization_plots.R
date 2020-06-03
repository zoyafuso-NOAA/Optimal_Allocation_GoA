###########################
## Tradeoff Plots for MS Optimization
###########################
setwd("/Users/zackoyafuso/Documents/GitHub/MS_OM_GoA/")

###########################
## Import Optimization results and species data
###########################
species_df = read.csv('/Users/zackoyafuso/Desktop/AK_BTS/data-raw/species.csv')
load('Optimization_GoA/optimization_results.RData')
species_names = c(as.character(species_df[species_df$SPECIES_CODE%in%species,'COMMON_NAME']),
                  'All Species')
species_names[2] = "rougheye and blackspotted rockfish"


library(rgdal); library(maptools); library(rgeos); library(raster)

test = readOGR('/Users/zackoyafuso/Desktop/survey_grids/goagrid_nolandsandman.shp')
test_dissolved = unionSpatialPolygons(SpP = test, IDs = test$STRATUM)


page = 1
{
 for(ispp in c(15, 1:14)){
  if(ispp %in% c(15,4,8,12) ){
   tiff(filename = paste0('Optimization_GoA/prelim_figures/tradeoff_page', 
                          page, '.tiff'), 
        width = 9, height = 8,
        units = 'in', compression = 'lzw', res = 200)
   page = page + 1
   par(mfrow = c(4,4), oma = c(0,0,0,0))
  }
  
  par(mar = c(5,4,1,1))
  plot(tot_mean ~ tot_var, data = res_df,
       subset = spp_scen == ispp, main = species_names[ispp], cex.main = 0.9,
       type = 'n', las = 1, pch = 16, lwd = 2, cex = 0.5, col = 'grey',
       xlab = 'Total Variance Criterion', ylab = "Total Mean Criterion")
  
  for(i in c(40, 30, 20)){
   temp_df = subset(res_df, n == i & spp_scen == ispp)
   temp_mat = matrix( res_mat[res_df$n == i & res_df$spp_scen == ispp,], 
                      ncol = nstrata )
   
   max_mean_idx = which.max(temp_df$tot_mean)
   max_var_idx = which.max(temp_df$tot_var)
   
   range_mean = diff(range(temp_df$tot_mean))
   range_var = diff(range(temp_df$tot_var))
   
   dist_mean = (temp_df$tot_mean - max(temp_df$tot_mean)) / range_mean
   dist_var =  (temp_df$tot_var - max(temp_df$tot_var)) / range_var
   
   comp_idx = which.min(sqrt(dist_mean^2 + dist_var^2))
   
   sol_idx = c(max_mean_idx, comp_idx, max_var_idx)
   
   points(tot_mean ~ tot_var, data = temp_df, 
          pch = 16, cex = 0.5, col = 'gray')
   points(tot_mean ~ tot_var, data = temp_df[sol_idx,], 
          pch = 16, cex = 1,
          col = c('black', 'darkblue', 'red'))   
   
   with(temp_df, 
        text(min(tot_var), min(tot_mean), paste('n =', i) ))
  }
  
  par(mar = c(0,0,0,0))
  
  for(i in c(20, 30, 40)){
   
   plot(1, type = 'n', axes = F, ann = F, asp = 1,
        xlim = extent(test_dissolved)[1:2],
        ylim = extent(test_dissolved)[3:4]+c(0,1e6))
   
   temp_df = subset(res_df, n == i & spp_scen == ispp)
   temp_mat = matrix( res_mat[res_df$n == i & res_df$spp_scen == ispp,], 
                      ncol = nstrata )
   
   max_mean_idx = which.max(temp_df$tot_mean)
   max_var_idx = which.max(temp_df$tot_var)
   
   range_mean = diff(range(temp_df$tot_mean))
   range_var = diff(range(temp_df$tot_var))
   
   dist_mean = (temp_df$tot_mean - max(temp_df$tot_mean)) / range_mean
   dist_var =  (temp_df$tot_var - max(temp_df$tot_var)) / range_var
   
   comp_idx = which.min(sqrt(dist_mean^2 + dist_var^2))
   
   sol_idx = c(max_mean_idx, comp_idx, max_var_idx)
   
   plot(test_dissolved, ann = F, axes = F, asp = 1, col = 'grey',
        add = T, border = F ) 
   #box()
   s = 3
   plot(test_dissolved[strata[temp_mat[sol_idx[s],] == 1],], add = T, 
        col = c('black', 'darkblue', 'red')[s], 
        border = c('black', 'darkblue', 'red')[s])
   
   test_dissolved2 = raster::shift(x = test_dissolved, dy = 5e5)
   
   plot(test_dissolved2, ann = F, axes = F, asp = 1, 
        col = 'grey', border = F, add = T); 
   
   s = 2
   plot(test_dissolved2[temp_mat[sol_idx[s],] == 1,], 
        add = T, col = c('black', 'darkblue', 'red')[s], 
        border = c('black', 'darkblue', 'red')[s])
   
   test_dissolved3 = raster::shift(x = test_dissolved2, dy = 5e5)
   
   plot(test_dissolved3, ann = F, axes = F, asp = 1, 
        col = 'grey', border = F, add = T); 
   
   s = 1
   plot(test_dissolved3[temp_mat[sol_idx[s],] == 1,], 
        add = T, col = c('black', 'darkblue', 'red')[s], 
        border = c('black', 'darkblue', 'red')[s])
   
   legend('top', legend = paste('n =', i), bty = 'n')
  }
  
  if(ispp %in% c(3,7,11,14) ) dev.off()
 }
}


# temp = test_dissolved
# 
# all_sol = array(data = 0, dim = c(ns+1,3,nstrata))
# 
# 
# for(ispp in 1:(ns+1)){
#    for(i in 5:49 ){
#       temp_df = subset(res_df, n == i & spp_scen == ispp)
#       temp_mat = matrix( res_mat[res_df$n == i & res_df$spp_scen == ispp,], 
#                          ncol = nstrata )
#       
#       max_mean_idx = which.max(temp_df$tot_mean)
#       max_var_idx = which.max(temp_df$tot_var)
#       
#       range_mean = diff(range(temp_df$tot_mean))
#       range_var = diff(range(temp_df$tot_var))
#       
#       dist_mean = (temp_df$tot_mean - max(temp_df$tot_mean)) / range_mean
#       dist_var =  (temp_df$tot_var - max(temp_df$tot_var)) / range_var
#       
#       comp_idx = which.min(sqrt(dist_mean^2 + dist_var^2))
#       
#       sol_idx = c(max_mean_idx, comp_idx, max_var_idx)
#       
#       all_sol[ispp,,] = all_sol[ispp,,] + temp_mat[sol_idx,]
#    }
#    
# }
# 
# all_sol[1,,]
# 
# plot(temp, col = grey.colors(max(all_sol[1,2,]))[all_sol[1,2,]], border = NA)

dev.off( )

>>>>>>> 31eac8fc54aae7c8725e74c6d672a64c8c49e92b
