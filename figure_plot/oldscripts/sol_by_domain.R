rm(list = ls())
library(raster); library(RColorBrewer); library(SamplingStrata)

which_machine = c('Zack_MAC' = 1)
result_wd = c('/Users/zackoyafuso/Documents/GitHub/MS_OM_GoA/Optimum_Allocation/', 
              'C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/Optimum_Allocation/')[which_machine]

setwd(result_wd)
load('../Extrapolation_depths.RData')
VAST_model = '6g'
load(paste0('model_', VAST_model, '/optimization.RData'))

which_cv = which(settings$cv[1:length(strata_list)] == 0.15)
settings = settings[which_cv,]
strata_list = strata_list[which_cv]
res_df = as.data.frame(res_df)
res_df = res_df[,paste0('V', which_cv + 2)] 
colnames(res_df) = paste0('sol_', 1:ncol(res_df))

sample_sizes = sapply(strata_list, FUN = function(x) sum(x$Allocation))

winner = which.min(sample_sizes)
ns = 15
goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')], 
                             data = data.frame(X1=res_df[,winner]) )

{tiff(paste0('model_', VAST_model, '/solution_map.tiff'), 
      res = 200, width = 8, height = 6, units = 'in', compression = 'lzw')
  par(mfrow = c(1,1), oma = rep(1,4), family = 'serif' )
  
  goa_ras = raster(goa, resolution = 5)
  goa_ras = rasterize(x = goa, y = goa_ras, field = 'X1')
  
  temp_df = strata_list[[winner]]
  nstrata = nrow(temp_df)
  
  par(mar = c(0,0,0,0))
  image(goa_ras, asp = 1, axes = F,
        col = brewer.pal(n = nstrata, name = 'Paired'))
  
  xrange = diff(par()$usr[1:2])
  yrange = diff(par()$usr[3:4])
  text(x = par()$usr[1]+xrange*0.725,
       y = par()$usr[3]+yrange*0.475,
       paste0('Optimal Sample Size: ', 
              sum(temp_df$Allocation), '\n',
              'CV constraint: ', settings$cv[winner]*100, '%'),
       cex = 1.5)
  
  dev.off()
}

frame$X2 = round(frame$X2)
temp_df = strata_list[[winner]]
nstrata = nrow(temp_df)
depth_cuts = unique(round(sort(c(temp_df$Lower_X1, temp_df$Upper_X1))) )
lon_cuts = unique(round(sort(c(temp_df$Lower_X2, temp_df$Upper_X2))))

matrix_space = matrix(data = nstrata + 1,
                      nrow = length(depth_cuts)-1, 
                      ncol = length(lon_cuts)-1,
                      dimnames = list(
                        paste0(depth_cuts[1:(length(depth_cuts)-1)], '-',
                               depth_cuts[2:length(depth_cuts)]),
                        paste0(lon_cuts[1:(length(lon_cuts)-1)], '-',
                               lon_cuts[2:length(lon_cuts)]) ) )

for(istrata in nstrata:1){
  depth_bounds = round(unlist(temp_df[istrata,c('Lower_X1','Upper_X1')]))
  col_idx = c(which(depth_cuts[1:(length(depth_cuts)-1)] == depth_bounds[1]),
              which(depth_cuts[2:length(depth_cuts)] == depth_bounds[2] ) )
  
  lon_bounds = round(unlist(temp_df[istrata,c('Lower_X2','Upper_X2')]))
  row_idx = c(which(lon_cuts[1:(length(lon_cuts)-1)] == lon_bounds[1]),
              which(lon_cuts[2:length(lon_cuts)] == lon_bounds[2] ))
  
  for(j in seq(row_idx[1],row_idx[2])){
    for(i in seq(col_idx[1],col_idx[2]) ){
      matrix_space[i,j] = ifelse(matrix_space[i,j] > istrata, 
                                 istrata, 
                                 matrix_space[i,j])
    }
  }
}

for(j in 1:(length(lon_cuts)-1) ){
  for(i in 1:(length(depth_cuts)-1)){
    matrix_space[i,j] = ifelse(any(frame$X1 >= depth_cuts[i] & 
                                     frame$X1 <= depth_cuts[i+1] &
                                     frame$X2 >= lon_cuts[j] & 
                                     frame$X2 <= lon_cuts[j+1]), 
                               matrix_space[i,j], 
                               nstrata+1)
  }
}

{tiff(paste0('model_', VAST_model, '/solution_space.tiff'), 
      res = 200, width = 5.85, height = 3.65, units = 'in', 
      compression = 'lzw')
  par(mar = c(4,4,1,1), oma = c(0,0.5,1,0.5))
  sol_ras = raster(matrix_space, xmn = 0, xmx = ncol(matrix_space),
                   ymn = 0, ymx = nrow(matrix_space) )
  image(sol_ras, axes = F, asp = 1, xlab = 'Eastings (km)', ylab = 'Depth (m)', 
        col = c(brewer.pal(n = nstrata, 'Paired'), 'black'))
  axis(side = 1, at = 0:(sol_ras@extent[2]) ,
       labels = lon_cuts, las = 2, cex.axis = 0.6)
  axis(side = 2, at = (sol_ras@extent[4]):0 , 
       labels = depth_cuts, las = 1)
  
  abline(h = 0:length(depth_cuts), v = 0:length(lon_cuts))
  legend(x = 0.25, y = 11.25, legend = 1:nstrata, horiz = T, xpd = NA, 
         cex = 0.75, bty = 'n', fill = brewer.pal(n = nstrata, 'Paired'))
  dev.off()
}

