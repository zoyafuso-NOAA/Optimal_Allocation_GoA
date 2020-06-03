#################################
## What do the stratifications look like?
#############################
rm(list = ls())

library(sp); library(raster); library(RColorBrewer)

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

load(paste0(dirname(dirname(output_wd)), '/Extrapolation_depths.RData' ))
load(paste0(output_wd, '/Stratified_RS_Simulation_Results.RData'))
load(paste0(output_wd, '/optimization_results.RData'))
load(paste0(output_wd, '/optimization_data_model_', VAST_model, '.RData'))

settings$id = 1:nrow(settings)
samples = c(820, 550, 280)
stratas = c(5,10,15)

yrange = diff(range(Extrapolation_depths[,c('N_km')]))

plot_random_sample = F

{
  png(file = paste0(PP_dir, 'sol_by_boat',
                   ifelse(plot_random_sample == T, '_withsamples', ''),
                   '_',optimization_type,
                   '.png'),
    width = 240, height = 90, units = 'mm', res = 1000)
par(mfrow = c(1,3), mar = c(0,0,3,0))
for(istrata in stratas){
 plot(1, type = 'n', axes = F, ann = F,
      xlim = range(Extrapolation_depths[,c('E_km')]),
      ylim = c(min(Extrapolation_depths[,c('N_km')])-1.5*yrange,
               max(Extrapolation_depths[,c('N_km')]))
 )
 mtext(side = 3, paste(istrata, 'Strata'))
 offset = 0
 
 for(isample in samples){
  sub_settings = subset(settings, nstrata == istrata)
  
  isol = sub_settings$id[which.min(abs(sub_settings$n - isample))]
  
  #Sample based on the stratification allocations
  sample_vec = c()
  for(i in 1:istrata ){
   available_cells = which(res_df[,isol+1] == i)
   
   if(length(available_cells) > 0){
    sample_cells = sample(x = available_cells, 
                          size = strata_list[[isol]]$Allocation[i], 
                          replace = F)
    sample_vec = c(sample_vec, sample_cells)
   }

  }
  
  #Organize sample set and total number of samples
  sample_vec = sort(sample_vec)
  sample_pts = Extrapolation_depths[sample_vec, c('E_km', 'N_km')]
  sample_pts[,'N_km'] = sample_pts[,'N_km'] + offset
  
  goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')],
                               data = data.frame(stratum = res_df[,isol+1]) )
  goa_ras = raster(goa, resolution = 5)
  goa_ras =rasterize(x = goa, y = goa_ras, field = 'stratum')
  
  goa_ras = raster::shift(goa_ras, dy = offset)
  offset = offset - yrange*.75
  
  image(goa_ras, asp = 1, axes = F, add = T,
        col = brewer.pal(n = istrata, name = 'Paired'))
  xrange = diff(goa_ras@extent[1:2])
  yrange = diff(goa_ras@extent[3:4])
  text(x = goa_ras@extent[1]+xrange*0.70,
       y = goa_ras@extent[3]+yrange*0.60,
       paste0('n = ', isample, '\n'),
       cex = 1.5)
  
  if(plot_random_sample) points(sample_pts, pch = 16, cex = 0.2)
  offset = offset + 1
 }
 
}

dev.off()
}

# frame$X2 = round(frame$X2)
# temp_df = strata_list[[isol]]
# depth_cuts = unique(round(sort(c(temp_df$Lower_X1, temp_df$Upper_X1))) )
# lon_cuts = unique(round(sort(c(temp_df$Lower_X2, temp_df$Upper_X2))))
# 
# matrix_space = matrix(data = istrata + 1,
#                       nrow = length(depth_cuts)-1, 
#                       ncol = length(lon_cuts)-1,
#                       dimnames = list(
#                        paste0(depth_cuts[1:(length(depth_cuts)-1)], '-',
#                               depth_cuts[2:length(depth_cuts)]),
#                        paste0(lon_cuts[1:(length(lon_cuts)-1)], '-',
#                               lon_cuts[2:length(lon_cuts)]) ) )
# 
# for(istrata in nstrata:1){
#  depth_bounds = round(unlist(temp_df[istrata,c('Lower_X1','Upper_X1')]))
#  col_idx = c(which(depth_cuts[1:(length(depth_cuts)-1)] == depth_bounds[1]),
#              which(depth_cuts[2:length(depth_cuts)] == depth_bounds[2] ) )
#  
#  lon_bounds = round(unlist(temp_df[istrata,c('Lower_X2','Upper_X2')]))
#  row_idx = c(which(lon_cuts[1:(length(lon_cuts)-1)] == lon_bounds[1]),
#              which(lon_cuts[2:length(lon_cuts)] == lon_bounds[2] ))
#  
#  for(j in seq(row_idx[1],row_idx[2])){
#   for(i in seq(col_idx[1],col_idx[2]) ){
#    matrix_space[i,j] = ifelse(matrix_space[i,j] > istrata, 
#                               istrata, 
#                               matrix_space[i,j])
#   }
#  }
# }
# 
# for(j in 1:(length(lon_cuts)-1) ){
#  for(i in 1:(length(depth_cuts)-1)){
#   matrix_space[i,j] = ifelse(any(frame$X1 >= depth_cuts[i] & 
#                                   frame$X1 <= depth_cuts[i+1] &
#                                   frame$X2 >= lon_cuts[j] & 
#                                   frame$X2 <= lon_cuts[j+1]), 
#                              matrix_space[i,j], 
#                              nstrata+1)
#  }
# }
