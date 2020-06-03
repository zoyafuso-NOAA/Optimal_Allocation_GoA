#######################################
## Plot Best solution
######################################
rm(list = ls())

######################################
## Import libraries
######################################
library(raster); library(RColorBrewer); library(SamplingStrata)

######################################
##Set up directories
######################################
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2)[1]
result_wd = c('/Users/zackoyafuso/Documents/GitHub/MS_OM_GoA/Optimum_Allocation/', 
              'C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/Optimum_Allocation/')[which_machine]

save_dir = paste0(c('/Users/zackoyafuso/Google Drive/',
                    'C:/Users/Zack Oyafuso/Google Drive/')[which_machine],
                  'MS_Optimizations/figure_plot/')

######################################
## load data
######################################
setwd(result_wd)
load('../Extrapolation_depths.RData')
VAST_model = '6g'

######################################
##Which strata to plot
######################################
nstrata = c('0.15' = 10, '0.2' = 8)

{
  tiff(paste0(save_dir, 'solution_ST.tiff'),
       res = 400, width =190, height = 100, units = 'mm', compression = 'lzw')
  
  #Plot settings
  par(oma = rep(0,4), family = 'serif' )
  layout(mat = matrix(c(1,2,3,4), ncol = 2), heights = c(.5,1))
  
  for(CV in c(0.15, 0.20)){ #for each CV constraint
    
    #Subset data depending on CV constraint
    load(paste0('model_', VAST_model, '/optimization_ST_master.RData'))
    which_cv = which(settings$cv[1:length(strata_list)] == CV 
                     & settings$nstata == nstrata[paste0(CV)])
    settings = settings[which_cv,]
    strata_list = strata_list[which_cv]
    res_df = as.data.frame(res_df)
    res_df = res_df[,which_cv] 
    
    #Find the solution wiht the best (lowest) sample size
    sample_sizes = sapply(strata_list, FUN = function(x) sum(x$Allocation))
    winner = names(sample_sizes)[which.min(sample_sizes)]
    
    #Set up spatial object
    ns = 15 #Number of species
    goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')], 
                                 data = data.frame(X1=res_df[,winner]) )
    
    
    goa_ras = raster(goa, resolution = 5)
    goa_ras = rasterize(x = goa, y = goa_ras, field = 'X1')
    
    temp_df = strata_list[[winner]]
    
    #Plot
    par(mar = c(0,0,0,0))
    image(goa_ras, asp = 1, axes = F,
          col =brewer.pal(n = nstrata[paste(CV)], 'Paired'))
    
    xrange = diff(par()$usr[1:2])
    yrange = diff(par()$usr[3:4])
    
    text(x = par()$usr[1]+xrange*0.72,
         y = par()$usr[3]+yrange*0.475,
         paste0(CV*100, '% CV Constraint\n',
                'Optimal Sample Size: ', sum(temp_df$Allocation) ),
         cex = 1, font = 2)

    #Plot stratification cuts
    depth_cuts = unique(ceiling(sort(c(temp_df$Lower_X1, temp_df$Upper_X1))))
    lon_cuts = unique(ceiling(sort(c(temp_df$Lower_X2, temp_df$Upper_X2))))
    
    matrix_space = matrix(data = nstrata[paste(CV)] + 1,
                          nrow = length(depth_cuts)-1, 
                          ncol = length(lon_cuts)-1,
                          dimnames = list(
                            paste0(depth_cuts[1:(length(depth_cuts)-1)], '-',
                                   depth_cuts[2:length(depth_cuts)]),
                            paste0(lon_cuts[1:(length(lon_cuts)-1)], '-',
                                   lon_cuts[2:length(lon_cuts)]) ) )
    
    for(istrata in nstrata[paste(CV)]:1){
      depth_bounds = ceiling(unlist(temp_df[istrata,c('Lower_X1','Upper_X1')]))
      col_idx = c(which(depth_cuts[1:(length(depth_cuts)-1)] == depth_bounds[1]),
                  which(depth_cuts[2:length(depth_cuts)] == depth_bounds[2] ) )
      
      lon_bounds = ceiling(unlist(temp_df[istrata,c('Lower_X2','Upper_X2')]))
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
        matrix_space[i,j] = ifelse(any(frame$X1 > depth_cuts[i] & 
                                         frame$X1 < depth_cuts[i+1] &
                                         frame$X2 > lon_cuts[j] & 
                                         frame$X2 < lon_cuts[j+1]), 
                                   matrix_space[i,j], 
                                   nstrata[paste(CV)]+1)
      }
    }
    
    par(mar = c(3,3,3,1.5))
    sol_ras = raster(matrix_space, xmn = 0, xmx = ncol(matrix_space),
                     ymn = 0, ymx = nrow(matrix_space) )
    image(sol_ras, axes = F, ann = F,
          col = c( brewer.pal(n = nstrata[paste(CV)], 'Paired'),
                   'black'))
    mtext(side = 1, 'Eastings (km)', line = 1.9, cex = 0.75)
    mtext(side = 2, 'Depth (m)', line = 2.2, cex = 0.75)
    axis(side = 1, at = 0:(sol_ras@extent[2]) ,
         labels = lon_cuts, las = 2, cex.axis = 0.75)
    axis(side = 2, at = (sol_ras@extent[4]):0 ,
         labels = depth_cuts, las = 1, cex.axis = 0.75)

    abline(h = 0:length(depth_cuts), v = 0:length(lon_cuts))
    
    #Plot legend for each stratum
    legend(x = 0, 
           y = c('0.15' = 12.75, '0.2' = 11)[paste(CV)], 
           horiz = F, xpd = NA, ncol = 2,
           legend = paste0('Stratum ', 1:nstrata[paste(CV)], ': ', 
                           temp_df$Population, ' units (n = ', 
                           temp_df$Allocation, ')'),
           cex = 0.7, bty = 'n', 
           fill = brewer.pal(n = nstrata[paste(CV)], 'Paired'))
  }

  dev.off()
}

