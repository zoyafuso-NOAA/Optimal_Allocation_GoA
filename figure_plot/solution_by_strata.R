#####################################
## Optimal Solutions from 5-20 strata
#####################################

library(sp); library(raster); library(RColorBrewer)

VAST_Model = '6g'
results_wd = paste0(c('C:/Users/Zack Oyafuso/Documents/'),
                    'GitHub/MS_OM_GoA/Optimum_Allocation/model_', 
                    VAST_Model, '/')

setwd(results_wd)

load('../../Extrapolation_depths.RData')
load('optimization_ST_master.RData')
strata_list = strata_list[settings$mut_change == 0.10 & settings$elitism_rate == 0.10 & settings$cv == 0.15]
res_df = res_df[,settings$mut_change == 0.10 & settings$elitism_rate == 0.10 & settings$cv == 0.15]
settings = subset(settings, (mut_change == 0.10 & elitism_rate == 0.10 & cv == 0.15))

settings$n = sapply(strata_list, FUN = function(x) sum(x$Allocation))

best_sol = aggregate(n ~ nstata + cv, data = settings, FUN = min)

yrange = diff(range(Extrapolation_depths[,c('N_km')]))

{tiff(paste0('../../figure_plot/solution_by_strata.tiff'),
      res = 400, width =190, height = 60, units = 'mm', compression = 'lzw')
  par(mfrow = c(1,4), mar = c(0,0,0,0), family = 'serif')
  for(istrata in 1:nrow(best_sol)) {
    if(istrata%%4 == 1) {
      plot(1, type = 'n', axes = F, ann = F,
           xlim = range(Extrapolation_depths[,c('E_km')]),
           ylim = c(min(Extrapolation_depths[,c('N_km')])-2.25*yrange,
                    max(Extrapolation_depths[,c('N_km')])))
      offset = 0
    }
    
    #rownames of settings with the strata number
    row_idx = row.names(settings)[settings$nstata == best_sol$nstata[istrata]]
    best_sol_idx = which.min(settings[row_idx, 'n'])
    idx = row_idx[best_sol_idx]
    
    goa = SpatialPointsDataFrame(coords=Extrapolation_depths[,c('E_km','N_km')], 
                                 data=data.frame(X1=res_df[,paste0('sol_',idx)]) )
    
    goa_ras = raster(goa, resolution = 5)
    goa_ras = rasterize(x = goa, y = goa_ras, field = 'X1')
    
    goa_ras = raster::shift(goa_ras, dy = -yrange*offset*0.75)
    
    image(goa_ras, asp = 1, axes = F, add = T,
          col = c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
                  '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a',
                  '#ffff99', '#8dd3c7', 'brown', '#bebada', '#fb8072',
                  '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9',
                  '#bc80bd', '#ccebc5')[1:best_sol$nstata[istrata]])
    text(x = goa_ras@extent[1] + 0.7*diff(goa_ras@extent[1:2]),
         y = goa_ras@extent[3]+ 0.6*diff(goa_ras@extent[3:4]),
         paste(best_sol$nstata[istrata], 'strata\n',
         'n =', settings[idx,'n']) ) 
    
    offset = offset + 1
  }
  dev.off()}
