######################################
## Where do 1%/5% of the highest catches occur for each species
######################################
rm(list = ls())

######################################
## Import Libraries
######################################
library(raster); library(RColorBrewer); library(SamplingStrata)

######################################
##Set up directories
######################################
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2)[1]
github_dir = paste0(c('', 
                      'C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/Optimum_Allocation/')

VAST_model = '6g'
VAST_dir = paste0(c('', 
                    'C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/',
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output', VAST_model)

######################################
## Import data
######################################
setwd(github_dir)

load('../Extrapolation_depths.RData')
load(paste0(VAST_dir, '/Spatial_Settings.RData'))
Data_Geostat[,c("E_km", "N_km")] = Spatial_List$loc_i

nstrata = c('0.15' = 10, '0.2' = 8)
yrange = diff(range(Extrapolation_depths[,c('N_km')]))

for(percentile in c(5,1)){
  {tiff(paste0('../figure_plot/top_', percentile, 'percent_CPUE.tiff'),
        res = 600, width =190, height = 100, units = 'mm', compression = 'lzw')
    par(mfrow = c(1,3), mar = c(0,0,0,0), oma = c(0,0,0,0), family = 'serif')
    for(ispp in 1:15){
      
      if(ispp%%5 == 1) {
        
        plot(1, type = 'n', axes = F, ann = F,
             xlim = range(Extrapolation_depths[,c('E_km')]),
             ylim = c(min(Extrapolation_depths[,c('N_km')])-3*yrange,
                      max(Extrapolation_depths[,c('N_km')])))
        offset = 0
      }
      
      for(CV in c(0.15, 0.20)[1]){

        load(paste0('model_', VAST_model, '/optimization_ST_master.RData'))
        which_cv = which(settings$cv[1:length(strata_list)] == CV 
                         & settings$nstata == nstrata[paste0(CV)])
        settings = settings[which_cv,]
        strata_list = strata_list[which_cv]
        res_df = as.data.frame(res_df)
        res_df = res_df[,which_cv] 
        
        sample_sizes = sapply(strata_list, FUN = function(x) sum(x$Allocation))
        
        winner = names(sample_sizes)[which.min(sample_sizes)]
        ns = 15
        
        ######################################
        ## Set up spatial object
        ######################################
        goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')], 
                                     data = data.frame(X1=res_df[,winner]) )
        
        
        goa_ras = raster(goa, resolution = 5)
        goa_ras = rasterize(x = goa, y = goa_ras, field = 'X1')
        
        goa_ras = raster::shift(goa_ras, dy = -yrange*offset*0.7)
        
        xrange = diff(goa_ras@extent[1:2])
        new_yrange = diff(goa_ras@extent[3:4])
        
        image(goa_ras, asp = 1, axes = F, ann = F, add = T,
              col =brewer.pal(n = nstrata[paste(CV)], 'Paired'))
        
        temp_df = subset(Data_Geostat, 
                         subset = spp == levels(Data_Geostat$spp)[ispp] )
        temp_df$N_km = temp_df$N_km -yrange*offset*0.7
        toohigh_cutoff = quantile(temp_df$Catch_KG/temp_df$AreaSwept_km2, 
                                  probs = c(1-percentile/100, 1),)[1]
        temp_df$toohigh = temp_df$Catch_KG /temp_df$AreaSwept_km2 >= toohigh_cutoff
        
        points(N_km ~ E_km, data = temp_df, subset = toohigh == T, pch = 16,
               cex = 0.25)
        
        text(x = goa_ras@extent[1] + 0.725*xrange,
             y = goa_ras@extent[3]+ 0.55*new_yrange,
             gsub(levels(Data_Geostat$spp)[ispp], pattern = ' ', 
                  replacement = '\n'), font = 3 ) 
        offset = offset + 1
      }
      
    }
    dev.off()}
  
}
