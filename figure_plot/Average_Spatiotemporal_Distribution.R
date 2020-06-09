#################################
## Plot the mean and temporal CV of density for each species
#################################

rm(list = ls())

#######################################
## Import Libraries
#######################################
library(sp); library(raster); library(RColorBrewer); library(plotrix)


#######################################
## Set up directories
#######################################
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[3]

VAST_model = "6g"
VAST_dir = paste0(c('/Users/zackoyafuso/Google Drive/', 
                    'C:/Users/Zack Oyafuso/Google Drive/',
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output', VAST_model, '/')

github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/',
                      'C:/Users/zack.oyafuso/Work/')[which_machine],
                    'GitHub/MS_OM_GoA/')

PP_dir = paste0(c('/Users/zackoyafuso/Google Drive/', 
                  'C:/Users/Zack Oyafuso/Google Drive/')[which_machine],
                'MS_Optimizations/powerpoint_plot/')

figure_dir = paste0(c('/Users/zackoyafuso/Google Drive/', 
                      'C:/Users/Zack Oyafuso/Google Drive/')[which_machine],
                    'MS_Optimizations/figure_plot/')

#######################################
## Load Fit, spatial data
#######################################
load(paste0(VAST_dir,'VAST_MS_GoA_Run.RData'))
load(paste0(VAST_dir,'Spatial_Settings.RData'))
load(paste0(github_dir, 'Extrapolation_depths.RData'))

Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

bbox_ = c('xmin' = min(Extrapolation_List$Data_Extrap[,c('E_km')]),
          'xmax' = max(Extrapolation_List$Data_Extrap[,c('E_km')]),
          'ymin' = min(Extrapolation_List$Data_Extrap[,c('N_km')]),
          'ymax' = max(Extrapolation_List$Data_Extrap[,c('N_km')]) )

yrange = diff(range(Extrapolation_List$Data_Extrap[,c('N_km')]))
xrange = diff(range(Extrapolation_List$Data_Extrap[,c('E_km')]))

ns = Save$TmbData$n_c
NTime = length(Years2Include)

#######################################
## Extract Average Spatial Distribution (Omega)
## 
#######################################
density_gct = Save$Report$D_gcy[,,Years2Include]
sd_density_gc = apply(density_gct, MARGIN = 1:2, sd)
mean_density_gc = apply(density_gct, MARGIN = 1:2, mean)

{
    png(paste0(figure_dir, 'Mean_CV.png'),
     width = 190, height = 190, units = 'mm', res = 500)
    par(mar = c(0,0,0,0), mfrow = c(5,3), oma = c(1,1,1,1))

    # png(paste0(figure_dir, 'Mean_CV.png'),
    #  width = 190, height = 190, units = 'mm', res = 500)
par(mar = c(0,0,0,0), mfrow = c(5,3))

for(ispp in 1:ns){
    
    #Base Plot
    plot(1, type = 'n', axes = F, ann = F,
         xlim = bbox_[c('xmin', 'xmax')],
         ylim = c(bbox_['ymin']-1*yrange, bbox_['ymax']))
    
    #Extract Data for a species
    temp = data.frame(mean = mean_density_gc[,ispp])#,
                      #cv = cv_density_gc[,ispp])
    
    temp_int = as.numeric(cut(x = temp$mean, 
                              breaks = quantile(temp$mean, 
                                                probs = seq(0,1,.1)) ))
    
    for(ispp in 1:ns){
        
        #Base Plot
        plot(1, type = 'n', axes = F, ann = F,
             xlim = bbox_[c('xmin', 'xmax')],
             ylim = c(bbox_['ymin'], bbox_['ymax']))
        box()
        
        #Extract Data for a species
        temp = data.frame(mean = mean_density_gc[,ispp])
        
        temp_int = as.integer(cut(x = temp$mean, 
                                  breaks = c(0,1,quantile(temp$mean[temp$mean > 1], 
                                                          probs = seq(0,1,length=10))[-1] )))
        
        #Spatial Object
        goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 
                                                                      'N_km')], 
                                     data = data.frame(val = temp_int))
        
        goa_ras = raster(goa, resolution = 20)
        goa_ras = rasterize(x = goa, y = goa_ras, field = 'val')
        image(goa_ras, col = c(brewer.pal(n = 9, name = 'YlOrRd'), 'black'),
              axes = F, asp = 1, ann = F, add = T)
        
        text(x = bbox_[1] + xrange*0.15, 
             y = bbox_[4] - yrange*0.15,
             gsub(Save$Spp[ispp], pattern = ' ', replacement = '\n'), font = 3)
        
        
        legend_vals = list(
            top = c("< 1", 
                    round(quantile(temp$mean[temp$mean > 1], 
                                   probs = seq(0,1,length=10)))[c(NA,3,NA, 5,NA, 7,NA, 9,NA)]),
            bottom = c(round(quantile(temp$mean[temp$mean > 1], 
                                      probs = seq(0,1,length=10)))[c(NA,2,NA,4,NA,6,NA,8,NA,10)]) )
        
        plotrix::color.legend(
            xl = bbox_[1] + xrange*0.4, 
            xr = bbox_[1] + xrange*1, 
            yb = bbox_[3] + yrange*0.05, 
            yt = bbox_[3] + yrange*0.15,
            legend = legend_vals[['top']], 
            rect.col = c(brewer.pal(n = 9, name = 'YlOrRd'),
                         'black'), cex = 0.7, align = 'rb' )
        plotrix::color.legend(
            xl = bbox_[1] + xrange*0.4, 
            xr = bbox_[1] + xrange*1, 
            yb = bbox_[3] + yrange*0.05, 
            yt = bbox_[3] + yrange*0.15,
            legend = legend_vals[['bottom']], 
            rect.col = c(brewer.pal(n = 9, name = 'YlOrRd'),
                         'black'), cex = 0.7, align = 'lt' )
        
    }
    
    dev.off()
    
}
