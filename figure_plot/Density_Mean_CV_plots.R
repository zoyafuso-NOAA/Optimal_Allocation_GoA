#################################
## Plot the mean and temporal CV of density for each species
#################################

rm(list = ls())

#######################################
## Import Libraries
#######################################
library(VAST); library(mvtnorm); library(SamplingStrata); 
library(sp); library(raster)

#######################################
## Set up directories
#######################################
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2)[2]

VAST_model = "6g"
VAST_dir = paste0(c('/Users/zackoyafuso/Google Drive/', 
                    'C:/Users/Zack Oyafuso/Google Drive/')[which_machine],
                  'VAST_Runs/VAST_output', VAST_model, '/')

github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/')[which_machine],
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

#######################################
## Extract Densities for each species across the extrapolation grid in time
## Calculate Mean density for each extrapolation grid for each species
## Calculate the temporal variance of density for each species across grids
## Calculate the temporal CV of density for each species across grids
#######################################
Ests = Save$Report$Index_gcyl[,,Years2Include,1]
TotVar = apply(X = Ests, MARGIN = c(1,2), FUN = function(x) var(as.vector(x)) )
Ests_means =apply(Ests, MARGIN = 1:2, mean)
Ests_CV = sqrt(TotVar)/ Ests_means


yrange = diff(range(Extrapolation_List$Data_Extrap[,c('N_km')]))
xrange = diff(range(Extrapolation_List$Data_Extrap[,c('E_km')]))


{
    png(paste0(figure_dir, 'Mean_CV.png'),
     width = 190, height = 190, units = 'mm', res = 500)
par(mar = c(0,0,0,0), mfrow = c(5,3))

for(i in 1:length(Save$Spp)){
    plot(1, type = 'n', axes = F, ann = F,
         xlim = range(Extrapolation_List$Data_Extrap[,c('E_km')]),
         ylim = c(min(Extrapolation_List$Data_Extrap[,c('N_km')])-1*yrange,
                  max(Extrapolation_List$Data_Extrap[,c('N_km')]))
    )
    
    temp = as.data.frame(cbind(Ests_means[,i], Ests_CV[,i]))
    names(temp) = c("Mean", 'CV')
    
    goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 
                                                                  'N_km')], 
                                 data = temp)
    
    for(itype in c('Mean', 'CV')){
        goa_ras = raster(goa, resolution = 5)
        goa_ras =rasterize(x = goa, y = goa_ras, field = itype)
        
        vals  = values(goa_ras)
        
        if(itype == 'CV') val_cuts = c(0,0.25,0.5,1,2,10)
        
        if(itype == 'Mean')     {
            vals = vals[vals > 1]
            val_cuts = c(0, quantile(vals , na.rm = T, 
                                     probs = seq(0,1, length = 5))[-1])
        } 
        
        values(goa_ras) = cut(x = values(goa_ras), breaks = val_cuts)
        
        if(itype == 'CV') goa_ras = raster::shift(x = goa_ras, dy = -yrange)
        
        image(goa_ras, asp = 1, axes = F, ann = F, add = T,
              col = list('Mean' = c(hcl.colors(n=4, "YlOrRd", rev = TRUE)[c(1,2,3)], 'black'), 
                         'CV' = hcl.colors(n=5, palette='viridis', rev=T))[[itype]])
        
        if(itype == 'Mean') val_cuts = round(val_cuts)
        if(itype == 'CV') val_cuts = c(0,0.25,0.5,1,2,10)
        legend(x = goa_ras@extent[1] + xrange*0.55,
               y = goa_ras@extent[3] + yrange*0.7,
               bty = 'n', cex = 0.75,
               fill = list('Mean' =c(hcl.colors(n=4, "YlOrRd", rev = TRUE)[c(1,2,3)], 'black'), 
                           'CV' = hcl.colors(n=5, palette='viridis', rev=T))[[itype]], 
               legend = paste0(val_cuts[1:(length(val_cuts)-1)], '-',
                               val_cuts[2:length(val_cuts)]) )
        
        if(itype == 'Mean') text(x = goa_ras@extent[1] + xrange*0.15, 
                                 y = goa_ras@extent[4] - yrange*0.15,
                                 gsub(Save$Spp[i], pattern = ' ', 
                                      replacement = '\n'),
                                 font = 3)
    }
    box()
}

dev.off()
}
