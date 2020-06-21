#################################
## Plot the mean and temporal CV of density for each species
#################################

rm(list = ls())

#######################################
## Import Libraries
#######################################
library(sp); library(raster); library(RColorBrewer); library(rworldxtra)
library(plotrix); library(rgdal); library(oce); library(ocedata)
library(rnaturalearth)

#######################################
## Set up directories
#######################################
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[1]

VAST_model = "6g"
VAST_dir = paste0(c('/Users/zackoyafuso/Google Drive/', 
                    'C:/Users/Zack Oyafuso/Google Drive/',
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output', VAST_model, '/')

github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/',
                      'C:/Users/zack.oyafuso/Work/')[which_machine],
                    'GitHub/Optimal_Allocation_GoA/')

PP_dir = paste0(c('/Users/zackoyafuso/Google Drive/', 
                  'C:/Users/Zack Oyafuso/Google Drive/')[which_machine],
                'MS_Optimizations/powerpoint_plot/')

figure_dir = paste0(c('/Users/zackoyafuso/Google Drive/', 
                      'C:/Users/Zack Oyafuso/Google Drive/')[which_machine],
                    'MS_Optimizations/figure_plot/')

#######################################
## Load Shapefiles
#######################################
#World Data for inset
data(coastlineWorld)
data(coastlineWorldFine, package="ocedata")

#Extract USA for AK land
# data(countriesHigh)
# AK = sp::spTransform(countriesHigh[c(39,236),], 
#                      crs('+proj=utm +zone=5 +units=km'))

AK <- rnaturalearth::ne_countries(continent = "north america", scale = "medium", returnclass = "sp")
AK <- sp::spTransform(AK, CRS = "+proj=utm +zone=5 +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")

#######################################
## Load Fit, spatial data
#######################################
load(paste0(VAST_dir,'VAST_MS_GoA_Run.RData'))
load(paste0(VAST_dir,'Spatial_Settings.RData'))
load(paste0(github_dir, 'data/Extrapolation_depths.RData'))

###############################################
## Spatial, Temporal, and Species Settings
###############################################

#Boundary Box of the extrapolation grid
bbox_ = c('xmin' = min(Extrapolation_List$Data_Extrap[,c('E_km')]),
          'xmax' = max(Extrapolation_List$Data_Extrap[,c('E_km')]),
          'ymin' = min(Extrapolation_List$Data_Extrap[,c('N_km')]),
          'ymax' = max(Extrapolation_List$Data_Extrap[,c('N_km')]) )

#x- and y-ranges of the extrapolation grid, used in plotting relative positions
yrange = diff(range(Extrapolation_List$Data_Extrap[,c('N_km')]))
xrange = diff(range(Extrapolation_List$Data_Extrap[,c('E_km')]))

#which years to include from the OM
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

ns = Save$TmbData$n_c #Number of Species
sci_names = c(Save$Spp, 'Bathymetry (meters)')

#######################################
## Create 200 m isobar
#######################################
bathy = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 
                                                                'N_km')], 
                               data = data.frame(val = Extrapolation_depths$depth))
bathy_ras = raster(bathy, resolution = 5)
bathy_ras = rasterize(x = bathy, y = bathy_ras, field = 'val')

bathy_iso_200 = raster::rasterToContour(x = bathy_ras, levels = c(200))

#######################################
## Extract Densities and Calculate Mean Density across years
#######################################
density_gct = Save$Report$D_gcy[,,Years2Include]
mean_density_gc = apply(density_gct, MARGIN = 1:2, mean)

#######################################
## Plot
#######################################
{
    #Set up file
    png(paste0(figure_dir, 'Mean_CV.png'),
        width = 190, height = 190, units = 'mm', res = 500)
    
    #Panel Layout: 6 rows by 3 columns, first 5 rows are spatial distributions
    #for the 15 species and the last row (2nd column) is for the inset.
    par(mar = c(0,0,0,0), oma = c(2,1,1,1))
    layout(mat = matrix(c(1:15,18,17,16), ncol = 3, byrow = T))
    
    for(ispp in 1:(ns+1)){     #For Each species
        
        #Base Plot
        plot(1, type = 'n', axes = F, ann = F,
             xlim = bbox_[c('xmin', 'xmax')],
             ylim = c(bbox_['ymin'], bbox_['ymax']))
        box()
        
        #Extract Data for a species
        if(ispp %in% 1:ns) temp = data.frame(mean = mean_density_gc[,ispp])
        if(ispp == 16) temp = data.frame(mean = Extrapolation_depths$depth)
        
        #Discretize by decile: the first bin is for densities < 1, then the 
        #densities > 1 are broken up into 9 quantile groups. In total, there 
        #are 10 density bins. This was done so that the color scales aren't 
        #skewed too much by very low and very high values.
        if(ispp %in% 1:ns){
            temp_int = as.integer(
                cut(x = temp$mean, 
                    breaks = c(0,1,quantile(temp$mean[temp$mean > 1], 
                                            probs = seq(0,1,length=10))[-1])))
        }
        
        if(ispp == 16){
            temp_int = as.integer(
                cut(x = temp$mean, include.lowest = T,
                    breaks = quantile(temp$mean, 
                                      probs = seq(0,1,length=11) )))
        }
        
        #Spatial Object: we create a raster by creating a spatial object from
        #the positions of the Extrapolation grid, then attributing the density
        #values across the grid. Then that object is rasterized.
        goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km',
                                                                      'N_km')],
                                     data = data.frame(val = temp_int))
        
        goa_ras = raster(goa, resolution = 5)
        goa_ras = rasterize(x = goa, y = goa_ras, field = 'val')
        
        #Plot density raster
        color_density = ifelse(ispp %in% 1:ns, 'Blues', 'OrRd')
        image(goa_ras, col = c('white', brewer.pal(n=9, name=color_density)),
              axes = F, asp = 1, ann = F, add = T)
        
        #Overlay 200 m isopleth
        lwd_density = ifelse(ispp %in% 1:ns, 0.5, 2)
        plot(bathy_iso_200, add = T, lty = 1, col = 'black', lwd = 0.5)
        
        #Isopleth legend
        legend(x = bbox_[1] + xrange*0.575,
               y = bbox_[3] + yrange*0.8,
               legend = '200 m\nisopleth', 
               cex = 0.75, lty = 1, bty = 'n', lwd = 2, col = 'black')
        
        #Overlay Land
        plot(AK, add = T, col = 'tan', border = F)
        
        #Overlay Species Name
        text(x = bbox_[1] + xrange*0.15, 
             y = bbox_[3] + yrange*0.55,
             gsub(sci_names[ispp], pattern = ' ', replacement = '\n'), 
             font = ifelse(ispp %in% 1:ns, 3, 1))
        
        #Overlay scale bar
        segments(x0 = bbox_[1] , x1 = bbox_[1] + 500,
                 y0 = bbox_[3] + yrange*0.9, lwd = 1)
        
        #Overlay scale bar ticks
        segments(x0 = seq(from = bbox_[1], to = bbox_[1] + 500, length = 6),
                 y0 = bbox_[3] + yrange*0.90,
                 y1 = bbox_[3] + yrange*0.93, lwd = 1)
        
        #Overlay scale bar label
        text(x = bbox_[1] + 250, y = bbox_[3] + yrange*0.85, 
             '500 km', cex = 0.8)
        
        #Create legend values
        if(ispp %in% 1:ns){
            legend_vals = list(
                top = c(
                    "< 1",
                    round(quantile(temp$mean[temp$mean > 1],
                                   probs = seq(0,1,length=10)))[c(NA,3,NA, 5,
                                                                  NA, 7,NA, 
                                                                  9,NA)]),
                bottom = c(round(quantile(temp$mean[temp$mean > 1],
                                          probs = seq(0,1,length=10)))[c(NA,2,NA,4,
                                                                         NA,6,NA,8,
                                                                         NA,10)]) )
        }
        
        if(ispp == 16) {
            legend_vals = list(
                top = round(quantile(temp$mean[temp$mean > 1],
                                     probs = seq(0,1,length=10)))[c(1,NA,3,NA, 5,
                                                                    NA, 7,NA, 
                                                                    9,NA)],
                bottom = round(quantile(temp$mean[temp$mean > 1],
                                          probs = seq(0,1,length=10)))[c(NA,2,NA,4,
                                                                         NA,6,NA,8,
                                                                         NA,10)])
        }
        
        #Overlay Legend
        plotrix::color.legend(
            xl = bbox_[1] + xrange*0.4,
            xr = bbox_[1] + xrange*1,
            yb = bbox_[3] + yrange*0.05,
            yt = bbox_[3] + yrange*0.15,
            legend = legend_vals[['top']],
            rect.col = c('white', brewer.pal(n = 9, name = color_density)), 
            cex = 0.5, align = 'rb' )
        plotrix::color.legend(
            xl = bbox_[1] + xrange*0.4,
            xr = bbox_[1] + xrange*1,
            yb = bbox_[3] + yrange*0.05,
            yt = bbox_[3] + yrange*0.15,
            legend = legend_vals[['bottom']],
            rect.col = c('white', brewer.pal(n = 9, name = color_density)), 
            cex = 0.5, align = 'lt' )
        
    }
    
    #Plot Inset
    par(mar = c(0,1,0,1))
    mapPlot(coastlineCut(coastlineWorldFine, -135),
            longitudelim=c(-180, -70), latitudelim=c(30, 60),
            proj="+proj=aea +lat_0=40 +lon_0=-135", col='tan',
            border = 'black', grid = c(15,15))
    rect(xleft = -2000000, xright = 100000,
         ybottom = 1600000, ytop = 2600000,
         lwd = 3)
    
    
    #Empty Plot
    plot(1, type = 'n', axes = F, ann = F)
    
    dev.off()
}
