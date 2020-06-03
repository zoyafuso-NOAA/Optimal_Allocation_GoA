##################################
## Assign bathymetry value to each extrapolation grid
## Gulf of Alaska
#################################
rm(list = ls())

##################################
## Import Libraries
#################################
library(marmap); library(sp); library(RANN); library(raster);

which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[2]

##################################
## Load Extrapolation Grid used in VAST
#################################
modelno = '8'

setwd( c('',
         'C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/',
         'C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/')[which_machine] )

load(paste0('VAST_output', modelno, '/', 'Spatial_Settings.RData'))
observed_depths = read.csv('data/data/GOA_multspp.csv')$BOTTOM_DEPTH 

##################################
## Extract fine-scale bathymetry map 
## Convert latlon to UTM zone 5
#################################
xmin <- -170
xmax <- -132
ymin <- 52
ymax <- 60.5
bathymap <- getNOAA.bathy(lon1 = xmin, lon2 = xmax,
                          lat1 = ymin, lat2 = ymax,
                          resolution = 1)
bathymap <- fortify.bathy(bathymap)

bathymap_coord = sp::SpatialPoints(coords = bathymap[,c('x', 'y')],
                                   proj4string = CRS('+proj=longlat') )
cord.UTM <- sp::spTransform(bathymap_coord, CRS("+proj=utm +zone=5N"))
bathymap[,c('E_km', 'N_km')] = cord.UTM@coords / 1000

#####################
## Asssign bathymetry values of extrapolation cells to the nearest
## value in the bathymetry map
#####################

bathy_idx = RANN::nn2(query=Extrapolation_List$Data_Extrap[,c('E_km','N_km')],
                      data = bathymap[,c('E_km', 'N_km')],
                      k = 1)$nn.idx

Extrapolation_List$Data_Extrap$depth = -bathymap$z[bathy_idx[,1]]

#####################
## Plot locations where depths are negative (land?) 
#####################
plot(Extrapolation_List$Data_Extrap[,c('E_km', 'N_km')], pch = '.')
points(Extrapolation_List$Data_Extrap[Extrapolation_List$Data_Extrap$depth <=0,c('E_km', 'N_km')], pch = 16, col = 'red')

#####################
## Assign negative bathymetry values (presumably land) to the shallowest
## bathymetry observed in the dataset
#####################
Extrapolation_List$Data_Extrap$depth[Extrapolation_List$Data_Extrap$depth <= min(observed_depths)] = min(observed_depths)

# neg_depths = sum(Extrapolation_List$Data_Extrap$depth <=0)
# k = 2
# while(neg_depths != 0){
#   idxs = which(Extrapolation_List$Data_Extrap$depth <= 0)
#   Extrapolation_List$Data_Extrap$depth[idxs] = -bathymap$z[bathy_idx[idxs,k]]
#   neg_depths = sum(Extrapolation_List$Data_Extrap$depth <=0)
#   k = k + 1
# }

#############################
## Center depth and calculate depth^2
#############################
Extrapolation_List$Data_Extrap$DEPTH = scale(log(Extrapolation_List$Data_Extrap$depth))
Extrapolation_List$Data_Extrap$DEPTH2 = Extrapolation_List$Data_Extrap$DEPTH^2

#############################
## Plot Bathyetry Field
#############################
test = raster::rasterize(
  x = Extrapolation_List$Data_Extrap[,c('E_km', 'N_km')],
  y = raster(nrows=100, ncols=320,
             xmn=min(Extrapolation_List$Data_Extrap$E_km),
             xmx=max(Extrapolation_List$Data_Extrap$E_km),
             ymn=min(Extrapolation_List$Data_Extrap$N_km),
             ymx=max(Extrapolation_List$Data_Extrap$N_km),
             crs = CRS("+proj=utm +zone=5N")),
  field = Extrapolation_List$Data_Extrap$depth)

par(mar = c(0,0,0,0))
plot(test, axes = F, legend = F)

##########################
## Update Spatial Settings
##########################
Extrapolation_depths = Extrapolation_List$Data_Extrap

save(list = c("Extrapolation_depths"), file = 'Extrapolation_depths.RData')

