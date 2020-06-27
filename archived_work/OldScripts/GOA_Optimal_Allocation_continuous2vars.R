#################################
## Sampling Strata 
## Method == "continous"
#################################
rm(list = ls())

library(VAST); 
library(mvtnorm); library(SamplingStrata); library(sp)
library(RColorBrewer); library(raster)

which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[3]

VAST_wd = c('/Users/zackoyafuso/Google Drive/VAST_Runs/',
            'C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/',
            'C:/Users/zack.oyafuso/Desktop/VAST_Runs/')[which_machine]

VAST_model = "6g"
output_wd = c(paste0('/Users/zackoyafuso/Documents/GitHub/MS_OM_GoA/',
                     'Optimum_Allocation/model_', VAST_model),
              paste0("C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/",
                     "Optimum_Allocation/model_", VAST_model),
              paste0("C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/",
                     "Optimum_ALlocation/model_", VAST_model))[which_machine]

setwd(VAST_wd)

if(!dir.exists(output_wd)) dir.create(output_wd)

load(paste0(VAST_wd, 'VAST_output',VAST_model,'/VAST_MS_GoA_Run.RData'))
load(paste0(VAST_wd, 'VAST_output',VAST_model,'/Spatial_Settings.RData'))


load(paste0(c('/Users/zackoyafuso/Documents/GitHub/MS_OM_GoA/',
              "C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/",
              'C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/'),
            'Extrapolation_depths.RData')[which_machine])

Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
NTime = length(Years2Include)

df = df_raw = NULL

df = cbind(
  data.frame(Domain = 1,
             x = 1:Save$TmbData$n_g,
             lat = Extrapolation_depths$N_km,
             lon = Extrapolation_depths$E_km - min(Extrapolation_depths$E_km),
             depth = Extrapolation_depths$depth),
  apply(X=Save$Report$Index_gcyl[,,Years2Include,], MARGIN = 1:2, FUN = mean ) )

names(df)[-(1:5)] = gsub(x = Save$Spp, pattern = ' ', replacement = '_')

frame <- buildFrameDF(df = df,
                      id = "x",
                      X = c("depth", 'lon'),
                      Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                      domainvalue = "Domain")

for(iT in 1:NTime){
  df_raw = rbind(df_raw, cbind(
    data.frame(Domain = 1,
               x = 1:Save$TmbData$n_g,
               year = iT,
               lat = Extrapolation_depths$N_km,
               lon = Extrapolation_depths$E_km - min(Extrapolation_depths$E_km),
               depth = Extrapolation_depths$depth),
    Save$Report$Index_gcyl[,,Years2Include[iT],] )
  )
}
names(df_raw)[-(1:6)] = gsub(x = Save$Spp, pattern = ' ', replacement = '_')
frame_raw <- buildFrameDF(df = df_raw,
                          id = "x",
                          X = c("depth", 'lon'),
                          Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                          domainvalue = "Domain")

#Settings for optimizer
settings = expand.grid(cv = c(0.2, 0.15),
                       mut_change = c(0.1, 0.01),
                       elitism_rate = c(0.2, 0.1),
                       nstata = c(5,7,10),
                       iter = 1:10)

ns = Save$TmbData$n_c

rm(list = c('Save', 'Spatial_List', 'spp_df', 'strata.limits', 'fine_scale',
            'Method', 'modelno', 'n_x', 'which_spp', 'Year_Set', 
            'Years2Include', 'Data_Geostat', 'df', 'Extrapolation_List',
            'gulf_of_alaska_grid'))

res_df = as.matrix(frame[,c('id', 'domainvalue')])
strata_list = list()

for(i in 1:nrow(settings)){
  par(mfrow = c(2,1))
  
  cv = list()
  for(spp in 1:ns) cv[[paste0('CV', spp)]] = settings$cv[i]
  cv[['DOM']] = 1
  cv[['domainvalue']] = 1
  cv <- as.data.frame(cv)
  
  set.seed(1234 + i)
  solution <- optimStrata(method = "continuous",
                          errors = cv, 
                          framesamp = frame,
                          iter = 50,
                          pops = 50,
                          elitism_rate = settings$elitism_rate[i],
                          mut_chance = settings$mut_change[i],
                          nStrata = settings$nstata[i],
                          showPlot = T,
                          writeFiles = F,
                          parallel = F)
  
  strata_list[[i]] =  summaryStrata(solution$framenew,
                                    solution$aggr_strata,
                                    progress=FALSE) 
  
  res_df = cbind(res_df, solution$framenew$STRATO)
  
  save(list = c('strata_list', 'res_df', 'settings', 'ns', 'frame'), 
       file = paste0(output_wd, '/optimization.RData'))
  
  #Plot
  goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')], 
                               data = cbind(solution$framenew[,paste0('Y',1:ns)],
                                            Str_no = solution$framenew$STRATO,
                                            depth = solution$framenew$X1,
                                            lon = solution$framenew$X2) )
  goa_ras = raster(goa, resolution = 5)
  goa_ras =rasterize(x = goa, y = goa_ras, field = 'Str_no')
  plot(goa_ras, col = terrain.colors(10)[-10], axes = F)
}

