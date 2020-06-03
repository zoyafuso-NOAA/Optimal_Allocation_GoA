#################################
## Optimization Data
#################################
rm(list = ls())

###############################
## Import required packages
###############################
library(VAST);  library(mvtnorm); library(sp); library(RColorBrewer); 
library(raster)
library(memoise); library(doParallel); library(foreach); library(iterators); 
library(parallel); library(pbapply); library(formattable)


###############################
## Set up directories
###############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[2]
optimization_type = c('_spatial', '')[1]
VAST_model = "6g"

SamplingStrata_dir = paste0(c('', 
                              'C:/Users/Zack Oyafuso',
                              'C:/Users/zack.oyafuso',
                              'C:/Users/zack.oyafuso')[which_machine],
                            '/Downloads/SamplingStrata-master/R')
github_dir = paste0(c('', 
                      'C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/MS_OM_GoA/Optimum_Allocation/')

VAST_dir = paste0(c('', 
                    'C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/',
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output', VAST_model)


output_wd = paste0(c('/Users/zackoyafuso/Documents/', 
                     'C:/Users/Zack Oyafuso/Documents/',
                     'C:/Users/zack.oyafuso/Work/', 
                     'C:/Users/zack.oyafuso/Work/' )[which_machine], 
                   "GitHub/MS_OM_GoA/Optimum_Allocation/model_", VAST_model,
                   optimization_type)

#########################
## Load functions from SamplingStrata packages into global environment
## Load modified buildStrataDF function
#########################
for(ifile in dir(SamplingStrata_dir, full.names = T)) source(ifile)
source(paste0(github_dir, '/buildStrataDF_Zack.R'))

#########################
## Load VAST products
#########################
load(paste0(VAST_dir, '/VAST_MS_GoA_Run.RData'))
load(paste0(VAST_dir, '/Spatial_Settings.RData'))
load(paste0(dirname(github_dir), '/Extrapolation_depths.RData'))

#########################
## Index years that had data
#########################
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
NTime = length(Years2Include)

##########################
## Create the data inputs to SamplingStrata
##########################
df = df_raw = NULL

##########################
## Mean density across years
#########################
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

############################
## Density for each observed year and cell
############################
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

ns = Save$TmbData$n_c

rm(list = c('Save', 'Spatial_List', 'spp_df', 'strata.limits', 'fine_scale',
            'Method', 'modelno', 'n_x', 'which_spp', 'Year_Set', 
            'Years2Include', 'Data_Geostat', 'df', 'Extrapolation_List',
            'gulf_of_alaska_grid', 'ifile', 'iT', 'df_raw'))

save(list = c('Extrapolation_depths', 'frame', 
              'frame_raw', 'ns', 'NTime'),
     file = paste0(output_wd, '/optimization_data_model_', 
            VAST_model, '.RData'))
