###############################################################################
## Project:       Data synthesis for stratified survey optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Create dataset used for all optimization runs based on a 
##                Gulf of Alaska Multispecies Groundfish VAST Spatiotemporal 
##                Operating Model specified on line 24
##
##                Calculate true mean density across years for each species
##
##                Set up other constants used in downstream processes
###############################################################################

##################################################
####  Import required packages   
##################################################
library(SamplingStrata)

##################################################
####    Set up directories
##################################################
rm(list = ls())

which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[2]
VAST_model = "10b" 

github_dir = paste0(c('/Users/zackoyafuso/Documents', 
                      'C:/Users/Zack Oyafuso/Documents',
                      'C:/Users/zack.oyafuso/Work',
                      'C:/Users/zack.oyafuso/Work')[which_machine],
                    '/GitHub/Optimal_Allocation_GoA/')

VAST_dir = paste0(c('/Users/zackoyafuso/Google Drive/', 
                    'C:/Users/Zack Oyafuso/Google Drive/', 
                    'C:/Users/zack.oyafuso/Desktop/',
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'GOA_VAST_Runs/VAST_output', VAST_model, '/')

##################################################
####   Set up Result Directories
##################################################
if(!dir.exists(paste0(github_dir, 'model_', VAST_model))){
  dir.create(paste0(github_dir, 'model_', VAST_model))
  dir.create(paste0(github_dir, 'model_', VAST_model, 
                    '/Survey_Comparison_Simulations/'))
  dir.create(paste0(github_dir, 'model_', VAST_model, 
                    '/Spatiotemporal_Optimization/'))
}

##################################################
####   Load VAST output and bathymetry data
##################################################
load(paste0(VAST_dir, '/fit.RData'))
load(paste0(github_dir, 'data/Extrapolation_depths.RData'))
spp_df = read.csv(paste0(github_dir, "data/spp_df.csv"), 
                  check.names=F, header = T, row.names = 'modelno')

##################################################
####   Constants
##################################################
## Index years that had data
Year_Set = seq(min(fit$data_frame[,'t_i']),
               max(fit$data_frame[,'t_i']))
Years2Include = which( Year_Set %in% sort(unique(fit$data_frame[,'t_i'])))
NTime = length(Years2Include)

#Number of sampling grids
N = nrow(Extrapolation_depths)

#Species names
which_spp_idx = unlist(spp_df[VAST_model,])
sci_names = sort( names(spp_df)[which_spp_idx] )

# common_names = c('arrowtooth flounder', 'walleye pollock', 'Pacific cod',
#                  'rex sole', 'flathead sole', 'Pacific halibut', 
#                  'southern rock sole', 'northern rock sole', 'yellowfin sole',
#                  'Dover sole', 'Pacific ocean perch', 
#                  'blackspotted/rougheye\nrockfishes',
#                  'northern rockfish', 'dusky rockfish', 
#                  'shortspine thornyhead')
ns = length(sci_names)

#Sample sizes
samples = c(280, 550, 820)
nboats = 3

#Number of times to simulate survey
Niters = 1000

##################################################
####   Create the data inputs to SamplingStrata
##################################################
df = df_raw = NULL

##################################################
####   Mean density across years
##################################################
df = cbind(
  data.frame(Domain = 1,
             x = 1:N,
             lat = Extrapolation_depths$N_km,
             lon = Extrapolation_depths$E_km - min(Extrapolation_depths$E_km),
             depth = Extrapolation_depths$depth),
  
  #Mean Density across years
  apply(X=fit$Report$D_gcy[,,Years2Include], MARGIN = 1:2, FUN = mean ))

names(df)[-(1:5)] = gsub(x = sci_names, pattern = ' ', replacement = '_')
frame <- SamplingStrata::buildFrameDF(df = df,
                                      id = "x",
                                      X = c("depth", 'lon'),
                                      Y = gsub(x = sci_names, 
                                               pattern = ' ', 
                                               replacement = '_'),
                                      domainvalue = "Domain")

##################################################
####   Predicted density for each observed year and cell
##################################################
for(iT in 1:NTime){
  df_raw = rbind(df_raw, cbind(
    data.frame(Domain = 1,
               x = 1:N,
               year = iT,
               lat = Extrapolation_depths$N_km,
               lon = Extrapolation_depths$E_km-min(Extrapolation_depths$E_km),
               depth = Extrapolation_depths$depth),
    fit$Report$D_gcy[,,Years2Include[iT]] )
  )
}
names(df_raw)[-(1:6)] = gsub(x = sci_names, pattern = ' ', replacement = '_')
frame_raw <- SamplingStrata::buildFrameDF(df = df_raw,
                                          id = "x",
                                          X = c("depth", 'lon'),
                                          Y = gsub(x = sci_names, 
                                                   pattern = ' ', 
                                                   replacement = '_'),
                                          domainvalue = "Domain")

##################################################
####   Calculate "true" mean density
##################################################
frame_raw$year = rep(1:NTime, each = N)
stmt = paste0('aggregate(cbind(',
              paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
              ") ~ year, data = frame_raw, FUN = mean)")
true_mean = eval(parse(text = stmt))[,-1]
colnames(true_mean) = sci_names

##################################################
####   Save Data
##################################################
save(list = c('frame', 'frame_raw', 'true_mean', 'ns', 'NTime', 'N',
              'sci_names', 'samples', 'nboats', 'Niters'),
     file = paste0(github_dir, 'model_', VAST_model, 
                   '/optimization_data.RData'))
