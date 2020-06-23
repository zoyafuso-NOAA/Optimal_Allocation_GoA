#################################
## Parallelize the Optimization Process
## Method == "continous"
#################################
rm(list = ls())
# library(SamplingStrata);
library(memoise); library(doParallel); library(foreach); library(iterators); library(parallel); library(pbapply); library(formattable)

for(ifile in dir('C:/Users/Zack Oyafuso/Downloads/SamplingStrata-master/R', full.names = T))
  source(ifile)
source('C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/Optimum_Allocation/buildStrataDF_Zack.R')

library(VAST); 
library(mvtnorm); 
# 
library(sp)
library(RColorBrewer); library(raster)

which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3)[2]

VAST_wd = c('/Users/zackoyafuso/Google Drive/VAST_Runs/',
            'C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/',
            'C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/')[which_machine]

VAST_model = "6d"
output_wd = c(paste0('/Users/zackoyafuso/Documents/GitHub/MS_OM_GoA/',
                     'Optimum_Allocation/model_', VAST_model),
              paste0("C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/",
                     "Optimum_Allocation/model_", VAST_model),
              paste0("C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/",
                     "Optimum_Allocation/model_", VAST_model))[which_machine]

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
  data.frame(Domain = cut(x = Extrapolation_depths$Lon, 
                          breaks = c(-171, -159, -154, -147, -140, -130), 
                          labels = c('Shumagin_1', 'Chirikof_2', 'Kodiak_3',
                                     'Yakutak_4', 'SE_5')),
             x = 1:Save$TmbData$n_g,
             lat = Extrapolation_depths$N_km,
             lon = Extrapolation_depths$E_km - min(Extrapolation_depths$E_km),
             depth = Extrapolation_depths$depth),
  apply(X=Save$Report$Index_gcyl[,,Years2Include,], MARGIN = 1:2, FUN = mean ) )

names(df)[-(1:5)] = gsub(x = Save$Spp, pattern = ' ', replacement = '_')

frame <- buildFrameDF(df = df,
                      id = "x",
                      X = c("depth"),#, 'lon'),
                      Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                      domainvalue = "Domain")

for(iT in 1:NTime){
  df_raw = rbind(df_raw, cbind(
    data.frame(Domain = cut(x = Extrapolation_depths$Lon, 
                            breaks = c(-171, -159, -154, -147, -140, -130), 
                            labels = paste(1:5)),
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
                          X = c("depth"),#, 'lon'),
                          Y = gsub(x = Save$Spp, pattern = ' ', replacement = '_'),
                          domainvalue = "Domain")

#Settings for optimizer
settings = expand.grid(cv = c(0.3),
                       mut_change = c(0.1),
                       elitism_rate = c(0.2),
                       dom = c(2:5))

ns = Save$TmbData$n_c
domains = unique(df$Domain)
ndom_ = length(unique(frame$domainvalue))


rm(list = c('Save', 'Spatial_List', 'spp_df', 'strata.limits', 'fine_scale',
            'Method', 'modelno', 'n_x', 'which_spp', 'Year_Set', 
            'Years2Include', 'Data_Geostat', 'df', 'Extrapolation_List',
            'gulf_of_alaska_grid', 'ifile', 'iT'))

res_df = as.matrix(frame[,c('id', 'domainvalue')])
strata_list = list()

for(ii in 2:nrow(settings)){
  
  par(mfrow = c(3,2))
  plot_this = (ii%%10 == 0)
  
  cv = list()
  for(spp in 1:ns) cv[[paste0('CV', spp)]] = rep(settings$cv[ii], ndom_)
  cv[['DOM']] = levels(domains)
  cv[['domainvalue']] = as.numeric(domains)
  cv <- as.data.frame(cv)
  
  set.seed(1234 + ii)
  solution <- optimStrata(method = "continuous",
                          errors = cv, 
                          framesamp = frame,
                          iter = 50,
                          pops = 10,
                          elitism_rate = settings$elitism_rate[ii],
                          mut_chance = settings$mut_change[ii],
                          # nStrata = rep(5, ndom_),
                          nStrata = rep(settings$dom[ii],ndom_),
                          showPlot = plot_this,
                          parallel = T,
                          cores = 2)
  
  strata_list[[ii]] =  summaryStrata(solution$framenew,
                                    solution$aggr_strata,
                                    progress=FALSE) 
  
  res_df = cbind(res_df, solution$framenew$STRATO)
  
  if(ii %% 1 == 0) {
    save(list = c('strata_list', 'res_df'), 
         file = paste0(output_wd, '/optimization_spatiotemporal.RData'))
  }
}

