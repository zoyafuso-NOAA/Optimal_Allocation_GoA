#################################
## Weird small strata
#############################
rm(list = ls())

library(VAST)
library(sp); library(raster); library(RColorBrewer)

############################
## Set up directories
#############################
which_machine = c('Zack_MAC'=1, 'Zack_PC' =2, 'Zack_GI_PC'=3, 'VM' = 4)[3]

VAST_model = "6g"
VAST_dir = paste0(c('/Users/zackoyafuso/Google Drive/', 
                    'C:/Users/Zack Oyafuso/Google Drive/',
                    'C:/Users/zack.oyafuso/Desktop/')[which_machine],
                  'VAST_Runs/VAST_output', VAST_model, '/')

github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/',
                      'C:/Users/zack.oyafuso/Work/')[which_machine],
                    'GitHub/Optimal_Allocation_GoA/')


load(paste0(github_dir, 'data/Extrapolation_depths.RData' ))
load(paste0(github_dir, 'Spatiotemporal_Optimization_Scheme2/spatiotemporal_Flexible_optimization_results.RData'))
load(paste0(github_dir, 'data/optimization_data.RData'))

settings$id = 1:nrow(settings)

########################################3
## Calculate Spatiotemporal Variance of each stratum
########################################
istrata = 10
isample = 550

sub_settings = subset(settings, nstrata == istrata)

isol = sub_settings$id[which.min(abs(sub_settings$n - isample))]

stratano = rep(res_df[,1+isol], NTime)
table(stratano)

#Calculate Stratum Mean Density and Variance
stmt = paste0('aggregate(cbind(',
              paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
              ") ~ stratano, data = frame_raw, FUN = var)")
sample_var = eval(parse(text = stmt))[,-1]

strata_list[isol]
goa = SpatialPointsDataFrame(coords = Extrapolation_depths[,c('E_km', 'N_km')],
                             data = data.frame(stratum = res_df[,isol+1]) )
goa_ras = raster(goa, resolution = 5)
goa_ras =rasterize(x = goa, y = goa_ras, field = 'stratum')

plot(goa_ras, col =  brewer.pal(n = istrata, name = 'Paired'))
values(goa_ras) = ifelse(values(goa_ras) %in% 4:7, 1, NA) 
plot(goa_ras, col = 'black', add = T, legend = F)
