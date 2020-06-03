##############################
## Bethel Algorithm on GoA current strata
##############################
rm(list = ls())

############################
## Import Libraries
############################
library(rgdal); library(raster); library(rgeos); library(tidyr)
library(SamplingStrata)

############################
## Set up directories
#############################
which_machine = c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI' = 3)[3]
modelno = '6g'
github_dir = paste0(c('/Users/zackoyafuso/Documents/', 
                      'C:/Users/Zack Oyafuso/Documents/',
                      'C:/Users/zack.oyafuso/Work/')[which_machine], 'GitHub/MS_OM_GoA/Optimum_Allocation_AI/')

############################
## Load Data
############################
# load(paste0(VAST_dir, '/Spatial_Settings.RData'))
# load(paste0(github_dir, '/Optimum_Allocation/model_', modelno, 
#             '/optimization_data_model_', modelno, '.RData'))

#Survey data
survey_data = readRDS(paste0(github_dir, "AI_BTS.rds"))
survey_data = subset(survey_data,
                     SPECIES_NAME %in% c('Pleurogrammus monopterygius',
                                         'Sebastes alutus', 
                                         'Sebastes aleutianus'),
                     select = c(YEAR, STRATUM, CPUE, SPECIES_NAME))

strata_allo = with(subset(survey_data, 
                          SPECIES_NAME == 'Sebastes alutus'),
                   table(YEAR, STRATUM))

#Grid data
survey_grid = read.csv(paste0(github_dir, "grid_AI.csv"))
stratapop = table(survey_grid$STRATUM)[colnames(strata_allo)]

## Constants
sci_names = sort(unique(survey_data$SPECIES_NAME))
ns = length(sci_names)
strata = colnames(strata_allo)
years = sort(unique(survey_data$YEAR))
NTime = length(years)
N = sum(stratapop)

######################
## Results Objects
######################

CVs = seq(0.20, 0.40, by = 0.005)
sample_allocation = array(data = 0, dim = c(length(CVs), length(strata)),
                          dimnames = list(CVs, strata))

spp_cvs = array(data = 0, dim = c(length(CVs), ns),
                dimnames = list(CVs, sci_names) )

##########################
## Create Input for bethel algorithm
## Stratum Means and Variances
##########################
sample_mean = spread(data = aggregate(CPUE ~ SPECIES_NAME + STRATUM,
                                      data = survey_data, FUN=mean, drop = F),
                     key = SPECIES_NAME, value = CPUE)[,-1]
names(sample_mean) = paste0('M', 1:ns)

sample_var = spread(data = aggregate(CPUE ~ SPECIES_NAME + STRATUM,
                                     data = survey_data, FUN = var, drop = F),
                    key = SPECIES_NAME, value = CPUE)[,-1]
names(sample_var) = paste0('S', 1:ns)

sample_mean[is.na(sample_var)] = 0
sample_var[is.na(sample_var)] = 0

df = cbind(data.frame(stratum = strata,
                      N = as.vector(stratapop)),#,
           #X1 = factor(1:length(temp_strata))),
           sample_mean, sqrt(sample_var),
           data.frame(cens = 0,
                      cost = 1,
                      DOM1 = 'tot'))

##################################
## Calculate optimal allocation across CV levels
##################################
for(icv in CVs){
 stmt = paste0('cbind(', paste0("CV", 1:ns, '=', icv, collapse = ', '), ')' )
 CV = eval(parse(text = stmt))
 errors = cbind(data.frame(DOM = 'DOM1'), CV, domainvalue = 1)
 
 n = bethel(stratif = df, errors = errors, printa=TRUE, 
            epsilon = 1e-11, maxiter = 200)
 sample_allocation[paste(icv),] = as.numeric(n)
 spp_cvs[paste(icv),] = as.numeric( attributes(n)$outcv[,'ACTUAL CV'] )
}

total_sample_size = apply(sample_allocation, MARGIN = 1, sum)

###############################
## Plots
###############################

{png(filename = paste0(github_dir, 'AI_tradeoff_ST.png'),
     width = 8, height = 5, units = 'in', res = 500)
 par(mfrow = c(1,1), mar = c(5,5,1,1))
 plot(1, pch = 16, las = 1, xlim = range(CVs), ylim = c(0, 800), 
      type = 'n', xlab = 'Upper CV Constraint', ylab = 'Total Sample Size')
 lines(CVs, total_sample_size)
 points(CVs, total_sample_size, pch = 16)
 
 abline(h = c(280, 550, 820), lty = 'dotted')
 text(x = .59, y = c(350, 620, 880), c('1 Boat', paste(2:3, 'Boats')))
 dev.off()}



{png(filename = paste0(github_dir,
                       '/AI_sppCV.png'),
     width = 8, height = 7, units = 'in', res = 500)
par(mfrow = c(1,1), mar = c(5,5,1,1))
matplot(CVs, spp_cvs, las = 1, ylab = 'Spatiotemporal CV',
        xlab = 'Upper Spatiotempral CV constraint',
        pch = 16)
legend('topleft', legend = sci_names, pch = 16, col = palette(), bty = 'n', text.font = 3)
abline(a = 0, b = 1)
dev.off()}

###############################
#
#Calculate Total Mean and Variance, calculate CV
sample_mean = aggregate(CPUE ~ SPECIES_NAME + STRATUM + YEAR,
                                      data = survey_data, FUN=mean, drop = F)
sample_mean = split(sample_mean, f = sample_mean$YEAR)

sample_var = aggregate(CPUE ~ SPECIES_NAME + STRATUM + YEAR,
                                     data = survey_data, FUN = var, drop = F)
sample_var = split(sample_var, f = sample_var$YEAR)

observed_meancv = array(dim = c(NTime, ns, 2),
                        dimnames = list(NULL, sci_names, c('mean', 'cv')))
for(iyear in 1:NTime){
   temp_sample_mean = spread(data = sample_mean[[iyear]], key = SPECIES_NAME, value = CPUE, fill = 0)[,-c(1:2)]
   temp_sample_var =  spread(data = sample_var[[iyear]], key = SPECIES_NAME, value = CPUE, fill = 0)[,-c(1:2)]
   temp_strata_allocation = strata_allo[iyear,strata]
   
   SRS_mean = colSums(sweep(x = temp_sample_mean, MARGIN = 1, 
                            STATS = stratapop / N,
                            FUN = '*'))
   
   SRS_var = colSums(sweep(x = temp_sample_var, MARGIN = 1, 
                           STATS = (stratapop/N)^2*(1 - temp_strata_allocation/stratapop)/temp_strata_allocation,
                           FUN = '*'))
   
   strata_CV = sqrt(SRS_var) / SRS_mean
   
   observed_meancv[iyear, , 'mean'] = SRS_mean
   observed_meancv[iyear, , 'cv'] = strata_CV
}

# par(mfrow = c(1,1), mar = c(5,5,1,1))
# plot(1, pch = 16, las = 1, xlim = c(0.1,0.6), ylim = c(0, 800), 
#      type = 'n', xlab = 'Upper CV Constraint', ylab = 'Total Sample Size')
# lines(CVs, total_sample_size)
# points(CVs, total_sample_size, pch = 16)
# 
# points(y = rowSums(strata_allo ), x = observed_meancv[,1,'cv'], pch = 16, col = 'red')
# points(y = rowSums(strata_allo ), x = observed_meancv[,2,'cv'], pch = 16, col = 'blue')
# points(y = rowSums(strata_allo ), x = observed_meancv[,3,'cv'], pch = 16, col = 'darkgreen')

