##########################
## Preliminary MS GoA 
## Survey Optimizations
## Zack Oyafuso 
## Created 9 February 2020
##########################

setwd("/Users/zackoyafuso/Documents/GitHub/MS_OM_GoA/")

library(tidyr)
library(glpkAPI)

#############################
## Load Optimization Functions
#############################
source('Optimization_GoA/optimization_functions.R')

# or from flat files exported from AFSC database
GOA = read.csv(paste0("/Users/zackoyafuso/Desktop/AK_BTS/data-raw/",
                      "cpue_GOA_selected_spp.csv"), 
               stringsAsFactors = FALSE) # CPUE is (num or kg / km^2)

# Filter Species: 
# Arrowtooth Flounder (Atherestes stomias, code 10110)
# Pacific Cod (Gadus macrocephalus, code 21720)
# Pacific Ocean Perch (Sebastes alutus, code 30060)
# Sablefish (Anoplopoma fimbria, code 20510)
# Walleye pollock (Gadus chalcogrammus, code 21740)
# Dover sole (Solea solea, code 10180)
# Pacific halibut (Hippoglossus stenolepis, code 10120)
# Flathead sole (Hippoglossoides elassodon, code 10130)
# Rex sole (Glyptocephalus zachirus, code 10200)
# Dusky rockfish (Sebastes variabilis, code 30152)
# Northern rockfish (Sebastes polyspinis, code 30420)
# Rougheye and blackspotted rockfishes (Sebastes aleutianus and Sebastes melanostictus, respectively, codes 30050,30051,30052)
# Northern and Southern rock sole (Lepidopsetta polyxystra and Lepidopseta bilineata, respectivity, codes 10260,10261,10262)

# data <- filter(GOA, SPECIES_CODE %in% c(10110, 21720, 30060, 20510, 21740, 
#                                         10180, 10120, 10130, 10200, 30152,
#                                         30420, 30050,30051,30052 ))

species = c(10110, 21720, 30060, 20510, 21740, 
            10180, 10120, 10130, 10200, 30152,
            30420, 30050,30051,30052 )
data = GOA[GOA$SPECIES_CODE %in% species,]


ns = length(species)
nstrata = length(unique(data$STRATUM))
strata = paste(sort(unique(data$STRATUM)) )

weight_scen = rbind(matrix(nrow=ns,ncol=ns,data=1/23)+
                      diag(9/23,nrow=ns,ncol=ns),
                    c(rep(1/ns, ns)))

res_df = res_mat = data.frame()

for(ispp in 1:(ns+1)){
  
  optim_df = calc_portfolio(weights = weight_scen[ispp,])

  for(i in 1:nstrata) {
    temp = 0.001; output_code = 0
    
    while(output_code %in% c(0,14)){
      x = do_optim(objvals = optim_df[, 'return'],
                   variances = optim_df[, 'TotalVar'],
                   number_of_stations = i,
                   var_constraint = temp)
      output_code = x$output_code
      
      if(output_code %in% c(0,14)){
        res_df = rbind(res_df, data.frame(spp_scen = ispp,
                                          n = i,
                                          tot_var = x$tot_var,
                                          rel_var = x$rel_var,
                                          tot_mean = x$objval) )
        
        res_mat = rbind(res_mat, as.integer(x$x))
        temp = x$rel_var + 0.001
      }
    }
  }
}

res_mat = as.matrix(res_mat)

save(list = c('res_df', 'res_mat', 'ns', 'nstrata', 'strata', 'species'),
     file = 'Optimization_GoA/optimization_results.RData')


