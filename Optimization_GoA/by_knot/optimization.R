
setwd('C:/Users/Zack Oyafuso/Google Drive/VAST_Runs/VAST_output2a')

library(VAST); library(mvtnorm);library(glpkAPI)

# setwd('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/')
load(paste0('VAST_MS_GoA_Run.RData'))
load(paste0('Spatial_Settings.RData'))

Opt = Save$Opt
Report = Save$Report
TmbData = Save$TmbData
Obj = Save$Obj
ns = length(Save$Spp)

Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))
weights = rep(1/ns, ns)
eps = Report$Epsilon2_gct

mean_crit = apply(X=log(Report$D_gcy), MARGIN=1:2, FUN=mean)
mean_crit_scale = apply(mean_crit, MARGIN = 2, FUN = scale)



weight_scen = rbind(matrix(nrow=ns,ncol=ns,data=1/23)+
                      diag(9/23,nrow=ns,ncol=ns),
                    c(rep(1/ns, ns)))

res_df = res_mat = data.frame()

for(ispp in 1:(ns+1)){
  
  optim_df = calc_portfolio(mean_crit_scale = mean_crit_scale,
                            eps = Report$Epsilon2_gct,
                            weights = weight_scen[ispp,],
                            n_x = n_x)
  
  for(i in 1:n_x) {
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
        temp = x$rel_var + 0.01
      }
    }
  }
}

res_mat = as.matrix(res_mat)

save(list = c('res_df', 'res_mat', 'ns', 'n_x'),
     file = paste0('C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/',
                   'Optimization_GoA/by_knot/optimization_results.RData'))




