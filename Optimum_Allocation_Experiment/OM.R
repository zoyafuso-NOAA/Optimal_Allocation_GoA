######################################
## Operating Model
## Create a domain with four examples of spatiotemporal variability
## 1. No spatial and spatiotemporal variability
## 2. Just Spatial Variability
## 3. Just Spatiotemporal Variability
## 4. Both Spatial and Spatiotemporal Variability
######################################
rm(list = ls())
which_machine = c('Zack_Mac' = 1)[1]
output_wd = paste0(c('/Users/zackoyafuso/Documents/')[which_machine],
                   'GitHub/MS_OM_GoA/Optimum_Allocation_Experiment/')

######################################
## Import Libaries
## Import OM function
######################################
library(RandomFields)

Sim_Fn = function( logmean=1, 
                   Scale=0.2, 
                   SD_omega=1, 
                   SD_epsilon=1, 
                   x_dim=50, 
                   y_dim=50,
                   n_years=10 ){
  loc_xy = expand.grid( "x"=seq(0,1,length=x_dim), 
                        "y"=seq(0,1,length=y_dim))
  
  RF_omega = RMgauss(var=SD_omega^2, scale=Scale)
  RF_epsilon = RMgauss(var=SD_epsilon^2, scale=Scale)
  
  log_d_rt = Epsilon_rt = matrix(NA, ncol=n_years, nrow=x_dim*y_dim)
  Omega_r = RFsimulate(model=RF_omega, 
                       x=loc_xy[,'x'], 
                       y=loc_xy[,'y'])@data[,1]
  for(t in 1:n_years){
    Epsilon_rt[,t] = RFsimulate(model=RF_epsilon, 
                                x=loc_xy[,'x'], 
                                y=loc_xy[,'y'])@data[,1]
    log_d_rt[,t] = logmean + Epsilon_rt[,t] + Omega_r
  }
  
  
  # Return stuff
  Return = list("log_d_rt"=log_d_rt, 
                "Epsilon_rt"=Epsilon_rt, 
                "Omega_r"=Omega_r, 
                "loc_xy"=loc_xy, 
                "n_per_year"=x_dim*y_dim, 
                "n_years"=n_years)
  return( Return )
}

###############################
## Settings for OM Scenarios
###############################
ST_settings = expand.grid(Scale = c(0.3),
                          SD_omega=c(0.5,0.1), 
                          SD_epsilon=c(0.5,0.1) )

ST_settings = ST_settings[-2,]

###############################
#Calculate Operation Models
###############################
OM = list()
for(irow in 1:nrow(ST_settings)){
  set.seed(100)
  OM[[irow]] = Sim_Fn(logmean=1, 
                      Scale = ST_settings$Scale[irow], 
                      SD_omega=ST_settings$SD_omega[irow], 
                      SD_epsilon=ST_settings$SD_epsilon[irow], 
                      x_dim=100, 
                      y_dim=100,
                      n_years=10)
}

###############################
## Save
###############################
save(list = c('OM', 'ST_settings'), file = paste0(output_wd, 'OM.RData') )
