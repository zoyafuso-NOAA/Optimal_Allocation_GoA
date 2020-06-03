####################################
## Multispecies Survey Optimization
####################################

setwd('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/Optimization_BeringSea')

#############################
## Import Libraries
#############################
library(glpkAPI)
library(tidyr)
library(reshape2)

#########################################
## Import Data
#########################################
EBS = read.csv('EBS_data_trimmed.csv')

########################################
## Calculate Mean and SD CPUE for each species at each station
########################################
spp_mean = spread(data = aggregate(WGTCPUE ~ SPECIES_CODE + STATIONID, data = EBS, FUN = mean), key = SPECIES_CODE, value = WGTCPUE)

spp_sd = spread(data = aggregate(WGTCPUE ~ SPECIES_CODE + STATIONID, data = EBS, FUN = sd), key = SPECIES_CODE, value = WGTCPUE)

########################################
## Species labels
########################################
spp = length(unique(EBS$SPECIES_CODE))
spp_labels =  c('Arrowtooth_Flounder', 'Flathead_Sole',
                'Yellowfin_Sole', 'AK_Plaice', 'Walleye_Pollock', 'All')

#############################
## Create a function that inputs a vector of species weights
## Calculates the weighted return and total variance across spp. for a station
#############################

calc_portfolio = function(weights = rep(1/5, 5)){
  
  port_ret = matrix(nrow = nrow(spp_mean), ncol = 2,
                    dimnames = list(spp_mean$STATIONID, 
                                    c('return', 'variance')))
  
  for(id in spp_mean$STATIONID){
    temp = spread(data = subset(x = EBS, subset = STATIONID == id), 
                  key = SPECIES_CODE,
                  value = WGTCPUE)
    
    scaled_mean = colMeans(temp[,-c(1:2)]) / apply(X = spp_mean[,-1], 
                                                   MARGIN = 2, 
                                                   FUN = max)
    
    scaled_sd = apply(temp[,-c(1:2)], 
                      MARGIN = 2, 
                      FUN = sd) / apply(X = spp_sd[,-1], 
                                        MARGIN = 2, 
                                        FUN = max)
    
    temp_cor = suppressWarnings(cor(x = temp[,-c(1:2)]))
    temp_cor[is.na(temp_cor)] = 0
    
    port_var = matrix(nrow = spp, ncol = spp)
    for(i in 1:spp){
      for(j in 1:spp){
        port_var[i,j] = weights[i] * weights[j] * scaled_sd[i] * scaled_sd[j] * temp_cor[i,j]
      }
    }
    
    port_ret[id, ] = c('return' = sum(weights * scaled_mean), 
                       'var' = sum(port_var) )
  }
  
  return(port_ret)
}

do_optim = function(objvals = optim_df$return,
                    variances = optim_df$variance,
                    number_of_stations = 100,
                    var_constraint = 0.1){
  
  ##############################
  ## Defining the Model
  ###############################
  
  # Initialize Model
  model <- glpkAPI::initProbGLPK()
  
  # Set objective function as a minimization funtion
  glpkAPI::setObjDirGLPK(lp = model, 
                         lpdir = glpkAPI::GLP_MAX)
  
  # Initialize decision variables (columns)
  glpkAPI::addColsGLPK(lp = model, 
                       ncols = length(objvals))
  
  # Set the objective function, specify no bounds on decision variables
  # GLP_FR means free variable
  glpkAPI::setColsBndsObjCoefsGLPK(lp = model, 
                                   j = seq_along(objvals),
                                   lb = NULL, ub = NULL,
                                   obj_coef = objvals,
                                   type = rep(glpkAPI::GLP_FR, length(objvals)))
  
  # Specify that decision variables are binary. GLP_BV means binary variable
  glpkAPI::setColsKindGLPK(lp = model, 
                           j = seq_along(objvals),
                           kind = rep(glpkAPI::GLP_BV, length(objvals)))
  
  # Initialize the structural constraints (rows)
  # There is 1 constraint for the total variance, another for the number of stations
  glpkAPI::addRowsGLPK(lp = model, 
                       nrows = 2)
  
  mr = rep(1:2, each = length(objvals))
  mc = rep(1:length(objvals), times = 2)
  mz = c(variances, rep(1, length(objvals)))
  
  # set non-zero elements of constraint matrix
  glpkAPI::loadMatrixGLPK(lp = model, 
                          ne = length(mz),
                          ia = mr, ja = mc, ra = mz)
  
  # Set the lower and upper bounds for the right-hand side
  # of the inequalities
  lower = c(sum(variances)*var_constraint, number_of_stations)
  
  upper = c(Inf, number_of_stations)
  
  # Specify the type of structural variable. Integers refer to different types.   
  # See ?glpkAPI::glpkConstants() for the set of variable types. 
  # "2" refers to a variable with a lower bound and 
  # "3" refers a variable with an upper bound
  # "5" refers to a fixed variable
  
  bound_type = c(2, #total variance has a lower bound
                 5) #Fixed variable
  
  glpkAPI::setRowsBndsGLPK(lp = model, 
                           i = seq_along(upper),
                           lb = lower, ub = upper,
                           type = bound_type )
  
  # Presolve and automatically calculate relaxed solution
  # otherwise glpkAPI::solveSimplexGLPK(model) must be called first
  glpkAPI::setMIPParmGLPK(parm = PRESOLVE, val = GLP_ON)
  glpkAPI::setMIPParmGLPK(parm = MSG_LEV, val = GLP_MSG_ALL)
  
  # Set the maximum optimality gap of the solution
  glpkAPI::setMIPParmGLPK(parm = MIP_GAP, val = 0.0001)
  
  # Stop after specified number of seconds, convert to milliseconds
  glpkAPI::setMIPParmGLPK(parm = TM_LIM, val = 1000 * 60) 
  
  #############################
  ## Solve Model
  #############################
  screen_out = glpkAPI::solveMIPGLPK(lp = model) 
  
  #############################
  ## Prepare return object
  #############################
  results <- list(output_code = screen_out,
                  status = glpkAPI::return_codeGLPK(screen_out),
                  objval = glpkAPI::mipObjValGLPK(model),
                  x = glpkAPI::mipColsValGLPK(model) )
  
  results$tot_var = sum(variances[results$x == 1])
  results$rel_var = results$tot_var / sum(variances)
  
  #Delete model
  glpkAPI::delProbGLPK(lp = model)
  
  return(results)
}

weight_scen = rbind(matrix(nrow=spp,ncol=spp,data=1/14)+
                      diag(9/14,nrow=spp,ncol=spp),
                    c(rep(1/spp, spp)))

res_df = res_mat = data.frame()

for(ispp in 1:(spp+1)){

  optim_df = calc_portfolio(weights = weight_scen[ispp,])
  
  # Plot distribution of weighted mean and total variance
  # par(mfrow = c(2,1), mar = c(0,0,0,0))
  # plot(lat ~ long, data = optim_df, cex = return * 20, axes = F); box()
  # plot(lat ~ long, data = optim_df, cex = variance * 100, axes = F); box()
  
  for(i in seq(from=50, to=300, by=25)){
    temp = 0.1; opt_res = 0
    
    while(opt_res %in% c(0, 14)){
      x = do_optim(objvals = optim_df[, 'return'], 
                   variances = optim_df[, 'variance'],
                   number_of_stations = i,
                   var_constraint = temp)
      
      opt_res = x$output_code
      
      if(x$output_code %in% c(0, 14)){
        res_df = rbind(res_df, data.frame(spp_scen = ispp,
                                          n = sum(x$x == 1),
                                          tot_var = x$tot_var,
                                          rel_var = x$rel_var,
                                          tot_mean = x$objval) )
        
        res_mat = rbind(res_mat, as.integer(x$x))
        
        temp = x$rel_var + 0.0001
        
      }
      
    }
  }
  
}

res_mat = as.matrix(res_mat)

save(list = c('res_df', 'res_mat', 'spp', 'spp_labels'),
     file = 'optimization_results.RData')
