#################################
## Functions in the Optimization Process
#################################

calc_portfolio = function(weights = rep(1/ns, ns),
                          df = data,
                          
                          decision_colname = 'STRATUM',
                          decision_labels = strata
){
  
  ##Result object: mean and variance criteria for each stratum/station
  ##given a vector of species weights
  port_ret = matrix(nrow = length(decision_labels), ncol = 3,
                    dimnames = list(decision_labels, 
                                    c('return', 'TotalVar', 'TotalCov')))
  
  for(id in decision_labels ){ #for each stratum/station
    
    #Subset observation at that stratum/station
    sub_df = subset(df, subset = get(decision_colname) == id)
    
    #Calculate the mean CPUE across years and species, then spread such that
    #each species is a column
    temp = spread(data = aggregate(WGTCPUE ~ YEAR + SPECIES_CODE,
                                   data = sub_df,
                                   FUN = mean), 
                  key = SPECIES_CODE, value = WGTCPUE)
    
    #Scale means years for a species by its respective observed maximum across
    #the time series
    max_means = apply(X = temp[,-1],  MARGIN = 2, FUN = max, na.rm = T)
    scaled_means = sweep(x = temp[,-1], MARGIN = 2, STATS = max_means, FUN='/')
    
    #Calculate the sd for a species on the scaled means
    scaled_sd = apply(X = scaled_means, MARGIN = 2, FUN = sd, na.rm = T)
    
    #Calculate the correlation among species means
    temp_cor = suppressWarnings(cor(x = temp[,-1]))
    
    #There are some strata/stations where certain species are not observed.
    #Set these pairwise correlations to zero. 
    temp_cor[is.na(temp_cor)] = 0
    
    #Weighted covariance among species
    port_var = matrix(nrow = ns, ncol = ns)
    for(i in 1:ns){
      for(j in 1:ns){
        port_var[i,j] = weights[i] * weights[j] * scaled_sd[i] * scaled_sd[j] * temp_cor[i,j]
      }
    }
    
    port_ret[id, ] = c('return' = sum(weights * colMeans(scaled_means), na.rm = T), 
                       'TotalVar' = sqrt(sum(port_var, na.rm = T)),
                       'TotalCov' = sum(port_var[upper.tri(port_var)], na.rm = T)*2 )
  }
  
  return(port_ret)
}

do_optim = function(objvals = optim_df$return,
                    variances = optim_df$TotalVar,
                    number_of_stations = 100,
                    var_constraint = 0.1,
                    covar_constraint = 0.1){
  
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
  # There is 1 constraint for the total variance,
  # 1 constratint for the number of stations
  
  glpkAPI::addRowsGLPK(lp = model, 
                       nrows = 2)
  
  mr = rep(1:2, each = length(objvals))
  mc = rep(1:length(objvals), times = 2)
  mz = c(variances, rep(1, length(objvals)) )
  
  # set non-zero elements of constraint matrix
  glpkAPI::loadMatrixGLPK(lp = model, 
                          ne = length(mz),
                          ia = mr, ja = mc, ra = mz)
  
  # Set the lower and upper bounds for the right-hand side
  # of the inequalities
  lower = c(sum(variances)*var_constraint,
            number_of_stations)
  
  upper = c(Inf,
            number_of_stations)
  
  # Specify the type of structural variable. Integers refer to different types.   
  # See ?glpkAPI::glpkConstants() for the set of variable types. 
  # "2" refers to a variable with a lower bound and 
  # "3" refers a variable with an upper bound
  # "5" refers to a fixed variable
  
  bound_type = c(2, #total variance has a lower bound
                 5) #Number of stations/strata is a fixed variable
  
  glpkAPI::setRowsBndsGLPK(lp = model, 
                           i = seq_along(upper),
                           lb = lower, ub = upper,
                           type = bound_type )
  
  # Presolve and automatically calculate relaxed solution
  # otherwise glpkAPI::solveSimplexGLPK(model) must be called first
  glpkAPI::setMIPParmGLPK(parm = PRESOLVE, val = GLP_ON)
  glpkAPI::setMIPParmGLPK(parm = MSG_LEV, val = GLP_MSG_ALL)
  
  # Set the maximum optimality gap of the solution
  # Set as an argument to generalize function
  glpkAPI::setMIPParmGLPK(parm = MIP_GAP, val = 0.0001)
  
  # Stop after specified number of seconds, convert to milliseconds
  # Set as an argument to generalize function
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