calc_portfolio = function(mean_crit_scale = mean_crit_scale,
                          eps = Report$Epsilon2_gct,
                          weights = weights,
                          n_x = n_x){
  #Weighted Sum
  wgt_mean_criteria = rowSums( sweep(x = mean_crit_scale, 
                                     MARGIN = 2, 
                                     STATS = weights, 
                                     FUN = '*') 
  )
  
  #Weighted Variance
  wgt_total_variance = vector(length = n_x)
  
  for(iknots in 1:n_x){
    temp = eps[iknots,,Years2Include] 
    temp_sds = diag( x = weights * apply(X = temp, MARGIN = 1, sd) )
    temp_cor = cor(t(temp))
    
    wgt_total_variance[iknots] = sum(temp_sds%*% temp_cor %*% temp_sds)
  }
  
  return(data.frame(return = wgt_mean_criteria,
                    TotalVar = wgt_total_variance))
}

do_optim = function(objvals = optim_df$return,
                    variances = optim_df$TotalVar,
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
