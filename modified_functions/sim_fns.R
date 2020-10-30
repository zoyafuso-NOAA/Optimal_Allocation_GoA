##################################
## Create function to calcualte a SRS and output mean, CV, and relative bias
##################################
do_SRS <- function(density = frame_raw[, paste0("Y", 1:ns)],
                   true_density = true_mean,
                   time = frame_raw$year,
                   n = c(280, 550, 820)[1],
                   cell_idx = 1:N){
  
  #Some constants extracted from the data
  n_time <- length(unique(time))
  n_spp <-  dim(density)[2]
  
  #Result objects
  mean_density <- cv <- rel_bias <- matrix(ncol = n_spp, nrow = n_time)
  
  for(iyear in 1:n_time) {
    # Take n random samples of the available cells
    sample_vec <- sample(x = cell_idx, 
                         size = n)
    
    #Subset the density df based on which cells were sampled and year
    sample_df <- density[time == iyear, ][sample_vec, ]
    
    #Calculate Mean and standard error
    temp_sim_mean <- colMeans(sample_df)
    temp_sim_var <- apply(sample_df, 
                          MARGIN = 2, 
                          FUN = var)
    temp_sim_se <- sqrt(temp_sim_var / n)
    
    #Save mean and cv of estimates across species
    mean_density[iyear, ] <- temp_sim_mean
    cv[iyear, ] <- temp_sim_se / temp_sim_mean
    
    #Calculate relative bias of mean estimate
    rel_bias[iyear, ] <- unlist(100 * (temp_sim_mean - true_density[iyear, ]) / 
                                  true_density[iyear, ])
  }
  
  return(list("mean_denisty" = round(mean_density, 2),
              "cv" = round(cv, 4),
              "rel_bias" = round(rel_bias, 2) ))
  
}

##################################
## Create function to calcualte a STRS and output mean, CV, and relative bias
##################################
do_STRS <- function(density = frame_raw[, paste0("Y", 1:ns)],
                    cell_idx = (1:N),
                    strata = Extrapolation_depths$stratum,
                    strata_to_use = allocations$Stratum, 
                    allocation = allocations$boat1,
                    true_density = true_mean,
                    time = frame_raw$year){
  
  #Some constants
  n_time <- length(unique(time))
  n_spp <-  dim(density)[2]
  
  n_cells <- length(cell_idx)
  Nh <- table(strata)[allocation > 0]
  Wh <- Nh / n_cells
  nh <- allocation[allocation > 0]
  wh <- nh / Nh
  
  strata_to_use <- strata_to_use[allocation > 0]

  #Result objects
  mean_density <- cv <- rel_bias <- matrix(ncol = n_spp, nrow = n_time)
  
  for (iyear in 1:n_time) {
    #Subset density df by year, then by available cells
    sub_df <- density[time == iyear, ]
    
    #Take a random sample based on the allocation and stratum
    sample_vec <- c()
    for(istrata in 1:length(strata_to_use)) {
      sample_vec <- c(sample_vec,
                      sample(x = which(strata == strata_to_use[istrata]),
                             size = nh[istrata]))
    }
    
    sampled_strata <- rep(x = strata_to_use, 
                          times = nh)
    
    #subset sub_df by which cells were chosen
    sample_df <- sub_df[sample_vec, ]
    
    #Calculate STRS mean density
    stmt <- paste0('aggregate(cbind(',
                   paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                   ") ~ sampled_strata, data = sample_df, FUN = mean)")
    sample_mean <- eval(parse(text = stmt))[, -1]
    STRS_mean <- colSums(sweep(x = sample_mean, 
                               MARGIN = 1, 
                               STATS = Wh,
                               FUN = '*'))
    
    #Calculate STRS variance
    stmt <- paste0('aggregate(cbind(',
                   paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                   ") ~ sampled_strata, data = sample_df, FUN = var)")
    sample_var <- eval(parse(text = stmt))[, -1]
    STRS_var <- colSums(sweep(x = sample_var, 
                              MARGIN = 1, 
                              STATS = Wh^2 * (1 - wh) / nh,
                              FUN = '*'))
    
    #Save mean and cv of estimates across species
    mean_density[iyear, ] <- STRS_mean
    cv[iyear, ] <- sqrt(STRS_var) / STRS_mean 
    
    #Calculate relative bias of mean estimate
    rel_bias[iyear, ] <- unlist(100 * (STRS_mean - true_density[iyear, ]) / 
                                  true_density[iyear, ])
  }
  
  return(list("mean_denisty" = round(mean_density, 2),
              "cv" = round(cv, 4),
              "rel_bias" = round(rel_bias, 2) ))
  
}

##################################
## Create function to calcualte a STRS for a single species 
## and output mean, CV, and relative bias
##################################
do_SS_STRS <- function(density = frame_raw[, paste0("Y15")],
                       cell_idx = (1:N)[Extrapolation_depths$stratum != 0],
                       strata = Extrapolation_depths$stratum[Extrapolation_depths$stratum != 0],
                       strata_to_use = allocations$Stratum, 
                       allocation = allocations$boat1,
                       true_density = true_mean[, 15],
                       time = frame_raw$year){
  
  #Do some input checks
  if(! length(cell_idx) == length(strata)) 
    stop("cell_idx must be same length as strata")
  
  if(! length(allocation) == length(strata_to_use)) 
    stop("cell_idx must be same length as strata")
  
  #Some constants
  n_time <- length(unique(time))
  n_spp <-  dim(density)[2]
  
  n_cells <- length(cell_idx)
  Nh <- table(strata)[allocation > 0]
  Wh <- Nh / n_cells
  nh <- allocation[allocation > 0]
  wh <- nh / Nh
  
  #Result objects
  mean_density <- cv <- rel_bias <- vector(length = n_time)
  
  for (iyear in 1:n_time) {
    #Subset density df by year, then by available cells
    sub_df <- density[time == iyear][cell_idx]
    
    #Take a random sample based on the allocation and stratum
    sample_vec <- c()
    for(istrata in 1:length(strata_to_use)) {
      sample_vec <- c(sample_vec,
                      sample(x = which(strata == strata_to_use[istrata]),
                             size = allocation[istrata]))
    }
    
    sampled_strata <- rep(x = strata_to_use, 
                          times = allocation)
    
    #subset sub_df by which cells were chosen
    sample_df <- sub_df[sample_vec]
    
    #Calculate STRS mean density
    sample_mean <- tapply(X = sample_df, 
                          INDEX = sampled_strata, 
                          FUN = mean)
    
    STRS_mean <- sum(sample_mean * Wh)
    
    #Calculate STRS variance
    sample_var <- tapply(X = sample_df, 
                         INDEX = sampled_strata, 
                         FUN = var)
    STRS_var <- sum( sample_var * Wh^2 * (1 - wh) / nh )
    
    #Save mean and cv of estimates across species
    mean_density[iyear] <- STRS_mean
    cv[iyear] <- sqrt(STRS_var) / STRS_mean 
    
    #Calculate relative bias of mean estimate
    rel_bias[iyear] <- unlist(100 * (STRS_mean - true_density[iyear]) / 
                                true_density[iyear])
  }
  
  return(list("mean_denisty" = round(mean_density, 2),
              "cv" = round(cv, 4),
              "rel_bias" = round(rel_bias, 2) ))
  
}
