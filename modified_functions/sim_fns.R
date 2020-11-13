##################################
## Create function to calcualte a STRS and output mean, CV, and relative bias
##################################
do_STRS <- function(input){
  
  #Some constants
  n_cells <- dim(input$density)[1]
  n_spp <-  dim(input$density)[2]
  n_time <- dim(input$density)[3]
  
  survey_detail <- 
    data.frame("Stratum" = as.integer(names(table(input$solution))),
               "Nh" = as.integer(table(input$solution)),
               "nh" = input$allocation)
  
  #Assume stratum weights include untrawlabe areas
  survey_detail$Wh <- survey_detail$Nh / n_cells
  survey_detail$wh <- with(survey_detail, nh/Nh)
  
  #Take strata with 0 effor allocation out
  survey_detail <- survey_detail[survey_detail$nh > 0, ]
  
  strata_to_use <- survey_detail$Stratum
  
  #Result objects
  mean_density <- cv <- rel_bias <- matrix(ncol = n_spp, 
                                           nrow = n_time)
  
  for (iyear in 1:n_time) {
    
    #Subset density df by year
    sub_df <- input$density[, , iyear]
    
    #Take a random sample based on the allocation and stratum
    sample_vec <- c()
    for(istrata in 1:nrow(survey_detail)) {
      sample_vec <- c(sample_vec,
                      sample(x = which(input$solution == survey_detail$Stratum[istrata]),
                             size = survey_detail$nh[istrata]))
    }
    
    sampled_strata <- rep(x = survey_detail$Stratum, 
                          times = survey_detail$nh)
    
    #subset sub_df by which cells were chosen
    sample_df <- as.data.frame(sub_df[sample_vec, ])
    colnames(sample_df) <- paste0("Y", 1:ns)
    
    #Calculate STRS mean density
    stmt <- paste0('aggregate(cbind(',
                   paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                   ") ~ sampled_strata, data = sample_df, FUN = mean)")
    sample_mean <- eval(parse(text = stmt))[, -1]
    STRS_mean <- colSums(sweep(x = sample_mean, 
                               MARGIN = 1, 
                               STATS = with(survey_detail, Wh),
                               FUN = '*'))
    
    #Calculate STRS variance of mean density
    stmt <- paste0('aggregate(cbind(',
                   paste0('Y', 1:(ns-1), sep = ',', collapse = ''), 'Y',ns, 
                   ") ~ sampled_strata, data = sample_df, FUN = var)")
    sample_var <- eval(parse(text = stmt))[, -1]
    STRS_var <- colSums(sweep(x = sample_var, 
                              MARGIN = 1, 
                              STATS = with(survey_detail, Wh^2 * (1 - wh) / nh),
                              FUN = '*'))
    
    #Save mean and cv of estimates across species
    mean_density[iyear, ] <- STRS_mean
    cv[iyear, ] <- sqrt(STRS_var) / STRS_mean 
    
    #Calculate relative bias of mean estimate
    rel_bias[iyear, ] <- unlist(
      100 * (STRS_mean - input$true_density[iyear, ]) / 
        input$true_density[iyear, ])
    
  }
  
  return(list("mean_denisty" = round(mean_density, 2),
              "cv" = round(cv, 4),
              "rel_bias" = round(rel_bias, 2) ))
  
}
