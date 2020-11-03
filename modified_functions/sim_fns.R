##################################
## Create function to calcualte a SRS and output mean, CV, and relative bias
##################################
# do_SRS <- function(density = frame_raw[, paste0("Y", 1:ns)],
#                    true_density = true_mean,
#                    time = frame_raw$year,
#                    n = c(280, 550, 820)[1],
#                    cell_idx = 1:N){
#   
#   #Some constants extracted from the data
#   n_time <- length(unique(time))
#   n_spp <-  dim(density)[2]
#   
#   #Result objects
#   mean_density <- cv <- rel_bias <- matrix(ncol = n_spp, nrow = n_time)
#   
#   for(iyear in 1:n_time) {
#     # Take n random samples of the available cells
#     sample_vec <- sample(x = cell_idx, 
#                          size = n)
#     
#     #Subset the density df based on which cells were sampled and year
#     sample_df <- density[time == iyear, ][sample_vec, ]
#     
#     #Calculate Mean and standard error
#     temp_sim_mean <- colMeans(sample_df)
#     temp_sim_var <- apply(sample_df, 
#                           MARGIN = 2, 
#                           FUN = var)
#     temp_sim_se <- sqrt(temp_sim_var / n)
#     
#     #Save mean and cv of estimates across species
#     mean_density[iyear, ] <- temp_sim_mean
#     cv[iyear, ] <- temp_sim_se / temp_sim_mean
#     
#     #Calculate relative bias of mean estimate
#     rel_bias[iyear, ] <- unlist(100 * (temp_sim_mean - true_density[iyear, ]) / 
#                                   true_density[iyear, ])
#   }
#   
#   return(list("mean_denisty" = round(mean_density, 2),
#               "cv" = round(cv, 4),
#               "rel_bias" = round(rel_bias, 2) ))
#   
# }

##################################
## Create function to calcualte a STRS and output mean, CV, and relative bias
##################################
do_STRS <- function(input){
  
  #Some constants
  n_cells <- length(input$domain_idx)
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
  
  strata_to_use <- survey_detail$Stratum#[survey_detail$nh > 0]
  input$solution <- input$solution[input$domain_idx]
  
  
  #Result objects
  mean_density <- cv <- rel_bias <- matrix(ncol = n_spp, 
                                           nrow = n_time)
  
  for (iyear in 1:n_time) {
    
    #Subset density df by year
    sub_df <- input$density[input$domain_idx, , iyear]
    
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
