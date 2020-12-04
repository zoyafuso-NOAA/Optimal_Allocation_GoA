##################################
## Create function to calcualte a STRS and output mean, CV, and relative bias
##################################
# input = list(
#   "density" = D_gct[, , Years2Include],
#   
#   "cell_areas" = Extrapolation_depths$Area_km2,
#   
#   "obs_CV" = obs_CV[ierror],
#   
#   "solution" = switch(
#     isurvey,
#     "Current" = Extrapolation_depths$stratum_new_label,
#     "STRS" = res_df[, iboat]),
#   
#   "allocation" = switch(
#     isurvey,
#     "Current" = allocations[, paste0("boat", iboat)],
#     "STRS" = strata_list[[iboat]]$Allocation),
#   
#   "true_density" = true_mean,
#   
#   "true_index_district" = true_index_district,
#   
#   "post_strata" = district_vals
# )

do_STRS <- function(input){
  
  #Some constants
  n_cells <- dim(input$density)[1]
  n_spp <-  dim(input$density)[2]
  n_time <- dim(input$density)[3]
  n_dom <- length(table(input$post_strata))
  
  survey_detail <- 
    data.frame("Stratum" = as.integer(names(table(input$solution))),
               "Nh" = as.integer(table(input$solution)),
               "nh" = input$allocation)
  
  #Assume stratum weights include untrawlabe areas
  survey_detail$Wh <- survey_detail$Nh / n_cells
  survey_detail$wh <- with(survey_detail, nh/Nh)
  
  #Take strata with 0 effor allocation out
  strata_to_use <- survey_detail$nh > 0
  survey_detail <- survey_detail[strata_to_use, ]
  
  #Strata Areas
  strata_areas <- aggregate(cell_areas ~ solution, 
                            FUN = sum,
                            data = with(input, data.frame(solution, cell_areas)))
  strata_areas <- subset(strata_areas, solution %in% survey_detail$Stratum)
  
  #Result objects
  mean_density <- cv <- rel_bias <- matrix(ncol = n_spp, 
                                           nrow = n_time)
  index_district <- array(dim = c(n_spp, n_time, n_dom))
  
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
    colnames(sample_df) <- paste0("Y", 1:n_spp)
    
    #add observation error proportional to predicted density
    sample_df <- apply(sample_df, 
                       MARGIN = 1:2,
                       FUN = function(x) {
                         m <- x
                         s <- x * input$obs_CV
                         location <- log(m^2 / sqrt(s^2 + m^2))
                         shape <- sqrt(log(1 + (s^2 / m^2)))
                         return(rlnorm(n = 1, location, shape))
                       })
    
    #Calculate STRS mean density
    stmt <- paste0('aggregate(cbind(',
                   paste0('Y', 1:(n_spp-1), sep = ',', collapse = ''), 
                   'Y', n_spp, 
                   ") ~ sampled_strata, data = sample_df, FUN = mean)")
    sample_mean <- eval(parse(text = stmt))[, -1]
    strata_mean <- sweep(x = sample_mean, 
                         MARGIN = 1, 
                         STATS = with(survey_detail, Wh),
                         FUN = '*')
    STRS_mean <- colSums(strata_mean)
    
    #Calculate STRS variance of mean density
    stmt <- paste0('aggregate(cbind(',
                   paste0('Y', 1:(n_spp-1), sep = ',', collapse = ''),
                   'Y',n_spp, 
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
      100 * (STRS_mean - input$true_density[, iyear]) / 
        input$true_density[, iyear])
    
    #Post-stratified indices of abundance
    index_df <- data.frame(Area_km2 = input$cell_areas,
                           stratum = input$solution,
                           district = input$post_strata)
    
    index_district[, iyear, ] <-  t(
      sapply(X = 1:n_spp, 
             FUN = function(x) {
               index_cell <- index_df$Area_km2 * 
                 sample_mean[index_df$stratum, paste0("Y", x)]
               tapply(X = index_cell,
                      INDEX = index_df$district,
                      FUN = sum,
                      na.rm = T)
             }) * 0.001 #Metric tonnes 
    )
  }
  
  return(list("mean_denisty" = round(mean_density, 2),
              "cv" = round(cv, 4),
              "rel_bias" = round(rel_bias, 2),
              "bias_index_district" = 100 * 
                (index_district - input$true_index_district) /
                input$true_index_district ))
}
