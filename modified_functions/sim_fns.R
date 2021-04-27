###############################################################################
## Project:       Simulate Surveys Function
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Create function to calcualte a STRS and output mean, 
##                CV, and relative bias
###############################################################################

do_STRS <- function(input){
  
  #Some constants
  n_cells <- dim(input$density)[1]
  n_time <-  dim(input$density)[2]
  n_dom <- length(table(input$post_strata))
  
  survey_detail <- 
    data.frame("Stratum" = as.integer(names(table(input$solution))),
               "Nh" = as.integer(table(input$solution)),
               "nh" = input$allocation)
  
  #Take strata with 0 effor allocation out
  strata_to_use <- survey_detail$nh > 0
  survey_detail <- survey_detail[strata_to_use, ]
  
  #Assume stratum weights include untrawlabe areas
  survey_detail$Wh <- survey_detail$Nh / n_cells
  survey_detail$wh <- with(survey_detail, nh/Nh)
  
  #Strata Areas
  strata_areas <- aggregate(cell_areas ~ solution, 
                            FUN = sum,
                            data = with(input, data.frame(solution, 
                                                          cell_areas)))
  strata_areas <- subset(strata_areas, solution %in% survey_detail$Stratum)
  
  #Result objects
  mean_density <- cv <- index <- rel_bias <- rel_log_bias <- c()
  index_district <- array(dim = c(n_time, n_dom))
  
  for (iyear in 1:n_time) {
    
    #Subset density df by year
    sub_df <- input$density[, iyear]
    
    #Take a random sample based on the allocation and stratum
    sample_vec <- c()
    for(istrata in 1:nrow(survey_detail)) {
      sample_vec <- c(sample_vec,
                      sample(x = which(input$solution == survey_detail$Stratum[istrata]),
                             size = survey_detail$nh[istrata]) )
    }
    
    sampled_strata <- rep(x = survey_detail$Stratum, 
                          times = survey_detail$nh)
    
    #subset sub_df by which cells were chosen
    sample_df <- as.vector(sub_df[sample_vec])
    
    #Calculate STRS mean density
    strata_mean <- tapply(X = sample_df, 
                          INDEX = sampled_strata,
                          FUN = mean)
    STRS_mean <- sum(strata_mean * survey_detail$Wh)
    
    #Calculate STRS variance of mean density
    strata_var <- tapply(X = sample_df, 
                         INDEX = sampled_strata,
                         FUN = var)
    STRS_var <- sum(strata_var * with(survey_detail, Wh^2 * (1 - wh) / nh) )
    
    #Save mean and cv of estimates across species
    cv[iyear] <- sqrt(STRS_var) / STRS_mean 
    
    #Calculate index of abundance by district
    index_df <- data.frame(Area_km2 = input$cell_areas,
                           stratum = input$solution,
                           district = input$post_strata,
                           mean_dens = strata_mean[paste(input$solution)])
    
    index_district[iyear, ] <- 
      tapply(X = index_df$mean_dens * index_df$Area_km2 * 0.001,
             INDEX = index_df$district,
             FUN = sum,
             na.rm = TRUE)
    
    # Calculate total index
    index[iyear] <- sum(strata_areas * strata_mean) * 0.001
    
  }
  
  #Calculate Relative bias of index over the entire domain and by districts
  rel_bias <-   
    100 * (index - rowSums(input$true_index_district)) /
    rowSums(input$true_index_district)
  rel_log_bias <- 
    log10(rowSums(index_district) / rowSums(input$true_index_district))
  
  bias_index_district <- 100 * (index_district - input$true_index_district) /
    input$true_index_district
  log_bias_index_district <- log10(index_district / 
                                     input$true_index_district)
  
  return(list("strs_mean" = STRS_mean,
              "strs_index" = index,
              "cv" = round(cv, 4),
              "rel_bias" = round(rel_bias, 2),
              "rel_log_bias" = round(rel_log_bias, 3),
              "bias_index_district" = round(bias_index_district, 2),
              "log_bias_index_district" = round(log_bias_index_district, 3) )
  )
}