###############################################################################
## Project:       Synthesize Optimization Results
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Synthesize all optimization results
##                sample sizes, expected CVs, solutions, allocations, etc
###############################################################################
rm(list = ls())

##################################################
####  Set up directories
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2, "Zack_GI_PC" = 3)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")

##################################################
####   Import Packages
##################################################
library(SamplingStrata)

##################################################
####   Load Data
##################################################
load(paste0(github_dir, "data/optimization_data.RData"))

scen <- data.frame(domain = rep(c("district", "full_domain"), each = 2),
                   strata = c(3, 5, 10, 15))

###########################
## Empty Result objects
###########################
master_res_df <- data.frame(id = 1:n_cells)
master_settings <- master_settings_district <- data.frame()
master_strata_list <- master_strata_stats_list <- list()

idx <- 0
for(irow in 1:nrow(scen)) {
   
   ###########################
   ## Indices
   ##########################
   idom <- scen$domain[irow]
   istrata <- scen$strata[irow]
   
   frame <- switch( idom,
                    "full_domain" = frame_all,
                    "district" = frame_district)
   
   n_dom <- length(unique(frame$domainvalue))
   
   for (iboat in 1:n_boats) {
      result_dir <- paste0(github_dir, "results/", idom, 
                           "/Multi_Species_Optimization/boat", 
                           iboat, "/Str", istrata, "/")
      
      runs <- grep(x = dir(result_dir), 
                   pattern = "Run", 
                   value = TRUE)
      
      for (irun in runs) {
         temp_dir <- paste0(result_dir,  irun, "/result_list.RData")
         
         if (file.exists(temp_dir)) {
            idx <- idx + 1
            load(temp_dir)
            
            ## Solution: which strata is assigned to each extrapolation cell
            solution <- 
               switch(idom,
                      "full_domain" = result_list$solution$indices$X1,
                      "district" = as.factor(paste0(
                         "DOM", result_list$solution$framenew$DOMAINVALUE,
                         " STR", result_list$solution$framenew$STRATO))
               )
            
            solution <- as.integer(solution)
            master_res_df <- cbind(master_res_df, 
                                   solution )
            
            ## Strata characteristics: sample size, population, sampling rate, 
            ## strata variable cuts
            master_strata_list <- c(master_strata_list, 
                                    list(result_list[[2]]))
            
            ## Strata statistics (mean and variance)
            temp_strata_stats_list <- result_list[[1]]$aggr_strata
            strata_order <- order(as.integer(temp_strata_stats_list$DOM1),
                                  as.integer(temp_strata_stats_list$STRATO))
            temp_strata_stats_list <- temp_strata_stats_list[strata_order, ]
            master_strata_stats_list <- c(master_strata_stats_list, 
                                          list(temp_strata_stats_list))
            
            ## High-level settings: total sample size and expected CV across
            ## species
            species_cv <- switch(
               idom,
               "full_domain" = {
                  cv_agg <- result_list[[3]]
                  attributes(cv_agg)$dimnames[[1]] <- ""
                  attributes(cv_agg)$dimnames[[2]] <- paste0("CV", 1:ns_opt)
                  cv_agg
               },
               "district" = {
                  agg_strata <- result_list$solution$aggr_strata
                  agg_strata$STRATO <- 1:nrow(agg_strata)
                  agg_strata$DOM1 <- 1
                  
                  cv_agg <- matrix(as.numeric(
                     SamplingStrata::expected_CV(strata = agg_strata)),
                     nrow = 1)
                  colnames(cv_agg) <- paste0("CV", 1:ns_opt)
                  cv_agg
               })
            
            n <- result_list$n
            
            master_settings <- rbind(master_settings, 
                                     cbind(data.frame(domain = idom,
                                                      id = idx,
                                                      strata = istrata, 
                                                      boat = iboat,
                                                      n), species_cv))
            
            
            ## High-level settings by district: total sample size and 
            ## expected CV across species 
            species_cv <- result_list[[3]]
            attributes(species_cv)$dimnames[[1]] <- rep("", n_dom)
            attributes(species_cv)$dimnames[[2]] <- paste0("CV", 1:ns_opt)
            
            master_settings_district <- rbind(
               master_settings_district,
               cbind(data.frame(id = idx,
                                boat = iboat,
                                strata = istrata,
                                domain = idom,
                                domain_no = 1:n_dom,
                                n = tapply(X = result_list$sum_stats$Allocation,
                                           INDEX = result_list$sum_stats$Domain,
                                           FUN = sum), species_cv))
            )
         }
      }
   }
}

##################################################
####   Subset to only the solutions closest to 280, 550, and 820
##################################################
idx <-  as.vector(
   #First split master_settings by domain and strata level
   sapply(X = split.data.frame(x = master_settings, 
                               f = list(master_settings$domain,
                                        master_settings$strata),
                               drop = TRUE,
                               sep = "_"),
          
          ## on each split, calculate which solution id is closest to the 
          ## 1, 2, and 3, boat (280, 550, 820 stations) scenarios
          FUN = function(x) {
             sapply(X = samples, 
                    FUN = function(y) 
                       x$id[which.min(abs(x$n - y))])
          }
   )
)

settings <- subset(x = master_settings, 
                   subset = id %in% idx,
                   select = -id)
settings_district <- subset(x = master_settings_district,
                            subset = id %in% idx,
                            select = -id)
res_df <- master_res_df[, 1 + idx]
strata_list <- master_strata_list[idx]
strata_stats_list <- master_strata_stats_list[idx]

names(res_df) <- names(strata_list) <- names(strata_stats_list) <- 
   paste0("sol_", 1:nrow(settings))

##################################################
####   Adjust sample size 
##################################################   
for(irow in 1:nrow(settings)) {
   
   temp_boat <- settings$boat[irow]
   temp_n <-  settings$n[irow]
   if((temp_n >= (samples[temp_boat] - 1) & 
       temp_n <= (samples[temp_boat] + 1)) ) next
   
   
   temp_domain <- settings$domain[irow]
   temp_strata <- settings$strata[irow]
   temp_cv <- subset(x = settings_district,
                     subset = boat == temp_boat &
                        strata == temp_strata &
                        domain == temp_domain,
                     select = paste0("CV", 1:ns_opt))

   
   temp_strata_stats_list <- subset(x = strata_stats_list[[irow]],
                                    select = -SOLUZ)
   n_dom <- ifelse(temp_domain == "full_domain", 1, 5)
   temp_n_by_strata <- actual_n_by_strata <- 
      vector(length = nrow(temp_strata_stats_list))
   names(temp_n_by_strata) <- names(actual_n_by_strata) <- 
      temp_strata_stats_list$DOM1
   
   while ( !(temp_n >= (samples[temp_boat] - 1) & 
             temp_n <= (samples[temp_boat] + 1)) ) {
      over_or_under <- temp_n > samples[temp_boat]
      
      temp_cv <- temp_cv * (1 + ifelse(test = over_or_under == T, 
                                       yes = 0.0005,
                                       no = -0.0005))
      temp_n_by_reg <- vector(length = n_dom)
      
      for (i in 1:n_dom) {
         error_df <- cbind(data.frame("DOM" = "DOM1",
                                      temp_cv[i, ],
                                      "domainvalue"  = 1))
         temp_bethel <- SamplingStrata::bethel(
            errors = error_df,
            stratif = subset(x = temp_strata_stats_list, subset = DOM1 == i), 
            realAllocation = T, 
            printa = T)
         
         temp_n_by_reg[i] <- sum(as.numeric(ceiling(temp_bethel)))
         
         temp_n_by_strata[names(temp_n_by_strata) == paste(i)] <- ceiling( temp_bethel)
         actual_n_by_strata[names(actual_n_by_strata) == paste(i)] <- temp_bethel
      }
      
      temp_n <- sum(temp_n_by_reg)
      print(temp_n)
   }
   
   temp_cv_agg <- SamplingStrata::expected_CV(
      strata = cbind(subset(temp_strata_stats_list,
                            select = -DOM1),
                     DOM1 = 1,
                     SOLUZ = actual_n_by_strata))
   
   ####################################
   ## Record new results
   ####################################
   settings[irow, paste0("CV", 1:ns_opt)] <- as.numeric(temp_cv_agg)
   settings[irow, "n"] <- temp_n
   strata_stats_list[[irow]]$SOLUZ <- actual_n_by_strata 
   strata_list[[irow]]$Allocation <- temp_n_by_strata
   
   settings_district[settings_district$boat == temp_boat &
                        settings_district$strata == temp_strata &
                        settings_district$domain == temp_domain, 
                     paste0("CV", 1:ns_opt)] <- temp_cv
}

##################################################
####   Save Objects
##################################################
save(list = c("res_df", "settings", "settings_district",
              "strata_list", "strata_stats_list"),
     file = paste0(github_dir,
                   "results/MS_optimization_knitted_results.RData"))
