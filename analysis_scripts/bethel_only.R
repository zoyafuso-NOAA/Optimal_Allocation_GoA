rm(list = ls())

##################################################
####  Set up directories  
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[3]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")
output_dir <- paste0(c("/Users/zackoyafuso/", 
                       "C:/Users/Zack Oyafuso/")[which_machine], 
                     "Google Drive/MS_Optimizations/TechMemo/figures/")

##################################################
####  Install a forked version of the SamplingStrata Package from 
####  zoyafuso-NOAA's Github page
####
####  Import other required packages
##################################################
library(devtools)
devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata")
library(SamplingStrata)

##################################################
####  Load Data
##################################################
load(paste0(github_dir, 
            "data/optimization_data.RData"))
load(paste0(github_dir,
            "results/MS_optimization_knitted_results.RData"))
load(paste0(github_dir, 
            "results/full_domain/srs_pop_cv.RData"))

##################################################
####  Isolate District-Level optimization, 3 strata per district
####  Set lower threshold values
##################################################
sol_idx <- 2
threshold <- c(0.1, 0.2, 0.25)

strata_df <- strata_stats_list[[sol_idx]]
strata_df$stratum = 1:nrow(strata_df)
strata_df$DOM1 = 1

##################################################
####  Result Object
##################################################
tradeoff_df <- matrix(nrow = ns_opt, ncol = length(threshold))

for (ithreshold in 1:length(threshold)) {
  ##################################################
  ####  Calculate allocation at Simple Random Sampling Population CV
  ####  and 2 boats, set lower CV threshold
  ##################################################
  error_df <- cbind(data.frame(DOM = "DOM1"),
                    t(srs_pop_cv_full_domain[spp_idx_opt, , 2]),
                    data.frame(domainvalue = 1))
  names(error_df)[names(error_df) %in% sci_names_opt] <- 
    paste0("CV", 1:ns_opt)
  
  error_df[, paste0("CV", 1:ns_opt)] <- 
    ifelse(test = error_df[, paste0("CV", 1:ns_opt)] <= threshold[ithreshold], 
           yes = threshold[ithreshold], 
           no = error_df[, paste0("CV", 1:ns_opt)])
  
  ##################################################
  ####  Run Bethel algorithm and calculate total sample size
  ##################################################
  bethel_allocation  <- bethel(errors = error_df, 
                               stratif = strata_df, 
                               printa=TRUE)
  pop_cvs <- as.numeric(attributes(bethel_allocation)$outcv[, "ACTUAL CV"])
  (n <- sum(bethel_allocation))
  
  ##################################################
  ####  Adjust CVs to be within 550 samples (desired under two boats)
  ##################################################
  over_or_under <- ifelse(n < 549, "under", 
                          ifelse(n > 551, "over", 
                                 "exact"))
  
  while(over_or_under %in% c("under", "over")) {
    change_rate <- switch(over_or_under, 
                          "under" = 0.985,
                          "over" = 1.015,
                          "exact" = 1)
    error_df[, paste0("CV", 1:ns_opt)] <- 
      pop_cvs * change_rate
    
    error_df[, paste0("CV", 1:ns_opt)] <- 
      ifelse(test = error_df[, paste0("CV", 1:ns_opt)] <= threshold[ithreshold], 
             yes = threshold[ithreshold], 
             no = error_df[, paste0("CV", 1:ns_opt)])
    
    bethel_allocation  <- bethel(errors = error_df, 
                                 stratif = strata_df, 
                                 printa=TRUE)
    
    pop_cvs <- as.numeric(attributes(bethel_allocation)$outcv[, "ACTUAL CV"])
    
    print(n <- sum(bethel_allocation))
    
    over_or_under <- ifelse(n < 550, "under", 
                            ifelse(n > 550, "over", 
                                   "exact"))
  }
  
  tradeoff_df[, ithreshold] <- pop_cvs
  
}

{
  png(filename = paste0(output_dir, "bethel_only_comparisons.png"), 
      width = 170, 
      height = 170, 
      units = "mm", 
      res = 500)
  
  par(mar = c(2, 1, 2, 1), oma = c(2, 10, 0, 0), mfrow = c(2, 2))
  for (i in 1:4) {
    plot(1,
         type = "n", 
         xlim = c(0, 0.4),
         ylim = c(1, 15),
         axes = F,
         ann = F)
    box()
    
    spp_order <- c(1:9,15,11, 10, 12:14)
    
    axis(side = 1)
    
    mtext(side = 3, 
          text = paste(LETTERS[i], ") ", 
                       c(0, threshold)[i],
                       " as Lower CV Threshold"),
          line = 0.25)
    
    if (i %in% c(1, 3)) axis(side = 2, 
                             labels = common_names_opt[spp_order], 
                             at = 1:ns_opt, 
                             las = 1)
    
    points(y = 1:ns_opt,
           x = settings[sol_idx, paste0("CV_", spp_order)],
           col = "black", 
           pch = 16)
    
    if (i > 1)   {
      matpoints(y = 1:ns_opt,
                x = tradeoff_df[spp_order, 1:(i - 1)],
                pch = 16,
                col = c("blue", "red", "green")[1:(i - 1)])
      
      abline(v = threshold[i - 1],
             lty = "dashed",
             col = c("blue", "red", "green")[i - 1])
    }
  }
  
  mtext(side = 1, 
        text = "Expected CV",
        outer = TRUE,
        line = 0.5)
  
  dev.off()
}
