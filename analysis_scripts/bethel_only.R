rm(list = ls())

##################################################
####  Set up directories  
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[1]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")

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
threshold <- c(0.1, 0.15, 0.20)

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
 over_or_under <- ifelse(n < 548, "under", 
                         ifelse(n > 553, "over", 
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

par(mar = c(3, 12, 1, 1))
plot(1,
     type = "n", 
     xlim = c(0, 0.4),
     ylim = c(1, 15),
     axes = F,
     ann = F)
box()

axis(side = 1)
axis(side = 2, 
     labels = common_names_opt, 
     at = 1:ns_opt, 
     las = 1)

points(y = 1:ns_opt,
       x = settings[sol_idx, paste0("CV_", 1:ns_opt)],
       col = "grey", 
       pch = 16)
matpoints(y = 1:ns_opt,
       x = tradeoff_df,
       pch = 16)
abline(v = threshold,
       lty = "dashed",
       col = palette())

legend("bottomright", 
       legend = c(0, threshold), 
       title = "Lower Population\nCV Thresholds",
       pch = 16, 
       col = c("grey", palette()), 
       cex = 1.5, 
       bty = "n")

