###############################################################################
## Project:         Plot Expected CVs of different sampling designs
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     
###############################################################################
rm(list = ls())

##################################################
####  Set up directories  
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/")[which_machine], 
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
            "data/Extrapolation_depths.RData"))
load(paste0(github_dir, 
            "data/fit_density.RData"))
load(paste0(github_dir,
            "results/MS_optimization_knitted_results.RData"))
load(paste0(github_dir, 
            "results/full_domain/Single_Species_Optimization",
            "/optimization_knitted_results.RData"))
load(paste0(github_dir, 
            "results/full_domain/srs_pop_cv.RData"))

GOA_allocations <- readxl::read_xlsx(
  path = paste0(github_dir, 
                '/data/GOA 2019 stations by stratum.xlsx'))

GOA_allocations3 <- readxl::read_xlsx(
  path = paste0(github_dir, 
                '/data/GOA2019_ 3 boat_825_RNDM_stations.xlsx')) 

##################################################
####  Current survey allocation
##################################################
## Rename Current Stratum labels
new_strata_labels = 1:length(unique(Extrapolation_depths$stratum))
names(new_strata_labels) <- sort(unique(Extrapolation_depths$stratum))

Extrapolation_depths$stratum_new_label <- 
  new_strata_labels[paste(Extrapolation_depths$stratum)]

allocations <- data.frame(Stratum = sort(unique(GOA_allocations3$stratum)),
                          boat3 = aggregate(id ~ stratum, 
                                            data = GOA_allocations3, 
                                            FUN = length)$id,
                          boat2 = c(GOA_allocations$`Number stations`,
                                    rep(0, 5)))
allocations$boat1 <- ceiling(allocations$boat2 / 2)

allocations$boat1 <- ifelse(allocations$boat1 == 0, 0, 
                            ifelse(allocations$boat1 == 1, 2, 
                                   allocations$boat1))

allocations <- rbind(data.frame(Stratum = 0, boat3 = 0, boat2 = 0, boat1 = 0),
                     allocations)

allocations$Stratum <- 1:nrow(allocations)

##################################################
####   General dataframe used to calculate CVs
##################################################
frame <- cbind(data.frame(domainvalue = 1,
                          id = 1:n_cells,
                          WEIGHT = n_years),
               
               matrix(data = apply(X = D_gct[, , years_included],
                                   MARGIN = c(1, 2), 
                                   FUN = sum),
                      ncol = ns_all,
                      dimnames = list(NULL, paste0("Y", 1:ns_all))),
               
               matrix(data = apply(X = D_gct[, , years_included],
                                   MARGIN = c(1, 2), 
                                   FUN = function(x) sum(x^2)),
                      ncol = ns_all,
                      dimnames = list(NULL, paste0("Y", 1:ns_all, "_SQ_SUM")))
)

##################################################
####  Calculate Expected CV of current design across ALL species
##################################################
current_pop_cv <- matrix(nrow = ns_all, 
                         ncol = n_boats)

for (isol in 1:n_boats) {
  frame$X1 <- Extrapolation_depths$stratum_new_label
  strata_stats <- SamplingStrata::buildStrataDF(dataset = frame)
  strata_stats <- strata_stats[order(as.numeric(strata_stats$STRATO)), ]
  strata_stats$SOLUZ <-  allocations[, paste0("boat", isol)]
  strata_stats <- subset(strata_stats, SOLUZ > 0)
  temp_cv <- unlist(SamplingStrata::expected_CV(strata = strata_stats))
  
  current_pop_cv[, isol] <- temp_cv
}

##################################################
####  Calculate expected CV for all species for the 
####  gulf-wide (full_domain) 15 strata solution (isol == 11)
##################################################
isol <- 14
frame$X1 <- res_df[, isol]
strata_stats <- SamplingStrata::buildStrataDF(dataset = frame)
# strata_stats <- strata_stats[order(as.numeric(strata_stats$STRATO)), ]
strata_stats$SOLUZ <- strata_stats_list[[isol]][, "SOLUZ"] 
strs_pop_cv <- unlist(SamplingStrata::expected_CV(strata = strata_stats))

##################################################
####  Plot
##################################################
{
  png(filename = paste0(output_dir, "Expected_CVs_Across_Design.png"),
      width = 190,
      height = 200,
      units = "mm",
      res = 500)
  
  ## Layout Plot
  par(mar = c(4,11,1,1))
  
  ## Base plot
  plot(1,
       type = "n",
       pch = 16,
       cex = 1.5,
       ylim = c(1, ns_all),
       xlim = c(0, 0.70),
       axes = F,
       ann = F)
  box()
  abline(h = 1:ns_all, 
         v = seq(from = 0, to = 1, by = 0.1),
         col = "lightgrey", 
         lty = "dashed")
  mtext(side = 1, 
        text = "Expected CV",
        line = 2.5,
        cex = 1.25,
        font = 2)
  axis(side = 1, at = seq(from = 0, to = 1, by = 0.05))
  axis(side = 2, 
       labels = common_names_all[c(spp_idx_opt)], 
       at = 1:ns_opt,
       col = "black",
       las = 1)
  axis(side = 2, 
       labels = common_names_all[spp_idx_eval], 
       at = (ns_opt + 1):ns_all,
       col.axis = "darkgrey",
       las = 1)
  
  ## Plot Expected CVs
  matpoints(y = 1:ns_all,
            x = cbind(srs_pop_cv_full_domain[c(spp_idx_opt, spp_idx_eval), 1, 2], 
                      strs_pop_cv[c(spp_idx_opt, spp_idx_eval)],
                      current_pop_cv[c(spp_idx_opt, spp_idx_eval), 2] ),
            pch = 16,
            col = c("black", "blue", "green"),
            cex = 1)
  points(y = 1:ns_all,
         x = subset(x = settings_agg_full_domain,
                    subset = iboat == 2)[c(spp_idx_opt, 
                                           spp_idx_eval), "cv"],
         pch = 16,
         col = "red",
         cex = 1)
  
  ## Legend
  legend("bottomright",
         pch = 16,
         col = c("red", "green", "blue", "black"),
         legend = c("SS STRS", "Current STRS", "Optimized STRS", "SRS"),
         cex = 1)
  
  dev.off()
}
