###############################################################################
## Project:         Expected CV constrints versus Simulated True CV
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     
###############################################################################
rm(list = ls())

##################################################
####  Install a forked version of the SamplingStrata Package from 
####  zoyafuso-NOAA's Github page
####
####  Import other required packages
##################################################
library(SamplingStrata)

##################################################
####   Set up directories based on whether the optimization is being conducted
####        on a multi-species or single-species level
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents", 
                       "C:/Users/Zack Oyafuso/Documents")[which_machine],
                     "/GitHub/Optimal_Allocation_GoA/")

output_dir <- paste0(c("/Users/zackoyafuso/Google Drive/",
                       "C:/Users/Zack Oyafuso/Google Drive/")[which_machine],
                     "MS_Optimizations/TechMemo/figures/")

##################################################
####   Import Data
##################################################
load(paste0(github_dir,  "data/fit_density.RData"))
load(paste0(github_dir,  "data/Extrapolation_depths.RData"))
load(paste0(github_dir,  "data/optimization_data.RData"))
load(paste0(github_dir, "results/MS_optimization_knitted_results.RData"))
load(paste0(github_dir, "results/SUR_opt_district_STR_5_simulation_results.RData"))

idx <- which(settings$domain == "district" & 
               settings$strata == 5 & 
               settings$boat == 2)

# expected_sample_cv <- matrix(ncol = ns_all, nrow = n_years)
# 
# for (iyear in 1:n_years) {
#   frame_district <- cbind(data.frame(
#     domainvalue = cut(x = Extrapolation_depths$Lon,
#                       breaks = c(-170, -159, -154, -147, -140, -132),
#                       labels = 1:5),
#     id = 1:n_cells,
#     X1 = res_df[, idx],
#     WEIGHT = 1),
# 
#     matrix(data = D_gct[, , years_included[iyear]],
#            ncol = ns_all,
#            dimnames = list(NULL, paste0("Y", 1:ns_all))),
# 
#     matrix(data = D_gct[, , years_included[iyear]]^2,
#            ncol = ns_all,
#            dimnames = list(NULL, paste0("Y", 1:ns_all, "_SQ_SUM")))
#   )
# 
#   temp_agg_strata <- buildStrataDF(dataset = frame_district)
#   temp_agg_strata$SOLUZ <- strata_stats_list[[idx]]$SOLUZ
#   temp_agg_strata$DOM1 <- 1
# 
#   expected_sample_cv[iyear, ] <-
#     as.numeric(expected_CV(strata= temp_agg_strata))
# 
# }

{
  
  png(filename = paste0(output_dir, "True_vs_Expected_CVs.png"), 
      width = 170, 
      height = 170, 
      units = "mm", 
      res = 500) 
  
  par(mar = c(3, 10, 3, 1), mfrow = c(1, 1))
  boxplot(SUR_opt_district_STR_5_true_cv["obsCV=0", 
                                         paste0("year_", 1:n_years),
                                         spp_idx_opt,
                                         "boat_2"],
          horizontal = T,
          las = 1, 
          names = NA,
          ylim = c(0, 0.5),
          col = "white",
          main = "Distribution of true CV vs expected CV (red points)")
  
  
  axis(side = 2, 
       labels = common_names_opt,
       at = 1:ns_opt, 
       las = 1)
  points(x = unlist(settings[idx, paste0("CV_", 1:ns_opt)]),
         y = 1:ns_opt,
         pch = 16,
         col = "red",
         cex = 1.5)
  
  dev.off()
}

# boxplot(100*(SUR_opt_district_STR_5_true_cv["obsCV=0",
#                                             paste0("year_", 1:n_years),
#                                             spp_idx_opt,
#                                             "boat_2"] -
#                expected_sample_cv[, spp_idx_opt]) /
#           expected_sample_cv[, spp_idx_opt],
#         horizontal = T,
#         las = 1,
#         pch = 16,
#         col = "white",
#         names = NA,
#         ylim = c(-12, 12),
#         main = "Year-specific Expected CV versus True CV")
# 
# axis(side = 2, labels = common_names_opt, at = 1:ns_opt, las = 1)
# abline(v = 0, lty = "dashed", col = "grey" )
