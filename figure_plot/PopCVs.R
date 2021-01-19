###############################################################################
## Project:         Plot Population CVs of different sampling designs
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     
###############################################################################
rm(list = ls())

##################################################
####  Set up directories  
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[1]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")
github_dir2 <- paste0(c("/Users/zackoyafuso/Documents/", 
                        "C:/Users/Zack Oyafuso/Documents/")[which_machine], 
                      "GitHub/MS_OM_GoA/data/")

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

dbe <- readRDS(paste0(dirname(github_dir), "/MS_OM_GoA/data/", 
                      "GOA_biomass_indices_wnames.rds"))
dbe <- subset(dbe, 
              SPECIES_NAME %in% sci_names_all &
                YEAR %in% (1996:2019)[years_included])

dbe$BIOMASS_CV <- sqrt(dbe$VAR_WGT_CPUE) / dbe$MEAN_WGT_CPUE
dbe <- dbe[order(dbe$YEAR),]

load(paste0(github_dir, 
            "data/vast_index.RData"))
vast_index <- subset(vast_index, 
                     spp %in% sci_names_opt)

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
####  Calculate Pop CV of optimized solutions across ALL species
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

current_pop_cv <- matrix(nrow = ns_all, 
                         ncol = n_boats, 
                         dimnames = list())

for (isol in 1:n_boats) {
  frame$X1 <- Extrapolation_depths$stratum_new_label
  
  strata_stats <- SamplingStrata::buildStrataDF(dataset = frame)
  nh <- allocations[, paste0("boat", isol)]
  Nh <- table(frame$X1)
  wh <- nh / Nh
  Wh <- Nh / sum(Nh)
  
  STRS_mean <- colSums(sweep(x = strata_stats[, paste0("M", 1:ns_all)], 
                             MARGIN = 1, 
                             STATS = Nh,
                             FUN = '*'))
  
  STRS_var <- as.matrix(sweep(x = strata_stats[, paste0("S", 1:ns_all)]^2, 
                              MARGIN = 1, 
                              STATS =  Nh^2 * (1 - wh) / nh,
                              FUN = '*'))
  
  STRS_var <- STRS_var[is.finite(rowSums(STRS_var)), ]
  
  STRS_var <- colSums(STRS_var)
  
  current_pop_cv[, isol] <- sqrt(STRS_var) / STRS_mean
}

rm(isol, strata_stats, nh, Nh, wh, Wh, STRS_mean, STRS_var)

##################################################
####  Calculate Pop CV of optimized solutions across ALL species
##################################################

strs_pop_cv <-  matrix(nrow = ns_all, ncol = nrow(settings), dimnames = list())

for (isol in 1:nrow(settings)) {
  frame$X1 <- res_df[, isol]
  
  strata_stats <- SamplingStrata::buildStrataDF(dataset = frame)
  strata_stats <- strata_stats[order(strata_stats$X1), ]
  
  nh <- strata_stats_list[[isol]]$SOLUZ
  Nh <- strata_stats_list[[isol]]$N
  wh <- nh / Nh
  # wh <- nh / table(frame$X1)
  Wh <- Nh / sum(Nh)
  
  STRS_mean <- colSums(sweep(x = strata_stats[, paste0("M", 1:ns_all)],
                             MARGIN = 1,
                             STATS = Nh,
                             FUN = '*'))
  
  STRS_var <- colSums(sweep(x = strata_stats[, paste0("S", 1:ns_all)]^2,
                            MARGIN = 1,
                            STATS =  Nh^2 * (1 - wh) / nh,
                            FUN = '*'))
  
  strs_pop_cv[, isol] <- sqrt(STRS_var) / STRS_mean
}

rm(isol, strata_stats, nh, Nh, wh, Wh, STRS_mean, STRS_var)

##################################################
####  Plot
##################################################
{
  png(filename = paste0(output_dir, "PopCVs.png"),
      width = 190,
      height = 200,
      units = "mm",
      res = 500)
  
  par(mar = c(4,11,1,1))
  plot(1,
       type = "n",
       pch = 16,
       cex = 1.5,
       ylim = c(1, ns_all),
       xlim = c(0, 0.90),
       axes = F,
       ann = F)
  
  box()
  abline(h = 1:ns_all, 
         v = seq(from = 0, to = 1, by = 0.1),
         col = "lightgrey", 
         lty = "dashed")
  mtext(side = 1, 
        text = "Population CV",
        line = 2.5,
        cex = 1.25,
        font = 2)
  
  matpoints(y = 1:ns_all,
            x = cbind(srs_pop_cv_full_domain[c(spp_idx_opt, spp_idx_eval), 1, 2], 
                      strs_pop_cv[c(spp_idx_opt, spp_idx_eval), 5],
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
  
  legend("bottomright",
         pch = 16,
         lwd = 2,
         col = c("red", 
                 "green", 
                 "blue", 
                 "black"),
         legend = c("SS STRS", 
                    "Current STRS", 
                    "Optimized STRS", 
                    "SRS"),
         cex = 1)
  
  dev.off()
}

###########################
## Distribution of observed DBE with the population CV 
###########################
# par(mfrow = c(5, 3), mar = c(2,3,2,1))
# for (ispp in c(1:10, 12:15)) {
#   temp_dbe <- subset( x = dbe,
#                       subset = SPECIES_NAME %in% sci_names_opt[ispp])
#   temp_vast <- subset( x = vast_index,
#                        subset = spp %in% sci_names_opt[ispp])
#   
#   boxplot(cbind(temp_dbe$BIOMASS_CV, temp_vast$cv),
#           ylim = c(0, max(cbind(temp_dbe$BIOMASS_CV, temp_vast$cv))),
#           at = 1:2,
#           xlim = c(0, 4),
#           col = "white",
#           las = 1,
#           names = NA)
#   
#   points(x = rep(3, 2), 
#          y = current_pop_cv[spp_idx_opt[ispp], c(2, 3)],
#          pch = 16,
#          col = c("black", "blue"),
#          cex = 2)
#   
#   mtext(side = 3, text = common_names_opt[ispp] )
# }
