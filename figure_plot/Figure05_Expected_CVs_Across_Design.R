###############################################################################
## Project:         Plot Expected CVs of different sampling designs
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:     Show only two-boat solutions
###############################################################################
rm(list = ls())

##################################################
####  Set up directories  
##################################################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[2]
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
library(readxl)

##################################################
####  Load Data
##################################################
load("data/processed/optimization_data.RData")
load("data/processed/grid_goa.RData")
load("data/processed/VAST_fit_D_gct.RData")

MS_objects <- load("results/MS_optimization_knitted_results.RData")

for (ivar in MS_objects) assign(x = paste0(ivar, "_MS"), value = get(ivar))

SS_objects <- load(paste0("results/full_domain/Single_Species_Optimization",
                          "/optimization_knitted_results.RData"))
for (ivar in SS_objects) assign(x = paste0(ivar, "_SS"), value = get(ivar))
rm(list = SS_objects)

goa_allocations <- 
  readxl::read_xlsx(path = 'data/GOA 2019 stations by stratum.xlsx')

##################################################
####  Current survey allocation
##################################################
## Rename Current Stratum labels
new_strata_labels = 1:length(unique(grid_goa$stratum))
names(new_strata_labels) <- sort(unique(grid_goa$stratum))

grid_goa$stratum_new_label <- 
  new_strata_labels[paste(grid_goa$stratum)]

allocations <- data.frame(stratum = sort(unique(goa_allocations$Stratum)),
                          nh = c(goa_allocations$`Number stations`))

allocations <- rbind(data.frame(stratum = 0, nh = 0),
                     allocations,
                     data.frame(stratum = c(510, 520, 530, 540, 550), 
                                nh = rep(0, 5)))

allocations$stratum <- 1:nrow(allocations)

##################################################
####   General dataframe used to calculate CVs
##################################################
frame <- cbind(data.frame(domainvalue = 1,
                          id = 1:n_cells,
                          WEIGHT = n_years),
               
               matrix(data = apply(X = D_gct,
                                   MARGIN = c(1, 2), 
                                   FUN = sum),
                      ncol = ns_all,
                      dimnames = list(NULL, paste0("Y", 1:ns_all))),
               
               matrix(data = apply(X = D_gct,
                                   MARGIN = c(1, 2), 
                                   FUN = function(x) sum(x^2)),
                      ncol = ns_all,
                      dimnames = list(NULL, paste0("Y", 1:ns_all, "_SQ_SUM")))
)

##################################################
####   Calculate Expected CV under Simple RS across ALL species
##################################################
srs_stats <- SamplingStrata::buildStrataDF( 
  dataset = cbind( frame, X1 = 1))

srs_n <- samples[2]
srs_var <- as.numeric((srs_stats[, paste0("S", 1:ns_all)])^2) * (1 - srs_n / n_cells) / srs_n
srs_mean <- as.numeric(srs_stats[, paste0("M", 1:ns_all)])

srs_cv <- as.numeric(sqrt(srs_var) / srs_mean)

##################################################
####  Calculate Expected CV of current design across ALL species
##################################################
frame$X1 <- grid_goa$stratum_new_label
strata_stats <- SamplingStrata::buildStrataDF(dataset = frame)
strata_stats <- strata_stats[order(as.numeric(strata_stats$STRATO)), ]
strata_stats$SOLUZ <-  allocations$nh
strata_stats <- subset(strata_stats, SOLUZ > 0)
temp_cv <- unlist(SamplingStrata::expected_CV(strata = strata_stats))

current_strs_cv <- as.numeric(temp_cv)


##################################################
####  Calculate expected CV for all species for the 
####  gulf-wide (full_domain) 15 strata solution (isol == 8)
##################################################
isol <- 2
frame$X1 <- res_df_MS[, isol]
strata_stats <- SamplingStrata::buildStrataDF(dataset = frame)
strata_stats <- strata_stats[order(as.numeric(strata_stats$STRATO)), ]
strata_stats$SOLUZ <- strata_stats_list_MS[[isol]][, "SOLUZ"] 
strs_pop_cv <- as.numeric(SamplingStrata::expected_CV(strata = strata_stats))

##################################################
####  Plot
##################################################
{
  png(filename = paste0(output_dir, "Figure05_Expected_CVs_Across_Design.png"),
      width = 180,
      height = 180,
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
       xlim = c(0, 0.8),
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
       at = ns_all:(ns_all-ns_opt+1),
       col = "black",
       las = 1)
  axis(side = 2, 
       labels = common_names_all[spp_idx_eval], 
       at = ns_eval:1,
       col.axis = "darkgrey",
       las = 1)
  
  ## Plot Expected CVs
  matpoints(y = ns_all:1,
            x = cbind(srs_cv[c(spp_idx_opt, spp_idx_eval)], 
                      strs_pop_cv[c(spp_idx_opt, spp_idx_eval)],
                      current_strs_cv[c(spp_idx_opt, spp_idx_eval)]),
            pch = 16,
            col = c("black", "orange", "green"),
            cex = 1)
  
  points(y = ns_all:1,
         x = subset(x = settings_SS,
                    subset = boat == 2)[c(spp_idx_opt, 
                                          spp_idx_eval), "cv"],
         pch = 16,
         col = "red",
         cex = 1)
  
  ## Legend
  legend("topright",
         pch = 16,
         col = c("red", "green", "orange", "black"),
         legend = c("SS STRS", "Current STRS", "Optimized STRS", "SRS"),
         cex = 1)
  
  dev.off()
}
