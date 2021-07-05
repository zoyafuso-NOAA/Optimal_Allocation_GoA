###############################################################################
## Project:       True and RRMSE of CV, Survey Comparison
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Compare True CV, and RRMSE of CV across species
##                for current and optimized STRS surveys
###############################################################################
rm(list = ls())

##################################################
####  Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2)[2]
output_dir <- paste0(c("/Users/zackoyafuso/",
                       "C:/Users/Zack Oyafuso/")[which_machine],
                     "Google Drive/MS_Optimizations/TechMemo/")

##################################################
####   Set up appendix output directories
##################################################
if(!dir.exists(paste0(output_dir, "/appendix/Appendix E plots/" )))
  dir.create(paste0(output_dir, "/appendix/Appendix E plots/" ))

##################################
## Import plotting percentile time series function
##################################
source("modified_functions/plot_percentiles.R")
library(SamplingStrata)
library(readxl)

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load("data/processed/optimization_data.RData")
load("data/processed/grid_goa.RData")
load("results/MS_optimization_knitted_results.RData")

GOA_allocations <- 
  readxl::read_xlsx(path = 'data/GOA 2019 stations by stratum.xlsx')

GOA_allocations3 <- 
  readxl::read_xlsx(path = 'data/GOA2019_ 3 boat_825_RNDM_stations.xlsx')

##################################
## Specify Management Districts
##################################
new_strata_labels = 1:length(unique(grid_goa$stratum))
names(new_strata_labels) <- sort(unique(grid_goa$stratum))

##################################
## Rename Current Stratum labels
##################################
grid_goa$stratum_new_label <- 
  new_strata_labels[paste(grid_goa$stratum)]

##################################################
####   Create dataframe of effort allocations across boats
##################################################
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
####   Import simulation results
##################################################
scen <- data.frame(survey_type = c("cur", rep("opt", 4) ),
                   strata = c("cur", 3, 5, 10, 15),
                   domain = c("full_domain",
                              rep(c("district", "full_domain"), each = 2)))

for (irow in 1:nrow(scen)) {
  scen_name <- paste0(scen$domain[irow], "_STR_", scen$strata[irow], "_")
  file_name <- paste0("results/survey_simulations/", 
                      scen_name, "simulation_results.RData")
  
  load(file_name)
  
}


##################################################
####   Calculate Expected CV of each species for each survey design
##################################################
expected_sample_cv <- matrix(ncol = ns_all, 
                             nrow = nrow(settings) + 3,
                             dimnames = list(NULL, common_names_all),
                             #add three more row for the current design
)

for (iscen in 1:(nrow(settings) + 3)) {
  
  #Temporary operating model data frame
  temp_frame <- subset(frame_all, select = -X2)
  
  if( !(iscen %in% (nrow(settings) + (1:3))) ){ #if looping through the 
    #optimized solutions
    
    #Assign to X1 the solution (stratum label)
    temp_frame$X1 <- res_df[, iscen]
    
    #Calculate strata means and variances
    temp_agg_strata <- buildStrataDF(dataset = temp_frame)
    # temp_agg_strata <- temp_agg_strata[order(as.integer(temp_agg_strata$STRATO)), ]
    temp_agg_strata$DOM1 <- 1
    
    #Clean up order of strata for the district-level solutions
    if(settings$domain[iscen] == "district") {
      idx <- order(as.numeric(temp_agg_strata$STRATO))
      temp_agg_strata <- temp_agg_strata[idx, ]
      temp_agg_strata$SOLUZ <- strata_stats_list[[iscen]]$SOLUZ
    }
    
    #Clean up order of strata for the district-level solutions
    if(settings$domain[iscen] == "full_domain") {
      idx <- order(as.numeric(temp_agg_strata$STRATO))
      temp_agg_strata <- temp_agg_strata[idx, ]
      temp_agg_strata$SOLUZ <- strata_stats_list[[iscen]]$SOLUZ[idx]
    }
    
  }
  
  #Clean up order of strata for the current survey design
  if(iscen %in% (nrow(settings) + (1:3)) ) {
    temp_frame$X1 <- grid_goa$stratum_new_label
    temp_agg_strata <- buildStrataDF(dataset = temp_frame)
    
    idx <- order(as.numeric(temp_agg_strata$STRATO))
    temp_agg_strata <- temp_agg_strata[idx, ]
    
    nh <-  allocations[, paste0("boat", iscen - nrow(settings))]
    temp_agg_strata <- temp_agg_strata[nh > 0, ]
    temp_agg_strata$SOLUZ <- nh[nh > 0]
  }
  
  #Calculate expected STRS CV
  expected_sample_cv[iscen, ] <- 
    as.numeric(expected_CV(strata = temp_agg_strata))
  
}

##################################################
####   Plot True CV and RRMSE
##################################################

plot_scen <- expand.grid(metric = c("true_cv", "rrmse_cv"),
                         boat = c(2, 1, 3),
                         stringsAsFactors = FALSE)
plot_scen$fig_name <- c("figures/Figure08_TrueCV.png", 
                        "figures/Figure09_RRMSECV.png",
                        paste0("appendix/Appendix E plots/Appendix E", 
                               1:4, '.png') )

for(irow in 1:nrow(plot_scen)){
  plot_filename <- paste0(output_dir, plot_scen$fig_name[irow])
  imetric <- plot_scen$metric[irow]
  iboat <- plot_scen$boat[irow]
  
  ## Open Device
  png(filename = plot_filename,
      width = 170, height = 200, units = "mm", res = 500)
  
  
  ## Plot layout
  layout(mat = matrix(c(1:ns_all, ns_all + 1, ns_all + 1), 
                      ncol = 4, 
                      byrow = TRUE))
  par(mar = c(0.5, 3, 1.5, 1),
      oma = c(0, 2.5, 0, 0))
  
  for (ispp in c(spp_idx_opt, spp_idx_eval)) { ## Loop over species -- start
    
    ## Calculate max y for plotting
    merged_metric <- 
      lapply(X = with(scen, paste0(domain, "_STR_", strata, "_", imetric )), 
             FUN = function(x) get(x)[1:n_years,
                                      common_names_all[ispp],
                                      paste0("boat_", iboat) ])
    ymax_ <- max(unlist(merged_metric)) * 1.1
    
    ## Plot metric
    boxplot(merged_metric,
            ylim = c(0, max(0.10, ymax_) ),
            las = 1,
            axes = F,
            pch = 16,
            col  = c("white", "cyan", "blue4", "chartreuse1", "darkgreen"),
            cex = 0.5,
            at = c(1, 3:4, 6:7))
    box()
    axis(side = 2, las = 1)
    
    ## Separate scenarios
    abline(v = c(2, 5),
           lty = "dotted", 
           col = "black")
    
    ## Species label
    mtext(side = 3, text = common_names_all[ispp], 
          col = ifelse(ispp %in% spp_idx_opt, "black", "darkgrey"),
          cex = 0.8)
    
  } ## Loop over species -- start
  
  ## Plot legend
  par(mar = c(0, 0, 0, 0))
  plot(x = 1, y = 1,
       xlim = c(0,1), ylim = c(0,1),
       type = "n", axes = F, ann = F)
  
  legend(x = 0,
         y = 0.9,
         legend = c("Current Survey",
                    "District-Level Optimized Survey, 3 Strata per District",
                    "District-Level Optimized Survey, 5 Strata per District",
                    "Gulf-Wide Optimized Survey, 10 Total Strata",
                    "Gulf-Wide Optimized Survey, 15 Total Strata"),
         fill = c("white", 
                  "cyan",  "blue4",
                  "chartreuse1", "darkgreen"),
         xpd = NA,
         cex = 1,
         bty = "n")
  
  ## y-axis label
  mtext(side = 2, 
        outer = T, 
        line = 0.5,
        font = 2,
        text = switch(imetric, 
                      "rrmse_cv" = "RRMSE of CV",
                      "true_cv" = "True CV") )
  
  ## Close Device
  dev.off()
}

