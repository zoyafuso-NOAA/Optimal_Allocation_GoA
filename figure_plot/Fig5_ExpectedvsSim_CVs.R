###############################################################################
## Project:       Expected versus Simulated CVs
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Figure 5: Comparison of the relative difference between 
##                expected and realized coefficient of variation (CV) of 
##                abundance. Specifically, this shows the percent difference of 
##                the true CV relative to the upper CV constraint associated 
##                with a ten-strata, two boat of survey effort scenario 
##                (n = 550) across all included species. The left and center 
##                plots show optimizations using the one-CV constraint approach.
##                The right plot shows an optimization using the species-
##                specific CV constraint approach (refer to the main text for 
##                how CVs were specified across species). For the species-
##                specific CV constraint approach, a value of 0.10 was chosen 
##                as the lowest a CV constraint could be specified (indicated 
##                by the blue borders). Positive values indicate that the 
##                observed True CV is greater than what is expected from the 
##                optimization. Negative or near-zero values indicate that the 
##                observed true CV meets the CV expectations of the 
##                optimization. Results were qualitatively consistent with 
##                other scenarios.
###############################################################################
rm(list = ls())

############################
## Set up directories
#############################
which_machine <- c("Zack_MAC" = 1, "Zack_PC" = 2)[2]
VAST_model <- "6g"
github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/")

figure_dir <- paste0(c('/Users/zackoyafuso/', 
                       'C:/Users/Zack Oyafuso/')[which_machine],
                     'Google Drive/MS_Optimizations/figure_plot/')

load(paste0(github_dir, "optimization_data.RData"))

################
## Plot Settings
################
plot_settings = data.frame(
  type = c("Spatial", "Spatiotemporal"),
  data_filename = c("spatial_only_", 
                    "spatiotemporal_"),
  subtitle = c("Spatial Only\n(One-CV Constraint)", 
               "Spatiotemporal\n(One-CV Constraint)")
)

####################
{
  png(filename = paste0(figure_dir, "Fig5_choke_spp.png"),
      width = 190, 
      height = 100, 
      units = "mm", 
      res = 500)
  
  par(mfrow = c(1, 3), 
      mar = c(3, 0, 3, 0), 
      oma = c(1, 12, 0, 1))
  
  for (irow in 1:2) {
    load(paste0(github_dir, plot_settings$type[irow], 
                "_Optimization/STRS_Sim_Res_",
                plot_settings$type[irow], ".RData"))
    load(paste0(github_dir,plot_settings$type[irow],"_Optimization/",
                plot_settings$data_filename[irow], 
                "optimization_results.RData"))
    
    sub_settings <- subset(settings, 
                           nstrata == 10)
    sample_idx <- which.min(abs(sub_settings$n - 550))
    
    #Subset out the 10-strata, 2 boat scenario
    abs_diff <- STRS_true_cv_array[,, 2, 2] - sub_settings$cv[sample_idx]
    rel_diff <- 100 * abs_diff / sub_settings$cv[sample_idx]
    
    boxplot(rel_diff, 
            horizontal = TRUE, 
            add = F, 
            axes = F,
            pch = 16, 
            cex = 0.5, 
            ylim = c(-110,160),
            main = plot_settings$subtitle[irow])
    box()
    abline(v = 0, 
           col = "darkgrey", 
           lty = "dashed")
    axis(side = 1)
    if (irow == 1) axis(side = 2, 
                        labels = sci_names, 
                        las = 1, 
                        font = 3, 
                        at = 1:ns)
    
  }
  mtext(side = 1, 
        text = paste0("Percent Difference of the True CV ",
                      "Relative to the Upper CV Constraint"), 
        outer = T, 
        line = 0)
  
  #######################
  ## Flexible Scheme
  #######################
  load(paste0(github_dir, "Spatiotemporal_Optimization_Scheme2/",
              "spatiotemporal_Flexible_optimization_results.RData"))
  load(paste0(github_dir, "Spatiotemporal_Optimization_Scheme2/",
              "STRS_Sim_Res_Spatiotemporal_Flexible.RData"))
  settings$id <- 1:nrow(settings)
  
  sub_settings <- subset(settings, 
                         nstrata == 10)
  idx <- which.min(abs(sub_settings$n - 550))
  
  expected_cv <- unlist(sub_settings[idx - 1, paste0("CV_", 1:ns)] )
  expected_cv <- sapply(expected_cv, function(x) max(x, 0.1))
  
  #Subset out the 10-strata, 2 boat scenario
  abs_diff <- sweep(x = STRS_true_cv_array[,, 2, 2], 
                    MARGIN = 2, 
                    STATS = expected_cv, 
                    FUN = "-" )
  rel_diff <- 100 * sweep(x = abs_diff, 
                          MARGIN = 2, 
                          STATS = expected_cv, 
                          FUN = "/" )
  
  boxplot(rel_diff, 
          horizontal = TRUE, 
          add = F, 
          axes = F,
          pch = 16, 
          cex = 0.5, 
          ylim = c(-110,160),
          main =  "Spatiotemporal\n(Spp-Specific CV Constraints)",
          border = ifelse(expected_cv == 0.1, "blue", "black" ))
  box()
  abline(v = 0, 
         col = "darkgrey", 
         lty = "dashed")
  axis(side = 1)
  
  dev.off()
  
}

############################################
## Supplementary Figure
############################################
stratas = c(5, 10, 15, 20, 30, 60)
which_strata = c(2,4:6)

{
  png(filename = paste0(figure_dir, "Supplemental_Figures/SFig3_choke_spp.png"),
      width = 190, 
      height = 200, 
      units = "mm", 
      res = 500)
  
  par(mfrow = c(4, 3), 
      mar = c(0, 0, 0, 0), 
      oma = c(3.5, 12, 3.5, 1))
  
  for (istrata in which_strata) {
    for (irow in 1:2){
      load(paste0(github_dir, plot_settings$type[irow], 
                  "_Optimization/STRS_Sim_Res_",
                  plot_settings$type[irow], ".RData"))
      load(paste0(github_dir,plot_settings$type[irow],"_Optimization/",
                  plot_settings$data_filename[irow], 
                  "optimization_results.RData"))
      
      sub_settings = subset(settings, 
                            strata == stratas[istrata])
      sample_idx = which.min(abs(sub_settings$n - 550))
      
      abs_diff = STRS_true_cv_array[,,istrata,2] - sub_settings$cv[sample_idx]
      rel_diff = 100 * abs_diff / sub_settings$cv[sample_idx]
      
      boxplot(rel_diff, 
              horizontal = TRUE, 
              add = F, 
              axes = F,
              pch = 16, 
              cex = 0.5, 
              ylim = c(-110,160) )
      if (istrata == 2) mtext(side = 3, 
                              text = plot_settings$subtitle[irow], 
                              line = 0.5)
      box()
      abline(v = 0, 
             col = "darkgrey", 
             lty = "dashed")
      if (istrata == 6) axis(side = 1)
      if (irow == 1) axis(side = 2, 
                          lebels = sci_names, 
                          las = 1, 
                          font = 3, 
                          at = 1:ns)
      
    }
    mtext(side = 1, 
          text = paste0("Percent Difference of the True CV ",
                        "Relative to the Upper CV Constraint"), 
          outer = T, 
          line = 2.5)
    
    #######################
    ## Flexible Scheme
    #######################
    load(paste0(github_dir, "Spatiotemporal_Optimization_Scheme2/",
                "spatiotemporal_Flexible_optimization_results.RData"))
    load(paste0(github_dir, "Spatiotemporal_Optimization_Scheme2/",
                "STRS_Sim_Res_Spatiotemporal_Flexible.RData"))
    
    sub_settings = subset(settings, 
                          strata == stratas[istrata])
    idx = which.min(abs(sub_settings$n - 550))
    expected_cv = unlist(sub_settings[idx-1, paste0("CV_", 1:ns)])
    expected_cv = sapply(X = expected_cv, 
                         FUN = function(x) max(x, 0.1))
    abs_diff = sweep(x = STRS_true_cv_array[,, istrata, 2], 
                     MARGIN = 2, 
                     STATS = expected_cv, 
                     FUN = "-" )
    rel_diff = 100 * sweep(x = abs_diff, 
                           MARGIN = 2, 
                           STATS = expected_cv, 
                           FUN = "/" )
    
    boxplot(rel_diff, 
            horizontal = TRUE, 
            add = F, 
            axes = F,
            pch = 16, 
            cex = 0.5, 
            ylim = c(-110,160),
            border = ifelse(expected_cv == 0.1, "blue", "black" ))
    box()
    abline(v = 0, 
           col = "darkgrey", 
           lty = "dashed")
    if (istrata == 6) 
      axis(side = 1) 
    if (istrata == 2) 
      mtext(side = 3, 
            text = "Spatiotemporal\n(Spp-Specific CV Constraint)",
            line = 0.5)
    
    legend("bottomright", 
           legend = paste(stratas[istrata], "Strata"),
           bty = "n", 
           cex = 1.75)
  }
  dev.off()
}

