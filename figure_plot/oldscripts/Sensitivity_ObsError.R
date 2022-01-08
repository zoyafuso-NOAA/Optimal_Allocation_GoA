###############################################################################
## Project:       Sensitivity: Observation Error
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Compare Relative Bias, True CV, and RRMSE of CV across species
##                for surveys with varying added observation error for 
##                two-boat optimized surveys
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2)[1]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")
result_dir <- paste0(c("/Users/zackoyafuso/", 
                       "C:/Users/Zack Oyafuso/")[which_machine],
                     "Google Drive/MS_Optimizations/TechMemo/figures")

##################################################
####   Load Data
##################################################
load(paste0(github_dir, "data/optimization_data.RData"))
load(paste0(github_dir, "results/MS_optimization_knitted_results.RData"))
source(paste0(github_dir, "modified_functions/plot_percentiles.R"))

scen <- data.frame(survey_type = c("cur", "opt", "opt"),
                   strata = c("cur", 5, 15),
                   domain = c("full_domain", "district", "full_domain"))

for (irow in 1:nrow(scen)) {
  scen_name <- paste0("SUR_", scen$survey_type[irow], "_",
                      scen$domain[irow], "_STR_", scen$strata[irow], "_")
  file_name <- paste0(github_dir, "results/",
                      scen_name, "simulation_results.RData")
  
  load(file_name)
}

##################################################
####   Plot
##################################################
spp_set = "opt"

irow = 2

rb <- with(scen[irow, ], 
           get(paste0("SUR_", survey_type, "_", domain, 
                      "_STR_", strata, "_rb_agg")))
true_cv <- with(scen[irow, ], 
                get(paste0("SUR_", survey_type, "_", domain, 
                           "_STR_", strata, "_true_cv")))
rrmse_cv <- with(scen[irow, ], 
                 get(paste0("SUR_", survey_type, "_", domain, 
                            "_STR_", strata, "_rrmse_cv")))



for (spp_set in  c("opt", "eval")) {
  
  ## Set up png file
  # png(filename = paste0(result_dir, 
  #                       "/Sensitivity_ObsError_", spp_set,
  #                       ".png"),
  #     width = 190,
  #     height = c("opt" = 220, "eval" = 120)[spp_set],
  #     units = "mm",
  #     res = 500)
  
  ## Set up plot layout based on which species to plot
  layout(mat = matrix(data = list(
    "opt" = c(1, 2, 3,       25, 26, 27,  
              4, 5, 6,       28, 29, 30,
              7, 8, 9,       31, 32, 33,
              10, 11, 12,    34, 35, 36,
              13, 14, 15,    37, 38, 39,
              16, 17, 18,    40, 41, 42,
              19, 20, 21,    43, 44, 45,
              22, 23, 24,    46, 46, 46),
    "eval" = c(1, 2, 3,       13, 14, 15,  
               4, 5, 6,       16, 17, 18,   
               7, 8, 9,       19, 20, 21,    
               10, 11, 12,    22, 22, 22))[[spp_set]],
    ncol = 6,
    byrow = T),
    widths = rep(x = c(1,0.5, 0.5), times = 3) )
  par(mar = c(0.25, 1.75, 0.25, 1.75), 
      oma = c(0, 1, 3, 0.25))
  
  ## Set up constants based on which species to plot
  ns = get(paste0("ns_", spp_set))
  spp_idx <- get(paste0("spp_idx_", spp_set))
  common_names <- get(paste0("common_names_", spp_set))
  
  for (ispp in 1:ns) {
    ## Relative Bias
    par(mar = c(0.25, 1.75, 0.25, 1.75))
    
    opt_ylim <- 1.0 * max(abs(apply( X = rb[paste0("obsCV=", c(0, 1)), 
                                            1:n_years, 
                                            spp_idx[ispp], 
                                            "boat_2",
                                            ],
                                     MARGIN = 1:2,
                                     FUN = quantile,
                                     probs = c(0.025, 0.975),
                                     na.rm = T))) 
    
    plot(1,  
         type = "n",
         ylim = c(-opt_ylim, opt_ylim),
         xlim = c(1, 23),
         axes = F,
         ann = F)
    
    # if (ispp %in% c(1, c("opt" = 9, "eval" = 5)[spp_set])) {
    #   mtext(side = 3, 
    #         text = "Percent\nBias", 
    #         font = 2)
    # }
    
    plot_percentiles(values = rb[paste0("obsCV=0"), 
                                 1:n_years, 
                                 spp_idx[ispp], 
                                 "boat_2",
                                 ],
                     xs = 1:11,
                     inner_color = "cadetblue1",
                     outer_color = "dodgerblue")
    
    plot_percentiles(values = rb["obsCV=1", 
                                 1:n_years, 
                                 spp_idx[ispp], 
                                 "boat_2",
                                 ],
                     xs = 13:23,
                     inner_color = "chartreuse3",
                     outer_color = "palegreen4")
    
    abline(h = 0)
    axis(side = 2, las = 1)
    box()
    
    ## Species Label
    mtext(side = 1,
          line = -1,
          text = common_names[ispp],
          cex = 0.7)
    
    ##  True CV
    par(mar = c(1.75, 1.5, 0.5, 1.25))
    boxplot(t(true_cv[paste0("obsCV=", c(0, 0.25, 0.5, 1)),
                      1:n_years,
                      spp_idx[ispp], 
                      "boat_2"]),
            xlim = c(-1, 8),
            ylim = c(0, 
                     1.25 * 
                       max(true_cv[paste0("obsCV=", c(0, 0.25, 0.5, 1)),
                                   1:n_years,
                                   spp_idx[ispp],
                                   "boat_2"])),
            las = 1,
            lwd = 0.5,
            col = c("dodgerblue", 
                    "white", "white",
                    "chartreuse3"),
            axes = F,
            pch = 16,
            cex = 0.25,
            at = seq(from = 0.5, to = 6.5, by = 2))
    box()
    axis(side = 2, 
         las = 1,
         cex.axis = 0.9)
    axis(side = 1,
         at = seq(from = 0.5, to = 6.5, by = 2),
         labels = NA, 
         tcl = -0.15)
    text(x = seq(from = 1.25, to = 7.25, by = 2),
         y = max(true_cv[paste0("obsCV=", c(0, 0.25, 0.5, 1)),
                         1:n_years,
                         spp_idx[ispp], 
                         "boat_2"]) * 
           -c(0.175, 0.35)[c(1, 2, 1, 2, 1)] ,
         labels = paste0(c(0, 25, 50, 100), "%"),
         cex = 0.7,
         col = c("blue", 
                 "black", "black",
                 "chartreuse3"),
         xpd = NA)
    
    # if (ispp %in% c(1, c("opt" = 9, "eval" = 5)[spp_set])) {
    #   mtext(side = 3, 
    #         "True\nCV", 
    #         font = 2, 
    #         line = 0.25)}
    
    ## RRMSE of CV
    
    boxplot(t(rrmse_cv[paste0("obsCV=", c(0, 0.25, 0.5, 1)),
                       1:n_years,
                       spp_idx[ispp],
                       "boat_2"]),
            xlim = c(-1, 8),
            ylim = c(0, 
                     1.25 * max(
                       rrmse_cv[paste0("obsCV=", 
                                       c(0, 0.25, 0.5, 1)),
                                1:n_years,
                                spp_idx[ispp],
                                "boat_2"])),
            las = 1,
            axes = F,
            lwd = 0.5,
            col = c("dodgerblue", 
                    "white", "white", 
                    "chartreuse3"),
            pch = 16,
            cex = 0.25,
            at = seq(from = 0.5, to = 6.5, by = 2))
    box()
    axis(side = 2, 
         las = 1, 
         cex.axis = 0.9)
    axis(side = 1,
         at = seq(from = 0.5, to = 6.5, by = 2),
         labels = NA, 
         tcl = -0.15)
    
    ## CV Label
    text(x = seq(from = 1.25, to = 7.25, by = 2),
         y = c(-0.175, -0.35)[c(1, 2, 1, 2, 1)] * 
           max(rrmse_cv[,
                        1:n_years,
                        spp_idx[ispp],
                        "boat_2"]),
         labels = paste0(c(0, 25, 50, 100), "%"),
         cex = 0.7,
         col = c("blue", 
                 "black", "black",
                 "darkgreen"),
         xpd = NA)
    
    # if (ispp %in% c(1, c("opt" = 9, "eval" = 5)[spp_set])) {
    #   mtext(side = 3, 
    #         text = "RRMSE\nCV", 
    #         font = 2, 
    #         line = 0.25)}
  }
  
  ## Legend
  par(mar = c(0, 1, 0, 2))
  plot(1,
       xlim = c(0, 1),
       ylim = c(0, 1),
       type = "n",
       axes = F,
       ann = F)
  
  ## Percentile legend
  rect(xleft = c(0.1, 0.1, 0.4, 0.4), 
       xright = c(0.2, 0.2, 0.5, 0.5),
       ybottom = c(0.0, 0.3, 0.0, 0.3), 
       ytop = c(0.8, 0.5, 0.8, 0.5),
       col = c("dodgerblue", "cadetblue1", "chartreuse3", "palegreen4"),
       border = F)
  points(x = c(0.15, 0.45), 
         y = c(0.4, 0.4), 
         pch = 15)
  segments(x0 = c(0.15, 0.45),
           y0 = 0.0,
           y1 = 0.8,
           lty = "dashed")
  text(x = 0.3,
       y = c(0.025, 0.275, 0.4, 0.525, 0.75),
       labels = paste0(c(5,25,50,75,95), "%"))
  text(x = c(0.15, 0.45),
       y = c(0.925, 0.925),
       labels = c("Pred.\nDensity", "+Obs Error\n(100% CV)") )
  
  ## Years Legend
  points(x = seq(from = 0.6, to = 1, length = 11),
         y = rep(0.5, 11),
         col = rev(grey.colors(11, start = 0, end = 0.9)),
         pch = 15,
         cex = 2)
  text(x = c(0.6, 1),
       y = 0.65,
       labels = c(1996, 2019) )
  
  # dev.off()
}

