###############################################################################
## Project:       Sensitivity: Observation Error
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Compare Relative Bias, True CV, and RRMSE of CV across species
##                for surveys where observation error was and was not included
##                when simulating surveys
##                2-boat survey only
###############################################################################
rm(list = ls())

##################################################
####   Import Libraries
##################################################
library(RColorBrewer)

##################################################
####   Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2)[2]

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
load(paste0(github_dir, "results/simulation_result.RData"))

##################################################
####   Relative Bias Plots:
####   differences in simulation type
####
####   Base Plot
##################################################
plot_percentiles <- function(
  values = NULL,
  xs = NULL,
  ispp = NULL,
  plot_years = F,
  inner_color = "cadetblue1",
  outer_color = "dodgerblue"
){
  
  if(is.null(values)) stop("Must supply values")
  
  temp_quants <- apply(X = values,
                       MARGIN = 1,
                       FUN = quantile,
                       probs = c(0.03, 0.05, 0.25, 0.50, 0.75, 0.95, 0.97),
                       na.rm = T)
  
  polygon(y = c(temp_quants["5%",], rev(temp_quants["95%",])),
          x = c(xs, rev(xs)),
          col = outer_color,
          border = F)
  polygon(y = c(temp_quants["25%",], rev(temp_quants["75%",])),
          x = c(xs, rev(xs)),
          col = inner_color,
          border = F)
  segments(y0 = temp_quants["5%",],
           y1 = temp_quants["95%",],
           x0 = xs,
           col = "black",
           lty = "dashed")
  abline(h = 0)
  
  points(y = temp_quants["50%", ],
         x = xs,
         pch = 15,
         col = rev(grey.colors(11, start = 0, end = 0.9)) )
  
  if(plot_years){
    text(x = seq(from = xs[1], to = xs[length(xs)], by = 2), 
         y = temp_quants["3%", seq(from = 1, to = 11, by = 2)],
         labels = c(1996, seq(from = 2003, to = 2019, by = 4)) )
  }
  
}

##################################################
####   Plot
##################################################
{
  png(filename = paste0(result_dir, "/Comparison_trawlable.png"),
      width = 190,
      height = 220,
      units = "mm",
      res = 500)
  
  layout(mat = matrix(c(1, 2, 3,       25, 26, 27,  
                        4, 5, 6,       28, 29, 30,
                        7, 8, 9,       31, 32, 33,
                        10, 11, 12,    34, 35, 36,
                        13, 14, 15,    37, 38, 39,
                        16, 17, 18,    40, 41, 42,
                        19, 20, 21,    43, 44, 45,
                        22, 23, 24,    46, 46, 46),
                      ncol = 6,
                      byrow = T),
         widths = rep(c(1,0.5, 0.5), times = 3) )
  par(mar = c(0.25, 1.75, 0.25, 1.75), oma = c(0, 1, 3, 1))
  
  for (isample in 2) {
    for (ispp_ in 1:ns) {
      
      #Relative Bias
      ylim_ <- max(abs(
        apply(X = cbind(STRS_rel_bias_est_trawl[1, , ispp_, isample, ], 
                        Current_rel_bias_est_trawl[1, , ispp_, isample, ]),
              MARGIN = 1,
              FUN = quantile,
              probs = c(0.025, 0.975),
              na.rm = T)
      ))
      
      par(mar = c(0.25, 1.75, 0.25, 1.75))
      plot(1,  
           type = "n",
           ylim = c(-ylim_, ylim_),
           xlim = c(1, 23),
           axes = F,
           ann = F)
      
      if(ispp_ %in% c(1, 9)) mtext(side = 3, "Percent\nBias", font = 2)
      
      plot_percentiles(values = STRS_rel_bias_est_trawl[1, , ispp_, isample,],
                       xs = 1:11, 
                       ispp = ispp_,
                       plot_years = F,
                       inner_color = "cadetblue1",
                       outer_color = "dodgerblue")
      
      plot_percentiles(values = Current_rel_bias_est_trawl[1, , ispp_, isample,],
                       xs = 13:23, 
                       ispp = ispp_,
                       plot_years = F,
                       inner_color = "tomato2",
                       outer_color =  "firebrick" )
      axis(side = 2, las = 1)
      box()
      
      #True CV
      par(mar = c(2, 1.75, 0.25, 1.5))
      ylim_ <- max(abs(cbind(STRS_true_cv_array_trawl[1, , ispp_, isample],
                             Current_true_cv_array_trawl[1, , ispp_, isample])))
      
      boxplot(cbind(STRS_true_cv_array_trawl[1, , ispp_, isample],
                    Current_true_cv_array_trawl[1, , ispp_, isample]),
              ylim = c(0, 1.25 * ylim_),
              las = 1,
              axes = F,
              pch = 16,
              col  = c("dodgerblue", "firebrick"))
      box()
      axis(side = 2, las = 1)
      if(ispp_ %in% c(1, 9)) mtext(side = 3, "True\nCV", font = 2)
      
      #RRMSE of CV
      ylim_ <- max(abs(cbind(STRS_rrmse_cv_array_trawl[1, , ispp_, isample],
                             Current_rrmse_cv_array_trawl[1, , ispp_, isample])))
      
      boxplot(cbind(STRS_rrmse_cv_array_trawl[1, , ispp_, isample],
                    Current_rrmse_cv_array_trawl[1, , ispp_, isample]),
              ylim = c(0, 1.25 * ylim_),
              las = 1,
              names = NA,
              axes = F,
              pch = 16,
              col = c("dodgerblue", "firebrick"))
      box()
      axis(side = 2, las = 1)
      if(ispp_ %in% c(1, 9)) mtext(side = 3, "RRMSE\nCV", font = 2)
      
      text(x = -1, 
           y = par("usr")[3] - diff(par("usr")[3:4])*0.15,
           labels = common_names[ispp_],
           xpd = NA,
           font = 3,
           cex = 1.25)
      
    }
  }
  
  ##Legend
  par(mar = c(0, 1, 0, 2))
  plot(1,
       xlim = c(0, 1),
       ylim = c(0, 1),
       type = "n",
       axes = F,
       ann = F)
  
  #Percentile legend
  rect(xleft = c(0.1, 0.1, 0.4, 0.4), xright = c(0.2, 0.2, 0.5, 0.5),
       ybottom = c(0.0, 0.3, 0.0, 0.3), ytop = c(0.8, 0.5, 0.8, 0.5),
       col = c("dodgerblue", "cadetblue1", "firebrick", "tomato2"),
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
       labels = c("Opt.\nSurvey", "Curr.\nSurvey") )
  
  #Years Legend
  points(x = seq(from = 0.6, to = 1, length = 11),
         y = rep(0.5, 11),
         col = rev(grey.colors(11, start = 0, end = 0.9)),
         pch = 15,
         cex = 2)
  text(x = c(0.6, 1),
       y = 0.65,
       labels = c(1996, 2019) )
  
  dev.off()
}

