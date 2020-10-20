###############################################################################
## Project:       Create Relatie bias plots
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Produce a function that takes in a distribution of relative
##                bias calculations for each year of data and produce a violin
##                plot
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[1]
VAST_model <- "11" 

github_dir <- paste0(c('/Users/zackoyafuso/Documents/', 
                       'C:/Users/Zack Oyafuso/Documents/',
                       'C:/Users/zack.oyafuso/Work/')[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/model_", VAST_model, "/")

##################################################
####   Load data
##################################################
load(paste0(github_dir, "/optimization_data.RData") )
load(paste0(github_dir, 
            "/Survey_Comparison_Simulations/",
            "Simple_RS_Simulation_Results.RData") )
load(paste0(github_dir, 
            "/Survey_Comparison_Simulations/",
            "Survey_Simulation_Results.RData") )

##################################################
####   Base Plot
##################################################
plot_percentiles <- function(
 values = NULL,
 xs = NULL,
 ispp = NULL,
 plot_years = F
){
 
 if(is.null(values)) stop("Must supply values")
 
 temp_rel_bias <- sweep(x = values,
                        MARGIN = 1,
                        STATS = true_mean[, ispp],
                        FUN = "-")
 temp_rel_bias <- sweep(x = temp_rel_bias,
                        MARGIN = 1,
                        STATS = true_mean[, ispp],
                        FUN = "/") * 100
 
 temp_quants <- apply(X = temp_rel_bias,
                      MARGIN = 1,
                      FUN = quantile,
                      probs = c(0.05, 0.25, 0.50, 0.75, 0.95))
 
 polygon(y = c(temp_quants["5%",], rev(temp_quants["95%",])),
         x = c(xs, rev(xs)),
         col = 'darkgrey',
         border = F)
 polygon(y = c(temp_quants["25%",], rev(temp_quants["75%",])),
         x = c(xs, rev(xs)),
         col = 'black')
 segments(y0 = temp_quants["5%",],
          y1 = temp_quants["95%",],
          x0 = xs,
          col = "lightgrey",
          lty = "dashed")
 abline(h = 0, col = 'purple', lty = 'dashed')
 
 matpoints(x = xs,
           y = t(temp_quants[c("5%", "25%", "75%", "95%"), ]),
           pch = 16,
           col = 'white',
           cex = 0.5)
 
 points(y = temp_quants["50%", ],
        x = xs,
        pch = 15,
        col = 'darkgrey')
 
 if(plot_years){
  text(x = seq(from = xs[1], to = xs[length(xs)], by = 2), 
       y = temp_quants["5%", seq(from = 1, to = 11, by = 2)],
       labels = c(1996, seq(from = 2003, to = 2019, by = 4)) )
  text(x = seq(from = xs[2], to = xs[length(xs)], by = 2), 
       y = temp_quants["95%", seq(from = 2, to = 10, by = 2)] ,
       labels = c(1999, seq(from = 2005, to = 2019, by = 4)) )
 }
 
}

par(mfrow = c(4,2), mar = c(1,4,2,0))
isample <- 2
for(ispp_ in 1:15){
 
 ylim_ <- sweep(x = SRS_sim_mean[, ispp_, isample, ],
               MARGIN = 1,
               STATS = true_mean[, ispp_],
               FUN = "-")
 ylim_ <- sweep(x = ylim_,
               MARGIN = 1,
               STATS = true_mean[, ispp_],
               FUN = "/") * 100
 ylim_ <- apply(ylim_, 
               MARGIN = 1,
               FUN = quantile,
               probs = c(0.05, 0.95))
 ylim_ <- max(abs(ylim_)) 
 
 plot(1,  
      type = "n",
      ylim = 1.25 * c(-ylim_, ylim_),
      xlim = c(1,23),
      axes = F,
      ann = F)
 box()
 mtext(side = 3, text = sci_names[ispp_], font = 3)
 axis(side = 2, 
      at = c(-100, -75, -50, -20, -10, 0, 10, 20, 50, 75, 100),
      las = 1)
 axis(side = 2, at = c(-5, 5), labels = NA)
 
 plot_percentiles(values = SRS_sim_mean[, ispp_, isample, ],
                  xs = 1:11, 
                  ispp = ispp_,
                  plot_years = T)
 
 plot_percentiles(values = Survey_sim_mean[, ispp_, , isample],
                  xs = 13:23, 
                  ispp = ispp_,
                  plot_years = T)
}

