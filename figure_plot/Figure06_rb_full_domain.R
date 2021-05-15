##############################################################################
# Project:       Sensitivity: Observation Error
# Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
# Description:   Compare Relative Bias across species
#                for current and optimized STRS surveys
##############################################################################
rm(list = ls())

##################################################
####  Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/",
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine],
                     "GitHub/Optimal_Allocation_GoA/")

output_dir <- paste0(c("/Users/zackoyafuso/",
                       "C:/Users/Zack Oyafuso/")[which_machine],
                     "Google Drive/MS_Optimizations/TechMemo/figures/")

##################################
## Import plotting percentile time series function
##################################
source(paste0(github_dir, "modified_functions/plot_percentiles.R"))

##################################
## Import Strata Allocations and spatial grid and predicted density
##################################
load(paste0(github_dir, '/data/optimization_data.RData'))
load(paste0(github_dir, '/data/Extrapolation_depths.RData'))

scen <- data.frame(survey_type = c("cur", rep("opt", 4) ),
                   strata = c("cur", 3, 5, 10, 15),
                   domain = c("full_domain",
                              rep(c("district", "full_domain"), each = 2)))

for (irow in 1:nrow(scen)) {
  scen_name <- paste0("SUR_", scen$survey_type[irow], "_",
                      scen$domain[irow], "_STR_", scen$strata[irow], "_")
  file_name <- paste0(github_dir, "results/sim_dens_surveys/",
                      scen_name, "simulation_results.RData")
  
  load(file_name)
}

##################################
## General layout of plots
##################################
gen_layout <- matrix(c(5, 1,2,3, 4,4,4,4), ncol = 2)

{
  for (spp_group in 1:2) { ## loop over species groupings (opt or eval) -- start
    
    ## Open Device
    png(filename = paste0(output_dir, "Figure6", LETTERS[spp_group],
                          "_RB_full_domain.png"),
        width = 190, height = c(220, 170)[spp_group], units = "mm",
        res = 500)
    
    ## Plot Layout
    par(mar = c(.5, 0, 0.25, 0), oma = c(1,5,0,0))
    plot_layout <- rbind(
      cbind(gen_layout + 5 * 0, gen_layout + 5 * 1, 
            gen_layout + 5 * 2, gen_layout + 5 * 3),
      cbind(gen_layout + 5 * 4, gen_layout + 5 * 5, 
            gen_layout + 5 * 6, gen_layout + 5 * 7), 
      cbind(gen_layout + 5 * 8, gen_layout + 5 * 9, 
            gen_layout + 5 * 10, gen_layout + 5 * 11),
      switch(paste0(spp_group),
             "1" = cbind(gen_layout + 5 * 12, gen_layout + 5 * 13, 
                         gen_layout + 5 * 14, gen_layout + 5 * 15),
             "2" = NULL)
    )
    layout(mat =  plot_layout, 
           widths = c(1,0.4, 1,0.4, 1,0.4, 1,0.1),
           heights = c(0.5, 1, 1, 1))
    
    ## Loop over species
    for (ispp in list(spp_idx_opt, spp_idx_eval)[[spp_group]] ) {
      ## Loop over three survey scenarios
      ## 1) Current
      ## 5) Gulf-wide optiization, 10 strata
      ## 2) District-level optimiztion, 3 strata per district 
      ## start
      
      ## Calculate y-max for the plot
      y_max <- max(abs(unlist(
        lapply(X = lapply(X = c(1, 5, 2),
                          FUN = function(x) 
                            get(paste0("SUR_", 
                                       scen$survey_type[irow], 
                                       "_",
                                       scen$domain[irow], 
                                       "_STR_", 
                                       scen$strata[irow], 
                                       "_rb_agg"))[,
                                                   ispp,
                                                   "boat_2" ,
                                                   ] ),
               FUN = function(x) plot_percentiles(values = x, plot = F) )))) 
      
      y_max <- max(y_max, 25)
      
      for (irow in c(1, 5, 2) ) {
        ## subset result object based on survey scenario
        scen_name <- paste0("SUR_", scen$survey_type[irow], "_",
                            scen$domain[irow], "_STR_", scen$strata[irow], "_")
        rb_agg <- get(paste0(scen_name, "rb_agg"))[ , ispp, "boat_2" , ]
        
        ## Base plot
        plot(x = 1, y = 1,
             type = "n", axes = F, ann = F,
             xlim = c(0, 12), ylim = c(-y_max, y_max))
        
        ## Color of the background corresponds to the  the survey type
        rect(xleft = -5, 
             xright = 15, 
             ybottom = -500, 
             ytop = 500, 
             col = ifelse(irow %in% 1, 
                          "white",
                          ifelse(irow %in% 5, "grey90", 
                                 "grey50")))
        axis(side = 2,
             las = 1,
             cex.axis = 0.75,
             at = pretty(x = c(-y_max, y_max), n = 3) )
        box()
        
        ## Time axis
        axis(side = 1,
             at = 1:11,
             labels = NA,
             tck = -0.05)
        
        if(irow == 2) {
          axis(side = 1, 
               at = c(1, 11), 
               labels = c("Yr 1", "Yr 11"),
               lwd = F, 
               tick = F, 
               line = -0.5, 
               cex.axis = 1)
        }
        
        ## Plot bias distributions
        plot_percentiles(values = rb_agg,
                         xs = 1:11, 
                         pt.cex = 0.5,
                         pt.colors = "black")
        abline(h = 0, lwd = 0.5, lty = "dotted")
      }
      
      plot(1,type = "n", axes = F, ann = F)  
      plot(1,type = "n", axes = F, ann = F, xlim = c(0, 1), ylim = c(0, 1))
      text(x = 0.5,
           y = 0.3, 
           labels = common_names_all[ispp], 
           cex = 1.25, 
           font = 2, 
           xpd = NA)
    } ## Loop over three survey scenarios -- end
    
    ## Plot Legend
    plot(x = 1, y = 1,
         type = "n", axes = F, ann = F, 
         xlim = c(0, 1), ylim = c(0, 1))
    rect(xleft = par("usr")[1], xright = par("usr")[3], 
         ybottom = par("usr")[2], ytop = par("usr")[4],
         col = "white")
    text(x = 0.5,
         y = 0.5, 
         labels = "Current STRS Design", 
         cex = 1, 
         font = 2, 
         xpd = NA)
    
    plot(x = 1, y = 1,
         type = "n", axes = F, ann = F, 
         xlim = c(0, 1), ylim = c(0, 1))
    rect(xleft = par("usr")[1], xright = par("usr")[3], 
         ybottom = par("usr")[2], ytop = par("usr")[4],
         col = "white")
    text(x = 0.5,
         y = 0.5, 
         labels = "Gulf-Wide\n(10 Strata)\nOptimized STRS Design",
         cex = 1, 
         font = 2, 
         xpd = NA)
    
    plot(x = 1, y = 1,
         type = "n", axes = F, ann = F, 
         xlim = c(0, 1), ylim = c(0, 1))
    rect(xleft = par("usr")[1], xright = par("usr")[3], 
         ybottom = par("usr")[2], ytop = par("usr")[4],
         col = "white")
    text(x = 0.5,
         y = 0.5, 
         labels = "District-Level\n(3 Strata per District)\nOptimized STRS Design", 
         cex = 1, 
         font = 2, 
         xpd = NA)
    
    mtext(side = 2, 
          outer = T, 
          text = "Percent Bias (100% (Sim - True) / True)", 
          line = 3, 
          font = 2)
    
    ## Close Device
    dev.off()
  } ## loop over species groupings (opt or eval) -- end
}
