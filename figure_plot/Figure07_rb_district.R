###############################################################################
## Project:       Sensitivity: Observation Error
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Compare Relative Bias across districts across species
##                for current and optimized STRS surveys
###############################################################################
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
  file_name <- paste0(github_dir, "results/pred_dens_surveys/",
                      scen_name, "simulation_results.RData")
  
  load(file_name)
}

##################################
## General layout of plots
##################################
gen_layout <- rbind(matrix(rep(17, 6), nrow = 1),
                    cbind(matrix(data = 1:15,
                                 ncol = 5,
                                 byrow = T),
                          matrix(rep(16, 3), ncol = 1)))

spp_idx_1 <- spp_idx_opt[1:7]
spp_idx_2 <- spp_idx_opt[8:14]
spp_idx_3 <- c(spp_idx_opt[15], spp_idx_eval[1:6])
spp_idx_4 <- spp_idx_eval[7:11]

for (which_spp in paste(1:4)) {
  
  ## Open png file
  png(filename = paste0(output_dir, "Figure07_RB_district_", which_spp, ".png"),
      width = 190,
      height = 220,
      units = "mm",
      res = 500)
  
  ## Set up plot layout
  par(mar = c(0, 0, 0,0), oma = c(1,5,0,0))
  plot_layout <- rbind( cbind(gen_layout + 17*0, gen_layout + 17*1),
                        cbind(gen_layout + 17*2, gen_layout + 17*3),
                        cbind(gen_layout + 17*4, gen_layout + 17*5),
                        cbind(gen_layout + 17*6, gen_layout + 17*7))
  
  layout(mat =  plot_layout, widths = c(rep(1, 11), 0.1))
  
  ## Loop over species
  for (ispp in get(paste0("spp_idx_", which_spp))) {
    ## Loop over three survey scenarios
    ## 1) Current
    ## 5) Gulf-wide optiization, 10 strata
    ## 2) District-level optimiztion, 3 strata per district
    for (irow in c(1, 4, 2) ) {
      ## subset result object based on survey scenario
      scen_name <- paste0("SUR_", scen$survey_type[irow], "_",
                          scen$domain[irow], "_STR_", scen$strata[irow], "_")
      # rb_district <- get(paste0(scen_name, "log_rb_district"))
      rb_district <- get(paste0(scen_name, "rb_district"))
      
      ## Loop over district
      for (idistrict in 1:5) {
        
        ## Base plot
        plot(1,
             type = "n",
             xlim = c(0, 12),
             ylim = c(-100, 100),
             # ylim = log10(c(0.05, 2000)),
             axes = F,
             ann = F)
        box()
        
        ## Color of the background determines the survey type
        # rect(xleft = -5, 
        #      xright = 15, 
        #      ybottom = -2, 
        #      ytop = 5, 
        #      col = ifelse(irow %in% 1, 
        #                   "white",
        #                   ifelse(irow %in% 4, "grey90", 
        #                          "grey50")))
        
        ## District Labels
        if(irow == 1) mtext(side = 3, 
                            line = -0.75,
                            text = districts$district[idistrict],
                            cex = 0.5)
        
        ## Time axis
        axis(side = 1, 
             at = 1:11, 
             labels = NA,
             tck = -0.05)
        
        if(idistrict %in% c(2, 4) ) {
          axis(side = 1, 
               at = c(1, 11), 
               labels = c("Yr 1", "Yr 11"),
               lwd = F, 
               tick = F, 
               line = -1, 
               cex.axis = 0.75)
          arrows(x0 = 3.5, x1 = 8.5, 
                 y0 = -2, 
                 xpd = NA, length = 0.02)
        }
        
        ## Plot time series
        plot_this <- rb_district[,
                                 ispp,
                                 "boat_2" ,
                                 idistrict,
                                 ]
        
        plot_percentiles(values = plot_this,
                         xs = 1:11, 
                         pt.cex = 0.25,
                         inner_color = "dodgerblue",
                         pt.colors = "black")
        
        ## 
        if(idistrict == 1) {
          axis(side = 2, 
               las = 1, 
               # at = log10(c(0.1, 1, 10, 100, 1000)),
               at = c(-100, -50, 0, 50, 100),
               labels = NA,
               cex.axis = 0.7,
               tck = -0.1)
          
          axis(side = 2, 
               las = 1, 
               # at = log10(c(0.1, 1, 10, 100, 1000)),
               at = c(-100, -50, 0, 50, 100),
               labels = c(-100, -50, 0, 50, 100),
               # labels = c(0.1, 1, 10, 100, 1000),
               cex.axis = 0.75,
               lwd = 0,
               line = -0.25)
          
          axis(side = 2, 
               las = 1, 
               at = c(-100, -50, 0, 50, 100),
               # at = log10(c(0.05, 0.5, 5, 50, 500)),
               labels = NA,
               cex.axis = 0.7,
               tck = -0.05)
        }
        
        abline(h = 0, lwd = 0.5, lty = "dotted")
      }
    }
    
    plot(1,type = "n", axes = F, ann = F)  
    plot(1,type = "n", axes = F, ann = F, xlim = c(0, 1), ylim = c(0, 1))
    text(x = 0.4, 
         y = 0.35, 
         labels = common_names_all[ispp], 
         cex = 1.4, 
         font = 2)
  }
  
  # mtext(side = 2, 
  #       outer = T, 
  #       text = "Bias Ratio (Simulated Value / True Value)", 
  #       line = 3, 
  #       font = 2)
  
  dev.off()
}
