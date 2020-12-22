rm(list = ls())

##################################################
####  Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")

output_dir <- paste0(c("C:/Users/Zack Oyafuso/"),
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

scen <- data.frame(survey_type = c("cur", rep("opt", 6) ),
                   strata = c("cur", 3, 5, 10, 10, 15, 20),
                   domain = c("full_domain", 
                              rep(c("district", "full_domain"), each = 3)))

for (irow in 1:7) {
  scen_name <- paste0("SUR_", scen$survey_type[irow], "_", 
                      scen$domain[irow], "_STR_", scen$strata[irow], "_")
  file_name <- paste0(github_dir, "results/", 
                      scen_name, "simulation_results.RData")
  
  load(file_name)
  
}

{png(filename = paste0(output_dir, "RB_full_domain.png"), 
     width = 190, 
     height = 200,
     units = "mm", 
     res = 500)
  
par(mar = c(0.5, 0, 0.5, 0),
    mfcol = c(6, 7),
    oma = c(0, 5, 5, 0))

for (irow in 1:7) {
  
  scen_name <- paste0("SUR_", scen$survey_type[irow], "_",
                      scen$domain[irow], "_STR_", scen$strata[irow], "_")
  rb_agg <- get(paste0(scen_name, "rb_agg"))
  
  for (ispp in c(2, 5, 14, 18, 15, 21)) {#, spp_idx_eval)) {
    
    ## Inner lapply: knits together all the relative bias arrays into a list
    ##               and extract the 0 obs CV, species ispp, 2-boat soln's
    ## Outer lapply: calculate 5% and 95% percentiles on each sub-list across
    ##               simulated surveys
    ## Then calculate the most extreme values and use as symmetrical y-limits
    y_range <- range(unlist(
      lapply(X = lapply(X = grep(x = ls(), 
                                 pattern = "_rb_agg", 
                                 value = TRUE), 
                        FUN = function(x) get(x)["obsCV=0", 
                                                 ,
                                                 ispp     ,
                                                 "boat_2" ,
                                                 ]),
             
             FUN = function(x) apply(x, MARGIN = 1, 
                                     FUN = quantile, 
                                     probs = c(0.05, 0.95), 
                                     na.rm = T))
    ))
    
    plot(1, 
         type = "n",
         xlim = c(1, 11),
         ylim = y_range,
         axes = F,
         ann = F)
    box()
    
    if (irow == 1) axis(side = 2, las = 1)
    
    if (ispp == 2)
      mtext(side = 3, 
            text = with(scen[irow,], 
                        paste0(survey_type, " Survey\n", domain, "\n",
                               strata, " Strata")))
    
    plot_percentiles(values = rb_agg["obsCV=0", 
                                     ,
                                     sci_names_all[ispp], 
                                     "boat_2" ,
                                     ],
                     xs = 1:11)
    abline(h = 0)
    
    if (irow == 4)
      mtext(side = 3, 
            text = gsub(x = sci_names_all[ispp], 
                        pattern = " ", 
                        replacement = "\n"), 
            line = -4)
  }
  
}

mtext(side = 2,
      outer = T,
      line = 3,
      text = "Relative Percent Bias (100% (Sim - True) / True)")

dev.off()
}



