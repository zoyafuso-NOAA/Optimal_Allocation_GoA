
###############################################################################
## Project:         Lower-CV constraint on species CV
## Author:          Zack Oyafuso (zack.oyafuso@noaa.go)
## Description:     Rerun the Bethel algorithm on a given suvery design with
##                  different lower-cv constraints
###############################################################################

rm(list = ls())

##################################################
####  Set up directories  
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")
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
library(sp)
library(raster)
library(RColorBrewer)
library(plotrix)

##################################################
####  Load Data
##################################################
load(paste0(github_dir, 
            "data/optimization_data.RData"))
load(paste0(github_dir, 
            "data/Extrapolation_depths.RData"))
load(paste0(github_dir,
            "results/MS_optimization_knitted_results.RData"))
source( paste0(github_dir, "/modified_functions/sim_fns.R") )
load(paste0(github_dir, '/data/fit_density.RData'))


##################################################
####  Isolate District-Level optimization, 3 strata per district
####  Set lower threshold values
##################################################
sol_idx <- 2
threshold <- c(0.15, 0.25)

strata_df <- strata_stats_list[[sol_idx]]
nh <- strata_list[[sol_idx]]$Allocation
strata_df$stratum <- 1:nrow(strata_df)
strata_df$DOM1 <- 1

spp_order <- c(1:9, 15, 11, 10, 12:14) #For plotting

plot_solution <- res_df[, sol_idx]
goa <- sp::SpatialPointsDataFrame(
  coords = Extrapolation_depths[,c("E_km", "N_km")],
  data = data.frame(Str_no = plot_solution) )
goa_ras <- raster::raster(x = goa, 
                          resolution = 5)
goa_ras <- raster::rasterize(x = goa, 
                             y = goa_ras, 
                             field = "Str_no")

xrange <- range(Extrapolation_depths[, "E_km"])
yrange <- range(Extrapolation_depths[, "N_km"])
xrange_diff <- diff(xrange)
yrange_diff <- diff(yrange)

##################################################
####  Result Object
##################################################
tradeoff_df <- matrix(nrow = ns_opt, 
                      ncol = length(threshold))
tradeoff_sample_allocations <- matrix(nrow = nrow(strata_df), 
                                      ncol = 1 + length(threshold))
tradeoff_sample_allocations[, 1] <- nh

for (ithreshold in 1:length(threshold)) {
  
  ## Initiate CVs to be those calculated under SRS
  srs_stats <- SamplingStrata::buildStrataDF( 
    dataset = cbind( subset(frame_all, select = -c(X1, X2)),
                     X1 = 1))
  
  srs_var <- as.matrix(srs_stats[, paste0("S", spp_idx_opt)])^2
  
  srs_var <- sweep(x = srs_var, 
                   MARGIN = 1, 
                   STATS = (1 - samples[2] / n_cells) / samples[2], 
                   FUN = "*")
  
  srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", spp_idx_opt)]
  names(srs_cv) <- paste0("CV", 1:ns_opt)
  
  ##################################################
  ####  Calculate allocation at Simple Random Sampling Population CV
  ####  and 2 boats, set lower CV threshold
  ##################################################
  error_df <- cbind(data.frame(DOM = "DOM1"),
                    srs_cv,
                    data.frame(domainvalue = 1))
  
  error_df[, paste0("CV", 1:ns_opt)] <- 
    ifelse(test = error_df[, paste0("CV", 1:ns_opt)] <= threshold[ithreshold], 
           yes = threshold[ithreshold], 
           no = error_df[, paste0("CV", 1:ns_opt)])
  
  ##################################################
  ####  Run Bethel algorithm and calculate total sample size
  ##################################################
  bethel_allocation  <- bethel(errors = error_df, 
                               stratif = strata_df, 
                               printa=TRUE)
  pop_cvs <- as.numeric(attributes(bethel_allocation)$outcv[, "ACTUAL CV"])
  (n <- sum(bethel_allocation))
  
  ##################################################
  ####  Adjust CVs to be within 550 samples (desired under two boats)
  ##################################################
  over_or_under <- ifelse(n < 549, "under", 
                          ifelse(n > 551, "over", 
                                 "exact"))
  
  while(over_or_under %in% c("under", "over")) {
    change_rate <- switch(over_or_under, 
                          "under" = 0.99,
                          "over" = 1.01,
                          "exact" = 1)
    error_df[, paste0("CV", 1:ns_opt)] <- 
      pop_cvs * change_rate
    
    error_df[, paste0("CV", 1:ns_opt)] <- 
      ifelse(test = error_df[, paste0("CV", 1:ns_opt)] <= threshold[ithreshold], 
             yes = threshold[ithreshold], 
             no = error_df[, paste0("CV", 1:ns_opt)])
    
    bethel_allocation  <- bethel(errors = error_df, 
                                 stratif = strata_df, 
                                 printa=TRUE)
    
    pop_cvs <- as.numeric(attributes(bethel_allocation)$outcv[, "ACTUAL CV"])
    
    print(n <- sum(bethel_allocation))
    
    over_or_under <- ifelse(n < 549, "under", 
                            ifelse(n > 551, "over", 
                                   "exact"))
  }
  
  tradeoff_df[, ithreshold] <- pop_cvs
  tradeoff_sample_allocations[, 1 + ithreshold] <- as.numeric(bethel_allocation)
}

##################################################
####   
##################################################

temp_res <- array(dim = c(length(threshold) + 1, n_years, ns_all, n_iters) )

for (isurvey in 1:(length(threshold)+1) ) {
  for (isim in 1:n_iters) {
    
    set.seed(isim + 23423)
    sim_survey <- 
      do_STRS( 
        input = list(
          "density" = D_gct[, , years_included],
          
          "cell_areas" = Extrapolation_depths$Area_km2,
          
          "obs_CV" = 0,
          
          "solution" = res_df[, sol_idx],
          
          "allocation" =  tradeoff_sample_allocations[, isurvey],
          
          "true_density" = true_mean,
          
          "true_index_district" = true_index_district,
          
          "post_strata" = district_vals
        )
      )
    temp_res[isurvey, , , isim] <- sim_survey$mean_denisty
    
    if (isim%%50 == 0) print(isim)
  }
}

par(mfrow = c(5, 3), mar = c(3,3,1,1))
for (ispp in 1:15) {
  plot_this <- t(sweep(x = apply(X = temp_res[, , spp_idx_opt[ispp], ],
                                 MARGIN = 1:2,
                                 FUN = sd),
                       MARGIN = 2,
                       STATS = true_mean[spp_idx_opt[ispp], ],
                       FUN = '/'))
  boxplot(plot_this,
          main = common_names_opt[ispp],
          ylim = c(0, max(plot_this)),
          las = 1)
}


##################################################
####   Plot 
##################################################

{
  png(filename = paste0(output_dir, "bethel_only_comparisons.png"),
      width = 170,
      height = 190,
      units = "mm",
      res = 500)
  
  ## Plot layout
  layout(mat = matrix(data = 1:8, ncol = 2, byrow = T))
  par(mar = c(3, 0, 0, 1), 
      oma = c(2, 7.5, 1, 0))
  
  for (i in 1:4) {
    
    ## Base plot
    plot(1,
         type = "n", 
         xlim = c(0, 0.425),
         ylim = c(1, 15),
         axes = F,
         ann = F)
    box()
    
    if(i == 4) mtext(side = 1,
                     text = "Expected CV",
                     line = 2.5)
    
    axis(side = 1, 
         at = seq(from = 0, to = 0.4, by = 0.1),
         labels = seq(from = 0, to = 0.4, by = 0.1))
    
    text(x = 0.35,
         y = 3, 
         labels = paste0(LETTERS[i], ") Lower CV\nThreshold: ", 
                         c(0, threshold)[i]),
         cex = 1)
    
    axis(side = 2, 
         labels = common_names_opt[spp_order], 
         at = 1:ns_opt, 
         las = 1,
         cex.axis = 0.75)
    
    ## Plot orignal solution's expected CV
    points(y = 1:ns_opt,
           x = settings[sol_idx, paste0("CV", spp_order)],
           col = "darkgrey", 
           pch = 16)
    
    ## Plot solution with different lower CV threshold
    if (i > 1)   {
      matpoints(y = 1:ns_opt,
                x = tradeoff_df[spp_order, i-1],
                pch = 16,
                col = "black")
      
      abline(v = threshold[i - 1],
             lty = "dashed",
             col = "black")
    }
    
    ## Calculate change in sample allocations across solutions 
    abs_diff <- apply(X = tradeoff_sample_allocations[, c(1, i)],
                      MARGIN = 1,
                      FUN = diff)
    
    ## Stratum colors based on change in sample size from original solution
    
    if (i == 1) strata_colors <- "white"
    if (i != 1) {
      pos_or_neg <- abs_diff > 0
      strata_colors <- ifelse(test = pos_or_neg,
                              yes = "blue", 
                              no = "red")
      
      for (icolor in 1:length(strata_colors)) {
        color_grad <- colorRampPalette(colors = c("grey90", 
                                                  strata_colors[icolor]))(100)
        
        strata_colors[icolor] <- color_grad[abs(abs_diff)[icolor] + 1]
      } 
    }
    
    
    ## Plot change in 
    image(goa_ras,
          asp = 1,
          axes = F,
          col = strata_colors,
          ann = F)
    
    if(i == 1){
      
      ## Add legend
      plotrix::color.legend(
        xl = xrange[1] + xrange_diff * 0.1,
        xr = xrange[1] + xrange_diff * 0.9,
        yb = yrange[1] + yrange_diff * 0.4,
        yt = yrange[1] + yrange_diff * 0.6,
        legend = c(-100, -50, 0, 50, 100) ,
        rect.col = colorRampPalette(colors = c("red", "grey90", "blue"))(1000) ,
        gradient = "x", 
        align = "rb",
        cex = 0.6)
      mtext(side = 3,
            line = -4, 
            cex = 0.9,
            text = "Absolute change in sample size")
    }
  }
  
  dev.off()
}




