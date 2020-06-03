#################################
## Optimization Plots
################################

setwd('C:/Users/zack.oyafuso/Work/GitHub/MS_OM_GoA/Optimization_BeringSea')

load('optimization_results.RData')
stations = read.csv('10110_stations_.csv')
stations = stations[stations$STATIONID != 'J-13',]
EBS = read.csv('EBS_data_trimmed.csv')

########################################
## Calculate Mean and SD CPUE for each species at each station
########################################
spp_mean = spread(data = aggregate(WGTCPUE ~ SPECIES_CODE + STATIONID, data = EBS, FUN = mean), key = SPECIES_CODE, value = WGTCPUE)

spp_sd = spread(data = aggregate(WGTCPUE ~ SPECIES_CODE + STATIONID, data = EBS, FUN = sd), key = SPECIES_CODE, value = WGTCPUE)


{tiff(filename = 'prelim_figures/Tradeoff.tiff', width = 6, height = 7,
      units = 'in', compression = 'lzw', res = 200)
  par(mfrow = c(3,2), mar = c(2,2,1,1), oma = c(2,2,0.5,0.5))
  for(ispp in 1:(spp + 1)){
    with(subset(x = res_df, spp_scen == ispp),
         plot(tot_mean ~ tot_var, 
              type = 'n', las = 1,
              xlim = range(tot_var)* c(1, 1.05), 
              ylim = range(tot_mean) * c(0.8, 1),
              ann = F)) 
    legend('topleft', spp_labels[ispp], bty = 'n', cex = 1, text.font = 2)
    
    for(i in c(75, 300, seq(from=50, to=200, by=50))){
      lines(tot_mean ~ tot_var, data = res_df, lwd = 2, 
            subset = (spp_scen == ispp) & (n == i) )
      points(tot_mean ~ tot_var, data = res_df, pch = 16, 
             subset = (spp_scen == ispp) & (n == i))
      with(subset(res_df, (spp_scen == ispp) & (n == i)),
           text(max(tot_var), min(tot_mean), paste('n =', i), pos = 1)
      )
    }
  }
  
  mtext(side = 1, 'Variance Criterion', outer = T, line = 0.5, font = 2)
  mtext(side = 2, 'Total Mean Criterion', outer = T, line = 0.5, font = 2)
  
  dev.off()
}

###########################
## Show specific solutions for a given effort level
## Solution with highest mean
## Solution with highest variance
## Compromise solution 
###########################
for(ispp in 1:(spp + 1)){
  
  {tiff(filename = paste0('prelim_figures/sol_placements_',
                          spp_labels[ispp], '.tiff'), 
        width = 8, height = 7, units = 'in', compression = 'lzw', res = 200)
  par(mfrow = c(4,4), oma = c(0,0,3,0))
  
  for(i in c(100, 150, 200, 250)){
    temp_df = subset(res_df, n == i & spp_scen == ispp)
    temp_mat = matrix( res_mat[res_df$n == i & res_df$spp_scen == ispp,], 
                       ncol = nrow(stations) )
    
    max_mean_idx = which.max(temp_df$tot_mean)
    max_var_idx = which.max(temp_df$tot_var)
    
    range_mean = diff(range(temp_df$tot_mean))
    range_var = diff(range(temp_df$tot_var))
    
    dist_mean = (temp_df$tot_mean - max(temp_df$tot_mean)) / range_mean
    dist_var =  (temp_df$tot_var - max(temp_df$tot_var)) / range_var
    
    comp_idx = which.min(sqrt(dist_mean^2 + dist_var^2))
    
    sol_idx = c(max_mean_idx, comp_idx, max_var_idx)
    
    par(mar = c(5,4,1,1))
    plot(tot_mean ~ tot_var, data = temp_df, 
         type = 'b', las = 1, pch = 16, lwd = 2, cex = 1, col = 'grey',
         xlim = range(temp_df$tot_var), ylim = range(temp_df$tot_mean),
         xlab = 'Total Variance', ylab = "Total Mean")
    legend('bottomleft', paste('n =', i), bty = 'n', cex = 1, text.font = 2)
    points(tot_mean ~ tot_var, data = temp_df[sol_idx,], pch = 16, cex = 1.5, 
           col = c('black', 'darkblue', 'red'))
    
    par(mar = rep(0.5, 4))
    for(s in 1:3){
      plot(lat ~ long, data = stations, ann = F, axes = F, asp = 1); box()
      points(lat ~ long, data = stations[temp_mat[sol_idx[s],] == 1,], 
             pch = 16, cex = 2, col = c('black', 'darkblue', 'red')[s])
    }
  }
  
  mtext(side = 3, spp_labels[ispp], outer = T)
  dev.off()}
}

library(tidyr)

# spp_labels = c('Arrowtooth_Flounder', 'Flathead_Sole', 'Yellowfin_Sole',
#                'AK_Plaice', 'Walleye_Pollock' )

#######################################
spp_tradeoff = list(spp_mean = data.frame(), 
                    spp_sd = data.frame())

for(i in c(100, 150, 200, 250, 300)){
  for(ispp in 1:6){
    temp_df = subset(res_df, n == i & spp_scen == ispp)
    temp_mat = matrix( res_mat[res_df$n == i & res_df$spp_scen == ispp,], 
                       ncol = nrow(stations) )
    
    range_mean = diff(range(temp_df$tot_mean))
    range_var = diff(range(temp_df$tot_var))
    
    dist_mean = (temp_df$tot_mean - max(temp_df$tot_mean)) / range_mean
    dist_var =  (temp_df$tot_var - max(temp_df$tot_var)) / range_var
    
    comp_idx = which.min(sqrt(dist_mean^2 + dist_var^2))
    
    rel_mean = colSums(spp_mean[temp_mat[comp_idx,] == 1,-1]) / colSums(spp_mean[,-1])
    names(rel_mean) = spp_labels
    rel_mean = c(ispp, i, rel_mean)
    
    spp_tradeoff$spp_mean = rbind(spp_tradeoff$spp_mean, rel_mean )
    
    rel_sd = colSums(spp_sd[temp_mat[comp_idx,] == 1,-1]) / colSums(spp_sd[,-1])
    names(rel_sd) = spp_labels
    rel_sd = c(ispp, i, rel_sd)
    spp_tradeoff$spp_sd = rbind(spp_tradeoff$spp_sd, rel_sd )
    
  }
  
}

names(spp_tradeoff$spp_mean) = names(spp_tradeoff$spp_sd) =  c('spp_scen', 'n', spp_labels)


{tiff(filename = 'prelim_figures/Tradeoff_heatmap.tiff', width = 6, height = 6,
      units = 'in', compression = 'lzw', res = 200)
  par(mar = c(5,5,2,0), mfrow = c(2,2), oma = c(1,1,0,0))
  for(i in c(100, 150, 200, 300)){
    image(z = t(as.matrix(subset(spp_tradeoff$spp_mean, 
                                 subset = n == i, 
                                 select = spp_labels) )) , 
          axes = F, col = hcl.colors(10, "YlOrRd", rev = TRUE),
          breaks = seq(from = 0, to = 1, by = 0.1), asp = 1)
    text(x = rep(seq(0, 1, length = 5), 6), 
         y = rep(seq(0, 1, length = 6), each = 5),
         
         labels = t(round(as.matrix(subset(spp_tradeoff$spp_mean, 
                                           subset = n == i, 
                                           select = spp_labels) ), 2)) )
    axis(side = 1, at = seq(0, 1, length = 5), las = 2,
         labels = c('ATF', 'Flat', 'YF_Sole', 'AK_Plaice', 'WeP'))
    axis(side = 2, at = seq(0, 1, length = 6), 
         labels = c('ATF', 'Flat', 'YF_Sole', 'AK_Plaice', 'WeP', 'Equal'), las = 1)
    box()
    mtext(side = 3, paste('n =', i))
  }
  
  mtext(side = 2, 'Weighting Scenarios', outer = T)
  mtext(side = 1, 'Species', outer = T)
  dev.off()}
