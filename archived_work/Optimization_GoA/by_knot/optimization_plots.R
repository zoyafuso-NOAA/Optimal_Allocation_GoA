i = 100
ispp = 2
Save$Spp[ispp]

temp_df = subset(res_df, subset = spp_scen == ispp & n == i)
temp_mat = res_mat[res_df$spp_scen == ispp & res_df$n == i,]

par(mfrow = c(1,1), mar = c(3,3,0,0))
plot(tot_var ~ tot_mean, data = temp_df)

max_spp = which.max(temp_df$tot_mean)
max_var = which.max(temp_df$tot_var)

par(mfrow = c(2,1), mar = c(0,0,0,0))
plot( Spatial_List$loc_g, pch = 16, axes = F )
points(Spatial_List$loc_g[temp_mat[max_spp,] == 1,], pch = 16, col = 'darkgreen')

plot( Spatial_List$loc_g, pch = 16, axes = F )
points(Spatial_List$loc_g[temp_mat[max_var,] == 1,], pch = 16, col = 'blue')

#################################################################
knot_locs = Spatial_List$loc_g

all_sol = array(data = 0, dim = c(ns+1,3,n_x))


for(ispp in 1:(ns+1)){
   for(i in 1:n_x ){
      temp_df = subset(res_df, n == i & spp_scen == ispp)
      
      if(nrow(temp_df) >= 3){
        temp_mat = matrix( res_mat[res_df$n == i & res_df$spp_scen == ispp,],
                           ncol = n_x )
        
        max_mean_idx = which.max(temp_df$tot_mean)
        max_var_idx = which.max(temp_df$tot_var)
        
        range_mean = diff(range(temp_df$tot_mean))
        range_var = diff(range(temp_df$tot_var))
        
        dist_mean = (temp_df$tot_mean - max(temp_df$tot_mean)) / range_mean
        dist_var =  (temp_df$tot_var - max(temp_df$tot_var)) / range_var
        
        comp_idx = which.min(sqrt(dist_mean^2 + dist_var^2))
        
        sol_idx = c(max_mean_idx, comp_idx, max_var_idx)
        
        all_sol[ispp,,] = all_sol[ispp,,] + temp_mat[sol_idx,]
      }
   }
}


par(mfrow = c(3,1), mar = c(0,0,0,0), family = 'serif')
for(ispp in 1:(ns+1)){
  for(isol in 1:3){
    
    col_mag = 100*all_sol[ispp,isol,]/max(all_sol[ispp,isol,])
    col_breaks = cut(col_mag, breaks = 4)
    colors = c('gold', 'orange', 'red', 'darkred')[col_breaks]
    
    plot(knot_locs, pch = 16, axes = F, ann = F, cex = 2, col = colors )
    
    if(isol == 3) legend('bottomright', 
                         legend = c(as.character(Save$Spp), 'all')[ispp],
                         cex = 2, bty = 'n')
  }
}



