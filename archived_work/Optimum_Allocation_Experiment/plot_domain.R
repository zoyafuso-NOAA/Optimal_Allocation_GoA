############################
## Plot Average Spatial Effect, Spatiotemporal Effects, and Log-Densities
## Of each OM scenario
############################
rm(list = ls())

######################################
## Import OM
######################################
which_machine = c('Zack_Mac' = 1)[1]
output_wd = paste0(c('/Users/zackoyafuso/Documents/')[which_machine],
                   'GitHub/MS_OM_GoA/Optimum_Allocation_Experiment/')
load(paste0(output_wd, '/OM.RData'))

######################################
## Import Libaries
######################################
library(sp); library(raster); library(RColorBrewer); library(plotrix)

{
  #Create the plot order of the panels
  plot_order = c()
  for(i in 1:2){ 
    plot_order = c(plot_order, 
                   8*(i-1) + 1:2,
                   8*2 + 4,
                   8*(i-1) + 3:5,
                   8*2 + 5,
                   8*(i-1) + 6:8)
  }
  plot_order = c(plot_order, rep(8*2 + 1, 2), rep(8*2 + 2, 4), rep(8*2 + 3, 4))
  
  par(mar = c(0,0,0,0), oma = c(0,1,3,1))
  layout(mat = matrix(plot_order, nrow = 3, byrow = T),
         widths = c(1,1, 0.1, 1,1,1, 0.1, 1,1,1),
         heights = c(rep(1,2,0.1)) )
  
  for(irow in 1:2){
    plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), axes = F, ann = F)
    #box()
    text(x = 0.45, y = 0.5, 
         paste0('Scale: ', ST_settings$Scale[irow], '\n',
                'Spatial SD: ', ST_settings$SD_omega[irow], '\n',
                'SpatTemp SD: ', ST_settings$SD_epsilon[irow]), xpd = NA)
    
    if(irow == 1) mtext(side = 3, 'Spatial\nSettings', line = 0.25)
    
    test = OM[[irow]]
    domain = SpatialPointsDataFrame(coords = test$loc_xy,
                                    data = data.frame(omega = test$Omega_r) )
    domain_ras = raster(domain, resolution = 1/100)
    domain_ras = rasterize(x = domain, y = domain_ras, field = 'omega')
    image(domain_ras, col = rev(brewer.pal(n = 11, name = 'RdBu')), 
          zlim = c(-2,2), axes = F)
    box()
    
    if(irow == 1) mtext(side = 3, 'Spatial\nEffect', line = 0.25)
    
    for(iyear in c(1,3,9)){
      epsilon = test$Epsilon_rt[,iyear]
      domain = SpatialPointsDataFrame( coords = test$loc_xy,
                                       data = data.frame(omega = epsilon) )
      
      domain_ras = raster(domain, resolution = 1/100)
      domain_ras = rasterize(x = domain, y = domain_ras, field = 'omega')
      
      image(domain_ras, axes = F, col = rev(brewer.pal(n = 11, name = 'RdBu')), 
            zlim = c(-2,2))
      box()
      
      if(irow == 1 & iyear == 3) 
        mtext(side = 3, 
              'Spatiotemporal Effect at\nYears 1, 3, and 9', line = 0.25)
    }
    
    for(iyear in c(1,3,9)){
      logden = test$log_d_rt[,iyear]
      domain = SpatialPointsDataFrame(coords = test$loc_xy,
                                      data = data.frame(omega = logden) )
      domain_ras = raster(domain, resolution = 1/100)
      domain_ras = rasterize(x = domain, y = domain_ras, field = 'omega')
      
      image(domain_ras, col = brewer.pal(n = 9, name = 'Purples') ,
            axes = F, zlim = c(-1,3.5) )
      box()
      if(irow == 1 & iyear == 3) 
        mtext(side = 3, 
              'Log-Density at\nYears 1, 3, and 9', line = 0.25)
    }
  }
  
  plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), axes = F, ann = F)
  
  plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), axes = F, ann = F)
  plotrix::color.legend(xl = 0.05, xr = 0.95, yb = 0.3, yt = 0.6,
                        rect.col = rev(brewer.pal(n = 11, name = 'RdBu')),
                        legend = -2:2)
  
  plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), axes = F, ann = F)
  plotrix::color.legend(xl = 0.05, xr = 0.95, yb = 0.3, yt = 0.6,
                        rect.col = brewer.pal(n = 9, name = 'Purples'),
                        legend = -1:3 )
  plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), axes = F, ann = F)
  plot(1, type = 'n', xlim = c(0,1), ylim = c(0,1), axes = F, ann = F)
  
}
