setwd("F:/")

res_df_spline <- res_df <- data.frame()

for(ispp in dir("VAST_Runs_splines/")){
  
  for (itype in 1:10) {
    temp_file <- paste0("VAST_Runs_splines/", ispp, "/CV_", itype, 
                        "/crossval_fit_performance.RData") 
    if(file.exists(temp_file)) {
      load(temp_file)
      res_df_spline <- rbind(res_df_spline, 
                             data.frame(spp_name = ispp,
                                        fold = itype,
                                        prednll = cv_performance$prednll))
      
    }
    
    temp_file <- paste0("VAST_Runs/", ispp, "/CV_", itype, 
                        "/crossval_fit_performance.RData") 
    if(file.exists(temp_file)) {
      load(temp_file)
      res_df <- rbind(res_df, data.frame(spp_name = ispp,
                                         fold = itype,
                                         prednll = cv_performance$prednll))
      
    }
    
  }
  
}

aggregate(prednll ~ spp_name, data = res_df, FUN = mean)
aggregate(prednll ~ spp_name, data = res_df_spline, FUN = mean)
