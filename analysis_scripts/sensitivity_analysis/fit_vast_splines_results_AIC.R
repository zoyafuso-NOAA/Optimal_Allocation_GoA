setwd("D:/VAST_Runs/")
library(tidyr)

AIC_df <- data.frame()

for(ispp in dir()) {
  load(paste0(ispp, "/fit.RData"))
  
  depth_in_model <- length(grep(x = ispp, pattern = "_depth", value = FALSE))
  
  spp_name <- ispp
  if(length(grep(pattern = "depth", x = spp_name)) == 1) 
    spp_name <- gsub(x = spp_name, pattern = "_depth", replacement = "")
  
  AIC_df <- 
    rbind(AIC_df,
          data.frame(species = spp_name,
                     depth = ifelse(test = depth_in_model == 1,
                                    yes = "depth", no = "none"),
                     AIC = round(fit$parameter_estimates$AIC),
                     BIC = round(2 * fit$parameter_estimates$objective + length(fit$parameter_estimates$par)*log(7900) )))
  
}

setwd("D:/VAST_Runs_splines/")
for(ispp in dir()[-c(1:3)]) {
  load(paste0(ispp, "/fit.RData"))
  
  depth_in_model <- length(grep(x = ispp, pattern = "_depth", value = FALSE))
  spp_name <- gsub(x = ispp, pattern = "_depth", replacement = "")
  
  AIC_df <- 
    rbind(AIC_df,
          data.frame(species = spp_name,
                     depth = "depth_splines",
                     AIC = round(fit$parameter_estimates$AIC),
                     BIC = round(2 * fit$parameter_estimates$objective + length(fit$parameter_estimates$par)*log(7900) )))
  
  
}

tidyr::spread(data = subset(AIC_df, select = -AIC), key = depth, value = BIC)
