# library(devtools)
# devtools::install_local('C:/Users/zack.oyafuso/Downloads/SamplingStrata-master/', force = T)
library(SamplingStrata)
library(sp)
data("meuse")# locations (155 observed points)
data("meuse.grid")# grid of points (3103)
library(gstat)
library(automap)

meuse.grid$id <- c(1:nrow(meuse.grid))
coordinates(meuse)<-c("x","y")
coordinates(meuse.grid)<-c("x","y")

#Spatiotemporal Model 
lm_lead <- lm(log(lead) ~ dist,data=meuse)
lm_zinc <- lm(log(zinc) ~ dist,data=meuse)

#Get values across sampling grid via kriging
kriging_lead = autoKrige(log(lead) ~ dist, meuse, meuse.grid)
kriging_zinc = autoKrige(log(zinc) ~ dist, meuse, meuse.grid)

#Create df with metal predictions and varainces across the domain
df <- NULL
df$id <- meuse.grid$id
df$lead.pred <- kriging_lead$krige_output@data$var1.pred
df$lead.var <- kriging_lead$krige_output@data$var1.var
df$zinc.pred <- kriging_zinc$krige_output@data$var1.pred
df$zinc.var <- kriging_zinc$krige_output@data$var1.var
df$lon <- meuse.grid$x
df$lat <- meuse.grid$y
df$dom1 <- 1
df <- as.data.frame(df)
head(df)

# Produce optimal stratification of the 3,103 points under a precision 
# constraint of 1% on the target estimates
frame <- buildFrameSpatial(df=df,
                           id="id",
                           X=c("lead.pred","zinc.pred"),
                           Y=c("lead.pred","zinc.pred"),
                           variance=c("lead.var","zinc.var"),
                           lon="lon",
                           lat="lat",
                           domainvalue = "dom1")
cv <- as.data.frame(list(DOM=rep("DOM1",1),
                         CV1=rep(0.01,1),
                         CV2=rep(0.01,1),
                         domainvalue=c(1:1) ))

set.seed(1234)
solution <- optimStrata (
  method = "spatial",
  errors=cv, 
  framesamp=frame,
  iter = 15,
  pops = 20,
  mut_chance = 50,
  nStrata = 5,
  fitting = c(summary(lm_lead)$r.square,summary(lm_zinc)$r.square),
  range = c(kriging_lead$var_model$range[2],kriging_zinc$var_model$range[2]),
  kappa=1,
  writeFiles = FALSE,
  showPlot = TRUE,
  parallel = TRUE)

expected_CV(solution$aggr_strata)

#Save Results
framenew <- solution$framenew
outstrata <- solution$aggr_strata

strataStructure <- summaryStrata(solution$framenew,
                                 solution$aggr_strata,
                                 progress=FALSE)
strataStructure

#Plot resuts
frameres <- SpatialPointsDataFrame(data=framenew, coords=cbind(framenew$LON,framenew$LAT) )
frameres2 <- SpatialPixelsDataFrame(points=frameres[c("LON","LAT")], data=framenew)
frameres2$LABEL <- as.factor(frameres2$LABEL)
spplot(frameres2,c("LABEL"), col.regions=bpy.colors(5))
