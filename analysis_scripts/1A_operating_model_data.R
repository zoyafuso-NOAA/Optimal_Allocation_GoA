###############################################################################
## Project:         GOA Groundfish CPUE Data Synthesis 
## Author:          Zack Oyafuso (zack.oyafuso@noaa.gov)
## Contributors:    Lewis Barnett (lewis.barnett@noaa.gov)
## Description:     Create CPUE dataset used for VAST for species of interest
###############################################################################
rm(list = ls())

##################################################
#### Require Packages
##################################################
library(dplyr)
library(raster)

##################################################
####   CRSs used
##################################################
lonlat_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
utm_crs <- "+proj=utm +zone=5N +units=km"

##################################################
#### Import CPUE survey data
#### Import Haul-level data
#### Import Species codes
##################################################
data <- read.csv(gzfile("data/raw_data/cpue_GOA_selected_spp.csv.gz"))

haul <- read.csv("data/raw_data/haul.csv", stringsAsFactors = FALSE)

species_codes <- read.csv("data/raw_data/species.csv", stringsAsFactors = FALSE)

##################################################
####   Merge together bathymetry rasters
##################################################
split_bathy <- list()
n_split_rasters <- length(dir("data/raw_data/split_goa_bathy_ras/")) / 2
for (i in 1:n_split_rasters) {
  split_bathy[[i]] <- raster::raster(paste0("data/raw_data/",
                                            "split_goa_bathy_ras",
                                            "/goa_bathy_processed_", 
                                            i, ".grd"))
}

bathy <- SpaDES.tools::mergeRaster(split_bathy)
rm(split_bathy)

##################################################
#### Join haul data to get coordinates, depth, bottom and surface temperature
##################################################
haul <- cbind(haul, 
              # get haul midpoints
              geosphere::midPoint(cbind(haul$START_LONGITUDE, 
                                        haul$START_LATITUDE), 
                                  cbind(haul$END_LONGITUDE, 
                                        haul$END_LATITUDE))) 
haul$DATE <- as.Date(haul$START_TIME, format = "%d-%b-%y")
haul$MONTH <- lubridate::month(haul$DATE)
haul$DAY <- lubridate::day(haul$DATE)
haul <- haul %>% dplyr::select(HAULJOIN, SURFACE_TEMPERATURE, 
                               GEAR_DEPTH, GEAR_TEMPERATURE, 
                               LATITUDE = lat, LONGITUDE = lon, 
                               DATE, DAY, MONTH)
data <- inner_join(data, haul)

##################################################
####   Join species names
##################################################
species_codes = dplyr::select(species_codes, -YEAR_ADDED)
data <- inner_join(data, species_codes)

##################################################
####  Select and rename columns, dropping rows with mising depths
##################################################
data <- data %>% dplyr::select(YEAR, SURVEY, 
                               BOTTOM_DEPTH = GEAR_DEPTH,
                               SURFACE_TEMPERATURE, GEAR_TEMPERATURE, 
                               # CPUE = WGTCPUE,
                               EFFORT, WEIGHT,
                               LATITUDE, LONGITUDE, DATE, DAY, MONTH, 
                               SPECIES_NAME, COMMON_NAME) %>%
  tidyr::drop_na(BOTTOM_DEPTH,LATITUDE,LONGITUDE) 

##################################################
####  Filter to GOA survey, remove tows with 0 bottom depth, and drop 2001,
####  the year when the survey was incomplete and years before 1996 when a 
####  different net/soak time was used
##################################################
data <- data %>% filter(SURVEY == "GOA", 
                        BOTTOM_DEPTH > 0, 
                        YEAR != 2001 & YEAR >= 1996)

##################################################
####  Sum catches of northern and southern rock sole with rock sole unid.
#### (not distinguished until 1996), rename species complex
##################################################
rock_soles <- data %>% 
  dplyr::filter(COMMON_NAME %in% c("rock sole unid.", 
                                   "southern rock sole", 
                                   "northern rock sole")) %>%
  group_by_at(vars(-WEIGHT, -COMMON_NAME, -SPECIES_NAME)) %>%
  summarise(WEIGHT = sum(WEIGHT)) %>%
  ungroup() %>%
  mutate(SPECIES_NAME = "Lepidopsetta spp.", COMMON_NAME = "rock soles")

data <- as.data.frame(rbind(data, rock_soles))

##################################################
####  Sum catches of blackspooted and rougheye rocks with rougheye and 
####  blackspotted rockfish unid.,  rename species complex
##################################################
B_R_rockfishes <- data %>% dplyr::filter(
  COMMON_NAME %in% c("blackspotted rockfish", 
                     "rougheye rockfish", 
                     "rougheye and blackspotted rockfish unid.")) %>%
  group_by_at(vars(-WEIGHT, -COMMON_NAME, -SPECIES_NAME)) %>%
  summarise(WEIGHT = sum(WEIGHT)) %>%
  ungroup() %>%
  mutate(SPECIES_NAME = "Sebastes B_R", COMMON_NAME = "BS and RE rockfishes")
data <- as.data.frame(rbind(data, B_R_rockfishes))

##################################################
####  Sum catches of sculpins
##################################################
sculpins <- data %>% dplyr::filter(
  COMMON_NAME %in% c("bigeye sculpin", "great sculpin", 
                     "plain sculpin", "yellow Irish lord")) %>%
  group_by_at(vars(-WEIGHT, -COMMON_NAME, -SPECIES_NAME)) %>%
  summarise(WEIGHT = sum(WEIGHT)) %>%
  ungroup() %>%
  mutate(SPECIES_NAME = "sculpins", COMMON_NAME = "sculpins")
data <- as.data.frame(rbind(data, sculpins))

##################################################
####   Change name of spiny dogfish to Pacific spiny dogifsh
##################################################
data$COMMON_NAME[data$COMMON_NAME == "spiny dogfish"] <- "Pacific spiny dogfish"


##################################################
##################################################
data <- subset(data,
               COMMON_NAME %in% c(
                 ## Species included in the survey optimization
                 "arrowtooth flounder", ## Atherestes stomias
                 "Pacific cod", ## Gadus macrocephalus
                 "walleye pollock", ## Gadus chalcogrammus
                 "rex sole", ## Glyptocephalus zachirus
                 "flathead sole", ## Hippoglossoides elassodon
                 "Pacific halibut", ## Hippoglossus stenolepis
                 "southern rock sole", ## Lepidopsetta bilineata
                 "northern rock sole", ## Lepidopsetta polyxystra
                 "Pacific ocean perch", ## Sebastes alutus
                 "silvergray rockfish", ## Sebastes brevispinis
                 "northern rockfish", ## Sebastes polyspinis
                 "dusky rockfish", ## Sebastes variabilis
                 "BS and RE rockfishes", ## Sebastes aleutianus and S. melanostictus
                 "Dover sole", ## Microstomus pacificus
                 "shortspine thornyhead", ## Sebastolobus alascanus
                 
                 ## Species not included in the survey optimization, but 
                 ## included when simulating surveys
                 "sablefish",
                 "Atka mackerel",
                 "shortraker rockfish",
                 "Pacific spiny dogfish",
                 "yelloweye rockfish",
                 "giant octopus",
                 "longnose skate",
                 "big skate",
                 "harlequin rockfish",
                 "giant grenadier",
                 "sculpins"
               ))

##################################################
####   Assign station depths from EFH layer
##################################################
cpue_shape = sp::SpatialPointsDataFrame(
  coords = data[, c("LONGITUDE", "LATITUDE")],
  data = data,
  proj4string = CRS(lonlat_crs))

cpue_shape_aea <- sp::spTransform(x = cpue_shape,
                                  CRSobj = crs(bathy))
cpue_shape_aea@data$depth =  raster::extract(x = bathy,
                                             y = cpue_shape_aea,
                                             method = "simple")

##################################################
####   Plot bathymetry and station locations along with stations without 
####   assigned depths
####
##################################################
plot(bathy)
plot(cpue_shape_aea, 
     add = T,
     pch = ".")
plot(cpue_shape_aea[is.na(cpue_shape_aea@data$depth),], 
     add = T,
     pch = 16,
     col = 'red')

mismatched_idx = which(is.na(cpue_shape_aea@data$depth))
summary(data[mismatched_idx, "BOTTOM_DEPTH"])

##################################################
####   Plot correlation between EFH depths and reported depth from BTS
##################################################
plot(depth ~ BOTTOM_DEPTH,
     data = cpue_shape_aea@data,
     subset = COMMON_NAME == "arrowtooth flounder",
     xlab = "Depth recorded by the BTS",
     ylab = "Depth extracted from EFH layer")
with(subset(cpue_shape_aea@data,
            subset = COMMON_NAME == "arrowtooth flounder"), 
     cor(depth, BOTTOM_DEPTH, use = "complete.obs"))
abline(a = 0, b = 1)

##################################################
####   Attach depths to dataset, scaled
##################################################
data$DEPTH_EFH = cpue_shape_aea@data$depth
data$DEPTH_EFH[is.na(data$DEPTH_EFH)] = data$BOTTOM_DEPTH[is.na(data$DEPTH_EFH)]
data$LOG_DEPTH_EFH = log(data$DEPTH_EFH)
data$LOG_DEPTH_EFH_CEN = scale(data$LOG_DEPTH_EFH)
data$LOG_DEPTH_EFH_CEN_SQ = data$LOG_DEPTH_EFH_CEN^2

##################################################
#### Save
##################################################
data <- data[order(data$YEAR, data$SPECIES_NAME),]

if(!dir.exists("data/processed/")) dir.create("data/processed/")
write.csv(x = data, 
          file = "data/processed/goa_vast_data_input.csv", 
          row.names = F)
