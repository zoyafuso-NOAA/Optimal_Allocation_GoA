## Get AK groundfish bottom trawl survey data for 3 primary surveys
# Result is cleaned cpue (kg/km^2) by haul, with zeros included

# Note: EBS includes all species, whereas GOA and AI are only a subset
# here, we will filter the data to only include species represented in all regions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)

# Extract directly from AFSC database query if have permissions
#install.packages("G:/Conner/R/Package Binaries/sumfish_3.1.25.zip", repos = NULL) #install specific version of sumfish if necessary
#library(sumfish)
#setUser()
#EBS = "SELECT * FROM ebsshelf.ebsshelf_cpue" %>% getSQL()
#GOA = "SELECT * FROM goa.cpue" %>% getSQL()
#AI = "SELECT * FROM ai.cpue" %>% getSQL()

# or from flat files exported from AFSC database
setwd("C:/Users/lewis.barnett/Work/AFSC/Data/AK_BTS/data-raw/")
GOA = read.csv("cpue_GOA_selected_spp.csv", stringsAsFactors = FALSE) # CPUE is (num or kg / km^2)

# filter to pollock
data <- filter(GOA, SPECIES_CODE == 21740)

# join haul data to get coordinates
haul <- read.csv("haul.csv", stringsAsFactors = FALSE)
haul <- cbind(haul, 
             geosphere::midPoint(cbind(haul$START_LONGITUDE, haul$START_LATITUDE), 
                                     cbind(haul$END_LONGITUDE, haul$END_LATITUDE))) # get haul midpoints
haul <- haul %>% select(HAULJOIN, LATITUDE = lat, LONGITUDE = lon)
data <- inner_join(data, haul)

# get biomass from cpue and effort
data <- data %>% mutate(BIOMASS = WGTCPUE * EFFORT) %>%
  select(HAULJOIN, YEAR, BIOMASS, EFFORT, LATITUDE, LONGITUDE) %>%
  tidyr::drop_na()

# pull in length data
# TO DO LATER
# length = read.csv("length.csv", stringsAsFactors = FALSE)

# subset data to western GOA by dropping all records east of extent of 2001 survey
#boundary = as.numeric(data %>% 
#  filter(YEAR == 2001) %>% 
#  summarise(max(LONGITUDE)))
# or subset to western and central GOA as defined in the assessment
boundary = -140

data = data %>% filter(LONGITUDE <= boundary)


write.csv(data, "C:/Users/lewis.barnett/Work/AFSC/Data/AK_BTS/data/AK_BTS_GOA_pollock_140.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# or use sumfish to get everything

#library(sumfish)
#EBS <- sumHaul(getRacebase(year=c(2018,1982), survey='EBS_SHELF')) #file size too large, query in smaller chunks of years and combine
#GOA2 <- sumHaul(getRacebase(year=c(2017,1984), survey='GOA'))
#AI <- sumHaul(getRacebase(year=c(2018,1980), survey='AI'))
#data <- rbind(EBS,GOA,AI)

# steps below are likely needed but havent been tested
#data$wCPUE = data$wCPUE * 100 # convert CPUE to kg/km2 from ha
#species_codes =  read.csv("C:/Users/lewis.barnett/Work/AFSC/Data/AK_BTS/species.csv", stringsAsFactors = FALSE) # or getRacebase()$species
#left_join(data, species_codes)
#data %>% dplyr::select(YEAR,REGION,START_LATITUDE,START_LONGITUDE,BOTTOM_DEPTH,SURFACE_TEMPERATURE,GEAR_TEMPERATURE,SPECIES_NAME,COMMON_NAME,wCPUE)

