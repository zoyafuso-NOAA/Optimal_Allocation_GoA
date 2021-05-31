###############################################################################
## Project:      Load and compare historical surveys
## Author:       Megsie Siple (margaret.siple@noaa.gov)
##               contributions from Lewis Barnett (lewis.barnett@noaa.gov)
##               and Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:  Load previous survey data from RACE_BASE (need to check with
##               Lewis). Compare basic info about station number, boat number,
##               station locations, and number of hauls
###############################################################################
rm(list = ls())

##################################################
####   Set up directories
##################################################
which_machine <- c('Zack_MAC' = 1, 'Zack_PC' = 2, 'Zack_GI_PC' = 3)[2]

github_dir <- paste0(c("/Users/zackoyafuso/Documents/", 
                       "C:/Users/Zack Oyafuso/Documents/",
                       "C:/Users/zack.oyafuso/Work/")[which_machine], 
                     "GitHub/Optimal_Allocation_GoA/")

##################################################
####   Import Libaries
##################################################
library(tidyverse)

##################################################
####   Import processed haul data
##################################################
## Load data from flat files exported from AFSC database
haul <- read.csv(paste0(github_dir, "data/haul.csv"))

## Calculate haul midpoints
haul <- cbind(
  haul,
  geosphere::midPoint(
    cbind(haul$START_LONGITUDE, haul$START_LATITUDE),
    cbind(haul$END_LONGITUDE, haul$END_LATITUDE)
  )
) 

## Temporal components of hauls
haul$DATE <- as.Date(haul$START_TIME, format = "%d-%b-%y")
haul$MONTH <- lubridate::month(haul$DATE)
haul$DAY <- lubridate::day(haul$DATE)
haul$YEAR <- lubridate::year(haul$DATE)

## Filter by region, year, and whether the haul was satisfactory for 
## inclusion in abundance index comp
haul <- haul %>%
  filter(
    REGION == "GOA",
    YEAR > 1989,
    YEAR < 2020,
    ABUNDANCE_HAUL == "Y"
  )

head(haul)

## Does the HAUL column contain unique identifiers for each haul within a year?
haul %>%
  group_by(YEAR, HAUL) %>%
  count()

haul %>%
  group_by(YEAR, HAUL, STATIONID) %>%
  count()

## The answer is no! A lot of haul numbers appear in the data more than once. 
## They show up a max of 3 times. I think this is the number of legs on a cruise.
## There are, for example, two "haul 3"s in the 1990 dataset
##
## Station ID corresponds to grid cell that the station is in (row X column)
## STATIONID x STRATUM = unique station ... some are sampled multiple times so
## look for hauls with performance >=0 to get the successful haul if there are
## multiple tries.

# EXAMPLE: 1990 survey
locations <- haul %>%
  dplyr::filter(YEAR == 1990) %>%
  dplyr::distinct(STATIONID, lat, lon, END_LATITUDE, END_LONGITUDE,
                  HAUL, DATE, HAULJOIN) %>%
  # HAULJOIN to order sites by when they were visited (multiple sites/day)
  dplyr::arrange(DATE, HAULJOIN) %>%
  tibble::rowid_to_column("order")

loc_sf <- sf::st_as_sf( x = locations,
                        coords = c("lon", "lat"),
                        crs = 4326, agr = "constant")

locations %>%
  ggplot() +
  geom_sf(data = loc_sf) +
  geom_point(data = locations, aes(x = lon, y = lat, colour = order)) +
  scale_colour_viridis_c("Sampling order", direction = 1) +
  xlab("Longitude") +
  ylab("Latitude")

##################################################
####   Plot the order of the survey
##################################################
get_map <- function(year, hauldata = haul) {
  
  locs <- hauldata %>%
    dplyr::filter(YEAR == year) %>%
    dplyr::distinct(STATIONID, lat, lon, END_LATITUDE, END_LONGITUDE,
             HAUL, DATE, HAULJOIN, VESSEL) %>%
    dplyr::arrange(DATE, HAULJOIN) %>%
    tibble::rowid_to_column("order")
  
  locs_sf <- sf::st_as_sf(x = locs,
                           coords = c("lon", "lat"),
                           crs = 4326, agr = "constant")
  
  mapplot <- locs %>%
    ggplot() +
    geom_sf(data = locs_sf) +
    geom_point(data = locs, aes(x = lon, y = lat, colour = order)) +
    scale_colour_viridis_c("Sampling order", direction = 1) +
    xlab("Longitude") +
    ylab("Latitude") +
    labs(title = paste("Year = ", year)) +
    facet_wrap(~VESSEL, ncol = 1)
  
  return(mapplot)
}

get_map(year = 1990, hauldata = haul)
get_map(year = 2019, hauldata = haul)

##################################################
####   Save updated haul information
##################################################
save(list = "haul", 
     file = paste0(github_dir, "data/processed_haul_data.RData"))
