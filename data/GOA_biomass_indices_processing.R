# read in GOA biomass indices and attach species names
library(dplyr)

# load flat file from database
GOA_indices <- read.csv("data/GOA_biomass_indices.csv")

# join species names
species_codes =  read.csv("data/species.csv", stringsAsFactors = FALSE)
species_codes = select(species_codes, -YEAR_ADDED)
GOA_indices <- inner_join(GOA_indices, species_codes)

# create lumped indices for the following species:
# sculpins, Rougheye and Blackspotted rockfishes, spiny dogfish
sculpins <- GOA_indices[grep("sculpin", GOA_indices$COMMON_NAME),] %>%
  group_by(YEAR) %>%
  summarise(TOTAL_BIOMASS = sum(TOTAL_BIOMASS), BIOMASS_VAR = sum(BIOMASS_VAR),
            MEAN_WGT_CPUE = sum(MEAN_WGT_CPUE), VAR_WGT_CPUE = sum(VAR_WGT_CPUE)) %>%
  mutate(SPECIES_NAME = "Cottoidea", COMMON_NAME = "sculpins") 

Sebastes_BR <- GOA_indices[grep(c("blackspotted|rougheye"), GOA_indices$COMMON_NAME),] %>%
  group_by(YEAR) %>%
  summarise(TOTAL_BIOMASS = sum(TOTAL_BIOMASS), BIOMASS_VAR = sum(BIOMASS_VAR),
            MEAN_WGT_CPUE = sum(MEAN_WGT_CPUE), VAR_WGT_CPUE = sum(VAR_WGT_CPUE)) %>%
  mutate(SPECIES_NAME = "Sebastes BR", COMMON_NAME = "rougheye and blackspotted rockfish") 

GOA_indices <- bind_rows(GOA_indices, sculpins, Sebastes_BR)
  
saveRDS(GOA_indices, "data/GOA_biomass_indices_wnames.rds")