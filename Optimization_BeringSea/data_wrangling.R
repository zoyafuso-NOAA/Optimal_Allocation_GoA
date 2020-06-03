#############################################
## Data wrangling for linear programming exercise
#############################################
setwd("C:/Users/Zack Oyafuso/Documents/GitHub/MS_OM_GoA/Optimization_BeringSea")

##############################
## Import Packages
##############################
library(tidyr); library(reshape2)

##############################
## Import EBS CPUE dataset
##############################
EBS = read.csv('C:/Users/Zack Oyafuso/Desktop/AK_BTS/data-raw/cpue_EBSshelf_all_spp.csv')
# EBS = read.csv('G:/Oyafuso/data/data-raw/cpue_EBSshelf_all_spp.csv')


##############################
## Subset species and years > 1986
#############################
spp = c('Arrowtooth_Flounder' = 10110, 'Flathead_Sole' = 10130, 
        'Yellowfin_Sole' = 10210, 'Walleye_Pollock' = 21740,
        'AK_Plaice' = 10285)

df = subset(x = EBS,
            subset = YEAR > 1986 
            &!(STATIONID %in% 
                 c('AZ0504', 'GF1918', 'GF2019', 'GF2120', 'GF2221', 'HG1918', 
                   'HG2019', 'HG2120','HG2221', 'IH1918', 'IH2019', 'IH2120', 
                   'IH2221', 'ON2524', 'ON2625', 'PO2423', 'PO2524', 'PO2625',
                   'PO2726', 'QP2423', 'QP2524', 'QP2625', 'QP2726', 'JI1918',
                   'JI2019', 'JI2120', 'JI2221', 'J-13'))
            & SPECIES_CODE %in% spp,  
            select = c('YEAR', 'STATIONID', 'SPECIES_CODE', 'WGTCPUE'))

write.csv(x = df, file = 'EBS_data_trimmed.csv', row.names = F)
