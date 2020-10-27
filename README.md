# Multispecies Stratified Survey Optimization for Gulf of Alaska Groundfishes
 
This repository is provides the code used for an In Prep NOAA Technical 
Memorandum manuscript by Zack S. Oyafuso, Lewis A.K. Barnett and Stan Kotwicki 
entitled "..." 

## Requirements

A handful of R packages are required. Some conventional ones:

```
library(sp)
library(raster)
library(RColorBrewer)
```

The bulk of the optimization is done within the SamplingStrata R Package (https://github.com/barcaroli/SamplingStrata). It is best to download the package and save the path of the directory. There is one function in the package, BuildStrataDF() that I modify for this analysis, so the functions in the SamplingStrata package, along with my modified BuildStrataDF() function need to be in the global environment. 

## Script Overview

Survey_Optimization.R : Conducts the survey optimization

knitting_runs.R : knits all the optimization runs into neat result outputs

Simulate_Opt_Survey.R: Simulates optimized surveys

## Species Included

The species set included in the manuscript are a complex of Gulf of Alaska
cods, flatfishes, and rockfishes:

| Scientific Name                     | Common Name                           |
|-------------------------------------|---------------------------------------|
| *Atheresthes stomias*               | arrowtooth flounder                   |
| *Gadus chalcogrammus*               | Alaska or walleye pollock             |
| *Gadus macrocephalus*               | Pacific cod                           |
| *Glyptocephalus zachirus*           | rex sole                              |
| *Hippoglossoides elassodon*         | flathead sole                         |
| *Hippoglossus stenolepis*           | Pacific halibut                       |
| *Lepidopsetta bilineata*            | southern rock sole                    |
| *Lepidopsetta polyxystra*           | northern rock sole                    |
| *Microstomus pacificus*             | Pacific Dover sole                    |
| *Sebastes alutus*                   | Pacific ocean perch                   |
| *Sebastes melanostictus/aleutianus* | blackspotted and rougheye rockfishes* |
| *Sebastes brevispinis*              | yellowfin sole                        |
| *Sebastes polyspinis*               | northern rockfish                     |
| *Sebastes variabilis*               | dusky rockfish                        |
| *Sebastolobus alascanus*            | shortspine thornyhead                 |

*Due to identification issues between two rockfishes these two species were 
combined into a species group we will refer as "Sebastes B_R" (blackspotted 
rockfish and rougheye rockfish, respectively) hereafter. 

## Input Data--Spatial Domain

The spatial domain of the survey optimization is the Gulf of Alaska 
divided into a X km resolution grid. The script used to create the survey grid
is contained in the MS_OM_GoA repo (https://github.com/zoyafuso-NOAA/MS_OM_GoA) 
in the using the script Extrapolation_Grid_Covariates.R in the data/ directory. 
That script produces an RData producted called Extrapolation_depths.RData that 
is contained within the data/ directory this repo.

In some scenarios, only trawlable areas were used as shown below:

![](graphics/domain.png)


Density of each species was predicted across the spatiotemporal domain using a 
vector autoregressive spatiotemporal model using the VAST package
(https://github.com/James-Thorson-NOAA/VAST). Gulf of Alaska bottom-trawl 
catch-per-unit area survey data were used from years 1996, 1999, and the odd
years from 2003-2019. Code in the repository zoyafuso-NOAA/MS_OM_GoA/ was used
to run the VAST models and the output was saved in this repo 
(model_11/fit_density.RData). 

Data for the optimization were synthesized in the optimization_data.R script. 
It's purpose is to take the VAST model density predictions and create an input 
dataset in the form that is used in the SamplingStrata R package. 

Extrapolation_depths.RData contains a variable called Extrapolation_depths 
which is a dataframe that contains the locations, Gulf of Alaska stratum ID,
area, and depths of each grid in the spatial domain. The depth and E_km fields
are used as strata variables. The output of the script is saved as 
optimization_data.RData and contains the following variables and constants. 
0