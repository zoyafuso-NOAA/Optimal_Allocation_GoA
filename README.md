# Multispecies Stratified Survey Optimization for Gulf of Alaska Groundfishes
 
This repository provides the code used for an In Prep NOAA Technical 
Memorandum manuscript by Zack Oyafuso, Lewis Barnett, Margaret Siple,
and Stan Kotwicki temporarily entitled "The expected performance and 
feasibility of a Gulf of Alaska groundfish bottom trawl survey optimized 
for abundance estimation." 

## Package Requirements

A handful of R packages are required. R version 4.0.2 was used for the 
analysis. Some conventional packages for plotting and manipulating data:

VAST version 3.6.1,
and FishStatsUtils version 2.8.0 (2020-09-22). 

```
library(tidyverse)
library(sp)
library(raster)
library(RColorBrewer)
```

The creation of the operating models was done using the [VAST R package](https://github.com/James-Thorson-NOAA/VAST)
version 3.6.1 and FishStatsUtils version 2.8.0.

```
library(VAST)
```

The bulk of the optimization is done within the SamplingStrata R Package 
(https://github.com/barcaroli/SamplingStrata). There is one function in 
the package, BuildStrataDF() that I modify for this analysis, so it is 
best to use a forked version of the package that I modified:

```
library(devtools)
devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata")
library(SamplingStrata)
```

Lastly, the TSP package was imported to calculate the shortest path
that visits each station.

```
library(TSP)
```

## Species Included

The species set included in the survey optimization is a complex of fifteen 
Gulf of Alaska cods, flatfishes, and rockfishes:

| Scientific Name                     | Common Name                           |
|-------------------------------------|---------------------------------------|
| *Gadus chalcogrammus*               | Alaska or walleye pollock             | 
| *Gadus macrocephalus*               | Pacific cod                           |   
| *Atheresthes stomias*               | arrowtooth flounder                   | 
| *Hippoglossoides elassodon*         | flathead sole                         |
| *Glyptocephalus zachirus*           | rex sole                              | 
| *Lepidopsetta polyxystra*           | northern rock sole                    |
| *Lepidopsetta bilineata*            | southern rock sole                    | 
| *Microstomus pacificus*             | Dover sole                            | 
| *Hippoglossus stenolepis*           | Pacific halibut                       |
| *Sebastes alutus*                   | Pacific ocean perch                   |
| *Sebastes melanostictus/aleutianus* | blackspotted and rougheye rockfishes* |
| *Sebastes brevispinis*              | silvergrey rockfish                   |
| *Sebastes variabilis*               | dusky rockfish                        |
| *Sebastes polyspinis*               | northern rockfish                     |
| *Sebastolobus alascanus*            | shortspine thornyhead                 |

*Due to identification issues between two rockfishes these two species were 
combined into a species group we will refer as "BE and RS rockfishes"
(blackspotted rockfish and rougheye rockfish, respectively). 

In addtion, eleven species/species groups were included in the survey evaluations.
These taxa were not included in the optimization but included when simulating
surveys:

| Scientific Name                     | Common Name                           |
|-------------------------------------|---------------------------------------|
| *Anoplopoma fimbria*                | sablefish                             |
| *Pleurogrammus monopterygius*       | Atka mackerel                         |
| *Sebastes borealis*                 | shortraker rockfish                   |
| *Sebastes variegatus*               | harlequin rockfish                    |
| *Sebastes ruberrimus*               | yelloweye rockfish                    |
| Species from Genera *Hemitripterus*, *Hemilepidotus*, and *Myoxocephalus* | sculpins (plain, great, and bigmouth  sculpins and yellow Irish lord) |
| *Beringraja binoculata*             | big skate                             |
| *Raja rhina*                        | longnose skate                        |
| *Albatrossia Pectoralis*            | giant grenadier                       |
| *Enteroctopus dofleini*             | giant octopus                         |
| *Squalus suckleyi*                  | Pacific spiny dogfish                 |

## Survey Optimization Scenarios
|     Survey Exploration                                              |     Scenario Option/Scenario Label        |     A    |     B    |     C    |     D    |     E    |     F    |     G    |     H    |     I    |     J    |     K    |     L    |     M    |
|---------------------------------------------------------------------|-------------------------------------------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|
|                                                                     |                                           |          |          |          |          |          |          |          |          |          |          |          |          |          |
|     Spatial scale of optimization                                   |     gulf-wide                             |     X    |          |          |          |          |          |          |          |          |          |          |          |          |
|                                                                     |     area-level                            |          |     X    |     X    |     X    |     X    |     X    |     X    |     X    |     X    |     X    |     X    |     X    |     X    |
|                                                                     |                                           |          |          |          |          |          |          |          |          |          |          |          |          |          |
|     Stratum variables                                               |     depth and longitude                   |     X    |     X    |          |     X    |          |     X    |     X    |          |          |     X    |          |          |          |
|                                                                     |     depth only                            |          |          |     X    |          |     X    |          |          |     X    |     X    |          |     X    |          |          |
|                                                                     |     existing STRS design                  |          |          |          |          |          |          |          |          |          |          |          |     X    |     X    |
|                                                                     |                                           |          |          |          |          |          |          |          |          |          |          |          |          |          |
|     Sensitivity of data inputs                                      |     MLE                                   |     X    |     X    |     X    |     X    |     X    |          |          |          |          |     X    |     X    |     X    |     X    |
|                                                                     |     measurement error simulated           |          |          |          |          |          |     X    |          |     X    |          |          |          |          |          |
|                                                                     |     fixed and random effects simulated    |          |          |          |          |          |          |     X    |          |     X    |          |          |          |          |
|                                                                     |                                           |          |          |          |          |          |          |          |          |          |          |          |          |          |
|     One deep stratum: cells deeper than 300   m are set to 300 m    |     no                                    |     X    |     X    |     X    |          |          |     X    |     X    |     X    |     X    |     X    |     X    |     X    |     X    |
|                                                                     |     yes                                   |          |          |          |     X    |     X    |          |          |          |          |          |          |          |          |
|                                                                     |                                           |          |          |          |          |          |          |          |          |          |          |          |          |          |
|     Maximum depth cutoff                                            |     1000 m                                |     X    |     X    |     X    |     X    |     X    |     X    |     X    |     X    |     X    |          |          |     X    |          |
|                                                                     |     700 m                                 |          |          |          |          |          |          |          |          |          |     X    |     X    |          |     X    |
## Script Overview 

Tasks are done modularly from data wrangling to VAST model fits to survey optimization to survey simulations to figure/table production. See the [wiki page](https://github.com/zoyafuso-NOAA/Optimal_Allocation_GoA/wiki) for a breakdown of the script workflow. 
