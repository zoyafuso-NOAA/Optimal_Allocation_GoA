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

## Script Overview (Optimal_Allocation_GoA/analysis_scripts/)

The survey optimization framework is modularized into separate scripts. These
are the scripts used below and the sections following are the order in which
the optimization is conducted. As of now, there are no high-level wrapper
functions that may ease wider general use. 

[optimization_data.R](https://github.com/zoyafuso-NOAA/Optimal_Allocation_GoA/blob/master/analysis_scripts/optimization_data.R):
Synthesizes data inputs and constants common to all subsequent scripts. 

[Survey_Optimization_SS.R](https://github.com/zoyafuso-NOAA/Optimal_Allocation_GoA/blob/master/analysis_scripts/Survey_Optimization_SS.R):
Conducts single-species survey optimization.

[knitting_runs_SS.R](https://github.com/zoyafuso-NOAA/Optimal_Allocation_GoA/blob/master/analysis_scripts/knitting_runs_SS.R):
Knits all the single-species optimization runs into neat result outputs.

[Survey_Optimization.R](https://github.com/zoyafuso-NOAA/Optimal_Allocation_GoA/blob/master/analysis_scripts/Survey_Optimization.R):
Conducts the multispecies survey optimization.

[knitting_runs.R](https://github.com/zoyafuso-NOAA/Optimal_Allocation_GoA/blob/master/analysis_scripts/knitting_runs.R): 
Knits all the multispecies optimization runs into neat result outputs.

[Simulate_Surveys.R](https://github.com/zoyafuso-NOAA/Optimal_Allocation_GoA/blob/master/analysis_scripts/Simulate_Surveys.R):
Simulates current and optimized stratified random surveys.

[survey_distance_travelled/haul_data.R](https://github.com/zoyafuso-NOAA/Optimal_Allocation_GoA/blob/master/analysis_scripts/survey_distance_travelled/haul_data.R): 
Synthesizes haul-level historical data

[survey_distance_travelled/historical_surveys.R](https://github.com/zoyafuso-NOAA/Optimal_Allocation_GoA/blob/master/analysis_scripts/survey_distance_travelled/historical_surveys.R): 
Calculates actual distance travelled and the approximate shortest 
total distance using the TSP package for the historically observed
stations for each survey year. The nearest and second closest stations 
are also calculated.

[survey_distance_travelled/survey_feasibility.R](https://github.com/zoyafuso-NOAA/Optimal_Allocation_GoA/blob/master/analysis_scripts/survey_distance_travelled/survey_feasibility.R): 
Simulates station locations under the current and optimized STRS designs
and calculates the approximate shortest total distance using the TSP package.
