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

The bulk of the optimization is done within the SamplingStrata R Package 
(https://github.com/barcaroli/SamplingStrata). It is best to download the 
package and save the path of the directory. There is one function in the 
package, BuildStrataDF() that I modify for this analysis, so the functions 
in the SamplingStrata package, along with my modified BuildStrataDF() 
function need to be in the global environment. 

## Script Overview (Optimal_Allocation_GoA/analysis_scripts/)

optimization_data.R : Synthesizes data inputs and constants common to 
all subsequent scripts. 

Calculate_Population_Variances.R : Calculates population variances of 
simple random, optimized single-species stratified random, and current 
stratified random surveys.

Survey_Optimization.R : Conducts the multi- and single-species survey 
optimization.

knitting_runs.R : knits all the optimization runs into neat result
outputs.

knitting_runs_SS.R : knits all the single-species optimization runs into
neat result outputs.

Simulate_Surveys.R : Simulates current and optimized stratified random 
surveys.

## Species Included

The species set (ns = 15) included in the manuscript are a complex of Gulf of 
Alaska cods, flatfishes, and rockfishes:

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

## Input Data -- Spatial Domain

The spatial domain of the survey optimization is the Gulf of Alaska 
divided into a X km resolution grid resulting in N = 22832 total survey cells.
The script used to create the survey grid is contained in the MS_OM_GoA repo 
(https://github.com/zoyafuso-NOAA/MS_OM_GoA) in the using the script 
Extrapolation_Grid_Covariates.R in the data/ directory. That script produces 
an RData producted called Extrapolation_depths.RData that is contained within 
the data/ directory this repo. Extrapolation_depths.RData contains a variable 
called Extrapolation_depths which is a dataframe of N rows. Useful fields for 
this analysis are stated in the table below:

| Field Name          | Description                                 |
|---------------------|---------------------------------------------|
| Area_km2            | num, Area of grid cell in square kilometers |
| Lon                 | num, Longitude                              |
| Lat                 | num, Latitude                               |
| Depth_EFH           | num, Depth in meters                        |
| E_km                | num, Eastings in kilometers, 5N UTM         |
| N_km                | num, Northings in kilometers, 5N UTM        |
| stratum             | int, Stratum ID in current STRS design      |
| trawlable           | logi, is the cell trawlable?                |
| shallower_than_700m | logi, is the cell < 700 m?                  |
| shallow_trawlable   | logi, is the cells trawlable and < 700 m    |

Two scenarios of the spatial domain are considered in this analysis and are 
shown below: 1) Full domain (black) and 2) domain when untrawlable and deep
(i.e., > 700 m ) survey grids are removed:

![](graphics/domain.png)

## Input Data -- Predicted denisity
Density of each species was predicted across the spatiotemporal domain using a 
vector autoregressive spatiotemporal model using the VAST package
(https://github.com/James-Thorson-NOAA/VAST). Gulf of Alaska bottom-trawl 
catch-per-unit area survey data were used from years 1996, 1999, and the odd
years from 2003-2019. Code in the repository zoyafuso-NOAA/MS_OM_GoA/ 
(https://github.com/zoyafuso-NOAA/MS_OM_GoA) was used to run the VAST models 
and the output was saved in this repo (model_11/fit_density.RData). This .RData
file contains a variable called "D_gct" which is a 3-D array of dimension 
(N, ns, 24). There are 24 total years (1996-2019), but only NTime = 11 years of
data used. 

## Input Data -- Putting it all together

Data for the optimization were synthesized in the optimization_data.R script. 
It's purpose is to take the VAST model density predictions and create an input 
dataset in the form that is used in the SamplingStrata package. The depth and 
E_km fields are used as strata variables. The script creates two subdirectories
model_11/full_domain/ and model_11/trawlable/, and an .RData file called 
optimization_data.RData is saved in each subdirectory. The output of 
optimization_data.RData contains the following variables and constants: 

| Variable Name | Description                                                                                                                        | Class Type and Dimensions                  |
|---------------|------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------|
| ns            | Number of species in optimization                                                                                                  | numeric vector, length 1                   |
| sci_names     | Scientific species names, used in plots                                                                                            | character vector, length ns                |
| nboats        | Total number of sample sizes of interest, (nboats = 3)                                                                             | numeric vector, length 1                   |
| samples       | Range of sample sizes of interest, corresponding to 1 (n = 280), 2 (n = 550), and 3 (n = 820) boats                                | numeric vector, length nboats              |
| NStrata       | Total number of strata scenarios, (NStrata = 6)                                                                                    | numeric vector, length 1                   |
| stratas       | Range of number of strata, (stratas = c(5, 10, 15, 20, 30, 60))                                                                    | numeric vector, length NStrata             |
| N             | Total number of grid cells in the spatial domain, (N = 23339 cells)                                                                | numeric vector, length 1                   |
| NTime         | Total number of years with data, (NTime = 11 years between 1996-2019)                                                              | numeric vector, length 1                   |
| Niters        | Total number of times a survey is simulated, (Niters = 1000)                                                                       | numeric vector, length 1                   |
| frame         | Annual mean densities for each species, longitude, and depth across grid cells                                                     | dataframe, N rows x 19 columns             |
| frame_raw     | Densities for each species across observed years, along with longitude and depth across cells                                      | dataframe, N*NTime rows x 20 columns       |
| true_mean     | True mean densities for each species and year. This is the "truth" that is used in the performance metrics when simulating surveys | dataframe, NTime rows x ns columns         |

The optimization_data.R script uses the indices in the Extrapolation_depths
dataframe to remove untrawlable/deep survey grids for the trawlable scenario.

## Survey Optimization

The SamplingStrata R package (https://github.com/barcaroli/SamplingStrata)
is used for the optimization. 

The optimization is run over a range of number of stratas from 5 to 60 on the
full domain (model_11/full_domain/Spatiotemporal_Optimization/) and the 
trawlable domain (model_11/trawlable/Spatiotemporal_Optimization/). 
Optimizations were conducted for each boat effort level (../boat1, ../boat2,
../boat3). Each run of the optimization is saved in its own directory with the 
code template of StrXRunY where X is the number of strata in the solution and Y
is the run number. Within each run folder contains:

| File Name            | Description                                                         |
|----------------------|---------------------------------------------------------------------|
| output/plotdom1.png  | Genetic algorithm results                                           |
| output/outstrata.txt | Stratum-level means and variances for each species                  |
| solution.png         | Low-quality snapshot of the solution mapped onto the spatial domain |
| result_list.RData    | Result workspace of the optimization                                |

The result_list.RData workspace contains a named list called result_list, which
consists of the elements:

| Variable Name                    | Description                                                                                                         | Class Type and Dimensions                      |
|----------------------------------|---------------------------------------------------------------------------------------------------------------------|------------------------------------------------|
| result_list$solution$indices     | Solution indexed by strata, contained in the X1 column                                                              | dataframe, N rows and 2 columns                |
| result_list$solution$aggr_strata | Stratum-level means and variances for each species                                                                  | dataframe, variable number of rows, 37 columns |
| result_list$solution$frame_new   | Original data, along with the solution in the STRATO column.                                                        | dataframe, N rows and 21 columns               |
| result_list$sum_stats            | Characteristics of the optimized strata, e.g., allocated sampling, population size, strata variable characteristics | dataframe, variable number of rows, 9 columns  |
| result_list$CV_constraints       | Expected CV across species                                                                                          | numeric vector, length ns                      |
| result_list$n                    | Optimized total sample size                                                                                         | numeric, length 1                              |

## Knitting Together Optimization Results

The results from each run are synthesized in the knitting_runs.R script. Four
variables are saved in the optimization_knitted_results.RData workspace:

| Variable Name     | Description                                                                 | Class Type and Dimensions                     |
|-------------------|-----------------------------------------------------------------------------|-----------------------------------------------|
| settings          | Optimized strata and expected CV for each species and number of strata      | dataframe, variable number of rows, 19 c      |
| res_df            | Solutions for each run                                                      | dataframe, N rows, variable number of columns |
| strata_list       | Collection of result_list$solution$aggr_strata from each run                | list of variable length                       |
| strata_stats_list | Collection of stratum-level means and variances across species for each run | list of variable length                       |

## Survey Simulation and Performance Metrics (work in progress)...

## Graphic Workflow

![](graphics/Workflow.png)

