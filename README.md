# Multispecies Stratified Survey Optimization (Still working on this main page...)

Code repository for the manuscript (In Prep) entitled "Incorporating spatiotemporal variability and prespecified uncertainty in the optimization of a multispecies survey design" by Zack Oyafuso, Lewis Barnett, and Stan Kotwicki.

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

