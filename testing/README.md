# Residuals and Diagnostics for Multinomial Regression Models

This repository provides the code and datasets used in the paper authored by Eric A. E. Gerber and Bruce A. Craig.

https://doi.org/10.1002/sam.11645

## Code Files:

* `sam_gerber_craig_code.R`: The main code file containing most of the relevant code to reproduce the majority of the manuscript.
* `mln_functions.R`: This file contains the necessary functions to fit the multinomial logistic-normal model, which serves as a candidate model throughout the paper.
* `ks_sim_study.R`: This file uses parallelized code to conduct the simulation studies of the Kolmogorov-Smirnov tests.

## Dataset Files:

* `bbs_63014_saxapahaw.xlsx`: This file contains one of the real datasets used in the paper.

Two other datasets are available via R packages:
* Forest pollen data: MM R package.
* Baseball data: Lahman R package.

### Dataset `bbs_63014_saxapahaw.xlsx` Description:

This dataset includes information about bird species and weather conditions:

* Column 1: Year (1998-2019)
* Columns 2-83: Bird Species
* Column 84: Maximum Temperature (Fahrenheit)
* Column 85: Minimum Temperature (Fahrenheit)
* Column 86: Total Precipitation (Inches)
* Column 87: Maximum Daily Precipitation (Inches)
* Column 88: Rainy Days in Period
* Column 89: Maximum Sustained Wind (Miles per Hour)

#### Data Sources:

Bird species count data (Columns 1-83) were obtained from the North American Breeding Bird Survey (BBS) for Saxapahaw (63014). Additional information:

* Start Time: 0532
* Valid Survey Window: May 22 - June 30
* Number of Years Run: 46
* First Year Run: 1966
* Most Recent Year Run: 2019
* Start Latitude: 35.9898149
* Start Longitude: -79.4688523
* Bird Conservation Region: Piedmont (29)

Weather data (Columns 84-89) were obtained from [Visual Crossing](https://www.visualcrossing.com/weather-history/35.9898149%2C-79.4688523/us/2002-05-22/2002-06-30).
