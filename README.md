---
editor_options: 
  markdown: 
    wrap: 72
---

# Shearwaters_plastic

# System Requirements

R version 4.4.1 (2024-06-14 ucrt)

Matrix products: default

attached base packages:

-   stats
-   graphics
-   grDevices
-   utils
-   datasets
-   methods
-   base

loaded via a namespace (and not attached):

-   compiler_4.4.1
-   tools_4.4.1
-   rstudioapi_0.16.0

Additional packages called within the code files

# Code tested on

Platform: x86_64-w64-mingw32/x64

Running under: Windows 10 x64 (build 19045)

With the above R installation

R package versions as tested can be loaded using command renv::restore()

# Setup

1.  Copy/clone the repo to a local directory
2.  Open Shearwaters.Rproj
3.  Working directory should automatically be recognised as the master
    folder

# Usage

All results in the manuscript can be reproduced by running the code
files detailed below

01_data_wrangling.R

-   Reads in the encounter information
-   Cleans data
-   Converts to a capture-recapture history format
-   Stores individual covariates

02_visualise_data.R

-   Exploratory visualisations of the observed data

03_exploratory_GLM_models.R

-   Fit, plot and obtain estimates from Binomial GLM models

04_HMM_models.R

-   Fit capture-recapture models
-   Produce model average results
-   Bootstrap for standard errors
-   Visualise model results
