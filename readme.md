# CRRM

## Introduction

This repository contains code for the manuscript:  
[Ahmed, Z., & Levesque, I. R. (2016). Increased robustness in reference region model analysis of DCE MRI using two-step constrained approaches. Magnetic Resonance in Medicine](http://doi.org/10.1002/mrm.26530).

Most of the code is written in [MATLAB](https://www.mathworks.com/products/matlab/).
There are some parts, specifically those that require non-linear least squares fitting, where [Julia](http://julialang.org/) is used with the [DCEMRI.jl](https://github.com/davidssmith/DCEMRI.jl) package.
Julia was used because it is much faster than Matlab when it comes to non-linear fitting.
There is also some [R](https://cloud.r-project.org/index.html) code, particularly for visualizing the results.
To identify these files by their extensions: .m is for Matlab, .jl is for Julia, and .R & .Rmd are for R.

## Using Julia

Julia is needed to do non-linear fitting. To get running, the steps are:

- Download [Julia](http://julialang.org/)
- Python is also required - if you don't have it, I suggest grabbing [Anaconda](https://www.continuum.io/downloads)
- If you wish to use an IDE with Julia, then check out [Juno](http://junolab.org/) which runs on [Atom](https://atom.io/). The IDE is practical, but is not required.
- Install pre-requisites by first running Julia, and entering the following in the Julia command prompt:
```
Pkg.add("DCEMRI")
Pkg.build("HDF5")
Pkg.build("PyCall")
Pkg.add("Hwloc")
Pkg.update()
```
- Now it should be possible to run Julia code by opening the Julia command prompt, navigating to the CRRM folder using `cd("path/to/CRRM")`, and then using the command: `include("scriptName.jl")` where scriptName is the name of the script you wish to run.

## Using R

R isn't necessary for running any of the code, but it can be used to visualizing the results.
We have included two user-interfaces for sifting through the results from the simulation.
To get running with R, the suggested steps are:

- Download & Install [R base](https://www.r-project.org/) and [RStudio](https://www.rstudio.com/products/rstudio/download/)
- Install the pre-requisite R packages by launching RStudio and entering the following in the command terminal:
`install.packages(c('ggplot2','Cairo','knitr','shiny'))`. There are actually more packages that will need to be installed, but that should happen automatically after the next step
- Open either of the files in `./analysis` in RStudio, and then click on 'Run Document' which is above the text-editor area in RStudio (default hotkey: Ctrl+Shift+K in windows)
- The interface should launch in new window and be usable

# Code

## Main functions of interest

For those only interested in the implemented models, the MATLAB implementation for the linear models is in `./mfiles/LRRM.m` and `/mfiles/CLRRM.m`.
The implementation in Julia is in `./jlfiles`.

A quick demo is included in the scripts starting with a01 to a03.

## Code for Simulations & In-vivo Data Analysis

This repository also contains the code used for the simulation and in-vivo data analysis, the results of which was presented in the manuscript.
The code is organized as a series of scripts, where each series represents a particular experiment (i.e. simulation, phantom, clinical) and is meant to be run in linear succession.

Code was tested on Matlab 2015b, Julia 0.4.5, R 3.32, all running on Windows 7 - 64 bit.
An estimated run time for each script is provided, based on how long it takes to run on our lab computer (i7-4770 / 3.40 GHz / 32 Gb RAM).
The simulations take approximately 500 mb of disk space and 1 Gb of RAM to do the linear fitting and 2.5 Gb of RAM to do the non-linear fitting.
The non-linear fit requires more RAM because of multi-threading - the 2.5 Gb estimate is from using 4 threads.
QIN Breast DCE-MRI dataset will use ~4.5 Gb of RAM for certain steps.
The Breast data itself takes 8 Gb of space originally, but it is preprocessed into a smaller size so that it only takes ~1Gb of disk space and the original 8Gb of data can be deleted afterwards.

Folders will be created as necessary, and this behaviour can be changed by modifying `./mfiles/DefaultFolders.m`

### Experiment 01 - Single Voxel Simulation

In this "experiment", the concentration-time curve for a reference tissue and a tumour voxel is simulated.
Noise is added to the tumour voxel to create 10,000 noisy copies of the tumour voxels' concentration-time curve.
These curves are then downsampled and the Reference Region Models are fitted to the simulated data.
Simulation properties can be set in `./mfiles/SimProperties.m`.

Files of interest:
```
Main:
e01_1_simMaker.m
  - creates the simulation data
e01_2a_simProcessNLSQ.jl
  - run by opening Julia command prompt, and using the command: include("path/to/e01_2a_simProcessNLSQ.jl")
  - fits the reference region models to simulated data using non-linear least-squares (i.e. NRRM and CNRRM)
  - results of fitting are saved as .mat files
e01_2b_simProcessLLSQ.m
  - fits the reference region models to simulated data using linear least squares fitting
  - results of fitting are saved as .mat files
e01_3a_collectResults.m
  - compiles a summary of the fits (e.g. mean, std.dev) from e01_2a and e01_2b
  - the summary is exported as a .csv file
e01_3b_collectInformation.m
  - collects additional information from the fitting results, including the run times and the error in the estimated kep,RR

Optional:
e01_x_simAnalyzerStatic.m
  - this script looks at the relationship between error in kepRR and error in CLRRM fit
./analysis/s01_analysisErrKepRR.Rmd
  - a shiny app for visualizing the results of the previous script
```

A summary of the results of the simulation are saved in `./dataResults/e01-simResults-refY.csv`.
These can be visualized by running `./analysis/e01_analysis.Rmd` in RStudio.

### Experiment 02 - QIBA Phantom

For convenience, the QIBA phantom data for the Tofts Model (version 6) is included in `./data/QIBA-ToftsV6/DICOM`.
Refer to https://sites.duke.edu/dblab/qibacontent/ for more details on this dataset.

Files of interest:
```
e02_1_qibaPrep.m
  - converts the phantom DCE-MRI data into concentration data
e02_2a_qibaProcessNLSQ.jl
  - fits the non-linear reference region models to the phantom data
e02_2b_qibaProcessLLSQ.m
  - fits the linear reference region models to the phantom data
e02_3_qibaDisplayResults.m
  - displays the results of the LRRM, NRRM and CRRM fits on the QIBA phantom data
```

### Experiment 03 - QIN Breast DCE-MRI Data

Files of interest:
```
e03_1_unzipData.m
  - unzips the downloaded data for future steps
e03_2_stripData.m
  - removes unnecessary data to reduce file size, and builds the tumour mask
e03_3_shipToJulia.m
  - pre-processed data so that Julia can analyze it with the NRRM and CNRRM
e03_4_qinProcessNLSQ.jl
  - fits the NRRM and CNRRM to the Breast data
e03_5_qinProcessLLSQ.m
  - fir the LRRM and CLRRM to the Breast data
e03_6_collectResults.m
  - collects the results from all the fits and summarizes them

Optional:
e03_x_qinCNR.m
  - Estimates the CNR in the breast datasets
e03_x_seeMaps.m
  - Displays the estimated maps from the fitting steps
  - Might not be fully functional
```
