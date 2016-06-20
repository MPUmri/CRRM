# CRRM

This repository contains a small example of the Linear Reference Region Model (LRRM), Non-linear Reference Region Model (NRRM), and their constrained variants (CLRRM and CNRRM).

## Description of main .m files

- e01_1_simMaker.m
    + Will make a quick simulation of 10,000 voxels with a CNR of 5 and temporal resolution of 10s
- e01_2a_simProcessNLSQ.jl
    + Process the simulated data using the NRRM and CNRRM
    + Requires Julia  and the [DCEMRI.jl package](https://github.com/davidssmith/DCEMRI.jl)
    + To run, launch julia, cd() to this directory, then run `include{"e01_2a_simProcessNLSQ.jl"}`
    + If you don't have Julia, you can skip it
- e01_2b_simProcessLLSQ.m
    + Process the simulated data using the LRRM and CLRRM
- e01_3a_collectResults.m
    + Plots boxplot of the fitting results

The implementation of the LRRM and CLRRM are in `./mfiles/LRRM.m` and `/mfiles/CLRRM.m`.

