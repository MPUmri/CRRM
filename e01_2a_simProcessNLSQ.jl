# This script process the simulation data using non-linear least squares fitting
# for both the original Reference Region Model (NRRM) and the Constrained approach (CNRRM)

# Run this script by launching julia then entering: include("path/to/e01_2a_simProcessNLSQ.jl")

# Estimated runtimes:
# ~2.6 hours for NRRM
# ~1.2 hours for CNRRM
# So total run time is almost 4 hours
# Note: These runtimes were obtained using 4 threads.
#       If only a single core is used, then the runtime could be 4 times longer

# Pick the choice of parameters for the reference tissue parameters
# (This should match the choice from step e01_1)
refName = "refY"

# The temporal resolution (in seconds) to process over
listTRes = [1,5,10,15,30,60]
# The choice of Contrast-Noise Ratios
# (This does not change noise, rather it is only here to remind the code what the CNRs are. Should match values from simMaker.)
listCNR = collect(5:5:50)

# First, figure out the current file path
curFilePath = Base.source_path() # Path to current file
workDir = dirname(curFilePath) # Our working directory (directory of current file)
## If the above doesn't work, then user has to manually enter the location of the CRRM directory
## i.e. wordDir = "/path/to/CLRRM/Directory"

import Hwloc
topology = Hwloc.topology_load()
counts = Hwloc.histmap(topology)
nCores = counts[:Core]
if nprocs()<nCores
  addprocs(nCores-nprocs()+1)
end

## Load packages
using DCEMRI  # Needed for the levenberg-marquardt fitting
using MAT      # Needed for loading/saving .mat files

# Some required directories/files
matDir = workDir * "/data/simData/$refName/rawData"
warmupFile = workDir * "/data/QIBA-ToftsV6/QIBAv6-Mini4Jl-Warmup.mat"
auxCode = workDir * "/jlfiles/auxCode.jl"
refRegCode = workDir * "/jlfiles/refRegionFunction.jl"

outDir = workDir * "/data/simData/$refName/NRRM"

# Make the output directory if it doesn't already exist
if !(isdir(outDir))
  mkdir(outDir)
end

# Load pre-requisite code
include(auxCode)

## Warmup run to make sure all is well
println("==== Warmup Run ====")
results = fitdata(datafile=warmupFile, models=[2])
rm("output.mat") # Delete output since it isn't needed

# If warmup succeeds, then proceed by loading Reference Region Model
include(refRegCode)

# Define the reference tissue parameters
# Default is 'refY'
refKTrans = 0.1;
refVe = 0.1;
if (refName == "refWS")
  refKTrans = 0.07;
  refVe = 0.14;
elseif (refName == "refP")
  refKTrans = 0.137;
  refVe = 0.115;
end
refKep = refKTrans/refVe

# Get the names of the mat files containing simulated data
matFiles = readdir(matDir)
numFiles = length(matFiles)

println("")
println("===============================")
println("==== Start processing data ====")
# Loop through the simulated datasets
for q=1:numFiles
  # Load the simulated data for each CNR
  curCNR = listCNR[q]
  println("")
  println("------------")
  println("Processing CNR = $curCNR")
  curFile = matDir * "/Sim-CNR-$curCNR.mat"
  matData = matread(curFile)
  cnrInd = listCNR[q]
  for i=1:length(listTRes)
    # Appropriately downsample simulated data to obtain desired temporal resolutions
    curTRes = listTRes[i]
    println("")
    println("---")
    println("Temporal Resolution = $curTRes s")
    # Downsample the data
    Ct = downsample(matData["CtNoisy"],curTRes)
    Crr = downsample(vec(matData["CrrClean"]),curTRes)
    t = downsample(vec(matData["T"]),curTRes)
    # Build a dummy mask (ones everywhere)
    (nT, nV) = size(Ct)
    # Non-Linear Reference Region Model
    # (Unconstrained and Constrained are both run if estKepRR is set to zero)
    (pkParamsCN, residCN, medianKepRR, pkParamsN, residN, runtimeN, runtimeC) = fitCNRRM(Ct, Crr, t, estKepRR=0.0, doTime=true)
    # Output results
    results = Dict()
    results["pkParamsN"] = pkParamsN
    results["residN"] = residN
    results["pkParamsCN"] = pkParamsCN
    results["residCN"] = residCN
    results["runtimeN"] = runtimeN
    results["runtimeC"] = runtimeC
    results["medianKepRR"] = medianKepRR
    matwrite(outDir * "/Sim-CNR-$curCNR-TRes-$curTRes.mat", results)
  end
end
println("==== Finished ====")
println("==================")
