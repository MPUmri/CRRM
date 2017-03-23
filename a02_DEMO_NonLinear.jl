# This script process the simulation data using non-linear least squares fitting
# for both the original Reference Region Model (NRRM) and the Constrained approach (CNRRM)

# Run this script by launching julia then entering: include("path/to/a02_DEMO_NonLinear.jl")

# Get current file path
curFilePath = Base.source_path() # Path to current file
workDir = dirname(curFilePath) # Our working directory (directory of current file)
## If the above doesn't work, then user has to manually enter the location of the CRRM directory
## i.e. wordDir = "/path/to/CLRRM/Directory"

# Figure out how many cores we have to work with
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
warmupFile = workDir * "/data/QIBA-ToftsV6/QIBAv6-Mini4Jl-Warmup.mat"
auxCode = workDir * "/jlfiles/auxCode.jl"
refRegCode = workDir * "/jlfiles/refRegionFunction.jl"

inFile = workDir * "/data/demoData4Julia.mat" # File to load data from
outFile = workDir * "/data/demoDataFromJulia.mat" # File to save results to

# Load pre-requisite code
include(auxCode)

## Warmup run to make sure all is well
println("==== Warmup Run ====")
results = fitdata(datafile=warmupFile, models=[2])
rm("output.mat") # Delete output since it isn't needed

# If warmup succeeds, then proceed by loading Reference Region Models
include(refRegCode)

println("==== Processing data ====")

# Load .mat data
matData = matread(inFile)
# Assign data from .mat file to variables in Julia, for simplicity
Ct = matData["Ct"]
Crr = vec(matData["Crr"])
t = vec(matData["t"])

# Do the fitting
# Non-Linear Reference Region Model
# (Unconstrained and Constrained are both run if estKepRR is set to zero)
(pkParamsCN, residCN, estimatedKepRR, pkParamsN, residN, runtimeN, runtimeCN) = fitCNRRM(Ct, Crr, t, estKepRR=0.0, doTime=true)

# Output results to a .mat file
results = Dict()
results["estParamsNRRM"] = pkParamsN
results["estParamsCNRRM"] = pkParamsCN
results["runtimeNRRM"] = runtimeN
results["runtimeCNRRM"] = runtimeCN
results["estKepRR_N"] = estimatedKepRR
matwrite(outFile, results)

println("==== Finished ====")
println("==================")
