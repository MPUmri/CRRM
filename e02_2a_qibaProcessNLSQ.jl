# This script processes the QIBA Phantom data using non-linear fitting
# for the original and constrained reference region models

# Run this script by launching julia then entering: include("path/to/e02_2a_qibaProcessNLSQ.jl")

# Estimated runtime: ~14 seconds

# First, figure out the current file path
curFilePath = Base.source_path() # Path to current file
workDir = dirname(curFilePath) # Our working directory (directory of current file)
## If the above doesn't work, then user has to manually enter the location of the CLRRM directory

## Load packages
import DCEMRI
using MAT

# Some required files
matFile = workDir * "/data/QIBA-ToftsV6/QIBAv6-Mini.mat"
warmupFile = workDir * "/data/QIBA-ToftsV6/QIBAv6-Mini4Jl-Warmup.mat"
auxCode = workDir * "/jlfiles/auxCode.jl"
refRegCode = workDir * "/jlfiles/refRegionFunction.jl"

outDir = workDir * "/data/QIBA-ToftsV6/Results"

if !(isdir(outDir))
  mkdir(outDir)
end

# Load pre-requisite code
include(auxCode)

## Warmup run to make sure all is well
results = DCEMRI.fitdata(datafile=warmupFile, models=[2])
rm("output.mat") # Delete output since it isn't needed

# If warmup succeeds, then proceed by loading Reference Region Model
include(refRegCode)

# Load the phantom data
mat = matread(matFile)
Ct = mat["concData"]
Ct = permutedims(Ct, [3, 1, 2])
t = vec(mat["T"])
Cp = vec(mat["aif"])
(sT, sX, sY) = size(Ct)
Ct = reshape(Ct,sT,sX*sY)
mask = ones(sX*sY) .> 0

## Using 'refY' parameters
refKTrans = 0.1;
refVe = 0.1;
refKep = refKTrans/refVe
Crr = DCEMRI.toftskety(t, [refKTrans,refKep], Cp)

# Fit the NRRM and CNRRM to the QIBA phantom data
(paramsCN, residCN, estKepRR, paramsN, residN, runtimeN, runtimeCN)=fitCNRRM(Ct, Crr, t, estKepRR=0.0, doTime=true)

# Export to .mat
results = Dict()
  results["pkParamsN"] = paramsN
  results["pkParamsCN"] = paramsCN
  results["residN"] = residN
  results["residCN"] = residCN
  results["estKepRRN"] = estKepRR
  results["refKt"] = refKTrans
  results["refVe"] = refVe
  results["refKep"] = refKep
  matwrite(outDir * "/QIBAv6-refY-NL.mat", results)
