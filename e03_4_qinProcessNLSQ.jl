# Process QIN Breast Data using NRRM and CNRRM

# Run this script by launching julia then entering: include("path/to/e03_4_qinProcessNLSQ.jl")

# Estimated runtime: ~500 seconds, on 4 cores

# Choice of parameters for the reference tissue. Should match e03_05
# It's best to leave this as it is
refName = "refY"

# First, figure out the current file path
curFilePath = Base.source_path() # Path to current file
workDir = dirname(curFilePath) # Our working directory (directory of current file)
## If the above doesn't work, then user should enter the location of the CLRRM directory
## by manually modifying this script

# Load pre-requisite packages
using DCEMRI
using MAT

# Get the directories and files of interest
matDir = workDir * "/data/QINBreast/ForJulia"
warmupFile = workDir * "/data/QIBA-ToftsV6/QIBAv6-Mini4Jl-warmup.mat"
auxCode = workDir * "/jlfiles/auxCode.jl"
refRegCode = workDir * "/jlfiles/refRegionFunction.jl"

outDir = workDir * "/data/QINBreast/FromJulia"

if !(isdir(outDir))
  mkdir(outDir)
end

############### Warmup run
results = fitdata(datafile=warmupFile, models=[2])
rm("output.mat") # Delete output since it isn't needed
##############

# This will load the functions for fitting the models
include(auxCode)
include(refRegCode)

# Define the reference tissue parameters
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

############## Main run
matFiles = readdir(matDir)
numFiles = length(matFiles)

# Loop through the .mat files
for q=1:numFiles
  curName = matFiles[q]
  curFile = matDir * "/$curName"

  # Load data
  matData = matread(curFile)
  Ct = matData["Ct"]
  (nT, nV) = size(Ct)
  Ct = reshape(Ct, nT, nV, 1)
  Cp = vec(matData["Cp"])
  t = vec(matData["t"])
  tResFactor = round(Int64,matData["tResFactor"])
  (nT, nV) = size(Ct)
  mask = ones(nV,1) .> 0

  # Simulate reference tissue curve
  Crr = DCEMRI.toftskety(t, [refKTrans, refKep], Cp)
  Ct = reshape(Ct,nT,nV)
  # Downsample so that Crr, Cp, and t have same temporal resolution as Ct
  Crr = downsample(Crr,tResFactor)
  Cp = downsample(Cp,tResFactor)
  t = downsample(t,tResFactor)

  # Non-Linear Tofts Model
  pkParamsT, residT = DCEMRI.fitdce(Ct, mask, t, Cp, models=[2])

  # Re-arrange the output from Tofts Model fit
  kt = pkParamsT[1,:,1]
  ve = pkParamsT[2,:,1]
  kep = kt./ve

  # Parameter mask using threshold on fitted Tofts Model parameters
  pMask = (kt.>0.001) & (kt.<1) & (ve.>0.001) & (ve.<1) & (kep.<35)

  kt = kt[pMask]
  ve = ve[pMask]
  kep = kep[pMask]

  # NRRM and CNRRM
  (pkParamsCN, residCN, estKepRR, pkParamsN, residN, runtimeN, runtimeCN)=fitCNRRM(Ct, Crr, t, estKepRR=0.0, doTime=true)

  # Output results
  results = Dict()
  results["pkParamsN"] = pkParamsN
  results["residN"] = residN
  results["pkParamsCN"] = pkParamsCN
  results["residCN"] = residCN
  results["runtimeN"] = runtimeN
  results["runtimeC"] = runtimeCN
  results["estKepRRN"] = estKepRR
  results["pkParamsT"] = pkParamsT
  results["residT"] = residT
  results["kt"] = kt
  results["ve"] = ve
  results["kep"] = kep
  results["Crr"] = Crr
  results["Cp"] = Cp
  results["t"] = t
  results["pMask"] = convert(Array{Bool,1}, pMask)
  matwrite(outDir * "/$curName", results)
end
