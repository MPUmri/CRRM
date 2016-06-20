# This script process the simulation data using non-linear least squares fitting
# for both the original Reference Region Model (NRRM) and the Constrained approach (CNRRM)

# Pick the choice of parameters for the reference tissue parameters
# (This should match the choice from step e01_1)
refName = "refY"

# The temporal resolution (in seconds) to process over
listTRes = [10]
# The choice of Contrast-Noise Ratios
# (This does not change noise, rather it is only here to remind the code what the CNRs are)
listCNR = [5]

# First, figure out the current file path
curFilePath = Base.source_path() # Path to current file
workDir = dirname(curFilePath) # Our working directory (directory of current file)
## If the above doesn't work, then user has to manually enter the location of the CLRRM directory
## i.e. wordDir = "/path/to/CLRRM/Directory"

## Load packages
import DCEMRI  # Needed for the levenberg-marquardt fitting
using MAT      # Needed for loading/saving .mat files

# Some required directories/files
matDir = workDir * "/data/simData/$refName/rawData"
warmupFile = workDir * "/data/QIBA-ToftsV6/QIBAv6-Mini4Jl-Warmup.mat"
auxCode = workDir * "/jlfiles/auxCode.jl"
refRegCode = workDir * "/jlfiles/refRegionFunction.jl"

outDir = workDir * "/data/simData/$refName/NRRM"

if !(isdir(outDir))
  mkdir(outDir)
end

# cd(outDir) # unnecessary

# Load pre-requisite code
include(auxCode)

## Warmup run to make sure all is well
println("==== Warmup Run ====")
results = DCEMRI.fitdata(datafile=warmupFile, models=[2])
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
  curCNR = listCNR[q]
  println("")
  println("------------")
  println("Processing CNR = $curCNR")
  curFile = matDir * "/Sim-CNR-$curCNR.mat"
  matData = matread(curFile)
  cnrInd = listCNR[q]
  for i=1:length(listTRes)
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
    mask = ones(nV) .> 0
    # Non-Linear Reference Region Model
    tic()
    pkParamsN, residN = doRRM(Ct, mask, t, Crr, constrain=false)
    runtimeN = toc()
    # Build a mask to filter out any voxel where any estimate is negative
    maskNegative = pkParamsN .<= 0
    maskNegative = sum(maskNegative,1)
    maskPositive = vec(maskNegative .== 0)
    # Apply the mask, and obtain the interquartile mean from the voxel-wise kepRR estimates
    pkParams = pkParamsN[:,maskPositive]
    estKepRR = vec(pkParams[4,:])
    meanKepRR = mean(iqrFilter(estKepRR,[.25, .75], true))
    # Use the interquartile mean-based kepRR with the Constrained Non-Linear Reference Region Model
    println("")
    tic()
    pkParamsCN, residCN = doRRM(Ct, mask, t, Crr, estKepRR=meanKepRR, constrain=true)
    runtimeC = toc()
    # Output results
    results = Dict()
    results["pkParamsN"] = pkParamsN
    results["residN"] = residN
    results["pkParamsCN"] = pkParamsCN
    results["residCN"] = residCN
    results["runtimeN"] = runtimeN
    results["runtimeC"] = runtimeC
    results["meanKepRR"] = meanKepRR
    results["stdKepRR"] = std(iqrFilter(estKepRR,[.25, .75], true))
    matwrite(outDir * "/Sim-CNR-$curCNR-TRes-$curTRes.mat", results)
  end
end
println("==== Finished ====")
println("==================")
