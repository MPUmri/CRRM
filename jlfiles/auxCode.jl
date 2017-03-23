# This file contains some auxiliary code that is needed to process the DCE-MRI data through Julia

# This function will do non-linear fitting for the constrained reference region model
function fitCNRRM{N}(Ct::Array{Float64,N}, Crr::Vector{Float64}, t::Vector{Float64}; estKepRR=0.0, doTime=false)
  # Options:
  # estKepRR - Estimate for kepRR to use with CNRRM. If estKepRR<=0, then it will be estimated by first using NRRM
  # doTime - Records how long it takes to do the fitting
  DCEMRI.@dprint "fitting DCE data"
  (nT,nV) = size(Ct)
  @assert nT == length(t)
  mask = ones(nV) .> 0
  idxs = find(mask)
  nidxs = length(idxs)
  didNRRM = false
  if (estKepRR<=0.0)
    # First do NRRM, then get kepRR estimate, then do CNRRM
    didNRRM = true
    DCEMRI.@dprint "attempting Unconstrained Reference Region Tofts-Kety model"

    paramsNRRM = zeros(4, nV)
    residNRRM = Inf*ones(nV)
    runtimeNRRM = 0.0
    (paramsNRRM, residNRRM, runtimeNRRM) = fitNRRM(Ct,Crr,t,doTime=true)
    maskNegative = paramsNRRM .<= 0
    maskNegative = sum(maskNegative,1)
    maskPositive = vec(maskNegative .== 0)
    # Apply the mask, and obtain the median from the voxel-wise kepRR estimates
    paramsTmp = paramsNRRM[:,maskPositive]
    estKepRR = median(vec(paramsTmp[4,:]))
  end

  DCEMRI.@dprint "attempting Constrained Reference Region Tofts-Kety model"
  if doTime
    tic()
  end
  paramsCNRRM = zeros(3,nV)
  residCNRRM = Inf*ones(nV)
  p0 = [0.01, 0.01]
  fr2(x,p) = CRRM(x, p, Crr) # Cp is reference region curve. NOT AIF.
  p, r, dof = DCEMRI.nlsfit(fr2, Ct, idxs, [t; [estKepRR]], p0)
  r = squeeze(sumabs2(r, 1), 1) / dof
  for k in idxs
    residCNRRM[k] = r[k]
    paramsCNRRM[1,k] = p[1,k] # Relative KTrans
    paramsCNRRM[2,k] = estKepRR*p[1,k]/p[2,k] # Relative ve
    paramsCNRRM[3,k] = p[2,k] # kep
	end
  if doTime
    runtimeCNRRM = toc()
  end
  if didNRRM
    if doTime
      return (paramsCNRRM, residCNRRM, estKepRR, paramsNRRM, residNRRM, runtimeNRRM, runtimeCNRRM)
    else
      return (paramsCNRRM, residCNRRM, estKepRR, paramsNRRM, residNRRM)
    end
  else
    if doTime
      return (paramsCNRRM, residCNRRM, estKepRR, runtimeCNRRM)
    else
      return (paramsCNRRM, residCNRRM, estKepRR)
    end
  end
end

# This function will do non-linear fitting for the unconstrained reference region model
function fitNRRM{N}(Ct::Array{Float64,N}, Crr::Vector{Float64}, t::Vector{Float64}; doTime=false)
  DCEMRI.@dprint "fitting DCE data"
  (nT,nV) = size(Ct)
  @assert nT == length(t)

  mask = ones(nV) .> 0
  idxs = find(mask)
  nidxs = length(idxs)

  resid = Inf*ones(nV)
  params = zeros(4, nV)

  DCEMRI.@dprint "attempting Unconstrained Reference Region Tofts-Kety model"
  if doTime
    tic()
    runTime=0.0
  end
  p0 = [0.01, 0.01, 0.01]
  fr1(x,p) = RRM(x, p, Crr)
  p, r, dof = DCEMRI.nlsfit(fr1, Ct, idxs, t, p0)
  r = squeeze(sumabs2(r, 1), 1) / dof
  for k in idxs
    resid[k] = r[k]
    params[1,k] = p[1,k] # Relative KTrans
    params[2,k] = p[2,k]*p[1,k]/p[3,k] # Relative ve
    params[3,k] = p[3,k] # kep or (KTrans/ve)
    params[4,k] = p[2,k] #kep,RR
  end
  runTime = toc()
  if doTime
    return (params, resid, runTime)
  else
    return (params, resid)
  end
end

function fitLRRM(Ct::Matrix{Float64}, Crr::Vector{Float64}, t::Vector{Float64}, getResid=0)
	# Linear reference region model (based on Tofts-Kety)
	# Reference: J Cardenas-Rodriguez et al. MRI 31 (2013) 497-507
	# Input:
	#		t = time points of DCE acquisition. Vector of length nT.
	#		Ct = cocentration-time data. Matrix of dimensions nT-by-N, where N is #voxels.
	#   Crr = concentration-time curve of reference region. Vector of length nT.
	#		getResid = whether to return residuals in results. (1=True, 0=False)
	# Output:
	#		p = Fit parameters. 3-by-N matrix, where N is number of voxels
	#			The 3 elements of p are: [relativeKTrans, relativeVe, kep]
	#   resid = Residuals for each voxel. Vector of length N.
  (nT, nV) = size(Ct)
  resid = zeros(nV)
  trapzCrr = zeros(nT)
  p=zeros(3,nV)
  trapzCrr=cumtrapz(t,Crr)
  for i in 1:nV
    M = zeros(nT, 3);
    curCt = Ct[:,i]
    M[:,1] = Crr
    M[:,2] = trapzCrr
    M[:,3] = -cumtrapz(t,curCt)
    p[:,i] = M\curCt
		resid[i] = norm(curCt-M*p[:,i])
		# p = [kTrans_TOI/kTrans_RR, kTrans_TOI/ve_RR, kTrans_TOI/ve_TOI]
		# desire: p = [kTrans_TOI/ktrans_RR, ve_TOI/ve_RR, kTrans_TOI/ve_TOI]
  end
  p[2,:] = p[2,:]./p[3,:]
  if getResid==0
    return p
  elseif getResid==1
    return p,resid
  end
end

function fitCLRRM(t::Vector{Float64}, Ct::Matrix{Float64}, Crr::Vector{Float64}; getResid=0, estKepRR=1.0)
	# Constrained Linear Reference Region Model
	# Reference: Ahmed & Levesque. (2016), MRM, http://doi.org/10.1002/mrm.26530
	# Input:
	#		t = time points of DCE acquisition. Vector of length nT.
	#		Ct = cocentration-time data. Matrix of dimensions nT-by-N, where N is #voxels.
	#   Crr = concentration-time curve of reference region. Vector of length nT.
	#		getResid = whether to return residuals in results. (1=True, 0=False)
  #   estKepRR = estimate for the kepRR which will be held fixed
	# Output:
	#		p = Fit parameters. 3-by-N matrix, where N is number of voxels
	#			The 3 elements of p are: [relativeKTrans, relativeVe, kep]
	#   resid = Residuals for each voxel. Vector of length N.
  (nT, nV) = size(Ct)
  resid = zeros(nV)
  trapzCrr = zeros(nT)
  p=zeros(2,nV)
  pOut = zeros(3,nV)
  trapzCrr=cumtrapz(t,Crr)
  for i in 1:nV
    M = zeros(nT, 2);
    curCt = Ct[:,i]
    M[:,1] = Crr + estKepRR* trapzCrr
    M[:,2] = -cumtrapz(t,curCt)
    p[1:2,i] = M\curCt
		resid[i] = norm(curCt-M*p[:,i])
  end
  pOut[3,:] = p[2,:]
  pOut[2,:] = estKepRR*p[1,:]
  pOut[1,:] = p[1,:]
  # right now, pOut = [kTrans_TOI/kTrans_RR, kTrans_TOI/ve_RR, kTrans_TOI/ve_TOI]
  # we desire: pOut = [kTrans_TOI/ktrans_RR, ve_TOI/ve_RR, kTrans_TOI/ve_TOI]
  pOut[2,:] = pOut[2,:]./pOut[3,:]
  if getResid==0
    return pOut
  elseif getResid==1
    return pOut,resid
  end
end


function trapz{T<:Real}(x::Vector{T}, y::Vector{T})
	# Trapezoidal integration
	# Reference: https://groups.google.com/d/msg/julia-users/CNYaDUYog8w/lveUAdjkjmIJ
	# Input:
	# 	x = x values
	# 	y = y values of f(x)
	# Output:
	#		r = integral of f(x) from x[1] to x[end] using trapezoidal rule
	if (length(y) != length(x))
	   error("Vectors must be of same length")
	end
	r = 0.0
	for i in 2:length(x)
	    r += (x[i] - x[i-1]) * (y[i] + y[i-1])
	end
	return r/2.0
end

function cumtrapz{T<:Real}(x::Vector{T}, y::Vector{T})
	# Cumulative Trapezoidal integration
	# Input:
	# 	x = x values
	# 	y = y values of f(x)
	# Output:
	#		r = vector, where r(k) is integral of f(x) from x[1] to x[k]
  if (length(y) != length(x))
	   error("Vectors must be of same length")
	end
  r = zeros(length(x))
  for k in 1:length(x)
    r[k]=trapz(x[1:k],y[1:k])
  end
  return r
end

# This is a simple way of downsampling by omitting data points
function downsample(x::Vector{Float64}, dx::Int)
	samplePoints = collect(1:dx:length(x))
	numPoints = length(samplePoints)
	newX = zeros(numPoints)
	for i in 1:numPoints
		newX[i] = x[samplePoints[i]]
	end
	return(newX)
end

function downsample(x::Array{Float64}, dx::Int)
	(nT,nV) = size(x)
	samplePoints = collect(1:dx:nT)
	numPoints = length(samplePoints)
	newX = zeros(numPoints,nV)
	for i in 1:numPoints
		newX[i,:] = x[samplePoints[i], :]
	end
	return(newX)
end
