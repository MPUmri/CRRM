# This file contains some auxiliary code that is needed to process the DCE-MRI data through Julia

# This function will do non-linear fitting for the reference region model
function doRRM{N}(Ct::Array{Float64,N}, mask::BitArray, t::Vector{Float64},
                   Crr::Vector{Float64}; estKepRR=0.0, constrain=false)
  DCEMRI.@dprint "fitting DCE data"
  (nT,nV) = size(Ct)
  @assert nT == length(t)

  idxs = find(mask)
  nidxs = length(idxs)

  resid = Inf*ones(nV)
  params = zeros(4, nV)
  if constrain
  	params = zeros(3,nV)
  end

  if !constrain
    DCEMRI.@dprint "attempting Unconstrained Reference Region Tofts-Kety model"
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
  else
    DCEMRI.@dprint "attempting Constrained Reference Region Tofts-Kety model"
    p0 = [0.01, 0.01]
    fr2(x,p) = CRRM(x, p, Crr) # Cp is reference region curve. NOT AIF.
    p, r, dof = DCEMRI.nlsfit(fr2, Ct, idxs, [t; [estKepRR]], p0)
    r = squeeze(sumabs2(r, 1), 1) / dof
    for k in idxs
	    resid[k] = r[k]
	    params[1,k] = p[1,k] # Relative KTrans
	    params[2,k] = estKepRR*p[1,k]/p[2,k] # Relative ve
	    params[3,k] = p[2,k] # kep
  	end
  end
  (params, resid)
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

function iqrFilter(x::Array{Float64,1}, qrValue::Vector, doPositive=false)
	newX = x
	if doPositive
		newX = x[x.>0]
	end
	qrRange = quantile(newX,qrValue)
	newX = newX[newX.>qrRange[1]]
	newX = newX[newX.<qrRange[2]]
	return(newX)
end
