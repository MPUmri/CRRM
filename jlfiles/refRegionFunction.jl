# This file contains the implementation of the Reference Region Model
# and Constrained Reference Region Model, which are to be used with non-linear
# least-squares fitting

@everywhere function RRM(t::Vector{Float64}, p::Vector{Float64}, Crr::Vector{Float64})
  # Reference Region Model
  # Tofts modle using a reference tissue curve (Crr) instead of AIF (Cp)
	# Ct = kR * Crr + kR * (kepRR-kep) * conv(Crr,exp(-kep*t))
	# Reference: Yankeelov et al. (2005), MRI, 23(4), 519–29. doi:10.1016/j.mri.2005.02.013
  kR = p[1] # relative kTrans
  kepRR = p[2] # kep of reference tissue
  kep = p[3] # kep of TOI
  Ct = zeros(t)
  for k in 1:length(Ct)
    tk = t[k]
    @simd for j = 1:k
      @inbounds y = exp(kep*(t[j] - tk)) * Crr[j]
      @inbounds Ct[k] += ifelse((j == 1 || j == k) && k > 1, 0.5*y, y)
    end
  end
  dtp = (t[2] - t[1])
  @simd for k = 1:length(Ct)
    @inbounds Ct[k] = kR * (Crr[k] + dtp*(kepRR-kep)*Ct[k])
  end
  Ct
end

@everywhere function CRRM(tmpT::Vector{Float64}, p::Vector{Float64}, Crr::Vector{Float64})
  # Constrained Reference Region Model
  # Tofts modle using a reference tissue curve (Crr) instead of AIF (Cp)
	# Ct = kR * Crr + kR * (kepRR-kep) * conv(Crr,exp(-kep*t))
	# Reference: Me?
  kR = p[1] # relative kTrans
  kep = p[2] # kep of TOI
	# This part is ugly. We needed to pass an extra fixed parameter (i.e. kepRR)
	# but I didn't want to modify the nlsfit() function...
	# ... So, I snuck kepRR as the last element of time.
  kepRR = tmpT[end] # use last element of time as kepRR (because that's how we smuggled it into this function)
  t=tmpT[1:end-1] # delete that last element from time because it actually doesn't represent a time point
  # The rest of the function is similar to RRM()
  Ct = zeros(t)
  for k in 1:length(Ct)
    tk = t[k]
    @simd for j = 1:k
      @inbounds y = exp(kep*(t[j] - tk)) * Crr[j]
      @inbounds Ct[k] += ifelse((j == 1 || j == k) && k > 1, 0.5*y, y)
    end
  end
  dtp = (t[2] - t[1])
  @simd for k = 1:length(Ct)
    @inbounds Ct[k] = kR * (Crr[k] + dtp*(kepRR-kep)*Ct[k])
  end
  Ct
end

function fitLRRM(t::Vector{Float64}, Ct::Matrix{Float64}, Crr::Vector{Float64}, getResid=0)
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

function fitCLRRM(t::Vector{Float64}, Ct::Matrix{Float64}, Crr::Vector{Float64}, getResid=0, estKepRR=1.0)
	# Constrained Linear Reference Region Model
	# Reference: Me?
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

function llsqRR(t::Vector{Float64}, Ct::Matrix{Float64}, cRR::Vector{Float64}, mode=1)
  (nT, nV) = size(Ct)
  resid = zeros(nV)
  trapzCrr = zeros(nT)
  p=zeros(3,nV)
  for k in 1:nT
    trapzCrr[k]=trapz(t[1:k],cRR[1:k])
  end
  for i in 1:nV
    M = zeros(nT, 3);
    curCt = Ct[:,i]
    for k in 1:nT
      M[k,1] = cRR[k]
      M[k,2] = trapzCrr[k]
      M[k,3] = -trapz(t[1:k],curCt[1:k]);
    end
    x = M\curCt
    resid[i] = norm(curCt-M*x)
    p[:,i]=x
  end
  if mode==0
    return p,resid
  elseif mode==1
    return p
  end
end
