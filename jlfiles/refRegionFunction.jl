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
	# Ct = kR * Crr + kR * (kepRR-kep) * conv(Crr,exp(-kep*t)) --- but kepRR is held fixed
	# Reference: Ahmed & Levesque. (2016), MRM, http://doi.org/10.1002/mrm.26530
  kR = p[1] # relative kTrans
  kep = p[2] # kep of TOI
	# This part is ugly. We needed to pass an extra fixed parameter (i.e. kepRR)
	# but I didn't want to modify the nlsfit() function...
	# ... So, I snuck kepRR as the last element of time.
  kepRR = tmpT[end] # use last element of time as kepRR (because that's how we smuggled it into this function)
  t=tmpT[1:end-1] # delete that last element from time because it doesn't represent an actual time point
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
