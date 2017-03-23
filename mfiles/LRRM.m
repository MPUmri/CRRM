function [pkParams, resid] = LRRM(Ct, Crr, t, doPure)
    % Fits the Linear Reference Region Model
  	% Reference: J Cardenas-Rodriguez et al. MRI 31 (2013) 497-507
  	% Inputs:
  	%		Ct = cocentration-time data. Matrix of dimensions nT-by-N, where nT is number of time points, and N is #voxels.
  	%   Crr = concentration-time curve of reference region. Vector of length nT.
  	%		t = time points of DCE acquisition. Vector of length nT.
  	% Output:
  	%		pkParams = Fitted parameters. 3-by-N matrix, where N is number of voxels
  	%		doPure = Boolean.
    %            If true, pkParams will return the fitted parameters,
    %            If false, pkParams will be modified to provide useful parameters.
    %            -- if false, the 3 elements of pkParams are: [relativeKTrans, relativeVe, kep]
  	%   resid = Norm of residuals for each voxel. Vector of length N.


    % Set doPure to false be default
    if nargin<4
        doPure = false;
    end

    stepSize = t(2)-t(1); % Assuming constant step size throughout acquisition

    % Initialize matrices
    [sT, sX] = size(Ct);

    M1 = zeros(sT,1);
    M2 = M1;

    M1 = Crr;
    M2 =  stepSize*cumtrapz(Crr);

    resid = zeros(sX,1);
    pkParams = zeros(sX,3);

    % Do the linear least-square fit for each voxel in Ct
    % Disabling warnings for speed (possible non-tissue regions in image give warnings)
    warning off

    for i=1:sX
            curCt = squeeze(Ct(:,i));
            M3 = -stepSize*cumtrapz(curCt);
            M = [M1, M2, M3];
            pkParams(i,:) = mldivide(M,curCt);
            resid(i) = norm(curCt-M*pkParams(i,:)');
    end
    % Turn warnings back on
    warning on

    % The fitted values of pkParams are: [kTrans_TOI/kTrans_RR, kTrans_TOI/ve_RR, kTrans_TOI/ve_TOI]
    % but it'd be more useful if we had: [kTrans_TOI/ktrans_RR, ve_TOI/ve_RR, kTrans_TOI/ve_TOI]
    if doPure==false
        pkParams(:,2) = pkParams(:,2)./pkParams(:,3);
    end
end
