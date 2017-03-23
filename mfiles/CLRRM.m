function [pkParams, resid, estKepRR, stdKepRR, p] = CLRRM(Ct, Crr, t, kepRef)
    % Constrained Linear Reference Region Model
  	% Reference: Ahmed & Levesque. (2016), MRM, http://doi.org/10.1002/mrm.26530
  	% Input:
  	%		Ct = cocentration-time data. Matrix of dimensions nT-by-N, where N is #voxels.
  	%   Crr = concentration-time curve of reference region. Vector of length nT.
    %		t = time points of DCE acquisition. Vector of length nT.
    %   kepRef = Estimated value for the reference region kep value
    %            If kepRef > 0, then this value will be used as kepRef
    %                      = 0, then median of positive values will be used as kepRef
    %                      = -1, then mean will be used
    %                      = -2, then interquartile mean will be used as kepRef
    % Output:
  	%		pkParams = Fitted parameters. 3-by-N matrix, where N is number of voxels.
    %              The 3 elements of pkParams are: [relativeKTrans, relativeVe, kep]
  	%   resid = Norm of residuals for each voxel. Vector of length N.
    %   estKepRR = estimated value for kepRR (aka kepRef)
    %              if kepRef > 0, then estKepRR = kepRef
    %              otherwise, estKepRR is the value estimated kepRef from LRRM
    %   stdKepRR = std deviation of estimate kepRR values
    %              if kepRef > 0, then stdKepRR = nan
    %              otherwise, stdKepRR is spread of estimated kepRef from LRRM
    %   p = Estimated parameters from LRRM, if LRRM fit was done to estimate kepRef. Otherwise, p=nan.


    % Use interquartile mean by default
    if nargin<4
        kepRef = 0;
    end

    if kepRef <=0
        % If kepRef is not known, then estimate it from the LRRM
        p = LRRM(Ct,Crr,t,1);
        p1 = p(:,1); % First fitted parameter from LRRM
        p2 = p(:,2); % Second fitted parameter from LRRM
        goodVals = (p1>0) & (p2>0) & p(:,3)>0; % Good values are when all parameters are positive
        % Get kepRR, which is ratio of p2/p1, for voxels where all fitted parameters are positive
        x= p2(goodVals)./p1(goodVals);
        % Use either the mean, interquartile mean, or median (default) to obtain a single value for kepRR
        if kepRef == -1
            % Using the mean estimated kepRR from LRRM
            estKepRR = mean(x);
            stdKepRR = std(x);
        elseif kepRef == -2
            % Using the interquartile mean of kepRR from LRRM
            qtRange = quantile(x,[.25 .75]);
            xMask = x>qtRange(1) & x<qtRange(2);
            estKepRR = mean(x(xMask));
            stdKepRR = std(x(xMask));
        else
            % Using the median of positive kepRR estimates
            estKepRR = median(x(x>0));
            qtRange = quantile(x,[.25 .75]);
            stdKepRR = (qtRange(2)-qtRange(1))/1.35; % Estimate stdKepR using interquartile range
        end
    else
      % If a positive value of kepRef is provided, then use that as kepRR
        estKepRR = kepRef;
        stdKepRR = nan;
        p = nan;
    end

    % Assuming constant step size throughout acquisition
    stepSize = t(2)-t(1);

    % Initialize matrices
    [sT, sX] = size(Ct);

    M1 = zeros(sT,1);
    M2 = M1;

    M1 = Crr;
    M2 = stepSize*cumtrapz(Crr);

    resid = zeros(sX,1);
    pkParams = zeros(sX,2);

    % Do the linear least-squares fit
    % Disabling warnings for speed (possible non-tissue regions in image give warnings)
    warning off

    for i=1:sX
        curCt = squeeze(Ct(:,i));
        M3 = -stepSize*cumtrapz(curCt);
        M = [M1+estKepRR*M2, M3];
        pkParams(i,:) = mldivide(M,curCt);
        resid(i) = norm(curCt-M*pkParams(i,:)');
    end
    warning on

    % Modify the CLRRM fit parameters so that they match the LRRM's parameters
    pkParams(:,3)=pkParams(:,2);
    pkParams(:,2)=estKepRR*pkParams(:,1);

    % The fitted values of pkParams are: [kTrans_TOI/kTrans_RR, kTrans_TOI/ve_RR, kTrans_TOI/ve_TOI]
    % but it'd be more useful if we had: pkParams = [kTrans_TOI/ktrans_RR, ve_TOI/ve_RR, kTrans_TOI/ve_TOI]
    pkParams(:,2) = pkParams(:,2)./pkParams(:,3);
end
