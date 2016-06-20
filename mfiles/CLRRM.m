function [pkParams, resid, estKepRR, stdKepRR, p] = CLRRM(Ct, Crr, t, kepRef)
  % Constrained Linear Reference Region fit
  % Inputs:
  % Ct - TxN Matrix where T is time and N is number of voxels
  %    - Concentration-time curve in tissue of interest
  % Crr - Tx1
  %     - Concentration-time curve in reference tissue
  % t  - Tx1
  %    - time, in minutes
  % kepRef - numeric value
  % If kepRef > 0, then this value will be used as kepRef
  %           = 0, then mean will be used as kepRef
  %           = -1, then interquartile mean will be used as kepRef
  %
  % Outputs:
  % pkParams - Nx3 matrix, [Ktrans/Ktrans_RR  v_e/v_e,rr  kep]
  % resid - residuals of the fit
  % estKepRR - estimated kepRR value
  % stdKepRR - std.Deviation of estimated kepRR value
  % p - pkParams from the initial LRRM fit

    % Use interquartile mean by default
    if nargin<4
        kepRef = -1;
    end
%%
    if kepRef <=0
        % Do first run of LRRM
        p = LRRM(Ct,Crr,t,1);
        if kepRef == 0
            % Using the mean of positive values
            p1 = p(:,1); % First fitted parameter
            p2 = p(:,2); % Second fitted parameter
            goodVals = (p1>0) & (p2>0) & p(:,3)>0; % Good values are when all parameters are positive
            kepRRs = p2(goodVals)./p1(goodVals); % Get kepRR, which is ratio of p2/p1
            estKepRR = mean(kepRRs);
            stdKepRR = std(kepRRs);
        elseif kepRef == -1
            % Using the interquartile mean
            p1 = p(:,1);
            p2 = p(:,2);
            goodVals = (p1>0) & (p2>0) & p(:,3)>0;
            kepRRs = p2(goodVals)./p1(goodVals);
            qtRange = quantile(kepRRs,[.25 .75]); % Get the 25th and 75th quantile
            xMask = kepRRs>qtRange(1) & kepRRs<qtRange(2);
            estKepRR = mean(kepRRs(xMask)); % Take mean of values within the 25th and 75th quantiles
            stdKepRR = std(kepRRs(xMask));
        end
    else
        % Otherwise, use the kepRef that is provided with function call
        estKepRR = kepRef;
        stdKepRR = nan;
        p = nan;
    end

    stepSize = t(2)-t(1); % Assuming constant step size throughout acquisition

    % Initialize matrices
    [sT, sX] = size(Ct);

    M1 = zeros(sT,1);
    M2 = M1;

    M1 = Crr;
    M2 = stepSize*cumtrapz(Crr);

    resid = zeros(sX,1);
    pkParams = zeros(sX,2);

    % Solve the linear problem
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

    pkParams(:,3)=pkParams(:,2);
    pkParams(:,2)=estKepRR*pkParams(:,1);

    % pkParams = [kTrans_TOI/kTrans_RR, kTrans_TOI/ve_RR, kTrans_TOI/ve_TOI]
    % desire: pkParams = [kTrans_TOI/ktrans_RR, ve_TOI/ve_RR, kTrans_TOI/ve_TOI]
    pkParams(:,2) = pkParams(:,2)./pkParams(:,3);
end
