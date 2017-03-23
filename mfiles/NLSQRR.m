function [pkParams, resid, exitFlag] = NLSQRR(Ct, Crr, t, doPure)
  % Non-Linear Fitting for Reference Region Model
%%
    if nargin<4
        doPure = 0;
    end
    % Initialize matrices
    [sT sX] = size(Ct);

    % Tofts Model
    p = zeros(sX,3);
    resid = zeros(sX,1);
    exitFlag = resid;

    options=optimset('Algorithm','levenberg-marquardt','display','off');
    pkGuess = [0.01,0.01,0.01];
    for i=1:sX
        [p(i,:), resid(i), grbg, exitFlag(i)] = lsqnonlin(@(x) Ct(:,i) - RefRegionModel(Crr,x,t),pkGuess, [], [],options);
        %[p(i,:)] = lsqnonlin(@(x) Ct(:,i) - TrapzKetyK(Cp,x,t)',pkGuess, [], [],options);
    end
    if doPure == 0
        pkParams(:,1) = p(:,1); % Relative KTrans
        pkParams(:,2) = p(:,2).*p(:,1)./p(:,3); % Relative ve
        pkParams(:,3) = p(:,3); % kep or (KTrans/ve)
    else
        pkParams = p;
    end

end
