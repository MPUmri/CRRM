function [pkParams, resid, exitFlag] = NLSQ(Ct, Cp, t, modType)
  % Non-Linear Fitting

  if nargin<4
    modType = 0;
  end

  % Initialize matrices
  datSize=size(Ct);
  [sT sX] = size(Ct);

  if modType == 0
    % Tofts Model
    p = zeros(sX,2);
    resid = zeros(sX,1);
    exitFlag = resid;

    options=optimset('Algorithm','levenberg-marquardt','display','off');
    pkGuess = [0.01,0.01];
    parfor i=1:sX
        [p(i,:), resid(i), grbg, exitFlag(i)] = lsqnonlin(@(x) Ct(:,i) - TrapzKety(Cp,x,t),pkGuess, [], [],options);
        %[p(i,:)] = lsqnonlin(@(x) Ct(:,i) - TrapzKetyK(Cp,x,t)',pkGuess, [], [],options);
    end
    pkParams = p;
  end
end
