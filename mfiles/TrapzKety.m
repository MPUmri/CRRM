function ct = TrapzKety(theAif, pkParams, t, doExt)
    % Tofts-Kety model using trapz() for integration
    % Inputs:
    % - theAif: [vector] The Arterial Input Function
    % - pkParams: [vector] Pharmacokinetic parameters, in format [kTrans,kep,vp]
    % - t: [vector] Time (in minutes) for each timepoint
    % - doExt: [boolean] Set to 1 for the Extended Tofts-Kety Model

    % Simple Tofts-Kety Model is used by default
		if nargin < 4
			doExt = 0;
    end

    % Obtain relevant parameters from input
    kTrans = pkParams(1);
    kep = pkParams(2);

    stepSize = t(2)-t(1);

    % Initialize vector
    ct = zeros(size(theAif));

    % Calculate concentration at each time point
    for k=1:length(t)
        ct(k) = kTrans*trapz(theAif(1:k).*exp(kep*(t(1:k)-t(k))));
    end

    % Scale by temporal resolution, and transpose (transpose is personal preference)
    ct = stepSize * ct';

    % Include the contribution from the plasma volume if Extended model is used
    if doExt == 1
    	ct = ct + pkParams(3)*theAif';
    end

end
