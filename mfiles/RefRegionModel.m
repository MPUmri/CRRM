function [ Ct ] = RefRegionModel( Crr, pkParams, t)
    kR = pkParams(1); % relative kTrans
    kepRR = pkParams(2); % kep of reference tissue
    kep = pkParams(3); % kep of TOI
    
    stepSize = t(2) - t(1);
    sT = length(t);
    Ct = zeros(sT,1);
    
    for k=1:sT
        Ct(k) = trapz(Crr(1:k).*exp(kep*(t(1:k)-t(k))));
    end
    
    Ct = kR * (Crr + stepSize*(kepRR-kep)*Ct);
end

