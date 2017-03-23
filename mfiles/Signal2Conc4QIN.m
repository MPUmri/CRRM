function [ Ct ] = Signal2Conc4QIN(S,R10,TR,alpha,r1)
% This converts MRI signal to concentration (made for QIN Breast data)
% Inputs:
%   - S: T1-weighted signal (Time-by-numVox)
%   - R10: R1 at pre-contrast frame
%   - TR: Repetition time, in seconds
%   - alpha: Flip angle, in degrees
%   - r: relaxitivity constant for contrast agent
%% Debug:
%     S = S';
%     R10 = dataInfo.R10;
%     TR = dataInfo.TR;
%     alpha = dataInfo.Alpha;
%     r1 = dataInfo.r1p; 

%%
    [nT nV] = size(S);
    cosAlpha = cos(alpha * pi / 180);
    TR = TR * 1E-3;
    S0 = mean(S(1:2,:),1); % Pre-contrast signal (NOT S0 from SPGR equation)
    E0 = exp(-R10*TR);
    %%
    s = S ./ repmat(S0,[nT 1]); % Signal normalized by pre-contrast signal
    %%
    E = (1 - s + s * E0 - E0 * cosAlpha) ./ ...
        (1 - s*cosAlpha + s*E0*cosAlpha - E0*cosAlpha);
    %%
    r1Data = zeros(size(E));
    Ct = r1Data;
    
    r1Data(E>0) = (-1/TR) * log(E(E>0));
    Ct(r1Data>0) = (r1Data(r1Data>0) - R10)/r1;
end

