function theAif = ParkerAif(t, t0)
% Generates an AIF using the model and parameters from:
% Parker et al. (2006) MRM, 56(5), 993-1000. doi:10.1002/mrm.21066
% Inputs:
%     t: [vector] time, in minutes, for each DCE frame
%     t0: [real value] time of injection, in minutes

    Hct = 0.42; % Common estimate for hematocritc from paper

    if nargin < 1
        t=0:499;
        t=t'/60;
    end
    if nargin < 2
        t0 = 0;
    end

    % Set all pre-injection time points to zero
    t = t-t0;
    t(t<0) = 0;

    % Population-averaged parameters from Parker et al. Paper
    A(1) = 0.809;
    A(2) = 0.330;
    T(1) = 0.17046;
    T(2) = 0.365;
    Sigma(1) = 0.0563;
    Sigma(2) = 0.132;
    Alpha = 1.050;
    Beta = 0.1685;
    s = 38.078;
    Tau = 0.483;


    C_b = (A(1)/(Sigma(1)*sqrt(2*pi))) * exp(-(t-T(1)).^2 / (2*(Sigma(1))^2)) ...
        + (A(2)/(Sigma(2)*sqrt(2*pi))) * exp(-(t-T(2)).^2 / (2*(Sigma(2))^2)) ...
        + Alpha * exp(-Beta*t) ./ (1+ exp(-s*(t-Tau)));

    % C_b is the concentration in the whole blood
    % Take hematocrit into account to get concentration in plasma (i.e. C_p)
    theAif = C_b / (1-Hct);
    theAif(t<=0)=0; % Ensure concentration is zero for all pre-injection time points
end
