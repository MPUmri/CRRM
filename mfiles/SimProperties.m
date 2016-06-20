function [simProp] = SimProperties(refName)
    % General properties for the simulation
    
    if nargin < 1
        % Use 'refY' if no input is provided
        refName = 'refY';
    end

    % Pharmacokinetics of reference tissue
    % Units are min^{-1} for KTrans, while EES is unitless
    if strcmp(refName,'refP')
        % Values from Tables 2 and 3 (Left Obturator)
        % from Padhani et al. (2002), NMR in Biomedicine, 15(2), 143–153. doi:10.1002/nbm.732
        simProp.name = refName;
        simProp.KtRR = 0.137;
        simProp.veRR = 0.115;
        simProp.kepRR = simProp.KtRR/simProp.veRR;
    elseif strcmp(refName,'refWS')
        % Walker-Samuel et al. (2007), PMB, 52(1), 75–89. doi:10.1088/0031-9155/52/1/006
        simProp.name = refName;
        simProp.KtRR = 0.07;
        simProp.veRR = 0.14;
        simProp.kepRR = simProp.KtRR/simProp.veRR;
    else
        % Values from Yankeelov et al. (2005), MRI, 23(4), 519–29. doi:10.1016/j.mri.2005.02.013
        simProp.name = 'refY';
        simProp.KtRR = 0.1;
        simProp.veRR = 0.1;
        simProp.kepRR = simProp.KtRR/simProp.veRR;
    end
    
    simProp.CNR = [5:5:50]; % Range of Contrast-Noise Ratio to simulate
    simProp.nVox = 10000; % Number of replications for each CNR
    simProp.TRes = [1,5,10,15,30,60]; % Temporal resolutions (in seconds)
    simProp.tDuration = 10; % Duration of DCE Acquisition (in minutes)

    % Pharmacokinetics of tumour tissue
    simProp.Kt = 0.25; % /min
    simProp.ve = 0.4; % unitless (fractional)
    simProp.kep = simProp.Kt/simProp.ve; % /min

end