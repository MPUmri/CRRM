function [simProp] = SimProperties(refName)
    % General properties for the simulation
    % Inputs:
    % refName: Choice of 'refY', 'refP' and 'refWS'.
    %          Default: 'refY'
    %          The choices determine the tracer-kinetic parameter of the simulated reference tissue.
    %          The paper uses 'refY'.
    %          'refY' : Reference Tissue KTrans = 0.1 /min, Reference Tissue EES = 0.1
    %          'refWS': Reference Tissue KTrans = 0.07 /min, Reference Tissue EES = 0.14
    %          'refP' : Reference Tissue KTrans = 0.137 /min, Reference Tissue EES = 0.115
    %          References:
    %           - Yankeelov et al. (2005), MRI, 23(4), 519-29. doi:10.1016/j.mri.2005.02.013
    %           - Walker-Samuel et al. (2007), PMB, 52(1), 75-89. doi:10.1088/0031-9155/52/1/006
    %           - Padhani et al. (2002), NMR in Biomedicine, 15(2), 143-153. doi:10.1002/nbm.732


    if nargin < 1
        % Use 'refY' if no input is provided
        refName = 'refY';
    end

    % Pharmacokinetics of reference tissue
    % Units are min^{-1} for KTrans & kep, while EES is unitless
    if strcmp(refName,'refP')
        % Values from Tables 2 and 3 (Left Obturator)
        % from Padhani et al. (2002), NMR in Biomedicine, 15(2), 143-153. doi:10.1002/nbm.732
        simProp.name = refName;
        simProp.KtRR = 0.137;
        simProp.veRR = 0.115;
        simProp.kepRR = simProp.KtRR/simProp.veRR;
    elseif strcmp(refName,'refWS')
        % Walker-Samuel et al. (2007), PMB, 52(1), 75-89. doi:10.1088/0031-9155/52/1/006
        simProp.name = refName;
        simProp.KtRR = 0.07;
        simProp.veRR = 0.14;
        simProp.kepRR = simProp.KtRR/simProp.veRR;
    else
        % Values from Yankeelov et al. (2005), MRI, 23(4), 519-29. doi:10.1016/j.mri.2005.02.013
        simProp.name = 'refY';
        simProp.KtRR = 0.1;
        simProp.veRR = 0.1;
        simProp.kepRR = simProp.KtRR/simProp.veRR;
    end

    % Pharmacokinetics of tumour tissue
    simProp.Kt = 0.25; % units: /min
    simProp.ve = 0.4; % unitless (fractional)
    simProp.kep = simProp.Kt/simProp.ve; % units: /min

    % Properties of the simulated data
    simProp.CNR = [5:5:50]; % Range of Contrast-Noise Ratios to simulate
    simProp.nVox = 10000; % Number of replications for each CNR
    simProp.TRes = [1,5,10,15,30,60]; % Temporal resolutions (in seconds) - Integers only. Can't have any value below 1.
    simProp.tDuration = 10; % Duration of DCE Acquisition (in minutes)

end
