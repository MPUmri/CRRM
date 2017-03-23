% Process QIN Breast DCE-MRI data using LRRM and CLRRM and variants.

% Estimated run time: ~45 seconds


%% General Setup
clearvars
addpath('./mfiles')

% Pick the choice of parameters for reference tissue
refName = 'refY'; % Choices: 'refY', 'refWS', 'refP'
% Make sure this matches the choice in e03_4_shipToJulia.jl
% Refer to './mfiles/SimProperties.m for further details

% Save output as .mat file
doSave = 1; % Should be 1 or future steps won't work

% This controls which kind of algorithms will be attemped
methodList = {'LRRM','CLRRM','CHRRM','LLTM'};
% Description:
% LRRM : Linear Reference Region Model (LRRM)
%       [Cardenas-Rodriguez et al. (2013) MRI, 31(4), 497-507.doi:10.1016/j.mri.2012.10.008]
% CLRRM : Contrained LRRM using interquartile mean from LRRM to estimate kepRef (paper version)
% CHRRM : Constrained LRRM using NRRM as first fit to estimate kepRef
% LLTM : Linear Tofts Model

%%%% End of Setup

%% Initialize input and output folders
dirLocation = DefaultFolders();
dataDir = fullfile(dirLocation.qin,'Stripped');
outDir = fullfile(dirLocation.qin,['MapsLL-' refName]);

if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% Initialize properties for reference tissue
refProp = SimProperties(refName);
ktRR = refProp.KtRR;
veRR = refProp.veRR;
kepRR = refProp.kepRR;


%% Process patient data
matFiles = dir([dataDir '/*.mat']);
mainTimer = tic;
for q=1:length(matFiles) % Cycle through all files
    curName = matFiles(q).name
    curFile = fullfile(dataDir,curName);
    load(curFile);
    %% Obtain frame times (upsampled)
    t=1:dataInfo.timeResolution/dataInfo.timeResolutionFactor:dataInfo.timePts*dataInfo.timeResolution;
    t=t/60; % Convert to minutes
    %% Simulate the reference tissue curve
    Crr = TrapzKety(aif,[ktRR, kepRR],t);
    Crr = downsample(Crr, dataInfo.timeResolutionFactor);
    Cp = downsample(aif,dataInfo.timeResolutionFactor);
    %% Mask and Reshape the DCE data to make future steps easier
    [nX nY nZ nT] = size(dynData);
    mask = reshape(mask,[nX*nY nZ]);
    S = reshape(dynData,[nX*nY*nZ nT]);
    S(mask<1,:) = [];
    clearvars dynData
    %% Convert DCE data to Concentration data
    Ct = Signal2Conc4QIN(S',dataInfo.R10,dataInfo.TR,dataInfo.Alpha,dataInfo.r1p);
    %% Initialize variables
    % Pharmacokinetic Parameters
    pkParamsLL = 0;
    pkParamsCL = 0;
    pkParamsCH = 0;
    pkParamsLT = 0;
    % Residual of fit
    residLL = 0;
    residCL = 0;
    residCH = 0;
    residLT = 0;
    % Runtimes
    runtimeLL = 0;
    runtimeCL = 0;
    runtimeCH = 0;
    runtimeLT = 0;
    % Estimated kep for reference region
    meanRefKepCL = 0;
    % StdDev of reference region kep estimate
    stdRefKepCL = 0;
    % Raw fit of LRRM that is used by CLRRM
    p=0;
    %% Fit the data using selected Algorithm/Model
    % We have to downsample 't' so it matches everything else
    t=downsample(t,dataInfo.timeResolutionFactor);
    % Cycle through all the selected algorithm/models
    for i=1:length(methodList)
        curM = methodList{i};
        if strcmp('LRRM',curM)
            % Linear Reference Region Model
            tic;
            [pkParamsLL, residLL] = LRRM(Ct,Crr,t,0);
            runtimeLL=toc;
        elseif strcmp('CLRRM',curM)
            % Constrained Linear Reference Region Model...
            % ...using interquartile mean as the kepRR estimate
            tic;
            [pkParamsCL, residCL,estKepRRL, stdKepRRL,p] = ...
                CLRRM(Ct,Crr,t);
            runtimeCL=toc;
        elseif strcmp('CHRRM',curM)
            % Constrained Hybrid Reference Region Model...
            % ...using the NRRM for kepRR estimate
            nrrmFile = fullfile(dirLocation.qin,'FromJulia',curName);
            load(nrrmFile,'estKepRRN');
            tic;
            [pkParamsCH, residCH] = CLRRM(Ct,Crr,t,estKepRRN);
            runtimeCH=toc;
        elseif strcmp('LLTM',curM)
            % Linear Tofts Model
            tic;
            [pkParamsLT, residLT] = LLSQ(Ct,Cp,t,0);
            runtimeLT=toc;
        end
    end
    nVox = sum(mask(:));
    % Save output
    if doSave == 1
        outFile = fullfile(outDir,curName);
        save(outFile, 'pkParamsLL','pkParamsCL','pkParamsCH','pkParamsLT',...
            'residLL','residCL','residCH','residLT','p','mask',...
            'estKepRRL','stdKepRRL','methodList','nVox',...
            'runtimeLL','runtimeCL','runtimeCH','runtimeLT');
    end

end
toc(mainTimer)
disp('Done');
