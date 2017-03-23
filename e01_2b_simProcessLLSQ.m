% Simulation Analyzer. This script:
% - reads in data made by Simulation Maker (e01_1_simMaker.m)
% - fits data using LRRM and variants of CRRM
% - saves results as .mat

% Total Runtimes:
% ~30 seconds for LRRM
% ~55 seconds for CLRRM
% ~30 seconds for CHRRM
% ~55 seconds for CLRRMm
% ~30 seconds for CLRRMt
% The actual run time will depend on which method(s) are entered in the
% variable 'methodList'
% If all variants are used, then total run time ~185 seconds

%% General setup
clearvars
fclose('all')
addpath('./mfiles');

% Match this with the choice in e01_1_simMaker.m, lest the code will get angry
refName = 'refY';
% Choices are: 'refY', 'refWS', 'refP'
% Refer to './mfiles/SimProperties.m for further details

% This controls which models/approaches will be attemped
methodList = {'LRRM','CLRRM','CHRRM','CLRRMm','CLRRMt'};
% Choices can include: 'LRRM', 'CLRRM', 'CLRRMn', 'CLRRMm', 'CLRRMt'
% Description:
% LRRM : Linear Reference Region Model (LRRM)
%       [Cardenas-Rodriguez et al. (2013) MRI, 31(4), 497-507.doi:10.1016/j.mri.2012.10.008]
% CLRRM : Contrained LRRM using median to estimate kepRR (paper version)
% CLRRMm : Constrained LRRM using mean over all voxels to estimate kepRR (crude method)
% CHRRM : Constrained Hybrid RRM using interquartile mean from non-linear estimates of kepRR
%         Requires e01_2a to be succesfully run before running this script.
% CLRRMt : Contrained LRRM using true value of kepRR (assuming kepRef was known a-priori)

% End of setup

%% Simulation Properties

simProp = SimProperties(refName);

% Copy variables to make code less verbose
listCNR = simProp.CNR; % Range of Contrast-Noise Ratio to simulate
nVox = simProp.nVox; % Number of replications for each CNR
listTRes = simProp.TRes; % Temporal resolutions (in seconds)
tDuration = simProp.tDuration; % Duration of DCE Acquisition (in minutes)

% Pharmacokinetics of reference tissue
KtRR = simProp.KtRR;
veRR = simProp.veRR;
kepRR = simProp.kepRR;

% Pharmacokinetics of tumour tissue
Kt = simProp.Kt;
ve = simProp.ve;
kep = simProp.kep;

% Relative parameters for tumour tissue
relKt = Kt/KtRR;
relVe = ve/veRR;

%% Initialize output folders

dirLocation = DefaultFolders();
dataDir = fullfile(dirLocation.sim,simProp.name,'rawData');
outDir = fullfile(dirLocation.sim,simProp.name);

% Check for inconsistencies
% Can't do CLRRMn if there are no pre-existing estimates from NRRM
if (strmatch('CLRRMn',methodList) & ~exist(fullfile(outDir,'NRRM')));
    error('The CLRRMn method needs estimates from non-linear fitting which should be in the NRRM folder, but no folder was found');
end

% Create output folders if needed
if (strmatch('LRRM',methodList) & ~exist(fullfile(outDir,'LRRM')));
    mkdir(fullfile(outDir,'LRRM'));
end
if (strmatch('CLRRM',methodList) & ~exist(fullfile(outDir,'CLRRM')));
    mkdir(fullfile(outDir,'CLRRM'));
end
if (strmatch('CHRRM',methodList) & ~exist(fullfile(outDir,'CHRRM')));
    mkdir(fullfile(outDir,'CHRRM'));
end
if (strmatch('CLRRMm',methodList) & ~exist(fullfile(outDir,'CLRRMm')));
    mkdir(fullfile(outDir,'CLRRMm'));
end
if (strmatch('CLRRMt',methodList) & ~exist(fullfile(outDir,'CLRRMt')));
    mkdir(fullfile(outDir,'CLRRMt'));
end

%% Analyze the simulation data
globalTimer = tic;

% Cycle through all concentration-noise ratios
for indCNR = 1:length(listCNR)
    inFile = ['Sim-CNR-' num2str(listCNR(indCNR)) '.mat'];
    inFile = fullfile(dataDir,inFile);
    load(inFile); % Load noisy simulation data
    % Cycle through all temporal resolutions
    for indTRes = 1:length(listTRes)
        curTRes = listTRes(indTRes);
        % Downsample the high temporal resolution (1s) data to desired resolution
        t = downsample(T, curTRes);
        Ct = downsample(CtNoisy, curTRes);
        Crr = downsample(CrrClean, curTRes);
        % Cycle through the method list and use each algorithm
        for i=1:length(methodList)
            curM = methodList{i};
            curFile = ['Sim-CNR-' num2str(curCNR) '-TRes-' num2str(curTRes) '.mat']; % Output file name
            if strcmp(curM,'LRRM')
                % The Linear Reference Region Model
                tic;
                [pkParamsLL, residLL] = LRRM(Ct, Crr, t, 0);
                runtimeLL = toc;
                % Save the raw results as a .mat file
                outFile = fullfile(outDir, 'LRRM', curFile);
                save(outFile, 'pkParamsLL', 'residLL','runtimeLL');
            elseif strcmp(curM,'CLRRM')
                % The Contrained Linear Reference Region Model
                % - using median from LRRM for kepRef estimate [paper's method]
                tic;
                [pkParamsCL, residCL, estKepRRCL, stdKepRRCL, pkParamsLLRaw] = CLRRM(Ct, Crr, t);
                runtimeCL=toc;
                outFile = fullfile(outDir, 'CLRRM', curFile);
                save(outFile, 'pkParamsCL', 'residCL','runtimeCL', 'estKepRRCL', 'stdKepRRCL', 'pkParamsLLRaw');
            elseif strcmp(curM,'CLRRMm')
                % The Contrained Linear Reference Region Model
                % - using mean from LRRM for all voxels for kepRef estimate
                [pkParamsCLm, residCLm, estKepRRCLm, stdKepRRCLm] = CLRRM(Ct, Crr, t, -1);
                outFile = fullfile(outDir, 'CLRRMm', curFile);
                save(outFile, 'pkParamsCLm', 'residCLm', 'estKepRRCLm', 'stdKepRRCLm');
            elseif strcmp(curM,'CHRRM')
                % The Contrained Linear Reference Region Model
                % - using estimates from the NRRM to get kepRR
                nrrmFile = fullfile(outDir,'NRRM',curFile);
                load(nrrmFile); % This file contains 'medianKepRR' from NRRM
                tic;
                [pkParamsCH, residCH] = CLRRM(Ct, Crr, t, medianKepRR);
                runtimeCH = toc;
                estKepRRCH = medianKepRR;
                outFile = fullfile(outDir, 'CHRRM', curFile);
                save(outFile, 'pkParamsCH', 'residCH', 'runtimeCH', 'estKepRRCH');
            elseif strcmp(curM,'CLRRMt')
                % The Contrained Linear Reference Region Model
                % - when the true kepRef value is known a-priori
                tic;
                [pkParamsCLt, residCLt] = CLRRM(Ct, Crr, t, kepRR);
                runtimeCLt = toc;
                outFile = fullfile(outDir, 'CLRRMt', curFile);
                save(outFile, 'pkParamsCLt', 'residCLt', 'runtimeCLt');
            end
        end
    end
end
toc(globalTimer)
disp('.')
disp('.')
disp('.')
disp('Done processing simulation data')
