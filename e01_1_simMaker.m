% Simulation Maker. This script:
% - simulates concentration-time curves at 0.5s temporal resolution [using Tofts Model]
% - adds gaussian noise at prescribed concentration-noise ratios
% - saves data as .mat file for further processing

% Run time: ~1 seconds

%% General setup

clearvars
addpath('./mfiles')

% Pick the choice of parameters for reference tissue
refName = 'refY';  % The paper uses 'refY'
% Choices are:
% 'refY' : Reference KTrans = 0.1 /min, Reference EES = 0.1
% 'refWS': Reference KTrans = 0.07 /min, Reference EES = 0.14
% 'refP' : Reference KTrans = 0.137 /min, Reference EES = 0.115
% References:
% - Yankeelov et al. (2005), MRI, 23(4), 519-29. doi:10.1016/j.mri.2005.02.013
% - Walker-Samuel et al. (2007), PMB, 52(1), 75-89. doi:10.1088/0031-9155/52/1/006
% - Padhani et al. (2002), NMR in Biomedicine, 15(2), 143-153. doi:10.1002/nbm.732

% For further settings, look in './mfiles/SimProperties.m'.

% Set a seed for reproducibility
rng(12345) % The paper arbitrarily uses 12345 as the seed

% End of Setup
disp('Simulation Maker started')
%% Obtain Simulation Properties

simProp = SimProperties(refName);

% Copying variables to make future code less verbose
listCNR = [5]; % List of Contrast-Noise Ratio to simulate
nVox = simProp.nVox; % Number of replications for each CNR
listTRes = [10]; % List Temporal resolutions (in seconds)
tDuration = simProp.tDuration; % Duration of DCE Acquisition (in minutes)

% Pharmacokinetics of tumour tissue
Kt = simProp.Kt;
ve = simProp.ve;
kep = simProp.kep;

% Pharmacokinetics of reference tissue
KtRR = simProp.KtRR;
veRR = simProp.veRR;
kepRR = simProp.kepRR;

%% Initialize output folder
dirLocation = DefaultFolders();
outDir = fullfile(dirLocation.sim,simProp.name,'rawData');

% Make the output folder if they don't already exist
if ~exist(dirLocation.sim)
    mkdir(dirLocation.sim)
end

if ~exist(outDir)
    mkdir(outDir)
end

%% Stimulate the simulation
tic

% Time points
T=0:0.5:60*simProp.tDuration; % Use a 0.5 second temporal resolution (initially)
T=T/60; % Convert to minutes

% Obtain AIF using model from Parker et al. (2006) MRM, 56(5), 993-1000. doi:10.1002/mrm.21066
Cp = ParkerAif(T,T(120)); % The injection occurs at 1 minute (frame 120 at 0.5s temporal resolution)

% Simulate the noiseless (clean) concentration-time curves at 0.5 temporal resolution
CtClean = TrapzKety(Cp,[Kt, kep],T,0);
CrrClean = TrapzKety(Cp,[KtRR, kepRR],T,0);
% Downsample to 1s temporal resolution (we don't need anything higher for simulation)
CtClean = downsample(CtClean, 2);
CrrClean = downsample(CrrClean, 2);
T = downsample(T,2);

% Add noise for each concentration-noise ratio (CNR) level and export
for indCNR = 1:length(listCNR)
    curCNR = listCNR(indCNR);
    curSigma = max(CtClean(:))/curCNR;
    % Noise is gaussian (can have negative concentration values)
    CtNoisy = repmat(CtClean,[1 nVox]) + curSigma * randn(length(CtClean),nVox);
    % Output noisy curves as .mat files
    curFile = ['Sim-CNR-' num2str(curCNR) '.mat'];
    outFile = fullfile(outDir, curFile);
    save(outFile, 'CtNoisy','CrrClean','T','curCNR');
end

disp('.')
disp('.')
disp('.')
disp ('Simulation Maker done')
toc
