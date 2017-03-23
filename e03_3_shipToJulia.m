% This script prepares the QIN data so that Julia can process it

% Estimated run time: ~25 seconds

clearvars
addpath('./mfiles')

%% Set the input/output folders
dirLocation = DefaultFolders();
dataDir = fullfile(dirLocation.qin,'Stripped');
outDir = fullfile(dirLocation.qin,'ForJulia');

if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% Loop through the .mat files
matFiles = dir([dataDir '/*.mat']);

tic
for q=1:length(matFiles)
    curName = matFiles(q).name
    curFile = fullfile(dataDir,curName);
    load(curFile);
    %% Set the DCE time points
    t=1:dataInfo.timeResolution/dataInfo.timeResolutionFactor:dataInfo.timePts*dataInfo.timeResolution;
    t=t/60; % Convert to minutes
    % Initially, t and Cp are at a higher temporal resolution than the DCE data
    % We'll note down the time resolution factor, and account for it in Julia
    tResFactor = dataInfo.timeResolutionFactor;
    Cp = aif;
    %% Reshape DCE data
    [nX nY nZ nT] = size(dynData);
    mask = reshape(mask,[nX*nY nZ]);
    S = reshape(dynData,[nX*nY*nZ nT]);
    S(mask<1,:) = [];
    clearvars dynData
    %% Convert DCE Data from Signal to Concentration
    Ct = Signal2Conc4QIN(S',dataInfo.R10,dataInfo.TR,dataInfo.Alpha,dataInfo.r1p);
    %% Output to mat file
    outFile = fullfile(outDir,curName);
    save(outFile, 'Ct', 'Cp', 't', 'mask','tResFactor');
end
toc
disp('Done');
