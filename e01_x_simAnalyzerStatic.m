% Simulation Analyzer using static kepRef estimates.
% Note: This script is UNNECESSARY. Comments in code will be minimal.
% The code here was not used in the paper, but (I think) it is interesting nonetheless.
% This script uses the CLRRM with a range of kepRef values in order to see
% the effect they have on the estimates.
% In other words, this script tries to answer the question:
% "How bad will the CLRRM estimates be if the kepRef was over/underestimated by +/-50%?"

% Run time: ~470 seconds

%% Setup
clearvars
fclose('all')
addpath('./mfiles');

refName = 'refY';
% Choices are: 'refY', 'refWS', 'refP'

% The range of kepRef estimates
% If 'refY' is used, then this range goes from -50% to +50% from the true kepRef (which is 1.0 /min)
refList = {0.5,0.6,0.7,0.8,0.85,0.9,0.925,0.95,0.975,1.0,1.025,1.05,1.075,1.1,1.15,1.2,1.3,1.4,1.5};

%% Simulation Properties
simProp = SimProperties(refName);
CNR = simProp.CNR; % Range of Contrast-Noise Ratio to simulate
nVox = simProp.nVox; % Number of replications for each CNR
TRes = simProp.TRes; % Temporal resolutions (in seconds)
tDuration = simProp.tDuration; % Duration of DCE Acquisition (in minutes)

% Pharmacokinetics of tumour tissue
Kt = simProp.Kt;
ve = simProp.ve;
kep = simProp.kep;

% Pharmacokinetics of reference tissue
KtRR = simProp.KtRR;
veRR = simProp.veRR;
kepRR = simProp.kepRR;

relKt = Kt/KtRR;
relVe = ve/veRR;

%% Folder setup
dirLocation = DefaultFolders();
dataDir = fullfile(dirLocation.sim,simProp.name,'rawData');
outDir = fullfile(dirLocation.sim,simProp.name);
outFile = fullfile(dirLocation.results,['e01-simResultsStatic-' simProp.name '.csv']);

if ~exist(dirLocation.results)
    mkdir(dirLocation.results)
end

%% Initialize CSV with header
hdr=['KepRef,AIF,CNR,TemporalRes,TrueRelKt,TrueRelVe,TrueKep,errKt,errVe,errKep,stdErrKt,stdErrVe,stdErrKep, meanResid, stdResid'];
outID = fopen(outFile, 'w+');
fprintf(outID, '%s\n', hdr);

%% Process the simulated data
tic

curAif = 'Default';

for indCNR = 1:length(CNR)
    inFile = ['Sim-CNR-' num2str(CNR(indCNR)) '.mat'];
    inFile = fullfile(dataDir,inFile);
    load(inFile); % Load noisy simulation data
    for indTRes = 1:length(TRes)
        curTRes = TRes(indTRes);
        % Downsample the high temporal resolution (1s) data to desired resolution
        t = downsample(T, curTRes);
        Ct = downsample(CtNoisy, curTRes);
        Crr = downsample(CrrClean, curTRes);
        for i=1:length(refList)
            curRef = refList{i};
            [pkParamsCL, residCL, meanRefKep, stdRefKep, p] = CLRRM(Ct, Crr, t, curRef);
            ktError = 100*(pkParamsCL(:,1)-relKt)/relKt;
            veError = 100*(pkParamsCL(:,2)-relVe)/relVe;
            kepError = 100*(pkParamsCL(:,3)-kep)/kep;
            meanKtErr = mean(ktError);
            meanVeErr = mean(veError);
            meanKepErr = mean(kepError);
            stdKtErr = std(ktError);
            stdVeErr = std(veError);
            stdKepErr = std(kepError);
            meanResid = mean(residCL(:));
            stdResid = std(residCL(:));
            outLine = {curRef,curAif,curCNR,curTRes,relKt,relVe,kep,meanKtErr,meanVeErr,...
                meanKepErr,stdKtErr,stdVeErr,stdKepErr,meanResid,stdResid};
            fprintf(outID,'%f,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
        end
    end
end
toc
fclose(outID);
fclose('all')
