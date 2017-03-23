% Process QIBA Phantom Data using LRRM and CLRRM
% Uses two combinations of reference region parameters:
% - 'refY' : Reference KTrans = 0.1 /min, Reference EES = 0.1
% - 'refWS': Reference KTrans = 0.07 /min, Reference EES = 0.14

% The input data is 30 voxels with 1321 time points.
% Run times are displayed in console. 

% Reminder:
% LRRM is the original model
% CLRRM is the modified model, using median for estimating kepRef

% Estimated total run time: < 1 sec

%% Initialize

clearvars 
addpath('./mfiles')

dirLocation = DefaultFolders();
qibaFile = fullfile(dirLocation.qiba,'QIBAv6-Mini.mat');
outDir = fullfile(dirLocation.qiba,'Results');

if ~exist(outDir)
    mkdir(outDir)
end

totalTimer = tic;

%% Load phantom data
load(qibaFile);

[nX nY nT] = size(concData);
ctData = reshape(concData,[nX*nY nT]);
ctData = ctData';

%% Process data using 'refY' parameters

simProp = SimProperties('refY');
refKt = simProp.KtRR;
refVe = simProp.veRR;
refKep = refKt/refVe;
Crr = TrapzKety(aif,[refKt, refKep],T,0);

disp('Using refY...')
tic
[pkParamsLL, residLL] = LRRM(ctData,Crr,T,0);
runTime = toc;
disp(['LRRM took ' num2str(runTime) ' seconds'])
tic
[pkParamsCL, residCL, refKepEst, refKepStd, pLL] = CLRRM(ctData,Crr,T);
runTime = toc;
disp(['CLRRM took ' num2str(runTime) ' seconds'])
tic

outFile = fullfile(outDir,'QIBAv6-refY-LL.mat');
save(outFile,'pkParamsLL','residLL',...
    'pkParamsCL','residCL','refKepEst', 'refKepStd', 'pLL');
%% Process data using 'refWS' parameters

simProp = SimProperties('refWS');
refKt = simProp.KtRR;
refVe = simProp.veRR;
refKep = refKt/refVe;
Crr = TrapzKety(aif,[refKt, refKep],T,0);

disp('Using refWS...')
tic
[pkParamsLL, residLL] = LRRM(ctData,Crr,T,0);
runTime = toc;
disp(['LRRM took ' num2str(runTime) ' seconds'])
tic
[pkParamsCL, residCL, refKepEst, refKepStd, pLL] = CLRRM(ctData,Crr,T);
runTime = toc;
disp(['CLRRM took ' num2str(runTime) ' seconds'])
outFile = fullfile(outDir,'QIBAv6-refWS-LL.mat');
save(outFile,'pkParamsLL','residLL',...
    'pkParamsCL','residCL','refKepEst', 'refKepStd', 'pLL');
%%
disp('.')
disp('..')
disp(['Total run time: ' num2str(toc(totalTimer))]);