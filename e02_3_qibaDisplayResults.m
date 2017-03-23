% This script calculates the percent errors in the QIBA phantom data and
% displays the results

% Estimated run time: < 1 s

%% General setup
clearvars
addpath('./mfiles')

% Choice of reference region parameters
refName = 'refY';
% Choices are: 'refY', 'refWS'

% Do we have results for NRRM (Nonlinear Reference Region Model)
% This requires running the julia script.
% Set to 0 to only compare linear models only.
doNL = 1;

%%
% Set doNL to false if refName is not refY
doNL = doNL && strcmp(refName,'refY');

%% Load data
dirLocation = DefaultFolders();
dataDir = fullfile(dirLocation.qiba,'Results');
outFile = fullfile(dirLocation.results,['e02-results-' refName '.mat']);

fileLL = fullfile(dataDir,['QIBAv6-' refName '-LL.mat']);
fileNL = fullfile(dataDir,['QIBAv6-' refName '-NL.mat']);
load(fileLL);
if doNL
    load(fileNL);
end

simProp = SimProperties(refName);
refKt = simProp.KtRR;
refVe = simProp.veRR;

%% Build true map
nX = 60;
nY = 50;

valKt = [0.01,0.02,0.05,0.1,0.2,0.35];
valVe = [0.01,0.05,0.1,0.2,0.5];
trueKt = zeros(nX,nY);
trueVe = trueKt;

for i=1:nX/10
    trueKt(1+(i-1)*10:i*10,:) = valKt(i);
end
for j=1:nY/10
    trueVe(:,1+(j-1)*10:j*10) = valVe(j);
end

% Subsample
trueKt = trueKt(5:10:(end-5), 5:10:(end-5),:);
trueVe = trueVe(5:10:(end-5), 5:10:(end-5),:);
trueKep = trueKt./trueVe;
nX = 6;
nY = 5;

%% Identify variables and convert relative estimates to absolute estimates
ktLL = refKt * pkParamsLL(:,1);
ktCL = refKt * pkParamsCL(:,1);

veLL = refVe * pkParamsLL(:,2);
veCL = refVe * pkParamsCL(:,2);

kepLL = pkParamsLL(:,3);
kepCL = pkParamsCL(:,3);

%% Reshape into 2D - matching the original phantom 
ktLL = reshape(ktLL,[nX nY]);
ktCL = reshape(ktCL,[nX nY]);

veLL = reshape(veLL,[nX nY]);
veCL = reshape(veCL,[nX nY]);

kepLL = reshape(kepLL,[nX nY]);
kepCL = reshape(kepCL,[nX nY]);

residLL = reshape(residLL,[nX nY]);
residCL = reshape(residCL,[nX nY]);

%% Calculate percent errors
errKtLL = 100*(ktLL-trueKt)./trueKt;
errKtCL = 100*(ktCL-trueKt)./trueKt;

errVeLL = 100*(veLL-trueVe)./trueVe;
errVeCL = 100*(veCL-trueVe)./trueVe;

errKepLL = 100*(kepLL-trueKep)./trueKep;
errKepCL = 100*(kepCL-trueKep)./trueKep;

%% Do all the above for NonLinear results too, if needed
if doNL
    pkParamsN = pkParamsN';
    residNL = residN';
    ktNL = refKt * pkParamsN(:,1);
    veNL = refVe * pkParamsN(:,2);
    kepNL = pkParamsN(:,3);
    
    ktNL = reshape(ktNL,[nX nY]);
    veNL = reshape(veNL,[nX nY]);
    kepNL = reshape(kepNL,[nX nY]);
    residNL = reshape(residNL,[nX nY]);
    
    errKtNL = 100*(ktNL-trueKt)./trueKt;
    errVeNL = 100*(veNL-trueVe)./trueVe;
    errKepNL = 100*(kepNL-trueKep)./trueKep;
    
    pkParamsCN = pkParamsCN';
    residCN = residCN';
    ktCN = refKt * pkParamsCN(:,1);
    veCN = refVe * pkParamsCN(:,2);
    kepCN = pkParamsCN(:,3);
    
    ktCN = reshape(ktCN,[nX nY]);
    veCN = reshape(veCN,[nX nY]);
    kepCN = reshape(kepCN,[nX nY]);
    residCN = reshape(residCN,[nX nY]);
    
    errKtCN = 100*(ktCN-trueKt)./trueKt;
    errVeCN = 100*(veCN-trueVe)./trueVe;
    errKepCN = 100*(kepCN-trueKep)./trueKep;
end
%% Display results
%% LRRM
figure
imagesc(abs(errKepLL))
set(gca,'YTick',[1:6], 'YTickLabel',valKt)
set(gca,'XTick',[1:5], 'XTickLabel',valVe)
xlabel('Ve')
ylabel('KTrans')
title('Absolute Percent Error in k_{ep} using LRRM')
colorbar
%% CLRRM using median-based kepRef estimate
figure
imagesc(abs(errKepCL))
set(gca,'YTick',[1:6], 'YTickLabel',valKt)
set(gca,'XTick',[1:5], 'XTickLabel',valVe)
xlabel('Ve')
ylabel('KTrans')
title('Absolute Percent Error in k_{ep} using CLRRM')
colorbar
%% Nonlinear Reference Region Model
if doNL == 1
    figure
    imagesc(abs(errKepNL))
    set(gca,'YTick',[1:6], 'YTickLabel',valKt)
    set(gca,'XTick',[1:5], 'XTickLabel',valVe)
    xlabel('Ve')
    ylabel('KTrans')
    title('Absolute Percent Error in k_{ep} using NRRM')
    colorbar
    
    figure
    imagesc(abs(errKepCN))
    set(gca,'YTick',[1:6], 'YTickLabel',valKt)
    set(gca,'XTick',[1:5], 'XTickLabel',valVe)
    xlabel('Ve')
    ylabel('KTrans')
    title('Absolute Percent Error in k_{ep} using CNRRM')
    colorbar
end
%%
