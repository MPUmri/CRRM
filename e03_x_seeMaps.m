clearvars
addpath('./mfiles')

refName = 'refY'; % Choices: 'refY', 'refWS', 'refP'
% Depends on choice used in e03_5

% Patient of Interest - Which patient to process?
PoI = 1;

%% Get the reference region properties
refProp = SimProperties(refName);
ktRR = refProp.KtRR;
veRR = refProp.veRR;
kepRR = refProp.kepRR;

%% Load Data

dirLocation = DefaultFolders();
llDir = fullfile(dirLocation.qin,['MapsLL-' refName]);
matFiles = dir([llDir '/*.mat']);

q=PoI;
curFile = matFiles(q).name;
llFile = fullfile(llDir, curFile);

load(llFile)

%% Re-arrange Tofts Model map
pkParamsT = pkParamsLT;
pkParamsT(:,3) = pkParamsT(:,2);
pkParamsT(:,2) = pkParamsT(:,1)./pkParamsT(:,3);

%% Get KTrans
ktLL = ktRR*pkParamsLL(:,1);
ktCL = ktRR*pkParamsCL(:,1);
ktT = pkParamsT(:,1);
%% Get Ve
veLL = veRR*pkParamsLL(:,2);
veCL = veRR*pkParamsCL(:,2);
veT = pkParamsT(:,2);
%% Get Kep
kepLL = pkParamsLL(:,3);
kepCL = pkParamsCL(:,3);
kepT = pkParamsT(:,3);
%% Get Mask
nZ = size(mask,2);
mask = reshape(mask,[320 320 nZ]);

%% Make the maps

% Initialize matrices
imKtLL = zeros(size(mask));
imKtCL = imKtLL;
imKtT = imKtLL;
imKepLL = imKtLL;
imKepCL = imKtLL;
imKepT = imKtLL;
imVeLL = imKtLL;
imVeCL = imKtLL;
imVeT = imKtLL;

% KTrans
imKtLL(mask>0) = ktLL;
imKtCL(mask>0) = ktCL;
imKtT(mask>0) = ktT;

% kep
imKepLL(mask>0) = kepLL;
imKepCL(mask>0) = kepCL;
imKepT(mask>0) = kepT;

% ve
imVeLL(mask>0) = veLL;
imVeCL(mask>0) = veCL;
imVeT(mask>0) = veT;

% Reshape the maps into 3D
imKtLL = reshape(imKtLL,[320 320 nZ]);
imKtCL = reshape(imKtCL,[320 320 nZ]);
imKtT = reshape(imKtT,[320 320 nZ]);

imKepLL = reshape(imKepLL,[320 320 nZ]);
imKepCL = reshape(imKepCL,[320 320 nZ]);
imKepT = reshape(imKepT,[320 320 nZ]);

imVeLL = reshape(imVeLL,[320 320 nZ]);
imVeCL = reshape(imVeCL,[320 320 nZ]);
imVeT = reshape(imVeT,[320 320 nZ]);

% Crop the masks to focus on tumour region
mask = reshape(mask,[320 320 nZ]);
imKtLL = AutoCrop(imKtLL,mask);
imKtCL = AutoCrop(imKtCL,mask);
imKtT = AutoCrop(imKtT,mask);

imKepLL = AutoCrop(imKepLL,mask);
imKepCL = AutoCrop(imKepCL,mask);
imKepT = AutoCrop(imKepT,mask);

imVeLL = AutoCrop(imVeLL,mask);
imVeCL = AutoCrop(imVeCL,mask);
imVeT = AutoCrop(imVeT,mask);
return
%% Pick which map to show
showMe(imKepLL)
