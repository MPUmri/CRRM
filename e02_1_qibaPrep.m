% This script pre-processes QIBA Virtual Phantom Data
% More specifically, it converts the signal data to concentration data

% The code was adapted from the 'test/q6/prep6.py' script of DCEMRI.jl
% https://github.com/davidssmith/DCEMRI.jl

% Estimated run time: ~4 seconds

%% Initialize
clearvars
addpath('./mfiles');

dirLocation = DefaultFolders();

dicomDir = fullfile(dirLocation.qiba,'DICOM');
outDir = dirLocation.qiba;

%% Load data

tic

% Get list of dicom files from dicom folder
dcmFiles = dir([dicomDir '/*.dcm']);

% Read single dicom file and get some basic info from dicom header
% (This could be replaced with hardcoded values from phantom documentation)
aFile = fullfile(dicomDir,dcmFiles(1).name);
x = dicominfo(aFile);
flipAngle = x.FlipAngle * pi / 180;
TR = x.RepetitionTime * 1E-3;
nX = x.Rows;
nY = x.Columns;

% Some values are just known (from phantom documentation)
nT = length(dcmFiles); 	% Number of frames
deltaT = 0.5; 					% Time between frames
r = 4.5;								% Contrast agent relaxivity
hct = 0.45; 						% Hematocrit
plasmaT1 = 1.440;				% Blood plasma T1
plasmaR1 = 1 / plasmaT1;

% Load all dicom data
dcmData = zeros(nX, nY, nT);
for i=1:nT
	aFile = fullfile(dicomDir, dcmFiles(i).name);
	dcmData(:,:,i) = double(dicomread(aFile));
end

%% Clean up the phantom

% Take signal from the AIF-only compartment
aifData = squeeze(dcmData(71:end,:,:));

% Remove the AIF compartment from our dceData so we have tissue-only data,
% and also remove the 'timer' compartment
dceData = dcmData(11:70,:,:);
[nX, nY, nT] = size(dceData);

%% Obtain AIF

% AIF signal
aifS = reshape(aifData, [10*nY nT]);
aifS = mean(aifS);

% Convert signal to concentration
aifS0 = mean(aifS(1,1:120),2); % Pre-injection signal. A mean isn't necessary since it is noiseless.
aifSN = aifS / aifS0;
E0 = exp(-plasmaR1 * TR);
E = (1 - aifSN + aifSN*E0 - E0*cos(flipAngle)) ./ ...
	(1 - aifSN*cos(flipAngle) + aifSN*E0*cos(flipAngle) - E0*cos(flipAngle));
bloodAif = ((-1/TR)*log(E) - plasmaR1) / r;

% Blood AIF is concentration in whole blood.
% It has to be converted to concentration in plasma.
aif = bloodAif / (1-hct);
aif(aif<0) = 0;

%% Convert the main phantom data to concentration data and export
% (This part is poorly commented, but there's not much to say)

% %% Export large version (unnecessary, but might be useful in future) %%
T = 0:0.5:0.5*(nT-1);
T = T/60;
[sX, sY, sT] = size(dceData);
mapT1 = ones(sX, sY);	% T1 is 1000 ms
mapR1 = 1 ./ mapT1;
mapS0 = ones(sX, sY, sT) * 50000; % Equlibrium magnetization in tissue

r1Data = zeros(sX, sY, sT);
S0 = squeeze(mean(dceData(:,:,1:120),3));
E0 = exp(-mapR1*TR);
E0 = repmat(E0, [1 1 sT]);
A = dceData ./ repmat(S0,[1 1 sT]);
E = (1 - A + A.* E0 - E0 * cos(flipAngle)) ./ ...
    (1 - A*cos(flipAngle) + A.*E0*cos(flipAngle) - E0*cos(flipAngle));
r1Data = (-1/TR) * log(E);
concData = (r1Data - repmat(mapR1,[1 1 sT]))/r;

valKt = [0.01,0.02,0.05,0.1,0.2,0.35];
valVe = [0.01,0.05,0.1,0.2,0.5];
injectionTime = T(120);
precontrastSignal = squeeze(dceData(1,1,1));
save([outDir '/QIBAv6-Full.mat'],'dceData','concData','mapR1','TR','flipAngle','r',...
    'T','aif','valKt','valVe','injectionTime','precontrastSignal');

% %% Export a subsampled version %%
T = 0:0.5:0.5*(nT-1);
T = T/60;

% Pick a single voxel form each compartment
% There is no need to average the compartment voxels since data is noiseless
dceData = dceData(5:10:(end-5), 5:10:(end-5),:);
[sX, sY, sT] = size(dceData);
% % Alternative method if averaging is required:
% blockFilter = ones(10,10,1)/100;
% smallData = convn(dceData, blockFilter);
% smallImage = smallData(6:10:(end-6), 6:10:(end-6),:);

mapT1 = ones(sX, sY);	% T1 is 1000 ms
mapR1 = 1 ./ mapT1;
mapS0 = ones(sX, sY, sT) * 50000; % Equlibrium magnetization in tissue

r1Data = zeros(sX, sY, sT);
S0 = squeeze(mean(dceData(:,:,1:120),3));
E0 = exp(-mapR1*TR);
E0 = repmat(E0, [1 1 sT]);
A = dceData ./ repmat(S0,[1 1 sT]);
E = (1 - A + A.* E0 - E0 * cos(flipAngle)) ./ ...
    (1 - A*cos(flipAngle) + A.*E0*cos(flipAngle) - E0*cos(flipAngle));
r1Data = (-1/TR) * log(E);
concData = (r1Data - repmat(mapR1,[1 1 sT]))/r;

valKt = [0.01,0.02,0.05,0.1,0.2,0.35];
valVe = [0.01,0.05,0.1,0.2,0.5];
injectionTime = T(120);
precontrastSignal = squeeze(dceData(1,1,1));
save([outDir '/QIBAv6-Mini.mat'],'dceData','concData','mapR1','TR','flipAngle','r',...
    'T','aif','valKt','valVe','injectionTime','precontrastSignal');
%% Julia will be used for some steps, and requires its own version
DCEdata = permute(dceData,[3 1 2]);
relaxivity = r;
DCEflip = flipAngle*180/pi;
t = T*60;
Cp = aif;
R10 = mapR1;
S0 = mapS0(:,:,1);
mask = ones(size(R10));
save([outDir '/QIBAv6-Mini4Jl.mat'],'DCEdata','R10','S0','TR','DCEflip',...
    'relaxivity','t','Cp','mask');
%% Save another version for Julia which is downsampled - for warmup runs
t=downsample(t,100);
DCEdata = downsample(DCEdata,100);
Cp = downsample(Cp,100);
save([outDir '/QIBAv6-Mini4Jl-Warmup.mat'],'DCEdata','R10','S0','TR','DCEflip',...
    'relaxivity','t','Cp','mask');
%% Finish
toc
disp('...')
disp('...')
disp('Done preparing QIBA phantom data');
