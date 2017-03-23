% Quick demo using the MATLAB implementation of the CLRRM.
% Some data is simulated and then the LRRM and CLRRM are fitted to it

% When running this script, please select the "change folder" option if prompted
%% 
clearvars

addpath('./mfiles')

% Directory to save data which will be used by Julia for CNRRM
outDir = './data';
if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% Simulate some data

rng(12345); % Set arbitrary seed for reproducibility

% Simulated reference region properties
KTransRR = 0.1;
kepRR = 1.0;
veRR = KTransRR ./ kepRR;

% Tissue of interest properties
KTrans = 0.25;
kep = 0.625;
ve = KTrans ./ kep;

% Temporal Resolution (in seconds)
tRes = 5;

% Concentration-Noise Ratio (lower value -> more noise)
CNR = 10;
nRep = 1000; % Number of replications (copies) in which noise is added

% Define the time-span
tEnd = 10; % Duration of acquisition (in minutes - most other variables are in seconds)
tInject = 60; % Time of injection (in seconds)
initialTRes = 0.1; % Simulate at 0.1 s initially
t = initialTRes:initialTRes:tEnd*60;
t = t'/60; % Convert time from seconds into minutes

% Generate the AIF - use the injection time defined earlier
Cp = ParkerAif(t,t(tInject./initialTRes)); 

% Simulate the reference regions concentration-time curve
Crr = TrapzKety(Cp,[KTransRR kepRR],t)';

% Simulate the tissue of interest's concentration-time curve
Ct = TrapzKety(Cp,[KTrans kep],t)';

% Add noise to Ct
sigmaC = max(Ct(:))/CNR; % Determine sigma for the noise
Ct = repmat(Ct,[1 nRep]); % Make the copies of Ct
Ct = Ct + sigmaC * randn(size(Ct)); % Add gaussian noise to each Ct copy

% Downsample
t = downsample(t,tRes/initialTRes);
Ct = downsample(Ct,tRes/initialTRes);
Crr = downsample(Crr,tRes/initialTRes);

%% Do the fitting

% The implemented LRRM and CLRRM require that Ct be nT-by-nV matrix, where
% nT is number of frames, and nV is number of voxels
% The other two inputs are Crr (concentration in RefRegion) and t (time, in minutes)

tic
[estParamsLRRM, resnormLRRM] = LRRM(Ct,Crr,t);
runtimeLRRM=toc;
tic
[estParamsCLRRM, resnormCLRRM, estKepRR, stdKepRR] = CLRRM(Ct,Crr,t);
runtimeCLRRM=toc;

% The format of estParams is:
% [ KTrans/KTransRR    ,    ve/veRR    ,    kep    ] 
% So, if we knew KTransRR and veRR, then:
%     estKTrans = KTransRR * estParams(:,1);
%     estVe = veRR * estParams(:,2);
% Otherwise, the relative parameters themselves could be useful on their own

%% Compute percent errors

% Percent error in relative KTrans
errEstKTrans_LRRM = PercentError(estParamsLRRM(:,1), KTrans/KTransRR);
% Percent error in relative Ve
errEstVe_LRRM = PercentError(estParamsLRRM(:,2), ve/veRR);
% Percent error in kep
errEstKep_LRRM = PercentError(estParamsLRRM(:,3), kep);

errEstKTrans_CLRRM = PercentError(estParamsCLRRM(:,1), KTrans/KTransRR);
errEstVe_CLRRM = PercentError(estParamsCLRRM(:,2), ve/veRR);
errEstKep_CLRRM = PercentError(estParamsCLRRM(:,3), kep);

%% Plot the percent errors to compare the two models

modelOrder = {'LRRM','CLRRM'};

figure
boxplot([errEstKTrans_LRRM errEstKTrans_CLRRM])
rLine = refline([0 0]);
set(rLine,'LineStyle',':')
set(gca,'XTick',[1:2], 'XTickLabel',modelOrder)
ylim([-50 50])
ylabel('Percent error in relative K^{Trans}')
title('Percent error in relative K^{Trans}')

figure
boxplot([errEstVe_LRRM errEstVe_CLRRM])
rLine = refline([0 0]);
set(rLine,'LineStyle',':')
set(gca,'XTick',[1:2], 'XTickLabel',modelOrder)
ylim([-20 20])
ylabel('Percent error in relative v_e')
title('Percent error in relative v_e')

figure
boxplot([errEstKep_LRRM errEstKep_CLRRM])
rLine = refline([0 0]);
set(rLine,'LineStyle',':')
set(gca,'XTick',[1:2], 'XTickLabel',modelOrder)
ylim([-100 100])
ylabel('Percent error in k_{ep}')
title('Percent error in k_{ep}')

%% Save original data for Julia
save(fullfile(outDir,'demoData4Julia.mat'),'Ct','Crr','t','KTrans','KTransRR','ve','veRR','kep','kepRR');