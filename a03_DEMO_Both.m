% This script will collect the output from Julia (i.e. from a02)
% and plot the results for LRRM, CLRRM, NRRM, CNRRM and CHRRM

% When running this script, please select the "change folder" option if prompted
%%
clearvars

addpath('./mfiles')

load('./data/demoData4Julia.mat');
load('./data/demoDataFromJulia.mat');

%% Re-Do the linear fitting
% (it's fast anyways)

tic
[estParamsLRRM, resnormLRRM] = LRRM(Ct,Crr,t);
runtimeLRRM=toc;
tic
[estParamsCLRRM, resnormCLRRM, estKepRR, stdKepRR] = CLRRM(Ct,Crr,t);
runtimeCLRRM=toc;

% CHRRM uses the estimatedKepRR from NRRM
tic
[estParamsCHRRM, resnormCHRRM] = CLRRM(Ct,Crr,t,estKepRR_N);
runtimeCHRRM=toc;

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

% The NRRM and CNRRM estimates tend to be transposed from Julia
% So un-transpose them first
estParamsNRRM = estParamsNRRM';
estParamsCNRRM = estParamsCNRRM';

errEstKTrans_NRRM = PercentError(estParamsNRRM(:,1), KTrans/KTransRR);
errEstVe_NRRM = PercentError(estParamsNRRM(:,2), ve/veRR);
errEstKep_NRRM = PercentError(estParamsNRRM(:,3), kep);

errEstKTrans_CNRRM = PercentError(estParamsCNRRM(:,1), KTrans/KTransRR);
errEstVe_CNRRM = PercentError(estParamsCNRRM(:,2), ve/veRR);
errEstKep_CNRRM = PercentError(estParamsCNRRM(:,3), kep);

errEstKTrans_CHRRM = PercentError(estParamsCHRRM(:,1), KTrans/KTransRR);
errEstVe_CHRRM = PercentError(estParamsCHRRM(:,2), ve/veRR);
errEstKep_CHRRM = PercentError(estParamsCHRRM(:,3), kep);

%% Plot the percent errors to compare the two models

modelOrder = {'LRRM','CLRRM','NRRM','CNRRM','CHRRM'};

figure
boxplot([errEstKTrans_LRRM errEstKTrans_CLRRM errEstKTrans_NRRM errEstKTrans_CNRRM errEstKTrans_CHRRM])
rLine = refline([0 0]);
set(rLine,'LineStyle',':')
set(gca,'XTick',[1:length(modelOrder)], 'XTickLabel',modelOrder)
ylim([-50 50])
ylabel('Percent error in relative K^{Trans}')
title('Percent error in relative K^{Trans}')

figure
boxplot([errEstVe_LRRM errEstVe_CLRRM errEstVe_NRRM errEstVe_CNRRM errEstVe_CHRRM])
rLine = refline([0 0]);
set(rLine,'LineStyle',':')
set(gca,'XTick',[1:length(modelOrder)], 'XTickLabel',modelOrder)
ylim([-20 20])
ylabel('Percent error in relative v_e')
title('Percent error in relative v_e')

figure
boxplot([errEstKep_LRRM errEstKep_CLRRM errEstKep_NRRM errEstKep_CNRRM errEstKep_CHRRM])
rLine = refline([0 0]);
set(rLine,'LineStyle',':')
set(gca,'XTick',[1:length(modelOrder)], 'XTickLabel',modelOrder)
ylim([-100 100])
ylabel('Percent error in k_{ep}')
title('Percent error in k_{ep}')
