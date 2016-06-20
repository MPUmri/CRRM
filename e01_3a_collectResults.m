% Result collector. This script:
% - reads in raw results from e01_2a_simProcessLLSQ.m & e01_2b_simProcessNLSQ.m
% - calculates summary statistics

% Runtime: ~2 seconds

%% General setup
clearvars
fclose('all')
addpath('./mfiles');

% Make sure this matches the choice used in e01_2*.m, or else bad things may happen
refName = 'refY';
% Choices are: 'refY', 'refWS', 'refP'
% Refer to './mfiles/SimProperties.m for further details

% This controls which kind of models/approaches should be summarized
methodList = {'LRRM','CLRRM','NRRM'};
% Choices can include: 'LRRM', 'CLRRM', 'CLRRMn', 'CLRRMm', 'CLRRMt', 'NRRM'
% Refer to e01_2_simProcessLLSQ.m for details

% End of setup

%% Simulation Properties

simProp = SimProperties(refName);

% Copy variables to make code less verbose
listCNR = [5]; % Range of Contrast-Noise Ratio to simulate
nVox = simProp.nVox; % Number of replications for each CNR
listTRes = [10]; % Temporal resolutions (in seconds)
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
dataDir = fullfile(dirLocation.sim,simProp.name);

% Check for inconsistencies
% This is a crude check that only makes sure that the results folder exist
% Obviously, if they exist but are empty, then they pass this check but the
% script will probably choke somewhere else down the line.
if (strmatch('LRRM',methodList) & ~exist(fullfile(dataDir,'LRRM')));
    error('Results for LRRM not found');
end
if (strmatch('CLRRM',methodList) & ~exist(fullfile(dataDir,'CLRRM')));
    error('Results for LRRM not found');
end
if (strmatch('CLRRMm',methodList) & ~exist(fullfile(dataDir,'CLRRMm')));
    error('Results for LRRM not found');
end
if (strmatch('CLRRMn',methodList) & ~exist(fullfile(dataDir,'CLRRMn')));
    error('Results for LRRM not found');
end
if (strmatch('CLRRMt',methodList) & ~exist(fullfile(dataDir,'CLRRMt')));
    error('Results for LRRM not found');
end
if (strmatch('NRRM',methodList) & ~exist(fullfile(dataDir,'NRRM')));
    error('Results for NRRM not found');
end

%% Analyze the simulation data
tic

% Cycle through all concentration-noise ratios
for indCNR = 1:length(listCNR)
    curCNR = listCNR(indCNR);
    % Cycle through all temporal resolutions
    for indTRes = 1:length(listTRes)
        curTRes = listTRes(indTRes);
        % Cycle through the method list
        for i=1:length(methodList)
            curM = methodList{i};
            curFile = ['Sim-CNR-' num2str(curCNR) '-TRes-' num2str(curTRes) '.mat'];
            load(fullfile(dataDir,curM,curFile));
            if strcmp(curM,'LRRM')
                % Calculate percent error for relative KTrans, relative EES and kep
                [ktErrsLL, meanKtErrLL, stdKtErrLL] = PercentError(pkParamsLL(:,1),relKt);
                [veErrsLL, meanVeErrLL, stdVeErrLL] = PercentError(pkParamsLL(:,2),relVe);
                [kepErrsLL, meanKepErrLL, stdKepErrLL] = PercentError(pkParamsLL(:,3),kep);
                % Note down the mean and standard deviation of residual too
                meanResidLL = mean(residLL(:));
                stdResidLL = std(residLL(:));
                % Output the mean and stdDeviation to CSV file...
            elseif strcmp(curM,'CLRRM')
                [ktErrsCL, meanKtErrCL, stdKtErrCL] = PercentError(pkParamsCL(:,1),relKt);
                [veErrsCL, meanVeErrCL, stdVeErrCL] = PercentError(pkParamsCL(:,2),relVe);
                [kepErrsCL, meanKepErrCL, stdKepErrCL] = PercentError(pkParamsCL(:,3),kep);
                meanResidCL = mean(residCL(:));
                stdResidCL = std(residCL(:));
            elseif strcmp(curM,'CLRRMm')
                [ktErrsCLm, meanKtErrCLm, stdKtErrCLm] = PercentError(pkParamsCLm(:,1),relKt);
                [veErrsCLm, meanVeErrCLm, stdVeErrCLm] = PercentError(pkParamsCLm(:,2),relVe);
                [kepErrsCLm, meanKepErrCLm, stdKepErrCLm] = PercentError(pkParamsCLm(:,3),kep);
                meanResidCLm = mean(residCLm(:));
                stdResidCLm = std(residCLm(:));
            elseif strcmp(curM,'CLRRMn')
                [ktErrsCLn, meanKtErrCLn, stdKtErrCLn] = PercentError(pkParamsCLn(:,1),relKt);
                [veErrsCLn, meanVeErrCLn, stdVeErrCLn] = PercentError(pkParamsCLn(:,2),relVe);
                [kepErrsCLn, meanKepErrCLn, stdKepErrCLn] = PercentError(pkParamsCLn(:,3),kep);
                meanResidCLn = mean(residCLn(:));
                stdResidCLn = std(residCLn(:));
            elseif strcmp(curM,'CLRRMt')
                [ktErrsCLt, meanKtErrCLt, stdKtErrCLt] = PercentError(pkParamsCLt(:,1),relKt);
                [veErrsCLt, meanVeErrCLt, stdVeErrCLt] = PercentError(pkParamsCLt(:,2),relVe);
                [kepErrsCLt, meanKepErrCLt, stdKepErrCLt] = PercentError(pkParamsCLt(:,3),kep);
                meanResidCLt = mean(residCLt(:));
                stdResidCLt = std(residCLt(:));
            elseif strcmp(curM,'NRRM')
                % This option contains NRRM and CNRRM
                pkParamsN = pkParamsN';
                pkParamsCN = pkParamsCN';
                % Do NRRM first
                [ktErrsNL, meanKtErrNL, stdKtErrNL] = PercentError(pkParamsN(:,1),relKt);
                [veErrsNL, meanVeErrNL, stdVeErrNL] = PercentError(pkParamsN(:,2),relVe);
                [kepErrsNL, meanKepErrNL, stdKepErrNL] = PercentError(pkParamsN(:,3),kep);
                meanResidNL = mean(residN(:));
                stdResidNL = std(residN(:));
                % Then do CNRRM
                [ktErrsCN, meanKtErrCN, stdKtErrCN] = PercentError(pkParamsCN(:,1),relKt);
                [veErrsCN, meanVeErrCN, stdVeErrCN] = PercentError(pkParamsCN(:,2),relVe);
                [kepErrsCN, meanKepErrCN, stdKepErrCN] = PercentError(pkParamsCN(:,3),kep);
                meanResidCN = mean(residCN(:));
                stdResidCN = std(residCN(:));
            end
        end
    end
end
toc
disp('.')
disp('.')
disp('.')
disp('Done summarizing simulation data')

%%
if exist('ktErrsNL','var')
    xLabels = {'NRRM','CNRRM','LRRM','CLRRM'};
    %%
    figure
    boxplot([ktErrsNL ktErrsCN ktErrsLL ktErrsCL])
    set(gca,'XTick',[1:4], 'XTickLabel',xLabels)
    ylabel('Percent error in K^{Trans}')
    title(['Percent error in K^{Trans}'])
    %%
    figure
    boxplot([kepErrsNL kepErrsCN kepErrsLL kepErrsCL])
    set(gca,'XTick',[1:4], 'XTickLabel',xLabels)
    ylabel('Percent error in k_{ep}')
    title(['Percent error in k_{ep}'])
    %%
    figure
    boxplot([veErrsNL veErrsCN veErrsLL veErrsCL])
    set(gca,'XTick',[1:4], 'XTickLabel',xLabels)
    ylabel('Percent error in v_e')
    title(['Percent error in v_e'])
    ylim([-200 200])
else
    xLabels = {'LRRM','CLRRM'};
    %%
    figure
    boxplot([ktErrsLL ktErrsCL])
    set(gca,'XTick',[1:2], 'XTickLabel',xLabels)
    ylabel('Percent error in K^{Trans}')
    title(['Percent error in K^{Trans}'])
    %%
    figure
    boxplot([kepErrsLL kepErrsCL])
    set(gca,'XTick',[1:2], 'XTickLabel',xLabels)
    ylabel('Percent error in k_{ep}')
    title(['Percent error in k_{ep}'])
    %%
    figure
    boxplot([veErrsLL veErrsCL])
    set(gca,'XTick',[1:2], 'XTickLabel',xLabels)
    ylabel('Percent error in v_e')
    title(['Percent error in v_e'])
    ylim([-200 200])
end