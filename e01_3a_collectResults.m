% Result collector. This script:
% - reads in raw results from e01_2a_simProcessLLSQ.m & e01_2b_simProcessNLSQ.m
% - calculates summary statistics
% - saves summary as .csv

% Runtime: ~2 seconds

%% General setup
clearvars
fclose('all')
addpath('./mfiles');

% Make sure this matches the choice used in e01_2*.m, or else bad things may happen
refName = 'refY';
% Choices are: 'refY', 'refWS', 'refP'
% Refer to './mfiles/SimProperties.m for further details

% Name of the csv file which will contain the summary results
csvName = ['e01-simResults-' refName '.csv'];

% This controls which kind of models/approaches should be summarized
methodList = {'LRRM','CLRRM','CHRRM','CLRRMm','CLRRMt','NRRM'};
% Choices can include: 'LRRM', 'CLRRM', 'CLRRMn', 'CLRRMm', 'CLRRMt', 'NRRM'
% Refer to e01_2_simProcessLLSQ.m for details

% End of setup

%% Simulation Properties

simProp = SimProperties(refName);

% Copy variables to make code less verbose
listCNR = simProp.CNR; % Range of Contrast-Noise Ratio to simulate
nVox = simProp.nVox; % Number of replications for each CNR
listTRes = simProp.TRes; % Temporal resolutions (in seconds)
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
outFile = fullfile(dirLocation.results,csvName);

% Create folder for csv output (if it doesn't exist)
if ~exist(dirLocation.results)
    mkdir(dirLocation.results)
end

% Check for inconsistencies
% This is a crude check that only makes sure that the results folder exist
% Obviously, if they exist but are empty, then they pass this check but the
% script will probably choke somewhere else down the line.
if (strmatch('LRRM',methodList) & ~exist(fullfile(dataDir,'LRRM')));
    error('Results for LRRM not found');
end
if (strmatch('CLRRM',methodList) & ~exist(fullfile(dataDir,'CLRRM')));
    error('Results for CLRRM not found');
end
if (strmatch('CHRRM',methodList) & ~exist(fullfile(dataDir,'CHRRM')));
    error('Results for CHRRM not found');
end
if (strmatch('CLRRMm',methodList) & ~exist(fullfile(dataDir,'CLRRMm')));
    error('Results for CLRRMm not found');
end
if (strmatch('CLRRMt',methodList) & ~exist(fullfile(dataDir,'CLRRMt')));
    error('Results for CLRRMt not found');
end
if (strmatch('NRRM',methodList) & ~exist(fullfile(dataDir,'NRRM')));
    error('Results for NRRM not found');
end

%% Initialize output CSV file

% The CSV header
hdr=['FitMethod,CNR,TemporalRes,TrueRelKt,TrueRelVe,TrueKep,errKt,errVe,errKep,stdErrKt,stdErrVe,stdErrKep, meanResid, stdResid'];
outID = fopen(outFile, 'w+');
fprintf(outID, '%s\n', hdr); % Print header into csv file

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
                [ktErrs, meanKtErr, stdKtErr] = PercentError(pkParamsLL(:,1),relKt);
                [veErrs, meanVeErr, stdVeErr] = PercentError(pkParamsLL(:,2),relVe);
                [kepErrs, meanKepErr, stdKepErr] = PercentError(pkParamsLL(:,3),kep);
                % Note down the mean and standard deviation of residual too
                meanResid = mean(residLL(:));
                stdResid = std(residLL(:));
                % Output the mean and stdDeviation to CSV file...
                outLine = {'LRRM',curCNR,curTRes,relKt,relVe,kep,meanKtErr,meanVeErr,...
                    meanKepErr,stdKtErr,stdVeErr,stdKepErr,meanResid,stdResid};
                fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
            elseif strcmp(curM,'CLRRM')
                [ktErrs, meanKtErr, stdKtErr] = PercentError(pkParamsCL(:,1),relKt);
                [veErrs, meanVeErr, stdVeErr] = PercentError(pkParamsCL(:,2),relVe);
                [kepErrs, meanKepErr, stdKepErr] = PercentError(pkParamsCL(:,3),kep);
                meanResid = mean(residCL(:));
                stdResid = std(residCL(:));
                outLine = {'CLRRM',curCNR,curTRes,relKt,relVe,kep,meanKtErr,meanVeErr,...
                    meanKepErr,stdKtErr,stdVeErr,stdKepErr,meanResid,stdResid};
                fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
            elseif strcmp(curM,'CLRRMm')
                [ktErrs, meanKtErr, stdKtErr] = PercentError(pkParamsCLm(:,1),relKt);
                [veErrs, meanVeErr, stdVeErr] = PercentError(pkParamsCLm(:,2),relVe);
                [kepErrs, meanKepErr, stdKepErr] = PercentError(pkParamsCLm(:,3),kep);
                meanResid = mean(residCLm(:));
                stdResid = std(residCLm(:));
                outLine = {'CLRRMm',curCNR,curTRes,relKt,relVe,kep,meanKtErr,meanVeErr,...
                    meanKepErr,stdKtErr,stdVeErr,stdKepErr,meanResid,stdResid};
                fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
            elseif strcmp(curM,'CHRRM')
                [ktErrs, meanKtErr, stdKtErr] = PercentError(pkParamsCH(:,1),relKt);
                [veErrs, meanVeErr, stdVeErr] = PercentError(pkParamsCH(:,2),relVe);
                [kepErrs, meanKepErr, stdKepErr] = PercentError(pkParamsCH(:,3),kep);
                meanResid = mean(residCH(:));
                stdResid = std(residCH(:));
                outLine = {'CHRRM',curCNR,curTRes,relKt,relVe,kep,meanKtErr,meanVeErr,...
                    meanKepErr,stdKtErr,stdVeErr,stdKepErr,meanResid,stdResid};
                fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
            elseif strcmp(curM,'CLRRMt')
                [ktErrs, meanKtErr, stdKtErr] = PercentError(pkParamsCLt(:,1),relKt);
                [veErrs, meanVeErr, stdVeErr] = PercentError(pkParamsCLt(:,2),relVe);
                [kepErrs, meanKepErr, stdKepErr] = PercentError(pkParamsCLt(:,3),kep);
                meanResid = mean(residCLt(:));
                stdResid = std(residCLt(:));
                outLine = {'CLRRMt',curCNR,curTRes,relKt,relVe,kep,meanKtErr,meanVeErr,...
                    meanKepErr,stdKtErr,stdVeErr,stdKepErr,meanResid,stdResid};
                fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
            elseif strcmp(curM,'NRRM')
                % This option contains NRRM and CNRRM
                pkParamsN = pkParamsN';
                pkParamsCN = pkParamsCN';
                % Do NRRM first
                [ktErrs, meanKtErr, stdKtErr] = PercentError(pkParamsN(:,1),relKt);
                [veErrs, meanVeErr, stdVeErr] = PercentError(pkParamsN(:,2),relVe);
                [kepErrs, meanKepErr, stdKepErr] = PercentError(pkParamsN(:,3),kep);
                meanResid = mean(residN(:));
                stdResid = std(residN(:));
                outLine = {'NRRM',curCNR,curTRes,relKt,relVe,kep,meanKtErr,meanVeErr,...
                    meanKepErr,stdKtErr,stdVeErr,stdKepErr,meanResid,stdResid};
                fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
                % Then do CNRRM
                [ktErrs, meanKtErr, stdKtErr] = PercentError(pkParamsCN(:,1),relKt);
                [veErrs, meanVeErr, stdVeErr] = PercentError(pkParamsCN(:,2),relVe);
                [kepErrs, meanKepErr, stdKepErr] = PercentError(pkParamsCN(:,3),kep);
                meanResid = mean(residCN(:));
                stdResid = std(residCN(:));
                outLine = {'CNRRM',curCNR,curTRes,relKt,relVe,kep,meanKtErr,meanVeErr,...
                    meanKepErr,stdKtErr,stdVeErr,stdKepErr,meanResid,stdResid};
                fprintf(outID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n', outLine{:});
            end
        end
    end
end
toc
fclose(outID); % Close the CSV file
fclose('all') % Close everything, for safety incase the previous line doesn't do it
disp('.')
disp('.')
disp('.')
disp('Done summarizing simulation data')
disp(['Summary saved in: ' outFile]);
