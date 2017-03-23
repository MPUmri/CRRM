% Information collector. This script:
% - reads in raw results from e01_2a_simProcessLLSQ.m & e01_2b_simProcessNLSQ.m
% - collects the following features:
%   - the runtimes of each method
%   - the kepRR estimates from the Linear and Non-Linear approach

% Estimated run time: < 1s

clearvars
fclose('all')
addpath('./mfiles');

refName = 'refY'; % Match this with e01_2*
% Choices are:
% 'refY'
% 'refWS'
% 'refP'

%% Simulation Properties
simProp = SimProperties(refName);
listCNR = simProp.CNR; % Range of Contrast-Noise Ratio to simulate
nVox = simProp.nVox; % Number of replications for each CNR
listTRes = simProp.TRes; % Temporal resolutions (in seconds)
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

%%
dirLocation = DefaultFolders();
dataDir = fullfile(dirLocation.sim,simProp.name);

%%
tic
for i = 1:length(listCNR)
    for j = 1:length(listTRes)
        curFile = ['Sim-CNR-' num2str(listCNR(i)) ...
            '-TRes-' num2str(listTRes(j)) '.mat'];
        clrrmFile = fullfile(dataDir, 'CLRRM', curFile);
        lrrmFile = fullfile(dataDir, 'LRRM', curFile);
        nrrmFile = fullfile(dataDir, 'NRRM', curFile);
        load(clrrmFile);
        load(lrrmFile);
        load(nrrmFile);
        % Note down the run times
        rN(i,j) = runtimeN;
        rCN(i,j) = runtimeC;
        rL(i,j) = runtimeLL;
        rCL(i,j) = runtimeCL;
        % Note down the estimated kepRR from LRRM and NRRM
        estKepRRL(i,j) = estKepRRCL;
        estKepRRN(i,j) = medianKepRR;
        % Also, note down the mean kepRR (instead of the mean)
        % estKepRRLm(i,j) = mean(pkParamsLLRaw(:,2)./pkParamsLLRaw(:,1));
        % estKepRRNm(i,j) = mean(pkParamsN(4,:));
    end
end
toc
return
% Code ends here, the next sections look into different features

%% Runtime as a function of DCE time points
noiseInd = 1; % Choose the CNR (index in variable 'listCNR')

% Plot the run time and do a polynomial fit
numTime = 600./listTRes; % Number of time points
x = 1:max(numTime);
fitRL = polyfit(numTime,rL(noiseInd,:),1); % Linear fit for Linear RRM
fitRN = polyfit(numTime,rN(noiseInd,:),2); % Quadratic fir for Non-Linear RRM

yRL = polyval(fitRL,x);
yRN = polyval(fitRN,x);

figure
scatter(numTime,rN(noiseInd,:))
hold on
plot(x,yRN);
hold off
xlabel('Number of DCE Frames')
ylabel('Run time [s]')
title(['Run time vs number of DCE timepoints - NRRM'])
legend('Measured','Quadratic Fit','Location','southeast')
figure
scatter(numTime,rL(noiseInd,:))
hold on
plot(x,yRL);
hold off
xlabel('Number of DCE Frames')
ylabel('Run time [s]')
title(['Run time vs number of DCE timepoints - LRRM'])
legend('Measured','Linear Fit','Location','southeast')

%% Runtime as function of DCE time points - for CLRRM & CNRRM at ALL noise

% Plot the run time and do a polynomial fit
numTime = 600./listTRes; % Number of time points
x = 1:max(numTime);
tempX=repmat(numTime,[10 1]);
fitRL = polyfit(tempX(:),rL(:),1); % Linear fit for Linear RRM
fitRCL = polyfit(tempX(:),rCL(:),1);
fitRN = polyfit(tempX(:),rN(:),2); % Quadratic fir for Non-Linear RRM
fitRCN = polyfit(tempX(:),rCN(:)+rN(:),2);

yRL = polyval(fitRL,x);
yRCL = polyval(fitRCL,x);
yRN = polyval(fitRN,x);
yRCN = polyval(fitRCN,x);

figure
scatter(tempX(:),rCN(:)+rN(:),'o')
hold on
plot(x,yRCN);
hold off
xlabel('Number of DCE Frames')
ylabel('Run time [s]')
title(['Run time vs number of DCE timepoints - NRRM'])
legend('Measured','Quadratic Fit','Location','southeast')
figure
scatter(tempX(:),rCL(:))
hold on
plot(x,yRCL);
hold off
xlabel('Number of DCE Frames')
ylabel('Run time [s]')
title(['Run time vs number of DCE timepoints - LRRM'])
legend('Measured','Linear Fit','Location','southeast')

%% Runtime as a function of Noise
% This is actually the runtime ratio where the reference value is the
% runtime at the lowest noise level.
% The question being asked here is whether noise affects the run time
% The answer seems to be... not really?
figure
plot(listCNR,rCN./repmat(rCN(end,:),[10 1]))
xlabel('CNR')
ylabel('Runtime Ratio')
title(['Runtime Ratio vs Noise at different temporal resolutions - CNRRM'])
legend('1s','5s','10s','15s','30s','60s')
figure
plot(listCNR,rCL./repmat(rCL(end,:),[10 1]))
xlabel('CNR')
ylabel('Runtime Ratio')
title(['Runtime Ratio vs Noise at different temporal resolutions - CLRRM'])
legend('1s','5s','10s','15s','30s','60s')

%% EstKepRR as a function of Noise - NRRM vs LRRM
TResInd = 6; % Choose the TemporalResolution (index in variable: listTRes)

figure
scatter(listCNR,PercentError(estKepRRL(:,TResInd),1),'^')
hold on
scatter(listCNR,PercentError(estKepRRN(:,TResInd),1),'o');
hold off
xlabel('CNR')
ylabel('Percent Error in Estimated k_{ep,RR}')
title(['Percent Error in Estimated k_{ep,RR} for TRes = ' num2str(listTRes(TResInd)) 's'])
legend('LRRM','NRRM')
% Use absolute percent error
figure
scatter(listCNR,abs(PercentError(estKepRRL(:,TResInd),1)),'^')
hold on
scatter(listCNR,abs(PercentError(estKepRRN(:,TResInd),1)),'o');
hold off
xlabel('CNR')
ylabel('Percent Error in Estimated k_{ep,RR}')
title(['Absolute Percent Error in Estimated k_{ep,RR} for TRes = ' num2str(listTRes(TResInd)) 's'])
legend('LRRM','NRRM')

%% EstKepRR as a function of TRes - NRRM vs LRRM
noiseInd = 1; % Choose the CNR (index in variable: 'listCNR')

figure
scatter(listTRes,PercentError(estKepRRL(noiseInd,:),1),'^')
hold on
scatter(listTRes,PercentError(estKepRRN(noiseInd,:),1),'o');
hold off
xlabel('Temporal Resolution [s]')
ylabel('Percent Error in Estimated k_{ep,RR}')
title(['Percent Error in Estimated k_{ep,RR} for CNR = ' num2str(listCNR(noiseInd))])
legend('LRRM','NRRM')
% Use absolute percent error
figure
scatter(listTRes,abs(PercentError(estKepRRL(noiseInd,:),1)),'^')
hold on
scatter(listTRes,abs(PercentError(estKepRRN(noiseInd,:),1)),'o');
hold off
xlabel('Temporal Resolution [s]')
ylabel('Percent Error in Estimated k_{ep,RR}')
title(['Absolute Percent Error in Estimated k_{ep,RR} for CNR = ' num2str(listCNR(noiseInd))])
legend('LRRM','NRRM')

%% Estimated kep,RR from LRRM
figure
scatter(listCNR,estKepRRL(:,1),'^')
hold on
scatter(listCNR,estKepRRL(:,4),'o')
hold on
scatter(listCNR,estKepRRL(:,6),'v')
hold off
xlabel('CNR')
ylabel('Estimated k_{ep,RR}')
title('Estimated k_{ep,RR} for different noise levels and temporal resolutions - LRRM')
legend('1s','15s','60s','Location','southeast')
%% Estimated kep,RR from NRRM
figure
scatter(listCNR,estKepRRN(:,1),'^')
hold on
scatter(listCNR,estKepRRN(:,4),'o')
hold on
scatter(listCNR,estKepRRN(:,5),'v')
hold off
xlabel('CNR')
ylabel('Estimated k_{ep,RR}')
title('Estimated k_{ep,RR} for different noise levels and temporal resolutions - NRRM')
legend('1s','15s','30s')
%% Difference of absolute percent error in kepRR
figure
scatter(listCNR,abs(PercentError(estKepRRL(:,1),1))-abs(PercentError(estKepRRN(:,1),1)),'^')
hold on
scatter(listCNR,abs(PercentError(estKepRRL(:,3),1))-abs(PercentError(estKepRRN(:,3),1)),'o')
scatter(listCNR,abs(PercentError(estKepRRL(:,5),1))-abs(PercentError(estKepRRN(:,5),1)),'v')
hold off
xlabel('CNR')
ylabel('Percent Error Difference')
title(['Difference in absolute percent error for kepRR (LRRM-NRRM)'])
legend('1s','15s','30s')
