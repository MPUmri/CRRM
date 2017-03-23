% Calculate the SNR and CNR in the QIN Breast datatset

% Estimated run time: ~19 s

clearvars
addpath('./mfiles')

dirLocation = DefaultFolders();
dataDir = fullfile(dirLocation.qin,'Stripped');

matFiles = dir([dataDir '/*.mat']);

%%
tic
for q=1:length(matFiles)
    curName = matFiles(q).name
    %% Load data
    inFile = fullfile(dataDir, curName);
    load(inFile);

    [nX nY nZ nT] = size(dynData);
    S = reshape(dynData,[nX*nY*nZ nT]);
    S(mask<1,:)=[];
    signal = S(:,1);
    sigma = std(S(:,2)-S(:,1));
    tmpSNR = sqrt(2)*signal./sigma;
    tmpSNR(tmpSNR<0) = [];
    tmpSNR(tmpSNR>500) = [];
    SNR{q} = tmpSNR;

    %%
    Ct = Signal2Conc4QIN(S',dataInfo.R10,dataInfo.TR,dataInfo.Alpha,dataInfo.r1p);
    ct0 = max(Ct);
    sigma = std(Ct(2,:)-Ct(1,:));
    tmpCNR = sqrt(2)*ct0./sigma;
    tmpCNR(tmpCNR<0) = [];
    tmpCNR(tmpCNR>500) = [];
    CNR{q} = tmpCNR';
end
%%
toc
%%
allSNR = 0;
allCNR = 0;
for q=1:length(SNR)
    allSNR = [allSNR; SNR{q}];
    mSNR(q) = nanmean(SNR{q});
    stdSNR(q) = nanstd(SNR{q});
    allCNR = [allCNR; CNR{q}];
    mCNR(q) = nanmean(CNR{q});
    stdCNR(q) = nanstd(CNR{q});
end
allSNR(1)=[];
allCNR(1)=[];
%%
figure
errorbar(mSNR,stdSNR);
title('mean SNR')
figure
errorbar(mCNR,stdSNR)
title('mean CNR')
figure
boxplot(allSNR)
title('all SNR')
figure
boxplot(allCNR)
title('all CNR')
