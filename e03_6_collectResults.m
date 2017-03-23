% Collect results from the QIN Analysis

% Estimated run time: ~1s

clearvars
close all
addpath('./mfiles')

refName = 'refY'; % Choices: 'refY', 'refWS', 'refP'

% Patient of Interest
% Figures for these patients will be shown individually
% This can be a vector, i.e. PoI=[1,3,20] will show figures for the first,
% third and twentieth dataset
PoI = -1; % Use -1 to disable.
%%
dirLocation = DefaultFolders();
synthDir = fullfile(dirLocation.qin,'FromJulia');
llDir = fullfile(dirLocation.qin,['MapsLL-' refName]);

%%
refProp = SimProperties(refName);
ktRR = refProp.KtRR;
veRR = refProp.veRR;
kepRR = refProp.kepRR;

matFiles = dir([llDir '/*.mat']);
%%
tic
for q=1:length(matFiles)
    curName = matFiles(q).name;
    llFile = fullfile(llDir, curName);
    synthFile = fullfile(synthDir, curName);

    load(llFile)
    load(synthFile);
    
    % Delete the 'LLTM' entry from method list (if it exists)
    methodList(find(not(cellfun('isempty',strfind(methodList,'LLTM')))))=[];
    %% Get the parameters from the Tofts model (ground truth)
    % Re-arrange Tofts Model Map
    % Original: [KTrans, kep]
    % Desire: [KTrans, ve, kep]
    pkParamsT = pkParamsLT;
    pkParamsT(:,3) = pkParamsT(:,2);
    pkParamsT(:,2) = pkParamsT(:,1)./pkParamsT(:,3);
    pMask = (pkParamsT(:,1)>0) & (pkParamsT(:,1)<=5) & (pkParamsT(:,2)>0) & (pkParamsT(:,2)<=1);
    pkParamsT(~pMask,:)=[]; % pMask contains the voxels which have reasonable fits
    residLT(~pMask)=[];
    ktT = pkParamsT(:,1);
    veT = pkParamsT(:,2);
    kepT = pkParamsT(:,3);
    % This next step of calcualting errors is unnecessary since they'll all be 0
    % but it is done anyways because it helps initialize the matrices
    errKtT=PercentError(ktT,ktT);
    errKepT=PercentError(kepT,kepT);
    errVeT=PercentError(veT,veT);
    errKt = errKtT;
    errKep = errKepT;
    errVe = errVeT;
    %% Go through the different approaches
    if strmatch('LRRM',methodList)
        pkParamsLL(~pMask,:)=[];
        residLL(~pMask)=[];
        ktLL = ktRR*pkParamsLL(:,1);
        veLL = veRR*pkParamsLL(:,2);
        kepLL = pkParamsLL(:,3);
        errKtLL=PercentError(ktLL,ktT);
        errKepLL=PercentError(kepLL,kepT);
        errVeLL=PercentError(veLL,veT);
        errKt = [errKt errKtLL];
        errKep = [errKep errKepLL];
        errVe = [errVe errVeLL];
        corrKtLL(q) = corr(ktLL,ktT);
        corrKepLL(q) = corr(kepLL,kepT);
        corrKtLL_S(q) = corr(ktLL,ktT,'type','Spearman');
        corrKepLL_S(q) = corr(kepLL,kepT,'type','Spearman');
        cccKtLL(q) = CCC(ktLL,ktT);
        cccKepLL(q) = CCC(kepLL,kepT);
    end
    if strmatch('CLRRM',methodList)
        pkParamsCL(~pMask,:)=[];
        residCL(~pMask)=[];
        ktCL = ktRR*pkParamsCL(:,1);
        veCL = veRR*pkParamsCL(:,2);
        kepCL = pkParamsCL(:,3);
        errKtCL=PercentError(ktCL,ktT);
        errKepCL=PercentError(kepCL,kepT);
        errVeCL=PercentError(veCL,veT);
        errKt = [errKt errKtCL];
        errKep = [errKep errKepCL];
        errVe = [errVe errVeCL];
        corrKtCL(q) = corr(ktCL,ktT);
        corrKepCL(q) = corr(kepCL,kepT);
        corrKtCL_S(q) = corr(ktCL,ktT,'type','Spearman');
        corrKepCL_S(q) = corr(kepCL,kepT,'type','Spearman');
        cccKtCL(q) = CCC(ktCL,ktT);
        cccKepCL(q) = CCC(kepCL,kepT);
    end
    if strmatch('CHRRM',methodList)
        pkParamsCH(~pMask,:)=[];
        residCH(~pMask)=[];
        ktCH = ktRR*pkParamsCH(:,1);
        veCH = veRR*pkParamsCH(:,2);
        kepCH = pkParamsCH(:,3);
        errKtCH=PercentError(ktCH,ktT);
        errKepCH=PercentError(kepCH,kepT);
        errVeCH=PercentError(veCH,veT);
        errKt = [errKt errKtCH];
        errKep = [errKep errKepCH];
        errVe = [errVe errVeCH];
        corrKtCH(q) = corr(ktCH,ktT);
        corrKepCH(q) = corr(kepCH,kepT);
        corrKtCH_S(q) = corr(ktCH,ktT,'type','Spearman');
        corrKepCH_S(q) = corr(kepCH,kepT,'type','Spearman');
        cccKtCH(q) = CCC(ktCH,ktT);
        cccKepCH(q) = CCC(kepCH,kepT);
    end
    
    %% Include the estimates from NRRM and CNRRM
    pkParamsN = pkParamsN';
    pkParamsCN = pkParamsCN';
    
    pkParamsN(~pMask,:)=[];
    residN(~pMask)=[];
    ktN = ktRR*pkParamsN(:,1);
    veN = veRR*pkParamsN(:,2);
    kepN = pkParamsN(:,3);
    errKtN=PercentError(ktN,ktT);
    errKepN=PercentError(kepN,kepT);
    errVeN=PercentError(veN,veT);
    errKt = [errKt errKtN];
    errKep = [errKep errKepN];
    errVe = [errVe errVeN];
    corrKtNL(q) = corr(ktN,ktT);
    corrKepNL(q) = corr(kepN,kepT);
    corrKtNL_S(q) = corr(ktN,ktT,'type','Spearman');
    corrKepNL_S(q) = corr(kepN,kepT,'type','Spearman');
    cccKtNL(q) = CCC(ktN,ktT);
    cccKepNL(q) = CCC(kepN,kepT);
    
    pkParamsCN(~pMask,:)=[];
    residCN(~pMask)=[];
    ktCN = ktRR*pkParamsCN(:,1);
    veCN = veRR*pkParamsCN(:,2);
    kepCN = pkParamsCN(:,3);
    errKtCN=PercentError(ktCN,ktT);
    errKepCN=PercentError(kepCN,kepT);
    errVeCN=PercentError(veCN,veT);
    errKt = [errKt errKtCN];
    errKep = [errKep errKepCN];
    errVe = [errVe errVeCN];
    corrKtCN(q) = corr(ktCN,ktT);
    corrKepCN(q) = corr(kepCN,kepT);
    corrKtCN_S(q) = corr(ktCN,ktT,'type','Spearman');
    corrKepCN_S(q) = corr(kepCN,kepT,'type','Spearman');
    cccKtCN(q) = CCC(ktCN,ktT);
    cccKepCN(q) = CCC(kepCN,kepT);
    
    methodList{end+1} = 'NRRM';
    methodList{end+1} = 'CNRRM';
    
    %% Combine the errors into a single variable which collects voxels from
    % all patients
    if q==1
        errKtAll = errKt;
        errKepAll = errKep;
        errVeAll = errVe;
        allP = p;
    else
        errKtAll = [errKtAll; errKt];
        errKepAll = [errKepAll; errKep];
        errVeAll = [errVeAll; errVe];
        allP = [allP; p];
    end
    
    %% Display the figures for an individual patient - if necessary
%     nM = length(methodList);
%     if any(q==PoI)
%         x=-200:200;
%         figure(1)
%         plot(x,histc(errKt(:,3),x)./length(errKt(:,1)));
%         hold on
%         figure(2)
%         plot(x,histc(errKt(:,5),x)./length(errKt(:,1)));
%         hold on
%     end
    %% Get the kepRR estimate from the different approaches
    estKepRR_L(q) = estKepRRL;
    estKepRR_N(q) = estKepRRN;
    stdKepRR_L(q) = stdKepRRL;
    
    %% Get the runtimes from the different approaches
    rLL(q) = runtimeLL;
    rCL(q) = runtimeCL;
    rCH(q) = runtimeCH;
    rNL(q) = runtimeN;
    rCN(q) = runtimeC;
    %%
    numVox(q) = sum(mask(:));
    numGoodVox(q) = sum(pMask(:));
end
nM = length(methodList);
toc
return
%% EXPORTED VERSIONS
%% Initial setting - changing the order of the percent errors
mySeq = [4 1 5 3 2];
nS = length(mySeq);
%% KTrans
figure
boxplot(errKtAll(:,mySeq+1))
set(gca,'XTick',[1:nS], 'XTickLabel',methodList(mySeq))
ylim([-400 200])
ylabel('Percent error in K^{Trans}')
title(['Percent error in K^{Trans} for All Patients'])
%% kep
figure
boxplot(errKepAll(:,mySeq+1))
set(gca,'XTick',[1:nS], 'XTickLabel',methodList(mySeq))
ylim([-400 200])
ylabel('Percent error in k_{ep}')
title(['Percent error in k_{ep} for All Patients'])
%% ve
figure
boxplot(errVeAll(:,mySeq+1))
set(gca,'XTick',[1:nS], 'XTickLabel',methodList(mySeq))
ylim([-200 200])
ylabel('Percent error in v_e')
title(['Percent error in v_e for All Patients'])
%%
%% Aditional code that is unnecessary
%%
%% KTrans - all methods
figure
boxplot(errKtAll(:,2:end))
set(gca,'XTick',[1:nM], 'XTickLabel',methodList)
ylim([-200 200])
ylabel('Percent error in K^{Trans}')
title(['Percent error in K^{Trans} for All Patients'])
%% kep - all methods
figure
boxplot(errKepAll(:,2:end))
set(gca,'XTick',[1:nM], 'XTickLabel',methodList)
ylim([-400 200])
ylabel('Percent error in k_{ep}')
title(['Percent error in k_{ep} for All Patients'])
%% ve - all methods
figure
boxplot(errVeAll(:,2:end))
set(gca,'XTick',[1:nM], 'XTickLabel',methodList)
ylim([-300 200])
ylabel('Percent error in v_e')
title(['Percent error in v_e for All Patients'])
%% Summary statistics
meanErr = [mean(errKtAll(:,2:end)); mean(errKepAll(:,2:end)); nanmean(errVeAll(:,2:end))];
printmat(meanErr,'Mean Percent Error','KTrans kep ve',strjoin(methodList,' '))
stdErr = [std(errKtAll(:,2:end)); std(errKepAll(:,2:end)); nanstd(errVeAll(:,2:end))]./sqrt(length(errKtAll(:,2:end)));
printmat(stdErr,'StdError of Percent Error','KTrans kep ve',strjoin(methodList,' '))
medianErr = [median(errKtAll(:,2:end)); median(errKepAll(:,2:end)); nanmedian(errVeAll(:,2:end))];
printmat(medianErr,'Median Percent Error','KTrans kep ve',strjoin(methodList,' '))
%%
skewM = [skewness(errKtAll); skewness(errKepAll); skewness(errVeAll)];
printmat(skewM(:,2:end),'Skewness of Percent Errors', 'KTrans kep ve', strjoin(methodList,' '));
%% Export correlation results to csv [not used]
outFile = './dataResults/QinSummary.csv';
hdr=['FitType,Visit,CorrKt,CorrKep,CorrKtS,CorrKepS,CCCKt,CCCKep,kepRR'];
outID = fopen(outFile, 'w+');
fprintf(outID, '%s\n', hdr);

for i=1:20
    if mod(i,2) == 0
        curV = '2';
    else 
        curV = '1';
    end
    
   outLine = {'NRRM', curV, corrKtNL(i), corrKepNL(i), corrKtNL_S(i), corrKepNL_S(i), cccKtNL(i), cccKepNL(i), nan};
   fprintf(outID,'%s,%s,%f,%f,%f,%f,%f,%f,%f\n', outLine{:}); 
   
   outLine = {'LRRM', curV, corrKtLL(i), corrKepLL(i), corrKtLL_S(i), corrKepLL_S(i), cccKtLL(i), cccKepLL(i), nan};
   fprintf(outID,'%s,%s,%f,%f,%f,%f,%f,%f,%f\n', outLine{:}); 
   
   outLine = {'CLRRM', curV, corrKtCL(i), corrKepCL(i), corrKtCL_S(i), corrKepCL_S(i), cccKtCL(i), cccKepCL(i), estKepRR_L(i)};
   fprintf(outID,'%s,%s,%f,%f,%f,%f,%f,%f,%f\n', outLine{:}); 
   
   outLine = {'CNRRM', curV, corrKtCN(i), corrKepCN(i), corrKtCN_S(i), corrKepCN_S(i), cccKtCN(i), cccKepCN(i), estKepRR_N(i)};
   fprintf(outID,'%s,%s,%f,%f,%f,%f,%f,%f,%f\n', outLine{:}); 
   
   outLine = {'CHRRM', curV, corrKtCH(i), corrKepCH(i), corrKtCH_S(i), corrKepCH_S(i), cccKtCH(i), cccKepCH(i), estKepRR_N(i)};
   fprintf(outID,'%s,%s,%f,%f,%f,%f,%f,%f,%f\n', outLine{:}); 
   
end
fclose('all')