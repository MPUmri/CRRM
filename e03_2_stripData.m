% Pre-processing of QIN Breast DCE-MRI data.
% This script will strip away all non-tumour slices and build a tumour mask.
% This is done to save time. The original data covers the entire breast,
% but we are only interested in the slices containing the tumour.

% This step was done because during the research work, the data was processed
% multiple times for various reasons. Loading the data was the slowest part
% of the procedure, so the non-tumour slices were stripped away in order to
% decrease wait times and keep RAM usage low.

% After running the script, the './data/QINBreast' folder is no longer needed
% and can be deleted to free up 8 Gb of space

% Run time: ~415 seconds

%% General setup
clearvars
addpath('./mfiles')

doOverwrite = 0; % Overwrite any existing output

%% Initialize
dirLocation = DefaultFolders();
matDir = fullfile(dirLocation.qin,'Unzipped');
outDir = fullfile(dirLocation.qin,'Stripped');

% Make output directory if it doesn't exist
if ~exist(outDir,'dir')
    mkdir(outDir)
end

%% Do the job
tic
% Get list of .mat files
matFiles = dir([matDir '/*.mat']);

%Cycle through .mat files
for q=1:length(matFiles)
    % Clear variables to keep memory usage low
    clearvars -except matDir matFiles outDir q doOverwrite
    % Load data
    curName = matFiles(q).name;
    curFile = fullfile(matDir,curName);
    outFile = fullfile(outDir,curName);
    if ~doOverwrite && exist(outFile,'file')
        % Skip file if output already exists
        continue;
    end
    disp(['Processing ' curName]);
    load(curFile);
    % Grab the information struct, and add info about rows/column/slices
    dataInfo = da.para;
    dataInfo.nX = da.hdr.Rows;
    dataInfo.nY = da.hdr.Columns;
    dataInfo.nZ = numel(da.slcList);
    % Grab the data from only the tumour-containing slices & grab AIF
    dynData = da.data(:,:,da.slcList,:);
    aif=da.aif;
    % Build the mask for the tumour region
    mask = zeros(dataInfo.nX, dataInfo.nY, dataInfo.nZ);
    for r=1:numel(da.slcList)
        curSlice = da.slcList(r);
        curMask = zeros(320*320,1);
        curMask(da.slcdata(curSlice).Q)=1;
        mask(:,:,r)=reshape(curMask,[dataInfo.nX dataInfo.nY]);
    end
    % Patient BC05 seems to have extra precontrast frame? Remove it
    if ~isempty(findstr(curName,'BC05'))
        dynData(:,:,:,1)=[];
        dataInfo.BaseLineEnd = dataInfo.BaseLineEnd-1;
    end
    % Save to output directory, using the same filename as the input file
    save(outFile, 'dynData', 'mask', 'aif', 'dataInfo');
end
toc
disp('...')
disp('...')
disp('done')
