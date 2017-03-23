% This script will unzip the MATLAB version of the QIN Breast DCE-MRI Data
% Detail about the dataset is available at:
% https://wiki.cancerimagingarchive.net/display/Public/QIN+Breast+DCE-MRI
%
% If the matlab version of the data is downloaded, then the downloaded
% folder will be called 'QIN Breast DCE-MRI' which contains individual
% patients, and then the image data is buried deeper in the folder in a zip file
% This script will unzip the data and make it easier to access for future steps
%
% In order for this script to work, the downloaded 'QIN Breast DCE-MRI'
% folder should be placed in the 'data' directory.
% So the data for patient one will be in './data/QIN Breast DCE-MRI/QIN-Breast-DCE-MRI-BC01/1.3.6.1.....'
%
% After running this script, the 'QIN Breast DCE-MRI' folder is no longer 
% needed and can be deleted to free up 8 Gb of space
%
% Estimated run time: ~275 seconds 

%% Initialize
clearvars
addpath('./mfiles')

dirLocation = DefaultFolders();
dataDir = dirLocation.qinZip;
outDir = fullfile(dirLocation.qin,'Unzipped');

% Make output directory if it doesn't exist
if ~exist(outDir,'dir')
    mkdir(outDir)
end

%% Do the job
tic
% Get list of patient directories
patientDirs = dir(dataDir);
patientDirs(1:2)=[]; % Discard the first two dir ('.' and '..')
% Cycle through all patient directories
for q=1:length(patientDirs)
    curPatient = patientDirs(q).name;
    disp(['Processing ' curPatient]);
    curDir = fullfile(dataDir,curPatient);
    % Get the directories for the two visits
    visitDirs = dir(curDir);
    for r=3:4 % Skip the first two because they are '.' and '..'
        % We are here...
        curVisitDir = fullfile(curDir,visitDirs(r).name);
        % ...but we need to go one directory further.
        % So, get the list of directories in current folder (there is only one)...
        lastDirs = dir(curVisitDir);
        % ... and our zip file will be in that directory
        zipDir = fullfile(curVisitDir,lastDirs(3).name);
        % Find the zip file
        zipFile = dir([zipDir '/*.zip']);
        % Now unzip it, unzip it good
        unzip(fullfile(zipDir,zipFile(1).name),outDir);
    end
end
toc
disp('...')
disp('...')
disp('Done')
