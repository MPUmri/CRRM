function [dirLocation] = DefaultFolders()
    % This sets the default folder locations

    dirLocation.sim = './data/simData'; % Simulation data
    dirLocation.qiba = './data/QIBA-ToftsV6'; % QIBA Phantom data
    dirLocation.qinZip = './data/QIN Breast DCE-MRI'; % Downloaded QIN Breast DCE-MRI data
    dirLocation.qin = './data/QINBreast'; % Unzipped and processed QIN Breast data

    dirLocation.results = './dataResults'; % Results data
end
