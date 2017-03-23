function [simProp] = QibaSimProperties()
    
    simProp.SNR = [5:5:50]; % Range of Contrast-Noise Ratio to simulate
    simProp.nClones = 100; % Number of replications for each SNR
    simProp.tRes = [1,5,10,15,30,60]; % Temporal resolutions (in seconds)
    
end