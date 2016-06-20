function [pkParams, resid] = LRRM(Ct, Crr, t, doPure)
    % Linear Reference Region fit
    % Inputs:
    % Ct - TxN Matrix where T is time and N is number of voxels
    %    - Concentration-time curve in tissue of interest
    % Crr - Tx1
    %     - Concentration-time curve in reference tissue
    % t  - Tx1
    %    - time, in minutes
    % doPure - boolean (true or false)
    %        - if false, then fitted parameters will be re-arranged to be more user friendly
    % Outputs:
    % pkParams - Nx3 matrix,
    %          - if doPute = false, then pkParams = [Ktrans/Ktrans_RR  v_e/v_e,rr  kep]
    %           else pkParams = [kTrans_TOI/kTrans_RR, kTrans_TOI/ve_RR, kep]

    %% Cárdenas-Rodríguez, J., Howison, C. M., & Pagel, M. D. (2013). MRI 31(4), 497–507. doi:10.1016/j.mri.2012.10.008
    if nargin<4
        doPure = false;
    end
    %%
    stepSize = t(2)-t(1); % Assuming constant step size throughout acquisition

    % Initialize matrices
    [sT, sX] = size(Ct);

    M1 = zeros(sT,1);
    M2 = M1;

    M1 = Crr;
    M2 =  stepSize*cumtrapz(Crr);

    resid = zeros(sX,1);
    pkParams = zeros(sX,3);

    % Solve the linear problem
    % Disabling warnings for speed (possible non-tissue regions in image give warnings)
    warning off

    for i=1:sX
            curCt = squeeze(Ct(:,i));
            M3 = -stepSize*cumtrapz(curCt);
            M = [M1, M2, M3];
            pkParams(i,:) = mldivide(M,curCt);
            resid(i) = norm(curCt-M*pkParams(i,:)');
    end
    warning on

    % pkParams = [kTrans_TOI/kTrans_RR, kTrans_TOI/ve_RR, kTrans_TOI/ve_TOI]
    % desire: pkParams = [kTrans_TOI/ktrans_RR, ve_TOI/ve_RR, kTrans_TOI/ve_TOI]
    if doPure==false
        pkParams(:,2) = pkParams(:,2)./pkParams(:,3);
    end
end
