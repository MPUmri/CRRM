function [pErr, mErr, sErr] = PercentError(estVals, trueVals)

    if isrow(trueVals)
        trueVals = trueVals';
    end
    
    if size(estVals) ~= size(trueVals)
        [nX nY] = size(estVals);
        if nX == length(trueVals)
            trueVals = repmat(trueVals,[1 nY]);
        else
            trueVals = repmat(trueVals',[nX 1]);
        end
    end
    pErr = 100*(estVals-trueVals)./trueVals;
    mErr = mean(pErr(:));
    sErr = std(pErr(:));
end

