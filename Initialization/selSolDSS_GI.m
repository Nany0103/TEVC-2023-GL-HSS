function [selVal,time] = selSolDSS_GI(objVal,selNum)
    tic;
    selVal = [];
    
    [~,index] = max(objVal(:,1));
    selVal = [selVal;objVal(index,:)];
    objVal(index,:) = [];
    
    while size(selVal,1) < selNum
        % Distance calculation, other distances can be used.
        distance = pdist2(objVal,selVal);
        [~,index] = max(min(distance,[],2));
        selVal = [selVal;objVal(index,:)];
        objVal(index,:) = [];
    end    
    time = toc;
end
