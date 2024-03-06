function [selVal,time] = selSolRVSS_PD(objVal,selNum)
    tic;
    %selVal = [];
    [~,M] = size(objVal);
    [W,selNum] = UniformPoint(selNum,M);
    
    
    Cosine   = 1 - pdist2(objVal,W,'cosine');
    Distance = repmat(sqrt(sum(objVal.^2,2)),1,selNum).*sqrt(1-Cosine.^2);
    % Associate each solution with its nearest reference point
    [d,pi] = min(Distance,[],1);
    
    selVal = objVal(pi,:);
   
    time = toc;
end

function d = MSSU(M)
    d1 = UniformVector(10000,M,1);
    d = 
end

function [V,N] = UniformVector(N,M,seed)
    rng(1);
    mu = zeros(1,M);
    sigma = eye(M,M);
    R = mvnrnd(mu,sigma,N);
    V = abs(R./sqrt(sum(R.^2,2)));
end