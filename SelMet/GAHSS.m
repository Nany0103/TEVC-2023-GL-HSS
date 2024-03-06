function [Subset,runTime] = GAHSS(objVal,selNum,num_vec,seed,ref)
    tic;
    rng(seed);
    [solNum,M] = size(objVal);            
    [W,num_vec] = UniformVector(num_vec, M);
    tensor = zeros(solNum,num_vec);
    r = ref;
    for i=1:solNum
        s = objVal(i,:);
        temp1 = min(abs(s-r)./W,[],2)';        
        tensor(i,:) = temp1;     
    end
    mintensor = tensor;
    selVal = zeros(selNum,M);
    selInd = zeros(1,selNum);
    for num = 1:selNum
        mintensor = min(mintensor, tensor);
        r2hvc = sum(mintensor,2);
        [~,bestindex] = max(r2hvc);
        for i=1:solNum
            s = objVal(i,:);
            temp1 = max((objVal(bestindex,:)-s)./W,[],2)';        
            tensor(i,:) = temp1;     
        end  
        selVal(num,:) = objVal(bestindex,:);
        selInd(1,num) = bestindex;
    end    
    Subset = selVal;
    runTime = toc;
end
