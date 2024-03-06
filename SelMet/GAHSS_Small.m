function [Subset,selInd] = GAHSS_Small(objVal,selNum,ref,seed,num_vec)
    rng(seed);
    [solNum,M] = size(objVal);
    W = UniformVector(num_vec,M);
    tensor = inf(solNum,num_vec);
    r = ref;
    for i=1:num_vec
        temp = min(abs(objVal-r)./W(i,:),[],2);        
        tensor(:,i) = temp;     
    end
    mintensor = tensor;
    selVal = zeros(selNum,M);
    selInd = zeros(1,selNum);
    for num = 1:selNum
        mintensor = min(mintensor, tensor);
        r2hvc = sum(mintensor,2);
        [~,bestindex] = max(r2hvc);
        for i=1:num_vec
            temp =  max((objVal(bestindex,:)-objVal)./W(i,:),[],2);
            tensor(:,i) = temp;      
        end
        selVal(num,:) = objVal(bestindex,:);
        selInd(1,num) = bestindex;
    end    
    Subset = selVal;
end
