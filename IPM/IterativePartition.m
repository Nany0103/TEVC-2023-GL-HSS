function [indices,recorder] = IterativePartition(candidate, n, methods)
    recorder = {};
    niter = length(methods);
    recorder.num = zeros(1,niter);
    recorder.time = zeros(1,niter);
    P = candidate;
    index = 1:size(P,1);
    indices = [];
    st = tic;
    for k=1:niter
        sets = division(P,n,methods(k));
        sorted_sets = {};
        removed = [];
        for i=1:n
            fprintf('Objective-%d-th-Iteration-%d-th-Subset\n',k,i);
            I = index(sets{i});
            [frontNum,~,] = NDSort(candidate(I,:),1);
            sorted_sets{i} = I(frontNum==1);
            removed = [removed,I(frontNum~=1)];
        end
        % for some division methods with overlap solutions, unique is neeeded
        if contains(methods(k).name, "overlap")
            removed = unique(removed); 
        end
        index = setdiff(index,removed);
        P = candidate(index,:);
        recorder.num(k)=length(index);
        recorder.time(k)=toc(st);
        indices{k}=index;
    end
end