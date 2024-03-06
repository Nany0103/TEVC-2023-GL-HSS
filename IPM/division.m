function sets = division(candidate, n, method)
    sets = {};
    [N, M] = size(candidate);
    if method.name == "objective"
        [B,I] = sort(candidate(:,method.dim));
        candidate = I;
        siz = ceil(N/n);
        for i=1:n
            sets{i}=candidate((i-1)*siz+1:min(i*siz,N));
        end
    elseif method.name == "objective_overlap"
        [B,I] = sort(candidate(:,method.dim));
        candidate = I;
        siz = ceil(N/(n+1));
        for i=1:n
            sets{i}=candidate((i-1)*siz+1:min(i*siz+siz,N));
        end
    elseif method.name == "random"
        rng(method.seed);
        candidate = randperm(N);
        siz = ceil(N/n);
        for i=1:n
            sets{i}=candidate((i-1)*siz+1:min(i*siz,N));
        end
    elseif method.name == "cosine"
        rng(method.seed)
        sets = mini_batch_kmeans(candidate,n,1000);    
    else
        error("No such division method!")
    end
end
function res = normalize(val, PF) % nomrmalize the obecjtive values by PF
    N = size(val,1);
    fmin   = min(PF, [], 1);
    fmax   = max(PF,[],1);
    res = (val-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
end