function g = gradient(data,ind,M,ref)
    s = data(ind,:);
    data(ind,:) = [];
    g = zeros(1,M);
    for i=1:M
        lobj = [1:i-1,i+1:M];
        ind = data(:,i)<s(i);
        data_consider = data(ind,:);
        datal = data_consider(:,lobj);
        sl = s(lobj);
        rl = ref(lobj);
        datalp = max(datal,sl);
        g(i) = stk_dominatedhv(datalp(NDSort(datalp,1)==1,:),rl)-prod(rl-sl);
    end
end