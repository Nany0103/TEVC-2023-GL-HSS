function [selData,runTime] = TGIHSS(data,selNum,ref,seed)
    tic;
    alpha = 20;
    alpha1 = floor(alpha/2);
    alpha2 = alpha-alpha1;
    [~,M] = size(data);
    ind = zeros(1,alpha*selNum);
    for i = 1:alpha1
        [~,ind((i-1)*selNum+1:i*selNum)] = GAHSS_Small(data,selNum,ref,(seed-1)*alpha1+i,1);
    end
    [~,ind(alpha1*selNum+1:end)] = DSS(data,alpha2*selNum);
    selData = data(ind,:);
    selData = LGIHSSEmbed(selData,selNum,ref);
    runTime = toc;
end