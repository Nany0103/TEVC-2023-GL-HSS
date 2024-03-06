function [selData,runTime] = TGIHSSIni(data,selNum,ref,seed)
    tic;
    alpha = 20;
    [~,M] = size(data);
    ind = zeros(1,alpha*selNum);
    for i = 1:alpha
        [~,ind((i-1)*selNum+1:i*selNum)] = GAHSS_Small(data,selNum,ref,(seed-1)*alpha+i,1);
    end
    selData = data(ind,:);
    selData = LGIHSSEmbed(selData,selNum,ref);
    runTime = toc;
end