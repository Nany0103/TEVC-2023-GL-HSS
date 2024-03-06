function [dataSelSet,runTime] = DLHSS3(data,dataSel,selNum,refSel)
    % d-a no projection
    tic;
    runTime = [0];
    M = size(data,2);
    dataSelSet = zeros(selNum,M,2);
    dataSelSet(:,:,1) = dataSel;
    gen = 2;
    count = 0;
    while true
        isC = false(1,selNum);
        for j=1:selNum
            hvC = stk_dominatedhv(dataSel,refSel);
            while true
                g = gradient(dataSel,j,M,refSel);
                % Calculate the angle of all candidate solutions along the
                % direction
                dist = pdist2(data,dataSel(j,:));
                [~,indd] = mink(dist,2000);
                cos_angle = sum((data(indd,:)-dataSel(j,:)).*g,2)./(sqrt(sum((data(indd,:)-dataSel(j,:)).^2,2)).*norm(g));
                [~,inda] = max(cos_angle);
                hvN = stk_dominatedhv([dataSel([1:j-1,j+1:end],:);data(indd(inda),:)],refSel);
                if hvN>hvC
                    count = 0;
                    isC(j) = true;
                    % disp([gen,j]);
                    dataSel(j,:) = data(indd(inda),:);
                    hvC = hvN;
                else
                    break;
                end
                count = count+1;
                if count>=selNum
                    break
                end
            end
        end
        dataSelSet(:,:,gen) = dataSel;
        runTime = [runTime,runTime(end)+toc];
        tic;
        gen = gen+1;
        if sum(isC)==0 || count>=selNum
            break
        end
    end
end