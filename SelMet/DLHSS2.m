function [dataSelSet,runTime] = DLHSS2(data,dataSel,selNum,refSel,proInd)
    % d-a
    tic;
    runTime = [0];
    M = size(data,2);
    dataSelSet = zeros(selNum,M,2);
    dataSelSet(:,:,1) = dataSel;
    gen = 2;
    while true
        isC = false(1,selNum);
        for j=1:selNum
            hvC = stk_dominatedhv(dataSel,refSel);
            while true
                g = gradient(dataSel,j,M,refSel);
                % project the gradient to the PF plane
                if proInd==1 || proInd==4
                    n = (-1/sqrt(M))*ones(1,M);
                elseif proInd==2 || proInd==5
                    n = -1*dataSel(j,:);
                    n = (1/norm(n))*n;
                elseif proInd==3 || proInd==6
                    n = dataSel(j,:)-1;
                    n = (1/norm(n))*n;
                end
                cos_angle = sum(g.*n)./(norm(g)*norm(n));
                gpn = norm(g)*cos_angle;
                gp = n*gpn;
                gf = g-gp;
                % Calculate the angle of all candidate solutions along the
                % direction
                dist = pdist2(data,dataSel(j,:));
                [~,indd] = mink(dist,20);
                cos_angle = sum((data(indd,:)-dataSel(j,:)).*gf,2)./(sqrt(sum((data(indd,:)-dataSel(j,:)).^2,2)).*norm(gf));
                [~,inda] = max(cos_angle);
                hvN = stk_dominatedhv([dataSel([1:j-1,j+1:end],:);data(indd(inda),:)],refSel);
                if hvN>hvC
                    isC(j) = true;
                    % disp([gen,j]);
                    dataSel(j,:) = data(indd(inda),:);
                    hvC = hvN;
                else
                    break;
                end
            end
        end
        dataSelSet(:,:,gen) = dataSel;
        runTime = [runTime,runTime(end)+toc];
        tic;
        gen = gen+1;
        if sum(isC)==0
            break
        end
    end
end