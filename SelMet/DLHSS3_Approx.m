function [dataSelSet,runTime] = DLHSS3_Approx(data,dataSel,selNum,seed,refSel)
    % d-a no projection
    rng(seed);
    tic;
    runTime = [0];
    M = size(data,2);
    dataSelSet = zeros(selNum,M,2);
    dataSelSet(:,:,1) = dataSel;
    gen = 2;
    count = 0;
    W = UniformVector(10000,M);
    while true
        isC = false(1,selNum);
        for j=1:selNum
            hvC = R2HV(dataSel,W,refSel);
            while true
                g = gradient_Approx(dataSel,j,M,seed,refSel);
                % Calculate the angle of all candidate solutions along the
                % direction
                dist = pdist2(data,dataSel(j,:));
                [~,indd] = mink(dist,2000);
                cos_angle = sum((data(indd,:)-dataSel(j,:)).*g,2)./(sqrt(sum((data(indd,:)-dataSel(j,:)).^2,2)).*norm(g));
                [~,inda] = max(cos_angle);
                hvN = R2HV([dataSel([1:j-1,j+1:end],:);data(indd(inda),:)],W,refSel);
                if hvN>hvC
                    count = 0;
                    isC(j) = true;
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