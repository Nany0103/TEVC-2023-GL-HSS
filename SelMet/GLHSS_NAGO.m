function [dataSelSet,runTime] = GLHSS_NAGO(data,dataSel,selNum,refSel)
    % d-a no projection
    tic;
    A = AssociationIni(data,dataSel);
    runTime = [0];
    M = size(data,2);
    dataSelSet = zeros(selNum,M,2);
    dataSelSet(:,:,1) = dataSel;
    gen = 2;
    count = 0;
    while true
        isC = false(1,selNum);
        for j=1:selNum
            hvC = HVCE(dataSel(j,:),dataSel([1:j-1,j+1:end],:),refSel);
            while true
                g = gradient(dataSel,j,M,refSel);
                % Calculate the angle of all candidate solutions along the
                % direction
                indd = find(A(:,2)==j);
                cos_angle = sum((data(indd,:)-dataSel(j,:)).*g,2)./(sqrt(sum((data(indd,:)-dataSel(j,:)).^2,2)).*norm(g));
                [~,inda] = max(cos_angle);
                hvN = HVCE(data(indd(inda),:),dataSel([1:j-1,j+1:end],:),refSel);
                if hvN>hvC
                    count = 0;
                    isC(j) = true;
                    A = AssociationUpd(data,dataSel,j,data(indd(inda),:),A);
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
function hvc = HVCE(s,data,r)
    dataP = max(data,s);
    hvc = prod(r-s)-stk_dominatedhv(dataP,r);
end
function A = AssociationIni(data,dataSel)
    [n,~] = size(data);
    A = zeros(n,2);
    dist = pdist2(data,dataSel);
    % A(:,1) min distance, A(:,2) min index
    [A(:,1),A(:,2)] = min(dist,[],2);
end
function A = AssociationUpd(data,dataSel,i,newS,A)
    [n,~] = size(data);
    indRel = find(A(:,2)==i);
    indOth = find(A(:,2)~=i);
    % distance of other candidate solutions to new
    distToNew = pdist2(data(indOth,:),newS);
    % find the index of other candidate solutions whose association changes
    indChangeOth = indOth(distToNew<A(indOth,1));
    % include the special other index to related index
    indRel = [indRel;indChangeOth];
    % calculate the distance of related candiate solutions to all the
    % sleected solutions
    dataSel(i,:) = newS;
    distToAll = pdist2(data(indRel,:),dataSel);
    [A(indRel,1),A(indRel,2)] = min(distToAll,[],2);
end