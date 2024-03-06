function [Subset,runTime,numDist,count] = TGAHSS(objVal,selNum,num_vec1,num_sol2,seed,ref)
    % this version is from TGAHSS8, but in this version, for all number of
    % vec1, the first vector are specific designed
    tic;
    rng(seed);
    [solNum,M] = size(objVal);
    num_vec2 = 100;
    if num_vec1>0
        W1 = zeros(num_vec1,M);
        W1(1,:) = sqrt(1/M)*ones(1,M);
        W1(2:end,:) = UniformVector(num_vec1-1,M);
    end
    [W2,num_vec2] = UniformVector(num_vec2,M);
    % tensor is a matrix each row refers to each candidate solution and
    % each colum refers to the minimum distance of each weight vector. The left
    % part of the tensor is for the first stage and all the tensor is for
    % the second stage.
    if num_vec1>0
        tensor1 = inf(solNum,num_vec1);
    end
    tensor2 = inf(solNum,num_vec2);
    % count records the number of selected solutions which have been considered for
    % the distance calculation in the right part of the tensor
    count = zeros(1,solNum);
    r = ref;
    % When calculating the minimum distance of WV, we also need to consider
    % the reference point. This part is to calculate the distance of all DVs
    % from all solutions to ref
    numDist = 0;
    if num_vec1>0
        if num_vec1<10
            for i=1:num_vec1
                temp1 = min(abs(objVal-r)./W1(i,:),[],2);        
                tensor1(:,i) = temp1;     
            end
        else
            for i=1:solNum
                temp1 = min(abs(objVal(i,:)-r)./W1,[],2)';        
                tensor1(i,:) = temp1;     
            end
        end
    end
    if num_sol2>1
        for i=1:solNum
            temp2 = min(abs(objVal(i,:)-r)./W2,[],2)';        
            tensor2(i,:) = temp2;     
        end
        numDist = (num_vec1+num_vec2)*solNum;
    end
    % mintensor records the minimum distance of each WVs. tensor is only a
    % temptation matrix for calculating the WV for only the current
    % selected solution
    if num_vec1>0
        mintensor1 = tensor1;
    end
    mintensor2 = tensor2;
    selVal = zeros(selNum,M);
    selInd = zeros(1,selNum);
    canInd = true(1,solNum);
    for num = 1:selNum
        % update the minimum distance of the left part by comparing each
        % elements in the left part of mintensor with tensor in the left
        % part. This tensor records the distance of the left part of WV
        % considering the last selected solution
        if num_vec1>0
            mintensor1 = min(mintensor1, tensor1);
            % r2hvc1 is the approximated HVC of each candidate solution using
            % the left part of mintensor in the first stage
            r2hvc1 = sum(mintensor1,2);
            % candidateInd is the index of the selected solutions from
            % candidate set according to the roughly approximatedly HVC using
            % the left part of mintensor
            [~,candidateInd] = maxk(r2hvc1,num_sol2);
        else
            unselectedInd = find(canInd);
            candidateInd = unselectedInd(randperm(solNum-num+1,num_sol2));
        end
        % update the right part of mintensor of each solution of the selected solutions
        % from candidate set according the roughly approximatedly HVC
        if num_sol2>1
            for i = 1:num_sol2
                % s refers to each selected candidate solution
                s = objVal(candidateInd(i),:);
                % since we have calculate the distance of all WVs for the ref
                % (count=0), when num=1, we don't have to do any update for
                % mintensor for the right part of mintensor
                if num>1
                    % temp2 records all the distance of the right part WV for
                    % consideration of the selected solutions which are not
                    % considered for the right part WV
                    temp3 = zeros(num-1-count(candidateInd(i)),num_vec2);
                    % since the currently num-1 solutions has been selected, we
                    % caculated the distance of right part WV for all selected
                    % solution which are not considered for the right part WV
                    for j = count(candidateInd(i))+1:num-1
                        % sn is the unconsidered selected solution in the
                        % subset
                        sn = selVal(j,:);
                        % calculate the distance of the right WV for the
                        % selected candidate solution s and the selected
                        % unconsidered solution sn
                        temp3(j-count(candidateInd(i)),:) = max((sn-s)./W2,[],2)';
                        numDist = numDist+num_vec2;
                    end
                    % update the distance of the right part WV in the mintensor
                    % of the selected candidate solution s
                    mintensor2(candidateInd(i),:) = ...
                        min(mintensor2(candidateInd(i),:),min(temp3,[],1));
                end
            end
            % update count
            count(candidateInd) = num-1;
            % calculate the approximated HVC for the selected candidate
            % solutions using all WVs
            r2hvc2 = sum(mintensor2(candidateInd,:),2);
            % maxInd is the index of the selected candidate solution with the
            % largerst approximated HVC
            [~,maxInd] = max(r2hvc2);
            % bestindex is the index of the selected candiate solution in
            % the candidate set
            bestindex = candidateInd(maxInd);
        else
            bestindex = candidateInd;
        end
        if num_vec1>0
            % calculate the distance of each wv in the left part for the
            % selected solution in this iteration
            if num_vec1<10
                for i=1:num_vec1
                    temp1 =  max((objVal(bestindex,:)-objVal)./W1(i,:),[],2);
                    tensor1(:,i) = temp1;      
                end
            else
                for i=1:solNum
                    temp1 =  max((objVal(bestindex,:)-objVal(i,:))./W1,[],2)';
                    tensor1(i,:) = temp1;  
                end
            end
            numDist = numDist+num_vec1*solNum;
        end
        % update selVal and selInd
        selVal(num,:) = objVal(bestindex,:);
        selInd(1,num) = bestindex;
        canInd(1,bestindex) = false;
    end
    Subset = selVal;
    runTime = toc;
end
