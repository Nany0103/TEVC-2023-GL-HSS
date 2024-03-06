clc; clear;
obj = [8,10,6,9,9,10];
pro = {'WFG3','WFG3','RE61','RE91','DDMOP1','DDMOP4'};
DPlus = [17,19,3,7,11,13];
selNum = 100;
numSol2 = 5000;
runNum = 21;
alpha = 1000;
for proInd = 4:6
    M = obj(proInd);
    D = DPlus(proInd);
    data = [];
    for runInd = 1:runNum
        dataSpec = load(sprintf('./Data/EMOA/NSGAIII_%s_M%d_D%d_%d',pro{proInd},M,D,runInd)).UEA;
        data = [data;dataSpec];
    end
    % use objective partition method to obtain the non-dominated
    % solutions
    method = getMethods(3,1,M);
    index = IterativePartition(data, 1000, method);
    index = index{end};
    data = data(index,:);
    fileName = sprintf('./Data/Approximated PF/ApproximatedPF_%s_M%d.mat',pro{proInd},M);
    save(fileName,'data');
end