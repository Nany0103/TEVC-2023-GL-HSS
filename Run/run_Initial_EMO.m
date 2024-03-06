clc; clear;
obj = [8,10];
pro = {'DTLZ1','DTLZ2','MinusDTLZ1','MinusDTLZ2','WFG3','DTLZ7'};
DPlus = [4,9,4,9,9,19];
selNum = 100;
runNum = 21;
% parameters for TGAHSS
numSol2 = 5000;
for objInd = 1:2
    M = obj(objInd);
    H = getH(selNum,M);
    r = 1+1/H;
    for proInd = 1:6
        D = M+DPlus(proInd);
        for runInd = 1:runNum
            data = load(sprintf('./Data/PF_New/%s_N1000000_M%d_%d',pro{proInd},M,runInd)).data;
            fmin = min(data); fmax = max(data);
            % normalization
            data = (data-fmin)./(fmax-fmin);
            ref = ones(1,M)*r;
            warning('off');
            % run TGAHSS and GL-HSS algorithm: the proposed algorithm
            run_GLTGA_EMO(ref,M,data,selNum,pro{proInd},numSol2,runInd,fmin,fmax);
            % run DSS-GLHSS algorithm
            run_GLDSS_EMO(ref,M,data,selNum,pro{proInd},runInd,fmin,fmax);
            % run IDSS-GLHSS algorithm
            run_GLIDSS_EMO(ref,M,data,selNum,pro{proInd},runInd,fmin,fmax);
            % run Rand-GLHSS algorithm
            run_GLRand_EMO(ref,M,data,selNum,pro{proInd},runInd,fmin,fmax);
        end
    end
end

function run_GLTGA_EMO(ref,M,data,selNum,proDispName,numSol2,runInd,fmin,fmax)
    fileNameTGA = sprintf('./Result/Result New/TGAHSS_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat',M,proDispName,selNum,numSol2,runInd);
    if exist(fileNameTGA,'file')
        dataSelTGA = load(fileNameTGA).dataSel;
    else
        [dataSel,runTime] = TGAHSS(data,selNum,1,numSol2,1,ref);
        % denormalize the selected subset
        dataSel = dataSel.*(fmax-fmin)+fmin;
        save(fileNameTGA,'dataSel','runTime');
        dataSelTGA = dataSel;
    end
    % normalize dataSel
    dataSelTGA = (dataSelTGA-fmin)./(fmax-fmin);
    fprintf('M%d_%s_GL-NAGO_numSol2=%d_runInd=%d\n',M,proDispName,numSol2,runInd);
    fileName = sprintf('./Result/Result New/GLHSSNAGO_TGA_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat',M,proDispName,selNum,numSol2,runInd);
    [dataSelSet,runTime] = GLHSS_NAGO(data,dataSelTGA,selNum,ref);
    % denormalize dataSelSet
    for i=1:size(dataSelSet,3)
        dataSelSet(:,:,i) = dataSelSet(:,:,i).*(fmax-fmin)+fmin;
    end
    save(fileName,'dataSelSet','runTime');
end

function run_GLDSS_EMO(ref,M,data,selNum,proDispName,runInd,fmin,fmax)
    fileNameDSS = sprintf('./Result/Result New/DSS_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    if exist(fileNameDSS,'file')
        dataSelDSS = load(fileNameDSS).dataSel;
    else
        [dataSel,runTime] = selSolDSS_GI(data,selNum);
        % denormalize the selected subset
        dataSel = dataSel.*(fmax-fmin)+fmin;
        save(fileNameDSS,'dataSel','runTime');
        dataSelDSS = dataSel;
    end
    % normalize dataSel
    dataSelDSS = (dataSelDSS-fmin)./(fmax-fmin);
    fprintf('M%d_%s_DSS_GL-NAGO_runInd=%d\n',M,proDispName,runInd);
    fileName = sprintf('./Result/Result New/GLHSSNAGO_DSS_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    [dataSelSet,runTime] = GLHSS_NAGO(data,dataSelDSS,selNum,ref);
    % denormalize dataSelSet
    for i=1:size(dataSelSet,3)
        dataSelSet(:,:,i) = dataSelSet(:,:,i).*(fmax-fmin)+fmin;
    end
    save(fileName,'dataSelSet','runTime');
end

function run_GLIDSS_EMO(ref,M,data,selNum,proDispName,runInd,fmin,fmax)
    fileNameDSS = sprintf('./Result/Result New/IDSS_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    if exist(fileNameDSS,'file')
        dataSelIDSS = load(fileNameDSS).dataSel;
    else
        [dataSel,runTime] = selSolDSS_I(data,selNum);
        % denormalize the selected subset
        dataSel = dataSel.*(fmax-fmin)+fmin;
        save(fileNameDSS,'dataSel','runTime');
        dataSelIDSS = dataSel;
    end
    % normalize dataSel
    dataSelIDSS = (dataSelIDSS-fmin)./(fmax-fmin);
    fprintf('M%d_%s_IDSS_GL-NAGO_runInd=%d\n',M,proDispName,runInd);
    fileName = sprintf('./Result/Result New/GLHSSNAGO_IDSS_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    [dataSelSet,runTime] = GLHSS_NAGO(data,dataSelIDSS,selNum,ref);
    % denormalize dataSelSet
    for i=1:size(dataSelSet,3)
        dataSelSet(:,:,i) = dataSelSet(:,:,i).*(fmax-fmin)+fmin;
    end
    save(fileName,'dataSelSet','runTime');
end

function run_GLRand_EMO(ref,M,data,selNum,proDispName,runInd,fmin,fmax)
    fileNameRand = sprintf('./Result/Result New/Rand_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    if exist(fileNameRand,'file')
        dataSelRand = load(fileNameRand).dataSel;
    else
        rng(1);
        tic;
        dataSel = data(randperm(size(data,1),selNum),:);
        runTime = toc;
        % denormalize the selected subset
        dataSel = dataSel.*(fmax-fmin)+fmin;
        save(fileNameRand,'dataSel','runTime');
        dataSelRand = dataSel;
    end
    % normalize dataSel
    dataSelRand = (dataSelRand-fmin)./(fmax-fmin);
    fprintf('M%d_%s_Rand_GL-NAGO_runInd=%d\n',M,proDispName,runInd);
    fileName = sprintf('./Result/Result New/GLHSSNAGO_Rand_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    [dataSelSet,runTime] = GLHSS_NAGO(data,dataSelRand,selNum,ref);
    % denormalize dataSelSet
    for i=1:size(dataSelSet,3)
        dataSelSet(:,:,i) = dataSelSet(:,:,i).*(fmax-fmin)+fmin;
    end
    save(fileName,'dataSelSet','runTime');
end