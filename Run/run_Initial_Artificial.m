clc; clear;
obj = [8,10];
pro = {'Linear','Concave','Convex','I-Linear','I-Concave','I-Convex'};
selNum = 100;
runNum = 21;
% parameters for TGAHSS
numSol2 = 5000;
for objInd = 1:2
    M = obj(objInd);
    H = getH(selNum,M);
    r = 1+1/H;
    for proInd = 1:6
        for runInd = 1:runNum
            data = load(sprintf('./Data/PF_New/%s_N1000000_M%d_%d',pro{proInd},M,runInd)).data;
            ref = ones(1,M)*r;
            warning('off');
            % run TGAHSS and GL-HSS algorithm: the proposed algorithm
            run_GLTGA_Artificial(ref,M,data,selNum,pro{proInd},numSol2,runInd);
            % run DSS-GLHSS algorithm
            run_GLDSS_Artificial(ref,M,data,selNum,pro{proInd},runInd);
            % run IDSS-GLHSS algorithm
            run_GLIDSS_Artificial(ref,M,data,selNum,pro{proInd},runInd);
            % run Rand-GLHSS algorithm
            run_GLRand_Artificial(ref,M,data,selNum,pro{proInd},runInd);
        end
    end
end

function run_GLTGA_Artificial(ref,M,data,selNum,proDispName,numSol2,runInd)
    fileNameGA = sprintf('./Result/Result New/TGAHSS_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat',M,proDispName,selNum,numSol2,runInd);
    if exist(fileNameGA,'file')
        dataSelGA = load(fileNameGA).dataSel;
    else
        [dataSel,runTime] = TGAHSS(data,selNum,1,numSol2,1,ref);
        save(fileNameGA,'dataSel','runTime');
        dataSelGA = dataSel;
    end
    fprintf('M%d_%s_GL-NAGO_numSol2=%d_runInd=%d\n',M,proDispName,numSol2,runInd);
    fileName = sprintf('./Result/Result New/GLHSSNAGO_TGA_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat',M,proDispName,selNum,numSol2,runInd);
    [dataSelSet,runTime] = GLHSS_NAGO(data,dataSelGA,selNum,ref);
    save(fileName,'dataSelSet','runTime');
end

function run_GLDSS_Artificial(ref,M,data,selNum,proDispName,runInd)
    fileNameDSS = sprintf('./Result/Result New/DSS_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    if exist(fileNameDSS,'file')
        dataSelDSS = load(fileNameDSS).dataSel;
    else
        [dataSel,runTime] = selSolDSS_GI(data,selNum);
        save(fileNameDSS,'dataSel','runTime');
        dataSelDSS = dataSel;
    end
    fprintf('M%d_%s_DSS_GL-NAGO_runInd=%d\n',M,proDispName,runInd);
    fileName = sprintf('./Result/Result New/GLHSSNAGO_DSS_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    [dataSelSet,runTime] = GLHSS_NAGO(data,dataSelDSS,selNum,ref);
    save(fileName,'dataSelSet','runTime');
end

function run_GLIDSS_Artificial(ref,M,data,selNum,proDispName,runInd)
    fileNameDSS = sprintf('./Result/Result New/IDSS_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    if exist(fileNameDSS,'file')
        dataSelIDSS = load(fileNameDSS).dataSel;
    else
        [dataSel,runTime] = selSolDSS_I(data,selNum);
        save(fileNameDSS,'dataSel','runTime');
        dataSelIDSS = dataSel;
    end
    fprintf('M%d_%s_IDSS_GL-NAGO_runInd=%d\n',M,proDispName,runInd);
    fileName = sprintf('./Result/Result New/GLHSSNAGO_IDSS_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    [dataSelSet,runTime] = GLHSS_NAGO(data,dataSelIDSS,selNum,ref);
    save(fileName,'dataSelSet','runTime');
end

function run_GLRand_Artificial(ref,M,data,selNum,proDispName,runInd)
    fileNameRand = sprintf('./Result/Result New/Rand_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    if exist(fileNameRand,'file')
        dataSelRand = load(fileNameRand).dataSel;
    else
        rng(1);
        tic;
        dataSel = data(randperm(size(data,1),selNum),:);
        runTime = toc;
        save(fileNameRand,'dataSel','runTime');
        dataSelRand = dataSel;
    end
    fprintf('M%d_%s_Rand_GL-NAGO_runInd=%d\n',M,proDispName,runInd);
    fileName = sprintf('./Result/Result New/GLHSSNAGO_Rand_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    [dataSelSet,runTime] = GLHSS_NAGO(data,dataSelRand,selNum,ref);
    save(fileName,'dataSelSet','runTime');
end