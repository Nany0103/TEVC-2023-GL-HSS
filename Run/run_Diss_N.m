clc; clear;
obj = [8,10];
pro = {'Linear','Concave','Convex','I-Linear','I-Concave','I-Convex'};
selNum = 100;
numSol2 = 5000;
runNum = 21;
NSet = [1000000];
for objInd = 2
    M = obj(objInd);
    H = getH(selNum,M);
    r = 1+1/H;
    for NSetInd = 1:length(NSet)
        N = NSet(NSetInd);
        for proInd = 6
            parfor runInd = 1:runNum
                data = load(sprintf('./Data/PF_New/%s_N%d_M%d_%d',pro{proInd},N,M,runInd)).data;
                ref = ones(1,M)*r;
                warning('off');
                if N<=10000
                    run_GL_NAGO_Rand(ref,M,data,N,selNum,pro{proInd},runInd);
                else
                    run_GL_NAGO_TGA(ref,M,data,N,selNum,pro{proInd},numSol2,runInd);
                end
%                 run_LGI(ref,M,data,N,selNum,pro{proInd},runInd);
            end
        end
    end
end

function run_LGI(ref,M,data,N,selNum,proDispName,runInd)
    fileName = sprintf('./Result/Result New/LGIHSS_M%d_%s_N=%d_selNum=%d_runInd=%d.mat', ...
        M,proDispName,N,selNum,runInd);
    fprintf('LGIHSS_M%d_%s_N=%d_selNum=%d_runInd=%d\n',M,proDispName,N,selNum,runInd);
    if ~exist(fileName,'file') % 文件不存在
        [dataSel,runTime] = LGIHSS(data,selNum,ref);
        save(fileName,'dataSel','runTime');
    end
end

function run_GL_NAGO_Rand(ref,M,data,N,selNum,proDispName,runInd)
    fileNameRand = sprintf('./Result/Result New/Rand_M%d_%s_N=%d_selNum=%d_runInd=%d.mat',M,proDispName,N,selNum,runInd);
    if exist(fileNameRand,'file')
        dataSelGA = load(fileNameRand).dataSel;
    else
        rng(1);
        tic;
        dataSel = data(randperm(size(data,1),selNum),:);
        runTime = toc;
        save(fileNameRand,'dataSel','runTime');
        dataSelGA = dataSel;
    end
    fprintf('M%d_%s_%d_GL-NAGO-Rand_runInd=%d\n',M,proDispName,N,runInd);
    fileName = sprintf('./Result/Result New/GLHSSNAGO_Rand_M%d_%s_N=%d_selNum=%d_runInd=%d.mat',M,proDispName,N,selNum,runInd);
    [dataSelSet,runTime] = GLHSS_NAGO(data,dataSelGA,selNum,ref);
    save(fileName,'dataSelSet','runTime');
end


function run_GL_NAGO_TGA(ref,M,data,N,selNum,proDispName,numSol2,runInd)
    fileNameGA = sprintf('./Result/Result New/TGAHSS_M%d_%s_N=%d_selNum=%d_numSol2=%d_runInd=%d.mat',M,proDispName,N,selNum,numSol2,runInd);
    if exist(fileNameGA,'file')
        dataSelGA = load(fileNameGA).dataSel;
    else
        [dataSel,runTime] = TGAHSS(data,selNum,1,numSol2,1,ref);
        save(fileNameGA,'dataSel','runTime');
        dataSelGA = dataSel;
    end
    fprintf('M%d_%s_%d_GL-NAGO_numSol2=%d_runInd=%d\n',M,proDispName,N,numSol2,runInd);
    fileName = sprintf('./Result/Result New/GLHSSNAGO_TGA_M%d_%s_N=%d_selNum=%d_numSol2=%d_runInd=%d.mat',M,proDispName,N,selNum,numSol2,runInd);
    [dataSelSet,runTime] = GLHSS_NAGO(data,dataSelGA,selNum,ref);
    save(fileName,'dataSelSet','runTime');
end

function H = getH(N,M)
    H = 1;
    while nchoosek(H+M,M-1) <= N
        H = H + 1;
    end
end