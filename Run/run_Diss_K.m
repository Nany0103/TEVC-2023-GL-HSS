clc; clear;
obj = 10;
pro = {'Linear','Concave','Convex','I-Linear','I-Concave','I-Convex'};
selNumSet = [50,200];
numSol2 = 5000;
runNum = 21;
N = 1000000;
for objInd = 1
    M = obj(objInd);
    for kInd = 1:5
        selNum = selNumSet(kInd);
        H = getH(selNum,M);
        r = 1+1/H;
        for proInd = 1:6
            parfor runInd = 1:runNum
                data = load(sprintf('./Data/PF_New/%s_N%d_M%d_%d',pro{proInd},N,M,runInd)).data;
                ref = ones(1,M)*r;
                warning('off');
                run_GL_NAGO_TGA(ref,M,data,selNum,pro{proInd},numSol2,runInd);
                run_LGI(ref,M,data,selNum,pro{proInd},runInd);
            end
        end
    end
end

function run_LGI(ref,M,data,selNum,proDispName,runInd)
    fileName = sprintf('./Result/Result New/LGIHSS_M%d_%s_selNum=%d_runInd=%d.mat', ...
        M,proDispName,selNum,runInd);
    fprintf('LGIHSS_M%d_%s_selNum=%d_runInd=%d\n',M,proDispName,selNum,runInd);
    if ~exist(fileName,'file')
        [dataSel,runTime] = LGIHSS(data,selNum,ref);
        save(fileName,'dataSel','runTime');
    end
end

function run_GL_NAGO_TGA(ref,M,data,selNum,proDispName,numSol2,runInd)
    fileNameGA = sprintf('./Result/Result New/TGAHSS_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat',M,proDispName,selNum,numSol2,runInd);
    if ~exist(fileNameGA,'file')
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

function H = getH(N,M)
    H = 1;
    while nchoosek(H+M,M-1) <= N
        H = H + 1;
    end
end