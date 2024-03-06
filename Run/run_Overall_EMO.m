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
            % run APL-HSS algorithm
            run_APLHSS_EMO(ref,M,data,selNum,pro{proInd},runInd,fmin,fmax);
            % run LGI-HSS algorithm
            run_LGIHSS_EMO(ref,M,data,selNum,pro{proInd},runInd,fmin,fmax);
            % run TGI-HSS algorithm
            run_TGIHSS_EMO(ref,M,data,selNum,pro{proInd},runInd,fmin,fmax);
            % run TGAHSS and GL-HSS algorithm: the proposed algorithm
            run_GLTGA_EMO(ref,M,data,selNum,pro{proInd},numSol2,runInd,fmin,fmax);
            % run DSS algorithm
            run_DSS_EMO(M,data,selNum,pro{proInd},runInd,fmin,fmax);
            % run CSS-MEA algorithm
            run_KME_EMO(M,data,selNum,pro{proInd},runInd,fmin,fmax);
            % run IDSS algorithm
            run_IDSS_EMO(M,data,selNum,pro{proInd},runInd,fmin,fmax);
            % run GAHSS algorithm
            run_GA_EMO(ref,M,data,selNum,pro{proInd},runInd,fmin,fmax);
        end
    end
end

function run_APLHSS_EMO(ref,M,data,selNum,proDispName,runInd,fmin,fmax)
    fileName = sprintf('./Result/Result New/APLHSS_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    if ~exist(fileName,'file')
        [dataSelSet,runTime] = APLHSS(data,selNum,1,1e-3,ref);
         % denormalize the selected subset
        for i=1:size(dataSelSet,3)
            dataSelSet(:,:,i) = dataSelSet(:,:,i).*(fmax-fmin)+fmin;
        end
        save(fileName,'dataSelSet','runTime');
    end
    fprintf('M%d_%s_APLHSS_runInd=%d\n',M,proDispName,runInd);
end

function run_LGIHSS_EMO(ref,M,data,selNum,proDispName,runInd,fmin,fmax)
    fileName = sprintf('./Result/Result New/LGIHSS_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    if ~exist(fileName,'file')
        [dataSel,runTime] = LGIHSS(data,selNum,ref,1000000);
        % denormalize the selected subset
        dataSel = dataSel.*(fmax-fmin)+fmin;
        save(fileName,'dataSel','runTime');
    end
    fprintf('M%d_%s_LGIHSS_runInd=%d\n',M,proDispName,runInd);
end

function run_TGIHSS_EMO(ref,M,data,selNum,proDispName,runInd,fmin,fmax)
    fileNameDSS = sprintf('./Result/Result New/TGIHSS_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    if ~exist(fileNameDSS,'file')
        [dataSel,runTime] = TGIHSS(data,selNum,ref,1);
        % denormalize the selected subset
        dataSel = dataSel.*(fmax-fmin)+fmin;
        save(fileNameDSS,'dataSel','runTime');
    end
    fprintf('M%d_%s_DSS_TGIHSS_runInd=%d\n',M,proDispName,runInd);
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

function run_DSS_EMO(M,data,selNum,proDispName,runInd,fmin,fmax)
    fileNameDSS = sprintf('./Result/Result New/DSS_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    if ~exist(fileNameDSS,'file')
        [dataSel,runTime] = selSolDSS_GI(data,selNum);
        % denormalize the selected subset
        dataSel = dataSel.*(fmax-fmin)+fmin;
        save(fileNameDSS,'dataSel','runTime');
    end
end

function run_KME_EMO(M,data,selNum,proDispName,runInd,fmin,fmax)
    fileName = sprintf('./Result/Result New/KME_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    if ~exist(fileName,'file')
        [dataSel,runTime] = K_means_S(data,selNum);
        % denormalize the selected subset
        dataSel = dataSel.*(fmax-fmin)+fmin;
        save(fileName,'dataSel','runTime');
    end
    fprintf('M%d_%s_KME_runInd=%d\n',M,proDispName,runInd);
end

function run_IDSS_EMO(M,data,selNum,proDispName,runInd,fmin,fmax)
    fileName = sprintf('./Result/Result New/IDSS_M%d_%s_selNum=%d_runInd=%d.mat',M,proDispName,selNum,runInd);
    if ~exist(fileName,'file')
        [dataSel,runTime] = selSolDSS_I(data,selNum);
        % denormalize the selected subset
        dataSel = dataSel.*(fmax-fmin)+fmin;
        save(fileName,'dataSel','runTime');
    end
    fprintf('M%d_%s_IDSS_GL-NAGO_runInd=%d\n',M,proDispName,runInd);
end

function run_GA_EMO(ref,M,data,selNum,proDispName,runInd,fmin,fmax)
    j = 100;
    fprintf('GAHSS%d_M%d_%s_selNum=%d_runInd=%d\n',j,M,proDispName,selNum,runInd);
    fileName = sprintf('./Result/Result New/GAHSS%d_M%d_%s_selNum=%d_runInd=%d.mat',j,M,proDispName,selNum,runInd);
    if ~exist(fileName,'file')
        [dataSel,runTime] =  GAHSS(data,selNum,j,1,ref);
        % denormalize the selected subset
        dataSel = dataSel.*(fmax-fmin)+fmin;
        save(fileName,'dataSel','runTime');
    end
end