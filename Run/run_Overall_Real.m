clc; clear;
obj = [6,9,9,10];
pro = {'RE61','RE91','DDMOP1','DDMOP4'};
DPlus = [3,7,11,13];
selNum = 100;
runNum = 21;
% parameters for TGAHSS
numSol2 = 5000;
for proInd = 1:4
    M = obj(proInd);
    H = getH(selNum,M);
    r = 1+1/H;
    D = DPlus(proInd);
    for runInd = 1:runNum
        data = load(sprintf('./Data/EMOA/NSGAIII_%s_M%d_D%d_%d',pro{proInd},M,D,runInd)).UEA;
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
    fileNameGA = sprintf('./Result/Result New/TGAHSS_M%d_%s_selNum=%d_numSol2=%d_runInd=%d.mat',M,proDispName,selNum,numSol2,runInd);
    if exist(fileNameGA,'file')
        dataSelTGA = load(fileNameGA).dataSel;
    else
        [dataSel,runTime] = TGAHSS(data,selNum,1,numSol2,1,ref);
        % denormalize the selected subset
        dataSel = dataSel.*(fmax-fmin)+fmin;
        save(fileNameGA,'dataSel','runTime');
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